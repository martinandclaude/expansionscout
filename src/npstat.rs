//! numpy-compatible statistics, and the arithmetic they rest on.
//!
//! numpy's reductions are not a plain left-to-right sum: a contiguous float64
//! `sum` is pairwise with eight accumulators and a 128-element block. That
//! order changes the last bits of a result, so it is reproduced exactly here
//! rather than approximated. `median`, `quantile` and `var` follow numpy 2.x's
//! own formulas for the same reason.
//!
//! Transcendental functions go through the `libm` crate, not the platform's
//! maths library. numpy uses its own SIMD `exp` and `log` on CPUs with
//! AVX-512F and glibc's elsewhere, and the two disagree in the last bit for a
//! few percent of arguments, so the Python implementation's mixture fit is
//! CPU-dependent at that level. `libm` gives the same bits on every platform.

/// numpy's `pairwise_sum` for float64, as `np.add.reduce` applies it to a
/// contiguous 1-D array.
pub fn pairwise_sum(a: &[f64]) -> f64 {
    0.0 + pw(a)
}

fn pw(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut res = 0.0;
        for &v in a {
            res += v;
        }
        res
    } else if n <= 128 {
        let mut r = [a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]];
        let mut i = 8;
        while i < n - (n % 8) {
            for j in 0..8 {
                r[j] += a[i + j];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        while i < n {
            res += a[i];
            i += 1;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        pw(&a[..n2]) + pw(&a[n2..])
    }
}

/// `float(np.mean(a))`.
pub fn mean(a: &[f64]) -> f64 {
    pairwise_sum(a) / a.len() as f64
}

/// `float(np.var(a))` (ddof = 0).
pub fn var(a: &[f64]) -> f64 {
    let n = a.len() as f64;
    let m = pairwise_sum(a) / n;
    let sq: Vec<f64> = a
        .iter()
        .map(|&v| {
            let d = v - m;
            d * d
        })
        .collect();
    pairwise_sum(&sq) / n
}

fn sorted(a: &[f64]) -> Vec<f64> {
    let mut v = a.to_vec();
    v.sort_by(|x, y| x.partial_cmp(y).expect("no NaN in counts"));
    v
}

/// `float(np.median(a))`: the middle value, or the mean of the two middle
/// values, which numpy computes as their sum over two.
pub fn median(a: &[f64]) -> f64 {
    let s = sorted(a);
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        (0.0 + s[n / 2 - 1] + s[n / 2]) / 2.0
    }
}

/// numpy's `_lerp`, which interpolates from whichever end is nearer.
fn lerp(a: f64, b: f64, t: f64) -> f64 {
    let diff = b - a;
    if t >= 0.5 {
        b - diff * (1.0 - t)
    } else {
        a + diff * t
    }
}

/// `np.quantile(a, q)` with the default `method="linear"`, on sorted data.
fn quantile_sorted(s: &[f64], q: f64) -> f64 {
    let n = s.len();
    let vi = (n as f64 - 1.0) * q;
    let (prev, next);
    let prev_f;
    if vi >= n as f64 - 1.0 {
        prev = n - 1;
        next = n - 1;
        prev_f = -1.0;
    } else if vi < 0.0 {
        prev = 0;
        next = 0;
        prev_f = 0.0;
    } else {
        prev_f = vi.floor();
        prev = prev_f as usize;
        next = prev + 1;
    }
    let gamma = vi - prev_f;
    lerp(s[prev], s[next], gamma)
}

/// `np.quantile(a, qs)`.
pub fn quantiles(a: &[f64], qs: &[f64]) -> Vec<f64> {
    let s = sorted(a);
    qs.iter().map(|&q| quantile_sorted(&s, q)).collect()
}

/// `float(np.percentile(a, p))`.
pub fn percentile(a: &[f64], p: f64) -> f64 {
    quantile_sorted(&sorted(a), p / 100.0)
}

/// Index of the first minimum, as `np.argmin` returns it.
pub fn argmin(a: &[f64]) -> usize {
    let mut best = 0;
    for (i, &v) in a.iter().enumerate() {
        if v < a[best] {
            best = i;
        }
    }
    best
}

/// Python's `round()` for a float: half to even.
pub fn py_round(v: f64) -> f64 {
    v.round_ties_even()
}

/// `math.fsum`: the correctly rounded sum, so independent of order and of the
/// Python version. Ported from CPython's `math_fsum` (Shewchuk's algorithm).
pub fn fsum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut partials: Vec<f64> = Vec::new();
    for mut x in values {
        let mut i = 0;
        for j in 0..partials.len() {
            let mut y = partials[j];
            if x.abs() < y.abs() {
                std::mem::swap(&mut x, &mut y);
            }
            let hi = x + y;
            let lo = y - (hi - x);
            if lo != 0.0 {
                partials[i] = lo;
                i += 1;
            }
            x = hi;
        }
        partials.truncate(i);
        if x != 0.0 {
            partials.push(x);
        }
    }
    let mut n = partials.len();
    let mut hi = 0.0;
    if n > 0 {
        n -= 1;
        hi = partials[n];
        let mut lo = 0.0;
        while n > 0 {
            let x = hi;
            n -= 1;
            let y = partials[n];
            hi = x + y;
            let yr = hi - x;
            lo = y - yr;
            if lo != 0.0 {
                break;
            }
        }
        // Round half-even correctly when the remaining partials would push
        // the result across a tie.
        if n > 0 && ((lo < 0.0 && partials[n - 1] < 0.0) || (lo > 0.0 && partials[n - 1] > 0.0)) {
            let y = lo * 2.0;
            let x = hi + y;
            let yr = x - hi;
            if y == yr {
                hi = x;
            }
        }
    }
    hi
}

pub fn exp(x: f64) -> f64 {
    libm::exp(x)
}

pub fn ln(x: f64) -> f64 {
    libm::log(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantile_and_median() {
        let a = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0];
        assert_eq!(median(&a), 3.5);
        assert_eq!(quantiles(&a, &[0.25, 0.75]), vec![1.75, 5.25]);
        assert_eq!(percentile(&[17.0], 5.0), 17.0);
    }

    #[test]
    fn fsum_is_exact() {
        assert_eq!(fsum([0.1; 10]), 1.0);
        assert_eq!(fsum([1e100, 1.0, -1e100, 1e-100]), 1.0);
    }
}
