//! Python-compatible text formatting.
//!
//! The outputs of this tool are text, and the Python implementation is the
//! reference for what that text is. Rust's `{:.N}` already agrees with
//! Python's `:.Nf` (both round the exact binary value half-to-even), but
//! `repr`, `:g`, `json.dumps`, `html.escape` and `textwrap.wrap` have no Rust
//! equivalent with the same output, so they are written out here against
//! CPython's own rules.

use std::fmt::Write;

/// `f"{v:.{nd}f}"`.
pub fn fixed(v: f64, nd: usize) -> String {
    format!("{:.*}", nd, v)
}

/// Shortest round-trip digits and the decimal exponent, from Rust's `{:e}`,
/// which uses the same shortest-digits rule as CPython's `repr`.
fn shortest(v: f64) -> (bool, String, i32) {
    let s = format!("{:e}", v);
    let (mant, exp) = s.split_once('e').expect("{:e} always has an exponent");
    let neg = mant.starts_with('-');
    let digits: String = mant.chars().filter(|c| c.is_ascii_digit()).collect();
    (neg, digits, exp.parse().unwrap())
}

fn exp_suffix(e: i32) -> String {
    if e < 0 {
        format!("e-{:02}", -e)
    } else {
        format!("e+{:02}", e)
    }
}

/// `repr(v)` for a float, as CPython writes it (`float_repr_style == 'short'`).
pub fn repr(v: f64) -> String {
    if v.is_nan() {
        return "nan".into();
    }
    if v.is_infinite() {
        return if v > 0.0 { "inf".into() } else { "-inf".into() };
    }
    if v == 0.0 {
        return if v.is_sign_negative() {
            "-0.0".into()
        } else {
            "0.0".into()
        };
    }
    let (neg, digits, e) = shortest(v);
    let decpt = e + 1;
    let mut out = String::new();
    if neg {
        out.push('-');
    }
    if decpt <= -4 || decpt > 16 {
        out.push_str(&digits[..1]);
        if digits.len() > 1 {
            out.push('.');
            out.push_str(&digits[1..]);
        }
        out.push_str(&exp_suffix(decpt - 1));
    } else if decpt <= 0 {
        out.push_str("0.");
        for _ in 0..(-decpt) {
            out.push('0');
        }
        out.push_str(&digits);
    } else {
        let d = decpt as usize;
        if digits.len() <= d {
            out.push_str(&digits);
            for _ in digits.len()..d {
                out.push('0');
            }
            out.push_str(".0");
        } else {
            out.push_str(&digits[..d]);
            out.push('.');
            out.push_str(&digits[d..]);
        }
    }
    out
}

/// `f"{v:g}"`: six significant digits, trailing zeros removed.
pub fn general(v: f64) -> String {
    if v.is_nan() {
        return "nan".into();
    }
    if v.is_infinite() {
        return if v > 0.0 { "inf".into() } else { "-inf".into() };
    }
    if v == 0.0 {
        return if v.is_sign_negative() { "-0".into() } else { "0".into() };
    }
    const P: i32 = 6;
    let s = format!("{:.*e}", (P - 1) as usize, v);
    let (_, exp) = s.split_once('e').unwrap();
    let x: i32 = exp.parse().unwrap();
    let strip = |t: &str| -> String {
        if t.contains('.') {
            t.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            t.to_string()
        }
    };
    if (-4..P).contains(&x) {
        strip(&format!("{:.*}", (P - 1 - x) as usize, v))
    } else {
        let (m, _) = s.split_once('e').unwrap();
        format!("{}{}", strip(m), exp_suffix(x))
    }
}

/// `f"{v:.0%}"`.
pub fn percent0(v: f64) -> String {
    format!("{:.0}%", v * 100.0)
}

/// `html.escape(s)` with its default `quote=True`.
pub fn html_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#x27;"),
            _ => out.push(c),
        }
    }
    out
}

// ---------------------------------------------------------------------------
// json.dumps
// ---------------------------------------------------------------------------

/// A JSON value that keeps Python's int/float distinction and key order.
#[derive(Clone, Debug, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    Int(i64),
    Float(f64),
    Str(String),
    List(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

impl Json {
    pub fn opt_f(v: Option<f64>) -> Json {
        v.map(Json::Float).unwrap_or(Json::Null)
    }
    pub fn opt_s(v: Option<&str>) -> Json {
        v.map(|s| Json::Str(s.to_string())).unwrap_or(Json::Null)
    }
}

fn json_str(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            ' '..='~' => out.push(c),
            _ => {
                let mut buf = [0u16; 2];
                for unit in c.encode_utf16(&mut buf) {
                    write!(out, "\\u{:04x}", unit).unwrap();
                }
            }
        }
    }
    out.push('"');
}

fn json_float(v: f64) -> String {
    if v.is_nan() {
        "NaN".into()
    } else if v.is_infinite() {
        if v > 0.0 {
            "Infinity".into()
        } else {
            "-Infinity".into()
        }
    } else {
        repr(v)
    }
}

/// `json.dumps(v)` with Python's defaults: `", "` and `": "` separators,
/// `ensure_ascii=True`, and `repr` for floats.
pub fn json_dumps(v: &Json) -> String {
    let mut out = String::new();
    dump(&mut out, v);
    out
}

fn dump(out: &mut String, v: &Json) {
    match v {
        Json::Null => out.push_str("null"),
        Json::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
        Json::Int(i) => write!(out, "{}", i).unwrap(),
        Json::Float(f) => out.push_str(&json_float(*f)),
        Json::Str(s) => json_str(out, s),
        Json::List(items) => {
            out.push('[');
            for (i, it) in items.iter().enumerate() {
                if i > 0 {
                    out.push_str(", ");
                }
                dump(out, it);
            }
            out.push(']');
        }
        Json::Obj(items) => {
            out.push('{');
            for (i, (k, it)) in items.iter().enumerate() {
                if i > 0 {
                    out.push_str(", ");
                }
                json_str(out, k);
                out.push_str(": ");
                dump(out, it);
            }
            out.push('}');
        }
    }
}

// ---------------------------------------------------------------------------
// textwrap.wrap
// ---------------------------------------------------------------------------

const WS: &[char] = &['\t', '\n', '\x0b', '\x0c', '\r', ' '];

fn is_ws(c: char) -> bool {
    WS.contains(&c)
}

/// Python's `\w`: Unicode alphanumerics and underscore.
fn is_word(c: char) -> bool {
    c.is_alphanumeric() || c == '_'
}

/// `[^\d\W]`: a word character that is not a digit.
fn is_letter(c: char) -> bool {
    is_word(c) && !c.is_numeric()
}

/// `[\w!"'&.,?]`
fn is_word_punct(c: char) -> bool {
    is_word(c) || matches!(c, '!' | '"' | '\'' | '&' | '.' | ',' | '?')
}

/// The chunks `TextWrapper._split` produces with `break_on_hyphens=True`,
/// i.e. `wordsep_re.split(text)` with empty strings removed. The regex uses
/// look-behind, which the `regex` crate does not have, so it is hand-coded.
fn split_chunks(text: &str) -> Vec<String> {
    let t: Vec<char> = text.chars().collect();
    let n = t.len();
    let mut out: Vec<String> = Vec::new();
    let mut pending_start = 0usize; // start of text not yet emitted
    let mut i = 0usize;
    // Scan like re.split: at each position try the alternatives in order; on a
    // match, the text before it is a chunk, then the match itself.
    while i < n {
        if let Some(end) = match_at(&t, i) {
            if end > i {
                if pending_start < i {
                    out.push(t[pending_start..i].iter().collect());
                }
                out.push(t[i..end].iter().collect());
                pending_start = end;
                i = end;
                continue;
            }
        }
        i += 1;
    }
    if pending_start < n {
        out.push(t[pending_start..n].iter().collect());
    }
    out.into_iter().filter(|s| !s.is_empty()).collect()
}

/// Length of the `wordsep_re` match starting at `i`, if any.
fn match_at(t: &[char], i: usize) -> Option<usize> {
    let n = t.len();
    // 1. whitespace+
    if is_ws(t[i]) {
        let mut j = i;
        while j < n && is_ws(t[j]) {
            j += 1;
        }
        return Some(j);
    }
    // 2. (?<=wp) -{2,} (?=\w)
    if t[i] == '-' && i > 0 && is_word_punct(t[i - 1]) {
        let mut j = i;
        while j < n && t[j] == '-' {
            j += 1;
        }
        // -{2,} is greedy but backtracks: find the longest run ending before \w
        let mut k = j;
        while k >= i + 2 {
            if k < n && is_word(t[k]) {
                return Some(k);
            }
            k -= 1;
        }
    }
    // 3. nws+? (?: hyphenated | end-of-word | em-dash )  -- lazy, so the
    // shortest prefix for which one of the endings matches.
    let mut j = i; // j = end of the nws+? part
    loop {
        if j >= n || is_ws(t[j]) {
            break;
        }
        j += 1; // consumed one more non-whitespace char
                // (a) hyphenated word: -(?: (?<=lt{2}-) | (?<=lt-lt-)) (?= lt -? lt)
        if j < n && t[j] == '-' {
            let h = j; // position of '-'
            let lb1 = h >= 2 && is_letter(t[h - 1]) && is_letter(t[h - 2]);
            let lb2 = h >= 3 && is_letter(t[h - 1]) && t[h - 2] == '-' && is_letter(t[h - 3]);
            if lb1 || lb2 {
                let a = h + 1;
                let la = a < n
                    && is_letter(t[a])
                    && ((a + 1 < n && is_letter(t[a + 1])) || (a + 2 < n && t[a + 1] == '-' && is_letter(t[a + 2])));
                if la {
                    return Some(h + 1);
                }
            }
        }
        // (b) end of word: (?=ws|\Z)
        if j >= n || is_ws(t[j]) {
            return Some(j);
        }
        // (c) em-dash: (?<=wp) (?=-{2,}\w)
        if is_word_punct(t[j - 1]) && j + 1 < n && t[j] == '-' && t[j + 1] == '-' {
            let mut k = j;
            while k < n && t[k] == '-' {
                k += 1;
            }
            if k < n && is_word(t[k]) {
                return Some(j);
            }
        }
    }
    None
}

/// `textwrap.wrap(text, width)` with the default options.
pub fn wrap(text: &str, width: usize) -> Vec<String> {
    // expand_tabs (tabsize 8), then replace_whitespace
    let mut expanded = String::new();
    let mut col = 0usize;
    for c in text.chars() {
        if c == '\t' {
            let n = 8 - col % 8;
            for _ in 0..n {
                expanded.push(' ');
            }
            col += n;
        } else if c == '\n' || c == '\r' {
            expanded.push(c);
            col = 0;
        } else {
            expanded.push(c);
            col += 1;
        }
    }
    let text: String = expanded.chars().map(|c| if is_ws(c) { ' ' } else { c }).collect();
    let mut chunks: Vec<String> = split_chunks(&text);
    chunks.reverse();
    let mut lines: Vec<String> = Vec::new();
    let clen = |s: &String| s.chars().count();
    while !chunks.is_empty() {
        let mut cur: Vec<String> = Vec::new();
        let mut cur_len = 0usize;
        let w = width;
        if chunks.last().map(|c| c.trim().is_empty()).unwrap_or(false) && !lines.is_empty() {
            chunks.pop();
        }
        while let Some(last) = chunks.last() {
            let l = clen(last);
            if cur_len + l <= w {
                cur_len += l;
                cur.push(chunks.pop().unwrap());
            } else {
                break;
            }
        }
        if let Some(last) = chunks.last() {
            if clen(last) > w {
                // _handle_long_word with break_long_words=True, break_on_hyphens=True
                let space_left = if w < 1 { 1 } else { w - cur_len };
                let chunk: Vec<char> = last.chars().collect();
                let mut end = space_left;
                if chunk.len() > space_left {
                    // rfind('-', 0, space_left) with hyphen > 0 and any non-hyphen before it
                    if let Some(h) = chunk[..space_left.min(chunk.len())].iter().rposition(|&c| c == '-') {
                        if h > 0 && chunk[..h].iter().any(|&c| c != '-') {
                            end = h + 1;
                        }
                    }
                }
                let head: String = chunk[..end.min(chunk.len())].iter().collect();
                let tail: String = chunk[end.min(chunk.len())..].iter().collect();
                cur.push(head);
                cur_len += end.min(chunk.len());
                let idx = chunks.len() - 1;
                chunks[idx] = tail;
                if chunks[idx].is_empty() {
                    chunks.pop();
                }
            } else if cur.is_empty() {
                // nothing fits and it is not long: cannot happen with width >= 1
            }
        }
        if let Some(last) = cur.last() {
            if last.trim().is_empty() {
                cur_len -= clen(last);
                cur.pop();
            }
        }
        let _ = cur_len;
        if !cur.is_empty() {
            lines.push(cur.concat());
        }
    }
    lines
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repr_matches_python() {
        let cases = [
            (0.1, "0.1"),
            (17.0, "17.0"),
            (1e-05, "1e-05"),
            (0.0001, "0.0001"),
            (1e16, "1e+16"),
            (1234567890123456.0, "1234567890123456.0"),
            (1.5e-07, "1.5e-07"),
            (-2.5, "-2.5"),
            (0.0, "0.0"),
            (-0.0, "-0.0"),
            (123456789.125, "123456789.125"),
            (1e22, "1e+22"),
            (0.3333333333333333, "0.3333333333333333"),
        ];
        for (v, want) in cases {
            assert_eq!(repr(v), want, "{v}");
        }
    }

    #[test]
    fn general_matches_python() {
        let cases = [
            (17.0, "17"),
            (17.5, "17.5"),
            (2.0, "2"),
            (1e-05, "1e-05"),
            (0.0001, "0.0001"),
            (123456.0, "123456"),
            (1234567.0, "1.23457e+06"),
            (0.5, "0.5"),
            (100.25, "100.25"),
        ];
        for (v, want) in cases {
            assert_eq!(general(v), want, "{v}");
        }
    }

    #[test]
    fn wrap_breaks_on_hyphens_like_python() {
        // Expected values are Python 3's own output for these inputs.
        assert_eq!(
            wrap("the homozygote/heterozygote decision is per-allele and well-known", 20),
            [
                "the homozygote/heter",
                "ozygote decision is",
                "per-allele and well-",
                "known"
            ]
        );
        assert_eq!(
            wrap("a--b well--known x ---y long-long-long-long-long-word-here", 12),
            [
                "a--b well--",
                "known x ---y",
                "long-long-",
                "long-long-",
                "long-word-",
                "here"
            ]
        );
    }
}

// ---------------------------------------------------------------------------
// Navigating loaded YAML/JSON the way the Python code does with dicts
// ---------------------------------------------------------------------------

impl Json {
    /// `d.get(key)` on a mapping; `None` for a missing key or a non-mapping.
    pub fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Obj(items) => items.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    /// `d.get(key)` treating an explicit null like a missing key.
    pub fn get_some(&self, key: &str) -> Option<&Json> {
        self.get(key).filter(|v| !matches!(v, Json::Null))
    }

    pub fn has(&self, key: &str) -> bool {
        self.get(key).is_some()
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Json::Str(s) => Some(s),
            _ => None,
        }
    }

    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Json::Int(i) => Some(*i as f64),
            Json::Float(f) => Some(*f),
            Json::Bool(b) => Some(if *b { 1.0 } else { 0.0 }),
            _ => None,
        }
    }

    pub fn as_list(&self) -> &[Json] {
        match self {
            Json::List(v) => v,
            _ => &[],
        }
    }

    pub fn entries(&self) -> &[(String, Json)] {
        match self {
            Json::Obj(v) => v,
            _ => &[],
        }
    }

    /// Python truthiness.
    pub fn truthy(&self) -> bool {
        match self {
            Json::Null => false,
            Json::Bool(b) => *b,
            Json::Int(i) => *i != 0,
            Json::Float(f) => *f != 0.0,
            Json::Str(s) => !s.is_empty(),
            Json::List(v) => !v.is_empty(),
            Json::Obj(v) => !v.is_empty(),
        }
    }

    /// `str(v)` for a scalar.
    pub fn py_str(&self) -> String {
        match self {
            Json::Null => "None".into(),
            Json::Bool(b) => {
                if *b {
                    "True".into()
                } else {
                    "False".into()
                }
            }
            Json::Int(i) => i.to_string(),
            Json::Float(f) => repr(*f),
            Json::Str(s) => s.clone(),
            other => json_dumps(other),
        }
    }

    pub fn from_serde(v: &serde_json::Value) -> Json {
        match v {
            serde_json::Value::Null => Json::Null,
            serde_json::Value::Bool(b) => Json::Bool(*b),
            serde_json::Value::Number(n) => {
                if let Some(i) = n.as_i64() {
                    Json::Int(i)
                } else if n.is_f64() {
                    Json::Float(n.as_f64().unwrap())
                } else {
                    Json::Float(n.as_f64().unwrap_or(f64::NAN))
                }
            }
            serde_json::Value::String(s) => Json::Str(s.clone()),
            serde_json::Value::Array(a) => Json::List(a.iter().map(Json::from_serde).collect()),
            serde_json::Value::Object(o) => {
                Json::Obj(o.iter().map(|(k, v)| (k.clone(), Json::from_serde(v))).collect())
            }
        }
    }

    pub fn from_yaml(v: &yaml_rust2::Yaml) -> Json {
        use yaml_rust2::Yaml;
        match v {
            Yaml::Null | Yaml::BadValue | Yaml::Alias(_) => Json::Null,
            Yaml::Boolean(b) => Json::Bool(*b),
            Yaml::Integer(i) => Json::Int(*i),
            Yaml::Real(s) => Json::Float(parse_yaml_float(s)),
            Yaml::String(s) => Json::Str(s.clone()),
            Yaml::Array(a) => Json::List(a.iter().map(Json::from_yaml).collect()),
            Yaml::Hash(h) => Json::Obj(
                h.iter()
                    .map(|(k, v)| {
                        let key = match k {
                            Yaml::String(s) | Yaml::Real(s) => s.clone(),
                            Yaml::Integer(i) => i.to_string(),
                            Yaml::Boolean(b) => b.to_string(),
                            _ => String::new(),
                        };
                        (key, Json::from_yaml(v))
                    })
                    .collect(),
            ),
        }
    }
}

fn parse_yaml_float(s: &str) -> f64 {
    let t = s.replace('_', "");
    match t.to_ascii_lowercase().as_str() {
        ".inf" | "+.inf" => f64::INFINITY,
        "-.inf" => f64::NEG_INFINITY,
        ".nan" => f64::NAN,
        other => other.parse().unwrap_or(f64::NAN),
    }
}

/// A number that remembers whether Python would hold it as an int or a float,
/// because `str()` and `json.dumps` print the two differently.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Num {
    Int(i64),
    Float(f64),
}

impl Num {
    pub fn f(self) -> f64 {
        match self {
            Num::Int(i) => i as f64,
            Num::Float(f) => f,
        }
    }

    pub fn from_json(v: &Json) -> Option<Num> {
        match v {
            Json::Int(i) => Some(Num::Int(*i)),
            Json::Float(f) => Some(Num::Float(*f)),
            Json::Bool(b) => Some(Num::Int(*b as i64)),
            _ => None,
        }
    }

    pub fn json(self) -> Json {
        match self {
            Num::Int(i) => Json::Int(i),
            Num::Float(f) => Json::Float(f),
        }
    }

    /// `str(v)`.
    pub fn py_str(self) -> String {
        match self {
            Num::Int(i) => i.to_string(),
            Num::Float(f) => repr(f),
        }
    }
}

/// `str(v)` for an optional number: `None` prints as "None".
pub fn opt_num_str(v: Option<Num>) -> String {
    v.map(Num::py_str).unwrap_or_else(|| "None".into())
}
