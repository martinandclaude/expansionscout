//! A self-contained review page for one sample.
//!
//! Port of `expansionscout/report.py`, which explains the page. Like the
//! Python, the page is built from the tables just written, read back as
//! text, so it can only show what the output files show. The style and the
//! script are the same files the Python reads, compiled in.

use crate::catalog::{Catalog, Locus};
use crate::pyfmt::{fixed, general, html_escape, json_dumps, Json, Num};
use crate::qc::{f, get, Finding, TsvRow};

pub const CSS: &str = include_str!("../data/report.css");
pub const JS: &str = include_str!("../data/report.js");

const TIER_FINDING: &str = "finding";
const TIER_CHECK: &str = "check";
const TIER_UNINFORMATIVE: &str = "uninformative";
const TIER_NORMAL: &str = "normal";

const TIERS: [(&str, &str, &str); 4] = [
    (TIER_FINDING, "Findings",
     "An allele in the pathogenic range, an expansion seen only as a lower bound, or reads the tool would not summarise as one allele."),
    (TIER_CHECK, "Worth a look",
     "Intermediate or uncertain calls, and calls sitting close enough to a band boundary that the assignment could move."),
    (TIER_UNINFORMATIVE, "Not evidence of a normal result",
     "Called in the normal range, but the reads at this locus could not have shown an expansion. A normal call here is the absence of a measurement, not a negative."),
    (TIER_NORMAL, "Normal, and informative", ""),
];

/// `_s(row, key)`.
fn s(row: &TsvRow, key: &str) -> Option<String> {
    match get(row, key) {
        None | Some("") | Some(".") => None,
        Some(v) => Some(v.to_string()),
    }
}

fn tier(row: &TsvRow) -> &'static str {
    let cls = [
        s(row, "a1_cls").unwrap_or_default(),
        s(row, "a2_cls").unwrap_or_default(),
    ];
    if cls.iter().any(|c| c.contains("path")) {
        return TIER_FINDING;
    }
    if s(row, "evidence").as_deref() == Some("EXPANSION_LB") {
        return TIER_FINDING;
    }
    let truthy = |v: Option<f64>| v.is_some_and(|x| x != 0.0);
    if cls.iter().any(|c| c == "no_call") || truthy(f(row, "a1_dispersed")) || truthy(f(row, "a2_dispersed")) {
        return TIER_FINDING;
    }
    if cls.iter().any(|c| c == "intermediate" || c == "uncertain") {
        return TIER_CHECK;
    }
    if let (Some(a1), Some(a2)) = (f(row, "a1_median"), f(row, "a2_median")) {
        if a1 != a2 && (a2 - a1).abs() <= 1.0 && cls[0] != cls[1] {
            return TIER_CHECK;
        }
    }
    if f(row, "negative_reliable") == Some(0.0) {
        return TIER_UNINFORMATIVE;
    }
    TIER_NORMAL
}

fn call_string(row: &TsvRow) -> Option<String> {
    let a1 = f(row, "a1_median")?;
    Some(match f(row, "a2_median") {
        None => general(a1),
        Some(a2) => format!("{} / {}", general(a1), general(a2)),
    })
}

fn jf(v: Option<f64>) -> Json {
    Json::opt_f(v)
}

fn js(v: Option<String>) -> Json {
    v.map(Json::Str).unwrap_or(Json::Null)
}

fn obj(items: Vec<(&str, Json)>) -> Json {
    Json::Obj(items.into_iter().map(|(k, v)| (k.to_string(), v)).collect())
}

/// Python's `int(float)`, truncating.
fn int_of(v: Option<f64>) -> i64 {
    v.filter(|x| *x != 0.0).map(|x| x.trunc() as i64).unwrap_or(0)
}

fn locus_payload(row: &TsvRow, reads: &[&TsvRow], locus: Option<&Locus>, band_prov: Option<&Vec<Json>>) -> Json {
    let mut alleles = Vec::new();
    for a in ["a1", "a2"] {
        let Some(med) = f(row, &format!("{a}_median")) else {
            continue;
        };
        let g = |k: &str| f(row, &format!("{a}_{k}"));
        alleles.push(obj(vec![
            ("n", Json::Str(a.into())),
            ("median", Json::Float(med)),
            ("p5", jf(g("p5"))),
            ("p95", jf(g("p95"))),
            ("support", jf(g("support"))),
            ("cls", js(s(row, &format!("{a}_cls")))),
            ("dispersed", Json::Bool(g("dispersed").is_some_and(|v| v != 0.0))),
            ("meth_tract", jf(g("meth_tract"))),
            ("meth_up", jf(g("meth_up"))),
            ("n_int", jf(g("n_interruption"))),
            ("n_path", jf(g("n_pathogenic"))),
            ("size_bp", jf(g("size_bp"))),
        ]));
    }

    let mut per_read: Vec<(f64, String, Json)> = Vec::new();
    for r in reads {
        let mut n = f(r, "n_total_units");
        let cls = s(r, "read_class").unwrap_or_else(|| "spanning".into());
        if n.is_none() && (cls == "left_partial" || cls == "right_partial") {
            n = f(r, "lower_bound_units");
        }
        let Some(n) = n else { continue };
        let id: String = s(r, "read_id").unwrap_or_default().chars().take(8).collect();
        per_read.push((
            n,
            id.clone(),
            obj(vec![
                ("id", Json::Str(id)),
                ("n", Json::Float(n)),
                ("cls", Json::Str(cls)),
                ("hp", Json::Int(int_of(f(r, "hp")))),
                ("allele", Json::Int(int_of(f(r, "allele")))),
                ("strand", js(s(r, "strand"))),
                ("mapq", jf(f(r, "mapq"))),
                ("struct", Json::Str(s(r, "structure_string").unwrap_or_default())),
                ("cpg", Json::Str(s(r, "cpg_5mc").unwrap_or_default())),
                ("meth", jf(f(r, "frac_5mc"))),
                ("meth_up", jf(f(r, "frac_5mc_up"))),
                ("lb", jf(f(r, "lower_bound_units"))),
                ("n_path", jf(f(r, "n_pathogenic"))),
            ]),
        ));
    }
    per_read.sort_by(|a, b| (-a.0).partial_cmp(&-b.0).unwrap().then_with(|| a.1.cmp(&b.1)));

    let mut bands: Vec<Json> = Vec::new();
    for b in band_prov.map(|v| v.as_slice()).unwrap_or(&[]) {
        let sub = |k: &str, f: &str| {
            b.get_some(k)
                .filter(|x| x.truthy())
                .and_then(|x| x.get(f))
                .cloned()
                .unwrap_or(Json::Null)
        };
        bands.push(obj(vec![
            ("upper", b.get("upper").cloned().unwrap_or(Json::Null)),
            ("label", b.get("label").cloned().unwrap_or(Json::Null)),
            ("boundary", sub("boundary", "source")),
            ("boundary_quote", sub("boundary", "quote")),
            ("naming", sub("naming", "source")),
            ("naming_quote", sub("naming", "quote")),
            ("inhouse", Json::Bool(b.get("inhouse").is_some_and(Json::truthy))),
        ]));
    }
    if bands.is_empty() {
        if let Some(loc) = locus {
            let num = |v: Option<Num>| v.map(Num::json).unwrap_or(Json::Null);
            for (lo, hi, label) in [
                (loc.benign_min, loc.benign_max, "benign"),
                (loc.intermediate_min, loc.intermediate_max, "intermediate"),
                (loc.pathogenic_min, loc.pathogenic_max, "pathogenic"),
            ] {
                if lo.is_none() && hi.is_none() {
                    continue;
                }
                bands.push(obj(vec![
                    ("lower", num(lo)),
                    ("upper", num(hi)),
                    ("label", Json::Str(label.into())),
                    ("boundary", Json::Str("strchive".into())),
                ]));
            }
        }
    }

    obj(vec![
        ("id", js(s(row, "locus"))),
        ("gene", js(s(row, "gene"))),
        ("chrom", js(s(row, "chrom"))),
        ("start", jf(f(row, "start"))),
        ("end", jf(f(row, "end"))),
        ("motif", js(s(row, "motif_gene"))),
        ("ref_motif", js(s(row, "motif_ref"))),
        ("unit_len", jf(f(row, "unit_len"))),
        ("regime", js(s(row, "regime"))),
        (
            "composition",
            Json::Bool(f(row, "composition").is_some_and(|v| v != 0.0)),
        ),
        ("ploidy", jf(f(row, "ploidy"))),
        ("method", js(s(row, "method"))),
        ("call", js(call_string(row))),
        ("tier", Json::Str(tier(row).into())),
        ("n_spanning", jf(f(row, "n_spanning"))),
        ("n_split", jf(f(row, "n_split"))),
        ("n_partial", jf(f(row, "n_partial"))),
        ("lb", jf(f(row, "lower_bound_units"))),
        ("lb_support", jf(f(row, "lower_bound_support"))),
        ("detect", jf(f(row, "detectability"))),
        ("evidence", js(s(row, "evidence"))),
        ("negative_reliable", jf(f(row, "negative_reliable"))),
        ("tagged_fraction", jf(f(row, "tagged_fraction"))),
        ("notes", js(s(row, "notes"))),
        ("alleles", Json::List(alleles)),
        ("reads", Json::List(per_read.into_iter().map(|r| r.2).collect())),
        ("bands", Json::List(bands)),
    ])
}

fn opt_str(v: Option<&Json>) -> Option<String> {
    v.and_then(|x| x.as_str()).map(String::from)
}

fn row_html(loc: &Json) -> String {
    let alleles = loc.get("alleles").map(Json::as_list).unwrap_or(&[]);
    let chips: Vec<String> = alleles
        .iter()
        .filter_map(|a| {
            let c = opt_str(a.get("cls")).filter(|c| !c.is_empty())?;
            Some(format!(
                "<span class=\"cls {}\">{}</span>",
                html_escape(&c),
                html_escape(&c.replace('_', " "))
            ))
        })
        .collect();
    let num = |k: &str| loc.get(k).and_then(Json::as_f64);
    let chrom = opt_str(loc.get("chrom")).unwrap_or_default();
    let coord = match (num("start"), num("end")) {
        (Some(a), Some(b)) if !chrom.is_empty() && a != 0.0 && b != 0.0 => {
            format!("{chrom}:{}-{}", a.trunc() as i64, b.trunc() as i64)
        }
        _ => String::new(),
    };
    let lb = num("lb")
        .map(|v| format!("&ge;&nbsp;{}", fixed(v, 0)))
        .unwrap_or_default();
    let detect = num("detect")
        .map(|v| format!("{}&nbsp;%", fixed(100.0 * v, 0)))
        .unwrap_or_default();
    let mut meth = String::new();
    for a in alleles {
        let v = a
            .get("meth_up")
            .and_then(Json::as_f64)
            .or_else(|| a.get("meth_tract").and_then(Json::as_f64));
        if let Some(v) = v {
            if !meth.is_empty() {
                meth.push_str(" / ");
            }
            meth.push_str(&fixed(100.0 * v, 0));
        }
    }
    if !meth.is_empty() {
        meth.push_str("&nbsp;%");
    }
    let reads = (num("n_spanning").unwrap_or(0.0) + num("n_split").unwrap_or(0.0)).trunc() as i64;
    let dash = "<span class=muted>&mdash;</span>";
    let or_dash = |v: &str| if v.is_empty() { dash.to_string() } else { v.to_string() };
    let call = opt_str(loc.get("call")).filter(|c| !c.is_empty());
    format!(
        "<tr class=\"locus\" data-locus=\"{}\"><td class=\"gene\">{}</td><td class=\"muted\">{}</td>\
         <td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td>\
         <td><button class=\"coord\" data-coord=\"{coord}\">{coord}</button></td>\
         </tr><tr class=\"detail\" hidden><td colspan=\"9\"></td></tr>",
        html_escape(&opt_str(loc.get("id")).unwrap_or_default()),
        html_escape(&opt_str(loc.get("gene")).unwrap_or_default()),
        html_escape(&opt_str(loc.get("id")).unwrap_or_default()),
        call.map(|c| html_escape(&c)).unwrap_or_else(|| "&mdash;".into()),
        if chips.is_empty() {
            dash.to_string()
        } else {
            chips.join(" ")
        },
        reads,
        or_dash(&lb),
        or_dash(&detect),
        or_dash(&meth),
    )
}

/// The page as a string. `date` is the ISO date printed in the subtitle.
pub fn build_report(
    loci_rows: &[TsvRow],
    per_read_rows: &[TsvRow],
    findings: &[Finding],
    sample: &str,
    catalog: &Catalog,
    meta: &[(&str, String)],
    date: &str,
) -> String {
    let mut payload = Vec::new();
    for row in loci_rows {
        let id = get(row, "locus");
        let reads: Vec<&TsvRow> = per_read_rows.iter().filter(|r| get(r, "locus") == id).collect();
        let locus = id.and_then(|i| catalog.loci.iter().find(|l| l.id == i));
        let prov = id
            .and_then(|i| catalog.provenance.bands.iter().find(|b| b.0 == i))
            .map(|b| &b.1);
        payload.push(locus_payload(row, &reads, locus, prov));
    }
    let tier_of = |p: &Json| opt_str(p.get("tier")).unwrap_or_default();
    let mut out: Vec<String> = vec![
        "<!doctype html><html lang=en><head><meta charset=utf-8>".into(),
        "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">".into(),
        format!("<title>{} &middot; repeat review</title>", html_escape(sample)),
        format!("<style>{CSS}</style></head><body><div class=wrap>"),
        format!("<h1>{} &middot; repeat expansion review</h1>", html_escape(sample)),
    ];
    let n_find = payload.iter().filter(|p| tier_of(p) == TIER_FINDING).count();
    out.push(format!(
        "<p class=\"sub\">{} loci &middot; {n_find} needing review &middot; ExpansionScout {} &middot; {date}</p>",
        payload.len(),
        env!("CARGO_PKG_VERSION")
    ));
    if !findings.is_empty() {
        out.push("<div class=\"qc\"><h2>Checks this run made on itself</h2>".into());
        for x in findings {
            let klass = match x.status.as_str() {
                "warn" => "warn",
                "fail" => "fail",
                "pass" => "pass",
                _ => "",
            };
            out.push(format!(
                "<div class=\"f\"><span class=\"tag {klass}\">{}</span><span>{}</span><span class=\"basis\">{}</span></div>",
                html_escape(&x.status), html_escape(&x.message), html_escape(&x.basis)));
        }
        out.push("</div>".into());
    }
    let head = "<tr><th>gene</th><th>locus</th><th>call</th><th>class</th><th>reads</th>\
                <th>bound</th><th>detect</th><th>5mC</th><th>position</th></tr>";
    for (t, title, why) in TIERS {
        let rows: Vec<&Json> = payload.iter().filter(|p| tier_of(p) == t).collect();
        if rows.is_empty() {
            continue;
        }
        out.push(format!(
            "<section><div class=\"tierhead\"><h2>{title}</h2><span class=\"count\">{}</span></div>",
            rows.len()
        ));
        if !why.is_empty() {
            out.push(format!("<p class=\"tierwhy\">{why}</p>"));
        }
        out.push(format!("<table>{head}"));
        out.extend(rows.iter().map(|p| row_html(p)));
        out.push("</table></section>".into());
    }
    out.push("<h2>Run provenance</h2><details class=\"prov\"><summary>What produced this page</summary><dl>".into());
    for (k, v) in meta {
        out.push(format!("<dt>{}</dt><dd>{}</dd>", html_escape(k), html_escape(v)));
    }
    out.push("</dl></details>".into());
    out.push("<footer>Every locus in the catalogue is on this page, including those that could not be called \
              and those called normal on reads that could not have shown an expansion. Nothing is filtered away.</footer>".into());
    out.push("</div><div id=tip></div>".into());
    let data = Json::Obj(vec![
        ("loci".into(), Json::List(payload)),
        ("sources".into(), catalog.provenance.sources.clone()),
    ]);
    out.push(format!("<script>const DATA={};</script>", json_dumps(&data)));
    out.push(format!("<script>{JS}</script></body></html>"));
    out.join("\n")
}

/// Today's date, or the date of SOURCE_DATE_EPOCH when it is set, in UTC.
pub fn report_date() -> String {
    let secs = std::env::var("SOURCE_DATE_EPOCH")
        .ok()
        .and_then(|v| v.trim().parse::<i64>().ok())
        .unwrap_or_else(|| {
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_secs() as i64)
                .unwrap_or(0)
        });
    let days = secs.div_euclid(86_400);
    // Howard Hinnant's civil_from_days
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    format!("{:04}-{:02}-{:02}", if m <= 2 { y + 1 } else { y }, m, d)
}
