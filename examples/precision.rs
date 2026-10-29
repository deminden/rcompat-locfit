//! Report precision against the R 4.6.2 DESeq2 fixtures.
#[path = "../tests/support/mod.rs"]
mod support;

use std::collections::BTreeMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let fixture = args.next().unwrap_or_else(|| support::HARD_FIXTURE.into());
    if fixture.starts_with('-') || args.next().is_some() {
        return Err("usage: precision [fixture.csv]".into());
    }
    let rows = support::evaluate(&support::load(fixture))?;
    let mut maxima: BTreeMap<&str, (f64, f64)> = BTreeMap::new();
    for row in &rows {
        let entry = maxima.entry(&row.case).or_default();
        entry.0 = entry.0.max(row.log_error);
        entry.1 = entry.1.max(row.normal_log_error);
    }
    for (name, (log, normal)) in maxima {
        println!("{name}: max log error={log:.17e}, normal-input log error={normal:.17e}");
    }
    Ok(())
}
