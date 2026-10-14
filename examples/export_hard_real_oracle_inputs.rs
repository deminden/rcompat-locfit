//! Export the committed hard-real fixture as exact IEEE-754 input bits for R.
//! Run from the repository root. No R package or extra Rust dependency is needed.
use std::fs;

fn bits(value: &str) -> String {
    if value.is_empty() {
        String::new()
    } else {
        value
            .parse::<f64>()
            .expect("invalid fixture number")
            .to_bits()
            .to_le_bytes()
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect()
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let source = "fixtures/hard_real_2026_subset_cases.csv";
    let output = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "target/hard_real_oracle_inputs.csv".into());
    let text = fs::read_to_string(source)?;
    let mut rows = String::from("# input_encoding=IEEE754_binary64_little_endian_hex\n# source=fixtures/hard_real_2026_subset_cases.csv\ncase,kind,index,label,mean_bits,disp_bits,log_mean_bits,log_disp_bits\n");
    for line in text
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
    {
        let fields: Vec<_> = line.split(',').collect();
        assert_eq!(fields.len(), 10, "invalid fixture row");
        rows.push_str(&format!(
            "{},{},{},{},{},{},{},{}\n",
            fields[0],
            fields[1],
            fields[2],
            fields[3],
            bits(fields[4]),
            bits(fields[5]),
            bits(fields[6]),
            bits(fields[7])
        ));
    }
    fs::write(&output, rows)?;
    println!("Wrote {output}");
    Ok(())
}
