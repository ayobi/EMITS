//! Output formatting for abundance estimates.

use crate::taxonomy::AggregatedResult;
use anyhow::Result;
use std::collections::HashMap;
use std::io::Write;

/// Write raw per-accession abundance estimates to a TSV file.
///
/// Columns: taxon, relative_abundance, estimated_reads
pub fn write_abundance_tsv<W: Write>(
    writer: &mut W,
    abundances: &HashMap<String, f64>,
    total_reads: usize,
) -> Result<()> {
    writeln!(writer, "taxon\trelative_abundance\testimated_reads")?;

    let mut sorted: Vec<_> = abundances.iter().collect();
    sorted.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap_or(std::cmp::Ordering::Equal));

    for (taxon, abundance) in sorted {
        let estimated_reads = abundance * total_reads as f64;
        writeln!(
            writer,
            "{}\t{:.6}\t{:.1}",
            taxon, abundance, estimated_reads
        )?;
    }

    Ok(())
}

/// Write aggregated abundance estimates to a TSV file.
///
/// Columns: taxon, lineage, relative_abundance, estimated_reads, n_accessions
pub fn write_aggregated_tsv<W: Write>(
    writer: &mut W,
    result: &AggregatedResult,
    total_reads: usize,
) -> Result<()> {
    writeln!(
        writer,
        "taxon\tlineage\trelative_abundance\testimated_reads\tn_accessions"
    )?;

    let mut sorted: Vec<_> = result.abundances.iter().collect();
    sorted.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap_or(std::cmp::Ordering::Equal));

    for (taxon, abundance) in sorted {
        let estimated_reads = abundance * total_reads as f64;
        let lineage = result.lineages.get(taxon).map(|s| s.as_str()).unwrap_or("");
        let n_acc = result.accession_counts.get(taxon).copied().unwrap_or(0);
        writeln!(
            writer,
            "{}\t{}\t{:.6}\t{:.1}\t{}",
            taxon, lineage, abundance, estimated_reads, n_acc
        )?;
    }

    Ok(())
}

/// Write a comparison table (truth vs naive vs EM) to TSV.
pub fn write_comparison_tsv<W: Write>(
    writer: &mut W,
    truth: &HashMap<String, f64>,
    naive: &HashMap<String, f64>,
    em: &HashMap<String, f64>,
) -> Result<()> {
    writeln!(
        writer,
        "taxon\ttrue_abundance\tnaive_abundance\tem_abundance\tnaive_error\tem_error"
    )?;

    let all_taxa: std::collections::HashSet<&String> =
        truth.keys().chain(naive.keys()).chain(em.keys()).collect();
    let mut sorted_taxa: Vec<_> = all_taxa.into_iter().collect();
    sorted_taxa.sort();

    for taxon in sorted_taxa {
        let t = truth.get(taxon).copied().unwrap_or(0.0);
        let n = naive.get(taxon).copied().unwrap_or(0.0);
        let e = em.get(taxon).copied().unwrap_or(0.0);
        writeln!(
            writer,
            "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            taxon,
            t,
            n,
            e,
            (n - t).abs(),
            (e - t).abs()
        )?;
    }

    Ok(())
}

/// Write an aggregated comparison table to TSV.
pub fn write_aggregated_comparison_tsv<W: Write>(
    writer: &mut W,
    em_agg: &AggregatedResult,
    naive_agg: &AggregatedResult,
    total_reads: usize,
) -> Result<()> {
    writeln!(
        writer,
        "taxon\tlineage\tem_abundance\tnaive_abundance\tem_reads\tnaive_reads\tn_accessions"
    )?;

    let all_taxa: std::collections::HashSet<&String> = em_agg
        .abundances
        .keys()
        .chain(naive_agg.abundances.keys())
        .collect();
    let mut sorted: Vec<_> = all_taxa.into_iter().collect();
    sorted.sort_by(|a, b| {
        let ea = em_agg.abundances.get(*b).unwrap_or(&0.0);
        let eb = em_agg.abundances.get(*a).unwrap_or(&0.0);
        ea.partial_cmp(eb).unwrap_or(std::cmp::Ordering::Equal)
    });

    for taxon in sorted {
        let em_ab = em_agg.abundances.get(taxon).copied().unwrap_or(0.0);
        let naive_ab = naive_agg.abundances.get(taxon).copied().unwrap_or(0.0);
        let lineage = em_agg
            .lineages
            .get(taxon)
            .or_else(|| naive_agg.lineages.get(taxon))
            .map(|s| s.as_str())
            .unwrap_or("");
        let n_acc = em_agg.accession_counts.get(taxon).copied().unwrap_or(0);
        writeln!(
            writer,
            "{}\t{}\t{:.6}\t{:.6}\t{:.1}\t{:.1}\t{}",
            taxon,
            lineage,
            em_ab,
            naive_ab,
            em_ab * total_reads as f64,
            naive_ab * total_reads as f64,
            n_acc
        )?;
    }

    Ok(())
}

/// Print a summary table to stdout showing aggregated results.
pub fn print_aggregated_summary(
    em_agg: &AggregatedResult,
    naive_agg: Option<&AggregatedResult>,
    top_n: usize,
) {
    let mut sorted: Vec<_> = em_agg.abundances.iter().collect();
    sorted.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap_or(std::cmp::Ordering::Equal));

    let rank_name = match em_agg.rank {
        crate::taxonomy::TaxRank::Kingdom => "Kingdom",
        crate::taxonomy::TaxRank::Phylum => "Phylum",
        crate::taxonomy::TaxRank::Class => "Class",
        crate::taxonomy::TaxRank::Order => "Order",
        crate::taxonomy::TaxRank::Family => "Family",
        crate::taxonomy::TaxRank::Genus => "Genus",
        crate::taxonomy::TaxRank::Species => "Species",
    };

    if let Some(naive) = naive_agg {
        println!("\nTop {} ({}-level):", top_n, rank_name);
        println!(
            "  {:<40} {:>10} {:>10} {:>6}",
            rank_name, "EM", "Naive", "Refs"
        );
        println!("  {}", "-".repeat(68));
        for (taxon, em_ab) in sorted.iter().take(top_n) {
            let naive_ab = naive.abundances.get(*taxon).copied().unwrap_or(0.0);
            let n_acc = em_agg.accession_counts.get(*taxon).copied().unwrap_or(0);
            println!(
                "  {:<40} {:>9.2}% {:>9.2}% {:>6}",
                taxon,
                *em_ab * 100.0,
                naive_ab * 100.0,
                n_acc
            );
        }
    } else {
        println!("\nTop {} ({}-level):", top_n, rank_name);
        println!("  {:<40} {:>10} {:>6}", rank_name, "EM", "Refs");
        println!("  {}", "-".repeat(58));
        for (taxon, em_ab) in sorted.iter().take(top_n) {
            let n_acc = em_agg.accession_counts.get(*taxon).copied().unwrap_or(0);
            println!("  {:<40} {:>9.2}% {:>6}", taxon, *em_ab * 100.0, n_acc);
        }
    }
}
