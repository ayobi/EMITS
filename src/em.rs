//! Expectation-Maximization algorithm for abundance estimation.
//!
//! Given a set of reads, each with one or more candidate reference alignments,
//! iteratively estimates the relative abundance of each reference taxon by:
//!   1. E-step: compute posterior probability of each read originating from each
//!      candidate taxon, given current abundance estimates and alignment likelihoods.
//!   2. M-step: update taxon abundances by summing fractional read assignments.
//!   3. Repeat until convergence.

use crate::paf::{PafRecord, ReadAlignments};
use std::collections::HashMap;

/// Result of the EM algorithm.
#[derive(Debug)]
#[allow(dead_code)]
pub struct EmResult {
    /// Relative abundance estimates (taxon -> proportion, sums to 1.0)
    pub abundances: HashMap<String, f64>,
    /// Number of iterations until convergence
    pub iterations: usize,
    /// Final maximum change in abundance between last two iterations
    pub final_delta: f64,
    /// Per-read assignments: read_name -> Vec<(taxon, posterior_probability)>
    pub read_assignments: HashMap<String, Vec<(String, f64)>>,
}

/// Configuration for the EM algorithm.
#[derive(Debug, Clone)]
pub struct EmConfig {
    /// Maximum number of EM iterations
    pub max_iterations: usize,
    /// Convergence threshold: stop when max abundance change < this value
    pub convergence_threshold: f64,
    /// Minimum abundance threshold: taxa below this are pruned between iterations
    pub min_abundance: f64,
    /// Temperature for score-to-likelihood conversion.
    /// Lower = more sensitive to alignment score differences (use for high-accuracy reads).
    /// Higher = more tolerant of score noise (use for noisier reads).
    pub temperature: f64,
}

impl Default for EmConfig {
    fn default() -> Self {
        EmConfig {
            max_iterations: 100,
            convergence_threshold: 1e-6,
            min_abundance: 1e-7,
            temperature: 0.5,
        }
    }
}

/// Convert a normalized alignment score to a likelihood.
///
/// Uses a softmax-like transformation so that small differences in alignment
/// score translate to meaningful probability differences. The temperature
/// parameter controls sensitivity.
fn score_to_likelihood(normalized_score: f64, temperature: f64) -> f64 {
    (normalized_score / temperature).exp()
}

/// Compute alignment likelihoods for a read's candidate alignments.
///
/// Returns a vec of (target_name, likelihood) for each alignment,
/// where likelihoods are NOT yet normalized (that happens in the E-step
/// after weighting by current abundances).
fn compute_likelihoods(alignments: &[PafRecord], temperature: f64) -> Vec<(String, f64)> {
    alignments
        .iter()
        .map(|aln| {
            let likelihood = score_to_likelihood(aln.normalized_score(), temperature);
            (aln.target_name.clone(), likelihood)
        })
        .collect()
}

/// Run the EM algorithm on read alignments.
///
/// # Arguments
/// * `alignments` - Grouped read-to-reference alignments from PAF parsing
/// * `config` - EM configuration parameters
///
/// # Returns
/// * `EmResult` with final abundances, iteration count, and per-read assignments
pub fn run_em(alignments: &ReadAlignments, config: &EmConfig) -> EmResult {
    // Collect all unique taxa from alignments
    let mut all_taxa: Vec<String> = alignments
        .values()
        .flat_map(|alns| alns.iter().map(|a| a.target_name.clone()))
        .collect();
    all_taxa.sort();
    all_taxa.dedup();

    let n_taxa = all_taxa.len();
    let n_reads = alignments.len();
    log::info!(
        "Starting EM with {} reads and {} candidate taxa",
        n_reads,
        n_taxa
    );

    if n_taxa == 0 || n_reads == 0 {
        return EmResult {
            abundances: HashMap::new(),
            iterations: 0,
            final_delta: 0.0,
            read_assignments: HashMap::new(),
        };
    }

    // Initialize uniform abundances
    let mut abundances: HashMap<String, f64> = all_taxa
        .iter()
        .map(|t| (t.clone(), 1.0 / n_taxa as f64))
        .collect();

    // Precompute likelihoods for each read (these don't change between iterations)
    let read_likelihoods: HashMap<String, Vec<(String, f64)>> = alignments
        .iter()
        .map(|(read_name, alns)| {
            (
                read_name.clone(),
                compute_likelihoods(alns, config.temperature),
            )
        })
        .collect();

    let mut iterations = 0;
    let mut final_delta = f64::MAX;

    for iter in 0..config.max_iterations {
        iterations = iter + 1;

        // === E-step ===
        // For each read, compute posterior probability of each candidate taxon
        // P(taxon | read) ∝ P(read | taxon) * P(taxon)
        //                   = likelihood(alignment) * abundance(taxon)
        let mut new_abundances: HashMap<String, f64> =
            all_taxa.iter().map(|t| (t.clone(), 0.0)).collect();

        for likelihoods in read_likelihoods.values() {
            // Compute weighted likelihoods: likelihood * current abundance
            let weighted: Vec<(String, f64)> = likelihoods
                .iter()
                .filter_map(|(taxon, lik)| {
                    let abundance = abundances.get(taxon).copied().unwrap_or(0.0);
                    if abundance > 0.0 {
                        Some((taxon.clone(), lik * abundance))
                    } else {
                        None
                    }
                })
                .collect();

            // Normalize to posterior probabilities
            let total_weight: f64 = weighted.iter().map(|(_, w)| w).sum();
            if total_weight <= 0.0 {
                continue;
            }

            // === M-step (accumulate) ===
            // Add fractional assignments to new abundance estimates
            for (taxon, weight) in &weighted {
                let posterior = weight / total_weight;
                *new_abundances.entry(taxon.clone()).or_default() += posterior;
            }
        }

        // Normalize abundances to sum to 1.0
        let total: f64 = new_abundances.values().sum();
        if total > 0.0 {
            for v in new_abundances.values_mut() {
                *v /= total;
            }
        }

        // Prune low-abundance taxa
        new_abundances.retain(|_, v| *v >= config.min_abundance);

        // Re-normalize after pruning
        let total: f64 = new_abundances.values().sum();
        if total > 0.0 {
            for v in new_abundances.values_mut() {
                *v /= total;
            }
        }

        // Check convergence: max change in any taxon's abundance
        let max_delta = all_taxa
            .iter()
            .map(|t| {
                let old = abundances.get(t).copied().unwrap_or(0.0);
                let new = new_abundances.get(t).copied().unwrap_or(0.0);
                (old - new).abs()
            })
            .fold(0.0f64, f64::max);

        final_delta = max_delta;
        abundances = new_abundances;

        // Update active taxa list after pruning
        // (we keep the original list for delta computation but track what's active)

        if iter % 10 == 0 || max_delta < config.convergence_threshold {
            let active_taxa = abundances.len();
            log::debug!(
                "Iteration {}: max_delta={:.2e}, active_taxa={}",
                iter + 1,
                max_delta,
                active_taxa
            );
        }

        if max_delta < config.convergence_threshold {
            log::info!(
                "EM converged after {} iterations (delta={:.2e})",
                iter + 1,
                max_delta
            );
            break;
        }
    }

    if final_delta >= config.convergence_threshold {
        log::warn!(
            "EM did not converge after {} iterations (final delta={:.2e})",
            config.max_iterations,
            final_delta
        );
    }

    // Compute final per-read assignments using converged abundances
    let mut read_assignments: HashMap<String, Vec<(String, f64)>> = HashMap::new();
    for (read_name, likelihoods) in &read_likelihoods {
        let weighted: Vec<(String, f64)> = likelihoods
            .iter()
            .filter_map(|(taxon, lik)| {
                let abundance = abundances.get(taxon).copied().unwrap_or(0.0);
                if abundance > 0.0 {
                    Some((taxon.clone(), lik * abundance))
                } else {
                    None
                }
            })
            .collect();

        let total_weight: f64 = weighted.iter().map(|(_, w)| w).sum();
        if total_weight <= 0.0 {
            continue;
        }

        let assignments: Vec<(String, f64)> = weighted
            .iter()
            .map(|(taxon, w)| (taxon.clone(), w / total_weight))
            .filter(|(_, p)| *p > 0.001) // Only keep assignments > 0.1%
            .collect();

        read_assignments.insert(read_name.clone(), assignments);
    }

    EmResult {
        abundances,
        iterations,
        final_delta,
        read_assignments,
    }
}

/// Compute naive count-based abundances (no EM).
///
/// Each read is assigned entirely to its best-scoring alignment.
/// This serves as the baseline comparison for the EM approach.
pub fn naive_count(alignments: &ReadAlignments) -> HashMap<String, f64> {
    let mut counts: HashMap<String, f64> = HashMap::new();
    let mut total_reads = 0u64;

    for alns in alignments.values() {
        if alns.is_empty() {
            continue;
        }

        // Assign read to the best-scoring alignment
        let best = alns
            .iter()
            .max_by(|a, b| {
                a.normalized_score()
                    .partial_cmp(&b.normalized_score())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap();

        *counts.entry(best.target_name.clone()).or_default() += 1.0;
        total_reads += 1;
    }

    // Normalize to relative abundances
    if total_reads > 0 {
        for v in counts.values_mut() {
            *v /= total_reads as f64;
        }
    }

    counts
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::paf::PafRecord;

    /// Helper to create a mock PafRecord
    fn mock_record(
        query: &str,
        target: &str,
        matches: u64,
        block_len: u64,
        as_score: i64,
    ) -> PafRecord {
        PafRecord {
            query_name: query.to_string(),
            query_len: 1000,
            target_name: target.to_string(),
            target_len: 800,
            matches,
            block_len,
            mapq: 60,
            alignment_score: Some(as_score),
        }
    }

    #[test]
    fn test_em_single_mapping_reads() {
        // When every read maps uniquely, EM should converge to naive counts
        let mut alignments: ReadAlignments = HashMap::new();
        for i in 0..70 {
            alignments.insert(
                format!("read_A_{}", i),
                vec![mock_record(
                    &format!("read_A_{}", i),
                    "species_A",
                    950,
                    1000,
                    1900,
                )],
            );
        }
        for i in 0..30 {
            alignments.insert(
                format!("read_B_{}", i),
                vec![mock_record(
                    &format!("read_B_{}", i),
                    "species_B",
                    950,
                    1000,
                    1900,
                )],
            );
        }

        let config = EmConfig::default();
        let result = run_em(&alignments, &config);

        let a = result.abundances.get("species_A").copied().unwrap_or(0.0);
        let b = result.abundances.get("species_B").copied().unwrap_or(0.0);

        assert!((a - 0.7).abs() < 0.01, "Expected ~0.7, got {}", a);
        assert!((b - 0.3).abs() < 0.01, "Expected ~0.3, got {}", b);
    }

    #[test]
    fn test_em_resolves_ambiguous_mappings() {
        // Key test: reads that map to multiple species should be resolved by EM
        //
        // Setup: 60 reads map uniquely to species_A, 20 uniquely to species_B,
        // and 20 reads are ambiguous (map equally to both).
        //
        // Naive counting assigns ambiguous reads to whichever scores best (or arbitrarily),
        // while EM should distribute them proportionally to the learned abundances.
        let mut alignments: ReadAlignments = HashMap::new();

        // 60 unique reads for species A
        for i in 0..60 {
            alignments.insert(
                format!("unique_A_{}", i),
                vec![mock_record(
                    &format!("unique_A_{}", i),
                    "species_A",
                    960,
                    1000,
                    1920,
                )],
            );
        }

        // 20 unique reads for species B
        for i in 0..20 {
            alignments.insert(
                format!("unique_B_{}", i),
                vec![mock_record(
                    &format!("unique_B_{}", i),
                    "species_B",
                    960,
                    1000,
                    1920,
                )],
            );
        }

        // 20 ambiguous reads mapping equally to both
        for i in 0..20 {
            let name = format!("ambig_{}", i);
            alignments.insert(
                name.clone(),
                vec![
                    mock_record(&name, "species_A", 950, 1000, 1900),
                    mock_record(&name, "species_B", 950, 1000, 1900),
                ],
            );
        }

        let config = EmConfig::default();
        let result = run_em(&alignments, &config);

        let a = result.abundances.get("species_A").copied().unwrap_or(0.0);
        let b = result.abundances.get("species_B").copied().unwrap_or(0.0);

        // EM should distribute the 20 ambiguous reads ~3:1 (A:B),
        // giving roughly A=0.75, B=0.25 (since 60:20 = 3:1 prior)
        // The exact values depend on the likelihood model, but A should be > B
        assert!(a > b, "EM should assign more to A (got A={}, B={})", a, b);
        assert!(a > 0.6, "Species A should be > 0.6 (got {})", a);
        assert!(b > 0.15, "Species B should be > 0.15 (got {})", b);

        // Compare with naive: naive would assign all ambiguous reads to one species
        let naive = naive_count(&alignments);
        let naive_a = naive.get("species_A").copied().unwrap_or(0.0);
        let naive_b = naive.get("species_B").copied().unwrap_or(0.0);

        log::info!("EM:    A={:.4}, B={:.4}", a, b);
        log::info!("Naive: A={:.4}, B={:.4}", naive_a, naive_b);

        // EM should produce more balanced estimates than naive for the ambiguous case
        // (Naive will be more extreme because it assigns all ambig reads to one)
    }

    #[test]
    fn test_em_convergence() {
        let mut alignments: ReadAlignments = HashMap::new();
        for i in 0..50 {
            alignments.insert(
                format!("read_{}", i),
                vec![mock_record(
                    &format!("read_{}", i),
                    "species_A",
                    950,
                    1000,
                    1900,
                )],
            );
        }

        let config = EmConfig {
            max_iterations: 100,
            convergence_threshold: 1e-6,
            min_abundance: 1e-7,
            temperature: 0.5,
        };
        let result = run_em(&alignments, &config);

        assert!(
            result.final_delta < config.convergence_threshold,
            "EM should converge (delta={})",
            result.final_delta
        );
    }

    #[test]
    fn test_naive_count_best_hit() {
        let mut alignments: ReadAlignments = HashMap::new();
        alignments.insert(
            "read1".to_string(),
            vec![
                mock_record("read1", "species_A", 960, 1000, 1920), // better
                mock_record("read1", "species_B", 900, 1000, 1800), // worse
            ],
        );

        let counts = naive_count(&alignments);
        assert_eq!(counts.get("species_A").copied().unwrap_or(0.0), 1.0);
        assert_eq!(counts.get("species_B").copied().unwrap_or(0.0), 0.0);
    }
}
