//! Simulation module for generating synthetic PAF-like alignments.
//!
//! Creates in-silico test scenarios with known ground truth community
//! compositions to validate that the EM algorithm improves over naive counting.

use crate::paf::{PafRecord, ReadAlignments};
use std::collections::HashMap;

/// A simulated community with known composition.
#[derive(Debug)]
pub struct SimulatedCommunity {
    /// Ground truth relative abundances (taxon -> proportion)
    pub true_abundances: HashMap<String, f64>,
    /// Simulated read alignments
    pub alignments: ReadAlignments,
    /// Total number of reads generated
    pub n_reads: usize,
}

/// Configuration for community simulation.
#[derive(Debug, Clone)]
pub struct SimConfig {
    /// Total number of reads to simulate
    pub n_reads: usize,
    /// Probability that a read aligns to a closely related species in addition
    /// to its true source (simulates multi-mapping due to ITS similarity)
    pub cross_mapping_rate: f64,
    /// Base alignment score for correct alignments
    pub true_alignment_score: i64,
    /// Alignment score for cross-mapped (incorrect) alignments
    /// Should be slightly lower than true_alignment_score
    pub cross_alignment_score: i64,
    /// Random seed for reproducibility
    pub seed: u64,
}

impl Default for SimConfig {
    fn default() -> Self {
        SimConfig {
            n_reads: 1000,
            cross_mapping_rate: 0.3,
            true_alignment_score: 1900,
            cross_alignment_score: 1850,
            seed: 42,
        }
    }
}

/// Simple deterministic pseudo-random number generator (xorshift64).
/// We avoid external dependencies for the PoC.
struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        SimpleRng {
            state: if seed == 0 { 1 } else { seed },
        }
    }

    fn next_u64(&mut self) -> u64 {
        self.state ^= self.state << 13;
        self.state ^= self.state >> 7;
        self.state ^= self.state << 17;
        self.state
    }

    /// Returns a float in [0, 1)
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Choose an index from a weighted distribution.
    fn weighted_choice(&mut self, weights: &[f64]) -> usize {
        let total: f64 = weights.iter().sum();
        let mut r = self.next_f64() * total;
        for (i, w) in weights.iter().enumerate() {
            r -= w;
            if r <= 0.0 {
                return i;
            }
        }
        weights.len() - 1
    }
}

/// Define a group of closely related species that can cross-map.
/// This simulates the real-world scenario where ITS sequences from
/// related fungi are similar enough to cause alignment ambiguity.
#[derive(Debug, Clone)]
pub struct SpeciesGroup {
    /// Species names in this group
    pub species: Vec<String>,
    /// True relative abundances for each species
    pub abundances: Vec<f64>,
}

/// Simulate a fungal community with known composition and realistic
/// cross-mapping patterns.
///
/// # Design
///
/// Species are organized into groups of closely related taxa. Reads from
/// any species in a group have a chance of also aligning to other species
/// in the same group (simulating ITS sequence similarity). Reads from
/// species in different groups do NOT cross-map (simulating distinct ITS
/// sequences between distant taxa).
///
/// This is a simplification but captures the key dynamic that makes EM
/// useful: ambiguous alignments between closely related species.
pub fn simulate_community(groups: &[SpeciesGroup], config: &SimConfig) -> SimulatedCommunity {
    let mut rng = SimpleRng::new(config.seed);

    // Build true abundance map (normalize across all groups)
    let mut true_abundances: HashMap<String, f64> = HashMap::new();
    let mut all_species: Vec<(String, f64, usize)> = Vec::new(); // (name, abundance, group_idx)
    let total_abundance: f64 = groups.iter().flat_map(|g| g.abundances.iter()).sum();

    for (gi, group) in groups.iter().enumerate() {
        for (si, species) in group.species.iter().enumerate() {
            let norm_abundance = group.abundances[si] / total_abundance;
            true_abundances.insert(species.clone(), norm_abundance);
            all_species.push((species.clone(), norm_abundance, gi));
        }
    }

    // Generate reads
    let mut alignments: ReadAlignments = HashMap::new();
    let weights: Vec<f64> = all_species.iter().map(|(_, a, _)| *a).collect();

    for read_idx in 0..config.n_reads {
        let source_idx = rng.weighted_choice(&weights);
        let (ref source_species, _, group_idx) = all_species[source_idx];

        let read_name = format!("read_{:06}", read_idx);
        let mut read_alns = Vec::new();

        // Primary alignment to true source
        read_alns.push(PafRecord {
            query_name: read_name.clone(),
            query_len: 1000,
            target_name: source_species.clone(),
            target_len: 800,
            matches: 950,
            block_len: 1000,
            mapq: 60,
            alignment_score: Some(config.true_alignment_score),
        });

        // Cross-mappings to other species in the same group
        let group = &groups[group_idx];
        if group.species.len() > 1 {
            for other_species in &group.species {
                if other_species == source_species {
                    continue;
                }
                if rng.next_f64() < config.cross_mapping_rate {
                    read_alns.push(PafRecord {
                        query_name: read_name.clone(),
                        query_len: 1000,
                        target_name: other_species.clone(),
                        target_len: 800,
                        matches: 920,
                        block_len: 1000,
                        mapq: 50,
                        alignment_score: Some(config.cross_alignment_score),
                    });
                }
            }
        }

        alignments.insert(read_name, read_alns);
    }

    SimulatedCommunity {
        true_abundances,
        alignments,
        n_reads: config.n_reads,
    }
}

/// Compute L1 error between estimated and true abundances.
/// L1 = sum of |estimated - true| for all taxa.
pub fn l1_error(estimated: &HashMap<String, f64>, truth: &HashMap<String, f64>) -> f64 {
    let all_taxa: std::collections::HashSet<&String> =
        estimated.keys().chain(truth.keys()).collect();

    all_taxa
        .iter()
        .map(|t| {
            let est = estimated.get(*t).copied().unwrap_or(0.0);
            let tru = truth.get(*t).copied().unwrap_or(0.0);
            (est - tru).abs()
        })
        .sum()
}

/// Compute Bray-Curtis dissimilarity between two abundance profiles.
/// BC = sum|est_i - true_i| / sum(est_i + true_i)
pub fn bray_curtis(estimated: &HashMap<String, f64>, truth: &HashMap<String, f64>) -> f64 {
    let all_taxa: std::collections::HashSet<&String> =
        estimated.keys().chain(truth.keys()).collect();

    let mut numerator = 0.0;
    let mut denominator = 0.0;

    for t in &all_taxa {
        let est = estimated.get(*t).copied().unwrap_or(0.0);
        let tru = truth.get(*t).copied().unwrap_or(0.0);
        numerator += (est - tru).abs();
        denominator += est + tru;
    }

    if denominator == 0.0 {
        return 0.0;
    }

    numerator / denominator
}

/// Run a complete simulation experiment and print comparison results.
pub fn run_simulation_experiment(
    name: &str,
    groups: &[SpeciesGroup],
    config: &SimConfig,
    em_config: &crate::em::EmConfig,
) {
    println!("\n{}", "=".repeat(70));
    println!("Experiment: {}", name);
    println!("{}", "=".repeat(70));

    let community = simulate_community(groups, config);

    // Print true composition
    println!("\nTrue community composition:");
    let mut sorted_truth: Vec<_> = community.true_abundances.iter().collect();
    sorted_truth.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());
    for (taxon, abundance) in &sorted_truth {
        println!(
            "  {:<30} {:.4} ({:.1}%)",
            taxon,
            abundance,
            *abundance * 100.0
        );
    }

    // Count multi-mapped reads
    let multi_mapped = community
        .alignments
        .values()
        .filter(|alns| alns.len() > 1)
        .count();
    println!(
        "\nReads: {} total, {} multi-mapped ({:.1}%)",
        community.n_reads,
        multi_mapped,
        multi_mapped as f64 / community.n_reads as f64 * 100.0
    );

    // Run naive counting
    let naive = crate::em::naive_count(&community.alignments);

    // Run EM
    let em_result = crate::em::run_em(&community.alignments, em_config);

    // Compare results
    println!(
        "\n{:<30} {:>10} {:>10} {:>10}",
        "Taxon", "Truth", "Naive", "EM"
    );
    println!("{}", "-".repeat(62));
    for (taxon, true_abundance) in &sorted_truth {
        let naive_est = naive.get(*taxon).copied().unwrap_or(0.0);
        let em_est = em_result.abundances.get(*taxon).copied().unwrap_or(0.0);
        println!(
            "{:<30} {:>9.4} {:>9.4} {:>9.4}",
            taxon, true_abundance, naive_est, em_est
        );
    }

    // Check for false positives (taxa estimated but not in truth)
    for taxon in em_result.abundances.keys() {
        if !community.true_abundances.contains_key(taxon) {
            let em_est = em_result.abundances[taxon];
            println!(
                "{:<30} {:>9.4} {:>9.4} {:>9.4}  <-- false positive",
                taxon,
                0.0,
                naive.get(taxon).copied().unwrap_or(0.0),
                em_est
            );
        }
    }

    // Error metrics
    let naive_l1 = l1_error(&naive, &community.true_abundances);
    let em_l1 = l1_error(&em_result.abundances, &community.true_abundances);
    let naive_bc = bray_curtis(&naive, &community.true_abundances);
    let em_bc = bray_curtis(&em_result.abundances, &community.true_abundances);

    println!("\nError metrics:");
    println!("  {:<25} {:>10} {:>10}", "", "Naive", "EM");
    println!("  {:<25} {:>10.6} {:>10.6}", "L1 error", naive_l1, em_l1);
    println!("  {:<25} {:>10.6} {:>10.6}", "Bray-Curtis", naive_bc, em_bc);

    let l1_improvement = if naive_l1 > 0.0 {
        (1.0 - em_l1 / naive_l1) * 100.0
    } else {
        0.0
    };
    let bc_improvement = if naive_bc > 0.0 {
        (1.0 - em_bc / naive_bc) * 100.0
    } else {
        0.0
    };

    println!(
        "\n  EM improvement: L1 {:.1}% better, Bray-Curtis {:.1}% better",
        l1_improvement, bc_improvement
    );
    println!("  EM converged in {} iterations", em_result.iterations);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_l1_error() {
        let mut a = HashMap::new();
        let mut b = HashMap::new();
        a.insert("sp_A".to_string(), 0.7);
        a.insert("sp_B".to_string(), 0.3);
        b.insert("sp_A".to_string(), 0.6);
        b.insert("sp_B".to_string(), 0.4);

        let err = l1_error(&a, &b);
        assert!((err - 0.2).abs() < 1e-10);
    }

    #[test]
    fn test_bray_curtis_identical() {
        let mut a = HashMap::new();
        a.insert("sp_A".to_string(), 0.5);
        a.insert("sp_B".to_string(), 0.5);
        let bc = bray_curtis(&a, &a);
        assert!(bc.abs() < 1e-10);
    }

    #[test]
    fn test_simulation_generates_reads() {
        let groups = vec![SpeciesGroup {
            species: vec!["sp_A".to_string(), "sp_B".to_string()],
            abundances: vec![0.7, 0.3],
        }];
        let config = SimConfig {
            n_reads: 100,
            ..Default::default()
        };
        let community = simulate_community(&groups, &config);
        assert_eq!(community.alignments.len(), 100);
        assert!((community.true_abundances["sp_A"] - 0.7).abs() < 1e-10);
    }
}
