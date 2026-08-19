//! Sequencing platform presets for EMU-ITS.
//!
//! Different long-read sequencing platforms produce reads with distinct
//! error profiles, quality score distributions, and alignment characteristics.
//! These presets bundle optimized parameters for each platform so users
//! don't need to tune EM internals manually.
//!
//! # Supported platforms
//!
//! - **ONT R10.4.1**: Oxford Nanopore R10.4.1 chemistry (current generation).
//!   ~Q20 median accuracy, moderate homopolymer errors. Default preset.
//! - **ONT R9.4.1**: Oxford Nanopore R9.4.1 chemistry (legacy).
//!   ~Q12-Q15 accuracy, higher error rate requires relaxed thresholds.
//! - **PacBio HiFi**: PacBio Circular Consensus Sequencing (HiFi/Revio).
//!   ~Q30+ accuracy, very low error rate allows stringent thresholds.
//! - **ONT Duplex**: Oxford Nanopore Duplex reads.
//!   ~Q30 accuracy from complementary strand consensus.

use std::fmt;

/// Sequencing platform preset identifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Platform {
    /// Oxford Nanopore R10.4.1 (current generation, ~Q20)
    OntR10,
    /// Oxford Nanopore R9.4.1 (legacy, ~Q12-Q15)
    OntR9,
    /// PacBio HiFi / Revio (~Q30+)
    PacBioHifi,
    /// Oxford Nanopore Duplex (~Q30)
    OntDuplex,
}

impl fmt::Display for Platform {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Platform::OntR10 => write!(f, "ont-r10"),
            Platform::OntR9 => write!(f, "ont-r9"),
            Platform::PacBioHifi => write!(f, "pacbio-hifi"),
            Platform::OntDuplex => write!(f, "ont-duplex"),
        }
    }
}

impl Platform {
    /// Parse a preset name string into a Platform variant.
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "ont-r10" | "ont_r10" | "r10" | "ont-r10.4.1" | "r10.4.1" => Some(Platform::OntR10),
            "ont-r9" | "ont_r9" | "r9" | "ont-r9.4.1" | "r9.4.1" => Some(Platform::OntR9),
            "pacbio-hifi" | "pacbio_hifi" | "hifi" | "pacbio" | "revio" => {
                Some(Platform::PacBioHifi)
            }
            "ont-duplex" | "ont_duplex" | "duplex" => Some(Platform::OntDuplex),
            _ => None,
        }
    }

    /// List all available preset names for help text.
    pub fn available() -> &'static str {
        "ont-r10, ont-r9, pacbio-hifi, ont-duplex"
    }

    /// Short description for each platform.
    pub fn description(&self) -> &'static str {
        match self {
            Platform::OntR10 => "Oxford Nanopore R10.4.1 (~Q20, current chemistry)",
            Platform::OntR9 => "Oxford Nanopore R9.4.1 (~Q12-15, legacy chemistry)",
            Platform::PacBioHifi => "PacBio HiFi/Revio (~Q30+, high accuracy)",
            Platform::OntDuplex => "Oxford Nanopore Duplex (~Q30, duplex consensus)",
        }
    }
}

/// Complete parameter set for a sequencing platform.
///
/// Bundles alignment filtering, EM algorithm, and minimap2 parameters
/// optimized for each platform's error profile.
#[derive(Debug, Clone)]
pub struct PresetParams {
    /// Source platform
    pub platform: Platform,

    // ── Alignment filtering ──
    /// Minimum alignment identity to retain a hit (0.0-1.0)
    pub min_identity: f64,
    /// Minimum mapping quality
    pub min_mapq: u8,

    // ── EM parameters ──
    /// Temperature for score-to-likelihood conversion.
    /// Lower = more sensitive to score differences (better for high-accuracy reads).
    /// Higher = more tolerant of score noise (better for noisy reads).
    pub temperature: f64,
    /// Maximum EM iterations
    pub max_iterations: usize,
    /// Convergence threshold
    pub convergence_threshold: f64,
    /// Minimum abundance to retain a taxon between iterations
    pub min_abundance: f64,

    // ── Recommended minimap2 parameters ──
    /// Recommended minimap2 preset flag (e.g., "map-ont", "map-hifi")
    pub minimap2_preset: &'static str,
    /// Recommended max secondary alignments (-N)
    pub minimap2_secondary_n: u32,
    /// Recommended secondary score ratio (-p)
    pub minimap2_secondary_p: f64,
}

impl PresetParams {
    /// Get optimized parameters for a given platform.
    pub fn for_platform(platform: Platform) -> Self {
        match platform {
            Platform::OntR10 => PresetParams {
                platform,
                // R10.4.1: ~Q20 median, decent accuracy
                min_identity: 0.80,
                min_mapq: 0,
                temperature: 0.5,
                // Raised from 100, which was truncating EM before convergence on
                // both published datasets: the simulated 21-species community
                // needs 232 iterations and the ATCC ONT mock 175. Truncation was
                // conservative rather than wrong -- on the mock, estimates at 100
                // iterations already matched the converged values to within 0.007
                // percentage points -- but on the simulated community it left L1
                // at 7.48% where convergence gives 7.30%.
                max_iterations: 500,
                convergence_threshold: 1e-6,
                min_abundance: 1e-7,
                minimap2_preset: "map-ont",
                minimap2_secondary_n: 10,
                minimap2_secondary_p: 0.9,
            },
            Platform::OntR9 => PresetParams {
                platform,
                // R9.4.1: ~Q12-Q15, higher error rate
                // Lower identity threshold to retain more alignments
                // Higher temperature to be more tolerant of score noise
                min_identity: 0.70,
                min_mapq: 0,
                temperature: 0.8,
                // Raised for the same reason as ont-r10. R9 admits more
                // candidates at a lower identity threshold and uses a higher
                // temperature, both of which slow convergence, so it needs at
                // least as much headroom. Not measured directly.
                max_iterations: 500,
                convergence_threshold: 1e-6,
                min_abundance: 1e-7,
                minimap2_preset: "map-ont",
                minimap2_secondary_n: 15,
                minimap2_secondary_p: 0.85,
            },
            Platform::PacBioHifi => PresetParams {
                platform,
                // HiFi: ~Q30+, very accurate reads.
                // Strict identity threshold since reads are high quality.
                // Low temperature = more sensitive to score differences.
                //
                // min_mapq MUST stay 0. minimap2 assigns MAPQ 0 to reads with
                // multiple equally-scoring alignments -- precisely the ambiguous
                // reads EM exists to resolve. An earlier value of 5 discarded
                // them before EM saw them, collapsing the candidate set from 187
                // taxa to 31 and driving the EM improvement over naive counting
                // to exactly zero (L1 6.43% vs 5.03% once corrected). See
                // scripts/sweep_hifi_preset.py.
                min_identity: 0.95,
                min_mapq: 0,
                temperature: 0.15,
                // Raised from 100. With min_mapq corrected to 0 the candidate
                // set grows to ~170 taxa and EM needs 152 iterations to reach
                // convergence_threshold on the simulated HiFi community; at 100
                // it terminated early with delta 3.67e-6.
                //
                // Set well above the observed 152 rather than just above it.
                // The cost of a high ceiling is bounded -- EM stops as soon as
                // it converges -- whereas the cost of a low one is a silently
                // truncated estimate, which is exactly the failure this preset
                // shipped with. Iteration count is also sensitive to
                // temperature: at temperature 0.5 the same data needs 376.
                // Real communities are more complex than this simulation, so
                // leave margin.
                max_iterations: 500,
                convergence_threshold: 1e-7,
                min_abundance: 1e-8,
                minimap2_preset: "map-hifi",
                minimap2_secondary_n: 10,
                minimap2_secondary_p: 0.95,
            },
            Platform::OntDuplex => PresetParams {
                platform,
                // Duplex: ~Q30, comparable to HiFi but ONT error profile
                // (still has some homopolymer issues unlike HiFi).
                //
                // min_mapq set to 0 for the same reason as pacbio-hifi: MAPQ 0
                // marks the multi-mapping reads EM is designed to resolve.
                // Not benchmarked directly -- corrected by analogy with the
                // PacBio HiFi result, since both are high-accuracy presets that
                // shared the same mistuning.
                min_identity: 0.90,
                min_mapq: 0,
                temperature: 0.2,
                // Raised for the same reason as the other presets. Not measured
                // directly; corrected by analogy.
                max_iterations: 500,
                convergence_threshold: 1e-7,
                min_abundance: 1e-8,
                minimap2_preset: "map-ont",
                minimap2_secondary_n: 10,
                minimap2_secondary_p: 0.93,
            },
        }
    }

    /// Generate the recommended minimap2 command for this preset.
    pub fn minimap2_cmd(&self, db_path: &str, reads_path: &str) -> String {
        format!(
            "minimap2 -c --secondary=yes -N {} -p {} -x {} {} {}",
            self.minimap2_secondary_n,
            self.minimap2_secondary_p,
            self.minimap2_preset,
            db_path,
            reads_path
        )
    }

    /// Print a summary of the preset parameters.
    pub fn print_summary(&self) {
        println!(
            "Platform preset: {} ({})",
            self.platform,
            self.platform.description()
        );
        println!("  Alignment filters:");
        println!("    min_identity:    {:.2}", self.min_identity);
        println!("    min_mapq:        {}", self.min_mapq);
        println!("  EM parameters:");
        println!("    temperature:     {:.2}", self.temperature);
        println!("    max_iterations:  {}", self.max_iterations);
        println!("    convergence:     {:.0e}", self.convergence_threshold);
        println!("    min_abundance:   {:.0e}", self.min_abundance);
        println!("  Recommended minimap2:");
        println!(
            "    minimap2 -c --secondary=yes -N {} -p {} -x {} <db> <reads>",
            self.minimap2_secondary_n, self.minimap2_secondary_p, self.minimap2_preset
        );
    }
}

/// Convert preset parameters into an EmConfig for the EM engine.
impl From<&PresetParams> for crate::em::EmConfig {
    fn from(preset: &PresetParams) -> Self {
        crate::em::EmConfig {
            max_iterations: preset.max_iterations,
            convergence_threshold: preset.convergence_threshold,
            min_abundance: preset.min_abundance,
            temperature: preset.temperature,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_platform_parsing() {
        assert_eq!(Platform::from_str("ont-r10"), Some(Platform::OntR10));
        assert_eq!(Platform::from_str("r10"), Some(Platform::OntR10));
        assert_eq!(Platform::from_str("R10.4.1"), Some(Platform::OntR10));
        assert_eq!(Platform::from_str("ont-r9"), Some(Platform::OntR9));
        assert_eq!(Platform::from_str("hifi"), Some(Platform::PacBioHifi));
        assert_eq!(
            Platform::from_str("pacbio-hifi"),
            Some(Platform::PacBioHifi)
        );
        assert_eq!(Platform::from_str("revio"), Some(Platform::PacBioHifi));
        assert_eq!(Platform::from_str("duplex"), Some(Platform::OntDuplex));
        assert_eq!(Platform::from_str("illumina"), None);
    }

    #[test]
    fn test_preset_params() {
        let ont = PresetParams::for_platform(Platform::OntR10);
        let hifi = PresetParams::for_platform(Platform::PacBioHifi);

        // HiFi should have stricter identity threshold
        assert!(hifi.min_identity > ont.min_identity);
        // HiFi should have lower temperature (more score-sensitive)
        assert!(hifi.temperature < ont.temperature);
    }

    /// Regression guard. `max_iterations` is a safety valve, not a stopping
    /// rule -- `convergence_threshold` should be what halts EM. Observed
    /// requirements: 232 iterations (simulated community, ont-r10), 175 (ATCC
    /// ONT mock, ont-r10), 376 (HiFi at temperature 0.5).
    #[test]
    fn test_iteration_ceilings_exceed_observed_requirements() {
        const OBSERVED_MAX: usize = 376;
        for platform in [
            Platform::OntR10,
            Platform::OntR9,
            Platform::PacBioHifi,
            Platform::OntDuplex,
        ] {
            let p = PresetParams::for_platform(platform);
            assert!(
                p.max_iterations > OBSERVED_MAX,
                "{platform} allows only {} iterations; EM has been observed to \
                 need {OBSERVED_MAX}, so this ceiling would truncate before \
                 convergence and silently return a conservative estimate",
                p.max_iterations
            );
        }
    }

    /// Regression guard. Filtering on MAPQ removes multi-mapping reads, which
    /// are the reads EM is for. No preset may filter them out.
    #[test]
    fn test_no_preset_filters_multimapping_reads() {
        for platform in [
            Platform::OntR10,
            Platform::OntR9,
            Platform::PacBioHifi,
            Platform::OntDuplex,
        ] {
            let p = PresetParams::for_platform(platform);
            assert_eq!(
                p.min_mapq, 0,
                "{platform} sets min_mapq={}; MAPQ 0 marks the ambiguous \
                 alignments EM is designed to resolve, so filtering on it \
                 discards exactly the signal EMITS depends on",
                p.min_mapq
            );
        }
    }

    #[test]
    fn test_r9_more_permissive_than_r10() {
        let r9 = PresetParams::for_platform(Platform::OntR9);
        let r10 = PresetParams::for_platform(Platform::OntR10);

        assert!(r9.min_identity < r10.min_identity);
        assert!(r9.temperature > r10.temperature);
    }

    #[test]
    fn test_minimap2_cmd() {
        let preset = PresetParams::for_platform(Platform::OntR10);
        let cmd = preset.minimap2_cmd("unite.fasta", "reads.fastq");
        assert!(cmd.contains("map-ont"));
        assert!(cmd.contains("-N 10"));
    }
}
