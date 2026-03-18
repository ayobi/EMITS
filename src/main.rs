//! EMU-ITS: Expectation-Maximization abundance estimation for fungal ITS
//! communities from long-read sequencing.
//!
//! Built for ONT and PacBio long reads with platform-specific presets
//! that optimize alignment filtering and EM parameters for each
//! sequencing technology's error profile.

mod em;
mod output;
mod paf;
mod preset;
mod sim;
mod taxonomy;

use anyhow::Result;
use clap::{Parser, Subcommand};
use std::fs::File;
use std::io::{BufReader, BufWriter};

#[derive(Parser)]
#[command(name = "emu-its")]
#[command(about = "EM-based abundance estimation for fungal ITS from long reads")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run EM abundance estimation on a PAF alignment file
    Run {
        /// Input PAF file from minimap2 (use - for stdin)
        #[arg(short, long)]
        input: String,

        /// Output TSV file for abundance estimates
        #[arg(short, long, default_value = "abundances.tsv")]
        output: String,

        /// Sequencing platform preset [ont-r10, ont-r9, pacbio-hifi, ont-duplex]
        /// Sets optimized defaults for alignment filtering and EM parameters.
        /// Individual parameters below override the preset if specified.
        #[arg(short, long, default_value = "ont-r10")]
        preset: String,

        /// Minimum alignment identity (0.0-1.0) [overrides preset]
        #[arg(long)]
        min_identity: Option<f64>,

        /// Minimum mapping quality [overrides preset]
        #[arg(long)]
        min_mapq: Option<u8>,

        /// Maximum EM iterations [overrides preset]
        #[arg(long)]
        max_iter: Option<usize>,

        /// EM convergence threshold [overrides preset]
        #[arg(long)]
        threshold: Option<f64>,

        /// EM temperature parameter [overrides preset]
        #[arg(long)]
        temperature: Option<f64>,

        /// Also output naive count-based estimates for comparison
        #[arg(long)]
        compare: bool,

        /// Taxonomic rank for aggregated output [species, genus, family, order, class, phylum]
        /// Collapses multiple UNITE accessions into one row per taxon at this rank.
        /// Raw per-accession output is always written; aggregated output is additional.
        #[arg(long, default_value = "species")]
        rank: String,
    },

    /// Show details of a sequencing platform preset
    Preset {
        /// Platform name [ont-r10, ont-r9, pacbio-hifi, ont-duplex]
        #[arg(default_value = "ont-r10")]
        platform: String,

        /// Show recommended minimap2 command for this preset
        #[arg(long)]
        minimap2_cmd: bool,

        /// UNITE database path (for minimap2 command output)
        #[arg(long, default_value = "unite_db.fasta")]
        db: String,

        /// Reads file path (for minimap2 command output)
        #[arg(long, default_value = "reads.fastq")]
        reads: String,
    },

    /// Run simulation experiments to validate EM vs naive counting
    Simulate {
        /// Number of reads to simulate
        #[arg(short, long, default_value = "5000")]
        n_reads: usize,

        /// Cross-mapping rate (0.0-1.0)
        #[arg(short, long, default_value = "0.3")]
        cross_rate: f64,

        /// Output comparison TSV
        #[arg(short, long)]
        output: Option<String>,

        /// Run all predefined experiments
        #[arg(long)]
        all: bool,
    },
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Run {
            input,
            output,
            preset,
            min_identity,
            min_mapq,
            max_iter,
            threshold,
            temperature,
            compare,
            rank,
        } => {
            // Resolve taxonomic rank
            let tax_rank = taxonomy::TaxRank::from_str(&rank).unwrap_or_else(|| {
                eprintln!(
                    "Unknown rank '{}'. Available: {}",
                    rank,
                    taxonomy::TaxRank::available()
                );
                eprintln!("Using default: species");
                taxonomy::TaxRank::Species
            });
            // Resolve platform preset
            let platform = preset::Platform::from_str(&preset).unwrap_or_else(|| {
                eprintln!(
                    "Unknown preset '{}'. Available: {}",
                    preset,
                    preset::Platform::available()
                );
                eprintln!("Using default: ont-r10");
                preset::Platform::OntR10
            });
            let params = preset::PresetParams::for_platform(platform);

            // Apply preset defaults, with CLI overrides taking precedence
            let eff_min_identity = min_identity.unwrap_or(params.min_identity);
            let eff_min_mapq = min_mapq.unwrap_or(params.min_mapq);
            let eff_max_iter = max_iter.unwrap_or(params.max_iterations);
            let eff_threshold = threshold.unwrap_or(params.convergence_threshold);
            let eff_temperature = temperature.unwrap_or(params.temperature);

            // Print active configuration
            println!("\nEMU-ITS v{}", env!("CARGO_PKG_VERSION"));
            println!("{}", "=".repeat(60));
            println!("Platform:       {} ({})", platform, platform.description());
            println!(
                "Min identity:   {:.2}{}",
                eff_min_identity,
                if min_identity.is_some() {
                    " (override)"
                } else {
                    ""
                }
            );
            println!(
                "Min mapq:       {}{}",
                eff_min_mapq,
                if min_mapq.is_some() {
                    " (override)"
                } else {
                    ""
                }
            );
            println!(
                "Temperature:    {:.2}{}",
                eff_temperature,
                if temperature.is_some() {
                    " (override)"
                } else {
                    ""
                }
            );
            println!(
                "Max iterations: {}{}",
                eff_max_iter,
                if max_iter.is_some() {
                    " (override)"
                } else {
                    ""
                }
            );
            println!("{}", "=".repeat(60));

            // Parse PAF alignments
            let reader: Box<dyn std::io::BufRead> = if input == "-" {
                Box::new(BufReader::new(std::io::stdin()))
            } else {
                Box::new(BufReader::new(
                    File::open(&input)
                        .map_err(|e| anyhow::anyhow!("Cannot open {}: {}", input, e))?,
                ))
            };

            let alignments = paf::parse_paf(reader, eff_min_identity, eff_min_mapq)?;
            let total_reads = alignments.len();

            if total_reads == 0 {
                anyhow::bail!(
                    "No alignments passed filters. Check input file and filter parameters."
                );
            }

            // Run EM with preset-derived config
            let config = em::EmConfig {
                max_iterations: eff_max_iter,
                convergence_threshold: eff_threshold,
                min_abundance: params.min_abundance,
                temperature: eff_temperature,
            };
            let result = em::run_em(&alignments, &config);

            // Write raw per-accession output
            let mut writer = BufWriter::new(File::create(&output)?);
            output::write_abundance_tsv(&mut writer, &result.abundances, total_reads)?;
            log::info!("Wrote raw EM abundance estimates to {}", output);

            // Write aggregated output
            let em_agg = taxonomy::aggregate_abundances(&result.abundances, tax_rank);
            let agg_path = output.replace(".tsv", &format!("_{}.tsv", rank));
            let mut agg_writer = BufWriter::new(File::create(&agg_path)?);
            output::write_aggregated_tsv(&mut agg_writer, &em_agg, total_reads)?;
            log::info!(
                "Wrote aggregated EM estimates ({}-level, {} taxa) to {}",
                rank,
                em_agg.abundances.len(),
                agg_path
            );

            if compare {
                let naive = em::naive_count(&alignments);

                // Raw comparison
                let comparison_path = output.replace(".tsv", "_comparison.tsv");
                let mut cmp_writer = BufWriter::new(File::create(&comparison_path)?);
                output::write_comparison_tsv(
                    &mut cmp_writer,
                    &result.abundances,
                    &naive,
                    &result.abundances,
                )?;
                log::info!("Wrote raw comparison to {}", comparison_path);

                // Aggregated comparison
                let naive_agg = taxonomy::aggregate_abundances(&naive, tax_rank);
                let agg_cmp_path = output.replace(".tsv", &format!("_{}_comparison.tsv", rank));
                let mut agg_cmp_writer = BufWriter::new(File::create(&agg_cmp_path)?);
                output::write_aggregated_comparison_tsv(
                    &mut agg_cmp_writer,
                    &em_agg,
                    &naive_agg,
                    total_reads,
                )?;
                log::info!("Wrote aggregated comparison to {}", agg_cmp_path);

                // Print aggregated summary with comparison
                println!("\nResults");
                println!("{}", "=".repeat(60));
                println!("Total reads:      {}", total_reads);
                println!("Raw taxa:         {}", result.abundances.len());
                println!(
                    "Aggregated taxa:  {} ({}-level)",
                    em_agg.abundances.len(),
                    rank
                );
                println!("EM iterations:    {}", result.iterations);
                println!("Final delta:      {:.2e}", result.final_delta);

                output::print_aggregated_summary(&em_agg, Some(&naive_agg), 15);
            } else {
                // Print aggregated summary without comparison
                println!("\nResults");
                println!("{}", "=".repeat(60));
                println!("Total reads:      {}", total_reads);
                println!("Raw taxa:         {}", result.abundances.len());
                println!(
                    "Aggregated taxa:  {} ({}-level)",
                    em_agg.abundances.len(),
                    rank
                );
                println!("EM iterations:    {}", result.iterations);
                println!("Final delta:      {:.2e}", result.final_delta);

                output::print_aggregated_summary(&em_agg, None, 15);
            }
        }

        Commands::Preset {
            platform,
            minimap2_cmd,
            db,
            reads,
        } => {
            let plat = preset::Platform::from_str(&platform).unwrap_or_else(|| {
                eprintln!(
                    "Unknown platform '{}'. Available: {}",
                    platform,
                    preset::Platform::available()
                );
                std::process::exit(1);
            });
            let params = preset::PresetParams::for_platform(plat);

            if minimap2_cmd {
                println!("{}", params.minimap2_cmd(&db, &reads));
            } else {
                println!();
                params.print_summary();
                println!();
            }
        }

        Commands::Simulate {
            n_reads,
            cross_rate,
            output,
            all,
        } => {
            let em_config = em::EmConfig::default();

            if all {
                run_all_experiments(n_reads, &em_config);
            } else {
                // Default experiment: simple two-group community
                let groups = vec![
                    sim::SpeciesGroup {
                        species: vec![
                            "Fusarium_oxysporum".to_string(),
                            "Fusarium_solani".to_string(),
                            "Fusarium_graminearum".to_string(),
                        ],
                        abundances: vec![0.40, 0.15, 0.05],
                    },
                    sim::SpeciesGroup {
                        species: vec![
                            "Trichoderma_harzianum".to_string(),
                            "Trichoderma_viride".to_string(),
                        ],
                        abundances: vec![0.25, 0.10],
                    },
                    sim::SpeciesGroup {
                        species: vec!["Saccharomyces_cerevisiae".to_string()],
                        abundances: vec![0.05],
                    },
                ];

                let sim_config = sim::SimConfig {
                    n_reads,
                    cross_mapping_rate: cross_rate,
                    ..Default::default()
                };

                sim::run_simulation_experiment(
                    "Default ITS community",
                    &groups,
                    &sim_config,
                    &em_config,
                );

                // Write output if requested
                if let Some(output_path) = output {
                    let community = sim::simulate_community(&groups, &sim_config);
                    let naive = em::naive_count(&community.alignments);
                    let em_result = em::run_em(&community.alignments, &em_config);
                    let mut writer = BufWriter::new(File::create(&output_path)?);
                    output::write_comparison_tsv(
                        &mut writer,
                        &community.true_abundances,
                        &naive,
                        &em_result.abundances,
                    )?;
                    log::info!("Wrote comparison to {}", output_path);
                }
            }
        }
    }

    Ok(())
}

/// Run a battery of simulation experiments testing different scenarios.
fn run_all_experiments(n_reads: usize, em_config: &em::EmConfig) {
    println!("╔══════════════════════════════════════════════════════════════════════╗");
    println!("║              EMU-ITS Proof of Concept: EM vs Naive Counting         ║");
    println!("║              Simulated ITS Community Experiments                    ║");
    println!("╚══════════════════════════════════════════════════════════════════════╝");

    // Experiment 1: Simple community with closely related species
    {
        let groups = vec![
            sim::SpeciesGroup {
                species: vec![
                    "Fusarium_oxysporum".to_string(),
                    "Fusarium_solani".to_string(),
                    "Fusarium_graminearum".to_string(),
                ],
                abundances: vec![0.40, 0.15, 0.05],
            },
            sim::SpeciesGroup {
                species: vec![
                    "Trichoderma_harzianum".to_string(),
                    "Trichoderma_viride".to_string(),
                ],
                abundances: vec![0.25, 0.10],
            },
            sim::SpeciesGroup {
                species: vec!["Saccharomyces_cerevisiae".to_string()],
                abundances: vec![0.05],
            },
        ];

        let config = sim::SimConfig {
            n_reads,
            cross_mapping_rate: 0.3,
            true_alignment_score: 1900,
            cross_alignment_score: 1850,
            seed: 42,
        };

        sim::run_simulation_experiment(
            "1. Moderate complexity (6 spp, 30% cross-mapping)",
            &groups,
            &config,
            em_config,
        );
    }

    // Experiment 2: High cross-mapping rate (very similar ITS sequences)
    {
        let groups = vec![
            sim::SpeciesGroup {
                species: vec![
                    "Aspergillus_niger".to_string(),
                    "Aspergillus_tubingensis".to_string(),
                    "Aspergillus_brasiliensis".to_string(),
                    "Aspergillus_carbonarius".to_string(),
                ],
                abundances: vec![0.35, 0.25, 0.10, 0.05],
            },
            sim::SpeciesGroup {
                species: vec![
                    "Penicillium_chrysogenum".to_string(),
                    "Penicillium_rubens".to_string(),
                ],
                abundances: vec![0.15, 0.10],
            },
        ];

        let config = sim::SimConfig {
            n_reads,
            cross_mapping_rate: 0.6, // High cross-mapping
            true_alignment_score: 1900,
            cross_alignment_score: 1870, // Very similar scores
            seed: 123,
        };

        sim::run_simulation_experiment(
            "2. High similarity (Aspergillus section Nigri, 60% cross-mapping)",
            &groups,
            &config,
            em_config,
        );
    }

    // Experiment 3: Rare species detection
    {
        let groups = vec![
            sim::SpeciesGroup {
                species: vec![
                    "Trichoderma_harzianum".to_string(),
                    "Trichoderma_atroviride".to_string(),
                    "Trichoderma_virens".to_string(),
                ],
                abundances: vec![0.60, 0.01, 0.005], // Very rare species
            },
            sim::SpeciesGroup {
                species: vec![
                    "Metarhizium_anisopliae".to_string(),
                    "Metarhizium_robertsii".to_string(),
                ],
                abundances: vec![0.20, 0.02],
            },
            sim::SpeciesGroup {
                species: vec!["Beauveria_bassiana".to_string()],
                abundances: vec![0.105],
            },
        ];

        let config = sim::SimConfig {
            n_reads,
            cross_mapping_rate: 0.4,
            true_alignment_score: 1900,
            cross_alignment_score: 1840,
            seed: 456,
        };

        sim::run_simulation_experiment(
            "3. Rare species detection (biocontrol agents, includes 0.5-2% taxa)",
            &groups,
            &config,
            em_config,
        );
    }

    // Experiment 4: Even community (all species at similar abundance)
    {
        let groups = vec![
            sim::SpeciesGroup {
                species: vec![
                    "Candida_albicans".to_string(),
                    "Candida_tropicalis".to_string(),
                    "Candida_parapsilosis".to_string(),
                ],
                abundances: vec![0.20, 0.18, 0.15],
            },
            sim::SpeciesGroup {
                species: vec![
                    "Pichia_kudriavzevii".to_string(),
                    "Pichia_membranifaciens".to_string(),
                ],
                abundances: vec![0.17, 0.14],
            },
            sim::SpeciesGroup {
                species: vec!["Debaryomyces_hansenii".to_string()],
                abundances: vec![0.16],
            },
        ];

        let config = sim::SimConfig {
            n_reads,
            cross_mapping_rate: 0.35,
            true_alignment_score: 1900,
            cross_alignment_score: 1855,
            seed: 789,
        };

        sim::run_simulation_experiment(
            "4. Even community (6 yeasts, similar abundances)",
            &groups,
            &config,
            em_config,
        );
    }

    // Experiment 5: Varying cross-mapping rates (sensitivity analysis)
    {
        println!("\n{}", "=".repeat(70));
        println!("Experiment 5: Sensitivity to cross-mapping rate");
        println!("{}", "=".repeat(70));

        let groups = vec![
            sim::SpeciesGroup {
                species: vec![
                    "Fusarium_oxysporum".to_string(),
                    "Fusarium_solani".to_string(),
                ],
                abundances: vec![0.60, 0.20],
            },
            sim::SpeciesGroup {
                species: vec!["Trichoderma_harzianum".to_string()],
                abundances: vec![0.20],
            },
        ];

        println!(
            "\n{:<20} {:>12} {:>12} {:>12} {:>12}",
            "Cross-map rate", "Naive L1", "EM L1", "Naive BC", "EM BC"
        );
        println!("{}", "-".repeat(70));

        for rate_pct in [10, 20, 30, 40, 50, 60, 70, 80] {
            let rate = rate_pct as f64 / 100.0;
            let config = sim::SimConfig {
                n_reads,
                cross_mapping_rate: rate,
                true_alignment_score: 1900,
                cross_alignment_score: 1860,
                seed: 42 + rate_pct,
            };

            let community = sim::simulate_community(&groups, &config);
            let naive = em::naive_count(&community.alignments);
            let em_result = em::run_em(&community.alignments, em_config);

            let naive_l1 = sim::l1_error(&naive, &community.true_abundances);
            let em_l1 = sim::l1_error(&em_result.abundances, &community.true_abundances);
            let naive_bc = sim::bray_curtis(&naive, &community.true_abundances);
            let em_bc = sim::bray_curtis(&em_result.abundances, &community.true_abundances);

            println!(
                "{:<20} {:>12.6} {:>12.6} {:>12.6} {:>12.6}",
                format!("{}%", rate_pct),
                naive_l1,
                em_l1,
                naive_bc,
                em_bc
            );
        }
    }

    println!("\n{}", "=".repeat(70));
    println!("All experiments complete.");
    println!("{}", "=".repeat(70));
}
