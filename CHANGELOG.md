# Changelog

All notable changes to EMITS will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-08-25

### Fixed

- **`pacbio-hifi` and `ont-duplex` presets set `min_mapq` to 5, which discarded the
  multi-mapping alignments EM exists to resolve.** minimap2 assigns MAPQ 0 to reads with
  several equally scoring alignments; filtering them out removed 25.9% of alignments,
  reduced the candidate set from ~170 taxa to 31, and eliminated any improvement over
  naive best-hit counting. Corrected to 0. On a simulated PacBio HiFi community this
  lowers aggregate L1 error from 6.43% to 5.02% and restores a 20.8% improvement over
  naive counting.

  **Users who ran EMITS 0.1.0 with `--preset pacbio-hifi` or `--preset ont-duplex` should
  re-run their analyses.** ONT R10 and R9 results are unaffected.

- `max_iterations` was 100 across all presets, truncating EM before convergence. Observed
  requirements: 232 iterations (simulated 21-species community), 175 (ATCC ONT mock), 376
  (HiFi at temperature 0.5). Raised to 500 on all presets so that `convergence_threshold`
  rather than the iteration ceiling determines termination. Estimates from 0.1.0 were
  conservative rather than wrong: on the ATCC mock, converged values differ by at most
  0.007 percentage points; on the simulated community, aggregate L1 improves from 7.48%
  to 7.30%.

- `check-db` could loop indefinitely if stdin entered a persistent error state
  (`filter_map` over `io::Lines` skips errors and continues; now `map_while`).

- Corrected stale `emu-its` binary name in `scripts/run_real_data_validation.sh`.

### Added

- EUKARYOME reference database support via `--db-format {unite,eukaryome,auto}`. Defaults
  to `auto`, which infers the format from alignment target names, so existing UNITE
  workflows are unchanged. Taxonomy is anchored on the trailing seven colon-delimited
  fields rather than fixed offsets; species epithets are combined with the genus to give
  keys comparable to UNITE.
- `emits check-db` subcommand: reads FASTA headers on stdin and reports how they parse,
  exiting non-zero if a database release cannot be interpreted.
- Regression tests asserting that no preset filters on MAPQ, and that every iteration
  ceiling exceeds the highest observed requirement.

### Changed

- Figure-generation scripts now produce three-way EMITS/naive/EMU comparisons and include
  peer-review revisions to layout.
- `scripts/simulate_its_community.py` gains `--platform
  {ont-r10,pacbio-hifi,pacbio-hifi-random}` for simulating platform-specific read accuracy.
- Added `scripts/sweep_hifi_preset.py` for preset parameter sweeps.

## [0.1.0] - 2026-03-18

### Added
- Initial release of EMITS
- EM algorithm for probabilistic abundance estimation from minimap2 PAF output
- Platform presets: `ont-r10`, `ont-r9`, `pacbio-hifi`, `ont-duplex`
- UNITE header parsing and species-level taxonomic aggregation
- Naive best-hit counting for comparison (`--compare` flag)
- Built-in simulation framework (`emits simulate`)
- Species-level and raw per-accession output formats
- Validation against ONT ATCC mock community (10 species)
- Validation against synthetic UNITE community (21 species)
