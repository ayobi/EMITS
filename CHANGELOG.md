# Changelog

All notable changes to EMITS will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
