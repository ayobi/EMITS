//! PAF (Pairwise mApping Format) parser for minimap2 output.
//!
//! Extracts read-to-reference alignments with mapping quality and
//! alignment score information needed for the EM algorithm.

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// A single alignment record from a PAF file.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PafRecord {
    /// Read/query name
    pub query_name: String,
    /// Query sequence length
    pub query_len: u64,
    /// Reference/target name (e.g., UNITE accession)
    pub target_name: String,
    /// Target sequence length
    pub target_len: u64,
    /// Number of matching bases
    pub matches: u64,
    /// Alignment block length
    pub block_len: u64,
    /// Mapping quality (0-255)
    pub mapq: u8,
    /// Alignment score from minimap2 (AS:i: tag), if present
    pub alignment_score: Option<i64>,
}

impl PafRecord {
    /// Parse a single PAF line into a PafRecord.
    ///
    /// PAF format (tab-separated):
    /// 0: query name, 1: query len, 2: query start, 3: query end,
    /// 4: strand (+/-), 5: target name, 6: target len, 7: target start,
    /// 8: target end, 9: num matches, 10: block len, 11: mapq,
    /// 12+: optional tags (SAM-like)
    pub fn from_line(line: &str) -> Result<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            anyhow::bail!("PAF line has fewer than 12 fields: {}", line);
        }

        // Parse the alignment score from optional tags (AS:i:VALUE)
        let alignment_score = fields[12..]
            .iter()
            .find(|f| f.starts_with("AS:i:"))
            .and_then(|f| f[5..].parse::<i64>().ok());

        Ok(PafRecord {
            query_name: fields[0].to_string(),
            query_len: fields[1].parse().context("invalid query length")?,
            target_name: fields[5].to_string(),
            target_len: fields[6].parse().context("invalid target length")?,
            matches: fields[9].parse().context("invalid match count")?,
            block_len: fields[10].parse().context("invalid block length")?,
            mapq: fields[11].parse().context("invalid mapq")?,
            alignment_score,
        })
    }

    /// Compute alignment identity as matches / block_len.
    pub fn identity(&self) -> f64 {
        if self.block_len == 0 {
            return 0.0;
        }
        self.matches as f64 / self.block_len as f64
    }

    /// Compute a normalized alignment score.
    /// Uses alignment score / query length if AS tag is present,
    /// otherwise falls back to identity * (block_len / query_len).
    pub fn normalized_score(&self) -> f64 {
        if let Some(as_score) = self.alignment_score {
            if self.query_len == 0 {
                return 0.0;
            }
            // Normalize by query length so shorter/longer ITS regions
            // don't systematically bias the score
            as_score as f64 / self.query_len as f64
        } else {
            // Fallback: identity weighted by alignment coverage
            if self.query_len == 0 {
                return 0.0;
            }
            self.identity() * (self.block_len as f64 / self.query_len as f64)
        }
    }
}

/// Grouped alignments: for each read, all its candidate reference hits.
pub type ReadAlignments = HashMap<String, Vec<PafRecord>>;

/// Parse a PAF file and group alignments by read name.
///
/// # Arguments
/// * `reader` - Any buffered reader (file, stdin, etc.)
/// * `min_identity` - Minimum alignment identity threshold (0.0-1.0)
/// * `min_mapq` - Minimum mapping quality threshold
pub fn parse_paf<R: BufRead>(reader: R, min_identity: f64, min_mapq: u8) -> Result<ReadAlignments> {
    let mut alignments: ReadAlignments = HashMap::new();
    let mut total = 0u64;
    let mut passed = 0u64;

    for line in reader.lines() {
        let line = line.context("failed to read PAF line")?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        total += 1;
        let record = PafRecord::from_line(&line)?;

        // Apply quality filters
        if record.identity() < min_identity {
            continue;
        }
        if record.mapq < min_mapq {
            continue;
        }

        passed += 1;
        alignments
            .entry(record.query_name.clone())
            .or_default()
            .push(record);
    }

    log::info!(
        "Parsed {} alignments, {} passed filters ({:.1}%), {} unique reads",
        total,
        passed,
        if total > 0 {
            passed as f64 / total as f64 * 100.0
        } else {
            0.0
        },
        alignments.len()
    );

    Ok(alignments)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_paf_line(query: &str, target: &str, matches: u64, block_len: u64, mapq: u8) -> String {
        // Minimal PAF line with required 12 fields + AS tag
        format!(
            "{}\t1000\t0\t{}\t+\t{}\t800\t0\t{}\t{}\t{}\t{}\tAS:i:{}",
            query,
            block_len,
            target,
            block_len,
            matches,
            block_len,
            mapq,
            matches * 2 // simple score
        )
    }

    #[test]
    fn test_parse_paf_record() {
        let line = make_paf_line("read1", "species_A", 900, 1000, 60);
        let record = PafRecord::from_line(&line).unwrap();
        assert_eq!(record.query_name, "read1");
        assert_eq!(record.target_name, "species_A");
        assert!((record.identity() - 0.9).abs() < 1e-6);
    }

    #[test]
    fn test_identity_filter() {
        let lines = vec![
            make_paf_line("read1", "sp_A", 950, 1000, 60), // 0.95 identity
            make_paf_line("read1", "sp_B", 700, 1000, 60), // 0.70 identity
        ];
        let input = lines.join("\n");
        let reader = std::io::BufReader::new(input.as_bytes());
        let alignments = parse_paf(reader, 0.8, 0).unwrap();

        let read1 = &alignments["read1"];
        assert_eq!(read1.len(), 1);
        assert_eq!(read1[0].target_name, "sp_A");
    }

    #[test]
    fn test_multi_mapping() {
        let lines = vec![
            make_paf_line("read1", "sp_A", 950, 1000, 60),
            make_paf_line("read1", "sp_B", 940, 1000, 60),
            make_paf_line("read2", "sp_A", 960, 1000, 60),
        ];
        let input = lines.join("\n");
        let reader = std::io::BufReader::new(input.as_bytes());
        let alignments = parse_paf(reader, 0.8, 0).unwrap();

        assert_eq!(alignments.len(), 2); // 2 reads
        assert_eq!(alignments["read1"].len(), 2); // read1 maps to 2 species
        assert_eq!(alignments["read2"].len(), 1); // read2 maps to 1 species
    }
}
