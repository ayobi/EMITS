//! Reference-database taxonomy parsing and abundance aggregation.
//!
//! Parses reference sequence headers to extract taxonomic information at all
//! ranks (kingdom through species), and provides aggregation of per-accession
//! abundances to any taxonomic level.
//!
//! Two header formats are supported:
//!
//! * [`DbFormat::Unite`] — UNITE, pipe-delimited with `k__`/`p__`/... prefixes.
//! * [`DbFormat::Eukaryome`] — EUKARYOME general FASTA, colon-delimited.
//!
//! The EUKARYOME general FASTA release is described in Tedersoo et al. (2024)
//! *Database* baae043: fields are separated by colons, the nomenclature and
//! sequence-origin fields are dropped relative to the spreadsheet release, and
//! taxonomy is carried at kingdom, phylum, class, order, family, genus and the
//! *species epithet only*. Ranks below the level of identification are written
//! as a single dot.
//!
//! Because the number of leading metadata fields is not guaranteed stable
//! across releases, the EUKARYOME parser anchors on the **trailing** seven
//! fields rather than on fixed offsets, and refuses to guess when a header does
//! not have the expected shape (see [`Taxonomy::try_from_eukaryome_header`]).

use std::collections::HashMap;

/// Reference database header format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DbFormat {
    /// UNITE: `Species_name|accession|SH_id|type|k__X;p__X;...;s__X`
    Unite,
    /// EUKARYOME general FASTA: colon-delimited, trailing 7 taxonomy fields.
    Eukaryome,
}

/// Read-coverage values used by EUKARYOME to mark which rRNA markers a
/// reference sequence spans.
const EUKARYOME_MARKERS: [&str; 6] = ["SSU", "ITS", "LSU", "SSU-ITS", "ITS-LSU", "all"];

/// Number of trailing colon-separated fields carrying taxonomy in EUKARYOME.
const EUKARYOME_TAX_FIELDS: usize = 7;

impl DbFormat {
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "unite" => Some(DbFormat::Unite),
            "eukaryome" | "euk" => Some(DbFormat::Eukaryome),
            _ => None,
        }
    }

    pub fn available() -> &'static str {
        "unite, eukaryome, auto"
    }

    /// Infer the header format from a single header string.
    ///
    /// UNITE headers are pipe-delimited and carry `k__`-style rank prefixes;
    /// either signal alone is taken as conclusive. Anything else with enough
    /// colon-separated fields to carry a full lineage is treated as EUKARYOME.
    /// Ambiguous input falls back to UNITE, preserving prior behaviour.
    pub fn detect(header: &str) -> Self {
        if header.contains("k__") || header.contains('|') {
            return DbFormat::Unite;
        }
        if header.matches(':').count() >= EUKARYOME_TAX_FIELDS {
            return DbFormat::Eukaryome;
        }
        DbFormat::Unite
    }

    /// Infer the format from a batch of headers by majority vote.
    ///
    /// More robust than [`DbFormat::detect`] on a single header, which may be
    /// unusually short or malformed.
    pub fn detect_many<'a, I: IntoIterator<Item = &'a str>>(headers: I) -> Self {
        let mut unite = 0usize;
        let mut euk = 0usize;
        for h in headers.into_iter().take(1000) {
            match DbFormat::detect(h) {
                DbFormat::Unite => unite += 1,
                DbFormat::Eukaryome => euk += 1,
            }
        }
        if euk > unite {
            DbFormat::Eukaryome
        } else {
            DbFormat::Unite
        }
    }
}

/// Parsed taxonomy from a reference sequence header.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct Taxonomy {
    /// Original full header string
    pub raw: String,
    /// Database accession (e.g., "JF910285" for UNITE, "EUK1210054" for EUKARYOME)
    pub accession: String,
    /// UNITE Species Hypothesis ID (e.g., "SH1061784.10FU"); empty for EUKARYOME
    pub sh_id: String,
    /// EUKARYOME read coverage (SSU, ITS, LSU, SSU-ITS, ITS-LSU, all); empty for UNITE
    pub marker: String,
    /// Taxonomic ranks
    pub kingdom: String,
    pub phylum: String,
    pub class: String,
    pub order: String,
    pub family: String,
    pub genus: String,
    pub species: String,
}

/// Normalise an EUKARYOME taxonomy field.
///
/// Ranks below the level of identification are written as a single dot; a run
/// of dots is also treated as absent. Underscores are converted to spaces so
/// that species keys match the UNITE code path.
fn euk_field(raw: &str) -> String {
    let t = raw.trim();
    if t.is_empty() || t.chars().all(|c| c == '.') {
        return String::new();
    }
    t.replace('_', " ")
}

impl Taxonomy {
    /// Parse a header using an explicit format.
    pub fn from_header(header: &str, format: DbFormat) -> Self {
        match format {
            DbFormat::Unite => Taxonomy::from_unite_header(header),
            DbFormat::Eukaryome => Taxonomy::from_eukaryome_header(header),
        }
    }

    /// Parse a UNITE-formatted header string.
    ///
    /// Expected format:
    /// `Species_name|accession|SH_id|type|k__X;p__X;c__X;o__X;f__X;g__X;s__X`
    pub fn from_unite_header(header: &str) -> Self {
        let parts: Vec<&str> = header.split('|').collect();

        let accession = parts.get(1).unwrap_or(&"").to_string();
        let sh_id = parts.get(2).unwrap_or(&"").to_string();

        // Parse taxonomy string (last pipe-separated field containing semicolons)
        let tax_str = parts.iter().rfind(|p| p.contains(';')).unwrap_or(&"");

        let mut kingdom = String::new();
        let mut phylum = String::new();
        let mut class = String::new();
        let mut order = String::new();
        let mut family = String::new();
        let mut genus = String::new();
        let mut species = String::new();

        for field in tax_str.split(';') {
            let field = field.trim();
            if let Some(val) = field.strip_prefix("k__") {
                kingdom = val.to_string();
            } else if let Some(val) = field.strip_prefix("p__") {
                phylum = val.to_string();
            } else if let Some(val) = field.strip_prefix("c__") {
                class = val.to_string();
            } else if let Some(val) = field.strip_prefix("o__") {
                order = val.to_string();
            } else if let Some(val) = field.strip_prefix("f__") {
                family = val.to_string();
            } else if let Some(val) = field.strip_prefix("g__") {
                genus = val.to_string();
            } else if let Some(val) = field.strip_prefix("s__") {
                species = val.replace('_', " ");
            }
        }

        // If species is empty, try to get it from the first field
        if species.is_empty() {
            if let Some(name) = parts.first() {
                species = name.replace('_', " ");
            }
        }

        Taxonomy {
            raw: header.to_string(),
            accession,
            sh_id,
            marker: String::new(),
            kingdom,
            phylum,
            class,
            order,
            family,
            genus,
            species,
        }
    }

    /// Parse an EUKARYOME general-FASTA header, returning `None` if the header
    /// does not have the expected shape.
    ///
    /// Callers that need a total function should use
    /// [`Taxonomy::from_eukaryome_header`], which degrades to an unclassified
    /// record instead of failing. Prefer this variant when validating a
    /// database up front, so that a format change surfaces as an error rather
    /// than as silently unclassified abundance.
    ///
    /// # Format
    ///
    /// Fields are colon-delimited. Taxonomy occupies the final seven fields:
    /// kingdom, phylum, class, order, family, genus, species **epithet**. The
    /// epithet is combined with the genus to produce a binomial, since
    /// EUKARYOME stores `sapiens` where UNITE stores `Homo_sapiens`.
    pub fn try_from_eukaryome_header(header: &str) -> Option<Self> {
        let parts: Vec<&str> = header.split(':').collect();
        if parts.len() < EUKARYOME_TAX_FIELDS + 1 {
            return None;
        }

        let split_at = parts.len() - EUKARYOME_TAX_FIELDS;
        let (meta, tax) = parts.split_at(split_at);

        let accession = meta.first().unwrap_or(&"").trim().to_string();
        if accession.is_empty() {
            return None;
        }

        // Read coverage is not at a guaranteed offset; locate it by vocabulary.
        let marker = meta
            .iter()
            .find(|f| {
                EUKARYOME_MARKERS
                    .iter()
                    .any(|m| m.eq_ignore_ascii_case(f.trim()))
            })
            .map(|f| f.trim().to_string())
            .unwrap_or_default();

        let kingdom = euk_field(tax[0]);
        let phylum = euk_field(tax[1]);
        let class = euk_field(tax[2]);
        let order = euk_field(tax[3]);
        let family = euk_field(tax[4]);
        let genus = euk_field(tax[5]);
        let epithet = euk_field(tax[6]);

        // EUKARYOME carries the epithet alone; UNITE carries the binomial.
        // Reconstruct a binomial so species keys are comparable across formats.
        let species = if epithet.is_empty() {
            String::new()
        } else if genus.is_empty() || epithet.starts_with(&genus) {
            epithet
        } else {
            format!("{} {}", genus, epithet)
        };

        Some(Taxonomy {
            raw: header.to_string(),
            accession,
            sh_id: String::new(),
            marker,
            kingdom,
            phylum,
            class,
            order,
            family,
            genus,
            species,
        })
    }

    /// Parse an EUKARYOME header, degrading to an unclassified record when the
    /// header cannot be interpreted.
    pub fn from_eukaryome_header(header: &str) -> Self {
        Taxonomy::try_from_eukaryome_header(header).unwrap_or_else(|| Taxonomy {
            raw: header.to_string(),
            accession: header.split(':').next().unwrap_or("").trim().to_string(),
            ..Default::default()
        })
    }

    /// Get the taxon name at a given rank level.
    pub fn at_rank(&self, rank: TaxRank) -> &str {
        match rank {
            TaxRank::Kingdom => &self.kingdom,
            TaxRank::Phylum => &self.phylum,
            TaxRank::Class => &self.class,
            TaxRank::Order => &self.order,
            TaxRank::Family => &self.family,
            TaxRank::Genus => &self.genus,
            TaxRank::Species => &self.species,
        }
    }

    /// Format the full lineage as a semicolon-separated string.
    pub fn lineage(&self) -> String {
        format!(
            "k__{};p__{};c__{};o__{};f__{};g__{};s__{}",
            self.kingdom,
            self.phylum,
            self.class,
            self.order,
            self.family,
            self.genus,
            self.species.replace(' ', "_")
        )
    }
}

/// Taxonomic rank for aggregation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TaxRank {
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
}

impl TaxRank {
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "kingdom" | "k" => Some(TaxRank::Kingdom),
            "phylum" | "p" => Some(TaxRank::Phylum),
            "class" | "c" => Some(TaxRank::Class),
            "order" | "o" => Some(TaxRank::Order),
            "family" | "f" => Some(TaxRank::Family),
            "genus" | "g" => Some(TaxRank::Genus),
            "species" | "s" => Some(TaxRank::Species),
            _ => None,
        }
    }

    pub fn available() -> &'static str {
        "kingdom, phylum, class, order, family, genus, species"
    }
}

/// Aggregate per-accession abundances to a given taxonomic rank.
///
/// Takes a map of reference header -> abundance and returns a map of
/// taxon_name_at_rank -> summed abundance, along with the lineage
/// for each aggregated taxon.
pub fn aggregate_abundances(
    abundances: &HashMap<String, f64>,
    rank: TaxRank,
    format: DbFormat,
) -> AggregatedResult {
    let mut agg: HashMap<String, f64> = HashMap::new();
    let mut lineages: HashMap<String, String> = HashMap::new();
    let mut accession_counts: HashMap<String, usize> = HashMap::new();

    for (header, abundance) in abundances {
        let tax = Taxonomy::from_header(header, format);
        let key = tax.at_rank(rank).to_string();

        if key.is_empty() {
            // Unclassified at this rank — use "Unclassified"
            *agg.entry("Unclassified".to_string()).or_default() += abundance;
            *accession_counts
                .entry("Unclassified".to_string())
                .or_default() += 1;
        } else {
            *agg.entry(key.clone()).or_default() += abundance;
            *accession_counts.entry(key.clone()).or_default() += 1;
            // Keep the lineage from the highest-abundance accession
            lineages.entry(key).or_insert_with(|| tax.lineage());
        }
    }

    AggregatedResult {
        abundances: agg,
        lineages,
        accession_counts,
        rank,
    }
}

/// Result of taxonomic aggregation.
#[derive(Debug)]
pub struct AggregatedResult {
    /// Aggregated abundances (taxon_name -> summed abundance)
    pub abundances: HashMap<String, f64>,
    /// Representative lineage for each taxon
    pub lineages: HashMap<String, String>,
    /// Number of database accessions collapsed into each taxon
    pub accession_counts: HashMap<String, usize>,
    /// Rank at which aggregation was performed
    pub rank: TaxRank,
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── UNITE ────────────────────────────────────────────────────────────

    #[test]
    fn test_parse_unite_header() {
        let header = "Fusarium_oxysporum|JF910285|SH1061784.10FU|refs|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusarium;s__Fusarium_oxysporum";
        let tax = Taxonomy::from_unite_header(header);

        assert_eq!(tax.accession, "JF910285");
        assert_eq!(tax.sh_id, "SH1061784.10FU");
        assert_eq!(tax.kingdom, "Fungi");
        assert_eq!(tax.phylum, "Ascomycota");
        assert_eq!(tax.class, "Sordariomycetes");
        assert_eq!(tax.order, "Hypocreales");
        assert_eq!(tax.family, "Nectriaceae");
        assert_eq!(tax.genus, "Fusarium");
        assert_eq!(tax.species, "Fusarium oxysporum");
    }

    #[test]
    fn test_parse_incertae_sedis() {
        let header = "Agaricomycetes_sp|KM065545|SH1187164.10FU|reps|k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricomycetes_ord_Incertae_sedis;f__Agaricomycetes_fam_Incertae_sedis;g__Agaricomycetes_gen_Incertae_sedis;s__Agaricomycetes_sp";
        let tax = Taxonomy::from_unite_header(header);

        assert_eq!(tax.genus, "Agaricomycetes_gen_Incertae_sedis");
        assert_eq!(tax.species, "Agaricomycetes sp");
    }

    // ── EUKARYOME ────────────────────────────────────────────────────────

    #[test]
    fn test_parse_eukaryome_header() {
        let header = "EUK1210054:ITS:Fungi:Ascomycota:Sordariomycetes:Hypocreales:Nectriaceae:Fusarium:oxysporum";
        let tax = Taxonomy::try_from_eukaryome_header(header).expect("should parse");

        assert_eq!(tax.accession, "EUK1210054");
        assert_eq!(tax.marker, "ITS");
        assert_eq!(tax.kingdom, "Fungi");
        assert_eq!(tax.phylum, "Ascomycota");
        assert_eq!(tax.class, "Sordariomycetes");
        assert_eq!(tax.order, "Hypocreales");
        assert_eq!(tax.family, "Nectriaceae");
        assert_eq!(tax.genus, "Fusarium");
        // Epithet is promoted to a binomial to match the UNITE code path.
        assert_eq!(tax.species, "Fusarium oxysporum");
    }

    #[test]
    fn test_eukaryome_dots_are_unclassified() {
        // Ranks below the level of identification are written as a single dot.
        let header = "EUK0000001:ITS:Fungi:Ascomycota:Sordariomycetes:.:.:.:.";
        let tax = Taxonomy::try_from_eukaryome_header(header).expect("should parse");

        assert_eq!(tax.class, "Sordariomycetes");
        assert_eq!(tax.order, "");
        assert_eq!(tax.family, "");
        assert_eq!(tax.genus, "");
        assert_eq!(tax.species, "");
    }

    #[test]
    fn test_eukaryome_extra_leading_fields() {
        // Taxonomy is anchored on the trailing seven fields, so additional
        // leading metadata must not shift the lineage.
        let header = "EUK1210054:INSDc:curated_2024:ITS:Fungi:Ascomycota:Eurotiomycetes:Eurotiales:Aspergillaceae:Aspergillus:fumigatus";
        let tax = Taxonomy::try_from_eukaryome_header(header).expect("should parse");

        assert_eq!(tax.accession, "EUK1210054");
        assert_eq!(tax.marker, "ITS");
        assert_eq!(tax.genus, "Aspergillus");
        assert_eq!(tax.species, "Aspergillus fumigatus");
    }

    #[test]
    fn test_eukaryome_epithet_already_binomial() {
        // Defensive: if a release ever carries the full binomial, don't
        // produce "Fusarium Fusarium oxysporum".
        let header = "EUK1:ITS:Fungi:Ascomycota:Sordariomycetes:Hypocreales:Nectriaceae:Fusarium:Fusarium_oxysporum";
        let tax = Taxonomy::try_from_eukaryome_header(header).expect("should parse");
        assert_eq!(tax.species, "Fusarium oxysporum");
    }

    #[test]
    fn test_eukaryome_rejects_malformed() {
        assert!(Taxonomy::try_from_eukaryome_header("EUK1:ITS:Fungi").is_none());
        assert!(Taxonomy::try_from_eukaryome_header("").is_none());
        // Degrading variant must not panic and must not invent taxonomy.
        let tax = Taxonomy::from_eukaryome_header("EUK1:ITS:Fungi");
        assert_eq!(tax.accession, "EUK1");
        assert_eq!(tax.species, "");
    }

    // ── Format detection ─────────────────────────────────────────────────

    #[test]
    fn test_detect_format() {
        let unite = "Fusarium_oxysporum|JF910285|SH1061784.10FU|refs|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusarium;s__Fusarium_oxysporum";
        let euk = "EUK1210054:ITS:Fungi:Ascomycota:Sordariomycetes:Hypocreales:Nectriaceae:Fusarium:oxysporum";

        assert_eq!(DbFormat::detect(unite), DbFormat::Unite);
        assert_eq!(DbFormat::detect(euk), DbFormat::Eukaryome);
        assert_eq!(
            DbFormat::detect_many(vec![euk, euk, unite]),
            DbFormat::Eukaryome
        );
        assert_eq!(
            DbFormat::detect_many(vec![unite, unite, euk]),
            DbFormat::Unite
        );
        // Ambiguous input falls back to UNITE (prior behaviour).
        assert_eq!(DbFormat::detect("something_odd"), DbFormat::Unite);
    }

    // ── Aggregation ──────────────────────────────────────────────────────

    #[test]
    fn test_aggregate_species() {
        let mut abundances = HashMap::new();
        // Two accessions of the same species
        abundances.insert(
            "Nakaseomyces_glabratus|MF767833|SH001|refs|k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Nakaseomyces;s__Nakaseomyces_glabratus".to_string(),
            0.119,
        );
        abundances.insert(
            "Nakaseomyces_glabratus|KP674599|SH002|refs|k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Nakaseomyces;s__Nakaseomyces_glabratus".to_string(),
            0.003,
        );
        // Different species
        abundances.insert(
            "Candida_albicans|AB001|SH003|refs|k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Debaryomycetaceae;g__Candida;s__Candida_albicans".to_string(),
            0.144,
        );

        let result = aggregate_abundances(&abundances, TaxRank::Species, DbFormat::Unite);
        assert!((result.abundances["Nakaseomyces glabratus"] - 0.122).abs() < 0.001);
        assert!((result.abundances["Candida albicans"] - 0.144).abs() < 0.001);
        assert_eq!(result.accession_counts["Nakaseomyces glabratus"], 2);
        assert_eq!(result.accession_counts["Candida albicans"], 1);
    }

    #[test]
    fn test_aggregate_genus() {
        let mut abundances = HashMap::new();
        abundances.insert(
            "Fusarium_oxysporum|A|SH1|refs|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusarium;s__Fusarium_oxysporum".to_string(),
            0.10,
        );
        abundances.insert(
            "Fusarium_solani|B|SH2|refs|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusarium;s__Fusarium_solani".to_string(),
            0.05,
        );

        let result = aggregate_abundances(&abundances, TaxRank::Genus, DbFormat::Unite);
        assert!((result.abundances["Fusarium"] - 0.15).abs() < 0.001);
        assert_eq!(result.accession_counts["Fusarium"], 2);
    }

    #[test]
    fn test_aggregate_eukaryome_consolidates_accessions() {
        // The accession-redundancy behaviour that motivates EMITS must hold
        // for EUKARYOME too: several accessions, one species key.
        let mut abundances = HashMap::new();
        abundances.insert(
            "EUK0000001:ITS:Fungi:Ascomycota:Saccharomycetes:Saccharomycetales:Saccharomycetaceae:Nakaseomyces:glabratus".to_string(),
            0.110,
        );
        abundances.insert(
            "EUK0000002:ITS:Fungi:Ascomycota:Saccharomycetes:Saccharomycetales:Saccharomycetaceae:Nakaseomyces:glabratus".to_string(),
            0.013,
        );
        abundances.insert(
            "EUK0000003:ITS:Fungi:Ascomycota:Saccharomycetes:Saccharomycetales:Debaryomycetaceae:Candida:albicans".to_string(),
            0.144,
        );

        let result = aggregate_abundances(&abundances, TaxRank::Species, DbFormat::Eukaryome);
        assert!((result.abundances["Nakaseomyces glabratus"] - 0.123).abs() < 0.001);
        assert_eq!(result.accession_counts["Nakaseomyces glabratus"], 2);
        assert!((result.abundances["Candida albicans"] - 0.144).abs() < 0.001);
    }

    #[test]
    fn test_aggregate_eukaryome_unclassified_not_silently_dropped() {
        let mut abundances = HashMap::new();
        abundances.insert(
            "EUK0000004:ITS:Fungi:Ascomycota:Sordariomycetes:.:.:.:.".to_string(),
            0.25,
        );
        let result = aggregate_abundances(&abundances, TaxRank::Species, DbFormat::Eukaryome);
        assert!((result.abundances["Unclassified"] - 0.25).abs() < 1e-9);
    }
}
