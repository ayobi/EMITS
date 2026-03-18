//! UNITE taxonomy parsing and abundance aggregation.
//!
//! Parses UNITE-formatted reference sequence headers to extract taxonomic
//! information at all ranks (kingdom through species), and provides
//! aggregation of per-accession abundances to any taxonomic level.

use std::collections::HashMap;

/// Parsed taxonomy from a UNITE sequence header.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct Taxonomy {
    /// Original full header string
    pub raw: String,
    /// UNITE accession (e.g., "JF910285")
    pub accession: String,
    /// Species Hypothesis ID (e.g., "SH1061784.10FU")
    pub sh_id: String,
    /// Taxonomic ranks
    pub kingdom: String,
    pub phylum: String,
    pub class: String,
    pub order: String,
    pub family: String,
    pub genus: String,
    pub species: String,
}

impl Taxonomy {
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
            kingdom,
            phylum,
            class,
            order,
            family,
            genus,
            species,
        }
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
/// Takes a map of UNITE header -> abundance and returns a map of
/// taxon_name_at_rank -> summed abundance, along with the lineage
/// for each aggregated taxon.
pub fn aggregate_abundances(abundances: &HashMap<String, f64>, rank: TaxRank) -> AggregatedResult {
    let mut agg: HashMap<String, f64> = HashMap::new();
    let mut lineages: HashMap<String, String> = HashMap::new();
    let mut accession_counts: HashMap<String, usize> = HashMap::new();

    for (header, abundance) in abundances {
        let tax = Taxonomy::from_unite_header(header);
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
    /// Number of UNITE accessions collapsed into each taxon
    pub accession_counts: HashMap<String, usize>,
    /// Rank at which aggregation was performed
    pub rank: TaxRank,
}

#[cfg(test)]
mod tests {
    use super::*;

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

        let result = aggregate_abundances(&abundances, TaxRank::Species);
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

        let result = aggregate_abundances(&abundances, TaxRank::Genus);
        assert!((result.abundances["Fusarium"] - 0.15).abs() < 0.001);
        assert_eq!(result.accession_counts["Fusarium"], 2);
    }
}
