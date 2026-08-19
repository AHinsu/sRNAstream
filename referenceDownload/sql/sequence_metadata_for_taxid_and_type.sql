-- sequence_metadata_for_taxid_and_type.sql
-- -------------------------------------------------------
-- Fetch core metadata for active RNA sequences in
-- RNAcentral for a given NCBI taxonomy ID and one or more
-- RNA types.  Output is a tab-separated table (TSV) with
-- a header row, suitable for inspection in a spreadsheet
-- or downstream analysis scripts.
--
-- This is the fallback query used when the WITH-features
-- query (sequence_metadata_for_taxid_and_type_with_features.sql)
-- returns an empty result (e.g. due to a Cartesian product
-- explosion for RNA types like tRNA that have many feature
-- records per sequence).
--
-- Performance: a CTE (filtered_rna) pre-filters
-- rnc_rna_precomputed to only the requested TAXID and
-- rna_type rows BEFORE the joins, so the expensive join
-- operations work on a small dataset rather than the full
-- rnc_rna_precomputed table.
--
-- Columns returned (17 total):
--   rnacentral_id      -- species-specific URS, e.g. URS0000000001_9913
--   rna_type           -- RNAcentral RNA type (miRNA, tRNA, rRNA, …)
--   description        -- human-readable sequence description
--   short_description  -- abbreviated sequence description
--   databases          -- contributing databases (pipe-separated)
--   accession_id       -- accession from rnc_accession_active
--   assembly_id        -- genome assembly identifier
--   region_id          -- genomic region identifier
--   gene_status        -- gene status from rnc_gene_status
--   region_name        -- name of the genomic region
--   chromosome         -- chromosome / scaffold name
--   strand             -- genomic strand (+ / -)
--   region_start       -- start coordinate of the region
--   region_stop        -- end coordinate of the region
--   was_mapped         -- whether the region was computationally mapped
--   providing_databases -- databases that provided the region annotation
--   exon_count         -- number of exons in the region
--
-- Schema notes:
--   * rnc_rna_precomputed.id is the full urs_taxid
--     (e.g. URS0000000001_9913); all urs_taxid columns in
--     joined tables reference this value.
--
-- Placeholders (replaced by the wrapper scripts):
--   :TAXID     -- NCBI taxonomy ID (integer), e.g. 9913
--   :RNA_TYPES -- Comma-separated quoted type list, e.g.
--                 'miRNA' for a single type, or
--                 'miRNA','tRNA','rRNA' for multiple types.
--
-- Run with psql flags:  -A -F $'\t' -q
-- to get clean TSV output (tab-separated, with header,
-- no row-count footer).
-- -------------------------------------------------------

WITH filtered_rna AS (
    SELECT * 
    FROM rnc_rna_precomputed
    WHERE taxid = :TAXID
      AND rna_type IN (:RNA_TYPES)
      AND is_active = true
)
SELECT DISTINCT
    pre.id                                  AS rnacentral_id,
    pre.rna_type,
    pre.description,
    pre.short_description,
    pre.databases,
    acc.accession                           AS accession_id,
    gs.assembly_id,
    gs.region_id,
    gs.status                               AS gene_status,
    sr.region_name,
    sr.chromosome,
    sr.strand,
    sr.region_start,
    sr.region_stop,
    sr.was_mapped,
    sr.providing_databases,
    sr.exon_count
FROM filtered_rna pre
LEFT JOIN rnc_accession_active acc
    ON acc.urs_taxid = pre.id
LEFT JOIN rnc_gene_status gs
    ON gs.urs_taxid = pre.id
LEFT JOIN rnc_sequence_regions sr
    ON sr.urs_taxid = pre.id
ORDER BY rnacentral_id, accession_id;
