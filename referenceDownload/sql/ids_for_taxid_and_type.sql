-- ids_for_taxid_and_type.sql
-- -------------------------------------------------------
-- Return a plain list of RNAcentral species-specific IDs
-- (format: URS<accession>_<taxid>) for a given NCBI
-- taxonomy ID and one or more RNA types.
--
-- This output is used as the input file for:
--   seqkit grep -f <output> rnacentral_species_specific_ids.fasta.gz
--
-- Schema note: rnc_rna_precomputed.id is the primary key
-- and already contains the full urs_taxid (e.g.
-- URS0000000001_9913), so it is selected directly.
--
-- Placeholders (replaced by the wrapper scripts):
--   :TAXID     -- NCBI taxonomy ID (integer), e.g. 9913
--   :RNA_TYPES -- Comma-separated quoted type list substituted
--                 directly into the IN clause below, e.g.
--                 'miRNA' for a single type, or
--                 'miRNA','tRNA','rRNA' for multiple types.
--
-- Run with psql flags:  -t -A  (tuples-only, unaligned)
-- so output is one ID per line with no extra formatting.
-- -------------------------------------------------------

SELECT DISTINCT
    pre.id AS rnacentral_id
FROM rnc_rna_precomputed pre
WHERE pre.taxid = :TAXID
  AND pre.rna_type IN (:RNA_TYPES)
  AND pre.is_active = true
ORDER BY rnacentral_id;
