-- list_types_for_taxid.sql
-- -------------------------------------------------------
-- List all RNA types (and counts) present in RNAcentral
-- for a given NCBI taxonomy ID.
--
-- Usage: replace :TAXID with your NCBI taxid (e.g. 9913
--        for Bos taurus) before running, OR use the
--        psql \set command / the wrapper scripts which
--        substitute TAXID automatically.
--
-- Example result columns:
--   rna_type | n
-- -------------------------------------------------------

SELECT
    pre.rna_type,
    COUNT(*) AS n
FROM rnc_rna_precomputed pre
WHERE pre.taxid = :TAXID
  AND pre.is_active = true
GROUP BY pre.rna_type
ORDER BY n DESC;
