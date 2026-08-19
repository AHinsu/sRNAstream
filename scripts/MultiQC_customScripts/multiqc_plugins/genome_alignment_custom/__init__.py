"""
MultiQC custom module: genome_alignment_custom
==============================================
Reads alignment/all_samples_alignment_stats.tsv produced by
3.3_alignment_stats.sh and renders:
  1. A stacked bar chart (mapped + unmapped reads per sample)
  2. A summary table with mapping percentages

This module is named "Mapping to Genome" in the report to match the
pipeline phase labelling requested by the user.

Search-pattern key (defined in multiqc_config.yaml):
    genome_alignment_custom/stats  →  all_samples_alignment_stats.tsv
"""

import logging
from collections import OrderedDict

log = logging.getLogger(__name__)

try:
    from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
except ImportError:
    from multiqc.modules.base_module import BaseMultiqcModule

    class ModuleNoSamplesFound(Exception):
        pass

from multiqc.plots import bargraph, table

ALIGN_CATS = OrderedDict(
    [
        ("mapped_reads",   {"name": "Mapped reads",   "color": "#4e79a7"}),
        ("unmapped_reads", {"name": "Unmapped reads", "color": "#cccccc"}),
    ]
)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="Mapping to Genome",
            anchor="genome_alignment_custom",
            info=(
                "Bowtie1 alignment statistics for each sample mapped to the Ensembl "
                "reference genome (Phase III.2). Mapped reads are passed to miRDeep2 "
                "(Phase IV) and miRNA quantification (Phase V)."
            ),
        )

        self.align_data = {}

        for f in self.find_log_files("genome_alignment_custom/stats"):
            parsed = _parse_alignment_tsv(f["f"])
            if parsed:
                self.align_data.update(parsed)
                self.add_data_source(f)

        if not self.align_data:
            raise ModuleNoSamplesFound

        log.info(f"Found genome alignment stats for {len(self.align_data)} samples")

        # ── Plot 1: stacked bargraph (mapped / unmapped) ──────────────────────
        bar_data = {
            s: {
                "mapped_reads":   v["mapped_reads"],
                "unmapped_reads": v["unmapped_reads"],
            }
            for s, v in self.align_data.items()
        }

        bar_pconfig = {
            "id": "genome_alignment_bargraph",
            "title": "Genome Alignment: Mapped vs Unmapped Reads per Sample",
            "ylab": "Read count",
            "stacking": "normal",
            "hide_zero_cats": False,
        }

        # ── Plot 2: summary table ─────────────────────────────────────────────
        tbl_data = {
            s: {
                "total_reads":    v["total_reads"],
                "mapped_reads":   v["mapped_reads"],
                "unmapped_reads": v["unmapped_reads"],
                "pct_mapped":     v["pct_mapped"],
            }
            for s, v in self.align_data.items()
        }

        tbl_headers = OrderedDict(
            [
                ("total_reads",    {"title": "Total reads",     "format": "{:,.0f}", "scale": "Blues"}),
                ("mapped_reads",   {"title": "Mapped reads",    "format": "{:,.0f}", "scale": "Greens"}),
                ("unmapped_reads", {"title": "Unmapped reads",  "format": "{:,.0f}", "scale": "Greys"}),
                ("pct_mapped",     {"title": "% Mapped",        "format": "{:.2f}%", "scale": "RdYlGn", "min": 0, "max": 100}),
            ]
        )

        tbl_pconfig = {
            "id": "genome_alignment_table",
            "title": "Genome Alignment Summary",
        }

        self.add_section(
            name="Alignment statistics",
            anchor="genome_alignment_stats",
            description=(
                "Per-sample Bowtie1 alignment rates to the Ensembl reference genome. "
                "Mapped reads (blue) are used in all downstream phases."
            ),
            plot=bargraph.plot(bar_data, ALIGN_CATS, bar_pconfig),
        )

        self.add_section(
            name="Alignment summary table",
            anchor="genome_alignment_table",
            description="Per-sample read counts and mapping percentages.",
            plot=table.plot(tbl_data, tbl_headers, tbl_pconfig),
        )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _parse_alignment_tsv(content):
    """
    Parse all_samples_alignment_stats.tsv.
    Expected columns: sample  total_reads  mapped_reads  unmapped_reads  pct_mapped  reported_alignments
    Returns dict of {sample: {...}}.
    """
    data = {}
    lines = content.splitlines()
    if not lines:
        return data

    header = None
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if header is None:
            header = line.split("\t")
            continue
        parts = line.split("\t")
        row = dict(zip(header, parts))
        sample = row.get("sample", "").strip()
        if not sample:
            continue
        try:
            data[sample] = {
                "total_reads":    int(row.get("total_reads", 0)),
                "mapped_reads":   int(row.get("mapped_reads", 0)),
                "unmapped_reads": int(row.get("unmapped_reads", 0)),
                "pct_mapped":     float(row.get("pct_mapped", 0.0)),
            }
        except (ValueError, TypeError):
            pass
    return data
