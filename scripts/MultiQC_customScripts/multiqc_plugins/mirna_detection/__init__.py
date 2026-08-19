"""
MultiQC custom module: mirna_detection
=======================================
Reads mirna_counts/*_mirna_stats.txt produced by 5.2_mirna_quantification.sh
and renders:
  1. Stacked bar chart – mapped vs unmapped reads per sample
  2. Bar chart (toggling between ≥1 / ≥10 / ≥100 read views) – distinct
     miRNAs detected per sample at each read-count threshold

Search-pattern key (defined in multiqc_config.yaml):
    mirna_detection/stats  →  *_mirna_stats.txt
"""

import logging
import re
from collections import OrderedDict

log = logging.getLogger(__name__)

try:
    from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
except ImportError:
    from multiqc.modules.base_module import BaseMultiqcModule

    class ModuleNoSamplesFound(Exception):
        pass

from multiqc.plots import bargraph

# Add the multiqc_plugins/ directory to sys.path so that _utils can be
# imported as a top-level module when loaded via MultiQC's --module-dir flag
# (in that context modules are not part of the multiqc_plugins package and
# relative imports would fail).
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from _utils import parse_rna_stats as _parse_stats_util

# ── Section metadata ─────────────────────────────────────────────────────────
MODULE_NAME     = "miRNA Detection"
MODULE_ANCHOR   = "mirna_detection"
MODULE_INFO     = (
    "Per-sample miRNA mapping statistics from Bowtie1 alignment to the unified "
    "mature miRNA reference (Phase V.2). Distinct miRNA counts are shown at "
    "three read-count thresholds: ≥1, ≥10, and ≥100 reads."
)
SP_KEY          = "mirna_detection/stats"
FILENAME_SUFFIX = "_mirna_stats.txt"
RNA_LABEL       = "miRNAs"

# ── Plot colours ─────────────────────────────────────────────────────────────
MAPPED_COLOR   = "#4e79a7"
UNMAPPED_COLOR = "#cccccc"
TH1_COLOR      = "#4e79a7"
TH10_COLOR     = "#f28e2b"
TH100_COLOR    = "#e15759"

# Standalone PNG produced by plot_mirna_detection.py
STANDALONE_PNG = "mirna_distinct_detected.png"


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name=MODULE_NAME,
            anchor=MODULE_ANCHOR,
            info=MODULE_INFO,
        )

        self.rna_data = {}

        for f in self.find_log_files(SP_KEY):
            fn = f["fn"]
            m = re.match(r"^(.+)" + re.escape(FILENAME_SUFFIX) + r"$", fn)
            s_name = self.clean_s_name(m.group(1) if m else f["s_name"], f)
            parsed = _parse_stats(f["f"], rna_label=RNA_LABEL)
            if parsed:
                self.rna_data[s_name] = parsed
                self.add_data_source(f, s_name)

        if not self.rna_data:
            raise ModuleNoSamplesFound

        log.info(f"Found {MODULE_NAME} stats for {len(self.rna_data)} samples")

        # ── Plot 1: mapped vs unmapped (stacked) ─────────────────────────────
        bar_data = {
            s: {"mapped": v["mapped"], "unmapped": v["unmapped"]}
            for s, v in self.rna_data.items()
        }
        map_cats = OrderedDict(
            [
                ("mapped",   {"name": "Mapped reads",   "color": MAPPED_COLOR}),
                ("unmapped", {"name": "Unmapped reads", "color": UNMAPPED_COLOR}),
            ]
        )
        map_pconfig = {
            "id": f"{MODULE_ANCHOR}_mapping_bargraph",
            "title": f"{MODULE_NAME}: Mapped vs Unmapped Reads per Sample",
            "ylab": "Read count",
            "stacking": "normal",
            "tt_percentage": True,
        }

        self.add_section(
            name="Alignment statistics",
            anchor=f"{MODULE_ANCHOR}_alignment",
            description=(
                f"Stacked bar chart of reads mapped to the {RNA_LABEL} reference "
                "versus reads that did not align. Input reads are the genome-mapped "
                "reads from Phase III."
            ),
            plot=bargraph.plot(bar_data, map_cats, map_pconfig),
        )

        # ── Plot 2: distinct detected per sample (toggle ≥1 / ≥10 / ≥100) ───
        d1_data   = {s: {"detected": v["d1"]}   for s, v in self.rna_data.items()}
        d10_data  = {s: {"detected": v["d10"]}  for s, v in self.rna_data.items()}
        d100_data = {s: {"detected": v["d100"]} for s, v in self.rna_data.items()}

        cats_d1   = OrderedDict([("detected", {"name": f"Distinct {RNA_LABEL} (\u22651 read)",   "color": TH1_COLOR})])
        cats_d10  = OrderedDict([("detected", {"name": f"Distinct {RNA_LABEL} (\u226510 reads)", "color": TH10_COLOR})])
        cats_d100 = OrderedDict([("detected", {"name": f"Distinct {RNA_LABEL} (\u2265100 reads)","color": TH100_COLOR})])

        distinct_pconfig = {
            "id": f"{MODULE_ANCHOR}_distinct_bargraph",
            "title": f"{MODULE_NAME}: Distinct {RNA_LABEL} Detected per Sample",
            "ylab": f"Distinct {RNA_LABEL} detected",
            "data_labels": [
                {"name": "\u22651 read",   "ylab": f"Distinct {RNA_LABEL}"},
                {"name": "\u226510 reads", "ylab": f"Distinct {RNA_LABEL}"},
                {"name": "\u2265100 reads","ylab": f"Distinct {RNA_LABEL}"},
            ],
        }

        self.add_section(
            name=f"Distinct {RNA_LABEL} detected",
            anchor=f"{MODULE_ANCHOR}_distinct",
            description=(
                f"Number of distinct {RNA_LABEL} with at least 1, 10, or 100 "
                "mapped reads per sample. Use the buttons above the plot to toggle "
                "between thresholds."
            ),
            plot=bargraph.plot(
                [d1_data, d10_data, d100_data],
                [cats_d1, cats_d10, cats_d100],
                distinct_pconfig,
            ),
        )


# ---------------------------------------------------------------------------
# Helpers – thin wrappers around shared utilities in multiqc_plugins/_utils.py
# ---------------------------------------------------------------------------
def _parse_stats(content, rna_label="miRNAs"):
    return _parse_stats_util(content)
