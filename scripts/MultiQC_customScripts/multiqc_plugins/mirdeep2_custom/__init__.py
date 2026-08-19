"""
MultiQC custom module: mirdeep2_custom
======================================
Reads mirdeep2/parsed/mirdeep2_all.tsv produced by 4.1_mirdeep2_analysis.sh
and renders a stacked bar chart showing the number of known and novel miRNA
predictions surviving each score threshold ≥1 through ≥10.

miRDeep2 runs on combined reads from all samples pooled together, so there is
one result file rather than one file per sample.  The x-axis groups represent
score thresholds, not individual samples.

Search-pattern key (defined in multiqc_config.yaml):
    mirdeep2_custom/all  →  mirdeep2_all.tsv
"""

import logging
import os
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
# imported as a top-level module when loaded via MultiQC's --module-dir flag.
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

SCORE_THRESHOLDS = list(range(1, 11))  # ≥1 … ≥10

MD2_CATS = OrderedDict(
    [
        ("known", {"name": "Known miRNA", "color": "#4e79a7"}),
        ("novel", {"name": "Novel miRNA", "color": "#f28e2b"}),
    ]
)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="miRDeep2 Analysis",
            anchor="mirdeep2_custom",
            info=(
                "miRDeep2 predictions from combined reads of all samples (Phase IV). "
                "The chart shows how many known and novel miRNA predictions survive "
                "each successive score threshold from ≥1 to ≥10."
            ),
        )

        rows = []
        found_file = False

        for f in self.find_log_files("mirdeep2_custom/all"):
            found_file = True
            parsed = _parse_mirdeep2_tsv(f["f"])
            rows.extend(parsed)
            self.add_data_source(f)

        if not found_file or not rows:
            raise ModuleNoSamplesFound

        log.info(f"Parsed {len(rows)} miRDeep2 predictions")

        # Build bargraph data: key = threshold label, value = {known, novel}
        plot_data = OrderedDict()
        for thresh in SCORE_THRESHOLDS:
            label = f"\u2265{thresh}"
            above = [r for r in rows if r["score"] >= thresh]
            plot_data[label] = {
                "known": sum(1 for r in above if r["category"] == "known"),
                "novel": sum(1 for r in above if r["category"] == "novel"),
            }

        pconfig = {
            "id": "mirdeep2_score_bargraph",
            "title": "miRDeep2: Known and Novel miRNA Predictions by Score Threshold",
            "ylab": "Number of miRNA predictions",
            "xlab": "miRDeep2 score threshold",
            "stacking": "normal",
            "tt_percentage": False,
        }

        self.add_section(
            name="Known and novel miRNA predictions by score threshold",
            anchor="mirdeep2_score_section",
            description=(
                "Stacked bar chart of miRDeep2 predictions at each score threshold. "
                "Known miRNAs (blue) have a match in the reference miRNA database; "
                "novel miRNAs (orange) are de-novo predictions. "
                "Higher score thresholds yield more confident (but fewer) predictions."
            ),
            plot=bargraph.plot(plot_data, MD2_CATS, pconfig),
        )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _parse_mirdeep2_tsv(content):
    """Parse mirdeep2_all.tsv content, return list of {score, category} dicts."""
    rows = []
    lines = content.splitlines()
    if not lines:
        return rows

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
        try:
            score = float(row.get("score", 0))
        except (ValueError, TypeError):
            score = 0.0
        category = row.get("category", "").strip().lower()
        if category in ("known", "novel"):
            rows.append({"score": score, "category": category})
    return rows
