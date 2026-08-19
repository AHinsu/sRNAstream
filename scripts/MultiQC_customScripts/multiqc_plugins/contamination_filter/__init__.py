"""
MultiQC custom module: contamination_filter
===========================================
Reads *_contamination_stats.txt files produced by 1.5b_contamination_filtering.sh
and renders a per-sample stacked bar chart showing reads removed at each
sequential contamination filtering step (rRNA → tRNA → snRNA → snoRNA →
piRNA → Other RNA) and the final clean reads that proceed to genome alignment.

Search-pattern key (defined in multiqc_config.yaml):
    contamination_filter/stats  →  *_contamination_stats.txt
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
# imported as a top-level module regardless of how MultiQC loads this plugin.
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

# ---------------------------------------------------------------------------
# Category definitions (order = stacking order in bar, bottom to top)
# ---------------------------------------------------------------------------
CATEGORIES = OrderedDict(
    [
        ("rRNA",   {"name": "rRNA",       "color": "#e15759"}),
        ("tRNA",   {"name": "tRNA",       "color": "#f28e2b"}),
        ("snRNA",  {"name": "snRNA",      "color": "#76b7b2"}),
        ("snoRNA", {"name": "snoRNA",     "color": "#edc948"}),
        ("piRNA",  {"name": "piRNA",      "color": "#b07aa1"}),
        ("other",  {"name": "Other RNA",  "color": "#bab0ac"}),
        ("clean",  {"name": "Clean reads","color": "#59a14f"}),
    ]
)

RNA_KEYS = {"rRNA", "tRNA", "snRNA", "snoRNA", "piRNA", "other"}


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="Contamination Filtering",
            anchor="contamination_filter",
            info=(
                "Reads removed at each sequential contamination filtering step "
                "(Phase I.5b). Filtering order: rRNA → tRNA → snRNA → snoRNA → "
                "piRNA → Other RNA. Clean reads (green) proceed to genome alignment."
            ),
        )

        self.contam_data = {}

        for f in self.find_log_files("contamination_filter/stats"):
            fn = f["fn"]
            m = re.match(r"^(.+)_contamination_stats\.txt$", fn)
            s_name = self.clean_s_name(m.group(1) if m else f["s_name"], f)
            parsed = _parse_contamination_stats(f["f"])
            if parsed:
                self.contam_data[s_name] = parsed
                self.add_data_source(f, s_name)

        if not self.contam_data:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.contam_data)} samples")

        plot_data = {
            s: {cat: counts.get(cat, 0) for cat in CATEGORIES}
            for s, counts in self.contam_data.items()
        }

        pconfig = {
            "id": "contamination_filter_bargraph",
            "title": "Contamination Filtering: Read Fate per Sample",
            "ylab": "Read count",
            "stacking": "normal",
            "hide_zero_cats": True,
        }

        self.add_section(
            name="Read fate per sample",
            anchor="contamination_filter_reads",
            description=(
                "Stacked bar chart of read counts at each contamination filtering step. "
                "Reads are removed sequentially; only unmapped reads advance to the next step."
            ),
            plot=bargraph.plot(plot_data, CATEGORIES, pconfig),
        )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _parse_contamination_stats(content):
    """Parse file content of *_contamination_stats.txt → dict of counts."""
    counts = {}
    for line in content.splitlines():
        line = line.strip()
        if line.startswith("Output Reads (Clean):"):
            val = line.split(":", 1)[1].strip().replace(",", "")
            try:
                counts["clean"] = int(val)
            except ValueError:
                pass
        else:
            m = re.match(r"^(\w+):\s+([\d,]+)\s*\(", line)
            if m and m.group(1) in RNA_KEYS:
                try:
                    counts[m.group(1)] = int(m.group(2).replace(",", ""))
                except ValueError:
                    pass
    return counts or None
