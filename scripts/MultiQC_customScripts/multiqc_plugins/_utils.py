"""
_utils.py – Shared helpers for the bovine-microRNA MultiQC plugin modules.

Provides:
    parse_rna_stats(content)  – parse *_stats.txt files produced by phases V/VI
    embed_png(png_filename)   – load a PNG from the analysis dir and return base64 HTML
"""

import base64
import os
import re


def parse_rna_stats(content):
    """
    Parse a per-sample stats file produced by 5.2_mirna_quantification.sh or
    6.2_snrna_snorna_quantification.sh.

    Expected file sections (UTF-8):
        Mapped to <type> reference : <N>  (M.MM%)
        Unmapped                   : <N>
        ...
        ≥1 read                    : <N>
        ≥10 reads                  : <N>
        ≥100 reads                 : <N>

    Returns a dict with keys: mapped, unmapped, d1, d10, d100
    or None if no values were found.
    """
    result = {"mapped": 0, "unmapped": 0, "d1": 0, "d10": 0, "d100": 0}
    for line in content.splitlines():
        line_s = line.strip()
        # Mapped reads: works for both "Mapped to miRNA reference" and
        # "Mapped to snRNA reference" etc.
        m = re.search(r"Mapped to \S.*?reference\s*:\s*([\d,]+)", line_s)
        if m:
            try:
                result["mapped"] = int(m.group(1).replace(",", ""))
            except ValueError:
                pass
        if line_s.startswith("Unmapped"):
            m = re.search(r":\s*([\d,]+)", line_s)
            if m:
                try:
                    result["unmapped"] = int(m.group(1).replace(",", ""))
                except ValueError:
                    pass
        # Unicode \u2265 = ≥
        # Check longest prefix first (≥100 > ≥10 > ≥1) so that e.g. the
        # "≥100 reads" line does not also trigger the "≥1" or "≥10" branches.
        if "\u2265100" in line_s and "reads" in line_s:
            m = re.search(r":\s*([\d,]+)", line_s)
            if m:
                try:
                    result["d100"] = int(m.group(1).replace(",", ""))
                except ValueError:
                    pass
        elif "\u226510" in line_s and "reads" in line_s:
            m = re.search(r":\s*([\d,]+)", line_s)
            if m:
                try:
                    result["d10"] = int(m.group(1).replace(",", ""))
                except ValueError:
                    pass
        elif "\u22651" in line_s and "read" in line_s:
            m = re.search(r":\s*([\d,]+)", line_s)
            if m:
                try:
                    result["d1"] = int(m.group(1).replace(",", ""))
                except ValueError:
                    pass
    return result if any(result.values()) else None


def embed_png(png_filename):
    """
    Look for *png_filename* inside the MultiQC analysis directory
    (``<analysis_dir>/multiqc/<png_filename>``) and return an HTML ``<img>``
    element with the image base64-encoded inline.

    Returns an empty string if the file cannot be found or read.
    """
    try:
        from multiqc import config as mqc_config

        for adir in getattr(mqc_config, "analysis_dir", []):
            png_path = os.path.join(adir, "multiqc", png_filename)
            if os.path.exists(png_path):
                with open(png_path, "rb") as fh:
                    b64 = base64.b64encode(fh.read()).decode()
                return (
                    '<div style="margin-top:15px;">'
                    "<p><strong>Standalone export (matplotlib):</strong></p>"
                    f'<img src="data:image/png;base64,{b64}" '
                    'style="max-width:100%;border:1px solid #ddd;border-radius:4px;" />'
                    "</div>"
                )
    except Exception:
        pass
    return ""
