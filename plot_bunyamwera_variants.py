#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Variant distribution by segment from VCFs (AF from DP4), with:
 - Main plot (no legend)
 - Legend-only image

Requirements: matplotlib, pandas (standard Python + pip install if needed)
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# 1) Inputs
# -----------------------------
# Update these paths if your files live elsewhere.
vcf_files = {
    "SRR34843736": "SRR34843736.bunyamwera.vcf",
    "SRR34843737": "SRR34843737.bunyamwera.vcf",
    "SRR34843738": "SRR34843738.bunyamwera.vcf",
    "SRR34843739": "SRR34843739.bunyamwera.vcf",
}

# Segment order for concatenation (must match contig names present in the VCFs)
chrom_order = ["MF926354.1", "KP795084.1", "KP795093.1"]

# Deep color mapping (you can change these if desired)
colors = {
    "SRR34843736": "crimson",      # deep red
    "SRR34843737": "navy",         # deep blue
    "SRR34843738": "indigo",       # deep purple
    "SRR34843739": "darkgoldenrod" # deep yellow
}

# -----------------------------
# 2) Helpers
# -----------------------------
def parse_vcf_dp4(file_path: str, sample_id: str):
    """
    Parse a VCF and compute AF from DP4 (ref_f, ref_r, alt_f, alt_r).
    Returns a list of dicts: {Sample, Chrom, Pos, AF}
    """
    out = []
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"VCF not found: {file_path}")

    with open(file_path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\t\r\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, flt, info = parts[:8]
            try:
                pos = int(pos)
            except ValueError:
                continue

            # Compute AF from DP4 if present
            af = None
            if "DP4=" in info:
                for kv in info.split(";"):
                    if kv.startswith("DP4="):
                        try:
                            ref_f, ref_r, alt_f, alt_r = map(int, kv.split("=")[1].split(","))
                            total = ref_f + ref_r + alt_f + alt_r
                            if total > 0:
                                af = (alt_f + alt_r) / total
                        except Exception:
                            pass

            if af is not None:
                out.append({"Sample": sample_id, "Chrom": chrom, "Pos": pos, "AF": af})
    return out

# -----------------------------
# 3) Load variants from all VCFs
# -----------------------------
records = []
for sid, path in vcf_files.items():
    records.extend(parse_vcf_dp4(path, sid))

df = pd.DataFrame(records)
if df.empty:
    raise RuntimeError("No variants parsed. Check that DP4 is present in INFO fields of the VCFs.")

# -----------------------------
# 4) Build cumulative coordinates by segment
# -----------------------------
offsets = {}
lengths = {}
cum = 0
for chrom in chrom_order:
    seg = df[df["Chrom"] == chrom]
    seg_len = int(seg["Pos"].max()) if not seg.empty else 0
    offsets[chrom] = cum
    lengths[chrom] = seg_len
    cum += seg_len

# Map each variant to cumulative x coordinate
df["Coord"] = df.apply(lambda r: r["Pos"] + offsets.get(r["Chrom"], 0), axis=1)

# -----------------------------
# 5) Plot: main chart WITHOUT legend
# -----------------------------
fig, ax = plt.subplots(figsize=(15, 6))

# Scatter per sample; marker size proportional to AF
for sid in df["Sample"].unique():
    sub = df[df["Sample"] == sid]
    ax.scatter(
        sub["Coord"],
        [sid] * len(sub),
        s=sub["AF"] * 220,
        alpha=0.65,
        color=colors.get(sid, None),
        marker="x",
        linewidths=1.2,
    )

# Segment bands, boundary lines, and labels with exact spans
ranges = []
for i, chrom in enumerate(chrom_order):
    start = offsets[chrom]
    end = start + lengths[chrom]
    ranges.append((chrom, start, end))

    # light band behind points
    ax.axvspan(start, end, alpha=0.04, color="gray")

    # central label below axis: segment name + span
    mid = (start + end) / 2 if end > start else start
    ax.annotate(
        f"{chrom}\n{start:,}–{end:,}",
        xy=(mid, -0.16),
        xycoords=("data", "axes fraction"),
        ha="center",
        va="top",
        fontsize=10,
    )

    # vertical dashed boundary + numeric marker at boundary
    if i > 0:
        ax.axvline(start, color="goldenrod", linestyle="--", alpha=0.9)
        ax.annotate(
            f"{start:,}",
            xy=(start, -0.02),
            xycoords=("data", "axes fraction"),
            ha="center",
            va="top",
            fontsize=9,
        )

# rightmost end marker
total_end = ranges[-1][2]
ax.annotate(
    f"{total_end:,}",
    xy=(total_end, -0.02),
    xycoords=("data", "axes fraction"),
    ha="center",
    va="top",
    fontsize=9,
)

# axis formatting
ax.set_xlabel("Cumulative coordinate (concatenated segments)")
ax.set_ylabel("")    # y-ticks already show SRR IDs
ax.grid(axis="x", linestyle=":", alpha=0.3)

# Leave legend OFF on the main figure
plt.subplots_adjust(bottom=0.28, left=0.13, right=0.97, top=0.92)

out_main = "variant_distribution_bunyamwera_segments_nolegend.jpg"
plt.savefig(out_main, dpi=300)
plt.close(fig)

print(f"Saved main plot: {out_main}")

# -----------------------------
# 6) Legend-only image
# -----------------------------
fig_leg = plt.figure(figsize=(4, 2))
ax_leg = fig_leg.add_subplot(111)

handles = []
labels = []
for sid in df["Sample"].unique():
    h = ax_leg.scatter([], [], color=colors.get(sid, None), marker="x", label=sid)
    handles.append(h)
    labels.append(sid)

leg = ax_leg.legend(handles=handles, labels=labels, title="Sample", loc="center")
ax_leg.axis("off")

out_leg = "variant_distribution_legend.jpg"
plt.savefig(out_leg, dpi=300, bbox_inches="tight")
plt.close(fig_leg)

print(f"Saved legend image: {out_leg}")
