#!/usr/bin/env python3
"""
STEP 3: Dual-Pol Decomposition + Visualization + Cleanup
==========================================================

Runs 5 decomposition indices on each scene's C2 folder using polsartools,
then generates:
  - ONE master chart per scene with ALL maps (VV, VH, + 7 decomposition params)
  - All-scenes mega chart (rows=dates, cols=parameters)
  - Time series chart across all scenes
  - Cleans up intermediate SNAP files

Usage:
    conda activate intnotebook
    python step3_decompose.py
"""

import os
import sys
import glob
import shutil
import numpy as np
import warnings
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import rasterio

WORK_DIR = os.path.expanduser("~/Desktop/Muddasir")
PREPROCESSED_DIR = os.path.join(WORK_DIR, "outputs", "preprocessed")
FIGURES_DIR = os.path.join(WORK_DIR, "outputs", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)
AOI_SHP = os.path.join(WORK_DIR, "aoi", "Battambang.shp")
WINDOW_SIZE = 5


def find_c2_folders():
    folders = []
    if not os.path.isdir(PREPROCESSED_DIR): return folders
    for d in sorted(os.listdir(PREPROCESSED_DIR)):
        c2 = os.path.join(PREPROCESSED_DIR, d, "C2")
        req = ["C11.tif", "C12_real.tif", "C12_imag.tif", "C22.tif"]
        if all(os.path.isfile(os.path.join(c2, f)) for f in req):
            folders.append((d, c2))
    return folders


def get_date(name):
    for p in name.split("_"):
        if len(p) >= 8 and p[:8].isdigit():
            return f"{p[:4]}-{p[4:6]}-{p[6:8]}"
    return name[:20]


def cleanup(scene_dir):
    freed = 0
    for iw in ["IW1", "IW2", "IW3"]:
        for sub in [f"intermediate_{iw}.data", f"intermediate_{iw}.dim",
                    f"geocoded_{iw}", f"graph_{iw}.xml"]:
            p = os.path.join(scene_dir, sub)
            if os.path.isdir(p):
                for f in os.listdir(p):
                    fp = os.path.join(p, f)
                    if os.path.isfile(fp): freed += os.path.getsize(fp)
                shutil.rmtree(p, ignore_errors=True)
            elif os.path.isfile(p):
                freed += os.path.getsize(p); os.remove(p)
    # Also old single-IW intermediates
    for sub in ["intermediate_c2.data", "intermediate_c2.dim"]:
        p = os.path.join(scene_dir, sub)
        if os.path.isdir(p):
            for f in os.listdir(p):
                fp = os.path.join(p, f)
                if os.path.isfile(fp): freed += os.path.getsize(fp)
            shutil.rmtree(p, ignore_errors=True)
        elif os.path.isfile(p):
            freed += os.path.getsize(p); os.remove(p)
    for x in glob.glob(os.path.join(scene_dir, "*.xml")): os.remove(x)
    return freed


def run_decompositions(c2_dir, scene_name):
    import polsartools as pst
    date = get_date(scene_name)
    print(f"\n  Scene: {date}")
    print(f"  C2:    {c2_dir}")
    results = {}

    print("  [1/5] H/Alpha...")
    try:
        pst.h_alpha_dp(c2_dir, win=WINDOW_SIZE, fmt="tif")
        results["Hdp"] = os.path.join(c2_dir, "Hdp.tif")
        results["alphadp"] = os.path.join(c2_dir, "alphadp.tif")
        print("        ✅ Hdp.tif, alphadp.tif")
    except Exception as e: print(f"        ❌ {e}")

    print("  [2/5] DpRVI...")
    try:
        pst.dprvi(c2_dir, win=WINDOW_SIZE, fmt="tif")
        results["dprvi"] = os.path.join(c2_dir, "dprvi.tif")
        print("        ✅ dprvi.tif")
    except Exception as e: print(f"        ❌ {e}")

    print("  [3/5] DOP...")
    try:
        pst.dop_dp(c2_dir, win=WINDOW_SIZE, fmt="tif")
        results["dopdp"] = os.path.join(c2_dir, "dopdp.tif")
        print("        ✅ dopdp.tif")
    except Exception as e: print(f"        ❌ {e}")

    print("  [4/5] RVI...")
    try:
        pst.rvi_dp(c2_dir, win=WINDOW_SIZE, fmt="tif")
        results["rvidp"] = os.path.join(c2_dir, "rvidp.tif")
        print("        ✅ rvidp.tif")
    except Exception as e: print(f"        ❌ {e}")

    print("  [5/5] PRVI...")
    try:
        pst.prvi_dp(c2_dir, win=WINDOW_SIZE, fmt="tif")
        results["prvidp"] = os.path.join(c2_dir, "prvidp.tif")
        print("        ✅ prvidp.tif")
    except Exception as e: print(f"        ❌ {e}")

    return results


def load_raster(path):
    if not path or not os.path.isfile(path): return None, None
    ds = rasterio.open(path)
    data = ds.read(1).astype(np.float32)
    data[data == 0] = np.nan
    return data, ds


def to_db(data):
    """Convert linear power to dB."""
    with np.errstate(divide='ignore', invalid='ignore'):
        db = 10 * np.log10(np.abs(data))
    db[~np.isfinite(db)] = np.nan
    return db


def make_master_chart(c2_dir, scene_name, results):
    """
    ONE single image with ALL maps for this scene.
    Layout: 3 rows × 3 columns = 9 panels
      Row 1: VV σ° (dB), VH σ° (dB), VH/VV ratio (dB)
      Row 2: Entropy (H), Alpha (α°), DpRVI
      Row 3: DOP, RVI, H/α scatter
    """
    date = get_date(scene_name)
    scene_fig_dir = os.path.join(FIGURES_DIR, scene_name)
    os.makedirs(scene_fig_dir, exist_ok=True)

    aoi = None
    try:
        import geopandas as gpd
        if os.path.isfile(AOI_SHP):
            aoi = gpd.read_file(AOI_SHP)
            if aoi.crs and aoi.crs.to_epsg() != 4326: aoi = aoi.to_crs(epsg=4326)
    except: pass

    # Load backscatter: C11 = VV power, C22 = VH power
    c11_data, c11_ds = load_raster(os.path.join(c2_dir, "C11.tif"))
    c22_data, c22_ds = load_raster(os.path.join(c2_dir, "C22.tif"))

    rasters = {}
    if c11_data is not None:
        rasters["VV_dB"] = (to_db(c11_data), c11_ds)
    if c22_data is not None:
        rasters["VH_dB"] = (to_db(c22_data), c22_ds)
    if c11_data is not None and c22_data is not None:
        # VH/VV cross-pol ratio — key indicator for vegetation
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio_db = 10 * np.log10(np.abs(c22_data) / np.abs(c11_data))
        ratio_db[~np.isfinite(ratio_db)] = np.nan
        rasters["VH_VV_dB"] = (ratio_db, c11_ds)

    for key, path in results.items():
        d, ds = load_raster(path)
        if d is not None: rasters[key] = (d, ds)

    # 3 rows × 3 columns
    panels = [
        # Row 1: Backscatter
        ("VV_dB",    "VV Backscatter σ° (dB)",  "gray",    None, None, "dB"),
        ("VH_dB",    "VH Backscatter σ° (dB)",  "gray",    None, None, "dB"),
        ("VH_VV_dB", "VH/VV Ratio (dB)",        "RdBu_r",  -18, -6,   "dB"),
        # Row 2: Decomposition
        ("Hdp",      "Entropy (H)",              "viridis", 0, 1,      "H"),
        ("alphadp",  "Alpha Angle (α°)",         "jet",     0, 90,     "Degrees"),
        ("dprvi",    "DpRVI (Vegetation Index)",  "RdYlGn", 0, 1,      "Index"),
        # Row 3: More indices + scatter
        ("dopdp",    "Degree of Polarization",   "plasma",  0, 1,      "DOP"),
        ("rvidp",    "RVI (Radar Veg. Index)",   "YlGn",    0, 1,      "Index"),
        # Last panel will be H/Alpha scatter (handled separately)
    ]

    nrows, ncols = 3, 3
    fig, axes = plt.subplots(nrows, ncols, figsize=(24, 20))
    fig.suptitle(
        f"Sentinel-1 Dual-Pol Analysis — {date}\n{scene_name[:80]}",
        fontsize=16, fontweight='bold', y=0.995)

    for idx, (key, title, cmap, vmin, vmax, cbar_label) in enumerate(panels):
        ax = axes.flat[idx]
        if key not in rasters:
            ax.text(0.5, 0.5, f"{title}\nN/A", ha='center', va='center',
                    transform=ax.transAxes, fontsize=14, color='gray')
            ax.set_facecolor('#f0f0f0'); ax.set_title(title, fontsize=12); continue

        data, ds = rasters[key]
        extent = [ds.bounds.left, ds.bounds.right, ds.bounds.bottom, ds.bounds.top]
        if vmin is None: vmin = np.nanpercentile(data, 2)
        if vmax is None: vmax = np.nanpercentile(data, 98)

        im = ax.imshow(data, extent=extent, cmap=cmap, vmin=vmin, vmax=vmax,
                       interpolation='nearest', origin='upper')
        if aoi is not None:
            aoi.boundary.plot(ax=ax, color='red', linewidth=1.5, linestyle='--')
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_xlabel("Lon", fontsize=9); ax.set_ylabel("Lat", fontsize=9)
        ax.tick_params(labelsize=8)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.05)
        plt.colorbar(im, cax=cax, label=cbar_label)

    # Last panel (row 2, col 2): H/Alpha scatter plot
    ax = axes.flat[8]
    if "Hdp" in rasters and "alphadp" in rasters:
        h_d = rasters["Hdp"][0]; a_d = rasters["alphadp"][0]
        h_flat = h_d[~np.isnan(h_d)].flatten()
        a_flat = a_d[~np.isnan(a_d)].flatten()
        n = min(len(h_flat), len(a_flat))
        h_flat = h_flat[:n]; a_flat = a_flat[:n]
        if len(h_flat) > 300000:
            idx_s = np.random.choice(len(h_flat), 300000, replace=False)
            h_flat, a_flat = h_flat[idx_s], a_flat[idx_s]
        hb = ax.hexbin(h_flat, a_flat, gridsize=150, cmap='jet', mincnt=1, norm=LogNorm())
        ax.set_xlabel("Entropy (H)", fontsize=10)
        ax.set_ylabel("Alpha (°)", fontsize=10)
        ax.set_title("H/α Feature Space", fontsize=12, fontweight='bold')
        ax.set_xlim(0, 1); ax.set_ylim(0, 90)
        ax.axhline(42.5, color='w', ls='--', alpha=0.6, lw=1)
        ax.axvline(0.5, color='w', ls='--', alpha=0.6, lw=1)
        # Zone labels
        ax.text(0.25, 20, "Surface", ha='center', fontsize=8, color='white', fontweight='bold')
        ax.text(0.25, 65, "Dihedral", ha='center', fontsize=8, color='white', fontweight='bold')
        ax.text(0.75, 20, "Mixed\nSurface", ha='center', fontsize=8, color='white', fontweight='bold')
        ax.text(0.75, 65, "Volume\nScattering", ha='center', fontsize=8, color='white', fontweight='bold')
        plt.colorbar(hb, ax=ax, label="Count", shrink=0.8)
    else:
        ax.text(0.5, 0.5, "H/α Scatter\nN/A", ha='center', va='center',
                transform=ax.transAxes, fontsize=14, color='gray')
        ax.set_facecolor('#f0f0f0')
        ax.set_title("H/α Feature Space", fontsize=12)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    out_path = os.path.join(scene_fig_dir, "master_chart.png")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  📊 Master chart → {out_path}")

    # Save statistics
    stats_path = os.path.join(scene_fig_dir, "statistics.txt")
    with open(stats_path, 'w') as f:
        f.write(f"STATISTICS — {date}\n{'='*50}\n")
        stat_keys = [
            ("VV_dB",    "VV Backscatter (dB)"),
            ("VH_dB",    "VH Backscatter (dB)"),
            ("VH_VV_dB", "VH/VV Ratio (dB)"),
            ("Hdp",      "Entropy (H)"),
            ("alphadp",  "Alpha (°)"),
            ("dprvi",    "DpRVI"),
            ("dopdp",    "DOP"),
            ("rvidp",    "RVI"),
            ("prvidp",   "PRVI"),
        ]
        for key, label in stat_keys:
            if key in rasters:
                d = rasters[key][0]
                v = d[~np.isnan(d)]
                if len(v) > 0:
                    f.write(f"\n{label}:\n")
                    f.write(f"  Mean:   {np.mean(v):.4f}\n")
                    f.write(f"  Std:    {np.std(v):.4f}\n")
                    f.write(f"  Min:    {np.min(v):.4f}\n")
                    f.write(f"  Max:    {np.max(v):.4f}\n")
                    f.write(f"  Median: {np.median(v):.4f}\n")
    print(f"  📝 Stats → {stats_path}")

    for k in list(rasters.keys()):
        try: rasters[k][1].close()
        except: pass


def make_all_scenes_chart(all_results):
    """
    Giant chart: rows = scenes (by date), columns = all parameters.
    Includes VV, VH, VH/VV ratio, H, α, DpRVI, DOP, RVI.
    """
    if len(all_results) < 1: return

    param_keys =   ["VV_dB", "VH_dB", "VH_VV_dB", "Hdp", "alphadp", "dprvi", "dopdp", "rvidp"]
    param_titles = ["VV σ°(dB)", "VH σ°(dB)", "VH/VV(dB)", "Entropy", "Alpha(°)", "DpRVI", "DOP", "RVI"]
    param_cmaps =  ["gray", "gray", "RdBu_r", "viridis", "jet", "RdYlGn", "plasma", "YlGn"]
    param_ranges = [(None,None), (None,None), (-18,-6), (0,1), (0,90), (0,1), (0,1), (0,1)]

    # Map keys to source files
    def get_source(c2_dir, key):
        if key == "VV_dB": return os.path.join(c2_dir, "C11.tif"), True
        elif key == "VH_dB": return os.path.join(c2_dir, "C22.tif"), True
        elif key == "VH_VV_dB": return "ratio", False  # special
        else: return os.path.join(c2_dir, f"{key}.tif"), False

    n_scenes = len(all_results)
    n_params = len(param_keys)

    fig, axes = plt.subplots(n_scenes, n_params, figsize=(3.5*n_params, 3*n_scenes))
    if n_scenes == 1: axes = axes.reshape(1, -1)

    fig.suptitle(
        f"Sentinel-1 Dual-Pol — All {n_scenes} Scenes × All Parameters\nBattambang, Cambodia",
        fontsize=16, fontweight='bold', y=1.01)

    for row, (scene_name, c2_dir, _) in enumerate(all_results):
        date = get_date(scene_name)

        # Pre-load C11 and C22 for this scene (needed for ratio)
        c11_data, c11_ds = load_raster(os.path.join(c2_dir, "C11.tif"))
        c22_data, c22_ds = load_raster(os.path.join(c2_dir, "C22.tif"))

        for col in range(n_params):
            ax = axes[row, col]
            key = param_keys[col]

            data = None; ds = None
            if key == "VV_dB" and c11_data is not None:
                data = to_db(c11_data); ds = c11_ds
            elif key == "VH_dB" and c22_data is not None:
                data = to_db(c22_data); ds = c22_ds
            elif key == "VH_VV_dB" and c11_data is not None and c22_data is not None:
                with np.errstate(divide='ignore', invalid='ignore'):
                    data = 10 * np.log10(np.abs(c22_data) / np.abs(c11_data))
                data[~np.isfinite(data)] = np.nan
                ds = c11_ds
            else:
                fpath = os.path.join(c2_dir, f"{key}.tif")
                data, ds = load_raster(fpath)

            if data is None or ds is None:
                ax.text(0.5, 0.5, "N/A", ha='center', va='center',
                        transform=ax.transAxes, fontsize=9, color='gray')
                ax.set_facecolor('#f5f5f5')
            else:
                extent = [ds.bounds.left, ds.bounds.right, ds.bounds.bottom, ds.bounds.top]
                vmin, vmax = param_ranges[col]
                if vmin is None: vmin = np.nanpercentile(data, 2)
                if vmax is None: vmax = np.nanpercentile(data, 98)
                im = ax.imshow(data, extent=extent, cmap=param_cmaps[col],
                              vmin=vmin, vmax=vmax, interpolation='nearest', origin='upper')
                if col == n_params - 1:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="5%", pad=0.05)
                    plt.colorbar(im, cax=cax)

            if col == 0:
                ax.set_ylabel(date, fontsize=10, fontweight='bold', rotation=0,
                             labelpad=65, va='center')
            else:
                ax.set_ylabel("")
            if row == 0:
                ax.set_title(param_titles[col], fontsize=11, fontweight='bold')
            ax.tick_params(labelsize=5)
            if row < n_scenes - 1: ax.set_xticklabels([])
            if col > 0: ax.set_yticklabels([])

        # Close rasterio handles
        if c11_ds: c11_ds.close()
        if c22_ds: c22_ds.close()

    plt.tight_layout()
    out = os.path.join(FIGURES_DIR, "ALL_SCENES_ALL_PARAMS.png")
    plt.savefig(out, dpi=120, bbox_inches='tight')
    plt.close()
    print(f"\n  📊 ALL-IN-ONE chart → {out}")


def plot_time_series(all_results):
    from datetime import datetime
    if len(all_results) < 2:
        print("  ⚠️ Need ≥2 scenes for time series."); return

    dates = []
    metrics = {k: [] for k in ["VV_dB", "VH_dB", "VH_VV_dB",
                                 "Hdp", "alphadp", "dprvi", "dopdp", "rvidp", "prvidp"]}

    for sn, c2, _ in all_results:
        ds_str = get_date(sn)
        try: dt = datetime.strptime(ds_str, "%Y-%m-%d")
        except: continue
        dates.append(dt)

        # Backscatter
        c11_d, c11_r = load_raster(os.path.join(c2, "C11.tif"))
        c22_d, c22_r = load_raster(os.path.join(c2, "C22.tif"))

        if c11_d is not None:
            vv_db = to_db(c11_d)
            metrics["VV_dB"].append(np.nanmean(vv_db))
        else: metrics["VV_dB"].append(np.nan)

        if c22_d is not None:
            vh_db = to_db(c22_d)
            metrics["VH_dB"].append(np.nanmean(vh_db))
        else: metrics["VH_dB"].append(np.nan)

        if c11_d is not None and c22_d is not None:
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = 10 * np.log10(np.abs(c22_d) / np.abs(c11_d))
            ratio[~np.isfinite(ratio)] = np.nan
            metrics["VH_VV_dB"].append(np.nanmean(ratio))
        else: metrics["VH_VV_dB"].append(np.nan)

        if c11_r: c11_r.close()
        if c22_r: c22_r.close()

        # Decomposition indices
        for k in ["Hdp", "alphadp", "dprvi", "dopdp", "rvidp", "prvidp"]:
            fp = os.path.join(c2, f"{k}.tif")
            if os.path.isfile(fp):
                d, r = load_raster(fp)
                metrics[k].append(np.nanmean(d) if d is not None else np.nan)
                if r: r.close()
            else: metrics[k].append(np.nan)

    if len(dates) < 2: return
    si = np.argsort(dates); dates = [dates[i] for i in si]
    for k in metrics: metrics[k] = [metrics[k][i] for i in si]

    fig, axes = plt.subplots(4, 1, figsize=(14, 16), sharex=True)
    fig.suptitle("Dual-Pol Time Series — Battambang, Cambodia",
                 fontsize=14, fontweight='bold')

    # Panel 1: Backscatter (VV, VH, ratio)
    ax = axes[0]
    ax.plot(dates, metrics["VV_dB"], '-o', color='navy', label='VV σ° (dB)', lw=2, ms=6)
    ax.plot(dates, metrics["VH_dB"], '-s', color='darkred', label='VH σ° (dB)', lw=2, ms=6)
    ax2 = ax.twinx()
    ax2.plot(dates, metrics["VH_VV_dB"], '-^', color='green', label='VH/VV (dB)', lw=2, ms=6)
    ax.set_ylabel("σ° (dB)"); ax2.set_ylabel("VH/VV (dB)", color='green')
    ax.set_title("Backscatter"); ax.grid(alpha=0.3)
    l1, lb1 = ax.get_legend_handles_labels()
    l2, lb2 = ax2.get_legend_handles_labels()
    ax.legend(l1+l2, lb1+lb2, loc='best')

    # Panel 2: Vegetation indices
    ax = axes[1]
    for k, l, c, m in [("dprvi","DpRVI","forestgreen","o"),
                        ("rvidp","RVI","limegreen","s"),
                        ("prvidp","PRVI","teal","^")]:
        if any(not np.isnan(v) for v in metrics[k]):
            ax.plot(dates, metrics[k], f'-{m}', color=c, label=l, lw=2, ms=6)
    ax.set_ylabel("Index"); ax.set_title("Vegetation Indices")
    ax.legend(); ax.grid(alpha=0.3); ax.set_ylim(0, 1)

    # Panel 3: H/Alpha
    ax = axes[2]
    ax.plot(dates, metrics["Hdp"], '-o', color='steelblue', label='Entropy (H)', lw=2, ms=6)
    ax3 = ax.twinx()
    ax3.plot(dates, metrics["alphadp"], '-s', color='darkorange', label='Alpha (°)', lw=2, ms=6)
    ax.set_ylabel("H", color='steelblue'); ax3.set_ylabel("α (°)", color='darkorange')
    ax.set_title("H/Alpha Decomposition")
    l1, lb1 = ax.get_legend_handles_labels()
    l2, lb2 = ax3.get_legend_handles_labels()
    ax.legend(l1+l2, lb1+lb2); ax.grid(alpha=0.3)

    # Panel 4: DOP
    ax = axes[3]
    ax.plot(dates, metrics["dopdp"], '-o', color='purple', label='DOP', lw=2, ms=6)
    ax.set_ylabel("DOP"); ax.set_xlabel("Acquisition Date")
    ax.set_title("Degree of Polarization")
    ax.legend(); ax.grid(alpha=0.3); ax.set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "time_series.png"), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  📊 Time series → {os.path.join(FIGURES_DIR, 'time_series.png')}")

    # CSV with all values
    csv = os.path.join(FIGURES_DIR, "time_series.csv")
    with open(csv, 'w') as f:
        f.write("date,VV_dB,VH_dB,VH_VV_dB,H,alpha,DpRVI,DOP,RVI,PRVI\n")
        for i, dt in enumerate(dates):
            vals = [f"{metrics[k][i]:.4f}" if not np.isnan(metrics[k][i]) else ""
                    for k in ["VV_dB","VH_dB","VH_VV_dB","Hdp","alphadp","dprvi","dopdp","rvidp","prvidp"]]
            f.write(f"{dt.strftime('%Y-%m-%d')},{','.join(vals)}\n")
    print(f"  📝 CSV → {csv}")


def main():
    print("=" * 70)
    print("  STEP 3: DECOMPOSITION + VISUALIZATION + CLEANUP")
    print("=" * 70)

    scenes = find_c2_folders()
    if not scenes:
        print(f"\n❌ No C2 folders in {PREPROCESSED_DIR}"); sys.exit(1)

    print(f"\n{len(scenes)} scene(s):")
    for s, c in scenes: print(f"  📅 {get_date(s)} → {c}")

    all_results = []; total_freed = 0

    for i, (sn, c2) in enumerate(scenes, 1):
        print(f"\n{'='*70}\n  SCENE {i}/{len(scenes)}\n{'='*70}")
        res = run_decompositions(c2, sn)
        if res:
            print(f"\n  Generating master chart...")
            make_master_chart(c2, sn, res)
        freed = cleanup(os.path.join(PREPROCESSED_DIR, sn))
        total_freed += freed
        if freed > 0: print(f"\n  🧹 Freed {freed/(1024**3):.2f} GB")
        all_results.append((sn, c2, res))

    # All-scenes mega chart
    print(f"\n{'='*70}\n  GENERATING ALL-SCENES CHART\n{'='*70}")
    make_all_scenes_chart(all_results)
    plot_time_series(all_results)

    print(f"\n{'='*70}")
    print(f"  🎉 ALL DONE!")
    print(f"{'='*70}")
    print(f"\n  Disk freed: {total_freed/(1024**3):.2f} GB")
    print(f"""
  OUTPUTS:
    📁 {PREPROCESSED_DIR}/<scene>/C2/  → C2 + decomposition TIFs
    📁 {FIGURES_DIR}/<scene>/          → master_chart.png, statistics.txt
    📊 {FIGURES_DIR}/ALL_SCENES_ALL_PARAMS.png  ← THE BIG CHART
    📊 {FIGURES_DIR}/time_series.png
    📝 {FIGURES_DIR}/time_series.csv
""")


if __name__ == "__main__":
    main()


# ==========================================================================
# DECOMPOSITION TECHNIQUES — DOCUMENTATION
# ==========================================================================
#
# ══════════════════════════════════════════════════════════════════════════
# WHAT DECOMPOSITION TECHNIQUES DID WE USE?
# ══════════════════════════════════════════════════════════════════════════
#
# 5 dual-pol decomposition methods from polsartools, plus 3 backscatter
# products derived directly from the C2 matrix:
#
# BACKSCATTER (from C2 matrix, no decomposition needed):
# ──────────────────────────────────────────────────────
#   VV σ° (dB)   = 10·log10(C11)  — Co-pol backscatter
#                   Sensitive to surface roughness, soil moisture,
#                   and structure geometry. Dominant scattering.
#
#   VH σ° (dB)   = 10·log10(C22)  — Cross-pol backscatter
#                   Sensitive to VOLUME scattering → vegetation biomass,
#                   canopy structure, crop height. THE key channel for
#                   agriculture monitoring.
#
#   VH/VV (dB)   = 10·log10(C22/C11) — Cross-pol ratio
#                   Normalizes for incidence angle effects.
#                   More negative = surface-dominated (bare soil ~ -12 to -15 dB)
#                   Less negative = volume-dominated (dense veg ~ -6 to -10 dB)
#
# DECOMPOSITION INDICES:
# ──────────────────────
# 1. H/Alpha (Cloude-Pottier dual-pol)
#    - Entropy H: 0=single mechanism, 1=random
#    - Alpha α: 0°=surface, 45°=volume, 90°=dihedral
#    - THE most informative for physical interpretation
#
# 2. DpRVI — Dual-pol Radar Vegetation Index
#    - 0 (bare) → 1 (dense vegetation)
#    - BEST single index for crop monitoring with S1
#    - Reference: Mandal et al. (2020)
#
# 3. DOP — Barakat Degree of Polarization
#    - 0 (unpolarized) → 1 (fully polarized)
#    - Low = volume scattering, High = surface
#
# 4. RVI — Radar Vegetation Index
#    - 4·σ°_VH / (σ°_VV + σ°_VH), simple ratio
#
# 5. PRVI — Polarimetric RVI
#    - (1-DOP)·RVI, improved noise rejection
#
# ══════════════════════════════════════════════════════════════════════════
# WHY VH IS CRITICAL (and was missing before)
# ══════════════════════════════════════════════════════════════════════════
#
# VV (co-pol): Dominated by surface/double-bounce scattering.
#   → Responds to soil moisture, surface roughness, water bodies
#
# VH (cross-pol): Requires depolarization → VOLUME scattering
#   → Directly proportional to vegetation biomass and canopy density
#   → THE primary channel for crop monitoring, forest mapping
#   → Changes much more dramatically during crop growth cycle
#
# VH/VV ratio: Removes incidence angle dependency
#   → More robust for comparing across scenes
#   → Standard metric in agricultural SAR studies
#
# ══════════════════════════════════════════════════════════════════════════
# OTHER DUAL-POL OPTIONS IN POLSARTOOLS (not used):
#   shannon_h_dp, dprvic, dp_desc, dprbi/dprbic, dprsi/dprsic,
#   powers_dp, powers_dp_grd
#
# FULL-POL (if available): Yamaguchi, Freeman-Durden, H/A/Alpha, etc.
# COMPACT-POL: m-chi, m-delta, CpRVI, MF3CC
# ══════════════════════════════════════════════════════════════════════════