#!/usr/bin/env python3
"""
STEP 1: Verify Setup - Environment, Packages, Data, AOI
========================================================
Run this FIRST to make sure everything is in place before processing.

Usage (from terminal with intnotebook activated):
    conda activate intnotebook
    python step1_verify_setup.py
"""

import os
import sys
import glob

print("=" * 70)
print("  SENTINEL-1 SLC DUAL-POL DECOMPOSITION - SETUP VERIFICATION")
print("=" * 70)

# -------------------------------------------------------------------------
# 1. Check Python & Conda environment
# -------------------------------------------------------------------------
print("\n[1/6] PYTHON & ENVIRONMENT")
print(f"  Python path   : {sys.executable}")
print(f"  Python version: {sys.version}")

conda_env = os.environ.get("CONDA_DEFAULT_ENV", "UNKNOWN")
print(f"  Conda env     : {conda_env}")
if conda_env != "intnotebook":
    print("  ⚠️  WARNING: You are NOT in 'intnotebook' environment!")
    print("  Run: conda activate intnotebook")
else:
    print("  ✅ Correct environment active.")

# -------------------------------------------------------------------------
# 2. Check critical packages
# -------------------------------------------------------------------------
print("\n[2/6] PACKAGE CHECKS")
packages = {
    "pyroSAR": "pyroSAR",
    "polsartools": "polsartools",
    "GDAL (osgeo)": "osgeo",
    "numpy": "numpy",
    "matplotlib": "matplotlib",
    "geopandas": "geopandas",
    "rasterio": "rasterio",
    "fiona": "fiona",
}

all_ok = True
for display_name, import_name in packages.items():
    try:
        mod = __import__(import_name)
        ver = getattr(mod, "__version__", "version unknown")
        print(f"  ✅ {display_name:20s} → {ver}")
    except ImportError:
        print(f"  ❌ {display_name:20s} → NOT INSTALLED")
        all_ok = False

if not all_ok:
    print("\n  ⚠️  Some packages are missing. Install them with:")
    print("    conda install -c conda-forge gdal rasterio geopandas fiona matplotlib")
    print("    pip install polsartools pyroSAR")

# -------------------------------------------------------------------------
# 3. Check SNAP / GPT availability (required for SLC preprocessing)
# -------------------------------------------------------------------------
print("\n[3/6] SNAP GPT CHECK")
import shutil

gpt_path = shutil.which("gpt")
if gpt_path:
    print(f"  ✅ GPT found at: {gpt_path}")
else:
    # Common install locations
    common_paths = [
        os.path.expanduser("~/snap/bin/gpt"),
        os.path.expanduser("~/esa-snap/bin/gpt"),
        "/usr/local/snap/bin/gpt",
        "/opt/snap/bin/gpt",
        "C:\\Program Files\\snap\\bin\\gpt.exe",
    ]
    found = False
    for p in common_paths:
        if os.path.isfile(p):
            print(f"  ✅ GPT found at: {p}")
            print(f"     (Not on PATH — will use this path directly)")
            gpt_path = p
            found = True
            break
    if not found:
        print("  ❌ GPT (SNAP Graph Processing Tool) NOT FOUND!")
        print("     SNAP is REQUIRED to preprocess Sentinel-1 SLC data.")
        print("     Download: https://step.esa.int/main/download/snap-download/")
        print("     After install, add to PATH or note the install location.")

# -------------------------------------------------------------------------
# 4. Check Sentinel-1 SLC zip files in ./downloads/
# -------------------------------------------------------------------------
print("\n[4/6] SENTINEL-1 SLC DATA")
downloads_dir = os.path.join(os.getcwd(), "downloads")

if not os.path.isdir(downloads_dir):
    print(f"  ❌ 'downloads' folder NOT FOUND at: {downloads_dir}")
    print("     Create it and place your S1 SLC .zip files inside.")
else:
    slc_zips = sorted(glob.glob(os.path.join(downloads_dir, "S1*.zip")))
    if len(slc_zips) == 0:
        # Also check for any .zip files
        slc_zips = sorted(glob.glob(os.path.join(downloads_dir, "*.zip")))

    if len(slc_zips) == 0:
        print(f"  ❌ No .zip files found in: {downloads_dir}")
    else:
        print(f"  ✅ Found {len(slc_zips)} SLC zip file(s):")
        for i, z in enumerate(slc_zips, 1):
            size_mb = os.path.getsize(z) / (1024 * 1024)
            print(f"     {i}. {os.path.basename(z)}  ({size_mb:.1f} MB)")

# -------------------------------------------------------------------------
# 5. Check AOI shapefile
# -------------------------------------------------------------------------
print("\n[5/6] AOI SHAPEFILE")
aoi_dir = os.path.join(os.getcwd(), "aoi")
shp_path = os.path.join(aoi_dir, "Battambang.shp")

if not os.path.isdir(aoi_dir):
    print(f"  ❌ 'aoi' folder NOT FOUND at: {aoi_dir}")
else:
    if os.path.isfile(shp_path):
        print(f"  ✅ Shapefile found: {shp_path}")
        # Check companion files
        companions = [".shx", ".dbf", ".prj"]
        for ext in companions:
            comp_file = shp_path.replace(".shp", ext)
            if os.path.isfile(comp_file):
                print(f"     ✅ {os.path.basename(comp_file)}")
            else:
                print(f"     ⚠️  Missing: {os.path.basename(comp_file)}")
        
        # Try to read and display AOI bounds
        try:
            import geopandas as gpd
            aoi = gpd.read_file(shp_path)
            print(f"     CRS     : {aoi.crs}")
            print(f"     Bounds  : {aoi.total_bounds}")
            print(f"     Features: {len(aoi)}")
        except Exception as e:
            print(f"     ⚠️  Could not read shapefile: {e}")
    else:
        shp_files = glob.glob(os.path.join(aoi_dir, "*.shp"))
        if shp_files:
            print(f"  ⚠️  Expected 'Battambang.shp' but found: {[os.path.basename(f) for f in shp_files]}")
        else:
            print(f"  ❌ No .shp files found in: {aoi_dir}")

# -------------------------------------------------------------------------
# 6. Create output directories
# -------------------------------------------------------------------------
print("\n[6/6] OUTPUT DIRECTORIES")
dirs_to_create = [
    os.path.join(os.getcwd(), "outputs"),
    os.path.join(os.getcwd(), "outputs", "preprocessed"),
    os.path.join(os.getcwd(), "outputs", "decomposition"),
    os.path.join(os.getcwd(), "outputs", "figures"),
]
for d in dirs_to_create:
    os.makedirs(d, exist_ok=True)
    print(f"  ✅ {d}")

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------
print("\n" + "=" * 70)
print("  VERIFICATION COMPLETE")
print("=" * 70)
print("""
NEXT STEP:
  If all checks passed, run:
    python step2_preprocess_slc.py

  If SNAP/GPT is missing, install it first from:
    https://step.esa.int/main/download/snap-download/
""")