# 🛰️ Sentinel-1 Dual-Pol Polarimetric Decomposition Pipeline

**Automated SAR polarimetric analysis for agricultural monitoring using Sentinel-1 SLC data**

This pipeline downloads Sentinel-1 IW SLC (Single Look Complex) data, processes all 3 sub-swaths (IW1+IW2+IW3) into geocoded C2 covariance matrices, runs 5 dual-pol polarimetric decomposition indices, and generates publication-ready visualizations and time series charts.

Developed for **crop monitoring** 

![Master Chart Example](images/master_chart.png)

---

## 📋 Table of Contents

- [What This Does](#-what-this-does)
- [Pipeline Overview](#-pipeline-overview)
- [Installation](#-installation)
- [Data Preparation](#-data-preparation)
- [Usage](#-usage)
- [Outputs](#-outputs)
- [Understanding the Results](#-understanding-the-results)
- [Troubleshooting](#-troubleshooting)

---

## 🔍 What This Does

This pipeline takes **raw Sentinel-1 radar satellite images** and extracts meaningful information about what's on the ground — crops, bare soil, water, forests, buildings — using **polarimetric decomposition**.

In simple terms:
1. The satellite sends a microwave pulse (VV) and listens for echoes in two orientations (VV and VH)
2. We capture the **power** of each echo AND the **phase relationship** between them
3. We store this in a 2×2 table called the **C2 covariance matrix** (4 numbers per pixel)
4. We mathematically decompose this table into physically meaningful parameters:
   - **How mixed is the scattering?** (Entropy)
   - **What type of scattering dominates?** (Alpha angle)
   - **How much vegetation is there?** (DpRVI, RVI, PRVI)
   - **How organized is the signal?** (Degree of Polarization)

---

## 🔄 Pipeline Overview

```
Step 1: Download Sentinel-1 SLC data (ASF Vertex / Copernicus Hub)
           ↓
Step 2: Preprocess SLC → C2 matrix (SNAP + GDAL)
           │
           ├── IW1: Split → Calibrate(complex) → Deburst → C2 → Multilook → Speckle Filter
           ├── IW2: Split → Calibrate(complex) → Deburst → C2 → Multilook → Speckle Filter
           ├── IW3: Split → Calibrate(complex) → Deburst → C2 → Multilook → Speckle Filter
           │         ↓
           ├── Geocode each sub-swath with GDAL (tie-point grids → GCPs → gdalwarp TPS)
           │         ↓
           └── Mosaic IW1+IW2+IW3 → clip to AOI → C11.tif, C12_real.tif, C12_imag.tif, C22.tif
           ↓
Step 3: Decomposition + Visualization (polsartools + matplotlib)
           │
           ├── H/Alpha decomposition (Entropy + Alpha angle)
           ├── DpRVI (Dual-pol Radar Vegetation Index)
           ├── DOP (Degree of Polarization)
           ├── RVI (Radar Vegetation Index)
           ├── PRVI (Polarimetric RVI)
           │         ↓
           ├── Master chart per scene (VV, VH, VH/VV, H, α, DpRVI, DOP, RVI, H/α scatter)
           ├── All-scenes mega chart (rows = dates, columns = parameters)
           ├── Time series plot + CSV
           └── Cleanup intermediate files
```

---

## 🛠️ Installation

### 1. Create Conda Environment

```bash
conda create -n sarpolsar -c conda-forge python=3.11 gdal "numpy<2" rasterio h5py geopandas cartopy bokeh matplotlib jupyter eccodes pygrib -y
```

```bash
conda activate sarpolsar
```

### 2. Install ESA SNAP (Sentinel Application Platform)

SNAP provides the SAR-specific processing operators (orbit correction, calibration, deburst, polarimetric matrices, speckle filtering).

```bash
# Download SNAP 13.0 installer
wget https://download.esa.int/step/snap/13.0/installers/esa-snap_all_linux-13.0.0.sh

# Make it executable
chmod +x esa-snap_all_linux-13.0.0.sh

# Run installer (follow prompts, install to ~/esa-snap/)
./esa-snap_all_linux-13.0.0.sh
```

After installation, **build the Graph Processing Tool (GPT)**:

```bash
# Update SNAP modules (first run takes a few minutes)
~/esa-snap/bin/snap --nosplash --nogui --modules --update-all
```

Verify GPT is working:

```bash
~/esa-snap/bin/gpt -h
```

You should see the help text for the Graph Processing Tool.

### 3. Install Python Packages

```bash
# MintPy (InSAR time series — optional but useful)
pip install git+https://github.com/insarlab/MintPy.git

# Jupyter widgets
pip install opensarlab_lib ipyfilechooser ipympl ipywidgets jupyterlab-widgets widgetsnbextension

# pyroSAR (SAR data handling utilities)
pip install pyroSAR

# polsartools (polarimetric decomposition — THE key package)
pip install polsartools
```

### 4. Verify Installation

```bash
conda activate sarpolsar
python -c "
import polsartools as pst
from osgeo import gdal
import rasterio
import geopandas
print('✅ All packages installed successfully!')
print(f'   polsartools version: {pst.__version__}')
print(f'   GDAL version: {gdal.__version__}')
"
```

---

## 📁 Data Preparation

### Folder Structure

Create this folder structure before running:

```
~/Desktop/Muddasir/          (or your working directory)
├── downloads/               ← Put Sentinel-1 SLC .zip files here
│   ├── S1A_IW_SLC__1SDV_20251107T112026_....zip
│   ├── S1A_IW_SLC__1SDV_20251119T112026_....zip
│   └── ...
├── aoi/                     ← Put your AOI shapefile here
│   ├── Battambang.shp
│   ├── Battambang.shx
│   ├── Battambang.dbf
│   └── Battambang.prj
├── step2_preprocess_slc.py
└── step3_decompose.py
```

### Download Sentinel-1 SLC Data

1. Go to [ASF Vertex](https://search.asf.alaska.edu/) or [Copernicus Browser](https://dataspace.copernicus.eu/)
2. Search for your AOI
3. Filter:
   - **Mission:** Sentinel-1
   - **Product Type:** SLC (Single Look Complex) — ⚠️ NOT GRD!
   - **Beam Mode:** IW (Interferometric Wide)
   - **Polarization:** VV+VH (Dual-pol)
4. Download .zip files into the `downloads/` folder

> **⚠️ Important:** You MUST use **SLC** data, not GRD. GRD data loses the complex phase information needed for polarimetric decomposition. The C12 (cross-correlation) component requires the phase relationship between VV and VH channels.

### AOI Shapefile

Your AOI shapefile should be in **EPSG:4326 (WGS84)** or any CRS (the script reprojects automatically). The script buffers the AOI by 0.05° on all sides to ensure complete coverage.

---

## 🚀 Usage

### Step 1: Download Data

Download Sentinel-1 SLC .zip files manually from ASF Vertex or Copernicus and place them in the `downloads/` folder.

### Step 2: Preprocess SLC → C2 Matrix

```bash
conda activate sarpolsar
python step2_preprocess_slc.py
```

**What it does:**
- Processes each of the 3 sub-swaths (IW1, IW2, IW3) separately through SNAP
- Geocodes each sub-swath using GDAL (reads tie-point grids, builds GCPs, warps with TPS)
- Mosaics the 3 sub-swaths into one seamless GeoTIFF covering the full ~250km swath
- Clips to your AOI
- Auto-skips already completed scenes on re-run

**Time:** ~30-90 minutes per scene (3 × 10-30 min for SNAP + 5 min for geocoding/mosaic)

**RAM:** 8GB minimum

### Step 3: Decomposition + Visualization

```bash
python step3_decompose.py
```

**What it does:**
- Runs 5 decomposition indices using polsartools
- Generates master chart per scene (3×3 grid with all parameters)
- Generates all-scenes mega chart
- Generates time series plot and CSV
- Cleans up intermediate SNAP files to free disk space

**Time:** ~5 minutes per scene

---

## 📦 Outputs

### Per Scene

```
outputs/preprocessed/<scene_name>/C2/
├── C11.tif           ← VV power (co-pol intensity)
├── C12_real.tif      ← VV-VH cross-correlation (real part)
├── C12_imag.tif      ← VV-VH cross-correlation (imaginary part)
├── C22.tif           ← VH power (cross-pol intensity)
├── Hdp.tif           ← Entropy (0 = pure scattering, 1 = random)
├── alphadp.tif       ← Alpha angle (0° = surface, 45° = volume, 90° = dihedral)
├── e1_norm.tif       ← Normalized eigenvalue 1
├── e2_norm.tif       ← Normalized eigenvalue 2
├── dprvi.tif         ← Dual-pol Radar Vegetation Index (0-1)
├── dopdp.tif         ← Degree of Polarization (0-1)
├── rvidp.tif         ← Radar Vegetation Index (0-1)
└── prvidp.tif        ← Polarimetric RVI (0-1)
```

### Figures

```
outputs/figures/
├── <scene_name>/
│   ├── master_chart.png     ← 9-panel overview (VV, VH, ratio, H, α, DpRVI, DOP, RVI, H/α scatter)
│   └── statistics.txt       ← Mean, std, min, max for all parameters
├── ALL_SCENES_ALL_PARAMS.png ← Giant grid (rows = dates, cols = all parameters)
├── time_series.png           ← 4-panel temporal analysis
└── time_series.csv           ← Machine-readable time series data
```

---

## 🧠 Understanding the Results

### Backscatter (from C2 matrix directly)

| Product | Source | What It Measures |
|---------|--------|-----------------|
| **VV σ° (dB)** | C11 | Co-pol backscatter — surface roughness, soil moisture, structure |
| **VH σ° (dB)** | C22 | Cross-pol backscatter — **vegetation biomass, canopy structure** |
| **VH/VV (dB)** | C22/C11 | Cross-pol ratio — vegetation density, normalized for incidence angle |

### Decomposition Products (from eigenvalue decomposition of C2)

| Product | Range | Low Value Means | High Value Means |
|---------|-------|----------------|-----------------|
| **Entropy (H)** | 0 → 1 | Single scattering mechanism (bare soil, water) | Random mixed scattering (dense vegetation, forest) |
| **Alpha (α)** | 0° → 90° | Surface scattering (flat ground, water) | Double-bounce (buildings, flooded forest) |
| **DpRVI** | 0 → 1 | No vegetation (bare soil, water) | Dense vegetation (crops, forest) |
| **DOP** | 0 → 1 | Unpolarized / chaotic (vegetation) | Fully polarized / organized (surface, man-made) |
| **RVI** | 0 → 1 | Low vegetation | High vegetation |
| **PRVI** | 0 → 1 | Low vegetation (DOP-weighted) | High vegetation (DOP-weighted) |

### Reading the H/α Scatter Plot

```
Alpha (°)
  90° ┌─────────────────┬─────────────────┐
      │   Dihedral /    │    Random       │
      │   Double-bounce │    Dihedral     │
  45° ├─────────────────┼─────────────────┤
      │   Surface       │    Volume       │
      │   Scattering    │    Scattering   │
   0° └─────────────────┴─────────────────┘
      0                0.5                1
                  Entropy (H)
```

### Agricultural Interpretation

For crop monitoring (like rice paddies in Cambodia):

- **Growing season:** DpRVI rises, Entropy rises, Alpha moves toward 45°, VH gets brighter
- **Harvest / fallow:** DpRVI drops, Entropy drops, Alpha moves toward 0°, VH gets darker
- **Flooded fields:** Very low VH, high VV, possible double-bounce (α → 90°) if crops standing in water

---

## 🔧 Troubleshooting

### SNAP GPT "/ by zero" Error
This is a known SNAP bug with the Terrain-Correction operator on polarimetric data. **Our pipeline avoids this entirely** by using GDAL for geocoding instead of SNAP's TC operator.

### "Cannot construct DataBuffer" Error
SNAP runs out of memory when processing all 3 sub-swaths in a single graph. **Our pipeline processes each sub-swath separately** then mosaics with GDAL.

### GDAL "no GeoTransform and no GCPs" Error
Individual `.img` files in BEAM-DIMAP have no georeferencing. **Our pipeline reads tie-point grids** from the `.data/tie_point_grids/` folder and builds GCPs from the `.dim` XML metadata.

### Out of Memory
Increase SNAP's memory allocation in `step2_preprocess_slc.py`:
```python
GPT_MEMORY = "12G"  # or "16G" if you have enough RAM
```

### Missing Orbit Files
SNAP auto-downloads precise orbit files. If it fails, set `continueOnFail` to `true` (already set in our graph) to fall back to predicted orbits.

### Scene Doesn't Cover AOI
Check which sub-swaths overlap your AOI. The pipeline is fault-tolerant — if one sub-swath fails, it still mosaics the other two.

---

## 📚 References

- Cloude, S.R. & Pottier, E. (1997). An entropy based classification scheme for land applications of polarimetric SAR. IEEE TGRS.
- Mandal, D. et al. (2020). Dual polarimetric radar vegetation index for crop growth monitoring. IEEE GRSL.
- Kim, Y. & van Zyl, J. (2009). A time-series approach to estimate soil moisture using polarimetric radar data. IEEE TGRS.
- [polsartools documentation](https://github.com/polsartools/polsartools)
- [ESA SNAP Toolbox](https://step.esa.int/main/toolboxes/snap/)

---

## 📄 License

This project is for research and educational purposes.
Sentinel-1 data is freely available from ESA/Copernicus.
