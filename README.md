# BPGA — A Shiny App for Basic Population Genetic Analysis
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16953021.svg)](https://doi.org/10.5281/zenodo.16953021)
<br>BPGA is an interactive **R/Shiny** application that lets you run core population-genetic workflows on **PLINK** datasets and explore results visually. It supports **PCA**, **FST**, and **ADMIXTURE**, can **merge user data with worldwide reference panels**, and exports publication-ready figures and an interactive Manhattan plot.
<br><br>This Shiny app provides an interactive platform for visualizing population genetic structure. Users can upload their own datasets or use example files to explore population relationships, genetic diversity, and ancestral components. The app dynamically reads PLINK binary files and integrates them with worldwide reference populations from the **1000 Genomes Project (1000G)** and the **Human Genome Diversity Project (HGDP)**, enabling basic population analyses such as **Principal Component Analysis (PCA)**, **ADMIXTURE**, and **FST** analyses. The app performs PCA for dimensionality reduction, producing scatterplots that reveal population clusters. ADMIXTURE analysis is employed to assess population ancestry components, with a geographic visualization feature that displays these components as pie charts on a world map—each slice representing an ancestral proportion, and chart sizes scaled according to population size. Finally, FST analysis enables users to quantify genetic differentiation between pairs of populations. Designed for both researchers and students in population genetics, the app offers a user-friendly interface for exploring genetic structure through multiple methods. It facilitates interactive data exploration and produces publication-ready plots. Overall, the app integrates PCA, ADMIXTURE, and FST analyses, providing a comprehensive view of population dynamics and evolutionary relationships.

---

## Table of Contents
- [Highlights](#highlights)
- [What the app does](#what-the-app-does)
- [Input data & formats](#input-data--formats)
- [Outputs](#outputs)
- [Quick start](#quick-start)
- [External binaries](#external-binaries)
- [Folder structure](#folder-structure)
- [How to use (workflow)](#how-to-use-workflow)
- [Deployment notes](#deployment-notes)
- [Troubleshooting](#troubleshooting)
- [Acknowledgments](#acknowledgments)
- [References](#references)
- [License](#license)

---

## Highlights
- ✔️ Load **example** datasets or your **own PLINK** data.
- ✔️ Optional **merge** with a worldwide reference dataset (1000G/HGDP-based) that is downloaded and decompressed automatically.
- ✔️ **PCA** (with pruning, MAF filter), **FST** (pairwise or multi-group), and **ADMIXTURE** (K=1..10, with CV errors).
- ✔️ Interactive **Manhattan plot** (FST) with **dbSNP** and **UCSC** links.
- ✔️ **Leaflet map** of ancestry components (pie charts per population) + PNG export.
- ✔️ **Download all plots** as a single **.zip**.
- ✔️ Works with **temporary session folders** (no persistent writes unless you download).

---

## What the app does
- **Data preparation**
  - Accepts user `.zip` (or `.gz`) bundles containing PLINK **.bed/.bim/.fam**.
  - Filters to autosomes and standardizes sample labels.
  - Optionally **merges** with a bundled **Reference** panel (downloaded on demand).

- **Analyses**
  - **PCA via PLINK** with user-set **MAF** and LD-pruning (`--indep-pairwise`), 10 PCs.
  - **FST via GCTA** on selected (sub)populations, plus an **interactive Manhattan plot**.
  - **ADMIXTURE** runs at single K or K-range, with **CV error** summaries.

- **Visualization**
  - PCA scatter (colored by population or subpopulation).
  - ADMIXTURE barplot and boxplot (sortable by component).
  - **Leaflet** map with per-population **pie charts** of ancestry fractions.

---

## Input data & formats
- **User input (compressed)**: `.zip` or `.gz` containing **.bed/.bim/.fam**.
- **FAM population coding** (first column) — choose one when loading (see examples):
  1. **First column preformatted as `POP1_POP2`** ("pop1pop2") → you’ll provide a coordinate file
  2. **First column describes one single population** (“unic”) → you’ll provide one 3-letter code (e.g., `USR`) and optional **LAT/LON**
  3. **First column describes several populations** (“subpop”) → you’ll provide one 3-letter **POP1** code and a coordinate file
- **Coordinate file (needed for mapping when using `pop1pop2` or `subpop`)**  
  Tab-separated text (.txt) with **five columns** (no header):
  
  POP1  POP2  Pop_name  LAT(POP2)  LON(POP2)
---
<h4>Examples of FAM file and assigned coordinate file</h4>

<table>
  <tr>
    <td valign="top">
      <strong>First column preformatted as <code>POP1_POP2</code></strong>
      <table border="1">
        <thead>
          <tr>
            <th>Family ID</th>
            <th>Individual ID</th>
            <th>Father</th>
            <th>Mother</th>
            <th>Sex</th>
            <th>Phenotype</th>
          </tr>
        </thead>
        <tbody>
          <tr><td>AFR_ASW</td><td>NA19916</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
          <tr><td>AFR_ASW</td><td>NA19703</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
          <tr><td>EUR_CEU</td><td>NA12341</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
          <tr><td>EUR_CEU</td><td>NA06984</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
          <tr><td>EAS_CHB</td><td>NA18532</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
          <tr><td>EAS_CHB</td><td>NA18561</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
        </tbody>
      </table>
    </td>
    <td valign="top">
      <strong>Assigned coordinate file (.txt)</strong>
      <table border="1">
        <thead>
          <tr>
            <th>POP1</th>
            <th>POP2</th>
            <th>Pop_name</th>
            <th>LAT</th>
            <th>LON</th>
          </tr>
        </thead>
        <tbody>
          <tr><td>AFR</td><td>ASW</td><td>African</td><td>-3.82</td><td>12.93</td></tr>
          <tr><td>EAS</td><td>CHD</td><td>Chinese</td><td>40.00</td><td>115.00</td></tr>
          <tr><td>EUR</td><td>CEU</td><td>European</td><td>51.30</td><td>11.66</td></tr>
        </tbody>
      </table>
    </td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top">
      <strong>First column describes one single population</strong>
      <table border="1">
        <thead>
          <tr>
            <th>Family ID</th>
            <th>Individual ID</th>
            <th>Father</th>
            <th>Mother</th>
            <th>Sex</th>
            <th>Phenotype</th>
          </tr>
        </thead>
        <tbody>
          <tr><td>AFR</td><td>NA19818</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
          <tr><td>AFR</td><td>NA20346</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
          <tr><td>AFR</td><td>NA19921</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
          <tr><td>AFR</td><td>NA20281</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
          <tr><td>AFR</td><td>NA20301</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
          <tr><td>AFR</td><td>NA20294</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
          <tr><td>AFR</td><td>NA20357</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
        </tbody>
      </table>
    </td>
    <td valign="top">
      <strong>Assigned coordinate file</strong>
      <table border="1">
        <thead>
          <tr>
            <th>Instructions</th>
          </tr>
        </thead>
        <tbody>
          <tr><td>User assigned population code</td></tr>
          <tr><td>User assigned LON, LAT coordinates</td></tr>
        </tbody>
      </table>
    </td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top">
      <strong>First column describes several populations</strong>
      <table border="1">
        <thead>
          <tr>
            <th>Family ID</th>
            <th>Individual ID</th>
            <th>Father</th>
            <th>Mother</th>
            <th>Sex</th>
            <th>Phenotype</th>
          </tr>
        </thead>
        <tbody>
          <tr><td>ASW</td><td>NA19916</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
          <tr><td>ASW</td><td>NA19703</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
          <tr><td>CEU</td><td>NA12341</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
          <tr><td>CEU</td><td>NA06984</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
          <tr><td>CHB</td><td>NA18532</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
          <tr><td>CHB</td><td>NA18561</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
        </tbody>
      </table>
    </td>
    <td valign="top">
      <strong>Assigned coordinate file (.txt)</strong>
      <table border="1">
        <thead>
          <tr>
            <th>POP1</th>
            <th>POP2</th>
            <th>Pop_name</th>
            <th>LAT</th>
            <th>LON</th>
          </tr>
        </thead>
        <tbody>
          <tr><td>ASW</td><td>ASW</td><td>Africa_sample</td><td>9.3</td><td>19.3</td></tr>
          <tr><td>CEU</td><td>CEU</td><td>Europe_sample</td><td>50</td><td>15</td></tr>
          <tr><td>CHB</td><td>CHB</td><td>Asia_sample</td><td>38</td><td>83</td></tr>
        </tbody>
      </table>
    </td>
  </tr>
</table>

---
 
## Outputs
All figures are written to the session’s temp folder (shown in the sidebar) under **`plots/`**:
- `PCA.png`
- `ADM.png` (ADMIXTURE barplot)
- `box_ADM.png` (ADMIXTURE boxplot)
- `map_ADM.png` (Leaflet snapshot via mapshot)
- `FST_interactive.html` (interactive Plotly Manhattan plot)

A **Download** button bundles everything into a **.zip**.

> Note: Session folders are temporary. **Download your results** before you reload/close the session.

---

## Quick start

### 1) Install R packages
```r
install.packages(c(
  "shiny","qqman","readr","data.table","ggplot2","shinyjs","dplyr","tidyr",
  "rmarkdown","knitr","ggforce","mapplots","maps","ggh4x","conflicted",
  "shinycssloaders","sf","leaflet","htmlwidgets","leaflet.minicharts",
  "mapview","plotly","Cairo"
))
```
> `sf` requires GDAL/GEOS/PROJ on your system (install via apt, brew, etc.). For map export, also install `webshot2` if needed.

### 2) Place required binaries and assets
See [External binaries](#external-binaries) and [Folder structure](#folder-structure). Ensure executables in `www/` are **present** and **executable**.

### 3) Run the app
```r
# from the project root
shiny::runApp(".")
```

**Upload size limit:** default is **600 MB** (set in `options(shiny.maxRequestSize)`); adjust if needed.

---

## External binaries
The app calls command-line tools via `system()/system2()`:

- **PLINK 1.9** → `www/plink19`
- **GCTA** → `www/gcta64`
- **ADMIXTURE** → `www/admixture`

**You must** provide platform-appropriate binaries in the `www/` folder and make them executable:
```bash
chmod +x www/plink19 www/gcta64 www/admixture
```

Also ensure `unzip` and `gzip`/`gunzip` are available on your system path.

---

## Folder structure
Expected key files and folders:
```
.
├─ app.R
├─ example/
│  ├─ example_single.zip
│  └─ example_multi.zip
├─ www/
│  ├─ plink19*         # executable (not provided, see Note)
│  ├─ gcta64*          # executable (not provided, see Note)
│  ├─ admixture*       # executable (not provided, see Note)
│  ├─ coord.txt        # reference POP coordinates for mapping
│  ├─ fons.jpg         # background image
│  ├─ workflow.jpg     # splash/workflow image
│  └─ infograph.png    # splash image
```
Temporary session folders like `user_output_xxx/` (with `pca/`, `admx/`, `fst/`, `merged/`, `plots/`, etc.) are created under your R session’s `tempdir()`.

---

## How to use (workflow)

1. **Start**
   - Click **Start analysis**. A temp workspace is created and shown in the sidebar.
   - Optional: **Reload session** to start fresh (clears temp workspace).

2. **Load data**
   - Choose **Example single/multi** or **User file**.
   - For **User file**, upload a `.zip` or `.gz` with `.bed/.bim/.fam`.
   - Choose your **FAM format**:
     - `POP1_POP2` preformatted, **or**
     - **unic** (single population): set `POP1/POP2` 3-letter code and (optionally) **LAT/LON**, **or**
     - **subpop** (multiple): supply a **coordinate file**.

3. **Process input & (optional) Merge**
   - **Only loaded population** (`usr`) → click **Decompress file**.
   - **Merge with Worldwide populations** (`wwp`) → the app downloads & decompresses the reference panel, then click **Merge files**.

4. **PCA**
   - Select POP1 groups to analyze.
   - Set **MAF** for pruning and click **Run PCA**.
   - A PCA plot will appear; the log snippet shows sample/variant counts.

5. **FST**
   - Select one or more **POP2** groups and click **Run FST**.
   - Explore the **interactive Manhattan plot**; click to open **dbSNP** and **UCSC** for the nearest SNP.

6. **ADMIXTURE**
   - Set **MAF**, choose **K** (single or range), click **Run ADMIXTURE**.
   - Review **CV errors** printed in the sidebar.
   - See **barplot**, **boxplot**, and **map** of ancestry components.

7. **Export**
   - Click **Download plots (.zip)** to get all figures (and the interactive HTML).

---

## Deployment notes

### Shiny Server (Linux)
- Copy the project folder to `/srv/shiny-server/BPGA` (example path) and set permissions:
  ```bash
  sudo chown -R shiny:shiny /srv/shiny-server/BPGA
  sudo chmod +x /srv/shiny-server/BPGA/www/plink19 /srv/shiny-server/BPGA/www/gcta64 /srv/shiny-server/BPGA/www/admixture
  ```
- Ensure `unzip`, `gzip`, and system libs for `sf` (GDAL/GEOS/PROJ) are installed.
- Adjust upload limit in `app.R` if needed:
  ```r
  options(shiny.maxRequestSize = 600 * 1024^2) # bytes
  ```

### RStudio / Local
- Same requirements as above. Run with `shiny::runApp(".")`.

### Hosted demo

You can try an interactive instance of **BPGA** here:

**http://15.188.54.171:3838/bpga_app/**

This server is meant for functionality testing and demos. Performance and availability may vary depending on server load. The demo uses the app’s default upload cap (configured in the code) and includes example datasets to explore the workflow quickly.

---

## Troubleshooting
- **“Executable not found / permission denied”**  
  Binaries must exist in `www/` and be executable (`chmod +x`).

- **“unzip: not found” or **“.gz” not unpacking**  
  Install `unzip` and `gzip`/`gunzip` on your system.

- **`sf` installation fails**  
  Install system libs (GDAL/GEOS/PROJ). Check your OS docs (e.g., `apt install libgdal-dev libgeos-dev libproj-dev`).

- **Map export doesn’t save**  
  Install `webshot2` and ensure a headless browser is available.

- **No plots in the zip**  
  Make sure you actually ran the analyses; figures are created on demand under `plots/`.

- **Upload too big**  
  Increase `options(shiny.maxRequestSize=...)` in `app.R`.

---

## Acknowledgments
- **Reference panel**: the app integrates worldwide populations (1000 Genomes Project & HGDP-derived).  
- Tools used: **PLINK**, **GCTA**, **ADMIXTURE**, **R/Shiny**, **Plotly**, **Leaflet**.

---

## References 
- Alexander, D. H., Novembre, J. & Lange, K. Fast model-based estimation of ancestry in unrelated individuals. Genome Res 19, 1655–1664 (2009). 
- Auton, A. et al. A global reference for human genetic variation. Nature 526, 68–74 (2015). 
- Bergström, A. et al. Insights into human genetic variation and population history from 929 diverse genomes. Science 367, (2020). 
- Purcell, S. et al. PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses. The American Journal of Human Genetics 81, 559–575 (2007).
- Yang, J., Lee, S. H., Goddard, M. E. & Visscher, P. M. GCTA: A Tool for Genome-wide Complex Trait Analysis. Am J Hum Genetics 88, 76–82 (2011).

---

**Note**: Third‑party tools **plink**, **gcta** and **admixture** are not provided. Download them from provided links and save these app dependencies at www folder. 
- **Download pages**: 
- **admixture** - http://dalexander.github.io/admixture/download.html
- **plink19** - https://www.cog-genomics.org/plink2/
- **gcta** - https://yanglab.westlake.edu.cn/software/gcta/#Download

---

## License
This project is licensed under the **Apache License, Version 2.0** — see the [`LICENSE`](./LICENSE) file for details.  
A [`NOTICE`](./NOTICE) file is included for attribution. If you distribute modified versions, keep relevant notices as required by Section 4(d) of the license.

> Third‑party tools and datasets (e.g., PLINK, GCTA, ADMIXTURE, 1000G/HGDP) are subject to their own licenses/terms. This repository licenses **the app code only**.

