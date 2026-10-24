# Storage Capacity

A Python-based computational screening tool for estimating the theoretical potassium-ion storage capacity of candidate materials from CIF crystal structures. The script evaluates framework composition, estimates potassium storage limits, searches for possible K-ion insertion sites, and exports modified structures for further investigation.

## Overview

The workflow processes CIF files and performs the following steps:

1. **Structure loading:** Reads CIF crystal structures using `pymatgen`, with a fallback for certain malformed numeric fields.
2. **Host framework preparation:** Removes existing potassium species from the structure before evaluating the host framework.
3. **Redox estimation:** Estimates the number of potassium ions that can be accommodated based on framework composition and predefined element-specific assumptions.
4. **Geometric screening:** Estimates insertion limits using BVSE energy-grid analysis through `bvlain`, with a geometric fallback when BVSE analysis fails.
5. **Insertion-site search:** Searches for candidate potassium positions while checking minimum distances to framework atoms and previously accepted potassium ions.
6. **Capacity estimation:** Calculates theoretical capacity using the estimated potassium-ion count and host molar mass.
7. **Output generation:** Saves capacity results to a CSV file and writes K-inserted CIF structures when insertion sites are found.

## Repository Structure

```text
KMOF-Storage-Capacity/
├── main.py
├── requirements.txt
├── README.md
├── cif/
│   └── *.cif
├── output_cif/
│   └── *.cif
└── results.csv
```

- `main.py` — Main storage-capacity screening script.
- `requirements.txt` — Python dependencies.
- `cif/` — Input CIF structures.
- `output_cif/` — Generated structures containing inserted potassium ions.
- `results.csv` — Calculated capacity results.

The script automatically creates the `cif/` and `output_cif/` directories if they do not exist.

## Requirements

- Python 3.10 recommended
- NumPy
- pymatgen
- bvlain for BVSE-based screening

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/grp27-kmof/Storage-Capacity.git
cd Storage-Capacity
```

Replace the repository URL if you use a different GitHub repository name.

### 2. Create a Conda environment

```bash
conda create -n storage-capacity python=3.10 -y
conda activate storage-capacity
```

### 3. Install dependencies

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

If `bvlain` is unavailable through your configured package sources, install it using the appropriate instructions for the version used by your project. The script can fall back to geometric screening if the BVSE import or calculation fails.

## Input Data

Place the CIF files to be screened in the `cif/` directory.

Example:

```text
cif/
├── mof_001.cif
├── mof_002.cif
└── mof_003.cif
```

The current script processes `.cif` files in this directory sequentially. It does not recursively search subdirectories.

## Configuration

The main parameters are defined near the beginning of `main.py`:

| Parameter | Default | Purpose |
|---|---:|---|
| `INPUT_DIR` | `cif` | Input CIF directory |
| `OUTPUT_CIF_DIR` | `output_cif` | Destination for K-inserted CIF files |
| `OUTPUT_CSV` | `results.csv` | Capacity results file |
| `FARADAY` | `96485.33` | Faraday constant in C/mol |
| `MIN_K_HOST_DIST` | `2.00 Å` | Minimum allowed K-to-host distance |
| `MIN_K_K_DIST` | `3.50 Å` | Minimum allowed distance between inserted K ions |
| `GRAPHITE_K` | `279.0 mAh/g` | Graphite capacity benchmark |
| `K_VOL_EST` | `11.0 Å³` | Estimated volume occupied per K ion |

These values can be modified in `main.py` if the screening methodology or assumptions change.

## Storage Capacity Calculation

The script calculates theoretical specific capacity using:

$$
C_{\mathrm{th}} = \frac{nF}{3.6M}
$$

Where:

- **Theoretical capacity (C_th):** The theoretical capacity in mAh/g.
- **Number of electrons (n):** The estimated number of transferable electrons, corresponding to the estimated potassium-ion count used by the script.
- **Faraday constant (F):** Approximately 96,485.33 C/mol.
- **Host molar mass (M):** The host framework molar mass in g/mol.

The estimated capacity is compared against the graphite benchmark of 279 mAh/g to calculate the percentage improvement:

$$
\mathrm{Improvement} =
\left(\frac{C_{\mathrm{th}}}{279} - 1\right)\times 100
$$

A positive value indicates a higher estimated theoretical capacity than the benchmark; a negative value indicates a lower one.

**Important:** These are screening estimates based on the assumptions and heuristics implemented in the script. They are not experimental capacities or definitive proof of electrochemical stability.

## Running the Script

After placing the CIF files in `cif/`, run:

```bash
python main.py
```

The script processes each CIF, reports progress and warnings in the terminal, and appends successful calculations to `results.csv`. At the start of a non-empty run, an existing `results.csv` is removed so that the new run starts with fresh results.

## Output Files

### 1. `results.csv`

The CSV contains the following columns:

| Column | Description |
|---|---|
| `cif_file` | Input CIF filename |
| `Formula` | Reduced chemical formula of the host framework |
| `Molar_Mass_Host` | Host framework molar mass in g/mol |
| `Max_K_Ions_n` | Estimated potassium-ion count used in the capacity calculation |
| `Theoretical_Capacity_mAh_g` | Estimated theoretical capacity in mAh/g |
| `Improvement_vs_KGr_perc` | Percentage improvement relative to the 279 mAh/g graphite benchmark |

Only successfully completed calculations are written to the CSV. Files that raise processing exceptions are reported in the terminal.

### 2. `output_cif/`

When the script finds acceptable insertion sites, it writes a modified CIF containing the host framework and the inserted potassium ions. The output filename is based on the corresponding input filename.

A capacity result may still be generated when no insertion sites are found, so check the terminal warnings and generated CIF files when interpreting results.