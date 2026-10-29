# cotopaxi-so2

Code and aggregated data for: **Association between Cotopaxi volcano SO₂ and respiratory hospital discharges, Cotopaxi Province, Ecuador (2022 vs 2024)**.

This repository supports the Data Availability Statement of the manuscript. It contains the analysis code and the aggregated (canton-by-week) data used for the tables and figures. It does **not** contain individual hospital records.

## Contents

| Path | What it is |
|---|---|
| `code/GeoSalud_Table1_Verification.ipynb` | Reproduces the annual respiratory discharges by canton (Table 1, ICD-10 J00–J99, by canton of residence) from the raw INEC microdata. |
| `code/GeoSalud_supplementary_tables_colab.py` | Builds Supplementary Tables S2–S7 (leave-one-out, residual autocorrelation, MAD cleaning counts, cleaned-series sensitivity, wind check, quasi-Poisson models adjusted for epidemiological week and its square) from `data/weekly_panel.csv`. |
| `data/weekly_panel.csv` | Weekly panel: canton, year, ISO week, mean TROPOMI/Sentinel-5P SO₂ column density (mol/m², raw, missing weeks left empty), and respiratory discharges (by discharge date). |
| `data/population_INEC_rev2024_cotopaxi.csv` | Population by canton, age group (0–5, 6–17, 18–59, ≥60, and total) for 2022 and 2024, summed from the INEC *Estimaciones y Proyecciones de Población, revisión 2024* (cantonal tables by single year of age, Table 3.1). Denominators of the hospitalization rates. |
| `data/ERA5_wind_highSO2_weeks.csv` | ERA5-Land 10 m wind at the crater for the five highest-SO₂ weeks of each year (post hoc wind check, Supplementary Table S6). |
| `results/` | Output tables used in the manuscript: Table 2 (raw series), Table 3 and Table 4 (diagnostic and age groups, four core cantons), and the canton-mean SO₂ values of Figure 2. |

## Data sources

* **SO₂:** `COPERNICUS/S5P/OFFL/L3_SO2`, band `SO2_column_number_density`, extracted over FAO GAUL level-2 canton polygons with Google Earth Engine (weekly means).
* **Hospital discharges:** INEC, *Registro Estadístico de Camas y Egresos Hospitalarios* (https://www.ecuadorencifras.gob.ec/camas-y-egresos-hospitalarios/). The raw microdata are public at INEC and are not redistributed here.
* **Wind:** ERA5-Land (Copernicus Climate Change Service).

## How to reproduce

1. Download the INEC discharge files for 2022 (`.sav`) and 2024 (`.csv`) and run `code/GeoSalud_Table1_Verification.ipynb` to obtain the annual totals of Table 1.
2. Put `data/weekly_panel.csv` in the working folder and run `code/GeoSalud_supplementary_tables_colab.py` (Python 3, `pandas`, `numpy`, `scipy`, `statsmodels`, `python-docx`). It writes `Table_S2.csv` … `Table_S7.csv` and `Supplementary_Material.docx`.

## Notes

* The primary analysis uses the raw (unreplaced) SO₂ series; the MAD-cleaned series (median + 2.5 × MAD, 0.0002 mol/m², weeks 40–52 of 2022 exempt) is a secondary sensitivity analysis.
* Excess-risk estimates are from quasi-Poisson models (Poisson standard errors multiplied by the square root of the dispersion ratio), per interquartile-range increase of SO₂.
* Week assignment follows ISO 8601 using the hospital discharge date.

## License

MIT (see `LICENSE`).
