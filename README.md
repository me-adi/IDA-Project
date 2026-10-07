# IDA-Project

# Regression Analysis for Establishing a Relation Between Weather Parameters

Group 19 project for **Introduction to Data Analytics (IDA)**. It uses the Daily
Delhi Climate dataset and asks how mean temperature, humidity, wind speed and
mean pressure are related to each other. Each of the six pairs of parameters is
fitted with simple linear regression and simple non-linear regression, and every
model is compared by R².

| Name | Roll No. | Email |
|---|---|---|
| Aditya Palapati (group leader) | S20230010176 | aditya.p23@iiits.in |
| Meghana Kuruva | S20230010134 | meghana.k23@iiits.in |
| Rohin Sai Bogadi | S20230010047 | rohinsai.b23@iiits.in |
| Srimeenakshi K S | S20262010002 | srimeenakshi.ks@iiits.in |

## Contents

```
data/
  DailyDelhiClimateTrain.csv   2013-01-01 .. 2017-01-01, 1462 days (raw, unchanged)
  DailyDelhiClimateTest.csv    2017-01-01 .. 2017-04-24, 114 days (raw, unchanged)
  processed/clean.csv          cleaned data written by the script
R/
  Project_code.R               the complete analysis (one script)
outputs/
  tables/                      results as .csv and as LaTeX .tex tables
  figures/                     plots as .png (300 dpi) and .pdf
report/
  main.tex                     LaTeX report
  main.pdf                     compiled report
```

## Running the analysis

You need R 4.x. No add-on packages are used.

From the project folder (the one that contains `data/`), run:

```
Rscript R/Project_code.R
```

The run takes about 1–2 minutes. It rebuilds `data/processed/clean.csv` and
everything in `outputs/tables/` and `outputs/figures/`, creating those folders
if they are missing, and prints a step-by-step log. The script uses relative
paths only and stops with a message if it is not started from the project
folder.

The script is written in **base R only**. It uses no packages and none of R's
ready-made statistics or modelling functions (`lm`, `nls`, `cor`, `var`, `sd`,
`quantile`, `median`, `approx`, `solve`, …). Every estimator is implemented by
hand, with its formula written in the comments. The only built-in statistical
function it calls is `pt()`, the Student-t distribution function, for p-values.

The script has six sections:

1. Setup and helper functions
2. Pre-processing (2a)
3. Relation analysis (2b)
4. Simple linear regression (2c-i)
5. Simple non-linear regression (2c-ii): quadratic, cubic, exponential, logarithmic, power
6. R² comparison and conclusions (2d)

## Building the report

The report reads its tables, figures and quoted numbers straight from
`outputs/`, so every number in it matches the code. Run the R script first,
then compile from the `report/` folder:

```
cd report
pdflatex main.tex
pdflatex main.tex
```

Run `pdflatex` twice so the table of contents and cross-references are filled in.
The report only needs standard LaTeX packages: graphicx, booktabs, amsmath,
hyperref, geometry, float, caption and listings. Keep `report/`, `outputs/` and
`R/` side by side, because `main.tex` refers to `../outputs/` and `../R/`.

## Main results

These are the full-data R² values, as reported in `outputs/tables/r2_comparison.csv`.

| Relation | Linear R² | Best curve | Best R² |
|---|---|---|---|
| meantemp ~ meanpressure | 0.776 | cubic | 0.819 |
| meantemp ~ humidity | 0.329 | cubic | 0.331 |
| humidity ~ meanpressure | 0.115 | quadratic | 0.187 |
| humidity ~ wind_speed | 0.160 | exponential | 0.165 |
| meantemp ~ wind_speed | 0.098 | exponential | 0.094 |
| wind_speed ~ meanpressure | 0.090 | cubic | 0.091 |

- **Temperature and pressure** have by far the strongest relation (r = −0.881).
  About 76% of that R² comes from the annual cycle the two share.
- **Temperature and humidity** have a moderate relation (r = −0.574). It is a
  genuine day-to-day link, and it survives removing the seasons.
- **Non-linear curves** give a meaningful gain (ΔR² > 0.05) for humidity and
  pressure only. For the other relations a straight line is adequate.
- **Pre-processing** fixed 10 impossible pressure values by interpolation and
  capped 33 extreme wind speeds at 17.61 km/h.

The report explains these results in detail, including the checks for
overfitting, outliers, autocorrelation and seasonality.
