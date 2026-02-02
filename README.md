# VAR Analysis – Germany (Quarterly Data)
This repository contains an R script for estimating a Vector Autoregression (VAR)
model using German quarterly macroeconomic and fiscal data.

## Methodology
- Structural VAR with Cholesky identification
- Quarterly frequency (ts object, frequency = 4)
- Lag length selected via AIC, BIC, HQIC, and FPE
- Constant term included

## Variables (Cholesky order)
1. Government expenditure growth
2. Tax revenue growth
3. GDP growth
4. Inflation
5. Short-term interest rate (Euribor 3M)

## Requirements
R packages:
- readr
- lubridate
- zoo
- vars
- ggplot2
- dplyr

## Usage
1. Adjust the file path to the quarterly CSV data.
2. Run `VAR_Rcode.R` to estimate the model and generate outputs.

## Notes
The script is fully reproducible and designed for academic use (e.g. bachelor thesis).
