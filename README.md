# Replication Package for: Robust forecast evaluation under extreme shocks

**August 04, 2026**

This replication package accompanies Iacone, F. and A. Viselli (2026). *Robust forecast evaluation under extreme shocks*.

## Data availability and provenance statements

### Statement about rights

The author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript.

### Summary of availability

All data are publicly available.

### Details on each data source

Data on forecast errors for US real Gross Domestic Product (GDP) growth nowcasts from the Survey of Professional Forecasters (SPF) were downloaded from the Federal Reserve Bank of Philadelphia on November 28, 2025. The analysis uses the median SPF nowcasts and the corresponding realizations.

Data can be downloaded from https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/error-statistics.

A copy of the data is provided as part of this archive. The data are in the public domain. Datafile: `SPF_Error_Statistics.mat` (MATLAB, proprietary).

### Details on the dataset

The data file contains the real GDP realization series (`Realiz1`) and the median one-step-ahead SPF forecast series (`SPFfor_Step1`). The empirical application uses quarterly observations from 2000:Q1 through 2024:Q4. A guide to the error-statistics dataset and Stark (2010) are available from the data-source page linked above.

## Description of programs/code

Run the following programs/code in MATLAB R2023a (recommended).
The replication files are compatible with macOS, Windows, and Linux, provided that a compatible version of MATLAB is installed.
The code does not include any OS- or Windows-specific paths or functions, so no platform-dependent adjustments should be necessary.

*   `ReplicationCode_MC_Table_1_2.m` generates Tables 1 and 2 in Section 3 (MONTE CARLO EXERCISE).
*   `run_spf_mc.m` is the helper function called by `ReplicationCode_MC_Table_1_2.m` to run each Monte Carlo design.
*   `ReplicationCode_Fig_1_Table_3.m` generates Figure 1 and Table 3 in Section 4 (EMPIRICAL APPLICATION TO US REAL GDP GROWTH) and reports the empirical calibration parameters used in the Monte Carlo exercise.

The Monte Carlo replication code uses 10,000 replications and a fixed random seed. The empirical script requires the MATLAB Econometrics Toolbox for the AR(1) calibration.

## License for Code

The code is licensed under a MIT license. See license.txt for details.

## List of tables and figures

The provided code reproduces all tables and figures in the paper.

| Figure/Table # | Program | Page Number |
| :--- | :--- | :--- |
| Table 1 | `ReplicationCode_MC_Table_1_2.m` | 8 |
| Table 2 | `ReplicationCode_MC_Table_1_2.m` | 9 |
| Table 3 | `ReplicationCode_Fig_1_Table_3.m` | 11 |
| Figure 1 | `ReplicationCode_Fig_1_Table_3.m` | 11 |

## Instructions to replicators

*   Set the working directory to the folder containing all MATLAB scripts and data files.
*   Execute `ReplicationCode_MC_Table_1_2.m` to reproduce Tables 1 and 2. The script calls `run_spf_mc.m` automatically.
*   Execute `ReplicationCode_Fig_1_Table_3.m` to reproduce Figure 1 and Table 3.
*   Tables are displayed in the **MATLAB Command Window**; Figure 1 is displayed in a MATLAB figure window.
*   The expected running time is approximately a few seconds for both the Monte Carlo simulations and `ReplicationCode_Fig_1_Table_3.m`.

## References

*   Coroneo, L. and F. Iacone (2020). Comparing predictive accuracy in small samples using fixed-smoothing asymptotics. *Journal of Applied Econometrics* 35, 391–409.
*   Diebold, F. X. and R. S. Mariano (1995). Comparing predictive accuracy. *Journal of Business & Economic Statistics* 20, 134–144.
*   Iacone, F., L. Rossini and A. Viselli (2026). Comparing predictive ability in the presence of instability over a very short time. *The Econometrics Journal* 29, 143–166.
*   Stark, T. (2010). Realistic evaluation of real-time forecasts in the Survey of Professional Forecasters. *Research Rap*, Special Report 1.
*   Survey of Professional Forecasters error-statistics documentation.

## Authors

*   Fabrizio Iacone
*   Andrea Viselli
