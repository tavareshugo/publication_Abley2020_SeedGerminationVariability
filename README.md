# Code for Abley & Formosa et al 2020

This repository holds the data analysis code for [Abley & Formosa et al 2020](). 
For the code used in simulations of the gene regulatory network please see the accompanying [GitLab repository](https://gitlab.com/slcu/teamJL/abaga). 

All of the code (and data) is deposited as a static snapshot in the Cambridge Apollo repository (_doi link to be added once published_), corresponding to release (_release to be made once published_) of the GitHub repository.

The code in the GitHub repository is non-static and may be updated if we find any bugs (see "NEWS" section at the bottom).

## Setup

1. Clone the repository (or [download as a zip file]() and extract the files) - you should then have a project directory called _"publication_Abley2020_SeedGerminationVariability"_.
2. Download and extract the zip files with the [data]() from the Appolo repository. Place them inside the project directory. 

You should then have this directory structure: 

```
publication_Abley2020_SeedGerminationVariability
          ├── data
          └── scripts
```

## Running the analysis

Most of the data analysis was done using `R`. 

To recreate the analysis and figures for the paper, open the `2020-Abley_Germination.Rproj` RStudio project file (which will ensure you have the correct working directory set - R scripts use relative paths from the project directory).

Some details about the files on `scripts/R`:

- All scripts named with prefix `Fig` recreate each figure in the paper. Some of the graphs were slightly edited in inkscape for the published version.
- The scripts in `data_processing` do different data processing steps and produce files already provided in the `data` folder in our repository. However, they can be used to recreate those files. 
  - `bsa_run_scans.R` - run the bulk segregant mapping analysis using [QTLSeqr](https://github.com/bmansfeld/QTLseqr).
  - `qtl_run_scans.R` - fit the QTL models on the MAGIC lines using the [MagicHelpR](https://github.com/tavareshugo/MagicHelpR) package. 
  - `GerminationData_MakingLongTables.R` - takes the raw phenotype data and reshapes it into a suitable format for figures. 
- The scripts in `functions` contain some custom functions used in some of the analysis (they are loaded from other scripts).

----

Whole-genome sequence data for bulk-segregant mapping of the _No x Col-0_ F2 was analysed using stand-alone command line tools. 

The code for this analysis is in `scripts/shell`. The scripts are numbered in the order they should be run. We appologise that some of the paths are hard-coded, but the analysis should be fully documented in the code. 

The main output of these scripts is the file `ColxNo_all_naive.tsv` (SNP calls using [freebayes](https://github.com/ekg/freebayes)' naive variant calling method), which we provide in the data repository. 


## Further details

The following data files are provided in our repository:

```
Data/
├── BSAmapping
│   ├── ColxNo_all_naive.tsv
│   ├── QTLseqr_input.tsv
│   └── QTLseqr_results.csv
├── ColxNo_F2germ.csv
├── ColxNo_F3germ.csv
├── cyp707a1_ga3oxMutantsGerm.csv
├── GAandABAdoseResponse.csv
├── germ_summaryForExpComparison.csv
├── germ_summaryForExpComparisonCV.csv
├── HeatShockExperimentGerm.csv
├── MAGICParents.csv
├── MAGICParentsRawData.csv
├── MAGICs.csv
├── MAGICsRawData.csv
├── QTLmapping
│   ├── germ_summaryPerLineForQTLMapping.tsv
│   ├── magic_gen_object.rds
│   ├── qtl_scan_all.csv
│   └── qtl_scan_no_outlier.csv
├── SingleSiliquesGermRawData.csv
├── Soil_vs_plates.csv
├── soil_vs_platesRawData.csv
├── SpanishAccessions.csv
├── SpanishAccessionsRawData.csv
├── statsVsSurvivalCombined234.csv
└── SummaryCVSowingTimes.csv
```

#### Figures 1 & 3

Data files `MAGICsRawData.csv`, `MAGICParentsRawDat.csv` and `SpanishAccessionsRawData.csv` contain germination time data for the MAGIC lines, the MAGIC parental accessions and the Spanish accessions respectively. `soil_vs_platesRawData.csv` contains germination time data for MAGIC lines and accessions used in the soil vs plate comparison in Fig.1 Figure Supplement 1G.

These files were read into the `GerminationData_MakingLongTables.R` script, which converts the daily germination counts into a long format for further analysis. This script outputs the following files: `MAGICs.csv`, `MAGICParents.csv`, `SpanishAccessions.csv`, `Soil_vs_plates.csv`.

These files were read into `Fig01_03_MAGICdataAnalysis.R` which was used to summarise the data and perform the analyses used to generate Fig. 1 and its supplements and Fig. 3 and its supplements. The `Fig01_03_MAGICdataAnalysis.R` script outputs the file `germ_summaryPerLineForQTLMapping.tsv` (file found in `data/QTLmapping`) which contains summary statistics for each MAGIC line used for QTL mapping. This script also uses the files `germ_summaryForExpComparison.csv`, `germ_summaryForExpComparisonCV.csv` and `SummaryCVSowingTimes.csv` for intermediate steps in the analysis used to generate Fig. 1C (see comments in script `Fig01_03_MAGICdataAnalysis.R` for full details).

#### Figure 2

`SingleSiliquesGermRawData.csv` contains germination time data for single and half siliques and was analysed with the script `Fig02_SingleandHalfSiliques.R` to generate the plots for Fig.2 and Fig.2, Figure Supplement 1.

#### Figure 4

Figure 4 was generated using files `QTLmapping/qtl_scan_no_outlier.csv`, `QTLmapping/qtl_scan_all.csv` and `BSAmapping/QTLseqr_results.csv` using `Fig04-qtl_bsa_mapping.R`.

Figure 4, Figure Supplement 3 was generated by using data file `ColxNo_F2germ.csv` with R script `Fig04_Sup3_Col x No F2germ.R`.

Figure 4, Figure Supplement 4 was generated by using data file `ColxNo_F3germ.csv` with R script `Fig04_Sup4_ColxNoF3_germination.R`


#### Figure 5 & 6

See the accompanying [simulations repository](https://gitlab.com/slcu/teamJL/abaga)


#### Figure 7

Figure 7 and its supplement were generated by using data file `GAandABAdoseResponse.csv` with script `Fig07_ABAandGA_doseResponse.R`

#### Figure 8

Figure 8 was generated by using data file `cyp707a1_ga3oxMutantsGerm.csv` with the script `Fig08_GA_ABAmutants_germ.R`

#### Figure 9

Figure 9 and its supplement was generated by using data file `HeatShockExperimentGerm.csv` with script `Fig09_HSexperimentSurvival.R`. The file `statsVsSurvivalCombined234.csv` is used as an intermediate step in this analysis (see R script for details).
