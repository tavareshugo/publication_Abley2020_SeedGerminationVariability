# Code for Abley & Formosa et al 2020

This repository holds the data analysis code for [Abley & Formosa et al 2020](). For the code used in simulations of the gene regulatory network please see the partner [GitLab repository](https://gitlab.com/slcu/teamJL/abaga). 

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
