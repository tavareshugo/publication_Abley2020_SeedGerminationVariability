#####
# Script to run QTL scans
# for Abley et al. 2019
#####

#################################################

# Run as: Rscript qtl_run_scans.R
# Adjust working directory and how many cores are available

# Ran on our HPC as:
# sbatch -n 48 --job-name=seed_qtl \
#   -D ~/projects/20161101_seed_qtl
#   -o ./logs/qtl_run_scans/qtl_run_scans.stdout \
#   --wrap "Rscript ./scripts/R/data_processing/qtl_run_scans.R"


# Define how many cores to use in this run (depending where the script is run)
ncores <- 48

# Define how many (if any) permutations
nperm <- 1000

#################################################

# Load libraries
library(methods)
library(MagicHelpR) # devtools::install_github("tavareshugo/MagicHelpR@v0.1")
library(tidyverse)
library(stringr)


#
# Read and tidy data ----
#
# Read MAGIC genotype object
magic_gen <- read_rds("./data/QTLmapping/magic_gen_object.rds")

# Read the tsv file with phenotypes
phen <- read_tsv("./data/QTLmapping/germ_summaryPerLineForQTLMapping.tsv",
                 col_types = cols())

# correct magic line names with "MAGIC." prefix
# select columns with relevant traits
# simplify column names
phen <- phen %>%
  mutate(magic_id = gsub("M", "MAGIC.", Magic_line)) %>%
  rename_all(list(~ str_replace(., "mean_", ""))) %>%
  select(magic_id, mean, mode, cv, percent)

# Remove MAGIC.345 (which has * in its name)
phen <- phen %>%
  filter(!str_detect(magic_id, "\\*"))


# Create subset without outlier lines
phen_filter <- phen %>% filter(min_rank(desc(cv)) > 8)


# Apply log and rank transformation to all variables as well as arcsine transform to percentage
phen <- phen %>%
  mutate_if(is.numeric, funs(log = log, rank = GenABEL::rntransform)) %>%
  mutate(percent_asin = asin(sqrt(percent/100)))

phen_filter <- phen_filter %>%
  mutate_if(is.numeric, funs(log = log, rank = GenABEL::rntransform)) %>%
  mutate(percent_asin = asin(sqrt(percent/100)))



#
# QTL scans excluding outliers ----
#
# Add phenotypes
magic_phen <- addPhenotypes(magic_gen, phen_filter, "magic_id")

# Phenotypes of interest
phenotypes <- names(phen_filter)[-1]

# QTL scan for all phenotypes
qtl_scan <- lapply(phenotypes, function(p, magic_phen){
  scanQtl(magic_phen, p, cores = ncores, perm = nperm)
}, magic_phen)

# Bind into a single table
names(qtl_scan) <- phenotypes
qtl_scan_outliers <- bind_rows(qtl_scan, .id = "trait") %>%
  mutate(marker_cov = "none")

# QTL scan using the Chr3 marker as covariate
qtl_scan <- lapply(phenotypes, function(p, magic_phen){
  scanQtl(magic_phen, p, cores = ncores, perm = nperm, marker_cov = "MN3_15679400")
}, magic_phen)

# Bind into a single table
names(qtl_scan) <- phenotypes
qtl_scan_outliers <- qtl_scan %>%
  bind_rows(.id = "trait") %>%
  mutate(marker_cov = "MN3_15679400") %>%
  bind_rows(qtl_scan_outliers) %>%
  separate(trait, c("trait", "transform"), sep = "_") %>%
  mutate(transform = ifelse(is.na(transform), "raw", transform))

write_csv(qtl_scan_outliers, "./data/QTLmapping/qtl_scan_no_outlier.csv")


#
# QTL scans all lines ----
#
# Add phenotypes
magic_phen <- addPhenotypes(magic_gen, phen, "magic_id")

# Phenotypes of interest
phenotypes <- names(phen_filter)[-1]

# QTL scan for all phenotypes
qtl_scan <- lapply(phenotypes, function(p, magic_phen){
  scanQtl(magic_phen, p, cores = ncores, perm = nperm)
}, magic_phen)

# Bind into a single table
names(qtl_scan) <- phenotypes
qtl_scan_all <- bind_rows(qtl_scan, .id = "trait") %>%
  mutate(marker_cov = "none")

# QTL scan using the Chr3 marker as covariate
qtl_scan <- lapply(phenotypes, function(p, magic_phen){
  scanQtl(magic_phen, p, cores = ncores, perm = nperm, marker_cov = "MN3_15679400")
}, magic_phen)

# Bind into a single table
names(qtl_scan) <- phenotypes
qtl_scan_all <- qtl_scan %>%
  bind_rows(.id = "trait") %>%
  mutate(marker_cov = "MN3_15679400") %>%
  bind_rows(qtl_scan_all) %>%
  separate(trait, c("trait", "transform"), sep = "_") %>%
  mutate(transform = ifelse(is.na(transform), "raw", transform))


write_csv(qtl_scan_all, "./data/QTLmapping/qtl_scan_all.csv")


