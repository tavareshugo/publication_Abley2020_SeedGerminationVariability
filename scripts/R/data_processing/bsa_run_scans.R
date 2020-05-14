#####
# Script to run bulk-segregant-mapping scans
# for Abley et al. 2019
#####

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(QTLseqr)  # devtools::install_github("bmansfeld/QTLseqr")

# custom function to run QTLseqr pipeline
source("./scripts/R/functions/pairwise_seqr.R")


# read data ---------------------------------------------------------------

# table of SNPs
snps <- read_tsv("./data/BSAmapping/ColxNo_all_naive.tsv",
                 col_types = cols(
                   CHROM = col_character(),
                   POS = col_integer(),
                   REF = col_character(),
                   ALT = col_character(),
                   QUAL = col_double(),
                   ColxNo_bulk_E1.AD = col_character(),
                   ColxNo_bulk_E1.DP = col_double(),
                   ColxNo_bulk_E1.GQ = col_integer(),
                   ColxNo_bulk_E1.GT = col_character(),
                   ColxNo_bulk_L1.AD = col_character(),
                   ColxNo_bulk_L1.DP = col_double(),
                   ColxNo_bulk_L1.GQ = col_integer(),
                   ColxNo_bulk_L1.GT = col_character(),
                   ColxNo_bulk_L2.AD = col_character(),
                   ColxNo_bulk_L2.DP = col_double(),
                   ColxNo_bulk_L2.GQ = col_integer(),
                   ColxNo_bulk_L2.GT = col_character()
                 ))

# Remove genotype quality columns (genotype quality was not called)
snps <- snps %>% select(-matches("GQ$"))

# Tidy up the columns
snps <- snps %>%
  # reshape to long format
  gather("key", "value", matches("^Col")) %>%
  # separate name of the bulk and the metric
  separate(key, c("bulk", "metric"), sep = "\\.") %>%
  # spread each metric to its own column
  spread(metric, value, convert = TRUE) %>%
  # retain only reference allele count and tidy name of bulks
  mutate(AD = as.integer(str_remove(AD, ",.*")),
         bulk = str_remove(bulk, "ColxNo_bulk_")) %>%
  # tidy column names
  rename_all(tolower) %>%
  # add SNP id and calculate frequency
  mutate(snp = paste(chrom, pos, sep = "_"),
         alt_freq = 1 - ad/dp)


# filter ------------------------------------------------------------------

snps_filter <- snps %>%
  # remove chloroplast and mitochondrial genomes
  filter(chrom %in% 1:5) %>%
  # retain only quality SNPs
  filter(qual > 20 & dp > 20 & dp < 200) %>%
  # remove SNPs with two alternative alleles
  filter(!str_detect(alt, ",")) %>%
  # retain only those that appear in all 3 bulks
  group_by(chrom, pos) %>%
  filter(n() == 3) %>%
  ungroup()



# QTLseqr analysis --------------------------------------------------------

# Write file in the right format
snps_filter %>%
  select(CHROM = chrom, POS = pos, AD = ad, DP = dp, bulk) %>%
  gather("key", "value", AD, DP) %>%
  mutate(bulk = paste(bulk, key, sep = ".")) %>%
  select(-key) %>%
  spread(bulk, value) %>%
  drop_na() %>%
  write_tsv("./data/BSAmapping/QTLseqr_input.tsv")


# Run pairwise comparisons of pools (using custom function that calls QTLseqr)
bsa_qtl <- pairwise_seqr("./data/BSAmapping/QTLseqr_input.tsv",
                         windowSize = 1e6, replications = 1e4,
                         pool_sizes = c(E1 = 152, L1 = 321, L2 = 213),
                         filter = 0)

# make SNP id
bsa_qtl$snp <- paste(bsa_qtl$CHROM, bsa_qtl$POS, sep = "_")

# write results
bsa_qtl %>%
  write_csv("./data/BSAmapping/QTLseqr_results.csv")
