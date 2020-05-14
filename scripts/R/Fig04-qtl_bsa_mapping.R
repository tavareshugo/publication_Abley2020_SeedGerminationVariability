################################
# Analysis of Germination data #
# Abley et al 2020             #
################################

library(tidyverse)
library(patchwork)

# Change ggplot2 default aesthetics
theme_set(theme_classic() +
            theme(panel.grid = element_blank(),
                  text = element_text(size = 12)))


#### read data ####
# QTL scans were pre-computed with ./scripts/R/data_processing/qtl_run_scans.R
# BSA scans were pre-computed with ./scripts/R/data_processing/bsa_run_scans.R

# results without outlier MAGICs
qtl_outliers <- read_csv("./data/QTLmapping/qtl_scan_no_outlier.csv",
                         col_types = cols())

# Calculate adjusted p-values (false discovery rate, fdr)
# highlight those that have a permuted p-value < 0.05
qtl_outliers <- qtl_outliers %>%
  drop_na(p, bp) %>%
  group_by(trait, transform, marker_cov) %>%
  mutate(fdr = p.adjust(p, method = "fdr"),
         perm_sig = ifelse(p_perm < 0.05, -log10(fdr), NA)) %>%
  ungroup()

# Results for all lines
qtl_all <- read_csv("./data/QTLmapping/qtl_scan_all.csv",
                    col_types = cols())

# Calculate adjusted p-values (false discovery rate, fdr)
# highlight those that have a permuted p-value < 0.05
qtl_all <- qtl_all %>%
  drop_na(p, bp) %>%
  group_by(trait, transform, marker_cov) %>%
  mutate(fdr = p.adjust(p, method = "fdr"),
         perm_sig = ifelse(p_perm < 0.05, -log10(fdr), NA)) %>%
  ungroup()

# Results from bulk-segregant analysis
bsa <- read_csv("./data/BSAmapping/QTLseqr_results.csv",
                col_types = cols(.default = col_double(),
                                 comparison = col_character(),
                                 snp = col_character())) %>%
  rename(chromosome = CHROM)

# location of DOG1 gene (AT5G45830)
dog1_gene <- tibble(chromosome = 5, start = 18589418, end = 18591506)


#### make plot ####

# panel A
p1 <- qtl_all %>%
  filter(transform == "rank" & trait == "cv") %>%
  mutate(cov_marker = ifelse(marker_cov == "none" & marker == "MN3_15679400", -log10(p), NA)) %>%
  ggplot(aes(bp/1e6, -log10(p))) +
  geom_line(aes(colour = marker_cov)) +
  geom_hline(yintercept = 3.5, linetype = "dashed", colour = "grey48") +
  geom_point(aes(y = cov_marker), colour = "orange2") +
  geom_vline(data = dog1_gene, aes(xintercept = start/1e6), linetype = 2) +
  facet_grid(. ~ chromosome, scale = "free_x", space = "free_x") +
  scale_colour_manual(values = c("orange2", "black")) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  theme(legend.position = "none", panel.border = element_rect(fill = NA)) +
  labs(x = "Mb", tag = "A)")

# panel B
p2 <- qtl_outliers %>%
  filter(transform == "rank" & trait == "cv") %>%
  mutate(cov_marker = ifelse(marker_cov == "none" & marker == "MN3_15679400", -log10(p), NA)) %>%
  ggplot(aes(bp/1e6, -log10(p))) +
  geom_line(aes(colour = marker_cov)) +
  geom_hline(yintercept = 3.5, linetype = "dashed", colour = "grey48") +
  geom_point(aes(y = cov_marker), colour = "orange2") +
  geom_vline(data = dog1_gene, aes(xintercept = start/1e6), linetype = 2) +
  facet_grid(. ~ chromosome, scale = "free_x", space = "free_x") +
  scale_colour_manual(values = c("orange2", "black")) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  theme(legend.position = "none", panel.border = element_rect(fill = NA)) +
  labs(x = "Mb", tag = "B)")

# Panel C - sample down SNPs for the graph
set.seed(99383)
p3 <- bsa %>%
  sample_n(10e3) %>%
  ggplot(aes(x = POS/1e6)) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_line(aes(y = CI_95, colour = comparison), linetype = 2) +
  geom_line(aes(y = -CI_95, colour = comparison), linetype = 2) +
  # geom_line(aes(y = CI_99, colour = comparison)) +
  # geom_line(aes(y = -CI_99, colour = comparison)) +
  geom_line(aes(y = tricubeDeltaSNP, colour = comparison), size = 1) +
  geom_vline(data = dog1_gene, aes(xintercept = start/1e6), linetype = 2) +
  facet_grid(~ chromosome, scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  labs(x = "Mb", y = "No-0 allele frequency difference",
       colour = "", tag = "C)") +
  scale_colour_brewer(palette = "Dark2") +
  theme(legend.position = c(0, 0), legend.justification = c(0, 0),
        legend.background = element_blank())


pdf("./figures/Fig04.pdf", width = 7.5, height = 7.5)
p1 + p2 + p3 + plot_layout(ncol = 1, heights = c(1, 1, 2.5))
dev.off()



