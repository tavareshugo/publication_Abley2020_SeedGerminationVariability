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



#### make plot ####

# panel A
p1 <- qtl_all %>%
  filter(transform == "rank") %>%
  mutate(cov_marker = ifelse(marker_cov == "none" & marker == "MN3_15679400", -log10(p), NA)) %>%
  ggplot(aes(bp/1e6, -log10(p))) +
  geom_line(aes(colour = marker_cov)) +
  geom_hline(yintercept = 3.5, linetype = "dashed", colour = "grey48") +
  geom_point(aes(y = cov_marker), colour = "orange2") +
  facet_grid(trait ~ chromosome, scale = "free_x", space = "free_x") +
  scale_colour_manual(values = c("orange2", "black")) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  theme(legend.position = "none", panel.border = element_rect(fill = NA)) +
  labs(x = "Mb", tag = "A)")

# panel B
p2 <- qtl_outliers %>%
  filter(transform == "rank" ) %>%
  mutate(cov_marker = ifelse(marker_cov == "none" & marker == "MN3_15679400", -log10(p), NA)) %>%
  ggplot(aes(bp/1e6, -log10(p))) +
  geom_line(aes(colour = marker_cov)) +
  geom_hline(yintercept = 3.5, linetype = "dashed", colour = "grey48") +
  geom_point(aes(y = cov_marker), colour = "orange2") +
  facet_grid(trait ~ chromosome, scale = "free_x", space = "free_x") +
  scale_colour_manual(values = c("orange2", "black")) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  theme(legend.position = "none", panel.border = element_rect(fill = NA)) +
  labs(x = "Mb", tag = "B)")



pdf("./figures/Fig04_S01.pdf", width = 7.5, height = 7.5)
p1 + p2 + plot_layout(ncol = 1)
dev.off()




