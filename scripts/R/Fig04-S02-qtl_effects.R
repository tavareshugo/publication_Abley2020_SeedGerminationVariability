################################
# Analysis of Germination data #
# Abley et al 2020             #
################################

library(tidyverse)
library(patchwork)
library(MagicHelpR)

# Change ggplot2 default aesthetics
theme_set(theme_classic() +
            theme(panel.grid = element_blank(),
                  text = element_text(size = 12)))


#### read data ####

# QTL scans
qtl <- read_csv("./data/QTLmapping/qtl_scan_no_outlier.csv")

# Read the tsv file with phenotypes
phen <- read_tsv("./data/QTLmapping/germ_summaryPerLineForQTLMapping.tsv",
                 col_types = cols()) %>%
  # correct magic line names with "MAGIC." prefix
  mutate(magic_id = gsub("M", "MAGIC.", Magic_line)) %>%
  # select columns with relevant traits
  select(magic_id, mean_median:mean_percent) %>%
  # simplify column names
  rename_all(list(~ str_replace(., "mean_", ""))) %>%
  rename(relrange = rel_range) %>%
  # Remove MAGIC.345 (which has * in its name)
  filter(!str_detect(magic_id, "\\*")) %>%
  # keep only four traits
  select(magic_id, mean, mode, cv, percent)

# Remove 8 most variable lines
phen_filter <- phen %>% filter(min_rank(desc(cv)) > 8)

# Object with MAGIC genotypes
magic_outliers <- read_rds("./data/QTLmapping/magic_gen_object.rds")

# Add phenotypes to the MagicGen object
magic_outliers <- addPhenotypes(magic_outliers, phenotypes = phen_filter, id = "magic_id")


#### Estimate QTL effects ####

# Get markers of interest
qtl %>%
  filter(trait == "cv" & transform == "rank") %>%
  group_by(marker_cov) %>%
  arrange(p) %>%
  slice(1) %>%
  pull(marker)

# The two markers of focus (from the scan for CV)
markers <- c(chr3 = "MN3_15679400", chr5 = "MN5_19823529")

# Calculate effects for both markers and for all 4 traits
set.seed(20190827) # for reproducibility
eff <- markers %>%
  map(function(x){

    # Estimate effect for CV
    eff <- estimateFounderEffect(magic_outliers, phenotype = "cv", marker = x,
                                 n_samples = 1e3) %>%
      rename_at(vars(-accession, -marker, -allele), list(~ paste0(., "_cv")))

    # Estimate effect for the other traits
    for(i in c("mode", "percent", "mean")){
      eff <- estimateFounderEffect(magic_outliers, phenotype = i, marker = x, 1e3) %>%
        rename_at(vars(-accession, -marker, -allele), list(~ paste0(., "_", i))) %>%
        full_join(eff, by = c("accession", "marker", "allele"))
    }
    return(eff)

  }) %>%
  bind_rows() %>%
  mutate(qtl = ifelse(marker == "MN3_15679400", "Chr3", "Chr5"))

# Add founder averages
eff <- read_tsv("./data/MAGICandAccessions/germ_summaryAccessions.tsv", col_types = cols()) %>%
  rename(accession = Magic_line) %>%
  # retain only MAGIC founder accessions
  filter(accession %in% c("Bur-0", "Can-0", "Col-0", "Ct-1",  "Edi-1", "Hi-0",  "Kn-0",
                           "Ler-0", "Mt-0",  "No-0", "Oy-0", "Po-0", "Rsch-4", "Sf-2",
                           "Tsu-0", "Wil-2", "Ws-0",  "Wu-0",  "Zu-0")) %>%
  # select traits of interest and tidy a bit
  select(accession, mean_mean, mean_mode, mean_percent, mean_cv) %>%
  rename_all(list(~ str_replace(., "mean_", "founder_"))) %>%
  mutate(accession = str_remove(accession, "-.$")) %>%
  # scale trait values
  mutate_at(vars(-accession), list(~ as.vector(scale(.)))) %>%
  # join with effect size table
  full_join(eff, by = "accession")


# Calculate correlation between effects
cor_effects <- eff %>%
  group_by(qtl) %>%
  summarise(founder_vs_magic = corLabel(founder_cv, effect_mean_cv, TRUE),
            percent_vs_cv = corLabel(effect_mean_percent, effect_mean_cv, TRUE),
            mode_vs_cv = corLabel(effect_mean_mode, effect_mean_cv, TRUE))


#### Make figure ####

p1 <- eff %>%
  ggplot(aes(effect_mean_cv, founder_cv)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey", linetype = "dashed") +
  geom_errorbarh(aes(xmin = effect_lo_cv, xmax = effect_up_cv), height = 0) +
  #geom_label(aes(label = accession), size = 2) +
  geom_point() +
  geom_text(data = cor_effects, x = 1.1, y = -2.3,
            aes(label = founder_vs_magic),
            hjust = 1, vjust = 0, size = 2.7) +
  ggrepel::geom_text_repel(aes(label = accession), size = 2, alpha = 0.5, direction = "y") +
  facet_grid( ~ qtl) +
  labs(y = "CV (founders)", x = "CV (QTL effect)", tag = "A)") +
  theme(panel.border = element_rect(fill = NA))

p2 <- eff %>%
  ggplot(aes(effect_mean_cv, effect_mean_percent)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey", linetype = "dashed") +
  geom_errorbarh(aes(xmin = effect_lo_cv, xmax = effect_up_cv), alpha = 0) +
  geom_point() +
  geom_text(data = cor_effects, x = -1, y = -0.6,
            aes(label = percent_vs_cv),
            hjust = 0, vjust = 0, size = 2.7) +
  ggrepel::geom_text_repel(aes(label = accession), size = 2, alpha = 0.5) +
  facet_grid( ~ qtl) +
  theme(panel.border = element_rect(fill = NA)) +
  labs(y = "Percent (QTL effect)", x = "CV (QTL effect)", tag = "B)") +
  scale_y_continuous(breaks = seq(-1, 1, 0.5))

p3 <- eff %>%
  ggplot(aes(effect_mean_cv, effect_mean_mode)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey", linetype = "dashed") +
  geom_errorbarh(aes(xmin = effect_lo_cv, xmax = effect_up_cv), alpha = 0) +
  geom_point() +
  geom_text(data = cor_effects, x = -1, y = 1,
            aes(label = mode_vs_cv),
            hjust = 0, vjust = 1, size = 2.7) +
  ggrepel::geom_text_repel(aes(label = accession), size = 2, alpha = 0.5) +
  facet_grid( ~ qtl) +
  theme(panel.border = element_rect(fill = NA)) +
  labs(y = "Mode (QTL effect)", x = "CV (QTL effect)", tag = "C)")

# pdf("./figures/FigS06.pdf", width = 7.5, height = 8)
p1 + p2 + p3 + plot_layout(ncol = 1)
# dev.off()
