################################
# Analysis of Germination data #
# Abley et al 2020             #
################################

library(tidyverse)
library(lme4)

# MAGIC data
germ_magic <- read_csv('./data/MAGICandAccessions/MAGICgerm_tidy.csv',
                       col_types = cols())

# Calculate germination statistics for each plant
germ_summary <- germ_magic %>%
  # restrict analysis to 20-45 range of days after rippening
  filter(MeanDAR > 20 & MeanDAR < 45) %>%
  # calculate stats for each plant
  group_by(Magic_line, Mother.plant.sowing.number, MeanDAR, Plant, id) %>%
  summarise(
    median = median(rep(day, number), na.rm=TRUE),
    mean = mean(rep(day, number), na.rm=TRUE),
    mode = mean(day[which(number == max(number))]),
    var =  var(rep(day, number), na.rm=TRUE),
    sd = sqrt(var),
    cv = sd/mean,
    ls = mean(abs(rep(day, number) - mean)/mean, na.rm = TRUE),
    range = diff(range(rep(day, number), na.rm=TRUE)) + 1,
    rel_range =range /mean,
    iqr = quantile(rep(day, number), 0.75, na.rm = TRUE)- quantile(rep(day, number), 0.25, na.rm = TRUE),
    pct_germ = (Total.germinated[1] /Total[1])*100
  ) %>%
  ungroup()

# transform some of the data
germ_summary <- germ_summary %>%
  mutate_at(vars(mean, mode, pct_germ, cv), list(log = log, rank = GenABEL::rntransform))


#### Fit models ####

# Mixed model to help partition variance into batch (Mother.plant.sowing.number)
# and genotype (Magic_line)

fit_trait_lmm <- function(trait, data){
  if(length(trait) != 1) stop("trait has to be length 1")
  if(!(trait %in% names(data))) stop("trait is not a column of data")

  # formula
  lmm_formula <- as.formula(paste(trait, "~", "(1|Magic_line) + (1|Mother.plant.sowing.number)"))

  # fit model
  lmm_fit <- lmer(lmm_formula, data = data)

  # extract variance and calculate percentage
  lmm_vars <- lmm_fit %>% VarCorr() %>%
    as.data.frame() %>%
    mutate(pct = vcov/sum(vcov)*100,
           trait = trait)

  return(lmm_vars)
}

# Fit model for traits of interest
germ_variances <- map_df(c("cv_rank", "mode_rank", "mean_rank", "pct_germ_log"),
                         fit_trait_lmm, data = germ_summary)

germ_variances %>%
  mutate(pct = paste0(round(pct), "%"),
         grp = case_when(grp == "Magic_line" ~ "MAGIC line",
                         grp == "Mother.plant.sowing.number" ~ "Batch",
                         TRUE ~ grp)) %>%
  select(grp, trait, pct) %>%
  spread(grp, pct) %>%
  gridExtra::grid.table()

# visualisation
germ_variances %>%
  mutate(grp = case_when(grp == "Magic_line" ~ "MAGIC line",
                         grp == "Mother.plant.sowing.number" ~ "Batch",
                         TRUE ~ grp)) %>%
  mutate(grp = factor(grp, levels = c("Residual", "Batch", "MAGIC line"))) %>%
  ggplot(aes(trait, pct, fill = grp)) +
  geom_col() +
  scale_fill_manual(values = c("MAGIC line" = "black",
                               "Residual" = "grey",
                               "Batch" = "grey48"),
                    name = "Group")


# Note:
# transformations improve residuals plot substantially
# for example compare mean with and without rank-transformation
lmer(mean ~ (1|Magic_line) + (1|Mother.plant.sowing.number), data = germ_summary) %>%
  plot()
lmer(mean_rank ~ (1|Magic_line) + (1|Mother.plant.sowing.number), data = germ_summary) %>%
  plot()



