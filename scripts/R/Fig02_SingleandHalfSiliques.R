################################
# Analysis of Germination data #
# Abley et al 2020             #
################################

# Load packages
library(gridExtra)
library(tidyverse)
theme_set(theme_bw())

# source function needed to bootstrap germination times
source("./scripts/R/functions/bootstrap_germ.R")


###
# Read and tidy data ####
###
# Raw data
germ <- read.csv('./data/SingleSiliquesGermRawData.csv',
                 na.strings = "")

# Add unique ID to each silique
germ$id <- 1:nrow(germ)

# Some plants had only one silique sampled. It is given an ID for consistency
germ$Silique <- ifelse(is.na(germ$Silique), "1", germ$Silique)


# Tidy data
## gather days to germination to have long format table
## convert days to numeric
## remove days with no scored germination
## remove siliques for which no seeds germinated

germ_long <- germ %>%
	gather(day, number, X1:X87) %>%
	mutate(day = as.numeric(gsub("X", "", day))) %>%
  filter(!is.na(number)) %>%
  filter(Total.germinated > 0) %>%
  mutate(name = paste(Magic_line, Plant, Silique, Sowing_date, sep = "_")) %>%
  arrange(Sowing_date, Magic_line, Plant, Silique, Top_or_bottom, day)

# Calculate proportion germinating on each day
germ_long <- germ_long %>%
	mutate(perc = number/Total*100,
	       perc_total_ger = number/Total.germinated*100,
	       germ_frac = number/Total.germinated)

# Calculate cumulative germination fraction for each silique
germ_long <- germ_long %>%
  group_by(id) %>%
  arrange(day) %>%
  mutate(cumgerm_frac = cumsum(germ_frac)) %>%
  ungroup()

# Calculate mean and variance for each plant's seeds, this has data for half siliques
germ_summary <- germ_long %>%
  group_by(id, Sowing_date, Magic_line, Plant, Silique, Top_or_bottom) %>%
  summarise(
            median = median(rep(day, number), na.rm=T),
            mean = mean(rep(day, number), na.rm=T),
            mode = mean(day[which(number == max(number))]),
            var =  var(rep(day, number), na.rm=T),
            sd = sqrt(var),
            cv = sd/mean,
            range = diff(range(rep(day, number), na.rm=T)) + 1,
            rel_range = range/mean,
            iqr = IQR(rep(day, number, na.rm = T)),
            pct_germ = (Total.germinated[1] / Total[1])*100,
            totalGerm = Total.germinated[1],
            maxDay=max(day),
            germAfterD20Fraction = ifelse(length(which(day == 20)) == 0, NA, 1-cumgerm_frac[which(day == 20)]),
            totalSeeds = Total[1]) %>%
  ungroup()


### Fig.2 figure supplement 1D
####
## Assess likelihood of distribution difference between bottom and top halves ####
####

# filter to just include the two experiments with top and bottom silique halves
germ_longTopBottom <- germ_long %>%
  filter((Sowing_date == "24-Nov" | Sowing_date == "19-Oct") &
           (Magic_line == "182" | Magic_line=="178"))

germ_summaryTopBottom <- germ_summary %>%
  filter((Sowing_date == "24-Nov" | Sowing_date == "19-Oct") &
           (Magic_line == "182" | Magic_line=="178"))

# We use a bootstrap-based approach to compare top and bottom halves
## In summary:
## 1. We pool data from all siliques to obtain a proxy for seed germination distributions of each genotype in these experiments
## 2. We sample germination times from this distribution using sample sizes equal to each silique assessed
## 3. Step 2 is repeated 1000 times for each silique half to get null distributions of CV differences
## 4. An empirical p-value is estimated by seeing how many of the 1000 simulations are above the observed CV difference

# Calculate fraction germination across entire experiment per genotype
pooled_germ_counts <- germ_longTopBottom %>%
  group_by(Magic_line, day) %>%
  summarise(number_germ = sum(number),
            total_germ = sum(Total.germinated)) %>%
  mutate(prop = number_germ/total_germ) %>%
  ungroup()

# Get sample sizes of each silique's top and bottom halves
silique_sample_sizes <- germ_summaryTopBottom %>%
  select(Sowing_date, Magic_line, Plant, Silique, Top_or_bottom, totalGerm) %>%
  spread(Top_or_bottom, totalGerm)

# Join this table to the pooled table
pooled_germ_counts <- left_join(pooled_germ_counts, silique_sample_sizes,
                                by = "Magic_line")

# Run our sampling bootstrap for each silique
set.seed(123409) # for reproducibility
null_cv <- pooled_germ_counts %>%
  filter(!is.na(bottom) & !is.na(top)) %>%
  group_by(Sowing_date, Magic_line, Plant, Silique) %>%
  do(boot_germ(.$day, .$prop, unique(.$bottom), unique(.$top))) %>%
  ungroup()

# Join up the observed CVs
null_cv <- germ_summaryTopBottom %>%
  select(Sowing_date, Magic_line, Plant, Silique, Top_or_bottom, cv,
         d20 = germAfterD20Fraction) %>%
  gather(trait, value, cv, d20) %>%
  unite("key", trait, Top_or_bottom) %>%
  spread(key, value) %>%
  full_join(null_cv, by = c("Sowing_date", "Magic_line", "Plant", "Silique")) %>%
  mutate(obs_cv_dif = cv_top - cv_bottom,
         obs_d20_dif = d20_top - d20_bottom,
         sim_cv_dif = sim_cv_top - sim_cv_bottom,
         sim_d20_dif = sim_d20_top - sim_d20_bottom)

# Visualisation
null_cv %>%
  filter((Sowing_date == "24-Nov" | Sowing_date == "19-Oct") &
           (Magic_line == "182" | Magic_line=="178")) %>%
  ggplot(aes(interaction(Plant, Silique), sim_cv_dif)) +
  geom_violin() +
  geom_point(aes(y = obs_cv_dif), colour = "brown") +
  facet_wrap(~ interaction(Sowing_date, Magic_line), scales = "free_x")

null_cv %>%
  filter((Sowing_date == "24-Nov" | Sowing_date == "19-Oct") &
           (Magic_line == "182" | Magic_line=="178")) %>%
  ggplot(aes(interaction(Plant, Silique), sim_d20_dif)) +
  geom_violin() +
  geom_point(aes(y = obs_d20_dif), colour = "brown") +
  facet_wrap(~ interaction(Sowing_date, Magic_line), scales = "free_x")

# Convert to "p-values"
silique_tests <- null_cv %>%
  group_by(Sowing_date, Magic_line, Plant, Silique, obs_cv_dif, obs_d20_dif) %>%
  summarise(p_cv = (1 + sum(abs(sim_cv_dif) > abs(obs_cv_dif)))/(n() + 1),
            p_d20 = (1 + sum(abs(sim_d20_dif) > abs(obs_d20_dif)))/(n() + 1)) %>%
  ungroup() %>%
  mutate(padj_cv = p.adjust(p_cv, method = "fdr"),
         padj_d20 = p.adjust(p_d20, method = "fdr")) %>%
  # join with germination summary
  full_join(germ_summaryTopBottom,
            by = c("Sowing_date", "Magic_line", "Plant", "Silique"))



### Fig.2 plotting

##M178 single siliques for Fig.2

#png("./plots/PNGs/SingleSiliques/M178SingleSiliquesJul2018.png", width = 9, height = 9, units = "cm", res = 600)
germ_long %>%
  filter(Magic_line=="178" & (Sowing_date=="14-Jun"| Sowing_date =="22-Sep")) %>%
  ggplot(aes(x = day, y = name, size=ifelse(perc==0, NA, perc))) + geom_point(stat = "identity") +
  xlim(0, 60)+
  theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Days to germination")+
  ylab("Replicate")+
  theme(legend.position="none")+
  ggtitle("M178 single siliques")
#dev.off()

##outputting source data for Fig 2B
# M178SingleSiliquesSourceData<-germ_long %>%  filter(Magic_line=="178" & (Sowing_date=="14-Jun"| Sowing_date =="22-Sep"))
# write.table(M178SingleSiliquesSourceData, './source_data/Figure2_source_data2.tsv',
#             quote = F, row.names = F, sep = "\t")

##M178 half silique plot for Fig 2
# png("./plots/M178SingleSiliques19OctTopBottom.png", width = 14, height = 10, units = "cm", res = 600)
germ_longTopBottom %>%
  filter(Magic_line=="178"& Sowing_date =="19-Oct") %>%
  ggplot(aes(x = day, y = name, size=ifelse(perc==0, NA, perc))) + geom_point(stat = "identity") +
  xlim(0, 60)+
  theme(text = element_text(size=12))+
  facet_grid(~Top_or_bottom)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Days to germination")+
  ylab("Replicate")+
  theme(legend.position="none")+
  ggtitle("M178 half siliques")
# dev.off()

##outputting source data for Fig 2C
# M178halfSiliquesSourceData<-germ_longTopBottom %>%  filter(Magic_line=="178"& Sowing_date =="19-Oct")
# write.table(M178halfSiliquesSourceData, './source_data/Figure2_source_data3.tsv',
#             quote = F, row.names = F, sep = "\t")


###Fig.2 Figure supplement 1

##M182 half silique plot

# png("./plots/M182SingleSiliques19OctTopBottom.png", width = 14, height = 10, units = "cm", res = 600)
germ_long%>%
  filter(Magic_line=="182" & Sowing_date == "19-Oct" ) %>%
  ggplot(aes(x = day, y = name, size=ifelse(perc==0, NA, perc))) + geom_point(stat = "identity") +
  xlim(0, 60)+
  theme(text = element_text(size=12))+
  facet_grid(~Top_or_bottom)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Days to germination")+
  ylab("Replicate")+
  theme(legend.position="none")+
  ggtitle("M182 half siliques")
# dev.off()


##M53 half silique plots

#png("./plots/M53SingleSiliquesTopBottom.png", width = 14, height = 10, units = "cm", res = 600)
germ_long %>%
  filter(Magic_line=="53" & Sowing_date=="19-Oct") %>%
  ggplot(aes(x = day, y = name, size=ifelse(perc==0, NA, perc))) + geom_point(stat = "identity") +
  xlim(0, 60)+
  facet_grid(. ~Top_or_bottom)+
  scale_colour_gradient2(low = "grey", mid = "orange", high = "red3", midpoint = 30) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Days to germination")+
  ylab("Replicate")+
  ggtitle("M53 half siliques")+
  theme(legend.position="none")
#dev.off()


##output source data for M182 and M53 half siliques
# M182_M53HalfSiliquesPlotting<-germ_long %>%  filter ((Magic_line=="182" & Sowing_date == "19-Oct" )| (Magic_line=="53" & Sowing_date=="19-Oct"))
# write.table(M182_M53HalfSiliquesPlotting, './source_data/Figure2_FigureSupplement1_source_data2.tsv',
#             quote = F, row.names = F, sep = "\t")


###M4 single silique

#png("./plots/M4SingleSiliques.png", width = 15, height = 10, units = "cm", res = 600)
germ_long %>%
  filter(Magic_line=="4" & (Sowing_date=="14-Jun" | Sowing_date=="22-Sep")) %>%
  ggplot(aes(x = day, y = name, size=ifelse(perc==0, NA, perc))) + geom_point(stat = "identity") +
  xlim(0, 60)+
  scale_colour_gradient2(low = "grey", mid = "orange", high = "red3", midpoint = 30) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Days to germination")+
  ylab("Replicate")+
  ggtitle("M4 half siliques")+
  theme(legend.position="none")
#dev.off()

##output source data for M4 single siliques
# M4SingleSiliquesPlotting<-germ_long %>%  filter(Magic_line=="4" & (Sowing_date=="14-Jun" | Sowing_date=="22-Sep"))
#
# write.table(M4SingleSiliquesPlotting, './source_data/Figure2_FigureSupplement1_source_data3.tsv',
#             quote = F, row.names = F, sep = "\t")



# Panel D-i - CV
# png("./plots/CVcomparisonSingleSiliques.png", width = 9, height = 7.5, units = "cm", res = 600)
silique_tests  %>%
  ggplot(aes(Top_or_bottom, cv, group = interaction(Plant, Silique))) +
  geom_point() +
  geom_line(aes(colour = padj_cv < 0.05), show.legend = FALSE) +
  facet_grid(Sowing_date ~ paste0("M", Magic_line)) +
  scale_colour_manual(values = c("black", "brown")) +
  theme(text = element_text(size=12))+
  labs(x = "Top or bottom of silique", y = "CV")
# dev.off()

# Panel D-ii - d20
# png("./plots/FractionAfterD20comparisonSingleSiliques.png", width = 9, height = 7.5, units = "cm", res = 600)
silique_tests   %>%
  ggplot(aes(Top_or_bottom, germAfterD20Fraction,
             group = interaction(Plant, Silique))) +
  geom_point() +
  geom_line(aes(colour = padj_d20 < 0.05), show.legend = FALSE) +
  facet_grid(Sowing_date ~ paste0("M", Magic_line)) +
  scale_colour_manual(values = c("black", "brown")) +
  labs(x = "Top or bottom of silique", y = "Fraction after D20")
# dev.off()

##output source data for D)
# siliqueTestsSourceData <- silique_tests %>%
#   select(Sowing_date, Magic_line, Plant, Silique, padj_cv, padj_d20, id, Top_or_bottom, cv, germAfterD20Fraction)
# write.table(siliqueTestsSourceData, './source_data/Figure2_FigureSupplement1_source_data4.tsv',
#             quote = F, row.names = F, sep = "\t")




