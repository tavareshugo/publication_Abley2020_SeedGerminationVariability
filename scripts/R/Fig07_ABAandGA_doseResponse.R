################################
# Analysis of Germination data #
# Abley et al 2020             #
################################

library(dplyr) # install.packages("dplyr")
library(tidyr) # install.packages("tidyr")
library(ggplot2)
theme_set(theme_bw())

germ <- read.csv('./data/GAandABAdoseResponse.csv')
germ$id <- 1:nrow(germ)

# Make data long format and tidy it
germ_long <- germ %>%
	gather(day, number, X1:X20) %>%
	arrange(Magic_line, Hormone, conc_uM, Experiment, High.or.low.variability) %>%
	mutate(day = as.numeric(gsub("X", "", day))) %>%
  filter(!is.na(number))


# Calculate proportion on each day
germ_long <- germ_long %>%
	mutate(perc = number/Total*100) %>%
  mutate(final_percent = Total.germinated_day20/ Total *100)


# Calculate mean and variance for each plant's seeds
germ_summary <- germ_long %>%
	group_by(Magic_line, Hormone, conc_uM, Experiment, High.or.low.variability) %>%
	summarise(
	          median = median(rep(day, number), na.rm=T),
						mean = mean(rep(day, number), na.rm=T),
						mode = mean(day[which(number == max(number))]),
						var =  var(rep(day, number), na.rm=T),
						final_percent= final_percent[1],
						sd = sqrt(var),
						cv = sd/mean,
						range = diff(range(rep(day, number), na.rm=T)) + 1,
						rel_range =range /mean,
						iqr = quantile(rep(day, number), 0.75, na.rm = T)- quantile(rep(day, number), 0.25, na.rm = T),
						pct_germ = (Total.germinated_day20[1] /Total[1])*100,
						quartile_coeff_disp= iqr / (quantile(rep(day, number), 0.75, na.rm = T)+ quantile(rep(day, number), 0.25, na.rm = T))
						)

##averaging over the replicates to get summaries per line
germ_summaryPerLinePerTreatment<-germ_summary %>%
  group_by(Magic_line, Hormone, conc_uM, High.or.low.variability) %>%
  summarise (
    mean_median = (mean (median, na.rm =T)),
    mean_mean= mean(mean, na.rm=T),
    mean_mode =mean(mode, na.rm=T),
    mean_cv = mean(cv, na.rm=T),
    mean_percent = mean(final_percent))%>%
    ungroup()

##averaging over the control treatments to get mean stats for a given line across the experiments, just for ordering in terms of average CV
germ_summaryPerMagicLine<-germ_summary %>%
  filter(conc_uM =="0")%>%
  group_by(Magic_line) %>%
  summarise (
    mean_median = (mean (median, na.rm =T)),
    mean_mean= mean(mean, na.rm=T),
    mean_mode =mean(mode, na.rm=T),
    mean_cv = mean(cv, na.rm=T),
    mean_percent = mean(final_percent))%>%
  ungroup()



#### plotting ABA and GA separately, facetting for high and low variability


##GA, % GERM
# png("./GA Percent vs concentrationNoLegend.png", width = 10, height = 5, units = "cm", res = 600)
germ_summaryPerLinePerTreatment%>%
  filter(Hormone =="GA")%>%
  mutate(Magic_line= factor(Magic_line, levels = c("151", "108", "213", "123", "Col", "467", "143", "393", "285", "178", "182", "53"))) %>%
  group_by(Magic_line, Hormone, High.or.low.variability)%>%
  ggplot(., aes(x = as.factor(conc_uM), y =mean_percent, colour=Magic_line, group = (Magic_line))) +
  geom_point()+
  geom_line()+
  facet_grid(.~High.or.low.variability)+
  ylab("Percentage germination")+
  xlab("GA concentration_uM")+
  ylim(0, 100)+
  theme(text = element_text(size=12))+
  theme(legend.position="none")
# dev.off()

##GA, CV

# png("./GA CV vs concentrationNoLegend.png", width = 10, height = 5, units = "cm", res = 600)
germ_summaryPerLinePerTreatment%>%
  filter(Hormone =="GA")%>%
  mutate(Magic_line= factor(Magic_line, levels = c("151", "108", "213", "123", "Col", "467", "143", "393", "285", "178", "182", "53"))) %>%
  group_by(Magic_line, Hormone, High.or.low.variability)%>%
  ggplot(., aes(x = as.factor(conc_uM), y =mean_cv, colour=Magic_line, group = (Magic_line))) +
  geom_point()+
  geom_line()+
  facet_grid(.~High.or.low.variability)+
  ylab("CV")+
  xlab("GA concentration_uM")+
  ylim(0, 0.75)+
  theme(text = element_text(size=12))+
  theme(legend.position="none")
# dev.off()

##GA, mode
# png("./GA Mode vs concentrationNoLegend.png", width = 10, height = 5, units = "cm", res = 600)
germ_summaryPerLinePerTreatment%>%
  filter(Hormone =="GA")%>%
  mutate(Magic_line= factor(Magic_line, levels = c("151", "108", "213", "123", "Col", "467", "143", "393", "285", "178", "182", "53"))) %>%
  group_by(Magic_line, Hormone, High.or.low.variability)%>%
  ggplot(., aes(x = as.factor(conc_uM), y =mean_mode, colour=Magic_line, group = (Magic_line))) +
  geom_point()+
  geom_line()+
  facet_grid(.~High.or.low.variability)+
  ylab("Mode")+
  xlab("GA concentration_uM")+
  ylim(0, 20)+
  theme(text = element_text(size=12))+
  theme(legend.position="none")
# dev.off()

##GA, mean
# png("./GA mean vs concentrationNoLegend.png", width = 10, height = 5, units = "cm", res = 600)
germ_summaryPerLinePerTreatment%>%
  filter(Hormone =="GA")%>%
  mutate(Magic_line= factor(Magic_line, levels = c("151", "108", "213", "123", "Col", "467", "143", "393", "285", "178", "182", "53"))) %>%
  group_by(Magic_line, Hormone, High.or.low.variability)%>%
  ggplot(., aes(x = as.factor(conc_uM), y =mean_mean, colour=Magic_line, group = (Magic_line))) +
  geom_point()+
  geom_line()+
  facet_grid(.~High.or.low.variability)+
  ylab("Mean")+
  xlab("GA concentration_uM")+
  ylim(0, 16)+
  theme(text = element_text(size=12))+
  theme(legend.position="none")
# dev.off()


##ABA, % germ
# png("./ABA Percent vs concentrationNoLegend.png", width = 10, height = 5, units = "cm", res = 600)
germ_summaryPerLinePerTreatment%>%
  filter(Hormone =="ABA")%>%
  mutate(Magic_line= factor(Magic_line, levels = c("151", "108", "213", "123", "Col", "467", "143", "393", "285", "178", "182", "53"))) %>%
  group_by(Magic_line, Hormone, High.or.low.variability)%>%
  ggplot(., aes(x = as.factor(conc_uM), y =mean_percent, colour=Magic_line, group = (Magic_line))) +
  geom_point()+
  geom_line()+
  facet_grid(.~High.or.low.variability)+
  ylab("Percentage germination")+
  xlab("ABA concentration_uM")+
  ylim(0, 100)+
  theme(text = element_text(size=12))+
  theme(legend.position="none")
# dev.off()

##ABA, CV

# png("./ABA CV vs concentrationNoLegend.png", width = 10, height = 5, units = "cm", res = 600)
germ_summaryPerLinePerTreatment%>%
  filter(Hormone =="ABA")%>%
  mutate(Magic_line= factor(Magic_line, levels = c("151", "108", "213", "123", "Col", "467", "143", "393", "285", "178", "182", "53"))) %>%
  group_by(Magic_line, Hormone, High.or.low.variability)%>%
  ggplot(., aes(x = as.factor(conc_uM), y =mean_cv, colour=Magic_line, group = (Magic_line))) +
  geom_point()+
  geom_line()+
  facet_grid(.~High.or.low.variability)+
  ylab("CV")+
  xlab("ABA concentration_uM")+
  ylim(0, 0.75)+
  theme(text = element_text(size=12))+
  theme(legend.position="none")
# dev.off()

##ABA, mode

# png("./ABA Mode vs concentrationNoLegend.png", width = 10, height = 5, units = "cm", res = 600)
germ_summaryPerLinePerTreatment%>%
  filter(Hormone =="ABA")%>%
  mutate(Magic_line= factor(Magic_line, levels = c("151", "108", "213", "123", "Col", "467", "143", "393", "285", "178", "182", "53"))) %>%
  group_by(Magic_line, Hormone, High.or.low.variability)%>%
  ggplot(., aes(x = as.factor(conc_uM), y =mean_mode, colour=Magic_line, group = (Magic_line))) +
  geom_point()+
  geom_line()+
  facet_grid(.~High.or.low.variability)+
  ylab("Mode")+
  xlab("ABA concentration_uM")+
  ylim(0, 20)+
  theme(text = element_text(size=12))+
  theme(legend.position="none")
# dev.off()

##ABA, mean

# png("./ABA mean vs concentrationNoLegend.png", width = 10, height = 5, units = "cm", res = 600)
germ_summaryPerLinePerTreatment%>%
  filter(Hormone =="ABA")%>%
  mutate(Magic_line= factor(Magic_line, levels = c("151", "108", "213", "123", "Col", "467", "143", "393", "285", "178", "182", "53"))) %>%
  group_by(Magic_line, Hormone, High.or.low.variability)%>%
  ggplot(., aes(x = as.factor(conc_uM), y =mean_mean, colour=Magic_line, group = (Magic_line))) +
  geom_point()+
  geom_line()+
  facet_grid(.~High.or.low.variability)+
  ylab("Mean")+
  xlab("ABA concentration_uM")+
  ylim(0, 16)+
  theme(text = element_text(size=12))+
  theme(legend.position="none")
# dev.off()


##output source data for Fig. 7
# write.table(germ_summaryPerLinePerTreatment, './source_data/Figure7_source_data1.tsv',
#             quote = F, row.names = F, sep = "\t")


###making figures for Figure 7, figure supplement 1

germ_longM182Col<- germ_long%>% filter(Magic_line == "Col"|Magic_line == "182")

# png("./Col M182 ABA dose response.png", width = 13, height = 10, units = "cm", res = 600)
germ_longM182Col %>%
  group_by(Magic_line, Experiment)%>%
  filter(Hormone =="ABA")%>%
  ggplot(., aes(x = day, y = perc, fill=as.factor(Experiment))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_grey(start = 0.7, end = 0.3)+
  ylab("Percentage of seeds")+
  xlab("Days to germination")+
  facet_grid(conc_uM ~ Magic_line) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size=12))
# dev.off()


# png("./Col M182 GA dose response.png", width = 13, height = 10, units = "cm", res = 600)
germ_longM182Col %>%
  group_by(Magic_line, Experiment)%>%
  filter(Hormone =="GA")%>%
  ggplot(., aes(x = day, y = perc, fill=as.factor(Experiment))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_grey(start = 0.7, end = 0.3)+
  ylab("Percentage of seeds")+
  xlab("Days to germination")+
  facet_grid(conc_uM ~ Magic_line) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size=12))
# dev.off()

##outputting source data for Figure 7, Figure supplement1
# write.table(germ_longM182Col, './source_data/Figure7_FigureSupplement1_source_data1.tsv',
#             quote = F, row.names = F, sep = "\t")


