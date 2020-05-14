################################
# Analysis of Germination data #
# Abley et al 2020             #
################################

library(dplyr) # install.packages("dplyr")
library(tidyr) # install.packages("tidyr")
library(ggplot2)
theme_set(theme_bw())

##read in the data
germ <- read.csv('./data/cyp707a1_ga3oxMutantsGerm.csv')
germ$id <- 1:nrow(germ)

# Make data long format and tidy it
germ_long <- germ %>%
	gather(day, number, X1:X48) %>%
	arrange(Line, Rep, Experiment) %>%
	mutate(day = as.numeric(gsub("X", " ", day))) %>%
  mutate(name = paste( Line, Rep, Experiment, sep = "_")) %>%
  filter(!is.na(number))

# Calculate proportion on each day
germ_long <- germ_long %>%
	mutate(perc = number/Total*100)

# Calculate mean and variance for each plant's seeds
germ_summary <- germ_long %>%
	group_by(Line, Rep, Sowing_date, Experiment) %>%
	summarise(total_germ = sum(number, na.rm = T),
						mean_day = ifelse(total_germ == 0, 0, sum(day*number, na.rm = T)/total_germ),
						var_day = ifelse (total_germ < 3, NA, sum(number * (day - mean_day)^2, na.rm = T )/total_germ),
						sd_day = sqrt(var_day),
						cv_day = sd_day/mean_day,
						pct_germ = total_germ[1]/Total[1]*100)

germ_summaryPerLinePerExp<-germ_summary %>%
  group_by(Line, Experiment) %>%
  summarise (
    mean_cv = mean(cv_day, na.rm=T),
    mean_percent = mean(pct_germ))%>%
  ungroup()




#plotting all lines across all experiments
germ_long%>%
  ggplot(aes(x = day, y = name, size=ifelse(perc==0, NA, perc), colour=Experiment)) +
  geom_point(stat = "identity") +
  xlim(0,18)+
  theme(
    axis.ticks.y = element_blank(),
    text = element_text(size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position="none") +
  theme(panel.spacing = unit(0, "lines"))+
  xlab("Days to germination")+
  ylab("Percentage germination")+
  facet_grid(Line ~ ., scales = "free_y", space = "free_y", switch ="y")


##plot for Fig 8.
# png("./ga3ox1ga3ox2_cyps_Experiment2.png", width = 12, height =18, units = "cm", res = 600)
germ_long%>%
  filter(Experiment =="2")%>%
  filter(Line != "cyp707a2-1")%>%
  ggplot(aes(x = day, y = perc, fill=as.factor(Rep))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(breaks=seq(0, 20,2), limits = c(0,  20))+
  scale_fill_grey(start = 0.7, end = 0.3)+
  theme(
    axis.ticks.y = element_blank(),
    text = element_text(size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position="none") +
  theme(panel.spacing = unit(0, "lines"))+
  xlab("Days to germination")+
  ylab("Percentage germination")+
  facet_grid(Line ~ .)
# dev.off()

##output source data for Fig. 8
# germLongExp2<-germ_long%>%  filter(Experiment =="2")%>%  filter(Line != "cyp707a2-1")
# write.table(germLongExp2, './source_data/Figure8_source_data1.tsv',
#             quote = F, row.names = F, sep = "\t")


