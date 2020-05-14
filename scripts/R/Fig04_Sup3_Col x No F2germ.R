################################
# Analysis of Germination data #
# Abley et al 2020             #
################################

library(dplyr) # install.packages("dplyr")
library(tidyr) # install.packages("tidyr")
library(ggplot2)
theme_set(theme_bw())
library(gridExtra)
library(moments)

germ <- read.csv('./data/ColxNo_F2germ.csv')

# Make data long format and tidy it
germ_long <- germ %>%
	gather(day, number, X1:X60) %>%
	arrange(Seed_sample, Rep) %>%
	mutate(day = as.numeric(gsub("X", "", day))) %>%
  filter(!is.na(number)) %>%
  mutate(name = paste(Seed_sample, Rep, sep = "_"))

# Calculate proportion on each day
germ_long <- germ_long %>%
	mutate(perc = number/Total_germ*100)

##output source data
# ColNoSourceData <- germ_long %>%
#   filter(Seed_sample== "Average F2" |Seed_sample== "Col-0"|Seed_sample== "No-0" )
# write.table(ColNoSourceData, './source_data/Figure4_FigureSupplement3_source_data1.tsv',
#             quote = F, row.names = F, sep = "\t")


##make plot for figure 4 figure supplement 3, A
# png("./F2 Col No_withlegendBig.png", width = 18, height = 12, units = "cm", res = 600)
germ_long %>%
  filter(Seed_sample== "Average F2" |Seed_sample== "Col-0"|Seed_sample== "No-0" )%>%
  ggplot(., aes(x = day, y = perc, fill=as.factor(Seed_sample))) + geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette="Set1")+
  xlab("Days to germination") + ylab("Percentage of seeds")+
  theme (text = element_text(size=12))+
  scale_x_continuous(breaks=seq(0, 60, 5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right")
# dev.off()

# png("./F2 Col No_nolegendEarly.png", width = 12, height = 12, units = "cm", res = 600)
germ_long %>%
  filter(Seed_sample== "Average F2" |Seed_sample== "Col-0"|Seed_sample== "No-0" )%>%
  ggplot(., aes(x = day, y = perc, fill=as.factor(Seed_sample))) + geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette="Set1")+
  xlab("Days to germination") + ylab("Percentage of seeds")+
  theme (text = element_text(size=12))+
  xlim(0, 10)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none")
# dev.off()

##make plot for figure 4 figure supplement 3, B
# png("./F2 Col No_nolegendLate.png", width = 12, height = 12, units = "cm", res = 600)
germ_long %>%
  filter(Seed_sample== "Average F2" |Seed_sample== "Col-0"|Seed_sample== "No-0" )%>%
  ggplot(., aes(x = day, y = perc, fill=as.factor(Seed_sample))) + geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette="Set1")+
  xlab("Days to germination") + ylab("Percentage of seeds")+
  theme (text = element_text(size=12))+
  xlim(10, 60)+
  ylim(0, 2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none")
# dev.off()



