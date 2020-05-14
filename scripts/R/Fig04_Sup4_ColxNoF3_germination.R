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


## read in data
germ <- read.csv('./data/ColxNo_F3germ.csv')
germ$id <- 1:nrow(germ)

# Make data long format and tidy it
germ_long <- germ %>%
	gather(day, number, D1:D10) %>%
	arrange(F3_line, Day_parent_germ, DAR, Early.or.late.germinator) %>%
	mutate(day = as.numeric(gsub("D", "", day))) %>%
  filter(!is.na(number)) %>%
  mutate(name = paste(F3_line, Day_parent_germ, sep = "_"))


# Calculate proportion on each day
germ_long <- germ_long %>%
	mutate(perc = number/Total*100, perc_total_ger = number/Total_germ*100)


# Calculate mean and variance for each plant's seeds
germ_summary <- germ_long %>%
	group_by(F3_line, Day_parent_germ, name) %>%
	summarise(
	          Early.or.late.germinator= Early.or.late.germinator[1],
	          median = median(rep(day, number), na.rm=T),
						mean = mean(rep(day, number), na.rm=T),
						mode = mean(day[which(number == max(number))]),
						var =  var(rep(day, number), na.rm=T),
						sd = sqrt(var),
						cv = sd/mean,
						range = diff(range(rep(day, number), na.rm=T)) + 1,
						rel_range =range /mean,
						iqr = quantile(rep(day, number), 0.75, na.rm = T)- quantile(rep(day, number), 0.25, na.rm = T),
						quartile_coeff_disp= iqr / (quantile(rep(day, number), 0.75, na.rm = T)+ quantile(rep(day, number), 0.25, na.rm = T)),
						pct_germ = (Total_germ[1] /Total[1])*100
						)

##summary stats averaged over all early vs all late germinators
germ_summary_EarlyVsLate<- germ_summary %>%
  group_by(Early.or.late.germinator) %>%
  summarise(
    mean_cv = mean(cv), na.rm=T,
    mean_mean = mean(mean), na.rm=T,
    mean_mode = mean(mode),
    mean_pct =  mean(pct_germ)
   )

##make plots for Figure4, figure supplement 4, all panels

# png("./F3 CV vs parent day germ.png", width = 8.9, height = 7.6, units = "cm", res = 600)
qplot(Day_parent_germ, cv, data = germ_summary, size=I(2)) + xlim(0, 60) + ylim(0, 0.5)+
theme(text = element_text(size=12))+
xlab("Day of F2 parent's germination") +
ylab("F3 CV of germination time")
# dev.off()

# png("./F3 Mean vs parent day germ.png", width = 8.9, height = 7.6, units = "cm", res = 600)
qplot(Day_parent_germ, mean, data = germ_summary, size=I(2)) + xlim(0, 60) + ylim(0, 5)+
theme(text = element_text(size=12))+
xlab("Day of F2 parent's germination") +
ylab("F3 Mean germination time")
# dev.off()


# png("./F3 Mode vs parent day germ.png", width = 8.9, height = 7.6, units = "cm", res = 600)
germ_summary%>%
ggplot(., aes(x = Day_parent_germ, y=mode)) +
  geom_jitter(height =0.01) +
  theme(text = element_text(size=12))+
  xlab("Day of F2 parent's germination") +
  ylab("F3 Mode germination time")+
  ylim(1, 5)
# dev.off()


# png("./F3 Percent vs parent day germ.png", width = 8.9, height = 7.6, units = "cm", res = 600)
germ_summary%>%
  ggplot(., aes(x = Day_parent_germ, y=pct_germ)) +
  geom_jitter(height =0.00, width=0) +
  theme(text = element_text(size=12))+
  xlab("Day of F2 parent's germination") +
  ylab("F3 Percent germination")+
  ylim(0, 100)
# dev.off()

##output source data
# write.table(germ_summary, './source_data/Figure4_FigureSupplement4_source_data1.tsv',
#             quote = F, row.names = F, sep = "\t")


####do tests

germ_summaryCV <- germ_summary %>% ungroup() %>%
  select(name, Early.or.late.germinator, cv) %>%
  spread(Early.or.late.germinator, cv)

res <- wilcox.test(germ_summaryCV$Early, germ_summaryCV$Late)
res

germ_summaryMode<-germ_summary%>%ungroup()%>%select (name, Early.or.late.germinator, mode)%>%spread(Early.or.late.germinator, mode)
res2 <- wilcox.test(germ_summaryMode$Early, germ_summaryMode$Late)
res2

germ_summaryMean<-germ_summary%>%ungroup()%>%select (name, Early.or.late.germinator, mean)%>%spread(Early.or.late.germinator, mean)
res3 <- wilcox.test(germ_summaryMean$Early, germ_summaryMean$Late)
res3

germ_summaryPercent<-germ_summary%>%ungroup()%>%select (name, Early.or.late.germinator, pct_germ)%>%spread(Early.or.late.germinator, pct_germ)
res4 <- wilcox.test(germ_summaryPercent$Early, germ_summaryPercent$Late)
res4

