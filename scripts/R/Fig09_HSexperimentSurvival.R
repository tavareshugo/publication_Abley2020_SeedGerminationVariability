################################
# Analysis of Germination data #
# Abley et al 2020             #
################################

library(dplyr) # install.packages("dplyr")
library(tidyr) # install.packages("tidyr")
library(ggplot2)
theme_set(theme_bw())

germ_long <- read.csv('./data/HeatShockExperimentGerm.csv')

##note that the data contain experiment numbers 2,3,4. In the text of the paper, they are reffered to as experiments 1,2,3 (so 2 here is 1 in the main text and so on)


# Calculate proportion on each day


germ_long <- germ_long %>%
  group_by(Magic_line, HS, Experiment) %>%
  mutate(Total.germinated = sum(number_germ)) %>%
  mutate(Total =Total.germinated + NonGer) %>%
  mutate(percentage_survived = (sum(number_survived_HS)/ sum(number_germ))*100) %>%
  ungroup()

my_magics <- unique(germ_long$Magic_line)

# Calculate mean and variance for each plant's seeds
germ_summaryStats <- germ_long %>%
  group_by(Magic_line,HS, High.or.low.variabiltiy, High_or_low_mode, Experiment) %>%
  summarise(
    median = median(rep(day, number_germ), na.rm=T),
    mean = mean(rep(day, number_germ), na.rm=T),
    mode = mean(day[which(number_germ == max(number_germ))]),
    var =  var(rep(day, number_germ), na.rm=T),
    percent_seedlings_survived= percentage_survived[1],
    sd = sqrt(var),
    cv = sd/mean,
    range = diff(range(rep(day, number_germ), na.rm=T)) + 1,
    rel_range =range /mean,
    iqr = quantile(rep(day, number_germ), 0.75, na.rm = T)- quantile(rep(day, number_germ), 0.25, na.rm = T),
    pct_germ = (Total.germinated[1] /Total[1])*100,
    quartile_coeff_disp= iqr / (quantile(rep(day, number_germ), 0.75, na.rm = T)+ quantile(rep(day, number_germ), 0.25, na.rm = T))
  )

##get summary stats for a given MAGIC line, in a given experiment
germ_summaryPerline <- germ_summaryStats %>%
  group_by(Magic_line,  Experiment) %>%
  summarise(
    mean_cv = mean(cv),
    mean_mode = mean(mode),
    mean_range = mean(range),
    mean_pecent = mean(pct_germ)
  )

##get summary stats for a given MAGIC line, across all experiments
germ_summaryPerlineAverage <- germ_summaryStats %>%
  group_by(Magic_line) %>%
  summarise(
    mean_cv = mean(cv),
    mean_mode = mean(mode),
    mean_range = mean(range),
    mean_pecent = mean(pct_germ)
  )



##get summary stats for a particular line across treatments in a particular experiment. This is used to correlate CV or mode with percent survival in a given experiment, for figure 9 figure supplement 1
germ_summaryPerlineExp3 <- germ_summaryStats %>%
  filter(Experiment ==3)%>%
  group_by(Magic_line, High.or.low.variabiltiy, Experiment)%>%
  summarise(
    mean_cv = mean(cv),
    mean_mode = mean(mode),
    mean_range = mean(range),
    mean_percent = mean(pct_germ)
  )

germ_summaryPerlineExp2 <- germ_summaryStats %>%
  filter(Experiment ==2)%>%
  group_by(Magic_line, High.or.low.variabiltiy, Experiment)%>%
  summarise(
    mean_cv = mean(cv),
    mean_mode = mean(mode),
    mean_range = mean(range),
    mean_percent = mean(pct_germ)
  )

germ_summaryPerlineExp4 <- germ_summaryStats %>%
  filter(Experiment ==4)%>%
  group_by(Magic_line, High.or.low.variabiltiy, Experiment)%>%
  summarise(
    mean_cv = mean(cv),
    mean_mode = mean(mode),
    mean_range = mean(range),
    mean_percent = mean(pct_germ)
  )

##getting table of cv, mode, range, percent, vs survival to HS, for each experiment,
#filtering experiment 2 to just include lines that then were repeated in experiments 3 and 4, this is for figure9 figure supplement 1

statsVsSurvivalExp2 <- germ_summaryStats %>%
  filter(Experiment == "2") %>%
  filter(Magic_line %in% c("108", "188", "461", "458", "Col-0", "151", "200",
                           "393", "4", "178", "305", "285", "351", "304", "182")) %>%
  select (Magic_line, HS, percent_seedlings_survived, Experiment)%>%
  spread(HS, percent_seedlings_survived) %>%
  select (Experiment, Magic_line, Day_0, Day_5, Day_8)%>%
  left_join(germ_summaryPerlineExp2, statsVsSurvivalExp2,
            by= c("Magic_line", "Experiment", "High.or.low.variabiltiy"))


statsVsSurvivalExp3<-germ_summaryStats %>%
  filter(Experiment =="3") %>%
  select (Magic_line, HS, percent_seedlings_survived, Experiment)%>%
  spread(HS, percent_seedlings_survived) %>%
  select (Experiment, Magic_line, Day_0, Day_5, Day_8)%>%
  left_join(germ_summaryPerlineExp3, statsVsSurvivalExp3,
            by= c("Magic_line", "Experiment", "High.or.low.variabiltiy"))


statsVsSurvivalExp4<-germ_summaryStats %>%
  filter(Experiment =="4") %>%
  select (Magic_line, HS, percent_seedlings_survived, Experiment)%>%
  spread(HS, percent_seedlings_survived) %>%
  select (Experiment, Magic_line, Day_0, Day_5, Day_8) %>%
  left_join(germ_summaryPerlineExp4, statsVsSurvivalExp4,
            by= c("Magic_line", "Experiment", "High.or.low.variabiltiy"))

##outputting these tables to the manually combine them in excel, with a composite name for each magic line of experiment_magic line
# write.table(statsVsSurvivalExp2, 'statsVsSurvivalExp2.tsv',
#             quote = F, row.names = F, sep = "\t")

# write.table(statsVsSurvivalExp3, 'statsVsSurvivalExp3.tsv',
#             quote = F, row.names = F, sep = "\t")

# write.table(statsVsSurvivalExp4, 'statsVsSurvivalExp4.tsv',
#             quote = F, row.names = F, sep = "\t")

##input the manually combined table with composite names
compositeSummary<- read.csv('./data/statsVsSurvivalCombined234.csv')


####plots


##box plots summarising survival for Fig 9
##here excluding 492 from exp2 because the % germ was too low. So the lines included here are the same as those in the scatter plots in the supplementary figure.
##the filtering for particular lines is done to make sure that a consistent set of lines is displayed in the graph where possible (see legened of figure 9)


# png("./PlotsForPaper/SurvivalBoxplot_allExperiments.png", width = 11, height = 11, units = "cm", res = 600)
germ_long%>% filter(Magic_line =="108" |Magic_line =="188"|Magic_line =="461"|Magic_line =="458"| Magic_line=="Col-0"| Magic_line=="151" | Magic_line=="200"| Magic_line=="393"| Magic_line=="4"| (Experiment =="3" & Magic_line=="492")| (Experiment =="4" & Magic_line=="492")| Magic_line=="178"| Magic_line=="305"| Magic_line=="285"| Magic_line=="351"| Magic_line=="304"| Magic_line=="182")%>%
  filter(HS != "Day_11")%>%
  ggplot(aes(HS, percentage_survived)) +
  geom_boxplot()+
  facet_grid(High.or.low.variabiltiy ~ Experiment)+
  xlab("Heat shock treatment")+
  ylab("Percentage of seedlings that survived")+
  theme(text = element_text(size=12))
# dev.off()

##outputting the source data for this
# percentageSurvivedSelectedLines<-germ_long%>% filter(Magic_line =="108" |Magic_line =="188"|Magic_line =="461"|Magic_line =="458"| Magic_line=="Col-0"| Magic_line=="151" | Magic_line=="200"| Magic_line=="393"| Magic_line=="4"| (Experiment =="3" & Magic_line=="492")| (Experiment =="4" & Magic_line=="492")| Magic_line=="178"| Magic_line=="305"| Magic_line=="285"| Magic_line=="351"| Magic_line=="304"| Magic_line=="182")%>% filter(HS != "Day_11")
# write.table(percentageSurvivedSelectedLines, './source_data/Figure9_source_data1.tsv',
#             quote = F, row.names = F, sep = "\t")



## example of histograms coloured according to survival for Figure 9b

#png("./PlotsForPaper/Col-0vsM812.png", width = 11, height = 11, units = "cm", res = 600)
germ_long %>%
  filter(Experiment==3 & (Magic_line =="Col-0" | Magic_line =="182"))%>%
  select(Magic_line, HS, Experiment, day, number_survived_HS, number_died_HS, Total) %>%
  gather("survived", "count", number_died_HS, number_survived_HS) %>%
  mutate(pct = (count/Total)*100,
         survived = ifelse(survived == "number_died_HS", "no", "yes")) %>%
  mutate(survived = factor(survived, level = c("yes", "no"))) %>%
  mutate(Magic_line = factor(Magic_line, level = c("Col-0", "182"))) %>%
  ggplot(aes(day, pct, fill = survived)) +
  geom_col() +
  facet_grid(HS ~Magic_line) +
  scale_fill_manual(values = c("blue", "red")) +
  xlab("Day")+
  ylab("Percent germinated")+
  theme(text = element_text(size=12))
#dev.off()

##output source data file for this plot
# ColM182ColouredSurvivalPlot<-germ_long %>%
#   filter(Experiment==3 & (Magic_line =="Col-0" | Magic_line =="182"))%>%
#   select(Magic_line, HS, Experiment, day, number_survived_HS, number_died_HS, Total) %>%
#   gather("survived", "count", number_died_HS, number_survived_HS) %>%
#   mutate(pct = (count/Total)*100,
#          survived = ifelse(survived == "number_died_HS", "no", "yes")) %>%
#   mutate(survived = factor(survived, level = c("yes", "no"))) %>%
#   mutate(Magic_line = factor(Magic_line, level = c("Col-0", "182")))
# write.table(ColM182ColouredSurvivalPlot, './source_data/Figure9_source_data2.tsv',
#             quote = F, row.names = F, sep = "\t")


##Figure 9 figure supplement 1
###plots of survival vs summary stats, for individual experiments. here the summary stats are summaries of a given line across all the treatments in a given experiment.

##for experiment 2, note that experiment 2 is labelled as experiment 1 in the sup fig


##Panels A
##CV vs % survival on D5
##filtering to only include the lines that are included in Fig9 box plots (this excludes lines with lowest percentage germination)

statsVsSurvivalExp2<-statsVsSurvivalExp2 %>%filter(Magic_line != "336" & Magic_line != "311"& Magic_line != "53"& Magic_line != "492"& Magic_line != "461"& Magic_line != "188"& Magic_line != "458")

# png("./PlotsForPaper/Experiment2_CV vs percent survivedD5.png", width = 5, height = 5, units = "cm", res = 600)
statsVsSurvivalExp2%>%
  ggplot(., aes(x = mean_cv, y =Day_5)) +
  geom_point(size=2)+
  theme(legend.position="none")+
  ylab("Percentage seedlings alive")+
  xlab("CV of germination time")+
  theme(text = element_text(size=9))
# dev.off()

cor.test(statsVsSurvivalExp2$mean_cv, statsVsSurvivalExp2$Day_5)

##output source data file for this

# write.table(statsVsSurvivalExp2, './source_data/Figure9_Figure_Supplement1_source_data1.tsv',
#             quote = F, row.names = F, sep = "\t")

##Mode vs % survival on D5

# png("./PlotsForPaper/Experiment2_Mode vs percent survivedD5.png", width = 5, height = 5, units = "cm", res = 600)
statsVsSurvivalExp2%>%
  ggplot(., aes(x = mean_mode, y =Day_5)) +
  geom_point(size=2)+
  theme(legend.position="none")+
  ylab("Percentage seedlings alive")+
  xlab("Mode germination time")+
  theme(text = element_text(size=9))
# dev.off()

cor.test(statsVsSurvivalExp2$mean_mode, statsVsSurvivalExp2$Day_5)

##panels B

##experiment 2, day 8##labelled as exp1, day 8 in the sup fig
##CV vs % survival on D8

# png("./PlotsForPaper/Experiment2_CV vs percent survivedD8.png", width = 5, height = 5, units = "cm", res = 600)
statsVsSurvivalExp2%>%
  ggplot(., aes(x = mean_cv, y =Day_8)) +
  geom_point(size=2)+
  theme(legend.position="none")+
  ylab("Percentage seedlings alive")+
  xlab("CV of germination time")+
  theme(text = element_text(size=9))
# dev.off()

cor.test(statsVsSurvivalExp2$mean_cv, statsVsSurvivalExp2$Day_8)


##Mode vs % survival on D8

# png("./PlotsForPaper/Experiment2_Mode vs percent survivedD8.png", width = 5, height = 5, units = "cm", res = 600)
statsVsSurvivalExp2%>%
  ggplot(., aes(x = mean_mode, y =Day_8)) +
  geom_point(size=2)+
  theme(legend.position="none")+
  ylab("Percentage seedlings alive")+
  xlab("Mode germination time")+
  theme(text = element_text(size=9))
# dev.off()

cor.test(statsVsSurvivalExp2$mean_mode, statsVsSurvivalExp2$Day_8)


##Panels C
##for experiment 3 ##exp3 is labelled as exp 2 in the sup fig

##write the source data
# write.table(statsVsSurvivalExp3, './source_data/Figure9_Figure_Supplement1_source_data2.tsv',
#             quote = F, row.names = F, sep = "\t")


##CV vs % survival on D5
# png("./PlotsForPaper/Experiment3_CV vs percent survivedD5.png", width = 5, height = 5, units = "cm", res = 600)
statsVsSurvivalExp3%>%
  ggplot(., aes(x = mean_cv, y =Day_5)) +
  geom_point(size=2)+
  theme(legend.position="none")+
  ylab("Percentage seedlings alive")+
  xlab("CV of germination time")+
  theme(text = element_text(size=9))+
  ylim(0,100)
# dev.off()

cor.test(statsVsSurvivalExp3$mean_cv, statsVsSurvivalExp3$Day_5)

##Mode vs % survival on D5

# png("./PlotsForPaper/Experiment3_Mode vs percent survivedD5.png", width = 5, height = 5, units = "cm", res = 600)
statsVsSurvivalExp3%>%
  ggplot(., aes(x = mean_mode, y =Day_5)) +
  geom_point(size=2)+
  theme(legend.position="none")+
  ylab("Percentage seedlings alive")+
  xlab("Mode germination time")+
  theme(text = element_text(size=9))+
  ylim(0,100)
# dev.off()

cor.test(statsVsSurvivalExp3$mean_mode, statsVsSurvivalExp3$Day_5)


##for experiment 3 day 8##exp3 is labelled as exp 2 in the sup fig
##panels D

##CV vs % survival on D8
# png("./PlotsForPaper/Experiment3_CV vs percent survivedD8.png", width = 5, height = 5, units = "cm", res = 600)
statsVsSurvivalExp3%>%
  ggplot(., aes(x = mean_cv, y =Day_8)) +
  geom_point(size=2)+
  theme(legend.position="none")+
  ylab("Percentage seedlings alive")+
  xlab("CV of germination time")+
  theme(text = element_text(size=9))+
  ylim(0,75)
# dev.off()

cor.test(statsVsSurvivalExp3$mean_cv, statsVsSurvivalExp3$Day_8)

##Mode vs % survival on D8

# png("./PlotsForPaper/Experiment3_Mode vs percent survivedD8.png", width = 5, height = 5, units = "cm", res = 600)
statsVsSurvivalExp3%>%
  ggplot(., aes(x = mean_mode, y =Day_8)) +
  geom_point(size=2)+
  theme(legend.position="none")+
  ylab("Percentage seedlings alive")+
  xlab("Mode germination time")+
  theme(text = element_text(size=9))+
  ylim(0,75)
# dev.off()


cor.test(statsVsSurvivalExp3$mean_mode, statsVsSurvivalExp3$Day_8)


##for experiment 4, which is labelled as exp3 in the sup fig

##write the source data
# write.table(statsVsSurvivalExp4, './source_data/Figure9_Figure_Supplement1_source_data3.tsv',
#             quote = F, row.names = F, sep = "\t")

##Panels E

##CV vs % survival on D5
# png("./PlotsForPaper/Experiment4_CV vs percent survivedD5.png", width = 5, height = 5, units = "cm", res = 600)
statsVsSurvivalExp4%>%
  ggplot(., aes(x = mean_cv, y =Day_5)) +
  geom_point(size=2)+
  theme(legend.position="none")+
  ylab("Percentage seedlings alive")+
  xlab("CV of germination time")+
  theme(text = element_text(size=9))+
  ylim(0,100)
# dev.off()

cor.test(statsVsSurvivalExp4$mean_cv, statsVsSurvivalExp4$Day_5)


##Mode vs % survival on D5

# png("./PlotsForPaper/Experiment4_Mode vs percent survivedD5.png", width = 5, height = 5, units = "cm", res = 600)
statsVsSurvivalExp4%>%
  ggplot(., aes(x = mean_mode, y =Day_5)) +
  geom_point(size=2)+
  theme(legend.position="none")+
  ylab("Percentage seedlings alive")+
  xlab("Mode germination time")+
  theme(text = element_text(size=9))+
  ylim(0,100)
# dev.off()

cor.test(statsVsSurvivalExp4$mean_mode, statsVsSurvivalExp4$Day_5)

##panels F
##CV vs % survival on D8
# png("./PlotsForPaper/Experiment4_CV vs percent survivedD8.png", width = 5, height = 5, units = "cm", res = 600)
statsVsSurvivalExp4%>%
  ggplot(., aes(x = mean_cv, y =Day_8)) +
  geom_point(size=2)+
  theme(legend.position="none")+
  ylab("Percentage seedlings alive")+
  xlab("CV of germination time")+
  theme(text = element_text(size=9))+
  ylim(0,20)
# dev.off()

cor.test(statsVsSurvivalExp4$mean_cv, statsVsSurvivalExp4$Day_8)

##Mode vs % survival on D8

# png("./PlotsForPaper/Experiment4_Mode vs percent survivedD8.png", width = 5, height = 5, units = "cm", res = 600)
statsVsSurvivalExp4%>%
  ggplot(., aes(x = mean_mode, y =Day_8)) +
  geom_point(size=2)+
  theme(legend.position="none")+
  ylab("Percentage seedlings alive")+
  xlab("Mode germination time")+
  theme(text = element_text(size=9))+
  ylim(0,20)
# dev.off()

cor.test(statsVsSurvivalExp4$mean_mode, statsVsSurvivalExp4$Day_8)


