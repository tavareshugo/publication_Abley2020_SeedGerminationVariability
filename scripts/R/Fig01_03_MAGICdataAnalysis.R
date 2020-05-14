################################
# Analysis of Germination data #
# Abley et al 2020             #
################################

#
# Setup ----
#
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_bw())

#
# Read data ----
#
# MAGIC data
germ_long <- read.csv('./data/MAGICs.csv', stringsAsFactors = FALSE)

# MAGIC parents data
germ_longParents<- read.csv('./data/MAGICParents.csv', stringsAsFactors = FALSE)

# Accessions data
germ_longAccessions<- read.csv('./data/SpanishAccessions.csv', stringsAsFactors = FALSE)

##soil vs plates comparison (Figure1 figure supplement 1 G)

germ_longSoilvsPlates<-read.csv('./data/Soil_vs_plates.csv', stringsAsFactors = FALSE)

####MAGIC lines#####

#filter out plates where there was less than 10 seeds germinated
germ_long<-germ_long %>%
  filter(Total.germinated > 10)

# summarise MAGIC germination data ----
#
# Calculate summary stats for each plant's seeds
germ_summary <- germ_long %>%
  group_by(Magic_line, Plant, MeanDAR, id, Mother.plant.sowing.number) %>%
  summarise(
    median = median(rep(day, number), na.rm=T),
    mean = mean(rep(day, number), na.rm=T),
    mode = mean(day[which(number == max(number))]),
    var =  var(rep(day, number), na.rm=T),
    sd = sqrt(var),
    cv = sd/mean,
    ls = mean(abs(rep(day, number) - mean)/mean, na.rm = TRUE),
    range = diff(range(rep(day, number), na.rm=T)) + 1,
    rel_range =range /mean,
    iqr = quantile(rep(day, number), 0.75, na.rm = T)- quantile(rep(day, number), 0.25, na.rm = T),
    pct_germ = (Total.germinated[1] /Total[1])*100
  ) %>%
  ungroup()


# Calculate the mean for each MAGIC line in a given experiment
# (mean over three batches of seeds from separate parent plants)
germ_summaryPerExp <- germ_summary %>%
  group_by(Magic_line, MeanDAR, Mother.plant.sowing.number) %>%
  summarise (
    mean_median = (mean (median, na.rm =T)),
    mean_mean= mean(mean, na.rm=T),
    mean_mode =mean(mode, na.rm=T),
    mean_var = mean(var, na.rm=T),
    mean_cv = mean(cv, na.rm=T),
    mean_ls=mean(ls, na.rm=T),
    mean_sd = mean(sd, na.rm=T),
    mean_range = mean(range, na.rm=T),
    mean_rel_range = mean(rel_range, na.rm=T),
    mean_iqr= mean(iqr, na.rm=T),
    mean_percent = mean(pct_germ, na.rm=T),
    n = n(),
    variance_cv = ifelse (n <= 2, NA, sum ((cv - mean_cv)^2, na.rm = T )/n),
    sd_cv= sqrt(variance_cv),
    sd_meanDay = ifelse (n <= 2, NA, sqrt(sum ((mean- mean_mean)^2, na.rm = T )/n)),
    sd_meanPercent = ifelse (n <= 2, NA, sqrt(sum ((pct_germ - mean_percent)^2, na.rm = T )/n))
  ) %>%
  ungroup()



# Generate a table with one entry for the summary stats for each line,
# to make scatter plots with one point per line and for QTL mapping.
# For some lines in the above table, there is more than one experiment
# So this does an average over different experiments (different Mother plant sowing numbers)
germ_summaryPerLine <- germ_summaryPerExp %>%
  group_by(Magic_line) %>%
  summarise (
    meanDAR = (mean (MeanDAR, na.rm =T)),
    mean_median = (mean (mean_median, na.rm =T)),
    mean_mean= mean(mean_mean, na.rm=T),
    mean_mode =mean(mean_mode, na.rm=T),
    mean_var = mean(mean_var, na.rm=T),
    mean_cv = mean(mean_cv, na.rm=T),
    mean_ls = mean(mean_ls, na.rm=T),
    mean_sd = mean(mean_sd, na.rm=T),
    mean_range = mean(mean_range, na.rm=T),
    mean_rel_range = mean(mean_rel_range, na.rm=T),
    mean_iqr= mean(mean_iqr, na.rm=T),
    n = n(),
    mean_percent = mean(mean_percent, na.rm=T)
  ) %>%
  ungroup()

# Add column that indicates these are the MAGIC lines
germ_summaryPerLine <- germ_summaryPerLine %>% mutate (Group="MAGIC")


# # output for QTL mapping
# write.table(germ_summaryPerLine, './data/QTLmapping/germ_summaryPerLineForQTLMapping.tsv',
#             quote = FALSE, row.names = FALSE, sep = "\t")



#
# summarise MAGIC parents ----

#filter out plates where there was less than 10 seeds germinated
germ_longParents<-germ_longParents %>%
  filter(Total.germinated > 10)

# Calculate summary stats for each plant's seeds
germ_summaryParents <- germ_longParents %>%
  group_by(Magic_Parent_line, Plant, Sowing_round, Parent_vernalised, Mean_DAR) %>%
  summarise(
    median = median(rep(day, number), na.rm=T),
    mean = mean(rep(day, number), na.rm=T),
    mode = mean(day[which(number == max(number))]),
    var =  var(rep(day, number), na.rm=T),
    sd = sqrt(var),
    cv = sd/mean,
    ls = mean(abs(rep(day, number) - mean)/mean, na.rm = TRUE),
    range = diff(range(rep(day, number), na.rm=T)) + 1,
    rel_range =range /mean,
    iqr = quantile(rep(day, number), 0.75, na.rm = T)- quantile(rep(day, number), 0.25, na.rm = T),
    pct_germ = (Total.germinated[1] /Total[1])*100
  ) %>%
  ungroup()

# Calculate the mean for each MAGIC parent line in a given experiment
# (mean over three batches of seeds from separate parent plants)
germ_summaryPerExpParents <- germ_summaryParents %>%
  group_by(Magic_Parent_line, Sowing_round, Parent_vernalised, Mean_DAR) %>%
  summarise(
    meanDAR = (mean (Mean_DAR, na.rm =T)),
    mean_median = (mean (median, na.rm =T)),
    mean_mean= mean(mean, na.rm=T),
    mean_mode =mean(mode, na.rm=T),
    mean_var = mean(var, na.rm=T),
    mean_cv = mean(cv, na.rm=T),
    mean_ls=mean(ls, na.rm=T),
    mean_sd = mean(sd, na.rm=T),
    mean_range = mean(range, na.rm=T),
    mean_rel_range = mean(rel_range, na.rm=T),
    mean_iqr= mean(iqr, na.rm=T),
    mean_percent = mean(pct_germ, na.rm=T)
  )  %>%
  ungroup()


# write.table(germ_summaryPerExpParents, './Data/germ_summaryPerExpParents.tsv',
#             quote = F, row.names = F, sep = "\t")

# Add column indicating that these are the MAGIC parental lines
germ_summaryPerExpParents <- germ_summaryPerExpParents %>%
  mutate (Group ="Parent") %>%
  rename(Magic_line = Magic_Parent_line)

# Select only a few columns
germ_summaryPerExpParents <- germ_summaryPerExpParents %>%
  select(Magic_line, Sowing_round, meanDAR, mean_median, mean_mean, mean_mode, mean_var, mean_cv,
         mean_ls, mean_sd, mean_range, mean_rel_range, mean_iqr, mean_percent, Group)


#
# summarising accessions ----
#
germ_longAccessions<-germ_longAccessions%>%filter(Total.germinated > 10)

germ_summaryAccessions <- germ_longAccessions %>%
  group_by(Accession, Plant, Sowing_date, MeanDAR) %>%
  summarise(
    median = median(rep(day, number), na.rm=T),
    mean = mean(rep(day, number), na.rm=T),
    mode = mean(day[which(number == max(number))]),
    var =  var(rep(day, number), na.rm=T),
    sd = sqrt(var),
    cv = sd/mean,
    ls = mean(abs(rep(day, number) - mean)/mean, na.rm = TRUE),
    range = diff(range(rep(day, number), na.rm=T)) + 1,
    rel_range =range /mean,
    iqr = quantile(rep(day, number), 0.75, na.rm = T)- quantile(rep(day, number), 0.25, na.rm = T),
    pct_germ = (Total.germinated[1] /Total[1])*100
  ) %>%
  ungroup()

# Calculate the mean for each Accession in a given experiment
# (mean over three batches of seeds from separate parent plants)
germ_summaryPerExpAccessions<-germ_summaryAccessions %>%
  group_by(Accession) %>%
  summarise (
    meanDAR = (mean (MeanDAR, na.rm =T)),
    mean_median = (mean (median, na.rm =T)),
    mean_mean= mean(mean, na.rm=T),
    mean_mode =mean(mode, na.rm=T),
    mean_var = mean(var, na.rm=T),
    mean_cv = mean(cv, na.rm=T),
    mean_ls = mean(ls, na.rm=T),
    mean_sd = mean(sd, na.rm=T),
    mean_range = mean(range, na.rm=T),
    mean_rel_range = mean(rel_range, na.rm=T),
    mean_iqr= mean(iqr, na.rm=T),
    mean_percent = mean(pct_germ, na.rm=T)
  )  %>%
  ungroup()

# Add column indicating these are Accessions
germ_summaryPerExpAccessions <- germ_summaryPerExpAccessions %>%
  mutate (Group ="Accession") %>%
  rename(Magic_line = Accession)



#
# Combine data ----
#
# combine the MAGIC parents data, spanish accessions and the MAGIC data into one data frame for histogram for Fig.1B. Only using MAGIC parents with DAR of 30, so just sowing round 1

germ_summaryPerExpParentsDAR30<- germ_summaryPerExpParents %>%
  filter(Sowing_round != "2" & Sowing_round != "3")

germ_summary_All <- bind_rows(germ_summaryPerLine,
                              germ_summaryPerExpParentsDAR30,
                              germ_summaryPerExpAccessions)

# write.table(germ_summary_All, './source_data/Figure1_source_data2.tsv',
#                      quote = F, row.names = F, sep = "\t")

# Combine accession data only (MAGIC parents and Spanish accessions)
germ_summary_accessionsOnly <- bind_rows(germ_summaryPerExpParentsDAR30, germ_summaryPerExpAccessions)



#
# Plotting ----
#
# Fig.1B histograms - just including MAGICs and parents
#png("./Plots/MAGICandParentsHistogram.png", width = 8.9, height = 9, units = "cm", res = 600)
germ_summary_All%>%
  filter (Group != "Accession") %>%
  ggplot(., aes(mean_cv, fill=Group)) +
  geom_histogram(colour ="black")+
  theme_bw()+
  theme(text = element_text(size=12))+
  ylab("Number of lines")+
  xlab("CV")+
  scale_fill_manual(values=c("#56B4E9", "#E69F00"))+
  theme(legend.position="top")
#dev.off()


# including, MAGIC, accessions and parents
# png("./plots/MAGICandParentsandAccessionsHistogram.png", width = 8.9, height = 9, units = "cm", res = 600)
germ_summary_All %>%
  mutate(Group=str_replace(Group, "Accession", "Spanish\nAccession"))%>%
  mutate(Group=str_replace(Group, "Parent", "Parental\nAccession"))%>%
  mutate(Group=str_replace(Group, "MAGIC", "MAGIC\nline"))%>%
  mutate(Group = factor(Group, levels = c("MAGIC\nline", "Parental\nAccession", "Spanish\nAccession"))) %>%
  ggplot(., aes(mean_cv, fill=Group)) +
  geom_histogram(colour ="black", bins = 50)+
  theme_bw()+
  theme(text = element_text(size=12))+
  ylab("Number of lines")+
  xlab("CV")+
  facet_grid(Group~., scales = "free_y")+
  scale_fill_manual(values=c("#56B4E9", "#E69F00", "purple"))+
  theme(legend.position="none")
# dev.off()


#
# Distribution of selected lines ----
#
# Filter data for plotting distributions of selected lines for Fig 1A
# MAGIC parents
GermLongParentsPlotting <- germ_longParents %>%
  filter((Magic_Parent_line=="Col-0" & Mean_DAR==33) |
           (Magic_Parent_line=="Bur-0" & Mean_DAR==30) |
           (Magic_Parent_line=="Ws-0" & Mean_DAR==29) |
           (Magic_Parent_line=="Sf-2" & Mean_DAR==29)) %>%
  rename(Line = Magic_Parent_line) %>%
  select(Line, day, perc, name, group, RepNumber)

# Spannish accessions
GermLongAccessionPlotting <- germ_longAccessions %>%
  filter((Accession =="Aul-0" & MeanDAR ==27) |
           (Accession =="Mad-0" & MeanDAR ==27) |
           (Accession =="Svi-0" & MeanDAR ==37)) %>%
  rename( Line = Accession) %>%
  select(Line, day, perc, name, group,  RepNumber)

# MAGIC lines
GermLongPlotting <- germ_long %>%
  filter((Magic_line=="M108" & (Plant=="1"|Plant=="2"| Plant=="4")) |
           (Magic_line=="M203" & MeanDAR==33) |
           (Magic_line=="M393") |
           (Magic_line=="M285") |
           (Magic_line=="M182" & MeanDAR==31) |
           (Magic_line=="M178" & MeanDAR==34 & (Plant=="1"|Plant=="2"| Plant=="3"))) %>%
  rename(Line = Magic_line) %>%
  select(Line, day, perc, name, group,  RepNumber)

# Combine all tables
ForPlottingCombined <- bind_rows(GermLongParentsPlotting, GermLongAccessionPlotting, GermLongPlotting)

# ##output this table as source data for Fig1A
# write.table(ForPlottingCombined, './source_data/Figure1_source_data1.tsv',
#             quote = F, row.names = F, sep = "\t")


# Make Fig. 1A
#png("./Plots/PNGs/CV distribution and comparison/AccessionsAndMAGICsCombinedDistributionsNoLegendJul2018.png", width = 18, height =15, units = "cm", res = 600)
ForPlottingCombined %>%
  mutate(Line = factor(Line,
                       levels = c("Col-0", "Bur-0", "Sf-2", "Ws-0", "Svi-0",
                                  "Mad-0", "Aul-0", "M108", "M203", "M305", "M393",
                                  "M285", "M182", "M178"))) %>%
  mutate(group = factor(group, levels = c("accessions", "magic parents", "magic"))) %>%
  ggplot(aes(x = day, y = name, size=ifelse(perc==0, NA, perc), colour=group)) +
  geom_point(stat = "identity") +
  scale_x_continuous(breaks=seq(0, 50,5), limits = c(0,  50)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="right") +
  theme(panel.spacing = unit(0, "lines"))+
  xlab("Days to germination")+
  ylab("Line")+
  theme(legend.position="none")+
  facet_grid(group + Line ~ ., scales = "free_y", space = "free_y", switch ="y")
#dev.off()



#
# Experiment comparison ----
#
# Make table with summary stats for each experiment for each line,
# for lines where there are repeats.
# This is to make the scatter plot for Fig.1C.
germ_summaryForExpComparison <- germ_summaryPerExp

# write.table(germ_summaryForExpComparison, './data/germ_summaryForExpComparison.tsv', quote = F, row.names = F, sep = "\t")


# Outputted table and then manually deleted lines where there is only one experiment.
# Manually edited to include experiment order, for plotting exp1 vs exp2.
# Choosing sowings where DAR is closest to 30.
# If there are more than 2 experiments for a line, where one is experiment 1,
# I remove experiment 1 as lines in this experiment had overall lower % germ.
germ_summaryForExpComparisonManuallyFiltered <- read.csv('./data/germ_summaryForExpComparison.csv',
                                                         stringsAsFactors = FALSE)

# now want to get different experiments in columns, rather than rows.
germ_summaryForExpComparison <- germ_summaryForExpComparisonManuallyFiltered %>%
  spread(Exp_order, mean_cv) %>%
  rename(CV_exp1 = `1`, CV_exp2 = `2`, CV_exp3 = `3`)

#
germ_summaryForExpComparisonCV <- germ_summaryForExpComparison %>%
  select(Magic_line, CV_exp1, CV_exp2, CV_exp3)

#now want to manually filter this to just have 1 row for each MAGIC line
# write.table(germ_summaryForExpComparisonCV, './data/germ_summaryForExpComparisonCV.tsv', quote = F, row.names = F, sep = "\t")

#this is ready for scatter plot of CV exp1 vs exp2
germ_summaryForExpComparisonCVManuallyEdited <-read.csv('./data/germ_summaryForExpComparisonCV.csv', stringsAsFactors = FALSE)

##write this to source data folder, it is the source data for Fig1C
# write.table(germ_summaryForExpComparisonCVManuallyEdited, './source_data/Figure1_source_data3.tsv', quote = F, row.names = F, sep = "\t")

# make plot for Fig. 1C
#png("./plots/CVexp1vsexp2CorrectedStraightLiney=x.png", width = 8.9, height = 7.6, units = "cm", res = 600)
ggplot(germ_summaryForExpComparisonCVManuallyEdited, aes(CV_exp1, CV_exp2)) +
  geom_point() +
  xlim(0, 1.5) + ylim(0, 1.5) +
  theme(text = element_text(size=12))+
  xlab("CV experiment 1") +
  ylab("CV experiment 2") +
  geom_abline(intercept=0.00, slope=1)
#dev.off()

# doing Pearson's rank correlation between two experiments, for either just lines with CV less than 0.5, or for all lines for which repeats are available
germ_summaryOnlyUnimodal <- germ_summaryForExpComparisonCVManuallyEdited %>%
  filter(CV_exp1 < 0.5 & CV_exp2 <0.5)

cor.test(germ_summaryForExpComparisonCVManuallyEdited$CV_exp1, germ_summaryForExpComparisonCVManuallyEdited$CV_exp2)

cor.test(germ_summaryOnlyUnimodal$CV_exp1, germ_summaryOnlyUnimodal$CV_exp2)


##Fig1_FigureSupplement1 Magic parents DAR comparison

germ_summaryPerExpParentsForDARcomparison<-germ_summaryPerExpParents%>%spread(Sowing_round, mean_cv)

##tidying data in excel to make columns with CV for sowing 1 (~30DAR), sowing 2 (~45DAR) and sowing 3 (~60DAR)
# write.table(germ_summaryPerExpParentsForDARcomparison, './data/MAGICParents_ResultsSummarySowingTimes.tsv', quote = F, row.names = F, sep = "\t")

SummaryCVSowingTimes<-read.csv('./data/SummaryCVSowingTimes.csv')

# write this table as source data for figure Figure1_Figure_Supplement1
# write.table(SummaryCVSowingTimes, './source_data/Figure1_FigureSupplement1_source_data1.tsv', quote = F, row.names = F, sep = "\t")


# png("./figures/new_parts_Feb2020/CVMagicParents 30 days vs 60 days.png", width = 10, height =10, units = "cm", res = 600)
SummaryCVSowingTimes%>%
  ggplot(., aes(x = Sowing1, y=Sowing3)) +
  geom_point(size=3.5) +
  scale_colour_manual(values = rainbow(25)) +
  xlab("~30 days") + ylab("~60 days") +
  theme (text = element_text(size=12))+
  theme(legend.position = "none") +
  ylim(0, 0.35) +
  xlim (0, 0.35)+
  geom_abline(a=0, b=1)
# dev.off()

cor.test(SummaryCVSowingTimes$Sowing1, SummaryCVSowingTimes$Sowing3)

# png("./figures/new_parts_Feb2020/CVMagicParents 30 days vs 45 days.png", width = 10, height =10, units = "cm", res = 600)
SummaryCVSowingTimes%>%
  ggplot(., aes(x = Sowing1, y=Sowing2)) +
  geom_point(size=3.5) +
  scale_colour_manual(values = rainbow(25)) +
  xlab("~30 days") + ylab("~45 days") +
  theme (text = element_text(size=12))+
  theme(legend.position = "none") +
  ylim(0, 0.35) +
  xlim (0, 0.35)+
  geom_abline(a=0, b=1)
# dev.off()

cor.test(SummaryCVSowingTimes$Sowing1, SummaryCVSowingTimes$Sowing2)




##Fig1 Figure supplement 1 C-F
# plot with one row for each MAGIC line, with the size and colour of the dot reflecting the percent germinaiton on each day
#png("./plots/M101DistributionsNoLegend.png", width = 10, height = 10, units = "cm", res = 600)
germ_long %>%
  filter(Magic_line=="M101" & MeanDAR > 30 & MeanDAR< 34) %>%
  ggplot(aes(x = day, y = name, size=ifelse(perc==0, NA, perc), colour =as.factor(Mother.plant.sowing.number))) + geom_point(stat = "identity") +
  xlim(0, 60)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Days to germination")+
  ylab("Replicate")+
  theme(legend.position="none")+
  ggtitle("M101 Whole plants")
#dev.off()

#png("./plots/M174DistributionsNoLegend.png", width = 10, height = 10, units = "cm", res = 600)
germ_long %>%
  filter(Magic_line=="M174") %>%
  ggplot(aes(x = day, y = name, size=ifelse(perc==0, NA, perc), colour =as.factor(Mother.plant.sowing.number))) + geom_point(stat = "identity") +
  xlim(0, 30)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Days to germination")+
  ylab("Replicate")+
  theme(legend.position="none")+
  ggtitle("M174 Whole plants")
#dev.off()

#png("./plots/M182DistributionsNoLegendBlack.png", width = 10, height = 10, units = "cm", res = 600)
germ_long %>%
  filter(Magic_line=="M182" & MeanDAR > 27 & MeanDAR < 34) %>%
  ggplot(aes(x = day, y = name, size=ifelse(perc==0, NA, perc), colour =as.factor(Mother.plant.sowing.number))) + geom_point(stat = "identity") +
  xlim(0, 60)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Days to germination")+
  ylab("Replicate")+
  theme(legend.position="none")+
  ggtitle("M182 Whole plants")
#dev.off()

#png("./plots/M178DistributionsNoLegend.png", width = 10, height = 10, units = "cm", res = 600)
germ_long %>%
  filter(Magic_line=="M178" & MeanDAR > 24 & MeanDAR < 36) %>%
  ggplot(aes(x = day, y = name, size=ifelse(perc==0, NA, perc), colour =as.factor(Mother.plant.sowing.number))) + geom_point(stat = "identity") +
  xlim(0, 60)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Days to germination")+
  ylab("Replicate")+
  theme(legend.position="none")+
  ggtitle("M178 Whole plants")
#dev.off()

##to output source data as one table for each of these lines, will filter germ long as it has been filtered for these graphs above

Fig1FigSup1CtoFSourceData<-germ_long %>%  filter((Magic_line=="M101" & MeanDAR > 30 & MeanDAR< 34)| Magic_line=="M174"| (Magic_line=="M182" & MeanDAR > 27 & MeanDAR < 34) | (Magic_line=="M178" & MeanDAR > 24 & MeanDAR < 36) )

# write this table as source data for  Figure1_Figure_Supplement1C-F
# write.table(Fig1FigSup1CtoFSourceData, './source_data/Figure1_FigureSupplement1_source_data2.tsv', quote = F, row.names = F, sep = "\t")


### plots for figure1 figure supplement 1G

# png("./SoilGermExp1.png", width = 12, height = 9, units = "cm", res = 600)
germ_longSoilvsPlates%>%
  mutate(Magic_line=factor(Magic_line, levels = c("194", "203", "311", "182", "178")))%>%
  filter(experiment=="1")%>%
  ggplot(., aes(x = day, y = Magic_line, size=ifelse(perc==0, NA, perc), colour = High_or_low_var)) + geom_point(stat="identity")+
  xlab("Days to germination") + ylab("Magic")+
  theme (text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none")
# dev.off()

# png("./PlotWithLegend.png", width = 12, height = 9, units = "cm", res = 600)
germ_longSoilvsPlates%>%
  mutate(Magic_line=factor(Magic_line, levels = c("194", "203", "311", "182", "178")))%>%
  filter(experiment=="1")%>%
  ggplot(., aes(x = day, y = Magic_line, size=ifelse(perc==0, NA, perc), colour = High_or_low_var)) + geom_point(stat="identity")+
  xlab("Days to germination") + ylab("Magic")+
  theme (text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# dev.off()


# png("./SoilGermExp2Col-0No.png", width = 12, height = 4, units = "cm", res = 600)
germ_longSoilvsPlates%>%
  filter(experiment=="2")%>%
  ggplot(., aes(x = day, y = Magic_line, size=ifelse(perc==0, NA, perc), colour = High_or_low_var)) + geom_point(stat="identity")+
  xlab("Days to germination") + ylab("Magic")+
  theme (text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none")
# dev.off()

# write the table as source data for  Figure1_Figure_Supplement1 G
# write.table(germ_longSoilvsPlates, './source_data/Figure1_FigureSupplement1_source_data3.tsv', quote = F, row.names = F, sep = "\t")


#
# Data for comparisons with single silique data ----
#
# Generate plots used  to compare whole plant germination time distributions with
# single silique and half silique - here just doing whole plant distributions
# (used for Fig. 2 and Fig 2 supplement 1).
# Drawing points whose size represents percent germination on each day.
# The half silique data is processed using a different script (02_SingleandHalfSiliques.R)

#png("./plots/PNGs/SingleSiliques/M178DistributionsNoLegendBlack.png", width = 9, height = 9, units = "cm", res = 600)
germ_long %>%
  filter(Magic_line=="M178" & MeanDAR > 24 & MeanDAR < 36) %>%
  ggplot(aes(x = day, y = name, size = ifelse(perc==0, NA, perc))) +
  geom_point() +
  xlim(0, 60)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Days to germination")+
  ylab("Replicate")+
  theme(legend.position="none")+
  ggtitle("M178 Whole plants")
#dev.off()

# output source data for Figure 2A
germ_longM178<-germ_long %>%
  filter(Magic_line=="M178" & MeanDAR > 24 & MeanDAR < 36)

# write.table(germ_longM178, './source_data/Figure2_source_data1.tsv', quote = F, row.names = F, sep = "\t")


#png("./plots/PNGs/SingleSiliques/M182DistributionsNoLegendBlack.png", width = 10, height = 10, units = "cm", res = 600)
germ_long %>%
  filter(Magic_line=="M182" & MeanDAR > 27 & MeanDAR < 34) %>%
  ggplot(aes(x = day, y = name, size = ifelse(perc==0, NA, perc))) +
  geom_point() +
  xlim(0, 60)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Days to germination")+
  ylab("Replicate")+
  theme(legend.position="none")+
  ggtitle("M182 Whole plants")
#dev.off()


#png("./plots/PNGs/SingleSiliques/M53DistributionsNoLegend.png", width = 10, height = 10, units = "cm", res = 600)
germ_long %>%
  filter(Magic_line=="M53" & MeanDAR ==27) %>%
  ggplot(aes(x = day, y = name, size=ifelse(perc==0, NA, perc))) +
  geom_point() +
  xlim(0, 60)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Days to germination")+
  ylab("Replicate")+
  theme(legend.position="none")+
  ggtitle("M53 Whole plants")
#dev.off()


#png("./plots/PNGs/SingleSiliques/M4DistributionsNoLegend.png", width = 10, height = 10, units = "cm", res = 600)
germ_long %>%
  filter(Magic_line=="M4") %>%
  ggplot(aes(x = day, y = name, size = ifelse(perc==0, NA, perc))) +
  geom_point() +
  xlim(0, 60)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Days to germination")+
  ylab("Replicate")+
  theme(legend.position="none")+
  ggtitle("M4 Whole plants")
#dev.off()


# output source data for Figure 2 figure supplement 1, M182, M53 and M4 whole plants.
germ_longFigure2FigSup1<-germ_long%>%filter((Magic_line=="M182" & MeanDAR > 27 & MeanDAR < 34)| (Magic_line=="M53" & MeanDAR ==27) |(Magic_line=="M4"))

# write.table(germ_longFigure2FigSup1, './source_data/Figure2_Figure_supplement1_source_data1.tsv', quote = F, row.names = F, sep = "\t")


#
# Mode vs CV ----
#
# Scatterplot of mode vs cv coloured by percent germ using mean values across
# replicate samples or experiments. Plotting all lines for Fig 3

# outputting source data for this section
# write.table(germ_summaryPerLine, './source_data/Figure3_source_data1.tsv', quote = F, row.names = F, sep = "\t")


#png("./plots/ModevsCVColouredPercentNolegendBLUE.png", width = 8.9, height = 7.6, units = "cm", res = 600)
ggplot(germ_summaryPerLine, aes(mean_mode, mean_cv, colour = mean_percent)) +
  geom_point(size = 2) +
  xlim(0, 8) + ylim(0, 1.75) +
  theme(text = element_text(size=12))+
  xlab("Mode days to germination")+
  ylab("CV of germination time") +
  theme(legend.position="none")+
  stat_smooth(method = "lm", se = FALSE, colour="black") +
  scale_colour_gradient2(low = "#ece7f2", mid = "#a6bddb", high = "#2b8cbe", midpoint = 50)
#dev.off()


#plotting only those with lower CV
#png("./plots/ModevsCVColouredPercentWithoutExtremesBLUEWithLegend.png", width = 8.9, height = 7.6, units = "cm", res = 600)
germ_summaryPerLine %>%
  filter (mean_cv < 0.58) %>%
  ggplot(aes(mean_mode, mean_cv, colour = mean_percent)) +
  geom_point(size = 2) +
  xlab("Mode days to germination")+
  ylab("CV of germination time") +
  xlim(0, 8) + ylim(0, 0.6)+
  theme(text = element_text(size=12)) +
  stat_smooth(method = "lm", se = FALSE, colour="black") +
  scale_colour_gradient2(low = "#ece7f2", mid = "#a6bddb", high = "#2b8cbe", midpoint = 50)
#dev.off()


##Now doing correlations between CV and percent for Supplmentary Fig 3
#png("./Plots/PNGs/ModevsCV/PercentvsCV.png", width = 8.9, height = 7.6, units = "cm", res = 600)
ggplot(germ_summaryPerLine, aes(mean_percent, mean_cv)) +
  geom_point(size = 2) +
  xlim(0, 100) + ylim(0, 1.75) +
  theme(text = element_text(size=12))+
  xlab("Percentage germination")+
  ylab("CV of germination time") +
  theme(legend.position="none") +
  stat_smooth(method = "lm", se = FALSE, colour="blue")
#dev.off()



#plotting only those with lower CV
#png("./Plots/PNGs/ModevsCV/PercentvsCVZoomNoTrend.png", width = 8.9, height = 7.6, units = "cm", res = 600)
germ_summaryPerLine %>%
  filter (mean_cv < 0.58) %>%
  ggplot(aes(mean_percent, mean_cv)) +
  geom_point(size = 2) +
  xlim(0, 100) + ylim(0, 0.6) +
  theme(text = element_text(size=12))+
  xlab("Percentage germination")+
  ylab("CV of germination time") +
  theme(legend.position="none")+
  stat_smooth(method = "lm", se = FALSE, colour="blue")
#dev.off()


##Now doing correlations between percent and mode (not shown in paper)
#png("./Plots/PNGs/ModevsCV/PercentvsMode.png", width = 8.9, height = 7.6, units = "cm", res = 600)
ggplot(germ_summaryPerLine, aes(mean_mode, mean_percent)) +
  geom_point(size = 2) +
  xlim(0, 10) + ylim(0, 100) +
  theme(text = element_text(size=12))+
  xlab("Mode days to germination")+
  ylab("Percentage germination") +
  theme(legend.position="none")+
  stat_smooth(method = "lm", se = FALSE, colour="blue")
#dev.off()

#png("./Plots/PNGs/ModevsCV/PercentvsModeNoTrend.png", width = 8.9, height = 7.6, units = "cm", res = 600)
ggplot(germ_summaryPerLine, aes(mean_mode, mean_percent)) +
  geom_point(size = 2) +
  xlim(1, 7.5) + ylim(0, 100) +
  theme(text = element_text(size=12), legend.position="none")+
  xlab("Mode days to germination")+
  ylab("Percentage germination")
#dev.off()



#doing Pearson's correlation between the different traits.
cor.test(germ_summaryPerLine$mean_mode, germ_summaryPerLine$mean_cv)

cor.test(germ_summaryPerLine$mean_percent, germ_summaryPerLine$mean_cv)

cor.test(germ_summaryPerLine$mean_mode, germ_summaryPerLine$mean_percent)

germSummaryPerLineUnimodal <- germ_summaryPerLine%>%filter(mean_cv <0.58)
cor.test(germSummaryPerLineUnimodal$mean_mode, germSummaryPerLineUnimodal$mean_cv)



###plotting pairs of lines' germination distributions for Fig3 and Fig 3 Fig supplement1


#comparing lines with different modes and different CVs
#png("./plots/M123andM141.png", width = 8.9, height = 7.6, units = "cm", res = 600)
germ_long %>%
  filter((Magic_line == "M123" & MeanDAR == 30) | (Magic_line == "M141" & MeanDAR == 33)) %>%
  ggplot(., aes(x = day, y = perc, fill=as.factor(RepNumber))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Magic_line ~.)+
  scale_fill_grey(start = 0.7, end = 0.3)+
  xlab("Days to germination") + ylab("Percentage of seeds")+
  scale_x_continuous(breaks=seq(0, 16,2), limits = c(0,  16)) +
  theme(text = element_text(size=12)) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#dev.off()

#comparing lines with mode of around 4
#png("./plots/M178andM192.png", width = 8.9, height = 7.6, units = "cm", res = 600)
germ_long %>%
  filter((Magic_line == "M192" ) |(Magic_line == "M178" & MeanDAR == 34) ) %>%
  ggplot(., aes(x = day, y = perc, fill= as.factor(RepNumber))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Magic_line ~.)+
  scale_fill_grey(start = 0.7, end = 0.3)+
  xlab("Days to germination") + ylab("Percentage of seeds")+
  theme (text = element_text(size=12))+
  scale_x_continuous(breaks=seq(0, 55, 5), limits = c(0,  55))+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#dev.off()


#comparing lines with mode of 2, M123 and M355

#png("./plots/M123andM355.png", width = 8.9, height = 7.6, units = "cm", res = 600)
germ_long %>%
  filter((Magic_line == "M123" ) | (Magic_line == "M355") ) %>%
  ggplot(., aes(x = day, y = perc, fill= as.factor(RepNumber))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Magic_line ~.)+
  scale_fill_grey(start = 0.7, end = 0.3)+
  xlab("Days to germination") + ylab("Percentage of seeds")+
  theme (text = element_text(size=12))+
  scale_x_continuous(breaks=seq(0, 15, 2), limits = c(0,  15))+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#dev.off()

#comparing lines with mode of 3, low cv and mode of 4, high cv

#png("./M108andM492.png", width = 8.9, height = 7.6, units = "cm", res = 600)
germ_long %>%
  filter((Magic_line == "M108" ) | (Magic_line == "M492") ) %>%
  ggplot(., aes(x = day, y = perc, fill= as.factor(RepNumber))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Magic_line ~.)+
  scale_fill_grey(start = 0.7, end = 0.3)+
  xlab("Days to germination") + ylab("Percentage of seeds")+
  theme (text = element_text(size=12))+
  scale_x_continuous(breaks=seq(0, 15, 2), limits = c(0,  15))+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#dev.off()

##comparing lines with similar modes and % germ but different cvs for Figure 3 figure supplement 1

#png("./426andM311.png", width = 8.9, height = 7.6, units = "cm", res = 600)
germ_long %>%
  filter((Magic_line == "M426" ) | (Magic_line == "M311" & MeanDAR==31) ) %>%
  ggplot(., aes(x = day, y = perc, fill= as.factor(RepNumber))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Magic_line ~.)+
  scale_fill_grey(start = 0.7, end = 0.3)+
  xlab("Days to germination") + ylab("Percentage of seeds")+
  theme (text = element_text(size=12))+
  scale_x_continuous(breaks=seq(0, 75, 5), limits = c(0,  75))+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#dev.off()

# png("./108andM393.png", width = 8.9, height = 7.6, units = "cm", res = 600)
germ_long %>%
  filter((Magic_line == "M108" ) | (Magic_line == "M393") ) %>%
  ggplot(., aes(x = day, y = perc, fill= as.factor(RepNumber))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Magic_line ~.)+
  scale_fill_grey(start = 0.7, end = 0.3)+
  xlab("Days to germination") + ylab("Percentage of seeds")+
  theme (text = element_text(size=12))+
  scale_x_continuous(breaks=seq(0, 20, 2), limits = c(0,  20))+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#dev.off()

##output source data with just the MAGIC lines used for these plots (For Fig 3 B-E and Fig 3, Fig supplement 1, B & C).
germLongPairComparison<- germ_long%>%filter((Magic_line == "M123" & MeanDAR == 30) | (Magic_line == "M141" & MeanDAR == 33) | (Magic_line == "M192" ) |(Magic_line == "M178" & MeanDAR == 34) | (Magic_line == "M123" ) | (Magic_line == "M355") | (Magic_line == "M108" ) | (Magic_line == "M492")  |(Magic_line == "M426" ) | (Magic_line == "M311" & MeanDAR==31) | (Magic_line == "M393") )

# write.table(germLongPairComparison, './source_data/Figure3_source_data2.tsv', quote = F, row.names = F, sep = "\t")



##doing scatter plots of mode vs cv and percent vs CV just for accessions for Figure 3, figure supplement 2

##output source data for this
# write.table(germ_summary_accessionsOnly, './source_data/Figure3_FigureSupplement2_source_data1.tsv', quote = F, row.names = F, sep = "\t")

#png("./Plots/PNGs/ModevsCV/ModevsCVColouredPercentNolegendBLUEAccessionsNoLine.png", width = 8.9, height = 7.6, units = "cm", res = 600)
ggplot(germ_summary_accessionsOnly, aes(mean_mode, mean_cv, colour = mean_percent)) +
  geom_point(size = 2) +
  xlim(0, 7.5) + ylim(0, 0.75)+   theme(text = element_text(size=12))+
  xlab("Mode days to germination")+
  ylab("CV of germination time") +
  theme(legend.position="none")+
  scale_colour_gradient2(low = "#ece7f2", mid = "#a6bddb", high = "#2b8cbe", midpoint = 50) +
  stat_smooth(method = "lm", se = FALSE, colour="black")
#dev.off()

#png("./Plots/PNGs/ModevsCV/PercentvsCVAccessionsNoLine.png", width = 8.9, height = 7.6, units = "cm", res = 600)
ggplot(germ_summary_accessionsOnly, aes(mean_percent, mean_cv)) +
  geom_point(size = 2) +
  xlim(0, 100) + ylim(0, 0.75) +
  theme(text = element_text(size=12))+
  xlab("Percentage germination")+
  ylab("CV of germination time") +
  theme(legend.position="none")+
  stat_smooth(method = "lm", se = FALSE, colour="blue")
#dev.off()

##making a long table for Ct-1 and Mad-0
GermLongCt <- germ_longParents %>%
  filter(Magic_Parent_line=="Ct-1" & Mean_DAR==28) %>%
  rename(Line = Magic_Parent_line) %>%
  select(Line, day, perc, name, group, RepNumber)

# Spannish accessions
GermLongMad <- germ_longAccessions %>%
  filter((Accession =="Mad-0" & MeanDAR ==27)) %>%
  rename( Line = Accession) %>%
  select(Line, day, perc, name, group,  RepNumber)

# Combine both tables
ForPlottingCombinedCtMad <- bind_rows(GermLongCt , GermLongMad)

##output this table as source data for Figure 3 figure supplement 2
# write.table(ForPlottingCombinedCtMad, './source_data/Figure3_FigureSupplement2_source_data2.tsv', quote = F, row.names = F, sep = "\t")



##comparing distributions of Ct-1 and Mad-0 for Supplmentary Fig 4c
#png("./Plots/PNGs/ModevsCV/CtvsMad.png", width = 8.9, height = 7.6, units = "cm", res = 600)
ForPlottingCombinedCtMad %>%
  mutate(Line = factor(Line, levels = c("Ct-1", "Mad-0")),
         RepNumber = factor(RepNumber)) %>%
  ggplot(aes(x = day, y = perc, fill= RepNumber)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Line ~.)+
  scale_fill_grey(start = 0.7, end = 0.3)+
  xlab("Days to germination") + ylab("Percentage of seeds")+
  theme (text = element_text(size=12))+
  scale_x_continuous(breaks=seq(0, 10, 1), limits = c(0,  10))+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#dev.off()


##testing Pearson's correlation for the traits.

cor.test(germ_summary_accessionsOnly$mean_mode, germ_summary_accessionsOnly$mean_cv)

cor.test(germ_summary_accessionsOnly$mean_percent, germ_summary_accessionsOnly$mean_cv)



