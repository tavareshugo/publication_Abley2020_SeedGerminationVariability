library(dplyr) # install.packages("dplyr")
library(tidyr) # install.packages("tidyr")
library(ggplot2)
theme_set(theme_bw())

##---------------------------For Magic lines---------------------------------

germ <- read.csv('./data/MAGICandAccessions/ForSubmission/MAGICsRawData.csv') # read in MAGIC all data that has experiments with different DAR
germ$id <- 1:nrow(germ)
germ<-germ%>% group_by(Magic_line, Above_or_below_agar, MeanDAR, Mother.plant.sowing.number)%>%
  mutate(RepNumber = as.character(1: n()))

# For MAGIC data, make data long format and tidy it
germ_long <- germ %>%
	gather(day, number, X1:X87) %>%
	arrange(Magic_line, Plant, Above_or_below_agar, DAR, MeanDAR, Mother.plant.sowing.number, RepNumber) %>%
	mutate(day = as.numeric(gsub("X", "", day))) %>%
  mutate(name = paste( Magic_line,Sowing_date, Plant, sep = "_")) %>%
  filter(!is.na(number))

# Calculate proportion on each day, filter for MeanDAR between 20 and 45
germ_long <- germ_long %>%
	mutate(perc = number/Total*100, perc_total_ger = number/Total.germinated*100, group="magic") %>%
  filter(Above_or_below_agar == "B")%>%
  filter(MeanDAR > 20 & MeanDAR < 45)

##this is outputting a tsv file.
write.table(germ_long, './data/MAGICandAccessions/ForSubmission/MAGICs.tsv', quote = F, row.names = F, sep = "\t")

##From this file, I manually delete columns RepNumber and Above_or_below_agar as they don't contain any info. Also manually deleted line 345* as there was a problem with this one. This it is then saved as Supplemental_File4_MAGICs.csv. Will be filtered for number germinated <10 in the later script used for summarising and making plots.



##---------------------------For Magic parent accessions---------------------------------


# read in MAGIC parent data
germParents <- read.csv('./data/MAGICandAccessions/ForSubmission/MAGICParentsRawData.csv') # read in MAGIC parent data
germParents$id <- 1:nrow(germParents)
germParents<-germParents%>% group_by(Magic_Parent_line, Sowing_round, Sowing_date)%>%
  mutate(RepNumber = as.character(1: n()))

# For MAGIC Parent data, make data long format and tidy it
germ_longParents <- germParents %>%
  gather(day, number, X1:X14) %>%
  arrange(Magic_Parent_line, Plant, DAR, Mean_DAR, Sowing_round, Parent_vernalised, RepNumber) %>%
  mutate(day = as.numeric(gsub("X", "", day))) %>%
  mutate(name = paste(Magic_Parent_line, Sowing_date, Plant, Mean_DAR, sep = "_")) %>%
  filter(!is.na(number))


# Calculate proportion on each day
germ_longParents <- germ_longParents %>%
  mutate(perc = number/Total*100, perc_total_ger = number/Total.germinated*100, group="magic parents")

write.table(germ_longParents, './data/MAGICandAccessions/ForSubmission/MAGICParents.tsv', quote = F, row.names = F, sep = "\t")

##I edited this table to change the numbers of the MAGIC parents to the accession names. Then saved as Supplemental_File2_MAGICparents
##Will be filtered for number germinated <10 in the later script used for summarising and making plots.


##---------------------------For spanish accessions---------------------------------


# Read in the Spanish accessions data (accession that are not MAGIC parents)

germAccessions <- read.csv('./data/MAGICandAccessions/ForSubmission/SpanishAccessionsRawData.csv') # read in Accessions data
germAccessions$id <- 1:nrow(germAccessions)
germAccessions<-germAccessions%>% group_by(Accession, Sowing_date)%>%
  mutate(RepNumber = as.character(1: n()))

germ_longAccessions <- germAccessions %>%
  gather(day, number, X1:X87) %>%
  arrange(Accession, Plant, Sowing_date, MeanDAR) %>%
  mutate(day = as.numeric(gsub("X", "", day))) %>%
  mutate(name = paste( Accession, Plant, Sowing_date, MeanDAR, sep = "_")) %>%
  filter(!is.na(number))

germ_longAccessions <- germ_longAccessions %>%
  mutate(perc = number/Total*100, perc_total_ger = number/Total.germinated*100, group="accessions")


write.table(germ_longAccessions, './data/MAGICandAccessions/ForSubmission/SpanishAccessions.tsv', quote = F, row.names = F, sep = "\t")

##outputted this file, then manually added length of time of vernalisation
##Will be filtered for number germinated <10 in the later script used for summarising and making plots. This removes 5 of the 15 lines



##---------------------------For soil vs plates experiment---------------------------------

germSoil <- read.csv('./data/MAGICandAccessions/ForSubmission/soil_vs_platesRawData.csv')

# Make data long format and tidy it
germ_longSoil <- germSoil %>%
  gather(day, number, X1:X60) %>%
  arrange(Magic_line, tray) %>%
  mutate(day = as.numeric(gsub("X", "", day))) %>%
  filter(!is.na(number)) %>%
  mutate(name = paste(Magic_line, tray, sep = "_"))

# Calculate proportion on each day
germ_longSoil <- germ_longSoil %>%
  mutate(perc = number/Total_germ*100)

write.table(germ_longSoil, './data/MAGICandAccessions/ForSubmission/Soil_vs_plates.tsv', quote = F, row.names = F, sep = "\t")

