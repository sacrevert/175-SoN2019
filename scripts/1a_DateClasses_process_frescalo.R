## O.L. Pescott
## started 19 Dec 2018
## Script: 1a_DateClasses_process_frescalo.R

## Script to process and run Frescalo on BSBI hectad data (vcs 1-112)
## using adjusted MOH taxa (v4) and BSBI date classes from 1970 - 2018
## ultimately for plant trends for State of Nature 2019 (Daniel Hayhow is RSPB link)

#rm(list=ls())
library(plyr)
library(reshape2)
library(lubridate)

#### Read and process data ####
## update (csv file v2) is an updated BSBI extract excluding the unmapped occurrences
#datDC <- read.csv(file = "data/bsbi/dateClassesDat_1970_2018.csv", header = T, stringsAsFactors = F)
#save(datDC, file = "data/bsbi/dateClassesDat_1970_2018.Rdata")
#load(file = "data/bsbi/dateClassesDat_1970_2018.Rdata")

#datDC <- read.csv(file = "data/bsbi/dateClassesDat_1970_2018_v2.csv", header = T, stringsAsFactors = F)
#save(datDC, file = "data/bsbi/dateClassesDat_1970_2018_v2.Rdata")
load(file = "data/bsbi/dateClassesDat_1970_2018_v2.Rdata")

head(datDC)
names(datDC) <- c("taxon", "qualifier", "hectad", "1970", "1987", "2000", "2010", "freq")
datDC <- datDC[,-8] # remove summary column
# re-concat name and qualifier
datDC$taxon <- ifelse(datDC$qualifier == "", datDC$taxon, paste(datDC$taxon, datDC$qualifier))
datDC <- datDC[,-2]
head(datDC)
# remove blank hectads
datDC <- datDC[!(datDC$hectad==""),]
# replace all record counts with presence/absence
str(datDC)
presAbs <- function(x){
  ifelse(x > 1, 1, x)
}
datDC[, 3:6] <- apply(datDC[,3:6], 2, FUN = presAbs)
head(datDC) # need to make data info into long-form for Frescalo

LdatDC <- melt(datDC, measure.vars = c("1970", "1987", "2000", "2010"))
LdatDC <- LdatDC[LdatDC$value==1,] # delete absences
# delete erroneous taxon (P. aviculare sensu Sell and Murrell) -- BSBI ddb bug?
LdatDC <- LdatDC[!(LdatDC$taxon == "Polygonum aviculare sensu Sell & Murrell"),] # delete taxon (only two rows)

names(LdatDC)[3] <- "year"
LdatDC$year <- as.Date(LdatDC$year, "%Y")
LdatDC$year <- lubridate::year(LdatDC$year)
rownames(LdatDC) <- 1:nrow(LdatDC)
head(LdatDC)

#### Run Frescalo on date-class data ####
library(sparta)
sinkdir <- "outputs/fres_out_DC" #frescalo folder will be dated 24 12 2018 for BSBI Extract v2
myFresPath <- 'C:\\Frescalo_3a_windows.exe'

fres_out_DC <- frescalo(Data = LdatDC,
                     time_periods = data.frame(start=c(1970,1987,2000,2010),end=c(1986,1999,2009,2018)),
                     sinkdir = sinkdir,
                     Fres_weights = 'LCGB',
                     site_col = 'hectad',
                     sp_col = 'taxon',
                     year_col = 'year',
                     frespath = myFresPath, phi = 0.87, ## Your value of phi (0.74) is smaller than the 98.5 percentile of input phi (0.87). 
                     alpha = 0.27)
