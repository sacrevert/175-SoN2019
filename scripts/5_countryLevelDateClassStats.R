## O.L. Pescott
## started 24 Dec 2018
## Script: 5_countryLevelDateClassStats.R

## Split BSBI data by country and run multiple Frescalo runs
library(reshape2)

countryHecs <- read.csv(file = "data/hectads/hec_mon_vc_country.csv", header = T, stringsAsFactors = F)
countryHecs <- countryHecs[countryHecs$CN_NAME %in% c("ENGLAND","WALES","SCOTLAND","NORTHERN IRELAND"),]

#### get and process BSBI DC data again ####
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

#### Subset data into countries based on hectads ####
engHecs <- unique((countryHecs[countryHecs$CN_NAME=="ENGLAND",])[,c("CN_NAME", "hectad")])
nirHecs <- unique((countryHecs[countryHecs$CN_NAME=="NORTHERN IRELAND",])[,c("CN_NAME", "hectad")]) ## no data currently
scoHecs <- unique((countryHecs[countryHecs$CN_NAME=="SCOTLAND",])[,c("CN_NAME", "hectad")])
walHecs <- unique((countryHecs[countryHecs$CN_NAME=="WALES",])[,c("CN_NAME", "hectad")])

## Extract data
engDat <- merge(engHecs, LdatDC, by.x = "hectad", by.y = "hectad", all.x = T, all.y = F)
scoDat <- merge(scoHecs, LdatDC, by.x = "hectad", by.y = "hectad", all.x = T, all.y = F)
walDat <- merge(walHecs, LdatDC, by.x = "hectad", by.y = "hectad", all.x = T, all.y = F)

#### England Frescalo #####
#### Run Frescalo on date-class data ####
library(sparta)
sinkdir <- "outputs/fres_out_DC_Eng" #frescalo folder will be dated 24 12 2018 for BSBI Extract v2
myFresPath <- 'C:\\Frescalo_3a_windows.exe'

fres_out_DC_E <- frescalo(Data = engDat,
                        time_periods = data.frame(start=c(1970,1987,2000,2010),end=c(1986,1999,2009,2018)),
                        sinkdir = sinkdir,
                        Fres_weights = 'LCGB',
                        site_col = 'hectad',
                        sp_col = 'taxon',
                        year_col = 'year',
                        frespath = myFresPath, phi = NULL, ## Your value of phi (0.74) is smaller than the 98.5 percentile of input phi (0.87). 
                        alpha = 0.27)
#####

#### Scotland Frescalo #####
#### Run Frescalo on date-class data ####
sinkdir <- "outputs/fres_out_DC_Sco" #frescalo folder will be dated 24 12 2018 for BSBI Extract v2

fres_out_DC_S <- frescalo(Data = scoDat,
                        time_periods = data.frame(start=c(1970,1987,2000,2010),end=c(1986,1999,2009,2018)),
                        sinkdir = sinkdir,
                        Fres_weights = 'LCGB',
                        site_col = 'hectad',
                        sp_col = 'taxon',
                        year_col = 'year',
                        frespath = myFresPath, phi = NULL, ## Your value of phi (0.74) is smaller than the 98.5 percentile of input phi (0.87). 
                        alpha = 0.27)

#####

#### Wales Frescalo #####
#### Run Frescalo on date-class data ####
sinkdir <- "outputs/fres_out_DC_Wal" #frescalo folder will be dated 24 12 2018 for BSBI Extract v2

fres_out_DC_W <- frescalo(Data = walDat,
                        time_periods = data.frame(start=c(1970,1987,2000,2010),end=c(1986,1999,2009,2018)),
                        sinkdir = sinkdir,
                        Fres_weights = 'LCGB',
                        site_col = 'hectad',
                        sp_col = 'taxon',
                        year_col = 'year',
                        frespath = myFresPath, phi = NULL, ## Your value of phi (0.74) is smaller than the 98.5 percentile of input phi (0.87). 
                        alpha = 0.27)

# (1) Also do above using Northern Irish data
library(reshape2)
#datDC_NI <- read.csv(file = "data/bsbi/dateClassesDat_1970_2018_v2_NI.csv", header = T, stringsAsFactors = F)
#save(datDC_NI, file = "data/bsbi/dateClassesDat_1970_2018_v2_NI.Rdata")
load(file = "data/bsbi/dateClassesDat_1970_2018_v2_NI.Rdata")

head(datDC_NI)
names(datDC_NI) <- c("taxon", "qualifier", "hectad", "1970", "1987", "2000", "2010", "freq")
datDC_NI <- datDC_NI[,-8] # remove summary column
# re-concat name and qualifier
datDC_NI$taxon <- ifelse(datDC_NI$qualifier == "", datDC_NI$taxon, paste(datDC_NI$taxon, datDC_NI$qualifier))
datDC_NI <- datDC_NI[,-2]
head(datDC_NI)
# remove blank hectads
datDC_NI <- datDC_NI[!(datDC_NI$hectad==""),]
# replace all record counts with presence/absence
str(datDC_NI)
presAbs <- function(x){
  ifelse(x > 1, 1, x)
}
datDC_NI[, 3:6] <- apply(datDC_NI[,3:6], 2, FUN = presAbs)
head(datDC_NI) # need to make data info into long-form for Frescalo

LdatDC_NI <- melt(datDC_NI, measure.vars = c("1970", "1987", "2000", "2010"))
LdatDC_NI <- LdatDC_NI[LdatDC_NI$value==1,] # delete absences
# delete erroneous taxon (P. aviculare sensu Sell and Murrell) -- BSBI ddb bug?
LdatDC_NI <- LdatDC_NI[!(LdatDC_NI$taxon == "Polygonum aviculare sensu Sell & Murrell"),] # delete taxon (only two rows)

names(LdatDC_NI)[3] <- "year"
LdatDC_NI$year <- as.Date(LdatDC_NI$year, "%Y")
LdatDC_NI$year <- lubridate::year(LdatDC_NI$year)
rownames(LdatDC_NI) <- 1:nrow(LdatDC_NI)
head(LdatDC_NI)

nirDat <- merge(nirHecs, LdatDC_NI, by.x = "hectad", by.y = "hectad", all.x = T, all.y = F) ## no data currently

library(sparta)
sinkdir <- "outputs/fres_out_DC_NI" #frescalo folder will be dated 24 12 2018 for BSBI Extract v2
myFresPath <- 'C:\\Frescalo_3a_windows.exe'

fres_out_DC_NI <- frescalo(Data = nirDat,
                          time_periods = data.frame(start=c(1970,1987,2000,2010),end=c(1986,1999,2009,2018)),
                          sinkdir = sinkdir,
                          Fres_weights = 'LCNI', # note change here
                          site_col = 'hectad',
                          sp_col = 'taxon',
                          year_col = 'year',
                          frespath = myFresPath, phi = NULL, ## Your value of phi (0.74) is smaller than the 98.5 percentile of input phi (0.87). 
                          alpha = 0.27)

##### TO DO ######
# (2) Process all country-level Frescalo outputs as per script 4.
# (3) Weighted average for UK (Britain + NI)?