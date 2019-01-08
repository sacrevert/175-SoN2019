## O.L. Pescott
## started 7th Jan 2019.
## Script: 5b_countryLevelDateClassStats.R
## update this process based on the fact that we now want to include the 1930-1969 date class
## (because the 1970-86 date-class is quite biased as a baseline)

## Split BSBI data by country and run multiple Frescalo runs
library(reshape2)

countryHecs <- read.csv(file = "data/hectads/hec_mon_vc_country.csv", header = T, stringsAsFactors = F)
# pare down
countryHecs <- countryHecs[countryHecs$CN_NAME %in% c("ENGLAND","WALES","SCOTLAND","NORTHERN IRELAND"),]

#### get and process BSBI DC data again ####
## now starting with 1930-69 DC
#datDC2 <- read.csv(file = "data/bsbi/dateClassesDat_1930_2018.csv", header = T, stringsAsFactors = F)
#save(datDC2, file = "data/bsbi/dateClassesDat_1930_2018.Rdata")
load(file = "data/bsbi/dateClassesDat_1930_2018.Rdata")

head(datDC2)
names(datDC2) <- c("taxon", "qualifier", "hectad", "1930", "1970", "1987", "2000", "2010", "freq")
datDC2 <- datDC2[,-9] # remove summary column
# re-concat name and qualifier
datDC2$taxon <- ifelse(datDC2$qualifier == "", datDC2$taxon, paste(datDC2$taxon, datDC2$qualifier))
datDC2 <- datDC2[,-2]
head(datDC2)
# remove blank hectads
datDC2 <- datDC2[!(datDC2$hectad==""),]
# replace all record counts with presence/absence
str(datDC2)
presAbs <- function(x){
  ifelse(x > 1, 1, x)
}
datDC2[, 3:7] <- apply(datDC2[,3:7], 2, FUN = presAbs)
head(datDC2) # need to make data info into long-form for Frescalo

LdatDC <- melt(datDC2, measure.vars = c("1930", "1970", "1987", "2000", "2010"))
LdatDC <- LdatDC[LdatDC$value==1,] # delete absences
# delete erroneous taxon (P. aviculare sensu Sell and Murrell) -- BSBI ddb bug?
LdatDC <- LdatDC[!(LdatDC$taxon == "Polygonum aviculare sensu Sell & Murrell"),] # delete taxon (only two rows)
head(LdatDC)
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
britHecs <- data.frame(hectad = unique(rbind(engHecs, scoHecs, walHecs)[,c("hectad")]))


## Extract data
engDat <- merge(engHecs, LdatDC, by.x = "hectad", by.y = "hectad", all.x = T, all.y = F)
scoDat <- merge(scoHecs, LdatDC, by.x = "hectad", by.y = "hectad", all.x = T, all.y = F)
walDat <- merge(walHecs, LdatDC, by.x = "hectad", by.y = "hectad", all.x = T, all.y = F)
nirDat <- merge(nirHecs, LdatDC, by.x = "hectad", by.y = "hectad", all.x = T, all.y = F)
britDat <- merge(britHecs, LdatDC, by.x = "hectad", by.y = "hectad", all.x = T, all.y = F)

### Do Britain first
library(sparta); myFresPath <- 'C:\\Frescalo_3a_windows.exe'
sinkdir <- "outputs/fres_out_DC1930_Brit"

fres_out_DC30_B <- frescalo(Data = britDat,
                          time_periods = data.frame(start=c(1930,1970,1987,2000,2010),end=c(1969,1986,1999,2009,2018)),
                          sinkdir = sinkdir,
                          Fres_weights = 'LCGB',
                          site_col = 'hectad',
                          sp_col = 'taxon',
                          year_col = 'year',
                          frespath = myFresPath, phi = NULL,
                          alpha = 0.27)

#### England Frescalo #####
#### Run Frescalo on date-class data ####
sinkdir <- "outputs/fres_out_DC1930_Eng"

fres_out_DC30_E <- frescalo(Data = engDat,
                        time_periods = data.frame(start=c(1930,1970,1987,2000,2010),end=c(1969,1986,1999,2009,2018)),
                        sinkdir = sinkdir,
                        Fres_weights = 'LCGB',
                        site_col = 'hectad',
                        sp_col = 'taxon',
                        year_col = 'year',
                        frespath = myFresPath, phi = NULL,
                        alpha = 0.27)

#### Scotland Frescalo #####
#### Run Frescalo on date-class data ####
sinkdir <- "outputs/fres_out_DC1930_Sco" 

fres_out_DC30_S <- frescalo(Data = scoDat,
                        time_periods = data.frame(start=c(1930,1970,1987,2000,2010),end=c(1969,1986,1999,2009,2018)),
                        sinkdir = sinkdir,
                        Fres_weights = 'LCGB',
                        site_col = 'hectad',
                        sp_col = 'taxon',
                        year_col = 'year',
                        frespath = myFresPath, phi = NULL, ## Your value of phi (0.74) is smaller than the 98.5 percentile of input phi (0.87). 
                        alpha = 0.27)

#### Wales Frescalo #####
#### Run Frescalo on date-class data ####
sinkdir <- "outputs/fres_out_DC1930_Wal" 

fres_out_DC30_W <- frescalo(Data = walDat,
                        time_periods = data.frame(start=c(1930,1970,1987,2000,2010),end=c(1969,1986,1999,2009,2018)),
                        sinkdir = sinkdir,
                        Fres_weights = 'LCGB',
                        site_col = 'hectad',
                        sp_col = 'taxon',
                        year_col = 'year',
                        frespath = myFresPath, phi = NULL,
                        alpha = 0.27)

#### Northern Ireland Frescalo #####
#### Run Frescalo on date-class data ####
sinkdir <- "outputs/fres_out_DC1930_NI" 

fres_out_DC30_NI <- frescalo(Data = nirDat,
                          time_periods = data.frame(start=c(1930,1970,1987,2000,2010),end=c(1969,1986,1999,2009,2018)),
                          sinkdir = sinkdir,
                          Fres_weights = 'LCNI', # note change here for NI land cover
                          site_col = 'hectad',
                          sp_col = 'taxon',
                          year_col = 'year',
                          frespath = myFresPath, phi = NULL,
                          alpha = 0.27)

##### Next tasks ######
# (2) Process all country-level Frescalo outputs as per script 4b.
# (3) Weighted average for UK (Britain + NI)?