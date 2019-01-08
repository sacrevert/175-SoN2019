## O.L. Pescott
## started 19 Dec 2018
## Script: 1b_FiveYear_process_frescalo.R

## Script to process and run Frescalo on BSBI hectad data (vcs 1-112)
## using adjusted MOH taxa (v4) and five-year windows from 1970 - 2018
## ultimately for plant trends for State of Nature 2019 (Daniel Hayhow is RSPB link)

#rm(list=ls())
library(plyr)
library(reshape2)
library(lubridate)

#### Read and process data ####
#datFY <- read.csv(file = "data/bsbi/fiveYear_1970_2018_v2.csv", header = T, stringsAsFactors = F) # v2 includes correct assigned vaguely dated records and excludes unmapped occurrences
#save(datFY, file = "data/bsbi/fiveYear_1970_2018_v2.Rdata")
load(file = "data/bsbi/fiveYear_1970_2018_v2.Rdata")
head(datFY)
names(datFY) <- c("taxon", "qualifier", "hectad", "1970","1975","1980","1985","1990","1995","2000","2005","2010","2015", "freq")
datFY <- datFY[,-(14:15)] # remove summary column and 2019+ column (side effect of custom date grouping on BSBI DDb)
# re-concat name and qualifier
datFY$taxon <- ifelse(datFY$qualifier == "", datFY$taxon, paste(datFY$taxon, datFY$qualifier))
datFY <- datFY[,-2]
head(datFY)
# remove blank hectads
datFY <- datFY[!(datFY$hectad==""),]
# replace all record counts with presence/absence
str(datFY)
presAbs <- function(x){
  ifelse(x > 1, 1, x)
}
datFY[, 3:12] <- apply(datFY[,3:12], 2, FUN = presAbs)
head(datFY) # need to make data info into long-form for Frescalo

LdatFY <- melt(datFY, measure.vars = c("1970","1975","1980","1985","1990","1995","2000","2005","2010","2015"))
LdatFY <- LdatFY[LdatFY$value==1,] # delete absences
names(LdatFY)[3] <- "year"
LdatFY$year <- as.Date(LdatFY$year, "%Y")
LdatFY$year <- lubridate::year(LdatFY$year)
rownames(LdatFY) <- 1:nrow(LdatFY)
head(LdatFY)

#### Run Frescalo on five-year data ####
library(sparta)
sinkdir <- "outputs/fres_out_FY"
myFresPath <- 'C:\\Frescalo_3a_windows.exe'

fres_out_FY <- frescalo(Data = LdatFY,
                        time_periods = data.frame(start=c(1970,1975,1980,1985,1990,1995,2000,2005,2010,2015),
                                                    end=c(1974,1979,1984,1989,1994,1999,2004,2009,2014,2018)),
                        sinkdir = sinkdir,
                        Fres_weights = 'LCGB',
                        site_col = 'hectad',
                        sp_col = 'taxon',
                        year_col = 'year',
                        frespath = myFresPath, phi = 0.87,
                        alpha = 0.27)