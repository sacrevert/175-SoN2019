## O.L. Pescott
## started 19 Dec 2018
## Script: 2a_DateClasses_process_multipleGAM.R

## Script to fit multiple GAMs to Frescalo BSBI date-class stats to capture
## trend uncertainty for plotting and trend fitting.
## "trend uncertainty" here is only in the form of standard errors for the GAMs (which are simulation dependent)

#rm(list=ls())
library(ggplot2)
library(mgcv)

#### Helper functions ####
#iterate function across rows of dataframe
f_lapply_row <- function(df) {
  lapply(seq_len(nrow(df)), function(i) as.list(df[i,,drop=FALSE]))
} ## end function


##### Read in data and create simulated tfactors for GAM ####
#tfDC <- read.csv(file = "outputs/fres_out_DC/frescalo_181219/Output/Trend.csv", header = T, stringsAsFactors = F)
tfDC <- read.csv(file = "outputs/fres_out_DC/frescalo_181224/Output/Trend.csv", header = T, stringsAsFactors = F) # v2 BSBI data, exclude unmapped occurrences
head(tfDC)
tfDC <- tfDC[order(tfDC$Species, tfDC$Time),]
rownames(tfDC) <- 1:nrow(tfDC)

# empty list for species data to go in
spLs <- spLsSim <- vector("list", length(unique(tfDC$Species)))
#names(spLsSim) <- unique(tfDC$Species)

# create list of species specific Frescalo trend data
spLs <- lapply(unique(tfDC$Species), function(x) tfDC[tfDC$Species == x,])
names(spLs) <- unique(tfDC$Species) # add names to list elements

# function designed to be applied to a single set of species trend info
# e.g. test <- tfDC[tfDC$Species == "Acer campestre",]
nsites <- nrow(read.csv(file = "outputs/fres_out_DC/frescalo_181224/Output/Stats.csv", header = T, stringsAsFactors = F))
simRnormSp <- function(x) { temp <- f_lapply_row(x) # convert data frame to list of lists
                            # sim new data by time period mean and sd
                            temp2 <- lapply(temp, function(x) rnorm(nsites, mean = x$TFactor, sd = x$StDev))
                            temp3 <- data.frame(do.call(cbind, temp2))
                            colnames(temp3) <- c("1978", "1993", "2004.5", "2014")
                            temp4 <- reshape2::melt(temp3, measure.vars = c("1978", "1993", "2004.5", "2014"))
                            temp4$Var1 <- temp[[1]][1]
                            temp4 <- tidyr::unnest(temp4,Var1)
                            names(temp4) <- c("year", "tfactor", "taxon")
                            temp4$year <- lubridate::year(as.Date(temp4$year, "%Y"))
                            return(temp4)
}
# run function, collecting sim data for every taxon
spLsSimDC_v2 <- lapply(spLs, FUN = simRnormSp)
save(spLsSimDC_v2, file = "outputs/sims/dateclassSims_v2.Rdata")

##### Run GAMs and plot results ####
# empty list for GAM results
#spGAMs <- vector("list", length(unique(tfDC$Species)))


# apply GAM over species, 3 knots
#spGAMs <- lapply(spLsSim, function(x) mgcv::gam(x$tfactor ~ s(x$year, k = 3), method = "REML"))

## plotting function
#spGamPlots_DC <- vector("list", length(unique(tfDC$Species)))
# produce gam plot for every species
#for (x in seq_along(spGAMs)){
#  spGamPlots_DC[[x]] <- plot(spGAMs[[x]], shade = T, rug = T, residuals = T, ylab = c("Tfactor"), xlab = c("Year"), 
#                          main = names(spGAMs)[x], seWithMean = T)
#}

## Ultimately I would like 5 year and DC plots side by side, so, do that in another script (#3) when I have an equivalent list for the 5 year GAMs


# see script 3 for continuation: 3_plotGraphsTogether.R
