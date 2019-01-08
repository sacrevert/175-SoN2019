## O.L. Pescott
## started 19 Dec 2018
## Script: 2b_FiveYear_process_multipleGAM.R

## Script to fit multiple GAMs to Frescalo BSBI five-year stats to capture
## trend uncertainty for plotting.

#rm(list=ls())
library(ggplot2)
library(mgcv)

#### Helper functions ####
#iterate function across rows of dataframe
f_lapply_row <- function(df) {
  lapply(seq_len(nrow(df)), function(i) as.list(df[i,,drop=FALSE]))
} ## end function


##### Read in data and create simulated tfactors for GAM ####
#tfFY <- read.csv(file = "outputs/fres_out_FY/frescalo_181219/Output/Trend.csv", header = T, stringsAsFactors = F)
tfFY <- read.csv(file = "outputs/fres_out_FY/frescalo_181224/Output/Trend.csv", header = T, stringsAsFactors = F)
head(tfFY)
tfFY <- tfFY[order(tfFY$Species, tfFY$Time),]
rownames(tfFY) <- 1:nrow(tfFY)

# empty list for species data to go in
spLs <- spLsSim <- vector("list", length(unique(tfFY$Species)))
#names(spLsSim) <- unique(tfDC$Species)

# create list of species specific Frescalo trend data
spLs <- lapply(unique(tfFY$Species), function(x) tfFY[tfFY$Species == x,])
names(spLs) <- unique(tfFY$Species) # add names to list elements

# function designed to be applied to a single set of species trend info
# e.g. test <- tfDC[tfDC$Species == "Acer campestre",]
nsites <- nrow(read.csv(file = "outputs/fres_out_FY/frescalo_181224/Output/Stats.csv", header = T, stringsAsFactors = F))
simRnormSp_FY <- function(x) { temp <- f_lapply_row(x) # convert data frame to list of lists
# sim new data by time period mean and sd
  temp2 <- lapply(temp, function(x) rnorm(nsites, mean = x$TFactor, sd = x$StDev))
  temp3 <- data.frame(do.call(cbind, temp2))
  colnames(temp3) <- c("1972.0", "1977.0", "1982.0", "1987.0", "1992.0", "1997.0", "2002.0", "2007.0", "2012.0", "2016.5")
  temp4 <- reshape2::melt(temp3, measure.vars = c("1972.0", "1977.0", "1982.0", "1987.0", "1992.0", "1997.0", "2002.0", "2007.0", "2012.0", "2016.5"))
  temp4$Var1 <- temp[[1]][1]
  temp4 <- tidyr::unnest(temp4,Var1)
  names(temp4) <- c("year", "tfactor", "taxon")
  temp4$year <- lubridate::year(as.Date(temp4$year, "%Y"))
  return(temp4)
}
# run function, collecting sim data for every taxon
spLsSimFY_v2 <- lapply(spLs, FUN = simRnormSp_FY)
save(spLsSimFY_v2, file = "outputs/sims/fiveyearSims_v2.Rdata")

##### Run GAMs and plot results ####
##############################
# SKip this part, as easier to just fit gams in ggplot
#
# empty list for GAM results
#spGAMs <- vector("list", length(unique(tfFY$Species)))
# apply GAM over species, 3 knots
#spGAMs <- lapply(spLsSim, function(x) mgcv::gam(x$tfactor ~ s(x$year, k = 3), method = "REML"))
##############################

## plotting function
# produce gam plot for every species
#for (x in seq_along(spGAMs)){
#  spGamPlots_FY[[x]] <- plot(spGAMs[[x]], shade = T, rug = T, residuals = T, ylab = c("Tfactor"), xlab = c("Year"), 
#                             main = names(spGAMs)[x], seWithMean = T)
#}
## OK, retreat to using ggplot for the following reason...
#https://stats.stackexchange.com/questions/7795/how-to-obtain-the-values-used-in-plot-gam-in-mgcv
# it's a faff reconstructing plots from listed components of plotted gams...

#### plotting GAMs in ggplot ####
## Ultimately I would like 5 year and DC plots side by side, so, do that in another script (#3) when I have an equivalent list for the 5 year GAMs

# see script 3 for continuation: 3_plotGraphsTogether.R
