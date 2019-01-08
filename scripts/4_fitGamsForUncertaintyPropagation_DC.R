## O.L. Pescott
## started 21 Dec 2018
## Script: 4_fitGamsForUncertaintyPropagation_DC.R

## Run GAMs on multiple single draws from Frescalo date-class results
## resulting in better propagation of uncertainty
## here we fit GAMs to n single draws from the mean/s.d. distribution returned from Frescalo
## and then extract the mean and some quantiles from this distribution of GAMs
##
## This seems a more honest way of capturing uncertainty compared to standard error of a single GAM fitted to larger draw from the mean/sd
## Frescalo time factor. ALternative might be to fit a single GAM to the means only of the Frescalo time factor
## but that, on the other hand, seems to overplay uncertainty, as we don't just have a single sample of the Frescalo time factor,
## the entire population of estimates across every hectad (although these are of course not independent)

#rm(list=ls())
library(ggplot2)
library(mgcv)
#### Helper functions ####
#iterate function across rows of dataframe
f_lapply_row <- function(df) {
  lapply(seq_len(nrow(df)), function(i) as.list(df[i,,drop=FALSE]))
} ## end function


##### Read in data and create simulated tfactors for GAM ####
tfDC <- read.csv(file = "outputs/fres_out_DC/frescalo_181224/Output/Trend.csv", header = T, stringsAsFactors = F)
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
times <- data.frame(year = c(1978, 1993, 2004.5, 2014))

fitGamsSp <- function(x) { 
  t1 <- f_lapply_row(x) # convert data frame to list of lists
  t2 <- lapply(t1, function(x) rnorm(100, mean = x$TFactor, sd = x$StDev))
  t3 <- data.frame(do.call(cbind, t2))
  names(t3) <- c("1978", "1993", "2004.5", "2014")
  t4 <- f_lapply_row(t3)
  t5 <- list()
  t5 <- lapply(t4, function(x) {tmp <- data.frame(tfactor = unlist(x), year = times)
  (mgcv::gam(tmp$tfactor ~ s(tmp$year, k = 3), seWithMean = T))$fitted.values
  })
}

### SLOW! ### 
# run function, extracting fitted values for 100 model runs for every taxon (this is 100 x 1353 GAMs = 135,300 models)
spGAMs <- lapply(spLs, FUN = fitGamsSp)
#save(spGAMs, file = "outputs/gams/spGAMs_v1.Rdata")
load(file = "outputs/gams/spGAMs_v1.Rdata")

##### Example of plotting all 100 fitted mean GAMs for one species ####
plot(1,type='n',ylim=c(0.35,0.50),xlim=c(1977,2015),xlab='year', ylab='tfactor')
for (z in 1:100){
lines( y = spGAMs[["Acer campestre"]][[z]][1:4], x = c(1978,1993,2004.5,2014), type = "l")
}
# compare to 
load(file = "outputs/plots/spGamPlots_DC_v2.Rdata")
spGamPlots_DC["Allium ursinum"] # pretty similar

#### Calculate medians and quantiles for plotting/sending to Dan Hayhow ####
extractDf <- function(x) {do.call(rbind, lapply(x, function(x) t(as.data.frame(x))))
}
spGAMsDfs <- lapply(spGAMs, FUN = extractDf)

spGAMsQuants <- lapply(spGAMsDfs, function(x){
  y <- data.frame(x)
  apply(y, MARGIN = 2, function(y) quantile(y, probs = c(0.025, 0.5, 0.975)))
})
#save(spGAMsQuants, file = "outputs/gams/spGAMsQuants.Rdata")
load(file = "outputs/gams/spGAMsQuants.Rdata")


#### Base plot example ####
test <- data.frame(spGAMsQuants[["Allium ursinum"]])
names(test) <- c(1978,1993,2004.5,2014)
plot(x = names(test), y = test[2,], ylim = c(min(test[1,]), max(test[3,]) ))
arrows(as.numeric(names(test)), as.numeric(test[1,]), as.numeric(names(test)), as.numeric(test[3,]), length=0.05, angle=90, code=3)
####
