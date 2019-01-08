## O.L. Pescott
## started 7 Jan 2019
## Script: 4b_fitGamsForUncertaintyPropagation_DC30_Wal.R
## Wales

## Run GAMs on multiple single draws from Frescalo date-class results for Wales
## resulting in better propagation of uncertainty
## here we fit GAMs to n single draws from the mean/s.d. distribution returned from Frescalo
## and then extract the mean and some quantiles from this distribution of GAMs
##
## This seems a more honest way of capturing uncertainty compared to standard error of a single GAM fitted to larger draw from the mean/sd
## Frescalo time factor. ALternative might be to fit a single GAM to the means only of the Frescalo time factor
## but that, on the other hand, seems to overplay uncertainty, as we don't just have a single sample of the Frescalo time factor,
## the entire population of estimates across every hectad (although these are of course not independent)
##
## 4b: update to include GAM predictions for 1970+

#rm(list=ls())
library(ggplot2)
library(mgcv)
#### Helper functions ####
#iterate function across rows of dataframe
f_lapply_row <- function(df) {
  lapply(seq_len(nrow(df)), function(i) as.list(df[i,,drop=FALSE]))
} ## end function


##### Read in data for Britain and create simulated tfactors for GAM ####
tfDC30_Wal <- read.csv(file = "outputs/fres_out_DC1930_Wal/frescalo_190107/Output/Trend.csv", header = T, stringsAsFactors = F)
head(tfDC30_Wal)
tfDC30_Wal <- tfDC30_Wal[order(tfDC30_Wal$Species, tfDC30_Wal$Time),]
rownames(tfDC30_Wal) <- 1:nrow(tfDC30_Wal)

# empty list for species data to go in
spLs <- spLsSim <- vector("list", length(unique(tfDC30_Wal$Species)))
#names(spLsSim) <- unique(tfDC$Species)

# create list of species specific Frescalo trend data
spLs <- lapply(unique(tfDC30_Wal$Species), function(x) tfDC30_Wal[tfDC30_Wal$Species == x,])
names(spLs) <- unique(tfDC30_Wal$Species) # add names to list elements

#### CHECK SPECIES THAT HAVE ANY ZERO Tfactors ####
spLsCHK <- lapply(spLs, function(x) length(which(x$TFactor==0)))
spLsCHK_TRUE <- which(unlist(spLsCHK) >= 4) # 4 or more Tfactor zeros
attr(spLsCHK_TRUE, "names")

# function designed to be applied to a single set of species trend info
# e.g. test <- tfDC30_Brit[tfDC30_Brit$Species == "Acer campestre",]
times <- data.frame(year = c(1949.5, 1993, 2004.5, 2014))

fitGamsSp <- function(x) { 
  t1 <- f_lapply_row(x) # convert data frame to list of lists (one list for each time period of Fres results)
  t2 <- lapply(t1, function(x) rnorm(100, mean = x$TFactor, sd = x$StDev))
  t3 <- data.frame(do.call(cbind, t2))
  t3 <- t3[,c(1,3,4,5)] # drop 1970-86 DC
  names(t3) <- c("1949.5", "1993", "2004.5", "2014")
  t4 <- f_lapply_row(t3)
  t5 <- list()
  predDat = data.frame(year = seq(1970,2018,1))
  ## Extract fitted values at modelled time points only
  #t5 <- lapply(t4, function(x) {tmp <- data.frame(tfactor = unlist(x), year = times)
  #(mgcv::gam(tmp$tfactor ~ s(tmp$year, k = 3), seWithMean = T))$fitted.values
  #})
  ## Extract predicted values from model fit for all years 1970-2018
    t5 <- lapply(t4, function(x) {
                tmp <- data.frame(tfactor = unlist(x), year = times)
                predict.gam(object = (mgcv::gam(tfactor ~ s(year, k = 3), seWithMean = T, data = tmp)), 
                  newdata = predDat, se.fit = T, newdata.guaranteed = T)
    }
  )
}

### SLOW! ### 
# run function, extracting fitted values for 100 model runs for every taxon (this is 100 x 1353 GAMs = 135,300 models)
#spGAMs <- lapply(spLs, FUN = fitGamsSp)
#save(spGAMs, file = "outputs/gams/spGAMs_v1.Rdata")
#load(file = "outputs/gams/spGAMs_v1.Rdata")

# run function, extracting PREDICTED values for 100 model runs for every taxon (this is 100 x 1353 GAMs = 135,300 models)
spGAMs <- lapply(spLs, FUN = fitGamsSp)
#save(spGAMs, file = "outputs/gams/spGAMs_v2_Wales_Preds.Rdata")
load(file = "outputs/gams/spGAMs_v2_Wales_Preds.Rdata")

## now we have a list of lists, with 100 predicted GAMs for every species
# we now just need to average these predictions
#test <- spGAMs["Acer campestre"][[1]]
#testDf <- spGAMs[1:2]
#as.data.frame(t((test[["Acer campestre"]][[1]][1]))
extractDf <- function(x) {do.call(rbind, lapply(x, function(x) t(as.data.frame(x)[[1]]))) # extract 100 fitted values (= rows) per species (49 years = columns)
}
#test2 <- lapply(testDf, FUN = extractDf)
spGAMsDfs <- lapply(spGAMs, FUN = extractDf)

# extract medians and quantiles per species
spGAMsPredsQs <- lapply(spGAMsDfs, function(x){
  y <- data.frame(x)
  tmp <- apply(y, MARGIN = 2, function(y) quantile(y, probs = c(0.025, 0.5, 0.975)))
  #print( dimnames(tmp))
  dimnames(tmp)[[2]] <- seq(1970,2018,1)
  return(tmp)
})

##### Normalise ####
#test <- spGAMsPredsQs[["Anacamptis laxiflora"]]
#((test/test[2,1])-1)+100
spGAMsPredsQs_FIN <- lapply(spGAMsPredsQs, function(x) (x/x[2,1]))
#### This is broken by Lizard Orchid (spGAMsPredsQs[["Himantoglossum hircinum]] == all zeros) - remove first
spGAMsPredsQs_FIN <- spGAMsPredsQs_FIN[!(names(spGAMsPredsQs_FIN) %in% c("Himantoglossum hircinum"))] ## See NI file for more general approach

save(spGAMsPredsQs_FIN, file = "outputs/DC1930_Wales_normGams.Rdata")



pdf("outputs/plots/Rplots_v3_1930Gams_Wal.pdf",
    width = 4, height = 4, onefile = T)
#x = seq(1970,2018,1)
for (i in names(spGAMsPredsQs_FIN)){
  plot(1,type='n',
       ylim= if( max(spGAMsPredsQs_FIN[[i]]) > 2.0 | min(spGAMsPredsQs_FIN[[i]]) < 0 ) { 
         c(min(spGAMsPredsQs_FIN[[i]]), max(spGAMsPredsQs_FIN[[i]])) } else {
         c(0,2.0)},
       xlim=c(1970,2018),xlab='year', ylab='index', main = names(spGAMsPredsQs_FIN[i]))
for (z in 1:3){
#lines(y = spGAMsPredsQs_FIN[["Acer campestre"]][z,], x = x, type = "l")
#lines(y = spGAMsPredsQs_FIN[["Ajuga chamaepitys"]][z,], x = x, type = "l")
#lines(y = spGAMsPredsQs_FIN[["Agrostis canina s.l."]][z,], x = x, type = "l")
#lines(y = spGAMsPredsQs_FIN[["Anemone nemorosa"]][z,], x = x, type = "l")
#lines(y = spGAMsPredsQs_FIN[["Achillea maritima"]][z,], x = x, type = "l")
#lines(y = spGAMsPredsQs_FIN[["Chrysosplenium oppositifolium"]][z,], x = x, type = "l")
#lines(y = spGAMsPredsQs_FIN[["Gagea lutea"]][z,], x = x, type = "l")
#lines(y = spGAMsPredsQs_FIN[["Galeopsis angustifolia"]][z,], x = x, type = "l")
#lines(y = spGAMsPredsQs_FIN[["Gnaphalium sylvaticum"]][z,], x = x, type = "l")
lines(y = spGAMsPredsQs_FIN[[i]][z,], x = seq(1970,2018,1), type = "l")
#lines(y = spGAMsPredsQs_FIN[["Anacamptis laxiflora"]][z,], x = x, type = "l")
  }
}
dev.off()
##

#### 100 species selection from script 3_plotGraphsTogether.R ####
load(file = "outputs/spGamPlots_DC_100.Rdata")
load(file = "outputs/100SpeciesSample.Rdata")
spGAMsPredsQs_FIN100 <- spGAMsPredsQs_FIN[names(spGAMsPredsQs_FIN) %in% allChkSps]

pdf("outputs/plots/Rplots_v3_1930Gams_Wal100.pdf",
    width = 4, height = 4, onefile = T)
#x = seq(1970,2018,1)
for (i in names(spGAMsPredsQs_FIN100)){
  plot(1,type='n',
       ylim= if( max(spGAMsPredsQs_FIN100[[i]]) > 2.0 | min(spGAMsPredsQs_FIN100[[i]]) < 0 ) { 
         c(min(spGAMsPredsQs_FIN100[[i]]), max(spGAMsPredsQs_FIN100[[i]])) } else {
           c(0,2.0)},
       xlim=c(1970,2018),xlab='year', ylab='index', main = names(spGAMsPredsQs_FIN100[i]))
  for (z in 1:3){
    #lines(y = spGAMsPredsQs_FIN[["Acer campestre"]][z,], x = x, type = "l")
    #lines(y = spGAMsPredsQs_FIN[["Ajuga chamaepitys"]][z,], x = x, type = "l")
    #lines(y = spGAMsPredsQs_FIN[["Agrostis canina s.l."]][z,], x = x, type = "l")
    #lines(y = spGAMsPredsQs_FIN[["Anemone nemorosa"]][z,], x = x, type = "l")
    #lines(y = spGAMsPredsQs_FIN[["Achillea maritima"]][z,], x = x, type = "l")
    #lines(y = spGAMsPredsQs_FIN[["Chrysosplenium oppositifolium"]][z,], x = x, type = "l")
    #lines(y = spGAMsPredsQs_FIN[["Gagea lutea"]][z,], x = x, type = "l")
    #lines(y = spGAMsPredsQs_FIN[["Galeopsis angustifolia"]][z,], x = x, type = "l")
    #lines(y = spGAMsPredsQs_FIN[["Gnaphalium sylvaticum"]][z,], x = x, type = "l")
    lines(y = spGAMsPredsQs_FIN100[[i]][z,], x = seq(1970,2018,1), type = "l")
    #lines(y = spGAMsPredsQs_FIN[["Anacamptis laxiflora"]][z,], x = x, type = "l")
  }
}
dev.off()
##


#### what follows below was for the fitted 4 time point data only ####
##### Example of plotting all 100 fitted mean GAMs for one species ####
#plot(1,type='n',ylim=c(0.35,0.50),xlim=c(1930,2020),xlab='year', ylab='tfactor')
#for (z in 1:100){
#lines(y = spGAMs[["Acer campestre"]][[z]][1:4], x = c(1949.5,1978,1993,2004.5,2014), type = "l")
#}
# compare to 
#load(file = "outputs/plots/spGamPlots_DC_v2.Rdata")
#spGamPlots_DC["Allium ursinum"] # pretty similar

#### Calculate medians and quantiles for plotting/sending to Dan Hayhow ####
#extractDf <- function(x) {do.call(rbind, lapply(x, function(x) t(as.data.frame(x))))
#}
#spGAMsDfs <- lapply(spGAMs, FUN = extractDf)

#spGAMsQuants <- lapply(spGAMsDfs, function(x){
#  y <- data.frame(x)
#  apply(y, MARGIN = 2, function(y) quantile(y, probs = c(0.025, 0.5, 0.975)))
#})
#save(spGAMsQuants, file = "outputs/gams/spGAMsQuants.Rdata")
#load(file = "outputs/gams/spGAMsQuants.Rdata")


#### Base plot example ####
#test <- data.frame(spGAMsQuants[["Allium ursinum"]])
#names(test) <- c(1978,1993,2004.5,2014)
#plot(x = names(test), y = test[2,], ylim = c(min(test[1,]), max(test[3,]) ))
#arrows(as.numeric(names(test)), as.numeric(test[1,]), as.numeric(names(test)), as.numeric(test[3,]), length=0.05, angle=90, code=3)

#test <- data.frame(spGAMsQuants[["Achillea maritima"]])
#names(test) <- c(1978,1993,2004.5,2014)
#plot(x = names(test), y = test[2,], ylim = c(min(test[1,]), max(test[3,]) ))
#arrows(as.numeric(names(test)), as.numeric(test[1,]), as.numeric(names(test)), as.numeric(test[3,]), length=0.05, angle=90, code=3)
####

## Quick check on missing data for Wales
for (i in names(spGAMsPredsQs_FIN100)){
if( max(spGAMsPredsQs_FIN[[i]])=="NaN" ) { 
  print(i)} else {
    NULL}
}
