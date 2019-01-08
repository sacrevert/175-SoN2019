## O.L. Pescott
## started 20 Dec 2018
## Script: 3_plotGraphsTogether.R

## Plot date class-based and five-year window-based GAMs together for comparison
## Also, extract 100 random species, stratified by sqrt(range size), for checking (with Kevin Walker and Pete Stroh, BSBI)

#rm(list=ls())
library(ggplot2)
library(gridExtra)

#### Load data, make graphs ####
## load data
load(file = "outputs/sims/fiveyearSims_v2.Rdata")
load(file = "outputs/sims/dateclassSims_v2.Rdata")

spGamPlots_FY <- vector("list", length(seq_along(spLsSimFY_v2)))
spGamPlots_DC <- vector("list", length(seq_along(spLsSimDC_v2)))

# apply GAM over species, 3 knots
for (x in seq_along(spLsSimFY_v2)){
  spGamPlots_FY[[x]] <- ggplot(spLsSimFY_v2[[x]], aes(x = year, y = tfactor)) +
    geom_point(alpha = 0.2) + 
    geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) +
    #geom_smooth(method = "gam", formula = y ~ s(x, k = 3), method.args=(list(seWithMean = T, unconditional = T))) +
    ggtitle(label = names(spLsSimFY_v2[x]))
  names(spGamPlots_FY) <- names(spLsSimFY_v2)
}

for (x in seq_along(spLsSimDC_v2)){
  spGamPlots_DC[[x]] <- ggplot(spLsSimDC_v2[[x]], aes(x = year, y = tfactor)) +
    geom_point(alpha = 0.2) + 
    geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) +
    ggtitle(label = names(spLsSimDC_v2[x]))
  names(spGamPlots_DC) <- names(spLsSimDC_v2)
}

###
#save(spGamPlots_FY, file = "outputs/plots/spGamPlots_FY_v2.Rdata")
#save(spGamPlots_DC, file = "outputs/plots/spGamPlots_DC_v2.Rdata")
load(file = "outputs/plots/spGamPlots_FY_v2.Rdata")
load(file = "outputs/plots/spGamPlots_DC_v2.Rdata")
###

#### select 100 species for checking! #####
# stratified random set
# get hectad counts first
hecDat <- read.csv(file = "data/bsbi/hecCount_1970_2018.csv", header = T, stringsAsFactors = F)
hecDat$taxon <- ifelse(hecDat$qualifier == "", hecDat$ï..group, paste(hecDat$ï..group, hecDat$qualifier))
hist(sqrt(hecDat$freq)) # flatten distribution
quantile(sqrt(hecDat$freq), probs = seq(0, 1, 0.25))
hecDat$srFreq <- sqrt(hecDat$freq)
set.seed(1001)
sampSp1 <- sample(hecDat[(hecDat$srFreq > 0) & (hecDat$srFreq < 11.83),4], 25, replace = F)
sampSp2 <- sample(hecDat[(hecDat$srFreq > 11.83) & (hecDat$srFreq < 27.65),4], 25, replace = F)
sampSp3 <- sample(hecDat[(hecDat$srFreq > 27.65) & (hecDat$srFreq < 42.59),4], 25, replace = F)
sampSp4 <- sample(hecDat[(hecDat$srFreq > 42.59) & (hecDat$srFreq < 53.24),4], 25, replace = F)
allChkSps <- c(sampSp1, sampSp2, sampSp3, sampSp4)
#save(allChkSps, file = "outputs/100SpeciesSample.Rdata")

spGamPlots_FY_100 <- spGamPlots_FY[names(spGamPlots_FY) %in% allChkSps]
spGamPlots_DC_100 <- spGamPlots_DC[names(spGamPlots_DC) %in% allChkSps]
#save(spGamPlots_DC_100, file = "outputs/spGamPlots_DC_100.Rdata")

#### Select chosen species to print graphs for ####
pdf("outputs/plots/Rplots_v2.1.pdf",
    width = 6, height = 4, onefile = T)
for (z in names(spGamPlots_DC_100)){
      do.call(grid.arrange, c(spGamPlots_DC_100[z], spGamPlots_FY_100[z], ncol = 2))
  }
dev.off()




