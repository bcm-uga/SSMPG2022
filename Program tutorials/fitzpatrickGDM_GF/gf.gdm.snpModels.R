################################################################################
# Scripts to model & map genetic turnover using gradient forests and generalized
# dissimilarity modeling described in Fitzpatrick & Keller (2015) Ecology Letters
################################################################################

# Mapping spatial genetic variation --------------------------------------------
###### functions to support mapping #####
# builds RGB raster from transformed environment from  GF or GDM
# snpPreds = dataframe of transformed variables from GF or GDM model
# rast = a raster mask to which RGB values are to be mapped
# mapCells = cell IDs to which RGB values should be assigned
pcaToRaster <- function(snpPreds, rast, mapCells){
  require(raster)
  
  pca <- prcomp(snpPreds, center=TRUE, scale.=FALSE)
  
  ##assigns to colors, edit as needed to maximize color contrast, etc.
  a1 <- pca$x[,1]; a2 <- pca$x[,2]; a3 <- pca$x[,3]
  r <- a1+a2; g <- -a2; b <- a3+a2-a1
  
  ##scales colors
  scalR <- (r-min(r))/(max(r)-min(r))*255
  scalG <- (g-min(g))/(max(g)-min(g))*255
  scalB <- (b-min(b))/(max(b)-min(b))*255
  
  ##assigns color to raster
  rast1 <- rast2 <- rast3 <- rast
  rast1[mapCells] <- scalR
  rast2[mapCells] <- scalG
  rast3[mapCells] <- scalB
  ##stacks color rasters
  outRast <- stack(rast1, rast2, rast3)
  return(outRast)
}

# Function to map difference between spatial genetic predictions
# predMap1 = dataframe of transformed variables from gf or gdm model for first set of SNPs
# predMap2 = dataframe of transformed variables from gf or gdm model for second set of SNPs
# rast = a raster mask to which Procrustes residuals are to be mapped
# mapCells = cell IDs to which Procrustes residuals values should be assigned
RGBdiffMap <- function(predMap1, predMap2, rast, mapCells){
  require(vegan)
  PCA1 <- prcomp(predMap1, center=TRUE, scale.=FALSE)
  PCA2 <- prcomp(predMap2, center=TRUE, scale.=FALSE)
  diffProcrust <- procrustes(PCA1, PCA2, scale=TRUE, symmetrical=FALSE)
  residMap <- residuals(diffProcrust)
  rast[mapCells] <- residMap
  return(list(max(residMap), rast))
}



################################################################################
# BEGIN SCRIPTS FOR GRADIENT FOREST MODELING
################################################################################
# load libraries, read in data, etc. -------------------------------------------
library(gradientForest)
library(raster)
library(dplyr)
library(colorRamps)

# read in data file with minor allele freqs & env/space variables
gfData <- read.csv("~/poplarSNP.ENV.data.4.GF.csv")
# in this tutorial, we will be using bioclimatic variables only
envGF <- gfData %>% select(contains("bio")) # bioclimatic covariates

# build different SNP datasets for modeling / comparision
SNPs_ref <- gfData %>% select(contains("REFERENCE")) # reference
SNPs_cand <- gfData %>% select(contains("GI5")) # GI5
################################################################################


# GRADIENT FOREST MODELING -----------------------------------------------------
maxLevel <- log2(0.368*nrow(envGF)/2) #account for correlations, see ?gradientForest 

# Fit gf models for reference SNPs 
gfRef <- gradientForest(data=cbind(envGF, SNPs_ref), 
                        predictor.vars=colnames(envGF),
                        response.vars=colnames(SNPs_ref), 
                        ntree=500, 
                        maxLevel=maxLevel, 
                        trace=T, 
                        corr.threshold=0.50)

# Fit gf models for candidate SNPs
gfCand <- gradientForest(data=cbind(envGF, SNPs_cand), 
                        predictor.vars=colnames(envGF),
                        response.vars=colnames(SNPs_cand), 
                        ntree=500, 
                        maxLevel=maxLevel, 
                        trace=T, 
                        corr.threshold=0.50)

# plot output, see ?plot.gradientForest
types <- c("Overall.Importance", 
           "Split.Density", 
           "Cumulative.Importance", 
           "Performance")
sel <- 1
plot(gfRef, plot.type=types[sel])
plot(gfCand, plot.type=types[sel])

# Cumulative importance plot
plot(gfCand, plot.type="C", 
     show.species=F, 
     common.scale=T, 
     imp.vars=names(sort(gfCand$overall.imp2, decreasing = T)), 
     lwd=2, 
     col="blue",
     cex.lab=1.5,
     cex.axis=0.6, 
     line.ylab=0.9,
     par.args=list(mgp=c(1.5, 0.5, 0), 
                   mar=c(2.5,1.0,1,1), 
                   omi=c(0,0.3,0,0)))
################################################################################

# prep spatial data
library(gtools)
# polygon of Balsam poplar geographic range
balm <- shapefile("~/popubals_Dissolve.shp")

# get climate rasters
# uncomment to download from online database
#library(geodata)
#climRasts <- worldclim_global(var="bio",
#                              res=5,
#                              path="./alpsSummerSchool")

# otherwise, load from file
rastFiles <- list.files(path="~/wc2.1_5m",
                        full.names=T,
                        pattern="_bio_")
#get files in correct order
rastFiles <- mixedsort(rastFiles)

# stack, subset, and clip to species range polygon
climRasts <- stack(rastFiles)
names(climRasts) <- paste("bio", 1:19, sep="_")
climRasts <- climRasts[[which(names(climRasts) %in% colnames(envGF))]]
climRasts <- raster::mask(stack(climRasts), balm)
climRasts <- raster::trim(climRasts)
plot(climRasts, col=rgb.tables(1000))

# create mask raster for clipping & mapping below
# trim() did not work on the future climate rasters
# for some reason...
# empty raster mask is used for mapping 
mask <- climRasts[[1]]>-100 

# download / prep future climate data
#cmip6Rasts <- stack(cmip6_world(model="UKESM1-0-LL",
#                                ssp="585",
#                                time="2061-2080",
#                                var="bioc",
#                                res=5,
#                                path="./alpsSummerSchool"))

# same as above, but for a future climate scenario
cmip6Rasts <- stack("~/wc2.1_5m_bioc_UKESM1-0-LL_ssp585_2061-2080.tif")
names(cmip6Rasts) <- paste("bio", 1:19, sep="_")
cmip6Rasts <- cmip6Rasts[[which(names(cmip6Rasts) %in% colnames(envGF))]]
cmip6Names <- names(cmip6Rasts)
cmip6Rasts <- cmip6Rasts*mask
names(cmip6Rasts) <- cmip6Names
plot(cmip6Rasts, col=rgb.tables(1000))

# Mapping spatial patterns - Script assumes:
# (1) a dataframe named env_trns containing extracted raster data (w/ cell IDs)
# and env. variables used in the models & with columns as follows: cell, bio1, bio2, etc.
#
# (2) a raster mask of the study region to which the RGB data will be written
mask[] <- NA # set all values to NA for better mapping 

# fix problems due to current & future data having different
# NA cells, often an issue with these types of data
# different numbers of NA cells, trouble...
sum(is.na(climRasts[[1]][]))
sum(is.na(cmip6Rasts[[1]][]))

# make both raster stacks have the same number of
# NAs in the same cells
climRasts[[6]][which(is.na(cmip6Rasts[[1]][]))] <- NA
climRasts[[5]][which(is.na(cmip6Rasts[[1]][]))] <- NA
climRasts[[4]][which(is.na(cmip6Rasts[[1]][]))] <- NA
climRasts[[3]][which(is.na(cmip6Rasts[[1]][]))] <- NA
climRasts[[2]][which(is.na(cmip6Rasts[[1]][]))] <- NA
climRasts[[1]][which(is.na(cmip6Rasts[[1]][]))] <- NA

# GF functions will not take rasters, so need to extract data from
# the rasters (cell ID as well)
env_trns <- raster::extract(climRasts, 1:ncell(climRasts[[1]]))
env_trns <- data.frame(cell=1:ncell(climRasts[[1]]), env_trns)
env_trns <- na.omit(env_trns)

# same for future
env_trns_future <- raster::extract(cmip6Rasts, 1:ncell(cmip6Rasts[[1]]))
env_trns_future <- data.frame(cell=1:ncell(cmip6Rasts[[1]]), env_trns_future)
env_trns_future <- na.omit(env_trns_future)

# check congruence
dim(env_trns)
dim(env_trns_future)

# transform env using gf models, see ?predict.gradientForest
predRef <- predict(gfRef, env_trns[,-1]) # remove cell column before transforming
predCand <- predict(gfCand, env_trns[,-1])

# map continuous variation - reference SNPs
refRGBmap <- pcaToRaster(predRef, mask, env_trns$cell)
plotRGB(refRGBmap)
#writeRaster(refRGBmap, "/.../refSNPs_map.tif", format="GTiff", overwrite=TRUE)

# map continuous variation - Candidate SNPs
candRGBmap <- pcaToRaster(predCand, mask, env_trns$cell)
plotRGB(candRGBmap)
#writeRaster(candRGBmap, "/.../candSNPs_map.tif", format="GTiff", overwrite=TRUE)

# Difference between maps (Candidate and Reference spatial patterns) 
diffCand <- RGBdiffMap(predRef, predCand, rast=mask, mapCells=env_trns$cell)
plot(diffCand[[2]])
#writeRaster(diffCand[[2]], "/.../diffRef_Cand.tif", format="GTiff", overwrite=TRUE)
################################################################################


# Calculate and map "local genetic offset" under climate change ----------------------
# Script assumes:
  # (1) a dataframe of transformed env. variables for CURRENT climate 
  # (e.g., predCand from above).
  #
  # (2) a dataframe named 'env_trns_future' containing extracted raster data of 
  # env. variables for FUTURE a climate scenario, same structure as env_trns

# first transform FUTURE env. variables
transCand <- predict(gfCand, env_trns_future[,-1])

# calculate Euclidean distance between current and future 
# genetic spaces  
genOffsetCand <- sqrt(rowSums((transCand-predCand)^2))
 
# assign values to raster
genOffset <- mask
genOffset[env_trns_future$cell] <- genOffsetCand
plot(genOffset, col=heat.colors(1000)[1000:1])
################################################################################
# END SCRIPTS FOR GRADIENT FOREST MODELING
################################################################################

#-------------------------------------------------------------------------------

################################################################################
# BEGIN SCRIPTS FOR GENERALIZED DISSIMILARITY MODELING
################################################################################
# load libraries, read in data, etc. -------------------------------------------
library(gdm)
library(raster)
# read in data file with Fst & env variables (in GDM site-pair format)
# note Fst values were scaled 0-1 to facilitate model convergence
gdmData <- read.csv("~/poplarFst.ENV.data.4.GDM.csv")

# build individual SNP datasets
SNPs_ref <- gdmData[,c(1,6:22)] # reference
SNPs_cand <- gdmData[,c(2,6:22)] # candidate
################################################################################


# GENERALIZED DISSIMILARITY MODELING -----------------------------------------------------
dev.off()
GEO <- T # use geographic distance as a predictor?
# reference SNPs
# can ignore warning here
gdmRef <- gdm(SNPs_ref, geo=GEO)
gdmRef$explained # % deviance explained
plot(gdmRef, plot.layout=c(3,4))
#refSplines <- isplineExtract(gdmRef) # extract spline data for custom plotting

# candidate SNPs
gdmCand <- gdm(SNPs_cand, geo=GEO)
gdmCand$explained
plot(gdmCand, plot.layout=c(3,4))
#candSplines <- isplineExtract(gdmCand) # extract spline data for custom plotting
################################################################################


# Mapping spatial genetic variation --------------------------------------------
dev.off()
# Follows identical procedure as above for GF, but uses ?gdm.transform instead
# Also note that if geo=T when fitting model, 
# x & y rasters can be supplied in addition to env. layers.
transRef <- gdm.transform(gdmRef, env_trns[,-1]) # remove cell column
transCand <- gdm.transform(gdmCand, env_trns[,-1])

# map continuous variation - reference SNPs
refRGBmap <- pcaToRaster(transRef, mask, env_trns$cell)
plotRGB(refRGBmap)
#writeRaster(refRGBmap, "/.../refSNPs_map.tif", format="GTiff", overwrite=TRUE)

# map continuous variation - candidate SNPs
candRGBmap <- pcaToRaster(transCand, mask, env_trns$cell)
plotRGB(candRGBmap)
#writeRaster(candRGBmap, "/.../candSNPs_map.tif", format="GTiff", overwrite=TRUE)

# Difference between maps (candidate and reference) 
diffCand <- RGBdiffMap(transRef, transCand, rast=mask, mapCells=env_trns$cell)
plot(diffCand[[2]])
#writeRaster(diffCand[[2]], "/.../diffRef_GI5.tif", format="GTiff", overwrite=TRUE)
################################################################################

# Calculate and map "local genetic offset" under climate change ----------------------
# Unlike GF, we do not calculate Euclidean distances
# between current & future transformed envirnments. Instead
# we predict offsets through time directly from the 
# fitted model.
# Note that the predict.gdm function can take rasters directly,
# which saves some work. Here, for simplicity, I fit 
# a new model with geo=F.

gdmCand <- gdm(SNPs_cand, geo=FALSE)
genOffsetCand <- predict(object=gdmCand, 
                         data=climRasts, 
                         time=TRUE, 
                         predRasts=cmip6Rasts)

plot(genOffsetCand, col=heat.colors(1000)[1000:1])
################################################################################
# END SCRIPTS FOR GENERALIZED DISSIMILARITY MODELING
################################################################################
