###################################
### Getting environmental data  ###
###-----------------------------###
#### Install and load packages ####
# install.packages("XML")
library(raster)
library(rgeos)
library(rgdal)
library(gdalUtils)
library(XML)
library(maps)
library(sf)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggpubr)
sessionInfo()

#### Clear workspace ####
rm(list=ls())
ls()

#### Set paths and working directory ####
path1 <- "/Users/dauphin/Documents/GitLab/prj/EnvironmentalData/res/"
path2 <- "/Users/dauphin/Documents/GitLab/prj/EnvironmentalData/dat/"
path3 <- "/Users/dauphin/Documents/GitLab/prj/EnvironmentalData/doc/"
path5 <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/"
setwd(path1)

#### Define variables ####
eu.extent <- c(-15,35,32,72) # define the European extent
bio.nb1 <- 1:19 # define the bioclim numbering
bio.nb2 <- sprintf("%02d", 1:19)

#### Get a background map and project it ####
# Download manually the raster layer and paste the zip file into the dat/ folder
browseURL("http://www.naturalearthdata.com/downloads/10m-raster-data/10m-natural-earth-2/",
          browser=getOption("browser"), encodeIfNeeded=F)
unzip(zipfile=paste(path2,"NE2_LR_LC.zip",sep=""), exdir=paste(path2,"NE2_LR_LC",sep=""))

# Load it
back.map <- stack(paste(path2,"NE2_LR_LC/NE2_LR_LC.tif",sep=""))

# Define geographic coordinate system
browseURL("https://epsg.io", browser=getOption("browser"), encodeIfNeeded=F) # here are the description of coordinate reference systems
crs(back.map) <- CRS("+init=EPSG:4326")

# Global map
plotRGB(back.map, axes=F)

# Europe map
plotRGB(back.map, ext=eu.extent, axes=F)

#### Get occurrence data and species range ####
# Download manually tree occurrences and extent (shapefile) and paste it into the dat/ folder
browseURL("https://figshare.com/collections/A_high-resolution_pan-European_tree_occurrence_dataset/3288407/1",
          browser=getOption("browser"), encodeIfNeeded=F)

# Download manually tree species range (shapefile) and paste it into the dat/ folder
browseURL("https://www.euforgen.org/species/fagus-sylvatica/",
          browser=getOption("browser"), encodeIfNeeded=F)

# Load tree occurrence data
occ <- read.table(paste(path2,"EUForestTreeOccurrences/EUForestspecies.csv",sep=""), header=T, sep=",")
occ.fs <- data.frame(occ[occ$SPECIES.NAME=="Fagus sylvatica",])
head(occ.fs); dim(occ.fs)
plot(occ.fs$X, occ.fs$Y, cex=0.2)

# Transform coordinates from ETRS89 to WGS84
browseURL("https://epsg.io", browser=getOption("browser"), encodeIfNeeded=F)
d <- data.frame(lon=occ.fs$X, lat=occ.fs$Y)
coordinates(d) <- c("lon", "lat")
proj4string(d) <- CRS("+init=epsg:3035") # ETRS89 / ETRS-LAEA
CRS.new <- CRS("+init=epsg:4326") # WGS84
d.wgs84 <- spTransform(d, CRS.new)
str(d.wgs84)
cd <- data.frame(d.wgs84@coords)
occ.fs$lonwgs84 <- cd$lon
occ.fs$latwgs84 <- cd$lat

# Load the species' occurrence (rg.fs1) and natural range (rg.fs2)
rg.fs1 <- readOGR(dsn=paste(path2,"EUForestTreeRanges/",sep=""), 
                  layer="Fagus_sylvatica", verbose=F)
rg.fs2 <- readOGR(dsn=paste(path2,"FsylSpeciesRangeEuforgen/shapefiles/",sep=""), 
                  layer="Fagus_sylvatica_sylvatica_plg", verbose=F)
crs(rg.fs1) <- CRS('+init=EPSG:4326')
crs(rg.fs2) <- CRS('+init=EPSG:4326')
rg.fs1 <- raster::crop(rg.fs1, eu.extent)
rg.fs2 <- raster::crop(rg.fs2, eu.extent)

# Visualise occurrence data and species' natural range
pdf(paste(path1,"FS_OccurrencesAndRange.pdf",sep=""), width=7.5, height=6.15)
par(mfrow=c(1,1), oma=c(4,4,0.5,0.5), mai=c(0,0,0,0), mar=c(0,0,0,0), mgp=c(0,1,0))
plotRGB(back.map, ext=eu.extent, alpha=250, axes=F)
map.scale(20, 36, relwidth=0.15, metric=T, ratio=F)
points(occ.fs$lonwgs84, occ.fs$latwgs84, pch=16, cex=0.2, col="gold")
plot(rg.fs1, col=alpha("red",0.0), lwd=1.5, lty=1, border="red", add=T)
plot(rg.fs2, col=alpha("blue",0.0), lwd=1.5, lty=1, border="blue", add=T)
axis(side=2, at=seq(32,72,10), lwd=0.5, cex.axis=0.8)
axis(side=1, at=seq(-15,35,10), lwd=0.5, cex.axis=0.8)
text(-20, 52, label="Latitude [°]", cex=1.2, srt=90, xpd=NA, pos=3)
text(15, 28, label="Longitude [°]", cex=1.2, xpd=NA, pos=2)
dev.off()

#### Download environmental data ####
# Climate -- download rasters from CHELSA
browseURL("https://chelsa-climate.org/downloads/",
          browser=getOption("browser"), encodeIfNeeded=F)
getOption("timeout") # check the time out setting
options(timeout=1000) # set a longer time if needed
for(i in 1:length(bio.nb1)){ # length(bio.nb1)
  download.file(paste(path5,"CHELSA_bio",bio.nb1[i],"_1981-2010_V.2.1.tif",sep=""),
                paste(path2,"CHELSA/CHELSA_bio",bio.nb2[i],"_1981-2010_V.2.1.tif",sep=""))
} # download climatic layers

download.file("https://chelsa-climate.org/wp-admin/download-page/CHELSA_tech_specification_V2.pdf",
              paste(path3,"CHELSA/CHELSA_tech_specification_V2.pdf",sep="")) # download the documentation too

# Topography -- download raster from GMTED2010
browseURL("https://topotools.cr.usgs.gov/gmted_viewer/gmted2010_global_grids.php",
          browser=getOption("browser"), encodeIfNeeded=F)
download.file("https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/topo/downloads/GMTED/Grid_ZipFiles/mn30_grd.zip",
              paste(path2,"mn30_grd.zip",sep=""))
unzip(zipfile=paste(path2,"mn30_grd.zip",sep=""), exdir=paste(path2,"mn30_grd",sep=""))

# Soil -- download rasters from ISRIC
browseURL("https://data.isric.org/geonetwork/srv/eng/catalog.search#/metadata/2fdf1a33-487f-4aa7-a6b6-b2b49483ad48",
          browser=getOption("browser"), encodeIfNeeded=F)
download.file("https://files.isric.org/soilgrids/latest/data_aggregated/1000m/phh2o/phh2o_5-15cm_mean_1000.tif",
              paste(path2,"SoilGrid/phh2o_5-15cm_mean_1000.tif",sep=""))
download.file("https://files.isric.org/soilgrids/latest/data_aggregated/1000m/phh2o/phh2o_60-100cm_mean_1000.tif",
              paste(path2,"SoilGrid/phh2o_60-100cm_mean_1000.tif",sep=""))

#### Extract environmental data from occurrences ####
head(occ.fs); dim(occ.fs)
env.fs <- occ.fs
for(i in 1:length(bio.nb2)){
  print(paste("Extraction of Bio", bio.nb2[i], sep=""))
  tmp.ras <- NULL; tmp.res <- NULL
  tmp.ras <- raster::raster(paste(path2, "CHELSA/CHELSA_bio", bio.nb2[i],"_1981-2010_V.2.1.tif", sep=""))
  tmp.res <- raster::extract(tmp.ras, env.fs[, c("lonwgs84","latwgs84")], "bilinear")
  env.fs[,i+12] <- tmp.res
  colnames(env.fs)[i+12] <- paste("Bio", bio.nb2[i], sep="")
  assign(paste("bio", bio.nb2[i], sep=""), tmp.ras)
}

gmted <- raster::raster(paste(path2, "mn30_grd/mn30_grd/w001000.adf", sep=""))
env.fs$gmted <- raster::extract(gmted, env.fs[, c("lonwgs84","latwgs84")], "bilinear")

pH.top <- raster::raster(paste(path2, "SoilGrid/phh2o_5-15cm_mean_1000.tif", sep="")) # ISRIC layers use the Homolosine projection
pH.top <- projectRaster(pH.top, crs=CRS("+init=EPSG:4326"))
raster::plot(pH.top)
env.fs$pH.top <- raster::extract(pH.top, env.fs[, c("lonwgs84","latwgs84")], "bilinear") / 10
pH.dee <- raster::raster(paste(path2, "SoilGrid/phh2o_60-100cm_mean_1000.tif", sep="")) # ISRIC layers use the Homolosine projection
pH.dee <- projectRaster(pH.dee, crs=CRS("+init=EPSG:4326"))
raster::plot(pH.dee)
env.fs$pH.dee <- raster::extract(pH.dee, env.fs[, c("lonwgs84","latwgs84")], "bilinear") / 10

head(env.fs)
env.fs <- env.fs[,-c(3:10)]
write.table(env.fs, paste(path1, "FS_OccurrencesEnvDatExtracted.csv", sep=""),
            row.names=F, col.names=T, quote=F, sep=";") # save table with extracted environmental data at sampling locations

#### Crop layers for EU extent ####
bio01.c <- crop(bio01, eu.extent)
bio02.c <- crop(bio02, eu.extent)
bio03.c <- crop(bio03, eu.extent)
gmted.c <- crop(gmted, eu.extent)
pH.top.c <- crop(pH.top, eu.extent)
pH.dee.c <- crop(pH.dee, eu.extent)

#### Visualise environmental data and explore ecological gradients ####
pdf(paste(path1,"FS_OccurrencesEnvDatExtracted1.pdf",sep=""), width=7.5*1.2, height=6.15*1.2)
par(mfrow=c(2,3))
plot(bio01.c, main=paste("Bio",bio.nb2[1],sep=""), alpha=150)
plot(bio02.c, main=paste("Bio",bio.nb2[2],sep=""), alpha=150)
plot(bio03.c, main=paste("Bio",bio.nb2[3],sep=""), alpha=150)
smoothScatter(env.fs$Bio01, env.fs$Bio02, main="Bio1 & Bio2")
smoothScatter(env.fs$Bio01, env.fs$Bio03, main="Bio1 & Bio3")
smoothScatter(env.fs$Bio02, env.fs$Bio03, main="Bio2 & Bio3")
dev.off()

pdf(paste(path1,"FS_OccurrencesEnvDatExtracted2.pdf",sep=""), width=7.5*1.2, height=6.15*1.2)
par(mfrow=c(2,3))
plot(bio01.c, main=paste("Bio",bio.nb2[1],sep=""), alpha=150)
plot(gmted.c, main=paste("GMTED",sep=""), alpha=150)
plot(pH.top.c, main=paste("pH",sep=""), alpha=150)
smoothScatter(env.fs$Bio01, env.fs$gmted, main="Bio1 & GMTED")
smoothScatter(env.fs$Bio01, env.fs$pH.top, main="Bio1 & pH")
smoothScatter(env.fs$gmted, env.fs$pH.top, main="GMTED & pH")
dev.off()

#### Sampling locations, data extraction and visualisation ####
sites <- read.table(paste(path2,"GenTree/FS_geo_coordinates.csv",sep=""), header=T, sep=",")
head(sites); dim(sites)

for(i in 1:length(bio.nb2)){
  print(paste("Extraction of Bio", bio.nb2[i], sep=""))
  tmp.ras <- NULL; tmp.res <- NULL
  tmp.ras <- raster::raster(paste(path2, "CHELSA/CHELSA_bio", bio.nb2[i],"_1981-2010_V.2.1.tif", sep=""))
  tmp.res <- raster::extract(tmp.ras, sites[, c("lonwgs84","latwgs84")], "bilinear")
  sites[,i+8] <- tmp.res
  colnames(sites)[i+8] <- paste("Bio", bio.nb2[i], sep="")
}

sites$gmted <- raster::extract(gmted.c, sites[, c("lonwgs84","latwgs84")], "bilinear")
sites$pH.top <- raster::extract(pH.top.c, sites[, c("lonwgs84","latwgs84")], "bilinear") / 10
sites$pH.dee <- raster::extract(pH.dee.c, sites[, c("lonwgs84","latwgs84")], "bilinear") / 10
head(sites); dim(sites)
write.table(sites, paste(path1, "FS_SamplingLocationsEnvDatExtracted.csv", sep=""),
            row.names=F, col.names=T, quote=F, sep=";")

pdf(paste(path1,"FS_OccurrencesAndSitesEnvDatExtracted1.pdf",sep=""), width=7.5*1.2, height=6.15*1.2)
par(mfrow=c(2,3))
plot(bio01.c, main=paste("Bio",bio.nb2[1],sep=""), alpha=150)
points(sites$lonwgs84, sites$latwgs84, cex=1.4, pch=21, bg="red", col="white")
plot(bio02.c, main=paste("Bio",bio.nb2[2],sep=""), alpha=150)
points(sites$lonwgs84, sites$latwgs84, cex=1.4, pch=21, bg="red", col="white")
plot(bio03.c, main=paste("Bio",bio.nb2[3],sep=""), alpha=150)
points(sites$lonwgs84, sites$latwgs84, cex=1.4, pch=21, bg="red", col="white")
smoothScatter(env.fs$Bio01, env.fs$Bio02, main="Bio1 & Bio2")
points(sites$Bio01, sites$Bio02, cex=1.0, pch=21, bg="red", col="white")
smoothScatter(env.fs$Bio01, env.fs$Bio03, main="Bio1 & Bio3")
points(sites$Bio01, sites$Bio03, cex=1.0, pch=21, bg="red", col="white")
smoothScatter(env.fs$Bio02, env.fs$Bio03, main="Bio2 & Bio3")
points(sites$Bio02, sites$Bio03, cex=1.0, pch=21, bg="red", col="white")
dev.off()

pdf(paste(path1,"FS_OccurrencesAndSitesEnvDatExtracted2.pdf",sep=""), width=7.5*1.2, height=6.15*1.2)
par(mfrow=c(2,3))
plot(bio01.c, main=paste("Bio",bio.nb2[1],sep=""), alpha=150)
points(sites$lonwgs84, sites$latwgs84, cex=1.4, pch=21, bg="red", col="white")
plot(gmted.c, main=paste("GMTED",sep=""), alpha=150)
points(sites$lonwgs84, sites$latwgs84, cex=1.4, pch=21, bg="red", col="white")
plot(pH.top.c, main=paste("pH",sep=""), alpha=150)
points(sites$lonwgs84, sites$latwgs84, cex=1.4, pch=21, bg="red", col="white")
smoothScatter(env.fs$Bio01, env.fs$gmted, main="Bio1 & GMTED")
points(sites$Bio01, sites$gmted, cex=1.0, pch=21, bg="red", col="white")
smoothScatter(env.fs$Bio01, env.fs$pH.top, main="Bio1 & pH")
points(sites$Bio01, sites$pH.top, cex=1.0, pch=21, bg="red", col="white")
smoothScatter(env.fs$gmted, env.fs$pH.top, main="GMTED & pH")
points(sites$gmted, sites$pH.top, cex=1.0, pch=21, bg="red", col="white")
dev.off()

#### Principal component analysis of environmental data at the sampling locations ####
env.sites <- read.table(paste(path1, "FS_SamplingLocationsEnvDatExtracted.csv",sep=""), header=T, sep=";")
pca <- PCA(env.sites[,9:NCOL(env.sites)], scale.unit=T, ncp=5, graph=F)
eig.val <- get_eigenvalue(pca)
var <- get_pca_var(pca)
corrplot(var$cos2, is.corr=F, mar=c(0.5, 0, 3, 0),
         title="Most contributing variables for each dimension")
var.plot1 <- fviz_pca_var(pca, axes=c(1,2), col.var="cos2",
                          title="Variable correlations from PCA axes 1-2",
                          gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"), repel=T)
var.plot2 <- fviz_pca_var(pca, axes=c(1,3), col.var="cos2",
                          title="Variable correlations from PCA axes 1-3",
                          gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"), repel=T)
var.plot3 <- fviz_eig(pca, addlabels=T, ylim=c(0,50),
                      main="Explained variance along PCs")
var.plot4 <- fviz_cos2(pca, choice="var", axes=1:2)
pdf(paste(path1, "FS_SamplingLocationsEnvDatExtracted_PCA_Plots.pdf", sep=""),
    width=8.27, height=8.27)
ggarrange(var.plot1, var.plot2, var.plot3, var.plot4, nrow=2, ncol=2,
          labels=c("A", "B", "C", "D")) ## Export in 10x10inches
dev.off()

#### Correlation plot of environmental data at the sampling locations ####
pdf(paste(path1, "FS_SamplingLocationsEnvDatExtracted_Corr_Plot.pdf", sep=""), width=8.27, height=8.27)
M <- cor(env.sites[,9:NCOL(env.sites)], method="pearson", use="pairwise.complete.obs")  # the correlation matrix of the original data
sigs <- cor.mtest(env.sites[,9:NCOL(env.sites)], conf.level=0.95)  # a test for correlation (if needed)
corrplot(M, type="upper", tl.pos="d", sig.level=0.05, tl.col="black", font=2, tl.cex=0.6,
         tl.offset=0.7, cl.pos="b", mar=c(0.5,0.5,0.5,0.5), mgp=c(0,0,0), oma=c(0,0,0,0))
#corrplot(M, type="upper", tl.pos="d", p.mat=sigs$p, sig.level=0.05, tl.col="black", font=2, tl.cex=0.6, tl.offset=0.7, cl.pos="b", mar=c(0.5,0.5,0.5,0.5), mgp=c(0,0,0), oma=c(0,0,0,0))
corrplot(M, add=T, type="lower", method="number", diag=F, tl.pos="n", cl.pos="n",
         col="black", number.cex=0.75, number.font=1,
         mar=c(0.5,0.5,0.5,0.5), mgp=c(0,0,0), oma=c(0,0,0,0))
dev.off()

#### Pearson correlation matrix of environmental data at the sampling locations ####
corr.pear <- cor(env.sites[,9:NCOL(env.sites)], method="pearson", use="pairwise.complete.obs")
corr.pear
write.table(corr.pear, paste(path1, "FS_SamplingLocationsEnvDatExtracted_Corr_Pearson.csv", sep=""),
            col.names=T, sep=",")

#### Retain uncorrelated variables based on Pearson"s coef < 0.75 and make a PCA ####
uncor.var <- c("Bio01","Bio04","Bio08","Bio12","Bio15","gmted","pH.top")
env.sites.uncor <- env.sites[,uncor.var]
pca <- PCA(env.sites.uncor, scale.unit=T, ncp=5, graph=F)
eig.val <- get_eigenvalue(pca)
var <- get_pca_var(pca)
corrplot(var$cos2, is.corr=F, mar=c(0.5, 0, 3, 0),
         title="Most contributing variables for each dimension")
var.plot1 <- fviz_pca_var(pca, axes=c(1,2), col.var="cos2",
                          title="Variable correlations from PCA axes 1-2",
                          gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"), repel=T)
var.plot2 <- fviz_pca_var(pca, axes=c(1,3), col.var="cos2",
                          title="Variable correlations from PCA axes 1-3",
                          gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"), repel=T)
var.plot3 <- fviz_eig(pca, addlabels=T, ylim=c(0,50), main="Explained variance along PCs")
var.plot4 <- fviz_cos2(pca, choice="var", axes=1:2)
pdf(paste(path1, "FS_SamplingLocationsEnvDatExtracted_PCA_UncorrVar_Plot.pdf", sep=""),
    width=8.27, height=8.27)
ggarrange(var.plot1, var.plot2, var.plot3, var.plot4, nrow=2, ncol=2,
          labels=c("A","B","C","D")) ## Export in 10x10inches
dev.off()

