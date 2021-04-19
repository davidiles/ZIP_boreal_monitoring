#------------------------------------------------
# Load/install packages and set graphical themes / working directory
#------------------------------------------------
my_packs <- c('tidyverse','readxl',        # Data formatting
              'rgeos','raster','sp','sf',  # Spatial packages
              "mgcv","jagsUI","glmmTMB"              # Analysis
) 
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

theme_set(theme_bw())
setwd("~/ECCC Monitoring/Boreal_SK/analysis/scripts")

#------------------------------------------------
# Read landscape rasters that are required for prediction
#------------------------------------------------
covar_rasters <- list(
  AnnTemp = raster("../data/AnnTemp.tif"),
  CMI = raster("../data/CMI.tif"),
  residual_closure = raster("../data/residual.closure.tif"),
  st_std_age = raster("../data/st_std_age.tif"),
  st_std_cc = raster("../data/st_std_cc.tif"),
  st_std_ht = raster("../data/st_std_ht.tif"),
  st_std_nl = raster("../data/st_std_nl.tif"),
  topoWet = raster("../data/topoWet.tif"),
  TREE = raster("../data/TREE.tif")
)


# Initialize the dataframe
covar_matrix <- xyFromCell(covar_rasters[[1]],1:ncell(covar_rasters[[1]]))

# Add columns for each covariate
for (i in 1:length(covar_rasters)) covar_matrix <- cbind(covar_matrix,values(covar_rasters[[i]]))

# Specify column names
for (i in 3:ncol(covar_matrix)) colnames(covar_matrix)[i] <- names(covar_rasters[[i-2]])

# Save in data folder
save(covar_matrix, file = "../data/covar_matrix.RData")

water_mask_raster <- raster("../data/TREE.tif") > 0
writeRaster(water_mask_raster, "../data/water_mask_raster.tif")