#------------------------------------------------
# Load/install packages and set graphical themes / working directory
#------------------------------------------------
my_packs <- c('tidyverse','readxl',        # Data formatting
              'rgeos','raster','sp','sf',  # Spatial packages
              "mgcv","jagsUI","glmmTMB",   # Analysis
              "viridis","MCMCvis","scales"          # Plotting
) 
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

theme_set(theme_bw())
setwd("~/ECCC Monitoring/Boreal_SK/analysis/scripts")

#------------------------------------------------
# Define study extent
#------------------------------------------------
ex_raster <- raster("../data/AnnTemp.tif") # Example raster

# Create polygon from extent of covariate raster
e <- extent(ex_raster)
study_area <- as(e, 'SpatialPolygons')
projection(study_area) <- projection(ex_raster)

#------------------------------------------------
# Model run specs
#------------------------------------------------
focal.species.vector <- c("WTSP","OSFL","OVEN","SAVS","YRWA","HETH","AMRO","DEJU")
model.name.vector <- c("m1","BAM_raster_as_covariate","log_BAM_raster_as_covariate") # Choose a name for this model

results.df <- data.frame()
for (focal.species in focal.species.vector){
  
  BAM_raster_filename <- paste0("../data/BAM_prediction_rasters/pred-",focal.species,"-CAN-Mean.tif")
  sp_raster <- raster(BAM_raster_filename)
  study_area <- st_as_sf(study_area) %>% st_transform(projection(sp_raster)) %>% as('Spatial')
  
  # Crop to study area
  sp_raster <- sp_raster %>%
    crop(study_area) %>%
    mask(study_area)
  
  # Calculate BAM prediction of total abundance
  BAM_total <- sum(values(sp_raster),na.rm = TRUE) * (res(sp_raster)[1]* res(sp_raster)[2])/10000
  
  results.df <- rbind(results.df, data.frame(Species = focal.species,
                                             Model.Type = "BAM Prediction",
                                             Bayesian.p.val = NA,
                                             estimate.mean = BAM_total,
                                             estimate.q025 = NA,
                                             estimate.q975 = NA))
  
  # Extract Bayesian model estimates of total abundance
  for (model.name in model.name.vector){
    
    if (file.exists(paste0("../output/empirical/",focal.species,"_",model.name,"_popsize_mcmc.RData"))){
    load(file = paste0("../output/empirical/",focal.species,"_",model.name,"_out.RData"))
    load(file = paste0("../output/empirical/",focal.species,"_",model.name,"_popsize_mcmc.RData"))
    
    Bayesian.p.val <- mean(out$sims.list$SSE.corrected.observed > out$sims.list$SSE.corrected.sim)
    
    results.df <- rbind(results.df, data.frame(Species = focal.species,
                                               Model.Type = model.name,
                                               Bayesian.p.val = Bayesian.p.val,
                                               estimate.mean = mean(popsize_mcmc),
                                               estimate.q025 = quantile(popsize_mcmc,0.025),
                                               estimate.q975 = quantile(popsize_mcmc,0.975)))
    
    }
  }
}

species.order <- subset(results.df, Model.Type == "BAM Prediction") %>% arrange(estimate.mean)

results.df$BP <- abs(results.df$Bayesian.p.val - 0.5)
results.df$Species <- factor(results.df$Species, levels = species.order$Species)

results.df <- subset(results.df, !(Species == "OVEN" & Model.Type == "log_BAM_raster_as_covariate"))

comparison_plot <- ggplot(results.df, aes(x = Species, col = Model.Type, shape = Model.Type, alpha = BP))+
  geom_errorbar(aes(ymin = estimate.q025, ymax = estimate.q975), position = position_dodge(width = 0.3), width = 0)+
  geom_point(aes(y = estimate.mean), position = position_dodge(width = 0.3))+
  scale_color_manual(values = RColorBrewer::brewer.pal(length(unique(model.name.vector))+1, "Set1"))+
  scale_alpha_continuous(limits = c(0,0.5), range = c(1,.3))+
  
  scale_y_continuous(label = comma)+
  ylab("SK population estimate")+
  theme_bw()
comparison_plot

pdf("../output/empirical/xx_model_comparison.pdf", width = 8, height = 4)
comparison_plot
dev.off()


log_comparison_plot <- ggplot(results.df, aes(x = Species, col = Model.Type, shape = Model.Type, alpha = BP))+
  geom_errorbar(aes(ymin = estimate.q025, ymax = estimate.q975), position = position_dodge(width = 0.25), width = 0)+
  geom_point(aes(y = estimate.mean), position = position_dodge(width = 0.25))+
  scale_color_manual(values = RColorBrewer::brewer.pal(length(unique(model.name.vector))+1, "Set1"))+
  scale_y_continuous(label = comma,trans="log10")+
  scale_alpha_continuous(limits = c(0,0.5), range = c(1,0.3))+
  
  ylab("SK population estimate")+
  theme_bw()
log_comparison_plot

pdf("../output/empirical/xx_model_comparison_logscale.pdf", width = 8, height = 4)
log_comparison_plot
dev.off()