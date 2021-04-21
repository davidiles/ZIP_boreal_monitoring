#------------------------------------------------
# Load/install packages and set graphical themes / working directory
#------------------------------------------------
my_packs <- c('tidyverse','readxl',        # Data formatting
              'rgeos','raster','sp','sf',  # Spatial packages
              "mgcv","jagsUI","glmmTMB",   # Analysis
              "viridis","MCMCvis"          # Plotting
) 
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

theme_set(theme_bw())
setwd("~/ECCC Monitoring/Boreal_SK/analysis/scripts") # MODIFY AS NEEDED

#------------------------------------------------
# Load BOSS data (used for sampling locations and covariate values)
#------------------------------------------------
dat <- read.csv("../data/Zerofilled_Data_with_offsets_incl_HASP.csv")
dat$residual.closure <- residuals(lm(st_std_cc~poly(st_std_ht,deg=2),dat))
dat$site_year <- paste0(dat$UniqueID,"_",dat$Year) %>% as.factor() %>% as.numeric()
dat$PSU_number <- dat$PSU %>% as.factor %>% as.numeric()

#------------------------------------------------
# Summarize site-level covariates
#------------------------------------------------
site_covar <- dat %>% dplyr::select(site_year,st_std_nl,st_std_age,TREE,AnnTemp,CMI,residual.closure) %>% unique()
#z-standardize site covariates
site_covar_Zstand <- site_covar
for (i in 2:ncol(site_covar_Zstand)) site_covar_Zstand[,i] <- scale(site_covar_Zstand[,i])

#------------------------------------------------
# Create extent object
#------------------------------------------------
ex_raster <- raster("../data/AnnTemp.tif") # Example raster (used for extent of study area / projection info)

# Create polygon from extent of covariate raster
e <- extent(ex_raster)
study_area <- as(e, 'SpatialPolygons')
projection(study_area) <- projection(ex_raster)

#------------------------------------------------
# Vector of species to simulate
#------------------------------------------------
species.counts <- sort(colSums(dat[,42:240]),decreasing = TRUE)
species.counts <- species.counts[species.counts > 100]
species.vector <- names(species.counts)

species.vector <- species.vector[-which(species.vector %in% c("HAWO","WTSP"))]

for (focal.species in species.vector){
  
  #------------------------------------------------
  # Load / prepare species prediction raster
  # Downloaded from: https://drive.google.com/drive/folders/1exWa6vfhGo1DNUL4ei2baDz77as7jYzY
  #------------------------------------------------
  
  # Load, crop, and reproject density raster from BAM
  BAM_raster_filename <- paste0("../data/BAM_prediction_rasters/pred-",focal.species,"-CAN-Mean.tif")
  if (!file.exists(BAM_raster_filename)) next
  
  sp_raster <- raster(BAM_raster_filename)
  study_area <- st_as_sf(study_area) %>% st_transform(projection(sp_raster)) %>% as('Spatial')
  
  # Crop to study area
  sp_raster <- sp_raster %>%
    crop(study_area) %>%
    mask(study_area)
  
  # Reproject species distribution raster to same crs/resolution as covariate map
  sp_raster_reproject <- projectRaster(from = sp_raster,to = ex_raster, method = "bilinear")
  
  # ***************************************************************
  # ***************************************************************
  # USER INPUT: Specify covariate relationships to be used in analysis below
  # ***************************************************************
  # ***************************************************************
  model.name <- "m2"
  log.offsets <- dat[,paste0(focal.species,".off")]
  count.model <- formula(site_year ~ st_std_nl + st_std_age + TREE + residual.closure) # site_year as response is just a placeholder
  zi.model <- formula(site_year ~ AnnTemp + I(AnnTemp^2) + CMI)                                       # site_year as response is just a placeholder
  
  #**** NOTE: code will "break" if either model is set to "intercept only".  Need to fix this eventually.
  
  #------------------------------------------------
  # Simulate counts at sampling locations
  #------------------------------------------------
  
  # Create spatial points object
  dat_sp <- dat
  coordinates(dat_sp) <- ~ Longitude + Latitude
  projection(dat_sp) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  # Extract density at each location
  dat$density <- extract(sp_raster_reproject, dat_sp)
  
  # Generate simulated counts
  dat$expected <- dat$density*exp(log.offsets)
  dat$sim.counts <- rpois(nrow(dat),dat$expected)
  
  #------------------------------------------------
  # Create design matrices 
  #------------------------------------------------
  # Design matrix generated for this model
  count.X <- model.matrix(count.model, data = site_covar_Zstand)
  zi.X <- model.matrix(zi.model, data = site_covar_Zstand)
  
  #------------------------------------------------
  # Fit model using JAGS
  #------------------------------------------------
  
  sink("zip.jags")
  cat("
    model {
      
      # ----------------------------
      # Count priors
      # ----------------------------
      
      # Random PSU effects
      PSU.sd ~ dunif(0,5)
      PSU.tau <- pow(PSU.sd,-2)
      for ( i in 1:nPSU ){ PSU.effect[i] ~ dnorm(0,PSU.tau) }
      
      # Log-scale coefficients
      for ( i in 1:k.count ){ log.alpha[i] ~ dnorm(0,0.01) }
      
      # ----------------------------
      # Zero inflation priors
      # ----------------------------
      beta[1] ~ dunif(0,1)
      logit.beta[1] <- log( beta[1]/(1-beta[1]) )
      for ( i in 2:k.zi ){ logit.beta[i] ~ dnorm(0,0.01) }
      
      # ----------------------------
      # Likelihood
      # ----------------------------
      
      # Site-level expected counts and 'occupancy' probabilities
      log.mu   <- count.X %*% log.alpha # Omits PSU-level random effects
      logit.p.z <- zi.X %*% logit.beta
      
      for ( j in 1:nSite ){
        p.z[j] <- 1/( 1+exp(-logit.p.z[j]) )   # Probability of site being 'occupied'
        z[j] ~ dbern( p.z[j] )
        z.sim[j] ~ dbern (p.z[j])
        }
      
      for ( i in 1:nObs ){
        log( lambda[i] ) <- log.mu[site[i]] + offset[i] + PSU.effect[PSU[i]]
        count[i] ~ dpois( lambda[i] * z[site[i]] )
      }
      
      # ----------------------------
      # Predictions for each data point (note that entire landscape predictions are done outside of mcmc engine)
      # Also Posterior Predictive Check
      # ----------------------------
      for ( i in 1:nObs ){
      
        # Predictions / simulated data without lognormal correction
        pred[i] <- exp( log.mu[site[i]] + offset[i]) * p.z[site[i]]
        
        # Predictions / simulated data with lognormal correction
        pred.corrected[i] <- exp( log.mu[site[i]] + offset[i] + 0.5 * PSU.sd * PSU.sd) * p.z[site[i]] # Including lognormal correction
        
      }
      
      # ----------------------------
      # Posterior Predictive Checks
      # ----------------------------
      for (  i in 1:nObs ){
      
        sim.count[i] ~ dpois(exp( log.mu[site[i]] + offset[i]) * z.sim[site[i]])
        sim.count.corrected[i] ~ dpois(exp( log.mu[site[i]] + offset[i] + 0.5 * PSU.sd * PSU.sd) * z.sim[site[i]]) # Including lognormal correction
      
        # Discrepancy measures for actual data
        sqError[i] <- pow(count[i] - pred[i],2)
        sqError.corrected[i] <- pow(count[i] - pred.corrected[i],2)
        
        # Discrepancy measures for simulated
        sqError.sim[i] <- pow(sim.count[i] - pred[i],2)
        sqError.corrected.sim[i] <- pow(sim.count.corrected[i] - pred.corrected[i],2)
        
      }
      
      SSE.observed <- sum(sqError[])
      SSE.sim <- sum(sqError.sim[])
      
      SSE.corrected.observed <- sum(sqError.corrected[])
      SSE.corrected.sim <- sum(sqError.corrected.sim[])
      

    }
    ",fill = TRUE)
  sink()
  
  jags.data <- list(nObs = nrow(dat),           # Number of "survey events"
                    nSite = nrow(site_covar),   # Number of survey locations
                    nPSU = max(dat$PSU_number), # Number of PSUs
                    
                    # Design matrices for count and zero-inflation covariates
                    count.X = count.X,
                    zi.X = zi.X,
                    
                    k.count = ncol(count.X), # Number of parameters for count component
                    k.zi = ncol(zi.X),       # Number of parameters in zero-inflation component
                    
                    # Counts
                    count = dat$sim.counts,        # Observed counts
                    offset = log.offsets,   # Log-scale offsets
                    PSU = dat$PSU_number,    # PSU identifier
                    site = dat$site_year)    # Survey location identifier
  
  parameters.to.save = c("PSU.sd","log.alpha","logit.beta",
                         "pred","pred.corrected",
                         
                         "SSE.observed",
                         "SSE.sim",
                         
                         "SSE.corrected.observed",
                         "SSE.corrected.sim")
  
  inits <- function()list(z = rep(1,jags.data$nSite))
  
  out <- jags(data = jags.data,
              model.file = "zip.jags",
              parameters.to.save = parameters.to.save,
              inits = inits,
              n.chains = 2,
              n.thin = 20,      # Heavy thinning to minimize storage of autocorrelated mcmc samples
              n.iter = 8000,
              n.burnin = 4000,
              parallel = TRUE)
  
  #save(out, file = paste0("../output/simulation/",focal.species,"_",model.name,"_out.RData"))
  out$mcmc.info$elapsed.mins
  max(unlist(out$Rhat),na.rm = TRUE)
  
  #------------------------------------------------
  # Traceplots (evaluate model convergence)
  #------------------------------------------------
  
  # pdf(paste0("../output/simulation/",focal.species,"_",model.name,"_traceplot_params.pdf"), width = 7, height = 5)
  # MCMCtrace(out, params = c("PSU.sd","log.alpha","logit.beta"),pdf = FALSE)
  # dev.off()
  
  #------------------------------------------------
  # Compare sums of predicted and observed counts; confirms that log-normal correction is necessary
  #------------------------------------------------
  
  # # Calculate the sum (with 95% credible intervals) of predictions for each observed value
  # sum_pred <- out$sims.list$pred %>% apply(., 1, sum)
  # sum_pred_corrected <- out$sims.list$pred.corrected %>% apply(., 1, sum)
  # 
  # # Evaluate whether log-normal correction results in correct sum
  # pdf(paste0("../output/simulation/",focal.species,"_",model.name,"_SumPredicted_vs_SumObserved.pdf"), width = 5, height = 5)
  # limits <- range(c(sum_pred,sum_pred_corrected))
  # par(mfrow = c(2,1))
  # hist(sum_pred,breaks = seq(limits[1],limits[2],length.out = 50), xlab = "Sum", main = "Predicted Sum\n")
  # abline(v = sum(species.counts), lwd = 2, col = "blue")
  # 
  # hist(sum_pred_corrected,breaks = seq(limits[1],limits[2],length.out = 50), xlab = "Sum", main = "Predicted Sum\n(with lognormal correction)")
  # abline(v = sum(species.counts), lwd = 2, col = "blue")
  # par(mfrow=c(1,1))
  # dev.off()
  
  #------------------------------------------------
  # Posterior predictive checks
  #------------------------------------------------
  # Bayesian.p.value <- mean(out$sims.list$SSE.observed > out$sims.list$SSE.sim)
  # pdf(paste0("../output/simulation/",focal.species,"_",model.name,"_PosteriorPredictiveCheck.pdf"), width = 5, height = 10)
  # limits <- range(c(out$sims.list$SSE.observed,out$sims.list$SSE.sim,out$sims.list$SSE.corrected.observed,out$sims.list$SSE.corrected.sim))
  # par(mfrow=c(2,1))
  # plot(out$sims.list$SSE.observed ~ out$sims.list$SSE.sim, pch = 19, xlim = limits, ylim = limits,
  #      xlab = "Sum Squared Error Simulated", ylab = "Sum Squared Error Observed",
  #      main = paste0("Bayesian p-value = ", round(Bayesian.p.value,3)))
  # abline(a = 0, b = 1)
  # 
  # Bayesian.p.value.corrected <- mean(out$sims.list$SSE.corrected.observed > out$sims.list$SSE.corrected.sim)
  # plot(out$sims.list$SSE.corrected.observed ~ out$sims.list$SSE.corrected.sim, pch = 19, xlim = limits, ylim = limits,
  #      xlab = "Sum Squared Error Simulated", ylab = "Sum Squared Error Observed",
  #      main = paste0("Bayesian p-value = ", round(Bayesian.p.value.corrected,3),"\n(with lognormal correction)"))
  # abline(a = 0, b = 1)
  # par(mfrow=c(1,1))
  # dev.off()
  
  #------------------------------------------------
  # Prepare province-wide covariate raster / design matrix for prediction
  #------------------------------------------------
  ex_raster <- raster("../data/AnnTemp.tif") # Example raster
  
  # Load covariate matrix
  load("../data/covar_matrix.RData") # this file is created by script0_prepareSKspatialdataframe.R
  
  # Remove columns not used in either count or zi model
  covars.in.analysis <- unique(c(colnames(jags.data$count.X),colnames(jags.data$zi.X)))[-1]
  covar_matrix_2 <- covar_matrix[,which(colnames(covar_matrix) %in% covars.in.analysis)]
  
  # Standardize columns (using mean and SD in empirical dataset; note that mean and SD won't be exactly 0 and 1)
  for (variable in colnames(covar_matrix_2)) covar_matrix_2[,variable] <- (covar_matrix_2[,variable] - mean(site_covar[,variable])) / sd(site_covar[,variable])
  
  # Convert to data.frame
  covar_df <- as.data.frame(covar_matrix_2)
  covar_df$site_year <- 1:nrow(covar_df)
  covar_df[is.na(covar_df)] <- 999 # Temporarily fill NAs (so design matrix has full number of rows - otherwise it removes NAs automatically)
  
  # Create a new design matrix for prediction
  count.X.pred <- model.matrix(count.model, data = covar_df)
  zi.X.pred <- model.matrix(zi.model,  data = covar_df)
  nrow(count.X.pred) == nrow(zi.X.pred) # Ensure this is true
  
  # Re-set the NA values in the design matrix
  count.X.pred[count.X.pred == 999] <- NA
  zi.X.pred[zi.X.pred == 999] <- NA
  water_mask <- covar_matrix[,which(colnames(covar_matrix)=="TREE")] > 0 # "Water" cells will have density fixed to zero
  
  # ----------------------------------------------
  # Compare maps of estimated mean count to actual mean count
  # ----------------------------------------------
  
  log.mu.cell <- (count.X.pred %*% apply(out$sims.list$log.alpha,2,mean))[,1] + 0.5 * mean(out$sims.list$PSU.sd) ^ 2
  logit.p.z.cell <- (zi.X.pred %*% apply(out$sims.list$logit.beta,2,mean))[,1]
  lambda <- exp(log.mu.cell) * plogis(logit.p.z.cell)
  
  lambda_fit_raster <- rasterFromXYZ(cbind(covar_matrix[,1:2],lambda*water_mask), res = res(ex_raster), crs = crs(ex_raster))
  
  bias_raster <- lambda_fit_raster - sp_raster_reproject
  limit <- max(abs(values(bias_raster)),na.rm = TRUE)
  
  pdf(paste0("../output/simulation/",focal.species,"_",model.name,"_Cell_Bias.pdf"), width = 12, height = 5)
  par(mfrow=c(1,3))
  plot(sp_raster_reproject * water_mask, col = viridis(10), colNA = "black", main = "True Density")
  plot(lambda_fit_raster, col = viridis(10), colNA = "black", main = "Estimated Density (mean)")
  plot(bias_raster, col = colorRampPalette(c("red","white","blue"))(11), zlim = c(-limit,limit), main = "Bias",colNA = "black")
  par(mfrow=c(1,1))
  dev.off()
  
  #------------------------------------------------
  # Posterior estimate of total population size (takes ~ 1 hour to run)
  #------------------------------------------------
  cell.area.ha <- (res(sp_raster_reproject)[1] * res(sp_raster_reproject)[2])/10000 # Average area of each cell in ha (needed to convert density to abundance)
  
  # Iterate through mcmc samples and calculate total population size (can't create a single giant matrix of predictions because of RAM limitations)
  nu <- out$mcmc.info$n.samples # By default use all samples for prediction - could use fewer to speed up processing if desired
  samps <- round(seq(1,out$mcmc.info$n.samples,length.out = nu))
  
  popsize_mcmc <- rep(NA,nu)
  
  for (i in 1:nu){
    log.mu <- (count.X.pred %*% out$sims.list$log.alpha[samps[i],])[,1] + 0.5 * out$sims.list$PSU.sd[samps[i]]^2
    logit.p.z <- (zi.X.pred %*% out$sims.list$logit.beta[samps[i],])[,1]
    lambda <- exp(log.mu) * plogis(logit.p.z)
    lambda <- lambda * water_mask
    
    lambda.noNA <- na.omit(lambda)
    popsize_mcmc[i] <- sum(lambda.noNA)*cell.area.ha
    print(i)
    
  }
  
  mean.pop <- round(mean(popsize_mcmc,na.rm = TRUE))
  lcl.pop <- round(quantile(popsize_mcmc,0.025,na.rm = TRUE))
  ucl.pop <- round(quantile(popsize_mcmc,0.975,na.rm = TRUE)) 
  
  #------------------------------------------------
  # Compare *actual* sum abundance to predicted sum abundance
  #------------------------------------------------
  
  actual.sum <- sum(values(sp_raster_reproject),na.rm = TRUE) * cell.area.ha * water_mask
  
  # pdf(paste0("../output/simulation/",focal.species,"_",model.name,"_posterior_popsize.pdf"),width = 5,height = 4)
  # hist(popsize_mcmc, breaks = 50, main = paste0("Estimated Total Abundance:\n",mean.pop, "\n(",lcl.pop," to ",ucl.pop,")"),
  #      xlab = "SK total abundance")
  # abline(v = actual.sum, lwd = 2, col = "blue")
  # text(x = actual.sum, y = 1, pos = 4, label = paste0("Truth = ",round(actual.sum)), col = "blue")
  # dev.off()
  
  
  # ----------------------------------------------
  # Mean and uncertainty for each cell in raster (takes ~3.5 hours to run)
  # Note: could write a more efficient loop that "chunks" the data into pieces and uses matrix algebra rather than looping through each cell individually
  # ----------------------------------------------
  # lambda.mean <- lambda.se <- lambda.lcl <- lambda.ucl <- rep(NA,nrow(count.X.pred))
  # 
  # start <- Sys.time()
  # for (i in 1:nrow(count.X.pred)){
  #   
  #   # Predictions for this cell
  #   log.mu.cell <- (count.X.pred[i,] %*% t(out$sims.list$log.alpha))[1,] + 0.5 * out$sims.list$PSU.sd ^ 2
  #   logit.p.z.cell <- (zi.X.pred[i,] %*% t(out$sims.list$logit.beta))[1,]
  #   lambda <- exp(log.mu.cell) * plogis(logit.p.z.cell)
  #   
  #   lambda.mean[i] <- mean(lambda)
  #   lambda.se[i] <- sd(lambda)
  #   lambda.ci <- quantile(lambda,c(0.025,0.975),na.rm = TRUE)
  #   lambda.lcl[i] <- lambda.ci[1]
  #   lambda.ucl[i] <- lambda.ci[2]
  #   
  # }
  # 
  # end <- Sys.time()
  # 
  # lambda_list <- list(lambda.mean = lambda.mean,
  #                     lambda.se = lambda.se,
  #                     lambda.lcl = lambda.lcl,
  #                     lambda.ucl = lambda.ucl)
  #
  # save(lambda_list, file = paste0("../output/simulation/",focal.species,"_",model.name,"_lambda_list.RData"))
  #
  # # ----------------------------------------------
  # # Generate maps
  # # ----------------------------------------------
  # 
  # lambda_mean <- rasterFromXYZ(cbind(covar_matrix[,1:2],lambda_list$lambda.mean*water_mask), res =res(ex_raster), crs = crs(ex_raster))
  # pdf(paste0("../output/simulation/",focal.species,"_",model.name,"_posterior_map_lambda_MEAN.pdf"),width = 5,height = 6)
  # plot(lambda_mean, col = viridis(10), main = "Predicted males/ha (mean prediction)")
  # dev.off()
  # 
  # lambda_se <- rasterFromXYZ(cbind(covar_matrix[,1:2],lambda_list$lambda.se*water_mask), res =res(ex_raster), crs = crs(ex_raster))
  # pdf(paste0("../output/simulation/",focal.species,"_",model.name,"_posterior_map_lambda_SE.pdf"),width = 5,height = 6)
  # plot(lambda_mean, col = viridis(10), main = "SE of predicted males/ha") 
  # dev.off()
  # 
  # lambda_CV <- lambda_se / lambda_mean
  # pdf(paste0("../output/simulation/",focal.species,"_",model.name,"_posterior_map_lambda_CV.pdf"),width = 5,height = 6)
  # plot(lambda_CV, col = viridis(10), main = "CV of predicted males/ha", colNA = "black") 
  # dev.off()
  # 
  # lambda_lcl <- rasterFromXYZ(cbind(covar_matrix[,1:2],lambda_list$lambda.lcl*water_mask), res =res(ex_raster), crs = crs(ex_raster))
  # pdf(paste0("../output/simulation/",focal.species,"_",model.name,"_posterior_map_lambda_LCL.pdf"),width = 5,height = 6)
  # plot(lambda_lcl, col = viridis(10), main = "Predicted males/ha (0.025 quantile)") 
  # dev.off()
  # 
  # lambda_ucl <- rasterFromXYZ(cbind(covar_matrix[,1:2],lambda_list$lambda.ucl*water_mask), res =res(ex_raster), crs = crs(ex_raster))
  # pdf(paste0("../output/simulation/",focal.species,"_",model.name,"_posterior_map_lambda_UCL.pdf"),width = 5,height = 6)
  # plot(lambda_ucl, col = viridis(10), main = "Predicted males/ha (0.975 quantile)") 
  # dev.off()
  # 
  # lambda_ciwidth <- lambda_ucl - lambda_lcl
  # pdf(paste0("../output/simulation/",focal.species,"_",model.name,"_posterior_map_lambda_CIwidth.pdf"),width = 5,height = 6)
  # plot(lambda_ciwidth, col = inferno(10),main = "Width of 95% credible intervals")
  # dev.off()
  
  #------------------------------------------------
  # Save Relevant Output
  #------------------------------------------------
  output_list <- list(focal.species = focal.species,
                      model.name = model.name,
                      sp_raster_reproject = sp_raster_reproject,
                      cell.area.ha = cell.area.ha,
                      jags.data = jags.data,
                      out = out,
                      popsize_mcmc = popsize_mcmc,
                      actual.sum = actual.sum,
                      lambda_fit_raster = lambda_fit_raster,
                      bias_raster = bias_raster)
  save(output_list,file = paste0("../output/simulation/",focal.species,"_",model.name,"_output_list.RData"))
}

#------------------------------------------------
# Comparison of results for all species with an "output_list"
#------------------------------------------------
dat <- read.csv("../data/Zerofilled_Data_with_offsets_incl_HASP.csv")

water_mask_raster <- raster("../data/water_mask_raster.tif")

SK_pop_results <- data.frame()
for (focal.species in names(sort(colSums(dat[,42:240]),decreasing = TRUE))){
  
  model.name <- "m2"
  filename <- paste0("../output/simulation/",focal.species,"_",model.name,"_output_list.RData")
  if (!file.exists(filename)) next
  
  load(file = filename)
  cell.area.ha <- output_list$cell.area.ha
  max.Rhat <- max(unlist(output_list$out$Rhat)[1:10],na.rm = TRUE)
  prop.Rhat.1.1 <- mean(unlist(output_list$out$Rhat)>1.1,na.rm = TRUE)
  Bayesian.p.value.corrected <- mean(output_list$out$sims.list$SSE.corrected.observed > output_list$out$sims.list$SSE.corrected.sim)
    
  actual.sum <- sum(values(output_list$sp_raster_reproject * water_mask_raster),na.rm = TRUE) * cell.area.ha
  
  
  species_results <- data.frame(Species = focal.species,
                                max.Rhat = max.Rhat,
                                prop.Rhat.1.1 = prop.Rhat.1.1,
                                Bayesian.p.value.corrected = Bayesian.p.value.corrected,
                                True_SK_pop = actual.sum,
                                Est_SK_pop_mean = mean(output_list$popsize_mcmc),
                                Est_SK_pop_q025 = quantile(output_list$popsize_mcmc,0.025),
                                Est_SK_pop_q500 = quantile(output_list$popsize_mcmc,0.500),
                                Est_SK_pop_q975 = quantile(output_list$popsize_mcmc,0.975))
  SK_pop_results <- rbind(SK_pop_results, species_results )
}


mx <- max(SK_pop_results[,c(4:8)])*1.2
SK_total_plot <- ggplot(SK_pop_results) +
  geom_abline(intercept = 0, slope = 1)+
  
  geom_errorbar(aes(x = True_SK_pop, ymin = Est_SK_pop_q025, ymax = Est_SK_pop_q975, col = Bayesian.p.value.corrected, alpha = max.Rhat), width = 0)+
  geom_point(aes(x = True_SK_pop, y = Est_SK_pop_mean, col = Bayesian.p.value.corrected, alpha = max.Rhat))+
  geom_text(aes(x = True_SK_pop, y = Est_SK_pop_mean, label = Species, col = Bayesian.p.value.corrected, alpha = max.Rhat), hjust = 0)+

  scale_color_gradientn(colors = c("dodgerblue","darkblue","black","darkred","orangered"), limits = c(0,1),name = "Bayesian P value\n(Goodness-of-fit)")+
  coord_cartesian(xlim = c(0,mx), ylim = c(0,mx))+
  xlab("True Pop Size")+
  ylab("Estimated Pop Size")+
  theme_bw()
SK_total_plot

pdf("../output/simulation/SK_species_simulation.pdf", width = 7, height = 5)
SK_total_plot
dev.off()
