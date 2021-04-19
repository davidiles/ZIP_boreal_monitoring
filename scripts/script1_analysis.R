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
setwd("~/ECCC Monitoring/Boreal_SK/analysis/scripts")

#------------------------------------------------
# Load BOSS data
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

# ***************************************************************
# ***************************************************************
# USER INPUT: Specify species and covariate relationships to be used in analysis below
# ***************************************************************
# ***************************************************************

for (focal.species in c("OSFL","OVEN","WTSP","SAVS","YRWA","HETH","AMRO","DEJU")){
  species.counts <- dat[,focal.species]
  log.offsets <- dat[,paste0(focal.species,".off")]
  
  model.name <- "m1" # Choose a name for this model
  count.model <- formula(site_year ~ st_std_nl + st_std_age + TREE + residual.closure) # site_year as response is just a placeholder
  zi.model <- formula(site_year ~ AnnTemp + I(AnnTemp^2) + CMI)                                       # site_year as response is just a placeholder
  
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
        pred.corrected[i] <- exp( log.mu[site[i]] + offset[i] + 0.5*PSU.sd*PSU.sd) * p.z[site[i]] # Including lognormal correction
        
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
                    count = species.counts,        # Observed counts
                    offset = log.offsets,   # Log-scale offsets
                    PSU = dat$PSU_number,    # PSU identifier
                    site = dat$site_year)    # Survey location identifier
  
  parameters.to.save = c("PSU.sd","log.alpha","logit.beta",
                         "pred","pred.corrected",
                         
                         "SSE.observed",
                         "SSE.sim",
                         
                         "SSE.corrected.observed",
                         "SSE.corrected.sim"
  )
  
  inits <- function()list(z = rep(1,jags.data$nSite))
  
  out <- jags(data = jags.data,
              model.file = "zip.jags",
              parameters.to.save = parameters.to.save,
              inits = inits,
              n.chains = 2,
              n.thin = 20,      # Heavy thinning to minimize storage of autocorrelated mcmc samples
              n.iter = 20000,
              n.burnin = 10000,
              parallel = TRUE)
  
  save(out, file = paste0("../output/empirical/",focal.species,"_",model.name,"_out.RData"))
  out$mcmc.info$elapsed.mins
  max(unlist(out$Rhat),na.rm = TRUE)
  
  #------------------------------------------------
  # Traceplots (evaluate model convergence)
  #------------------------------------------------
  
  pdf(paste0("../output/empirical/",focal.species,"_",model.name,"_traceplot_params.pdf"), width = 7, height = 5)
  MCMCtrace(out, params = c("PSU.sd","log.alpha","logit.beta"),pdf = FALSE)
  dev.off()
  
  #------------------------------------------------
  # Compare sums of predicted and observed counts
  #------------------------------------------------
  
  # Calculate the sum (with 95% credible intervals) of predictions for each observed value
  sum_pred <- out$sims.list$pred %>% apply(., 1, sum)
  sum_pred_corrected <- out$sims.list$pred.corrected %>% apply(., 1, sum)
  
  # Evaluate whether log-normal correction results in correct sum
  pdf(paste0("../output/empirical/",focal.species,"_",model.name,"_SumPredicted_vs_SumObserved.pdf"), width = 5, height = 5)
  limits <- range(c(sum_pred,sum_pred_corrected))
  par(mfrow = c(2,1))
  hist(sum_pred,breaks = seq(limits[1],limits[2],length.out = 50), xlab = "Sum", main = "Predicted Sum\n")
  abline(v = sum(species.counts), lwd = 2, col = "blue")
  
  hist(sum_pred_corrected,breaks = seq(limits[1],limits[2],length.out = 50), xlab = "Sum", main = "Predicted Sum\n(with lognormal correction)")
  abline(v = sum(species.counts), lwd = 2, col = "blue")
  par(mfrow=c(1,1))
  dev.off()
  
  # Confirms that log-normal correction is necessary
  
  #------------------------------------------------
  # Posterior predictive checks
  #------------------------------------------------
  Bayesian.p.value <- mean(out$sims.list$SSE.observed > out$sims.list$SSE.sim)
  pdf(paste0("../output/empirical/",focal.species,"_",model.name,"_PosteriorPredictiveCheck.pdf"), width = 5, height = 5)
  limits <- range(c(out$sims.list$SSE.observed,out$sims.list$SSE.sim))
  plot(out$sims.list$SSE.observed ~ out$sims.list$SSE.sim, pch = 19, xlim = limits, ylim = limits,
       xlab = "Sum Squared Error Simulated", ylab = "Sum Squared Error Observed",
       main = paste0("Bayesian p-value = ", round(Bayesian.p.value,3)))
  abline(a = 0, b = 1)
  dev.off()
  
  
  Bayesian.p.value.corrected <- mean(out$sims.list$SSE.corrected.observed > out$sims.list$SSE.corrected.sim)
  pdf(paste0("../output/empirical/",focal.species,"_",model.name,"_PosteriorPredictiveCheckcorrected.pdf"), width = 5, height = 5)
  limits <- range(c(out$sims.list$SSE.corrected.observed,out$sims.list$SSE.corrected.sim))
  plot(out$sims.list$SSE.corrected.observed ~ out$sims.list$SSE.corrected.sim, pch = 19, xlim = limits, ylim = limits,
       xlab = "Sum Squared Error Simulated", ylab = "Sum Squared Error Observed",
       main = paste0("Bayesian p-value = ", round(Bayesian.p.value.corrected,3),"\n(with lognormal correction)"))
  abline(a = 0, b = 1)
  dev.off()
  
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
  
  #------------------------------------------------
  # Posterior estimate of total population size (takes ~ 1 hour to run)
  #------------------------------------------------
  cell.area.ha <- 3.03 # Average area of each cell (needed to convert density to abundance)
  water_mask <- covar_matrix[,which(colnames(covar_matrix)=="TREE")] > 0 # "Water" cells will have density fixed to zero
  
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
  
  save(popsize_mcmc, file = paste0("../output/empirical/",focal.species,"_",model.name,"_popsize_mcmc.RData"))
  
  mean.pop <- round(mean(popsize_mcmc,na.rm = TRUE))
  lcl.pop <- round(quantile(popsize_mcmc,0.025,na.rm = TRUE))
  ucl.pop <- round(quantile(popsize_mcmc,0.975,na.rm = TRUE)) 
  
  pdf(paste0("../output/empirical/",focal.species,"_",model.name,"_posterior_popsize.pdf"),width = 5,height = 4)
  hist(popsize_mcmc, breaks = 50, main = paste0(mean.pop, "\n(",lcl.pop," to ",ucl.pop,")"),
       xlab = "SK total abundance")
  dev.off()
  
  
  # ----------------------------------------------
  # Mean and uncertainty for each cell in raster (takes ~3.5 hours to run)
  # Note: could write a more efficient loop that "chunks" the data into pieces and uses matrix algebra rather than looping through each cell individually
  # ----------------------------------------------
  lambda.mean <- lambda.se <- lambda.lcl <- lambda.ucl <- rep(NA,nrow(count.X.pred))
  
  start <- Sys.time()
  for (i in 1:nrow(count.X.pred)){
    
    # Predictions for this cell
    log.mu.cell <- (count.X.pred[i,] %*% t(out$sims.list$log.alpha))[1,] + 0.5 * out$sims.list$PSU.sd ^ 2
    logit.p.z.cell <- (zi.X.pred[i,] %*% t(out$sims.list$logit.beta))[1,]
    lambda <- exp(log.mu.cell) * plogis(logit.p.z.cell)
    
    lambda.mean[i] <- mean(lambda)
    lambda.se[i] <- sd(lambda)
    lambda.ci <- quantile(lambda,c(0.025,0.975),na.rm = TRUE)
    lambda.lcl[i] <- lambda.ci[1]
    lambda.ucl[i] <- lambda.ci[2]
    
  }
  
  end <- Sys.time()
  
  lambda_list <- list(lambda.mean = lambda.mean,
                      lambda.se = lambda.se,
                      lambda.lcl = lambda.lcl,
                      lambda.ucl = lambda.ucl)
  
  save(lambda_list, file = paste0("../output/empirical/",focal.species,"_",model.name,"_lambda_list.RData"))
  
  # ----------------------------------------------
  # Generate maps
  # ----------------------------------------------
  
  lambda_mean <- rasterFromXYZ(cbind(covar_matrix[,1:2],lambda_list$lambda.mean*water_mask), res =res(ex_raster), crs = crs(ex_raster))
  pdf(paste0("../output/empirical/",focal.species,"_",model.name,"_posterior_map_lambda_MEAN.pdf"),width = 5,height = 6)
  plot(lambda_mean, col = viridis(10), main = "Predicted males/ha (mean prediction)")
  dev.off()
  
  lambda_se <- rasterFromXYZ(cbind(covar_matrix[,1:2],lambda_list$lambda.se*water_mask), res =res(ex_raster), crs = crs(ex_raster))
  pdf(paste0("../output/empirical/",focal.species,"_",model.name,"_posterior_map_lambda_SE.pdf"),width = 5,height = 6)
  plot(lambda_mean, col = viridis(10), main = "SE of predicted males/ha")
  dev.off()
  
  lambda_CV <- lambda_se / lambda_mean
  pdf(paste0("../output/empirical/",focal.species,"_",model.name,"_posterior_map_lambda_CV.pdf"),width = 5,height = 6)
  plot(lambda_CV, col = viridis(10), main = "CV of predicted males/ha", colNA = "black")
  dev.off()
  
  lambda_lcl <- rasterFromXYZ(cbind(covar_matrix[,1:2],lambda_list$lambda.lcl*water_mask), res =res(ex_raster), crs = crs(ex_raster))
  pdf(paste0("../output/empirical/",focal.species,"_",model.name,"_posterior_map_lambda_LCL.pdf"),width = 5,height = 6)
  plot(lambda_lcl, col = viridis(10), main = "Predicted males/ha (0.025 quantile)")
  dev.off()
  
  lambda_ucl <- rasterFromXYZ(cbind(covar_matrix[,1:2],lambda_list$lambda.ucl*water_mask), res =res(ex_raster), crs = crs(ex_raster))
  pdf(paste0("../output/empirical/",focal.species,"_",model.name,"_posterior_map_lambda_UCL.pdf"),width = 5,height = 6)
  plot(lambda_ucl, col = viridis(10), main = "Predicted males/ha (0.975 quantile)")
  dev.off()
  
  lambda_ciwidth <- lambda_ucl - lambda_lcl
  pdf(paste0("../output/empirical/",focal.species,"_",model.name,"_posterior_map_lambda_CIwidth.pdf"),width = 5,height = 6)
  plot(lambda_ciwidth, col = inferno(10),main = "Width of 95% credible intervals")
  dev.off()
}

