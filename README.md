## R Code--------------------------------------------------------------------
## Title: Analysis of Negative Power Transformation-GEV Distribution (NPT-GEVD) for Estimation of Minimum IAT Data
## Authors: Prahadchai T., Yoon S., and Busababodhin P.
## Version 1 (1 Dec 2024)
## --------------------------------------------------------------------------
## Short Description:
## There are a total of 8steps to obtain the parameter and return level (RL) estimates.
## Please follow steps 1)–8). Note that steps 3)–7) require setting the directory, providing minimum data, and obtaining the results.
## Step 8) involves saving the final output in Excel format and storing it in your folder directory.
## --------------------------------------------------------------------------

## 1) Packages ---
library(SpatialExtremes)
library(ismev)
library(extRemes)
library(evd)
library(lmom)
library(goftest)
library(dplyr)
library(stringr)
library(openxlsx)

## 2) Important functions --------

### Calculation function of the RL for CT-GEVD using the MLE method and the standard error (SE) of RL using the delta method.
getReturnLevels_ct <- function(f, rp=c(25,50,100), conf=0.95){
  #### make empty list
  dn <- list(NULL,c('period','2.5% ','Estimate','97.5%','se'))
  rl <- matrix(NA,nrow=length(rp),ncol=5,dimnames=dn)
  rl[,1] <- rp
  
  #### Extract the estimated parameters (location, scale, shape)
  params <- f$mle
  names(params) <- c("location", "scale", "shape")
  
  ##### Manually compute the covariance matrix using the Fisher Information
  Cov_mat <- f$cov
  if(any(is.na(Cov_mat))){
    resam_x <- sample(f$data, size=(length(f$data)-5), replace = TRUE)
    refit_x <- gev.fit(resam_x, show=FALSE)
    cov_matrix <- refit_x$cov
  }else{
    cov_matrix <- f$cov
  }
  
  ##### Define the return period and calculate the return level for the transformed data
  return_level_transformed <- params["location"] - (params["scale"] / params["shape"]) * 
    (1 - (-log(1 - 1/rp))^(-params["shape"]))
  
  ##### Transform the return level back to the original scale by multiplying by -1
  ##### return_level_minima <- -return_level_transformed
  rl[,3] <-  -return_level_transformed
  
  ##### Delta method: Calculate the gradient of the rl with respect to the parameters (for maxima data)
  gradient <- NULL
  for (l in 1:length(rp)) {
    location_grad <- 1
    scale_grad <- -((1 - (-log(1 - 1/rp[l]))^(-params["shape"])) / params["shape"])
    shape_grad <- (params["scale"] / params["shape"]^2) * 
      (1 - (-log(1 - 1/rp[l]))^(-params["shape"])) - 
      (params["scale"] / params["shape"]) * 
      (-log(1 - 1/rp[l]))^(-params["shape"]) * log(-log(1 - 1/rp[l]))
    
    # Combine the gradients into a numeric vector
    gradient <- rbind(gradient, c(location_grad, scale_grad, shape_grad)) 
  }
  
  ##### Compute the standard error of the return level for the transformed data (maxima)
  v <- gradient%*%cov_matrix%*%t(gradient)
  
  ##### Conf with rl
  rl[,2] <- rl[,3]-qnorm(1-(1-conf)/2)*sqrt(diag(v))
  rl[,4] <- rl[,3]+qnorm(1-(1-conf)/2)*sqrt(diag(v))
  rl[,5] <- sqrt(diag(v))
  
  return(rl)
}

#### Calculation function of the RL for RT-GEVD using the MLE method and the SE of RL using the delta method.
getReturnLevels_rt <- function(f, rp=c(25,50,100), conf=0.95, sd_scale=NULL, trans=TRUE){
  # make empty list
  dn <- list(NULL,c('period','2.5% ','Estimate','97.5%','se'))
  rl <- matrix(NA,nrow=length(rp),ncol=5,dimnames=dn)
  rl[,1] <- rp
  
  # Extract the estimated parameters (location, scale, shape)
  params <- f$mle
  names(params) <- c("location", "scale", "shape")
  
  # Manually compute the covariance matrix using the Fisher Information
  Cov_mat <- f$cov
  if(any(is.na(Cov_mat))){
    resam_x <- sample(f$data, replace = TRUE)
    refit_x <- gev.fit(resam_x, show=FALSE)
    cov_matrix <- refit_x$cov
  }else{
    cov_matrix <- f$cov
  }
  
  # Define the return period and calculate the return level for the transformed data
  rl_transformed <- 
    params["location"] - (params["scale"] / params["shape"]) * (1 - (-log(1 - 1/rp))^(-params["shape"]))
  
  if(trans){
    rl_transformed_scaled <- (rl_transformed^(-1) - 1) * (sd(f$min.data) / sd_scale) + mean(f$min.data)
    rl[,3] <-  rl_transformed_scaled
  }else{
    rl[,3] <-  rl_transformed^(-1)
  }
  
  # Delta method: Calculate the gradient of the rl with respect to the parameters (for maxima data)
  gradient <- NULL
  for (l in 1:length(rp)) {
    location_grad <- 1
    scale_grad <- -((1 - (-log(1 - 1/rp[l]))^(-params["shape"])) / params["shape"])
    shape_grad <- (params["scale"] / params["shape"]^2) * 
      (1 - (-log(1 - 1/rp[l]))^(-params["shape"])) - 
      (params["scale"] / params["shape"]) * 
      (-log(1 - 1/rp[l]))^(-params["shape"]) * log(-log(1 - 1/rp[l]))
    
    # Combine the gradients into a numeric vector
    gradient <- rbind(gradient, c(location_grad, scale_grad, shape_grad)) 
  }
  
  # Compute the standard error of the return level for the transformed data (maxima)
  v <- gradient%*%cov_matrix%*%t(gradient)
  
  # Conf with rl
  rl[,2] <- rl[,3]-qnorm(1-(1-conf)/2)*sqrt(diag(v))
  rl[,4] <- rl[,3]+qnorm(1-(1-conf)/2)*sqrt(diag(v))
  rl[,5] <- sqrt(diag(v))
  
  return(rl)
}

# Calculation function of the RL for NPT-GEVD using the MLE method and the SE of RL using the delta method.
getReturnLevels_npt <- function(f, a, rp=c(25,50,100), conf=0.95, sd_scale=NULL, trans=TRUE){
  # make empty list
  dn <- list(NULL,c('period','2.5% ','Estimate','97.5%','se'))
  rl <- matrix(NA,nrow=length(rp),ncol=5,dimnames=dn)
  rl[,1] <- rp
  
  # Extract the estimated parameters (location, scale, shape)
  params <- f$mle
  names(params) <- c("location", "scale", "shape")
  
  # Manually compute the covariance matrix using the Fisher Information
  Cov_mat <- f$cov
  if(any(is.na(Cov_mat))){
    resam_x <- sample(f$data, replace = TRUE)
    refit_x <- gev.fit(resam_x, show=FALSE)
    cov_matrix <- refit_x$cov
  }else{
    cov_matrix <- f$cov
  }
  
  # Define the return period and calculate the return level for the transformed data
  rl_transformed <- 
    params["location"] - (params["scale"] / params["shape"]) * (1 - (-log(1 - 1/rp))^(-params["shape"]))
  
  
  if(trans){
    rl_transformed_scaled <- (rl_transformed^(-1/a) - 1) * (sd(f$min.data) / sd_scale) + mean(f$min.data)
    rl[,3] <-  rl_transformed_scaled
  }else{
    rl[,3] <-  rl_transformed^(-1/a)
  }
  
  # Delta method: Calculate the gradient of the rl with respect to the parameters (for maxima data)
  gradient <- NULL
  for (l in 1:length(rp)) {
    location_grad <- 1
    scale_grad <- -((1 - (-log(1 - 1/rp[l]))^(-params["shape"])) / params["shape"])
    shape_grad <- (params["scale"] / params["shape"]^2) * 
      (1 - (-log(1 - 1/rp[l]))^(-params["shape"])) - 
      (params["scale"] / params["shape"]) * 
      (-log(1 - 1/rp[l]))^(-params["shape"]) * log(-log(1 - 1/rp[l]))
    
    # Combine the gradients into a numeric vector
    gradient <- rbind(gradient, c(location_grad, scale_grad, shape_grad)) 
  }
  
  # Compute the standard error of the return level for the transformed data (maxima)
  v <- gradient%*%cov_matrix%*%t(gradient)
  
  # Conf with rl
  rl[,2] <- rl[,3]-qnorm(1-(1-conf)/2)*sqrt(diag(v))
  rl[,4] <- rl[,3]+qnorm(1-(1-conf)/2)*sqrt(diag(v))
  rl[,5] <- sqrt(diag(v))
  
  return(rl)
}

# Bootstrap sampling technique for refining the SE of parameter estimates in cases where the Fisher information matrix returns NA values.
boots.lmom_par <- function(data, nboot=1000, method="l-moment"){
  # Initialize a matrix to store bootstrap parameters
  boot_params <- matrix(NA, nrow = nboot, ncol = 3)
  
  # Bootstrap loop
  for (i in 1:nboot) {
    sample_data <- sample_fit  <- NULL
    # Resample data with replacement
    sample_data <- sample(data, replace = TRUE, size = length(data))
    
    # Refit distribution
    if(method=="l-moment"){
      sample_fit <- tryCatch({
        fevd(sample_data, method="Lmoments", type="GEV")$results
      }, error = function(e) rep(NA, 3))
    }else if(method=="mle"){
      sample_fit <- tryCatch({
        gev.fit(sample_data, show=FALSE)$mle
      }, error = function(e) rep(NA, 3))
    }
    # Store fitted parameters
    boot_params[i, ] <- sample_fit
  }
  # Calculate standard errors for each parameter
  param_se <- apply(boot_params, 2, sd, na.rm = TRUE)
  return(param_se)
}

# Calculation function of the RL for CT-GEVD using the LM method and the SE of RL using the bootstrap sampling method.
boots.lmom_rl_ct <- function(data, nboot=1000, rp=c(25,50,100), conf=0.95, method0="Lmoments"){
  # Set results matrix
  dn <- list(NULL,c('period','2.5% ','Estimate','97.5%','se'))
  rl <- matrix(NA,nrow=length(rp),ncol=5,dimnames=dn)
  rl[,1] <- rp
  
  # True RL
  true_fit <- fevd(data, method=method0, type="GEV")
  true_rl <- tryCatch({ -return.level(true_fit, return.period=rp) %>% as.numeric() 
  }, error = function(e) rep(NA, length(rp)))
  
  # Initialize a matrix to store bootstrap parameters
  boot_params <- matrix(NA, nrow = nboot, ncol = length(rp))
  
  # Bootstrap loop
  for (i in 1:nboot) { 
    sample_data <- sample_fit <- sample_rl <- NULL
    # Resample data with replacement
    sample_data <- sample(data, replace = TRUE, size = length(data))
    
    # Refit distribution
    sample_fit <- fevd(sample_data, method=method0, type="GEV")
    
    # Recal RL
    sample_rl <- tryCatch({ -return.level(sample_fit, return.period=rp) %>% as.numeric() 
    }, error = function(e) rep(NA, length(rp)))
    
    # Store fitted parameters
    boot_params[i, ] <- sample_rl
  }
  # Calculate standard errors for each parameter
  param_se <- apply(boot_params, 2, sd, na.rm = TRUE)
  rl[,3] <- true_rl
  rl[,2] <- rl[,3]-qnorm(1-(1-conf)/2)*param_se
  rl[,4] <- rl[,3]+qnorm(1-(1-conf)/2)*param_se
  rl[,5] <- param_se
  return(rl)
}

# Calculation function of the RL for RT-GEVD using the LM method and the SE of RL using the bootstrap sampling method.
boots.lmom_rl_rt <- function(data, min.data=NULL, nboot=1000, rp=c(25,50,100), conf=0.95, sd_scale=NULL, method0="Lmoments", trans=TRUE){
  # Set results matrix
  dn <- list(NULL,c('period','2.5% ','Estimate','97.5%','se'))
  rl <- matrix(NA,nrow=length(rp),ncol=5,dimnames=dn)
  rl[,1] <- rp
  
  if(trans){
    # True RL
    true_fit <- fevd(data, method=method0, type="GEV")
    true_rl <- tryCatch({ 
      ((return.level(true_fit, return.period=rp)^(-1) - 1) * (sd(min.data) / sd_scale) + mean(min.data)) %>% as.numeric() 
    }, error = function(e) rep(NA, length(rp)))
    
    # Initialize a matrix to store bootstrap parameters
    boot_params <- matrix(NA, nrow = nboot, ncol = length(rp))
    
    # Bootstrap loop
    for (i in 1:nboot) { 
      sample_data <- sample_fit <- sample_rl <- NULL
      # Resample data with replacement
      sample_data <- sample(data, replace = TRUE, size = length(data))
      
      # Refit distribution
      sample_fit <- fevd(sample_data, method=method0, type="GEV")
      
      # Recal RL
      sample_rl <- tryCatch({ 
        ((return.level(sample_fit, return.period=rp)^(-1) - 1) * (sd(min.data) / sd_scale) + mean(min.data)) %>% as.numeric() 
      }, error = function(e) rep(NA, length(rp)))
      
      # Store fitted parameters
      boot_params[i, ] <- sample_rl
    }
    
  }else{
    # True RL
    true_fit <- fevd(data, method=method0, type="GEV")
    true_rl <- tryCatch({ return.level(true_fit, return.period=rp)^(-1) %>% as.numeric() 
    }, error = function(e) rep(NA, length(rp)))
    
    # Initialize a matrix to store bootstrap parameters
    boot_params <- matrix(NA, nrow = nboot, ncol = length(rp))
    
    # Bootstrap loop
    for (i in 1:nboot) { 
      sample_data <- sample_fit <- sample_rl <- NULL
      # Resample data with replacement
      sample_data <- sample(data, replace = TRUE, size = length(data))
      
      # Refit distribution
      sample_fit <- fevd(sample_data, method=method0, type="GEV")
      
      # Recal RL
      sample_rl <- tryCatch({ return.level(sample_fit, return.period=rp)^(-1) %>% as.numeric() 
      }, error = function(e) rep(NA, length(rp)))
      
      # Store fitted parameters
      boot_params[i, ] <- sample_rl
    }
  }
  
  # Calculate standard errors for each parameter
  param_se <- apply(boot_params, 2, sd, na.rm = TRUE)
  rl[,3] <- true_rl
  rl[,2] <- rl[,3]-qnorm(1-(1-conf)/2)*param_se
  rl[,4] <- rl[,3]+qnorm(1-(1-conf)/2)*param_se
  rl[,5] <- param_se
  return(rl)
}

# Calculation function of the RL for NPT-GEVD using the LM method and the SE of RL using the bootstrap sampling method.
boots.lmom_rl_npt <- function(data, a, min.data=NULL, nboot=1000, rp=c(25,50,100), conf=0.95, sd_scale=NULL, method0="Lmoments", trans=TRUE){
  # Set results matrix
  dn <- list(NULL,c('period','2.5% ','Estimate','97.5%','se'))
  rl <- matrix(NA,nrow=length(rp),ncol=5,dimnames=dn)
  rl[,1] <- rp
  
  if(trans){
    # True RL
    true_fit <- fevd(data, method=method0, type="GEV")
    true_rl <- tryCatch({ 
      ((return.level(true_fit, return.period=rp)^(-1/a) - 1) * (sd(min.data) / sd_scale) + mean(min.data)) %>% as.numeric() 
    }, error = function(e) rep(NA, length(rp)))
    
    # Initialize a matrix to store bootstrap parameters
    boot_params <- matrix(NA, nrow = nboot, ncol = length(rp))
    
    # Bootstrap loop
    for (i in 1:nboot) { 
      sample_data <- sample_fit <- sample_rl <- NULL
      # Resample data with replacement
      sample_data <- sample(data, replace = TRUE, size = length(data))
      
      # Refit distribution
      sample_fit <- fevd(sample_data, method=method0, type="GEV")
      
      # Recal RL
      sample_rl <- tryCatch({ 
        ((return.level(sample_fit, return.period=rp)^(-1/a) - 1) * (sd(min.data) / sd_scale) + mean(min.data)) %>% as.numeric() 
      }, error = function(e) rep(NA, length(rp)))
      
      # Store fitted parameters
      boot_params[i, ] <- sample_rl
    }
    
  }else{
    # True RL
    true_fit <- fevd(data, method=method0, type="GEV")
    true_rl <- tryCatch({ return.level(true_fit, return.period=rp)^(-1/a) %>% as.numeric() 
    }, error = function(e) rep(NA, length(rp)))
    
    # Initialize a matrix to store bootstrap parameters
    boot_params <- matrix(NA, nrow = nboot, ncol = length(rp))
    
    # Bootstrap loop
    for (i in 1:nboot) { 
      sample_data <- sample_fit <- sample_rl <- NULL
      # Resample data with replacement
      sample_data <- sample(data, replace = TRUE, size = length(data))
      
      # Refit distribution
      sample_fit <- fevd(sample_data, method=method0, type="GEV")
      
      # Recal RL
      sample_rl <- tryCatch({ return.level(sample_fit, return.period=rp)^(-1/a) %>% as.numeric() 
      }, error = function(e) rep(NA, length(rp)))
      
      # Store fitted parameters
      boot_params[i, ] <- sample_rl
    }
  }
  
  # Calculate standard errors for each parameter
  param_se <- apply(boot_params, 2, sd, na.rm = TRUE)
  rl[,3] <- true_rl
  rl[,2] <- rl[,3]-qnorm(1-(1-conf)/2)*param_se
  rl[,4] <- rl[,3]+qnorm(1-(1-conf)/2)*param_se
  rl[,5] <- param_se
  return(rl)
}

# Optimization function based on the L-BFGS-B method by minimizing the Cramér-von Mises (CvM) statistical test with LM method
cal.f.cvm.lm <- function(a, sam.min.data, scale_min.data){
  # transform data
  trans.data <- sam.min.data^(-a)
  
  # gev distribution fit
  gev.data <- fevd(trans.data, method="Lmoments", type="GEV")
  params <- gev.data$results %>% as.numeric()
  
  # test gof
  gof.test <- goftest::cvm.test(scale_min.data,'pgev', params[1], params[2], params[3])
  return(gof.test$statistic)
}

# Optimization function based on the L-BFGS-B method by minimizing the CvM statistical test with MLE method
cal.f.cvm.mle <- function(a, sam.min.data, scale_min.data){
  # transform data
  trans.data <- sam.min.data^(-a)
  
  # gev distribution fit
  # gev.data <- fevd(trans.data, method="MLE", type="GEV")
  # params <- gev.data$results$par %>% as.numeric()
  params <- gev.fit(trans.data, show=FALSE)$mle
  
  # test gof
  gof.test <- goftest::cvm.test(scale_min.data,'pgev', params[1], params[2], params[3])
  return(gof.test$statistic)
}

# Hyperparameter estimation function using resampling with an added small noise technique.
opt_a0_addnoise <- function(target_func=NULL, scale_min.data=NULL, a0=NULL, n0=NULL){
  opt.a <- NULL
  rep_n <- 50
  rep_max <- 100
  r.vec <- NULL
  j <- 0
  count <- 0
  while ((j < rep_n) && (count < rep_max)) {
    # red_n <- n0*per_n0
    resam.min.data <- scale_min.data + rnorm(n0, mean=0, sd=0.01)
    op <- NULL
    
    try({op <- optim(par=1, fn=get(paste0(target_func)), sam.min.data=resam.min.data, scale_min.data=scale_min.data,
                     lower = 0, upper = 100, method="L-BFGS-B"
    )}, silent = TRUE)
    
    if(!is.null(op) && (op$par>0) && (round(op$par,4)!=a0) &&
       (round(op$par,4)!=0) && (round(op$par,4)!=100)){
      j <- j+1
      r.vec <- c(r.vec, op$par)
      cat("rep::",j,"est.a::",op$par,"\n")
    }
    count <- count + 1
  }
  # attention here!!!!
  if(!is.null(r.vec)){
    if(length(r.vec)>2){
      opt.a = median(r.vec)
    }
  }
  return(opt.a)
}


# 3) Setting path --------
setwd("/Users/thanawanp/Desktop/R Code - NPT_GEVD/IAT of hourly rainfall data")

# 4) Read a list of files in a folder --------
f.list <- list.files()

# 5) The main function estimates parameters and return levels (RL) --------
# 'min.data' is the input for minimum IAT data.
# 'rp0' specifies the RL estimate corresponding to 25-, 50-, and 100-year return periods.
# 'conf0' represents the 95% confidence interval for parameter estimates and RL.
# 'trans0' is require if you want to rescale data to mu->1 and sig->0.1

cal.g <- function(min.data=NULL, rp0=c(25,50,100), conf0=0.95, trans0=TRUE){ 
  # Initialize an empty output to store results. 
  # Vector 'A' is for storing output from CT-GEVD
  A1 <- A1_par <- A1_se <- A1_nllh <- A1_rl <- A1_cvm <- A1_ad <- NULL
  A2 <- A2_par <- A2_se <- A2_rl <- A2_cvm <- A2_ad <- NULL
  # Vector 'B' is for storing output from RT-GEVD
  B1 <- B1_par <- B1_se <- B1_nllh <- B1_rl <- B1_cvm <- B1_ad <- NULL
  B2 <- B2_par <- B2_se <- B2_rl <- B2_cvm <- B2_ad <- NULL
  # Vector 'C' is for storing output from NPT-GEVD
  C1 <- C1_par <- C1_se <- C1_nllh <- C1_rl <- C1_rl2 <- C1_cvm <- C1_ad <- NULL
  C2 <- C2_par <- C2_se <- C2_rl <- C2_rl2 <- C2_cvm <- C2_ad <- NULL
  # Empty vector for hyper-parameter estimates through MLE and LM method
  mle_a <- lm_a <- NULL
  
  # CT-GEVD ----
  m1.data = -min.data
  # MLE method ---
  # Parameters estimate
  A1 <- gev.fit(m1.data, show=FALSE)
  A1_cov <- A1$cov
  if(any(is.na(A1_cov)) || any(abs(A1_cov) < 1e-04)){ # Checking the Hessian matrix to verify if it works.
    A1_se <- boots.lmom_par(data=m1.data, method="mle") # Using the bootstrap sampling technique to obtain the standard error
  }else{
    A1_se <- A1$se
  }
  A1_par <- data.frame(est=A1$mle, se=A1_se)
  A1_nllh<- A1$nllh
  # RLs estimate
  A1_rl  <- getReturnLevels_ct(f=A1, rp=rp0, conf=conf0) # SE of RL calculated by delta method.
  if(any(is.na(A1_rl[,2]))){ # If the delta method does not work, the bootstrap sampling technique can be applied.
    A1_rl <- NULL
    A1_rl <- boots.lmom_rl_ct(data=m1.data, rp=rp0, conf=conf0, method0="MLE")
  }
  A1_cvm <- cvm.test(A1$data, "pgev", A1$mle[1], A1$mle[2], A1$mle[3])
  A1_ad  <- ad.test(A1$data, "pgev", A1$mle[1], A1$mle[2], A1$mle[3])
  
  # LM method ---
  A2 <- fevd(x=m1.data, method="Lmoments",type="GEV")
  A2_se <- boots.lmom_par(data=m1.data)
  A2_par<- data.frame(est=as.numeric(A2$results), se=A2_se)
  A2_rl <- boots.lmom_rl_ct(data=m1.data, rp=rp0, conf=conf0)
  A2_cvm<- cvm.test(A2$x, "pgev", A2$results[1], A2$results[2], A2$results[3])
  A2_ad <- ad.test(A2$x, "pgev", A2$results[1], A2$results[2], A2$results[3])
  
  # RT-GEVD ----
  if(trans0){ # if 'trans0' is TRUE
    # min.data need to be rescaled with mu->1, sig->0.1 for RT-GEVD and NPT-GEVD.
    sd_scale <- 0.1
    scale_min.data <- (min.data - mean(min.data)) * (sd_scale / sd(min.data)) + 1
    m2.data <- scale_min.data^(-1)
  }else{
    m2.data <- min.data^(-1)
  }
  
  # MLE method ---
  # Parameters estimate
  B1 <- gev.fit(m2.data, show=FALSE)
  B1_cov <- B1$cov
  if(any(is.na(B1_cov)) || any(abs(B1_cov) < 1e-04)){
    B1_se <- boots.lmom_par(data=m2.data, method="mle")
  }else{
    B1_se <- B1$se
  }
  B1_par <- data.frame(est=B1$mle, se=B1_se)
  B1_nllh<- B1$nllh
  B1$min.data <- min.data
  # RLs estimate
  B1_rl <- getReturnLevels_rt(f=B1, rp=rp0, conf=conf0, sd_scale=sd_scale, trans=trans0) # SE of RL calculated by delta method.
  if(any(is.na(B1_rl[,2]))){ # If the delta method does not work, the bootstrap sampling technique can be applied.
    B1_rl <- NULL
    B1_rl <- boots.lmom_rl_rt(data=m2.data, min.data=min.data, rp=rp0, conf=conf0, sd_scale=sd_scale, method0 = "MLE", trans=trans0)
  }
  B1_cvm <- cvm.test(B1$data, "pgev", B1$mle[1], B1$mle[2], B1$mle[3])
  B1_ad  <- ad.test(B1$data, "pgev", B1$mle[1], B1$mle[2], B1$mle[3])
  
  # LM method ---
  B2 <- fevd(x=m2.data, method="Lmoments",type="GEV")
  B2_se <- boots.lmom_par(data=m2.data, method="l-moment")
  B2_par<- data.frame(est=as.numeric(B2$results), se=B2_se)
  B2_rl <- boots.lmom_rl_rt(data=m2.data, min.data=min.data, rp=rp0, conf=conf0, sd_scale=sd_scale, trans=trans0)
  B2_cvm<- cvm.test(B2$x, "pgev", B2$results[1], B2$results[2], B2$results[3])
  B2_ad <- ad.test(B2$x, "pgev", B2$results[1], B2$results[2], B2$results[3])
  
  # NPT-GEVD ---
  # Finding the optimal hyper-parameter using an optimization algorithm based on the L-BFGS-B method.
  if(trans0){
    scale_min.data <- (min.data - mean(min.data)) * (sd_scale / sd(min.data)) + 1
  }else{
    scale_min.data <- min.data
  }
  # The output is the median of the estimates of 50 hyperparameters obtained through the resampling technique and MLE method.
  mle_a <- opt_a0_addnoise(target_func="cal.f.cvm.mle", scale_min.data=scale_min.data, a0=1, n0=length(scale_min.data))
  
  # Checking the validation of hyperparameter estimate is returned as output.
  if(is.null(mle_a)){ # If TRUE, then change the estimation method to LM
    mle_a_using_lm <- TRUE
    mle_a <- opt_a0_addnoise(target_func="cal.f.cvm.lm",  scale_min.data=scale_min.data, a0=1, n0=length(scale_min.data))
    
    recal_mle_a <- NULL
    recal_mle_a <- opt_a0_addnoise(target_func="cal.f.cvm.lm",  scale_min.data=scale_min.data, a0=1, n0=length(scale_min.data)) 
    if(is.null(recal_mle_a)){ # If the hyperparameter estimate is still not valid, then resample the data by reducing the sample size once.
      recal_mle_a2 <- NULL
      recal_mle_a2 <- opt_a0_addnoise(target_func="cal.f.cvm.lm",
                                      scale_min.data=sample(scale_min.data, size=(length(scale_min.data)*0.8), replace = TRUE), 
                                      a0=1, n0=length(scale_min.data)) 
      if(is.null(recal_mle_a2)){
        mle_a = 1
      }else{
        mle_a = recal_mle_a2
      }
    }else{
      mle_a = recal_mle_a
    }
  }else{
    mle_a_using_lm <- FALSE
  }
  
  # Transform minimum data to maximum.
  m3.data.mle <- scale_min.data^(-mle_a)
  
  # MLE method ---
  # Parameter estimate
  C1 <- gev.fit(m3.data.mle, show=FALSE)
  C1_cov <- C1$cov
  if(any(is.na(C1_cov)) || any(abs(C1_cov) < 1e-04)){
    C1_se <- boots.lmom_par(data=m3.data.mle, method="mle")
  }else{
    C1_se <- C1$se
  }
  C1_par <- data.frame(est=C1$mle, se=C1_se)
  C1_nllh<- C1$nllh
  C1$min.data <- min.data
  # RLs estimate
  C1_rl <- getReturnLevels_npt(f=C1, a=mle_a, rp=rp0, conf=conf0, sd_scale=sd_scale, trans=trans0) # SE of RL calculated by delta method.
  if(any(is.na(C1_rl[,2]))){ # If the delta method does not work, the bootstrap sampling technique can be applied.
    C1_rl <- NULL
    C1_rl <- boots.lmom_rl_npt(data=m3.data.mle, a=mle_a, min.data=min.data, rp=rp0, conf=conf0, sd_scale=sd_scale, method0 = "MLE", trans=trans0)
  }
  C1_cvm <- cvm.test(C1$data, "pgev", C1$mle[1], C1$mle[2], C1$mle[3])
  C1_ad  <- ad.test(C1$data, "pgev", C1$mle[1], C1$mle[2], C1$mle[3])
  
  # LM method ---
  # The output is the median of the estimates of 50 hyperparameters obtained through the resampling technique and LM method.
  # And checking the validation of hyperparameter estimate is returned as output.
  if(mle_a_using_lm){
    lm_a <- mle_a
  }else{
    recal_lm_a <- NULL
    recal_lm_a <- opt_a0_addnoise(target_func="cal.f.cvm.lm", scale_min.data=scale_min.data, a0=1, n0=length(scale_min.data))
    if(is.null(recal_lm_a)){
      recal_lm_a2 <- NULL
      recal_lm_a2 <- opt_a0_addnoise(target_func="cal.f.cvm.lm",
                                     scale_min.data=sample(scale_min.data, size=(length(scale_min.data)*0.8), replace = TRUE), 
                                     a0=1, n0=length(scale_min.data)) 
      if(is.null(recal_lm_a2)){
        lm_a = 1
      }else{
        lm_a = recal_lm_a2
      }
    }else{
      lm_a = recal_lm_a
    }
  }
  
  # Transform minimum data to maximum.
  m3.data.lm <- scale_min.data^(-lm_a)
  
  C2 <- fevd(m3.data.lm, method="Lmoments", type="GEV")
  C2_se <- boots.lmom_par(data=m3.data.lm, method="l-moment")
  C2_par<- data.frame(est=C2$results, se=C2_se)
  C2_rl <- boots.lmom_rl_npt(data=m3.data.lm, a=lm_a, min.data=min.data, rp=rp0, conf=conf0, sd_scale=sd_scale, trans=trans0)
  C2_cvm<- cvm.test(C2$x, "pgev", C2$results[1], C2$results[2], C2$results[3])
  C2_ad <- ad.test(C2$x, "pgev", C2$results[1], C2$results[2], C2$results[3])
  
  # Combine all results----
  result <- list()
  result$mle <- data.frame(CT=A1_par, RT=B1_par, NPT=C1_par)
  result$lmom <- data.frame(CT=A2_par, RT=B2_par, NPT=C2_par)
  result$nllh <- data.frame(CT=A1_nllh, RT=B1_nllh, NPT=C1_nllh)
  result$mle_rl <- data.frame(CT=A1_rl, RT=B1_rl[,-1], NPT=C1_rl[,-1])
  result$lmom_rl <- data.frame(CT=A2_rl, RT=B2_rl[,-1], NPT=C2_rl[,-1])
  
  result$a.mle <- mle_a
  result$a.lm <- lm_a
  
  # GOF - MLE results--
  result$cvm_test_mle <- data.frame(cvm.stat = c(A1_cvm$statistic, B1_cvm$statistic, C1_cvm$statistic),
                                    cvm.pval = c(A1_cvm$p.value, B1_cvm$p.value, C1_cvm$p.value))
  
  result$ad_test_mle <- data.frame(ad.stat = c(A1_ad$statistic, B1_ad$statistic, C1_ad$statistic),
                                   ad.pval = c(A1_ad$p.value, B1_ad$p.value, C1_ad$p.value))
  
  # GOF - LM results--
  result$cvm_test_lmom <- data.frame(cvm.stat = c(A2_cvm$statistic, B2_cvm$statistic, C2_cvm$statistic),
                                     cvm.pval = c(A2_cvm$p.value, B2_cvm$p.value, C2_cvm$p.value))
  
  result$ad_test_lmom <- data.frame(ad.stat = c(A2_ad$statistic, B2_ad$statistic, C2_ad$statistic),
                                    ad.pval = c(A2_ad$p.value, B2_ad$p.value, C2_ad$p.value))
  
  return(result)
}


# 6) This function is required for converting the output from the function 'cal.g' to adjust the results table into a dataframe format.
save_output <- function(id=NULL, x.list=NULL, idx=NULL, set.past0=NULL){
  # set colname rl ---
  colnames(x.list[[paste0(set.past0,"_rl")]]) <- c("CT.period","CT.L", "CT.Est", "CT.U", "CT.se",
                                                   "RT.L", "RT.Est", "RT.U", "RT.se",
                                                   "NPT.L", "NPT.Est", "NPT.U", "NPT.se") 
  
  # separate output
  if(idx==1){
    out_method <- c(id=id, 
                    mu_se = as.numeric(x.list[[paste0(set.past0)]][1,c("CT.est","CT.se")]),
                    sig_se= as.numeric(x.list[[paste0(set.past0)]][2,c("CT.est","CT.se")]),
                    xi_se = as.numeric(x.list[[paste0(set.past0)]][3,c("CT.est","CT.se")]),
                    nllh = if(set.past0=="lmom"){
                      NA
                    }else{
                      x.list$nllh[,idx]
                    },
                    cvm_test = as.numeric(x.list[[paste0('cvm_test_',set.past0)]][idx,]),
                    ad_test = as.numeric(x.list[[paste0('ad_test_',set.past0)]][idx,]),
                    rl25 = x.list[[paste0(set.past0,"_rl")]][,"CT.Est"][1],
                    rl25.LU = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("CT.L","CT.U")][1,]),
                    rl25.se = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("CT.se")][1]),
                    
                    rl50 = x.list[[paste0(set.past0,"_rl")]][,"CT.Est"][2],
                    rl50.LU = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("CT.L","CT.U")][2,]), 
                    rl50.se = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("CT.se")][2]),
                    
                    rl100 = x.list[[paste0(set.past0,"_rl")]][,"CT.Est"][3],
                    rl100.LU = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("CT.L","CT.U")][3,]),
                    rl100.se = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("CT.se")][3])
    )
    
  }else if(idx==2){
    out_method <- c(id=id, 
                    mu_se = as.numeric(x.list[[paste0(set.past0)]][1,c("RT.est","RT.se")]),
                    sig_se= as.numeric(x.list[[paste0(set.past0)]][2,c("RT.est","RT.se")]),
                    xi_se = as.numeric(x.list[[paste0(set.past0)]][3,c("RT.est","RT.se")]),
                    nllh = if(set.past0=="lmom"){
                      NA
                    }else{
                      x.list$nllh[,idx]
                    },
                    cvm_test = as.numeric(x.list[[paste0('cvm_test_',set.past0)]][idx,]),
                    ad_test = as.numeric(x.list[[paste0('ad_test_',set.past0)]][idx,]),
                    rl25 = x.list[[paste0(set.past0,"_rl")]][,"RT.Est"][1],
                    rl25.LU = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("RT.L","RT.U")][1,]), 
                    rl25.se = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("RT.se")][1]),
                    
                    rl50 = x.list[[paste0(set.past0,"_rl")]][,"RT.Est"][2],
                    rl50.LU = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("RT.L","RT.U")][2,]), 
                    rl50.se = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("RT.se")][2]),
                    
                    rl100 = x.list[[paste0(set.past0,"_rl")]][,"RT.Est"][3],
                    rl100.LU = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("RT.L","RT.U")][3,]),
                    rl100.se = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("RT.se")][3])
    )
  }else if(idx==3){
    out_method <- c(id=id, 
                    mu_se = as.numeric(x.list[[paste0(set.past0)]][1,c("NPT.est","NPT.se")]),
                    sig_se= as.numeric(x.list[[paste0(set.past0)]][2,c("NPT.est","NPT.se")]),
                    xi_se = as.numeric(x.list[[paste0(set.past0)]][3,c("NPT.est","NPT.se")]),
                    opt_a.mle = x.list$a.mle,
                    opt_a.lm = x.list$a.lm,
                    
                    nllh = if(set.past0=="lmom"){
                      NA
                    }else{
                      x.list$nllh[,idx]
                    },
                    cvm_test = as.numeric(x.list[[paste0('cvm_test_',set.past0)]][idx,]),
                    ad_test = as.numeric(x.list[[paste0('ad_test_',set.past0)]][idx,]),
                    rl25 = x.list[[paste0(set.past0,"_rl")]][,"NPT.Est"][1],
                    rl25.LU = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("NPT.L","NPT.U")][1,]), 
                    rl25.se = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("NPT.se")][1]),
                    
                    rl50 = x.list[[paste0(set.past0,"_rl")]][,"NPT.Est"][2],
                    rl50.LU = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("NPT.L","NPT.U")][2,]), 
                    rl50.se = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("NPT.se")][2]),
                    
                    rl100 = x.list[[paste0(set.past0,"_rl")]][,"NPT.Est"][3],
                    rl100.LU = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("NPT.L","NPT.U")][3,]), 
                    rl100.se = as.numeric(x.list[[paste0(set.past0,"_rl")]][,c("NPT.se")][3])
    )
  }
  return(out_method)
}


# 7) Loop for obtaining output, where i,…,n represent the station indices.
# Initialize an empty output or dataframe to store results for each method: CT-MLE, RT-MLE, NPT-MLE, CT-LM, RT-LM, and NPT-LM
cb_output_df1 <- NULL
cb_output_df2 <- NULL
cb_output_df3 <- NULL
cb_output_df4 <- NULL
cb_output_df5 <- NULL
cb_output_df6 <- NULL

for (i in 1:length(f.list)){ # i=1
  out <- NULL
  mle_method1  <- mle_method2  <- mle_method3 <-  NULL
  lmom_method1 <- lmom_method2 <- lmom_method3 <- NULL
  
  # station id
  id <- str_extract(f.list[i],"\\d+")
  
  # read data 
  dat <- read.csv(paste0(id,".csv"))
  
  # extract min data
  min.data <- dat$min
  cat('Station ID',id,':::', "minimum of data:", min(dat$min), "\n")
  
  # calculate est. para and RL 
  out <- NULL
  out <- cal.g(min.data=min.data)

  # MLE ----
  # CT-GEVD
  mle_method1 <- save_output(id=id, x.list=out, idx=1, set.past0="mle")
  
  # RT-GEVD
  mle_method2 <- save_output(id=id, x.list=out, idx=2, set.past0="mle")
  
  # NPT-GEVD
  mle_method3 <- save_output(id=id, x.list=out, idx=3, set.past0="mle")
  
  # L-moments ----
  # CT-GEVD
  lmom_method1 <- save_output(id=id, x.list=out, idx=1, set.past0="lmom")
  
  # RT-GEVD
  lmom_method2 <- save_output(id=id, x.list=out, idx=2, set.past0="lmom")
  
  # NPT-GEVD
  lmom_method3 <- save_output(id=id, x.list=out, idx=3, set.past0="lmom")
  
  # Combine result
  cb_output_df1 <- rbind(cb_output_df1, mle_method1)
  cb_output_df2 <- rbind(cb_output_df2, mle_method2)
  cb_output_df3 <- rbind(cb_output_df3, mle_method3)
  cb_output_df4 <- rbind(cb_output_df4, lmom_method1)
  cb_output_df5 <- rbind(cb_output_df5, lmom_method2)
  cb_output_df6 <- rbind(cb_output_df6, lmom_method3)
  
  # Print id number
  cat('Station ID:',id,':::',paste(i, "/",length(f.list)), "\n")
}

# 8) Save output as an excel format --------
# Set the desired path where you want to save your files.
Save_PATH <- paste0("/Users/thanawanp/Desktop/R Code - NPT_GEVD/Output.xlsx")

# Create a new workbook
wb <- createWorkbook()
addWorksheet(wb, "CT-GEVD_MLE");  writeData(wb, "CT-GEVD_MLE", cb_output_df1)
addWorksheet(wb, "RT-GEVD_MLE");  writeData(wb, "RT-GEVD_MLE", cb_output_df2)
addWorksheet(wb, "NPT-GEVD_MLE"); writeData(wb, "NPT-GEVD_MLE", cb_output_df3)
addWorksheet(wb, "CT-GEVD_LM");   writeData(wb, "CT-GEVD_LM", cb_output_df4)
addWorksheet(wb, "RT-GEVD_LM");   writeData(wb, "RT-GEVD_LM", cb_output_df5)
addWorksheet(wb, "NPT-GEVD_LM");  writeData(wb, "NPT-GEVD_LM", cb_output_df6)
saveWorkbook(wb, Save_PATH, overwrite = TRUE)

