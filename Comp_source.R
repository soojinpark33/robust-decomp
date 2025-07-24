decompose_effects <- function(dat, 
                              Xs = NULL,  # exposure covariates, e.g. c("X1","X2")
                              Cs = NULL,  # baseline covariates, e.g. c("Age","Gender")
                              B = 200, 
                              n_sample = 1000, 
                              seed = 1234) {
  
  # helper to paste covariate lists
  paste_if <- function(vars) {
    if (is.null(vars) || length(vars)==0) "" else paste(vars, collapse = " + ")
  }
  
  Xpart <- paste_if(Xs)
  Cpart <- paste_if(Cs)
  adjust_part <- paste(c(Xpart, Cpart)[c(nzchar(Xpart), nzchar(Cpart))],
                       collapse = " + ")
  
  # 1) Difference-in-coefficients (DIC) ---------------------------------------
  dic_fun <- function(samp) {
    # total effect: Y ~ R + adjust
    if (nzchar(adjust_part)) {
      f_tot <-as.formula(paste("Y ~ R +", adjust_part)) 
    } else {
      f_tot <-as.formula(paste("Y ~ R"))
    }
    
    # direct effect: Y ~ R + M + adjust
    if (nzchar(adjust_part)) {
      f_dir <-as.formula(paste("Y ~ R + M +", adjust_part)) 
    } else {
      f_dir <-as.formula(paste("Y ~ R + M"))
    }
    
    
    tot <- lm(f_tot, data=samp)
    dir <- lm(f_dir, data=samp)
    tau   <- coef(tot)["R"]
    zeta  <- coef(dir)["R"]
    delta <- tau - zeta
    c(delta=delta, zeta=zeta, tau=tau)
  }
  
  set.seed(seed)
  dic_boot <- replicate(B, dic_fun(dat %>% sample_n(n_sample)), simplify=FALSE) %>% 
    bind_rows()
  dic_boot <- dic_boot %>% 
    rename(
      tau   = tau.R,
      delta = delta.R,
      zeta  = zeta.R
    )
  dic_sum <- tibble(
    Method = "DIC",
    Effect = c("Red","Rem","Ini"),
    Mean   = c(mean(dic_boot$delta), mean(dic_boot$zeta), mean(dic_boot$tau)),
    `2.5%` = c(quantile(dic_boot$delta,.025),
               quantile(dic_boot$zeta,.025),
               quantile(dic_boot$tau,.025)),
    `97.5%`= c(quantile(dic_boot$delta,.975),
               quantile(dic_boot$zeta,.975),
               quantile(dic_boot$tau,.975))
  )
  
  # 2) Kitagawa-Oaxaca-Blinder (KOB) ----------------------------------------
  kob_fun <- function(samp) {
    samp$R <- factor(samp$R)
    # build right-hand side for oaxaca:  Y ~ [Xs]+M+[Cs]  | R
    rhs_parts <- c(Xs, "M", Cs)
    rhs       <- paste(rhs_parts[rhs_parts!=""], collapse=" + ")
    form_ox   <- as.formula(paste("Y ~", rhs, "| R"))
    
    ox <- oaxaca(form_ox, data=samp)
    raw  <- ox$y$y.diff
    expM <- ox$twofold$variables[[6]]["M","coef(explained)"]
    c(delta=expM, zeta=raw-expM, tau=raw)
  }
  
  set.seed(seed+1)
  kob_boot <- replicate(B, kob_fun(dat %>% sample_n(n_sample)), simplify=FALSE) %>% 
    bind_rows()
  kob_sum <- tibble(
    Method = "KOB",
    Effect = c("Red","Rem","Ini"),
    Mean   = c(-mean(kob_boot$delta), -mean(kob_boot$zeta), -mean(kob_boot$tau)),
    `2.5%` = c(-quantile(kob_boot$delta,.025),
               -quantile(kob_boot$zeta,.025),
               -quantile(kob_boot$tau,.025)),
    `97.5%`= c(-quantile(kob_boot$delta,.975),
               -quantile(kob_boot$zeta,.975),
               -quantile(kob_boot$tau,.975))
  )
  
  # 3) Causal Decomposition Analysis (CDA) ----------------------------------
  cda_fun <- function(samp) {
    samp$R <- factor(samp$R)
    
    # mediator model: M ~ R + Cs
    if (nzchar(Cpart)) {
      form_m <-as.formula(paste("M ~ R +", Cpart)) 
    } else {
      form_m <-as.formula(paste("M ~ R"))
    } 
    
    
    # outcome model: Y ~ R + M + adjust
    if (nzchar(Cpart) & nzchar(Xpart)) {
      fit.y <- lm(Y ~ R + X + M + C, data=samp)  
    } else if (nzchar(Cpart)) {
      fit.y <- lm(Y ~ R + M + C, data=samp) 
    } else if (nzchar(Xpart)) {
      fit.y <- lm(Y ~ R + M + X, data=samp) 
    } else {
      fit.y <- lm(Y ~ R + M, data=samp) 
    }
    
    
    fit.m <- lm(form_m, data=samp)
    
    if (nzchar(Cpart)){
      fit_cda <- smi(fit.m=fit.m, fit.y=fit.y, sims=200, conf.level = .95,
                     covariates =  c("C"),
                     treat="R")
      
      c(tau   = fit_cda$result[1,"estimate"],
        zeta  = fit_cda$result[2,"estimate"],
        delta = fit_cda$result[3,"estimate"])
    } else {
      tau   <- coef( lm(Y ~ R, data = samp)       )[2]
      delta <- coef(fit.y)[3] * coef(fit.m)[2]
      
      c(tau = tau,
        delta = delta, 
        zeta = tau-delta)
    }
    
  }
  
  set.seed(seed+2)
  cda_boot <- replicate(B, cda_fun(dat %>% sample_n(n_sample)), simplify=FALSE) %>% 
    bind_rows()
  if (nzchar(Cpart)) {
    cda_sum <- tibble(
      Method = "CDA",
      Effect = c("Red","Rem","Ini"),
      Mean   = c(mean(cda_boot$delta), mean(cda_boot$zeta), mean(cda_boot$tau)),
      `2.5%` = c(quantile(cda_boot$delta,.025),
                 quantile(cda_boot$zeta,.025),
                 quantile(cda_boot$tau,.025)),
      `97.5%`= c(quantile(cda_boot$delta,.975),
                 quantile(cda_boot$zeta,.975),
                 quantile(cda_boot$tau,.975))
    )} else {
      cda_sum <- tibble(
        Method = "CDA",
        Effect = c("Red","Rem","Ini"),
        Mean   = c(mean(cda_boot$delta.M), mean(cda_boot$zeta.R1), mean(cda_boot$tau.R1)),
        `2.5%` = c(quantile(cda_boot$delta.M,.025),
                   quantile(cda_boot$zeta.R1,.025),
                   quantile(cda_boot$tau.R1,.025)),
        `97.5%`= c(quantile(cda_boot$delta.M,.975),
                   quantile(cda_boot$zeta.R1,.975),
                   quantile(cda_boot$tau.R1,.975))
      ) }
  
  # bind and return all three
  bind_rows(dic_sum, kob_sum, cda_sum)
}

sens.for.se <- function(boot.res, fit.y, fit.m = NULL, mediators = NULL, covariates, treat, sel.lev.treat,
                        conf.level = 0.95, ry, rm){
  
  # Match arguments
  call <- match.call()
  data <- model.frame(fit.y)
  
  # Warning for inappropriate settings
  if(!(sel.lev.treat %in% levels(data[, treat]))){
    stop("'sel.lev.treat' must be one of the levels of treatment")
  }
  if(sel.lev.treat == levels(data[, treat])[1]){
    stop("'sel.lev.treat' must not be the reference level of treatment")
  }
  
  # Extract outcome model and data from bootstrapping results
  outcome <- all.vars(formula(fit.y))[[1]]
  len.treat <- length(levels(data[, treat]))
  # new data with releveled treatment
  data.new <- data
  data.new[, treat] <- relevel(data.new[, treat], ref = sel.lev.treat)
  # new outcome model with releveled treatment
  fit.y.new <- update(fit.y, data = data.new)
  lev <- which(levels(data[, treat]) == sel.lev.treat)
  
  if(is.null(fit.m) & !is.null(mediators)){
    
    # Possibility of multiple mediators
    if(length(class(fit.m)) == 1 && class(fit.m) == "list"){
      isMultiMediators <- TRUE
      num.meds <- length(fit.m)
    } else {
      isMultiMediators <- FALSE
      num.meds <- 1
    }
    ## Fit mediator model(s)
    # make an empty list to save mediator model(s)
    fit.meds <- vector(mode = "list", length = num.meds)
    for (i in 1:num.meds) {
      # The formula for each mediator model: lm(M1 ~ R + C); lm(M2 ~ R + C);...
      med.formula <- as.formula(paste(mediators[i], paste(c(treat, covariates), collapse = " + "), sep = " ~ "))
      if(is.factor(data[, mediators[i]])){
        fit.meds[[i]] <- lm(med.formula, data = data)
      } else {
        fit.meds[[i]] <- glm(med.formula, data = data)
      }
    }
    #fit.meds.new1 <- update(fit.meds[[1]], data = data.new)
    #vcov(fit.meds.new1)["racesex41", "racesex41"]
    
  } else if (!is.null(fit.m) & is.null(mediators)){
    
    # Possibility of multiple mediators
    if(length(class(fit.m)) == 1 && class(fit.m) == "list"){
      isMultiMediators <- TRUE
      num.meds <- length(fit.m)
    } else {
      isMultiMediators <- FALSE
      num.meds <- 1
    }
    
    fit.meds <- vector(mode = "list", length = num.meds)
    mediators <- rep(NA, num.meds)
    for (i in 1:num.meds) {
      fit.meds[[i]] <- fit.m[[i]]
      mediators[i] <- all.vars(formula(fit.m[[i]]))[[1]]
    }
    
  } else if (!is.null(fit.m) & !is.null(mediators)) {
    
    # Possibility of multiple mediators
    if(length(class(fit.m)) == 1 && class(fit.m) == "list"){
      isMultiMediators <- TRUE
      num.meds <- length(fit.m)
    } else {
      isMultiMediators <- FALSE
      num.meds <- 1
    }
    
    fit.meds <- vector(mode = "list", length = num.meds)
    mediators0 <- rep(NA, num.meds)
    for (i in 1:num.meds) {
      if(!isMultiMediators){
        fit.meds[[i]] <- fit.m
        mediators0[i] <- all.vars(formula(fit.m))[[1]]
      } else {
        fit.meds[[i]] <- fit.m[[i]]
        mediators0[i] <- all.vars(formula(fit.m[[i]]))[[1]]
      }
      if(mediators0[i] != mediators[i]){
        stop("Response variable(s) of 'fit.m' and 'mediators' must be match.")
      }
    }
    
  } else if (is.null(fit.m) & is.null(mediators)) {
    stop("Either 'fit.m' or 'mediators' must not be NULL.")
  }
  
  ##1. Calculate the SE of gamma_dm and df from the coefficient of M for a comparison group in fit.y.
  if (num.meds == 1) {
    meds <- all.vars(formula(fit.m))[[1]]
    if(is.factor(data[, mediators])){
      meds <- which(meds == substring(colnames(vcov(fit.y.new)), 1, nchar(meds)))[1]
    }
    var_gamma <- vcov(fit.y.new)[meds, meds]
    beta.resm.est <- coef(fit.y.new)[meds]
  } else if (num.meds == 2) {
    meds1 <- all.vars(formula(fit.m[[1]]))[[1]]
    meds2 <- all.vars(formula(fit.m[[2]]))[[1]]
    if(is.factor(data[, mediators[1]])){
      meds1 <- which(meds1 == substring(colnames(vcov(fit.y.new)), 1, nchar(meds1)))[1]
    }
    if(is.factor(data[, mediators[2]])){
      meds2 <- which(meds2 == substring(colnames(vcov(fit.y.new)), 1, nchar(meds2)))[1]
    }
    var_gamma <- vcov(fit.y.new)[meds1, meds1] + vcov(fit.y.new)[meds2, meds2] +
      2 * vcov(fit.y.new)[meds1, meds2]
    beta.resm.est <- sum(coef(fit.y.new)[c(meds1, meds2)]) ##??? when there are two Ms.
  } else if (num.meds > 2) {
    # will be added.
    stop("The case of three or more mediators is not supported yet.")
  }
  se_gamma <- sqrt(var_gamma)
  df <- fit.y.new$df.residual
  
  ##2. Calculate the effect of R on D and M (the coefficient for R in fit.m).
  treat.lev <- paste(treat, sel.lev.treat, sep = "")
  # sample covariance of initial disparity and disparity reduction
  cov.ini.red <- cov(boot.res$all.result[3 * (lev - 1) - 2, ],
                     boot.res$all.result[3 * (lev - 1), ])
  # sample covariance of initial disparity and alpha_r
  cov.ini.alphar <- cov(boot.res$all.result[3 * (lev - 1) - 2, ],
                        boot.res$alpha.r[, lev - 1])
  # sample covariance of initial disparity and se_gamma
  cov.ini.segamma <- cov(boot.res$all.result[3 * (lev - 1) - 2, ],
                         boot.res$se.gamma[, lev - 1])
  # sample variance of alpha_r
  var_alphahat.r <- var(boot.res$alpha.r[, lev - 1]) # 09/05/2022, 2/28/2023
  
  if (num.meds == 1) {
    alpha.r <- coef(fit.meds[[1]])[treat.lev] # coef of R[lev] on M = abs(alpha_r)
    ab <- abs(alpha.r)
    var_ab <- vcov(fit.meds[[1]])[treat.lev, treat.lev]
  } else if (num.meds == 2) {
    alpha.r <- coef(fit.meds[[1]])[treat.lev] + coef(fit.meds[[2]])[treat.lev]  # 09/05/2022
    ab <- abs(alpha.r)
    var_ab <- vcov(fit.meds[[1]])[treat.lev, treat.lev] +
      vcov(fit.meds[[2]])[treat.lev, treat.lev] - 2 * cov.ini.red
  } else if (num.meds > 2) {
    # will be added.
    stop("The case of three or more mediators is not supported yet.")
  }
  
  ## Calculate the variance of initial disparity # Use the one from the estimator function. 
  # Y ~ R + C
  ini.formula <- as.formula(paste(outcome, paste(c(treat, covariates), collapse = " + "), sep = " ~ "))
  ini <- lm(ini.formula, data = data)
  #var_tau <- vcov(ini)[treat.lev, treat.lev]
  var_tau <- var(boot.res$all.result[3 * (lev - 1) - 2, ])  # 09/05/2022
  
  ##3.
  # Point Estimate of Disparity Reduction for ref vs lev
  res.red <- boot.res$result[3 * (lev - 1), 1]
  # Point Estimate of Disparity Remaining for ref vs lev
  res.rem <- boot.res$result[3 * (lev - 1) - 1, 1]
  red <- rem <- lower_red <- lower_rem <- upper_red <- upper_rem <- var_ngamma <- matrix(NA, length(ry), length(rm))
  
  for (i in 1:length(ry)){
    for (j in 1:length(rm)){
      
      # true disp reduction (= est - bias)
      red[i, j] <- res.red - ab * se_gamma * sqrt(ry[i] * rm[j] * df)/sqrt(1 - rm[j]) #res.red
      beta.m <- beta.resm.est + ifelse(beta.resm.est > 0, - 1, 1) * se_gamma * sqrt(ry[i] * rm[j] * df)/sqrt(1 - rm[j]) # 09/05/2022
      # new variance of gamma
      var_ngamma[i, j] <- var_gamma * (1 - ry[i])/(1 - rm[j]) * (df/(df - 1)) 
      # new CI for disp reduction
      qn <- qnorm(1/2 + conf.level/2)
      se_red <- sqrt(alpha.r^2 * var_ngamma[i, j] + beta.m^2 * var_alphahat.r)  # Eq (15)
      lower_red[i, j] <- red[i, j] - qn * se_red
      upper_red[i, j] <- red[i, j] + qn * se_red
      
      # true disp remaining
      rem[i, j] <- res.rem + ab * se_gamma * sqrt(ry[i] * rm[j] * df)/sqrt(1 - rm[j]) #res.rem
      # new CI for disp remaining
      qn <- qnorm(1/2 + conf.level/2)
      k <- 1 * sqrt(ry[i] * rm[j] * df)/sqrt(1 - rm[j])
      se_rem <- sqrt(var_tau + se_red^2 - 2 * cov.ini.red +
                       2 * k * mean(boot.res$se.gamma[, lev - 1]) * cov.ini.alphar + # se_gamma
                       2 * k * mean(boot.res$alpha.r[, lev - 1]) * cov.ini.segamma)
      lower_rem[i, j] <- rem[i, j] - qn * se_rem
      upper_rem[i, j] <- rem[i, j] + qn * se_rem
      
    }
  }
  
  # output
  output <- c(lower_rem, upper_rem, lower_red, upper_red, red, rem)
  names(output) <- c("lower_rem", "upper_rem", "lower_red", "upper_red","red","rem")
  return(output)
  
}

