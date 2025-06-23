library(oaxaca)   
library(truncnorm)
library(dplyr)
library(causal.decomp)
source("U:/Comp_Decomp/Comp_source.R")


set.seed(123)          # reproducible
n  <- 1000000            # sample size
R <- rbinom(n,1,0.5)
alpha = c(-1, -0.3)
beta = c(0.6, -0.3, 0.5, 0.5)
gamma = c(1, 0.5, 0.1, -0.5, 0.8)      
delta = c(0.2, 0.7, -0.5, 0.4, 0.7, 0.8)

#:::::::::::::::::::::
# Scenario 1 (no C and X)
#:::::::::::::::::::::
M <- gamma[1] + gamma[2]*R + rnorm(n)
Y <- delta[1] + delta[2]*R + delta[3]*M + rnorm(n)
podat_noCX <- data.frame(R , M, Y)
# #true value of cda
tau_cda <- coef(lm(Y~R , data=podat_noCX))[2]   
exp_cda <- delta[3]*coef(lm(M~R, data=podat_noCX))[2] 
unexp_cda <- tau_cda-exp_cda 

#:::::::::::::::::::::
# Scenario 2 (C only)
#:::::::::::::::::::::
C <- rnorm(n,mean=1,sd=1)
p_R1   <- plogis(alpha[1] + alpha[1] * C)   # P(R = 1 | C)
R      <- rbinom(n, size = 1, prob = p_R1)
M <- gamma[1] + gamma[2]*R + gamma[4]*C + rnorm(n)
Y <- delta[1] + delta[2]*R + delta[3]*M + delta[5]*C + rnorm(n)
podat_C <- data.frame(R , M, Y, C)
# #true value of cda
tau_cda <- coef(lm(Y~R +C, data=podat_C))[2]   
exp_cda <- delta[3]*coef(lm(M~R+C, data=podat_C))[2] 
unexp_cda <- tau_cda-exp_cda 

#::::::::::::::::::::::
# Scenario 3 (X only)
#::::::::::::::::::::::
R <- rbinom(n,1,0.5)
X     <- beta[1] + beta[2] * R + rnorm(n, sd = 1)
M <- gamma[1] + gamma[2]*R + gamma[3]*X + rnorm(n)
Y <- delta[1] + delta[2]*R + delta[3]*M + delta[4]*X + rnorm(n)
podat_X <- data.frame(R , M, Y, X)
# #true value of cda
tau_cda <- coef(lm(Y~R , data=podat_X))[2]   
exp_cda <- delta[3]*coef(lm(M~R, data=podat_X))[2] 
unexp_cda <- tau_cda-exp_cda #-0.1072181

#:::::::::::::::::::::::::::::
# Scenario 4 (X and C)
#:::::::::::::::::::::::::::::
p_R1   <- plogis(alpha[1] + alpha[2] * C)   # P(R = 1 | C)
R      <- rbinom(n, size = 1, prob = p_R1)
X     <- beta[1] + beta[2] * R + beta[3] * C + rnorm(n, sd = 1)
M      <- gamma[1] + gamma[2] * R + gamma[3] * X + gamma[4] * C + rnorm(n, sd = 1)
Y      <- delta[1] + delta[2] * R  + delta[3] * M + delta[4] * X+ delta[5] * C +
  rnorm(n, sd = 1)
popdat <- data.frame(C, R, X, M, Y)
#true value of dic
tau_dic <- coef(lm(Y ~ R + X + C , data=popdat))[2]
unexp_dic <- coef(lm(Y ~ R + X + C + M , data=popdat))[2]
exp_dic <- tau_dic- unexp_dic 
# true value of cda
tau_cda <- coef(lm(Y~R +C, data=popdat))[2]   
exp_cda <- delta[3]*coef(lm(M~R + C, data=popdat))[2] 
unexp_cda <- tau_cda-exp_cda 

#:::::::::::::::::::::::::::::
# Scenario 5 (X, C, and U)
#:::::::::::::::::::::::::::::
U <-  rnorm(n, mean = 0, sd = 1)
p_R1   <- plogis(alpha[1] + alpha[2] * C)   # P(R = 1 | C)
R      <- rbinom(n, size = 1, prob = p_R1)
X     <- beta[1] + beta[2] * R + beta[3] * C + beta[4]* U + rnorm(n, sd = 1)
M      <- gamma[1] + gamma[2] * R + gamma[3] * X + gamma[4] * C + gamma[5]*U + rnorm(n, sd = 1)
Y      <- delta[1] + delta[2] * R + delta[3] * M + delta[4] * X + delta[5] * C +
  rnorm(n, sd = 1)
## ---- bundle into a data frame -----------------------------------------
popdat1 <- data.frame(C, R, X, M, Y,U)

#true value of kob
tau_ox <- coef(lm(Y~R , data=popdat1))[2]  
exp_ox <- delta[3]*coef(lm(M~R , data=popdat1))[2]
unexp_ox <- tau_ox-exp_ox 

#true value of cda
tau_cda <- coef(lm(Y~R + C, data=popdat1))[2]   
exp_cda <- delta[3]*coef(lm(M~R +C, data=popdat1))[2] 
unexp_cda <- tau_cda-exp_cda 

#:::::::::::::::::::::::::::::
# Scenario 6 (X, C, and U)
#:::::::::::::::::::::::::::::
M      <- gamma[1] + gamma[2] * R + gamma[3] * X + gamma[4] * C + gamma[5]* U + rnorm(n, sd = 1)
Y      <- delta[1] + delta[2] * R + delta[3] * M + delta[4] * X + delta[5] * C + delta[6]* U +
  rnorm(n, sd = 1)
## ---- bundle into a data frame -----------------------------------------
popdat2 <- data.frame(C, R, X, M, Y,U)

#_________________

 n_iter<-200
# Scenario 1: neither C nor X
res_no <- decompose_effects(dat=podat_noCX, Xs = NULL, Cs = NULL, B=n_iter, n_sample=2000)
print(res_no)

# Scenario 2: only C
res_C <- decompose_effects(dat=podat_C, Xs = NULL, Cs = "C", B=n_iter, n_sample=2000)
print(res_C)

# Scenario 3: only X
res_X <- decompose_effects(dat=podat_X, Xs = "X", Cs = NULL, B=n_iter, n_sample=2000)
print(res_X)

#Scenario 4
res1 <- decompose_effects(dat=popdat, Xs = "X", Cs = "C", B=n_iter, n_sample=2000)
print(res1)

#Scenario 5
res2 <- decompose_effects(dat=popdat1, Xs = "X", Cs = "C", B=n_iter, n_sample=2000)
print(res2)

#Scenario 6
res3 <- decompose_effects(dat=popdat2, Xs = "X", Cs = "C", B=n_iter, n_sample=2000)
print(res3)

#::::::::::::::::::::::::::::::::
## Adjusted CDA
resM <- lm(M ~ R + X + C, data = popdat2)$residuals
# regress U on R,X,C
resU <- lm(U ~ R + X + C, data = popdat2)$residuals
# their correlation
pcor_MU <- cor(resM, resU)

# 3) partial correlation Corr(Y, U | R,X,M,C)
# regress Y on R,X,M,C
resY <- lm(Y ~ R + X + M + C, data = popdat2)$residuals
# regress U on R,X,M,C
resU2 <- lm(U ~ R + X + M + C, data = popdat2)$residuals
# their correlation
pcor_YU <- cor(resY, resU2)

ini <- rem <- red <- matrix(NA, n_iter, 1)
n_sample <- 2000
for (i in seq_len(n_iter)) {
  
  samp <- popdat2 %>%                         # simple-random sample
    sample_n(n_sample, replace = FALSE)
  samp$R <- as.factor(samp$R)
  fit.m <- lm(M ~ R  + C , data = samp)
  fit.y <- lm(Y ~ R  + X  + M + C  , data = samp)
  fit <- smi(fit.m = fit.m, fit.y = fit.y, sims = 200, conf.level = .95,
             covariates = c("C"), treat = "R")
  
  # Store point estimate
  sens.res <- sens.for.se(boot.res = fit, fit.y = fit.y, fit.m = fit.m, mediators = "M",
                          covariates = "C", treat = "R", sel.lev.treat = "1", ry = pcor_YU^2, rm = pcor_MU^2)
  
  rem[i] <- sens.res[6]
  red[i] <- sens.res[5]
  ini[i] <- rem[i] + red[i]
}

summary_df <- data.frame(
  Effect      = c("Red", "Rem", "Ini"),
  Mean        = c(mean(red), mean(rem), mean(ini)),
  SD          = c(sd(red),   sd(rem), sd(ini)),
  `2.5%`      = c(quantile(red, .025), quantile(rem, .025), quantile(ini, .025)),
  `97.5%`     = c(quantile(red, .975), quantile(rem, .975), quantile(ini, .975))
)


save.image(file="sce2.RData")
