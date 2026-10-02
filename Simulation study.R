library(oaxaca)   
library(truncnorm)
library(dplyr)
library(causal.decomp)
# Local working copy (original is U:/Comp_Decomp/Comp_source.R — do not edit)
source("C:/Users/sykim/Desktop/R/Robust/comp_source_rv.R")


set.seed(123)          # reproducible
n  <- 1000000            # sample size
Cs_all <- c("C1", "C2")
Xs_all <- c("X1", "X2", "X3")
b_YM   <- -0.5           # structural M -> Y coefficient (CDA/KOB product-of-coefficients)

# Baseline C: different distributions (generated once, reused in Scenarios 2/4–6)
#   C1 age-like integer 18–80, stored centered
#   C2 binary (sex-like)
gen_C <- function(n) {
  C1_raw <- round(rtruncnorm(n, a = 18, b = 80, mean = 50, sd = 12))
  data.frame(
    C1 = as.numeric(scale(C1_raw, center = TRUE, scale = FALSE)),
    C2 = rbinom(n, 1, 0.55)
  )
}

# R | C  (logistic). use_C = FALSE keeps P(R=1)=0.5
gen_R <- function(n, C = NULL, use_C = TRUE) {
  if (use_C) rbinom(n, 1, plogis(-0.35 - 0.010 * C$C1 + 0.40 * C$C2))
  else       rbinom(n, 1, 0.5)
}

# Intermediate X: different distributions
#   X1 continuous, Gaussian noise (SES-like)
#   X2 continuous, smaller noise (health-like)
#   X3 binary logistic (event-like)
gen_X <- function(R, C, U = NULL, use_C = TRUE) {
  n <- length(R)
  u   <- if (is.null(U)) 0 else U
  c1  <- if (use_C) C$C1 else 0
  c2  <- if (use_C) C$C2 else 0
  data.frame(
    X1 =  0.40 - 0.35 * R + 0.20 * c1 + 0.15 * c2 + 0.80 * u + rnorm(n, sd = 1.00),
    X2 =  0.10 - 0.20 * R + 0.10 * c1 + 0.25 * c2 + 0.80 * u + rnorm(n, sd = 0.60),
    X3 = rbinom(n, 1, plogis(-1.20 + 0.40 * R + 0.08 * c1 + 0.30 * c2 + 0.80 * u))
  )
}

gen_M <- function(R, C, X, U = NULL, use_C = TRUE, use_X = TRUE) {
  n <- length(R)
  u  <- if (is.null(U)) 0 else U
  c1 <- if (use_C) C$C1 else 0
  c2 <- if (use_C) C$C2 else 0
  x1 <- if (use_X) X$X1 else 0
  x2 <- if (use_X) X$X2 else 0
  x3 <- if (use_X) X$X3 else 0
  1.0 + 0.50 * R + 0.12 * x1 + 0.18 * x2 + 0.22 * x3 -
    0.35 * c1 - 0.20 * c2 + 1.20 * u + rnorm(n)
}

gen_Y <- function(R, C, X, M, U = NULL, use_C = TRUE, use_X = TRUE) {
  n <- length(R)
  u  <- if (is.null(U)) 0 else U
  c1 <- if (use_C) C$C1 else 0
  c2 <- if (use_C) C$C2 else 0
  x1 <- if (use_X) X$X1 else 0
  x2 <- if (use_X) X$X2 else 0
  x3 <- if (use_X) X$X3 else 0
  0.2 + 0.70 * R + b_YM * M + 0.25 * x1 + 0.15 * x2 + 0.20 * x3 +
    0.04 * c1 + 0.30 * c2 + 1.20 * u + rnorm(n)
}

est_R <- function(f, dat) unname(coef(lm(f, data = dat))["R"])
fml <- function(y, ...) {
  parts <- unlist(list(...), use.names = FALSE)
  parts <- parts[!is.na(parts) & nzchar(as.character(parts))]
  as.formula(paste(y, "~", paste(parts, collapse = " + ")))
}
truth_vec <- function(ini, red) c(Red = unname(red), Rem = unname(ini - red), Ini = unname(ini))
truth_cda <- function(dat, Cvars = NULL) {
  if (length(Cvars)) truth_vec(est_R(fml("Y", "R", Cvars), dat),
                               b_YM * est_R(fml("M", "R", Cvars), dat))
  else               truth_vec(est_R(Y ~ R, dat), b_YM * est_R(M ~ R, dat))
}
truth_kob <- function(dat) {
  truth_vec(est_R(Y ~ R, dat), b_YM * est_R(M ~ R, dat))
}
truth_dic <- function(dat, Cvars = NULL, Xvars = NULL, rem_U = FALSE, ini_U = FALSE) {
  ini <- est_R(fml("Y", "R", Xvars, Cvars, if (ini_U) "U"), dat)
  rem <- est_R(fml("Y", "R", Xvars, "M", Cvars, if (rem_U) "U"), dat)
  c(Red = unname(ini - rem), Rem = unname(rem), Ini = unname(ini))
}
expand_truth <- function(scenario, method, tv) {
  data.frame(Scenario = scenario, Method = method,
             Effect = names(tv), True = as.numeric(tv),
             row.names = NULL)
}

#:::::::::::::::::::::
# Scenario 1 (no C and X)
#:::::::::::::::::::::
R <- gen_R(n, use_C = FALSE)
M <- gen_M(R, NULL, NULL, use_C = FALSE, use_X = FALSE)
Y <- gen_Y(R, NULL, NULL, M, use_C = FALSE, use_X = FALSE)
podat_noCX <- data.frame(R, M, Y)
truth_s1_cda <- truth_cda(podat_noCX)
truth_s1_dic <- truth_dic(podat_noCX)
truth_s1_kob <- truth_kob(podat_noCX)

#:::::::::::::::::::::
# Scenario 2 (C only)
#:::::::::::::::::::::
Cdat <- gen_C(n)
R <- gen_R(n, Cdat, use_C = TRUE)
M <- gen_M(R, Cdat, NULL, use_C = TRUE, use_X = FALSE)
Y <- gen_Y(R, Cdat, NULL, M, use_C = TRUE, use_X = FALSE)
podat_C <- data.frame(R, M, Y, Cdat)
truth_s2_cda <- truth_cda(podat_C, Cs_all)
truth_s2_dic <- truth_dic(podat_C, Cvars = Cs_all)
truth_s2_kob <- truth_kob(podat_C)

#::::::::::::::::::::::
# Scenario 3 (X only)
#::::::::::::::::::::::
R <- gen_R(n, use_C = FALSE)
Xdat <- gen_X(R, NULL, use_C = FALSE)
M <- gen_M(R, NULL, Xdat, use_C = FALSE, use_X = TRUE)
Y <- gen_Y(R, NULL, Xdat, M, use_C = FALSE, use_X = TRUE)
podat_X <- data.frame(R, M, Y, Xdat)
truth_s3_cda <- truth_cda(podat_X)
truth_s3_dic <- truth_dic(podat_X, Xvars = Xs_all)
truth_s3_kob <- truth_kob(podat_X)

#:::::::::::::::::::::::::::::
# Scenario 4 (X and C)
#:::::::::::::::::::::::::::::
R <- gen_R(n, Cdat, use_C = TRUE)
Xdat <- gen_X(R, Cdat, use_C = TRUE)
M <- gen_M(R, Cdat, Xdat, use_C = TRUE, use_X = TRUE)
Y <- gen_Y(R, Cdat, Xdat, M, use_C = TRUE, use_X = TRUE)
popdat <- data.frame(Cdat, R, Xdat, M, Y)
truth_s4_cda <- truth_cda(popdat, Cs_all)
truth_s4_dic <- truth_dic(popdat, Cvars = Cs_all, Xvars = Xs_all)
truth_s4_kob <- truth_kob(popdat)

#:::::::::::::::::::::::::::::
# Scenario 5 (X, C, and U): U -> X and U -> M; not Y
#:::::::::::::::::::::::::::::
U <- rnorm(n)
R <- gen_R(n, Cdat, use_C = TRUE)
Xdat <- gen_X(R, Cdat, U = U, use_C = TRUE)
M <- gen_M(R, Cdat, Xdat, U = U, use_C = TRUE, use_X = TRUE)
Y <- gen_Y(R, Cdat, Xdat, M, U = NULL, use_C = TRUE, use_X = TRUE)
popdat1 <- data.frame(Cdat, R, Xdat, M, Y, U)
truth_s5_cda <- truth_cda(popdat1, Cs_all)
truth_s5_dic <- truth_dic(popdat1, Cvars = Cs_all, Xvars = Xs_all,
                          ini_U = TRUE, rem_U = TRUE)
truth_s5_kob <- truth_kob(popdat1)

#:::::::::::::::::::::::::::::
# Scenario 6 (X, C, and U): U -> M and U -> Y; not X
#:::::::::::::::::::::::::::::
Xdat <- gen_X(R, Cdat, U = NULL, use_C = TRUE)
M <- gen_M(R, Cdat, Xdat, U = U, use_C = TRUE, use_X = TRUE)
Y <- gen_Y(R, Cdat, Xdat, M, U = U, use_C = TRUE, use_X = TRUE)
popdat2 <- data.frame(Cdat, R, Xdat, M, Y, U)
truth_s6_cda <- truth_cda(popdat2, Cs_all)
truth_s6_dic <- truth_dic(popdat2, Cvars = Cs_all, Xvars = Xs_all, rem_U = TRUE)
truth_s6_kob <- truth_kob(popdat2)

# Table 1 uses one CDA true column for all methods.
# Tables 2–3 use method-specific trues; DIC.mod / KOB.mod share the CDA true.
truths <- rbind(
  expand_truth(1, "CDA",     truth_s1_cda), expand_truth(1, "DIC",     truth_s1_cda),
  expand_truth(1, "KOB",     truth_s1_cda), expand_truth(1, "DIC.mod", truth_s1_cda),
  expand_truth(1, "KOB.mod", truth_s1_cda),
  expand_truth(2, "CDA",     truth_s2_cda), expand_truth(2, "DIC",     truth_s2_cda),
  expand_truth(2, "KOB",     truth_s2_cda), expand_truth(2, "DIC.mod", truth_s2_cda),
  expand_truth(2, "KOB.mod", truth_s2_cda),
  expand_truth(3, "CDA",     truth_s3_cda), expand_truth(3, "DIC",     truth_s3_cda),
  expand_truth(3, "KOB",     truth_s3_cda), expand_truth(3, "DIC.mod", truth_s3_cda),
  expand_truth(3, "KOB.mod", truth_s3_cda),
  expand_truth(4, "CDA",     truth_s4_cda), expand_truth(4, "DIC",     truth_s4_cda),
  expand_truth(4, "KOB",     truth_s4_cda), expand_truth(4, "DIC.mod", truth_s4_cda),
  expand_truth(4, "KOB.mod", truth_s4_cda),
  expand_truth(5, "CDA",     truth_s5_cda), expand_truth(5, "DIC",     truth_s5_dic),
  expand_truth(5, "KOB",     truth_s5_kob), expand_truth(5, "DIC.mod", truth_s5_cda),
  expand_truth(5, "KOB.mod", truth_s5_cda),
  expand_truth(6, "CDA",     truth_s6_cda), expand_truth(6, "DIC",     truth_s6_dic),
  expand_truth(6, "KOB",     truth_s6_kob), expand_truth(6, "DIC.mod", truth_s6_cda),
  expand_truth(6, "KOB.mod", truth_s6_cda),
  expand_truth(6, "CDA.adj", truth_s6_cda)
)
truths_method <- rbind(
  expand_truth(1, "CDA", truth_s1_cda), expand_truth(1, "DIC", truth_s1_dic),
  expand_truth(1, "KOB", truth_s1_kob),
  expand_truth(2, "CDA", truth_s2_cda), expand_truth(2, "DIC", truth_s2_dic),
  expand_truth(2, "KOB", truth_s2_kob),
  expand_truth(3, "CDA", truth_s3_cda), expand_truth(3, "DIC", truth_s3_dic),
  expand_truth(3, "KOB", truth_s3_kob),
  expand_truth(4, "CDA", truth_s4_cda), expand_truth(4, "DIC", truth_s4_dic),
  expand_truth(4, "KOB", truth_s4_kob),
  expand_truth(5, "CDA", truth_s5_cda), expand_truth(5, "DIC", truth_s5_dic),
  expand_truth(5, "KOB", truth_s5_kob),
  expand_truth(6, "CDA", truth_s6_cda), expand_truth(6, "DIC", truth_s6_dic),
  expand_truth(6, "KOB", truth_s6_kob)
)
#_________________

 n_iter<-500
out_dir <- "C:/Users/sykim/Desktop/R/Robust/sim_rv_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_tbl <- function(obj, name) {
  path <- file.path(out_dir, paste0(name, ".csv"))
  write.csv(as.data.frame(obj), path, row.names = FALSE)
  cat("Wrote", path, "\n")
}
attach_true <- function(res, scenario, truth_df = truths) {
  res <- as.data.frame(res)
  res$Scenario <- scenario
  tr <- truth_df[truth_df$Scenario == scenario, c("Method", "Effect", "True")]
  merge(res, tr, by = c("Method", "Effect"), all.x = TRUE, sort = FALSE)
}

# Each decompose_effects() call returns all 5 methods:
#   DIC, KOB, CDA, DIC.mod, KOB.mod
# Scenario 1: neither C nor X
res_no <- decompose_effects(dat=podat_noCX, Xs = NULL, Cs = NULL, B=n_iter, n_sample=2000)
res_no <- attach_true(res_no, 1)
print(res_no)
write_tbl(res_no, "scenario1_noCX")

# Scenario 2: only C
res_C <- decompose_effects(dat=podat_C, Xs = NULL, Cs = Cs_all, B=n_iter, n_sample=2000)
res_C <- attach_true(res_C, 2)
print(res_C)
write_tbl(res_C, "scenario2_C_only")

# Scenario 3: only X
res_X <- decompose_effects(dat=podat_X, Xs = Xs_all, Cs = NULL, B=n_iter, n_sample=2000)
res_X <- attach_true(res_X, 3)
print(res_X)
write_tbl(res_X, "scenario3_X_only")

# Scenario 4 (C and X)
res1 <- decompose_effects(dat=popdat, Xs = Xs_all, Cs = Cs_all, B=n_iter, n_sample=2000)
res1 <- attach_true(res1, 4)
print(res1)
write_tbl(res1, "scenario4_CX")

# Scenario 5
res2 <- decompose_effects(dat=popdat1, Xs = Xs_all, Cs = Cs_all, B=n_iter, n_sample=2000)
res2 <- attach_true(res2, 5)
print(res2)
write_tbl(res2, "scenario5_U_on_XM")

# Scenario 6
res3 <- decompose_effects(dat=popdat2, Xs = Xs_all, Cs = Cs_all, B=n_iter, n_sample=2000)
res3 <- attach_true(res3, 6)
print(res3)
write_tbl(res3, "scenario6_U_on_MY")

#::::::::::::::::::::::::::::::::
## Adjusted CDA
resM <- lm(M ~ R + X1 + X2 + X3 + C1 + C2, data = popdat2)$residuals
resU <- lm(U ~ R + X1 + X2 + X3 + C1 + C2, data = popdat2)$residuals
pcor_MU <- cor(resM, resU)

resY  <- lm(Y ~ R + X1 + X2 + X3 + M + C1 + C2, data = popdat2)$residuals
resU2 <- lm(U ~ R + X1 + X2 + X3 + M + C1 + C2, data = popdat2)$residuals
pcor_YU <- cor(resY, resU2)

ini <- rem <- red <- matrix(NA, n_iter, 1)
n_sample <- 2000
for (i in seq_len(n_iter)) {
  
  samp <- popdat2 %>%
    sample_n(n_sample, replace = FALSE)
  samp$R <- as.factor(samp$R)
  fit.m <- lm(M ~ R + C1 + C2, data = samp)
  fit.y <- lm(Y ~ R + X1 + X2 + X3 + M + C1 + C2, data = samp)
  fit <- smi(fit.m = fit.m, fit.y = fit.y, sims = 200, conf.level = .95,
             covariates = Cs_all, group = "R")
  
  sens.res <- sens.for.se(boot.res = fit, fit.y = fit.y, fit.m = fit.m, mediators = "M",
                          covariates = Cs_all, treat = "R", sel.lev.treat = "1",
                          ry = pcor_YU^2, rm = pcor_MU^2)
  
  rem[i] <- sens.res[6]
  red[i] <- sens.res[5]
  ini[i] <- rem[i] + red[i]
}

summary_df <- data.frame(
  Method      = "CDA.adj",
  Effect      = c("Red", "Rem", "Ini"),
  Mean        = c(mean(red), mean(rem), mean(ini)),
  `2.5%`      = unname(c(quantile(red, .025), quantile(rem, .025), quantile(ini, .025))),
  `97.5%`     = unname(c(quantile(red, .975), quantile(rem, .975), quantile(ini, .975))),
  Scenario    = 6,
  True        = unname(truth_s6_cda[c("Red", "Rem", "Ini")]),
  SD          = c(sd(red), sd(rem), sd(ini)),
  check.names = FALSE
)
print(summary_df)
write_tbl(summary_df, "scenario6_adjusted_CDA")

table1 <- rbind(res_no, res_C, res_X, res1)
table2 <- res2
table3 <- rbind(res3, summary_df[, names(res3)])
write_tbl(truths, "truths_for_tables")
write_tbl(truths_method, "truths_method_specific")
write_tbl(table1, "table1_scenarios1to4")
write_tbl(table2, "table2_scenario5")
write_tbl(table3, "table3_scenario6")

save(res_no, res_C, res_X, res1, res2, res3, summary_df,
     truths, truths_method, table1, table2, table3,
     truth_s1_cda, truth_s1_dic, truth_s1_kob,
     truth_s2_cda, truth_s2_dic, truth_s2_kob,
     truth_s3_cda, truth_s3_dic, truth_s3_kob,
     truth_s4_cda, truth_s4_dic, truth_s4_kob,
     truth_s5_cda, truth_s5_dic, truth_s5_kob,
     truth_s6_cda, truth_s6_dic, truth_s6_kob,
     pcor_MU, pcor_YU, red, rem, ini,
     file = file.path(out_dir, "sim_rv_tables.RData"))

save.image(file = file.path(out_dir, "sim_rv.RData"))
