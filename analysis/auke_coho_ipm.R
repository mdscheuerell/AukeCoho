#===========================================================================
# SETUP
#===========================================================================

## @knitr setup
options(device = ifelse(.Platform$OS.type == "windows", "windows", "quartz"))
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)

if(!require(rstan)) install.packages("rstan")
library(rstan) 
if(!require(salmonIPM)) devtools::install_github("ebuhle/salmonIPM")
library(salmonIPM)
if(!require(here)) install.packages("here")
library(here) 
if(!require(matrixStats)) install.packages("matrixStats")
library(matrixStats)
if(!require(yarrr)) install.packages("yarrr")
library(yarrr)
if(!require(Hmisc)) install.packages("Hmisc")
library(Hmisc)
if(!require(viridis)) install.packages("viridis")
library(viridis)

if(file.exists(here("analysis","results","auke_coho_ipm.RData")))
  load(here("analysis","results","auke_coho_ipm.RData"))
## @knitr

#===========================================================================
# DATA
#===========================================================================

## @knitr data

# Coho population data
fish_data <- read.csv(here("data","auke_coho_data_1980-2019.csv"))

# Environmental covariate data
cov_raw <- read.csv(here("data","covariates_1980-2019.csv"))

# Standardize covariates for modeling
env_data <- cov_raw
env_data[,-1] <- scale(env_data[,-1])

## @knitr

#===========================================================================
# FIT MODELS
#===========================================================================

# Density-independent
## @knitr fit_exp
fit_exp <- salmonIPM(fish_data, stan_model = "IPM_SMaS_np", SR_fun = "exp", 
                     pars = "Rmax", include = FALSE, log_lik = TRUE, 
                     chains = 3, cores = 3, iter = 1500, warmup = 500,
                     control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_exp
print(fit_exp,  probs = c(0.025,0.5,0.975),
      pars = c("p_M","q_M","s_MS","p_MS","q_MS","q_GR","M","S","R","B_rate_all","LL"), 
      include = FALSE)
## @knitr

launch_shinystan(fit_exp)

# Beverton-Holt
## @knitr fit_BH
fit_BH <- salmonIPM(fish_data, stan_model = "IPM_SMaS_np", SR_fun = "BH", 
                    pars = c(stan_pars("IPM_SMaS_np"), "epsilon_M"), log_lik = TRUE, 
                    chains = 3, cores = 3, iter = 1500, warmup = 500,
                    control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_BH
print(fit_BH, probs = c(0.025,0.5,0.975),
      pars = c("epsilon_M","p_M","q_M","s_MS","p_MS","q_MS","q_GR","M","S","R","B_rate_all","LL"), 
      include = FALSE)
## @knitr

launch_shinystan(fit_BH)

# Ricker
## @knitr fit_Ricker
fit_Ricker <- salmonIPM(fish_data, stan_model = "IPM_SMaS_np", SR_fun = "Ricker", 
                        pars = c(stan_pars("IPM_SMaS_np"), "epsilon_M"), log_lik = TRUE, 
                        chains = 3, cores = 3, iter = 1500, warmup = 500,
                        control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_Ricker
print(fit_Ricker, probs = c(0.025,0.5,0.975),
      pars = c("epsilon_M","p_M","q_M","s_MS","p_MS","q_MS","q_GR","M","S","R","B_rate_all","LL"), 
      include = FALSE)
## @knitr

launch_shinystan(fit_Ricker)

#--------------------------------------------------------------
# Model selection using LOO
#--------------------------------------------------------------

# Observationwise log-likelihood of each fitted model
# Here an observation is a row of fish_data, and the total likelihood includes 
# components for smolt abundance, spawner abundance, smolt age-frequency,
# ocean age-frequency, and adult (Gilbert-Rich) age-frequency
## @knitr loo
LL <- lapply(list(exp = fit_exp, BH = fit_BH, Ricker = fit_Ricker),
             loo::extract_log_lik, parameter_name = "LL", merge_chains = FALSE)

# Relative ESS of posterior draws of observationwise likelihood 
r_eff <- lapply(LL, function(x) relative_eff(exp(x)))

# PSIS-LOO
LOO <- lapply(1:length(LL), function(i) loo(LL[[i]], r_eff = r_eff[[i]]))
names(LOO) <- names(LL)

## Compare all three models
loo_compare(LOO)

## Exponential vs. Ricker
loo_compare(LOO[c("exp","Ricker")])

## Exponential vs. Beverton-Holt
loo_compare(LOO[c("exp","BH")])

## Beverton-Holt vs. Ricker
loo_compare(LOO[c("BH","Ricker")])
## @knitr


#--------------------------------------------------------------
# Save stanfit objects
#--------------------------------------------------------------

save(list = ls()[sapply(ls(), function(x) do.call(class, list(as.name(x)))) == "stanfit"], 
     file = here("analysis","results","auke_coho_ipm.RData"))


#===========================================================================
# FIGURES
#===========================================================================

#-------------------------------------------------------------------------
# S-R curve overlaid with observations and states
#-------------------------------------------------------------------------

mod_name <- "fit_Ricker"

# dev.new(width = 7, height = 7)
png(filename=here("analysis","results",paste0("SR_",mod_name,".png")),
    width=7, height=7, units="in", res=200, type="cairo-png")

## @knitr plot_SR
# observations
S_obs <- fish_data$S_obs
M_obs <- fish_data$M_obs/1000
n_Mage_obs <- stan_data(fish_data, stan_model = "IPM_SMaS_np")$n_Mage_obs
q_M_obs <- sweep(n_Mage_obs, 1, rowSums(n_Mage_obs), "/")
N <- nrow(fish_data)
M0_obs <- rep(NA,N)
for(i in 1:N)
  M0_obs[i] <- ifelse((i + 2) <= N, M_obs[i+2]*q_M_obs[i+2,1], NA) + 
               ifelse((i + 3) <= N, M_obs[i+3]*q_M_obs[i+3,2], NA)
# states
S <- do.call(extract1, list(as.name(mod_name), "S"))
M <- do.call(extract1, list(as.name(mod_name), "M"))/1000
M0 <- matrix(NA, nrow(M), ncol(M))
q_M <- do.call(extract1, list(as.name(mod_name), "q_M"))
for(i in 1:N) {
  if((i + 2) <= N) m2 <- M[,i+2]*q_M[,i+2,1] else m2 <- NA
  if((i + 3) <= N) m3 <- M[,i+3]*q_M[,i+3,2] else m3 <- NA
  M0[,i] <- m2 + m3
}
# predicted states
SR <- function(alpha, Rmax, S, A, SR_fun) 
{
  switch(SR_fun, 
         exp = alpha*S/A,
         BH = alpha*S/(A + alpha*S/Rmax),
         Ricker = alpha*(S/A)*exp(-alpha*S/(A*exp(1)*Rmax)))
}
SR_fun <- unlist(strsplit(mod_name, "_"))[2]
alpha <- as.vector(do.call(extract1, list(as.name(mod_name), "alpha")))
Rmax <- as.vector(do.call(extract1, list(as.name(mod_name), "Rmax")))
Smat <- matrix(seq(0, max(S_obs, apply(S, 2, quantile, 0.975), na.rm = T)*1.02, length = 500),
               nrow = length(alpha), ncol = 500, byrow = TRUE)
M0hat <- SR(alpha = alpha, Rmax = Rmax, S = Smat, A = 1, SR_fun)/1000

c_obs <- transparent("orangered3", trans.val = 0.3)
c_sr <- "blue4"
c_est <- transparent(c_sr, trans.val = 0.5)
c_srci <- transparent(c_sr, trans.val = 0.8)
c_arr <- "darkgray"

par(mar = c(5,4.5,2,1))
plot(S_obs, M0_obs, pch = 16, col = c_obs, las = 1,
     cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5, xaxs = "i", yaxs = "i",
     xlab = "Spawners", ylab = "Smolts (thousands)", xlim = c(0, max(Smat)),
     ylim = range(0, M0_obs, colQuantiles(M0, probs = 0.975, na.rm = TRUE), na.rm = TRUE)*1.02)

points(colMedians(S), colMedians(M0), pch = 16, col = c_est)
arrows(S_obs, M0_obs, colMedians(S), colMedians(M0, na.rm = TRUE), 
       col = c_arr, length = 0.1)
segments(x0 = colQuantiles(S, probs = 0.025, na.rm = TRUE), 
         y0 = colMedians(M0, na.rm = TRUE), 
         x1 = colQuantiles(S, probs = 0.975), col = c_est)
segments(x0 = colMedians(S), 
         y0 = colQuantiles(M0, probs = 0.025, na.rm = TRUE), 
         y1 = colQuantiles(M0, probs = 0.975, na.rm = TRUE), col = c_est)

lines(Smat[1,], colMedians(M0hat, na.rm = TRUE), lwd = 3, col = c_sr)
polygon(c(Smat[1,], rev(Smat[1,])), 
        c(colQuantiles(M0hat, probs = 0.025, na.rm = TRUE), 
          rev(colQuantiles(M0hat, probs = 0.975, na.rm = TRUE))), 
        col = c_srci, border = NA)
legend("topright", legend = c("observations", "states (95% CI)", "fit (95% CI)"), 
       pch = c(16,16,NA), col = c(c_obs, c_est, c_sr), lty = c(NA,1,1), lwd = c(NA,1,3),
       bty = "n")
legend("topright", legend = c("","",""), pch = c(NA,3,NA), col = c(NA, c_est, c_srci),
       lty = c(NA,NA,1), lwd = c(NA,NA,20), bty = "n", inset = c(0.19,0))

rm(list = c("mod_name","S","alpha","Rmax","M0_obs","M0","S_obs","Smat","n_Mage_obs",
            "q_M_obs","q_M","M0hat","c_obs","c_sr","c_est","c_srci","c_arr",
            "m2","m3"))
## @knitr
dev.off()


#-------------------------------------------------------------------------
# Posterior distributions (and priors) of S-R parameters
#-------------------------------------------------------------------------

mod_name <- "fit_Ricker"

# dev.new(width = 10, height = 5)
png(filename=here("analysis","results",paste0("SR_params_",mod_name,".png")),
    width=10, height=5, units="in", res=200, type="cairo-png")

## @knitr plot_SR_params
alpha <- do.call(extract1, list(as.name(mod_name),"alpha"))
Rmax <- do.call(extract1, list(as.name(mod_name),"Rmax"))/1000

par(mfrow = c(1,2), mar = c(5.1,2,3,1), oma = c(0,1,0,0))

# Posterior of log(alpha)
hist(log(alpha), 15, prob = TRUE,  border = "white",
     las = 1, cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5, yaxs = "i", yaxt = "n",
     xlab = bquote(log(alpha)), ylab = "", main = "Intrinsic productivity")
curve(dnorm(x,2,2), add = TRUE)
box(bty = "l")
mtext("Probability density", side = 2, line = 1, cex = par("cex")*1.5)

# Posterior of Rmax
hist(Rmax, 15, prob = TRUE,  border = "white",
     las = 1, cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5, yaxs = "i", yaxt = "n",
     xlab = bquote(italic(R)[max] * " (thousands)"), ylab = "", main = "Maximum smolts")
curve(dlnorm(x,2,3), add = TRUE)
box(bty = "l")

rm(list = c("mod_name","alpha","Rmax"))

## @knitr
dev.off()


#-------------------------------------------------------------------------
# Time series of observed and estimated smolts and spawners
#-------------------------------------------------------------------------

mod_name <- "fit_Ricker"

# dev.new(width = 7, height = 7)
png(filename=here("analysis","results",paste0("M_S_",mod_name,".png")),
    width=7, height=7, units="in", res=200, type="cairo-png")

## @knitr plot_smolt_spawner_timeseries
year <- fish_data$year
M_obs <- fish_data$M_obs/1000
S_obs <- fish_data$S_obs
M <- do.call(extract1, list(as.name(mod_name),"M"))/1000
tau_M <- do.call(extract1, list(as.name(mod_name),"tau_M"))
M_obs_IPM <- M*rlnorm(length(M), 0, tau_M)
S <- do.call(extract1, list(as.name(mod_name),"S"))
tau_S <- do.call(extract1, list(as.name(mod_name),"tau_S"))
S_obs_IPM <- S*rlnorm(length(S), 0, tau_S)

c_obs <- transparent("orangered3", trans.val = 0.3)
c_st <- "blue4"
c_stci <- transparent(c_st, trans.val = 0.8)
c_obsci <- transparent(c_st, trans.val = 0.8)

par(mfrow = c(2,1), mar = c(4.5, 5.1, 0.5, 0.5))

# Smolts
plot(year, colMedians(M), type = "l", lwd = 3, col = c_st, 
     las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n", 
     ylim = range(0, M_obs, colQuantiles(M_obs_IPM, probs = 0.975), na.rm = TRUE),
     xlab = "", ylab = "Smolts (thousands)")
polygon(c(year, rev(year)), 
        c(colQuantiles(M, probs = 0.025), rev(colQuantiles(M, probs = 0.975))),
        col = c_stci, border = NA)
polygon(c(year, rev(year)), 
        c(colQuantiles(M_obs_IPM, probs = 0.025), rev(colQuantiles(M_obs_IPM, probs = 0.975))),
        col = c_obsci, border = NA)
points(year, M_obs, pch = 16, cex = 1.5, col = c_obs)
axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.2)
rug(year[year %% 10 != 0], ticksize = -0.01)
rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.04)
legend("top", horiz = TRUE, text.width = c(5.5,4,7,7),
       legend = c("observations","states","process error","observation error"),
       cex = 0.8, pch = c(16,NA,NA,NA), pt.cex = 1.5, lty = c(NA,1,1,1), lwd = c(NA,3,10,10),
       col = c(c_obs, c_st, transparent(c_st, trans.val = 0.6), c_obsci), bty = "n")

# Spawners
plot(year, colMedians(S), type = "l", lwd = 3, col = c_st, 
     las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n",
     ylim = range(0, S_obs, colQuantiles(S_obs_IPM, probs = 0.975), na.rm = TRUE),
     xlab = "Year", ylab = "")
mtext("Spawners", side = 2, line = 3.5, cex = par("cex")*1.5)
polygon(c(year, rev(year)), 
        c(colQuantiles(S, probs = 0.025), rev(colQuantiles(S, probs = 0.975))),
        col = c_stci, border = NA)
polygon(c(year, rev(year)), 
        c(colQuantiles(S_obs_IPM, probs = 0.025), rev(colQuantiles(S_obs_IPM, probs = 0.975))),
        col = c_obsci, border = NA)
points(year, S_obs, pch = 16, cex = 1.5, col = c_obs)
axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.2)
rug(year[year %% 10 != 0], ticksize = -0.01)
rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.04)

rm(list = c("mod_name","year","S_obs","M_obs","S","M","c_obs","c_st","c_stci","c_obsci"))
## @knitr
dev.off()


#------------------------------------------------------------------------
# Time series of observed and estimated smolt and spawner age structure
#------------------------------------------------------------------------

mod_name <- "fit_Ricker"

# dev.new(width = 5, height = 7)
png(filename=here("analysis","results",paste0("q_",mod_name,".png")),
    width=5, height=7, units="in", res=200, type="cairo-png")

## @knitr plot_smolt_spawner_age_timeseries
year <- fish_data$year
dat <- stan_data(fish_data, stan_model = "IPM_SMaS_np")
n_Mage_obs <- dat$n_Mage_obs
q_Mage2_obs <- binconf(n_Mage_obs[,"n_Mage2_obs"], rowSums(n_Mage_obs), alpha = 0.05)
n_MSage_obs <- dat$n_MSage_obs
q_MSage0_obs <- binconf(n_MSage_obs[,"n_MSage0_obs"], rowSums(n_MSage_obs), alpha = 0.05)
n_GRage_obs <- dat$n_GRage_obs

q_M <- do.call(extract1, list(as.name(mod_name), "q_M"))
q_MS <- do.call(extract1, list(as.name(mod_name), "q_MS"))
q_GR <- do.call(extract1, list(as.name(mod_name), "q_GR"))

cobs <- transparent("orangered3", trans.val = 0.3)
cst <- "blue4"
cstci <- transparent(cst, trans.val = 0.8)
cGR <- viridis(ncol(n_GRage_obs), end = 0.9) 
cGRt <- transparent(cGR, trans.val = 0.2)
cGRtt <- transparent(cGR, trans.val = 0.7)

par(mfrow = c(3,1), mar = c(4.5, 5.1, 0.5, 1))

# proportion age-2 smolts
plot(year, colMedians(q_M[,,1]), type = "l", lwd = 2, col = cst, ylim = c(0,1), 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, xaxt = "n", 
     xlab = "", ylab = "Proportion age-2 smolts")
axis(side = 1, at = year[year %% 5 == 0], cex.axis = 1.2)
rug(year[year %% 5 != 0], ticksize = -0.02)
polygon(c(year, rev(year)),
        c(colQuantiles(q_M[,,1], probs = 0.025), 
          rev(colQuantiles(q_M[,,1], probs = 0.975))),
        col = cstci, border = NA)
points(year, q_Mage2_obs[,"PointEst"], pch = 16, col = cobs, cex = 1.5)
segments(x0 = year, y0 = q_Mage2_obs[,"Lower"], y1 = q_Mage2_obs[,"Upper"], col = cobs)
legend("top", horiz = TRUE, legend = c("observations (95% CI)   ","states (95% CI)"),
       text.col = "white", pch = c(16,NA), pt.cex = 1.5, lwd = c(1,2), col = c(cobs, cst), 
       bty = "n")
legend("top", horiz = TRUE, legend = c("observations (95% CI)   ","states (95% CI)"),
       pch = NA, lty = c(NA,1), lwd = c(NA,10), col = c(NA, cstci), bty = "n")

# proportion jacks
plot(year, colMedians(q_MS[,,1]), type = "l", lwd = 2, col = cst, 
     ylim = range(0, q_MSage0_obs[,"Upper"], colQuantiles(q_MS[,,1], probs = 0.975)), 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, xaxt = "n", 
     xlab = "", ylab = "Proportion jacks")
axis(side = 1, at = year[year %% 5 == 0], cex.axis = 1.2)
rug(year[year %% 5 != 0], ticksize = -0.02)
polygon(c(year, rev(year)),
        c(colQuantiles(q_MS[,,1], probs = 0.025), 
          rev(colQuantiles(q_MS[,,1], probs = 0.975))),
        col = cstci, border = NA)
points(year, q_MSage0_obs[,"PointEst"], pch = 16, col = cobs, cex = 1.5)
segments(x0 = year, y0 = q_MSage0_obs[,"Lower"], y1 = q_MSage0_obs[,"Upper"], col = cobs)

# Gilbert-Rich age proportions
plot(year, rep(0.5, length(year)), type = "n", ylim = c(0,0.9), 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, xaxt = "n", 
     xlab = "", ylab = "Proportion at age")
axis(side = 1, at = year[year %% 5 == 0], cex.axis = 1.2)
rug(year[year %% 5 != 0], ticksize = -0.02)
for(a in 1:ncol(n_GRage_obs))
{
  q_obs <- binconf(n_GRage_obs[,a], rowSums(n_GRage_obs), alpha = 0.05)
  lines(year, colMedians(q_GR[,,a]), col = cGRt[a], lwd = 2)
  polygon(c(year, rev(year)),
          c(colQuantiles(q_GR[,,a], probs = 0.025),
            rev(colQuantiles(q_GR[,,a], probs = 0.975))),
          col = cGRtt[a], border = NA)
  points(year, q_obs[,"PointEst"], pch = 16, col = cGRt[a], cex = 1.5)
  segments(x0 = year, y0 = q_obs[,"Lower"], y1 = q_obs[,"Upper"], col = cGRt[a])
}
legend("top", horiz = TRUE, title = "Gilbert-Rich age", x.intersp = 0.5,
       legend = parse(text = paste0(substring(colnames(n_GRage_obs), 8, 8), "[",
                                    substring(colnames(n_GRage_obs), 10, 10), "]")),
       col = cGRt, pch = 16, pt.cex = 1.5, lwd = 1, xjust = 0.5, bty = "n")

rm(list = c("mod_name","year","n_Mage_obs","q_Mage2_obs","n_MSage_obs","q_MSage0_obs",
            "n_GRage_obs","q_M","q_MS","q_GR","q_obs","cobs","cst","cstci","cGR","cGRt","cGRtt"))
## @knitr
dev.off()


#-----------------------------------------------------------------------------------
# Time series of smolt recruitment process errors and age proportions by brood year
#-----------------------------------------------------------------------------------

mod_name <- "fit_Ricker"

# dev.new(width = 7, height = 7)
png(filename=here("analysis","results",paste0("p_M_",mod_name,".png")),
    width=7, height=7, units="in", res=200, type="cairo-png")

## @knitr plot_smolt_recruitment_timeseries
year <- fish_data$year
# X_M <- env_data
# beta_M <- do.call(extract1, list(as.name(mod_name),"beta_M"))
epsilon_M <- do.call(extract1, list(as.name(mod_name),"epsilon_M"))
p_M <- do.call(extract1, list(as.name(mod_name),"p_M"))

c_st <- "blue4"
c_stci <- transparent(c_st, trans.val = 0.8)

par(mfrow = c(2,1), mar = c(4.5, 5.1, 0.5, 0.5))

# smolt recruitment process errors
plot(year, colMedians(epsilon_M), type = "n", 
     las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n", 
     ylim = range(colQuantiles(epsilon_M, probs = c(0.025,0.975))),
     xlab = "", ylab = "Recruitment anomaly")
abline(h = 0, lty = 2, lwd = 0.5)
lines(year, colMedians(epsilon_M), lwd = 3, col = c_st)
polygon(c(year, rev(year)), 
        c(colQuantiles(epsilon_M, probs = 0.025), 
          rev(colQuantiles(epsilon_M, probs = 0.975))),
        col = c_stci, border = NA)
axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.2)
rug(year[year %% 10 != 0], ticksize = -0.01)
rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.04)

# smolt age composition
plot(year, rep(0.5, length(year)), type = "n", 
     las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n", ylim = range(0,1), 
     xlab = "Brood year", ylab = "Proportion age-2 smolts")
abline(h = 0.5, lty = 2, lwd = 0.5)
lines(year, colMedians(p_M[,,1]), lwd = 3, col = c_st)
polygon(c(year, rev(year)), 
        c(colQuantiles(p_M[,,1], probs = 0.025), rev(colQuantiles(p_M[,,1], probs = 0.975))),
        col = c_stci, border = NA)
axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.2)
rug(year[year %% 10 != 0], ticksize = -0.01)
rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.04)

rm(list = c("mod_name","year","epsilon_M","c_st","c_stci","p_M"))
## @knitr
dev.off()


#--------------------------------------------------------------------------------------
# Time series of SAR and ocean-age proportions for each smolt age, by outmigration year
#--------------------------------------------------------------------------------------

mod_name <- "fit_Ricker"

# dev.new(width = 7, height = 7)
png(filename=here("analysis","results",paste0("SAR-p_MS_",mod_name,".png")),
    width=7, height=7, units="in", res=200, type="cairo-png")

## @knitr plot_SAR_jack_timeseries
year <- fish_data$year
s_MS <- do.call(extract1, list(as.name(mod_name),"s_MS"))
p_MS <- do.call(extract1, list(as.name(mod_name),"p_MS"))


c2 <- viridis(5)[2]
c2t <- transparent(c2, trans.val = 0.6)
c3 <- viridis(5)[4]
c3t <- transparent(c3, trans.val = 0.6)

par(mfrow = c(2,1), mar = c(4.5, 5.1, 0.5, 0.5))

# SAR
plot(year, colMedians(s_MS[,,1]), type = "l", lwd = 3, col = c2, 
     las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n", 
     ylim = range(0, apply(s_MS, 2:3, quantile, 0.975)),
     xlab = "", ylab = "Smolt-to-adult survival")
polygon(c(year, rev(year)), 
        c(colQuantiles(s_MS[,,1], probs = 0.025), 
          rev(colQuantiles(s_MS[,,1], probs = 0.975))),
        col = c2t, border = NA)
lines(year, colMedians(s_MS[,,2]), lwd = 3, col = c3)
polygon(c(year, rev(year)), 
        c(colQuantiles(s_MS[,,2], probs = 0.025), 
          rev(colQuantiles(s_MS[,,2], probs = 0.975))),
        col = c3t, border = NA)
axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.2)
rug(year[year %% 10 != 0], ticksize = -0.01)
rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.04)
legend("topright", title = "smolt age", legend = 2:3, 
       lty = 1, lwd = 3, col = c(c2,c3), bty = "n")

# Jack proportions
plot(year, colMedians(p_MS[,,1,1]), type = "l", lwd = 3, col = c2, 
     las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n", 
     ylim = range(0, apply(p_MS[,,,1], 2:3, quantile, 0.975)), 
     xlab = "Outmigration year", ylab = "Proportion jacks")
polygon(c(year, rev(year)), 
        c(colQuantiles(p_MS[,,1,1], probs = 0.025), 
          rev(colQuantiles(p_MS[,,1,1], probs = 0.975))),
        col = c2t, border = NA)
lines(year, colMedians(p_MS[,,2,1]), lwd = 3, col = c3)
polygon(c(year, rev(year)), 
        c(colQuantiles(p_MS[,,2,1], probs = 0.025), 
          rev(colQuantiles(p_MS[,,2,1], probs = 0.975))),
        col = c3t, border = NA)
axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.2)
rug(year[year %% 10 != 0], ticksize = -0.01)
rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.04)

rm(list = c("mod_name","year","s_MS","p_MS","c2","c2t","c3","c3t"))
## @knitr
dev.off()


#--------------------------------------------------------
# Fitted vs. observed catch (modeled as broodstock take)
#--------------------------------------------------------

mod_name <- "fit_Ricker"

S <- do.call(extract1, list(as.name(mod_name), "S"))
q_MS <- do.call(extract1, list(as.name(mod_name), "q_MS"))
B_rate <- do.call(extract1, list(as.name(mod_name), "B_rate_all"))
B_take <- B_rate*S*q_MS[,,2]/(1 - B_rate) 

dev.new()
par(mar = c(5,5,1,1))
plot(colMedians(B_take), fish_data$B_take_obs, pch = 1, cex = 1.5, las = 1,
     xlim = range(colQuantiles(B_take, probs = c(0.025, 0.975))),
     cex.axis = 1.2, cex.lab = 1.5, xlab = "Observed catch", ylab = "Estimated catch")
segments(x0 = colQuantiles(B_take, probs = 0.025), x1 = colQuantiles(B_take, probs = 0.975),
         y0 = fish_data$B_take_obs)
abline(0,1)

rm(list = c('mod_name','S','q_MS','B_rate','B_take'))


#---------------------------------------------------------------------
# Smolt recruitment process errors vs. jack proportion of spawners
#---------------------------------------------------------------------

mod_name <- "fit_Ricker"

epsilon_M <- do.call(extract1, list(as.name(mod_name),"epsilon_M"))
q_MS <- do.call(extract1, list(as.name(mod_name),"q_MS"))
dat <- stan_data(fish_data, stan_model = "IPM_SMaS_np")
n_MSage_obs <- dat$n_MSage_obs
q_MSage0_obs <- binconf(n_MSage_obs[,"n_MSage0_obs"], rowSums(n_MSage_obs), alpha = 0.05)

dev.new(width = 10, height = 5)
par(mfrow = c(1,2), mar = c(5,5,1,1))

# use estimated jack fraction (states)
plot(colMedians(q_MS[,,1]), colMedians(epsilon_M), pch = 16, cex = 1.5, 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, 
     xlim = range(colQuantiles(q_MS[,,1], probs = c(0.025,0.975))),
     ylim = range(colQuantiles(epsilon_M, probs = c(0.025,0.975))),
     xlab = "Estimated proportion jacks", ylab = "Smolt recruitment anomaly")
segments(x0 = colQuantiles(q_MS[,,1], probs = 0.025), 
         x1 = colQuantiles(q_MS[,,1], probs = 0.975),
         y0 = colMedians(epsilon_M))
segments(x0 = colMedians(q_MS[,,1]), 
         y0 = colQuantiles(epsilon_M, probs = 0.025), 
         y1 = colQuantiles(epsilon_M, probs = 0.975))

# use observed jack fraction
plot(q_MSage0_obs[,"PointEst"], colMedians(epsilon_M), pch = 16, cex = 1.5, 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, 
     xlim = range(q_MSage0_obs[,"Lower"], q_MSage0_obs[,"Upper"]),
     ylim = range(colQuantiles(epsilon_M, probs = c(0.025,0.975))),
     xlab = "Observed proportion jacks", ylab = "Smolt recruitment anomaly")
segments(x0 = q_MSage0_obs[,"Lower"], x1 = q_MSage0_obs[,"Upper"], 
         y0 = colMedians(epsilon_M))
segments(x0 = q_MSage0_obs[,"PointEst"], 
         y0 = colQuantiles(epsilon_M, probs = 0.025), 
         y1 = colQuantiles(epsilon_M, probs = 0.975))

rm(list = c('epsilon_M','q_MS','dat','mod_name','n_MSage_obs','q_MSage0_obs'))


