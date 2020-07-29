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
if(!require(viridis)) install.packages("viridis")
library(viridis)
if(!require(yarrr)) install.packages("yarrr")
library(yarrr)

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
                    log_lik = TRUE, chains = 3, cores = 3, iter = 1500, warmup = 500,
                    control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_BH
print(fit_BH, probs = c(0.025,0.5,0.975),
      pars = c("p_M","q_M","s_MS","p_MS","q_MS","q_GR","M","S","R","B_rate_all","LL"), 
      include = FALSE)
## @knitr

launch_shinystan(fit_BH)

# Ricker
## @knitr fit_Ricker
fit_Ricker <- salmonIPM(fish_data, stan_model = "IPM_SMaS_np", SR_fun = "Ricker", 
                        log_lik = TRUE, chains = 3, cores = 3, iter = 1500, warmup = 500,
                        control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_Ricker
print(fit_Ricker, probs = c(0.025,0.5,0.975),
      pars = c("p_M","q_M","s_MS","p_MS","q_MS","q_GR","M","S","R","B_rate_all","LL"), 
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
M0 <- matrix(NA, nrow(M), ncol(M))/1000
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
     ylim = range(0, M0_obs, apply(M0, 2, quantile, 0.975, na.rm = TRUE), na.rm = T)*1.02)

points(apply(S, 2, median), apply(M0, 2, median), pch = 16, col = c_est)
arrows(S_obs, M0_obs, apply(S, 2, median), apply(M0, 2, median, na.rm = TRUE), 
       col = c_arr, length = 0.1)
segments(x0 = apply(S, 2, quantile, 0.025, na.rm = TRUE), 
         y0 = apply(M0, 2, median, na.rm = TRUE), 
         x1 = apply(S, 2, quantile, 0.975), col = c_est)
segments(x0 = apply(S, 2, median), 
         y0 = apply(M0, 2, quantile, 0.025, na.rm = TRUE), 
         y1 = apply(M0, 2, quantile, 0.975, na.rm = TRUE), col = c_est)

lines(Smat[1,], apply(M0hat, 2, median, na.rm = TRUE), lwd = 3, col = c_sr)
polygon(c(Smat[1,], rev(Smat[1,])), 
        c(apply(M0hat, 2, quantile, 0.025, na.rm = TRUE), 
          rev(apply(M0hat, 2, quantile, 0.975, na.rm = TRUE))), 
        col = c_srci, border = NA)
lg <- legend("topright", legend = c("observations", "states (95% CI)", "fit (95% CI)"), 
             pch = c(16,16,NA), col = c(c_obs, c_est, c_sr), lty = c(NA,1,1), lwd = c(NA,1,3),
             bty = "n")
legend("topright", legend = c("","",""), pch = c(NA,3,NA), col = c(NA, c_est, c_srci),
       lty = c(NA,NA,1), lwd = c(NA,NA,20), bty = "n", inset = c(0.19,0))

rm(list = c("mod_name","S","alpha","Rmax","M0_obs","M0","S_obs","Smat","n_Mage_obs",
            "q_M_obs","q_M","M0hat","c_obs","c_sr","c_est","c_srci","c_arr",
            "m2","m3","lg"))
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
plot(year, apply(M, 2, median), type = "l", lwd = 3, col = c_st, 
     las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n", 
     ylim = range(0, M_obs, apply(M_obs_IPM, 2, quantile, 0.975), na.rm = TRUE),
     xlab = "", ylab = "Smolts (thousands)")
polygon(c(year, rev(year)), 
        c(apply(M, 2, quantile, 0.025), rev(apply(M, 2, quantile, 0.975))),
        col = c_stci, border = NA)
polygon(c(year, rev(year)), 
        c(apply(M_obs_IPM, 2, quantile, 0.025), rev(apply(M_obs_IPM, 2, quantile, 0.975))),
        col = c_obsci, border = NA)
points(year, M_obs, pch = 16, cex = 1.2, col = c_obs)
axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.2)
rug(year[year %% 10 != 0], ticksize = -0.01)
rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.04)
legend("top", horiz = TRUE, text.width = c(5.5,4,7,7),
       legend = c("observations","states","process error","observation error"),
       cex = 0.8, pch = c(16,NA,NA,NA), pt.cex = 1.2, lty = c(NA,1,1,1), lwd = c(NA,3,10,10),
       col = c(c_obs, c_st, transparent(c_st, trans.val = 0.6), c_obsci), bty = "n")

# Spawners
plot(year, apply(S, 2, median), type = "l", lwd = 3, col = c_st, 
     las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n",
     ylim = range(0, S_obs, apply(S_obs_IPM, 2, quantile, 0.975), na.rm = TRUE),
     xlab = "Year", ylab = "")
mtext("Spawners", side = 2, line = 3.5, cex = par("cex")*1.5)
polygon(c(year, rev(year)), 
        c(apply(S, 2, quantile, 0.025), rev(apply(S, 2, quantile, 0.975))),
        col = c_stci, border = NA)
polygon(c(year, rev(year)), 
        c(apply(S_obs_IPM, 2, quantile, 0.025), rev(apply(S_obs_IPM, 2, quantile, 0.975))),
        col = c_obsci, border = NA)
points(year, S_obs, pch = 16, cex = 1.2, col = c_obs)
axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.2)
rug(year[year %% 10 != 0], ticksize = -0.01)
rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.04)

rm(list = c("mod_name","year","S_obs","M_obs","S","M","c_obs","c_st","c_stci","c_obsci"))
## @knitr
dev.off()


# #-------------------------------------------------------------------------
# # Time series of observed and estimated proportion jacks
# #-------------------------------------------------------------------------
# 
# mod_name <- "fit_Ricker"
# 
# dev.new(width = 7, height = 5)
# # png(filename=here("analysis","results",paste0("p_jack",mod_name,".png")),
# #     width=7, height=5, units="in", res=200, type="cairo-png")
# 
# ## @knitr plot_p_jack_timeseries
# year <- fish_data$year
# p_jack_obs <- run_recon(fish_data)$p_age2_obs
# p_jack_IPM <- do.call(extract1, list(as.name(mod_name),"q_MS"))[,,1]
# 
# c_obs <- transparent("orangered3", trans.val = 0.3)
# c_sr <- "blue4"
# c_srci <- transparent(c_sr, trans.val = 0.8)
# 
# par(mar = c(4.5, 4.5, 0.5, 0.5))
# plot(year, apply(p_jack_IPM, 2, median), type = "l", lwd = 3, col = c_sr, 
#      las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n",
#      ylim = range(p_jack_obs, apply(p_jack_IPM, 2, quantile, 0.975)),
#      xlab = "Year", ylab = "Proportion jacks")
# polygon(c(year, rev(year)), 
#         c(apply(p_jack_IPM, 2, quantile, 0.025), rev(apply(p_jack_IPM, 2, quantile, 0.975))),
#         col = c_srci, border = NA)
# points(year[fish_data$obs_type=="past"], p_jack_obs[fish_data$obs_type=="past"], 
#        pch = 16, cex = 1.2, col = c_obs)
# axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.2)
# rug(year[year %% 10 != 0], ticksize = -0.01)
# rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.02)
# 
# rm(list = c("mod_name","year","p_jack_obs","p_jack_IPM","c_obs","c_sr","c_srci"))
# ## @knitr
# # dev.off()

