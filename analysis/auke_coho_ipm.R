## ----set_options, echo=FALSE, cache=TRUE---------------------------------
options(width = 100)
if(Sys.info()["sysname"] == "Windows") options(device = windows)

## ----load_pkgs, message=FALSE, warning=FALSE-----------------------------
## for inference
#unleash below to update master branch version of IPM
detach(package:salmonIPM, unload = TRUE)
devtools::install_github("ebuhle/salmonIPM@ICchinook-models")
library("salmonIPM")
## for dir management
if(!require("here")) {
  install.packages("here")
  library("here")
}
## for data munging
if(!require("readr")) {
  install.packages("readr")
  library("readr")
}
if(!require("dplyr")) {
  install.packages("dplyr")
  library("dplyr")
}
## for output/plotting
if(!require("viridis")) {
  install.packages("viridis")
  library("viridis")
}
if(!require("yarrr")) {
  install.packages("yarrr")
  library("yarrr")
}
options(tibble.print_max = Inf, tibble.width = Inf)

## ----get_data------------------------------------------------------------
## set working directory
setwd(here("analysis"))
## directory
datadir <- here("data")
## read data
datafile <- dir(datadir)[grep("auke_coho_data", dir(datadir))]
fishdata <- read_csv(file.path(datadir, datafile))
fishdata <- fishdata[fishdata$year <= 2019,]  # truncate future years
# read covariate data
pdofile <- dir(datadir)[grep("pdo", dir(datadir))]
pdodat <- read_csv(file.path(datadir, pdofile))
covfile <- dir(datadir)[grep("covariates", dir(datadir))]
cov_raw <- read_csv(file.path(datadir, covfile))
## Subset PDO to Nov-Feb and index to outmigration year
pdo <- pdodat[pdodat$month %in% c(11,12,1,2),]
pdo$winter_year <- ifelse(pdo$month %in% 11:12, pdo$year, pdo$year - 1)
pdo <- aggregate(pdo ~ winter_year, data = pdo, mean)

## Add PDO to covariate data
env_data <- cov_raw
env_data$winter_pdo <- pdo$pdo[match(env_data$year, pdo$winter_year)]
# HPC release in 100s of millions of fish
env_data$hpc_release <- (env_data$hpc_release - mean(env_data$hpc_release))/1e8
env_data$mean_stream_temp <- scale(env_data$mean_stream_temp, scale = FALSE)
env_data$winter_pdo <- scale(env_data$winter_pdo)

## Truncate fish data to non-missing covariate data
fishdata_env <- fishdata[fishdata$year %in% env_data$year,]

## ----fit_density_independent, eval=TRUE------------------------------------------------
fit_exp <- salmonIPM(fishdata, stan_model = "IPM_SMaS_np", SR_fun = "exp", 
                    pars = c(stan_pars("IPM_SMaS_np"), "LL"),
                    chains = 3, cores = 3, iter = 1500, warmup = 500,
                    control = list(adapt_delta = 0.99, max_treedepth = 13))

## ----print_density_independent---------------------------------------------------
print(fit_exp,  probs = c(0.025,0.5,0.975),
      pars = c("p_M","q_M","s_MS","p_MS","q_MS","q_GR","M","S","R","LL"), include = FALSE)

## ----shinystan_density_indepndent------------------------------------------------------------
launch_shinystan(fit_exp)


## ----fit_BH, eval=TRUE------------------------------------------------
fit_BH <- salmonIPM(fishdata, stan_model = "IPM_SMaS_np", SR_fun = "BH", 
                    pars = c(stan_pars("IPM_SMaS_np"), "LL"),
                    chains = 3, cores = 3, iter = 1500, warmup = 500,
                    control = list(adapt_delta = 0.99, max_treedepth = 13))

## ----print_fit_BH---------------------------------------------------
print(fit_BH, probs = c(0.025,0.5,0.975),
      pars = c("p_M","q_M","s_MS","p_MS","q_MS","q_GR","M","S","R","LL"), include = FALSE)

## ----shinystan_fit_BH--------------------------------------------------------
launch_shinystan(fit_BH)


## ----fit_BH, eval=TRUE------------------------------------------------
fit_BH_env <- salmonIPM(fishdata_env, 
                        env_data = list(M = env_data[,"mean_stream_temp", drop = FALSE], 
                                        MS = env_data[,c("hpc_release","winter_pdo")]), 
                        stan_model = "IPM_SMaS_np", SR_fun = "BH", 
                        pars = c(stan_pars("IPM_SMaS_np"), "LL"),
                        chains = 3, cores = 3, iter = 1500, warmup = 500,
                        control = list(adapt_delta = 0.99, max_treedepth = 13))

## ----print_fit_BH---------------------------------------------------
print(fit_BH_env, probs = c(0.025,0.5,0.975),
      pars = c("p_M","q_M","s_MS","p_MS","q_MS","q_GR","M","S","R","LL"), include = FALSE)

## ----shinystan_fit_BH--------------------------------------------------------
launch_shinystan(fit_BH_env)


## ----fit_Ricker, eval=TRUE------------------------------------------------
fit_Ricker <- salmonIPM(fishdata, stan_model = "IPM_SMaS_np", SR_fun = "Ricker", 
                    pars = c(stan_pars("IPM_SMaS_np"), "LL"),
                    chains = 3, cores = 3, iter = 1500, warmup = 500,
                    control = list(adapt_delta = 0.99, max_treedepth = 13))

## ----print_Ricker---------------------------------------------------
print(fit_Ricker, probs = c(0.025,0.5,0.975),
      pars = c("p_M","q_M","s_MS","p_MS","q_MS","q_GR","M","S","R","LL"), include = FALSE)

## ----shinystan_Ricker------------------------------------------------------------
launch_shinystan(fit_Ricker)


##---MODEL SELECTION-------------------------------------------------------
##---compare_spawner_recruit_models-----------------------------------
loo_exp <- loo(extract1(fit_exp,"LL"))
loo_BH <- loo(extract1(fit_BH,"LL"))
loo_Ricker <- loo(extract1(fit_Ricker,"LL"))
# compare all three S-R models
compare(loo_exp, loo_BH, loo_Ricker)
# pairwise comparisons
compare(loo_exp, loo_BH)
compare(loo_exp, loo_Ricker)
compare(loo_BH, loo_Ricker)



## ----FIGURES--------------------------------------------------------------

#-------------------------------------------------------------------------
# S-R curve overlaid with observations and states
#-------------------------------------------------------------------------

dev.new(width = 7, height = 7)
# png(filename="SR.png", width=7, height=7, units="in", res=200, type="cairo-png")
par(mar = c(5,5,2,1))
  
SR_fun <- "BH"
SR <- function(alpha, Rmax, S, A, SR_fun) 
{
  switch(SR_fun, 
         BH = alpha*S/(A + alpha*S/Rmax),
         Ricker = alpha*(S/A)*exp(-alpha*S/(A*exp(1)*Rmax)))
}

S_obs <- fishdata$S_obs
M_obs <- fishdata$M_obs
n_Mage_obs <- stan_data(fishdata, stan_model = "IPM_SMaS_np")$n_Mage_obs
q_M_obs <- sweep(n_Mage_obs, 1, rowSums(n_Mage_obs), "/")
N <- nrow(fishdata)
M0_obs <- rep(NA,N)
for(i in 1:N)
  M0_obs[i] <- ifelse((i + 2) <= N, M_obs[i+2]*q_M_obs[i+2,1], NA) + 
               ifelse((i + 3) <= N, M_obs[i+3]*q_M_obs[i+3,2], NA)
S_IPM <- extract1(fit_BH,"S")
M_IPM <- extract1(fit_BH,"M")
M0_IPM <- matrix(NA, nrow(M_IPM), ncol(M_IPM))
q_M_IPM <- extract1(fit_BH,"q_M")
for(i in 1:N) {
  if((i + 2) <= N) m2 <- M_IPM[,i+2]*q_M_IPM[,i+2,1] else m2 <- NA
  if((i + 3) <= N) m3 <- M_IPM[,i+3]*q_M_IPM[,i+3,2] else m3 <- NA
  M0_IPM[,i] <- m2 + m3
}

S <- matrix(seq(0, max(S_obs, apply(S_IPM, 2, quantile, 0.975), na.rm = T)*1.02, length = 500),
            nrow = sum(fit_BH@sim$n_save - fit_BH@sim$warmup2), ncol = 500, byrow = T)
alpha <- as.vector(extract1(fit_BH,"alpha"))
Rmax <- as.vector(extract1(fit_BH,"Rmax"))
M0hat_IPM <- SR(alpha = alpha, Rmax = Rmax, S = S, A = 1, SR_fun)

c_obs <- transparent("orangered3", trans.val = 0.3)
c_sr <- "blue4"
c_est <- transparent(c_sr, trans.val = 0.5)
c_srci <- transparent(c_sr, trans.val = 0.8)
c_arr <- "darkgray"

plot(S_obs, M0_obs, pch = 16, col = c_obs, las = 1,
     cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5, xaxs = "i", yaxs = "i",
     xlab = "Spawners", ylab = "", xlim = c(0, max(S)),
     ylim = range(0, M0_obs, apply(M0_IPM, 2, quantile, 0.975, na.rm = T), na.rm = T)*1.02)
mtext("Smolts", side = 2, line = 3.5, cex = 1.5)

points(apply(S_IPM, 2, median), apply(M0_IPM, 2, median), pch = 16, col = c_est)
arrows(S_obs, M0_obs, apply(S_IPM, 2, median), apply(M0_IPM, 2, median, na.rm = T), 
       col = c_arr, length = 0.1)
segments(x0 = apply(S_IPM, 2, quantile, 0.025, na.rm = T), 
         y0 = apply(M0_IPM, 2, median, na.rm = T), 
         x1 = apply(S_IPM, 2, quantile, 0.975), col = c_est)
segments(x0 = apply(S_IPM, 2, median), 
         y0 = apply(M0_IPM, 2, quantile, 0.025, na.rm = T), 
         y1 = apply(M0_IPM, 2, quantile, 0.975, na.rm = T), col = c_est)

lines(S[1,], apply(M0hat_IPM, 2, median, na.rm = T), lwd = 3, col = c_sr)
polygon(c(S[1,], rev(S[1,])), 
        c(apply(M0hat_IPM, 2, quantile, 0.025, na.rm = T), 
          rev(apply(M0hat_IPM, 2, quantile, 0.975, na.rm = T))), 
        col = c_srci, border = NA)
legend("topright", legend = c("observations", "states (95% CI)"), 
       pch = 16, col = c(c_obs, c_est), lty = c(NA,1))

rm(list = c("S","alpha","Rmax","M0_obs","M0_IPM","S_obs","S_IPM","n_Mage_obs","q_M_obs",
            "q_M_IPM","M0hat_IPM","c_obs","c_sr","c_est","c_srci","c_arr","m2","m3"))

# dev.off()


#-------------------------------------------------------------------------
# Posterior distributions (and priors) of S-R parameters
#-------------------------------------------------------------------------

dev.new(width = 10, height = 5)
# png(filename="SR_params.png", width=10, height=5, units="in", res=200, type="cairo-png")
par(mfrow = c(1,2), mar = c(5.1,4.5,1,0))

c1 <- transparent("blue4", trans.val = 0.3)

# Posterior of log(alpha)
hist(log(extract1(fit_BH,"alpha")), 15, prob = TRUE, col = c1, las = 1, 
     cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5,
     xlab = bquote(log(alpha)), ylab = "Probability density", main = "Intrinsic productivity")
curve(dnorm(x,2,2), col = "gray", lwd = 2, add = TRUE)

# Posterior of log(Rmax)
hist(extract1(fit_BH,"Rmax"), 15, prob = TRUE, col = c1, las = 1, 
     cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5,
     xlab = bquote(italic(R)[max]), ylab = "", main = "Maximum smolts")
curve(dlnorm(x,2,3), col = "gray", lwd = 2, add = TRUE)

rm(c1)
# dev.off()


#-------------------------------------------------------------------------
# Time series of observed and estimated smolts and spawners
#-------------------------------------------------------------------------

dev.new(width = 7, height = 7)
# png(filename="M_S_timeseries.png", width=7, height=7, units="in", res=200, type="cairo-png")
par(mfrow = c(2,1), mar = c(4.5, 5.1, 0.5, 0.5))

year <- fishdata$year
M_obs <- fishdata$M_obs/1000
S_obs <- fishdata$S_obs
M_IPM <- extract1(fit_BH,"M")/1000
S_IPM <- extract1(fit_BH,"S")

c_obs <- transparent("orangered3", trans.val = 0.3)
c_sr <- "blue4"
c_srci <- transparent(c_sr, trans.val = 0.8)

# Smolts
plot(year, apply(M_IPM, 2, median), type = "l", lwd = 3, col = c_sr, 
     las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n", 
     ylim = range(0, M_obs, apply(M_IPM, 2, quantile, 0.975), na.rm = TRUE),
     xlab = "", ylab = "Smolts (thousands)")
abline(v = max(fishdata$year[fishdata$obs_type == "past"]), col = "darkgray", lwd = 2)
polygon(c(year, rev(year)), 
        c(apply(M_IPM, 2, quantile, 0.025), rev(apply(M_IPM, 2, quantile, 0.975))),
        col = c_srci, border = NA)
points(year, M_obs, pch = 16, cex = 1.2, col = c_obs)
axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.5)
rug(year[year %% 10 != 0], ticksize = -0.01)
rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.04)

# Spawners
plot(year, apply(S_IPM, 2, median), type = "l", lwd = 3, col = c_sr, 
     las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n",
     ylim = range(0, S_obs, apply(S_IPM, 2, quantile, 0.975), na.rm = TRUE),
     xlab = "Year", ylab = "")
mtext("Spawners", side = 2, line = 3.5, cex = par("cex")*1.5)
abline(v = max(fishdata$year[fishdata$obs_type == "past"]), col = "darkgray", lwd = 2)
polygon(c(year, rev(year)), 
        c(apply(S_IPM, 2, quantile, 0.025), rev(apply(S_IPM, 2, quantile, 0.975))),
        col = c_srci, border = NA)
points(year, S_obs, pch = 16, cex = 1.2, col = c_obs)
axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.5)
rug(year[year %% 10 != 0], ticksize = -0.01)
rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.04)

rm(list = c("year","S_obs","M_obs","S_IPM","M_IPM","c_obs","c_sr","c_srci"))

# dev.off()


#-------------------------------------------------------------------------
# Posterior distributions of 2019 spawner forecast
#-------------------------------------------------------------------------

dev.new(width = 7, height = 7)
# png(filename="S_forecast_2019.png", width=7, height=7, units="in", res=200, type="cairo-png")
par(mar = c(5.1,6,1,0))

c1 <- transparent("blue4", trans.val = 0.3)

hist(extract1(fit_BH,"S"), 15, prob = TRUE, col = c1, las = 1, 
     cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5,
     xlab = "Spawners", ylab = "", main = "2019 Forecast")
mtext("Probability density", side = 2, line = 4.5, cex = par("cex")*1.5)

rm(c1)
# dev.off()


#-------------------------------------------------------------------------
# Time series of observed and estimated proportion jacks
#-------------------------------------------------------------------------

dev.new(width = 7, height = 5)
# png(filename="p_jack_timeseries.png", width=7, height=5, units="in", res=200, type="cairo-png")
par(mar = c(4.5, 4.5, 0.5, 0.5))

year <- fishdata$year
p_jack_obs <- run_recon(fishdata)$p_age2_obs
p_jack_IPM <- extract1(fit_BH,"q")[,,1]

c_obs <- transparent("orangered3", trans.val = 0.3)
c_sr <- "blue4"
c_srci <- transparent(c_sr, trans.val = 0.8)

plot(year, apply(p_jack_IPM, 2, median), type = "l", lwd = 3, col = c_sr, 
     las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n",
     ylim = range(p_jack_obs, apply(p_jack_IPM, 2, quantile, 0.975)),
     xlab = "Year", ylab = "Proportion jacks")
abline(v = max(fishdata$year[fishdata$obs_type == "past"]), col = "darkgray", lwd = 2)
polygon(c(year, rev(year)), 
        c(apply(p_jack_IPM, 2, quantile, 0.025), rev(apply(p_jack_IPM, 2, quantile, 0.975))),
        col = c_srci, border = NA)
points(year[fishdata$obs_type=="past"], p_jack_obs[fishdata$obs_type=="past"], 
       pch = 16, cex = 1.2, col = c_obs)
axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.5)
rug(year[year %% 10 != 0], ticksize = -0.01)
rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.02)

rm(list = c("year","p_jack_obs","p_jack_IPM","c_obs","c_sr","c_srci"))

# dev.off()

