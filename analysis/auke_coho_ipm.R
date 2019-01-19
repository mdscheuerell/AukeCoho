## ----set_options, echo=FALSE, cache=TRUE---------------------------------
options(width = 120)
if(Sys.info()["sysname"] == "Windows") options(device = windows)

## ----load_pkgs, message=FALSE, warning=FALSE-----------------------------
## for inference
if(!require("rstan")) {
  install.packages("rstan")
  library("rstan")
}
if(!require("loo")) {
  install.packages("loo")
  library("loo")
}
if(!require("salmonIPM")) {
  if(file.exists("github_auth_token.txt")) 
    github_auth_token <- scan("github_auth_token.txt", what = "character")
  devtools::install_github("ebuhle/salmonIPM", ref = "spawner-smolt-models", 
                           auth_token = github_auth_token)
  library("salmonIPM")
}
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
# fishdata <- fishdata[fishdata$year <= 2019,]  # truncate future years
## read covariate data
pdofile <- dir(datadir)[grep("pdo", dir(datadir))]
pdodat <- read_csv(file.path(datadir, pdofile))
covfile <- dir(datadir)[grep("covariates", dir(datadir))]
cov_raw <- read_csv(file.path(datadir, covfile))
## Subset PDO to Nov-Feb and index to outmigration year
pdo <- pdodat[pdodat$month %in% c(11,12,1,2),]
pdo$winter_year <- ifelse(pdo$month %in% 11:12, pdo$year, pdo$year - 1)
pdo <- aggregate(pdo ~ winter_year, data = pdo, mean)

## Add PDO to covariate data
cov_raw$winter_pdo <- pdo$pdo[match(cov_raw$year, pdo$winter_year)]
covdata <- data.frame(brood_year = fishdata$year[fishdata$obs_type == "past"])

covdata$parr0_year <- covdata$brood_year + 1
covdata$parr1_year <- covdata$brood_year + 2
covdata$smolt1_year <- covdata$brood_year + 2
covdata$smolt2_year <- covdata$brood_year + 3

covdata$parr0_stream_temp <- cov_raw$mean_stream_temp[match(covdata$parr0_year, cov_raw$year)]
covdata$parr1_stream_temp <- cov_raw$mean_stream_temp[match(covdata$parr1_year, cov_raw$year)]
covdata$mean_stream_temp <- rowMeans(covdata[,c("parr0_stream_temp","parr1_stream_temp")])

covdata$smolt1_hpc_release <- cov_raw$hpc_release[match(covdata$smolt1_year, cov_raw$year)]
covdata$smolt2_hpc_release <- cov_raw$hpc_release[match(covdata$smolt2_year, cov_raw$year)]
covdata$mean_hpc_release <- rowMeans(covdata[,c("smolt1_hpc_release","smolt2_hpc_release")])/1e6

covdata$smolt1_pdo <- pdo$pdo[match(covdata$smolt1_year, cov_raw$year)]
covdata$smolt2_pdo <- pdo$pdo[match(covdata$smolt2_year, cov_raw$year)]
covdata$mean_pdo <- rowMeans(covdata[,c("smolt1_pdo","smolt2_pdo")])

env_data <- covdata[,c("brood_year","mean_stream_temp","mean_hpc_release","mean_pdo")]
env_data[,-1] <- scale(env_data[,-1])
env_data <- na.omit(env_data)

## Truncate fish data to non-missing covariate data
fishdata_env <- fishdata[fishdata$year %in% na.omit(env_data)$brood_year,]

## ----fit_model, eval=TRUE------------------------------------------------
fit_ipm <- salmonIPM(fishdata, model = "IPM", SR_fun = "Ricker", pool_pops = FALSE, 
                     chains = 3, iter = 1500, warmup = 1000,
                     control = list(adapt_delta = 0.99))

## ----print_fitted_model---------------------------------------------------
print(fit_ipm, pars = c("B_rate_all","p","q","S","R"), include = FALSE)

## ----shinystan------------------------------------------------------------
launch_shinystan(fit_ipm)

## ----fit_model_with_covariates, eval=TRUE------------------------------------------------
fit_ipm_env <- salmonIPM(fishdata_env, env_data = env_data[,-1,drop = FALSE], 
                         model = "IPM", SR_fun = "Ricker", pool_pops = FALSE, 
                         chains = 3, iter = 1500, warmup = 1000,
                         control = list(adapt_delta = 0.99))

## ----print_fitted_model---------------------------------------------------
print(fit_ipm_env, pars = c("B_rate_all","p","q","S","R"), include = FALSE)

## ----shinystan------------------------------------------------------------
launch_shinystan(fit_ipm_env)

## ----FIGURES--------------------------------------------------------------

#-------------------------------------------------------------------------
# S-R curve overlaid with observations and states
#-------------------------------------------------------------------------

dev.new(width = 7, height = 7)
# png(filename="SR.png", width=7, height=7, units="in", res=200, type="cairo-png")
SR_fun <- "Ricker"
SR <- function(alpha, Rmax, S, A, SR_fun) 
{
  switch(SR_fun, 
         BH = alpha*S/(A + alpha*S/Rmax),
         Ricker = alpha*(S/A)*exp(-alpha*S/(A*exp(1)*Rmax)))
}

yy <- stan_data(fishdata, model = "RR")$year
S_obs <- fishdata$S_obs
R_obs <- run_recon(fishdata)$R
S_IPM <- extract1(fit_ipm,"S")
R_IPM <- extract1(fit_ipm,"R")

S <- matrix(seq(0, max(S_obs, apply(S_IPM, 2, quantile, 0.975), na.rm = T)*1.02, length = 500),
            nrow = sum(fit_ipm@sim$n_save - fit_ipm@sim$warmup2), ncol = 500, byrow = T)
alpha <- as.vector(extract1(fit_ipm,"alpha"))
Rmax <- as.vector(extract1(fit_ipm,"Rmax"))
Rhat_IPM <- SR(alpha = alpha, Rmax = Rmax, S = S, A = 1, SR_fun)

c_obs <- transparent("orangered3", trans.val = 0.3)
c_sr <- "blue4"
c_est <- transparent(c_sr, trans.val = 0.5)
c_srci <- transparent(c_sr, trans.val = 0.8)
c_arr <- "darkgray"

plot(S_obs, R_obs, pch = 16, col = c_obs, las = 1,
     cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5, xaxs = "i", yaxs = "i",
     xlab = "Spawners", ylab = "Recruits", xlim = c(0, max(S)),
     ylim = range(0, R_obs, apply(R_IPM, 2, quantile, 0.975), na.rm = T)*1.02)

points(apply(S_IPM, 2, median), apply(R_IPM, 2, median), pch = 16, col = c_est)
arrows(S_obs, R_obs, apply(S_IPM, 2, median), apply(R_IPM, 2, median), col = c_arr, length = 0.1)
segments(x0 = apply(S_IPM, 2, quantile, 0.025), y0 = apply(R_IPM, 2, median), 
         x1 = apply(S_IPM, 2, quantile, 0.975), col = c_est)
segments(x0 = apply(S_IPM, 2, median), y0 = apply(R_IPM, 2, quantile, 0.025), 
         y1 = apply(R_IPM, 2, quantile, 0.975), col = c_est)

lines(S[1,], apply(Rhat_IPM, 2, median), lwd = 3, col = c_sr)
polygon(c(S[1,], rev(S[1,])), 
        c(apply(Rhat_IPM, 2, quantile, 0.025), rev(apply(Rhat_IPM, 2, quantile, 0.975))), 
        col = c_srci, border = NA)
legend("topright", legend = c("observations", "states (95% CI)"), 
       pch = 16, col = c(c_obs, c_est), lty = c(NA,1))

rm(list = c("yy","S","alpha","Rmax","R_obs","R_IPM","S_obs","S_IPM",
            "Rhat_IPM","c_obs","c_sr","c_est","c_srci","c_arr"))

# dev.off()


#-------------------------------------------------------------------------
# Posterior distributions of S-R parameters
#-------------------------------------------------------------------------

dev.new(width = 10, height = 5)
# png(filename="SR_params.png", width=10, height=5, units="in", res=200, type="cairo-png")
par(mfrow = c(1,2), mar = c(5.1,4.5,1,0))

c1 <- transparent("blue4", trans.val = 0.3)

# Posterior of log(alpha)
hist(log(extract1(fit_ipm,"alpha")), 15, prob = TRUE, col = c1, las = 1, 
     cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5,
     xlab = bquote(log(alpha)), ylab = "Probability density", main = "Intrinsic productivity")

# Posterior of log(Rmax)
hist(extract1(fit_ipm,"Rmax"), 15, prob = TRUE, col = c1, las = 1, 
     cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5,
     xlab = bquote(italic(R)[max]), ylab = "", main = "Maximum recruits")

rm(c1)
# dev.off()


#-------------------------------------------------------------------------
# Time series of observed and estimated spawners and R/S
#-------------------------------------------------------------------------

dev.new(width = 7, height = 7)
# png(filename="S_RS_timeseries.png", width=7, height=7, units="in", res=200, type="cairo-png")
par(mfrow = c(2,1), mar = c(4.5, 5.1, 0.5, 0.5))

year <- fishdata$year
S_obs <- fishdata$S_obs
R_obs <- run_recon(fishdata)$R
RS_obs <- R_obs/S_obs
S_IPM <- extract1(fit_ipm,"S")
R_IPM <- extract1(fit_ipm,"R")
RS_IPM <- R_IPM/S_IPM

c_obs <- transparent("orangered3", trans.val = 0.3)
c_sr <- "blue4"
c_srci <- transparent(c_sr, trans.val = 0.8)

# Spawners
plot(year, apply(S_IPM, 2, median), type = "l", lwd = 3, col = c_sr, 
     las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n",
     ylim = range(0, S_obs, apply(S_IPM, 2, quantile, 0.975), na.rm = TRUE),
     xlab = "", ylab = "")
mtext("Spawners", side = 2, line = 3.5, cex = par("cex")*1.5)
abline(v = max(fishdata$year[fishdata$obs_type == "past"]), col = "darkgray", lwd = 2)
polygon(c(year, rev(year)), 
        c(apply(S_IPM, 2, quantile, 0.025), rev(apply(S_IPM, 2, quantile, 0.975))),
        col = c_srci, border = NA)
points(year, S_obs, pch = 16, cex = 1.2, col = c_obs)
axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.5)
rug(year[year %% 10 != 0], ticksize = -0.01)
rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.04)

# Recruits per spawner
plot(year, apply(RS_IPM, 2, median), type = "l", lwd = 3, col = c_sr, 
     las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n", 
     ylim = range(0, RS_obs, apply(RS_IPM, 2, quantile, 0.975), na.rm = TRUE),
     xlab = "Year", ylab = "Recruits per spawner")
abline(h = 1, lty = 2, lwd = 2)
abline(v = max(fishdata$year[fishdata$obs_type == "past"]), col = "darkgray", lwd = 2)
polygon(c(year, rev(year)), 
        c(apply(RS_IPM, 2, quantile, 0.025), rev(apply(RS_IPM, 2, quantile, 0.975))),
        col = c_srci, border = NA)
points(year, RS_obs, pch = 16, cex = 1.2, col = c_obs)
axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.5)
rug(year[year %% 10 != 0], ticksize = -0.01)
rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.04)

rm(list = c("year","S_obs","R_obs","RS_obs","S_IPM","R_IPM","RS_IPM",
            "c_obs","c_sr","c_srci"))

# dev.off()


#-------------------------------------------------------------------------
# Posterior distributions of 2019 spawner forecast
#-------------------------------------------------------------------------

dev.new(width = 7, height = 7)
# png(filename="S_forecast_2019.png", width=7, height=7, units="in", res=200, type="cairo-png")
par(mar = c(5.1,6,1,0))

c1 <- transparent("blue4", trans.val = 0.3)

hist(extract1(fit_ipm,"S"), 15, prob = TRUE, col = c1, las = 1, 
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
p_jack_IPM <- extract1(fit_ipm,"q")[,,1]

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

