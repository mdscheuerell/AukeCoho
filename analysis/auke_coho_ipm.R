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
  devtools::install_github("ebuhle/salmonIPM", 
                           auth_token = "")
  library("salmonIPM")
}
## for data munging
if(!require("readr")) {
  install.packages("readr")
  library("readr")
}
if(!require("dplyr")) {
  install.packages("dplyr")
  library("dplyr")
  options(tibble.print_max = Inf, tibble.width = Inf)
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
## for dir management
if(!require("here")) {
  install.packages("here")
  library("here")
}

## ----get_data------------------------------------------------------------
## directory
datadir <- here("data")
## filename
datafile <- dir(datadir)[grep("auke_coho_data", dir(datadir))]
## read data
fishdat <- read_csv(file.path(datadir, datafile))

## ----fit_model, eval=TRUE------------------------------------------------
fit_ipm <- salmonIPM(fishdat, model = "IPM", SR_fun = "Ricker", pool_pops = FALSE, 
                     chains = 3, iter = 1500, warmup = 1000,
                     control = list(adapt_delta = 0.999, stepsize = 0.01, max_treedepth = 13))

## ----print_fitted_model---------------------------------------------------
print(fit_ipm, pars = c("B_rate_all","p","q","S_tot","R_tot"), include = FALSE)

## ----shinystan------------------------------------------------------------
launch_shinystan(fit_ipm)

## ----FIGURES--------------------------------------------------------------

#-------------------------------------------------------------------------
# S-R curve overlaid with observations and states
#-------------------------------------------------------------------------

dev.new(width = 7, height = 7)
# png(filename="SR.png", width=7, height=7, units="in", res=200, type="cairo-png")
SR_fun <- "Ricker"
SR <- function(a, Rmax, S, A, SR_fun) 
{
  switch(SR_fun, 
         BH = a*S/(A + a*S/Rmax),
         Ricker = a*(S/A)*exp(-a*S/(A*exp(1)*Rmax)))
}

yy <- stan_data(fishdat, model = "RR")$year
S_tot_obs <- fishdat$S_tot_obs
R_tot_obs <- run_recon(fishdat)$R
S_tot_IPM <- extract1(fit_ipm,"S_tot")
R_tot_IPM <- extract1(fit_ipm,"R_tot")

S <- matrix(seq(0, max(S_tot_obs, apply(S_tot_IPM, 2, quantile, 0.975), na.rm = T)*1.02, length = 500),
            nrow = sum(fit_ipm@sim$n_save - fit_ipm@sim$warmup2), ncol = 500, byrow = T)
a <- as.vector(extract1(fit_ipm,"a"))
Rmax <- as.vector(extract1(fit_ipm,"Rmax"))
R_IPM <- SR(a = a, Rmax = Rmax, S = S, A = 1, SR_fun)

c_obs <- transparent("orangered3", trans.val = 0.4)
c_sr <- "blue4"
c_est <- transparent(c_sr, trans.val = 0.5)
c_srci <- transparent(c_sr, trans.val = 0.8)
c_arr <- "darkgray"

plot(S_tot_obs, R_tot_obs, pch = 16, col = c_obs, las = 1,
     cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5, xaxs = "i", yaxs = "i",
     xlab = "Spawners", ylab = "Recruits", xlim = c(0, max(S)),
     ylim = range(0, R_tot_obs, apply(R_tot_IPM, 2, quantile, 0.975), na.rm = T)*1.02)

points(apply(S_tot_IPM, 2, median), apply(R_tot_IPM, 2, median), pch = 16, col = c_est)
arrows(S_tot_obs, R_tot_obs, apply(S_tot_IPM, 2, median), apply(R_tot_IPM, 2, median), col = c_arr, length = 0.1)
segments(x0 = apply(S_tot_IPM, 2, quantile, 0.025), y0 = apply(R_tot_IPM, 2, median), 
         x1 = apply(S_tot_IPM, 2, quantile, 0.975), col = c_est)
segments(x0 = apply(S_tot_IPM, 2, median), y0 = apply(R_tot_IPM, 2, quantile, 0.025), 
         y1 = apply(R_tot_IPM, 2, quantile, 0.975), col = c_est)

lines(S[1,], apply(R_IPM, 2, median), lwd = 3, col = c_sr)
polygon(c(S[1,], rev(S[1,])), 
        c(apply(R_IPM, 2, quantile, 0.025), rev(apply(R_IPM, 2, quantile, 0.975))), 
        col = c_srci, border = NA)
legend("topright", legend = c("observations", "states (95% CI)"), 
       pch = 16, col = c(c_obs, c_est), lty = c(NA,1))

rm(list = c("yy","S","a","Rmax","R_tot_obs","R_tot_IPM","S_tot_obs","S_tot_IPM",
            "R_IPM","c_obs","c_sr","c_est","c_srci","c_arr"))

# dev.off()


#-------------------------------------------------------------------------
# Posterior distributions of S-R parameters
#-------------------------------------------------------------------------

dev.new(width = 10, height = 5)
# png(filename="SR_params.png", width=7, height=3.5, units="in", res=200, type="cairo-png")
par(mfrow = c(1,2), mar = c(5.1,4.5,1,0))

c1 <- transparent("blue4", trans.val = 0.3)

# Posterior of log(a)
hist(log(extract1(fit_ipm,"a")), prob = TRUE, col = c1, las = 1, 
     cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5,
     xlab = bquote(log(alpha)), ylab = "Probability density", main = "Intrinsic productivity")

# Posterior of log(Rmax)
hist(extract1(fit_ipm,"Rmax"), prob = TRUE, col = c1, las = 1, 
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

year <- fishdat$year
S_tot_obs <- fishdat$S_tot_obs
R_tot_obs <- run_recon(fishdat)$R
RS_obs <- R_tot_obs/S_tot_obs
S_tot_IPM <- extract1(fit_ipm,"S_tot")
R_tot_IPM <- extract1(fit_ipm,"R_tot")
RS_IPM <- R_tot_IPM/S_tot_IPM

c_obs <- transparent("orangered3", trans.val = 0.3)
c_sr <- "blue4"
c_srci <- transparent(c_sr, trans.val = 0.8)

# Spawners
plot(year, apply(S_tot_IPM, 2, median), type = "l", lwd = 3, col = c_sr, 
     las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n",
     ylim = range(0, S_tot_obs, apply(S_tot_IPM, 2, quantile, 0.975)),
     xlab = "", ylab = "")
mtext("Spawners", side = 2, line = 3.5, cex = par("cex")*1.5)
polygon(c(year, rev(year)), 
        c(apply(S_tot_IPM, 2, quantile, 0.025), rev(apply(S_tot_IPM, 2, quantile, 0.975))),
        col = c_srci, border = NA)
points(year, S_tot_obs, pch = 16, cex = 1.2, col = c_obs)
axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.5)
rug(year[year %% 10 != 0], ticksize = -0.01)
rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.04)

# Recruits per spawner
plot(year, apply(RS_IPM, 2, median), type = "l", lwd = 3, col = c_sr, 
     las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n", 
     ylim = range(0, RS_obs, apply(RS_IPM, 2, quantile, 0.975), na.rm = TRUE),
     xlab = "Year", ylab = "Recruits per spawner")
abline(h = 1, lty = 2, lwd = 2)
polygon(c(year, rev(year)), 
        c(apply(RS_IPM, 2, quantile, 0.025), rev(apply(RS_IPM, 2, quantile, 0.975))),
        col = c_srci, border = NA)
points(year, RS_obs, pch = 16, cex = 1.2, col = c_obs)
axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.5)
rug(year[year %% 10 != 0], ticksize = -0.01)
rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.04)

rm(list = c("year","S_tot_obs","R_tot_obs","RS_obs","S_tot_IPM","R_tot_IPM","RS_IPM",
            "c_obs","c_sr","c_srci"))

# dev.off()


#-------------------------------------------------------------------------
# Time series of observed and estimated proportion jacks
#-------------------------------------------------------------------------

dev.new(width = 7, height = 5)
# png(filename="p_jack_timeseries.png", width=7, height=5, units="in", res=200, type="cairo-png")
par(mar = c(4.5, 4.5, 0.5, 0.5))

year <- fishdat$year
p_jack_obs <- run_recon(fishdat)$p_age2_obs
p_jack_IPM <- extract1(fit_ipm,"q")[,,1]

c_obs <- transparent("orangered3", trans.val = 0.3)
c_sr <- "blue4"
c_srci <- transparent(c_sr, trans.val = 0.8)

plot(year, apply(p_jack_IPM, 2, median), type = "l", lwd = 3, col = c_sr, 
     las = 1, cex.lab = 1.5, cex.axis = 1.2, xaxt = "n",
     ylim = range(p_jack_obs, apply(p_jack_IPM, 2, quantile, 0.975)),
     xlab = "Year", ylab = "Proportion jacks")
polygon(c(year, rev(year)), 
        c(apply(p_jack_IPM, 2, quantile, 0.025), rev(apply(p_jack_IPM, 2, quantile, 0.975))),
        col = c_srci, border = NA)
points(year, p_jack_obs, pch = 16, cex = 1.2, col = c_obs)
axis(side = 1, at = year[year %% 10 == 0], cex.axis = 1.5)
rug(year[year %% 10 != 0], ticksize = -0.01)
rug(year[year %% 10 != 0 & year %% 5 == 0], ticksize = -0.02)

rm(list = c("year","p_jack_obs","p_jack_IPM","c_obs","c_sr","c_srci"))

# dev.off()

