## ----set_options, echo=FALSE, cache=TRUE---------------------------------
options(width = 120)

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
  devtools::install_github("ebuhle/salmonIPM")
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
}
## for output/plotting
if(!require("viridis")) {
  install.packages("viridis")
  library("viridis")
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
fit_ipm <- salmonIPM(fishdat, model = "IPM", pool_pops = FALSE, 
                     chains = 3, iter = 1500, warmup = 1000,
                     control = list(adapt_delta = 0.999, stepsize = 0.01, max_treedepth = 13))

## ----print_fitted_model---------------------------------------------------
print(fit_ipm, pars = c("B_rate_all","p","q","S_tot","R_tot"), include = FALSE)

## ----shinystan------------------------------------------------------------
launch_shinystan(fit_ipm)


