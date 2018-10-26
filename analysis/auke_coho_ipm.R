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
datafile <- dir(datadir)[grep("ipm", dir(datadir))]
## read data
fishdat <- read_csv(file.path(datadir, datafile))

## ----trim_years----------------------------------------------------------
## trim early years
fishdat <- filter(fishdat, year >= 1994)

## ----set_area------------------------------------------------------------
# set spawning area
fishdat$A <- rep(1, dim(fishdat)[1])

## ----set_harvest---------------------------------------------------------
## remove years with no F estimates
fishdat <- filter(fishdat, !is.na(F_rate))

## ----fit_model, eval=TRUE------------------------------------------------
fit_ipm <- salmonIPM(fishdat, model = "IPM",
                     chains = 3, iter = 1000, warmup = 500,
                     control = list(adapt_delta = 0.95, stepsize = 0.01, max_treedepth = 13))

