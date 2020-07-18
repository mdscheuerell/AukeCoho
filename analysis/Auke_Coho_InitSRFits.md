Auke\_Coho\_InitSRFits
================
DT
7/17/2020

-----

[**Mark D. Scheuerell**](https://faculty.washington.edu/scheuerl/)  
*Fish Ecology Division, Northwest Fisheries Science Center, National
Marine Fisheries Service, National Oceanic and Atmospheric
Administration, Seattle, WA, <mark.scheuerell@noaa.gov>*

**Eric R. Buhle**  
*Quantitative Science Inc., Northwest Fisheries Science Center, National
Marine Fisheries Service, National Oceanic and Atmospheric
Administration, Seattle, WA USA, <eric.buhle@noaa.gov>*

**David Tallmon, Scott Vulstek**

-----

This is version 0.20.07.17.

Note that author list and order subject to change.

-----

# Requirements

All analyses require the [R software](https://cran.r-project.org/)
(v3.4+) for data retrieval, data processing, and summarizing model
results, and the [Stan software](http://mc-stan.org/) (v2.17.0) for
Hamiltonian Monte Carlo (HMC) sampling.

We also need a few packages that are not included with the base
installation of R, so we begin by installing them (if necessary) and
then loading them.

``` r
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
```

# Load fish data

We begin by reading in the fish data file

``` r
setwd(here("analysis"))
datadir <- here("data")
datafile <- dir(datadir)[grep("auke_coho_data_1980-2019", dir(datadir))]
fishdata <- read_csv(file.path(datadir, datafile))
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   obs_type = col_character()
    ## )

    ## See spec(...) for full column specifications.

``` r
fishdata <- fishdata[fishdata$year <= 2020,]  # truncate future years
```

Read in covariates that will be used in future analyses. Need to discuss
indexing covariates before using.

``` r
covfile <- dir(datadir)[grep("covariates_1980-2019", dir(datadir))]
cov_raw <- read_csv(file.path(datadir, covfile))
```

    ## Parsed with column specification:
    ## cols(
    ##   year = col_double(),
    ##   hpc_release = col_double(),
    ##   pdo_nov_jan = col_double(),
    ##   gauge_spring = col_double()
    ## )

``` r
env_data <- cov_raw
# HPC release in 100s of millions of fish
env_data$hpc_release <- (env_data$hpc_release - mean(env_data$hpc_release))/1e8
# winter PDO
env_data$pdo_nov_jan <- scale(env_data$pdo_nov_jan)
# spring gauge
env_data$gauge_spring <- scale(env_data$gauge_spring)
```

Fit the first SR function: Ricker. Print the output.

``` r
fit_Ricker <- salmonIPM(fishdata, stan_model = "IPM_SMaS_np", SR_fun = "Ricker",
pars = c(stan_pars("IPM_SMaS_np"), "LL"),chains = 3, cores = 3, iter = 1000, warmup = 500,
control = list(adapt_delta = 0.99, max_treedepth = 13))
```

    ## Warning: The largest R-hat is NA, indicating chains have not mixed.
    ## Running the chains for more iterations may help. See
    ## http://mc-stan.org/misc/warnings.html#r-hat

    ## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
    ## Running the chains for more iterations may help. See
    ## http://mc-stan.org/misc/warnings.html#bulk-ess

    ## Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
    ## Running the chains for more iterations may help. See
    ## http://mc-stan.org/misc/warnings.html#tail-ess

``` r
print(fit_Ricker, pars = c("p_M","q_M","s_MS","p_MS","q_MS","q_GR","M","S","R","LL"), include = FALSE, probs = c(0.025,0.5,0.975))
```

    ## Inference for Stan model: IPM_SMaS_np.
    ## 3 chains, each with iter=1000; warmup=500; thin=1; 
    ## post-warmup draws per chain=500, total post-warmup draws=1500.
    ## 
    ##                        mean se_mean     sd      2.5%       50%     97.5% n_eff Rhat
    ## alpha[1]              16.95    0.15   2.92     11.97     16.75     23.12   372 1.00
    ## Rmax[1]             5935.14   31.11 582.00   4889.75   5902.86   7107.85   350 1.01
    ## rho_M[1]               0.40    0.01   0.17      0.05      0.40      0.72   305 1.00
    ## sigma_M[1]             0.29    0.00   0.04      0.23      0.29      0.38   541 1.00
    ## tau_M[1]               0.05    0.00   0.01      0.03      0.05      0.07   264 1.01
    ## mu_p_M[1,1]            0.55    0.00   0.02      0.50      0.55      0.60   174 1.01
    ## mu_p_M[1,2]            0.45    0.00   0.02      0.40      0.45      0.50   174 1.01
    ## sigma_p_M[1,1]         0.53    0.00   0.07      0.42      0.53      0.69   354 1.00
    ## R_p_M[1,1,1]           1.00     NaN   0.00      1.00      1.00      1.00   NaN  NaN
    ## mu_MS[1,1]             0.15    0.00   0.04      0.09      0.15      0.22   304 1.00
    ## mu_MS[1,2]             0.32    0.00   0.07      0.19      0.32      0.46   427 1.00
    ## rho_MS[1,1]            0.58    0.01   0.14      0.29      0.59      0.81   405 1.00
    ## rho_MS[1,2]            0.63    0.01   0.13      0.34      0.65      0.84   564 1.00
    ## sigma_MS[1,1]          0.63    0.00   0.09      0.48      0.62      0.82   412 1.01
    ## sigma_MS[1,2]          0.67    0.00   0.11      0.49      0.65      0.90   529 1.01
    ## R_MS[1,1,1]            1.00     NaN   0.00      1.00      1.00      1.00   NaN  NaN
    ## R_MS[1,1,2]            0.02    0.01   0.18     -0.32      0.02      0.36   427 1.01
    ## R_MS[1,2,1]            0.02    0.01   0.18     -0.32      0.02      0.36   427 1.01
    ## R_MS[1,2,2]            1.00    0.00   0.00      1.00      1.00      1.00  1279 1.00
    ## mu_p_MS[1,1,1]         0.03    0.00   0.01      0.02      0.03      0.04   601 1.00
    ## mu_p_MS[1,1,2]         0.97    0.00   0.01      0.96      0.97      0.98   601 1.00
    ## mu_p_MS[1,2,1]         0.23    0.00   0.02      0.19      0.23      0.27   375 1.01
    ## mu_p_MS[1,2,2]         0.77    0.00   0.02      0.73      0.77      0.81   375 1.01
    ## sigma_p_MS[1,1,1]      0.93    0.01   0.19      0.62      0.91      1.36   671 1.00
    ## sigma_p_MS[1,2,1]      0.57    0.00   0.08      0.43      0.56      0.77   489 1.01
    ## R_p_MS[1,1,1]          1.00     NaN   0.00      1.00      1.00      1.00   NaN  NaN
    ## R_p_MS[1,1,2]          0.23    0.01   0.21     -0.18      0.24      0.63   231 1.02
    ## R_p_MS[1,2,1]          0.23    0.01   0.21     -0.18      0.24      0.63   231 1.02
    ## R_p_MS[1,2,2]          1.00    0.00   0.00      1.00      1.00      1.00  1324 1.00
    ## tau_S[1]               0.05    0.00   0.01      0.03      0.05      0.07   611 1.01
    ## lp__              -41124.46    1.35  21.03 -41165.49 -41123.61 -41083.34   241 1.01
    ## 
    ## Samples were drawn using NUTS(diag_e) at Fri Jul 17 15:40:38 2020.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

Fit a Beverton-Holt function. Print the output.

``` r
fit_BH <- salmonIPM(fishdata, stan_model = "IPM_SMaS_np", SR_fun = "BH",
pars = c(stan_pars("IPM_SMaS_np"), "LL"),chains = 3, cores = 3, iter = 1000, warmup = 500,
control = list(adapt_delta = 0.99, max_treedepth = 13))
```

    ## Warning: There were 142 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See
    ## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

    ## Warning: There were 1338 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 13. See
    ## http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded

    ## Warning: Examine the pairs() plot to diagnose sampling problems

    ## Warning: The largest R-hat is NA, indicating chains have not mixed.
    ## Running the chains for more iterations may help. See
    ## http://mc-stan.org/misc/warnings.html#r-hat

    ## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
    ## Running the chains for more iterations may help. See
    ## http://mc-stan.org/misc/warnings.html#bulk-ess

    ## Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
    ## Running the chains for more iterations may help. See
    ## http://mc-stan.org/misc/warnings.html#tail-ess

``` r
print(fit_BH, pars = c("p_M","q_M","s_MS","p_MS","q_MS","q_GR","M","S","R","LL"), include = FALSE, probs = c(0.025,0.5,0.975))
```

    ## Inference for Stan model: IPM_SMaS_np.
    ## 3 chains, each with iter=1000; warmup=500; thin=1; 
    ## post-warmup draws per chain=500, total post-warmup draws=1500.
    ## 
    ##                             mean se_mean            sd       2.5%           50%          97.5% n_eff Rhat
    ## alpha[1]           3.748432e+152     NaN 4.264749e+153 1339979.30  1.789952e+85  2.902373e+151   NaN 1.00
    ## Rmax[1]             5.473780e+03   55.83  5.161300e+02    4364.07  5.477730e+03   6.450760e+03    85 1.03
    ## rho_M[1]            4.700000e-01    0.01  1.700000e-01       0.15  4.700000e-01   7.800000e-01   132 1.01
    ## sigma_M[1]          2.800000e-01    0.00  3.000000e-02       0.22  2.700000e-01   3.500000e-01   225 1.01
    ## tau_M[1]            5.000000e-02    0.00  1.000000e-02       0.03  5.000000e-02   8.000000e-02   253 1.01
    ## mu_p_M[1,1]         5.500000e-01    0.00  2.000000e-02       0.50  5.500000e-01   5.900000e-01    50 1.03
    ## mu_p_M[1,2]         4.500000e-01    0.00  2.000000e-02       0.41  4.500000e-01   5.000000e-01    50 1.03
    ## sigma_p_M[1,1]      5.400000e-01    0.01  8.000000e-02       0.42  5.400000e-01   7.000000e-01   123 1.01
    ## R_p_M[1,1,1]        1.000000e+00     NaN  0.000000e+00       1.00  1.000000e+00   1.000000e+00   NaN  NaN
    ## mu_MS[1,1]          1.700000e-01    0.00  5.000000e-02       0.10  1.600000e-01   2.900000e-01   102 1.03
    ## mu_MS[1,2]          3.200000e-01    0.01  7.000000e-02       0.19  3.100000e-01   4.700000e-01    93 1.01
    ## rho_MS[1,1]         5.900000e-01    0.01  1.300000e-01       0.31  6.000000e-01   8.200000e-01   107 1.04
    ## rho_MS[1,2]         6.100000e-01    0.01  1.300000e-01       0.33  6.100000e-01   8.300000e-01   111 1.02
    ## sigma_MS[1,1]       6.300000e-01    0.01  9.000000e-02       0.48  6.300000e-01   8.400000e-01   121 1.01
    ## sigma_MS[1,2]       6.600000e-01    0.01  1.100000e-01       0.49  6.500000e-01   9.100000e-01   178 1.00
    ## R_MS[1,1,1]         1.000000e+00     NaN  0.000000e+00       1.00  1.000000e+00   1.000000e+00   NaN  NaN
    ## R_MS[1,1,2]         0.000000e+00    0.02  1.700000e-01      -0.35  1.000000e-02   3.200000e-01    64 1.04
    ## R_MS[1,2,1]         0.000000e+00    0.02  1.700000e-01      -0.35  1.000000e-02   3.200000e-01    64 1.04
    ## R_MS[1,2,2]         1.000000e+00    0.00  0.000000e+00       1.00  1.000000e+00   1.000000e+00  1657 1.00
    ## mu_p_MS[1,1,1]      3.000000e-02    0.00  1.000000e-02       0.02  3.000000e-02   4.000000e-02   189 1.01
    ## mu_p_MS[1,1,2]      9.700000e-01    0.00  1.000000e-02       0.96  9.700000e-01   9.800000e-01   189 1.01
    ## mu_p_MS[1,2,1]      2.300000e-01    0.00  2.000000e-02       0.19  2.300000e-01   2.600000e-01   134 1.02
    ## mu_p_MS[1,2,2]      7.700000e-01    0.00  2.000000e-02       0.74  7.700000e-01   8.100000e-01   134 1.02
    ## sigma_p_MS[1,1,1]   9.100000e-01    0.01  1.800000e-01       0.62  8.900000e-01   1.300000e+00   252 1.02
    ## sigma_p_MS[1,2,1]   5.700000e-01    0.01  9.000000e-02       0.42  5.700000e-01   7.500000e-01   167 1.02
    ## R_p_MS[1,1,1]       1.000000e+00     NaN  0.000000e+00       1.00  1.000000e+00   1.000000e+00   NaN  NaN
    ## R_p_MS[1,1,2]       2.400000e-01    0.02  2.100000e-01      -0.20  2.500000e-01   6.400000e-01    74 1.02
    ## R_p_MS[1,2,1]       2.400000e-01    0.02  2.100000e-01      -0.20  2.500000e-01   6.400000e-01    74 1.02
    ## R_p_MS[1,2,2]       1.000000e+00    0.00  0.000000e+00       1.00  1.000000e+00   1.000000e+00  1537 1.00
    ## tau_S[1]            5.000000e-02    0.00  1.000000e-02       0.03  5.000000e-02   7.000000e-02   360 1.00
    ## lp__               -4.112445e+04    1.60  2.014000e+01  -41167.02 -4.112358e+04  -4.108621e+04   158 1.02
    ## 
    ## Samples were drawn using NUTS(diag_e) at Fri Jul 17 17:06:58 2020.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

Fit an exponential function. Print the output.

``` r
fit_exp <- salmonIPM(fishdata, stan_model = "IPM_SMaS_np", SR_fun = "exp",
pars = c(stan_pars("IPM_SMaS_np"), "LL"),chains = 3, cores = 3, iter = 1000, warmup = 500,
control = list(adapt_delta = 0.99, max_treedepth = 13))
```

    ## Warning: There were 6 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 13. See
    ## http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded

    ## Warning: Examine the pairs() plot to diagnose sampling problems

    ## Warning: The largest R-hat is NA, indicating chains have not mixed.
    ## Running the chains for more iterations may help. See
    ## http://mc-stan.org/misc/warnings.html#r-hat

    ## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
    ## Running the chains for more iterations may help. See
    ## http://mc-stan.org/misc/warnings.html#bulk-ess

    ## Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
    ## Running the chains for more iterations may help. See
    ## http://mc-stan.org/misc/warnings.html#tail-ess

``` r
print(fit_exp, pars = c("p_M","q_M","s_MS","p_MS","q_MS","q_GR","M","S","R","LL"), include = FALSE, probs = c(0.025,0.5,0.975))
```

    ## Inference for Stan model: IPM_SMaS_np.
    ## 3 chains, each with iter=1000; warmup=500; thin=1; 
    ## post-warmup draws per chain=500, total post-warmup draws=1500.
    ## 
    ##                        mean se_mean      sd      2.5%       50%     97.5% n_eff Rhat
    ## alpha[1]               6.54    0.05    0.84      5.05      6.46      8.44   265 1.01
    ## Rmax[1]             4899.18  181.55 5974.02    493.16   2966.67  20724.42  1083 1.00
    ## rho_M[1]               0.37    0.01    0.17      0.04      0.37      0.70   274 1.00
    ## sigma_M[1]             0.44    0.00    0.05      0.35      0.43      0.56   516 1.00
    ## tau_M[1]               0.05    0.00    0.01      0.03      0.05      0.08   502 1.01
    ## mu_p_M[1,1]            0.55    0.00    0.02      0.50      0.55      0.60   206 1.02
    ## mu_p_M[1,2]            0.45    0.00    0.02      0.40      0.45      0.50   206 1.02
    ## sigma_p_M[1,1]         0.53    0.00    0.07      0.42      0.52      0.70   419 1.00
    ## R_p_M[1,1,1]           1.00     NaN    0.00      1.00      1.00      1.00   NaN  NaN
    ## mu_MS[1,1]             0.15    0.00    0.04      0.09      0.15      0.24   405 1.01
    ## mu_MS[1,2]             0.32    0.00    0.07      0.19      0.31      0.46   423 1.01
    ## rho_MS[1,1]            0.58    0.01    0.14      0.29      0.59      0.82   445 1.00
    ## rho_MS[1,2]            0.62    0.01    0.14      0.30      0.63      0.83   596 1.00
    ## sigma_MS[1,1]          0.63    0.00    0.09      0.49      0.62      0.84   518 1.01
    ## sigma_MS[1,2]          0.67    0.00    0.10      0.49      0.66      0.89   717 1.00
    ## R_MS[1,1,1]            1.00     NaN    0.00      1.00      1.00      1.00   NaN  NaN
    ## R_MS[1,1,2]            0.00    0.01    0.18     -0.36      0.00      0.33   489 1.00
    ## R_MS[1,2,1]            0.00    0.01    0.18     -0.36      0.00      0.33   489 1.00
    ## R_MS[1,2,2]            1.00    0.00    0.00      1.00      1.00      1.00  1459 1.00
    ## mu_p_MS[1,1,1]         0.03    0.00    0.01      0.02      0.03      0.04   616 1.01
    ## mu_p_MS[1,1,2]         0.97    0.00    0.01      0.96      0.97      0.98   616 1.01
    ## mu_p_MS[1,2,1]         0.23    0.00    0.02      0.19      0.23      0.26   543 1.01
    ## mu_p_MS[1,2,2]         0.77    0.00    0.02      0.74      0.77      0.81   543 1.01
    ## sigma_p_MS[1,1,1]      0.92    0.01    0.19      0.62      0.90      1.34   675 1.00
    ## sigma_p_MS[1,2,1]      0.56    0.00    0.09      0.42      0.55      0.76   460 1.00
    ## R_p_MS[1,1,1]          1.00     NaN    0.00      1.00      1.00      1.00   NaN  NaN
    ## R_p_MS[1,1,2]          0.26    0.01    0.20     -0.17      0.27      0.60   302 1.00
    ## R_p_MS[1,2,1]          0.26    0.01    0.20     -0.17      0.27      0.60   302 1.00
    ## R_p_MS[1,2,2]          1.00    0.00    0.00      1.00      1.00      1.00  1372 1.00
    ## tau_S[1]               0.05    0.00    0.01      0.04      0.05      0.08   593 1.00
    ## lp__              -41126.21    1.31   20.82 -41169.45 -41125.43 -41089.09   254 1.01
    ## 
    ## Samples were drawn using NUTS(diag_e) at Fri Jul 17 18:02:55 2020.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).
