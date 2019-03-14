detach(package:salmonIPM, unload = TRUE)
library(salmonIPM)
set.seed(666)
dat <- fishdata[1:36,]
fit_BH <- salmonIPM(dat, stan_model = "IPM_SMaS_np", SR_fun = "BH", 
                    env_data = list(M = matrix(0,nrow(dat),1),
                                    MS = matrix(0,nrow(dat),1)),
                    chains = 3, cores = 3, iter = 100, warmup = 50,
                    control = list(adapt_delta = 0.99))

showConnections(); closeAllConnections()

print(fit_BH, pars = c("p_M","q_M","s_MS","p_MS","q_MS","q_GR","M","S","R"), include = FALSE)
fit_BH@par_dims[c("beta_M","beta_MS")]

