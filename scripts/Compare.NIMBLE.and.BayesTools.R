rm(list = ls())

library(nimble)
library(igraph)
library(coda)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(BayesianTools)
library(rrtm)
library(Rprospect)
library(pracma)

#################################################################################################################
# MCMC parameters
Nchains = 2     # Number of MCMC chains
Niter = 5000
Nburnin = 1000
thin = 5

# Wavelength parameters
WLa <- 520                      # Min. Wavelength
WLb <- 2500                     # Max. WL
Delta_WL <- 1                   # WL step
WLs <- seq(WLa,WLb,Delta_WL)    # all wavelengths
pos <- which(400:2500 %in% WLs) # used wavelengths
Nwl <- length(pos)              # Number of used wavelengths

# Data
data <- readRDS("./data/example.spectrum.RDS")
Nleaves <- length(unique(data$leaf.id))

#################################################################################################################
# Nimble

# Approximation of the exponential integral NIMBLE function
e1_approx <- nimbleFunction(
  run = function(x = double(1)) {
    returnType(double(1))

    A <- log((0.56146 / x + 0.65) * (1 + x))
    B <- x^4 * exp(7.7 * x) * (2 + x)^3.7

    return((A^-7.7 + B)^-0.13)
  })

# PROSPECT5 NIMBLE function
NIMprospect5 <- nimbleFunction(
  run = function(N = double(0),Cab = double(0),Car = double(0), Cw = double(0), Cm = double(0),
                 dataspec_p5 = double(2),talf = double(1),t12 = double(1),t21 = double(1), Nwl = double(0)) {

    cc <- matrix(NA,nrow = 5,ncol = 1)
    k <- numeric(length = Nwl)

    Cbrown <- 0

    cc[1,1] <- Cab / N
    cc[2,1] <- Car / N
    cc[3,1] <- Cbrown / N
    cc[4,1] <- Cw / N
    cc[5,1] <- Cm / N

    k[] <- dataspec_p5[,] %*% cc[,]

    trans <- (1 - k)*exp(-k) + k^2 *e1_approx(k)
    trans[trans < 0] <- 0

    ralf <- 1 - talf
    r12 <- 1 - t12
    r21 <- 1 - t21

    denom <- 1 - (r21 ^ 2) * (trans ^ 2)
    Ta <- talf * trans * t21 / denom
    Ra <- ralf + r21 * trans * Ta

    tt <- t12 * trans * t21 / denom
    rr <- r12 + r21 * trans * tt

    gt1 <- rr + tt >= 1
    tgt1 <- tt[gt1]

    Tsub <- 0*tt
    Rsub <- 0*tt
    r <- 0*tt
    t <- 0*tt

    Tsub[gt1] <- tgt1 / (tgt1 + (1 - tgt1) * (N - 1))
    Rsub[gt1] <- 1 - Tsub[gt1]

    inf <- rr == 0 | tt == 0
    Tsub[inf] <- 0
    Rsub[inf] <- 0

    r <- rr[!gt1 & !inf]
    t <- tt[!gt1 & !inf]

    D <- sqrt((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t))
    r2 <- r ^ 2
    t2 <- t ^ 2
    va <- (1 + r2 - t2 + D) / (2 * r)
    vb <- (1 - r2 + t2 + D) / (2 * t)

    vbNN <- vb ^ (N - 1)
    vbNN2 <- vbNN ^ 2
    va2 <- va ^ 2
    denomx <- va2 * vbNN2 - 1
    Rsub[!gt1 & !inf] <- va * (vbNN2 - 1) / denomx
    Tsub[!gt1 & !inf] <- vbNN * (va2 - 1) / denomx

    denomy <- 1 - Rsub * rr

    reflectance <- Ra + Ta * Rsub * tt / denomy

    returnType(double(1))

    if (N < 1.1 | Car < 0 | Cab < 0 | Cm <= 0 | Cw <0) {
      return(-9999*(reflectance**0))
    } else{
      return(reflectance)
    }
  })

# NIMBLE main model
run_prospect5 <- nimbleCode({

  # Run NIMprospect5 with mean parameters
  reflectance[1:Nwl] <- NIMprospect5(N,
                                     Cab,
                                     Car,
                                     Cw,
                                     Cm,
                                     dataspec_p5[,], talf[],t12[],t21[], Nwl)

  for (i in 1:Nleaves){
    for (j in 1:Nwl){
      obs_reflectance[j,i] ~ dnorm(reflectance[j], sd = Standard.Dev)
    }
  }


  # Priors
  ## Mean (Tree/PNM)
  N ~ dunif(1.,5)
  Cab ~ dunif(0,100)
  Car ~ dunif(0,50)
  Cw ~ dunif(0.,0.1)
  Cm ~ dunif(0.,0.1)
  Standard.Dev ~ dunif(0,1)

})

# List of data (obs_reflectance)
Data <- list(obs_reflectance = t(matrix(data %>% filter(wv %in% WLs) %>% pull(value),ncol = Nwl)))
colnames(Data$obs_reflectance) <- NULL

# Plot the mean data for liana (blue) and tree (green) for PNM (solid) and FTS (dash) sites
par(mfrow = c(1,1))
matplot(WLs,Data$obs_reflectance,type = 'l',col = "darkgreen")

# Define constants for NIMBLE
Constants <- list(Nwl = Nwl,
                  Nleaves = Nleaves,
                  talf = rrtm:::p45_talf[pos],
                  t12 = rrtm:::p45_t12[pos],
                  t21 = rrtm:::p45_t21[pos],
                  dataspec_p5 = rrtm:::dataspec_p5[pos,1:5])

# To test the model, give initial values
Inits <- list(Nmean = 2,
              Cabmean = 40,
              Carmean = 10,
              Cwmean = 0.01,
              Cmmean = 0.01,
              Standard.Dev = 0.05)

# Build the PROSPECT5-NIMBLE model
P5model <- nimbleModel(run_prospect5,
                       dimensions = list(dataspec_p5 = c(Nwl,5),
                                         talf = Nwl,
                                         t12 = Nwl,
                                         t21 = Nwl,
                                         reflectance = c(Nwl)),
                       data = Data,
                       constants = Constants,
                       debug = FALSE,
                       inits = Inits)

# Check if initialized properly
P5model$initializeInfo()

mcmc.out <- nimbleMCMC(code = P5model,
                       constants = Constants,
                       monitors = c("N","Cab","Car","Cw","Cm","Standard.Dev"),       # which nodes to monitor?
                       data = Data,
                       inits = Inits,
                       nburnin = Nburnin,
                       nchains = Nchains,
                       niter = Niter,
                       thin = thin,
                       summary = TRUE,
                       WAIC = TRUE,
                       samplesAsCodaMCMC = TRUE)

MCMCsamples <- mcmc.out$samples
param <- MCMCsamples[,]

if (Nchains == 1){
  param_all <- MCMCsamples[,1:6]
} else {
  param_all <- do.call(rbind,lapply(1:Nchains,function(i) MCMCsamples[[i]]))
}


Simu.Nimble <- matrix(unlist(lapply(1:nrow(param_all),function(isimu){
  rrtm::prospect5(N = param_all[isimu,5],
                  Cab = param_all[isimu,1],
                  Car = param_all[isimu,2],
                  Cbrown = 0,
                  Cw = param_all[isimu,4],
                  Cm = param_all[isimu,3])[["reflectance"]][pos]})),ncol = nrow(param_all))

param.Nimble <- as.data.frame(as.matrix(param)) %>%
  pivot_longer(cols = c("N","Cab","Car","Cw","Cm","Standard.Dev"),
               names_to = "param",
               values_to = "value")


#################################################################################################################
# Bayesian Tools

run_prospect <- function(params,waves){

  # PROSPECT parameters
  Nlayers <- params[1]
  Cab <- params[2]
  Car <- params[3]
  Cw <- params[4]
  Cm <- params[5]


  # Call RTM
  result <- tryCatch(
    Rprospect::prospect5(Nlayers,Cab,Car,Cw,Cm),
    error = function(e) NULL)
  if (is.null(result)) return(-1e20)
  reflectance <- result[["Reflectance"]]
  if (any(!is.finite(reflectance))) return(-1e20)
  if (any(reflectance < 0)) return(-1e20)
  if (any(reflectance > 1)) return(-1e20)

  reflectance_waves <- interp1(x = result[["Wavelength"]],
                               y = reflectance,
                               waves)
  return(reflectance_waves)
}

create_likelihood <- function(observed, waves) {
  function(params) {

    ssigma <- params[6]
    reflectance_waves <- run_prospect(params,waves)

    # Calculate likelihood
    ll <- 0

    for (i in seq(1,ncol(observed))){
      ll <- ll + sum(dnorm(reflectance_waves, observed[,i], ssigma, log = TRUE))
    }

    return(ll)
  }
}

Prospect_param_names <- c("N","Cab","Car","Cw","Cm","Standard.Dev")
pft_lowers <- c(N = 1, Cab = 0 , Car = 0,Cw = 0, Cm = 0, Standard.Dev = 0)
pft_uppers <-  c(N = 5, Cab = 100, Car = 50,Cw = 0.1, Cm = 0.1, Standard.Dev = 1)
params <- (pft_lowers + pft_uppers)/2

prior <- createUniformPrior(pft_lowers, pft_uppers)
likelihood <- create_likelihood(observed = Data$obs_reflectance,waves = WLs)
settings_MCMC <- list(iterations = Niter, nrChains = Nchains)

# Run inversion
setup <- BayesianTools::createBayesianSetup(likelihood, prior, parallel = FALSE)
samples <- BayesianTools::runMCMC(setup,settings = settings_MCMC)
# samples <- BayesianTools::runMCMC(samples,settings = settings_MCMC)

param_all_BT <- getSample(samples,start = Nburnin,thin = thin)
colnames(param_all_BT) <- Prospect_param_names

param_BT <- as.data.frame(param_all_BT) %>% pivot_longer(cols = c(Prospect_param_names),
                                                     names_to = "param",
                                                     values_to = "value")

Simu.BT <- matrix(unlist(lapply(1:nrow(param_all_BT),function(isimu){
  rrtm::prospect5(N = param_all_BT[isimu,1],
                  Cab = param_all_BT[isimu,2],
                  Car = param_all_BT[isimu,3],
                  Cbrown = 0,
                  Cw = param_all_BT[isimu,4],
                  Cm = param_all_BT[isimu,5])[["reflectance"]][pos]})),ncol = nrow(param_all_BT))

#################################################################################################################################
# Plotting outputs

# Parameter distributions
df_all <- bind_rows(list(param.Nimble %>% mutate(type = "NIMBLE"),
                         param_BT %>% mutate(type = "Bayes.Tools")))

ggplot(df_all) +
  geom_density_ridges(aes(x = value, y = type, fill = param)) +
  facet_wrap(~ param, scales = "free") +
  theme_bw() +
  labs(x = "", y = "") +
  guides(fill = FALSE)


# Spectra

matplot(WLs,Data$obs_reflectance,type = 'l',col = "darkgrey",ylim = c(0,0.7))
matlines(WLs,rowMeans(Simu.BT),col = "red",lty = 1)
matlines(WLs,rowMeans(Simu.Nimble),col = "blue",lty = 2)

