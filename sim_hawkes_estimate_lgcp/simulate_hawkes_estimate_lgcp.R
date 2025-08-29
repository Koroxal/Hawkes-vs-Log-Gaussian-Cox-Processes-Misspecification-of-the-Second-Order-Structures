########################################################################################
#  Simulating from a Hawkes process
########################################################################################

# Load required packages

library(ggplot2)
library(MASS)
library(dplyr)
library(knitr)
library(stpp)
library(lgcp)
library(spatstat)
library(inlabru)
library(INLA)
library(parallel)
library(ggplot2)
library(ggnewscale)
library(profvis)
library(viridis)
library(fields)
library(geoR)

# Code to simulate a generic Hawkes process using the acceptance–rejection method

# Load functions for generating Hawkes process data (Model 1)

source("functions_generate_data_model1_hawkes.R")

# Background intensity (baseline rate of events)
mu <- 3

# Define the kernel parameters for the Hawkes process
# C: excitation constant, w: spatial decay parameter, beta: temporal decay parameter
kernel <- list(C = .75, w = 1, beta = .02)

# Combine the baseline intensity and kernel parameters, and set a maximum allowed intensity
lam <- list(mu = mu, kernel = kernel, maximum = 1e3)

# Simulate the Hawkes process:
# T: time interval [0, 100]
# S: spatial domain (square [0,1] x [0,1])
# batch_size: number of simulations 
# min_n_points: minimum number of points to generate
# verbose: display progress in the console
pp <- generate(lam, T = c(0, 100), S = matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE),
               batch_size = 1, min_n_points = 5, verbose = TRUE)

# Remove all objects from the environment (clean-up)
rm(list=ls())

########################################################################################
# 3.3 Case study 1: Simulating from a Hawkes process, estimating with an LGCP
########################################################################################

# We load the dataset used in our simulation study

load("data_hawkes_1_002.RData")

# Plot reference data set (Figure 1(a))

ggplot() + 
  geom_point(data = pp_data_1_002, aes(x = x, y = y, color = times), size = 2.3) +  
  scale_color_gradient(
    name = "time", 
    breaks = c(0, 25, 50, 75, 100), 
    low = "#56B1F7", 
    high = "#132B43", 
    limits = c(0, 100),  
    guide = guide_colorbar(
      barwidth = 2,     
      barheight = 15,   
      title.theme = element_text(size = 20),  
      label.theme = element_text(size = 18)   
    )
  ) + 
  new_scale_color() + 
  labs(x = "x", y = "y") +
  theme(
    legend.title = element_text(size = 28),   
    legend.text = element_text(size = 18),    
    axis.title = element_text(size = 20),     
    axis.text = element_text(size = 18)        
  )

data_hawkes <- data.frame(x=pp_data_1_002$x, y=pp_data_1_002$y, t=pp_data_1_002$times)

# Spatial domain (unit square)

spatial_window <- owin(c(0, 1), c(0, 1))
win <- spatial_window

# Temporal domain

tlim <- c(0, 100)

# Data stppp

xyt_hawkes <- stppp(list(data = data_hawkes, tlim = tlim, window = win))
xyt_hawkes <- integerise(xyt_hawkes)

# Spatial intensity

den <- density.ppp(xyt_hawkes)
plot(den)
sar <- spatialAtRisk(den)
plot(sar)

# Temporal intensity

mut1 <- muEst(xyt_hawkes)
tar <- mut1
plot(tar)

# Compute the k-function from data

kin_hawk_original <- KinhomAverage(xyt_hawkes, sar, tar)
plot(kin_hawk_original)

set.seed(4353)

# Estimate parameters with an LGCP model (using an exponential covariance function for both time and space)

# Load required LGCP functions

source("lgcp_functions.R")

# Parameter estimation using the minimum contrast ad-hoc method
# Model 3 (exponential covariance function)

parameters1 <- minimum.contrast.spatiotemporal.adhoc(
  data = xyt_hawkes,           # Spatio-temporal point pattern data (from Hawkes process)
  model = "exponential",       # Exponential covariance model for time and space
  method = "g",                # Use the g-function (pair correlation) for estimation
  spatial.dens = sar,          # Estimated spatial intensity
  temporal.intens = tar,       # Estimated temporal intensity
  power = 0.25,                 # Weighting exponent for the contrast criterion
  transform = NULL,            # No transformation applied
  temporal.interval = c(.1, 5), # Interval of temporal lags to consider
  verbose = TRUE               # Show estimation progress
)

parameters1
phi <- parameters1$estimates[1]
phi_pr <- practicalRange(cov.model = "matern", kappa = 0.5, phi = round(parameters1$estimates[1], 3))
theta <- parameters1$estimates[3]
sigma <- parameters1$estimates[2]

rbind(sigma=sigma, phi=phi, phi_pr=phi_pr, theta=1/theta) %>% kable(., caption ="Estimates Model 3 (Exponential)")

# phi_pr is the practical range derived from phi, calculated from the phi parameter, using the practicalRange function with the "matern" covariance model and kappa = 0.5 (Exponential covariance function).
# theta is computed as the inverse of the third estimated parameter (1 / theta), to ensure consistency with the parameterisation used in the paper.

# Model 4 (Matérn covariance function)

parameters2 <- minimum.contrast.spatiotemporal.adhoc2(
  data = xyt_hawkes,           # Spatio-temporal point pattern data (from Hawkes process)
  model = "matern",            # Matern covariance model for space
  additionalparameters=list(nu=3/2), #set v=3/2
  method = "g",                # Use the g-function (pair correlation) for estimation
  spatial.dens = sar,          # Estimated spatial intensity
  temporal.intens = tar,       # Estimated temporal intensity
  power = 0.25,                 # Weighting exponent for the contrast criterion
  transform = NULL,            # No transformation applied
  temporal.interval = c(.1, 5), # Interval of temporal lags to consider
  verbose = TRUE               # Show estimation progress
)

parameters2
phi2 <- parameters2$estimates[1]
phi2_pr <- practicalRange(cov.model = "matern", kappa = 3/2, phi = round(parameters2$estimates[1], 3))
theta2 <- parameters2$estimates[3]
sigma2 <- parameters2$estimates[2]

rbind(sigma=sigma2, phi=phi2, phi_pr=phi2_pr, theta=1/theta2) %>% kable(., caption ="Estimates Model 4 (Matérn)")

# Table estimates for Model 3 and 4

combined_table <- data.frame(
  "Model 3" = c(sigma, phi, phi_pr, 1/theta),
  "Model 4" = c(sigma2, phi2, phi2_pr, 1/theta2),
  row.names = c("sigma", "phi", "phi_pr", "theta")
)
kable(combined_table, caption = "Estimates: Model 3 (Exponential) vs Model 4 (Matérn)")

# Bayesian MCMC estimation is performed below, yielding the same results; 
# however, the computation time is significantly longer.

# Spatial window and temporal limits
win <- owin(c(0, 1), c(0, 1))
tlim <- c(0, 100)

# Select cell width
chooseCellwidth(xyt_hawkes, cwinit = 0.05)
CELLWIDTH <- 0.05
  
# Select extension
EXT <- 2
  
# Perform polygon overlay operations and compute computational grid
polyolay <- getpolyol(data = xyt_hawkes, regionalcovariates = NULL, cellwidth = CELLWIDTH, ext = EXT)
plot(polyolay$fftpoly)
  
# Define a formula without covariates
FORM <- X ~ 1
  
# Define priors for LGCP parameters 
lgprior <- PriorSpec(LogGaussianPrior(mean = log(c(2, .03,1/3)), variance = diag(c(0.2, .02, .2))))
lgprior <- PriorSpec(LogGaussianPrior(mean = log(c(2, .03,3)), variance = diag(c(0.5, .5, 1))))
gprior <- PriorSpec(GaussianPrior(mean = rep(0, 1), variance = diag(5, 1))) # Adjusted for no covariates
priors <- lgcpPrior(etaprior = lgprior, betaprior = gprior)
  
# Set initial values for MCMC
INITS <- lgcpInits(etainit = log(c(2, 0.03, 1/3)), betainit = NULL)
  
# Choose the covariance function (maternCovFct15 and maternCovFct25 are also available)
CF <- CovFunction(exponentialCovFct)
DIRNAME <- getwd()

# Perform Bayesian MCMC estimation
# lg <- lgcpPredictSpatioTemporalPlusPars(
#     formula = FORM, xyt = xyt_hawkes, T = 100, laglength = 10,
#     ZmatList = NULL, model.priors = priors, model.inits = NULL,
#     spatial.covmodel = CF, cellwidth = CELLWIDTH,
#     mcmc.control = mcmcpars(
#       mala.length = 500000, burnin = 50000, retain = 2000, # It is necessary to increase the number of iterations and parameters of the MCMC
#       adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, targetacceptance = 0.574)
#     ),
#     output.control = setoutput(
#       gridfunction = dump2dir(dirname = file.path(DIRNAME, paste0("aegiss_iter__", 1)), forceSave = TRUE)
#     ),
#     ext = EXT
#   )


lg2 <- lgcpPredictSpatioTemporalPlusPars(
  formula = FORM, xyt = xyt_hawkes, T = 100, laglength = 10,
  ZmatList = NULL, model.priors = priors, model.inits = NULL,
  spatial.covmodel = CF, cellwidth = CELLWIDTH,
  mcmc.control = mcmcpars(
    mala.length = 800000, burnin = 80000, retain = 200, # It is necessary to increase the number of iterations and parameters of the MCMC
    adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, targetacceptance = 0.574)
  ),
  output.control = setoutput(
    gridfunction = dump2dir(dirname = file.path(DIRNAME, paste0("aegiss_iter__", 1)), forceSave = TRUE)
  ),
  ext = EXT
)


# Diagnostic plots for MCMC chains

par(mfrow=c(2,2))

# produce traceplots
traceplots(lg2, ask = FALSE)

# produce autocorrelation plots
parautocorr(lg2, ask = FALSE)

# summary table of posterior samples
parsummary(lg2)


########################################################################################
# K functions
########################################################################################

# Load required LGCP functions

source("lgcp_functions.R")

data_lgcps <- list()
set.seed(4353)

# Simulate 300 LGCP datasets with estimated parameters (Model 3)
# For Model 4, the procedure is the same but changing the covariance function

for(n_sim in 1:300){
  print(n_sim)
  A <- lgcpSim(
    owin = win,
    tlim = attr(tar, "tlim"),
    spatial.intensity = sar,
    temporal.intensity = tar,
    cellwidth = 0.04,
    model.parameters = lgcppars(sigma = sigma, phi = phi, theta = theta),
    spatial.covmodel = "exponential",
    covpars = c(),
    returnintensities = TRUE,
    progressbar = TRUE,
    ext = 2,
    plot = FALSE,
    ratepow = 0.25,sleeptime = 0,
    inclusion = "touching"
  )
  
  data <- data.frame(x=A$x, y=A$y, t=A$t)
  data_lgcps[[n_sim]] <- data
  print(nrow(data))
}

# output: 
load("data_k_functions_model3_1_002.RData")

# Plot Hawkes process Kinhom function

plot(kin_hawk_original$r, kin_hawk_original$iso, type="l", ylim=c(0, 0.3), xlab="r", ylab=expression(g[inhom](r)))
lines(kin_hawk_original$r, kin_hawk_original$theo, lty=2, col="black")

# Compute Kinhom for each LGCP simulation and plot (evelopes)

kas_lgcp <- list()
tlim <- c(0, 100)
win <- owin(c(0,1), c(0,1))
for(n_bd in 1:300){
  data_p <- data_lgcps[[n_bd]]
  
  data_lgcp <- data.frame(x=data_p$x, y=data_p$y, t=data_p$t)
  xyt_lgcp <- stppp(list(data = data_lgcp, tlim = tlim, window = win))
  xyt_lgcp <- integerise(xyt_lgcp)
  
  den_h <- density(xyt_lgcp)
  sar_h <- spatialAtRisk(den_h)
  mut1_h <- muEst(xyt_lgcp)
  tar_h <- mut1_h

  kin_lgcp <- KinhomAverage(xyt_lgcp, sar_h, tar_h)
  
  kas_lgcp[[n_bd]] <- kin_lgcp
  lines(kin_lgcp$r, kin_lgcp$iso, col="grey")
}

lines(kin_hawk_original$r, kin_hawk_original$iso, type="l", ylim=c(0, 0.6))
lines(kin_hawk_original$r, kin_hawk_original$theo, lty=2, col="black")

legend("topleft",c("iso (Hawkes)", "theo", "iso(LGCP simul)"), col = c("black","black",  "grey"),lty = c(1, 6, 2, 2, 3, 1), cex=.8)

# Remove extreme curves (top and bottom 5% based on mean iso value)

# Compute mean iso values for each LGCP Kinhom
mean_iso_values <- sapply(kas_lgcp, function(kin_hawk) {
  mean(kin_hawk$iso, na.rm = TRUE) # Calculamos la media de los valores "iso"
})

# Sort by mean iso
sorted_indices <- order(mean_iso_values)

# Identify lower and upper 5% (15 curves each out of 300)
lower_indices <- sorted_indices[1:15]
upper_indices <- sorted_indices[(length(sorted_indices) - 14):length(sorted_indices)]

# Keep remaining curves
remaining_indices <- sorted_indices[16:(length(sorted_indices) - 3)]
kas_lgcp_filtered <- kas_lgcp[remaining_indices]

# Plot evelopes bands and k functuon from Hawkes process
plot(kin_hawk_original$r, kin_hawk_original$iso, type = "l", ylim = c(0, 0.3),
     xlab = "r", ylab = expression(K[inhom](r)),
     cex.axis = 1.2,
     cex.lab = 1.5,
     yaxt = "n" 
)

axis(2, at = c(0.0, 0.1, 0.2, 0.3), labels = c("0.0", "0.1", "0.2", "0.3"), cex.axis = 1.2)    
lines(kin_hawk_original$r, kin_hawk_original$theo, lty=2, col="black")

for (i in remaining_indices) {
  kin_lgcp <- kas_lgcp[[i]]
  lines(kin_lgcp$r, kin_lgcp$iso, col="lightgrey")
}

lines(kin_hawk_original$r, kin_hawk_original$iso, type="l", ylim=c(0, 0.3))
lines(kin_hawk_original$r, kin_hawk_original$theo, lty=2, col="black")

# Plot together K-functions from Model 3 and Model 4 (Figure 2(a))

load("k_m3_1002.RData"); kas_lgcp_m3 <- kas_hawk
load("k_m4_1002.RData"); kas_lgcp_m4 <- kas_hawk

# Function to filter curves 
filter_curves <- function(kas_hawk) {
  mean_iso_values <- sapply(kas_hawk, function(kin_hawk) mean(kin_hawk$iso, na.rm = TRUE))
  sorted_indices <- order(mean_iso_values)
  remaining_indices <- sorted_indices[16:(length(sorted_indices) - 1)]
  kas_hawk[remaining_indices]
}

# Filter both datasets
kas_lgcp_filtered_m3 <- filter_curves(kas_lgcp_m3)
kas_lgcp_filtered_m4 <- filter_curves(kas_lgcp_m4)

# Set margins: bottom, left, top, right
par(mar = c(5, 5, 2, 2))  
# Plot the original K-function
plot(kin_hawk_original$r, kin_hawk_original$iso, type = "l", ylim = c(0, 0.3),
     xlab = "r", ylab = expression(K[inhom](r)),
     cex.axis = 1.2,
     cex.lab = 1.5,
     yaxt = "n"  # Disable y-axis ticks and labels
)
# Add custom y-axis labels
axis(2, at = c(0.0, 0.1, 0.2, 0.3), labels = c("0.0", "0.1", "0.2", "0.3"), cex.axis = 1.2)    

# Function to calculate min and max range for shading
calc_minmax <- function(kas_hawk_filtered) {
  r_vals <- kas_hawk_filtered[[1]]$r
  iso_mat <- sapply(kas_hawk_filtered, function(kh) kh$iso)
  
  lower <- apply(iso_mat, 1, min, na.rm = TRUE)
  upper <- apply(iso_mat, 1, max, na.rm = TRUE)
  
  list(r = r_vals, lower = lower, upper = upper)
}

# Calculate ranges for Model 3 and Model 4
range_m3 <- calc_minmax(kas_lgcp_filtered_m3)
range_m4 <- calc_minmax(kas_lgcp_filtered_m4)

# Shade area for Model 3 
polygon(c(range_m3$r, rev(range_m3$r)),
        c(range_m3$lower, rev(range_m3$upper)),
        col = adjustcolor("red", alpha.f = 0.4), border = NA)

# Shade area for Model 4
polygon(c(range_m4$r, rev(range_m4$r)),
        c(range_m4$lower, rev(range_m4$upper)),
        col = adjustcolor("blue", alpha.f = 0.4), border = NA)

# Plot original curve in black
lines(kin_hawk_original$r, kin_hawk_original$iso, col = "black", lwd = 2)

# Plot theoretical K-function as dashed line
lines(kin_hawk_original$r, kin_hawk_original$theo, lty = 2, col = "black")

# Add legend
legend("topleft",
       legend = c("Model 3", "Model 4"),
       col = c(adjustcolor("red", alpha.f = 0.4),
               adjustcolor("blue", alpha.f = 0.4)),
       lty = c(1, 1), lwd = c(5, 5), bty = "n", cex = 1.5)


########################################################################################
# Forecast
########################################################################################

set.seed(4353)

# Load observed Hawkes process dataset and preapre data

load("data_hawkes_1_002.RData")
data <- pp_data_1_002
data$t <- data$times
win <- owin(c(0,1), c(0,1))
tlim <- c(0, 100)
data <- data %>% dplyr::select(x, y, t)
xyt <- stppp(list(data = data, tlim = tlim, window = win))
xyt <- integerise(xyt)

# Estimate spatial and temporal intensities

den <- density(xyt)
plot(den)
sar <- spatialAtRisk(den)
mut1 <- muEst(xyt)
mut1
mut <- constantInTime(xyt)
tar <- mut1

# Define exceedance probability thresholds

exceed <- exceedProbs(c(1.5, 2, 3))
fun2 <- function(x){return(x^2)} # computes E(X^2). Can be used with the
gfun <- function(Y){return(Y)} # gives the mean
exceed <- exceedProbs(1.000017)
tmpdr <- tempdir()
data$t <- floor(data$t)
n.t <- max(data$t)
n.t.pred <- 10
n.lag <- 20

# Run LGCP prediction for Model 3 with MCMC (For Model 4, the procedure is the same but changing the covariance function)

lg <- lgcpPredict(
  xyt = xyt,
  T = as.integer(90),  
  laglength = as.integer(n.lag),  
  model.parameters = lgcppars(sigma = sigma, phi = phi, theta = theta),
  cellwidth = 0.05,
  spatial.intensity = sar,
  temporal.intensity = tar,
  mcmc.control = mcmcpars(
    mala.length = 12000,
    burnin = 2000,
    retain = 100,
    adaptivescheme = andrieuthomsh(
      inith = 1,
      alpha = 0.5,
      C = 1,
      targetacceptance = 0.574
    )
  ),
  output.control = setoutput(
    gridfunction = dump2dir(dirname = tmpdr),
    gridmeans = MonteCarloAverage(c("gfun", "fun2", "exceed"))
  )
)

# output
load("forecast_lg_m3_1002.RData")
lg

# Prepare colour scales for visualisation

ncol = 100
colores <- plasma(ncol)
posiciones <- seq(0, 1, length.out = ncol)
posiciones_modificadas <- posiciones^(1/1.5)
colores_modificados <- colores[as.integer(posiciones_modificadas * (ncol - 1)) + 1]

# Create binary mask for plotting

ff <- discreteWindow(lg)
CIm <- function(riesgolgcp, ff, numeropred){ 
  X <- riesgolgcp$grid[[numeropred]]
  X[!ff] <- NA
  return(im(t(X), xcol = riesgolgcp$xvals, yrow = riesgolgcp$yvals))
}

# Forecasting over future time steps

fcast <- lgcpForecast(lg, 91:100, spatial.intensity = sar, temporal.intensity = tar)
Famba <- fcast$intensity

# Sum predicted intensities for the first 9 forecast days

sum_intens <- (Famba$grid[[1]] + Famba$grid[[2]] + Famba$grid[[3]] + Famba$grid[[4]] + Famba$grid[[5]] +
           Famba$grid[[6]] + Famba$grid[[7]] + Famba$grid[[8]] + Famba$grid[[9]])

# Determine z-limits for plotting

zlim_90 <- range(as.vector(sum_intens), na.rm = TRUE)
zlim_9091 <- range(as.vector(CIm(Famba, ff, 1)$v), na.rm = TRUE)
zlim_9192 <- range(as.vector(CIm(Famba, ff, 2)$v), na.rm = TRUE)

# Plot cumulative forecast (days 91–99)

par(mar = c(1, 1, 1, 1))
sum_intens <- (Famba$grid[[1]] + Famba$grid[[2]] + Famba$grid[[3]] + Famba$grid[[4]] + Famba$grid[[5]] +
           Famba$grid[[6]] + Famba$grid[[7]] + Famba$grid[[8]] + Famba$grid[[9]])
X <- sum_intens
X[!ff] <- NA
image(im(t(X), xcol = Famba$xvals, yrow = Famba$yvals), col = colores_modificados, xlab = "", ylab = "", main = "", zlim = zlim_90)
plot(unmark(xyt[xyt$t >= 90 & xyt$t <= 100]), add = TRUE, chars = 16, pch=21, col = "black", bg = "white",cex = 1.2)

# Plot daily forecasts with observed events

par(mar = c(1, 1, 1, 1))
i <- 1
kk <- 89 + i
image(CIm(Famba, ff, i), col = colores_modificados, xlab = "", ylab = "", main = "", zlim = zlim_9091)
plot(unmark(xyt[xyt$t == kk]), add = TRUE, chars = 16, , pch=21, col = "black", bg = "white",cex = 1.2)

par(mar = c(1, 1, 1, 1))
i <- 2
kk <- 89 + i
image(CIm(Famba, ff, i), col = colores_modificados, xlab = "", ylab = "", main = "", zlim = zlim_9192)
plot(unmark(xyt[xyt$t == kk]), add = TRUE, chars = 16, pch=21, col = "black", bg = "white",cex = 1.2)

# The legends in the figures of the associated paper are slightly different because they were designed to display two models being compared, both of which needed to be represented uniformly. 
# However, this script focuses on reproducing results for a single model, to avoid overloading the code and to keep the process clear, rather than complicating it with additional details and replicates. 
# If more information are required, please contact the authors (bernebeu@uji.es).

