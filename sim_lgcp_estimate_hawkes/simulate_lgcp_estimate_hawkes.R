########################################################################################
#  Simulating from an LGCP
########################################################################################

# Load required packages
library(lgcp)
library(spatstat.geom)
library(ggplot2)
library(ggnewscale)
library(spatstat)
library(stpp)
library(inlabru)
library(INLA)
library(dplyr)
library(parallel)
library(viridisLite)
library(fields)

# ------------------------------------------------------------------
# Model: Simulating from a spatio-temporal LGCP with:
#   λ(x, t) = exp(S(x, t)),  where
#   S(x, t) ~ GP(μ, Σ) with separable covariance in space and time
#
# The expected number of points is:
#   E[N] = |D|exp(μ + σ^2 / 2)
# where |D| represents the volume of the space-time observation domain
# it is the product of the spatial area and the temporal length
# ------------------------------------------------------------------

# Desired expected number of points

target_E_N <- 1000

# Spatial domain (unit square by default)

spatial_window <- owin(c(0, 1), c(0, 1))
win <- spatial_window

# Temporal domain

tlim <- c(0, 100)

# Parameters of the spatio-temporal Gaussian process

sigma <- 2   # marginal standard deviation
phi <- 0.03     # spatial range
theta <- 3    # temporal correlation parameter

# Grid resolution (adjust if needed)

cellwidth <- 1/32

# Compute mu to achieve the desired expected count

volume <- area.owin(spatial_window) * diff(tlim)
sigma2 <- sigma^2
mu <- log(target_E_N / volume) - sigma2 / 2

# Define constant temporal intensity

temporal.intensity <- constantInTime(1, tlim = tlim)

# Simulate one LGCP realization

sim <- lgcpSim(
  owin = spatial_window,
  tlim = tlim,
  cellwidth = cellwidth,
  model.parameters = lgcppars(sigma = sigma, phi = phi, theta = theta, mu = mu),
  spatial.covmodel = "exponential",
  spatial.intensity = NULL,        
  temporal.intensity = temporal.intensity,
  ext = 2,
  returnintensities = TRUE,
  plot = FALSE,
  progressbar = TRUE
)

# Number of simulated points

sim$n

# Plot generated data

data <- data.frame(x=sim$x, y=sim$y, t=sim$t)

ggplot() +
  geom_point(data = data, aes(x = x, y = y, color = t), size = 2.3) +  # Aumenta el tamaño de los puntos
  scale_color_gradient(
    name = "time",
    breaks = c(0, 25, 50, 75, 100),
    low = "#56B1F7",
    high = "#132B43",
    limits = c(0, 100),  # Especificamos el rango de la leyenda
    guide = guide_colorbar(
      barwidth = 2,     # Ancho de la barra
      barheight = 15,   # Alto de la barra (aumentar este valor hace la barra más grande)
      title.theme = element_text(size = 20),  # Tamaño del título de la barra
      label.theme = element_text(size = 18)   # Tamaño del texto de la barra
    )
  ) +
  new_scale_color() +
  labs(x = "x", y = "y") +
  theme(
    legend.title = element_text(size = 28),    # Aumenta el tamaño del título de la leyenda
    legend.text = element_text(size = 18),     # Aumenta el tamaño del texto de la leyenda
    axis.title = element_text(size = 20),      # Aumenta el tamaño de los títulos de los ejes
    axis.text = element_text(size = 18)        # Aumenta el tamaño de los números de los ejes
  )


########################################################################################
# 3.2 Model misspecification: Simulating from an LGCP, estimating with a Hawkes process
########################################################################################

# We load the dataset used in our simulation study

load("data_from_lgcp_3003.RData")
data_RAW <- data

# Plot reference data set

ggplot() +
  geom_point(data = data, aes(x = x, y = y, color = t), size = 2.3) +  # Aumenta el tamaño de los puntos
  scale_color_gradient(
    name = "time",
    breaks = c(0, 25, 50, 75, 100),
    low = "#56B1F7",
    high = "#132B43",
    limits = c(0, 100),  # Especificamos el rango de la leyenda
    guide = guide_colorbar(
      barwidth = 2,     # Ancho de la barra
      barheight = 15,   # Alto de la barra (aumentar este valor hace la barra más grande)
      title.theme = element_text(size = 20),  # Tamaño del título de la barra
      label.theme = element_text(size = 18)   # Tamaño del texto de la barra
    )
  ) +
  new_scale_color() +
  labs(x = "x", y = "y") +
  theme(
    legend.title = element_text(size = 28),    # Aumenta el tamaño del título de la leyenda
    legend.text = element_text(size = 18),     # Aumenta el tamaño del texto de la leyenda
    axis.title = element_text(size = 20),      # Aumenta el tamaño de los títulos de los ejes
    axis.text = element_text(size = 18)        # Aumenta el tamaño de los números de los ejes
  )

data_lgcp <- data.frame(x=data$x, y=data$y, t=data$t)

# Spatial domain (unit square by default)

spatial_window <- owin(c(0, 1), c(0, 1))
win <- spatial_window

# Temporal domain

tlim <- c(0, 100)

# Data stppp

xyt_lgcp <- stppp(list(data = data_lgcp, tlim = tlim, window = win))
xyt_lgcp <- integerise(xyt_lgcp)

# Spatial intensity
den <- density.ppp(xyt_lgcp)
plot(den)
sar <- spatialAtRisk(den)
plot(sar)

# Temporal intensity
mut1 <- muEst(xyt_lgcp)
tar <- mut1
plot(tar)

# Compute the k-function from data
kin_lgcp_original <- KinhomAverage(xyt_lgcp, sar, tar)
plot(kin_lgcp_original)

set.seed(4353)

# Estimate with Hawkes model (gaussian triggering for space and temporal for time)

source("estim_model1_hawkes.R")

#inlabru convergence

ggplot(fit_etas$bru_iinla$track, aes(x = iteration, y = mode)) +
  geom_line() +
  facet_wrap(facets = vars(effect), scales = 'free')

inlabru:::make_track_plots(fit_etas)$default

#results

rbind(
  inlabru_0.025 = internal_to_natural(param = fit_etas$summary.fixed$"0.025quant", link.functions = link.functions),
  inlabru = internal_to_natural(param = fit_etas$summary.fixed$mean, link.functions = link.functions),
  inlabru_0.975 = internal_to_natural(param = fit_etas$summary.fixed$"0.975quant", link.functions = link.functions)
)

########################################################################################
# K functions
########################################################################################

param_inla <- internal_to_natural(param = fit_etas$summary.fixed$mean, link.functions = link.functions)
mu <-  param_inla[[1]]
kernel <- list(C = param_inla[[2]], w = param_inla[[3]], beta = param_inla[[4]])

lam <- list(mu = mu, kernel = kernel, maximum = 1e3)

# We simulated 300 datasets to subsequently compute their K functions

source("functions_generate_data_model1_hawkes.R")
set.seed(4353)
pp <- generate(lam, T = c(0, 100), S = matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE),
               batch_size = 300, min_n_points = 5, verbose = TRUE)


load("hawkes_data_k_functions.RData")

points <- pp$data
sizes <- pp$sizes
kas_hawk <- list()

# K fucntion from lgcp data

plot(kin_lgcp_original$r, kin_lgcp_original$iso, type="l")
lines(kin_lgcp_original$r, kin_lgcp_original$theo, lty=2, col="black")

# calculate envelope bands (K functions from Hawkes process)

win <- owin(c(0, 1), c(0, 1))
tlim <- c(0, 100)
for(n_bd in 1:300){
  data_p <- points[n_bd,,]
  data_p <- as.data.frame(data_p)
  colnames(data_p) <- c("t", "x", "y")
  data_p <- data_p %>% filter(t != 0)
  data_p$times <- data_p$t
  data_hawk <- data.frame(x=data_p$x, y=data_p$y, t=data_p$times)
  xyt_hawk <- stppp(list(data = data_hawk, tlim = tlim, window = win))
  xyt_hawk <- integerise(xyt_hawk)
  
  den_h <- density(xyt_hawk)
  sar_h <- spatialAtRisk(den_h)
  mut1_h <- muEst(xyt_hawk)
  tar_h <- mut1_h

  kin_hawk <- KinhomAverage(xyt_hawk, sar_h, tar_h)
  kas_hawk[[n_bd]] <- kin_hawk
  lines(kin_hawk$r, kin_hawk$iso, col="grey")
}

lines(kin_lgcp_original$r, kin_lgcp_original$iso, type="l")
lines(kin_lgcp_original$r, kin_lgcp_original$theo, lty=2, col="black")

legend("topleft",c("iso (LGCP)", "theo", "iso(hwk_simul)"), col = c("black","black",  "grey"),lty = c(1, 6, 2, 2, 3, 1), cex=.8)

# We remove the top and bottom 5%. Five percent of 300 is 15. Therefore, we remove the lowest 15 and the highest 15.

#Step 1: We calculate the mean of the "iso" values

mean_iso_values <- sapply(kas_hawk, function(kin_hawk) {
  mean(kin_hawk$iso, na.rm = TRUE) 
})

# Step 2: Sort the curves by the mean of "iso"

sorted_indices <- order(mean_iso_values)

# The 5% lowest are the first 15 indices

lower_indices <- sorted_indices[1:15]

# The 5% highest are the last 15 indices

upper_indices <- sorted_indices[(length(sorted_indices) - 14):length(sorted_indices)]

# Step 3: Select all curves except the 15 lowest and the 15 highest

remaining_indices <- sorted_indices[16:(length(sorted_indices) - 15)]
kas_hawk_filtered <- kas_hawk[remaining_indices]

# plot k functions with evelopes

par(mar = c(5, 5, 2, 2))  

plot(kin_lgcp_original$r, kin_lgcp_original$iso, type = "l", ylim = c(0, 0.25),
     xlab = "r", ylab = expression(K[inhom](r)),
     cex.axis = 1.2,
     cex.lab = 1.5,
     yaxt = "n" 
)

axis(2, at = c(0.0, 0.1, 0.2, 0.3), labels = c("0.0", "0.1", "0.2", "0.3"), cex.axis = 1.2)    
lines(kin_lgcp_original$r, kin_lgcp_original$theo, lty=2, col="black")

for (i in remaining_indices) {
  kin_hawk <- kas_hawk[[i]]
  lines(kin_hawk$r, kin_hawk$iso, col="lightgrey")
}

lines(kin_lgcp_original$r, kin_lgcp_original$iso, type="l", ylim=c(0, 0.6))
lines(kin_lgcp_original$r, kin_lgcp_original$theo, lty=2, col="black")


########################################################################################
# Forecast
########################################################################################

data <- data_RAW
history <- data %>% filter(t <= 90)  # Filter data to include only times up to 90

source("functions_generate_data_model1_hawkes.R")  # Load custom function definitions
set.seed(4353)  

# Generate historical point process data using specified parameters
pp <- generate_history(lam, T = c(90, 100), S = matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE),
                       batch_size = 100, min_n_points = 5, verbose = TRUE, retained_points = as.matrix(history))

load("forecast_data_model1.RData")  # Load forecast data for model 1

ncol <- 100  # Number of colours for the palette
colores <- plasma(ncol)  # Generate plasma colour palette with 100 colours

posiciones <- seq(0, 1, length.out = ncol)  # Create sequence from 0 to 1 with 100 points
posiciones_modificadas_4 <- posiciones^(1/1.5)  # Apply a non-linear transformation to positions
# Select modified colours based on transformed positions (adjusting colour intensity)
colores_modificados_4 <- colores[as.integer(posiciones_modificadas_4 * (ncol - 1)) + 1]

library(splancs)  # Load spatial analysis package
# Define polygon vertices for spatial window (unit square)
poly <- matrix(c(0, 0,
                 1, 0,
                 1, 1,
                 0, 1,
                 0, 0), 
               ncol = 2, byrow = TRUE)

h0_silverman <- 0.05  # Bandwidth parameter for kernel density estimation (Silverman’s rule)

points <- pp$data
sizes <- pp$sizes

# Create a list to store data frames of simulated points
list_simulated_points <- vector("list", dim(points)[1])

# Loop through each row of points to convert and filter data
for(k in 1:dim(points)[1]){
  kk <- as.data.frame(points[k,,])      # Convert array slice to data frame
  colnames(kk) <- c("t", "x", "y")      # Assign column names
  kk <- kk %>% filter(t != 0)            # Filter out rows where time equals zero
  list_simulated_points[[k]] <- kk       # Store the filtered data frame in the list
}

# Combine all data frames into one
predict_simulated_points <- do.call(rbind, list_simulated_points)
predict_simulated_points <- as.data.frame(predict_simulated_points)
summary(predict_simulated_points)

# Filter points for specific time ranges
pred_90 <- predict_simulated_points %>% filter(t >= 90)
pred_90_91 <- predict_simulated_points %>% filter(t >= 90 & t < 91)
pred_91_92 <- predict_simulated_points %>% filter(t >= 91 & t < 92)

# Compute kernel density estimates for each time interval
denisty_90 <- kernel2d(as.points(pred_90), poly, h0_silverman, nx=100, ny=100)
denisty_90_91 <- kernel2d(as.points(pred_90_91), poly, h0_silverman, nx=100, ny=100)
denisty_91_92 <- kernel2d(as.points(pred_91_92), poly, h0_silverman, nx=100, ny=100)

# Define color limits for plots based on density ranges
zlim_90 <- range(denisty_90$z, na.rm = TRUE)
zlim_9091 <- range(denisty_90_91$z, na.rm = TRUE)
zlim_9192 <- range(denisty_91_92$z, na.rm = TRUE)

# Plot density for t >= 90
par(mar = c(5, 4, 4, 8))  # Adjust margins to avoid overlap
image(denisty_90, col = colores_modificados_4, main = "", xlab = "", ylab = "", zlim = zlim_90, axes=FALSE)
axis(1, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), cex.axis = 1.8)
axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), cex.axis = 1.8)
points(data_RAW %>% filter(t >= 90 ) %>% dplyr::select(x, y), pch = 19, cex=1.3)
image.plot(zlim = zlim_90, col = colores_modificados_4, legend.only = TRUE, legend.shrink = 1, legend.width = 1.5, 
           horizontal = FALSE, legend.mar = 5,
           axis.args = list(at=seq(zlim_90[1], zlim_90[2], length.out = 6)[-c(1,6)],
                            labels=floor(seq(zlim_90[1], zlim_90[2], length.out = 6)[-c(1,6)]), cex.axis=2, las=3))

# Plot density for 90 <= t < 91
par(mar = c(5, 4, 4, 8))
image(denisty_90_91, col = colores_modificados_4, main = "", xlab = "", ylab = "", zlim = zlim_9091, axes=FALSE)
axis(1, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), cex.axis = 1.8)
axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), cex.axis = 1.8)
points(data_RAW %>% filter(t >= 90 & t < 91) %>% dplyr::select(x, y), pch = 19, cex=1.3)
image.plot(zlim = zlim_9091, col = colores_modificados_4, legend.only = TRUE, legend.shrink = 1, legend.width = 1.5,
           horizontal = FALSE, legend.mar = 5,
           axis.args = list(at=seq(zlim_9091[1], zlim_9091[2], length.out = 6)[-c(1,6)],
                            labels=floor(seq(zlim_9091[1], zlim_9091[2], length.out = 6)[-c(1,6)]), cex.axis=2, las=3))

# Plot density for 91 <= t < 92
par(mar = c(5, 4, 4, 8))
image(denisty_91_92, col = colores_modificados_4, main = "", xlab = "", ylab = "", zlim = zlim_9192, axes=FALSE)
axis(1, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), cex.axis = 1.8)
axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), cex.axis = 1.8)
points(data_RAW %>% filter(t >= 91 & t < 92) %>% dplyr::select(x, y), pch = 19, cex=1.3)
image.plot(zlim = zlim_9192, col = colores_modificados_4, legend.only = TRUE, legend.shrink = 1, legend.width = 1.5,
           horizontal = FALSE, legend.mar = 5,
           axis.args = list(at=seq(zlim_9192[1], zlim_9192[2], length.out = 6)[-c(1,6)],
                            labels=floor(seq(zlim_9192[1], zlim_9192[2], length.out = 6)[-c(1,6)]), cex.axis=2, las=3))

# The legends in the figures of the associated paper are slightly different because they were designed to display two models being compared, both of which needed to be represented uniformly. 
# However, this script focuses on reproducing results for a single model, to avoid overloading the code and to keep the process clear, rather than complicating it with additional details and replicates. 
# If more information are required, please contact the authors (bernebeu@uji.es).


