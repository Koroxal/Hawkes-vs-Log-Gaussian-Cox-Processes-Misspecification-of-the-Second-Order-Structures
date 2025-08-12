
# gamma copula transformation
gamma.t <- function(x, a, b) {
  bru_forward_transformation(qgamma, x, a, b)
}
# uniform copula transformation
unif.t <- function(x, a, b) {
  bru_forward_transformation(qunif, x, min = a, max = b)
}
# log-gaussian copula transformation
loggaus.t <- function(x, m, s) {
  bru_forward_transformation(qlnorm, x, meanlog = m, sdlog = s)
}

# exponential copula transformation
exp.t <- function(x, e) {
  bru_forward_transformation(qexp, x, rate = e)
}

# inversions

# gamma copula transformation
gamma.t.inv <- function(x, a, b) {
  bru_inverse_transformation(pgamma, x, a, b)
}
# uniform copula transformation
unif.t.inv <- function(x, a, b) {
  bru_inverse_transformation(punif, x, min = a, max = b)
}
# log-gaussian copula transformation
loggaus.t.inv <- function(x, m, s) {
  bru_inverse_transformation(plnorm, x, meanlog = m, sdlog = s)
}
# exponential copula transformation
exp.t.inv <- function(x, e) {
  bru_inverse_transformation(pexp, x, rate = e)
}

# scale change

internal_to_natural <- function(param, link.functions, inverse = FALSE) {
  if (inverse) {
    values <- c(
      link.functions$mu.inv(param[1]),
      link.functions$K0.inv(param[2]),
      link.functions$w.inv(param[3]),
      link.functions$sig.inv(param[4])
    )
  } else {
    values <- c(
      link.functions$mu(param[1]),
      link.functions$K0(param[2]),
      link.functions$w(param[3]),
      link.functions$sig(param[4])
    )
  }
  if (inverse) {
    names(values) <- c("theta_mu", "theta_k0", "theta_w", "theta_sigma")
  } else {
    names(values) <- c("mu", "k0", "w", "sigma")
  }
  values
}

# priors

link.f.be <- list(
  mu = \(x) gamma.t(x, 100, 8),
  K0 = \(x) loggaus.t(x, -0.5, 0.5),
  w = \(x)  gamma.t(x, 15, 5),
  sig = \(x) loggaus.t(x, -5, 2),
  mu.inv = \(x) gamma.t.inv(x, 100, 8),
  K0.inv = \(x) loggaus.t.inv(x, -0.5, 0.5),
  w.inv = \(x)  gamma.t.inv(x, 15, 5),
  sig.inv = \(x) loggaus.t.inv(x, -5, 2)
)

# initial values

k0 <- .5
w <- 1
mu <- 10
sigma <- .01
param <-
  internal_to_natural(c(mu, k0, w, sigma),
                      link.functions = link.f.be,
                      inverse = TRUE)

param <- c(.1, -1, .1, .01) 

# functions to create spatial bins

within_region <- function(row) {
  xx <- row["xx"]
  yy <- row["yy"]
  radio <- row["xy_end"]
  if (xx - radio < 0 | xx + radio > 1 | yy - radio < 0 | yy + radio > 1) {
    return("No")
  } else {
    return("Yes")
  }
}
point_inside_band_and_square <- function(x, y, x_centre, y_centre, r) {
  sq_norm <- (x - x_centre)^2 + (y - y_centre)^2
  return(sq_norm >= r[1]^2 & sq_norm <= r[2]^2 & x <= 1 & y <= 1 & x >= 0 & y >= 0)
}

weight_integral2 <- function(row) {
  # Number of points on the grid (higher number for greater accuracy)
  n_points <- 10000
  r <- c(as.numeric(row["xy_start"]), as.numeric(row["xy_end"]))
  x_centre <- as.numeric(row["xx"])
  y_centre <- as.numeric(row["yy"])
  angulos = runif(n_points, 0, 2*pi)
  radiis = sqrt(runif(n_points, r[1]^2, r[2]^2))
  
  # Generate random points within the circular band.
  x_coords <- x_centre + radiis * cos(angulos)
  y_coords <- y_centre + radiis * sin(angulos)
  
  # Count how many points are inside the unit square.
  points_inside_band <-
    sum(point_inside_band_and_square(x_coords, y_coords, x_centre, y_centre, r))
  
  # Calculate the proportion of points that fall inside the square.
  area_band_approx <- points_inside_band / n_points
  
  return(area_band_approx)
}

# find points defining the temporal bins for an observed point

breaks_exp <- function(tt_, T2_, coef_, delta_, N_exp_ = 10) {
  tt_breaks <- tt_ + delta_ * ((1 + coef_)^(0:N_exp_))
  tt_breaks <- tt_breaks[tt_breaks < T2]
  if (T2_ - tt_ < delta_) {
    return(c(tt_, T2_))
  }
  if (T2 - tt_breaks[length(tt_breaks)] < delta_) {
    tt_breaks[length(tt_breaks)] <- T2_
  }
  if (tt_breaks[length(tt_breaks)] < T2_) {
    tt_breaks <- c(tt_breaks, T2_)
  }
  return(c(tt_, tt_breaks))
}

# create the time grid

time.grid <- function(data.point, coef.t, delta.t,
                      T2., displaygrid = FALSE, N.exp.) {
  tt. <- data.point$ts
  idx.p <- data.point$idx.p
  # time bins
  # find bins break points
  t_b <- breaks_exp(tt., T2., coef_ = coef.t, delta_ = delta.t, N_exp_ = N.exp.)
  
  time.bins <- data.frame(
    t.start = t_b[-length(t_b)],
    t.end = t_b[-1]
  ) %>%
    mutate(t.bin.name = paste0(round(t.start, 3), "-", round(t.end, 3)))
  
  if (nrow(time.bins) - 1 == 0) {
    time.bins$t.ref_layer <- paste0("last-", idx.p)
  } else {
    time.bins$t.ref_layer <- c(1:(nrow(time.bins) - 1), paste0("last-", idx.p))
  }
  time.bins <- cbind(time.bins, data.point, row.names = NULL)
  time.bins
}

# functions to compute the integral

# temporal integral
It_df <- function(param_, time.df) {
  tth <- as.numeric(time.df$ts)
  T1b <- as.numeric(time.df$t.start)
  T2b <- as.numeric(time.df$t.end)
  param_w <- param_[3]
  T.l <- pmax(tth, T1b)
  1 / param_w * (exp(-param_w * (T.l - tth)) - exp(-param_w * (T2b - tth)))
}

# spatial integral
Is_df <- function(param_, space.df) {
  r1 <- as.numeric(space.df$xy_start)
  r2 <- as.numeric(space.df$xy_end)
  param_sig <- param_[4]
  int <- NULL
  for (i in seq_along(r1)) {
    # print(i)
    int[i] <- 2 * pi * param_sig^2 * (exp(-r1[i]^2 / (2 * param_sig^2)) - exp(-r2[i]^2 / (2 * param_sig^2)))
  }
  int
}

# compute temporal integral

compute.grid <- function(param., list.input_) {
  It.vec <- It_df(param_ = param., time.df = list.input_$time.sel)
  It.vec[list.input_$Imapping]
}

# compute the weighted spatial integral

compute.grid_s <- function(param., list.input_) {
  Is.vec <- Is_df(param_ = param., space.df = list.input_$space.sel)
  Is.vec_weights <- Is.vec[list.input_$Imapping_s] * list.input_$weights
  Is.vec_weights
}

# triggering function

gt <- function(th, t, ti, x, xi, y, yi) {
  output <- rep(0, length(ti))
  t.diff <- t - ti
  x.diff <- x - xi
  y.diff <- y - yi
  neg <- t.diff <= 0
  if (sum(!neg) > 0) {
    log.out <- log(th[2]) + log(th[3]) - th[3] * t.diff[!neg] - log(2 * pi * th[4]^2) - (x.diff^2 + y.diff^2) / (2 * th[4]^2)
    output[!neg] <- exp(log.out)
  } else {
    output
  }
  output
}

# Hawkes process conditional intensity

lambda_ <- function(th, t, ti.v, x, xi.v, y, yi.v, k) {
  if (is.null(ti.v) | all(ti.v > t)) {
    a <- th[1]
  } else {
    a <- th[1] + sum(gt(th, t, ti.v, x, xi.v, y, yi.v))
  }
  a
}

# function to calclulate the integral of the triggering

logLambda.i.inla <- function(th.K0, th.w, th.sig, list.input_, link.functions) {
  theta_ <- c(
    0,
    link.functions$K0(th.K0[1]),
    link.functions$w(th.w[1]),
    link.functions$sig(th.sig[1])
  )
  
  # compute the integral efficiently for each bin of each observation
  comp. <- compute.grid(param. = theta_, list.input_ = list.input_)
  comp.s <- compute.grid_s(param. = theta_, list.input_ = list.input_)
  
  out <- log(theta_[2]) + log(theta_[3]) - log(2 * pi * theta_[4]^2) + log(comp. + 1e-10) + log(comp.s + 1e-10)
  out
}

# function to calculate the Hawkes process conditional log-intensity for a set of observations given the history of the process. 

loglambda.inla <- function(th.mu, th.K0, th.w, th.sig, tt, xx, yy, th, link.functions) {
  # if no link.functions are provided
  if (is.null(link.functions)) {
    th.p <- c(th.mu[1], th.K0[1], th.w[1], th.sig[1])
  } else {
    th.p <- c(
      link.functions$mu(th.mu[1]),
      link.functions$K0(th.K0[1]),
      link.functions$w(th.w[1]),
      link.functions$sig(th.sig[1])
    )
  }
  out <- mean(unlist(mclapply(tt, function(x) {
    th_x <- th < x
    id <- which(x==tt)
    log(lambda_(th = th.p, t = x, ti.v = th[th_x], x = xx[id], xi.v=xx[th_x], y=yy[id], yi.v=yy[th_x]))
  }, mc.cores = 5)))
  out
}

# predictor function for the surrogate Poisson model

predictor.fun <- function(th.mu, th.K0, th.w, th.sig,
                          list.input, T1, T2, X1, X2, Y1, Y2,
                          link.functions = NULL) {
  out <- rep(0, list.input$n)
  out[list.input$idx.bkg] <-
    log(link.functions$mu(th.mu[1])) + log(T2 - T1) + log(X2 - X1) + log(Y2 - Y1)
  out[list.input$idx.trig] <- logLambda.i.inla(
    th.K0 = th.K0, th.w = th.w,
    th.sig = th.sig,
    list.input_ = list.input,
    link.functions = link.functions
  )
  
  out[list.input$idx.sl] <- loglambda.inla(
    th.mu = th.mu, th.K0 = th.K0,
    th.w = th.w,
    th.sig = th.sig,
    tt = list.input$sample.s$ts,
    th=list.input$sample.s$ts,
    xx = list.input$sample.s$x,
    yy = list.input$sample.s$y,
    link.functions = link.functions
  )
  out
}

# initial values

th.init <- list(
  th.mu = param[1],
  th.K0 = param[2],
  th.w = param[3],
  th.sig = param[4]
)

param_init <- param

domain_area=1
T1 <- 0
T2 <- 100
X1 <- 0
X2 <- 1
Y1 <- 0
Y2 <- 1
dd.ama <- data
dd.ama <- dd.ama %>% mutate(idx.p = 1:nrow(dd.ama))
dd.ama <- dd.ama %>% rename("ts" = t)

data.bru <- dd.ama %>% dplyr::select(idx.p, x, y, ts)
link.functions <- link.f.be # link functions
data.bru %>% summary

sample.s <- data.bru
bru.opt.list <- list(
  bru_verbose = 4, # type of visual output
  bru_max_iter = 100, # maximum number of iterations
  bru_method = list(max_step = 1.5, rel_tol = 0.1), 
  bru_initial = th.init,
  control.mode = list(x = param_init)
) 
bru.opt <- bru.opt.list
link.functions <- link.f.be 

# settings for temporal bins

coef.t. <- .8 # binning parameter (delta)
delta.t. <- .2 # binning parameter (delta)
N.max. <- 11 # binning parameter (n.max)
delta.t. * ((1 + coef.t.)^(0:N.max.))

df.time <- NULL
for (idx in 1:nrow(data.bru)) {
  # print(idx)
  result <- time.grid(
    data.point = data.bru[idx, ],
    coef.t = coef.t.,
    delta.t = delta.t.,
    T2. = T2,
    N.exp. = N.max.,
  )
  df.time <- rbind(df.time, result)
}

# settings for circular spatial bins

coef_ <- 0.31
delta_ <- 0.02
N_exp_ <- 16
radio <- delta_ * ((1 + coef_)^(0:N_exp_))
radio
xy_start <- c(0, radio[-(N_exp_ + 1)])
xy_end <- radio

medir_t <- proc.time()
df.space <- list()
for (i in 1:dim(data)[1]) {
  xx <- data$x[i]
  yy <- data$y[i]
  df.space[[i]] <- as.data.frame(cbind(
    idx.p = i, xx, yy, xy_start, xy_end,
    s.ref_layer = seq_along(xy_start)
  ))
  df.space[[i]]$within <- apply(df.space[[i]], 1, within_region)
  weight_values <- rep(1, nrow(df.space[[i]])) #Create a vector of NA with the length of df
  if (any(df.space[[i]]$within == "No")) {
    # If there are such rows, assign the values calculated by apply() only to those rows
    weight_values[df.space[[i]]$within == "No"] <- apply(df.space[[i]][df.space[[i]]$within == "No", ], 1, weight_integral2)
  }
  df.space[[i]]$weight <- weight_values
}
df.space <- do.call(rbind, df.space)
t_comp <- (proc.time() - medir_t)
t_comp

# combine the two (temporal and spatial bins)

df.j <- merge(df.time, df.space, by = "idx.p")
df.j <- df.j %>% dplyr::select(t.start, t.end, t.bin.name, t.ref_layer, idx.p, x, y, ts, xy_start, xy_end, weight, s.ref_layer)

# input that will be stored in list.input and needed for the efficient calculation of the expected number of triggered events

t.names <- unique(df.j$t.ref_layer)
time.sel <-
  df.j[vapply(t.names, \(bname) match(TRUE, df.j$t.ref_layer == bname), 0L), , drop = FALSE]
Imapping <- match(df.j$t.ref_layer, t.names)

s.names <- unique(df.j$s.ref_layer)
space.sel <-
  df.j[vapply(s.names, \(bname) match(TRUE, df.j$s.ref_layer == bname), 0L), , drop = FALSE]
Imapping_s <- match(df.j$s.ref_layer, s.names)

list.input <- list(
  df_grid = df.j,
  Imapping = Imapping,
  Imapping_s = Imapping_s,
  time.sel = time.sel,
  space.sel = space.sel,
  weights = df.j$weight
)

df.0 <- data.frame(counts = 0, exposures = 1, part = "background")

df.j$counts <- 0
df.j$exposures <- 1
df.j$part <- "triggered"

df.j <- df.j %>% dplyr::select(-x, -y)

# input that will be stored in list.input and needed for the efficient calculation of the expected number of triggered events

t.names <- unique(df.j$t.ref_layer)
time.sel <-
  df.j[vapply(t.names, \(bname) match(TRUE, df.j$t.ref_layer == bname), 0L), , drop = FALSE]
Imapping <- match(df.j$t.ref_layer, t.names)

s.names <- unique(df.j$s.ref_layer)
space.sel <-
  df.j[vapply(s.names, \(bname) match(TRUE, df.j$s.ref_layer == bname), 0L), , drop = FALSE]
Imapping_s <- match(df.j$s.ref_layer, s.names)

list.input <-
  list(
    Imapping = Imapping,
    time.sel = time.sel,
    Imapping.s = Imapping_s,
    space.sel = space.sel,
    weights = df.j$weight
  )

# third is for the sum of the log intensities

df.s <- data.frame(counts = nrow(sample.s), exposures = 0, part = "SL")

# bind data.frames together

data.input <- bind_rows(df.0, df.s, df.j)

# create list.input

list.input <- list(
  n = nrow(data.input),
  df_grid = df.j,
  Imapping = Imapping,
  Imapping_s = Imapping_s,
  time.sel = time.sel,
  space.sel = space.sel,
  sample.s = sample.s,
  weights = df.j$weight,
  idx.bkg = data.input$part == "background",
  idx.trig = data.input$part == "triggered",
  idx.sl = data.input$part == "SL"
)

# create formula representing the logintensity of the surrogate Poisson counts model

merged.form <- counts ~ predictor.fun(
  th.mu = th.mu, th.K0 = th.K0,
  th.w = th.w,
  th.sig = th.sig,
  list.input = list.input,
  T1 = T1, T2 = T2, X1 = X1,
  X2 = X2, Y1 = Y1, Y2 = Y2,
  link.functions = link.functions
)

# create components representing the parameters in the internal scale

cmp.part <- counts ~ -1 +
  th.mu(1, model = "linear", mean.linear = 0, prec.linear = 1) +
  th.K0(1, model = "linear", mean.linear = 0, prec.linear = 1) +
  th.w(1, model = "linear", mean.linear = 0, prec.linear = 1) +
  th.sig(1, model = "linear", mean.linear = 0, prec.linear = 1)

options(warn = 0)

# Run inlabru

  fit_etas <- bru(
    components = cmp.part,
    like(
      formula = merged.form,
      data = data.input,
      family = "poisson",
      E = exposures
    ),
    options = bru.opt
  )



