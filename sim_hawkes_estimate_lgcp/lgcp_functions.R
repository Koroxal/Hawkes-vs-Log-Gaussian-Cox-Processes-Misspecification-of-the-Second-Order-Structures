
g.diff.single <- function(ps,ghat,useq,model,transform,power,...){
  if(any(ps<=0)) return(NA)
  g.parametric <- exp(CovarianceFct(useq,model=model,sigma=ps[2], phi=ps[1]))
  return(sum((transform(ghat)^power-transform(g.parametric)^power)^2))
}


C.diff.single.im <- function(theta,data,ps,Chat,vseq,spat,model){
  myC <- Cvb(data,spat,100,model,covpars=NULL)
  theo <- sapply(vseq,myC,sigma=sqrt(ps[2]),phi=ps[1],theta=theta)
  return(sum((theo*(Chat[1]/theo[1])-Chat)^2))
}


minimum.contrast.spatiotemporal.adhoc <- function(data, model, method, 
                                                  spatial.dens, temporal.intens, 
                                                  power, transform = NULL, 
                                                  temporal.interval = NULL, 
                                                  verbose = TRUE, ...) 
{
  W <- data$window
  uni.t <- unique(data$t)
  rec <- as.rectangle(W)
  xs <- spatial.dens$X
  xc <- xs[1] + 0.5 * (xs[2] - xs[1]) + 0:(length(xs) - 2) * (xs[2] - xs[1])
  ys <- spatial.dens$Y
  yc <- ys[1] + 0.5 * (ys[2] - ys[1]) + 0:(length(ys) - 2) * (ys[2] - ys[1])
  
  if (verbose) cat("[Univariate spatio-temporal minimum contrast]\n")
  spatial.startvals <- c(min(diff(rec$xrange), diff(rec$yrange))/100, log(data$n)/2)
  
  if (method == "g") {
    if (verbose) 
      cat("Spatial: Time-averaged PCF estimation...\n")
    nonpar <- ginhomAverage(data, spatial.dens, temporal.intens, 
                            suppresswarnings = T, verbose = verbose)
    if (verbose) 
      cat(paste("Spatial: Starting values are (", round(spatial.startvals[1], 
                                                        2), ", ", round(spatial.startvals[2], 2), "); optimising ", 
                model, " correlation function...", sep = ""))
    MINCON.SPATIAL <- optim(par = spatial.startvals, 
                            fn = g.diff.single, ghat = nonpar$iso[-1], useq = nonpar$r[-1], 
                            model = model, transform = transform, power = power, 
                            ...)$par
    parametric <- exp(CovarianceFct(nonpar$r[-1], sigma =MINCON.SPATIAL[2], phi=MINCON.SPATIAL[1],  model = model))
    disc.vec <- sum((transform(nonpar$iso[-1])^power - 
                       transform(parametric)^power)^2) * nonpar$r[2]
    if (verbose) 
      cat("done.\n")
  }
  else if (method == "K") {
    if (verbose) cat("Spatial: Time-averaged K-function estimation...\n")
    nonpar <- KinhomAverage(data, spatial.dens, temporal.intens, suppresswarnings = T, verbose = verbose)
    if (nonpar$iso[2] == 0) {nonpar$iso[2] = nonpar$theo[2]}
    
    if (verbose) 
      cat(paste("Spatial: Starting values are (", round(spatial.startvals[1], 
                                                        2), ", ", round(spatial.startvals[2], 2), "); optimising ", 
                model, " correlation function...", sep = ""))
    MINCON.SPATIAL <- optim(par = spatial.startvals, 
                            fn = K.diff.single, khat = nonpar$iso[-1], useq = nonpar$r[-1], 
                            model = model, transform = transform, power = power, 
                            ...)$par
    parametric <- K.u(nonpar$r[-1], MINCON.SPATIAL[1], MINCON.SPATIAL[2], model = model, ...)
    disc.vec <- sum((transform(nonpar$iso[-1])^power - transform(parametric)^power)^2) * nonpar$r[2]
    if (verbose) 
      cat("done.\n")
  }
  if (verbose) 
    cat("Temporal: Covariance function estimation...")
  uqt <- as.numeric(names(table(as.integer(data$t))))
  tvals <- c()
  for (i in 1:length(uni.t)) tvals[i] <- temporal.intens(uni.t[i])
  autocov <- acf(table(as.integer(data$t)) - tvals, type = "covariance", plot = F)
  v <- 0:5
  blen <- 10
  index <- 1
  vseq <- Chat <- rep(NA, (blen - 1) * max(v) + 1)
  for (i in 1:(length(v) - 1)) {
    xs <- seq(v[i], v[i + 1], length = blen)
    if (i > 1) 
      xs <- xs[-1]
    vseq[index:(index + length(xs) - 1)] <- xs
    Chat[index:(index + length(xs) - 1)] <- approx(v[c(i, i + 1)], autocov$acf[c(i, i + 1)], xout = xs)$y
    index <- index + length(xs)
  }
  if (is.null(temporal.interval)) 
    temporal.interval <- c(0.1, 10)
  if (verbose) 
    cat("done.\n")
  if (verbose) 
    cat(paste("Temporal: Starting interval is (", round(temporal.interval[1], 
                                                        2), ", ", round(temporal.interval[2], 2), "); 
                optimising for exponential dependence...", sep = ""))
  MINCON.TEMPORAL <- optimise(f = C.diff.single.im, interval = temporal.interval, 
                              data = data, ps = MINCON.SPATIAL, Chat = Chat, vseq = vseq, 
                              spat = im(t(spatial.dens$Zm), 
                                        xcol = spatial.dens$X, yrow = spatial.dens$Y), 
                              model = model)$minimum
  if (verbose) 
    cat("done.\n")
  result <- matrix(c(MINCON.SPATIAL, MINCON.TEMPORAL), 
                   1, 3, byrow = T, dimnames = list(NULL, c("scale (spatial)", 
                                                            "variance (spatial)", "scale (temporal)")))
  return(list(estimates = result, discrepancy = matrix(disc.vec, 
                                                       1, 1, dimnames = list(NULL, "Squared discrepancy (spatial)"))))}


gu <- lgcp::CovarianceFct





blockcircbase2 <- function (x, y, sigma, phi, model, additionalparameters, inverse = FALSE) 
{
  M <- length(x)
  N <- length(y)
  xidx <- rep(1:M, N)
  yidx <- rep(1:N, each = M)
  dxidx <- pmin(abs(xidx - xidx[1]), M - abs(xidx - xidx[1]))
  dyidx <- pmin(abs(yidx - yidx[1]), N - abs(yidx - yidx[1]))
  d <- sqrt(((x[2] - x[1]) * dxidx)^2 + ((y[2] - y[1]) * dyidx)^2)
  covbase <- matrix(gu2(d, sigma = sigma, phi = phi, model = model, 
                        additionalparameters = additionalparameters), M, N)
  if (!inverse) {
    return(covbase)
  }
  else {
    return(inversebase(covbase))
  }
}

maternCov <- function (h, sigma, phi, nu) 
{
  if (nu == 1/2) {
    C <- sigma^2 * exp(-phi * abs(h))
  }
  else {
    C <- (sigma^2/(2^(nu - 1) * gamma(nu))) * ((sqrt(2*nu)*abs(h)/phi)^nu) * 
      besselK(sqrt(2*nu)*abs(h)/phi, nu)
  }
  C[h == 0] <- sigma^2
  return(as.matrix(C))
}


g.diff.single2 <- function(ps, ghat, useq, model, transform, power, nu = NULL, ...) {
  # Verificar que los parámetros sean válidos (ps > 0)
  if (any(ps <= 0)) return(NA)
  
  # Seleccionar el modelo de covarianza
  if (model == "exponential") {
    # Covarianza exponencial
    g.parametric <- exp(CovarianceFct(useq,model=model,sigma=ps[2], phi=ps[1]))
  } else if (model == "matern") {
    # Covarianza de Matérn (se requiere nu)
    if (is.null(nu)) stop("El parámetro 'nu' es necesario para el modelo Matérn")
    g.parametric <- exp(maternCov(useq, sigma = ps[2], phi = ps[1], nu = nu))
  } else {
    stop("Modelo desconocido: elige 'exponential' o 'matern'")
  }
  
  # Calcular la suma de las diferencias elevadas a la potencia
  return(sum((transform(ghat)^power - transform(g.parametric)^power)^2))
}


C.diff.single.im2 <- function(theta, data, ps, Chat, vseq, spat, model, additionalparameters) {
  # Asegúrate de que nu es válido antes de llamar a Cvb2
  nu <- additionalparameters$nu
  if (is.na(nu) || nu <= 0) {
    stop("nu debe ser mayor que 0.")
  }
  
  # Llamar a Cvb2
  myC <- Cvb2(data, spat, 100, model, additionalparameters)
  
  # Ajuste: aplicar vseq a la función myC con los otros parámetros
  theo <- sapply(vseq, myC,sigma = sqrt(ps[2]), phi = ps[1], theta = theta, additionalparameters = additionalparameters)
  
  # Retornar el error cuadrático
  return(sum((theo * (Chat[1] / theo[1]) - Chat)^2))
}


Cvb2 <- function (xyt, spatial.intensity, N = 100, spatial.covmodel, 
                  additionalparameters) {
  verifyclass(spatial.intensity, "im")
  sar <- spatialAtRisk(list(X = spatial.intensity$xcol, Y = spatial.intensity$yrow, 
                            Zm = t(spatial.intensity$v)))
  gsx <- length(xvals(sar))
  gsy <- length(yvals(sar))
  xy <- cbind(rep(xvals(sar), gsy), rep(yvals(sar), each = gsx))
  wt <- as.vector(zvals(sar))
  wt[is.na(wt)] <- 0
  sidx <- sample(1:(gsx * gsy), N, prob = wt)
  xy <- xy[sidx, ]
  pd <- as.vector(pairdist(xy))
  cvb <- function(nu, sigma, phi, theta, additionalparameters) {
    return(mean(exp(exp(-nu * theta) * gu2(pd, sigma = sigma, 
                                           phi = phi, model = spatial.covmodel, additionalparameters = additionalparameters)) - 
                  1))
  }
  return(cvb)
}


minimum.contrast.spatiotemporal.adhoc2 <- function(data, model, method, 
                                                   spatial.dens, temporal.intens, 
                                                   power, transform = NULL, 
                                                   temporal.interval = NULL, 
                                                   additionalparameters=NULL,
                                                   verbose = TRUE, ...) 
{
  W <- data$window
  uni.t <- unique(data$t)
  rec <- as.rectangle(W)
  xs <- spatial.dens$X
  xc <- xs[1] + 0.5 * (xs[2] - xs[1]) + 0:(length(xs) - 2) * (xs[2] - xs[1])
  ys <- spatial.dens$Y
  yc <- ys[1] + 0.5 * (ys[2] - ys[1]) + 0:(length(ys) - 2) * (ys[2] - ys[1])
  
  if (verbose) cat("[Univariate spatio-temporal minimum contrast]\n")
  spatial.startvals <- c(min(diff(rec$xrange), diff(rec$yrange))/100, log(data$n)/2)
  
  if (method == "g") {
    if (verbose) 
      cat("Spatial: Time-averaged PCF estimation...\n")
    nonpar <- ginhomAverage(data, spatial.dens, temporal.intens, 
                            suppresswarnings = T, verbose = verbose)
    if (verbose) 
      cat(paste("Spatial: Starting values are (", round(spatial.startvals[1], 
                                                        2), ", ", round(spatial.startvals[2], 2), "); optimising ", 
                model, " correlation function...", sep = ""))
    MINCON.SPATIAL <- optim(par = spatial.startvals, 
                            fn = g.diff.single2, ghat = nonpar$iso[-1], useq = nonpar$r[-1], 
                            model = model, transform = transform, power = power, nu=additionalparameters$nu,
                            ...)$par
    parametric <- exp(maternCov(nonpar$r[-1], sigma =MINCON.SPATIAL[2], phi=MINCON.SPATIAL[1],additionalparameters$nu))
    disc.vec <- sum((transform(nonpar$iso[-1])^power - 
                       transform(parametric)^power)^2) * nonpar$r[2]
    if (verbose) 
      cat("done.\n")
  }
  else if (method == "K") {
    if (verbose) cat("Spatial: Time-averaged K-function estimation...\n")
    nonpar <- KinhomAverage(data, spatial.dens, temporal.intens, suppresswarnings = T, verbose = verbose)
    if (nonpar$iso[2] == 0) {nonpar$iso[2] = nonpar$theo[2]}
    
    if (verbose) 
      cat(paste("Spatial: Starting values are (", round(spatial.startvals[1], 
                                                        2), ", ", round(spatial.startvals[2], 2), "); optimising ", 
                model, " correlation function...", sep = ""))
    MINCON.SPATIAL <- optim(par = spatial.startvals, 
                            fn = K.diff.single, khat = nonpar$iso[-1], useq = nonpar$r[-1], 
                            model = model, transform = transform, power = power, 
                            ...)$par
    parametric <- K.u(nonpar$r[-1], MINCON.SPATIAL[1], MINCON.SPATIAL[2], model = model, ...)
    disc.vec <- sum((transform(nonpar$iso[-1])^power - transform(parametric)^power)^2) * nonpar$r[2]
    if (verbose) 
      cat("done.\n")
  }
  if (verbose) 
    cat("Temporal: Covariance function estimation...")
  uqt <- as.numeric(names(table(as.integer(data$t))))
  tvals <- c()
  for (i in 1:length(uni.t)) tvals[i] <- temporal.intens(uni.t[i])
  autocov <- acf(table(as.integer(data$t)) - tvals, type = "covariance", plot = F)
  v <- 0:5
  blen <- 10
  index <- 1
  vseq <- Chat <- rep(NA, (blen - 1) * max(v) + 1)
  for (i in 1:(length(v) - 1)) {
    xs <- seq(v[i], v[i + 1], length = blen)
    if (i > 1) 
      xs <- xs[-1]
    vseq[index:(index + length(xs) - 1)] <- xs
    Chat[index:(index + length(xs) - 1)] <- approx(v[c(i, i + 1)], autocov$acf[c(i, i + 1)], xout = xs)$y
    index <- index + length(xs)
  }
  if (is.null(temporal.interval)) 
    temporal.interval <- c(0.1, 10)
  if (verbose) 
    cat("done.\n")
  if (verbose) 
    cat(paste("Temporal: Starting interval is (", round(temporal.interval[1], 
                                                        2), ", ", round(temporal.interval[2], 2), "); 
                optimising for exponential dependence...", sep = ""))
  # MINCON.TEMPORAL <- optimise(f = C.diff.single.im2, interval = temporal.interval, 
  #                             data = data, ps = MINCON.SPATIAL, Chat = Chat, vseq = vseq, 
  #                             spat = im(t(spatial.dens$Zm), 
  #                                       xcol = spatial.dens$X, yrow = spatial.dens$Y), 
  #                             model = model, nu=additionalparameters$nu)$minimum
  print("hola")
  MINCON.TEMPORAL <- optimise(
    f = C.diff.single.im2, 
    interval = temporal.interval, 
    data = data, 
    ps = MINCON.SPATIAL, 
    Chat = Chat, 
    vseq = vseq, 
    spat = im(t(spatial.dens$Zm), xcol = spatial.dens$X, yrow = spatial.dens$Y), 
    model = model, 
    additionalparameters = additionalparameters
  )$minimum
  print("adios")
  if (verbose) 
    cat("done.\n")
  result <- matrix(c(MINCON.SPATIAL, MINCON.TEMPORAL), 
                   1, 3, byrow = T, dimnames = list(NULL, c("scale (spatial)", 
                                                            "variance (spatial)", "scale (temporal)")))
  return(list(estimates = result, discrepancy = matrix(disc.vec, 
                                                       1, 1, dimnames = list(NULL, "Squared discrepancy (spatial)"))))}


#rectificación parameteros inicilaes:

minimum.contrast.spatiotemporal.adhoc2 <- function(data, model, method, 
                                                   spatial.dens, temporal.intens, 
                                                   power, transform = NULL, 
                                                   temporal.interval = NULL, 
                                                   additionalparameters=NULL,
                                                   spatial.startvals=NULL,
                                                   verbose = TRUE, ...) 
{
  W <- data$window
  uni.t <- unique(data$t)
  rec <- as.rectangle(W)
  xs <- spatial.dens$X
  xc <- xs[1] + 0.5 * (xs[2] - xs[1]) + 0:(length(xs) - 2) * (xs[2] - xs[1])
  ys <- spatial.dens$Y
  yc <- ys[1] + 0.5 * (ys[2] - ys[1]) + 0:(length(ys) - 2) * (ys[2] - ys[1])
  
  if (verbose) cat("[Univariate spatio-temporal minimum contrast]\n")
  
  if(is.null(spatial.startvals)){
    spatial.startvals <- c(min(diff(rec$xrange), diff(rec$yrange))/100, log(data$n)/2)
  } else {
    spatial.startvals=spatial.startvals
  }
  
  if (method == "g") {
    if (verbose) 
      cat("Spatial: Time-averaged PCF estimation...\n")
    nonpar <- ginhomAverage(data, spatial.dens, temporal.intens, 
                            suppresswarnings = T, verbose = verbose)
    if (verbose) 
      cat(paste("Spatial: Starting values are (", round(spatial.startvals[1], 
                                                        2), ", ", round(spatial.startvals[2], 2), "); optimising ", 
                model, " correlation function...", sep = ""))
    MINCON.SPATIAL <- optim(par = spatial.startvals, 
                            fn = g.diff.single2, ghat = nonpar$iso[-1], useq = nonpar$r[-1], 
                            model = model, transform = transform, power = power, nu=additionalparameters$nu,
                            ...)$par
    parametric <- exp(maternCov(nonpar$r[-1], sigma =MINCON.SPATIAL[2], phi=MINCON.SPATIAL[1],additionalparameters$nu))
    disc.vec <- sum((transform(nonpar$iso[-1])^power - 
                       transform(parametric)^power)^2) * nonpar$r[2]
    if (verbose) 
      cat("done.\n")
  }
  else if (method == "K") {
    if (verbose) cat("Spatial: Time-averaged K-function estimation...\n")
    nonpar <- KinhomAverage(data, spatial.dens, temporal.intens, suppresswarnings = T, verbose = verbose)
    if (nonpar$iso[2] == 0) {nonpar$iso[2] = nonpar$theo[2]}
    
    if (verbose) 
      cat(paste("Spatial: Starting values are (", round(spatial.startvals[1], 
                                                        2), ", ", round(spatial.startvals[2], 2), "); optimising ", 
                model, " correlation function...", sep = ""))
    MINCON.SPATIAL <- optim(par = spatial.startvals, 
                            fn = K.diff.single, khat = nonpar$iso[-1], useq = nonpar$r[-1], 
                            model = model, transform = transform, power = power, 
                            ...)$par
    parametric <- K.u(nonpar$r[-1], MINCON.SPATIAL[1], MINCON.SPATIAL[2], model = model, ...)
    disc.vec <- sum((transform(nonpar$iso[-1])^power - transform(parametric)^power)^2) * nonpar$r[2]
    if (verbose) 
      cat("done.\n")
  }
  if (verbose) 
    cat("Temporal: Covariance function estimation...")
  uqt <- as.numeric(names(table(as.integer(data$t))))
  tvals <- c()
  for (i in 1:length(uni.t)) tvals[i] <- temporal.intens(uni.t[i])
  autocov <- acf(table(as.integer(data$t)) - tvals, type = "covariance", plot = F)
  v <- 0:5
  blen <- 10
  index <- 1
  vseq <- Chat <- rep(NA, (blen - 1) * max(v) + 1)
  for (i in 1:(length(v) - 1)) {
    xs <- seq(v[i], v[i + 1], length = blen)
    if (i > 1) 
      xs <- xs[-1]
    vseq[index:(index + length(xs) - 1)] <- xs
    Chat[index:(index + length(xs) - 1)] <- approx(v[c(i, i + 1)], autocov$acf[c(i, i + 1)], xout = xs)$y
    index <- index + length(xs)
  }
  if (is.null(temporal.interval)) 
    temporal.interval <- c(0.1, 10)
  if (verbose) 
    cat("done.\n")
  if (verbose) 
    cat(paste("Temporal: Starting interval is (", round(temporal.interval[1], 
                                                        2), ", ", round(temporal.interval[2], 2), "); 
                optimising for exponential dependence...", sep = ""))
  # MINCON.TEMPORAL <- optimise(f = C.diff.single.im2, interval = temporal.interval, 
  #                             data = data, ps = MINCON.SPATIAL, Chat = Chat, vseq = vseq, 
  #                             spat = im(t(spatial.dens$Zm), 
  #                                       xcol = spatial.dens$X, yrow = spatial.dens$Y), 
  #                             model = model, nu=additionalparameters$nu)$minimum
  print("hola")
  MINCON.TEMPORAL <- optimise(
    f = C.diff.single.im2, 
    interval = temporal.interval, 
    data = data, 
    ps = MINCON.SPATIAL, 
    Chat = Chat, 
    vseq = vseq, 
    spat = im(t(spatial.dens$Zm), xcol = spatial.dens$X, yrow = spatial.dens$Y), 
    model = model, 
    additionalparameters = additionalparameters
  )$minimum
  print("adios")
  if (verbose) 
    cat("done.\n")
  result <- matrix(c(MINCON.SPATIAL, MINCON.TEMPORAL), 
                   1, 3, byrow = T, dimnames = list(NULL, c("scale (spatial)", 
                                                            "variance (spatial)", "scale (temporal)")))
  return(list(estimates = result, discrepancy = matrix(disc.vec, 
                                                       1, 1, dimnames = list(NULL, "Squared discrepancy (spatial)"))))}

gu2 <- function (u, sigma, phi, model, additionalparameters=NULL) 
{
  if (model == "exponential") {
    return(sigma^2 * exp(-u/phi))
  }
  if (model == "matern"){
    return(maternCov(u,  sigma, phi, additionalparameters$nu))
  }
  else {
    stop("Functionality is temporarily unavailable due to deprecation of RandomFields package. Model 'exponential' is still available")
  }
}



lgcpSim2 <- function (owin = NULL, tlim = as.integer(c(0, 10)), spatial.intensity = NULL, 
                      temporal.intensity = NULL, cellwidth = 0.05, model.parameters = lgcppars(sigma = 2, 
                                                                                               phi = 0.2, theta = 1), spatial.covmodel = "exponential", 
                      covpars = c(), returnintensities = FALSE, progressbar = TRUE, 
                      ext = 2, plot = FALSE, ratepow = 0.25, sleeptime = 0, inclusion = "touching") 
{
  if (!inherits(tlim, "integer")) {
    warning("Converting tlim into integer values, see ?as.integer")
    tlim <- as.integer(tlim)
  }
  tlim <- sort(tlim)
  if (tlim[1] == tlim[2]) {
    stop("Length of time interval given by as.integer(tlim) must be >= 1")
  }
  toffset <- tlim[1]
  maxt <- tlim[2] - toffset
  sigma <- model.parameters$sigma
  phi <- model.parameters$phi
  mu <- model.parameters$mu
  theta <- model.parameters$theta
  if (is.null(owin)) {
    owin <- owin()
  }
  if (is.null(temporal.intensity)) {
    temporal.intensity <- constantInTime(100, tlim)
  }
  else {
    if (!inherits(temporal.intensity, "temporalAtRisk")) {
      temporal.intensity <- temporalAtRisk(temporal.intensity, 
                                           tlim)
    }
    if (!all(tlim == attr(temporal.intensity, "tlim"))) {
      stop("Incompatible temporal.intensity, integer time limits (tlim and temporal.intensity$tlim) do not match")
    }
  }
  ndivs <- diff(tlim)
  tdiff = maxt/ndivs
  times <- tdiff/2 + tdiff * (0:(ndivs - 1))
  mut <- sapply(times + toffset, temporal.intensity)
  if (progressbar) {
    pb <- txtProgressBar(min = 1, max = ndivs, style = 3)
  }
  const0 <- 0.05
  c2 <- -phi * log(const0)
  if (cellwidth > c2/2) {
    warning(paste("cellwidth should be at least", c2/2, "to get accurate results."))
  }
  xyt <- ppp(window = owin)
  ow <- selectObsWindow(xyt, cellwidth)
  xyt <- ow$xyt
  M <- ow$M
  N <- ow$N
  cat(paste("FFT Grid size: [", ext * M, " , ", ext * N, "]\n", 
            sep = ""))
  if (is.null(spatial.intensity)) {
    spatial <- spatialAtRisk(list(X = seq(xyt$window$xrange[1], 
                                          xyt$window$xrange[2], length.out = M), Y = seq(xyt$window$yrange[1], 
                                                                                         xyt$window$yrange[2], length.out = N), Zm = matrix(1/(M * 
                                                                                                                                                 N), M, N)))
  }
  else {
    if (!any(class(spatial.intensity) == "spatialAtRisk")) {
      spatial <- spatialAtRisk(spatial.intensity)
    }
    else {
      spatial <- spatial.intensity
    }
  }
  study.region <- xyt$window
  del1 <- (study.region$xrange[2] - study.region$xrange[1])/M
  del2 <- (study.region$yrange[2] - study.region$yrange[1])/N
  Mext <- ext * M
  Next <- ext * N
  mcens <- study.region$xrange[1] + 0.5 * del1 + (0:(Mext - 
                                                       1)) * del1
  ncens <- study.region$yrange[1] + 0.5 * del2 + (0:(Next - 
                                                       1)) * del2
  xg <- mcens[1:M]
  yg <- ncens[1:N]
  cellarea <- del1 * del2
  if (inclusion == "centroid") {
    cellInside <- inside.owin(x = rep(mcens, Next), y = rep(ncens, 
                                                            each = Mext), w = study.region)
  }
  else if (inclusion == "touching") {
    cellInside <- touchingowin(x = mcens, y = ncens, w = study.region)
  }
  else {
    stop("Invlaid choice for argument 'inclusion'.")
  }
  cellInside <- as.numeric(matrix(as.logical(cellInside), Mext, 
                                  Next)[1:M, 1:N])
  spatialvals <- fftinterpolate(spatial, mcens, ncens, ext = ext)
  spatialvals <- spatialvals[1:M, 1:N]
  spatialvals <- spatialvals * cellInside
  spatialvals <- spatialvals/(cellarea * sum(spatialvals))
  bcb <- blockcircbase2(x = mcens, y = ncens, sigma = sigma, 
                        phi = phi, model = spatial.covmodel, additionalparameters = covpars)
  Qeigs <- eigenfrombase(inversebase(bcb))
  rqe <- sqrt(Qeigs)
  irqe <- 1/rqe
  if (returnintensities) {
    intensities <- array(NA, c(M, N, ndivs))
    truefield <- array(NA, c(M, N, ndivs))
  }
  else {
    intensities <- NULL
    truefield <- NULL
  }
  cases <- NULL
  t <- NULL
  Y <- YfromGamma(matrix(rnorm(Mext * Next), Mext, Next), invrootQeigs = irqe, 
                  mu = mu)[1:M, 1:N]
  rate <- as.vector(mut[1] * spatialvals * cellarea * exp(Y))
  if (returnintensities) {
    intensities[, , 1] <- rate
    truefield[, , 1] <- Y
  }
  cmat <- matrix(rpois(M * N, rate), M, N)
  ncases <- sum(cmat)
  if (ncases > 0) {
    caseidx <- which(cmat > 0)
    caseidx <- unlist(sapply(caseidx, function(x) {
      rep(x, cmat[x])
    }))
    cases <- cbind(rep(xg, length(yg)), rep(yg, each = length(xg)))[caseidx, 
    ] + cbind(runif(ncases, -del1/2, del1/2), runif(ncases, 
                                                    -del2/2, del2/2))
    t <- sort(runif(ncases, times[1] - tdiff/2, times[1] + 
                      tdiff/2))
  }
  if (plot) {
    rate[rate == 0] <- NA
    image.plot(xg, yg, matrix(rate, M, N)^ratepow)
    points(cases, pch = "+", cex = 0.5)
    Sys.sleep(sleeptime)
  }
  for (i in 2:ndivs) {
    Y <- mu * (1 - exp(-theta)) + exp(-theta) * Y + sqrt(1 - 
                                                           exp(-2 * theta)) * YfromGamma(matrix(rnorm(Mext * 
                                                                                                        Next), Mext, Next), invrootQeigs = irqe, mu = 0)[1:M, 
                                                                                                                                                         1:N]
    rate <- as.vector(mut[i] * spatialvals * cellarea * exp(Y))
    cmat <- matrix(rpois(M * N, rate), M, N)
    ncases <- sum(cmat)
    if (ncases > 0) {
      caseidx <- which(cmat > 0)
      caseidx <- unlist(sapply(caseidx, function(x) {
        rep(x, cmat[x])
      }))
      newcases <- cbind(rep(xg, length(yg)), rep(yg, each = length(xg)))[caseidx, 
      ] + cbind(runif(ncases, -del1/2, del1/2), runif(ncases, 
                                                      -del2/2, del2/2))
      cases <- rbind(cases, newcases)
      t <- c(t, sort(runif(ncases, times[i] - tdiff/2, 
                           times[i] + tdiff/2)))
      if (plot) {
        rate[rate == 0] <- NA
        image.plot(xg, yg, matrix(rate, M, N)^ratepow)
        points(newcases, pch = "+", cex = 0.5)
        Sys.sleep(sleeptime)
      }
    }
    if (returnintensities) {
      intensities[, , i] <- rate
      truefield[, , i] <- Y
    }
    if (progressbar) {
      setTxtProgressBar(pb, i)
    }
  }
  if (progressbar) {
    close(pb)
  }
  if (is.null(t)) {
    stop("No data generated for chosen parameters")
  }
  if (!all(inside.owin(cases[, 1], cases[, 2], owin))) {
    remidx <- which(!inside.owin(cases[, 1], cases[, 2], 
                                 owin))
    cases <- cases[-remidx, ]
    t <- t[-remidx]
  }
  xyt <- stppp(ppp(x = cases[, 1], y = cases[, 2], window = owin), 
               t = (t + toffset), tlim = tlim)
  attr(xyt, "rejects") <- NULL
  attr(xyt, "spatialatrisk") <- spatial
  attr(xyt, "temporalfitted") <- mut
  attr(xyt, "cellwidth") <- cellwidth
  attr(xyt, "sigma") <- sigma
  attr(xyt, "phi") <- phi
  attr(xyt, "theta") <- theta
  attr(xyt, "temporalintensity") <- temporal.intensity
  attr(xyt, "temporalfitted") <- mut
  attr(xyt, "spatialcovmodel") <- spatial.covmodel
  attr(xyt, "covpars") <- covpars
  attr(xyt, "ext") <- ext
  attr(xyt, "xvals") <- xg
  attr(xyt, "yvals") <- yg
  attr(xyt, "intensities") <- intensities
  attr(xyt, "truefield") <- truefield
  attr(xyt, "inclusion") <- inclusion
  return(xyt)
}


lgcpPredict2 <- function (xyt, T, laglength, model.parameters = lgcppars(), spatial.covmodel, 
                          covpars = c(), cellwidth = NULL, gridsize = NULL, spatial.intensity, 
                          temporal.intensity, mcmc.control, output.control = setoutput(), 
                          missing.data.areas = NULL, autorotate = FALSE, gradtrunc = Inf, 
                          ext = 2, inclusion = "touching") 
{
  starttime <- Sys.time()
  if (!inherits(T, "integer")) {
    warning("Converting T into integer value, see ?as.integer", 
            immediate. = TRUE)
    T <- as.integer(T)
  }
  if (!inherits(laglength, "integer")) {
    warning("Converting laglength into integer values, see ?as.integer", 
            immediate. = TRUE)
    laglength <- as.integer(laglength)
  }
  if (!inherits(xyt$tlim, "integer")) {
    warning("Converting xyt$tlim into integer values, see ?as.integer", 
            immediate. = TRUE)
    xyt$tlim <- as.integer(xyt$tlim)
  }
  if (!inherits(xyt$t, "integer")) {
    warning("Converting xyt$t into integer values, see ?as.integer", 
            immediate. = TRUE)
    xyt$t <- as.integer(xyt$t)
  }
  if (xyt$window$type == "rectangle") {
    xyt$window <- as.polygonal(xyt$window)
  }
  if (is.null(cellwidth) & is.null(gridsize)) {
    stop("Either cell width OR grid size must be specified")
  }
  if (!is.null(cellwidth) & !is.null(gridsize)) {
    stop("Either cell width OR grid size must be specified")
  }
  if (!all(sapply(gridsize, is.pow2))) {
    stop("All elements of gridsize must be a power of 2")
  }
  if (!is.null(gridsize) & autorotate == TRUE) {
    warning("In order to use autorotate, you must specify a cell width instead ... SETTING autorotate=FALSE.", 
            immediate. = TRUE)
    autorotate <- FALSE
  }
  if (!is.null(gridsize)) {
    approxcw <- diff(xyt$window$xrange)/gridsize[1]
    cwseq <- seq(approxcw/2, 2 * approxcw, length.out = 500)
    cwfun <- function(cw) {
      ow <- selectObsWindow(xyt, cw)
      return(c(ow$M, ow$N))
    }
    gsmat <- t(sapply(cwseq, cwfun))
    tf <- apply(gsmat, 1, function(x) {
      return(all(x == gridsize))
    })
    if (sum(tf) == 0) {
      stop("No sensible observation window found: either change gridsize, or specify cellwidth instead")
    }
    else {
      cellwidth <- cwseq[min(which(tf))]
    }
  }
  if (!is.null(gradtrunc)) {
    if (gradtrunc < 0) {
      stop("gradtrunc must be non-negative")
    }
  }
  if (!inherits(temporal.intensity, "temporalAtRisk")) {
    temporal.intensity <- temporalAtRisk(temporal.intensity, 
                                         tlim = xyt$tlim, xyt = xyt)
  }
  else {
    if (!all(as.integer(xyt$tlim) == attr(temporal.intensity, 
                                          "tlim"))) {
      stop("Incompatible temporal.intensity, integer time limits (xyt$tlim and temporal.intensity$tlim) do not match")
    }
  }
  if (laglength == 0) {
    stop("laglength must be >= 1")
  }
  if (mcmc.control$burnin > mcmc.control$mala.length) {
    stop("Number of burnin iterations must be less than the total number of iterations")
  }
  aggtimes <- T - laglength:0
  nobser <- 0
  for (i in 1:(laglength + 1)) {
    nobser <- nobser + sum(xyt$t == aggtimes[i])
  }
  if (nobser == 0) {
    cat("NOTE: time data should be integer-valued.\n")
    stop("No data in chosen time interval")
  }
  temporalfit <- sapply(aggtimes, temporal.intensity)
  if (any(is.na(temporalfit))) {
    stop("Missing temporal fitted values")
  }
  sigma <- model.parameters$sigma
  phi <- model.parameters$phi
  mu <- model.parameters$mu
  theta <- model.parameters$theta
  tdiff <- c(Inf, diff(aggtimes))
  numt <- length(tdiff)
  bt <- exp(-theta * tdiff)
  if (!autorotate) {
    test <- roteffgain(xyt, cellwidth)
  }
  else {
    test <- roteffgain(xyt, cellwidth)
    if (!test) {
      warning("There is no gain in efficiency by rotating, see ?roteffgain", 
              immediate. = TRUE)
      cat("Not rotating window.\n")
    }
    else {
      rotmat <- getRotation(xyt)$rotation
      xyt <- affine(xyt, mat = rotmat)
    }
  }
  ow <- selectObsWindow(xyt, cellwidth)
  xyt <- ow$xyt
  M <- ow$M
  N <- ow$N
  if (!is.null(missing.data.areas)) {
    if (autorotate) {
      if (test) {
        missing.data.areas <- lapply(missing.data.areas, 
                                     affine, mat = rotmat)
        lapply(1:numt, function(i) {
          missing.data.areas[[i]]$xrange <<- xyt$window$xrange
          missing.data.areas[[i]]$yrange <<- xyt$window$yrange
        })
      }
    }
  }
  tst <- mget("lgcpPredictstapptriggertestvalue", envir = parent.frame(), 
              ifnotfound = FALSE)$lgcpPredictstapptriggertestvalue
  if (tst) {
    del1 <- (xyt$window$xrange[2] - xyt$window$xrange[1])/M
    del2 <- (xyt$window$yrange[2] - xyt$window$yrange[1])/N
    mcens <- xyt$window$xrange[1] + 0.5 * del1 + (0:(M - 
                                                       1)) * del1
    ncens <- xyt$window$yrange[1] + 0.5 * del2 + (0:(N - 
                                                       1)) * del2
    xls <- rep(mcens, N)
    yls <- rep(ncens, each = M)
    spdf <- get("app", envir = parent.frame())$spdf
    olay <- over(SpatialPoints(cbind(xls, yls)), geometry(spdf))
    if (length(table(olay)) != length(spdf)) {
      cat("\n")
      warning("With chosen cell width, will not be able to generate aggregated inference for all regions.", 
              .immediate = TRUE)
      cat("\n")
    }
    olay[is.na(olay)] <- 0
    olay <- matrix(olay, M, N)
  }
  if (M * N >= (256^2)) {
    Sys.sleep(1)
    cat("\n")
    warning("USING LARGE FFT GRID: COMPUTATION MAY BE SLOW ON SOME MACHINES ...", 
            .immediate = TRUE)
    cat("\n")
  }
  cat(paste("FFT Grid size: [", ext * M, " , ", ext * N, "]\n", 
            sep = ""))
  Sys.sleep(1)
  rm(ow)
  if (!any(class(spatial.intensity) == "spatialAtRisk")) {
    spatial <- spatialAtRisk(spatial.intensity)
  }
  else {
    spatial <- spatial.intensity
  }
  if (autorotate) {
    if (test) {
      spatial <- affine(spatial, mat = rotmat)
    }
  }
  study.region <- xyt$window
  del1 <- (study.region$xrange[2] - study.region$xrange[1])/M
  del2 <- (study.region$yrange[2] - study.region$yrange[1])/N
  Mext <- ext * M
  Next <- ext * N
  mcens <- study.region$xrange[1] + 0.5 * del1 + (0:(Mext - 
                                                       1)) * del1
  ncens <- study.region$yrange[1] + 0.5 * del2 + (0:(Next - 
                                                       1)) * del2
  cellarea <- del1 * del2
  if (is.null(missing.data.areas)) {
    if (inclusion == "centroid") {
      cellInside <- inside.owin(x = rep(mcens, Next), y = rep(ncens, 
                                                              each = Mext), w = study.region)
    }
    else if (inclusion == "touching") {
      cellInside <- touchingowin(x = mcens, y = ncens, 
                                 w = study.region)
    }
    else {
      stop("Invlaid choice for argument 'inclusion'.")
    }
    cellInside <- matrix(as.numeric(cellInside), Mext, Next)
    cellInside <- rep(list(cellInside), numt)
  }
  else {
    cellInside <- list()
    for (i in 1:numt) {
      if (inclusion == "centroid") {
        cellInside[[i]] <- inside.owin(x = rep(mcens, 
                                               Next), y = rep(ncens, each = Mext), w = missing.data.areas[[i]])
      }
      else if (inclusion == "touching") {
        cellInside[[i]] <- touchingowin(x = mcens, y = ncens, 
                                        w = missing.data.areas[[i]])
      }
      else {
        stop("Invlaid choice for argument 'inclusion'.")
      }
      cellInside[[i]] <- matrix(as.numeric(cellInside[[i]]), 
                                Mext, Next)
    }
  }
  if (is.null(missing.data.areas)) {
    spatialvals <- fftinterpolate(spatial, mcens, ncens, 
                                  ext = ext)
    spatialvals <- spatialvals * cellInside[[1]]
    spatialvals <- spatialvals/(cellarea * sum(spatialvals))
    spatialvals <- rep(list(spatialvals), numt)
  }
  else {
    if (inclusion == "centroid") {
      cellIns <- inside.owin(x = rep(mcens, Next), y = rep(ncens, 
                                                           each = Mext), w = study.region)
    }
    else if (inclusion == "touching") {
      cellIns <- touchingowin(x = mcens, y = ncens, w = study.region)
    }
    else {
      stop("Invlaid choice for argument 'inclusion'.")
    }
    cellIns <- matrix(as.numeric(cellIns), Mext, Next)
    spatialinterp <- fftinterpolate(spatial, mcens, ncens, 
                                    ext = ext)
    tempinterp <- spatialinterp * cellIns
    NC <- cellarea * sum(tempinterp)
    spatialvals <- list()
    for (i in 1:numt) {
      spatialvals[[i]] <- spatialinterp * cellInside[[i]]
      spatialvals[[i]] <- spatialvals[[i]]/NC
    }
  }
  bcb <- blockcircbase2(x = mcens, y = ncens, sigma = sigma, 
                        phi = phi, model = spatial.covmodel, additionalparameters = covpars)
  Qeigs <- eigenfrombase(inversebase(bcb))
  rqe <- sqrt(Qeigs)
  irqe <- 1/rqe
  mLoop = mcmcLoop(N = mcmc.control$mala.length, burnin = mcmc.control$burnin, 
                   thin = mcmc.control$retain, progressor = mcmcProgressTextBar)
  nsamp <- floor((mLoop$N - mLoop$burnin)/mLoop$thin)
  if (!is.null(output.control$gridfunction) & class(output.control$gridfunction)[1] == 
      "dump2dir") {
    cat("WARNING: disk space required for saving is approximately ", 
        round(nsamp * object.size(array(runif(M * N), dim = c(M, 
                                                              N)))/1024^2, 2), " Mb, ", sep = "")
    if (!output.control$gridfunction$forceSave) {
      m <- menu(c("yes", "no"), title = "continue?")
      if (m == 1) {
        cat("Note: to bypass this menu, set forceSave=TRUE in dump2dir\n")
        Sys.sleep(2)
      }
      else {
        stop("Stopped")
      }
    }
  }
  nis <- list()
  for (i in 1:numt) {
    if (sum(xyt$t == aggtimes[i]) > 0) {
      nis[[i]] <- getCounts(xyt = xyt, subset = (xyt$t == 
                                                   aggtimes[i]), M = M, N = N, ext = ext)
    }
    else {
      nis[[i]] <- matrix(0, ext * M, ext * N)
    }
    ct1 <- sum(nis[[i]])
    nis[[i]] <- nis[[i]] * (spatialvals[[i]] > 0)
    ct2 <- sum(nis[[i]])
    if (ct2 < ct1) {
      warning(paste("Time ", aggtimes[i], ": ", ct1 - ct2, 
                    " data points lost due to discretisation.", sep = ""), 
              immediate. = TRUE)
    }
  }
  if (is.null(gradtrunc)) {
    gradtrunc <- computeGradtruncSpatioTemporal(nsims = 100, 
                                                scale = 1, nis = nis, mu = mu, rootQeigs = rqe, invrootQeigs = irqe, 
                                                spatial = spatialvals, temporal = temporalfit, bt = bt, 
                                                cellarea = cellarea)
  }
  gridfun <- output.control$gridfunction
  if (is.null(gridfun)) {
    gridfun <- nullFunction()
  }
  gridav <- output.control$gridmeans
  if (is.null(gridav)) {
    gridav <- nullAverage()
  }
  lg <- MALAlgcp(mcmcloop = mLoop, inits = mcmc.control$inits, 
                 adaptivescheme = mcmc.control$adaptivescheme, M = M, 
                 N = N, Mext = Mext, Next = Next, sigma = sigma, phi = phi, 
                 theta = theta, mu = mu, nis = nis, cellarea = cellarea, 
                 spatialvals = spatialvals, temporal.fitted = temporalfit, 
                 tdiff = tdiff, rootQeigs = rqe, invrootQeigs = irqe, 
                 cellInside = cellInside, MCMCdiag = mcmc.control$MCMCdiag, 
                 gradtrunc = gradtrunc, gridfun = gridfun, gridav = gridav, 
                 mcens = mcens, ncens = ncens, aggtimes = aggtimes)
  endtime <- Sys.time()
  timetaken <- endtime - starttime
  lg$xyt <- xyt
  lg$M <- M
  lg$N <- N
  lg$aggtimes <- aggtimes
  lg$tdiffs <- tdiff
  lg$vars <- bt
  lg$spatial <- spatial
  lg$temporal <- temporalfit
  lg$grid <- spatialvals
  lg$nis <- lgcpgrid(nis, xvals = mcens[1:M], yvals = ncens[1:N], 
                     zvals = aggtimes)
  lg$mcens <- mcens[1:M]
  lg$ncens <- ncens[1:N]
  lg$cellarea <- diff(mcens[1:2]) * diff(ncens[1:2])
  lg$sigma <- sigma
  lg$phi <- phi
  lg$mu <- mu
  lg$theta <- theta
  lg$mcmcpars <- mcmc.control
  lg$timetaken <- timetaken
  lg$ext <- ext
  lg$cellInside <- lapply(cellInside, function(x) {
    x[1:M, 1:N]
  })
  lg$spatialonly <- FALSE
  lg$inclusion <- inclusion
  if (tst) {
    lg$overlay <- olay
  }
  class(lg) <- c("lgcpPredict", "lgcpobject")
  return(lg)
}

plot_autocovariance <- function(xyt, spatial.intensity = NULL, temporal.intensity = NULL,
                                sigma, phi, theta, spatial.covmodel = "exponential",
                                covpars = c(), N = 100) {
  # Preparar intensidad espacial
  if (inherits(spatial.intensity, "spatialAtRisk")) {
    spatial.intensity <- as.im(spatial.intensity)
  }
  if (is.null(spatial.intensity)) {
    spatial.intensity <- density(xyt)  # Convertir xyt a intensidad espacial
  }
  
  # Preparar intensidad temporal
  if (is.null(temporal.intensity)) {
    temporal.intensity <- muEst(xyt)
  } else if (!inherits(temporal.intensity, "temporalAtRisk")) {
    temporal.intensity <- temporalAtRisk(temporal.intensity, tlim = xyt$tlim, xyt = xyt)
  }
  
  # Verificar compatibilidad temporal
  if (!all(xyt$tlim == attr(temporal.intensity, "tlim"))) {
    stop("Incompatible temporal.intensity: tlim values do not match.")
  }
  
  # Calcular autocovarianza empírica
  uqt <- as.numeric(names(table(as.integer(xyt$t))))
  tfit <- sapply(uqt, temporal.intensity)
  autocov <- acf(table(as.integer(xyt$t)) - tfit, type = "covariance", plot = FALSE)
  
  # Rango de desfases (lags)
  r <- 0:(length(autocov$acf) - 1)
  
  # Graficar la autocovarianza empírica
  plot(r, autocov$acf, type = "l", main = "Autocovariance",
       xlab = "Lag", ylab = "Autocovariance", col = "black")
  abline(h = 0)
  
  # Calcular autocovarianza teórica
  cvb <- Cvb(xyt = xyt, spatial.intensity = spatial.intensity,
             N = N, spatial.covmodel = spatial.covmodel, covpars = covpars)
  rseq <- seq(min(r), max(r), length.out = 100)
  theo <- sapply(rseq, cvb, sigma = sigma, phi = phi, theta = theta)
  
  # Escalar la teórica para que coincida con la empírica en lag 0
  theo <- theo * (autocov$acf[1] / theo[1])
  
  # Graficar la autocovarianza teórica
  lines(rseq, theo, type = "l", col = "red")
  legend("topright", legend = c("Theoretical","Empirical"),
         col = c("black", "red"), lty = 1)
}
