# functions to generate Hawkes data (gaussian triggering for space and exponential triggering for time)
# using the acceptance-rejection method

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

# triggering function

nu <- function(kernel, t, s, his_t, his_s) {
  
  delta_t <- t - his_t
  if(is.vector(his_s)){
    delta_s <- s - his_s
    delta_x <- delta_s[1]
    delta_y <- delta_s[2]
  }else{
    delta_s <- matrix(, nrow = nrow(his_s), ncol = ncol(his_s))
    delta_s[,1] <- s[1]-his_s[,1]
    delta_s[,2] <- s[2]-his_s[,2]
    delta_x <- delta_s[, 1]
    delta_y <- delta_s[, 2]
  }
  return(
    kernel$C * kernel$w * exp(- kernel$w * delta_t) *
      (1 / (2*pi*kernel$beta^2)) * exp(- (1 / (2*kernel$beta^2)) * (delta_x^2 + delta_y^2))
  )
}

# lambda (conditional intensity)

value <- function(lam, t, his_t, s, his_s) {
  if (length(his_t) > 0) {
    val <- lam$mu + sum(nu(kernel=lam$kernel, t=t, s=s, his_t=his_t, his_s=his_s))
  } else {
    val <- lam$mu
  }
  return(val)
}

# upper bound of lambda

upper_bound <- function(lam) {
  return(lam$maximum)
}

# lebesgue measure function

lebesgue_measure <- function(subspace){
  sub_lebesgue_ms <- c()
  for(i in 1:length(subspace[,1])){
    sub_lebesgue_ms[i] = subspace[i, 2] - subspace[i, 1]
  }
  for (j in 1: length(sub_lebesgue_ms)){
    if (sub_lebesgue_ms[j] == 0){sub_lebesgue_ms[j]=NA}
  }
  producto <- prod(sub_lebesgue_ms, na.rm = TRUE)
  return(producto)
}

# homogeneous poisson sampling function

homogeneous_poisson_sampling <- function(lam, T, S) {
  S_ <- rbind(T, S)
  n <- lebesgue_measure(S_)
  N <- rpois(1, lambda = upper_bound(lam) * n)
  
  points_prov <- c()
  points_apend <- c()
  for(i in 1:length(S_[,1])){
    points_prov <-  runif(n=N, min= S_[i,1], max=S_[i,2])
    points_apend <- rbind(points_apend, points_prov)
  }
  points <- t(points_apend)
  
  points <- points[order(points[, 1]), ]
  return(points)
}

# thinning function 

inhomogeneous_poisson_thinning <- function(lam, homo_points, verbose) {
  retained_points <- matrix(numeric(0), ncol = ncol(homo_points))
  if (verbose) {
    cat(sprintf("[%s] generate %s samples from homogeneous poisson point process\n", 
                Sys.time(), dim(homo_points)))
  }
  
  for (i in 1:nrow(homo_points)) {
    t <- homo_points[i, 1]
    s <- homo_points[i, -1]
    his_t <- retained_points[, 1]
    his_s <- retained_points[, -1]
    
    lam_value <- value(lam=lam, t=t, his_t=his_t, s=s, his_s=his_s)
    lam_bar <- upper_bound(lam)
    D <- runif(1)
    
    if (lam_value > lam_bar) {
      cat(sprintf("intensity %f is greater than upper bound %f.\n", lam_value, lam_bar))
      return(NULL)
    }
    
    if (lam_value >= D * lam_bar) {
      retained_points <- rbind(retained_points, homo_points[i, , drop = FALSE])
    }
    
    if (verbose && i != 1 && i %% (nrow(homo_points) %/% 10) == 0) {
      cat(sprintf("[%s] %d raw samples have been checked. %d samples have been retained.\n", 
                  Sys.time(), i, nrow(retained_points)))
    }
  }
  
  if (verbose) {
    cat(sprintf("[%s] thinning samples %s based on %s.\n", 
                Sys.time(), dim(retained_points), deparse(substitute(lam))))
  }
  
  return(retained_points)
}

# function to generate data

generate <- function(lam, T, S, batch_size, min_n_points, verbose) {
  points_list <- list()
  sizes <- numeric()
  max_len <- 0
  b <- 0
  
  while (b < batch_size) {
    homo_points <- homogeneous_poisson_sampling(lam, T, S)
    points <- inhomogeneous_poisson_thinning(lam, homo_points, verbose)
    
    if (is.null(points) || nrow(points) < min_n_points) {
      next
    }
    
    max_len <- ifelse(max_len < nrow(points), nrow(points), max_len)
    points_list[[b + 1]] <- points
    sizes <- c(sizes, nrow(points))
    cat(sprintf("[%s] %d-th sequence is generated.\n", Sys.time(), b+1))
    b <- b + 1
  }
  
  # fit the data into a tensor
  data <- array(0, dim = c(batch_size, max_len, ncol(points)))
  for (b in seq_along(points_list)) {
    data[b, 1:nrow(points_list[[b]]), ] <- points_list[[b]]
  }
  
  return(list(data = data, sizes = sizes))
}

# Thinning function with history included

inhomogeneous_poisson_thinning_history <- function(lam, homo_points, verbose, retained_points) {
  #retained_points <- matrix(numeric(0), ncol = ncol(homo_points))
  if (verbose) {
    cat(sprintf("[%s] generate %s samples from homogeneous poisson point process\n", 
                Sys.time(), dim(homo_points)))
  }
  
  for (i in 1:nrow(homo_points)) {
    t <- homo_points[i, 1]
    s <- homo_points[i, -1]
    his_t <- retained_points[, 1]
    his_s <- retained_points[, -1]
    
    lam_value <- value(lam=lam, t=t, his_t=his_t, s=s, his_s=his_s)
    lam_bar <- upper_bound(lam)
    D <- runif(1)
    
    if (lam_value > lam_bar) {
      cat(sprintf("intensity %f is greater than upper bound %f.\n", lam_value, lam_bar))
      return(NULL)
    }
    
    if (lam_value >= D * lam_bar) {
      retained_points <- rbind(retained_points, homo_points[i, , drop = FALSE])
    }
    
    if (verbose && i != 1 && i %% (nrow(homo_points) %/% 10) == 0) {
      cat(sprintf("[%s] %d raw samples have been checked. %d samples have been retained.\n", 
                  Sys.time(), i, nrow(retained_points)))
    }
  }
  
  if (verbose) {
    cat(sprintf("[%s] thinning samples %s based on %s.\n", 
                Sys.time(), dim(retained_points), deparse(substitute(lam))))
  }
  
  return(retained_points)
}

# function to generate data with history included

generate_history <- function(lam, T, S, batch_size, min_n_points, verbose, retained_points) {
  points_list <- list()
  sizes <- numeric()
  max_len <- 0
  b <- 0
  
  while (b < batch_size) {
    homo_points <- homogeneous_poisson_sampling(lam, T, S)
    points <- inhomogeneous_poisson_thinning_history(lam, homo_points, verbose, retained_points)
    
    if (is.null(points) || nrow(points) < min_n_points) {
      next
    }
    
    max_len <- ifelse(max_len < nrow(points), nrow(points), max_len)
    points_list[[b + 1]] <- points
    sizes <- c(sizes, nrow(points))
    cat(sprintf("[%s] %d-th sequence is generated.\n", Sys.time(), b+1))
    b <- b + 1
  }
  
  # fit the data into a tensor
  data <- array(0, dim = c(batch_size, max_len, ncol(points)))
  for (b in seq_along(points_list)) {
    data[b, 1:nrow(points_list[[b]]), ] <- points_list[[b]]
  }
  
  return(list(data = data, sizes = sizes))
}
