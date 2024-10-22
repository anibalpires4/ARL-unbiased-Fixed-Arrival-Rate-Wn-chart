################
### M/Ek/1,  ###
################

library(spc)

cdf.VU1_m <- Vectorize(function(x, rho, k) {
  lambda <- sqrt(rho)
  mu <- 1 / sqrt(rho)
  
  if (x <= 0) {
    result <- exp(lambda * x) * (k * mu / (k * mu + lambda))^k
  } else {
    result <- pgamma(x, k, k * mu) + exp(lambda * x) * (k * mu / (k * mu + lambda))^k * (1 - pgamma(x, k, k * mu + lambda))
  }
  result
})

pdf.VU1_m <- function(x, rho, k) {
  lambda <- sqrt(rho)
  mu <- 1 / sqrt(rho)
  
  as.numeric(x <= 0) * lambda * cdf.VU1_m(x, rho, k) + 
    as.numeric(x > 0) * (dgamma(x, k, k * mu) + 
                           lambda * exp(lambda * x) * (k * mu / (k * mu + lambda))^k * (1 - pgamma(x, k, k * mu + lambda)) - 
                           exp(lambda * x) * (k * mu / (k * mu + lambda))^k * dgamma(x, k, k * mu + lambda))
}

qf.VU1_m <- function(y, rho, k) {
  lambda <- sqrt(rho)
  mu <- 1 / sqrt(rho)
  
  ff <- function(x) y - cdf.VU1_m(x, rho, k)
  x <- uniroot(ff, c(-50, 50))$root
  x
}



# Chebyshev polynomials
Tn <- function(z, n) {
  if ( n==0 ) result <- 1
  if ( n==1 ) result <- z
  if ( n==2 ) result <- 2*z^2 - 1
  if ( n==3 ) result <- 4*z^3 - 3*z
  if ( n==4 ) result <- 8*z^4 - 8*z^2 + 1
  if ( n==5 ) result <- 16*z^5 - 20*z^3 + 5*z
  if ( n>5 )  result <- cos( n*acos(z) ) 
  result
}


# derivatives of Chebyshev polynomials
dTn <- function(z, n) {
  if ( n==0 ) result <- 0
  if ( n==1 ) result <- 1
  if ( n==2 ) result <- 4*z
  if ( n==3 ) result <- 12*z^2 - 3
  if ( n==4 ) result <- 32*z^3 - 16*z
  if ( n==5 ) result <- 80*z^4 - 60*z^2 + 5
  if ( n>5  ) result <- n * ( Tn(z,n-1) - z*Tn(z,n) ) / (1-z^2)
  result
}


coll.arl_m <- function(h, gamma, rho, k, z0 = 0, r = 40, qm = 30) {
  # Substitute lambda and mu
  lambda <- sqrt(rho)
  mu <- 1 / sqrt(rho)
  
  GQ <- quadrature.nodes.weights(qm, type = "GL", x1 = -1, x2 = 1)
  z <- GQ$nodes
  w <- GQ$weights
  zch <- h * (1 + cos(pi * (2 * (r:1) - 1) / 2 / r)) / 2
  A <- matrix(NA, nrow = r, ncol = r)
  
  for (i in 1:r) {
    zi <- zch[i]
    xm1 <- zi / 2
    h_ <- qf.VU1_m(1 - 1e-9, rho, k) + zi
    if (h_ > h) h_ <- h
    xm2 <- (zi + h_) / 2
    xw1 <- zi / 2
    xw2 <- (h_ - zi) / 2
    w1 <- xw1 * w
    w2 <- xw2 * w
    z1 <- xm1 + xw1 * z
    z2 <- xm2 + xw2 * z
    
    for (j in 1:r) {
      integral <- sum(w1 * pdf.VU1_m(z1 - zi, rho, k) * Tn(2 * z1 / h - 1, j - 1)) +
        sum(w2 * pdf.VU1_m(z2 - zi, rho, k) * Tn(2 * z2 / h - 1, j - 1))
      A[i, j] <- Tn(2 * zi / h - 1, j - 1) - (1 - gamma) * cdf.VU1_m(-zi, rho, k) * Tn(-1, j - 1) - integral
    }
  }
  
  one <- rep(1, r) 
  g <- solve(A, one)
  result <- 0
  for (i in 1:r) {
    result <- result + g[i] * Tn(2 * z0 / h - 1, i - 1)
  }
  result
}



coll2.arl_m <- function(h, gamma, rho, k, z0 = 0, r = 40, qm = 30) {
  # Substitute lambda and mu
  lambda <- sqrt(rho)
  mu <- 1 / sqrt(rho)
  
  GQ <- quadrature.nodes.weights(qm, type = "GL", x1 = -1, x2 = 1)
  z <- GQ$nodes
  w <- GQ$weights
  zch <- h * (1 + cos(pi * (2 * (r:1) - 1) / 2 / r)) / 2
  A <- matrix(NA, nrow = r, ncol = r)
  
  for (i in 1:r) {
    zi <- zch[i]
    xm1 <- zi / 2
    xm2 <- (zi + h) / 2
    xw1 <- zi / 2
    xw2 <- (h - zi) / 2
    w1 <- xw1 * w
    w2 <- xw2 * w
    z1 <- xm1 + xw1 * z
    z2 <- xm2 + xw2 * z
    
    for (j in 1:r) {
      integral <- cdf.VU1_m(h - zi, rho, k) - 
        cdf.VU1_m(-zi, rho, k) * (-1)^(j - 1) - 
        sum(w1 * cdf.VU1_m(z1 - zi, rho, k) * dTn(2 * z1 / h - 1, j - 1) * 2 / h) - 
        sum(w2 * cdf.VU1_m(z2 - zi, rho, k) * dTn(2 * z2 / h - 1, j - 1) * 2 / h)
      A[i, j] <- Tn(2 * zi / h - 1, j - 1) - (1 - gamma) * cdf.VU1_m(-zi, rho, k) * (-1)^(j - 1) - integral
    }
  }
  
  one <- rep(1, r) 
  g <- solve(A, one)
  result <- 0
  for (i in 1:r) {
    result <- result + g[i] * Tn(2 * z0 / h - 1, i - 1)
  }
  result
}



get.h_m <- function(L0, rho, k, z0 = 0, r = 40, qm = 30, OUTPUT = FALSE) {
  # Substitute lambda and mu
  lambda <- sqrt(rho)
  mu <- 1 / sqrt(rho)
  
  h1 <- 1
  L1 <- coll.arl_m(h1, 0, rho, k, z0 = z0, r = r, qm = qm)
  
  if (L1 > L0) {
    while (L1 > L0) {
      h2 <- h1
      L2 <- L1
      h1 <- h1 * 0.8
      L1 <- coll.arl_m(h1, 0, rho, k, z0 = z0, r = r, qm = qm)
      if (OUTPUT) cat(paste("1:\t", h1, "\t", L1, "\n"))
    }
  } else {
    while (L1 <= L0) {
      h2 <- h1
      L2 <- L1
      h1 <- h1 * 1.2
      L1 <- coll.arl_m(h1, 0, rho, k, z0 = z0, r = r, qm = qm)
      if (OUTPUT) cat(paste("1:\t", h1, "\t", L1, "\n"))
    }
  }
  
  h.error <- 1
  L.error <- 1
  while (abs(h.error) > 1e-9 & abs(L.error) > 1e-9) {
    h3 <- h1 + (L0 - L1) / (L2 - L1) * (h2 - h1)
    L3 <- coll.arl_m(h3, 0, rho, k, z0 = z0, r = r, qm = qm)
    if (OUTPUT) cat(paste("2:\t", h3, "\t", L3, "\n"))
    h1 <- h2
    h2 <- h3
    L1 <- L2
    L2 <- L3
    h.error <- h2 - h1
    L.error <- L2 - L1
  }
  h3
}



get.gamma_m <- function(h, L0, rho, k, z0 = 0, r = 40, qm = 30, OUTPUT = FALSE) {
  # Substitute lambda and mu
  lambda <- sqrt(rho)
  mu <- 1 / sqrt(rho)
  
  L0.full_m <- coll.arl_m(h, 0, rho, k, z0 = z0, r = r, qm = qm)
  if (OUTPUT) cat(paste("0:\t", 0, "\t", L0.full_m, "\n"))
  if (L0.full_m <= L0) stop("h too small or L0 too large")
  
  g1 <- 0.05
  L1 <- coll.arl_m(h, g1, rho, k, z0 = z0, r = r, qm = qm)
  
  if (L1 > L0) {
    while (L1 > L0) {
      g2 <- g1
      L2 <- L1
      g1 <- g1 * 1.2
      L1 <- coll.arl_m(h, g1, rho, k, z0 = z0, r = r, qm = qm)
      if (OUTPUT) cat(paste("1:\t", g1, "\t", L1, "\n"))
    }
  } else {
    while (L1 <= L0) {
      g2 <- g1
      L2 <- L1
      g1 <- g1 * 0.8
      L1 <- coll.arl_m(h, g1, rho, k, z0 = z0, r = r, qm = qm)
      if (OUTPUT) cat(paste("1:\t", g1, "\t", L1, "\n"))
    }
  }
  
  g.error <- 1
  L.error <- 1
  while (abs(g.error) > 1e-9 & abs(L.error) > 1e-9) {
    g3 <- g1 + (L0 - L1) / (L2 - L1) * (g2 - g1)
    L3 <- coll.arl_m(h, g3, rho, k, z0 = z0, r = r, qm = qm)
    if (OUTPUT) cat(paste("2:\t", g3, "\t", L3, "\n"))
    g1 <- g2
    g2 <- g3
    L1 <- L2
    L2 <- L3
    g.error <- g2 - g1
    L.error <- L2 - L1
  }
  g3
}



get.gammaAh_m <- function(L0, rho, k, param_type = "rho", z0 = 0, r = 40, qm = 30, eps = 1e-6, OUTPUT = FALSE) {
  # Substitute lambda and mu
  lambda <- sqrt(rho)
  mu <- 1 / sqrt(rho)
  
  # Determine the parameter bounds based on the param_type
  if (param_type == "rho") {
    rho_m <- rho - eps
    rho_p <- rho + eps
  } else {
    stop("Invalid param_type. Use 'rho'.")
  }
  
  # Initial calculations
  h1 <- get.h_m(2 * L0, rho, k, z0 = z0, r = r, qm = qm)
  g1 <- get.gamma_m(h1, L0, rho, k, z0 = z0, r = r, qm = qm)
  
  # Compute ARL values based on param_type
  ARLm <- coll.arl_m(h1, g1, rho_m, k, z0 = z0, r = r, qm = qm)
  ARLp <- coll.arl_m(h1, g1, rho_p, k, z0 = z0, r = r, qm = qm)
  
  dratio1 <- (ARLm - ARLp) / (2 * eps)
  if (OUTPUT) cat(paste("0\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
  
  # Iterative process to adjust h1
  if (dratio1 > 0) {
    while (dratio1 > 0) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 1.2
      g1 <- get.gamma_m(h1, L0, rho, k, z0 = z0, r = r, qm = qm)
      
      ARLm <- coll.arl_m(h1, g1, rho_m, k, z0 = z0, r = r, qm = qm)
      ARLp <- coll.arl_m(h1, g1, rho_p, k, z0 = z0, r = r, qm = qm)
      
      dratio1 <- (ARLm - ARLp) / (2 * eps)
      if (OUTPUT) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  } else {
    while (dratio1 <= 0) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 0.8
      g1 <- get.gamma_m(h1, L0, rho, k, z0 = z0, r = r, qm = qm)
      
      ARLm <- coll.arl_m(h1, g1, rho_m, k, z0 = z0, r = r, qm = qm)
      ARLp <- coll.arl_m(h1, g1, rho_p, k, z0 = z0, r = r, qm = qm)
      
      dratio1 <- (ARLm - ARLp) / (2 * eps)
      if (OUTPUT) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  }
  
  # Refine h1 using a secant-like method
  h.error <- 1
  dr.error <- 1
  while (abs(h.error) > 1e-10 & abs(dr.error) > 1e-8) {
    h3 <- h1 - dratio1 / (dratio2 - dratio1) * (h2 - h1)
    g3 <- get.gamma_m(h3, L0, rho, k, z0 = z0, r = r, qm = qm)
    
    ARLm <- coll.arl_m(h3, g3, rho_m, k, z0 = z0, r = r, qm = qm)
    ARLp <- coll.arl_m(h3, g3, rho_p, k, z0 = z0, r = r, qm = qm)
    
    dratio3 <- (ARLm - ARLp) / (2 * eps)
    if (OUTPUT) cat(paste("2\th3 =", h3, "\tg3 =", g3, "\tdratio1 =", dratio3, "\n"))
    h1 <- h2
    h2 <- h3
    dratio1 <- dratio2
    dratio2 <- dratio3
    h.error <- h2 - h1
    dr.error <- dratio2 - dratio1
  }
  
  data.frame(h = h3, gamma = g3)
}



##############
### Ek/M/1 ###
##############


cdf.VU1 <- Vectorize(function(x, rho, k) {
  # Substitute lambda and mu
  lambda <- sqrt(rho)
  mu <- 1 / sqrt(rho)
  
  if (x <= 0) {
    result <- 1 - pgamma(-x, k, k * lambda) - exp(-mu * x) * (k * lambda / (k * lambda + mu))^k * (1 - pgamma(-x, k, k * lambda + mu))
  } else {
    result <- 1 - exp(-mu * x) * (k * lambda / (k * lambda + mu))^k
  }
  result
})  

pdf.VU1 <- function(x, rho, k) {
  # Substitute lambda and mu
  lambda <- sqrt(rho)
  mu <- 1 / sqrt(rho)
  
  as.numeric(x <= 0) * (dgamma(-x, k, k * lambda) + 
                          mu * exp(-mu * x) * (k * lambda / (k * lambda + mu))^k * (1 - pgamma(-x, k, k * lambda + mu)) - 
                          exp(-mu * x) * (k * lambda / (k * lambda + mu))^k * dgamma(-x, k, k * lambda + mu)) + 
    as.numeric(x > 0) * mu * exp(-mu * x) * (k * lambda / (k * lambda + mu))^k
}

qf.VU1 <- function(y, rho, k) {
  # Substitute lambda and mu
  lambda <- sqrt(rho)
  mu <- 1 / sqrt(rho)
  
  ff <- function(x) y - cdf.VU1(x, rho, k)
  x <- uniroot(ff, c(-50, 50))$root
  x
}



# Chebyshev polynomials
Tn <- function(z, n) {
  if ( n==0 ) result <- 1
  if ( n==1 ) result <- z
  if ( n==2 ) result <- 2*z^2 - 1
  if ( n==3 ) result <- 4*z^3 - 3*z
  if ( n==4 ) result <- 8*z^4 - 8*z^2 + 1
  if ( n==5 ) result <- 16*z^5 - 20*z^3 + 5*z
  if ( n>5 )  result <- cos( n*acos(z) ) 
  result
}


# derivatives of Chebyshev polynomials
dTn <- function(z, n) {
  if ( n==0 ) result <- 0
  if ( n==1 ) result <- 1
  if ( n==2 ) result <- 4*z
  if ( n==3 ) result <- 12*z^2 - 3
  if ( n==4 ) result <- 32*z^3 - 16*z
  if ( n==5 ) result <- 80*z^4 - 60*z^2 + 5
  if ( n>5  ) result <- n * ( Tn(z,n-1) - z*Tn(z,n) ) / (1-z^2)
  result
}


# solve ARL integral equation with collocation
coll.arl <- function(h, gamma, rho, k, z0 = 0, r = 40, qm = 30) {
  # Substitute lambda and mu
  lambda <- sqrt(rho)
  mu <- 1 / sqrt(rho)
  
  GQ <- quadrature.nodes.weights(qm, type = "GL", x1 = -1, x2 = 1)
  z <- GQ$nodes
  w <- GQ$weights
  zch <- h * (1 + cos(pi * (2 * (r:1) - 1) / 2 / r)) / 2
  A <- matrix(NA, nrow = r, ncol = r)
  
  for (i in 1:r) {
    zi <- zch[i]
    xm1 <- zi / 2
    h_ <- qf.VU1(1 - 1e-9, rho, k) + zi
    if (h_ > h) h_ <- h
    xm2 <- (zi + h_) / 2
    xw1 <- zi / 2
    xw2 <- (h_ - zi) / 2
    w1 <- xw1 * w
    w2 <- xw2 * w
    z1 <- xm1 + xw1 * z
    z2 <- xm2 + xw2 * z
    
    for (j in 1:r) {
      integral <- sum(w1 * pdf.VU1(z1 - zi, rho, k) * Tn(2 * z1 / h - 1, j - 1)) + 
        sum(w2 * pdf.VU1(z2 - zi, rho, k) * Tn(2 * z2 / h - 1, j - 1))
      A[i, j] <- Tn(2 * zi / h - 1, j - 1) - (1 - gamma) * cdf.VU1(-zi, rho, k) * Tn(-1, j - 1) - integral
    }
  }
  
  one <- rep(1, r)
  g <- solve(A, one)
  result <- 0
  for (i in 1:r) {
    result <- result + g[i] * Tn(2 * z0 / h - 1, i - 1)
  }
  result
}



# solve ARL integral equation with collocation
coll2.arl <- function(h, gamma, rho, k, z0 = 0, r = 40, qm = 30) {
  # Substitute lambda and mu
  lambda <- sqrt(rho)
  mu <- 1 / sqrt(rho)
  
  GQ <- quadrature.nodes.weights(qm, type = "GL", x1 = -1, x2 = 1)
  z <- GQ$nodes
  w <- GQ$weights
  zch <- h * (1 + cos(pi * (2 * (r:1) - 1) / 2 / r)) / 2
  A <- matrix(NA, nrow = r, ncol = r)
  
  for (i in 1:r) {
    zi <- zch[i]
    xm1 <- zi / 2
    xm2 <- (zi + h) / 2
    xw1 <- zi / 2
    xw2 <- (h - zi) / 2
    w1 <- xw1 * w
    w2 <- xw2 * w
    z1 <- xm1 + xw1 * z
    z2 <- xm2 + xw2 * z
    
    for (j in 1:r) {
      integral <- cdf.VU1(h - zi, rho, k) - 
        cdf.VU1(-zi, rho, k) * (-1)^(j - 1) - 
        sum(w1 * cdf.VU1(z1 - zi, rho, k) * dTn(2 * z1 / h - 1, j - 1) * 2 / h) - 
        sum(w2 * cdf.VU1(z2 - zi, rho, k) * dTn(2 * z2 / h - 1, j - 1) * 2 / h)
      A[i, j] <- Tn(2 * zi / h - 1, j - 1) - (1 - gamma) * cdf.VU1(-zi, rho, k) * (-1)^(j - 1) - integral
    }
  }
  
  one <- rep(1, r)
  g <- solve(A, one)
  result <- 0
  for (i in 1:r) {
    result <- result + g[i] * Tn(2 * z0 / h - 1, i - 1)
  }
  result
}



# calibrate (determine threshold h) without randomization (aka gamma=0), collocation based
get.h <- function(L0, rho, k, z0 = 0, r = 40, qm = 30, OUTPUT = FALSE) {
  # Substitute lambda and mu
  lambda <- sqrt(rho)
  mu <- 1 / sqrt(rho)
  
  h1 <- 1
  L1 <- coll.arl(h1, 0, rho, k, z0 = z0, r = r, qm = qm)
  
  if (L1 > L0) {
    while (L1 > L0) {
      h2 <- h1
      L2 <- L1
      h1 <- h1 * 0.8
      L1 <- coll.arl(h1, 0, rho, k, z0 = z0, r = r, qm = qm)
      if (OUTPUT) cat(paste("1:\t", h1, "\t", L1, "\n"))
    }
  } else {
    while (L1 <= L0) {
      h2 <- h1
      L2 <- L1
      h1 <- h1 * 1.2
      L1 <- coll.arl(h1, 0, rho, k, z0 = z0, r = r, qm = qm)
      if (OUTPUT) cat(paste("1:\t", h1, "\t", L1, "\n"))
    }
  }
  
  h.error <- 1
  L.error <- 1
  while (abs(h.error) > 1e-9 & abs(L.error) > 1e-9) {
    h3 <- h1 + (L0 - L1) / (L2 - L1) * (h2 - h1)
    L3 <- coll.arl(h3, 0, rho, k, z0 = z0, r = r, qm = qm)
    if (OUTPUT) cat(paste("2:\t", h3, "\t", L3, "\n"))
    h1 <- h2
    h2 <- h3
    L1 <- L2
    L2 <- L3
    h.error <- h2 - h1
    L.error <- L2 - L1
  }
  h3
}




# search for gamma so that for fixed threshold the ARL becomes pre-defined L0, collocation based
get.gamma <- function(h, L0, rho, k, z0 = 0, r = 40, qm = 30, OUTPUT = FALSE) {
  # Substitute lambda and mu
  lambda <- sqrt(rho)
  mu <- 1 / sqrt(rho)
  
  L0.full <- coll.arl(h, 0, rho, k, z0 = z0, r = r, qm = qm)
  if (OUTPUT) cat(paste("0:\t", 0, "\t", L0.full, "\n"))
  if (L0.full <= L0) stop("h too small or L0 too large")
  
  g1 <- 0.05
  L1 <- coll.arl(h, g1, rho, k, z0 = z0, r = r, qm = qm)
  
  if (L1 > L0) {
    while (L1 > L0) {
      g2 <- g1
      L2 <- L1
      g1 <- g1 * 1.2
      L1 <- coll.arl(h, g1, rho, k, z0 = z0, r = r, qm = qm)
      if (OUTPUT) cat(paste("1:\t", g1, "\t", L1, "\n"))
    }
  } else {
    while (L1 <= L0) {
      g2 <- g1
      L2 <- L1
      g1 <- g1 * 0.8
      L1 <- coll.arl(h, g1, rho, k, z0 = z0, r = r, qm = qm)
      if (OUTPUT) cat(paste("1:\t", g1, "\t", L1, "\n"))
    }
  }
  
  g.error <- 1
  L.error <- 1
  while (abs(g.error) > 1e-9 & abs(L.error) > 1e-9) {
    g3 <- g1 + (L0 - L1) / (L2 - L1) * (g2 - g1)
    L3 <- coll.arl(h, g3, rho, k, z0 = z0, r = r, qm = qm)
    if (OUTPUT) cat(paste("2:\t", g3, "\t", L3, "\n"))
    g1 <- g2
    g2 <- g3
    L1 <- L2
    L2 <- L3
    g.error <- g2 - g1
    L.error <- L2 - L1
  }
  g3
}



# Combined function to calculate (h, gamma) for ARL = L0 and ARL maximization based on lambda or mu
get.gammaAh <- function(L0, rho, k, param_type = "rho", z0 = 0, r = 40, qm = 30, eps = 1e-6, OUTPUT = FALSE) {
  # Substitute lambda and mu
  lambda <- sqrt(rho)
  mu <- 1 / sqrt(rho)
  
  # Determine the parameter bounds based on the param_type
  if (param_type == "rho") {
    rho_m <- rho - eps
    rho_p <- rho + eps
  } else {
    stop("Invalid param_type. Use 'rho'.")
  }
  
  # Initial calculations
  h1 <- get.h(2 * L0, rho, k, z0 = z0, r = r, qm = qm)
  g1 <- get.gamma(h1, L0, rho, k, z0 = z0, r = r, qm = qm)
  
  # Compute ARL values based on rho
  ARLm <- coll.arl(h1, g1, rho_m, k, z0 = z0, r = r, qm = qm)
  ARLp <- coll.arl(h1, g1, rho_p, k, z0 = z0, r = r, qm = qm)
  
  dratio1 <- (ARLm - ARLp) / (2 * eps)
  if (OUTPUT) cat(paste("0\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
  
  # Iterative process to adjust h1
  if (dratio1 > 0) {
    while (dratio1 > 0) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 1.2
      g1 <- get.gamma(h1, L0, rho, k, z0 = z0, r = r, qm = qm)
      
      ARLm <- coll.arl(h1, g1, rho_m, k, z0 = z0, r = r, qm = qm)
      ARLp <- coll.arl(h1, g1, rho_p, k, z0 = z0, r = r, qm = qm)
      
      dratio1 <- (ARLm - ARLp) / (2 * eps)
      if (OUTPUT) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  } else {
    while (dratio1 <= 0) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 0.8
      g1 <- get.gamma(h1, L0, rho, k, z0 = z0, r = r, qm = qm)
      
      ARLm <- coll.arl(h1, g1, rho_m, k, z0 = z0, r = r, qm = qm)
      ARLp <- coll.arl(h1, g1, rho_p, k, z0 = z0, r = r, qm = qm)
      
      dratio1 <- (ARLm - ARLp) / (2 * eps)
      if (OUTPUT) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  }
  
  # Refine h1 using a secant-like method
  h.error <- 1
  dr.error <- 1
  while (abs(h.error) > 1e-10 & abs(dr.error) > 1e-8) {
    h3 <- h1 - dratio1 / (dratio2 - dratio1) * (h2 - h1)
    g3 <- get.gamma(h3, L0, rho, k, z0 = z0, r = r, qm = qm)
    
    ARLm <- coll.arl(h3, g3, rho_m, k, z0 = z0, r = r, qm = qm)
    ARLp <- coll.arl(h3, g3, rho_p, k, z0 = z0, r = r, qm = qm)
    
    dratio3 <- (ARLm - ARLp) / (2 * eps)
    if (OUTPUT) cat(paste("2\th3 =", h3, "\tg3 =", g3, "\tdratio1 =", dratio3, "\n"))
    h1 <- h2
    h2 <- h3
    dratio1 <- dratio2
    dratio2 <- dratio3
    h.error <- h2 - h1
    dr.error <- dratio2 - dratio1
  }
  
  data.frame(h = h3, gamma = g3)
}







#########################################
#########################################
#########################################

# Load necessary libraries
library(shiny)
library(ggplot2)
library(shinycssloaders)


# Define UI
ui <- fluidPage(
  titlePanel("ARL-unbiased Wn chart - Proportional changes in arrival and service rates"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons("distribution", "Select queue:",
                   choices = c("M/Ek/1", "Ek/M/1")),
      conditionalPanel(
        condition = "input.distribution == 'M/Ek/1'",
        sliderInput("rho", HTML("Target Traffic Intensity (&rho;0):"), min = 0.1, max = 0.9, value = 0.1, step = 0.1),
        sliderInput("L0", "Target ARL:", min = 100, max = 600, value = 200, step = 1),
        sliderInput("k", "Erlang Parameter (k):", min = 1, max = 10, value = 1, step = 1)
      ),
      conditionalPanel(
        condition = "input.distribution == 'Ek/M/1'",
        sliderInput("rho", HTML("Target Traffic Intensity (&rho;0):"), min = 0.1, max = 0.9, value = 0.1, step = 0.1),
        sliderInput("L0", "Target ARL:", min = 100, max = 600, value = 200, step = 1),
        sliderInput("k", "Erlang Parameter (k):", min = 1, max = 10, value = 1, step = 1)
      ),
      downloadButton("downloadPlot", "Download Plot")
    ),
    
    mainPanel(
      withSpinner(plotOutput("ARLPlot", height = "600px"), type = 8, color = "#0D6EFD")
    )
  )
)


server <- function(input, output) {
  
  # Store the calculated plot data for re-use in downloadHandler
  rv <- reactiveValues()
  
  output$ARLPlot <- renderPlot({
    # Create a progress object
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    
    distribution <- input$distribution
    
    if (distribution == "M/Ek/1") {
      progress$set(message = "Calculating", value = 0)
      progress$inc(0.2, message = "Initializing variables")
      
      rho0 <- input$rho
      L0 <- input$L0
      k <- input$k
      
      progress$inc(0.3, message = "Calculating limits")
      calibration_results_m <- get.gammaAh_m(L0, rho0, k)
      h <- calibration_results_m$h
      gamma <- calibration_results_m$gamma
      
      progress$inc(0.3, message = "Generating plot data")
      rho_seq <- seq(0.01, 0.99, by = 0.01)
      ARL_values <- sapply(rho_seq, function(rho) coll.arl_m(h, gamma, rho, k))
      
      # Cache results in reactiveValues for downloading
      rv$rho_seq <- rho_seq
      rv$ARL_values <- ARL_values
      rv$h <- h
      rv$gamma <- gamma
      
      progress$inc(0.2, message = "Rendering plot")
      legend_pos <- if (rho0 > 0.5) "topleft" else "topright"
      
      plot(rho_seq, ARL_values, type = 'l', col = 'black', xlab = "ρ", ylab = "ARL")
      abline(v=rho0, h=L0, lty=4, col="grey")
      legend(legend_pos, legend = c(paste("UCL =", round(h, 6)), paste("γL =", round(gamma, 6))), cex = 1.2, inset = c(0.025, 0.1))
      
    } else if (distribution == "Ek/M/1") {
      progress$set(message = "Calculating", value = 0)
      progress$inc(0.2, message = "Initializing variables")
      
      rho0 <- input$rho
      L0 <- input$L0
      k <- input$k
      
      progress$inc(0.3, message = "Calculating limits")
      calibration_results <- get.gammaAh(L0, rho0, k)
      h <- calibration_results$h
      gamma <- calibration_results$gamma
      
      progress$inc(0.3, message = "Generating plot data")
      rho_seq <- seq(0.01, 0.99, by = 0.01)
      ARL_values <- sapply(rho_seq, function(rho) coll.arl(h, gamma, rho, k))
      
      # Cache results in reactiveValues for downloading
      rv$rho_seq <- rho_seq
      rv$ARL_values <- ARL_values
      rv$h <- h
      rv$gamma <- gamma
      
      progress$inc(0.2, message = "Rendering plot")
      legend_pos <- if (rho0 > 0.5) "topleft" else "topright"
      
      plot(rho_seq, ARL_values, type = 'l', col = 'black', xlab = "ρ", ylab = "ARL")
      abline(v=rho0, h=L0, lty=4, col="grey")
      legend(legend_pos, legend = c(paste("UCL =", round(h, 6)), paste("γL =", round(gamma, 6))), cex = 1.2, inset = c(0.025, 0.1))
    }
  })
  
  # Add downloadHandler for downloading the plot
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("ARLPlot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      png(file)
      # Reuse the cached plot data from reactiveValues
      rho_seq <- rv$rho_seq
      ARL_values <- rv$ARL_values
      h <- rv$h
      gamma <- rv$gamma
      
      legend_pos <- if (input$rho > 0.5) "topleft" else "topright"
      
      plot(rho_seq, ARL_values, type = 'l', col = 'black', xlab = "ρ", ylab = "ARL")
      abline(v=input$rho, h=input$L0, lty=4, col="grey")
      legend(legend_pos, legend = c(paste("UCL =", round(h, 6)), paste("γL =", round(gamma, 6))), cex = 1.2, inset = c(0.025, 0.1))
      
      dev.off()
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)

