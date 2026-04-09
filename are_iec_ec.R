###############################################################################
# are_iec_ec.R
# Asymptotic Relative Efficiency of the Integrated Empirical Copula Test
# (IECT) versus the Empirical Copula Test (ECT) for Independence
#
# Contents:
#   1. Bessel-zero computation
#   2. Eigenvalues and eigenfunctions (EC and IEC)
#   3. Copula score derivatives at independence
#   4. Projection coefficients (double numerical integration)
#   5. Gil-Pelaez inversion for null distributions
#   6. Local power slopes and ARE
#   7. Example: 9 copula families
#
# Platform note:
#   This script is single-threaded and runs on any platform (Unix, macOS,
#   Windows) without modification.
###############################################################################

rm(list = ls())

cat("============================================================\n")
cat("  ARE(IEC, EC): Genest EC + Proposition 3.5 IEC\n")
cat("  All critical values computed from their own truncated null laws\n")
cat("============================================================\n\n")

###############################################################################
# 1. Bessel zeros of J_{1/4}
###############################################################################

find_bessel_zeros <- function(nu, n_zeros, x_max = 250, dx = 0.0025) {
  x    <- seq(dx, x_max, by = dx)
  vals <- besselJ(x, nu)
  idx  <- which(diff(sign(vals)) != 0)
  out  <- numeric(0)

  for (i in idx) {
    if (length(out) >= n_zeros) break
    root <- uniroot(
      function(t) besselJ(t, nu),
      interval = c(x[i], x[i + 1]),
      tol = 1e-14
    )$root
    out <- c(out, root)
  }
  out[seq_len(min(length(out), n_zeros))]
}

###############################################################################
# 2. Eigenvalues and eigenfunctions
###############################################################################

ec_lambda_mat <- function(K) {
  outer(1:K, 1:K, function(k, l) 1 / (pi^4 * k^2 * l^2))
}

iec_lambda_mat <- function(a) {
  outer(a, a, function(ak, al) 1 / (16 * ak^2 * al^2))
}

ec_phi <- function(u, k) {
  sqrt(2) * sin(k * pi * u)
}

iec_phi <- function(u, a_k) {
  2 * u^(3 / 2) * besselJ(a_k * u^2, 1 / 4) / besselJ(a_k, 5 / 4)
}

###############################################################################
# 3. Copula score derivatives  d/dtheta C_theta(u,v)  evaluated at theta_0
###############################################################################

get_cdot <- function(family) {
  eps <- 1e-15

  if (family == "gaussian") {
    return(function(u, v) dnorm(qnorm(u)) * dnorm(qnorm(v)))
  }
  if (family == "fgm") {
    return(function(u, v) u * (1 - u) * v * (1 - v))
  }
  if (family == "frank") {
    return(function(u, v) 0.5 * u * v * (1 - u) * (1 - v))
  }
  if (family == "clayton") {
    return(function(u, v) u * v * log(pmax(u, eps)) * log(pmax(v, eps)))
  }
  if (family == "clayton180") {
    return(function(u, v) {
      (1 - u) * (1 - v) * log(pmax(1 - u, eps)) * log(pmax(1 - v, eps))
    })
  }
  if (family == "gumbel") {
    return(function(u, v) {
      a <- -log(pmax(u, eps))
      b <- -log(pmax(v, eps))
      s <- a + b
      u * v * (a * log(s / a) + b * log(s / b))
    })
  }
  if (family == "plackett") {
    return(function(u, v) u * (1 - u) * v * (1 - v))
  }
  if (family == "amh") {
    return(function(u, v) u * (1 - u) * v * (1 - v))
  }
  if (family == "joe") {
    return(function(u, v) {
      a  <- pmax(1 - u, eps)
      b  <- pmax(1 - v, eps)
      g  <- a + b - a * b
      gp <- a * log(a) + b * log(b) - a * b * (log(a) + log(b))
      -g * (-log(g) + gp / g)
    })
  }

  stop("Unknown family: ", family)
}

###############################################################################
# 4. Projection coefficients (double numerical integration)
###############################################################################

compute_I_ec <- function(K, cdot_fun, rel.tol = 1e-9) {
  I <- matrix(0, K, K)
  for (k in 1:K) {
    for (l in 1:K) {
      outer_int <- function(u) {
        sapply(u, function(uu) {
          inner_val <- integrate(
            function(v) sin(l * pi * v) * cdot_fun(uu, v),
            lower = 1e-10, upper = 1 - 1e-10,
            subdivisions = 400L, rel.tol = rel.tol
          )$value
          sin(k * pi * uu) * inner_val
        })
      }
      val <- integrate(
        outer_int,
        lower = 1e-10, upper = 1 - 1e-10,
        subdivisions = 400L, rel.tol = rel.tol
      )$value
      I[k, l] <- 2 * k * l * pi^2 * val
    }
  }
  I
}

compute_I_iec <- function(a, cdot_fun, rel.tol = 1e-9) {
  K    <- length(a)
  Ibar <- matrix(0, K, K)
  for (k in 1:K) {
    for (l in 1:K) {
      outer_int <- function(u) {
        sapply(u, function(uu) {
          inner_val <- integrate(
            function(v) iec_phi(v, a[l]) * uu * v * cdot_fun(uu, v),
            lower = 1e-10, upper = 1 - 1e-10,
            subdivisions = 400L, rel.tol = rel.tol
          )$value
          iec_phi(uu, a[k]) * inner_val
        })
      }
      val <- integrate(
        outer_int,
        lower = 1e-10, upper = 1 - 1e-10,
        subdivisions = 400L, rel.tol = rel.tol
      )$value
      Ibar[k, l] <- 4 * a[k] * a[l] * val
    }
  }
  Ibar
}

###############################################################################
# 5. Gil-Pelaez inversion for weighted chi-squared null distributions
#    Q = sum_j lambda_j * Z_j^2
###############################################################################

gp_prepare <- function(lambda_vec, t_max = 1500, dt = 1 / 2500) {
  t       <- seq(dt, t_max, by = dt)
  log_mod <- rep(0, length(t))
  phase   <- rep(0, length(t))

  for (lam in lambda_vec) {
    log_mod <- log_mod - 0.25 * log1p(4 * lam^2 * t^2)
    phase   <- phase + 0.5 * atan2(2 * lam * t, 1)
  }

  list(t = t, dt = dt, log_mod = log_mod, phase = phase)
}

## Survival function P(Q >= x)
gp_surv <- function(x, gp) {
  with(gp, {
    0.5 + sum(exp(log_mod) * sin(phase - t * x) / t) * dt / pi
  })
}

## Density of Q + lambda0 * chi^2_2
hk_density <- function(x, lambda0, gp) {
  with(gp, {
    extra_log_mod <- -0.5 * log1p(4 * lambda0^2 * t^2)
    extra_phase   <- atan2(2 * lambda0 * t, 1)
    sum(exp(log_mod + extra_log_mod) *
          cos(phase + extra_phase - t * x)) * dt / pi
  })
}

## Critical value at level alpha
find_critical_value <- function(alpha, gp, lower = 0, upper = 0.2,
                                tol = 1e-10) {
  f <- function(x) gp_surv(x, gp) - alpha
  while (f(upper) > 0) upper <- 2 * upper
  uniroot(f, interval = c(lower, upper), tol = tol)$root
}

###############################################################################
# 6. Local power slopes and ARE
###############################################################################

compute_slope_ec <- function(K, I_ec, alpha = 0.05,
                             t_max = 1500, dt = 1 / 2500) {
  lam_mat <- ec_lambda_mat(K)
  lam_vec <- as.vector(lam_mat)

  gp      <- gp_prepare(lam_vec, t_max = t_max, dt = dt)
  c_alpha <- find_critical_value(alpha, gp)
  I_vec   <- as.vector(I_ec)

  slope <- sum(vapply(seq_along(lam_vec), function(j) {
    lam_vec[j] * I_vec[j]^2 * hk_density(c_alpha, lam_vec[j], gp)
  }, numeric(1)))

  list(slope = slope, critical = c_alpha, lambda = lam_vec)
}

compute_slope_iec <- function(a, I_iec, alpha = 0.05,
                              t_max = 1500, dt = 1 / 2500) {
  lam_mat <- iec_lambda_mat(a)
  lam_vec <- as.vector(lam_mat)

  gp      <- gp_prepare(lam_vec, t_max = t_max, dt = dt)
  c_alpha <- find_critical_value(alpha, gp)
  I_vec   <- as.vector(I_iec)

  slope <- sum(vapply(seq_along(lam_vec), function(j) {
    lam_vec[j] * I_vec[j]^2 * hk_density(c_alpha, lam_vec[j], gp)
  }, numeric(1)))

  list(slope = slope, critical = c_alpha, lambda = lam_vec)
}

###############################################################################
# 7. Main wrapper: compute ARE for a given copula family
###############################################################################

compute_are <- function(family, K = 10, alpha = 0.05,
                        t_max = 1500, dt = 1 / 2500) {
  cat("\n------------------------------------------------------------\n")
  cat("Family:", family, "\n")
  cat("------------------------------------------------------------\n")

  cdot_fun <- get_cdot(family)
  a        <- find_bessel_zeros(1 / 4, K)

  I_ec  <- compute_I_ec(K, cdot_fun)
  I_iec <- compute_I_iec(a, cdot_fun)

  ec_out  <- compute_slope_ec(K, I_ec, alpha = alpha,
                              t_max = t_max, dt = dt)
  iec_out <- compute_slope_iec(a, I_iec, alpha = alpha,
                               t_max = t_max, dt = dt)

  are <- iec_out$slope / ec_out$slope

  cat(sprintf("EC critical value:  %.10f\n", ec_out$critical))
  cat(sprintf("IEC critical value: %.10f\n", iec_out$critical))
  cat(sprintf("EC local slope:     %.10e\n", ec_out$slope))
  cat(sprintf("IEC local slope:    %.10e\n", iec_out$slope))
  cat(sprintf("ARE(IEC, EC):       %.10f\n", are))

  invisible(list(
    family       = family,
    ec_critical  = ec_out$critical,
    iec_critical = iec_out$critical,
    ec_slope     = ec_out$slope,
    iec_slope    = iec_out$slope,
    are          = are,
    I_ec         = I_ec,
    I_iec        = I_iec
  ))
}

###############################################################################
# 8. Run: 9 copula families
###############################################################################

res_joe        <- compute_are("joe")
res_gumbel     <- compute_are("gumbel")
res_gaussian   <- compute_are("gaussian")
res_clayton    <- compute_are("clayton")
res_clayton180 <- compute_are("clayton180")
res_frank      <- compute_are("frank")
res_fgm        <- compute_are("fgm")
res_amh        <- compute_are("amh")
res_plackett   <- compute_are("plackett")
