###############################################################################
# Simulation_IEC.R
# Monte Carlo Comparison of the Integrated Empirical Copula Test (IECT)
# and the Empirical Copula Test (ECT) for Independence
#
# Contents:
#   1. DGP definitions (8 dependence structures)
#   2. Test statistics (IECT: Rnbar, ECT: Bn)
#   3. Empirical level and power (Tables)
#   4. Joe copula power curves
#   5. Null p-value histograms
#   6. Real-data application: Uranium exploration (Li vs Ti)
#
# Platform note:
#   This script uses mclapply() from the {parallel} package, which relies on
#   POSIX fork() and runs natively on Unix / macOS only.
#
#   Windows users: replace every mclapply_progress() call with an equivalent
#   using parLapply() over a PSOCK cluster, e.g.
#
#     library(parallel)
#     cl <- makeCluster(detectCores() - 4)
#     clusterExport(cl, c("compute_stats", ...))    # export needed objects
#     clusterEvalQ(cl, library(copula))              # load packages on workers
#     result <- parLapply(cl, seq_len(N), function(i) { ... })
#     stopCluster(cl)
#
#   Alternatively, set  ncores <- 1  to run sequentially on any platform.
###############################################################################

rm(list = ls())

library(copula)
library(parallel)

set.seed(20240817)

###############################################################################
# 0. Utility: mclapply with progress bar and ETA
###############################################################################

mclapply_progress <- function(X, FUN, ...,
                              mc.cores   = getOption("mc.cores", 1L),
                              chunk_size = mc.cores * 5L,
                              label      = "") {
  n      <- length(X)
  start  <- proc.time()["elapsed"]
  chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  out    <- vector("list", n)

  for (k in seq_along(chunks)) {
    idx      <- chunks[[k]]
    out[idx] <- mclapply(X[idx], FUN, ..., mc.cores = mc.cores)
    done     <- max(idx)
    elapsed  <- proc.time()["elapsed"] - start
    eta      <- (n - done) / (done / elapsed)
    cat(sprintf("\r    %-28s %4d/%d (%3.0f%%)  elapsed: %5.1fs  ETA: %5.1fs",
                label, done, n, done / n * 100, elapsed, eta))
    flush.console()
  }

  elapsed <- proc.time()["elapsed"] - start
  cat(sprintf("\r    %-28s %4d/%d (100%%)  elapsed: %5.1fs  done.           \n",
              label, n, n, elapsed))
  out
}

###############################################################################
# 1. Global settings
###############################################################################

n_grid          <- c(50, 100, 200, 500)
num_simulations <- 5000   # null critical-value and empirical-level replicates
num_power       <- 1000   # power replicates per DGP
ncores          <- detectCores() - 4

###############################################################################
# 2. Data-generating processes (DGPs)
#    All return raw (n x 2) matrices; pobs() is applied at the test stage.
###############################################################################

gen_type1 <- function(n) {                        # Four Clouds
  cx <- sample(c(-1, 1), n, replace = TRUE)
  cy <- sample(c(-1, 1), n, replace = TRUE)
  cbind(cx + rnorm(n) / 3, cy + rnorm(n) / 3)
}

gen_type2 <- function(n) {                        # W-shape
  x <- seq(-1, 1, length.out = n)
  cbind(x + runif(n) / 3, 4 * ((x^2 - 0.5)^2 + runif(n) / 500))
}

gen_type3 <- function(n) {                        # Diamond (rotated square)
  x <- runif(n, -1, 1)
  y <- runif(n, -1, 1)
  theta <- -pi / 4
  R <- rbind(c(cos(theta), -sin(theta)),
             c(sin(theta),  cos(theta)))
  cbind(x, y) %*% R
}

gen_type4 <- function(n) {                        # Parabola
  x <- seq(-1, 1, length.out = n)
  cbind(x, (x^2 + runif(n)) / 2)
}

gen_type5 <- function(n) {                        # Two Parabolas
  x <- seq(-1, 1, length.out = n)
  y <- (x^2 + runif(n) / 2) * sample(c(-1, 1), n, replace = TRUE)
  cbind(x, y)
}

gen_type6 <- function(n) {                        # Circle
  x <- seq(-1, 1, length.out = n)
  cbind(sin(x * pi) + rnorm(n) / 8, cos(x * pi) + rnorm(n) / 8)
}

gen_type7 <- function(n) {                        # Hyperplane
  dx <- runif(n)
  cbind(dx, dx + runif(n))
}

gen_type8 <- function(n) {                        # Student-t copula (tau=0, df=2)
  rCopula(n, tCopula(param = 0, dim = 2, df = 2))
}

dgp_funs  <- list(gen_type1, gen_type2, gen_type3, gen_type4,
                  gen_type5, gen_type6, gen_type7, gen_type8)
dgp_names <- c("Four Clouds", "W-shape", "Diamond", "Parabola",
               "Two Parabolas", "Circle", "Hyperplane", "Student-t copula")

###############################################################################
# 3. Illustrative scatter plots (n = 1000)
###############################################################################

plot_titles <- c("(a) Four Clouds",   "(b) W-shape",
                 "(c) Diamond",       "(d) Parabola",
                 "(e) Two Parabolas", "(f) Circle",
                 "(g) Hyperplane",    "(h) Student-t copula")

par(mfrow = c(2, 4), mar = c(2, 2, 2, 2))
for (i in seq_along(dgp_funs)) {
  samp <- dgp_funs[[i]](1000)
  plot(samp[, 1], samp[, 2],
       xlab = "", ylab = "", axes = FALSE,
       main = plot_titles[i],
       col = "steelblue", pch = 16, cex = 0.4)
}
par(mfrow = c(1, 1))

###############################################################################
# 4. Test statistics: IECT (Rnbar) and ECT (Bn)
#    Input: (n x 2) matrix of pseudo-observations in [0,1]^2.
###############################################################################

compute_stats <- function(uvdat) {
  n    <- nrow(uvdat)
  U    <- uvdat[, 1]
  V    <- uvdat[, 2]
  Ubar <- 1 - U
  Vbar <- 1 - V

  ## --- S2, S3, S4 via O(n^2) sweep (IECT) ---
  ord        <- order(U)
  Ubar_s     <- Ubar[ord]
  Vbar_s     <- Vbar[ord]
  Ubar_breaks <- c(1, Ubar_s, 0)

  S2 <- S3 <- S4 <- 0
  vbar_set <- numeric(0)

  for (m in seq_len(n)) {
    pos      <- findInterval(Vbar_s[m], vbar_set)
    vbar_set <- append(vbar_set, Vbar_s[m], after = pos)
    du       <- Ubar_breaks[m + 1] - Ubar_breaks[m + 2]
    widths   <- diff(c(0, vbar_set, 1))
    Mvals    <- m:0
    S2 <- S2 + du * sum(widths * Mvals^2)
    S3 <- S3 + du * sum(widths * Mvals^3)
    S4 <- S4 + du * sum(widths * Mvals^4)
  }

  ## --- T1, T2 for IECT (O(n^2)) ---
  T1   <- sum((1 - U^3) * (1 - V^3))
  Umat <- outer(U, U, pmax)
  Vmat <- outer(V, V, pmax)
  T2   <- sum((1 - Umat^3) * (1 - Vmat^3))

  ## --- IECT statistic (Proposition 2.1) ---
  Rnbar <- n * ((S4 + 2 * S3 + S2) / (4 * n^4) -
                  (T2 + T1) / (18 * n^2) + 1 / 100)

  ## --- ECT statistic Bn (Genest, Quessy & Remillard 2006, Section 1) ---
  R  <- round(U * (n + 1))
  S  <- round(V * (n + 1))
  DR <- outer(R, R, function(s, t)
    (2 * n + 1) / (6 * n) +
      s * (s - 1) / (2 * n * (n + 1)) +
      t * (t - 1) / (2 * n * (n + 1)) -
      pmax(s, t) / (n + 1))
  DS <- outer(S, S, function(s, t)
    (2 * n + 1) / (6 * n) +
      s * (s - 1) / (2 * n * (n + 1)) +
      t * (t - 1) / (2 * n * (n + 1)) -
      pmax(s, t) / (n + 1))
  Bn <- sum(DR * DS) / n

  c(Rnbar = Rnbar, Bn = Bn)
}

###############################################################################
# 5. Main simulation: empirical level and power
###############################################################################

all_results <- list()
cv_list     <- list()

for (n in n_grid) {
  cat(sprintf("\n====== n = %d ======\n", n))

  ## --- Critical values under H0 ---
  cat(sprintf("  Critical values (%d simulations)...\n", num_simulations))
  h0_stats <- mclapply_progress(
    seq_len(num_simulations),
    function(i) compute_stats(pobs(cbind(runif(n), runif(n)))),
    mc.cores = ncores, label = "critical values"
  )
  h0_mat  <- do.call(rbind, h0_stats)
  cvRnbar <- quantile(h0_mat[, "Rnbar"], 0.95)
  cvBn    <- quantile(h0_mat[, "Bn"],    0.95)
  cv_list[[as.character(n)]] <- list(cvRnbar = cvRnbar, cvBn = cvBn)
  cat(sprintf("  Critical value -- IECT: %.4f | ECT: %.4f\n", cvRnbar, cvBn))

  ## --- Empirical level ---
  cat(sprintf("  Empirical level (%d simulations)...\n", num_simulations))
  level_stats <- mclapply_progress(
    seq_len(num_simulations),
    function(i) compute_stats(pobs(cbind(runif(n), runif(n)))),
    mc.cores = ncores, label = "empirical level"
  )
  level_mat <- do.call(rbind, level_stats)

  rows <- list()
  rows[[1]] <- data.frame(
    DGP        = "Level (H0)",
    n          = n,
    Power_IECT = mean(level_mat[, "Rnbar"] >= cvRnbar),
    Power_ECT  = mean(level_mat[, "Bn"]    >= cvBn)
  )

  ## --- Power for each DGP ---
  cat(sprintf("  Power for 8 DGPs (%d simulations each)...\n", num_power))
  for (i in seq_along(dgp_funs)) {
    pow_stats <- mclapply_progress(
      seq_len(num_power),
      function(j) compute_stats(pobs(dgp_funs[[i]](n))),
      mc.cores = ncores, label = dgp_names[i]
    )
    pow_mat <- do.call(rbind, pow_stats)
    rows[[i + 1]] <- data.frame(
      DGP        = dgp_names[i],
      n          = n,
      Power_IECT = mean(pow_mat[, "Rnbar"] >= cvRnbar),
      Power_ECT  = mean(pow_mat[, "Bn"]    >= cvBn)
    )
  }

  all_results[[as.character(n)]] <- do.call(rbind, rows)
}

## --- Final table ---
final_results <- do.call(rbind, all_results)
rownames(final_results) <- NULL

cat("\n=== Empirical level and power (alpha = 0.05) ===\n")
print(final_results, digits = 3, row.names = FALSE)

###############################################################################
# 6. Power curve: Joe copula (theta = 1 -> independence, theta > 1 -> dep.)
###############################################################################

theta_grid  <- seq(1, 2, by = 0.1)
joe_results <- list()

for (n in n_grid) {
  cat(sprintf("\n====== Joe copula power curve: n = %d ======\n", n))

  ## Fresh critical values under H0
  cat(sprintf("  Critical values (%d simulations)...\n", num_simulations))
  h0_joe <- mclapply_progress(
    seq_len(num_simulations),
    function(i) compute_stats(pobs(cbind(runif(n), runif(n)))),
    mc.cores = ncores, label = "critical values"
  )
  h0_joe_mat <- do.call(rbind, h0_joe)
  cvRnbar    <- quantile(h0_joe_mat[, "Rnbar"], 0.95)
  cvBn       <- quantile(h0_joe_mat[, "Bn"],    0.95)
  cat(sprintf("  Critical value -- IECT: %.4f | ECT: %.4f\n", cvRnbar, cvBn))

  rows <- list()
  for (k in seq_along(theta_grid)) {
    theta <- theta_grid[k]
    joe   <- joeCopula(param = theta, dim = 2)
    pow_stats <- mclapply_progress(
      seq_len(num_power),
      function(j) compute_stats(pobs(rCopula(n, joe))),
      mc.cores = ncores,
      label = sprintf("Joe(theta=%.2f)", theta)
    )
    pow_mat <- do.call(rbind, pow_stats)
    rows[[k]] <- data.frame(
      n          = n,
      theta      = theta,
      Power_IECT = mean(pow_mat[, "Rnbar"] >= cvRnbar),
      Power_ECT  = mean(pow_mat[, "Bn"]    >= cvBn)
    )
  }
  joe_results[[as.character(n)]] <- do.call(rbind, rows)
}

joe_final <- do.call(rbind, joe_results)
rownames(joe_final) <- NULL

cat("\n=== Joe copula power curve (alpha = 0.05) ===\n")
print(joe_final, digits = 3, row.names = FALSE)

## --- Power-curve plots ---
for (n in n_grid) {
  sub <- joe_final[joe_final$n == n, ]
  par(cex.lab = 1.7, cex.axis = 1.7, mgp = c(2.8, 0.5, 0))
  plot(sub$theta, sub$Power_IECT,
       type = "b", pch = 16, col = "steelblue", lwd = 2,
       ylim = c(0, 1),
       xlab = expression(theta), ylab = "Power",
       main = sprintf("Joe copula  (n = %d)", n),
       las = 1, cex.main = 1.9)
  lines(sub$theta, sub$Power_ECT,
        type = "b", pch = 17, col = "darkorange", lwd = 2)
  abline(h = 0.05, lty = 2, col = "gray60")
  legend("topleft",
         legend = c("IECT", "ECT"),
         col = c("steelblue", "darkorange"),
         pch = c(16, 17), lty = 1, lwd = 2, bty = "n", cex = 1.3)
}

###############################################################################
# 7. Null p-value histograms (IECT and ECT)
###############################################################################

par(cex.lab = 1.7, cex.axis = 1.7, mgp = c(2.9, 0.5, 0))

for (n in n_grid) {
  cat(sprintf("\n====== P-value histogram: n = %d ======\n", n))

  ## Null reference distribution
  cat(sprintf("  Null reference distribution (%d simulations)...\n", num_simulations))
  null_ref <- mclapply_progress(
    seq_len(num_simulations),
    function(i) compute_stats(pobs(cbind(runif(n), runif(n)))),
    mc.cores = ncores, label = sprintf("null ref (n=%d)", n)
  )
  null_mat <- do.call(rbind, null_ref)

  ## Separate H0 test statistics
  cat(sprintf("  H0 test statistics (%d simulations)...\n", num_power))
  h0_test <- mclapply_progress(
    seq_len(num_power),
    function(i) compute_stats(pobs(cbind(runif(n), runif(n)))),
    mc.cores = ncores, label = sprintf("H0 test (n=%d)", n)
  )
  test_mat <- do.call(rbind, h0_test)

  ## Compute p-values
  pval_IECT <- sapply(test_mat[, "Rnbar"], function(t)
    mean(null_mat[, "Rnbar"] >= t))
  pval_ECT  <- sapply(test_mat[, "Bn"], function(t)
    mean(null_mat[, "Bn"] >= t))

  ## Plot
  hist(pval_IECT, breaks = 20, freq = FALSE,
       main = bquote("IECT -- " * H[0] * " p-values  (n =" ~ .(n) * ")"),
       xlab = "p-value", col = "steelblue", border = "white",
       las = 1, cex.main = 1.9)
  abline(h = 1, lty = 2, col = "red", lwd = 2)

  hist(pval_ECT, breaks = 20, freq = FALSE,
       main = bquote("ECT -- " * H[0] * " p-values  (n =" ~ .(n) * ")"),
       xlab = "p-value", col = "darkorange", border = "white",
       las = 1, cex.main = 1.9)
  abline(h = 1, lty = 2, col = "red", lwd = 2)
}

###############################################################################
# 8. Real-data application: Uranium exploration (Li vs Ti)
###############################################################################

data(uranium, package = "copula")

v1 <- "Li"
v2 <- "Ti"
X  <- cbind(uranium[, v1], uranium[, v2])
U  <- pobs(X)
n  <- nrow(U)

## --- Summary statistics ---
tau_hat <- cor(X[, 1], X[, 2], method = "kendall")
rho_hat <- cor(X[, 1], X[, 2], method = "spearman")
cat(sprintf("Pair: %s vs %s  (n = %d)\n", v1, v2, n))
cat(sprintf("Kendall's tau:  %.4f\n", tau_hat))
cat(sprintf("Spearman's rho: %.4f\n", rho_hat))

## --- Empirical tail dependence ---
emp_upper_tail <- function(U, q = 0.9) {
  idx <- U[, 1] > q
  if (sum(idx) < 5) return(NA)
  mean(U[idx, 2] > q)
}

emp_lower_tail <- function(U, q = 0.1) {
  idx <- U[, 1] < q
  if (sum(idx) < 5) return(NA)
  mean(U[idx, 2] < q)
}

q_grid   <- seq(0.80, 0.95, by = 0.01)
lambda_U <- sapply(q_grid, function(q) emp_upper_tail(U, q))
lambda_L <- sapply(1 - q_grid, function(q) emp_lower_tail(U, q))

cat(sprintf("\nEmpirical upper tail dependence:\n"))
cat(sprintf("  lambda_U(0.85) = %.3f\n", emp_upper_tail(U, 0.85)))
cat(sprintf("  lambda_U(0.90) = %.3f\n", emp_upper_tail(U, 0.90)))
cat(sprintf("  lambda_U(0.95) = %.3f\n", emp_upper_tail(U, 0.95)))
cat(sprintf("Empirical lower tail dependence:\n"))
cat(sprintf("  lambda_L(0.15) = %.3f\n", emp_lower_tail(U, 0.15)))
cat(sprintf("  lambda_L(0.10) = %.3f\n", emp_lower_tail(U, 0.10)))
cat(sprintf("  lambda_L(0.05) = %.3f\n", emp_lower_tail(U, 0.05)))

## --- Exploratory scatter plots (2 x 2) ---
par(cex.lab = 1.7, cex.axis = 1.7, mgp = c(2.1, 0.5, 0))

## (a) Original scale
plot(X[, 1], X[, 2], pch = 16, cex = 0.4, col = "steelblue",
     xlab = "Li", ylab = "Ti",
     main = "(a) Original scale", cex.main = 1.9)

## (b) Pseudo-observations
plot(U[, 1], U[, 2], pch = 16, cex = 0.4, col = "steelblue",
     xlab = expression(hat(U)[1]), ylab = expression(hat(U)[2]),
     main = sprintf("(b) Pseudo observations (tau = %.3f)", tau_hat),
     cex.main = 1.9)

## (c) Upper tail zoom (U > 0.7)
upper_idx <- U[, 1] > 0.7 & U[, 2] > 0.7
plot(U[, 1], U[, 2], pch = 16, cex = 0.3, col = "steelblue",
     xlab = expression(hat(U)[1]), ylab = expression(hat(U)[2]),
     main = "(c) Upper tail region", cex.main = 1.9)
points(U[upper_idx, 1], U[upper_idx, 2], pch = 16, cex = 0.6, col = "red")
abline(h = 0.7, v = 0.7, lty = 2, col = "gray40")
legend("bottomright",
       legend = sprintf("%d / %d points", sum(upper_idx), n),
       text.col = "red", bty = "n", cex = 1.7)

## (d) Empirical tail dependence function
plot(q_grid, lambda_U,
     type = "b", pch = 16, col = "red", lwd = 2,
     ylim = c(0, max(0.3, max(lambda_U, lambda_L, na.rm = TRUE) + 0.05)),
     xlab = "Threshold q", ylab = expression(hat(lambda)(q)),
     main = "(d) Tail dependence", cex.main = 1.9)
lines(q_grid, lambda_L, type = "b", pch = 17, col = "blue", lwd = 2)
abline(h = 0, lty = 2, col = "gray60")
legend("topright",
       legend = c(expression(hat(lambda)[U](q)),
                  expression(hat(lambda)[L](1 - q))),
       col = c("red", "blue"), pch = c(16, 17),
       lty = 1, lwd = 2, bty = "n", cex = 1.7)

## --- Observed test statistics ---
obs_stats <- compute_stats(U)
cat(sprintf("\nObserved IECT statistic (Rnbar): %.6f\n", obs_stats["Rnbar"]))
cat(sprintf("Observed ECT  statistic (Bn):    %.6f\n", obs_stats["Bn"]))

## --- Monte Carlo null distribution ---
num_mc <- 5000
cat(sprintf("\nGenerating null distribution (%d replicates, n = %d)...\n",
            num_mc, n))

null_stats <- mclapply_progress(
  seq_len(num_mc),
  function(i) compute_stats(pobs(cbind(runif(n), runif(n)))),
  mc.cores = ncores, label = "null distribution"
)
null_mat <- do.call(rbind, null_stats)

## --- P-values and decisions ---
pval_IECT <- mean(null_mat[, "Rnbar"] >= obs_stats["Rnbar"])
pval_ECT  <- mean(null_mat[, "Bn"]    >= obs_stats["Bn"])

cat("\n====== Test Results: Uranium Li vs Ti ======\n")
cat(sprintf("  n = %d,  tau = %.4f,  lambda_U(0.90) = %.3f\n",
            n, tau_hat, emp_upper_tail(U, 0.90)))
cat(sprintf("  IECT:  Rnbar = %.6f,  p-value = %.4f  (%s)\n",
            obs_stats["Rnbar"], pval_IECT,
            ifelse(pval_IECT < 0.05, "Reject H0", "Fail to reject H0")))
cat(sprintf("  ECT:   Bn    = %.6f,  p-value = %.4f  (%s)\n",
            obs_stats["Bn"], pval_ECT,
            ifelse(pval_ECT < 0.05, "Reject H0", "Fail to reject H0")))

## --- Null distribution histograms with observed values ---
par(mfrow = c(1, 1), mar = c(5, 5, 3, 1),
    cex.lab = 1.7, cex.axis = 1.7, mgp = c(3.3, 0.8, 0))

## IECT
h1 <- hist(null_mat[, "Rnbar"], breaks = 40, plot = FALSE)
hist(null_mat[, "Rnbar"], breaks = 40, freq = FALSE,
     col = "steelblue", border = "white",
     main = "IECT null distribution",
     xlab = expression(bar(R)[n]), las = 1, cex.main = 1.9,
     xlim = c(0, max(obs_stats["Rnbar"] * 1.1, max(h1$breaks))))
abline(v = obs_stats["Rnbar"], col = "red", lwd = 2.5)
legend("topright",
       legend = c(bquote(bar(R)[n] == .(sprintf("%.5f", obs_stats["Rnbar"]))),
                  bquote(p * "-value" == .(sprintf("%.4f", pval_IECT)))),
       text.col = "red", bty = "n", cex = 1.7)

## ECT
h2 <- hist(null_mat[, "Bn"], breaks = 40, plot = FALSE)
hist(null_mat[, "Bn"], breaks = 40, freq = FALSE,
     col = "darkorange", border = "white",
     main = "ECT null distribution",
     xlab = expression(R[n]), las = 1, cex.main = 1.9,
     xlim = c(0, max(obs_stats["Bn"] * 1.1, max(h2$breaks))))
abline(v = obs_stats["Bn"], col = "red", lwd = 2.5)
legend("topright",
       legend = c(bquote(R[n] == .(sprintf("%.5f", obs_stats["Bn"]))),
                  bquote(p * "-value" == .(sprintf("%.4f", pval_ECT)))),
       text.col = "red", bty = "n", cex = 1.7)

## --- Copula fitting (for reference) ---
cat("\n====== Copula Fitting (for reference) ======\n")
fit_joe    <- fitCopula(joeCopula(),    data = U, method = "mpl")
fit_gumbel <- fitCopula(gumbelCopula(), data = U, method = "mpl")
fit_frank  <- fitCopula(frankCopula(),  data = U, method = "mpl")
fit_normal <- fitCopula(normalCopula(), data = U, method = "mpl")

fits <- list(Joe = fit_joe, Gumbel = fit_gumbel,
             Frank = fit_frank, Normal = fit_normal)
for (nm in names(fits)) {
  cat(sprintf("  %-8s  param = %6.3f   logLik = %8.2f   AIC = %8.2f  BIC = %8.2f\n",
              nm, coef(fits[[nm]])[1], logLik(fits[[nm]]),
              AIC(fits[[nm]]), BIC(fits[[nm]])))
}
