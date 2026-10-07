# =============================================================================
#  Title    : Regression analysis for establishing a relation between weather
#             parameters
#  Group ID : 19
#  Course   : Introduction to Data Analytics (IDA)
#  Members  : Aditya Palapati    S20230010176  (group leader)
#             Meghana Kuruva     S20230010134
#             Rohin Sai Bogadi   S20230010047
#             Srimeenakshi K S   S20262010002
#  Date     : 2026-10-07
#
#  Purpose
#    Find and quantify the relation between every pair of the daily weather
#    parameters of Delhi (mean temperature, humidity, wind speed, mean
#    pressure) with simple linear and simple non-linear regression, and
#    compare all models by R^2 (computed on the original scale).
#    (a) pre-processing, (b) pairwise relation analysis, (c) linear and
#    non-linear regression, (d) R^2 for every case and conclusions.
#
#  Data
#    data/DailyDelhiClimateTrain.csv  2013-01-01 .. 2017-01-01 (1462 days)
#    data/DailyDelhiClimateTest.csv   2017-01-01 .. 2017-04-24 (114 days)
#
#  Implementation
#    Base R only: no add-on packages and no ready-made statistics or modelling
#    functions (lm, nls, cor, var, sd, quantile, median, approx, solve, ...).
#    Every statistic is implemented by hand below, with the formulas in the
#    comments; pt() (Student-t CDF) is the only built-in statistical function
#    used.
#
#  How to run (from the root folder of the project, which contains data/):
#    Rscript R/Project_code.R
#  All paths are relative to that folder. Output: data/processed/clean.csv,
#  tables (.csv and LaTeX .tex) in outputs/tables/, figures (.png at 300 dpi
#  and .pdf) in outputs/figures/. Runtime about 1-2 minutes.
#
#  Sections
#    1  Setup & helper functions
#    2  Pre-processing
#    3  Relation analysis
#    4  Linear regression
#    5  Non-linear regression
#    6  R^2 comparison & conclusions
# =============================================================================


# =============================================================================
# 1  SETUP & HELPER FUNCTIONS
# =============================================================================

## ---- 1.1 Paths and global settings -----------------------------------------

set.seed(42)  # fixed seed for any step that uses randomness

if (!file.exists(file.path("data", "DailyDelhiClimateTrain.csv"))) {
  stop("Run this script from the project root: Rscript R/Project_code.R")
}

DATA_DIR      <- "data"
PROCESSED_DIR <- file.path("data", "processed")
FIG_DIR       <- file.path("outputs", "figures")
TAB_DIR       <- file.path("outputs", "tables")
for (d in c(PROCESSED_DIR, FIG_DIR, TAB_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

VARS <- c("meantemp", "humidity", "wind_speed", "meanpressure")
VAR_LABELS <- c(meantemp     = "Mean temperature (\u00B0C)",
                humidity     = "Humidity (%)",
                wind_speed   = "Wind speed (km/h)",
                meanpressure = "Mean pressure (hPa)")

# Plot colours: one main series colour, one highlight for corrected values,
# grey for raw (pre-cleaning) data.
COL_MAIN      <- "#2a78d6"
COL_HIGHLIGHT <- "#eb6834"
COL_RAW       <- "grey70"


## ---- 1.2 Descriptive statistics ---------------------------------------------

# my_mean: arithmetic mean
#   xbar = (1/n) * sum_{i=1..n} x_i
my_mean <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  sum(x) / length(x)
}

# my_var: unbiased sample variance (divisor n - 1, Bessel's correction)
#   s^2 = sum_{i=1..n} (x_i - xbar)^2 / (n - 1)
my_var <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  sum((x - my_mean(x))^2) / (n - 1)
}

# my_sd: sample standard deviation
#   s = sqrt(s^2)
my_sd <- function(x, na.rm = TRUE) {
  sqrt(my_var(x, na.rm))
}

# my_skewness: adjusted Fisher-Pearson sample skewness G1
#   m_k = (1/n) * sum (x_i - xbar)^k        (k-th central sample moment)
#   g1  = m_3 / m_2^(3/2)                     (moment coefficient of skewness)
#   G1  = g1 * sqrt(n (n - 1)) / (n - 2)      (small-sample bias adjustment)
#   G1 > 0: long right tail, G1 < 0: long left tail, |G1| < 0.5: ~symmetric.
my_skewness <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  n  <- length(x)
  d  <- x - my_mean(x)
  m2 <- sum(d^2) / n
  m3 <- sum(d^3) / n
  g1 <- m3 / m2^(3 / 2)
  g1 * sqrt(n * (n - 1)) / (n - 2)
}

# my_quantile: sample quantile, Hyndman & Fan type 7 (R's default)
#   Sort the data: x_(1) <= x_(2) <= ... <= x_(n). For probability p:
#     h = (n - 1) * p + 1,   j = floor(h),   g = h - j
#     Q(p) = x_(j) + g * (x_(j+1) - x_(j))        (with x_(n+1) := x_(n))
#   i.e. linear interpolation between the order statistics at positions
#   (k - 1) / (n - 1).
my_quantile <- function(x, probs, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  xs  <- sort(x)
  n   <- length(xs)
  out <- numeric(length(probs))
  for (k in seq_along(probs)) {
    h <- (n - 1) * probs[k] + 1
    j <- floor(h)
    g <- h - j
    upper  <- if (j < n) xs[j + 1] else xs[n]
    out[k] <- xs[j] + g * (upper - xs[j])
  }
  out
}

# my_median: middle value of the sorted data
#   n odd : median = x_((n+1)/2)
#   n even: median = (x_(n/2) + x_(n/2 + 1)) / 2
#   (identical to the type-7 quantile at p = 0.5)
my_median <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  xs <- sort(x)
  n  <- length(xs)
  if (n %% 2 == 1) xs[(n + 1) / 2] else (xs[n / 2] + xs[n / 2 + 1]) / 2
}

# my_zscore: standardised score of every observation
#   z_i = (x_i - xbar) / s
#   Under approximate normality |z| > 3 occurs with probability ~0.27 %.
#   Note: xbar and s are themselves inflated by outliers, so the z-score
#   rule is less robust than the IQR rule below.
my_zscore <- function(x) {
  (x - my_mean(x)) / my_sd(x)
}

# my_iqr_fences: Tukey's outlier fences
#   IQR   = Q3 - Q1
#   lower = Q1 - k * IQR,   upper = Q3 + k * IQR      (k = 1.5 by default)
#   Quartiles are robust (unaffected by the extreme values themselves).
my_iqr_fences <- function(x, k = 1.5) {
  q <- my_quantile(x, c(0.25, 0.75))
  iqr <- q[2] - q[1]
  list(q1 = q[1], q3 = q[2], iqr = iqr,
       lower = q[1] - k * iqr, upper = q[2] + k * iqr)
}

# my_rolling_median: centred moving median with window width k (k odd)
#   m_i = median{ x_j : |j - i| <= (k - 1) / 2 }   (window truncated at ends)
#   Used to judge whether an extreme value is an isolated spike relative to
#   its neighbouring days.
my_rolling_median <- function(x, k = 7) {
  n <- length(x)
  h <- (k - 1) %/% 2
  out <- numeric(n)
  for (i in seq_len(n)) {
    out[i] <- my_median(x[max(1, i - h):min(n, i + h)])
  }
  out
}

# my_boxplot_stats: components of a Tukey box plot from our own quartiles
#   box      = [Q1, Q3], centre line = median
#   whiskers = most extreme observations inside [Q1 - 1.5 IQR, Q3 + 1.5 IQR]
#   points   = observations outside the fences (potential outliers)
my_boxplot_stats <- function(x) {
  x <- x[!is.na(x)]
  q <- my_quantile(x, c(0.25, 0.5, 0.75))
  f <- my_iqr_fences(x)
  inside <- x[x >= f$lower & x <= f$upper]
  list(stats = c(min(inside), q[1], q[2], q[3], max(inside)),
       out   = x[x < f$lower | x > f$upper],
       n     = length(x))
}


# binned_means: a simple hand-made smoother for diagnostic plots
#   Sort the points by x, split them into k groups of (almost) equal size and
#   return the mean x and mean y of each group. A systematic trend in the
#   binned means of the residuals reveals curvature a model has missed.
binned_means <- function(x, y, k = 20) {
  o  <- order(x)
  g  <- ceiling(seq_along(o) * k / length(o))
  xs <- x[o]
  ys <- y[o]
  data.frame(x = sapply(1:k, function(j) my_mean(xs[g == j])),
             y = sapply(1:k, function(j) my_mean(ys[g == j])))
}


## ---- 1.3 Correlation --------------------------------------------------------

# my_rank: ranks 1..n, ties receive the average of the positions they occupy
#   If sorted values x_(i) = ... = x_(j) are tied, each gets rank (i + j) / 2.
#   (Spearman's rho = Pearson correlation of the ranks.)
my_rank <- function(x) {
  n  <- length(x)
  o  <- order(x)
  xs <- x[o]
  r  <- numeric(n)
  i  <- 1
  while (i <= n) {
    j <- i
    while (j < n && xs[j + 1] == xs[i]) j <- j + 1
    r[o[i:j]] <- (i + j) / 2
    i <- j + 1
  }
  r
}

# my_pearson: Pearson product-moment correlation coefficient
#   Sxx = sum (x_i - xbar)^2,  Syy = sum (y_i - ybar)^2
#   Sxy = sum (x_i - xbar)(y_i - ybar)
#   r   = Sxy / sqrt(Sxx * Syy),   -1 <= r <= 1
#   Only complete (x_i, y_i) pairs are used.
my_pearson <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  dx <- x - my_mean(x)
  dy <- y - my_mean(y)
  sum(dx * dy) / sqrt(sum(dx^2) * sum(dy^2))
}

# my_spearman: Spearman's rank correlation coefficient
#   rho = Pearson correlation of rank(x_i) and rank(y_i)   (average ranks for ties)
#   Measures MONOTONIC association; it equals Pearson's r only when the
#   relation is linear, so a clear gap |rho - r| hints at non-linearity
#   (or at the influence of extreme values on r).
my_spearman <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  my_pearson(my_rank(x[ok]), my_rank(y[ok]))
}

# cor_t_test: t-test of H0: population correlation = 0
#   t  = r * sqrt((n - 2) / (1 - r^2)),   df = n - 2
#   two-sided p-value = 2 * P(T_{n-2} <= -|t|)      (pt = Student-t CDF)
#   Assumes independent observations; with autocorrelated daily data the
#   effective sample size is smaller, so these p-values are optimistic.
cor_t_test <- function(r, n) {
  t <- r * sqrt((n - 2) / (1 - r^2))
  list(t = t, df = n - 2, p = 2 * pt(-abs(t), df = n - 2))
}

# fisher_ci: confidence interval for a correlation via Fisher's z-transformation
#   z  = atanh(r) = 0.5 * log((1 + r) / (1 - r))    (approximately normal)
#   SE = 1 / sqrt(n - 3)
#   [z_lo, z_hi] = z -/+ z_crit * SE   (z_crit = 1.96, the 97.5 % standard
#                                        normal quantile -> 95 % interval)
#   back-transform: r = tanh(z) = (exp(2 z) - 1) / (exp(2 z) + 1)
fisher_ci <- function(r, n, z_crit = 1.96) {
  z  <- 0.5 * log((1 + r) / (1 - r))
  se <- 1 / sqrt(n - 3)
  # back-transform: r = tanh(z) = (exp(2 z) - 1) / (exp(2 z) + 1)
  back <- function(z) (exp(2 * z) - 1) / (exp(2 * z) + 1)
  c(lower = back(z - z_crit * se), upper = back(z + z_crit * se))
}

# strength_label: verbal label for the size of a correlation coefficient
#   |r| < 0.1 negligible, < 0.3 weak, < 0.5 moderate, < 0.7 strong,
#   >= 0.7 very strong (common rule of thumb), plus the sign.
strength_label <- function(r) {
  a <- abs(r)
  size <- if (a < 0.1) "negligible" else if (a < 0.3) "weak" else
          if (a < 0.5) "moderate" else if (a < 0.7) "strong" else "very strong"
  if (a < 0.1) size else paste(size, if (r > 0) "positive" else "negative")
}

# my_qt: quantile of Student's t distribution, found by inverting pt()
#   Find q with pt(q, df) = p. pt() is strictly increasing in q, so bisection
#   works: keep a bracket [lo, hi] with pt(lo) < p <= pt(hi) and halve it
#   100 times (bracket width 100 / 2^100, far below machine precision).
#   df = Inf gives the standard normal distribution, so my_qt(p, Inf) = z_p.
#   (Vectorised over p.)
my_qt <- function(p, df) {
  lo <- rep(-50, length(p))
  hi <- rep(50, length(p))
  for (it in 1:100) {
    mid   <- (lo + hi) / 2
    below <- pt(mid, df) < p
    lo[below]  <- mid[below]
    hi[!below] <- mid[!below]
  }
  (lo + hi) / 2
}

# my_qnorm: standard normal quantile z_p = my_qt(p, df = Inf)
my_qnorm <- function(p) {
  my_qt(p, Inf)
}

# format_p: p-value for display: "< 0.001" or three decimals.
format_p <- function(p) {
  ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
}


## ---- 1.4 Interpolation and linear algebra -----------------------------------

# my_linear_interp: fill NA values by linear interpolation in time
#   For a missing value at time t with the nearest observed neighbours
#   (t0, v0) before and (t1, v1) after:
#     v(t) = v0 + (t - t0) / (t1 - t0) * (v1 - v0)
#   Time is the date as a day count, so unequal gaps are weighted correctly.
#   NAs before the first / after the last observation take the nearest
#   observed value (constant extrapolation).
my_linear_interp <- function(dates, values) {
  t   <- as.numeric(dates)
  obs <- which(!is.na(values))
  if (length(obs) < 2) stop("my_linear_interp: need at least two observed values")
  out <- values
  for (i in which(is.na(values))) {
    before <- obs[t[obs] < t[i]]
    after  <- obs[t[obs] > t[i]]
    if (length(before) == 0) {
      out[i] <- values[after[which.min(t[after])]]
    } else if (length(after) == 0) {
      out[i] <- values[before[which.max(t[before])]]
    } else {
      i0 <- before[which.max(t[before])]
      i1 <- after[which.min(t[after])]
      out[i] <- values[i0] + (t[i] - t[i0]) / (t[i1] - t[i0]) *
                             (values[i1] - values[i0])
    }
  }
  out
}

# solve_linear_system: solve A beta = b (A square) by Gaussian elimination
# with partial pivoting; used for the least-squares normal equations
# (X'X) beta = X'y.
#   Forward elimination, for column k = 1..n-1:
#     pivot row p = argmax_{i >= k} |A[i, k]|   (largest pivot limits
#                                                round-off error growth)
#     swap rows k and p of A and b
#     for each row i > k:  m = A[i, k] / A[k, k]
#                          A[i, ] <- A[i, ] - m * A[k, ]
#                          b[i]   <- b[i]   - m * b[k]
#   Back substitution, for i = n..1:
#     beta_i = (b_i - sum_{j > i} A[i, j] * beta_j) / A[i, i]
solve_linear_system <- function(A, b) {
  A <- as.matrix(A)
  b <- as.numeric(b)
  n <- nrow(A)
  if (ncol(A) != n || length(b) != n) stop("solve_linear_system: A must be n x n and b length n")
  tol <- 1e-12 * max(abs(A))
  for (k in seq_len(n)) {
    p <- (k - 1) + which.max(abs(A[k:n, k]))
    if (abs(A[p, k]) < tol) stop("solve_linear_system: matrix is (numerically) singular")
    if (p != k) {
      tmp <- A[k, ]; A[k, ] <- A[p, ]; A[p, ] <- tmp
      tmp <- b[k];   b[k]   <- b[p];   b[p]   <- tmp
    }
    if (k < n) {
      for (i in (k + 1):n) {
        m <- A[i, k] / A[k, k]
        A[i, ] <- A[i, ] - m * A[k, ]
        b[i]   <- b[i]   - m * b[k]
      }
    }
  }
  beta <- numeric(n)
  for (i in n:1) {
    s <- if (i < n) sum(A[i, (i + 1):n] * beta[(i + 1):n]) else 0
    beta[i] <- (b[i] - s) / A[i, i]
  }
  beta
}


# gauss_newton: non-linear least squares, minimise SSE(theta) = sum (y_i - f(x_i; theta))^2
#   Linearise f around the current estimate:  f(x; theta + delta) ~ f(x; theta) + J delta,
#   with Jacobian J_ij = d f(x_i; theta) / d theta_j. Least squares for delta
#   gives the normal equations
#     (J'J) delta = J'r,        r = y - f(x; theta)       (solved by solve_linear_system)
#   Step halving: try theta + lambda * delta for lambda = 1, 1/2, 1/4, ...
#   (at most 30 halvings) and accept the first step that lowers SSE, so SSE
#   decreases monotonically.
#   Converged when the relative decrease of SSE or the relative step size
#   max_j |lambda delta_j| / (|theta_j| + 1e-8) falls below tol. Not converged
#   after max_iter iterations, or when no step size lowers SSE while the step
#   is still large. f and jac are functions of (x, theta).
gauss_newton <- function(f, jac, theta, x, y, max_iter = 200, tol = 1e-10) {
  # SSE(theta) = sum (y_i - f(x_i; theta))^2; Inf if not finite (rejects the step)
  sse <- function(th) {
    v <- sum((y - f(x, th))^2)
    if (is.finite(v)) v else Inf
  }
  S <- sse(theta)
  if (!is.finite(S)) stop("SSE is not finite at the start values")
  for (iter in 1:max_iter) {
    r <- y - f(x, theta)
    J <- jac(x, theta)
    delta <- solve_linear_system(t(J) %*% J, t(J) %*% r)
    lambda <- 1
    accepted <- FALSE
    for (h in 0:30) {
      theta_new <- theta + lambda * delta
      S_new <- sse(theta_new)
      if (S_new < S) { accepted <- TRUE; break }
      lambda <- lambda / 2
    }
    rel_step <- max(abs(lambda * delta) / (abs(theta) + 1e-8))
    if (!accepted) {
      # no improvement possible: fine if we are already at the minimum
      ok <- max(abs(delta) / (abs(theta) + 1e-8)) < 1e-6
      return(list(theta = theta, SSE = S, iterations = iter, converged = ok,
                  message = if (ok) "converged" else "step halving failed"))
    }
    rel_dec <- (S - S_new) / S
    theta <- theta_new
    S <- S_new
    if (rel_dec < tol || rel_step < tol) {
      return(list(theta = theta, SSE = S, iterations = iter, converged = TRUE,
                  message = "converged"))
    }
  }
  list(theta = theta, SSE = S, iterations = max_iter, converged = FALSE,
       message = sprintf("no convergence in %d iterations", max_iter))
}


## ---- 1.5 Goodness-of-fit measures -------------------------------------------

# r2: coefficient of determination
#   SS_res = sum (y_i - yhat_i)^2,   SS_tot = sum (y_i - ybar)^2
#   R^2 = 1 - SS_res / SS_tot
#   (share of the variance of y explained by the model; computed on the
#   original scale of y, also for models fitted on a transformed scale)
r2 <- function(y, yhat) {
  1 - sum((y - yhat)^2) / sum((y - my_mean(y))^2)
}

# adj_r2: adjusted R^2, penalises the number of predictor terms
#   adj R^2 = 1 - (1 - R^2) * (n - 1) / (n - p - 1)
#   n = sample size, p = number of predictor terms excluding the intercept
#   (p = 1 for y = b0 + b1 x, p = 2 for a quadratic b0 + b1 x + b2 x^2).
adj_r2 <- function(y, yhat, n = length(y), p = 1) {
  1 - (1 - r2(y, yhat)) * (n - 1) / (n - p - 1)
}

# rmse: root mean squared error, in the units of y
#   RMSE = sqrt( (1/n) * sum (y_i - yhat_i)^2 )
rmse <- function(y, yhat) {
  sqrt(sum((y - yhat)^2) / length(y))
}

# mae: mean absolute error, in the units of y (less sensitive to outliers)
#   MAE = (1/n) * sum |y_i - yhat_i|
mae <- function(y, yhat) {
  sum(abs(y - yhat)) / length(y)
}


## ---- 1.6 Output helpers -----------------------------------------------------

# latex_escape: escape LaTeX special characters in table cells / headers.
# <, > and | become text commands, because with LaTeX's default (OT1) font
# encoding they would print as other glyphs.
latex_escape <- function(s) {
  s <- gsub("([&%$#_])", "\\\\\\1", s)
  s <- gsub("<", "\\textless{}", s, fixed = TRUE)
  s <- gsub(">", "\\textgreater{}", s, fixed = TRUE)
  s <- gsub("|", "\\textbar{}", s, fixed = TRUE)
  gsub("~", "$\\sim$", s, fixed = TRUE)   # a bare ~ is a non-breaking space in LaTeX
}

# format_column: numbers -> fixed decimals (whole-number columns without
# decimals), NA -> "--", text is escaped for LaTeX.
format_column <- function(x, digits) {
  if (is.numeric(x)) {
    x <- round(x, digits) + 0   # "+ 0" turns -0 into 0, so no "-0.000" is printed
    whole <- all(is.na(x) | abs(x - round(x)) < 1e-9)
    out <- if (whole) formatC(x, format = "d", big.mark = "")
           else formatC(x, format = "f", digits = digits)
    out[is.na(x)] <- "--"
    out
  } else {
    latex_escape(as.character(x))
  }
}

# write_latex_table: write a data frame as a booktabs-style LaTeX table
# (\toprule / \midrule / \bottomrule) inside a table float with caption and
# label, ready for \input{} in the report (needs \usepackage{booktabs}).
# Numeric columns are right-aligned, text columns left-aligned. The caption
# is written verbatim, so it may contain LaTeX (e.g. $R^2$).
# The tabular is wrapped in \resizebox{min(natural width, \linewidth)}{!}{...},
# so a wide table is scaled down to the text width and a narrow one is left
# unchanged (needs \usepackage{graphicx}).
write_latex_table <- function(df, file, caption, label, digits = 3) {
  align <- paste(ifelse(sapply(df, is.numeric), "r", "l"), collapse = "")
  cells <- sapply(df, format_column, digits = digits)
  if (is.null(dim(cells))) cells <- matrix(cells, nrow = 1)
  con <- file(file, open = "w", encoding = "UTF-8")
  on.exit(close(con))
  cat("% Generated by R/Project_code.R -- requires \\usepackage{booktabs,graphicx}\n",
      "\\begin{table}[htbp]\n",
      "  \\centering\n",
      "  \\caption{", caption, "}\n",
      "  \\label{", label, "}\n",
      "  \\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{%\n",
      "  \\begin{tabular}{", align, "}\n",
      "    \\toprule\n",
      "    ", paste(sapply(names(df), latex_header), collapse = " & "), " \\\\\n",
      "    \\midrule\n",
      sep = "", file = con)
  for (i in seq_len(nrow(cells))) {
    cat("    ", paste(cells[i, ], collapse = " & "), " \\\\\n", sep = "", file = con)
  }
  cat("    \\bottomrule\n",
      "  \\end{tabular}}\n",
      "\\end{table}\n",
      sep = "", file = con)
  invisible(file)
}

# latex_header: column name -> LaTeX header. Underscores become spaces
# (except in variable names), and R2 / Delta R2 / rho / n eff become math.
latex_header <- function(h) {
  if (!(h %in% VARS)) h <- gsub("_", " ", h, fixed = TRUE)
  s <- latex_escape(h)
  s <- gsub("Delta R2", "$\\Delta R^2$", s, fixed = TRUE)
  s <- gsub("(?<![A-Za-z^])R2(?![0-9])", "$R^2$", s, perl = TRUE)
  s <- gsub("rho", "$\\rho$", s, fixed = TRUE)
  gsub("n eff", "$n_{\\mathrm{eff}}$", s, fixed = TRUE)
}

# save_table: write a table to outputs/tables/ as <name>.csv and <name>.tex.
# tex_df (optional) is a display version for LaTeX, e.g. with p-values shown
# as "< 0.001"; the .csv always keeps the full-precision numbers.
save_table <- function(df, name, caption, label = paste0("tab:", name), digits = 3,
                       tex_df = df) {
  write.csv(df, file.path(TAB_DIR, paste0(name, ".csv")), row.names = FALSE)
  write_latex_table(tex_df, file.path(TAB_DIR, paste0(name, ".tex")),
                    caption = caption, label = label, digits = digits)
  cat(sprintf("  saved table  %s/%s.{csv,tex}\n", TAB_DIR, name))
}

# save_figure: draw a figure twice, as 300-dpi PNG and as vector PDF, in
# outputs/figures/. plot_fun is a function with no arguments that draws it.
save_figure <- function(name, plot_fun, width = 10, height = 6) {
  png(file.path(FIG_DIR, paste0(name, ".png")),
      width = width, height = height, units = "in", res = 300)
  plot_fun()
  dev.off()
  pdf(file.path(FIG_DIR, paste0(name, ".pdf")), width = width, height = height)
  plot_fun()
  dev.off()
  cat(sprintf("  saved figure %s/%s.{png,pdf}\n", FIG_DIR, name))
}

# draw_boxes: side-by-side box plots of one or more groups. The statistics
# come from my_boxplot_stats (our quartiles and Tukey fences); bxp() only draws.
draw_boxes <- function(groups, names, main, ylab, fill, xlab = "") {
  st    <- lapply(groups, my_boxplot_stats)
  out   <- unlist(lapply(st, function(s) s$out))
  group <- unlist(lapply(seq_along(st), function(k) rep(k, length(st[[k]]$out))))
  bxp(list(stats = sapply(st, function(s) s$stats), n = sapply(st, function(s) s$n),
           out = if (is.null(out)) numeric(0) else out,
           group = if (is.null(group)) numeric(0) else group, names = names),
      main = main, ylab = ylab, xlab = xlab, boxfill = fill, whisklty = 1,
      staplewex = 0.5, outpch = 16, outcex = 0.5, outcol = "grey30",
      frame.plot = FALSE, boxwex = 0.6, las = 1)
}

LOG <- character(0)
# log_line: print one line of a step-by-step log and keep it in LOG, so the
# whole log can also be written to a text file (sprintf syntax).
log_line <- function(fmt, ...) {
  txt <- sprintf(fmt, ...)
  cat(txt, "\n", sep = "")
  LOG <<- c(LOG, txt)
}


# =============================================================================
# 2  PRE-PROCESSING
# =============================================================================

cat("\n================ 2  PRE-PROCESSING ================\n\n")

## ---- 2.1 Read the raw files --------------------------------------------------

train <- read.csv(file.path(DATA_DIR, "DailyDelhiClimateTrain.csv"), stringsAsFactors = FALSE)
test  <- read.csv(file.path(DATA_DIR, "DailyDelhiClimateTest.csv"),  stringsAsFactors = FALSE)
train$source <- "train"
test$source  <- "test"

log_line("Step 1 - Read data")
log_line("  Train file: %d rows x %d columns, %s to %s",
         nrow(train), ncol(train) - 1, min(train$date), max(train$date))
log_line("  Test file : %d rows x %d columns, %s to %s",
         nrow(test), ncol(test) - 1, min(test$date), max(test$date))
log_line("  Missing values: train %d, test %d", sum(is.na(train)), sum(is.na(test)))


## ---- 2.2 Combine, parse dates, sort ------------------------------------------

dat <- rbind(train, test)
dat$date <- as.Date(dat$date, format = "%Y-%m-%d")
if (any(is.na(dat$date))) stop("Some dates could not be parsed")
dat <- dat[order(dat$date), ]
rownames(dat) <- NULL
n_combined <- nrow(dat)

log_line("Step 2 - Combine and sort")
log_line("  Bound train and test with a 'source' column, parsed dates with as.Date(),")
log_line("  sorted by date: %d rows.", n_combined)


## ---- 2.3 Overlapping date ----------------------------------------------------
# 2017-01-01 is the last day of the train file and the first day of the test
# file. The two rows disagree:
#   train: 10.0 degC, 100 %, 0.0 km/h, 1016.0 hPa -- round numbers, 100 %
#          humidity and zero wind, i.e. a single (partial-day) reading rather
#          than a daily mean;
#   test : 15.91 degC, 85.87 %, 2.74 km/h, 59.0 hPa -- fractional values that
#          are averages over many readings, but an impossible pressure.
# Decision: keep the test row (a proper daily mean for three variables, and
# its pressure is treated in step 6 like every other invalid pressure) and
# drop the train row. The dropped train pressure (1016.0 hPa) is used only as
# an independent cross-check of the interpolated value.

dup_dates <- unique(dat$date[duplicated(dat$date)])
log_line("Step 3 - Overlapping dates")
log_line("  Dates present in both files: %d (%s)", length(dup_dates),
         paste(format(dup_dates), collapse = ", "))
overlap_train_pressure <- NA
for (dd in as.list(dup_dates)) {
  rows <- dat[dat$date == dd, ]
  for (r in seq_len(nrow(rows))) {
    log_line("    %-5s row: meantemp %.2f, humidity %.2f, wind_speed %.2f, meanpressure %.2f",
             rows$source[r], rows$meantemp[r], rows$humidity[r],
             rows$wind_speed[r], rows$meanpressure[r])
  }
  overlap_train_pressure <- rows$meanpressure[rows$source == "train"][1]
}
drop <- dat$date %in% dup_dates & dat$source == "train"
dat  <- dat[!drop, ]
rownames(dat) <- NULL
log_line("  Kept the test row (averaged daily values), dropped the train row")
log_line("  (single-reading stub: round values, 100%% humidity, 0 wind). Rows: %d -> %d.",
         n_combined, nrow(dat))

# Snapshot of the data before any value is changed ("before" in figures/tables).
raw <- dat


## ---- 2.4 Continuity of the daily date index ----------------------------------

full_dates    <- seq(min(dat$date), max(dat$date), by = "day")
missing_dates <- full_dates[!(full_dates %in% dat$date)]
n_dup_dates   <- sum(duplicated(dat$date))
log_line("Step 4 - Daily date index")
log_line("  Expected %d consecutive days from %s to %s; found %d rows.",
         length(full_dates), format(min(full_dates)), format(max(full_dates)), nrow(dat))
log_line("  Missing dates: %d, duplicate dates: %d.", length(missing_dates), n_dup_dates)
if (n_dup_dates > 0) stop("Duplicate dates remain after removing the overlap")
if (length(missing_dates) > 0) {
  # Insert empty rows; their values are filled by interpolation in step 8.
  filler <- data.frame(date = missing_dates, meantemp = NA, humidity = NA,
                       wind_speed = NA, meanpressure = NA, source = "inserted")
  dat <- rbind(dat, filler)
  dat <- dat[order(dat$date), ]
  rownames(dat) <- NULL
  log_line("  Inserted %d empty rows for the missing dates.", length(missing_dates))
} else {
  log_line("  The series is a complete daily sequence; nothing to insert.")
}


## ---- 2.5 Units and range checks ---------------------------------------------
# The dataset documentation gives degC, %, km/h and hPa. We confirm these
# units from the observed ranges:
#   - temperature 6-39: degC for Delhi (in degF the same days would be 43-102);
#   - humidity within 0-100: relative humidity in percent;
#   - wind speed mean ~7: km/h (= ~1.9 m/s, typical for Delhi); in m/s the
#     mean would be ~25 km/h, implausibly windy for every day of the year;
#   - pressure median ~1009: hPa (= mbar), mean sea-level pressure.
# Physical bounds: humidity in [0, 100], wind speed >= 0. Values breaking
# them would be set to NA and interpolated in step 8.

log_line("Step 5 - Units and range checks")
for (v in VARS) {
  log_line("  %-12s min %9.2f  median %8.2f  mean %8.2f  max %9.2f",
           v, min(dat[[v]], na.rm = TRUE), my_median(dat[[v]]),
           my_mean(dat[[v]]), max(dat[[v]], na.rm = TRUE))
}
log_line("  Units consistent with degC, %%, km/h and hPa (see ranges above).")
bad_hum  <- which(!is.na(dat$humidity) & (dat$humidity < 0 | dat$humidity > 100))
bad_wind <- which(!is.na(dat$wind_speed) & dat$wind_speed < 0)
log_line("  humidity outside [0, 100]: %d values", length(bad_hum))
log_line("  negative wind_speed       : %d values", length(bad_wind))
dat$humidity[bad_hum]    <- NA
dat$wind_speed[bad_wind] <- NA
log_line("  wind_speed exactly 0      : %d days (kept: calm days are physically possible)",
         sum(dat$wind_speed == 0, na.rm = TRUE))


## ---- 2.6 Invalid mean pressure ----------------------------------------------
# Justification of the plausible range 990-1030 hPa:
#   - all credible observations in the data lie between ~991 and ~1023 hPa;
#   - Delhi's daily mean sea-level pressure moves smoothly between ~995 hPa
#     (monsoon low) and ~1020 hPa (winter high);
#   - the flagged values (-3, 12, 59, 310, 634, 938, 946, 1350, 1353,
#     7679 hPa) are physically impossible at the surface or would mean a
#     severe cyclone over Delhi, which did not happen; their neighbouring days
#     are all normal, so they are recording/transcription errors.
# Treatment: set to NA and fill by linear interpolation in time. Pressure is
# strongly autocorrelated from day to day, so the mean of the neighbouring
# days is an accurate estimate, and the day's other three variables are kept
# (deleting rows would lose valid data). Only 10 of 1575 values are affected.

P_MIN <- 990
P_MAX <- 1030
dat$meanpressure_raw <- dat$meanpressure
bad_p <- which(!is.na(dat$meanpressure) &
               (dat$meanpressure < P_MIN | dat$meanpressure > P_MAX))
dat$meanpressure[bad_p] <- NA
dat$meanpressure <- my_linear_interp(dat$date, dat$meanpressure)
dat$pressure_imputed <- as.integer(seq_len(nrow(dat)) %in% bad_p)

log_line("Step 6 - Invalid meanpressure (outside %d-%d hPa)", P_MIN, P_MAX)
log_line("  %d invalid values set to NA and filled by linear interpolation in time:",
         length(bad_p))
log_line("    %-10s  %-6s  %12s  %12s", "date", "source", "old (hPa)", "new (hPa)")
for (i in bad_p) {
  log_line("    %-10s  %-6s  %12.2f  %12.2f", format(dat$date[i]), dat$source[i],
           dat$meanpressure_raw[i], dat$meanpressure[i])
}
if (!is.na(overlap_train_pressure)) {
  i <- which(dat$date == dup_dates[1])
  log_line("  Cross-check: interpolated %s = %.2f hPa vs %.2f hPa in the dropped train row.",
           format(dat$date[i]), dat$meanpressure[i], overlap_train_pressure)
}
p_valid <- dat$meanpressure_raw[-bad_p]
log_line("  Valid values range %.2f-%.2f hPa; after correction %.2f-%.2f hPa.",
         min(p_valid), max(p_valid), min(dat$meanpressure), max(dat$meanpressure))


## ---- 2.7 Extreme wind speeds ------------------------------------------------
# Diagnostics: Tukey IQR fences (robust) and z-scores (|z| > 3), plus a
# comparison with the centred 7-day rolling median to see whether the
# extremes are isolated one-day spikes or part of a windy spell.
#
# Decision: CAP (winsorise) wind_speed at the upper Tukey fence Q3 + 1.5 IQR.
# Justification:
#   - The largest values (42.2, 34.5, 33.3, 30.7 km/h) are DAILY MEANS; a
#     24-hour mean of 30-42 km/h would be a strong breeze (Beaufort 5-6)
#     blowing all day and night in Delhi, where daily means are typically
#     2-15 km/h.
#   - Most of them are isolated spikes, 4-17 times the 7-day median with
#     ordinary values on the days either side, which points to
#     sensor/averaging errors rather than a genuine windy spell.
#   - A few (e.g. 2013-03-01, 2014-06-12: ~2x the 7-day median) sit inside
#     windy spells and may be real (18-25 km/h can occur on dust-storm or
#     western-disturbance days). As we cannot prove every high value is
#     wrong, we cap instead of deleting: the day stays in the data (its temperature, humidity and
#     pressure are valid), it is still among the windiest days, and a single
#     point can no longer dominate a least-squares fit (high leverage).
#   - The IQR rule is used rather than the z-score rule because the mean and
#     standard deviation are themselves inflated by the extremes, while the
#     quartiles are not.
# The original values are kept in wind_speed_raw for transparency.

w        <- dat$wind_speed
w_fences <- my_iqr_fences(w)
w_z      <- my_zscore(w)
w_rmed   <- my_rolling_median(w, k = 7)
above_fence <- which(w > w_fences$upper)
above_z3    <- which(abs(w_z) > 3)

log_line("Step 7 - Extreme wind_speed values")
log_line("  IQR rule : Q1 = %.2f, Q3 = %.2f, IQR = %.2f -> upper fence Q3 + 1.5 IQR = %.2f km/h; %d values above.",
         w_fences$q1, w_fences$q3, w_fences$iqr, w_fences$upper, length(above_fence))
log_line("  z-score  : mean = %.2f, sd = %.2f -> |z| > 3 means > %.2f km/h; %d values.",
         my_mean(w), my_sd(w), my_mean(w) + 3 * my_sd(w), length(above_z3))
log_line("  Values with |z| > 3 compared with the 7-day rolling median:")
log_line("    %-10s  %10s  %6s  %14s  %6s", "date", "km/h", "z", "7-day median", "ratio")
for (i in above_z3[order(-w[above_z3])]) {
  log_line("    %-10s  %10.2f  %6.2f  %14.2f  %6.1f", format(dat$date[i]), w[i], w_z[i],
           w_rmed[i], w[i] / w_rmed[i])
}
WIND_CAP <- w_fences$upper
dat$wind_speed_raw <- dat$wind_speed
dat$wind_capped    <- as.integer(!is.na(w) & w > WIND_CAP)
dat$wind_speed[dat$wind_capped == 1] <- WIND_CAP
log_line("  Decision: most extremes are isolated spikes (up to ~17x the 7-day median) that are")
log_line("  implausible as 24-hour means, but a few (ratio ~2) lie in windy spells and may be real.")
log_line("  So we cap at the upper IQR fence %.2f km/h (winsorising) instead of deleting: no", WIND_CAP)
log_line("  single day can dominate a least-squares fit, and the days' other variables are kept.")
log_line("  %d values capped.",
         sum(dat$wind_capped))
log_line("  wind_speed max %.2f -> %.2f km/h, sd %.2f -> %.2f, skewness %.2f -> %.2f.",
         max(w, na.rm = TRUE), max(dat$wind_speed), my_sd(w), my_sd(dat$wind_speed),
         my_skewness(w), my_skewness(dat$wind_speed))


## ---- 2.8 Fill any remaining gaps ----------------------------------------------
# Covers values set to NA in step 5 and rows inserted in step 4 (none in this
# data set, but the script stays correct if the input changes).

log_line("Step 8 - Remaining missing values")
n_filled <- integer(length(VARS))
names(n_filled) <- VARS
for (v in VARS) {
  n_filled[v] <- sum(is.na(dat[[v]]))
  if (n_filled[v] > 0) dat[[v]] <- my_linear_interp(dat$date, dat[[v]])
}
log_line("  Interpolated: %s", paste(sprintf("%s %d", VARS, n_filled), collapse = ", "))


## ---- 2.9 Calendar features ---------------------------------------------------
# Seasons follow the India Meteorological Department (IMD) definition:
#   Winter Dec-Feb, Summer (pre-monsoon) Mar-May, Monsoon Jun-Sep,
#   Post-monsoon Oct-Nov.

# season_of_month: IMD season of a month number (1-12)
season_of_month <- function(m) {
  ifelse(m %in% c(12, 1, 2), "Winter",
  ifelse(m %in% 3:5,         "Summer",
  ifelse(m %in% 6:9,         "Monsoon", "Post-monsoon")))
}
dat$month      <- as.integer(format(dat$date, "%m"))
dat$month_name <- month.abb[dat$month]
dat$season     <- season_of_month(dat$month)

log_line("Step 9 - Calendar features")
log_line("  Added month (1-12), month_name and season (IMD: Winter Dec-Feb, Summer Mar-May,")
log_line("  Monsoon Jun-Sep, Post-monsoon Oct-Nov). Days per season: %s.",
         paste(sprintf("%s %d", c("Winter", "Summer", "Monsoon", "Post-monsoon"),
                       sapply(c("Winter", "Summer", "Monsoon", "Post-monsoon"),
                              function(s) sum(dat$season == s))), collapse = ", "))


## ---- 2.10 Final validation and save -----------------------------------------

stopifnot(nrow(dat) == length(full_dates),
          all(dat$date == full_dates),
          !any(is.na(dat[VARS])),
          all(dat$humidity >= 0 & dat$humidity <= 100),
          all(dat$wind_speed >= 0),
          all(dat$meanpressure >= P_MIN & dat$meanpressure <= P_MAX))

clean <- data.frame(date = format(dat$date), source = dat$source,
                    dat[VARS],
                    month = dat$month, month_name = dat$month_name, season = dat$season,
                    meanpressure_raw = dat$meanpressure_raw, pressure_imputed = dat$pressure_imputed,
                    wind_speed_raw = dat$wind_speed_raw, wind_capped = dat$wind_capped)
# write.csv keeps only 15 significant digits. The raw data contain values that
# differ only in the 16th-17th digit, and rounding them would create artificial
# ties when clean.csv is read back (changing ranks and Spearman's rho). Numbers
# are therefore written with 17 significant digits, which reproduces every
# double exactly (lossless round trip).
num_cols  <- sapply(clean, is.numeric)
clean_out <- clean
clean_out[num_cols] <- lapply(clean[num_cols], function(v) sprintf("%.17g", as.numeric(v)))
write.csv(clean_out, file.path(PROCESSED_DIR, "clean.csv"), row.names = FALSE,
          quote = which(!num_cols))

log_line("Step 10 - Final checks and output")
log_line("  Final data: %d rows (%d consecutive days, %s to %s), no missing values,",
         nrow(dat), length(full_dates), format(min(dat$date)), format(max(dat$date)))
log_line("  all values within physical bounds. Saved to %s/clean.csv.", PROCESSED_DIR)
cat("\n")


## ---- 2.11 Summary tables -----------------------------------------------------

prep_summary <- data.frame(
  Step = c("Read", "Read", "Combine", "Overlap", "Date index", "Date index",
           "Range check", "Range check", "Pressure", "Wind", "Wind", "Final"),
  Check = c("Rows in train file", "Rows in test file", "Rows after combining",
            "Overlapping dates removed (train row dropped)", "Missing dates inserted",
            "Duplicate dates", "Humidity outside 0-100", "Negative wind speed",
            sprintf("Outside %d-%d hPa, interpolated", P_MIN, P_MAX),
            "Above Q3 + 1.5 IQR, capped", "Diagnostic: |z| > 3",
            "Rows after pre-processing"),
  Variable = c("all", "all", "all", "all", "all", "all", "humidity", "wind_speed",
               "meanpressure", "wind_speed", "wind_speed", "all"),
  Count = c(nrow(train), nrow(test), n_combined, length(dup_dates), length(missing_dates),
            n_dup_dates, length(bad_hum), length(bad_wind), length(bad_p),
            sum(dat$wind_capped), length(above_z3), nrow(dat)),
  stringsAsFactors = FALSE)

cat("Saving outputs\n")
save_table(prep_summary, "preprocessing_summary",
           caption = "Pre-processing steps: rows before/after, values corrected per variable and outliers handled.",
           label = "tab:preprocessing_summary")

# Descriptive statistics of every variable before (raw, after removing the
# overlapping row) and after cleaning.
describe <- function(x, variable, stage) {
  q <- my_quantile(x, c(0.25, 0.5, 0.75))
  data.frame(Variable = variable, Stage = stage, n = length(x),
             Mean = my_mean(x), SD = my_sd(x), Min = min(x), Q1 = q[1],
             Median = q[2], Q3 = q[3], Max = max(x), Skewness = my_skewness(x),
             stringsAsFactors = FALSE)
}
desc <- do.call(rbind, lapply(VARS, function(v) {
  rbind(describe(raw[[v]], v, "before"), describe(dat[[v]], v, "after"))
}))
save_table(desc, "preprocessing_descriptives",
           caption = "Descriptive statistics of each variable before and after pre-processing.",
           label = "tab:preprocessing_descriptives", digits = 2)

writeLines(LOG, file.path(TAB_DIR, "preprocessing_log.txt"))
cat(sprintf("  saved log    %s/preprocessing_log.txt\n", TAB_DIR))


## ---- 2.12 Figures ------------------------------------------------------------

# Figure: box plots before (top row) and after (bottom row) pre-processing.
# Each panel has its own y-axis so the cleaned distributions stay readable.
save_figure("boxplots_before_after", function() {
  op <- par(mfrow = c(2, 4), mar = c(1.5, 4.5, 3, 1), oma = c(0, 0, 2, 0), cex.main = 1)
  on.exit(par(op))
  for (v in VARS) draw_boxes(list(raw[[v]]), "", paste(v, "- before"), VAR_LABELS[v], COL_RAW)
  for (v in VARS) draw_boxes(list(dat[[v]]), "", paste(v, "- after"), VAR_LABELS[v], COL_MAIN)
  mtext("Distributions before (top) and after (bottom) pre-processing",
        outer = TRUE, cex = 1.1, font = 2)
}, width = 12, height = 7)

# Figure: cleaned daily time series, one panel per variable. Corrected values
# are marked: interpolated pressures and capped wind speeds.
save_figure("timeseries_clean", function() {
  op <- par(mfrow = c(4, 1), mar = c(2.5, 4.5, 2, 1), oma = c(0, 0, 2, 0), las = 1)
  on.exit(par(op))
  for (v in VARS) {
    # extra headroom in the wind panel so the legend sits above the cap line
    ylim <- range(dat[[v]])
    if (v == "wind_speed") ylim[2] <- ylim[2] * 1.35
    plot(dat$date, dat[[v]], type = "l", col = COL_MAIN, lwd = 0.8, ylim = ylim,
         xlab = "", ylab = "", main = VAR_LABELS[v], bty = "l")
    if (v == "meanpressure") {
      idx <- which(dat$pressure_imputed == 1)
      points(dat$date[idx], dat$meanpressure[idx], pch = 21, bg = COL_HIGHLIGHT, cex = 1.2)
      legend("bottomleft", legend = sprintf("interpolated (%d invalid values)", length(idx)),
             pch = 21, pt.bg = COL_HIGHLIGHT, pt.cex = 1.2, bty = "n", cex = 0.9)
    }
    if (v == "wind_speed") {
      idx <- which(dat$wind_capped == 1)
      abline(h = WIND_CAP, lty = 2, col = "grey40")
      points(dat$date[idx], dat$wind_speed[idx], pch = 21, bg = COL_HIGHLIGHT, cex = 1.2)
      legend("top", legend = c(sprintf("capped value (%d days)", length(idx)),
                               sprintf("cap = upper IQR fence, %.1f km/h", WIND_CAP)),
             pch = c(21, NA), pt.bg = COL_HIGHLIGHT, pt.cex = 1.2, lty = c(NA, 2),
             col = c("black", "grey40"), horiz = TRUE, bty = "n", cex = 0.9)
    }
  }
  mtext("Daily Delhi climate after pre-processing (2013-01-01 to 2017-04-24)",
        outer = TRUE, cex = 1.1, font = 2)
}, width = 11, height = 10)


# =============================================================================
# 3  RELATION ANALYSIS
# =============================================================================
# All statistics in this section use the cleaned data from section 2.

cat("\n================ 3  RELATION ANALYSIS ================\n\n")
LOG   <- character(0)   # new log: the plain-language relation summary
n_obs <- nrow(dat)
YEAR  <- as.integer(format(dat$date, "%Y"))


## ---- 3.1 Descriptive statistics ---------------------------------------------

desc_clean <- do.call(rbind, lapply(VARS, function(v) describe(dat[[v]], v, "clean")))
desc_clean$Stage <- NULL
cat("Descriptive statistics (cleaned data)\n")
print(cbind(desc_clean[1:2], round(desc_clean[-(1:2)], 2)), row.names = FALSE)
cat("\n")
save_table(desc_clean, "descriptive_stats",
           caption = sprintf("Descriptive statistics of the cleaned daily data ($n = %d$ days).", n_obs),
           label = "tab:descriptive_stats", digits = 2)


## ---- 3.2 Histograms ------------------------------------------------------------

save_figure("histograms", function() {
  op <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1), las = 1)
  on.exit(par(op))
  for (v in VARS) {
    x  <- dat[[v]]
    sk <- my_skewness(x)
    hist(x, breaks = 30, col = COL_MAIN, border = "white",
         main = sprintf("%s  (skewness %.2f)", VAR_LABELS[v], sk),
         xlab = VAR_LABELS[v], ylab = "Number of days")
    abline(v = my_mean(x),   lwd = 2, lty = 2)
    abline(v = my_median(x), lwd = 2, lty = 1)
    legend(if (sk < 0) "topleft" else "topright",
           legend = c(sprintf("mean %.1f", my_mean(x)), sprintf("median %.1f", my_median(x))),
           lty = c(2, 1), lwd = 2, bty = "n", cex = 0.85)
  }
}, width = 10, height = 7)


## ---- 3.3 Correlation matrices ----------------------------------------------------

P_MAT <- matrix(1, length(VARS), length(VARS), dimnames = list(VARS, VARS))
S_MAT <- P_MAT
for (i in 1:(length(VARS) - 1)) {
  for (j in (i + 1):length(VARS)) {
    P_MAT[i, j] <- P_MAT[j, i] <- my_pearson(dat[[VARS[i]]], dat[[VARS[j]]])
    S_MAT[i, j] <- S_MAT[j, i] <- my_spearman(dat[[VARS[i]]], dat[[VARS[j]]])
  }
}
cat("Pearson correlation matrix\n");  print(round(P_MAT, 3)); cat("\n")
cat("Spearman correlation matrix\n"); print(round(S_MAT, 3)); cat("\n")

# matrix_table: correlation matrix -> data frame with a Variable column (for save_table)
matrix_table <- function(M) {
  data.frame(Variable = rownames(M), M, check.names = FALSE, row.names = NULL)
}
save_table(matrix_table(P_MAT), "correlation_pearson",
           caption = "Pearson correlation matrix of the four weather parameters.",
           label = "tab:correlation_pearson")
save_table(matrix_table(S_MAT), "correlation_spearman",
           caption = "Spearman rank correlation matrix of the four weather parameters.",
           label = "tab:correlation_spearman")


## ---- 3.4 Pairwise inference, ranked by |r| -----------------------------------------
# For each of the 6 pairs: Pearson r with its t-test and 95 % Fisher-z CI,
# Spearman rho, and the gap rho - r. A gap of at least 0.10 is flagged as a
# hint of a non-linear (but monotonic) relation or of influential points.
#
# Extra column: daily weather is strongly autocorrelated (today resembles
# yesterday), so the n = 1575 days are not independent and the usual p-value
# is too optimistic. As a check we also report the effective sample size for
# two AR(1)-like series (Bartlett 1935; Bretherton et al. 1999):
#   n_eff = n * (1 - a_x * a_y) / (1 + a_x * a_y)
# where a_x, a_y are the lag-1 autocorrelations, and the p-value with n_eff.

# my_acf1: lag-1 autocorrelation, the Pearson correlation of consecutive days
#   a = r( (x_1, ..., x_{n-1}), (x_2, ..., x_n) )
my_acf1 <- function(x) {
  n <- length(x)
  my_pearson(x[1:(n - 1)], x[2:n])
}
ACF1 <- sapply(VARS, function(v) my_acf1(dat[[v]]))

DISAGREE <- 0.10
pair_rows <- list()
for (i in 1:(length(VARS) - 1)) {
  for (j in (i + 1):length(VARS)) {
    r     <- P_MAT[i, j]
    rho   <- S_MAT[i, j]
    tt    <- cor_t_test(r, n_obs)
    ci    <- fisher_ci(r, n_obs)
    aa    <- ACF1[i] * ACF1[j]
    n_eff <- n_obs * (1 - aa) / (1 + aa)
    pair_rows[[length(pair_rows) + 1]] <- data.frame(
      Pair = paste(VARS[i], "vs", VARS[j]), x = VARS[i], y = VARS[j],
      r = r, t = tt$t, df = tt$df, p = tt$p, CI_lower = ci["lower"], CI_upper = ci["upper"],
      rho = rho, rho_minus_r = rho - r, Strength = strength_label(r),
      n_eff = n_eff, p_eff = cor_t_test(r, n_eff)$p,
      stringsAsFactors = FALSE, row.names = NULL)
  }
}
pair_tab <- do.call(rbind, pair_rows)
pair_tab <- pair_tab[order(-abs(pair_tab$r)), ]
pair_tab <- data.frame(Rank = seq_len(nrow(pair_tab)), pair_tab, row.names = NULL)

cat("All pairs ranked by |Pearson r|\n")
print(data.frame(pair_tab[c("Rank", "Pair")], round(pair_tab[c("r", "CI_lower", "CI_upper", "rho")], 3),
                 t = round(pair_tab$t, 2), p = format_p(pair_tab$p),
                 n_eff = round(pair_tab$n_eff), p_eff = format_p(pair_tab$p_eff),
                 Strength = pair_tab$Strength), row.names = FALSE)
cat("\n")

pair_tex <- data.frame(
  Rank = pair_tab$Rank, Pair = pair_tab$Pair, `Pearson r` = pair_tab$r,
  `95% CI` = sprintf("[%.3f, %.3f]", pair_tab$CI_lower, pair_tab$CI_upper),
  t = sprintf("%.2f", pair_tab$t), p = format_p(pair_tab$p),
  `Spearman rho` = pair_tab$rho, `rho - r` = pair_tab$rho_minus_r,
  `n eff` = round(pair_tab$n_eff), `p (n eff)` = format_p(pair_tab$p_eff),
  check.names = FALSE, stringsAsFactors = FALSE)
save_table(pair_tab[c("Rank", "Pair", "x", "y", "r", "t", "df", "p", "CI_lower", "CI_upper",
                      "rho", "rho_minus_r", "Strength", "n_eff", "p_eff")],
           "correlation_pairs_ranked",
           caption = paste0("All six pairs ranked by $|r|$: Pearson $r$ with $t$-test ",
                            "($t = r\\sqrt{(n-2)/(1-r^2)}$, $df = n-2$) and 95\\% Fisher-$z$ ",
                            "interval, Spearman $\\rho$, and the $p$-value recomputed with the ",
                            "autocorrelation-adjusted effective sample size $n_{\\mathrm{eff}}$."),
           label = "tab:correlation_pairs_ranked", tex_df = pair_tex)


## ---- 3.5 Correlation heatmaps ----------------------------------------------------
# Diverging colour scale: blue (-1) -> light grey (0) -> orange (+1).

CORR_PAL <- colorRampPalette(c("#2a78d6", "#f0efec", "#eb6834"))(201)

# draw_corr_heatmap: correlation matrix as coloured cells (-1 blue .. +1 orange)
# with the coefficient printed in each cell
draw_corr_heatmap <- function(M, main) {
  p <- nrow(M)
  # image() puts z[i, j] at (x_i, y_j); flip rows so row 1 is drawn at the top
  image(1:p, 1:p, t(M[p:1, ]), zlim = c(-1, 1), col = CORR_PAL,
        axes = FALSE, xlab = "", ylab = "", main = main)
  axis(1, at = 1:p, labels = colnames(M), tick = FALSE, cex.axis = 0.9)
  axis(2, at = 1:p, labels = rev(rownames(M)), tick = FALSE, las = 1, cex.axis = 0.9)
  # diagonal (r = 1 by definition) in neutral grey so it does not dominate
  for (k in 1:p) rect(k - 0.5, p - k + 0.5, k + 0.5, p - k + 1.5, col = "grey88", border = NA)
  abline(h = 0.5 + 0:p, v = 0.5 + 0:p, col = "white", lwd = 2)
  for (i in 1:p) for (j in 1:p) {
    text(j, p - i + 1, sprintf("%.2f", M[i, j]), cex = 1.2,
         font = if (i == j) 1 else 2)
  }
}

save_figure("correlation_heatmaps", function() {
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 0.22))
  op <- par(mar = c(3, 7.5, 3, 1))
  on.exit({ par(op); layout(1) })
  draw_corr_heatmap(P_MAT, "Pearson r")
  draw_corr_heatmap(S_MAT, "Spearman rho")
  # colour key
  par(mar = c(3, 0.5, 3, 3.5))
  key <- seq(-1, 1, length.out = 201)
  image(1, key, matrix(key, nrow = 1), col = CORR_PAL, axes = FALSE, xlab = "", ylab = "")
  axis(4, at = seq(-1, 1, 0.5), las = 1)
}, width = 12, height = 5)


## ---- 3.6 Scatterplot matrix ---------------------------------------------------------
# Lower panels: scatter plots; diagonal: histograms; upper panels: Pearson r
# and Spearman rho (font size grows with |r|).

# panel_points: lower panels, plain scatter plot
panel_points <- function(x, y, ...) {
  points(x, y, pch = 16, cex = 0.35, col = adjustcolor(COL_MAIN, alpha.f = 0.35))
}
# panel_corr: upper panels, Pearson r (my_pearson) and Spearman rho (my_spearman)
panel_corr <- function(x, y, ...) {
  usr <- par("usr"); on.exit(par(usr = usr))
  par(usr = c(0, 1, 0, 1))
  r <- my_pearson(x, y)
  text(0.5, 0.6, sprintf("r = %.2f", r), cex = 1 + 1.2 * abs(r), font = 2)
  text(0.5, 0.3, sprintf("rho = %.2f", my_spearman(x, y)), cex = 1.1, col = "grey30")
}
# panel_hist: diagonal panels, histogram of the variable
panel_hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr = usr))
  h <- hist(x, breaks = 25, plot = FALSE)
  par(usr = c(usr[1:2], 0, max(h$counts) * 1.3))
  rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$counts, col = COL_RAW, border = "white")
}

save_figure("scatterplot_matrix", function() {
  pairs(dat[VARS], labels = c("meantemp\n(°C)", "humidity\n(%)",
                              "wind_speed\n(km/h)", "meanpressure\n(hPa)"),
        lower.panel = panel_points, upper.panel = panel_corr, diag.panel = panel_hist,
        gap = 0.4, las = 1, cex.labels = 1.2, main = "Scatterplot matrix of the cleaned daily data")
}, width = 10, height = 10)


## ---- 3.7 Seasonal profiles ---------------------------------------------------------
# Mean (and +/- 1 SD band) of every variable by calendar month, and the
# distribution per year. 2017 only covers January to April.

monthly <- data.frame(Month = 1:12, Name = month.abb)
for (v in VARS) {
  monthly[[v]]                 <- sapply(1:12, function(m) my_mean(dat[[v]][dat$month == m]))
  monthly[[paste0(v, "_sd")]]  <- sapply(1:12, function(m) my_sd(dat[[v]][dat$month == m]))
}
monthly$days <- sapply(1:12, function(m) sum(dat$month == m))
save_table(monthly[c("Month", "Name", VARS, "days")], "monthly_means",
           caption = "Mean of each variable by calendar month (2013--2017).",
           label = "tab:monthly_means", digits = 2)

save_figure("monthly_profiles", function() {
  op <- par(mfrow = c(2, 2), mar = c(3, 4.5, 3, 1), oma = c(0, 0, 2.5, 0), las = 1)
  on.exit(par(op))
  for (v in VARS) {
    m <- monthly[[v]]
    s <- monthly[[paste0(v, "_sd")]]
    plot(1:12, m, type = "n", ylim = range(m - s, m + s), xaxt = "n",
         xlab = "", ylab = VAR_LABELS[v], main = VAR_LABELS[v], bty = "l")
    axis(1, at = 1:12, labels = month.abb, cex.axis = 0.85)
    polygon(c(1:12, 12:1), c(m - s, rev(m + s)),
            col = adjustcolor(COL_MAIN, alpha.f = 0.18), border = NA)
    lines(1:12, m, col = COL_MAIN, lwd = 2)
    points(1:12, m, pch = 21, bg = COL_MAIN, col = "white", cex = 1.4)
  }
  mtext("Monthly mean (line) and +/- 1 SD (band), 2013-2017", side = 3, outer = TRUE,
        line = 0.8, font = 2, cex = 1.1)
}, width = 10, height = 7)

YEARS <- sort(unique(YEAR))
save_figure("boxplots_by_year", function() {
  op <- par(mfrow = c(2, 2), mar = c(4, 4.5, 3, 1), las = 1)
  on.exit(par(op))
  for (v in VARS) {
    draw_boxes(lapply(YEARS, function(yr) dat[[v]][YEAR == yr]),
               names = ifelse(YEARS == 2017, "2017*", as.character(YEARS)),
               main = VAR_LABELS[v], ylab = VAR_LABELS[v], fill = COL_MAIN)
  }
  mtext("* 2017: January-April only", side = 1, outer = TRUE, line = -1, adj = 0.98, cex = 0.8)
}, width = 10, height = 7)


## ---- 3.8 Plain-language summary --------------------------------------------------------

# describe_pair: one plain-language sentence for row k of the ranked pair table
describe_pair <- function(k) {
  row <- pair_tab[k, ]
  direction <- if (row$r > 0) "higher" else "lower"
  sprintf("%s: r = %.3f (95%% CI %.3f to %.3f), rho = %.3f -> %s. Days with higher %s tend to have %s %s.",
          row$Pair, row$r, row$CI_lower, row$CI_upper, row$rho, row$Strength,
          row$x, direction, row$y)
}

cat("\n")
log_line("Relation analysis summary (n = %d days)", n_obs)
k_pos <- which.max(pair_tab$r)
k_neg <- which.min(pair_tab$r)
log_line("Strongest positive relation:")
log_line("  %s", describe_pair(k_pos))
log_line("Strongest negative relation:")
log_line("  %s", describe_pair(k_neg))
log_line("All six pairs, ranked by |r|:")
for (k in seq_len(nrow(pair_tab))) {
  log_line("  %d. %s: r = %6.3f, rho = %6.3f, %s (t = %.2f, p %s)", k, pair_tab$Pair[k],
           pair_tab$r[k], pair_tab$rho[k], pair_tab$Strength[k], pair_tab$t[k],
           ifelse(pair_tab$p[k] < 0.001, "< 0.001", sprintf("= %.3f", pair_tab$p[k])))
}

flag <- which(abs(pair_tab$rho_minus_r) >= DISAGREE)
log_line("Pearson vs Spearman (|rho - r| >= %.2f suggests a non-linear monotonic relation", DISAGREE)
log_line("or influential points, to be examined with non-linear models in section 5):")
if (length(flag) == 0) {
  k_gap <- which.max(abs(pair_tab$rho_minus_r))
  log_line("  none: Pearson and Spearman agree for every pair (largest gap %+.3f for %s).",
           pair_tab$rho_minus_r[k_gap], pair_tab$Pair[k_gap])
  log_line("  Note: both coefficients only measure MONOTONIC association. A non-monotonic")
  log_line("  pattern (e.g. U-shaped, with high values at both ends of x) is invisible to both,")
  log_line("  so the scatterplot matrix and the non-linear fits in section 5 are still needed.")
} else {
  for (k in flag) {
    log_line("  %s: r = %.3f vs rho = %.3f (difference %+.3f) -> %s",
             pair_tab$Pair[k], pair_tab$r[k], pair_tab$rho[k], pair_tab$rho_minus_r[k],
             if (abs(pair_tab$rho[k]) > abs(pair_tab$r[k]))
               "rank correlation stronger: monotonic but curved relation likely"
             else "linear correlation stronger: likely driven by a few extreme points")
  }
}

log_line("Caution - autocorrelation: consecutive days are not independent (lag-1")
log_line("autocorrelation %s).", paste(sprintf("%s %.2f", VARS, ACF1), collapse = ", "))
log_line("With the effective sample size n_eff the p-values become:")
for (k in seq_len(nrow(pair_tab))) {
  log_line("  %s: n_eff = %4.0f, p %s", pair_tab$Pair[k], pair_tab$n_eff[k],
           ifelse(pair_tab$p_eff[k] < 0.001, "< 0.001", sprintf("= %.3f", pair_tab$p_eff[k])))
}
log_line("Pairs still significant at 5%% after this adjustment: %s.",
         if (any(pair_tab$p_eff < 0.05)) paste(pair_tab$Pair[pair_tab$p_eff < 0.05], collapse = ", ")
         else "none")
log_line("Much of the shared variation is seasonal (see monthly_profiles): correlation shows")
log_line("association through the common annual cycle, not that one variable causes the other.")

writeLines(LOG, file.path(TAB_DIR, "relation_summary.txt"))
cat(sprintf("  saved log    %s/relation_summary.txt\n", TAB_DIR))


# =============================================================================
# 4  LINEAR REGRESSION
# =============================================================================

cat("\n================ 4  LINEAR REGRESSION ================\n\n")
LOG <- character(0)   # new log: linear regression summary


## ---- 4.1 Simple linear regression by ordinary least squares -------------------

# my_slr: fit y_i = a + b * x_i + e_i by ordinary least squares (OLS)
#   xbar, ybar = sample means
#   Sxx = sum (x_i - xbar)^2
#   Sxy = sum (x_i - xbar) * (y_i - ybar)
#   b = Sxy / Sxx               slope (solves d SSE / d b = 0)
#   a = ybar - b * xbar         intercept (the line passes through (xbar, ybar))
#   yhat_i = a + b * x_i        fitted values
#   e_i    = y_i - yhat_i       residuals
#   SSE = sum e_i^2                     error (residual) sum of squares
#   SST = sum (y_i - ybar)^2            total sum of squares
#   SSR = sum (yhat_i - ybar)^2         regression sum of squares (SST = SSR + SSE)
#   R^2 = 1 - SSE / SST                 (= r_xy^2 in simple linear regression)
#   adj R^2 = 1 - (1 - R^2)(n - 1)/(n - 2)
#   RMSE = sqrt(SSE / n),  MAE = (1/n) sum |e_i|
#   s = sqrt(SSE / (n - 2))             residual standard error (2 parameters estimated)
#   SE(b) = s / sqrt(Sxx),   SE(a) = s * sqrt(1/n + xbar^2 / Sxx)
#   t_b = b / SE(b) ~ t_{n-2} under H0: slope = 0,   p = 2 * P(T_{n-2} <= -|t_b|)
#   95 % CI of the slope: b -/+ t_{0.975, n-2} * SE(b)
#   Durbin-Watson (residuals in time order):
#     DW = sum_{t=2..n} (e_t - e_{t-1})^2 / sum_{t=1..n} e_t^2  ~  2 (1 - rho_1)
#     DW ~ 2: no autocorrelation; DW < 2: positive autocorrelation (DW < 1.5
#     is a common warning level), which makes SE(b) and p-values too small.
my_slr <- function(x, y) {
  n    <- length(y)
  xbar <- my_mean(x)
  ybar <- my_mean(y)
  Sxx  <- sum((x - xbar)^2)
  Sxy  <- sum((x - xbar) * (y - ybar))
  b    <- Sxy / Sxx
  a    <- ybar - b * xbar
  yhat <- a + b * x
  e    <- y - yhat
  SSE  <- sum(e^2)
  SST  <- sum((y - ybar)^2)
  SSR  <- sum((yhat - ybar)^2)
  s    <- sqrt(SSE / (n - 2))
  se_b <- s / sqrt(Sxx)
  se_a <- s * sqrt(1 / n + xbar^2 / Sxx)
  t_b  <- b / se_b
  tcrit <- my_qt(0.975, n - 2)
  list(n = n, a = a, b = b, se_a = se_a, se_b = se_b, t_b = t_b,
       p_b = 2 * pt(-abs(t_b), df = n - 2),
       ci_b = c(b - tcrit * se_b, b + tcrit * se_b),
       fitted = yhat, resid = e, SSE = SSE, SST = SST, SSR = SSR,
       R2 = 1 - SSE / SST, adjR2 = adj_r2(y, yhat, n, p = 1),
       RMSE = rmse(y, yhat), MAE = mae(y, yhat), RSE = s,
       DW = sum((e[-1] - e[-n])^2) / sum(e^2))
}

# slr_predict: predictions of a fitted line for new x values
#   yhat = a + b * x
slr_predict <- function(fit, x) {
  fit$a + fit$b * x
}


## ---- 4.2 Model directions -------------------------------------------------------
# Regression is not symmetric: y ~ x minimises vertical errors in y, so the
# slope of y on x is not 1 / (slope of x on y). Only R^2 (= r^2) is the same
# in both directions. We therefore fix a "primary" direction per pair and
# report the reverse direction as well:
#   - meantemp is the response for every pair that contains it (it is the
#     variable of main interest and the best predicted one);
#   - for the other pairs we follow the chain large-scale state -> local
#     weather: meanpressure (synoptic pressure field) drives wind_speed, and
#     both wind (ventilation, advection of dry/moist air) and pressure act on
#     humidity. Hence humidity ~ wind_speed, humidity ~ meanpressure and
#     wind_speed ~ meanpressure.
#
# Evaluation: the headline numbers are the full-data fits (2013-01-01 to
# 2017-04-24). As an out-of-sample check each model is also fitted on the
# train period 2013-2016 and evaluated on the test period 2017 (Jan-Apr):
#   R^2_test = 1 - sum (y_t - yhat_t)^2 / sum (y_t - ybar_test)^2
# with yhat_t from the TRAIN coefficients. R^2_test < 0 means the line does
# worse than simply predicting the test-period mean.

LIN_PAIRS <- list(c(y = "meantemp",   x = "humidity"),
                  c(y = "meantemp",   x = "wind_speed"),
                  c(y = "meantemp",   x = "meanpressure"),
                  c(y = "humidity",   x = "wind_speed"),
                  c(y = "humidity",   x = "meanpressure"),
                  c(y = "wind_speed", x = "meanpressure"))

IS_TRAIN <- dat$date <= as.Date("2016-12-31")
IS_TEST  <- !IS_TRAIN
log_line("Train period: %s to %s (%d days); test period: %s to %s (%d days).",
         format(min(dat$date[IS_TRAIN])), format(max(dat$date[IS_TRAIN])), sum(IS_TRAIN),
         format(min(dat$date[IS_TEST])),  format(max(dat$date[IS_TEST])),  sum(IS_TEST))


## ---- 4.3 Fit all models ---------------------------------------------------------------

lin_rows  <- list()
lin_fits  <- list()   # full-data fits of the primary direction (for the figures)
lin_train <- list()   # train-period fits of the primary direction
for (k in seq_along(LIN_PAIRS)) {
  for (direction in c("primary", "reverse")) {
    yv <- if (direction == "primary") LIN_PAIRS[[k]]["y"] else LIN_PAIRS[[k]]["x"]
    xv <- if (direction == "primary") LIN_PAIRS[[k]]["x"] else LIN_PAIRS[[k]]["y"]
    x <- dat[[xv]]
    y <- dat[[yv]]
    f  <- my_slr(x, y)                         # full data (headline)
    ft <- my_slr(x[IS_TRAIN], y[IS_TRAIN])     # train period only
    yhat_test <- slr_predict(ft, x[IS_TEST])
    if (direction == "primary") {
      lin_fits[[k]]  <- f
      lin_train[[k]] <- ft
    }
    lin_rows[[length(lin_rows) + 1]] <- data.frame(
      Pair = k, Direction = direction, Model = paste(yv, "~", xv),
      Response = yv, Predictor = xv, n = f$n,
      a = f$a, SE_a = f$se_a, b = f$b, SE_b = f$se_b, t_b = f$t_b, p_b = f$p_b,
      CI_b_lower = f$ci_b[1], CI_b_upper = f$ci_b[2],
      SSE = f$SSE, SSR = f$SSR, SST = f$SST, R2 = f$R2, adj_R2 = f$adjR2,
      RMSE = f$RMSE, MAE = f$MAE, RSE = f$RSE, DW = f$DW,
      a_train = ft$a, b_train = ft$b, R2_train = ft$R2,
      R2_test = r2(y[IS_TEST], yhat_test), RMSE_test = rmse(y[IS_TEST], yhat_test),
      stringsAsFactors = FALSE, row.names = NULL)
  }
}
lin_tab <- do.call(rbind, lin_rows)

cat("Simple linear regression results (full data; train 2013-2016 -> test 2017)\n")
print(data.frame(Model = lin_tab$Model, a = signif(lin_tab$a, 5), b = signif(lin_tab$b, 4),
                 p_b = format_p(lin_tab$p_b), R2 = round(lin_tab$R2, 4),
                 RMSE = round(lin_tab$RMSE, 3), DW = round(lin_tab$DW, 3),
                 R2_train = round(lin_tab$R2_train, 4), R2_test = round(lin_tab$R2_test, 4)),
      row.names = FALSE)
cat("\n")

lin_tex <- data.frame(
  Model = ifelse(lin_tab$Direction == "reverse", paste(lin_tab$Model, "(rev.)"), lin_tab$Model),
  `a` = sprintf("%.3f", lin_tab$a), `b` = sprintf("%.4f", lin_tab$b),
  `95% CI of b` = sprintf("[%.4f, %.4f]", lin_tab$CI_b_lower, lin_tab$CI_b_upper),
  `t(b)` = sprintf("%.2f", lin_tab$t_b), p = format_p(lin_tab$p_b),
  `R2` = lin_tab$R2, `adj. R2` = lin_tab$adj_R2, RMSE = lin_tab$RMSE, DW = lin_tab$DW,
  `R2 train` = lin_tab$R2_train, `R2 test` = lin_tab$R2_test,
  check.names = FALSE, stringsAsFactors = FALSE)
save_table(lin_tab, "linear_results",
           caption = paste0("Simple linear regression $y = a + b\\,x$ for every pair, in the primary ",
                            "direction and reversed (rev.). $R^2$, adjusted $R^2$, RMSE and the ",
                            "Durbin--Watson statistic (DW) refer to the full data (2013--2017); ",
                            "$R^2$ train/test: fitted on 2013--2016, evaluated on 2017 (Jan--Apr)."),
           label = "tab:linear_results", tex_df = lin_tex)


## ---- 4.4 Figures ---------------------------------------------------------------------

# fmt_coef: coefficient with 4 significant digits for equations
fmt_coef <- function(v) formatC(v, digits = 4, format = "g")
# equation_text: fitted line as text, y = a + b x
equation_text <- function(f, yv, xv) {
  sprintf("%s = %s %s %s %s", yv, fmt_coef(f$a), if (f$b < 0) "-" else "+",
          fmt_coef(abs(f$b)), xv)
}

# Scatter plots with the full-data line (solid) and the train-period line
# (dashed); train days in blue, 2017 test days in orange.
save_figure("linear_fits", function() {
  op <- par(mfrow = c(2, 3), mar = c(4.5, 4.5, 4, 1), las = 1, cex.main = 0.95)
  on.exit(par(op))
  for (k in seq_along(LIN_PAIRS)) {
    yv <- LIN_PAIRS[[k]]["y"]; xv <- LIN_PAIRS[[k]]["x"]
    f <- lin_fits[[k]]; ft <- lin_train[[k]]
    plot(dat[[xv]][IS_TRAIN], dat[[yv]][IS_TRAIN], pch = 16, cex = 0.4,
         col = adjustcolor(COL_MAIN, alpha.f = 0.35), bty = "l",
         xlim = range(dat[[xv]]), ylim = range(dat[[yv]]),
         xlab = VAR_LABELS[xv], ylab = VAR_LABELS[yv],
         main = sprintf("%s\nR² = %.3f  (test 2017: %.3f)", equation_text(f, yv, xv),
                        f$R2, lin_tab$R2_test[lin_tab$Pair == k & lin_tab$Direction == "primary"]))
    points(dat[[xv]][IS_TEST], dat[[yv]][IS_TEST], pch = 16, cex = 0.5,
           col = adjustcolor(COL_HIGHLIGHT, alpha.f = 0.8))
    abline(a = f$a,  b = f$b,  lwd = 2.5, col = "black")
    abline(a = ft$a, b = ft$b, lwd = 1.5, lty = 2, col = "grey35")
    if (k == 1) {
      legend("topright", legend = c("2013-2016 (train)", "2017 (test)", "fit, all data",
                                    "fit, train only"),
             pch = c(16, 16, NA, NA), col = c(COL_MAIN, COL_HIGHLIGHT, "black", "grey35"),
             lty = c(NA, NA, 1, 2), lwd = c(NA, NA, 2.5, 1.5), bty = "n", cex = 0.8, bg = "white")
    }
  }
}, width = 13, height = 8.5)

# Residuals vs fitted values, with binned means (orange) as a hand-made
# smoother: a curved band of binned means means the straight line misses a
# non-linear pattern.
save_figure("linear_residuals_vs_fitted", function() {
  op <- par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1), las = 1, cex.main = 0.95)
  on.exit(par(op))
  for (k in seq_along(LIN_PAIRS)) {
    yv <- LIN_PAIRS[[k]]["y"]; xv <- LIN_PAIRS[[k]]["x"]
    f <- lin_fits[[k]]
    plot(f$fitted, f$resid, pch = 16, cex = 0.4, col = adjustcolor(COL_MAIN, alpha.f = 0.35),
         xlab = sprintf("Fitted %s", yv), ylab = "Residual", bty = "l",
         main = sprintf("%s ~ %s   (DW = %.2f)", yv, xv, f$DW))
    abline(h = 0, lty = 2, col = "grey30")
    bm <- binned_means(f$fitted, f$resid, k = 20)
    lines(bm$x, bm$y, col = COL_HIGHLIGHT, lwd = 2.5)
    points(bm$x, bm$y, pch = 21, bg = COL_HIGHLIGHT, col = "white", cex = 1.1)
    if (k == 1) {
      legend("topright", legend = "binned mean of residuals (20 groups)", lwd = 2.5,
             col = COL_HIGHLIGHT, pch = 21, pt.bg = COL_HIGHLIGHT, bty = "n", cex = 0.8)
    }
  }
}, width = 13, height = 8.5)

# Normal Q-Q plots of the residuals. Theoretical quantiles z_i = my_qnorm(p_i)
# at plotting positions p_i = (i - 0.5) / n; the reference line passes through
# the first and third quartiles (as in R's qqline), slope
# (Q3_e - Q1_e) / (z_0.75 - z_0.25).
QQ_Z <- my_qnorm((seq_len(n_obs) - 0.5) / n_obs)
Z_Q  <- my_qnorm(c(0.25, 0.75))
save_figure("linear_qq", function() {
  op <- par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1), las = 1, cex.main = 0.95)
  on.exit(par(op))
  for (k in seq_along(LIN_PAIRS)) {
    yv <- LIN_PAIRS[[k]]["y"]; xv <- LIN_PAIRS[[k]]["x"]
    e  <- sort(lin_fits[[k]]$resid)
    qe <- my_quantile(e, c(0.25, 0.75))
    slope <- (qe[2] - qe[1]) / (Z_Q[2] - Z_Q[1])
    plot(QQ_Z, e, pch = 16, cex = 0.4, col = adjustcolor(COL_MAIN, alpha.f = 0.5), bty = "l",
         xlab = "Theoretical normal quantile", ylab = "Sample quantile of residuals",
         main = sprintf("%s ~ %s  (skewness %.2f)", yv, xv, my_skewness(e)))
    abline(a = qe[1] - slope * Z_Q[1], b = slope, lwd = 2, col = "black")
  }
}, width = 13, height = 8.5)


## ---- 4.5 Summary --------------------------------------------------------------------

prim <- lin_tab[lin_tab$Direction == "primary", ]
prim <- prim[order(-prim$R2), ]
cat("\n")
log_line("Linear regression summary (primary directions, ranked by full-data R^2)")
for (i in seq_len(nrow(prim))) {
  log_line("  %-26s R^2 = %.3f: %s explains %4.1f%% of the variance of %s; slope %s per unit (95%% CI %s to %s), p %s.",
           prim$Model[i], prim$R2[i], prim$Predictor[i], 100 * prim$R2[i], prim$Response[i],
           fmt_coef(prim$b[i]), fmt_coef(prim$CI_b_lower[i]), fmt_coef(prim$CI_b_upper[i]),
           ifelse(prim$p_b[i] < 0.001, "< 0.001", sprintf("= %.3f", prim$p_b[i])))
}
log_line("Reverse direction: R^2 is identical (= r^2) but the slope is not the reciprocal;")
for (k in seq_along(LIN_PAIRS)) {
  rows <- lin_tab[lin_tab$Pair == k, ]
  log_line("  %-26s b = %9.4f | reverse b = %9.4f | product %.4f = r^2 = %.4f", rows$Model[1],
           rows$b[1], rows$b[2], rows$b[1] * rows$b[2], rows$R2[1])
}
log_line("Durbin-Watson: %s.", paste(sprintf("%s %.2f", prim$Model, prim$DW), collapse = ", "))
if (all(lin_tab$DW < 1.5)) {
  log_line("  All DW values are far below 2: the residuals are strongly positively autocorrelated")
  log_line("  (seasonal structure left in the residuals), so SE(b), t and p are optimistic. The")
  log_line("  slopes and R^2 remain valid descriptions of the fitted relation.")
}
log_line("Out-of-sample check (train 2013-2016 -> test 2017, Jan-Apr):")
for (i in seq_len(nrow(prim))) {
  log_line("  %-26s R^2 train %.3f -> test %6.3f (RMSE test %.3f vs full-data RMSE %.3f)",
           prim$Model[i], prim$R2_train[i], prim$R2_test[i], prim$RMSE_test[i], prim$RMSE[i])
}
log_line("  The test period covers only one season (winter to early summer), so its variance and")
log_line("  range are smaller than in the full data; a lower test R^2 is partly a consequence of")
log_line("  that narrower range, not only of a worse fit (compare the test RMSE, which is not")
log_line("  affected by the range of y).")

writeLines(LOG, file.path(TAB_DIR, "linear_summary.txt"))
cat(sprintf("  saved log    %s/linear_summary.txt\n", TAB_DIR))


# =============================================================================
# 5  NON-LINEAR REGRESSION
# =============================================================================
# Same 6 pairs and directions as section 4. Candidate curves (one predictor):
#   linear       y = a + b x                       (reference from section 4)
#   quadratic    y = c0 + c1 u + c2 u^2            u = (x - xbar) / s_x
#   cubic        y = c0 + c1 u + c2 u^2 + c3 u^3
#   exponential  y = a * exp(b x)
#   logarithmic  y = a + b ln(x)                   (needs x > 0)
#   power        y = a * x^b                       (needs x > 0)
# R^2 is ALWAYS computed on the original y scale, R^2 = 1 - SSE/SST with
# SSE = sum (y - yhat)^2, so all curves are compared on the same footing
# (the R^2 of a linearised fit such as ln y on x would not be comparable).

cat("\n================ 5  NON-LINEAR REGRESSION ================\n\n")
LOG <- character(0)   # new log: non-linear regression summary


## ---- 5.1 Model fitting functions -------------------------------------------------
# Every fitter returns the same structure: model, k (number of parameters),
# valid, converged, iterations, status (text), params, center, scale,
# equation (text) and curve(x_new), the fitted curve as a function.

MODEL_NAMES  <- c("linear", "poly2", "poly3", "exponential", "logarithmic", "power")
MODEL_K      <- c(linear = 2, poly2 = 3, poly3 = 4, exponential = 2, logarithmic = 2, power = 2)
MODEL_LABELS <- c(linear = "linear", poly2 = "quadratic", poly3 = "cubic",
                  exponential = "exponential", logarithmic = "logarithmic", power = "power")

# term: " + 1.23 u^2" style text for equations
term <- function(v, s) {
  sprintf(" %s %s%s", if (v < 0) "-" else "+", fmt_coef(abs(v)), s)
}

# invalid_fit: placeholder result for a model that is undefined or failed
# (keeps the reason in status, so the model is reported but not used)
invalid_fit <- function(model, status) {
  list(model = model, k = MODEL_K[[model]], valid = FALSE, converged = FALSE,
       iterations = NA, status = status, params = numeric(0), center = NA, scale = NA,
       equation = "", curve = NULL)
}

# Linear reference: closed-form OLS from section 4.
fit_linear_model <- function(x, y) {
  f <- my_slr(x, y)
  list(model = "linear", k = 2, valid = TRUE, converged = TRUE, iterations = 0,
       status = "ok (closed form)", params = c(a = f$a, b = f$b), center = NA, scale = NA,
       equation = paste0("y = ", fmt_coef(f$a), term(f$b, " x")),
       curve = function(xn) f$a + f$b * xn)
}

# poly_design: design matrix of a degree-d polynomial in the standardised predictor
#   u = (x - m) / s,   X = [1, u, u^2, ..., u^d]   (n x (d + 1))
poly_design <- function(x, m, s, d) {
  u <- (x - m) / s
  X <- matrix(1, nrow = length(x), ncol = d + 1)
  for (j in 1:d) X[, j + 1] <- u^j
  X
}

# fit_poly_model: polynomial least squares through the normal equations
#   minimise ||y - X beta||^2   =>   (X'X) beta = X'y
#   solved with our solve_linear_system (Gaussian elimination, partial pivoting).
#   Standardising x keeps X'X well conditioned: with raw pressure (~1000 hPa)
#   the columns would be ~1, 1e3, 1e6, 1e9 and X'X numerically singular.
#   The fitted curve is the same as with raw x; the coefficients refer to u.
fit_poly_model <- function(x, y, d) {
  m <- my_mean(x)
  s <- my_sd(x)
  X <- poly_design(x, m, s, d)
  beta <- solve_linear_system(t(X) %*% X, t(X) %*% y)
  names(beta) <- paste0("c", 0:d)
  eq <- paste0("y = ", fmt_coef(beta[1]),
               paste(sapply(1:d, function(j) term(beta[j + 1], if (j == 1) " u" else paste0(" u^", j))),
                     collapse = ""),
               sprintf(", u = (x - %s) / %s", fmt_coef(m), fmt_coef(s)))
  list(model = paste0("poly", d), k = d + 1, valid = TRUE, converged = TRUE, iterations = 0,
       status = "ok (normal equations)", params = beta, center = m, scale = s, equation = eq,
       curve = function(xn) as.vector(poly_design(xn, m, s, d) %*% beta))
}

# fit_exp_model: y = a * exp(b x)
#   Computed in the equivalent centred form y = A exp(b (x - c)), c = xbar,
#   a = A exp(-b c), so that A and b have similar magnitudes (with x ~ 1000 hPa,
#   a itself can be ~1e15, which makes Gauss-Newton badly conditioned).
#   1) Start values by linearising: ln y = ln A + b (x - c)  -> our my_slr:
#      A0 = exp(intercept), b0 = slope. Requires y > 0; if some y = 0 (wind
#      speed), the start uses only the days with y > 0.
#   2) Refine on the ORIGINAL scale with Gauss-Newton, Jacobian
#      d f / d A = exp(b (x - c)),   d f / d b = A (x - c) exp(b (x - c)).
fit_exp_model <- function(x, y) {
  pos <- y > 0
  if (sum(pos) < 3) return(invalid_fit("exponential", "skipped: fewer than 3 values with y > 0"))
  c0 <- my_mean(x)
  start <- my_slr(x[pos] - c0, log(y[pos]))
  # model f(x) = A exp(b (x - c)) and its Jacobian [df/dA, df/db]
  f   <- function(x, th) th[1] * exp(th[2] * (x - c0))
  jac <- function(x, th) {
    e <- exp(th[2] * (x - c0))
    cbind(e, th[1] * (x - c0) * e)
  }
  gn <- gauss_newton(f, jac, c(exp(start$a), start$b), x, y)
  A <- gn$theta[1]
  b <- gn$theta[2]
  a <- A * exp(-b * c0)
  note <- if (all(pos)) "" else sprintf("; start values from the %d values with y > 0", sum(pos))
  list(model = "exponential", k = 2, valid = TRUE, converged = gn$converged,
       iterations = gn$iterations,
       status = sprintf("%s (Gauss-Newton, %d iterations)%s",
                        if (gn$converged) "ok" else gn$message, gn$iterations, note),
       params = c(a = a, b = b), center = c0, scale = NA,
       equation = sprintf("y = %s * exp(%s x)", fmt_coef(a), fmt_coef(b)),
       curve = function(xn) A * exp(b * (xn - c0)))
}

# fit_log_model: y = a + b ln(x), an ordinary linear regression of y on ln(x).
#   Only defined for x > 0.
fit_log_model <- function(x, y) {
  if (any(x <= 0)) {
    return(invalid_fit("logarithmic",
                       sprintf("skipped: ln(x) undefined, %d predictor values <= 0", sum(x <= 0))))
  }
  f <- my_slr(log(x), y)
  list(model = "logarithmic", k = 2, valid = TRUE, converged = TRUE, iterations = 0,
       status = "ok (linear in ln x)", params = c(a = f$a, b = f$b), center = NA, scale = NA,
       equation = paste0("y = ", fmt_coef(f$a), term(f$b, " ln(x)")),
       curve = function(xn) f$a + f$b * log(xn))
}

# fit_power_model: y = a * x^b   (x > 0)
#   Computed in the equivalent form y = A (x / c)^b with c = geometric mean
#   of x, a = A c^(-b) (keeps A and b well scaled).
#   1) Start values by linearising: ln y = ln A + b ln(x / c)  -> my_slr
#      (days with y > 0 only).
#   2) Gauss-Newton on the original scale, Jacobian
#      d f / d A = (x / c)^b,   d f / d b = A (x / c)^b ln(x / c).
fit_power_model <- function(x, y) {
  if (any(x <= 0)) {
    return(invalid_fit("power",
                       sprintf("skipped: ln(x) undefined, %d predictor values <= 0", sum(x <= 0))))
  }
  pos <- y > 0
  if (sum(pos) < 3) return(invalid_fit("power", "skipped: fewer than 3 values with y > 0"))
  c0 <- exp(my_mean(log(x)))
  start <- my_slr(log(x[pos] / c0), log(y[pos]))
  # model f(x) = A (x / c)^b and its Jacobian [df/dA, df/db]
  f   <- function(x, th) th[1] * (x / c0)^th[2]
  jac <- function(x, th) {
    p <- (x / c0)^th[2]
    cbind(p, th[1] * p * log(x / c0))
  }
  gn <- gauss_newton(f, jac, c(exp(start$a), start$b), x, y)
  A <- gn$theta[1]
  b <- gn$theta[2]
  a <- A * c0^(-b)
  note <- if (all(pos)) "" else sprintf("; start values from the %d values with y > 0", sum(pos))
  list(model = "power", k = 2, valid = TRUE, converged = gn$converged,
       iterations = gn$iterations,
       status = sprintf("%s (Gauss-Newton, %d iterations)%s",
                        if (gn$converged) "ok" else gn$message, gn$iterations, note),
       params = c(a = a, b = b), center = c0, scale = NA,
       equation = sprintf("y = %s * x^%s", fmt_coef(a), fmt_coef(b)),
       curve = function(xn) A * (xn / c0)^b)
}

# safe_fit: run a fitter; any numerical error becomes an invalid fit with the
# error message as status instead of stopping the script.
safe_fit <- function(model, fitter, x, y) {
  tryCatch(fitter(x, y),
           error = function(e) invalid_fit(model, paste("failed:", conditionMessage(e))))
}

# fit_all_models: fit all six candidate curves to (x, y)
fit_all_models <- function(x, y) {
  list(linear      = safe_fit("linear",      fit_linear_model, x, y),
       poly2       = safe_fit("poly2",       function(x, y) fit_poly_model(x, y, 2), x, y),
       poly3       = safe_fit("poly3",       function(x, y) fit_poly_model(x, y, 3), x, y),
       exponential = safe_fit("exponential", fit_exp_model, x, y),
       logarithmic = safe_fit("logarithmic", fit_log_model, x, y),
       power       = safe_fit("power",       fit_power_model, x, y))
}


## ---- 5.2 Fit and evaluate every model -------------------------------------------
# For each pair, direction and curve:
#   full data   -> R^2, adjusted R^2 (p = k - 1 predictor terms), RMSE, MAE, DW
#   train 2013-2016 -> R^2_train; prediction of 2017 -> R^2_test, RMSE_test
# Model choice (per pair and direction), to avoid rewarding a polynomial
# simply for having more parameters:
#   score = (adjusted R^2 on the full data + R^2 on the 2017 test period) / 2
#   adjusted R^2 penalises extra parameters in-sample; R^2_test checks that
#   the curve still works on data not used for fitting. Only valid fits whose
#   Gauss-Newton iterations converged (full and train) are eligible.

nl_rows <- list()
NL_FITS <- list()   # full-data fits of the primary direction (for figures)
for (k in seq_along(LIN_PAIRS)) {
  for (direction in c("primary", "reverse")) {
    yv <- if (direction == "primary") LIN_PAIRS[[k]]["y"] else LIN_PAIRS[[k]]["x"]
    xv <- if (direction == "primary") LIN_PAIRS[[k]]["x"] else LIN_PAIRS[[k]]["y"]
    x <- dat[[xv]]
    y <- dat[[yv]]
    fits_full  <- fit_all_models(x, y)
    fits_train <- fit_all_models(x[IS_TRAIN], y[IS_TRAIN])
    if (direction == "primary") NL_FITS[[k]] <- fits_full
    for (m in MODEL_NAMES) {
      ff <- fits_full[[m]]
      ft <- fits_train[[m]]
      R2 <- adjR2 <- RMSE <- MAE <- DW <- R2_train <- R2_test <- RMSE_test <- NA
      if (ff$valid) {
        yhat  <- ff$curve(x)
        e     <- y - yhat
        R2    <- r2(y, yhat)
        adjR2 <- adj_r2(y, yhat, length(y), p = ff$k - 1)
        RMSE  <- rmse(y, yhat)
        MAE   <- mae(y, yhat)
        DW    <- sum((e[-1] - e[-length(e)])^2) / sum(e^2)
      }
      if (ft$valid) {
        R2_train <- r2(y[IS_TRAIN], ft$curve(x[IS_TRAIN]))
        yhat_te  <- ft$curve(x[IS_TEST])
        if (all(is.finite(yhat_te))) {
          R2_test   <- r2(y[IS_TEST], yhat_te)
          RMSE_test <- rmse(y[IS_TEST], yhat_te)
        }
      }
      status <- ff$status
      if (ff$valid && !(ft$valid && ft$converged)) status <- paste0(status, "; train fit: ", ft$status)
      eligible <- ff$valid && ff$converged && ft$valid && ft$converged && is.finite(R2_test)
      p <- c(ff$params, rep(NA, 4))[1:4]
      nl_rows[[length(nl_rows) + 1]] <- data.frame(
        Pair = k, Direction = direction, Relation = paste(yv, "~", xv),
        Response = yv, Predictor = xv, Model = m, Curve = MODEL_LABELS[[m]], k = MODEL_K[[m]],
        Valid = ff$valid, Converged = ff$converged, Iterations = ff$iterations,
        Status = status, Equation = ff$equation,
        p0 = p[1], p1 = p[2], p2 = p[3], p3 = p[4], center = ff$center, scale = ff$scale,
        R2 = R2, adj_R2 = adjR2, RMSE = RMSE, MAE = MAE, DW = DW,
        R2_train = R2_train, R2_test = R2_test, RMSE_test = RMSE_test,
        Eligible = eligible, Score = if (eligible) (adjR2 + R2_test) / 2 else NA,
        stringsAsFactors = FALSE, row.names = NULL)
    }
  }
}
nl_tab <- do.call(rbind, nl_rows)

# best model and gain over the linear reference, per pair and direction
nl_tab$Best <- FALSE
nl_tab$R2_gain_vs_linear <- NA
for (k in seq_along(LIN_PAIRS)) {
  for (direction in c("primary", "reverse")) {
    idx  <- which(nl_tab$Pair == k & nl_tab$Direction == direction)
    lin  <- idx[nl_tab$Model[idx] == "linear"]
    nl_tab$R2_gain_vs_linear[idx] <- nl_tab$R2[idx] - nl_tab$R2[lin]
    cand <- idx[nl_tab$Eligible[idx]]
    nl_tab$Best[cand[which.max(nl_tab$Score[cand])]] <- TRUE
  }
}

cat("Non-linear regression results (primary directions; R^2 on the original y scale)\n")
show <- nl_tab[nl_tab$Direction == "primary", ]
print(data.frame(Relation = show$Relation, Curve = show$Curve,
                 R2 = round(show$R2, 4), adj_R2 = round(show$adj_R2, 4),
                 RMSE = round(show$RMSE, 3), R2_test = round(show$R2_test, 4),
                 Score = round(show$Score, 4), Best = ifelse(show$Best, "*", ""),
                 Status = substr(show$Status, 1, 45)), row.names = FALSE)
cat("\n")


## ---- 5.3 Tables -------------------------------------------------------------------------

# note_text: short remark for the LaTeX table (skipped / not converged)
note_text <- function(row) {
  if (!row$Valid) return(sub("^skipped: ", "", sub(", [0-9]+ predictor values <= 0", " (x <= 0)", row$Status)))
  if (!row$Converged) return("not converged")
  ""
}
nl_tex_rows <- nl_tab[nl_tab$Direction == "primary", ]
nl_tex <- data.frame(
  Relation = ifelse(nl_tex_rows$Model == "linear", nl_tex_rows$Relation, ""),
  Curve = paste0(nl_tex_rows$Curve, ifelse(nl_tex_rows$Best, " *", "")),
  k = nl_tex_rows$k, R2 = nl_tex_rows$R2, `adj. R2` = nl_tex_rows$adj_R2,
  RMSE = nl_tex_rows$RMSE, MAE = nl_tex_rows$MAE, `R2 test` = nl_tex_rows$R2_test,
  Score = nl_tex_rows$Score,
  Note = sapply(seq_len(nrow(nl_tex_rows)), function(i) note_text(nl_tex_rows[i, ])),
  check.names = FALSE, stringsAsFactors = FALSE)
save_table(nl_tab, "nonlinear_results",
           caption = paste0("Linear and non-linear simple regression for the six pairs (primary ",
                            "direction; the reverse direction is in the CSV). $R^2$ is computed on the ",
                            "original $y$ scale for every curve. Score = (adjusted $R^2$ + test $R^2$)/2; ",
                            "* marks the selected curve. $k$ = number of parameters."),
           label = "tab:nonlinear_results", tex_df = nl_tex)

best_tab <- do.call(rbind, lapply(split(nl_tab, paste(nl_tab$Pair, nl_tab$Direction)), function(d) {
  b <- d[d$Best, ]
  l <- d[d$Model == "linear", ]
  data.frame(Pair = b$Pair, Direction = b$Direction, Relation = b$Relation, Best = b$Curve,
             Equation = b$Equation, R2_linear = l$R2, R2_best = b$R2, Gain = b$R2 - l$R2,
             adj_R2_best = b$adj_R2, R2_test_linear = l$R2_test, R2_test_best = b$R2_test,
             RMSE_linear = l$RMSE, RMSE_best = b$RMSE, stringsAsFactors = FALSE)
}))
best_tab <- best_tab[order(best_tab$Pair, best_tab$Direction != "primary"), ]
rownames(best_tab) <- NULL
best_tex <- data.frame(
  Relation = paste0(best_tab$Relation, ifelse(best_tab$Direction == "reverse", " (rev.)", "")),
  `Best curve` = best_tab$Best, `R2 linear` = best_tab$R2_linear, `R2 best` = best_tab$R2_best,
  Gain = best_tab$Gain, `R2 test linear` = best_tab$R2_test_linear,
  `R2 test best` = best_tab$R2_test_best, `RMSE linear` = best_tab$RMSE_linear,
  `RMSE best` = best_tab$RMSE_best, check.names = FALSE, stringsAsFactors = FALSE)
save_table(best_tab, "nonlinear_best",
           caption = paste0("Selected curve per relation (both directions) compared with the linear ",
                            "fit: full-data $R^2$, gain in $R^2$, and $R^2$ on the 2017 test period."),
           label = "tab:nonlinear_best", tex_df = best_tex)


## ---- 5.4 Figures ------------------------------------------------------------------------

# Colour by curve (fixed order) plus a line type as a second cue.
CURVE_COL <- c(linear = "black", poly2 = "#2a78d6", poly3 = "#eb6834",
               exponential = "#1baf7a", logarithmic = "#eda100", power = "#4a3aa7")
CURVE_LTY <- c(linear = 2, poly2 = 1, poly3 = 1, exponential = 1, logarithmic = 4, power = 5)

# legend_corner: the plot corner with the fewest data points (for the legend)
legend_corner <- function(x, y) {
  u <- (x - min(x)) / (max(x) - min(x))
  v <- (y - min(y)) / (max(y) - min(y))
  counts <- c(topleft     = sum(u < 0.45 & v > 0.6),  topright    = sum(u > 0.55 & v > 0.6),
              bottomleft  = sum(u < 0.45 & v < 0.4),  bottomright = sum(u > 0.55 & v < 0.4))
  names(counts)[which.min(counts)]
}

save_figure("nonlinear_fits", function() {
  op <- par(mfrow = c(2, 3), mar = c(4.5, 4.5, 4, 1), las = 1, cex.main = 0.95)
  on.exit(par(op))
  for (k in seq_along(LIN_PAIRS)) {
    yv <- LIN_PAIRS[[k]]["y"]; xv <- LIN_PAIRS[[k]]["x"]
    x <- dat[[xv]]; y <- dat[[yv]]
    fits <- NL_FITS[[k]]
    rows <- nl_tab[nl_tab$Pair == k & nl_tab$Direction == "primary", ]
    best <- rows$Model[rows$Best]
    plot(x, y, pch = 16, cex = 0.35, col = adjustcolor("grey45", alpha.f = 0.35), bty = "l",
         xlab = VAR_LABELS[xv], ylab = VAR_LABELS[yv],
         main = sprintf("%s ~ %s\nbest: %s (R² %.3f, test %.3f)", yv, xv, MODEL_LABELS[[best]],
                        rows$R2[rows$Best], rows$R2_test[rows$Best]))
    grid <- seq(min(x), max(x), length.out = 400)
    for (m in c(setdiff(MODEL_NAMES, best), best)) {    # best curve drawn last (on top)
      if (fits[[m]]$valid) {
        lines(grid, fits[[m]]$curve(grid), col = CURVE_COL[[m]], lty = CURVE_LTY[[m]],
              lwd = if (m == best) 3.5 else 1.8)
      }
    }
    labs <- sapply(MODEL_NAMES, function(m) {
      r <- rows[rows$Model == m, ]
      if (!r$Valid) sprintf("%s: skipped (x <= 0)", MODEL_LABELS[[m]])
      else sprintf("%s  R² = %.3f%s", MODEL_LABELS[[m]], r$R2, if (m == best) "  (best)" else "")
    })
    valid <- sapply(MODEL_NAMES, function(m) fits[[m]]$valid)
    legend(legend_corner(x, y), legend = labs, col = ifelse(valid, CURVE_COL[MODEL_NAMES], "grey60"),
           lty = ifelse(valid, CURVE_LTY[MODEL_NAMES], 0), lwd = ifelse(MODEL_NAMES == best, 3.5, 1.8),
           text.col = ifelse(valid, "black", "grey50"), bg = adjustcolor("white", alpha.f = 0.85),
           box.col = NA, cex = 0.72)
  }
}, width = 13, height = 8.5)

# Residuals of the selected curve vs fitted values, with binned means.
save_figure("nonlinear_best_residuals", function() {
  op <- par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1), las = 1, cex.main = 0.95)
  on.exit(par(op))
  for (k in seq_along(LIN_PAIRS)) {
    yv <- LIN_PAIRS[[k]]["y"]; xv <- LIN_PAIRS[[k]]["x"]
    rows <- nl_tab[nl_tab$Pair == k & nl_tab$Direction == "primary", ]
    best <- rows$Model[rows$Best]
    yhat <- NL_FITS[[k]][[best]]$curve(dat[[xv]])
    e <- dat[[yv]] - yhat
    plot(yhat, e, pch = 16, cex = 0.4, col = adjustcolor(COL_MAIN, alpha.f = 0.35), bty = "l",
         xlab = sprintf("Fitted %s", yv), ylab = "Residual",
         main = sprintf("%s ~ %s, %s  (DW = %.2f)", yv, xv, MODEL_LABELS[[best]], rows$DW[rows$Best]))
    abline(h = 0, lty = 2, col = "grey30")
    bm <- binned_means(yhat, e, k = 20)
    lines(bm$x, bm$y, col = COL_HIGHLIGHT, lwd = 2.5)
    points(bm$x, bm$y, pch = 21, bg = COL_HIGHLIGHT, col = "white", cex = 1.1)
    if (k == 1) {
      legend("topright", legend = "binned mean of residuals (20 groups)", lwd = 2.5,
             col = COL_HIGHLIGHT, pch = 21, pt.bg = COL_HIGHLIGHT, bty = "n", cex = 0.8)
    }
  }
}, width = 13, height = 8.5)


## ---- 5.5 Summary ---------------------------------------------------------------------

# verdict: verbal judgement of the R^2 gain of the selected curve over the line
verdict <- function(gain, best) {
  if (best == "linear") "linear fit is best"
  else if (gain <= 0) "no in-sample gain; chosen only for a slightly better 2017 test R^2 -> essentially linear"
  else if (gain < 0.01) "essentially linear (gain < 0.01)"
  else if (gain < 0.05) "mild curvature (gain 0.01-0.05)"
  else "clearly non-linear (gain >= 0.05)"
}

cat("\n")
log_line("Non-linear regression summary (R^2 on the original scale; selection by")
log_line("score = (adjusted R^2 + test R^2) / 2)")
for (dirn in c("primary", "reverse")) {
  log_line("%s direction:", if (dirn == "primary") "Primary" else "Reverse")
  bt <- best_tab[best_tab$Direction == dirn, ]
  for (i in seq_len(nrow(bt))) {
    log_line("  %-26s best %-11s R^2 %.3f -> %.3f (gain %+.3f), test R^2 %6.3f -> %6.3f: %s",
             bt$Relation[i], bt$Best[i], bt$R2_linear[i], bt$R2_best[i], bt$Gain[i],
             bt$R2_test_linear[i], bt$R2_test_best[i], verdict(bt$Gain[i], bt$Best[i]))
  }
}

# Overfitting check: where the cubic beats the quadratic in-sample but not out of sample.
over <- c()
for (k in seq_along(LIN_PAIRS)) for (dirn in c("primary", "reverse")) {
  d  <- nl_tab[nl_tab$Pair == k & nl_tab$Direction == dirn, ]
  p2 <- d[d$Model == "poly2", ]; p3 <- d[d$Model == "poly3", ]
  if (p3$adj_R2 > p2$adj_R2 && p3$R2_test < p2$R2_test) over <- c(over, d$Relation[1])
}
log_line("Overfitting check: the cubic has a higher adjusted R^2 than the quadratic but a lower")
log_line("  test R^2 for %d of 12 relations%s.", length(over),
         if (length(over) > 0) paste0(" (", paste(over, collapse = ", "), ")") else "")

skipped <- nl_tab[!nl_tab$Valid, ]
log_line("Skipped models: %d (%s).", nrow(skipped),
         if (nrow(skipped) > 0) paste(unique(sprintf("%s with predictor %s: %s", skipped$Curve,
                                                       skipped$Predictor, sub("^skipped: ", "", skipped$Status))),
                                      collapse = "; ") else "none")
gn_rows <- nl_tab[nl_tab$Model %in% c("exponential", "power") & nl_tab$Valid, ]
log_line("Gauss-Newton fits (exponential, power): %d of %d converged on the full data (%s iterations).",
         sum(gn_rows$Converged), nrow(gn_rows),
         paste(range(gn_rows$Iterations), collapse = "-"))
zero_start <- unique(nl_tab$Relation[grepl("start values from", nl_tab$Status)])
if (length(zero_start) > 0) {
  log_line("  Start values from y > 0 only (y has zeros) for: %s.", paste(zero_start, collapse = ", "))
}

writeLines(LOG, file.path(TAB_DIR, "nonlinear_summary.txt"))
cat(sprintf("  saved log    %s/nonlinear_summary.txt\n", TAB_DIR))


# =============================================================================
# 6  R^2 COMPARISON & CONCLUSIONS
# =============================================================================
# Merges sections 4 and 5, adds robustness checks (influential points, wind
# capping, seasonality) and writes an auto-generated draft of the conclusions.

cat("\n================ 6  R^2 COMPARISON & CONCLUSIONS ================\n\n")
LOG <- character(0)   # new log: key findings / conclusions draft

# Strength of a relation by (full-data) R^2 of the selected curve:
#   >= 0.50 strong, 0.25-0.50 moderate, 0.10-0.25 weak, < 0.10 very weak
r2_strength <- function(R2) {
  if (R2 >= 0.5) "strong" else if (R2 >= 0.25) "moderate" else if (R2 >= 0.1) "weak" else "very weak"
}
GAIN_MEANINGFUL <- 0.05   # Delta R^2 above which a non-linear curve counts as a real improvement


## ---- 6.1 R^2 comparison table ------------------------------------------------------
# For every relation (both directions): R^2 of every curve, the selected
# curve (section 5 score), the best NON-linear curve by the same score, and
#   Delta R^2 = R^2(best non-linear) - R^2(linear)
# plus the largest in-sample gain of any non-linear curve (to show how much
# of the "improvement" is only in-sample).

comp_rows <- list()
for (k in seq_along(LIN_PAIRS)) {
  for (direction in c("primary", "reverse")) {
    d   <- nl_tab[nl_tab$Pair == k & nl_tab$Direction == direction, ]
    lin <- d[d$Model == "linear", ]
    nl  <- d[d$Model != "linear" & d$Eligible, ]
    bnl <- nl[which.max(nl$Score), ]
    mx  <- nl[which.max(nl$R2), ]
    sel <- d[d$Best, ]
    # R2_of: full-data R^2 of curve m for this relation
    R2_of <- function(m) d$R2[d$Model == m]
    comp_rows[[length(comp_rows) + 1]] <- data.frame(
      Pair = k, Direction = direction, Relation = lin$Relation,
      Response = lin$Response, Predictor = lin$Predictor,
      r = my_pearson(dat[[lin$Predictor]], dat[[lin$Response]]),
      Linear = R2_of("linear"), Quadratic = R2_of("poly2"), Cubic = R2_of("poly3"),
      Exponential = R2_of("exponential"), Logarithmic = R2_of("logarithmic"), Power = R2_of("power"),
      Selected = sel$Curve, R2_selected = sel$R2, R2_test_selected = sel$R2_test,
      Best_nonlinear = bnl$Curve, R2_best_nonlinear = bnl$R2, Delta_R2 = bnl$R2 - lin$R2,
      R2_test_linear = lin$R2_test, R2_test_best_nonlinear = bnl$R2_test,
      Max_gain_curve = mx$Curve, Max_gain_in_sample = mx$R2 - lin$R2,
      DW_linear = lin$DW, stringsAsFactors = FALSE, row.names = NULL)
  }
}
comp <- do.call(rbind, comp_rows)
comp$Strength <- sapply(comp$R2_selected, r2_strength)

comp_tex <- data.frame(
  Relation = paste0(comp$Relation, ifelse(comp$Direction == "reverse", " (rev.)", "")),
  Linear = comp$Linear, Quadr. = comp$Quadratic, Cubic = comp$Cubic, Exp. = comp$Exponential,
  Log. = comp$Logarithmic, Power = comp$Power, Selected = comp$Selected,
  `Delta R2` = sprintf("%+.3f", comp$Delta_R2), `R2 test sel.` = comp$R2_test_selected,
  check.names = FALSE, stringsAsFactors = FALSE)
save_table(comp, "r2_comparison",
           caption = paste0("$R^2$ (original scale, full data) of every curve for every relation. ",
                            "Selected = curve with the highest score (adjusted $R^2$ + test $R^2$)/2; ",
                            "$\\Delta R^2$ = best non-linear curve minus linear; -- = not defined ",
                            "($\\ln x$ with $x \\le 0$)."),
           label = "tab:r2_comparison", tex_df = comp_tex)


## ---- 6.2 Robustness: influential points and wind capping --------------------------------

# cooks_distance: influence of each observation on a simple linear fit
#   leverage    h_i = 1/n + (x_i - xbar)^2 / Sxx
#   Cook's D    D_i = e_i^2 / (p * s^2) * h_i / (1 - h_i)^2,   p = 2, s^2 = SSE / (n - 2)
#   D_i > 4/n is a common flag for an influential point.
cooks_distance <- function(x, y) {
  f <- my_slr(x, y)
  n <- length(x)
  h <- 1 / n + (x - my_mean(x))^2 / sum((x - my_mean(x))^2)
  f$resid^2 / (2 * f$RSE^2) * h / (1 - h)^2
}

robust_rows <- list()
for (k in seq_along(LIN_PAIRS)) {
  yv <- LIN_PAIRS[[k]][["y"]]; xv <- LIN_PAIRS[[k]][["x"]]
  x <- dat[[xv]]; y <- dat[[yv]]
  D    <- cooks_distance(x, y)
  infl <- D > 4 / length(x)
  R2_all  <- my_slr(x, y)$R2
  R2_wo   <- my_slr(x[!infl], y[!infl])$R2
  # wind capping: same model with the uncapped wind speed
  R2_raw_wind <- NA
  if ("wind_speed" %in% c(xv, yv)) {
    xr <- if (xv == "wind_speed") dat$wind_speed_raw else x
    yr <- if (yv == "wind_speed") dat$wind_speed_raw else y
    R2_raw_wind <- my_slr(xr, yr)$R2
  }
  robust_rows[[k]] <- data.frame(Relation = paste(yv, "~", xv), R2_linear = R2_all,
                                 n_influential = sum(infl), max_Cooks_D = max(D),
                                 R2_without_influential = R2_wo, Change = R2_wo - R2_all,
                                 R2_uncapped_wind = R2_raw_wind, stringsAsFactors = FALSE)
}
robust <- do.call(rbind, robust_rows)
save_table(robust, "robustness_check",
           caption = paste0("Robustness of the linear fits: number of influential days (Cook's ",
                            "$D > 4/n$), $R^2$ after removing them, and $R^2$ with the uncapped wind speed."),
           label = "tab:robustness_check")


## ---- 6.3 Seasonality check --------------------------------------------------------------
# anomaly: x_anom = x - (mean of x in the same calendar month, over all years)
#   Subtracting the monthly climatology removes the common annual cycle; the
#   correlation of the anomalies measures the day-to-day (weather) association
#   that remains. R^2_month = 1 - SSE/SST with the monthly means as predictions
#   is the share of each variable's variance that is just the seasonal cycle.
monthly_climatology <- function(v) {
  sapply(1:12, function(m) my_mean(dat[[v]][dat$month == m]))[dat$month]
}
R2_MONTH <- sapply(VARS, function(v) r2(dat[[v]], monthly_climatology(v)))
SEASONS  <- c("Winter", "Summer", "Monsoon", "Post-monsoon")

seas_rows <- list()
for (k in seq_along(LIN_PAIRS)) {
  yv <- LIN_PAIRS[[k]][["y"]]; xv <- LIN_PAIRS[[k]][["x"]]
  ya <- dat[[yv]] - monthly_climatology(yv)
  xa <- dat[[xv]] - monthly_climatology(xv)
  r_all  <- my_pearson(dat[[xv]], dat[[yv]])
  r_anom <- my_pearson(xa, ya)
  within <- sapply(SEASONS, function(s) my_pearson(dat[[xv]][dat$season == s], dat[[yv]][dat$season == s]))
  seas_rows[[k]] <- data.frame(Relation = paste(yv, "~", xv), r_all = r_all, R2_all = r_all^2,
                               r_anomaly = r_anom, R2_anomaly = r_anom^2,
                               Seasonal_share = 1 - r_anom^2 / r_all^2,
                               r_Winter = within[1], r_Summer = within[2],
                               r_Monsoon = within[3], r_Postmonsoon = within[4],
                               stringsAsFactors = FALSE, row.names = NULL)
}
seas <- do.call(rbind, seas_rows)
save_table(seas, "seasonality_check",
           caption = paste0("Seasonality check: correlation of the raw daily values and of the anomalies ",
                            "from the monthly mean (annual cycle removed), the share of $R^2$ lost by ",
                            "removing the annual cycle, and the correlation within each season."),
           label = "tab:seasonality_check")


## ---- 6.4 Figures --------------------------------------------------------------------------

R2_PAL <- colorRampPalette(c("#f3f7fd", "#9ec3ef", "#2a78d6", "#0b3d7a"))(101)

# draw_r2_heatmap: R^2 matrix (rows = response, columns = predictor) with values in the cells
draw_r2_heatmap <- function(M, main) {
  p <- nrow(M)
  Z <- M
  Z[is.na(Z)] <- 0
  image(1:p, 1:p, t(Z[p:1, ]), zlim = c(0, 1), col = R2_PAL, axes = FALSE,
        xlab = "predictor", ylab = "", main = main)
  axis(1, at = 1:p, labels = colnames(M), tick = FALSE, cex.axis = 0.9)
  axis(2, at = 1:p, labels = rev(rownames(M)), tick = FALSE, las = 1, cex.axis = 0.9)
  for (i in 1:p) for (j in 1:p) {
    if (i == j) {
      rect(j - 0.5, p - i + 0.5, j + 0.5, p - i + 1.5, col = "grey88", border = NA)
    } else {
      text(j, p - i + 1, sprintf("%.3f", M[i, j]), font = 2, cex = 1.1,
           col = if (M[i, j] > 0.45) "white" else "black")
    }
  }
  abline(h = 0.5 + 0:p, v = 0.5 + 0:p, col = "white", lwd = 2)
}

R2_LIN_MAT  <- P_MAT^2
diag(R2_LIN_MAT) <- NA
R2_BEST_MAT <- matrix(NA, length(VARS), length(VARS), dimnames = list(VARS, VARS))
for (i in seq_len(nrow(comp))) R2_BEST_MAT[comp$Response[i], comp$Predictor[i]] <- comp$R2_selected[i]

save_figure("r2_heatmaps", function() {
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 0.22))
  op <- par(mar = c(4, 7.5, 3.5, 1))
  on.exit({ par(op); layout(1) })
  draw_r2_heatmap(R2_LIN_MAT,  "Linear R² (= r², symmetric)")
  draw_r2_heatmap(R2_BEST_MAT, "R² of the selected curve (row = response)")
  par(mar = c(4, 0.5, 3.5, 3.5))
  key <- seq(0, 1, length.out = 101)
  image(1, key, matrix(key, nrow = 1), col = R2_PAL, axes = FALSE, xlab = "", ylab = "")
  axis(4, at = seq(0, 1, 0.25), las = 1)
}, width = 12, height = 5)

save_figure("r2_gain_barplot", function() {
  o <- order(comp$Delta_R2)
  g <- comp$Delta_R2[o]
  labs <- paste0(comp$Relation[o], ifelse(comp$Direction[o] == "reverse", " (rev.)", ""))
  cols <- ifelse(comp$Direction[o] == "primary", COL_MAIN, "#9ec3ef")
  op <- par(mar = c(4.5, 13, 3, 6), las = 1)
  on.exit(par(op))
  xr <- range(c(g, 0, GAIN_MEANINGFUL)) + c(-0.01, 0.02)
  bp <- barplot(g, horiz = TRUE, names.arg = labs, col = cols, border = NA, xlim = xr,
                cex.names = 0.85, xlab = "Gain in R² = R²(best non-linear curve) - R²(linear)",
                main = "Gain in R² from the best non-linear curve over the straight line")
  abline(v = 0, col = "grey30")
  abline(v = GAIN_MEANINGFUL, lty = 2, col = COL_HIGHLIGHT, lwd = 1.5)
  text(GAIN_MEANINGFUL + 0.001, max(bp) + 0.85, sprintf("meaningful gain = %.2f", GAIN_MEANINGFUL),
       col = COL_HIGHLIGHT, cex = 0.8, adj = 0, xpd = TRUE)
  # bar labels on a white backing so the threshold line does not run through them
  labs_val <- sprintf("%+.3f %s", g, comp$Best_nonlinear[o])
  x0 <- pmax(g, 0) + 0.002
  w  <- strwidth(labs_val, cex = 0.8)
  h  <- strheight("X", cex = 0.8)
  rect(x0 - 0.0005, bp - h, x0 + w + 0.0005, bp + h, col = "white", border = NA, xpd = TRUE)
  text(x0, bp, labs_val, adj = 0, cex = 0.8, xpd = TRUE)
  legend("bottomright", legend = c("primary direction", "reverse direction"),
         fill = c(COL_MAIN, "#9ec3ef"), border = NA, bty = "n", cex = 0.85)
}, width = 10, height = 6.5)


## ---- 6.5 Ranked list of the strongest relations ------------------------------------------

prim_comp <- comp[comp$Direction == "primary", ]
prim_comp <- prim_comp[order(-prim_comp$R2_selected), ]
cat("Strongest relations (primary direction, ranked by R^2 of the selected curve)\n")
for (i in seq_len(nrow(prim_comp))) {
  s <- seas[seas$Relation == prim_comp$Relation[i], ]
  cat(sprintf("  %d. %-26s r = %6.3f | R^2 linear %.3f | selected %-11s R^2 %.3f (test %6.3f) | %-9s | anomaly R^2 %.3f\n",
              i, prim_comp$Relation[i], prim_comp$r[i], prim_comp$Linear[i], prim_comp$Selected[i],
              prim_comp$R2_selected[i], prim_comp$R2_test_selected[i], prim_comp$Strength[i],
              s$R2_anomaly))
}
ranked <- data.frame(Rank = seq_len(nrow(prim_comp)), Relation = prim_comp$Relation,
                     r = prim_comp$r, R2_linear = prim_comp$Linear, Selected = prim_comp$Selected,
                     R2_selected = prim_comp$R2_selected, R2_test = prim_comp$R2_test_selected,
                     Strength = prim_comp$Strength,
                     R2_anomaly = seas$R2_anomaly[match(prim_comp$Relation, seas$Relation)],
                     stringsAsFactors = FALSE)
save_table(ranked, "relations_ranked",
           caption = paste0("Relations ranked by $R^2$ of the selected curve (primary direction), with ",
                            "the $R^2$ that remains after removing the annual cycle (anomaly $R^2$)."),
           label = "tab:relations_ranked")
cat("\n")


## ---- 6.6 Conclusions (auto-generated draft) ------------------------------------------------
# Each finding is one bullet. Tokens {R2}, {dR2} are rendered as R^2 / Delta R^2
# in the text file and as $R^2$ / $\Delta R^2$ in the LaTeX version.

BULLETS <- character(0)
# bullet: add one finding (sprintf syntax) to the conclusions draft
bullet <- function(fmt, ...) BULLETS <<- c(BULLETS, sprintf(fmt, ...))
# fmt_list: join items with '; ' or return 'none'
fmt_list <- function(v) if (length(v) == 0) "none" else paste(v, collapse = "; ")
# sel_desc: 'relation (R^2 x, curve)' text for a relation
sel_desc <- function(row) {
  sprintf("%s ({R2} %.3f, %s)", row$Relation, row$R2_selected, row$Selected)
}

top <- prim_comp[1, ]
bullet("Data: %d consecutive days (%s to %s) after cleaning; %d invalid pressures interpolated and %d extreme wind speeds capped at %.1f km/h.",
       n_obs, format(min(dat$date)), format(max(dat$date)), length(bad_p), sum(dat$wind_capped), WIND_CAP)
bullet("Strongest relation: %s, r = %.3f. A straight line explains %.1f%% of the variance of %s ({R2} = %.3f); the %s curve raises this to %.1f%% ({R2} = %.3f, 2017 test {R2} = %.3f).",
       top$Relation, top$r, 100 * top$Linear, top$Response, top$Linear, top$Selected,
       100 * top$R2_selected, top$R2_selected, top$R2_test_selected)
for (lvl in c("strong", "moderate", "weak", "very weak")) {
  rows <- prim_comp[prim_comp$Strength == lvl, ]
  bullet("%s relations ({R2} %s): %s.",
         paste0(toupper(substring(lvl, 1, 1)), substring(lvl, 2)),
         c(strong = ">= 0.50", moderate = "0.25-0.50", weak = "0.10-0.25", `very weak` = "< 0.10")[[lvl]],
         fmt_list(sapply(seq_len(nrow(rows)), function(i) sel_desc(rows[i, ]))))
}

meaningful <- comp[comp$Delta_R2 > GAIN_MEANINGFUL, ]
mild       <- comp[comp$Delta_R2 > 0.01 & comp$Delta_R2 <= GAIN_MEANINGFUL, ]
none       <- comp[comp$Delta_R2 <= 0.01, ]
# gain_desc: text for each relation with its Delta R^2 and the out-of-sample check
gain_desc <- function(d) {
  sapply(seq_len(nrow(d)), function(i) {
    sprintf("%s%s: %s {dR2} = %+.3f (test {R2} %.3f -> %.3f, %s)", d$Relation[i],
            if (d$Direction[i] == "reverse") " (rev.)" else "", d$Best_nonlinear[i], d$Delta_R2[i],
            d$R2_test_linear[i], d$R2_test_best_nonlinear[i],
            if (d$R2_test_best_nonlinear[i] >= d$R2_test_linear[i]) "confirmed out of sample"
            else "NOT confirmed out of sample")
  })
}
bullet("Non-linear curves improve meaningfully ({dR2} > %.2f) only for: %s.", GAIN_MEANINGFUL, fmt_list(gain_desc(meaningful)))
bullet("Mild improvement (0.01 < {dR2} <= %.2f): %s.", GAIN_MEANINGFUL, fmt_list(gain_desc(mild)))
bullet("No meaningful improvement ({dR2} <= 0.01; a straight line is adequate): %s.",
       fmt_list(paste0(none$Relation, ifelse(none$Direction == "reverse", " (rev.)", ""))))

# Overfitting flags
poly_flag <- comp[comp$Best_nonlinear %in% c("quadratic", "cubic") & comp$Delta_R2 > 0.01 &
                  comp$R2_test_best_nonlinear < comp$R2_test_linear, ]
insample_only <- comp[comp$Max_gain_in_sample > 0.01 & comp$Delta_R2 <= 0.01, ]
bullet("Possible overfitting: %s. Here the polynomial gains in-sample but predicts the 2017 test period worse than the straight line, so its higher {R2} must be read with caution.",
       fmt_list(sprintf("%s%s (%s, {dR2} = %+.3f, test {R2} %.3f vs %.3f linear)", poly_flag$Relation,
                        ifelse(poly_flag$Direction == "reverse", " (rev.)", ""), poly_flag$Best_nonlinear,
                        poly_flag$Delta_R2, poly_flag$R2_test_best_nonlinear, poly_flag$R2_test_linear)))
bullet("In-sample gains that the selection rule rejected (not supported on the test period): %s.",
       fmt_list(sprintf("%s%s (%s %+.3f)", insample_only$Relation,
                        ifelse(insample_only$Direction == "reverse", " (rev.)", ""),
                        insample_only$Max_gain_curve, insample_only$Max_gain_in_sample)))
tp_grid  <- seq(min(dat$meanpressure), max(dat$meanpressure), length.out = 1000)
tp_cubic <- NL_FITS[[3]]$poly3$curve(tp_grid)
bullet("Quadratic and cubic curves bend near the edges of the data, where few days lie (e.g. the cubic temperature-pressure curve peaks at %.1f hPa and falls again towards the lowest pressures); they describe the observed range only and must not be extrapolated.",
       tp_grid[which.max(tp_cubic)])
bullet("The 2017 test period covers only January-April (114 days, one season), so test {R2} values are noisy and can be negative even for real relations; they are a check, not a definitive ranking.")

# Outliers / influence. Removing high-influence points is a sensitivity check,
# not a correction: dropping badly fitted days tends to raise R^2 by itself.
inflated <- robust[robust$Change < -0.02, ]
bullet("Influential points: removing the days with Cook's D > 4/n (%d-%d days per model) changes the linear {R2} by %+.3f to %+.3f. %s",
       min(robust$n_influential), max(robust$n_influential), min(robust$Change), max(robust$Change),
       if (nrow(inflated) > 0)
         sprintf("For %s the {R2} drops, i.e. it is inflated by a few influential days.",
                 paste(inflated$Relation, collapse = ", "))
       else if (all(robust$Change >= 0))
         "R^2 rises in every case: the influential days are poorly fitted days that LOWER the fit, so no reported R^2 is inflated by a few outliers."
       else "No R^2 is inflated by a few influential days.")
wr <- robust[!is.na(robust$R2_uncapped_wind), ]
wd <- wr$R2_linear - wr$R2_uncapped_wind
bullet("Wind capping (pre-processing): for %s the linear {R2} is %s with capped wind but %s with the raw values (difference %s). %s",
       paste(wr$Relation, collapse = ", "),
       paste(sprintf("%.3f", wr$R2_linear), collapse = " / "),
       paste(sprintf("%.3f", wr$R2_uncapped_wind), collapse = " / "),
       paste(sprintf("%+.3f", wd), collapse = " / "),
       if (max(wd) > 0.01)
         "The capping of 33 extreme values raises these R^2 values slightly, so part of the (weak) wind relations depends on that pre-processing decision; the ranking and the conclusion 'weak' do not change."
       else "The capping does not affect the conclusions.")

# Seasonality: wording driven by the computed share of R^2 that is seasonal.
bullet("Seasonality: the calendar month alone explains %s of the daily variance (R^2 of the monthly means), so temperature and pressure are dominated by the annual cycle.",
       paste(sprintf("%s %.0f%%", VARS, 100 * R2_MONTH), collapse = ", "))
bullet("Share of each linear {R2} that disappears when the annual cycle is removed (anomaly correlation): %s.",
       paste(sprintf("%s %.0f%% (r %.2f -> %.2f)", seas$Relation, 100 * seas$Seasonal_share,
                     seas$r_all, seas$r_anomaly), collapse = "; "))
# season_text: within-season correlations as text
season_text <- function(s) {
  paste(sprintf("%.2f (%s)", unlist(s[c("r_Winter", "r_Summer", "r_Monsoon", "r_Postmonsoon")]), SEASONS),
        collapse = ", ")
}
s_tp <- seas[seas$Relation == "meantemp ~ meanpressure", ]
s_th <- seas[seas$Relation == "meantemp ~ humidity", ]
s_hp <- seas[seas$Relation == "humidity ~ meanpressure", ]
bullet("Temperature-pressure: r = %.3f on the raw data but %.3f on the anomalies, so %.0f%% of the {R2} is %s. Within seasons r = %s. %s",
       s_tp$r_all, s_tp$r_anomaly, 100 * s_tp$Seasonal_share,
       if (s_tp$Seasonal_share >= 0.5) "the shared annual cycle" else "due to the annual cycle",
       season_text(s_tp),
       if (s_tp$Seasonal_share >= 0.5)
         "The very strong overall relation is therefore largely seasonal (hot, low-pressure summer/monsoon vs cool, high-pressure winter); a real but weaker day-to-day link remains."
       else "The relation is mainly a day-to-day link, not a seasonal artefact.")
mons <- dat$season == "Monsoon"
bullet("Temperature-humidity: r = %.3f on the raw data and %.3f on the anomalies, so only %.0f%% of the {R2} is due to the annual cycle; within seasons r = %s. %s",
       s_th$r_all, s_th$r_anomaly, 100 * s_th$Seasonal_share, season_text(s_th),
       if (s_th$Seasonal_share < 0.5)
         sprintf("Unlike temperature-pressure, this relation is NOT mainly seasonal: %s Mixing the seasons actually weakens the overall correlation, because monsoon days are both hot and humid (on average %.1f degC, %.0f%% humidity), which puts them on a separate band of the scatter plot; that is why no single curve in humidity reaches a high R^2.",
                 if (all(unlist(s_th[c("r_Winter", "r_Summer", "r_Monsoon", "r_Postmonsoon")]) < 0))
                   "the correlation is negative within every season (hotter days are drier), strongly in summer and monsoon."
                 else "the sign is not the same in every season.",
                 my_mean(dat$meantemp[mons]), my_mean(dat$humidity[mons]))
       else "The relation is largely seasonal.")
bullet("Humidity-pressure: r = %.3f on the raw data but %.3f on the anomalies (%.0f%% of the {R2} is seasonal). The U-shaped quadratic fit reflects the order of the seasons along the pressure axis (humid low-pressure monsoon, dry pre-monsoon, humid high-pressure winter) rather than a physical law, so the curve must not be extrapolated.",
       s_hp$r_all, s_hp$r_anomaly, 100 * s_hp$Seasonal_share)
bullet("Correlation does not imply causation: the regressions describe how the parameters co-vary, not that one causes the other. Temperature, pressure and humidity are all driven by common factors (solar heating, the monsoon circulation), and the choice of response variable is a modelling convention, not a causal direction.")
bullet("Independence: Durbin-Watson statistics of the linear fits are %.2f-%.2f (2 = no autocorrelation), so the residuals are strongly autocorrelated. The reported p-values and confidence intervals are therefore too optimistic; the R^2 values remain valid as descriptive measures of fit.",
       min(comp$DW_linear), max(comp$DW_linear))

# render_text: tokens {R2}, {dR2} -> plain text R^2, Delta R^2
render_text  <- function(s) gsub("{dR2}", "Delta R^2", gsub("{R2}", "R^2", s, fixed = TRUE), fixed = TRUE)
# render_latex: plain text -> LaTeX. Symbols are first turned into tokens so
# that escaping cannot break them, then the tokens become math.
render_latex <- function(s) {
  s <- gsub("R^2", "{R2}", s, fixed = TRUE)
  s <- gsub("->", "{TO}", s, fixed = TRUE)
  s <- gsub(">=", "{GEQ}", s, fixed = TRUE)
  s <- gsub("<=", "{LEQ}", s, fixed = TRUE)
  s <- latex_escape(s)
  s <- gsub("{dR2}", "$\\Delta R^2$", s, fixed = TRUE)
  s <- gsub("{R2}", "$R^2$", s, fixed = TRUE)
  s <- gsub("{TO}", "$\\rightarrow$", s, fixed = TRUE)
  s <- gsub("{GEQ}", "$\\geq$", s, fixed = TRUE)
  s <- gsub("{LEQ}", "$\\leq$", s, fixed = TRUE)
  s
}

cat("KEY FINDINGS (auto-generated draft of the conclusions)\n\n")
for (b in BULLETS) {
  cat(paste(strwrap(render_text(b), width = 95, initial = "  * ", prefix = "    "), collapse = "\n"), "\n")
}
writeLines(paste("-", render_text(BULLETS)), file.path(TAB_DIR, "conclusions_draft.txt"))
writeLines(c("% Generated by R/Project_code.R -- auto-generated draft, edit before use",
             "\\begin{itemize}", paste("  \\item", render_latex(BULLETS)), "\\end{itemize}"),
           file.path(TAB_DIR, "conclusions_draft.tex"))
cat(sprintf("\n  saved draft  %s/conclusions_draft.{txt,tex}\n", TAB_DIR))


## ---- 6.7 Key numbers for the LaTeX report -------------------------------------------------
# The report quotes numbers through LaTeX macros defined in key_numbers.tex
# (\input in main.tex), so every number in the text is the one computed here.
# Relation codes: T = meantemp, H = humidity, W = wind_speed, P = meanpressure;
# e.g. TP = meantemp ~ meanpressure (first letter = response).

KEY <- character(0)
# macro: add \newcommand{\name}{value} to the key-number file
macro <- function(name, value) KEY <<- c(KEY, sprintf("\\newcommand{\\%s}{%s}", name, value))
# num: number in math mode (proper minus sign), with d decimals; sgn adds a "+"
num <- function(v, d = 3, sgn = FALSE) {
  sprintf("\\ensuremath{%s}", sprintf(if (sgn) paste0("%+.", d, "f") else paste0("%.", d, "f"), v))
}

REL_CODE <- c("TH", "TW", "TP", "HW", "HP", "WP")   # same order as LIN_PAIRS
VAR_CODE <- c(meantemp = "T", humidity = "H", wind_speed = "W", meanpressure = "P")

macro("nDays", n_obs);                   macro("nTrainFile", nrow(train))
macro("nTestFile", nrow(test));          macro("nTrainPeriod", sum(IS_TRAIN))
macro("nTestPeriod", sum(IS_TEST))
macro("dateStart", format(min(dat$date))); macro("dateEnd", format(max(dat$date)))
macro("nPressFixed", length(bad_p));      macro("pressMin", P_MIN); macro("pressMax", P_MAX)
macro("pressValidLow", num(min(p_valid), 2)); macro("pressValidHigh", num(max(p_valid), 2))
macro("overlapInterp", num(dat$meanpressure[dat$date == dup_dates[1]], 2))
macro("overlapTrainP", num(overlap_train_pressure, 2))
macro("nWindCapped", sum(dat$wind_capped)); macro("windCap", num(WIND_CAP, 2))
macro("windMaxRaw", num(max(dat$wind_speed_raw), 2)); macro("nWindZthree", length(above_z3))
macro("windQone", num(w_fences$q1, 2));   macro("windQthree", num(w_fences$q3, 2))
macro("windSkewRaw", num(my_skewness(dat$wind_speed_raw), 2))
macro("windSkewClean", num(my_skewness(dat$wind_speed), 2))
macro("nWindZero", sum(dat$wind_speed == 0))
macro("monsoonTemp", num(my_mean(dat$meantemp[mons]), 1))
macro("monsoonHum", num(my_mean(dat$humidity[mons]), 0))
macro("DWmin", num(min(comp$DW_linear), 2)); macro("DWmax", num(max(comp$DW_linear), 2))
macro("cookMin", num(min(robust$Change), 3, TRUE)); macro("cookMax", num(max(robust$Change), 3, TRUE))
macro("nInflMin", min(robust$n_influential)); macro("nInflMax", max(robust$n_influential))
macro("nOverfit", length(over));          macro("nSkipped", nrow(skipped))
month_full <- c("January", "February", "March", "April", "May", "June", "July",
                "August", "September", "October", "November", "December")
hot <- which.max(monthly$meantemp); cold <- which.min(monthly$meantemp)
dry <- which.min(monthly$humidity)
macro("hotMonth", month_full[hot]);   macro("hotMonthTemp", num(monthly$meantemp[hot], 1))
macro("hotMonthPress", num(monthly$meanpressure[hot], 1))
macro("coldMonth", month_full[cold]); macro("coldMonthTemp", num(monthly$meantemp[cold], 1))
macro("coldMonthPress", num(monthly$meanpressure[cold], 1))
macro("coldMonthHum", num(monthly$humidity[cold], 0))
macro("dryMonth", month_full[dry]);   macro("dryMonthHum", num(monthly$humidity[dry], 0))
macro("pressSDraw", num(desc$SD[desc$Variable == "meanpressure" & desc$Stage == "before"], 1))
macro("pressSDclean", num(desc$SD[desc$Variable == "meanpressure" & desc$Stage == "after"], 2))
macro("cookThreshold", num(4 / n_obs, 4))
macro("maxRhoGap", num(max(abs(pair_tab$rho_minus_r))))
macro("spikeRatioMax", num(max(w[above_z3] / w_rmed[above_z3]), 0))
macro("windCapGainMin", num(min(wd))); macro("windCapGainMax", num(max(wd)))
macro("nGNconv", sum(gn_rows$Converged)); macro("nGNtotal", nrow(gn_rows))
macro("GNitMin", min(gn_rows$Iterations)); macro("GNitMax", max(gn_rows$Iterations))
for (v in VARS) {
  macro(paste0("Rmonth", VAR_CODE[[v]]), num(100 * R2_MONTH[[v]], 0))
  macro(paste0("acf", VAR_CODE[[v]]), num(ACF1[[v]], 2))
}
for (k in seq_along(LIN_PAIRS)) {
  yv <- LIN_PAIRS[[k]][["y"]]; xv <- LIN_PAIRS[[k]][["x"]]; X <- REL_CODE[k]
  r  <- P_MAT[yv, xv]
  ci <- fisher_ci(r, n_obs)
  aa <- ACF1[[xv]] * ACF1[[yv]]
  li <- lin_tab[lin_tab$Pair == k & lin_tab$Direction == "primary", ]
  cp <- comp[comp$Pair == k & comp$Direction == "primary", ]
  sz <- seas[k, ]
  rb <- robust[k, ]
  macro(paste0("r", X), num(r));                macro(paste0("rho", X), num(S_MAT[yv, xv]))
  macro(paste0("rCIlo", X), num(ci[["lower"]])); macro(paste0("rCIhi", X), num(ci[["upper"]]))
  macro(paste0("neff", X), num(n_obs * (1 - aa) / (1 + aa), 0))
  macro(paste0("Rlin", X), num(li$R2));          macro(paste0("RlinTest", X), num(li$R2_test))
  macro(paste0("RlinTrain", X), num(li$R2_train))
  macro(paste0("intercept", X), num(li$a, 2));   macro(paste0("slope", X), num(li$b, 4))
  macro(paste0("slopeCIlo", X), num(li$CI_b_lower, 4)); macro(paste0("slopeCIhi", X), num(li$CI_b_upper, 4))
  macro(paste0("slopeRev", X), num(lin_tab$b[lin_tab$Pair == k & lin_tab$Direction == "reverse"], 4))
  macro(paste0("DW", X), num(li$DW, 2))
  macro(paste0("best", X), cp$Selected);         macro(paste0("Rbest", X), num(cp$R2_selected))
  macro(paste0("RbestTest", X), num(cp$R2_test_selected))
  macro(paste0("bestNL", X), cp$Best_nonlinear); macro(paste0("gain", X), num(cp$Delta_R2, 3, TRUE))
  macro(paste0("RbestNLTest", X), num(cp$R2_test_best_nonlinear))
  macro(paste0("maxGain", X), num(cp$Max_gain_in_sample, 3, TRUE))
  macro(paste0("maxGainCurve", X), cp$Max_gain_curve)
  macro(paste0("ranom", X), num(sz$r_anomaly)); macro(paste0("Ranom", X), num(sz$R2_anomaly))
  macro(paste0("seasShare", X), num(100 * sz$Seasonal_share, 0))
  macro(paste0("rWinter", X), num(sz$r_Winter, 2)); macro(paste0("rSummer", X), num(sz$r_Summer, 2))
  macro(paste0("rMonsoon", X), num(sz$r_Monsoon, 2)); macro(paste0("rPost", X), num(sz$r_Postmonsoon, 2))
  macro(paste0("cookChange", X), num(rb$Change, 3, TRUE)); macro(paste0("nInfl", X), rb$n_influential)
  if (!is.na(rb$R2_uncapped_wind)) macro(paste0("RrawWind", X), num(rb$R2_uncapped_wind))
}
# cubic temperature-pressure curve: location of its maximum
macro("cubicPeakTP", num(tp_grid[which.max(tp_cubic)], 1))

# Line ranges of the key functions in this file, so that the code appendix of
# the report (\lstinputlisting with firstline/lastline) always shows the
# current code: from the first line of the comment block above the function
# to its closing brace.
SRC <- readLines(file.path("R", "Project_code.R"), encoding = "UTF-8")
# fun_range: first comment line above "name <- function" .. closing "}" line
fun_range <- function(name) {
  s <- grep(paste0("^", name, " <- function"), SRC)[1]
  e <- s - 1 + which(SRC[s:length(SRC)] == "}")[1]
  while (s > 1 && grepl("^#", SRC[s - 1])) s <- s - 1
  c(s, e)
}
KEY_FUNS <- c(my_rank = "Rank", my_pearson = "Pearson", my_slr = "Slr",
              solve_linear_system = "Solve", gauss_newton = "GaussNewton",
              fit_poly_model = "Poly", fit_exp_model = "Exp")
for (fn in names(KEY_FUNS)) {
  rg <- fun_range(fn)
  macro(paste0("line", KEY_FUNS[[fn]], "First"), rg[1])
  macro(paste0("line", KEY_FUNS[[fn]], "Last"), rg[2])
}

writeLines(c("% Generated by R/Project_code.R -- key numbers quoted in the report", KEY),
           file.path(TAB_DIR, "key_numbers.tex"))
cat(sprintf("  saved macros %s/key_numbers.tex (%d numbers)\n", TAB_DIR, length(KEY)))

cat("\nDone. All results are in outputs/.\n")
