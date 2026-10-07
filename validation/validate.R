# =============================================================================
#  validation/validate.R  --  OUR OWN CHECK, NOT PART OF THE SUBMISSION
#
#  Compares the hand-written results of R/Project_code.R (read from its output
#  tables) with R's built-in functions (lm, cor.test, nls, quantile, approx).
#  Any failure must be fixed in R/Project_code.R, never here.
#
#  Usage, from the project root:
#    Rscript validation/validate.R            # check the existing outputs
#    Rscript validation/validate.R --rerun    # re-run Project_code.R first
#
#  Difference reported: |ours - reference| / max(1, |reference|), i.e. the
#  absolute difference, made relative only for values larger than 1 in size
#  (intercepts around 1000, exponential coefficients around 1e15).
# =============================================================================

TOL       <- 1e-6   # closed-form results
TOL_GN    <- 1e-4   # Gauss-Newton (ours) vs nls: both stop at a tolerance
TOL_FISHER <- 1e-4  # we hard-code z = 1.96; cor.test uses qnorm(0.975) = 1.959964

if (!file.exists("R/Project_code.R")) stop("Run from the project root: Rscript validation/validate.R")

if ("--rerun" %in% commandArgs(trailingOnly = TRUE)) {
  cat("Re-running R/Project_code.R ...\n")
  rscript <- file.path(R.home("bin"), "Rscript")
  status <- system2(rscript, "R/Project_code.R", stdout = FALSE, stderr = FALSE)
  if (status != 0) stop("R/Project_code.R failed (exit status ", status, ")")
}

TAB <- "outputs/tables"
needed <- file.path(TAB, c("correlation_pairs_ranked.csv", "linear_results.csv",
                           "nonlinear_results.csv", "descriptive_stats.csv"))
needed <- c(needed, "data/processed/clean.csv")
missing <- needed[!file.exists(needed)]
if (length(missing) > 0) stop("Missing outputs (use --rerun): ", paste(missing, collapse = ", "))
if (any(file.mtime(needed) < file.mtime("R/Project_code.R"))) {
  cat("WARNING: some outputs are older than R/Project_code.R -- consider --rerun\n\n")
}

d       <- read.csv("data/processed/clean.csv", stringsAsFactors = FALSE)
d$date  <- as.Date(d$date)
train   <- d$date <= as.Date("2016-12-31")
test    <- !train
cors    <- read.csv(file.path(TAB, "correlation_pairs_ranked.csv"), stringsAsFactors = FALSE)
lin     <- read.csv(file.path(TAB, "linear_results.csv"), stringsAsFactors = FALSE)
nl      <- read.csv(file.path(TAB, "nonlinear_results.csv"), stringsAsFactors = FALSE)
desc    <- read.csv(file.path(TAB, "descriptive_stats.csv"), stringsAsFactors = FALSE)
VARS    <- c("meantemp", "humidity", "wind_speed", "meanpressure")

RES <- list()
check <- function(section, item, ours, ref, tol = TOL) {
  diff <- abs(ours - ref) / max(1, abs(ref))
  status <- if (is.na(ours) || is.na(ref)) "NA" else if (diff <= tol) "ok" else "FAIL"
  RES[[length(RES) + 1]] <<- data.frame(Section = section, Item = item, Ours = ours, Reference = ref,
                                        Diff = diff, Tol = tol, Status = status,
                                        stringsAsFactors = FALSE)
}
note <- function(section, item, text) {
  RES[[length(RES) + 1]] <<- data.frame(Section = section, Item = item, Ours = NA, Reference = NA,
                                        Diff = NA, Tol = NA, Status = paste("SKIP:", text),
                                        stringsAsFactors = FALSE)
}
r2_of <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)


## ---- 1. Descriptive statistics -------------------------------------------------
for (v in VARS) {
  x  <- d[[v]]
  dv <- desc[desc$Variable == v, ]
  q  <- quantile(x, c(0.25, 0.5, 0.75), type = 7, names = FALSE)
  check("descriptive", paste(v, "mean"),   dv$Mean,   mean(x))
  check("descriptive", paste(v, "sd"),     dv$SD,     sd(x))
  check("descriptive", paste(v, "Q1"),     dv$Q1,     q[1])
  check("descriptive", paste(v, "median"), dv$Median, median(x))
  check("descriptive", paste(v, "Q3"),     dv$Q3,     q[3])
  check("descriptive", paste(v, "min"),    dv$Min,    min(x))
  check("descriptive", paste(v, "max"),    dv$Max,    max(x))
}


## ---- 2. Pressure interpolation vs approx() ------------------------------------------
imp <- d$pressure_imputed == 1
ref <- approx(as.numeric(d$date[!imp]), d$meanpressure[!imp], xout = as.numeric(d$date[imp]))$y
for (i in seq_along(ref)) {
  check("interpolation", paste("meanpressure", format(d$date[imp][i])), d$meanpressure[imp][i], ref[i])
}


## ---- 3. Correlations ---------------------------------------------------------------------
for (i in seq_len(nrow(cors))) {
  x  <- d[[cors$x[i]]]; y <- d[[cors$y[i]]]
  ct <- cor.test(x, y, method = "pearson")
  lab <- cors$Pair[i]
  check("correlation", paste(lab, "Pearson r"),  cors$r[i],        unname(ct$estimate))
  check("correlation", paste(lab, "t"),          cors$t[i],        unname(ct$statistic))
  check("correlation", paste(lab, "p"),          cors$p[i],        ct$p.value)
  check("correlation", paste(lab, "CI lower"),   cors$CI_lower[i], ct$conf.int[1], TOL_FISHER)
  check("correlation", paste(lab, "CI upper"),   cors$CI_upper[i], ct$conf.int[2], TOL_FISHER)
  check("correlation", paste(lab, "Spearman rho"), cors$rho[i],    cor(x, y, method = "spearman"))
}


## ---- 4. Simple linear regression vs lm() --------------------------------------------------
for (i in seq_len(nrow(lin))) {
  y <- d[[lin$Response[i]]]; x <- d[[lin$Predictor[i]]]
  f <- lm(y ~ x); s <- summary(f); ci <- confint(f)["x", ]
  lab <- lin$Model[i]
  check("linear", paste(lab, "intercept a"), lin$a[i],          coef(f)[[1]])
  check("linear", paste(lab, "slope b"),     lin$b[i],          coef(f)[[2]])
  check("linear", paste(lab, "SE(a)"),       lin$SE_a[i],       s$coefficients[1, 2])
  check("linear", paste(lab, "SE(b)"),       lin$SE_b[i],       s$coefficients[2, 2])
  check("linear", paste(lab, "t(b)"),        lin$t_b[i],        s$coefficients[2, 3])
  check("linear", paste(lab, "p(b)"),        lin$p_b[i],        s$coefficients[2, 4])
  check("linear", paste(lab, "CI(b) lower"), lin$CI_b_lower[i], ci[[1]])
  check("linear", paste(lab, "CI(b) upper"), lin$CI_b_upper[i], ci[[2]])
  check("linear", paste(lab, "R2"),          lin$R2[i],         s$r.squared)
  check("linear", paste(lab, "adj R2"),      lin$adj_R2[i],     s$adj.r.squared)
  check("linear", paste(lab, "RSE"),         lin$RSE[i],        s$sigma)
  check("linear", paste(lab, "RMSE"),        lin$RMSE[i],       sqrt(mean(resid(f)^2)))
  check("linear", paste(lab, "MAE"),         lin$MAE[i],        mean(abs(resid(f))))
  check("linear", paste(lab, "Durbin-Watson"), lin$DW[i],       sum(diff(resid(f))^2) / sum(resid(f)^2))
  ft <- lm(y ~ x, subset = train)
  check("linear", paste(lab, "R2 train"),    lin$R2_train[i],   summary(ft)$r.squared)
  check("linear", paste(lab, "R2 test"),     lin$R2_test[i],
        r2_of(y[test], predict(ft, newdata = data.frame(x = x[test]))))
}


## ---- 5. Non-linear models -------------------------------------------------------------------
# Polynomials via lm(y ~ poly(x, d)) (orthogonal basis: same fitted values, no
# ill-conditioning), logarithmic via lm(y ~ log(x)), exponential and power via
# nls with independent start values from lm on the log scale. The nls models
# are written in centred form (as nls fails on y = a*exp(b*x) with x ~ 1000 hPa);
# the curve and its R^2 are the same.

nls_fit <- function(formula, data, start) {
  f <- tryCatch(nls(formula, data, start = start, control = nls.control(maxiter = 500)),
                error = function(e) NULL)
  if (is.null(f)) {
    f <- tryCatch(nls(formula, data, start = start, algorithm = "port",
                      control = nls.control(maxiter = 500)), error = function(e) NULL)
  }
  f
}

fit_reference <- function(model, x, y) {
  df <- data.frame(x = x, y = y)
  pos <- y > 0
  if (model == "poly2")       return(list(fit = lm(y ~ poly(x, 2), df), k = 2))
  if (model == "poly3")       return(list(fit = lm(y ~ poly(x, 3), df), k = 3))
  if (model == "logarithmic") return(list(fit = lm(y ~ log(x), df), k = 1))
  if (model == "exponential") {
    cx <- mean(x)
    st <- coef(lm(log(y) ~ I(x - cx), df, subset = pos))
    f  <- nls_fit(y ~ A * exp(b * (x - cx)), df, list(A = exp(st[[1]]), b = st[[2]]))
    if (is.null(f)) return(NULL)
    return(list(fit = f, k = 1, a = coef(f)[["A"]] * exp(-coef(f)[["b"]] * cx), b = coef(f)[["b"]]))
  }
  if (model == "power") {
    cx <- exp(mean(log(x)))
    st <- coef(lm(log(y) ~ log(x / cx), df, subset = pos))
    f  <- nls_fit(y ~ A * (x / cx)^b, df, list(A = exp(st[[1]]), b = st[[2]]))
    if (is.null(f)) return(NULL)
    return(list(fit = f, k = 1, a = coef(f)[["A"]] * cx^(-coef(f)[["b"]]), b = coef(f)[["b"]]))
  }
}

for (i in which(nl$Model != "linear")) {
  m   <- nl$Model[i]
  lab <- paste(nl$Relation[i], m)
  y <- d[[nl$Response[i]]]; x <- d[[nl$Predictor[i]]]
  needs_pos_x <- m %in% c("logarithmic", "power")
  if (needs_pos_x && any(x <= 0)) {
    # our script must have skipped it as well
    check("nonlinear", paste(lab, "skipped (x <= 0)"), as.numeric(!nl$Valid[i]), 1)
    next
  }
  tol <- if (m %in% c("exponential", "power")) TOL_GN else TOL
  ref <- fit_reference(m, x, y)
  if (is.null(ref)) { note("nonlinear", lab, "nls did not converge"); next }
  yhat <- fitted(ref$fit)
  n <- length(y)
  R2  <- r2_of(y, yhat)
  check("nonlinear", paste(lab, "R2 (original scale)"), nl$R2[i], R2, tol)
  check("nonlinear", paste(lab, "adj R2"), nl$adj_R2[i], 1 - (1 - R2) * (n - 1) / (n - ref$k - 1), tol)
  check("nonlinear", paste(lab, "RMSE"),   nl$RMSE[i],   sqrt(mean((y - yhat)^2)), tol)
  if (m %in% c("poly2", "poly3")) {
    check("nonlinear", paste(lab, "R2 = lm r.squared"), nl$R2[i], summary(ref$fit)$r.squared)
  }
  if (m %in% c("exponential", "power")) {
    check("nonlinear", paste(lab, "a"), nl$p0[i], ref$a, tol)
    check("nonlinear", paste(lab, "b"), nl$p1[i], ref$b, tol)
  }
  if (m == "logarithmic") {
    check("nonlinear", paste(lab, "a"), nl$p0[i], coef(ref$fit)[[1]])
    check("nonlinear", paste(lab, "b"), nl$p1[i], coef(ref$fit)[[2]])
  }
  # train -> test
  rt <- fit_reference(m, x[train], y[train])
  if (is.null(rt)) { note("nonlinear", paste(lab, "R2 test"), "nls (train) did not converge"); next }
  yhat_te <- predict(rt$fit, newdata = data.frame(x = x[test]))
  check("nonlinear", paste(lab, "R2 test"), nl$R2_test[i], r2_of(y[test], yhat_te), tol)
}
# the linear rows of the non-linear table must equal section 4
for (i in which(nl$Model == "linear")) {
  li <- lin[lin$Model == nl$Relation[i], ]
  check("nonlinear", paste(nl$Relation[i], "linear R2 = section 4"), nl$R2[i], li$R2)
}


## ---- Report ---------------------------------------------------------------------------------
res <- do.call(rbind, RES)
write.csv(res, "validation/validation_results.csv", row.names = FALSE)

old <- options(width = 200)
show <- data.frame(Section = res$Section, Item = res$Item,
                   Ours = formatC(res$Ours, digits = 10, format = "g"),
                   Reference = formatC(res$Reference, digits = 10, format = "g"),
                   Diff = formatC(res$Diff, digits = 2, format = "e"),
                   Tol = formatC(res$Tol, digits = 0, format = "e"), Status = res$Status)
print(show, row.names = FALSE, right = FALSE)
options(old)

cat("\n================ SUMMARY ================\n")
for (s in unique(res$Section)) {
  r <- res[res$Section == s, ]
  cat(sprintf("  %-14s %3d checks: %3d ok, %d FAIL, %d skipped/NA, max diff %.2e\n", s, nrow(r),
              sum(r$Status == "ok"), sum(r$Status == "FAIL"),
              sum(!(r$Status %in% c("ok", "FAIL"))), max(r$Diff, na.rm = TRUE)))
}
fails <- res[res$Status == "FAIL", ]
if (nrow(fails) > 0) {
  cat("\nFAILURES (fix in R/Project_code.R):\n")
  print(fails, row.names = FALSE)
  quit(status = 1)
} else {
  cat(sprintf("\nAll %d checks passed (tolerance %.0e; %.0e for Gauss-Newton and the 1.96 Fisher interval).\n",
              sum(res$Status == "ok"), TOL, TOL_GN))
  cat("Full table: validation/validation_results.csv\n")
}
