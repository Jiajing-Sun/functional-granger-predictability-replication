#!/usr/bin/env Rscript
# Monte Carlo study for relevant functional Granger predictability.
# The statistic is the residual--predictor partial cross-covariance:
#   C_hat = n^{-1} sum u_t \otimes v_{t-1},
# where both u_t and v_{t-1} are residualized on the restricted Y history.

rm(list = ls())
options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
} else {
  getwd()
}
repo_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

source(file.path(script_dir, "FGC_cv_lambda20_fixed.R"))
source(file.path(script_dir, "FGC_simulation_functions_lambda20_plugin.R"))

alpha <- 0.05
m_grid <- 20L
D <- 21L
burn <- 200L
a_y <- 0.4
c_x <- 0.5
psi <- 0.7
theta0 <- 0.30
theta_grid <- c(0, 0.10, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50)
n_list <- c(100L, 200L, 400L)
error_types <- c("iid", "fMA1", "t5")
R <- as.integer(Sys.getenv("SN_MAIN_R", "1000"))
if (is.na(R) || R < 1L) stop("SN_MAIN_R must be a positive integer.")
lambda_fun <- function(n) n^(-1)
plugin_bw <- NULL

out_dir_env <- Sys.getenv("SN_OUTPUT_DIR", unset = "")
out_dir <- if (nzchar(out_dir_env)) {
  normalizePath(out_dir_env, mustWork = FALSE)
} else {
  file.path(repo_dir, "processed", "simulation")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

q_quad <- get_cv_lambda20(alpha = alpha, which = "quad")
q_range <- get_cv_lambda20(alpha = alpha, which = "range")
write.csv(cv_lambda20_fixed(), file.path(out_dir, "critical_values_lambda20_fixed.csv"),
          row.names = FALSE)

b_vec <- 1 / seq_len(D)
sigma_vec <- 1 / seq_len(D)
thresholds <- design_specific_duX_thresholds(
  theta0 = theta0, a_y = a_y, c_x = c_x,
  b_vec = b_vec, sigma_vec = sigma_vec, psi = psi
)
psi0_check_error <- attr(thresholds, "psi0_white_noise_error")
zero_signal_error <- attr(thresholds, "zero_signal_error")
Delta_by_error <- setNames(thresholds$Delta, thresholds$error)
write.csv(thresholds, file.path(out_dir, "design_specific_thresholds.csv"),
          row.names = FALSE)

cat(sprintf("Simulation study: R=%d, D=%d, theta0=%.2f, eta_n=n^{-1}\n",
            R, D, theta0))
for (i in seq_len(nrow(thresholds))) {
  cat(sprintf(
    "  threshold[%s]: Delta=%.12g, boundary_error=%+.3e, method=%s\n",
    thresholds$error[i], thresholds$Delta[i], thresholds$boundary_error[i],
    thresholds$second_moment_method[i]
  ))
}
cat(sprintf(
  "Self-check psi=0 augmented-state versus white-noise formula: error=%+.3e\n",
  psi0_check_error
))
cat(sprintf(
  "Self-check theta=0 population effect across white/fMA1 designs: max=%+.3e\n",
  zero_signal_error
))
cat(sprintf("Critical values: q_quad=%.3f, q_range=%.3f\n", q_quad, q_range))

one_rep <- function(rep_id, n, theta, err, base_seed = 730531L) {
  err_offset <- match(err, error_types) * 10000000L
  seed <- base_seed + err_offset + 100000L * n + 1000L * round(theta * 100) + rep_id
  dat <- sim_FARX1_diag(
    n = n, D = D, theta = theta,
    a_y = a_y, c_x = c_x,
    b_vec = b_vec, sigma_vec = sigma_vec,
    error_type = err, psi = psi,
    burn = burn, seed = seed
  )
  sn <- relevant_SN_test_duX(
    X = dat$X, Y = dat$Y,
    Delta = unname(Delta_by_error[[err]]),
    q_quad = q_quad,
    q_range = q_range,
    m = m_grid,
    lambda_fun = lambda_fun
  )
  pl <- plugin_test_duX_Delta0(
    X = dat$X, Y = dat$Y,
    alpha = alpha,
    lambda_fun = lambda_fun,
    L_bw = plugin_bw
  )
  c(rej_quad = sn$reject_quad,
    rej_range = sn$reject_range,
    rej_plugin = pl$reject)
}

suggest_ncores <- function(reserve = 2L, unix_cap = 10L, env_var = "SN_NCORES") {
  env <- Sys.getenv(env_var, unset = "")
  if (nzchar(env)) {
    nc <- suppressWarnings(as.integer(env))
    if (!is.na(nc) && nc >= 1L) return(nc)
  }
  nc <- parallel::detectCores()
  if (is.na(nc) || nc < 1L) nc <- 1L
  min(max(1L, nc - reserve), unix_cap)
}

use_parallel <- TRUE
n_cores <- suggest_ncores()
cat(sprintf("Parallel: %s, n_cores=%d\n", if (use_parallel) "ON" else "OFF", n_cores))

parallel_apply <- function(X, FUN) {
  if (!use_parallel || n_cores <= 1L) return(lapply(X, FUN))
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
    parallel::parLapply(cl, X, FUN)
  } else {
    parallel::mclapply(X, FUN, mc.cores = n_cores)
  }
}

grid <- expand.grid(n = n_list, theta = theta_grid, error = error_types,
                    stringsAsFactors = FALSE)
all_results <- vector("list", nrow(grid))
t0 <- proc.time()[3]

for (g in seq_len(nrow(grid))) {
  cfg <- grid[g, ]
  cat(sprintf("[%02d/%02d] n=%d, theta=%.2f, error=%s\n",
              g, nrow(grid), cfg$n, cfg$theta, cfg$error))
  sim_list <- parallel_apply(seq_len(R), function(r) {
    one_rep(rep_id = r, n = cfg$n, theta = cfg$theta, err = cfg$error)
  })
  sim_mat <- do.call(rbind, sim_list)
  if (any(!is.finite(sim_mat))) {
    stop(sprintf("Non-finite simulation result at n=%d, theta=%.2f, error=%s.",
                 cfg$n, cfg$theta, cfg$error))
  }
  rej_quad <- mean(sim_mat[, "rej_quad"])
  rej_range <- mean(sim_mat[, "rej_range"])
  rej_plugin <- mean(sim_mat[, "rej_plugin"])
  all_results[[g]] <- data.frame(
    n = cfg$n,
    theta = cfg$theta,
    error = cfg$error,
    Delta = unname(Delta_by_error[[cfg$error]]),
    R = R,
    rej_quad = rej_quad,
    rej_range = rej_range,
    rej_plugin = rej_plugin,
    se_quad = sqrt(rej_quad * (1 - rej_quad) / R),
    se_range = sqrt(rej_range * (1 - rej_range) / R),
    se_plugin = sqrt(rej_plugin * (1 - rej_plugin) / R),
    stringsAsFactors = FALSE
  )
}

res_df <- do.call(rbind, all_results)
res_df <- res_df[order(res_df$error, res_df$n, res_df$theta), ]
write.csv(res_df, file.path(out_dir, "rejection_probabilities.csv"), row.names = FALSE)
size_boundary <- subset(res_df, abs(theta - theta0) < 1e-12)
write.csv(size_boundary, file.path(out_dir, "size_boundary.csv"), row.names = FALSE)
cat(sprintf("Simulation finished in %.2f seconds\n", proc.time()[3] - t0))

fmt_prob <- function(p, se) sprintf("%.3f (%.3f)", p, se)

write_boundary_tex <- function(dt, thresholds, n_rep, path) {
  err_order <- c("iid", "fMA1", "t5")
  delta_iid <- thresholds$Delta[match("iid", thresholds$error)]
  delta_fma1 <- thresholds$Delta[match("fMA1", thresholds$error)]
  delta_t5 <- thresholds$Delta[match("t5", thresholds$error)]
  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\begin{threeparttable}",
    "\\caption{Boundary rejection probabilities at the design-specific relevance thresholds evaluated at $\\theta=\\theta_0=0.30$ (Monte Carlo s.e.\\ in parentheses).}",
    "\\label{tab:sizeBoundary}",
    "\\begin{spacing}{1}",
    "\\begin{tabular}{llccc}",
    "\\toprule",
    "Disturbance design & $n$ & Quadratic SN & Adjusted-range SN & Plug-in ($\\Delta=0$) \\\\",
    "\\midrule"
  )
  labels <- c(iid = "iid", fMA1 = "fMA(1)", t5 = "$t_5$")
  for (err in err_order) {
    sub <- dt[dt$error == err, ]
    sub <- sub[order(sub$n), ]
    for (i in seq_len(nrow(sub))) {
      left <- if (i == 1L) labels[[err]] else ""
      lines <- c(lines, sprintf(
        "%s & %d & %s & %s & %s \\\\",
        left, sub$n[i],
        fmt_prob(sub$rej_quad[i], sub$se_quad[i]),
        fmt_prob(sub$rej_range[i], sub$se_range[i]),
        fmt_prob(sub$rej_plugin[i], sub$se_plugin[i])
      ))
    }
  }
  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{spacing}",
    "\\begin{tablenotes}[flushleft]",
    "\\footnotesize",
    sprintf("\\item \\textit{Notes:} Nominal level $\\alpha=0.05$ and $R=%d$ Monte Carlo replications.", n_rep),
    sprintf(paste0(
      "The relevance boundary is design-specific: $\\Delta_{\\mathrm{iid}}=%.8g$, ",
      "$\\Delta_{\\mathrm{fMA(1)}}=%.8g$, and $\\Delta_{t_5}=%.8g$, ",
      "each equal to the corresponding population $d_{uX}(\\theta_0)$ at $\\theta_0=0.30$."
    ), delta_iid, delta_fma1, delta_t5),
    "The iid and $t_5$ thresholds use matched-variance white-disturbance second moments; the fMA(1) threshold uses stationary second moments from the augmented-state discrete Lyapunov equation.",
    "The fMA(1) case is a serially correlated disturbance stress design, not an MDS-innovation design.",
    "The SN normalizers use the 19-point grid $\\lambda_i=0.25+0.75i/20$ ($m=20$).",
    "The positive-threshold SN columns use the diagonal correction; the point-null column is an approximate diagonal-HAC/Satterthwaite benchmark.",
    "Both $Y_t$ and $X_{t-1}$ are residualized on the restricted $Y$ history; a prefix with $q$ usable lagged pairs uses ridge penalty $\\eta_q=q^{-1}$.",
    "Disturbance designs: iid = Gaussian i.i.d.; {fMA(1)} = serially correlated functional MA(1) disturbances $e_{t,k}=z_{t,k}+0.7z_{t-1,k}$ scaled to match $\\operatorname{Var}(e_{t,k})=\\sigma_k^2$; {$t_5$} = i.i.d. Student-$t$ disturbances with 5 degrees of freedom rescaled to variance $\\sigma_k^2$.",
    "Monte Carlo s.e.\\ are $\\sqrt{\\hat p(1-\\hat p)/R}$.",
    "The first two SN columns correspond to the proposed relevant test of $H_0(\\Delta): d_{uX}\\le \\Delta$; therefore the boundary row $\\theta=\\theta_0$ is a size diagnostic at the least favorable null point.",
    "The last column (Plug-in with $\\Delta=0$) tests the classical point null $H_0: d_{uX}=0$ and is evaluated at $\\theta=\\theta_0>0$; hence it is power, not a size diagnostic for the relevant null.",
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\end{table}"
  )
  writeLines(lines, path, useBytes = TRUE)
}
write_boundary_tex(size_boundary, thresholds, R,
                   file.path(out_dir, "tab_size_boundary.tex"))

error_panel_titles <- c(
  iid = "iid Gaussian",
  fMA1 = "serial fMA(1)",
  t5 = "iid t(5)"
)

simulation_pdf_mode <- NULL

open_simulation_pdf <- function(out_pdf) {
  if (is.null(simulation_pdf_mode) || identical(simulation_pdf_mode, "cairo")) {
    device_before <- dev.cur()
    cairo_ok <- tryCatch(
      {
        withCallingHandlers(
          grDevices::cairo_pdf(
            filename = out_pdf,
            width = 10.5,
            height = 4.2,
            pointsize = 14,
            family = "sans",
            bg = "white"
          ),
          warning = function(w) stop(conditionMessage(w), call. = FALSE)
        )
        TRUE
      },
      error = function(e) {
        if (dev.cur() != device_before) dev.off()
        FALSE
      }
    )
    if (cairo_ok) {
      simulation_pdf_mode <<- "cairo"
      return(list(mode = "cairo", out_pdf = out_pdf))
    }
    simulation_pdf_mode <<- "embed"
    message("Cairo PDF is unavailable; using pdf() plus Ghostscript font embedding.")
  }

  if (!nzchar(Sys.which("gs"))) {
    stop("Neither a working Cairo PDF device nor Ghostscript is available for font embedding.")
  }
  raw_pdf <- tempfile(pattern = "fgc_sim_", fileext = ".pdf")
  grDevices::pdf(
    file = raw_pdf,
    width = 10.5,
    height = 4.2,
    pointsize = 14,
    family = "Helvetica",
    bg = "white",
    useDingbats = FALSE
  )
  list(mode = "embed", raw_pdf = raw_pdf, out_pdf = out_pdf)
}

close_simulation_pdf <- function(device) {
  dev.off()
  if (identical(device$mode, "embed")) {
    on.exit(unlink(device$raw_pdf), add = TRUE)
    if (file.exists(device$out_pdf)) unlink(device$out_pdf)
    grDevices::embedFonts(device$raw_pdf, outfile = device$out_pdf)
  }
  invisible(device$out_pdf)
}

plot_curves_single <- function(df, ycol, main_prefix, out_pdf) {
  pdf_device <- open_simulation_pdf(out_pdf)
  op <- par(
    mfrow = c(1, 3),
    mar = c(4.6, 4.8, 4.0, 1.0),
    mgp = c(2.9, 0.9, 0),
    tcl = -0.3,
    cex.axis = 1.05,
    cex.lab = 1.12,
    cex.main = 0.95
  )
  on.exit({par(op); close_simulation_pdf(pdf_device)}, add = TRUE)
  errs <- c("iid", "fMA1", "t5")
  ns <- sort(unique(df$n))
  ltys <- seq_along(ns)
  for (j in seq_along(errs)) {
    err <- errs[j]
    sub <- df[df$error == err, ]
    plot(range(df$theta), c(0, 1), type = "n",
         xlab = "Signal coefficient (theta)", ylab = "Rejection probability",
         main = paste0(main_prefix, "\n", error_panel_titles[[err]]))
    abline(h = alpha, lty = 3)
    abline(v = theta0, lty = 2)
    for (i in seq_along(ns)) {
      n0 <- ns[i]
      s2 <- sub[sub$n == n0, ]
      s2 <- s2[order(s2$theta), ]
      lines(s2$theta, s2[[ycol]], lty = ltys[i], lwd = 2)
    }
    if (j == 1L) legend("topleft", legend = paste0("n=", ns),
                        lty = ltys, lwd = 2, bty = "n", cex = 0.95)
  }
}

plot_curves_both <- function(df, out_pdf) {
  pdf_device <- open_simulation_pdf(out_pdf)
  op <- par(
    mfrow = c(1, 3),
    mar = c(4.6, 4.8, 4.0, 1.0),
    mgp = c(2.9, 0.9, 0),
    tcl = -0.3,
    cex.axis = 1.05,
    cex.lab = 1.12,
    cex.main = 0.95
  )
  on.exit({par(op); close_simulation_pdf(pdf_device)}, add = TRUE)
  errs <- c("iid", "fMA1", "t5")
  ns <- sort(unique(df$n))
  ltys <- seq_along(ns)
  for (j in seq_along(errs)) {
    err <- errs[j]
    sub <- df[df$error == err, ]
    ymax <- c(iid = 0.9, fMA1 = 0.6, t5 = 0.6)[[err]]
    ymax <- max(ymax, ceiling(10 * max(sub$rej_quad, sub$rej_range)) / 10)
    plot(range(df$theta), c(0, ymax), type = "n", yaxs = "i", yaxt = "n",
         xlab = "Signal coefficient (theta)", ylab = "Rejection probability",
         main = paste0("Quadratic vs range\n", error_panel_titles[[err]]))
    axis(2, at = seq(0, ymax, by = 0.1))
    abline(h = alpha, lty = 3)
    abline(v = theta0, lty = 2)
    for (i in seq_along(ns)) {
      n0 <- ns[i]
      s2 <- sub[sub$n == n0, ]
      s2 <- s2[order(s2$theta), ]
      lines(s2$theta, s2$rej_quad, lty = ltys[i], lwd = 2, col = "black")
      lines(s2$theta, s2$rej_range, lty = ltys[i], lwd = 2, col = "red")
    }
    if (j == 1L) {
      legend("topleft",
             legend = c(paste0("quad n=", ns), paste0("range n=", ns)),
             lty = rep(ltys, 2), lwd = 2,
             col = c(rep("black", length(ns)), rep("red", length(ns))),
             bty = "n", cex = 0.90, ncol = 2)
    }
  }
}

plot_curves_single(res_df, "rej_quad", "Quadratic SN",
                   file.path(out_dir, "power_curves_quad.pdf"))
plot_curves_single(res_df, "rej_range", "Adjusted-range SN",
                   file.path(out_dir, "power_curves_range.pdf"))
plot_curves_single(res_df, "rej_plugin", "Plug-in Delta=0",
                   file.path(out_dir, "power_curves_plugin_Delta0.pdf"))
plot_curves_both(res_df, file.path(out_dir, "power_curves_both.pdf"))

cat("Saved simulation outputs in ", normalizePath(out_dir), "\n", sep = "")
