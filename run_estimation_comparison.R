#!/usr/bin/env Rscript

## Common-slope estimation after method-specific anomaly removal.
## The data-generating design matches Section 4.1, except that the signal
## levels are zeta = 0.10, 0.30, and 0.50 and rho = 0.20.
##
## Environment variables:
##   EST_REPS=100
##   EST_CORES=7
##   EST_OUT_DIR=retained_methods_estimation_final
##   EST_ZETA_VALUES=0.10,0.30,0.50
##   EST_MAX_CELLS=0       0 runs all cells; a positive value is for smoke tests
##   EST_FORCE=0           1 overwrites completed cell files
##   VSOM_BOOTSTRAP_R=1000 used if detection-stage cutoffs are unavailable

options(stringsAsFactors = FALSE, scipen = 999)

env_int_est <- function(name, default) {
  x <- Sys.getenv(name, unset = "")
  if (identical(x, "")) default else as.integer(x)
}

env_num_vec_est <- function(name, default) {
  x <- Sys.getenv(name, unset = "")
  if (identical(x, "")) default else as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
}

script_arg_est <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
this_file_est <- normalizePath(sub("^--file=", "", script_arg_est[1]),
                               winslash = "/", mustWork = TRUE)
code_dir_est <- dirname(this_file_est)
## Keep all generated results inside this self-contained submission directory.
study_dir_est <- code_dir_est
source_file_est <- file.path(code_dir_est, "simulation_comparison.R")

## Load method implementations without executing the original main block.
expr_est <- parse(source_file_est)
main_start_est <- which(vapply(expr_est, function(e) {
  grepl("^tasks <- make_tasks\\(mode\\)", paste(deparse(e), collapse = " "))
}, logical(1)))[1]
if (!is.finite(main_start_est)) stop("Could not locate the comparison main block")
for (j in seq_len(main_start_est - 1L)) eval(expr_est[[j]], envir = .GlobalEnv)

reps_est <- env_int_est("EST_REPS", 100L)
detected_est <- parallel::detectCores(logical = FALSE)
if (!is.finite(detected_est)) detected_est <- 2L
cores_est <- max(1L, min(env_int_est("EST_CORES", max(1L, detected_est - 1L)),
                         detected_est))
max_cells_est <- env_int_est("EST_MAX_CELLS", 0L)
force_est <- env_int_est("EST_FORCE", 0L) == 1L
vsom_R_est <- env_int_est("VSOM_BOOTSTRAP_R", 1000L)
out_name_est <- Sys.getenv("EST_OUT_DIR", unset = "retained_methods_estimation_final")
out_dir_est <- file.path(study_dir_est, out_name_est)
cell_dir_est <- file.path(out_dir_est, "cells")
dir.create(cell_dir_est, recursive = TRUE, showWarnings = FALSE)

predictors_est <- c("iid", "ar", "ma")
errors_est <- c("normal", "mixnorm")
T_est <- 100L
rho_est <- 0.20
zeta_est <- env_num_vec_est("EST_ZETA_VALUES", c(0.10, 0.30, 0.50))
method_order_est <- c(
  "Raw", "Proposed", "Oracle LM--BH", "Pooled LM--BH",
  "Batch Conformal--BH (Oracle)", "Batch Conformal--BH (Screened)",
  "VSOM", "Cook's distance", "DFFITS"
)

design_est <- expand.grid(
  predictor = predictors_est, error = errors_est, T = T_est, rho = rho_est,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
design_est <- design_est[order(match(design_est$error, errors_est),
                               match(design_est$predictor, predictors_est),
                               design_est$rho), ]
rownames(design_est) <- NULL
design_est$cell <- seq_len(nrow(design_est))
if (max_cells_est > 0L) design_est <- head(design_est, max_cells_est)

raw_estimation_row <- function(pre, dat) {
  beta_raw <- pooled_joint_beta(pre, seq_len(pre$n))
  err <- beta_raw - dat$beta
  data.frame(
    method = "Raw", selected = 0L, retained = pre$n, fallback = 0L,
    TP = 0L, FP = 0L, TPR = NA_real_, FPR = NA_real_, FDR = NA_real_,
    RMSE = sqrt(mean(err^2)), MAE = mean(abs(err)),
    reference_purity = NA_real_, reference_clean_recall = NA_real_,
    stringsAsFactors = FALSE
  )
}

run_estimation_replication <- function(task) {
  dat <- generate_panel(task$scenario, task$zeta, task$seed)
  pre <- prepare_panel(dat)
  dat_ref <- generate_panel(task$scenario, task$zeta, task$reference_seed)
  pre_ref <- prepare_panel(dat_ref)
  all_idx <- seq_len(dat$n)
  clean_idx <- which(!dat$is_out)

  prop <- proposed_method(pre)
  oracle <- lm_method(pre, clean_idx)
  pooled <- lm_method(pre, all_idx)
  bc_oracle <- batch_conformal_method(
    pre, pre_ref, dat_ref, reference = "oracle",
    seed = task$reference_seed + 1000L
  )
  bc_screened <- batch_conformal_method(
    pre, pre_ref, dat_ref, reference = "screened",
    seed = task$reference_seed + 2000L
  )
  vsom <- vsom_method(pre, task$vsom_cutoff)
  infl <- influence_methods(pre)

  rows <- rbind(
    raw_estimation_row(pre, dat),
    evaluate_method("Proposed", prop, pre, dat),
    evaluate_method("Oracle LM--BH", oracle, pre, dat),
    evaluate_method("Pooled LM--BH", pooled, pre, dat),
    evaluate_method("Batch Conformal--BH (Oracle)", bc_oracle, pre, dat),
    evaluate_method("Batch Conformal--BH (Screened)", bc_screened, pre, dat),
    evaluate_method("VSOM", vsom, pre, dat),
    evaluate_method("Cook's distance", infl$Cook, pre, dat),
    evaluate_method("DFFITS", infl$DFFITS, pre, dat)
  )
  rows$predictor <- task$scenario$predictor
  rows$error <- task$scenario$error
  rows$T <- task$scenario$T
  rows$rho <- task$scenario$rho
  rows$zeta <- task$zeta
  rows$replication <- task$replication
  rows$seed <- task$seed
  rows
}

summarize_estimation <- function(raw) {
  group_names <- c("predictor", "error", "T", "rho", "zeta", "method")
  keys <- unique(raw[group_names])
  metrics <- c("RMSE", "MAE", "selected", "retained", "fallback")
  out <- vector("list", nrow(keys))
  for (j in seq_len(nrow(keys))) {
    take <- rep(TRUE, nrow(raw))
    for (v in group_names) take <- take & raw[[v]] == keys[[v]][j]
    d <- raw[take, , drop = FALSE]
    row <- keys[j, , drop = FALSE]
    for (v in metrics) {
      x <- d[[v]]
      row[[paste0(v, "_mean")]] <- mean(x, na.rm = TRUE)
      row[[paste0(v, "_sd")]] <- sd(x, na.rm = TRUE)
    }
    row$replications <- length(unique(d$replication))
    row$failures <- sum(!is.finite(d$RMSE) | !is.finite(d$MAE))
    out[[j]] <- row
  }
  ans <- do.call(rbind, out)
  ans <- ans[order(ans$rho, ans$zeta, match(ans$method, method_order_est),
                   match(ans$error, errors_est),
                   match(ans$predictor, predictors_est)), ]
  rownames(ans) <- NULL
  ans
}

message(sprintf(
  "Aligned estimation comparison: %d cells, %d replications, %d workers",
  nrow(design_est), reps_est, cores_est
))
message("zeta = ", paste(format(zeta_est, nsmall = 2), collapse = ", "))

cutoff_file_est <- file.path(study_dir_est, "retained_methods_detection_final",
                             "vsom_cutoffs.csv")
if (file.exists(cutoff_file_est)) {
  cutoff_data_est <- read.csv(cutoff_file_est, stringsAsFactors = FALSE)
  vsom_cutoff_est <- setNames(cutoff_data_est$cutoff,
                              cutoff_data_est$predictor)
} else {
  vsom_cutoff_est <- setNames(numeric(length(predictors_est)), predictors_est)
  for (pp in predictors_est) {
    scenario_vsom <- list(
      name = paste0("VSOM-null-", pp), n = 200L, T = T_est, K = 4L, m = 4L,
      rho = 0, predictor = pp, error = "normal", direction = "aligned",
      signal_lower = 0.10
    )
    pp_code <- match(pp, predictors_est)
    dat_vsom <- generate_panel(scenario_vsom, zeta = 0,
                               seed = 202608220L + pp_code * 10000L)
    pre_vsom <- prepare_panel(dat_vsom)
    vsom_cutoff_est[[pp]] <- vsom_bootstrap_cutoff(
      pre_vsom, R = vsom_R_est, alpha = 0.05,
      seed = 202608221L + pp_code * 10000L
    )
  }
}
write.csv(data.frame(
  predictor = names(vsom_cutoff_est), cutoff = as.numeric(vsom_cutoff_est),
  bootstrap_replications = vsom_R_est
), file.path(out_dir_est, "vsom_cutoffs.csv"), row.names = FALSE)
message("VSOM cutoffs: ", paste(sprintf("%s=%.4f", names(vsom_cutoff_est),
                                        vsom_cutoff_est), collapse = ", "))

cl_est <- parallel::makeCluster(cores_est)
on.exit(parallel::stopCluster(cl_est), add = TRUE)
objects_est <- ls(envir = .GlobalEnv, all.names = TRUE)
functions_est <- objects_est[vapply(mget(objects_est, envir = .GlobalEnv),
                                    is.function, logical(1))]
parallel::clusterExport(cl_est, functions_est, envir = .GlobalEnv)

start_all_est <- Sys.time()
for (jj in seq_len(nrow(design_est))) {
  dd <- design_est[jj, ]
  stem <- sprintf("cell_%02d_%s_%s_T%d_rho%02d", dd$cell, dd$predictor,
                  dd$error, dd$T, round(100 * dd$rho))
  rds_path <- file.path(cell_dir_est, paste0(stem, ".rds"))
  csv_path <- file.path(cell_dir_est, paste0(stem, ".csv"))
  if (file.exists(rds_path) && !force_est) {
    message(sprintf("[%02d/%02d] reuse %s", jj, nrow(design_est), stem))
    next
  }

  scenario <- list(
    name = stem, n = 200L, T = as.integer(dd$T), K = 4L, m = 4L,
    rho = dd$rho, predictor = dd$predictor, error = dd$error,
    direction = "aligned", signal_lower = 0.10
  )
  tasks <- vector("list", reps_est * length(zeta_est))
  pos <- 1L
  for (zz in zeta_est) {
    for (rr in seq_len(reps_est)) {
      tasks[[pos]] <- list(
        scenario = scenario, zeta = zz, replication = rr,
        ## Common random numbers across zeta within a design cell.
        seed = 202608210L + dd$cell * 100000L + rr,
        reference_seed = 802608210L + dd$cell * 100000L + rr,
        vsom_cutoff = unname(vsom_cutoff_est[[dd$predictor]])
      )
      pos <- pos + 1L
    }
  }
  tic <- Sys.time()
  pieces <- parallel::parLapplyLB(cl_est, tasks, run_estimation_replication)
  raw_cell <- do.call(rbind, pieces)
  saveRDS(raw_cell, rds_path)
  write.csv(raw_cell, csv_path, row.names = FALSE)
  message(sprintf("[%02d/%02d] saved %s in %.1f minutes", jj,
                  nrow(design_est), stem,
                  as.numeric(difftime(Sys.time(), tic, units = "mins"))))
}

cell_files <- file.path(cell_dir_est, sprintf(
  "cell_%02d_%s_%s_T%d_rho%02d.rds", design_est$cell,
  design_est$predictor, design_est$error, design_est$T,
  round(100 * design_est$rho)
))
missing_est <- cell_files[!file.exists(cell_files)]
if (length(missing_est)) stop("Missing completed cell files: ", paste(missing_est, collapse = ", "))

raw_est <- do.call(rbind, lapply(cell_files, readRDS))
raw_est$method[raw_est$method == "Oracle-reference LM--BH"] <- "Oracle LM--BH"
summary_est <- summarize_estimation(raw_est)
saveRDS(list(raw = raw_est, summary = summary_est, design = design_est,
             zeta = zeta_est, replications = reps_est),
        file.path(out_dir_est, "estimation_comparison.rds"))
write.csv(raw_est, file.path(out_dir_est, "estimation_comparison_raw.csv"),
          row.names = FALSE)
write.csv(summary_est,
          file.path(out_dir_est, "estimation_comparison_summary.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir_est, "sessionInfo.txt"))

message(sprintf("Saved %d raw rows and %d summary rows under %s in %.1f minutes",
                nrow(raw_est), nrow(summary_est), out_dir_est,
                as.numeric(difftime(Sys.time(), start_all_est, units = "mins"))))
