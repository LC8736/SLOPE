#!/usr/bin/env Rscript

## Formal comparison under the aligned main simulation design in Section 4.1.
## The script is resumable at the predictor/error/T/rho cell level.
## Environment variables:
##   EXP_REPS=100
##   EXP_CORES=7
##   EXP_OUT_DIR=retained_methods_detection_final
##   EXP_MAX_CELLS=0       0 runs all cells; a positive value is useful for smoke tests
##   EXP_FORCE=0           1 overwrites completed cell files
##   VSOM_BOOTSTRAP_R=1000 number of parametric bootstrap samples for VSOM

options(stringsAsFactors = FALSE, scipen = 999)

env_int <- function(name, default) {
  x <- Sys.getenv(name, unset = "")
  if (identical(x, "")) default else as.integer(x)
}

script_arg_exp <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
this_file_exp <- normalizePath(sub("^--file=", "", script_arg_exp[1]),
                               winslash = "/", mustWork = TRUE)
code_dir_exp <- dirname(this_file_exp)
## Keep all generated results inside this self-contained submission directory.
study_dir_exp <- code_dir_exp
source_file_exp <- file.path(code_dir_exp, "simulation_comparison.R")

## Load the tested method implementations without executing their original main block.
expr_exp <- parse(source_file_exp)
main_start_exp <- which(vapply(expr_exp, function(e) {
  grepl("^tasks <- make_tasks\\(mode\\)", paste(deparse(e), collapse = " "))
}, logical(1)))[1]
if (!is.finite(main_start_exp)) stop("Could not locate the comparison main block")
for (j in seq_len(main_start_exp - 1L)) eval(expr_exp[[j]], envir = .GlobalEnv)

reps_exp <- env_int("EXP_REPS", 100L)
detected_exp <- parallel::detectCores(logical = FALSE)
if (!is.finite(detected_exp)) detected_exp <- 2L
cores_exp <- max(1L, min(env_int("EXP_CORES", max(1L, detected_exp - 1L)),
                         detected_exp))
max_cells_exp <- env_int("EXP_MAX_CELLS", 0L)
force_exp <- env_int("EXP_FORCE", 0L) == 1L
vsom_R_exp <- env_int("VSOM_BOOTSTRAP_R", 1000L)
out_name_exp <- Sys.getenv("EXP_OUT_DIR", unset = "retained_methods_detection_final")
out_dir_exp <- file.path(study_dir_exp, out_name_exp)
cell_dir_exp <- file.path(out_dir_exp, "cells")
dir.create(cell_dir_exp, recursive = TRUE, showWarnings = FALSE)

predictors_exp <- c("iid", "ar", "ma")
errors_exp <- c("normal", "mixnorm")
T_exp <- 100L
rho_exp <- c(0.10, 0.20)
zeta_exp <- seq(0.10, 0.25, by = 0.05)

design_exp <- expand.grid(
  predictor = predictors_exp, error = errors_exp, T = T_exp, rho = rho_exp,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
design_exp <- design_exp[order(match(design_exp$predictor, predictors_exp),
                               match(design_exp$error, errors_exp),
                               design_exp$T, design_exp$rho), ]
rownames(design_exp) <- NULL
design_exp$cell <- seq_len(nrow(design_exp))
if (max_cells_exp > 0L) {
  design_exp <- head(design_exp, max_cells_exp)
}

run_task_expanded <- function(task) {
  ans <- run_replication(task)
  ans$predictor <- task$scenario$predictor
  ans$error <- task$scenario$error
  ans$T <- task$scenario$T
  ans$rho <- task$scenario$rho
  ans
}

summarize_expanded <- function(raw) {
  group_names <- c("predictor", "error", "T", "rho", "zeta", "method")
  keys <- unique(raw[group_names])
  metrics <- c("TPR", "FPR", "FDR", "RMSE", "MAE", "selected",
               "reference_purity", "reference_clean_recall")
  out <- vector("list", nrow(keys))
  for (j in seq_len(nrow(keys))) {
    take <- rep(TRUE, nrow(raw))
    for (v in group_names) take <- take & raw[[v]] == keys[[v]][j]
    d <- raw[take, , drop = FALSE]
    row <- keys[j, , drop = FALSE]
    for (v in metrics) {
      x <- d[[v]]
      row[[paste0(v, "_mean")]] <- if (all(is.na(x))) NA_real_ else
        mean(x, na.rm = TRUE)
      row[[paste0(v, "_sd")]] <- if (sum(!is.na(x)) <= 1L) NA_real_ else
        sd(x, na.rm = TRUE)
    }
    row$replications <- length(unique(d$replication))
    row$failures <- sum(!is.finite(d$TPR) | !is.finite(d$FPR))
    out[[j]] <- row
  }
  ans <- do.call(rbind, out)
  method_order <- c(
    "Proposed", "Oracle LM--BH", "Pooled LM--BH",
    "Batch Conformal--BH (Oracle)", "Batch Conformal--BH (Screened)",
    "VSOM", "Cook's distance", "DFFITS"
  )
  ans <- ans[order(match(ans$predictor, predictors_exp),
                   match(ans$error, errors_exp), ans$T, ans$rho, ans$zeta,
                   match(ans$method, method_order)), ]
  rownames(ans) <- NULL
  ans
}

message(sprintf(
  "Aligned main-design comparison: %d cells, %d replications, %d workers",
  nrow(design_exp), reps_exp, cores_exp
))
message("zeta = ", paste(format(zeta_exp, nsmall = 2), collapse = ", "))

## VSOM follows the observation-level bootstrap rule of Ismadyaliana et al.
## For each predictor process, we generate an independent clean null panel and
## take the median of the bootstrap 95th percentiles of squared standardized
## residuals. The resulting scale-free cutoff is reused across cells having
## the same predictor process.
vsom_cutoff_exp <- setNames(numeric(length(predictors_exp)), predictors_exp)
for (pp in predictors_exp) {
  scenario_vsom <- list(
    name = paste0("VSOM-null-", pp), n = 200L, T = T_exp, K = 4L, m = 4L,
    rho = 0, predictor = pp, error = "normal", direction = "aligned",
    signal_lower = 0.10
  )
  pp_code <- match(pp, predictors_exp)
  dat_vsom <- generate_panel(scenario_vsom, zeta = 0,
                             seed = 202608220L + pp_code * 10000L)
  pre_vsom <- prepare_panel(dat_vsom)
  vsom_cutoff_exp[[pp]] <- vsom_bootstrap_cutoff(
    pre_vsom, R = vsom_R_exp, alpha = 0.05,
    seed = 202608221L + pp_code * 10000L
  )
}
vsom_cutoff_table_exp <- data.frame(
  predictor = names(vsom_cutoff_exp), cutoff = as.numeric(vsom_cutoff_exp),
  bootstrap_replications = vsom_R_exp
)
write.csv(vsom_cutoff_table_exp, file.path(out_dir_exp, "vsom_cutoffs.csv"),
          row.names = FALSE)
message("VSOM cutoffs: ", paste(sprintf("%s=%.4f", names(vsom_cutoff_exp),
                                        vsom_cutoff_exp), collapse = ", "))

cl_exp <- parallel::makeCluster(cores_exp)
on.exit(parallel::stopCluster(cl_exp), add = TRUE)
objects_exp <- ls(envir = .GlobalEnv, all.names = TRUE)
functions_exp <- objects_exp[vapply(mget(objects_exp, envir = .GlobalEnv),
                                    is.function, logical(1))]
parallel::clusterExport(cl_exp, functions_exp, envir = .GlobalEnv)

start_all_exp <- Sys.time()
for (jj in seq_len(nrow(design_exp))) {
  dd <- design_exp[jj, ]
  stem <- sprintf("cell_%02d_%s_%s_T%d_rho%02d", dd$cell, dd$predictor,
                  dd$error, dd$T, round(100 * dd$rho))
  rds_path <- file.path(cell_dir_exp, paste0(stem, ".rds"))
  csv_path <- file.path(cell_dir_exp, paste0(stem, ".csv"))
  if (file.exists(rds_path) && !force_exp) {
    message(sprintf("[%02d/%02d] reuse %s", jj, nrow(design_exp), stem))
    next
  }

  scenario_name <- sprintf("%s--%s--T%d--rho%.2f", toupper(dd$predictor),
                           dd$error, dd$T, dd$rho)
  scenario <- list(
    name = scenario_name, n = 200L, T = as.integer(dd$T), K = 4L, m = 4L,
    rho = dd$rho, predictor = dd$predictor, error = dd$error,
    direction = "aligned", signal_lower = 0.10
  )
  tasks <- vector("list", reps_exp * length(zeta_exp))
  pos <- 1L
  for (zz in zeta_exp) {
    for (rr in seq_len(reps_exp)) {
      tasks[[pos]] <- list(
        scenario = scenario, zeta = zz, replication = rr,
        ## Common random numbers over zeta within a design cell.
        seed = 202610010L + dd$cell * 100000L + rr,
        reference_seed = 802610010L + dd$cell * 100000L + rr,
        vsom_cutoff = unname(vsom_cutoff_exp[[dd$predictor]])
      )
      pos <- pos + 1L
    }
  }
  tic <- Sys.time()
  pieces <- parallel::parLapplyLB(cl_exp, tasks, run_task_expanded)
  raw_cell <- do.call(rbind, pieces)
  saveRDS(raw_cell, rds_path)
  write.csv(raw_cell, csv_path, row.names = FALSE)
  message(sprintf("[%02d/%02d] saved %s in %.1f minutes", jj,
                  nrow(design_exp), stem,
                  as.numeric(difftime(Sys.time(), tic, units = "mins"))))
}

cell_files <- file.path(cell_dir_exp, sprintf(
  "cell_%02d_%s_%s_T%d_rho%02d.rds", design_exp$cell,
  design_exp$predictor, design_exp$error, design_exp$T,
  round(100 * design_exp$rho)
))
missing_exp <- cell_files[!file.exists(cell_files)]
if (length(missing_exp)) stop("Missing completed cells: ", paste(missing_exp, collapse = ", "))
raw_exp <- do.call(rbind, lapply(cell_files, readRDS))
## Older cached cell files may contain the previously explored zeta=0.05
## setting. Retain only the grid declared above when rebuilding combined files.
raw_exp <- raw_exp[raw_exp$zeta %in% zeta_exp, , drop = FALSE]
raw_exp$method[raw_exp$method == "Oracle-reference LM--BH"] <- "Oracle LM--BH"
summary_exp <- summarize_expanded(raw_exp)

saveRDS(list(raw = raw_exp, summary = summary_exp, design = design_exp,
             zeta = zeta_exp, replications = reps_exp),
        file.path(out_dir_exp, "expanded_main_comparison.rds"))
write.csv(raw_exp, file.path(out_dir_exp, "expanded_main_comparison_raw.csv"),
          row.names = FALSE)
write.csv(summary_exp,
          file.path(out_dir_exp, "expanded_main_comparison_summary.csv"),
          row.names = FALSE)
batch_diag_exp <- summary_exp[grepl("Batch Conformal", summary_exp$method),
                              c("predictor", "error", "T", "rho", "zeta",
                                "method", "reference_purity_mean",
                                "reference_purity_sd",
                                "reference_clean_recall_mean",
                                "reference_clean_recall_sd")]
write.csv(batch_diag_exp,
          file.path(out_dir_exp, "batch_reference_diagnostics.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir_exp, "sessionInfo.txt"))

message(sprintf("Completed in %.1f minutes; raw rows=%d, summary rows=%d",
                as.numeric(difftime(Sys.time(), start_all_exp, units = "mins")),
                nrow(raw_exp), nrow(summary_exp)))
