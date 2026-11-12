# Dataset selection from the profile fits: for each dataset, the oxygen metric (model 13, 14 or 15) with the
# lowest AIC, its cutpoint and profile-likelihood interval, and whether the dataset passes the checks below.
# Sourced by 06 and 06b. Run after 05_aic_table.R. Writes output/<output_folder>/fig2_profile_all_datasets.csv.
#
# A dataset is kept if (the same checks as filter_reason_direct() in code/breakpoint_refit/21_ensemble_funs.R):
#   - a metabolic-index model has the lowest AIC of all models fitted to that dataset
#   - the profile interval is closed
#   - each side of the cutpoint has at least min_tows tows or min_prop of the tows
#   - no more than max_prop_below of the tows are below the cutpoint
#   - delta AIC against the best no-oxygen model is at least min_daic (skipped if neither model 7 nor 8 passed)

library(dplyr)
library(readxl)
library(here)

here()

if(!exists("output_folder")) output_folder <- "output"

# Selection rule settings
min_tows <- 100
min_prop <- 0.1
min_daic <- 4
max_prop_below <- 0.9

mi_models <- c("model13", "model14", "model15")
mi_var <- c(model13 = "mi1", model14 = "mi2", model15 = "mi3")
species_iphc <- c("sablefish", "pacific cod", "pacific halibut", "yelloweye rockfish", "longnose skate", "big skate", "spiny dogfish", "rougheye rockfish")

aic_sel <- as.data.frame(readRDS(paste0("output/", output_folder, "/aic_table.rds")))
aic_sel[, 1:(ncol(aic_sel) - 5)] <- sapply(aic_sel[, 1:(ncol(aic_sel) - 5)], as.numeric) # delta AIC; NA = failed checks

best_one <- function(i) {
  row <- aic_sel[i, ]
  ms <- mi_models[!is.na(unlist(row[mi_models]))]
  out <- data.frame(species = row$species, data = row$`data type`, reason = NA_character_)
  if (!length(ms)) { out$reason <- "no MI model passed checks"; return(out) }
  best <- ms[which.min(unlist(row[ms]))]
  out$model <- best
  f <- paste0("output/", output_folder, "/", row$species, "_", row$`data type`)
  pr <- try(readRDS(paste0(f, "_", best, ".rds")), silent = TRUE)
  dat <- try(readRDS(paste0(f, "_dat.rds")), silent = TRUE)
  if (inherits(pr, "try-error") || inherits(dat, "try-error") || is.null(pr$estimate)) {
    out$reason <- "model or data file missing"
    return(out)
  }
  e <- pr$estimate
  v <- dat[[mi_var[[best]]]]
  unscale <- function(x) x * sd(v, na.rm = TRUE) + mean(v, na.rm = TRUE) # same scaling as 04a/06
  # Interval: the contiguous stretch of the profile around the best cutpoint where the NLL is within
  # qchisq(level, 1) / 2 of the minimum. The fit's own lwr/upr span every piece of the allowed set, so
  # when the profile has a second dip (e.g., a steep step on a few low-oxygen tows) they cross a gap.
  # Edges are linear interpolations between profile points, as in sdmTMB; a stretch that reaches the
  # end of the profile is open (NA).
  prof <- pr$profile[!is.na(pr$profile$nll), ]
  prof <- prof[order(prof$cut), ]
  within <- stats::qchisq(if (is.null(pr$level)) 0.95 else pr$level, 1) / 2
  i0 <- which.min(abs(prof$cut - e$cut))
  lo <- hi <- i0
  while (lo > 1 && prof$dnll[lo - 1] < within) lo <- lo - 1
  while (hi < nrow(prof) && prof$dnll[hi + 1] < within) hi <- hi + 1
  cross <- function(i, j) prof$cut[i] + (within - prof$dnll[i]) * (prof$cut[j] - prof$cut[i]) / (prof$dnll[j] - prof$dnll[i])
  lwr_s <- if (lo > 1) cross(lo, lo - 1) else NA_real_
  upr_s <- if (hi < nrow(prof)) cross(hi, hi + 1) else NA_real_
  out$cut <- unscale(e$cut)
  out$lwr <- unscale(lwr_s)
  out$upr <- unscale(upr_s)
  out$lwr_fit <- unscale(e$lwr) # the fit's own interval, spanning all pieces
  out$upr_fit <- unscale(e$upr)
  out$pieces <- e$pieces
  out$n <- e$n_below + e$n_above
  out$n_below <- e$n_below
  out$prop_below <- e$prop_below
  # delta AIC is relative to the best model overall, so the best MI model's value is the gap to the best no-oxygen model
  no_o2 <- unlist(row[c("model7", "model8")])
  no_o2 <- no_o2[!is.na(no_o2)]
  out$daic_vs_no_o2 <- if (length(no_o2)) min(no_o2) - min(unlist(row[ms])) else NA_real_
  n_above <- out$n - out$n_below
  out$reason <-
    if (min(unlist(row[ms])) > 0) "MI model not best by AIC"
    else if (!is.finite(out$lwr) || !is.finite(out$upr)) "open profile interval"
    else if (out$n_below < min_tows && out$n_below / out$n < min_prop) "too few tows below cutpoint"
    else if (n_above < min_tows && n_above / out$n < min_prop) "too few tows above cutpoint"
    else if (out$prop_below > max_prop_below) "over 90% of tows below cutpoint"
    else if (!is.na(out$daic_vs_no_o2) && out$daic_vs_no_o2 < min_daic) "delta AIC below 4"
    else "pass"
  out
}

selection <- bind_rows(lapply(seq_len(nrow(aic_sel)), best_one))
selection$region <- gsub(" _iphc", "", selection$data)
selection$iphc <- grepl("iphc", selection$data)
write.csv(selection, paste0("output/", output_folder, "/fig2_profile_all_datasets.csv"), row.names = FALSE)
print(table(selection$reason))

