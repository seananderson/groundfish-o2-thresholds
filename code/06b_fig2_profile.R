# Alternative Fig 2: the profile likelihood of the cutpoint, as a violin per species, for the best-AIC oxygen
# metric of each dataset that passes the checks in select_datasets.R. No ensemble over metrics and no draws.
# Run after 05_aic_table.R. Works on its own or after 06.

library(dplyr)
library(readxl)
library(ggplot2)
library(ggimage)
library(stringr)
library(here)

here()

output_folder <- "output"
dir.create(paste0("output/", output_folder, "/plots"), recursive = TRUE, showWarnings = FALSE)
source("code/select_datasets.R") # defines `selection`, one row per dataset

keep <- filter(selection, reason == "pass")
keep_all <- keep # every passing dataset, for the all-regions version at the end

# One dataset per species, as in 06 (pick_dataset() in select_datasets.R): coastwide if it passed
# (IPHC version first for IPHC species), then BC, then the other regions
keep <- bind_rows(lapply(unique(keep$species), function(s) {
  d <- filter(keep, species == s)
  filter(d, data == pick_dataset(d$data, s))
}))
keep$img <- paste0("data/pictures/", gsub(" ", "_", keep$species), ".png")
keep$Species <- str_to_sentence(keep$species)
write.csv(keep, paste0("output/", output_folder, "/fig2_profile_plotted.csv"), row.names = FALSE)

# Each violin is the relative profile likelihood, exp(-delta NLL), across cutpoints (peak = 1), from the
# best-AIC metric's profile. Not a probability density. No point and no interval are drawn.
profile_rows <- function(k, max_dnll = 4) {
  f <- paste0("output/", output_folder, "/", k$species, "_", k$data)
  pr <- readRDS(paste0(f, "_", k$model, ".rds"))
  dat <- readRDS(paste0(f, "_dat.rds"))
  v <- dat[[mi_var[[k$model]]]]
  p <- pr$profile[!is.na(pr$profile$nll), ]
  p <- p[order(p$cut), ]
  x <- p$cut * sd(v, na.rm = TRUE) + mean(v, na.rm = TRUE)
  grid <- seq(min(x), max(x), length.out = 1000)
  dn <- approx(x, p$dnll, grid)$y # dnll is relative to the profile minimum
  # Drop cutpoints more than max_dnll NLL worse than the best (4: relative likelihood < 0.018), so no hairline tails
  data.frame(Species = k$Species, data = k$data, kpa = grid, rel = ifelse(dn > max_dnll, 0, exp(-dn)))
}
prof_df <- bind_rows(lapply(seq_len(nrow(keep)), function(i) profile_rows(keep[i, ])))
lev <- keep$Species[order(keep$cut, decreasing = TRUE)]
keep$y <- match(keep$Species, lev)
# Colour each violin by the dataset it comes from
grp_lab <- function(data) {
  ifelse(data == "coastwide _iphc", "Coastwide, trawl & IPHC",
    ifelse(data == "coastwide", "Coastwide, trawl",
      ifelse(data == "bc", "British Columbia, trawl", data))) # any other dataset keeps its own name
}
grp_cols <- c("Coastwide, trawl" = "#332288", "Coastwide, trawl & IPHC" = "#AA4499",
              "British Columbia, trawl" = "#44AA99")
prof_df$group <- factor(grp_lab(prof_df$data), levels = names(grp_cols))
prof_df$y <- match(prof_df$Species, lev)
prof_df$rel[!is.finite(prof_df$rel)] <- 0

xlim_v <- c(min(prof_df$kpa[prof_df$rel > 0]) - 0.3, max(prof_df$kpa[prof_df$rel > 0]) + 0.3)
prof_df <- prof_df[prof_df$kpa >= xlim_v[1] & prof_df$kpa <= xlim_v[2], ]
x_img <- xlim_v[1] - 0.1 * diff(xlim_v)
g2 <- ggplot() +
  scale_y_continuous(breaks = seq_along(lev), labels = lev, expand = expansion(add = 0.6)) +
    geom_image(data = keep, aes(y = y, x = x_img, image = img), color = "grey25", size = 0.07) +
      scale_x_continuous(breaks = seq(2.5, 10, by = 2.5)) +
        geom_ribbon(data = prof_df, aes(x = kpa, ymin = y - 0.42 * rel, ymax = y + 0.42 * rel, group = Species,
          fill = group), alpha = 0.75, colour = NA) +
        scale_fill_manual(values = grp_cols, drop = TRUE, name = NULL) +
        coord_cartesian(xlim = xlim_v, expand = FALSE, clip = "off") +
          theme_bw(base_size = 20) +
  theme(panel.grid.major = element_line(colour = "grey88", linewidth = 0.45),
        panel.border = element_rect(colour = "grey45", fill = NA, linewidth = 0.4),
        axis.ticks = element_line(colour = "grey45", linewidth = 0.4),
        panel.grid.minor = element_blank(), strip.background = element_blank(),
        axis.title.x = element_text(size = 16, hjust = 1), # right-aligned so it fits at the narrower width
        axis.text.y = element_text(margin = margin(r = 40)), # room for the silhouettes between ticks and names
        axis.text = element_text(color = "black"),
        legend.position = "top", legend.justification = "right", legend.direction = "vertical",
        legend.text = element_text(size = 13), legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(0.9, "lines")) +
  xlab("Temperature-Adjusted Oxygen Threshold (kPa)") +
  ylab("")

ggsave(paste0("output/", output_folder, "/plots/Fig2_profile_violin.png"), plot = g2,
       width = 7, height = 7.3, units = "in", bg = "white", dpi = 600)

# All-regions version: every dataset that passes, overlaid in its species' row, coloured by region
# (same colours as Fig S3 in 06), with a dashed outline where the dataset includes IPHC data.
# Each violin uses that dataset's own best-AIC metric. Rows are in the same order as above.
# Each violin is cut at the 95% profile cutoff, so it spans that dataset's 95% profile interval
# (two pieces where the interval is in two pieces). Square-root x axis to fit the GOA and EBS fits.
max_dnll_regions <- qchisq(0.95, 1) / 2
keep_all$Species <- str_to_sentence(keep_all$species)
prof_all <- bind_rows(lapply(seq_len(nrow(keep_all)), function(i) profile_rows(keep_all[i, ], max_dnll_regions)))
prof_all$rel[!is.finite(prof_all$rel)] <- 0
prof_all$y <- match(prof_all$Species, lev)
prof_all$region <- factor(gsub(" _iphc", "", prof_all$data), levels = c("ebs", "goa", "bc", "cc", "coastwide"))
prof_all$survey <- ifelse(grepl("iphc", prof_all$data), "Bottom trawl & IPHC", "Bottom trawl only")
reg_cols <- c(ebs = "#88CCEE", goa = "#999933", bc = "#44AA99", cc = "#CC6677", coastwide = "#332288")
reg_labs <- c(ebs = "Eastern Bering Sea", goa = "Gulf of Alaska", bc = "British Columbia",
              cc = "California Current", coastwide = "Coastwide")
write.csv(keep_all, paste0("output/", output_folder, "/fig2_profile_all_regions_plotted.csv"), row.names = FALSE)

xr <- range(prof_all$kpa[prof_all$rel > 0])
xlim_a <- c(0.95 * xr[1], 1.05 * xr[2])
prof_all <- prof_all[prof_all$kpa >= xlim_a[1] & prof_all$kpa <= xlim_a[2], ]
# NA (not 0) outside each violin, so the outlines break there instead of running along the row
prof_all$rel[prof_all$rel == 0] <- NA
# Silhouettes 10% of the panel width left of it, on the square-root scale
x_img_a <- max(sqrt(xlim_a[1]) - 0.1 * diff(sqrt(xlim_a)), 0)^2
g3 <- ggplot() +
  scale_y_continuous(breaks = seq_along(lev), labels = lev, expand = expansion(add = 0.6)) +
  geom_image(data = keep, aes(y = y, x = x_img_a, image = img), color = "grey25", size = 0.07) +
  scale_x_sqrt(breaks = c(1, 2, 5, 10, 15, 20, 30, 40), labels = function(x) format(x, drop0trailing = TRUE, trim = TRUE)) +
  geom_ribbon(data = prof_all, aes(x = kpa, ymin = y - 0.42 * rel, ymax = y + 0.42 * rel,
    group = interaction(Species, data), fill = region, colour = region, linetype = survey),
    alpha = 0.3, linewidth = 0.35) +
  scale_fill_manual(values = reg_cols, labels = reg_labs, drop = TRUE, name = NULL) +
  scale_colour_manual(values = reg_cols, labels = reg_labs, drop = TRUE, name = NULL) +
  scale_linetype_manual(values = c("Bottom trawl & IPHC" = "22", "Bottom trawl only" = "solid"), name = NULL) +
  coord_cartesian(xlim = xlim_a, expand = FALSE, clip = "off") +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_line(colour = "grey88", linewidth = 0.45),
        panel.border = element_rect(colour = "grey45", fill = NA, linewidth = 0.4),
        axis.ticks = element_line(colour = "grey45", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 16, hjust = 1),
        axis.text.y = element_text(margin = margin(r = 40)),
        axis.text = element_text(color = "black"),
        legend.position = "top", legend.box = "vertical", legend.justification = "right",
        legend.text = element_text(size = 12), legend.margin = margin(0, 0, 0, 0),
        legend.spacing.y = unit(2, "pt"), legend.key.size = unit(0.9, "lines")) +
  guides(fill = guide_legend(nrow = 2), colour = guide_legend(nrow = 2),
         linetype = guide_legend(override.aes = list(fill = NA, colour = "grey20"))) +
  xlab("Temperature-Adjusted Oxygen Threshold (kPa)") +
  ylab("")

ggsave(paste0("output/", output_folder, "/plots/Fig2_profile_violin_regions.png"), plot = g3,
       width = 7.5, height = 7.5, units = "in", bg = "white", dpi = 600)
