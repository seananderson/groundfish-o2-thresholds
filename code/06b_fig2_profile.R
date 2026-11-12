# Alternative Fig 2: cutpoint and profile-likelihood interval of the best-AIC oxygen metric for each
# dataset that passes the checks in select_datasets.R. No ensemble over metrics and no draws.
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

# One dataset per species, as in 06: coastwide if available, and the IPHC version for IPHC species.
# If a species' coastwide IPHC dataset did not pass, use its coastwide trawl dataset.
keep <- bind_rows(lapply(unique(keep$species), function(s) {
  d <- filter(keep, species == s)
  if (any(d$region == "coastwide")) {
    d <- filter(d, region == "coastwide")
    if (s %in% species_iphc && any(d$iphc)) d <- filter(d, iphc)
  }
  if (nrow(d) > 1) d <- filter(d, region == "bc") # then the BC dataset, as in 06
  d
}))
keep$img <- paste0("data/pictures/", gsub(" ", "_", keep$species), ".png")
keep$Species <- str_to_sentence(keep$species)
write.csv(keep, paste0("output/", output_folder, "/fig2_profile_plotted.csv"), row.names = FALSE)

g <- ggplot(keep, aes(y = reorder(Species, cut, decreasing = TRUE), x = cut)) +
  geom_image(aes(image = img), color = "blue", alpha = 0.3, size = 0.08) +
  geom_linerange(aes(xmin = lwr, xmax = upr), linewidth = 0.4) +
  geom_point(size = 1) +
  theme_bw(base_size = 20) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        axis.text = element_text(color = "black")) +
  xlab("Temperature-Adjusted Oxygen Threshold (kPa)") +
  ylab("")

ggsave(paste0("output/", output_folder, "/plots/Fig2_profile.png"), plot = g,
       width = 9.5, height = 8.5, units = "in", bg = "white", dpi = 600)

# Violin version: each violin is the relative profile likelihood, exp(-delta NLL), across cutpoints
# (peak = 1), from the best-AIC metric's profile. Not a probability density. Violins only, with no
# point or interval; the figure with the best cutpoint and 95% interval is Fig2_profile.png.
prof_df <- bind_rows(lapply(seq_len(nrow(keep)), function(i) {
  k <- keep[i, ]
  f <- paste0("output/", output_folder, "/", k$species, "_", k$data)
  pr <- readRDS(paste0(f, "_", k$model, ".rds"))
  dat <- readRDS(paste0(f, "_dat.rds"))
  v <- dat[[mi_var[[k$model]]]]
  p <- pr$profile[!is.na(pr$profile$nll), ]
  p <- p[order(p$cut), ]
  x <- p$cut * sd(v, na.rm = TRUE) + mean(v, na.rm = TRUE)
  grid <- seq(min(x), max(x), length.out = 1000)
  dn <- approx(x, p$dnll, grid)$y # dnll is relative to the profile minimum
  # Drop cutpoints more than 4 NLL worse than the best (relative likelihood < 0.018), so no hairline tails
  data.frame(Species = k$Species, kpa = grid, rel = ifelse(dn > 4, 0, exp(-dn)))
}))
lev <- keep$Species[order(keep$cut, decreasing = TRUE)]
keep$y <- match(keep$Species, lev)
prof_df$y <- match(prof_df$Species, lev)
prof_df$rel[!is.finite(prof_df$rel)] <- 0

xlim_v <- c(min(prof_df$kpa[prof_df$rel > 0]) - 0.3, max(keep$upr) + 0.5)
# With clipping off (so the silhouettes can sit outside the panel), cut the violins at the axis limits ourselves
prof_df <- prof_df[prof_df$kpa >= xlim_v[1] & prof_df$kpa <= xlim_v[2], ]
# Silhouettes sit just outside the left edge of the panel, between the tick marks and the species names
x_img <- xlim_v[1] - 0.075 * diff(xlim_v)
g2 <- ggplot() +
  geom_ribbon(data = prof_df, aes(x = kpa, ymin = y - 0.42 * rel, ymax = y + 0.42 * rel, group = Species),
              fill = "#6666cc", alpha = 0.65, colour = NA) +
  scale_y_continuous(breaks = seq_along(lev), labels = lev, expand = expansion(add = 0.6)) +
  geom_image(data = keep, aes(y = y, x = x_img, image = img), color = "grey25", size = 0.07) +
  scale_x_continuous(breaks = seq(2.5, 10, by = 2.5)) +
  coord_cartesian(xlim = xlim_v, expand = FALSE, clip = "off") +
  theme_bw(base_size = 20) +
  # Faint grid: one horizontal line per species (breaks are the species rows) and vertical lines at the ticks
  theme(panel.grid.major = element_line(colour = "grey88", linewidth = 0.45),
        panel.border = element_rect(colour = "grey45", fill = NA, linewidth = 0.4),
        axis.ticks = element_line(colour = "grey45", linewidth = 0.4),
        panel.grid.minor = element_blank(), strip.background = element_blank(),
        axis.title.x = element_text(size = 16, hjust = 1), # right-aligned so it fits at the narrower width
        axis.text.y = element_text(margin = margin(r = 48)), # room for the silhouettes between ticks and names
        axis.text = element_text(color = "black")) +
  xlab("Temperature-Adjusted Oxygen Threshold (kPa)") +
  ylab("")

ggsave(paste0("output/", output_folder, "/plots/Fig2_profile_violin.png"), plot = g2,
       width = 8, height = 8.5, units = "in", bg = "white", dpi = 600)
