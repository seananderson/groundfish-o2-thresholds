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

xlim_v <- c(min(prof_df$kpa[prof_df$rel > 0]) - 0.3, max(prof_df$kpa[prof_df$rel > 0]) + 0.3)
prof_df <- prof_df[prof_df$kpa >= xlim_v[1] & prof_df$kpa <= xlim_v[2], ]
x_img <- xlim_v[1] - 0.1 * diff(xlim_v)
g2 <- ggplot() +
  scale_y_continuous(breaks = seq_along(lev), labels = lev, expand = expansion(add = 0.6)) +
    geom_image(data = keep, aes(y = y, x = x_img, image = img), color = "grey25", size = 0.07) +
      scale_x_continuous(breaks = seq(2.5, 10, by = 2.5)) +
        geom_ribbon(data = prof_df, aes(x = kpa, ymin = y - 0.42 * rel, ymax = y + 0.42 * rel, group = Species),
        fill = "#6666cc", alpha = 0.65, colour = NA) +
        coord_cartesian(xlim = xlim_v, expand = FALSE, clip = "off") +
          theme_bw(base_size = 20) +
  theme(panel.grid.major = element_line(colour = "grey88", linewidth = 0.45),
        panel.border = element_rect(colour = "grey45", fill = NA, linewidth = 0.4),
        axis.ticks = element_line(colour = "grey45", linewidth = 0.4),
        panel.grid.minor = element_blank(), strip.background = element_blank(),
        axis.title.x = element_text(size = 16, hjust = 1), # right-aligned so it fits at the narrower width
        axis.text.y = element_text(margin = margin(r = 40)), # room for the silhouettes between ticks and names
        axis.text = element_text(color = "black")) +
  xlab("Temperature-Adjusted Oxygen Threshold (kPa)") +
  ylab("")

ggsave(paste0("output/", output_folder, "/plots/Fig2_profile_violin.png"), plot = g2,
       width = 7, height = 6.5, units = "in", bg = "white", dpi = 600)
