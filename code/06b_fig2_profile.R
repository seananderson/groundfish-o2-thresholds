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
  xlab("Temperature-Adjusted Oxygen Threshold (kPa)\nBest-AIC metric, profile-likelihood interval") +
  ylab("")

ggsave(paste0("output/", output_folder, "/plots/Fig2_profile.png"), plot = g,
       width = 9.5, height = 8.5, units = "in", bg = "white", dpi = 600)
