library(dplyr)
library(ggplot2)
library(celltrackR)
library(simpleboot)
library(boot)
source("align_tracks.R")
source("plots.R")

angle_to_negative_y_axis <- function(dx, dy) {
  atan2(dx, -dy) * 180 / pi
}

calculate_angles_to_CCL21 <- function(aligned_tracks) {
  aligned_tracks %>%
    group_by(donor, condition, id) %>%
    summarize(
      dx = last(x) - first(x),
      dy = last(y) - first(y),
      angle = atan2(dx, -dy) * 180 / pi
    ) %>%
    select(-dx, -dy, -id)
}

diffboot_KD_WT <- function(angles) {
  diff_name <- "gal9 KD - WT"
  dboots <- data.frame()
  confs <- data.frame()
  for (don in unique(angles$donor)) {
    angles_KD <- filter(angles, condition == "gal9 KD", donor == don)$angle
    angles_WT <- filter(angles, condition == "WT", donor == don)$angle
    diffmeanboot <- two.boot(angles_KD, angles_WT, mean, R = 100000)
    ci <- boot.ci(diffmeanboot, type = "perc")$percent[4:5]

    dboots <- rbind(
      dboots,
      data.frame(
        x = diff_name,
        y = diffmeanboot$t,
        estimate_WT = mean(angles_WT),
        donor = don,
        condition = diff_name
      )
    )

    confs <- rbind(confs, data.frame(
      ci_low = ci[1],
      ci_high = ci[2],
      estimate = mean(angles_KD) - mean(angles_WT),
      estimate_WT = mean(angles_WT),
      donor = don,
      condition = diff_name
    ))
  }
  list(dboots = dboots, confs = confs)
}


tracks <- read.csv("data/tracks.csv")
borders <- read.csv("data/borders_CCL21.csv")
aligned_tracks <- align_tracks(tracks, borders)

angles_to_CCL21 <- calculate_angles_to_CCL21(aligned_tracks)
angles_all_donors <- mutate(angles_to_CCL21, donor = "all")
abs_angles_to_CCL21 <- angles_to_CCL21 %>% mutate(angle = abs(angle))
abs_angles_all_donors <- mutate(abs_angles_to_CCL21, donor = "all")

diffboot_all <- diffboot_KD_WT(abs_angles_all_donors)
diffboot_per_donor <- diffboot_KD_WT(abs_angles_to_CCL21)


p_tracks <- plot_tracks(tracks, borders)
p_aligned_tracks <- plot_aligned_tracks(aligned_tracks)
p_angles_all_donors <- polar_hist(angles_all_donors)
p_diff_KD_WT <- (plot_difference_KD_WT_in_donor(abs_angles_all_donors, diffboot_all, "all") +
  theme(strip.text.x = element_blank()))
p_diff_KD_WT_per_donor <- plot_difference_KD_WT_all_donors(abs_angles_to_CCL21, diffboot_per_donor)

# Plots for README.md
ggsave(p_tracks,
  file = "figures/tracks.png", create.dir = TRUE,
  width = 60, height = 60, units = "mm", scale = 2
)
ggsave(p_aligned_tracks,
  file = "figures/tracks_aligned.png",
  width = 60, height = 60, units = "mm", scale = 2
)

ggsave(p_angles_all_donors,
  file = "figures/polar_hists.png",
  width = 124,
  height = 45, units = "mm", scale = 1.5
)


# Plots for paper
ggsave(p_angles_all_donors,
  file = "figures/polar_hists.pdf",
  width = 124,
  height = 80, units = "mm", useDingbats = FALSE
)
ggsave(p_diff_KD_WT,
  file = "figures/effect-size.pdf",
  width = 62,
  height = 80, units = "mm", useDingbats = FALSE, scale = 1.5
)
ggsave(p_diff_KD_WT_per_donor,
  file = "figures/effect-size-per-donor.pdf",
  width = 124,
  height = 80, units = "mm", useDingbats = FALSE, scale = 2
)
