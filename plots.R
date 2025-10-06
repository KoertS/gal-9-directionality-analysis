library(ggplot2)
library(dplyr)
library(patchwork)
library(ggbeeswarm)
library(see)

my_theme <- theme_bw(base_size = 9) +
  theme(
    strip.background = element_blank(),
    text = element_text(size = 9),
    strip.text = element_text(size = 9),
    panel.spacing = unit(1, "lines")
  )
theme_set(my_theme)

set_factor_levels <- function(df) {
  mutate(df, condition = factor(condition, levels = c("- CCL21", "WT", "gal9 KD", "gal9 KD-WT")))
}

plot_tracks <- function(tracks, borders) {
  tracks %>%
    set_factor_levels() %>%
    ggplot(aes(x = x, y = y, color = as.factor(id))) +
    geom_path(show.legend = FALSE) +
    geom_segment(data = set_factor_levels(borders), aes(
      x = x1, y = y1,
      xend = x2, yend = y2
    ), color = "black") +
    facet_grid(donor ~ condition, scales = "free") +
    theme(panel.grid = element_blank(), panel.border = element_blank())
}

plot_aligned_tracks <- function(aligned_tracks) {
  aligned_tracks %>%
    set_factor_levels() %>%
    ggplot(aes(x = x, y = y, color = as.factor(id))) +
    geom_path(show.legend = FALSE) +
    geom_hline(yintercept = 0, color = "black") +
    facet_grid(donor ~ condition) +
    theme(panel.grid = element_blank(), panel.border = element_blank())
}

plot_difference_KD_WT_in_donor <- function(angles, diffboot, donorName) {
  data <- filter(angles, donor == donorName)
  dboots_don <- filter(diffboot$dboots, donor == donorName)
  confs_don <- filter(diffboot$confs, donor == donorName)

  means <- data %>%
    group_by(condition) %>%
    summarize(mean_angle = mean(angle)) %>%
    cbind(xend = c(4.5, 4.5, 4)) %>%
    set_factor_levels()

  p <- data %>%
    set_factor_levels() %>%
    ggplot(aes(x = as.numeric(condition), y = angle, color = condition), show.legend = FALSE) +
    geom_quasirandom(show.legend = FALSE, alpha = .3) +
    geom_violinhalf(
      data = dboots_don,
      aes(x = 4, y = y + estimate_WT, fill = x),
      color = NA,
      alpha = .5, flip = TRUE, show.legend = FALSE
    ) +
    geom_pointrange(data = confs_don, x = 4, size = .75, color = "black", aes(
      y = estimate + estimate_WT,
      ymin = ci_low + estimate_WT,
      ymax = ci_high + estimate_WT
    )) +
    labs(y = "absolute angle to y-axis (degrees)", x = NULL) +
    scale_y_continuous(
      limits = c(0, 180), breaks = c(0, 90, 180), expand = c(0, 0),
      sec.axis = sec_axis(~ . - means$mean_angle[means$condition == "WT"],
        name = "difference in mean angle (KD - WT)",
        breaks = seq(-45, 45, 15)
      )
    ) +
    geom_segment(
      data = means, aes(y = mean_angle, xend = xend),
      lty = 2, show.legend = FALSE
    ) +
    scale_color_manual(values = c("WT" = "black", "gal9 KD" = "red")) +
    annotate("text",
      x = 4.4, y = means$mean_angle[means$condition == "WT"], angle = 90,
      hjust = -.1, label = "less directional", alpha = .5
    ) +
    annotate("text",
      x = 4.4, y = means$mean_angle[means$condition == "WT"], angle = 90,
      hjust = 1.1, label = "more directional", alpha = .5
    ) +
    scale_x_continuous(
      breaks = 1:4,
      labels = c("- CCL21", "WT", "gal9 KD", "gal9 KD-WT")
    ) +
    facet_wrap(~donor) +
    theme(panel.grid = element_blank())
  p
}

plot_difference_KD_WT_all_donors <- function(angles, diffboot) {
  donors <- unique(angles$donor)
  plots_donors <- lapply(donors, function(don) {
    plot_difference_KD_WT_in_donor(
      angles = angles,
      donorName = don,
      diffboot = diffboot
    ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  })
  wrap_plots(plots_donors, ncol = 4, axis_titles = "collect")
}

polar_hist <- function(angles) {
  angles %>%
    set_factor_levels() %>%
    ggplot(aes(x = angle)) +
    geom_histogram(
      breaks = seq(-180, 180, 15)
    ) +
    coord_radial(expand = FALSE) +
    scale_x_continuous(breaks = seq(-180, 180, 90)) +
    facet_wrap(~condition) +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      panel.border = element_blank()
    )
}
