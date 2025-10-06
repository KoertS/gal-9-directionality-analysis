library(dplyr)

rotate <- function(x, y, angle, ox, oy) {
  qx <- ox + cos(angle) * (x - ox) - sin(angle) * (y - oy)
  qy <- oy + sin(angle) * (x - ox) + cos(angle) * (y - oy)
  cbind(qx, qy)
}

align_tracks <- function(tracks, borders) {
  transform_df <- borders %>%
    mutate(y2 = y2 - y1, angle = atan2(y2, x2), rot_by = (2 * pi - angle)) %>%
    select(donor, condition, x1, y1, rot_by)
  tracks %>%
    merge(transform_df) %>%
    mutate(
      x = rotate(x, y, rot_by, x1, y1)[, 1] - x1,
      y = rotate(x, y, rot_by, x1, y1)[, 2] - y1,
    ) %>%
    select(id, t, donor, condition, x, y)
}
