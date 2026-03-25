library(tidyverse)
library(cowplot)

theme_set(
  theme_cowplot() +
    theme(
      plot.background = element_rect(fill = 'white', color = NA)
    )
)

powerpoint_half_width <- 5.68
powerpoint_height <- 4.76
powerpoint_full_width <- 11.49

input_path <- 'crystal_bowl_depths.csv'
output_path <- 'crystal_bowl_depth_histograms.png'

fd_binwidth <- function(x)
{
  x <- x[!is.na(x)]
  if (length(x) < 2) return(1)

  bw <- 2 * IQR(x) / (length(x)^(1 / 3))
  if (is.finite(bw) && bw > 0) return(bw)

  fallback <- diff(range(x)) / 30
  if (is.finite(fallback) && fallback > 0) return(fallback)

  return(1)
}

bowl_depths <- read_csv(
  input_path,
  col_types = cols(
    mol_id = col_character(),
    plane = col_character(),
    plane_label = col_character(),
    bowl_depth = col_double()
  )
) |>
  mutate(
    plane = factor(
      plane,
      levels = c(
        'outer_terminal_carbon',
        'imine_nitrogen',
        'pyrrole_nitrogen'
      )
    ),
    plane_label = factor(
      plane_label,
      levels = c(
        'outer terminal carbons',
        'three imine nitrogens',
        'three pyrrole nitrogens'
      )
    )
  )

bowl_depth_histograms_plt <- bowl_depths |>
  ggplot(aes(x = bowl_depth)) +
  geom_histogram(binwidth = fd_binwidth) +
  facet_wrap(~plane_label, scales = 'free') +
  labs(
    x = 'Boron-to-plane distance (Angstrom)',
    y = 'Count'
  )

ggsave(
  output_path,
  bowl_depth_histograms_plt,
  width = powerpoint_full_width,
  height = powerpoint_height,
  units = 'in'
)

# Top and bottom 5 molecules for each measure
ranked <- bowl_depths |>
  group_by(plane) |>
  arrange(bowl_depth) |>
  slice(c(1:5, (n() - 4):n())) |>
  ungroup()
write_csv(ranked, 'crystal_bowl_depths_ranked.csv')
