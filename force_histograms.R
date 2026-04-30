library(tidyverse)
library(cowplot)
theme_set(theme_cowplot() + theme(plot.background = element_rect(fill = 'white')))

bsubpc_cif_match_indices <- read_csv('bsubpc_cif_match_indices.csv', col_types = cols(
    mol_id = col_character(),
    cif_match_idx = col_double(),
    bsubpc_idx = col_double(),
    bsubpc_label = col_character(),
    cif_atom_idx = col_double(),
    atom_site_label = col_character(),
    element = col_character()
)) |>
    # Get a category by stripping the underscore and digit off the end of the
    # BsubPc template label
    # Only the boron category doesn't have this
    mutate(bsubpc_category = str_replace(bsubpc_label, "_[0-9]+$", ""))

forces <- read_csv('forces.csv', col_types = cols( 
    inpath = col_character(), 
    dirname = col_character(), 
    mol_id = col_character(), 
    frame = col_double(), 
    atom_id = col_double(), 
    element = col_character(), 
    x = col_double(), 
    y = col_double(), 
    z = col_double(), 
    fx = col_double(), 
    fy = col_double(), 
    fz = col_double() 
)) |>
    mutate(`|f|` = sqrt(fx^2 + fy^2 + fz^2)) |>
    left_join(bsubpc_cif_match_indices, by = c('mol_id', atom_id = 'cif_atom_idx', 'element'))

stopifnot(all(forces$frame == 0)) # Only one frame per molecule, so we can ignore the frame column



fd_binwidth <- function(x) {
    2 * IQR(x) / length(x)^(1/3)
}
# Units are eV/angstrom
force_labs <- labs(x = 'Force (eV/angstrom)', y = 'Count')
force_lims <- coord_cartesian(xlim = c(0, 10))
# Facet plot for each element
histograms <- forces |>
    filter(dirname == 'crystal_forces') |>
    filter(element %in% c('B', 'C', 'N', 'O', 'H', 'F')) |>
    ggplot(aes(x = `|f|`)) +
    facet_wrap(~ element, scales = 'free_y') +
    geom_histogram(binwidth = fd_binwidth) +
    force_labs + force_lims +
    ggtitle('Force histograms by element for crystal positions')

ggsave('force_histograms_by_element.png', histograms, width = 10, height = 6)

histograms <- forces |>
    filter(dirname == 'crystal_forces_h_optimized') |>
    filter(element %in% c('B', 'C', 'N', 'O', 'H', 'F')) |>
    ggplot(aes(x = `|f|`)) +
    facet_wrap(~ element, scales = 'free_y') +
    geom_histogram(binwidth = fd_binwidth) +
    force_labs + force_lims +
    ggtitle('Force histograms by element after optimizing H positions')

ggsave('force_histograms_by_element_h_optimized.png', histograms, width = 10, height = 6)

# Now instead of breaking it down by element, breaking it down by the special
# template positions in the bsubpc category column
histograms <- forces |>
    filter(!is.na(bsubpc_category)) |>
    filter(dirname == 'crystal_forces') |>
    ggplot(aes(x = `|f|`)) +
    facet_wrap(~ bsubpc_category, scales = 'free_y') +
    geom_histogram(binwidth = fd_binwidth) +
    force_labs + force_lims +
    ggtitle('Force histograms by template position for crystal positions')

ggsave('force_histograms_by_template_position.png', histograms, width = 10, height = 6)

histograms <- forces |>
    filter(!is.na(bsubpc_category)) |>
    filter(dirname == 'crystal_forces_h_optimized') |>
    ggplot(aes(x = `|f|`)) +
    facet_wrap(~ bsubpc_category, scales = 'free_y') +
    geom_histogram(binwidth = fd_binwidth) +
    force_labs + force_lims +
    ggtitle('Force histograms by template position after optimizing H positions')

ggsave('force_histograms_by_template_position_h_optimized.png', histograms, width = 10, height = 6)
