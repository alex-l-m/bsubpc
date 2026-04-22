library(tidyverse)
library(matsindf)
library(glue)
library(cowplot)
theme_set(theme_cowplot() + theme(plot.background = element_rect(fill = 'white')))

input_path <- 'atom_table.csv'
output_path <- 'crystal_bowl_pca_scores.csv'
atoms <- read_csv(
    input_path,
    col_types = cols(
        mol_id = col_character(),
        atom_id = col_double(),
        atom_name = col_character(),
        symbol = col_character(),
        formal_charge = col_double(),
        x = col_double(),
        y = col_double(),
        z = col_double(),
        bsubpc_idx = col_double(),
        bsubpc_label = col_character(),
        .default = col_guess()
    )
)

bowl_atoms <- atoms |>
    filter(!is.na(bsubpc_idx))

bowl_atom_pos <- bowl_atoms |>
    select(mol_id, bsubpc_idx, bsubpc_label, x, y, z) |>
    pivot_longer(cols = c(x, y, z), names_to = 'axis', values_to = 'coord')

long_pca_embeddings <- bowl_atom_pos |>
    filter(
           # Boron never changes position
           bsubpc_label != 'boron'
           # This I think shouldn't really be necessary... why don't pyrrole
           # nitrogen positions always vary? Symmetry?
#               ! str_detect('pyrrole_nitrogen', bsubpc_label)
    ) |>
    select(-bsubpc_label) |>
    # If I didn't standardize rotations, I need to look just at z
    # I did though... didn't I?
    filter(axis == 'z') |>
    mutate(idx_axis = glue('{bsubpc_idx}_{axis}')) |>
    select(mol_id, idx_axis, coord) |>
    collapse_to_matrices(
        rownames = "mol_id", colnames="idx_axis", matvals = "coord"
    ) |>
    transmute(projection = lapply(coord, function(mat)
        prcomp(mat, center=TRUE, scale=TRUE)$x
    )) |>
    expand_to_tidy(
        rownames = "mol_id", colnames = "component", matvals = "projection"
    )

pca_embeddings <- long_pca_embeddings |>
    pivot_wider(names_from = component, values_from = projection)

pc_plot <- pca_embeddings |>
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point()

ggsave(
    filename = 'crystal_bowl_pca_scores.png',
    plot = pc_plot,
    width = 5, height = 5, units = 'in', dpi = 300
)
