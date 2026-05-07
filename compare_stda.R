library(tidyverse)
library(tidymodels)
library(cowplot)
theme_set(theme_cowplot() + theme(plot.background = element_rect(fill = 'white')))
library(ggrepel)
library(robustbase)

paper_table <- read_csv('table1.csv', col_types = cols(
    mol_id = col_character(),
    film_homo = col_double(),
    avg_film_homo = col_double(),
    toluene_homo = col_double(),
    toluene_gap = col_double(),
    dcm_homo = col_double(),
    dcm_gap = col_double(),
    dmso_homo = col_double(),
    dmso_gap = col_double()
))

# Select whichever solvent isn't missing
experiment_homolumo <- paper_table |>
    group_by(mol_id) |>
    transmute(
        homo = ifelse(!is.na(toluene_homo), toluene_homo,
                      ifelse(!is.na(dcm_homo), dcm_homo, dmso_homo)),
        gap = ifelse(!is.na(toluene_gap), toluene_gap,
                      ifelse(!is.na(dcm_gap), dcm_gap, dmso_gap))
    ) |>
    ungroup()

csd_matches <- read_csv('csd_matches.csv', col_types = cols(
    mol_id = col_character(),
    entry = col_character()
))
all_crystal_stda <- read_csv('crystal_stda.csv.gz', col_types = cols(
    mol_id = col_character(),
    energy = col_double()
)) |>
    # Renaming because these are not the molecule IDs in the paper, they're CSD
    # codes
    rename(entry = mol_id)

# STDA predicted emission energy of the CSD entries matching each molecule in
# the dataset
crystal_geom_stda <- csd_matches |>
    inner_join(all_crystal_stda, by = 'entry') |>
    select(mol_id, energy)

# Combine into a comparison table
predicted_energies <- crystal_geom_stda |>
    rename(stda_excitation_energy = energy)
comparison_table <- experiment_homolumo |>
    inner_join(predicted_energies, by = 'mol_id')

stda_comparison_plot <- comparison_table |>
    ggplot(aes(y = stda_excitation_energy, x = gap, label = mol_id)) +
    geom_point() +
    geom_smooth(method = lmrob, se = FALSE) +
    geom_label_repel() +
    # Axis labels
    xlab('Experimental HOMO-LUMO gap (eV)') +
    ylab('sTDA predicted excitation energy (eV)')
    
ggsave('stda_comparison_plot.png', stda_comparison_plot, width = unit(12, 'in'), height = unit(6, 'in'))
