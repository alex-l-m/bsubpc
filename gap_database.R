library(tidyverse)

csd_matches <- read_csv('csd_matches.csv', col_types = cols(
    mol_id = col_character(),
    entry = col_character()
))

table1 <- read_csv('table1.csv', col_types = cols(
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

entry_gap <- csd_matches |>
    left_join(table1, by = 'mol_id') |>
    mutate(gap = ifelse(!is.na(toluene_gap), toluene_gap,
            ifelse(!is.na(dcm_gap), dcm_gap, dmso_gap))) |>
    select(entry, gap) |>
    rename(mol_id = entry)

write_csv(entry_gap, "entry_gap.csv")
