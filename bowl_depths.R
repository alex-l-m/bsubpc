library(tidyverse)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', args, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub('^--file=', '', file_arg)))
} else {
  getwd()
}

input_path <- file.path(script_dir, 'atom_table.csv')
output_path <- file.path(script_dir, 'crystal_bowl_depths.csv')

plane_definitions <- tribble(
  ~plane, ~plane_label, ~expected_count,
  'outer_terminal_carbon', 'outer terminal carbons', 6L,
  'imine_nitrogen', 'three imine nitrogens', 3L,
  'pyrrole_nitrogen', 'three pyrrole nitrogens', 3L
)

expected_counts <- tribble(
  ~category, ~expected_count,
  'boron', 1L,
  'outer_terminal_carbon', 6L,
  'imine_nitrogen', 3L,
  'pyrrole_nitrogen', 3L
)

xyz_matrix <- function(tbl)
{
  tbl |>
    select(x, y, z) |>
    as.matrix()
}

best_fit_plane <- function(points)
{
  if (nrow(points) < 3) {
    stop('Need at least three points to fit a plane.', call. = FALSE)
  }

  center <- colMeans(points)
  centered <- sweep(points, 2, center, '-')
  sv <- svd(centered)
  normal <- sv$v[, ncol(sv$v)]
  normal_norm <- sqrt(sum(normal^2))

  if (!is.finite(normal_norm) || normal_norm == 0) {
    stop('Could not determine a plane normal.', call. = FALSE)
  }

  list(
    center = center,
    normal = normal / normal_norm
  )
}

point_plane_distance <- function(point, points)
{
  plane <- best_fit_plane(points)
  as.numeric(abs(sum((point - plane$center) * plane$normal)))
}

normalize_categories <- function(atoms)
{
  if ('category' %in% names(atoms)) {
    atoms |>
      mutate(category = na_if(category, ''))
  } else if ('bsubpc_label' %in% names(atoms)) {
    atoms |>
      mutate(category = str_remove(bsubpc_label, '_[0-9]+$'))
  } else {
    stop(
      'Input table must contain either a category column or a bsubpc_label column.',
      call. = FALSE
    )
  }
}

validate_molecule <- function(mol_atoms, mol_id)
{
  observed_counts <- mol_atoms |>
    filter(!is.na(category)) |>
    count(category, name = 'n')

  count_check <- expected_counts |>
    left_join(observed_counts, by = 'category') |>
    mutate(n = replace_na(n, 0L))

  bad_counts <- count_check |>
    filter(n != expected_count)

  if (nrow(bad_counts) > 0) {
    bad_counts_text <- bad_counts |>
      transmute(text = paste0(category, '=', n, ' (expected ', expected_count, ')')) |>
      pull(text) |>
      paste(collapse = ', ')

    stop(
      sprintf('%s: unexpected category counts {%s}', mol_id, bad_counts_text),
      call. = FALSE
    )
  }
}

bowl_depths_for_molecule <- function(mol_atoms, mol_id)
{
  validate_molecule(mol_atoms, mol_id)

  boron <- mol_atoms |>
    filter(category == 'boron') |>
    xyz_matrix() |>
    as.numeric()

  plane_definitions |>
    mutate(
      mol_id = mol_id,
      bowl_depth = map_dbl(
        plane,
        function(plane_name)
        {
          plane_points <- mol_atoms |>
            filter(category == plane_name) |>
            xyz_matrix()

          point_plane_distance(boron, plane_points)
        }
      )
    ) |>
    select(mol_id, plane, plane_label, bowl_depth)
}

if (!file.exists(input_path)) {
  stop(sprintf('Input file not found: %s', input_path), call. = FALSE)
}

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
    .default = col_guess()
  ),
  show_col_types = FALSE
) |>
  normalize_categories()

required_columns <- c('mol_id', 'x', 'y', 'z', 'category')
missing_columns <- setdiff(required_columns, names(atoms))
if (length(missing_columns) > 0) {
  stop(
    sprintf(
      'Input table is missing required columns: %s',
      paste(missing_columns, collapse = ', ')
    ),
    call. = FALSE
  )
}

if ('atom_id' %in% names(atoms)) {
  atoms <- atoms |>
    arrange(mol_id, atom_id)
} else {
  atoms <- atoms |>
    arrange(mol_id)
}

bowl_depths <- split(atoms, atoms$mol_id) |>
  imap(bowl_depths_for_molecule) |>
  bind_rows() |>
  mutate(
    plane = factor(plane, levels = plane_definitions$plane),
    plane_label = factor(plane_label, levels = plane_definitions$plane_label)
  ) |>
  arrange(mol_id, plane) |>
  mutate(
    plane = as.character(plane),
    plane_label = as.character(plane_label),
    bowl_depth = formatC(bowl_depth, format = 'f', digits = 6)
  )

write_csv(bowl_depths, output_path)
