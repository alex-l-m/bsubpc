library(tidyverse)
library(cowplot)
theme_set(theme_cowplot() + theme(plot.background = element_rect(fill = 'white')))

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
)) 

fd_binwdith <- function(x) {
    2 * IQR(x) / length(x)^(1/3)
}
# Facet plot for each element
histograms <- forces |>
    filter(element %in% c('B', 'C', 'N', 'O')) |>
    mutate(`|f|` = sqrt(fx^2 + fy^2 + fz^2)) |>
    ggplot(aes(x = `|f|`)) +
    facet_wrap(~ element, scales = 'free') +
    geom_histogram(binwidth = fd_binwdith) +
    # Units are eV/angstrom
    labs(x = 'Force (eV/angstrom)', y = 'Count', title = 'Distribution of atomic forces')

ggsave('force_histograms.png', histograms, width = 10, height = 6)
