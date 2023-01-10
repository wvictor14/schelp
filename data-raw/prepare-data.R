## code to prepare `genelists` dataset goes here
library(readr)
library(dplyr)
library(here)
library(tidyr)
library(tibble)
library(purrr)

treg <- read_csv(here::here('data-raw', 'Ferraro2014_genelist.csv'))
genelists <- list(

  # Compartment genes From Smillie 2019
  'smillie2019_epithelial' = c('EPCAM', 'KRT8','KRT18'),

  'smillie2019_stromal' =
    c('COL1A1', 'COL1A2', 'COL6A1', 'COL6A2','VWF', 'PLVAP', 'CDH5', 'S100B'),

  'smillie2019_immune' =
    c('CD52', 'CD2', 'CD3D', 'CD3G', 'CD3E', 'CD79A', 'CD79B', 'CD14',
      'CD16','FCGR3A', # note that FCGR3A is an alias for CD16
      'CD68','CD83','CSF1R','FCER1G'),

  'ferraro2014_treg_tconv_up' = treg$`Treg up`,
  'ferraro2014_treg_tconv_down' = treg$`Treg down`[complete.cases(treg$`Treg down`)],

  # From Kanke 2022 and Smillie 2019
  'colon_stem' <- c("LGR5"),
  'colon_goblet' <- c("MUC2"),
  'colon_ta_cycling' <- c('CENPA'), #referred to as g2-m-g1
  'colon_ta_s' <- c('PCNA'),
  'colon_enteroendocrine' <- c("CHGA"),
  'colon_colonocytes' <- c('CA1', 'CAECAM7'),
  'colon_enterocytes_BEST4' <- c('BEST4', 'OTOP2', 'SPIB', 'CA7'),

  'cryptaxis_top' <- c('SELENOP'),

  # Immune cells kanke 2022 + Smillie + others
  'immune_tcells_tregs' <- c('FOXP3', 'HELIOS', 'IL2RA'),
  'immune_bcells' <- c('CD27', 'CD79A', 'IGHM'),
  'immune_bcells_memory' <- c('MS4A1', 'CD19'),
  'immune_bcells_plasma' <- c('SLAMF7', 'IGHG1'),

)

map(genelists, length)

usethis::use_data(genelists, overwrite = TRUE)
