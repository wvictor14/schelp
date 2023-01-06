## code to prepare `genelists` dataset goes here
library(readr)
library(dplyr)
library(here)
library(tidyr)
library(tibble)
library(purrr)

treg <- read_csv(here::here('data-raw', 'Ferraro2014_genelist.csv'))
genelists <- list(

  # From Smillie 2019
  'smillie2019_epithelial' = c('EPCAM', 'KRT8','KRT18'),

  'smillie2019_stromal' =
    c('COL1A1', 'COL1A2', 'COL6A1', 'COL6A2','VWF', 'PLVAP', 'CDH5', 'S100B'),

  'smillie2019_immune' =
    c('CD52', 'CD2', 'CD3D', 'CD3G', 'CD3E', 'CD79A', 'CD79B', 'CD14',
      'CD16','FCGR3A', # note that FCGR3A is an alias for CD16
      'CD68','CD83','CSF1R','FCER1G'),

  'ferraro2014_treg_tconv_up' = treg$`Treg up`,
  'ferraro2014_treg_tconv_down' = na.omit(treg$`Treg down`)
)

map(genelists, length)

usethis::use_data(genelists, overwrite = TRUE)
