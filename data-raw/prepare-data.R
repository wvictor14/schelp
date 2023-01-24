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


#############################################################################
### smillie 2019 counts and metadata
########## .rds data objects were copied from "output" folder

library(dplyr)
library(here)
library(forcats)
library(readr)
library(purrr)

# clean up metadata
metadata <- read_csv(here::here('data-raw', 'Cells metadata.csv')) %>%

  # recode and relevel health_group
  mutate(Health_group_pretty = case_when(
    Health_group == 'Healthy' ~ 'Healthy',
    Health_group == 'Non-inflamed' ~ 'UC Non-Infl.',
    Health_group == 'Inflamed' ~ 'UC Infl.'
  ) %>% fct(levels = c('Healthy', 'UC Non-Infl.', 'UC Infl.'))) %>%

  # get cell numbers
  group_by(ori.Cluster_group) %>%
  mutate(ori_cluster_id_ncell = paste0(ori.Cluster_group, ' (', n(), ')'))%>%
  ungroup()

# extract gene names
counts <- readRDS(here::here('data-raw', 'Count matrix - log TP10k.rds'))

# filer to quality cells
metadata <- metadata %>%
  filter(qc_status == 'Passed', scrublet_call %in% c('Singlet', 'Not evaluated'))
counts <- counts[,metadata$cellid]

# Save example data
set.seed(1)
metadata <- metadata %>%
  filter(Compartment %in% c('Epithelial', 'Immune')) %>%
  group_by(Compartment, Health_group, ori.Cluster_group) %>%
  dplyr::slice_sample(prop = 0.01) %>% ungroup()

counts <- rbind(
  counts[c('FOXP3', 'NOX1', 'NXPE1', 'MS4A10'), metadata$cellid],
  counts[sample(1:nrow(counts), 26), metadata$cellid]
)

usethis::use_data(metadata, counts, overwrite = TRUE)
