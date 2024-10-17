# an_EA_taxa_gllvm.R ####
# analysis of taxonomic zoop data using gllvm.
# replaces an.EA_zoops_offsets.R

# load packages ####
ld_pkgs <- c("tidyverse")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)


# set metadata & load data ####
source("R/set_meta.R")
source("R/imp.load_data_all_zoops.R") # load biological data
source("R/imp.load_data_wims.R") #load WIMS data

# Lifeforms ####
## Carbon ####
## generate carbon-content data for analysis ####
### widen taxon data & append WIMS data ####
dfl %>% #names()
  dplyr::select(.,c(Pot.Number:WB_lb,LF02,mn_carbTot_raw)) %>% 
  #### assign net volume value of 1 if net is empty
  dplyr::mutate(netVol_use = ifelse(is.na(`Net.volume.sampled.(m3)`),1,`Net.volume.sampled.(m3)`)) %>%
  ### remove 'old' net value
  dplyr::select(.,-c(`Net.volume.sampled.(m3)`,Aphia.ID, Taxa)) %>% 
  ### remove #NA values
  dplyr::filter(., !is.na(mn_carbTot_raw)) %>% 
  #### sum carbon across duplicate lifeforms
  group_by(across(c(!mn_carbTot_raw))) %>% 
  summarise(.,mn_carbTot_raw=sum(mn_carbTot_raw),.groups = "drop") %>%
  ### widen and fill gaps with 0 values
  group_by(across(c(!mn_carbTot_raw))) %>% 
  pivot_wider(.,names_from = LF02, values_from = mn_carbTot_raw,values_fill = 0) %>% 
  ### append WIMS data
  left_join(., df_wims_w, by = "PRN") %>% 
  ungroup() -> tmp_carb_lf
  
### split the data into meta, lifeform, and wims chunks ####
tmp_carb_lf %>%
  ## metadata
  dplyr::select("Pot.Number":"netVol_use") %>% 
  relocate(PRN) -> df_carb_meta

tmp_carb_lf %>% #names()
  dplyr::select(.,netVol_use:WIMS.Code.y) %>%
  dplyr::select(.,-c(netVol_use,WIMS.Code.y)) -> df_carb_lf
df_carb_lf$PRN <- df_carb_meta$PRN
df_carb_lf %>% relocate(PRN) -> df_carb_lf

tmp_carb_lf %>% #names()
  dplyr::select(WIMS.Code.y:ncol(.)) -> df_carb_wims
df_carb_wims$PRN <- df_carb_meta$PRN
df_carb_wims %>% relocate(PRN) -> df_carb_wims

### collapse into single list item ####
df_carb <- list(df_carb_meta=df_carb_meta,
                df_carb_lf=df_carb_lf,
                df_carb_wims=df_carb_wims)

rm(df_carb_meta, df_carb_lf, df_carb_wims, tmp_carb_lf)

### Prep GLLVM ####
### CURRENT CODE INCLUDES 3 seperate list elements.  Column 1 of each list is the PRN.
### this will allow for trimming of data

# tic("Model set up")
# dfw %>% 
#   dplyr::select(-c(1:21)) %>% ###remove metadata info
#   dplyr::select_if(~ !is.numeric(.) || sum(.) !=0) %>% ### drop cols summing to 0
#   dplyr::select_if(~ !is.numeric(.) || sum(. != 0) >= 100) %>%  # Drop numeric columns with <10 non-zero values
#   filter(rowSums(across(where(is.numeric)))!=0) -> dftmp ###remove 'empty' rows
