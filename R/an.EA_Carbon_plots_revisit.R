# an.EA_Carbon_plots_revisit.R ####
# Re-examiniation of carbon contents lifeforms data


# 00 set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "tictoc","gllvm","patchwork", "lubridate","ggpubr")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

## set universals ####
tictoc::tic.clearlog();tic("Set universals");print("set universals")
source("R/set_meta.R")
toc(log=TRUE)

tic("Load zoop data")
source("R/imp.load_data_all_zoops.R")
toc(log=TRUE)
rm(dfl0)

# create WIDE versions of data based on abundances and carbon contents per m3 ####
## taxonomic ####
### Abundance/m3 ####
dfl %>% 
  dplyr::select(.,c(Pot.Number:"Sample.comments",
                    WB_lb, Aphia.ID,DisplayName, Abund_m3)) %>%
  group_by(across(-c(Abund_m3))) %>% 
  summarise(Abund_m3=sum(Abund_m3),.groups = "drop") %>% 
  dplyr::select(.,-Aphia.ID) %>% 
  pivot_wider(.,names_from = DisplayName, values_from = Abund_m3,
              values_fill = 0) -> dfw_tx_abund_m3

### Carbon (mean) ####
dfl %>% #names()
  dplyr::select(.,c(Pot.Number:"Sample.comments",
                    WB_lb, Aphia.ID,DisplayName, mn_carbTot_m3)) %>%
  group_by(across(-c(mn_carbTot_m3))) %>% 
  summarise(mn_CarbTot_m3=sum(mn_carbTot_m3),.groups = "drop") %>% 
  dplyr::select(.,-Aphia.ID) %>% 
  pivot_wider(.,names_from = DisplayName, values_from = mn_CarbTot_m3,
              values_fill = 0) -> dfw_tx_Cmn_m3

### Carbon (median) ####
dfl %>% #names()
  dplyr::select(.,c(Pot.Number:"Sample.comments",
                    WB_lb, Aphia.ID,DisplayName, md_carbTot_m3)) %>%
  group_by(across(-c(md_carbTot_m3))) %>% 
  summarise(md_carbTot_m3=sum(md_carbTot_m3),.groups = "drop") %>% 
  dplyr::select(.,-Aphia.ID) %>% 
  pivot_wider(.,names_from = DisplayName, values_from = md_carbTot_m3,
              values_fill = 0) -> dfw_tx_Cmd_m3
