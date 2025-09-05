# an_EA_taxa_gllvm_no_offset.R ####
# analysis of taxonomic zoop data using gllvm.
# replaces an.EA_zoops_offsets.R

# This version uses the carbon per m3 data, removing the computationally-hungry
# need for an offset term

# load packages ####
ld_pkgs <- c("tidyverse", "tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)


tic("set metadata & load data")
# set metadata & load data ####
source("R/set_meta.R")
source("R/imp.load_data_all_zoops.R") # load biological data
source("R/imp.load_data_wims.R") #load WIMS data
toc(log=TRUE)

tic("Format & prep data")
# Lifeforms ####
## Carbon ####
## generate carbon-content data for analysis ####
### widen taxon data & append WIMS data ####
dfl %>% #names()
  dplyr::select(.,c(Pot.Number:WB_lb,LF02,mn_carbTot_raw)) %>% #names()
  #### assign net volume value of 1 if net is empty
  dplyr::mutate(netVol_use = ifelse(is.na(`Net.volume.sampled.(m3)`),1,`Net.volume.sampled.(m3)`)) %>%
  ### remove 'old' net value
  dplyr::select(.,-c(`Net.volume.sampled.(m3)`,Aphia.ID, Taxa)) %>% #names()
  ### remove #NA values
  dplyr::filter(., !is.na(mgC_per.m3_mn)) %>% 
  
  ## drop unnecessary columns which interfere with pivot_wider
  dplyr::select(-c(
    max_axis_length_mm, ugC_per_individ_mn, ugC_in_sample_mn,
    ugC_per.m3_mn, mgC_per_individ_mn, mgC_in_sample_mn,
    ugC_per_individ_md, ugC_in_sample_md, ugC_per.m3_md,
    mgC_per_individ_md, mgC_in_sample_md, mgC_per.m3_md,
    mn_carbTot_raw
    )) %>% 
  
  #### sum carbon across duplicate lifeforms
  group_by(across(c(!mgC_per.m3_mn))) %>% 
  summarise(.,mgC_per.m3_mn=sum(mgC_per.m3_mn),.groups = "drop") %>%
  ### widen and fill gaps with 0 values
  group_by(across(c(!mgC_per.m3_mn))) %>% 
  pivot_wider(.,names_from = LF02, values_from = mgC_per.m3_mn,values_fill = 0) %>% 
  ### append WIMS data
  left_join(., df_wims_w, by = "PRN") %>% 
  ungroup() -> tmp_carb_lf
  
### split the data into meta, lifeform, and wims chunks ####
tmp_carb_lf %>%
  ## metadata
  dplyr::select(
   Pot.Number,                                                                 
   date_site,                                                                  
   sample.date,                                                                
   DJF,                                                                        
   month,                                                                      
   yday,                                                                       
   site.name,                                                                  
   BIOSYS.Code,                                                                
   WIMS.Code.x,                                                                
   Analyst.Initial,                                                            
   Sample.Depth.m,                                                             
   `Time.of.sampling.(GMT)?`,
   `Is.a.replicate?`,
   any.other.comments.on.sample.label,                                         
   PRN,                                                                        
   Sample.comments,                                                            
   CEA.Notes,                                                                  
   Eastings,                                                                   
   Northings,                                                                  
   Region,                                                                     
   WBID,                                                                       
   WB,                                                                         
   WB_lb,                                                                      
   netVol_use
   ) %>% 
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
toc(log=TRUE)
