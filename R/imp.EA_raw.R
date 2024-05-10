#####################################################################
#####################################################################
#####################################################################
# TO DO: ####
# Check import and appending/correcting of taxon names: Some taxa reading
# as NA values since update on 05/04/24
#####################################################################
#####################################################################
#####################################################################


### imp.EA_raw.R ####
#### Import EA data, attach WIMS data & create a block for later analysis ####
#### analysed Jul 2023
### RH has amended the files returned from the lab to include
### sample volumes

### useful links for plotting MDS with ggplot
# https://chrischizinski.github.io/rstats/vegan-ggplot2/
# https://www.youtube.com/watch?v=Y0GI34S-ZMI

## set up ####
#### install required packages ####
#### set local package library 

#### check and install required packages ####
# ptm <- Sys.time()
# req_packages <- c("lubridate",# working with dates
#                   "patchwork", #arranging multiple plots
#                   "tidyverse",# general data manipulation
#                   "vegan",# NMDS, multivariate analysis of ecological data
#                   "vegan3d",# as above, but with 3D figs
#                   "mvabund", #multivariate abundance data analyses
#                   "gllvm", # latent variable models
#                   "ecoCopula",#model based ordination/graphical modelling
#                   "ggthemes",# sensible visualisation styles
#                   "openxlsx",# read data from xlsx
#                   "MASS",# fit a negative binomial glm
#                   "gclus",#clustering of data
#                   "corrplot",#correlation plots
#                   "performance",# model checking
#                   "seas"#converting dates to seasonal factors
# )
# 
# new_packages <- req_packages[!(req_packages %in% installed.packages()[,"Package"])]
# if(length(new_packages)) install.packages(new_packages,library=libfolder,type="binary")
# Sys.time() - ptm;rm(ptm,req_packages,new_packages)

#### load packages ####
ld_pkgs <- c("tidyverse","MASS","lubridate", "tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

#### set universals ####
tictoc::tic.clearlog();tic("set universals");print("set universals")
### set up folders & import functions ###
source("R/folder.links.R")

perms <- 9999 ### number of permutations to run for multivariate analyses
theme_set(ggthemes::theme_few())###set theme for all ggplot objects
nit <- 200 #number of iterations
ppi <- 300 #image resolution
# colourblind friendly colour palette (RGB values also commented)
cbPalette <- c("#999999", #153/153/153
               "#E69F00",#230/159/000
               "#56B4E9",#086/180/233
               "#CC79A7", #204/121/167
               "#009E73",#000/158/115
               "#F0E442",#240/228/066
               "#0072B2",#000/114/178
               "#D55E00",#213/094/000
               
               "#444444", 
               "#C34D55",
               "#33A2C4",
               "#554C31",
               "#C5C221",
               "#5531A1",
               "#B32C55",
               "#BB3593" 
               
)

cbPalette2 <- c("#646464", #100/100/100
                "#B46D00",#180/109/0
                "#2482BA",#036/130/186
                "#006C41",#000/108/065
                "#BEB210",#190/178/016
                "#004080",#000/064/128
                "#A32C00",#163/044/000
                "#9A4775"#154/071/117
)
toc(log=TRUE)

### load data ####
tic("Load taxon data and correct taxon names")
print("Load taxon data and correct taxon names")
### taxon data
df_tx0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,
                                               "processedData/MBA_Returns_Amalgamated_USE.xlsx"),
                                        sheet="outR04"))

### append updated taxon names
tx_chk0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,
                                               "processedData/MBA_Returns_Amalgamated_USE.xlsx"),
                                        sheet="TaxonomicRaw"))

tx_chk <- tx_chk0 %>% 
  rename(Taxa=ScientificName_accepted,
         Aphia.ID=AphiaID_accepted) %>% 
  #dplyr::select(.,-ScientificName) %>% 
  distinct()

tx_chk %>% dplyr::select(., Taxa, Aphia.ID) %>% 
  distinct() -> tx_chktrm

df_tx0 <- left_join(df_tx0, tx_chktrm, by="Aphia.ID")
df_tx0$Taxa.x <- df_tx0$Taxa.y;df_tx0$Taxa.y <- NULL
df_tx0 %>%
  rename(Taxa=Taxa.x) %>% 
  dplyr::select(.,-AbundanceRaw) %>% 
  group_by(across(c(!Abund_m3))) %>% 
  summarise(.,Abund_m3=sum(Abund_m3),.groups = "drop") %>%
  ungroup() %>% 
  as_tibble(.) -> df_tx
rm(tx_chk,tx_chk0,tx_chktrm)
toc(log=TRUE)

### WIMS chemical data
tic("Load WIMS data")
print("Load WIMS data")
# WIMS Extract based on:
# Materials = 2HZZ & 2IZZ; SMPT_TYPE = CD, CC, CE; Dates from 01/06/2022-present
# SMP_Code:
# 42100171, 42100174, 42100179, 45400826, 60510027, 73015085, 82510555,
# 82615055, 82615255, 88002837, 88007163, 88007172, 88025879, BE061099,
# E0001449, E0001450, E0004730, G0003532, G0003572, LC544405, LC560357,
# PTTR0026, WA560349, Y0004367, Y0017477, YC536426

df_wims0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,
                                                 "/WIMS_Extract_WaterQuality_Zoop_Samples_240405.xlsx"),
                                          sheet="allDat"))

### counts for pots 1&2 have been summed as to have those for pots 3&4
### removed whitespace from AphiaIDs & taxon names

### prep data ####
### format & widen WIMS data ####
df_wims <- df_wims0

df_wims$PRN <- df_wims$SAMP_SCHEDULE_SAMPLE_ID
df_wims %>%
  dplyr::mutate(det=paste0(DETE_DESC,"_",UNIT_SHORT_DESC)) %>% ##create new variable label
  dplyr::mutate(Result=ifelse(is.na(df_wims$MEAS_SIGN == "<"), MEAS_RESULT,
                              paste0(MEAS_SIGN,MEAS_RESULT))) %>% #add "<" to results
  dplyr::select(.,c(WIMS.Code,REGION,Biosys.ID,
                    SMPT_LONG_NAME, SAMP_SAMPLE_DATE, SAMP_SAMPLE_TIME,
                    SAMP_Notes,
                    PRN,det,Result)) %>% ##only keep variables of interest
  tidyr::pivot_wider(.,names_from=det, values_from = Result) -> df_wims_w###widen data

## write WIMS wide csv
# write.csv(df_wims_w,file=paste0(datfol,"processedData/WIMS_wide.csv"),row.names = FALSE)
toc(log=TRUE)

### remove odd data
#df_tx <- as_tibble(df_tx0) ### create new data (keep df0 as 'raw')
tic("Prep for export & export data");print("Prep for export & export data")

### convert dates
df_tx$sample.date <- as.Date(df_tx$sample.date, origin = "1899-12-30")

# Remove 100 µm data [OPTIONAL] ####
df_tx_100um <- df_tx %>% 
  filter(str_starts(Sample.comments,"100um"))

df_tx %>% 
  filter(!str_starts(Sample.comments,"100um")) -> df_tx

## export taxa data for upload to PLET ####
### this version does not retain 0's unlike the zoopLONG version below
write.csv(df_tx, file = paste0(datfol,
                               "processedData/",
                               "zooplong_no_zero.csv"),
          row.names = FALSE)

###widen data & fill NAs with 0s ####
df_tx %>% 
  dplyr::select(.,-c("Aphia.ID")) %>% 
  pivot_wider(names_from = "Taxa",values_from = "Abund_m3",
              values_fill = 0) -> df_tx_w

### join & save data ####  
dfw <- left_join(df_tx_w,df_wims_w,by="PRN")

write.csv(dfw,file=paste0(datfol,"processedData/","zoopWIDEAbund_m3_WIMS_USE.csv"),row.names = FALSE)
write.csv(df_tx_w,file=paste0(datfol,"processedData/","zoopWIDEAbund_m3_taxOnly_USE.csv"),row.names = FALSE)
saveRDS(dfw, file=paste0(datfol,"processedData/","zoopWIDEAbund_m3_WIMS_USE.RDat"))

### generate LONG version WITH zero values (for calculation of means)
df_tx_w %>% 
  pivot_longer(.,
               cols = Acartia:"Subeucalanus crassus" #check this
               ) -> df_tx_l
write.csv(df_tx_l,file=paste0(datfol,"processedData/","zoopLONG.csv"),row.names = FALSE)
toc(log=TRUE)

### tidy up ###
tic("tidy up")
# unload packages
# detach("package:MASS", unload=TRUE)
detach("package:lubridate", unload=TRUE)
detach("package:tidyverse", unload=TRUE)

# remove data
rm(list = ls(pattern = "^df"))
rm(list = ls(pattern = "^cbPalette"))
rm(datfol,nit,perms, ppi)
toc(log=TRUE);unlist(tictoc::tic.log())