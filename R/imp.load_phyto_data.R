# imp.load_phyto_data.R ####
# load phytoplankton data


#### load packages ####
ld_pkgs <- c("tidyverse","tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

tic("Load phyto data downloaded from gov.uk Open Data site")
print("Load phyto data downloaded from gov.uk Open Data site")

### load phyto taxonomy info
phy_tx <- read.csv(file = "data/OPEN_DATA_TAXON_INFO.csv")
### load phyto abundance data
#### Data DLed from OPEN DATA
phy_raw <- read.csv(file="data/PHYT_OPEN_DATA_TAXA.csv") #608279, 54
phy_raw$SAMPLE_DATE <- as.Date(phy_raw$SAMPLE_DATE, format = "%d/%m/%Y")

left_join(phy_raw,phy_tx, by="TAXON_LIST_ITEM_KEY") %>% as_tibble(.) %>% 
  dplyr::select(.,c(WIMS_SITE_ID,SAMPLE_ID,SAMPLE_DATE,CELLS_LITRE,
                    COLS_LITRE,PREFERRED_TAXON_NAME)) %>% 
  rename(TAXON_NAME=PREFERRED_TAXON_NAME) %>% 
  filter(.,!is.na(WIMS_SITE_ID)) %>% 
  mutate(Abund_litre = ifelse(is.na(.$COLS_LITRE),.$CELLS_LITRE,.$COLS_LITRE)) %>%
  dplyr::select(.,-c(CELLS_LITRE,COLS_LITRE)) %>% 
  mutate(codeDate = paste0(WIMS_SITE_ID,SAMPLE_DATE)) %>%
  group_by(across(-c(Abund_litre))) %>% 
  summarise(Abund_litre=sum(Abund_litre),.groups = "drop") %>% 
  {list(dfphy_l = .,
        dfphy_w = pivot_wider(.,names_from = TAXON_NAME, values_from = Abund_litre,
                              values_fill = 0)
  )
    } -> phy_results_tmp

dfphy_l <- phy_results_tmp$dfphy_l
dfphy_w <- phy_results_tmp$dfphy_w
rm(phy_results_tmp)

dfphy_w %>% 
  filter(., SAMPLE_DATE > "2022-05-30") %>% 
  # remove numeric variables that are empty
  dplyr::select(where(~ !is.numeric(.x)||sum(.x) !=0))-> dfphy_w_recent
toc(log=TRUE)

### load wims data and append PRN to phyto data
## WIMS data
tic("load data sets: WIMS data")
source("R/imp.load_data_wims.R")
rm(df_wims,df_wims0,df_wims_l)
toc(log=TRUE)

#### append PRN value to phyto data ###
df_wims_w %>% mutate(codeDate = paste0(WIMS.Code,SAMP_SAMPLE_DATE)) %>% 
  dplyr::select(.,c(codeDate,PRN))-> df_wims_PRN
rm(df_wims_w)

left_join(dfphy_w_recent,df_wims_PRN, by="codeDate") %>% 
  filter(.,!is.na(PRN)) %>% 
  dplyr::select(.,-codeDate) -> dfphy_w_recent_prn
  # write.csv(., file="tmp.csv",row.names = FALSE)
toc(log = TRUE)
