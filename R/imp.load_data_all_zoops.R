# imp.load_data_all_zoops.R ####
### load and combine taxonomic and lifeform data into 1 all-encompassing Zoop data set

#### load packages ####
ld_pkgs <- c("tidyverse","seas","tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

## set universals ####
tic("Set universals");print("set universals")
source("R/set_meta.R")
toc(log=TRUE)

tic("Load lifeforms data, format data and correct names")
print("Load lifeforms data and correct taxon names")

## load data
dfl0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,
                                            "processedData/MBA_Returns_Amalgamated_USE.xlsx"),
                                     sheet="outR04_LF"))

### import  updated taxon names
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

### append updated taxon names to abundance data
dfl0 <- left_join(dfl0, tx_chktrm, by="Aphia.ID")
dfl0$Taxa.x <- dfl0$Taxa.y;dfl0$Taxa.y <- NULL

dfl0 %>%
  rename(Taxa=Taxa.x) %>% # clear up taxon column name
  # remove repeated taxon names within samples: sum by updated taxon names
  group_by(across(-c(AbundanceRaw,Abund_m3))) %>% 
  summarise(.,
            AbundanceRaw=sum(AbundanceRaw, na.rm = TRUE),
            Abund_m3 = sum(Abund_m3, na.rm = TRUE),
            .groups = "drop") %>%
  ungroup() %>% 
  # remove 100um mesh samples [OPTIONAL]
  filter(!str_starts(Sample.comments,"100um")) %>%
  as_tibble(.) -> dfl
rm(tx_chk,tx_chk0,tx_chktrm)
toc(log=TRUE)

# Prep metadata formatting ####
tic("Prep metadata formatting");print("Prep metadata formatting")
### convert sample dates from Excel to Date format
dfl$sample.date <- as.Date(dfl$sample.date, origin = "1899-12-30")
### add day of year variable
dfl$yday <- lubridate::yday(dfl$sample.date)
dfl %>% relocate(.,yday, .after = sample.date) -> dfl
### add month variable
dfl$month <- lubridate::month(dfl$sample.date)
dfl %>% relocate(.,month, .after = sample.date) -> dfl
### add season variable
dfl$DJF <- as.factor(seas::mkseas(dfl$sample.date, width="DJF"))#convert dates to 3month seasonal block
dfl %>% relocate(.,DJF, .after = sample.date) -> dfl

# ### set Region as a factor & order in clockwise from NE to NW
dfl$Region <- factor(dfl$Region, levels = c("NEast","Anglian","Thames",
                                          "Southern","SWest","NWest"))

### generate labels
### Likely wil require updating as samples are gathered at new stations/WBs
LFRegion <- dfl$Region
LFWB <- dfl$WB

WB_lb1 <- ifelse(LFRegion == "Southern","Sth",
                 ifelse(LFRegion == "Thames","Thm",
                        ifelse(LFRegion == "Anglian","Ang",
                               ifelse(LFRegion == "NWest","NW",
                                      ifelse(LFRegion == "NEast","NE",
                                             ifelse(LFRegion == "SWest","SW",NA)
                                      )))))

WB_lb2 <- ifelse(LFWB == "Solent","Solent",
                 ifelse(LFWB == "SOUTHAMPTON WATER","SotonWtr",
                        ifelse(LFWB == "Solway Outer South","SolwOtr",
                               ifelse(LFWB == "THAMES LOWER","ThmLwr",
                                      ifelse(LFWB == "Blackwater Outer","BlckwOtr",
                                             ifelse(LFWB == "Cornwall North","CornwNth",
                                                    ifelse(LFWB == "Barnstaple Bay","BrnstpB",
                                                           ifelse(LFWB == "Kent South","KentSth",
                                                                  ifelse(LFWB == "Mersey Mouth","MerseyMth",
                                                                         ifelse(LFWB == "Wash Outer","WashOtr",
                                                                                ifelse(LFWB == "Lincolnshire","Lincs",
                                                                                       ifelse(LFWB == "Yorkshire South","YorksSth",
                                                                                              ifelse(LFWB == "TEES","Tees",
                                                                                                     ifelse(LFWB == "Northumberland North","NrthmbNth",
                                                                                                            ifelse(LFWB == "Farne Islands to Newton Haven","FarneIs",
                                                                                                                   ifelse(LFWB == "Bristol Channel Inner South","BristInSth",
                                                                                                                          ifelse(LFWB == "Lincs Offshore","LncsOffsh",
                                                                                                                                 ifelse(LFWB == "Isle of Wight East","IoWE",
                                                                                                                                        NA)))))))))))))
                                      )))))
WB_lb <- paste0(WB_lb1,"_",WB_lb2)

dfl$WB_lb <- WB_lb
dfl %>% relocate(WB_lb,.after = WB) -> dfl
dfl$WB_lb <- factor(dfl$WB_lb, levels = c(
  "NE_NrthmbNth",
  "NE_FarneIs",
  "NE_Tees",
  "Ang_YorksSth",
  "Ang_Lincs",
  "Ang_LncsOffsh",
  "Ang_WashOtr",
  "Ang_BlckwOtr",
  "Thm_ThmLwr",
  "Sth_KentSth",
  "Sth_IoWE",
  "Sth_Solent",
  "Sth_SotonWtr",
  "SW_CornwNth",
  "SW_BrnstpB",
  "SW_BristInSth",
  "NW_MerseyMth",
  "NW_SolwOtr"
))
rm(WB_lb,WB_lb1,WB_lb2,LFRegion,LFWB)
toc(log=TRUE)

# load and append carbon content data ####
## Carbon content data
tic("Carbon content data")
## this version has been tweaked.  See 'changeLog' worksheet in .xlsx for info
df_carb <- readxl::read_xlsx(
  paste0(datfol,"Lifeforms/ZOOPLANKTON carbon mass data_v7Oct 2024_AMENDED_USE.xlsx"),
  sheet = "FINAL SHEET TO USE") %>% ## load data
  dplyr::select(.,c(2:6)) %>% ## keep only taxon info, carbon and length info
  janitor::clean_names(.) %>% # tidy names
  rename(copNonCop = copepod_or_non_copepod,
         CPerIndiv_ug = body_mass_as_ug_c_per_individ,
         longMaxAxis_mm = starts_with("longest_max"),
         Aphia.ID = aphia_id) %>% 
  mutate(Aphia.ID=as.numeric(Aphia.ID)) %>% 
  group_by(Aphia.ID) %>% 
  ### calculate mean & median lengths & carbon contents by Aphia IDs
  mutate(mnlongMaxAxis_mm = mean(longMaxAxis_mm, na.rm = TRUE),
         mdlongMaxAxis_mm = median(longMaxAxis_mm, na.rm = TRUE),
         mnCPerIndiv_ug = mean(CPerIndiv_ug, na.rm = TRUE),
         mdCPerIndiv_ug = median(CPerIndiv_ug, na.rm = TRUE)) %>% 
  ungroup(.) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) %>% 
  as_tibble(.) %>% 
  dplyr::select(.,-c(taxon, longMaxAxis_mm,CPerIndiv_ug)) %>%
  distinct(.)

# Append carbon content to taxon data & multiply abundance by Carbon
dfl <- left_join(dfl, df_carb, by="Aphia.ID")
rm(df_carb)

dfl %>%
  mutate(mn_carbTot_raw = AbundanceRaw*mnCPerIndiv_ug,
         md_carbTot_raw = AbundanceRaw*mdCPerIndiv_ug,
         mn_carbTot_m3 = Abund_m3*mnCPerIndiv_ug,
         md_carbTot_m3 = Abund_m3*mdCPerIndiv_ug
  ) -> dfl

toc(log=TRUE)

# Write data [if needed]####
tic("Write data to csv");print("Write data to csv")
# write.csv(dfl, file=paste0(datfol,"processedData/zoopsAll.csv"),
#           row.names = FALSE)
toc(log = TRUE)
unlist(tictoc::tic.log())
