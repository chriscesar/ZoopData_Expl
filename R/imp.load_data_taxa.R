# imp.load_data_taxa.R ####
## load taxon data

tic("Load taxon data, format data and correct taxon names")
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
  dplyr::select(.,-Abund_m3) %>% 
  group_by(across(c(!AbundanceRaw))) %>% 
  summarise(.,AbundanceRaw=sum(AbundanceRaw),.groups = "drop") %>%
  ungroup() %>% 
  as_tibble(.) -> df_tx
rm(tx_chk,tx_chk0,tx_chktrm)
toc(log=TRUE)

tic("Prep for export & export data");print("Prep for export & export data")

# convert & format dates ####
df_tx$sample.date <- as.Date(df_tx$sample.date, origin = "1899-12-30")

### add day of year ###
df_tx$yday <- lubridate::yday(df_tx$sample.date)
df_tx %>% relocate(.,yday, .after = sample.date) -> df_tx
### add month ###
df_tx$month <- lubridate::month(df_tx$sample.date)
df_tx %>% relocate(.,month, .after = sample.date) -> df_tx
### add season ###
df_tx$DJF <- as.factor(seas::mkseas(df_tx$sample.date, width="DJF"))#convert dates to 3month seasonal block
df_tx %>% relocate(.,DJF, .after = sample.date) -> df_tx

# Remove 100 µm data [OPTIONAL] ###
df_tx_100um <- df_tx %>% 
  filter(str_starts(Sample.comments,"100um"))

df_tx %>% 
  filter(!str_starts(Sample.comments,"100um")) -> df_tx

# set Region as a factor ####
df_tx$Region <- factor(df_tx$Region, levels = c("NEast","Anglian","Thames",
                                                "Southern","SWest","NWest"))

# create labels ####
LFRegion <- df_tx$Region
LFWB <- df_tx$WB

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
                                                                                                                          ifelse(LFWB == "Isle of Wight East","IoWE",
                                                                                                                                 ifelse(LFWB == "Lincs Offshore","LncsOffsh",
                                                                                                                          NA)))))))))))))
                                      )))))
WB_lb <- paste0(WB_lb1,"_",WB_lb2)
df_tx$WB_lb <- WB_lb
df_tx %>% relocate(WB_lb,.after = WB) -> df_tx
df_tx$WB_lb <- factor(df_tx$WB_lb, levels = c(
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

# widen data & fill NAs with 0s ####
df_tx %>% 
  dplyr::select(.,-c("Aphia.ID")) %>% 
  pivot_wider(names_from = "Taxa",values_from = "AbundanceRaw",#Abund_m3
              values_fill = 0) %>% 
  filter(.,!is.na(`Net.volume.sampled.(m3)`)) %>% ##OPTIONAL: remove 'empty' net volumes
  ungroup() -> df_tx_w
toc(log=TRUE)
