# imp.load_data_lifeforms.R ####
## load lifeforms data

#### load packages ####
ld_pkgs <- c("tidyverse","seas","tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)


tic("Load lifeforms data, format data and correct names")
print("Load lifeforms data and correct taxon names")

### lifeforms data
### taxon data
df_tx0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,
                                               "processedData/MBA_Returns_Amalgamated_USE.xlsx"),
                                        sheet="outR04_LF"))
### prep taxon data ####
### remove odd data
df_tx <- as_tibble(df_tx0) ### create new data (keep df0 as 'raw')

### convert dates
df_tx$sample.date <- as.Date(df_tx$sample.date, origin = "1899-12-30")

# Remove 100 µm data [OPTIONAL] ####
df_tx_100um <- df_tx %>% 
  filter(str_starts(Sample.comments,"100um"))

df_tx %>% 
  filter(!str_starts(Sample.comments,"100um")) -> df_tx

toc(log=TRUE)

###############
tic("Prep for export & export data");print("Prep for export & export data")

### convert dates
df_tx$sample.date <- as.Date(df_tx$sample.date, origin = "1899-12-30")

### add month
df_tx$month <- lubridate::month(df_tx$sample.date)
df_tx %>% relocate(.,month, .after = sample.date) -> df_tx
### add season
df_tx$DJF <- as.factor(seas::mkseas(df_tx$sample.date, width="DJF"))#convert dates to 3month seasonal block
df_tx %>% relocate(.,DJF, .after = sample.date) -> df_tx

# Remove 100 µm data [OPTIONAL] ####
df_tx %>% 
  filter(!str_starts(Sample.comments,"100um")) -> df_tx

### set Region as a factor
df_tx$Region <- factor(df_tx$Region, levels = c("NEast","Anglian","Thames",
                                                "Southern","SWest","NWest"))

toc(log=TRUE)

###widen data & fill NAs with 0s ####
tic("widen data")
df_tx %>% 
  dplyr::select(.,-c("Aphia.ID","Abund_m3","Taxa","Category":"Unallocated",
                     "LF0","Kingdom":"Subspecies","DisplayName")) %>% #drop unneeded cols
  group_by(across(c(-AbundanceRaw))) %>% # group by everything except abundance
  summarise(AbundanceRaw=sum(AbundanceRaw), #sum abundances
            .groups="drop") %>% 
  pivot_wider(names_from = "LF02",values_from = "AbundanceRaw", #widen
              values_fill = 0) %>% 
  ### remove samples missing net volume data
  dplyr::filter(., !is.na(`Net.volume.sampled.(m3)`)) -> df_tx_w
toc(log=TRUE)