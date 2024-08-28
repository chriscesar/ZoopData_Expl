# imp.load_data_wims.R ####
## load WIMS data

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
                                                 #"/WIMS_Extract_WaterQuality_Zoop_Samples_240618.xlsx"),
                                                 "/WIMS_Extract_240618.xlsx"),
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
