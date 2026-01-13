# GenerateOSPAR.v2.R ####
# amalgamate data for OSPAR submission

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","sf","units")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

tic("TOTAL")
#source("R/helperFunctions.R")
source("R/set_meta.R")
rm(cbPalette,cbPalette2,ppi,nit,perms)

tic("Load zoop data")
df0 <- readxl::read_xlsx(
  paste0(datfol,
         "processedData/260113_MBA_Returns_Amalgamated_USE.xlsx"),
                         sheet = "outR04_LF")
toc(log=TRUE)

## Convert grid refs ####
tictoc::tic("Convert grid refs")

# Convert BNG to WGS84
dfl_sf <- st_as_sf(df0, coords = c("Eastings","Northings"),crs=27700)
dfl_sf_4326 <- dfl_sf %>% st_transform(crs = 4326)  

LATIT <-  data.frame(st_coordinates(dfl_sf_4326))[,2]
LONGI <- data.frame(st_coordinates(dfl_sf_4326))[,1]
rm(dfl_sf, dfl_sf_4326)
tictoc::toc(log = TRUE)

tictoc::tic("Rename taxon variables")
dfout <- df0

dfout %>% dplyr::rename(
  SPECI = `Aphia ID`,
  # SPECI = DisplayName,
  RLIST = TaxList) -> dfout

tictoc::toc(log=TRUE)

tictoc::tic("append latlong values, convert depths & append AphiaIDs")
## append latlong values
dfout$LATIT <- LATIT
dfout$LONGI <- LONGI
rm(LATIT,LONGI)
tictoc::toc(log = TRUE)

# create sample dates ####
dfout %>% 
  dplyr::mutate(
    SDATE = `sample date` |>
      as.POSIXct(tz = "UTC") |>
      format("%Y%m%d")
  ) %>% 
  dplyr::mutate(
    dep_clean = str_extract_all(`Sample Depth m`, "\\d+\\.?\\d*") %>%  # Extract all numbers
      lapply(as.numeric) %>%                              # Convert to numeric
      sapply(max)                                         # Take max if multiple
  ) %>% 
  dplyr::mutate(dep_clean = case_when(
    is.na(dep_clean) ~ 30,
    TRUE ~ dep_clean
  )) -> dfout

aphia_ids <- dfout %>% dplyr::select(SPECI,Taxa) %>% distinct() %>% 
  rename(aphia_id = SPECI,
         Taxon_recorded = Taxa)
tictoc::toc(log=TRUE)

tictoc::tic("Create data for export")
# Create data for export ####
dfout %>% 
  dplyr::rename(
    PotNum = `Pot Number`,
    SiteNam = `site name`
  ) %>% 
  ## remove a bit of clutter
  dplyr::select(-c(
    Notes_Comments,AFBI,BSH,CEFAS,CNRS,EA,IEO,IFREMER,LLUR,MBA,MSS,NLWKN,
    NOVANA,NRW,NU,PML,RWS,SEPA,SMHI,SAMS,VLIZ,Unallocated,LF0,LF02,Kingdom,
    Phylum,Class,Order,Family,Genus,Subgenus,Species,Subspecies,
  )) %>% 
  ### add/calculate variables
  dplyr::mutate(
    # Reporting laboratory; 32 = EA Head Office
    # RLABO = 32,
    RLABO = "AWUK",
    # Monitoring year
    MYEAR = lubridate::year(`sample date`),
    # Ship or platform code
    # SHIPC = "Coastal Survey Vessel",
    SHIPC = "AA31",  # changed from SHIPC = "AA00",
    # Cruise identifier (series of sampling occasions)
    # If CRUIS is empty set year_month as  default value
    CRUIS = paste0(sprintf("%02d", month(`sample date`)),"_",year(`sample date`)),
    # Station identification /Sampling event ID
    # using BIOSYS Site ID
    RLIST = "ERID", # all taxa have an aphia ID
    SDATE = `sample date` |>
      as.POSIXct(tz = "UTC") |>
      format("%Y%m%d"),
    DTYPE = "ZP", # Zooplankton data
    SMPNO = paste0("Pot_",PotNum), #Need to consider using SampleID in future exports
    # STNNO = `BIOSYS Code`,
    ## concatenate BIOSYS Code & SDATE
    # #############
    # ## FIX1 (FAILED): Append BIOSYS & POTNUMBER #####
    # #############
    # STNNO = paste0(`BIOSYS Code`,PotNum),
    #############
    ## FIX2 (TEST): Append BIOSYS & SAMPLE DATE #####
    #############
    STNNO = paste0(`BIOSYS Code`,"_",SDATE),
    PARAM = "ABUNDNR",
    # MUNIT = "nr",# change to num per m3
    MUNIT = "nr/m3",# change to num per m3 #DONE
    ALABO = "AWUK", #AMENDED #DONE
    SLABO = "AWUK", #AMENDED #DONE
    LATIT = round(LATIT,5),##rounded to 5 dp
    LONGI = round(LONGI,5),##rounded to 5 dp
    SMTYP = "NET", # update with appropriate net type
    PURPM = "B~S~T~E",
    ## Create 'empty' columns
    POSYS = NA_character_,
    # STATN = SiteNam,
    ## remove special characters from site names
    # STATN = str_replace_all(SiteNam,"[^A-Za-z0-9 ]", ""),
    STATN = `BIOSYS Code`,
    WADEP = NA_character_,
    EDATE = NA_character_,
    STIME = NA_character_,
    ATIME = NA_character_,
    ETIME = NA_character_,
    # MNDEP = NA_character_, # change to 0?  #DONE
    # MXDEP = NA_character_, # change to 30? #DONE
    MNDEP = 0, # change to 0? #DONE
    MXDEP = dep_clean, # change missing depths to 30? #DONE
    NOAGG = NA_character_,
    FNFLA = NA_character_,
    FINFL = NA_character_,
    SMVOL = `Net volume sampled (m3)`, #volume sampled
    WIRAN = NA_character_,
    CLMET = "H",
    FLVOL = NA_character_,
    NPORT = NA_character_,
    SFLAG = NA_character_,
    STRID = NA_character_,
    SIZCL = NA_character_,
    SIZRF = NA_character_,
    MAGNI = NA_character_,
    COEFF = NA_character_,
    TRPHY = NA_character_,
    VFLAG = NA_character_,
    QFLAG = NA_character_,
    # VALUE = AbundanceRaw, #change to density per m3
    VALUE = Abund_m3, #change to density per m3 #DONE
    CPORT = NA_character_,
    SDVOL = NA_character_,
    REFSK = NA_character_,
    METST = NA_character_,
    METFP = NA_character_,
    METPT = NA_character_,
    METCX = NA_character_,
    METPS = NA_character_,
    METOA = NA_character_,
    FORML = NA_character_,
    ACCRD = NA_character_,
    ACORG = NA_character_,
    MESHS = case_when(
      str_detect(`Sample comments`, "200um") ~ "200",
      str_detect(`Sample comments`, "100um") ~ "100",
      TRUE ~ "ERROR"
    ),
    SAREA = NA_character_,
    SPEED = NA_character_,
    PDMET = NA_character_,
    SPLIT = NA_character_,
    DURAT = NA_character_,
    ICCOD = NA_character_,
    SUBST = NA_character_,
    DEPOS = NA_character_,
    PCNAP = NA_character_,
    PRSUB = NA_character_,
    MATRX = NA_character_,
    WLTYP = NA_character_,
    MSTAT = NA_character_,
    MPROG = "CEMP~NATL", #MPROG = "NATL"
  ) %>% 
  dplyr::select(
    RLABO,
    MYEAR,
    SHIPC,
    CRUIS,
    STNNO,
    LATIT,
    LONGI,
    # POSYS,
    STATN,
    # WADEP,
    SDATE,
    # EDATE,
    # STIME,
    # ATIME,
    # ETIME,
    DTYPE,
    SMPNO,
    MNDEP,
    MXDEP,
    # NOAGG,
    # FNFLA,
    # FINFL,
    SMVOL,
    # WIRAN,
    CLMET,
    # FLVOL,
    # NPORT,
    SPECI,
    RLIST,
    # SFLAG,
    # STRID,
    # SIZCL,
    # SIZRF,
    # MAGNI,
    # COEFF,
    # TRPHY,
    STAGE,
    PARAM,
    MUNIT,
    # VFLAG,
    # QFLAG,
    VALUE,
    # CPORT,
    # SDVOL,
    ALABO,
    # REFSK,
    # METST,
    # METFP,
    # METPT,
    # METCX,
    # METPS,
    # METOA,
    # FORML,
    # ACCRD,
    # ACORG,
    SLABO,
    SMTYP,
    MESHS,
    # SAREA,
    # SPEED,
    # PDMET,
    # SPLIT,
    # DURAT,
    # ICCOD,
    # SUBST,
    # DEPOS,
    # PCNAP,
    # PRSUB,
    # MATRX,
    # WLTYP,
    # MSTAT,
    PURPM,
    MPROG
  ) -> dfout_exp
tictoc::toc(log = TRUE)

# remove NA sample volume values ####
dfout_exp %>% dplyr::filter(!is.na(SMVOL)) -> dfout_exp

tictoc::tic("Write data, appending today's date to filename")
write.csv(x = dfout_exp,
          file = paste0("outputs/OSPAR/",format(Sys.Date(), format="%Y%m%d"),"EA_OSPAR_Zoop_export_UTF8.csv"),
          row.names = FALSE,
          fileEncoding = "UTF-8" #ensure UTF-8 encoding
          )
tictoc::toc(log = TRUE)
tictoc::toc(log = TRUE)
unlist(tictoc::tic.log())
