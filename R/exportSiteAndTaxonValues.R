# exportSiteAndTaxonValues.R ####
# creates csv of abundance and carbon across:
# 1) each zoop taxon; and
# 2) summed values across all taxa in each sample

# 00 set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "tictoc")
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

### export taxon-based version
tic("export taxon-based version")
write.csv(dfl, file = "outputs/data_zoops.csv",row.names = FALSE)
toc(log=TRUE)

### generate summaries by sample
tic("export site-based version")
dfl %>% 
  dplyr::select(.,-c(date_site,
                     Analyst.Initial,
                     Sample.Depth.m,
                     "Time.of.sampling.(GMT)?",
                     "Is.a.replicate?",
                     any.other.comments.on.sample.label,
                     Sample.comments,
                     CEA.Notes,
                     Taxa,
                     Aphia.ID,
                     SizeClass:DisplayName,
                     copNonCop)) %>%
  # names() #%>% 
  group_by(across(-c(AbundanceRaw,
                     Abund_m3,
                     mnlongMaxAxis_mm,
                     mdlongMaxAxis_mm,
                     mnCPerIndiv_ug,
                     mdCPerIndiv_ug,
                     mn_carbTot_raw,
                     md_carbTot_raw,
                     mn_carbTot_m3,
                     md_carbTot_m3))) %>% 
  summarise(
    AbundanceRaw = sum(AbundanceRaw,na.rm=TRUE),
    Abund_m3 = sum(Abund_m3,na.rm=TRUE),
    mnlongMaxAxis_mm = mean(mnlongMaxAxis_mm,na.rm=TRUE),
    mdlongMaxAxis_mm = mean(mdlongMaxAxis_mm,na.rm=TRUE),
    mnCPerIndiv_ug = mean(mnCPerIndiv_ug,na.rm=TRUE),
    mdCPerIndiv_ug = mean(mdCPerIndiv_ug,na.rm=TRUE),
    mn_carbTot_raw = sum(mn_carbTot_raw,na.rm=TRUE),
    md_carbTot_raw = sum(md_carbTot_raw,na.rm=TRUE),
    mn_carbTot_m3 = sum(mn_carbTot_m3,na.rm=TRUE),
    md_carbTot_m3 = sum(md_carbTot_m3,na.rm=TRUE),
    .groups = "drop"
  ) %>% 
  ungroup() %>% 
  write.csv(., file="outputs/data_sites.csv",row.names = FALSE)
toc(log=TRUE)

unlist(tictoc::tic.log())
