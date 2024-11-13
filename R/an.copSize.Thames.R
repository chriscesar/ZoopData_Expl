# an.copSize.Thames.R ####
# comparing prevalence of BIG and SMALL copepods in the Thames with other regions
#### load packages ####
ld_pkgs <- c("tidyverse","MASS","lubridate", "tictoc","gllvm","purrr","patchwork")
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

# load data ####
source("R/imp.load_data_all_zoops.R")

## remove non-copep & plot abundance of copepods
dfl %>% 
  filter(.,!is.na(copNonCop)) %>% 
  filter(., copNonCop == "COPEPOD") %>% 
  #remove taxon IDs
  dplyr::select(., -c(Taxa,Aphia.ID, SizeClass:LF0,Kingdom:DisplayName)) %>%
  group_by(across(-c(AbundanceRaw,
                     Abund_m3,
                     mnlongMaxAxis_mm,
                     mdlongMaxAxis_mm,
                     mnCPerIndiv_ug,
                     mdCPerIndiv_ug,
                     mn_carbTot_raw,
                     md_carbTot_raw,
                     mn_carbTot_m3,
                     md_carbTot_m3
                     ))) %>% 
  summarise(.,
            AbundanceRaw = sum(AbundanceRaw, na.rm = TRUE),
            Abund_m3 = sum(Abund_m3, na.rm = TRUE),
            mnlongMaxAxis_mm = sum(mnlongMaxAxis_mm, na.rm = TRUE),
            mdlongMaxAxis_mm = sum(mdlongMaxAxis_mm, na.rm = TRUE),
            mnCPerIndiv_ug = mean(mnCPerIndiv_ug, na.rm = TRUE),
            mdCPerIndiv_ug = mean(mdCPerIndiv_ug, na.rm = TRUE),
            mn_carbTot_raw = sum(mn_carbTot_raw, na.rm = TRUE),
            md_carbTot_raw = sum(md_carbTot_raw, na.rm = TRUE),
            mn_carbTot_m3 = sum(mn_carbTot_m3, na.rm = TRUE),
            md_carbTot_m3 = sum(md_carbTot_m3, na.rm = TRUE),
            .groups = "drop") %>%
  filter(.,LF02 %in% c(
    "Cop_Lg",
    "Cop_Sm"
  )) %>% 
  filter(., Region %in% c(
    "Thames",
    "SWest",
    "NEast",
    "Anglian",
    "NWest",
    "Southern"
    )) %>% 
  ungroup() %>% 
  ggplot(., aes(y = log(mn_carbTot_m3), x = LF02))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(alpha=0.2)+
  facet_wrap(.~Region)


# total carbon content of each sample by Site/Region
