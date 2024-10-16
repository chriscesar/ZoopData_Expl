# an.len.carb.relation.R ####
# investigation into length-carbon relationship for zoops 

# load packages ####
ld_pkgs <- c("tidyverse","MASS","lubridate","vegan","mvabund","seas","patchwork",
             "ecoCopula","performance","gclus","corrplot","gllvm","tictoc","ggpubr")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)
tictoc::tic.clearlog();tic("set universals");print("set universals")

# set universals ####
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

# Load data ####
## load zoops ####
source("R/imp.load_data_all_zoops.R")

#####
# initial plot ####
dfl %>% 
  dplyr::select(Aphia.ID, DisplayName,LF02,
                mnlongMaxAxis_mm,mdlongMaxAxis_mm,
                mnCPerIndiv_ug,mdCPerIndiv_ug) %>% 
  distinct(.) %>% 
  ggplot(., aes(
    x=log(mnlongMaxAxis_mm),
    y=log(mnCPerIndiv_ug))) +
  geom_point()+
  geom_smooth(method = "lm")

# dfl %>% 
#   dplyr::select(Aphia.ID, DisplayName,LF02,
#                 mnlongMaxAxis_mm,mdlongMaxAxis_mm,
#                 mnCPerIndiv_ug,mdCPerIndiv_ug) %>% 
#   distinct(.) %>%# View()
#   filter(.,!is.na(mnlongMaxAxis_mm)) %>% 
#   filter(.,!is.na(mnCPerIndiv_ug)) %>% 
#   ggplot(., aes(
#     x=log(mnlongMaxAxis_mm),
#     y=log(mnCPerIndiv_ug))) +
#   geom_point()+
#   #facet_wrap(.~LF02)+
#   geom_smooth(method = "lm")+
#   stat_regline_equation(label.x = -1, label.y = 7.5)+
#   stat_cor(aes(label=..rr.label..), label.x=-1, label.y=6.5)+
#   labs(title = "Relationship between log(zooplankter length) & log(individual carbon content) in
# zooplankton assemblages recorded in English estuarine & coastal water bodies",
#        x = "log(mean maximum length (mm))",
#        y = "log(mean carbon content per individual (ug)",
#        caption = paste0("Based on length-carbon relationships of ",nrow(.), " taxa"))
# ggsave("figs/2410_zoopSize/len_carb_all.pdf",
#        width = 18,height = 12,units = "in")

dfl %>%
  dplyr::select(Aphia.ID, DisplayName, LF02,
                mnlongMaxAxis_mm, mdlongMaxAxis_mm,
                mnCPerIndiv_ug, mdCPerIndiv_ug) %>%
  distinct() %>%
  filter(!is.na(mnlongMaxAxis_mm)) %>%
  filter(!is.na(mnCPerIndiv_ug)) %>%
  { ggplot(., aes(
    x = log(mnlongMaxAxis_mm),
    y = log(mnCPerIndiv_ug))) +
      geom_point() +
      geom_smooth(method = "lm") +
      facet_wrap(.~LF02)+
      stat_regline_equation(label.x = -1, label.y = 7.5) +
      stat_cor(aes(label = ..rr.label..), label.x = -1, label.y = 6.5) +
      labs(
        title = "Relationship between log(zooplankter length) & log(individual carbon content) in zooplankton assemblages recorded in English estuarine & coastal water bodies",
        x = "log(mean maximum length (mm))",
        y = "log(mean carbon content per individual (ug))",
        caption = paste0("Based on length-carbon relationships of ", nrow(.), " taxa.
                         Values based on those produced by Atkinson et al. (in prep.) as amended by Cesar (unpublished)")
      )
  }
ggsave("figs/2410_zoopSize/len_carb_all_faceted.pdf",
       width = 18,height = 12,units = "in")
