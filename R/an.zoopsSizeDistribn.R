# an.zoopsSizeDistribn.R ####
# Initial analyses of size distributions of zooplankters

# 00 set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "tictoc","gllvm","patchwork", "lubridate","ggpubr",
             "cowplot","ggridges", "seas")
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

## initial look at a single sample
dfl %>% 
  filter(.,Pot.Number==5) -> dftmp

### generate new variable: yyy_mm
## create year_month variable
dfl$yyyy_mm <- paste0(year(dfl$sample.date),"_",sprintf("%02d", month(dfl$sample.date)))

# generate factor levels
fac_tmp22 <- paste0(rep(2022, 12), "_", sprintf("%02d", 1:12))
fac_tmp22 <- fac_tmp22[6:12]
fac_tmp23 <- paste0(rep(2023, 12), "_", sprintf("%02d", 1:12))
fac_tmp24 <- paste0(rep(2024, 12), "_", sprintf("%02d", 1:12))
fac_tmp24 <- fac_tmp24[1:8]

fac_tmp <- c(fac_tmp22,fac_tmp23,fac_tmp24)

# assign dates to factors
dfl$yyyy_mm <- factor(dfl$yyyy_mm, levels = fac_tmp)
rm(fac_tmp,fac_tmp22,fac_tmp23,fac_tmp24)

### make a 'complete' version so that 'missing' surveys will be included in the charts
dfl_complete <- dfl %>%
  group_by(BIOSYS.Code) %>%  # Group by site
  tidyr::complete(yyyy_mm, fill = list(mn_carbTot_m3 = 0))  # Fill missing months with 0 or NA

# dftmp_complete %>%
# # dftmp %>% 
#   ggplot(., aes(
#     x=mnlongMaxAxis_mm,
#     y=yyyy_mm,
#     # y=DJF,
#     weight = mn_carbTot_m3))+
#   geom_density_ridges(alpha=0.7, aes(fill=DJF),
#                       jittered_points=TRUE,
#                       position = position_points_jitter(width = 0.05, height = 0),
#                       point_shape= "|",
#                       point_size = 3,
#                       point_alpha=1
#                       )+
#   ggthemes::theme_few()+
#   xlim(-10,120)+
#   labs(x = "Mean taxon length (mm)",
#        subtitle = paste0(unique(dftmp$Region)[1],": ", unique(dftmp$WB)[1],
#                          " (",unique(dftmp$BIOSYS.Code)[1],")"))+
#   scale_fill_manual(values=c("skyblue","lightgreen","orange","sienna"))+
#   scale_y_discrete(limits=rev)+
#   theme(axis.title.y = element_blank(),
#         legend.position = "none")

# generate pdfs for individual BIOSYS sites ####
## MONTLY ####
# Define the unique levels of BIOSYS.Code (each site)
unique_sites <- unique(dfl_complete$BIOSYS.Code)
# unique_sites <- unique_sites[20:25]

# Loop through each site and create/export a unique plot
purrr::walk(unique_sites, function(site) {
  
  # Filter data for the current site using dplyr::filter()
  dfl_filtered <- dplyr::filter(dfl_complete, BIOSYS.Code == site) #%>% 
    # dplyr::filter(., !is.na(mnlongMaxAxis_mm)) %>% filter(., !is.na(mn_carbTot_m3))
  
  # Create the plot with ggplot2
  p <- dfl_filtered %>%
    ggplot2::ggplot(aes(x = mnlongMaxAxis_mm, y = yyyy_mm, weight = mn_carbTot_m3)) +
    ggridges::geom_density_ridges(alpha=0.7, aes(fill=DJF),
                                  jittered_points=TRUE,
                                  position = position_points_jitter(width = 0.05, height = 0),
                                  point_shape= "|",
                                  point_size = 3,
                                  point_alpha=1
    )+
    ggthemes::theme_few() +  # Theme from ggthemes
    xlim(-10,120)+
    labs(x = "Mean taxon length (mm)",
         subtitle = paste0("Contribution of zooplankton size classes to total carbon content\n",na.omit(dfl_filtered$Region)[1],": ", na.omit(dfl_filtered$WB)[1],
                           " (",site,")"),
         caption = "Monthly sampling program. Fill colours indicate season.")+
    scale_fill_manual(values=c("skyblue","lightgreen","orange","sienna"))+
    scale_y_discrete(limits=rev)+
    theme(axis.title.y = element_blank(),
          legend.position = "none")
  
  # Export the plot to a PDF file using ggplot2::ggsave()
  ggplot2::ggsave(filename = paste0("figs/2410_zoopSize/",
    "zoopSize_monthly_", na.omit(dfl_filtered$WB_lb)[1],"_",site, ".pdf"), plot = p, width = 8, height = 6)
})

## Seasonally ####
# Define the unique levels of BIOSYS.Code (each site)
unique_sites <- unique(dfl_complete$BIOSYS.Code)
#unique_sites <- unique_sites[20:25]

# Loop through each site and create/export a unique plot
purrr::walk(unique_sites, function(site) {
  
  # Filter data for the current site using dplyr::filter()
  dfl_filtered <- dplyr::filter(dfl_complete, BIOSYS.Code == site) %>% 
  dplyr::filter(.,!is.na(DJF))
  # dplyr::filter(., !is.na(mnlongMaxAxis_mm)) %>% filter(., !is.na(mn_carbTot_m3))
  
  # Create the plot with ggplot2
  p <- dfl_filtered %>%
    ggplot2::ggplot(aes(x = mnlongMaxAxis_mm, y = DJF, weight = mn_carbTot_m3)) +
    ggridges::geom_density_ridges(alpha=0.7, aes(fill=DJF),
                                  jittered_points=TRUE,
                                  position = position_points_jitter(width = 0.05, height = 0),
                                  point_shape= "|",
                                  point_size = 3,
                                  point_alpha=1
    )+
    ggthemes::theme_few() +  # Theme from ggthemes
    xlim(-10,120)+
    labs(x = "Mean taxon length (mm)",
         subtitle = paste0("Contribution of zooplankton size classes to total carbon content\n",na.omit(dfl_filtered$Region)[1],": ", na.omit(dfl_filtered$WB)[1],
                           " (",site,")"),
         caption = "Monthly sampling program. Values grouped by season")+
    scale_fill_manual(values=c("skyblue","lightgreen","orange","sienna"))+
    scale_y_discrete(limits=rev)+
    theme(axis.title.y = element_blank(),
          legend.position = "none")
  
  # Export the plot to a PDF file using ggplot2::ggsave()
  ggplot2::ggsave(filename = paste0("figs/2410_zoopSize/",
                                    "zoopSize_seasoanlly_", na.omit(dfl_filtered$WB_lb)[1],"_",site, ".pdf"), plot = p, width = 8, height = 6)
})
