# an.EA_Carbon_plots_for_Stats_Doc.R ####
# Re-examiniation of carbon contents lifeforms data

# 00 set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "tictoc","patchwork", "lubridate","ggpubr",
             "ggridges")
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

# Which WB do we want?
wb_label <- "NE_FarneIs"

# Carbon by zoop type by WB ####
## Median ####
dfl %>% #names(.)
  mutate(.,LF03 = if_else(Copepod == "Y", paste0("Cop_",CopSize),ZooType)) %>%
  mutate(.,LF03 = if_else(is.na(LF03) | LF03 == "NYA", LF02, LF03)) %>%
  dplyr::select(.,c(Pot.Number, sample.date, yday, BIOSYS.Code, WIMS.Code,
                    LF03, PRN, 
                    md_carbTot_m3, Region,WB_lb)) %>% 
  mutate(md_carbTot_m3 = replace_na(md_carbTot_m3, 0)) %>% 
  filter(.,md_carbTot_m3 != 0) %>% 
  group_by(across(c(!md_carbTot_m3))) %>% 
  summarise(md_carbTot_m3 = sum(md_carbTot_m3),.groups = "drop") %>%
  ungroup(.) %>%
  dplyr::filter(.,WB_lb == wb_label) %>% 
  ggplot(., aes(x=LF03, y=log(md_carbTot_m3), col=Region))+
  geom_boxplot(varwidth = FALSE,outlier.shape = NA)+
  coord_flip()+
  facet_wrap(.~WB_lb)+
  # facet_wrap(.~WB_lb)+
  geom_point(aes(group=LF03), position=position_jitter(width = 0.2),alpha=0.5,size=3)+
  scale_colour_manual(values = cbPalette2)+
  scale_x_discrete(limits=rev)+
  labs(
    title = "Total carbon content by zooplankton type\nNorth East: Farne Islands",
    y="log(Total carbon)",
    caption = paste0("Samples gathered between ",format(min(dfl$sample.date), "%d/%m/%Y")," & ",format(max(dfl$sample.date), "%d/%m/%Y"),
                     "\nCarbon content values are based on *median* estimates of carbon contents per taxon"))+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_text(face=2),
        axis.text.x = element_text(face=2),
        axis.text.y = element_text(face=2,
                                   size = 12)) -> pl_zoopC_Box

# Carbon-weighted ####
### generate new variable: yyy_mm
## create year_month variable
dfl$yyyy_mm <- paste0(year(dfl$sample.date),"_",sprintf("%02d", month(dfl$sample.date)))

# generate factor levels
fac_tmp22 <- paste0(rep(2022, 12), "_", sprintf("%02d", 1:12))
fac_tmp22 <- fac_tmp22[6:12]
fac_tmp23 <- paste0(rep(2023, 12), "_", sprintf("%02d", 1:12))
fac_tmp24 <- paste0(rep(2024, 12), "_", sprintf("%02d", 1:12))
fac_tmp24 <- fac_tmp24[1:10]

fac_tmp <- c(fac_tmp22,fac_tmp23,fac_tmp24)

# assign dates to factors
dfl$yyyy_mm <- factor(dfl$yyyy_mm, levels = fac_tmp)
rm(fac_tmp,fac_tmp22,fac_tmp23,fac_tmp24)

### make a 'complete' version so that 'missing' surveys will be included in the charts
dfl_complete <- dfl %>%
  group_by(BIOSYS.Code) %>%  # Group by site
  # tidyr::complete(yyyy_mm, fill = list(mn_carbTot_m3 = 0)) %>% # Fill missing months with 0 or NA
  # tidyr::complete(DJF, fill = list(mn_carbTot_m3 = 0)) #Fill missing seasons
  tidyr::complete(yyyy_mm, DJF, fill = list(mn_carbTot_m3 = 0))

# Define the unique levels of BIOSYS.Code (each site)
tic("Carbon-weighted plots of monthly data by BIOSYS site")

## create ridge plot ####
dfl_complete %>% 
  dplyr::filter(.,WB_lb == wb_label) %>%droplevels(.) %>% 
  ggplot2::ggplot(aes(x = mnlongMaxAxis_mm, y = yyyy_mm, weight = mn_carbTot_m3)) +
  ggridges::geom_density_ridges(alpha=0.7, aes(fill=DJF),
                                jittered_points=TRUE,
                                position = ggridges::position_points_jitter(width = 0.05, height = 0),
                                point_shape= "|",
                                point_size = 3,
                                point_alpha=1
                                )+
    ggthemes::theme_few() +  # Theme from ggthemes
    xlim(-10,120)+
    labs(x = "Mean taxon length (mm)",
         title = paste0("Contribution of zooplankton size classes to total carbon content\n","North East: Farne Islands"),
         caption = "Monthly sampling program. Fill colours indicate season.")+
    scale_fill_manual(values=c("skyblue","lightgreen","orange","sienna"),
                      limits=c("DJF", "MAM", "JJA", "SON")) +  # Ensure all seasons are present
    scale_y_discrete(limits=rev)+
    theme(axis.title.y = element_blank(),
          legend.position = "none",
          axis.text.y = element_text(face=2,
                                     size = 10)) -> pl_zoopC_size_ridge
toc(log = TRUE)

pl_both <- pl_zoopC_Box+pl_zoopC_size_ridge
pl_both + patchwork::plot_annotation(tag_levels = "A")

png(file = "figs/FarneIslExample.png",
    width=18*ppi, height=12*ppi, res=ppi)
pl_both + patchwork::plot_annotation(tag_levels = "A")
dev.off()

unlist(tictoc::tic.log())

# Tidy up ####
rm(list=ls(pattern = "^df"))
rm(list=ls(pattern = "^pl"))
rm(list=ls(pattern = "^cb"))
rm(datfol,nit,perms,ppi,wb_label)

detach("package:patchwork", unload=TRUE)
detach("package:lubridate", unload=TRUE)
detach("package:ggpubr", unload=TRUE)
detach("package:ggridges", unload=TRUE)
detach("package:tictoc", unload=TRUE)
detach("package:tidyverse", unload=TRUE)
