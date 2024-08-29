# an.EA_TSPlots_for_report.R ####
# Comparison of Joint Species Distribution Model approaches for zooplankton
# lifeforms data

# 00 set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "tictoc","gllvm","patchwork", "lubridate","ggpubr")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

## set universals ####
tictoc::tic.clearlog();tic("Set universals");print("set universals")
source("R/set_meta.R")
toc(log=TRUE)

### create subdirectory for figures
tic("Create subdirectory for figures")
check_and_create_dir <- function(parent_dir, sub_dir) {
  # Construct the full path to the sub-directory
  full_path <- file.path(parent_dir, sub_dir)
  
  # Check if the sub-directory exists
  if (dir.exists(full_path)) {
    return("Directory already exists")
  } else {
    # Create the sub-directory
    dir.create(full_path)
    return(paste("Directory", sub_dir, "created in", parent_dir))
  }
}

parent_directory <- "figs"
sub_directory <- "2407dd_timeseries"
message <- check_and_create_dir(parent_directory, sub_directory)
print(message)
rm(parent_directory,sub_directory,message,check_and_create_dir)
toc(log=TRUE)

## load data ####

## lifeforms data
tic("load data sets: lifeforms data")
source("R/imp.load_data_lifeforms.R")
df_lf_l <- df_tx
df_lf_w <- df_tx_w
rm(df_tx_w,df_tx,df_tx_100um,df_tx0)
toc(log = TRUE)

## taxon data
tic("load data sets: taxon data")
source("R/imp.load_data_taxa.R")
df_tx_l <- df_tx
rm(df_tx,df_tx_100um,df_tx0)
toc(log=TRUE)

## WIMS data
tic("load data sets: WIMS data")
source("R/imp.load_data_wims.R")
df_wims_l <- df_wims
rm(df_wims,df_wims0)
toc(log=TRUE)

## join data ####  
tic("Join taxon & WIMS data")
dfw_tx <- left_join(df_tx_w,df_wims_w,by="PRN")
toc(log=TRUE)

tic("Join lifeforms & WIMS data")
dfw_lf <- left_join(df_lf_w,df_wims_w,by="PRN")
toc(log=TRUE)
rm(df_wims_w,df_lf_w,df_tx_w)

## Carbon content data
tic("Carbon content data")
df_carb <- readxl::read_xlsx(paste0(datfol,"Lifeforms/ZOOPLANKTON carbon mass data_v5Aug 2024.xlsx"),
                  sheet = "FINAL SHEET TO USE") %>% 
  dplyr::select(.,c(2:6)) %>% 
  janitor::clean_names(.) %>% 
  rename(copNonCop = copepod_or_non_copepod,
         CPerIndiv_ug = body_mass_as_ug_c_per_individ,
         longMaxAxis_mm = starts_with("longest_max"),
         Aphia.ID = aphia_id) %>% 
  mutate(Aphia.ID=as.numeric(Aphia.ID)) %>% 
  group_by(Aphia.ID) %>% 
  mutate(mnlongMaxAxis_mm = mean(longMaxAxis_mm, na.rm = TRUE),
         mdlongMaxAxis_mm = median(longMaxAxis_mm, na.rm = TRUE),
         mnCPerIndiv_ug = mean(CPerIndiv_ug, na.rm = TRUE),
         mdCPerIndiv_ug = median(CPerIndiv_ug, na.rm = TRUE)) %>% 
  ungroup(.) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) %>% 
  as_tibble(.)

df_carb %>% 
  dplyr::select(.,-c(taxon, longMaxAxis_mm,CPerIndiv_ug)) %>%
  distinct(.) -> df_carb_summary
toc(log=TRUE)

### create Region_WB variable
# LFRegion <- dfw_lf$Region
# LFWB <- dfw_lf$WB
# 
# WB_lb1 <- ifelse(LFRegion == "Southern","Sth",
#                     ifelse(LFRegion == "Thames","Thm",
#                            ifelse(LFRegion == "Anglian","Ang",
#                                   ifelse(LFRegion == "NWest","NW",
#                                          ifelse(LFRegion == "NEast","NE",
#                                                 ifelse(LFRegion == "SWest","SW",NA)
#                                          )))))
# 
# WB_lb2 <- ifelse(LFWB == "Solent","Solent",
#                     ifelse(LFWB == "SOUTHAMPTON WATER","Soton Wtr",
#                            ifelse(LFWB == "Solway Outer South","Solway O",
#                                   ifelse(LFWB == "THAMES LOWER","Thm Low",
#                                          ifelse(LFWB == "Blackwater Outer","Blckw Out",
#                                                 ifelse(LFWB == "Cornwall North","Cornw Nth",
#                                                        ifelse(LFWB == "Barnstaple Bay","Brnstp B",
#                                                               ifelse(LFWB == "Kent South","Kent Sth",
#                                                                      ifelse(LFWB == "Mersey Mouth","Mersey Mth",
#                                                                             ifelse(LFWB == "Wash Outer","Wash Out",
#                                                                                    ifelse(LFWB == "Lincolnshire","Lincs",
#                                                                                           ifelse(LFWB == "Yorkshire South","Yorks Sth",
#                                                                                                  ifelse(LFWB == "TEES","Tees",
#                                                                                                         ifelse(LFWB == "Northumberland North","Nrthmb Nth",
#                                                                                                                ifelse(LFWB == "Farne Islands to Newton Haven","Farne Is",
#                                                                                                                       ifelse(LFWB == "Bristol Channel Inner South","Brist Ch In Sth",
#                                                                                                                              NA)))))))))))
#                                          )))))
# WB_lb <- paste0(WB_lb1,"_",WB_lb2)
# #rm(WB_lb1,WB_lb2)
# 
# dfw_lf$WB_lb <- WB_lb
# dfw_lf %>% relocate(WB_lb,.after = WB) -> dfw_lf
# dfw_lf$WB_lb <- factor(dfw_lf$WB_lb, levels = c(
#   "NE_Nrthmb Nth",
#   "NE_Farne Is",
#   "NE_Tees",
#   "Ang_Yorks Sth",
#   "Ang_Lincs",
#   "Ang_Wash Out",
#   "Ang_Blckw Out",
#   "Thm_Thm Low",
#   "Sth_Kent Sth",
#   "Sth_Solent",
#   "Sth_Soton Wtr",  
#   "SW_Cornw Nth",
#   "SW_Brnstp B",
#   "SW_Brist Ch In Sth",
#   "NW_Mersey Mth", 
#   "NW_Solway O"
#   ))
# 
# dfw_tx$WB_lb <- WB_lb
# dfw_tx %>% relocate(WB_lb,.after = WB) -> dfw_tx
# dfw_tx$WB_lb <- factor(dfw_tx$WB_lb, levels = c(
#   "NE_Nrthmb Nth",
#   "NE_Farne Is",
#   "NE_Tees",
#   "Ang_Yorks Sth",
#   "Ang_Lincs",
#   "Ang_Wash Out",
#   "Ang_Blckw Out",
#   "Thm_Thm Low",
#   "Sth_Kent Sth",
#   "Sth_Solent",
#   "Sth_Soton Wtr",  
#   "SW_Cornw Nth",
#   "SW_Brnstp B",
#   "SW_Brist Ch In Sth",
#   "NW_Mersey Mth", 
#   "NW_Solway O"
# ))
########

#####

# PLOTS! ####
## Lifeforms ####

## Copepod size classes ###
# by wb
dfw_lf %>%
  dplyr::select(., c(Region,WB,WB_lb,sample.date,yday,month,`Net.volume.sampled.(m3)`,Cop_Sm, Cop_Ambi, Cop_Lg,Cop_NYA)) %>%
  mutate(Cop_Sm=(`Net.volume.sampled.(m3)`*Cop_Sm)/max(`Net.volume.sampled.(m3)`*Cop_Sm),
         Cop_Ambi=(`Net.volume.sampled.(m3)`*Cop_Ambi)/max(`Net.volume.sampled.(m3)`*Cop_Ambi),
         Cop_Lg=(`Net.volume.sampled.(m3)`*Cop_Lg)/max(`Net.volume.sampled.(m3)`*Cop_Lg),
         Cop_NYA=(`Net.volume.sampled.(m3)`*Cop_NYA)/max(`Net.volume.sampled.(m3)`*Cop_NYA)) %>% 
  pivot_longer(.,cols = Cop_Sm:Cop_NYA,names_to = "lf",values_to = "abund") %>%
  ggplot(.,aes(
    # x=sample.date,
    x=yday,
    y=log(abund+1),
    colour=lf)) +
  geom_line()+
  facet_wrap(.~WB_lb, scales = "free_y")+
  labs(x="Day of year",
       y="Relative abundance",
       title="Relative abundance of copepod size classes recorded in EA water bodies",
       caption="Relative abundance of a given size class is the value recorded in a sample divided by the maximum recorded value for that size class across all samples")+
  theme(legend.title = element_blank(),
        strip.text = element_text(face=2),
        axis.title = element_text(face=2)) -> pl

ggsave(plot = pl, filename = "figs/2407dd_timeseries/copSizesTSByWB.pdf",
       width = 12,height = 8,units = "in")
rm(pl)

# by Region
dfw_lf %>%
  dplyr::select(., c(Region,WB,WB_lb,sample.date,month,DJF,yday,`Net.volume.sampled.(m3)`,Cop_Sm, Cop_Ambi, Cop_Lg,Cop_NYA)) %>%
  mutate(Cop_Sm=(`Net.volume.sampled.(m3)`*Cop_Sm)/max(`Net.volume.sampled.(m3)`*Cop_Sm),
         Cop_Ambi=(`Net.volume.sampled.(m3)`*Cop_Ambi)/max(`Net.volume.sampled.(m3)`*Cop_Ambi),
         Cop_Lg=(`Net.volume.sampled.(m3)`*Cop_Lg)/max(`Net.volume.sampled.(m3)`*Cop_Lg),
         Cop_NYA=(`Net.volume.sampled.(m3)`*Cop_NYA)/max(`Net.volume.sampled.(m3)`*Cop_NYA)) %>% 
  pivot_longer(.,cols = Cop_Sm:Cop_NYA,names_to = "lf",values_to = "abund") %>%
  ggplot(.,aes(
    x=sample.date,
    # x=yday,
    y=abund,
    colour=lf
    )) +
  # geom_line()+
  geom_smooth(se=FALSE,method="loess",span=0.9)+
  # geom_point()+
  facet_wrap(.~Region, scales = "free_y")+
  labs(x=NULL,
       y="Relative abundance",
       title="Modelled relative abundances of copepod size classes recorded in EA water bodies",
       subtitle="Curves represent loess smooths with span = 0.9",
       caption="Relative abundance of a given size class is the value recorded in a sample
       divided by the maximum recorded value for that size class across all samples")+
  theme(legend.title = element_blank(),
        axis.title.y= element_text(face=2),
        strip.text = element_text(face=2)) -> pl

ggsave(plot = pl, filename = "figs/2407dd_timeseries/copSizeSmoothTSByRgn.pdf",
       width = 12,height = 8,units = "in")
rm(pl)

## small copepods
# by WB
dfw_lf %>% 
  ggplot(.,aes(
    x=sample.date,
    # x=yday,
    y= log(Cop_Sm*`Net.volume.sampled.(m3)`+1)
    ))+
  geom_line(colour=1)+
  scale_fill_manual(values=c("slategray1","green","sienna2","brown"))+
  scale_shape_manual(values=c(21:24))+
  geom_point(aes(shape=DJF, fill=DJF), size=2)+
  facet_wrap(.~WB_lb, scales = "free_y")+
  labs(title = "Temporal trends in small copepod abundances in EA water bodies",
       y = "log(Abundance, n+1)",#bquote("Abudance (m"^3~")"),
       x=NULL,
       caption = "Colour and shape of points indicates sampling season:
       Winter (grey circle), Spring (green square), Summer (Orange diamond), Autumn (brown triangle)") +
  theme(legend.position = "none",
        axis.title.y = element_text(face=2),
        strip.text = element_text(face=2)) -> pl

ggsave(plot = pl, filename = "figs/2407dd_timeseries/smCopTSByWBLogN.pdf",
       width = 12,height = 8,units = "in")
rm(pl)

# by Region
dfw_lf %>% 
  ggplot(.,aes(
    x=sample.date,
    #x=yday,
    y= Cop_Sm*`Net.volume.sampled.(m3)`
  ))+
  # geom_line(colour=1)+
  scale_fill_manual(values=c("slategray1","green","sienna2","brown"))+
  scale_shape_manual(values=c(21:24))+
  geom_smooth(se=FALSE, span=0.9)+
  geom_point(aes(shape=DJF, fill=DJF), size=2)+
  facet_wrap(.~Region, scales = "free_y")+
  labs(title = "Temporal trends in small copepod abundances in EA Regions",
       y = expression(bold("Abudance (m"^-3~")")),
       x=NULL,
       caption = "Colour and shape of points indicates sampling season:
       Winter (grey circle), Spring (green square), Summer (Orange diamond), Autumn (brown triangle).
       Line represents loess smoother (span = 0.9)") +
  theme(legend.position = "none",
        axis.title.y = element_text(face=2),
        strip.text = element_text(face=2)) -> pl

ggsave(plot = pl, filename = "figs/2407dd_timeseries/smCopTSByRgn.pdf",
       width = 12,height = 8,units = "in")
rm(pl)

## Fish
# by WB
dfw_lf %>% 
  ggplot(.,aes(
    x=sample.date,
    y= Fish_Mero*`Net.volume.sampled.(m3)`
  ))+
  geom_line(colour=1)+
  scale_fill_manual(values=c("slategray1","green","sienna2","brown"))+
  scale_shape_manual(values=c(21:24))+
  geom_point(aes(shape=DJF, fill=DJF), size=2)+
  ylim(0,NA)+
  facet_wrap(.~WB_lb, scales = "free_y")+
  labs(title = "Temporal trends in fish abundances in EA water bodies",
       y = expression(bold("Abudance (m"^3~")")),
       x=NULL,
       caption = "Colour and shape of points indicates sampling season:
       Winter (grey circle), Spring (green square), Summer (Orange diamond), Autumn (brown triangle)") +
  theme(legend.position = "none",
        axis.title.y = element_text(face=2),
        strip.text = element_text(face=2)) -> pl

ggsave(plot = pl, filename = "figs/2407dd_timeseries/fishTSByWB.pdf",
       width = 12,height = 8,units = "in")
rm(pl)

# by Region
dfw_lf %>% 
  ggplot(.,aes(
    x=sample.date,
    y= Fish_Mero*`Net.volume.sampled.(m3)`
  ))+
  # geom_line(colour=1)+
  scale_fill_manual(values=c("slategray1","green","sienna2","brown"))+
  scale_shape_manual(values=c(21:24))+
  geom_smooth(se=FALSE)+
  geom_point(aes(shape=DJF, fill=DJF), size=2)+
  facet_wrap(.~Region,scales = "free_y")+
  labs(title = "Temporal trends in fish abundances in EA Regions",
       y = bquote("Abudance (m"^3~")"),
       x=NULL,
       caption = "Colour and shape of points indicates sampling season:
       Winter (grey circle), Spring (green square), Summer (Orange diamond), Autumn (brown triangle)") +
  theme(legend.position = "none",
        axis.title.y = element_text(face=2),
        strip.text = element_text(face=2)) -> pl

ggsave(plot = pl, filename = "figs/2407dd_timeseries/fishTSByRgn.pdf",
       width = 12,height = 8,units = "in")
rm(pl)

### all fish vs all copepods
dfw_lf %>% 
  mutate(Cop_all = Cop_NYA+Cop_Sm+Cop_Ambi+Cop_Lg) %>% 
  mutate(Cop_all=Cop_all*`Net.volume.sampled.(m3)`,
         Fish_Mero=Fish_Mero*`Net.volume.sampled.(m3)`) %>% 
  filter(.,WB != "Solway Outer South") %>% 
  ggplot(.,aes(x=log(Cop_all+1), y= log(Fish_Mero+1)))+
  # geom_smooth(method="loess", colour=2,fill="pink",span = 0.9)+
  geom_smooth(method="lm", fill="lightblue")+
  geom_point(alpha=0.3)+
  coord_cartesian(ylim=c(0, NA))+
  facet_wrap(.~WB_lb, scales = "free")+
  ggpubr::stat_cor(method="pearson")+
  labs(title = "Relationship between log total copepod and total fish abundances in zooplankton assemblages in EA water bodies",
       y = expression(bold("Log fish abudance (m"^-3~")")),
       #y = bquote("Fish abudance (m"^-3~")"),
       #x = bquote("Total copepod abudance (m"^-3~")"),
       x = expression(bold("Log total copepod abudance (m"^-3~")")),
       caption=paste0("Samples gathered between ",format(min(dfw_lf$sample.date), "%d/%m/%Y")," & ",format(max(dfw_lf$sample.date), "%d/%m/%Y"),
                      "\nBlue lines inidate linear model")) +
  theme(legend.position = "none",
        axis.title = element_text(face=2),
        strip.text = element_text(face=2)) -> pl

ggsave(plot = pl, filename = "figs/2407dd_timeseries/copepodsFishByWB_allSmooths.pdf",
       width = 12,height = 8,units = "in")
rm(pl)

### all copepods vs temperature
dfw_lf %>% 
  mutate(Cop_all = Cop_NYA+Cop_Sm+Cop_Ambi+Cop_Lg) %>% 
  mutate(Cop_all=Cop_all*`Net.volume.sampled.(m3)`) %>% 
  filter(.,WB != "Solway Outer South") %>%
  ggplot(.,
         aes(x=log(as.numeric(`Chlorophyll : Acetone Extract_ug/l`)+1),
             y= log(Cop_all+1)))+
  geom_point()+
  facet_wrap(.~WB_lb, scales = "free")+
  geom_smooth(method="lm", fill="lightblue")+
  ggpubr::stat_cor(method="pearson",
                   label.y.npc="top", label.x.npc = "left")+
  labs(title = "Relationship between total copepod abundance and chlorophyll concentrations in EA water bodies",
       #y = bquote("log Copepod abundance (m"^-3~")"),
       y = expression(bold("log Copepod abundance (m"^-3~")")),
       #x = bquote("log Chlorophyll concentration (ug/l)"),
       x = expression(bold("log Chlorophyll concentration (ug/l)")),
       caption=paste0("Samples gathered between ",format(min(dfw_lf$sample.date), "%d/%m/%Y")," & ",format(max(dfw_lf$sample.date), "%d/%m/%Y"),
                      "\nBlue lines inidate linear model.")) +
  theme(legend.position = "none",
        strip.text = element_text(face=2)) -> pl

ggsave(plot = pl, filename = "figs/2407dd_timeseries/copepodsChlByWB_lm.pdf",
       width = 12,height = 8,units = "in")
rm(pl)
#####

# Carbon explorations ####

## flag which taxa we don't have carbon data for
left_join(df_lf_l, df_carb_summary, by="Aphia.ID") %>%
  filter(.,is.na(mdCPerIndiv_ug)) %>% 
  dplyr::select(.,c("Aphia.ID","Taxa")) %>% 
  distinct(.) %>% write.csv(.,
                            file="outputs/ZoopsWithoutCEstimate.csv",
                            row.names = FALSE)

# Append carbon content to taxon data & multiply abundance by Carbon
df_lf_l <- left_join(df_lf_l, df_carb_summary, by="Aphia.ID")

df_lf_l %>% 
  mutate(mn_carbTot_raw = AbundanceRaw*mnCPerIndiv_ug,
         md_carbTot_raw = AbundanceRaw*mdCPerIndiv_ug,
         mn_carbTot_m3 = Abund_m3*mnCPerIndiv_ug,
         md_carbTot_m3 = Abund_m3*mdCPerIndiv_ug
  ) -> df_lf_l

df_lf_l$sample.date <- as.Date(df_lf_l$sample.date, origin = "1899-12-30")
df_lf_l$yday <- lubridate::yday(df_lf_l$sample.date)
df_lf_l %>% relocate(.,yday, .after = sample.date) -> df_lf_l

### widen lifeforms data: remove counts, replace with Carbon
df_lf_w_C <- df_lf_l %>% 
  dplyr::select(.,c(date_site, sample.date, yday, BIOSYS.Code, WIMS.Code,
                    "Net.volume.sampled.(m3)", PRN, LF02, md_carbTot_m3,
                    WB,Region,WB_lb)) %>%
  mutate(md_carbTot_m3 = replace_na(md_carbTot_m3, 0)) %>% 
  group_by(across(c(!md_carbTot_m3))) %>% 
  summarise(md_carbTot_m3 = sum(md_carbTot_m3),.groups = "drop") %>%
  ungroup(.) %>% 
  pivot_wider(.,names_from = LF02, values_from = md_carbTot_m3,
              values_fill = 0) %>%
  rowwise() %>% #names(.)
  mutate(SUM = sum(c_across(-(1:10)))) %>% 
  ungroup()

df_lf_w_C %>% 
  ggplot(., aes(x=sample.date, y = log(SUM)))+ 
  geom_point()+
  geom_smooth(se=FALSE)+
  facet_wrap(.~WB_lb)+#, scale="free_y")+
  ylim(0,NA)+
  labs(title = "Trend in total carbon content within zooplankton assemblages by date in EA water bodies",
       y = expression(bold("Log total carbon content (ug m"^-3~")")),
       # y = expression(bold("Square root of total carbon content (ug m"^-3~")")),
       # y = expression(bold("Total carbon content (ug m"^-3~")")),
       x = "Date",
       caption=paste0("Samples gathered between ",format(min(dfw_lf$sample.date), "%d/%m/%Y")," & ",format(max(dfw_lf$sample.date), "%d/%m/%Y"),
                      "\nBlue lines inidate loess smooths.")) +
  theme(legend.position = "none",
        axis.title = element_text(face=2),
        strip.text = element_text(face=2)) -> pl

ggsave(plot = pl, filename = "figs/2407dd_timeseries/logtotCByDateByWB_Fixed_Y.pdf",
       width = 12,height = 8,units = "in")
rm(pl)

glob.mean <- mean(df_lf_w_C$SUM)
glob.mean.log <- mean(log(df_lf_w_C$SUM))
glob.median <- median(df_lf_w_C$SUM)
glob.median.log <- median(log(df_lf_w_C$SUM))

df_lf_w_C %>% 
  ggplot(.,aes(x=WB_lb, y=log(SUM), colour=Region))+
  geom_hline(yintercept = glob.mean.log,lty=2)+
  geom_boxplot(varwidth = TRUE,outlier.shape = NA)+
  geom_point(aes(group=Region), position=position_jitterdodge(),alpha=0.3)+
  scale_colour_manual(values = cbPalette2)+
  labs(title = "Total carbon content in zooplankton assemblages by EA water body",
       y="log(total carbon)",
       caption = paste0("Dashed line indicates global mean log carbon content across all water bodies",
                        "\nBox widths are proportional to the number of observations"))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face=2),
        axis.text.x = element_text(face=2)) -> pl
ggsave(plot = pl, filename = "figs/2407dd_timeseries/carbonByWB.pdf",
       width = 20,height = 12,units = "in")

df_lf_w_C %>%
  group_by(WB_lb) %>% 
  summarise(min=min(SUM),
            q1 = quantile(SUM, 0.25),
            median=median(SUM),
            mean = mean(SUM),
            q3 = quantile(SUM, 0.75),
            max=max(SUM),
            globalMean=glob.mean,
            globalMedian=glob.median,
            .groups = "drop"
            ) -> descriptives
write.csv(descriptives,
          file="outputs/descriptWB_SUMCarbon.csv",
          row.names = FALSE)

df_lf_w_C %>%
  mutate(SUM=log(SUM)) %>% group_by(WB_lb) %>% 
  summarise(min=min(SUM),
            q1 = quantile(SUM, 0.25),
            median=median(SUM),
            mean = mean(SUM),
            q3 = quantile(SUM, 0.75),
            max=max(SUM),
            globalMean=glob.mean.log,
            globalMedian=glob.median.log,
            .groups = "drop"
  ) -> logdescriptives
write.csv(logdescriptives,
          file="outputs/descriptWB_logSUMCarbon.csv",
          row.names = FALSE)
rm(descriptives,logdescriptives)

#####
# df_lf_l %>% 
#   dplyr::select(.,c(Pot.Number,date_site, sample.date, yday, BIOSYS.Code, WIMS.Code,
#                     ZooType, "Net.volume.sampled.(m3)", PRN, 
#                     md_carbTot_m3, WB,Region, WB_lb)) %>%
#   mutate(md_carbTot_m3 = replace_na(md_carbTot_m3, 0)) %>% 
#   group_by(across(c(!md_carbTot_m3))) %>% 
#   summarise(md_carbTot_m3 = sum(md_carbTot_m3),.groups = "drop") %>%
#   ungroup(.) %>% 
#   pivot_wider(.,names_from = ZooType, values_from = md_carbTot_m3,
#               values_fill = 0) %>%
#   rowwise() %>% #names(.)
#   mutate(SUM = sum(c_across(-(1:11)))) %>% 
#   ungroup() -> df_lf_w_C_zooType

df_lf_l %>% 
  dplyr::select(.,c(Pot.Number,date_site, sample.date, yday, BIOSYS.Code, WIMS.Code,
                    ZooType, "Net.volume.sampled.(m3)", PRN, LF02,
                    md_carbTot_m3, WB,Region, WB_lb)) %>%
  mutate(md_carbTot_m3 = replace_na(md_carbTot_m3, 0)) %>% 
  mutate(ZooType = if_else(is.na(ZooType) | ZooType == "NYA", LF02, ZooType)) %>%
  dplyr::select(.,-LF02) %>% 
  group_by(across(c(!md_carbTot_m3))) %>% 
  summarise(md_carbTot_m3 = sum(md_carbTot_m3),.groups = "drop") %>%
  ungroup(.) %>% 
  pivot_wider(.,names_from = ZooType, values_from = md_carbTot_m3,
              values_fill = 0) %>%
  rowwise() %>% #names(.)
  mutate(SUM = sum(c_across(-(1:11)))) %>% 
  ungroup() -> df_lf_w_C_zooType

df_lf_w_C_zooType %>% 
  dplyr::select(.,-SUM) %>%
  pivot_longer(.,cols=12:ncol(.),names_to = "ZooType",values_to = "TotCarb") %>%
  filter(.,TotCarb != 0) %>% 
  ggplot(., aes(x=ZooType, y=log(TotCarb)))+
  geom_boxplot(varwidth = FALSE,outlier.shape = NA)+
  coord_flip()+
  facet_wrap(.~Region)+
  # facet_wrap(.~WB_lb)+
  geom_point(aes(group=ZooType), position=position_jitter(width = 0.2),alpha=0.2)+
  scale_colour_manual(values = cbPalette2)+
  scale_x_discrete(limits=rev)+
  labs(
    title = "Total carbon content by zooplankton type by EA region",
    y="log(Total carbon)"
       )+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        strip.text = element_text(face=2),
        axis.title.x = element_text(face=2),
        axis.text.x = element_text(face=2),
        axis.text.y = element_text(face=2)) -> pl

ggsave(plot = pl, filename = "figs/2407dd_timeseries/carbonByZooTypeRegion.pdf",
       width = 20,height = 12,units = "in")
rm(pl)

df_lf_l %>% #names(.)
  mutate(.,LF03 = if_else(Copepod == "Y", paste0("Cop_",CopSize),ZooType)) %>%
  mutate(.,LF03 = if_else(is.na(LF03) | LF03 == "NYA", LF02, LF03)) %>%
  dplyr::select(.,c(Pot.Number, sample.date, yday, BIOSYS.Code, WIMS.Code,
                    LF03, "Net.volume.sampled.(m3)", PRN, 
                    md_carbTot_m3, WB,Region,WB_lb)) %>% 
  mutate(md_carbTot_m3 = replace_na(md_carbTot_m3, 0)) %>% 
  filter(.,md_carbTot_m3 != 0) %>% 
  group_by(across(c(!md_carbTot_m3))) %>% 
  summarise(md_carbTot_m3 = sum(md_carbTot_m3),.groups = "drop") %>%
  ungroup(.) %>%
  ggplot(., aes(x=LF03, y=log(md_carbTot_m3)))+
  geom_boxplot(varwidth = FALSE,outlier.shape = NA)+
  coord_flip()+
  facet_wrap(.~WB_lb)+
  # facet_wrap(.~WB_lb)+
  geom_point(aes(group=LF03), position=position_jitter(width = 0.2),alpha=0.2)+
  scale_colour_manual(values = cbPalette2)+
  scale_x_discrete(limits=rev)+
  labs(
    title = "Total carbon content by zooplankton type by EA water body",
    y="log(Total carbon)"
  )+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        strip.text = element_text(face=2),
        axis.title.x = element_text(face=2),
        axis.text.x = element_text(face=2),
        axis.text.y = element_text(face=2,
                                   size = 7)) -> pl

ggsave(plot = pl, filename = "figs/2407dd_timeseries/carbonByZooTypeWB_v2.pdf",
       width = 20,height = 12,units = "in")
rm(pl)

# SEAHORSE FOOD ####
## total carbon
df_lf_w_C %>% 
  filter(.,WB_lb %in% c("Sth_Solent","Sth_SotonWtr")) %>% 
  filter(.,WIMS.Code != "Y00017477") %>% 
  mutate(.,WIMS.Code = factor(WIMS.Code, levels = c("G0003572",
                                                    "Y0017477",
                                                    "Y0004367",
                                                    "G0003532"))) %>% 
  ggplot(.,aes(x=BIOSYS.Code,y=SUM, colour=WB_lb))+
  geom_hline(yintercept = seq(from=0, to=100000, by=10000), colour = "grey",
             linetype=2)+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(width = 0.25)+
  labs(title = "Total carbon content of zooplankters recorded in the Solent and Southampton Water",
       x="BIOSYS site code",
       y= "Total carbon (ug C/m3)")+
  scale_y_continuous(breaks = c(0,50000,100000),labels = scales::comma_format())+
  coord_flip()+
  theme(legend.title = element_blank(),
        axis.title = element_text(face=2)) -> pl

ggsave(plot = pl, filename = "figs/2407dd_timeseries/carbonTot_Soton.pdf",
       width = 20,height = 12,units = "in");rm(pl)

## by taxon
dfw_lf %>% 
  filter(.,WB_lb %in% c("Sth_Solent","Sth_SotonWtr")) %>% names(.)

df_lf_l %>% 
  filter(.,WB_lb %in% c("Sth_Solent","Sth_SotonWtr")) %>%
  mutate(.,label = ifelse(LF02 == "Cop_Sm","Cop_Sm",
                          ifelse(LF02 == "Cop_Lg","Cop_Lg",
                                 ifelse(LF02 == "Cop_Ambi","Cop_Ambi",
                                        ifelse(LF02 == "Cop_NYA","Cop_NYA",
                                               ifelse(Order == "Mysida","Mysida","X")))))) %>% 
  filter(.,label != "X") %>% 
  dplyr::select(.,c("BIOSYS.Code","WIMS.Code","sample.date","PRN",
                    "WB_lb","label","Abund_m3","md_carbTot_m3")) %>% 
  group_by(BIOSYS.Code, WIMS.Code, sample.date, PRN, WB_lb, label) %>%
  summarise(
    Abund_m3 = sum(Abund_m3, na.rm = TRUE),
    md_carbTot_m3 = sum(md_carbTot_m3, na.rm = TRUE)
  ) %>%
  ungroup() %>% 
  ggplot(.,aes(x= BIOSYS.Code,y=log(md_carbTot_m3+1), colour=WB_lb))+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(.~label, scales = "free_y")+
  geom_jitter(width = 0.25)+
  labs(title = "Total carbon content of selected zooplankters recorded in the Solent and Southampton Water",
       x="BIOSYS site code",
       y= "Log (n+1) carbon content (ug C/m3)")+
  theme(legend.title = element_blank(),
        axis.title = element_text(face=2))+
  coord_flip() -> pl

ggsave(plot = pl, filename = "figs/2407dd_timeseries/log_CarbonTot_Soton.pdf",
       width = 20,height = 12,units = "in");rm(pl)

# how many taxa do we not have C data for? ####
left_join(df_tx_l, df_carb_summary, by="Aphia.ID") %>% 
  mutate(.,WB_lb = WB_lb) %>%
  mutate(.,cdat=ifelse(is.na(.$mdCPerIndiv_ug),"N","Y")) %>%
  mutate(Abund_m3 = ifelse(is.na(.$"Net.volume.sampled.(m3)"),AbundanceRaw,
                           .$AbundanceRaw*.$"Net.volume.sampled.(m3)")) %>% 
  dplyr::select(.,c(sample.date,BIOSYS.Code,WIMS.Code,PRN,
                    Abund_m3,Region,WB_lb,cdat)) %>% 
  group_by(sample.date,BIOSYS.Code,WIMS.Code,PRN,
           Region,WB_lb,cdat) %>% 
  summarise(totAbnd=sum(Abund_m3),.groups = "drop") %>% 
  ggplot(.,aes(cdat,
               y=log(totAbnd+1),
               # y=totAbnd,
               colour=cdat))+
  # geom_bar(stat = "identity")+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(width=.2)+
  labs(
    title="Total zooplankton abundances in samples gathered in EA water bodies",
    subtitle = "Abundances split between taxa for which we do (Y) and do not (N) have carbon content estimates",
    # y="Total zooplankton abundance (individuals/m3)",
    y="Log (n+1) zooplankton abundance",
    x=NULL
    )+
  facet_wrap(.~WB_lb, scales = "free_y")+
  theme(
    axis.title.y = element_text(face=2),
    axis.text.x = element_text(face=2),
    strip.text = element_text(face=2),
    legend.position = "none"
  ) -> pl
  
ggsave(plot = pl, filename = "figs/2407dd_timeseries/LogAbundTot_Missing_C.pdf",
       width = 20,height = 12,units = "in");rm(pl)

# Prevalence of missing Carbon taxa ####
left_join(df_tx_l, df_carb_summary, by="Aphia.ID") %>% 
  mutate(.,WB_lb = WB_lb) %>% 
  mutate(.,cdat=ifelse(is.na(.$mdCPerIndiv_ug),"N","Y")) %>%
  mutate(Abund_m3 = ifelse(is.na(.$"Net.volume.sampled.(m3)"),AbundanceRaw,
                           .$AbundanceRaw*.$"Net.volume.sampled.(m3)")) %>% 
  dplyr::select(.,c(sample.date,BIOSYS.Code,WIMS.Code,PRN,
                    Abund_m3,Region,WB_lb,cdat)) %>% 
  mutate(preval=1) %>% dplyr::select(.,-Abund_m3) %>% 
  group_by(sample.date,BIOSYS.Code,WIMS.Code,PRN,
           Region,WB_lb,cdat) %>% 
  summarise(preval=sum(preval),.groups = "drop") %>% 
  mutate(.,WB_lb = factor(WB_lb, levels = c("NE_Nrthmb Nth",
    "NE_Farne Is","NE_Tees", "Ang_Yorks Sth", "Ang_Lincs", "Ang_Wash Out",
    "Ang_Blckw Out", "Thm_Thm Low", "Sth_Kent Sth", "Sth_Solent", "Sth_Soton Wtr",
    "SW_Cornw Nth", "SW_Brnstp B","SW_Brist Ch In Sth","NW_Mersey Mth", "NW_Solway O"))) %>% 
  filter(.,cdat=="N") %>% 
  ggplot(., aes(x=preval))+
  geom_histogram(fill="grey",colour=1)+
  facet_wrap(.~WB_lb)+
  labs(title = "Prevalence of zooplankton taxa which do not have an estimate of carbon content",
       subtitle="Histograms show the number of taxa without estimate of carbon content within each sample across EA water bodies",
       x="Number of taxa within sample with missing carbon estimates",
       y="Frequency")+
  theme(
    axis.title.y = element_text(face=2),
    axis.title.x = element_text(face=2),
    axis.text.x = element_text(face=2),
    strip.text = element_text(face=2)
  ) -> pl
ggsave(plot = pl, filename = "figs/2407dd_timeseries/OccurencesMissing_C.pdf",
       width = 20,height = 12,units = "in");rm(pl)

df_tx_l %>% 
  dplyr::select(.,c(Pot.Number,PRN)) %>% distinct() %>% count() -> samples

left_join(df_tx_l,df_carb_summary, by="Aphia.ID") %>% #names() 
  mutate(Abund_m3 = ifelse(is.na(.$"Net.volume.sampled.(m3)"),AbundanceRaw,
                           .$AbundanceRaw*.$"Net.volume.sampled.(m3)")) %>% 
  dplyr::select(Pot.Number, PRN, sample.date, Region, Aphia.ID,
                Taxa,Abund_m3,mdCPerIndiv_ug) %>% 
  mutate(present = 1) %>% 
  mutate(carbon = ifelse(is.na(mdCPerIndiv_ug), "N", "Y")) %>% 
  dplyr::select(-mdCPerIndiv_ug, -Pot.Number, -PRN, -sample.date, -Region,-Abund_m3) %>%
  group_by(across(-present)) %>%  # Group by all columns except 'present'
  summarise(occurences = sum(present), 
            proportionOfSamples = as.numeric(occurences/samples), 
            .groups = "drop") %>% 
  ungroup() %>% 
  arrange(-proportionOfSamples) %>% 
  rename(.,carbonEstimateAvailable = carbon) %>% #View(.)
  write.csv(.,
            file="outputs/PrevalenceOfMissingCarbonTax0.csv",
            row.names = FALSE)
rm(samples)

left_join(df_tx_l,df_carb_summary, by="Aphia.ID") %>% #names() 
  mutate(Abund_m3 = ifelse(is.na(.$"Net.volume.sampled.(m3)"),AbundanceRaw,
                           .$AbundanceRaw*.$"Net.volume.sampled.(m3)")) %>% 
  dplyr::select(Pot.Number, PRN, sample.date, Region, Aphia.ID,
                Taxa,Abund_m3,mdCPerIndiv_ug) %>% 
  # mutate(present = 1) %>% 
  mutate(carbon = ifelse(is.na(mdCPerIndiv_ug), "N", "Y")) %>% 
  dplyr::select(-mdCPerIndiv_ug, -Pot.Number, -PRN, -sample.date, -Region) %>%
  group_by(across(-Abund_m3)) %>%  # Group by all columns except 'present'
  summarise(totalAbund = sum(Abund_m3), 
            .groups = "drop") %>% 
  ungroup() %>% 
  arrange(-totalAbund) %>% 
  rename(.,carbonEstimateAvailable = carbon) %>% 
  mutate(.,proportionOfIndividuals = totalAbund/sum(totalAbund)*100) %>% 
  write.csv(.,
            file="outputs/ProportionOfCounts_MissingCarbonTax0.csv",
            row.names = FALSE)
