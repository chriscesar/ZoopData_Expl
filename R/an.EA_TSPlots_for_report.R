# an.EA_TSPlots_for_report.R ####
# Comparison of Joint Species Distribution Model approaches for zooplankton
# lifeforms data

# 00 set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "tictoc","gllvm","patchwork", "lubridate")
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
df_lf_l <- df_tx0
df_lf_w <- df_tx_w
rm(df_tx_w,df_tx,df_tx_100um,df_tx0)
toc(log = TRUE)

## taxon data
tic("load data sets: taxon data")
source("R/imp.load_data_lifeforms.R")
df_tx_l <- df_tx0
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

### create Region_WB variable
LFRegion <- dfw_lf$Region
LFWB <- dfw_lf$WB

WB_lb1 <- ifelse(LFRegion == "Southern","Sth",
                    ifelse(LFRegion == "Thames","Thm",
                           ifelse(LFRegion == "Anglian","Ang",
                                  ifelse(LFRegion == "NWest","NW",
                                         ifelse(LFRegion == "NEast","NE",
                                                ifelse(LFRegion == "SWest","SW",NA)
                                         )))))

WB_lb2 <- ifelse(LFWB == "Solent","Solent",
                    ifelse(LFWB == "SOUTHAMPTON WATER","Soton Wtr",
                           ifelse(LFWB == "Solway Outer South","Solway O",
                                  ifelse(LFWB == "THAMES LOWER","Thm Low",
                                         ifelse(LFWB == "Blackwater Outer","Blckw Out",
                                                ifelse(LFWB == "Cornwall North","Cornw Nth",
                                                       ifelse(LFWB == "Barnstaple Bay","Brnstp B",
                                                              ifelse(LFWB == "Kent South","Kent Sth",
                                                                     ifelse(LFWB == "Mersey Mouth","Mersey Mth",
                                                                            ifelse(LFWB == "Wash Outer","Wash Out",
                                                                                   ifelse(LFWB == "Lincolnshire","Lincs",
                                                                                          ifelse(LFWB == "Yorkshire South","Yorks Sth",
                                                                                                 ifelse(LFWB == "TEES","Tees",
                                                                                                        ifelse(LFWB == "Northumberland North","Nrthmb Nth",
                                                                                                               ifelse(LFWB == "Farne Islands to Newton Haven","Farne Is",
                                                                                                                      ifelse(LFWB == "Bristol Channel Inner South","Brist Ch In Sth",
                                                                                                                             NA)))))))))))
                                         )))))
WB_lb <- paste0(WB_lb1,"_",WB_lb2)
rm(WB_lb1,WB_lb2)

dfw_lf$WB_lb <- WB_lb
dfw_lf %>% relocate(WB_lb,.after = WB) -> dfw_lf
dfw_tx$WB_lb <- WB_lb
dfw_tx %>% relocate(WB_lb,.after = WB) -> dfw_tx
########

#####

# PLOTS! ####
## Lifeforms ####

## Copepod size classes ###
# by wb
dfw_lf %>%
  dplyr::select(., c(Region,WB,WB_lb,sample.date,`Net.volume.sampled.(m3)`,Cop_Sm, Cop_Ambi, Cop_Lg,Cop_NYA)) %>%
  mutate(Cop_Sm=(`Net.volume.sampled.(m3)`*Cop_Sm)/max(`Net.volume.sampled.(m3)`*Cop_Sm),
         Cop_Ambi=(`Net.volume.sampled.(m3)`*Cop_Ambi)/max(`Net.volume.sampled.(m3)`*Cop_Ambi),
         Cop_Lg=(`Net.volume.sampled.(m3)`*Cop_Lg)/max(`Net.volume.sampled.(m3)`*Cop_Lg),
         Cop_NYA=(`Net.volume.sampled.(m3)`*Cop_NYA)/max(`Net.volume.sampled.(m3)`*Cop_NYA)) %>% 
  pivot_longer(.,cols = Cop_Sm:Cop_NYA,names_to = "lf",values_to = "abund") %>%
  ggplot(.,aes(x=sample.date, y=abund, colour=lf)) +
  geom_line()+
  facet_wrap(.~WB_lb)+
  labs(x=NULL,
       y="Relative abundance",
       title="Relative abundance of copepod size classes recorded in EA water bodies",
       caption="Relative abundance of a given size class is the value recorded in a sample
       divided by the maximum recorded value for that size class across all samples")+
  theme(legend.title = element_blank()) -> pl

ggsave(plot = pl, filename = "figs/2407dd_timeseries/copSizesTSByWB.pdf",
       width = 12,height = 8,units = "in")
rm(pl)

# by Region
dfw_lf %>%
  dplyr::select(., c(Region,WB,WB_lb,sample.date,`Net.volume.sampled.(m3)`,Cop_Sm, Cop_Ambi, Cop_Lg,Cop_NYA)) %>%
  mutate(Cop_Sm=(`Net.volume.sampled.(m3)`*Cop_Sm)/max(`Net.volume.sampled.(m3)`*Cop_Sm),
         Cop_Ambi=(`Net.volume.sampled.(m3)`*Cop_Ambi)/max(`Net.volume.sampled.(m3)`*Cop_Ambi),
         Cop_Lg=(`Net.volume.sampled.(m3)`*Cop_Lg)/max(`Net.volume.sampled.(m3)`*Cop_Lg),
         Cop_NYA=(`Net.volume.sampled.(m3)`*Cop_NYA)/max(`Net.volume.sampled.(m3)`*Cop_NYA)) %>% 
  pivot_longer(.,cols = Cop_Sm:Cop_NYA,names_to = "lf",values_to = "abund") %>%
  ggplot(.,aes(x=sample.date, y=abund, colour=lf)) +
  # geom_line()+
  geom_smooth(se=FALSE)+
  # geom_point()+
  facet_wrap(.~Region)+
  labs(x=NULL,
       y="Relative abundance",
       title="Relative abundance of copepod size classes recorded in EA water bodies",
       caption="Relative abundance of a given size class is the value recorded in a sample
       divided by the maximum recorded value for that size class across all samples")+
  theme(legend.title = element_blank()) -> pl

ggsave(plot = pl, filename = "figs/2407dd_timeseries/copSizeSmoothTSByRgn.pdf",
       width = 12,height = 8,units = "in")
rm(pl)

## small copepods
# by WB
dfw_lf %>% 
  ggplot(.,aes(
    x=sample.date,
    y= Cop_Sm*`Net.volume.sampled.(m3)`
    ))+
  geom_line(colour=1)+
  scale_fill_manual(values=c("slategray1","green","sienna2","brown"))+
  scale_shape_manual(values=c(21:24))+
  geom_point(aes(shape=DJF, fill=DJF), size=2)+
  facet_wrap(.~WB_lb)+
  labs(title = "Temporal trends in small copepod abundances in EA water bodies",
       y = bquote("Abudance (m"^3~")"),
       x=NULL,
       caption = "Colour and shape of points indicates sampling season:
       Winter (grey circle), Spring (green square), Summer (Orange diamond), Autumn (brown triangle)") +
  theme(legend.position = "none")-> pl
ggsave(plot = pl, filename = "figs/2407dd_timeseries/smCopTSByWB.pdf",
       width = 12,height = 8,units = "in")
rm(pl)

# by Region
dfw_lf %>% 
  ggplot(.,aes(
    x=sample.date,
    y= Cop_Sm*`Net.volume.sampled.(m3)`
  ))+
  # geom_line(colour=1)+
  scale_fill_manual(values=c("slategray1","green","sienna2","brown"))+
  scale_shape_manual(values=c(21:24))+
  geom_smooth(se=FALSE)+
  geom_point(aes(shape=DJF, fill=DJF), size=2)+
  facet_wrap(.~Region)+
  labs(title = "Temporal trends in small copepod abundances in EA Regions",
       y = bquote("Abudance (m"^3~")"),
       x=NULL,
       caption = "Colour and shape of points indicates sampling season:
       Winter (grey circle), Spring (green square), Summer (Orange diamond), Autumn (brown triangle)") +
  theme(legend.position = "none")-> pl
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
  facet_wrap(.~WB_lb)+
  labs(title = "Temporal trends in fish abundances in EA water bodies",
       y = bquote("Abudance (m"^3~")"),
       x=NULL,
       caption = "Colour and shape of points indicates sampling season:
       Winter (grey circle), Spring (green square), Summer (Orange diamond), Autumn (brown triangle)") +
  theme(legend.position = "none")-> pl
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
  facet_wrap(.~Region)+
  labs(title = "Temporal trends in fish abundances in EA Regions",
       y = bquote("Abudance (m"^3~")"),
       x=NULL,
       caption = "Colour and shape of points indicates sampling season:
       Winter (grey circle), Spring (green square), Summer (Orange diamond), Autumn (brown triangle)") +
  theme(legend.position = "none")-> pl
ggsave(plot = pl, filename = "figs/2407dd_timeseries/fishTSByRgn.pdf",
       width = 12,height = 8,units = "in")
rm(pl)
