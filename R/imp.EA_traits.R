### imp.EA_traits.R ####
#### Import EA trait data ####

### useful links for plotting MDS with ggplot
# https://chrischizinski.github.io/rstats/vegan-ggplot2/
# https://www.youtube.com/watch?v=Y0GI34S-ZMI

## set up ####
#### install required packages ####
#### set local package library 
#libfolder <- "U:/Rlibrary"
# libfolder <- "M:/R/Shared Library"
# .libPaths(libfolder)

#### check and install required packages ####
ptm <- Sys.time()
req_packages <- c("lubridate",# working with dates
                  "seas", #more dates
                  "tidyverse",# general data manipulation
                  "vegan",# NMDS, multivariate analysis of ecological data
                  "vegan3d",# as above, but with 3D figs
                  "mvabund", #multivariate abundance data analyses
                  "ecoCopula",#model based ordination/graphical modelling
                  "ggthemes",# sensible visualisation styles
                  "openxlsx",# read data from xlsx
                  "MASS",# fit a negative binomial glm
                  "gclus",#clustering of data
                  "corrplot",#correlation plots
                  "performance",# model checking
                  "patchwork"
)

new_packages <- req_packages[!(req_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages,library=libfolder,type="binary")
Sys.time() - ptm;rm(ptm,req_packages,new_packages)

#### load packages ####
ld_pkgs <- c("tidyverse","MASS","lubridate","vegan","mvabund","seas",
             "ecoCopula","performance","gclus","corrplot", "patchwork")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

#### set universals ####
source("R/folder.links.R") ## data folders
perms <- 999 ### number of permutations to run for multivariate analyses
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

### load data ####
### taxon data
df_tx0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,
                                               "processedData/MBA_Returns_Amalgamated_USE.xlsx"),
                                        sheet="outR04_LF"))
### WIMS data
### await update of WIMS data
df_wims0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,


                                                 "/WIMS_Extract_WaterQuality_Zoop_Samples_230809.xlsx"),
                                          sheet="allDat"))

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

### remove odd data
df_tx <- as_tibble(df_tx0) ### create new data (keep df0 as 'raw')

### convert dates
df_tx$sample.date <- as.Date(df_tx$sample.date, origin = "1899-12-30")

# Remove 100 µm data [OPTIONAL] ####
df_tx_100um <- df_tx %>% 
  filter(str_starts(Sample.comments,"100um"))

df_tx %>% 
  filter(!str_starts(Sample.comments,"100um")) -> df_tx

###widen data & fill NAs with 0s ####
df_tx %>% 
  dplyr::select(.,-c("Aphia.ID","AbundanceRaw","Taxa","Category":"Unallocated")) %>% #drop unneeded cols
  group_by(across(c(-Abund_m3))) %>% 
  summarise(Abund_m3=sum(Abund_m3),
            .groups="drop") %>% 
  pivot_wider(names_from = "LF0",values_from = "Abund_m3",
              values_fill = 0) -> df_tx_w

### join & save data ####  
dfw <- left_join(df_tx_w,df_wims_w,by="PRN")

# write.csv(dfw,file=paste0(datfol,"processedData/","zoopWideTraitAbund_m3_WIMS_USE.csv"),row.names = FALSE)
# write.csv(df_tx_w,file=paste0(datfol,"processedData/","zoopWideTraitAbund_m3_taxOnly_USE.csv"),row.names = FALSE)

###############################
## look at taxon data only ####
###############################

### tax number
df_tx_w$S <- vegan::specnumber(df_tx_w[,-c(1:20)])
df_tx_w$N <- rowSums(df_tx_w[,-c(1:20,length(df_tx_w))])

### how often do different taxa appear?
xx <- df_tx_w %>% 
  dplyr::select(-c(1:20)) %>% ###remove metadata info
  dplyr::select(-last_col()) %>% ###remove taxon abundance
  dplyr::select(-last_col()) %>% ###remove taxon richness
  dplyr::mutate_all(~ifelse(. != 0, 1, .))

xx$WB <- df_tx_w$WB
xx$Region <- df_tx_w$Region

xx$WB_lb1 <- ifelse(xx$Region == "Southern","Sth",
                    ifelse(xx$Region == "Thames","Thm",
                           ifelse(xx$Region == "Anglian","Ang",
                                  ifelse(xx$Region == "NWest","NW",
                                         ifelse(xx$Region == "NEast","NE",
                                                ifelse(xx$Region == "SWest","SW",NA)
                                         )))))

xx$WB_lb2 <- ifelse(xx$WB == "Solent","Solent",
                    ifelse(xx$WB == "SOUTHAMPTON WATER","Soton Wtr",
                           ifelse(xx$WB == "THAMES LOWER","Thm Low",
                                  ifelse(xx$WB == "Blackwater Outer","Blckw Out",
                                         ifelse(xx$WB == "Cornwall North","Cornw Nth",
                                                ifelse(xx$WB == "Barnstaple Bay","Brnstp B",
                                                       ifelse(xx$WB == "Kent South","Kent Sth",
                                                              ifelse(xx$WB == "Mersey Mouth","Mersey Mth",
                                                                     ifelse(xx$WB == "Wash Outer","Wash Out",
                                                                            ifelse(xx$WB == "Lincolnshire","Lincs",
                                                                                   ifelse(xx$WB == "Yorkshire South","Yorks Sth",
                                                                                          ifelse(xx$WB == "TEES","Tees",
                                                                                                 ifelse(xx$WB == "Northumberland North","Nrthmb Nth",
                                                                                                        ifelse(xx$WB == "Farne Islands to Newton Haven","Farne Is",
                                                                                                               ifelse(xx$WB == "Bristol Channel Inner South","Brist Ch In Sth",
                                                                                                                      NA))))))))))
                                         )))))
xx$WB_lb <- paste0(xx$WB_lb1,"_",xx$WB_lb2)
xx$WB_lb1 <- NULL; xx$WB_lb2 <- NULL

xx %>% 
  # dplyr::select(-c(WB,WB_lb)) %>% ## plot by Region
  dplyr::select(-c(WB,Region)) %>% ## plot by WB
  # group_by(Region) %>% ## plot by Region
  group_by(WB_lb) %>%  ## plot by WB
  summarise(across(everything(),sum)) %>% # -> df_tx_prev_Rgn
  ### need to pivot to longer for plotting
  # pivot_longer(!Region, names_to="taxon",values_to="preval") %>% #by Region
  pivot_longer(!WB_lb, names_to="taxon",values_to="preval") %>% #by WB
  dplyr::filter(preval >0) -> xxprev

xxprev %>% 
  ggplot(.,aes(x=preval))+
  geom_histogram(fill = "lightgrey", colour = 1,bins=40)+
  # facet_wrap(.~Region)+ #plot by Region
  facet_wrap(.~WB_lb)+ #plot by WB
  # labs(title="Prevalence of zooplankton traits recorded by Region", #by Region
  labs(title="Prevalence of zooplankton traits recorded by Region_Water body", # by WB
       # caption=paste0("Prevalence is the number of times that individual traits have been recorded in a given region (at any abundance)\nZero values excluded\n",#Region
       caption=paste0("Prevalence is the number of times that individual traits have been recorded in a given water body (at any abundance)\nZero values excluded\n",#WB
                      "Samples gathered between ",min(dfw$sample.date)," & ",max(dfw$sample.date)),
       y = "Count",
       x = "Prevalence")

# ggsave(filename = "figs/zoopTraitsPrevalence_by_Region.pdf",width = 12,height = 7,units = "in")
ggsave(filename = "figs/zoopTraitsPrevalence_by_WB.pdf",width = 12,height = 10,units = "in")

#### quick ordinations ####
df_tx_w %>% 
  dplyr::select(-c(1:20)) %>% ###remove metadata info
  dplyr::select(-last_col()) %>% ###remove taxon abundance
  dplyr::select(-last_col()) %>% ###remove taxon richness
  dplyr::select_if(~ !is.numeric(.) || sum(.) !=0) %>% ### drop cols summing to 0
  filter(rowSums(across(where(is.numeric)))!=0) -> dftmp ###remove 'empty' rows

### NMDS ####
ptm <- Sys.time()###
set.seed(pi);ord <-   vegan::metaMDS(dftmp,trymax = 500)
Sys.time() - ptm;rm(ptm)
plot(ord)

### extract site info
df_tx_w %>% 
  mutate_at(c("Net.volume.sampled.(m3)",
              "Time.of.sampling.(GMT)?"),
            as.character) %>% 
  dplyr::select(-last_col()) %>% ###remove taxon abundance
  dplyr::select(-last_col()) %>% ###remove taxon richness
  dplyr::select_if(~ !is.numeric(.) || sum(.) !=0) %>% ### drop cols summing to 0
  dplyr::filter(rowSums(across(where(is.numeric)))!=0)  %>% ###remove 'empty' rows
  dplyr::select(c(1:21)) %>% ###remove taxon info
  # dplyr::mutate(mesh=substr(Sample.comments,1,5)) ->   scores_site#check this###
  dplyr::mutate(mesh=ifelse(grepl("200um",Sample.comments),"200um",
                            ifelse(grepl("100um",Sample.comments),"100um",NA))) -> scores_site

#### plot in ggplot2 ####
tmp_sites <- as.data.frame(scores(ord,display = "site"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
scores_site$NMDS1 <- tmp_sites$NMDS1 #location of individual samples in NMDS space
scores_site$NMDS2 <- tmp_sites$NMDS2 #location of individual samples in NMDS space
scores_site$yymm <- as.factor(format(scores_site$sample.date, "%Y-%m")) 
scores_site$DJF <- as.factor(mkseas(scores_site$sample.date, width="DJF"))#convert dates to 3month seasonal block
rm(tmp_sites)

# Using the scores function from vegan to extract the species scores and
# convert to a data.frame
scores_species <- as.data.frame(scores(ord,display = "species"))
scores_species$lb0 <-  make.cepnames(row.names(scores_species))#shorten names
scores_species$lb <-  row.names(scores_species)

#### generate mean centroids by WB and Region ####
scores_site %>% 
  group_by(WB) %>%
  summarise(mn_ax1_WB=mean(NMDS1),mn_ax2_WB=mean(NMDS2)) %>%
  ungroup() -> centr

scores_site <- left_join(scores_site,centr,by="WB");rm(centr)

scores_site %>% 
  group_by(Region) %>% 
  summarise(mn_ax1_Rgn=mean(NMDS1),mn_ax2_Rgn=mean(NMDS2)) %>%
  ungroup() -> centr

scores_site <- left_join(scores_site,centr,by="Region");rm(centr)

pl <- ggplot()+
  geom_hline(colour="grey",yintercept = 0, lty=2)+
  geom_vline(colour="grey",xintercept = 0, lty=2)+
  geom_text(data=scores_species, aes(x = NMDS1, y=NMDS2, label=lb),
            size=3,
            alpha=0.2)+
  ## segment by WB (comment out if not using)##
  # geom_segment(data=scores_site,aes(x=NMDS1,y=NMDS2,
  #                                   colour=WB,
  #                                   # colour=Region,
  #                                   xend=mn_ax1_WB,yend=mn_ax2_WB),
  #              show.legend = FALSE)+
  # geom_point(data=scores_site,
  #            aes(x=NMDS1, y=NMDS2,
  #                fill = WB,
  #                shape=WB),
  #            size=3)+
# scale_fill_manual(values=c(cbPalette))+
# scale_colour_manual(values=c(cbPalette))+
# scale_shape_manual(values = rep(c(25:21),3))+
# segment by Region (comment out if not using) ##
geom_segment(data=scores_site,aes(x=NMDS1,y=NMDS2,
                                  colour=Region,
                                  xend=mn_ax1_Rgn,yend=mn_ax2_Rgn),
             show.legend = FALSE)+
  geom_point(data=scores_site, show.legend=TRUE,
             aes(x=NMDS1, y=NMDS2,
                 fill = Region,
                 shape = Region),
             size=3)+
  scale_fill_manual(values=c(cbPalette))+
  scale_colour_manual(values=c(cbPalette))+
  # scale_shape_manual(values = rep(c(25:21),each=2))+
  scale_shape_manual(values = rep(c(24,25,23),each=2))+
  coord_equal()+
  labs(title="Non-metric Multidimensional Scaling of zooplankton lifeforms",
       subtitle="Colours & shapes indicate region",
       caption=paste0("Stress = ",round(ord$stress,3),
                      "\nOrdination based on Lifeform data",
                      "\nSamples gathered between ",
                      min(dfw$sample.date)," & ",max(dfw$sample.date)))+
  theme(legend.title = element_blank(),
        axis.title = element_text(face="bold"));pl

ggsave(filename = "figs/nmds_by_Region_Traits.pdf",width = 12,height = 12,units = "in",
       plot=pl);rm(pl)

#### plot by YY_MM ####
# Create an empty list to store individual plots
# Create an empty list to store individual plots
plot_list <- list()

# Define a consistent color and shape palette
consistent_palette <- c("Anglian" = "#999999",
                        "NEast" = "#E69F00",
                        "NWest"="#56B4E9",
                        "Southern" = "#CC79A7",
                        "SWest"="#009E73",
                        "Thames"="#F0E442")
consistent_shapes <- c("Anglian" = 24,
                       "NEast" = 24,
                       "NWest"=25,
                       "Southern" = 25,
                       "SWest"= 23,
                       "Thames"= 23)

# Iterate over each level of the 'yymm' factor
for (level in levels(scores_site$yymm)) {
  # Subset the data for the current 'yymm' level
  subset_data <- scores_site[scores_site$yymm == level, ]
  
  # Create a plot for the current level
  current_plot <- ggplot() +
    geom_hline(colour="grey",yintercept = 0, lty=2)+
    geom_vline(colour="grey",xintercept = 0, lty=2)+
    geom_text(data = scores_species, aes(x = NMDS1, y = NMDS2, label = lb), size = 3, alpha = 0.2) +
    #geom_segment(data = subset_data, aes(x = NMDS1, y = NMDS2, colour = Region, xend = mn_ax1_Rgn, yend = mn_ax2_Rgn), show.legend = FALSE) +
    geom_point(data = subset_data,
               aes(x = NMDS1, y = NMDS2, fill = Region, shape = Region),
               show.legend=FALSE,
               size = 3) +
    scale_fill_manual(values = consistent_palette) +
    scale_colour_manual(values = consistent_palette) +
    scale_shape_manual(values = consistent_shapes) +
    coord_equal() +
    labs(title = level)+  # Set the title to the 'yymm' level
    theme(axis.title=element_blank(),
          axis.text=element_blank())
  
  # Add the current plot to the list
  plot_list[[as.character(level)]] <- current_plot
}

# Combine all the individual plots into a single plot
(final_plot <- wrap_plots(plotlist = plot_list, ncol = 4)+  # Adjust the number of columns as needed
    plot_annotation(title="Non-metric Multidimensional Scaling Plots by sampling Year-Month",
                    subtitle = "Based on zooplankton taxon lifeforms",
                    caption = paste0("Stress = ",round(ord$stress,3),
                                     "\nSamples gathered between ",
                                     min(dfw$sample.date)," & ",
                                     max(dfw$sample.date)),
                    theme = theme(plot.title = element_text(size = 16, face="bold"))))

ggsave(filename = "figs/nmds_by_Region&YYMM_Traits.pdf",width = 12,height = 12,units = "in",
       plot=final_plot);rm(final_plot,plot_list)

#### plot by season ####
# Create an empty list to store individual plots
plot_list <- list()

# Iterate over each level of the 'yymm' factor
for (level in levels(scores_site$DJF)) {
  # Subset the data for the current 'yymm' level
  subset_data <- scores_site[scores_site$DJF == level, ]
  
  # Create a plot for the current level
  current_plot <- ggplot() +
    geom_hline(colour="grey",yintercept = 0, lty=2)+
    geom_vline(colour="grey",xintercept = 0, lty=2)+
    geom_text(data = scores_species, aes(x = NMDS1, y = NMDS2, label = lb), size = 3, alpha = 0.2) +
    # geom_segment(data = subset_data, aes(x = NMDS1, y = NMDS2, colour = Region, xend = mn_ax1_Rgn, yend = mn_ax2_Rgn), show.legend = FALSE) +
    geom_point(data = subset_data,
               aes(x = NMDS1, y = NMDS2, fill = Region, shape = Region),
               show.legend=FALSE,
               size = 3) +
    scale_fill_manual(values = consistent_palette) +
    scale_colour_manual(values = consistent_palette) +
    scale_shape_manual(values = consistent_shapes) +
    coord_equal() +
    labs(title = level)+  # Set the title to the 'DJF' level
    theme(axis.title=element_blank(),
          axis.text=element_blank())
  
  # Add the current plot to the list
  plot_list[[as.character(level)]] <- current_plot
}

# # Combine all the individual plots into a single plot
# final_plot <- wrap_plots(plotlist = plot_list, ncol = 2)+  # Adjust the number of columns as needed
#   plot_annotation(title="Non-metric Multidimensional Scaling Plots by sampling season",
#                   subtitle = "Based on zooplankton taxon lifeforms",
#                   caption = "DJF = December-February; MAM = March-May,\nJJA=June-August,SON=September-November",
#                   theme = theme(plot.title = element_text(size = 16, face="bold")))
# 
# ggsave(filename = "figs/nmds_by_Region&season_Traits.pdf",width = 12,height = 12,units = "in",
#        plot=final_plot);rm(final_plot)

# Combine all the individual plots into a single plot
(final_plot <- wrap_plots(plotlist = plot_list, ncol = 2)+  # Adjust the number of columns as needed
    plot_annotation(title="Non-metric Multidimensional Scaling Plots by sampling season",
                    subtitle = "Based on zooplankton taxon lifeforms",
                    caption = paste0("Stress = ",round(ord$stress,3),
                                     "\nSamples gathered between ",
                                     min(dfw$sample.date)," & ",
                                     max(dfw$sample.date),
                                     "\nDJF = December-February; MAM = March-May,\nJJA=June-August,SON=September-November"),
                    theme = theme(plot.title = element_text(size = 16, face="bold"))))

ggsave(filename = "figs/nmds_by_Region&season_Traits.pdf",width = 12,height = 12,units = "in",
       plot=final_plot);rm(final_plot, plot_list)
