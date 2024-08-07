### an.EA_traits.R ####
#### Initial explorations of data imported in imp.EA_traits.R

## set up ####
#### load packages ####
ld_pkgs <- c("tidyverse","MASS","lubridate","vegan","mvabund","seas","patchwork",
             "ecoCopula","performance","gclus","corrplot","gllvm","tictoc", "GGally")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

tictoc::tic.clearlog()
tic("set universals");print("set universals")
#### set universals ####
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
toc(log = TRUE)

tic("Load & format data")
#### load data ####
dfw <- readRDS(paste0(datfol,"processedData/","zoopWideTraitAbund_m3_taxOnly_USE.RDat"))

df_tx_w <- dfw %>% # remove WIMS data from object
  dplyr::select(-c(WIMS.Code.y:"Zinc, Dissolved_ug/l"))

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
                           ifelse(xx$WB == "Solway Outer South","Solway O",
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
                                                                                                                             NA)))))))))))
                                         )))))
xx$WB_lb <- paste0(xx$WB_lb1,"_",xx$WB_lb2)
xx$WB_lb1 <- NULL; xx$WB_lb2 <- NULL

toc(log = TRUE)

tic("Initial plots")
### plot by WB####
xx %>% 
  dplyr::select(-c(WB,Region)) %>% ## plot by WB
  group_by(WB_lb) %>%  ## plot by WB
  summarise(across(everything(),sum)) %>%
  ### need to pivot to longer for plotting
  pivot_longer(!WB_lb, names_to="taxon",values_to="preval") %>% #by WB
  dplyr::filter(preval >0) -> xxprev

xxprev %>% 
  ggplot(.,aes(x=preval))+
  geom_histogram(fill = "lightgrey", colour = 1,bins=40)+
  facet_wrap(.~WB_lb)+ #plot by WB
  labs(title="Prevalence of zooplankton lifeforms recorded by Region_Water body", # by WB
       caption=paste0("Prevalence is the number of times that individual lifeforms have been recorded in a given water body (at any abundance)\nZero values excluded\n",#WB
                      "Samples gathered between ",min(dfw$sample.date)," & ",max(dfw$sample.date)),
       y = "Count",
       x = "Prevalence")

ggsave(filename = "figs/zoopTraitsPrevalence_by_WB.pdf",width = 12,height = 10,units = "in")

### plot by Region ####
xx %>% 
  dplyr::select(-c(WB,WB_lb)) %>% ## plot by Region
  group_by(Region) %>% ## plot by Region
  summarise(across(everything(),sum)) %>% 
  ### need to pivot to longer for plotting
  pivot_longer(!Region, names_to="taxon",values_to="preval") %>% #by Region
  dplyr::filter(preval >0) -> xxprev

xxprev %>% 
  ggplot(.,aes(x=preval))+
  geom_histogram(fill = "lightgrey", colour = 1,bins=40)+
  facet_wrap(.~Region)+ #plot by Region
  labs(title="Prevalence of zooplankton lifeforms recorded by Region", #by Region
       caption=paste0("Prevalence is the number of times that individual lifeforms have been recorded in a given region (at any abundance)\nZero values excluded\n",#Region
                      "Samples gathered between ",min(dfw$sample.date)," & ",max(dfw$sample.date)),
       y = "Count",
       x = "Prevalence")

ggsave(filename = "figs/zoopTraitsPrevalence_by_Region.pdf",width = 12,height = 7,units = "in")
rm(xx,xxprev)

# Quick ordination ####
df_tx_w %>% 
  dplyr::select(-c(1:20)) %>% ###remove metadata info
  dplyr::select(-last_col()) %>% ###remove taxon abundance
  dplyr::select(-last_col()) %>% ###remove taxon richness
  dplyr::select_if(~ !is.numeric(.) || sum(.) !=0) %>% ### drop cols summing to 0
  filter(rowSums(across(where(is.numeric)))!=0) -> dftmp ###remove 'empty' rows

### NMDS ####
# ptm <- Sys.time()###
# set.seed(pi);ord <-   vegan::metaMDS(dftmp,trymax = 500)
# saveRDS(ord, file="figs/nmds_trt.Rdat")
# Sys.time() - ptm;rm(ptm)
ord <- readRDS("figs/nmds_trt.Rdat")
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

rm(ord, scores_site, scores_species, subset_data,current_plot)
toc(log=TRUE)

#############
# GLLVMs ####
#############
tic("Set up gllvm model")
# choose interesting environmental variables
keep <- c("Ammoniacal Nitrogen, Filtered as N_mg/l",
          "Chlorophyll : Acetone Extract_ug/l",
          "NGR : Easting_NGR",
          "NGR : Northing_NGR",
          "Nitrate, Filtered as N_mg/l",
          "Nitrite, Filtered as N_mg/l",
          "Nitrogen, Dissolved Inorganic : as N_mg/l",
          "Nitrogen, Total Oxidised, Filtered as N_mg/l",
          "Orthophosphate, Filtered as P_mg/l",
          "Oxygen, Dissolved as O2_mg/l",
          "Oxygen, Dissolved, % Saturation_%",
          "Salinity : In Situ_ppt",
          "Silicate, Filtered as SiO2_mg/l",
          "Temperature of Water_CEL",
          "Turbidity : In Situ_FTU",
          "Water Depth_m")

# keep only interesting variables
kp <- names(dfw) %in% keep
df_wims_w_trim <- dfw[,kp]
rm(kp, keep)

### replace LESS THAN values with value*0.5
# define function to replace "less than" values by  half their value
replace_values <- function(x) {
  if_else(str_detect(x, "^<"), as.numeric(sub("^<", "", x))/2, as.numeric(x))
}

# amend < values in wims data
df_wims_w_trim %>% 
  mutate_all(.,replace_values) -> df_wims_w_trim

### append region & WB
df_wims_w_trim$Region <- dfw$Region
df_wims_w_trim$WB <- dfw$WB

# extract taxon density data
dfw %>% 
  dplyr::select(-c(Pot.Number:WB,
                   "WIMS.Code.y":"Zinc, Dissolved_ug/l")
  ) -> df_tx_w_trm

#OPTIONAL: remove lifeforms which only appear n>=4 times ####
n <- 4
df_tx_w_trm <- df_tx_w_trm %>% 
  dplyr::select(which(apply(., 2, function(x) sum(x != 0) > n)))

## replace NA values with mean values for respective column.
## Leave non-numeric cols unchanged
df_wims_w_trim %>% 
  mutate(across(where(is.numeric),
                ~replace(.,
                         is.na(.),
                         mean(.,
                              na.rm = TRUE)
                )
  )
  ) -> df_wims_w_trim0

## rename colums
df_wims_w_trim0 <- df_wims_w_trim0 %>% 
  rename(
    nh4 = "Ammoniacal Nitrogen, Filtered as N_mg/l",
    chla = "Chlorophyll : Acetone Extract_ug/l",
    ngr_e = "NGR : Easting_NGR",
    ngr_n = "NGR : Northing_NGR",
    no3 = "Nitrate, Filtered as N_mg/l",
    no2 = "Nitrite, Filtered as N_mg/l",
    din = "Nitrogen, Dissolved Inorganic : as N_mg/l",
    ton = "Nitrogen, Total Oxidised, Filtered as N_mg/l",
    po4 = "Orthophosphate, Filtered as P_mg/l",
    o2_dis_mgl = "Oxygen, Dissolved as O2_mg/l",
    o2_dis_sat = "Oxygen, Dissolved, % Saturation_%",
    sal_ppt = "Salinity : In Situ_ppt",
    si = "Silicate, Filtered as SiO2_mg/l",
    tempC = "Temperature of Water_CEL",
    turb = "Turbidity : In Situ_FTU",
    depth = "Water Depth_m"
  )

df_wims_w_trim0$Region <- ifelse(df_wims_w_trim0$Region == "Southern", "Sth",
                                 ifelse(df_wims_w_trim0$Region == "NEast", "NE",
                                        ifelse(df_wims_w_trim0$Region == "NWest", "NW",
                                               ifelse(df_wims_w_trim0$Region == "Anglian", "An",
                                                      ifelse(df_wims_w_trim0$Region == "Thames", "Th",
                                                             ifelse(df_wims_w_trim0$Region == "SWest", "SW",NA)
                                                      )))))

### create scaled version for comparison of effects on model
df_wims_w_trim0 %>% 
  mutate_if(is.numeric,scale) -> df_wims_w_trim0_scale
toc(log=TRUE)

### fit models ####

# distribution families (from ?gllvm & https://cran.r-project.org/web/packages/gllvm/vignettes/vignette1.html)
# COUNTS:
# "negative.binomial" (with log link), poisson(link = "log"),
# zero-inflated poisson ("ZIP"), zero-inflated negative-binomial ("ZINB"),
# BINARY: 
# binomial(link = "probit") (and also with link = "logit" when method = "LA" or method = "EVA"),
# NORMAL:
# gaussian(link = "identity"),
# BIOMASS:
# Tweedie ("tweedie") (with log link, for "LA" and "EVA"-method),
# PERCENTAGE COVER:
# beta ("beta") (with logit and probit link, for "LA" and "EVA"-method)
# ORDINAL DATA:
# "ordinal" (only with "VA"-method).
# POSITIVE CONTINUOUS:
# "gamma" (with log link),
# NON-NEGATIVE CONTINUOUS:
# "exponential" (with log link),

### Fit models ####
# Unconstrained w/ Random ####
#### Tweedie ####
# tic("Unconstrained Tweedie model, nested by Region"")
# sDsn <- data.frame(Region = df_wims_w_trim0$Region)
# m_lvm_0 <- gllvm(df_tx_w_trm, # unconstrained model
#                  studyDesign = sDsn, row.eff = ~(1|Region),
#                  family = "tweedie"
#                  )
# saveRDS(m_lvm_0, file="figs/gllvm_traits_uncon_tweed.Rdat") #3.265326 mins
# toc(log=TRUE)
m_lvm_0 <- readRDS("figs/gllvm_traits_uncon_tweed.Rdat")

###################
#### Gaussian ####
# ptm <- Sys.time()
# sDsn <- data.frame(Region = df_wims_w_trim0$Region)
# m_lvm_0 <- gllvm(df_tx_w_trm, # unconstrained model
#                  studyDesign = sDsn, row.eff = ~(1|Region),
#                  family = gaussian()
#                  )
# saveRDS(m_lvm_0, file="figs/gllvm_traits_uncon_gauss.Rdat")
# Sys.time() - ptm;rm(ptm)
# m_lvm_0 <- readRDS("figs/gllvm_traits_uncon_gauss.Rdat")

##########################
# Constrained w/ Random ####
#### Nested by Region ####
##### Tweedie ####
# tic("Constrained Tweedie, nested by Region")
# sDsn <- data.frame(Region = df_wims_w_trim0$Region)
# m_lvm_4 <- gllvm(y=df_tx_w_trm, # model with environmental parameters
#                  # X=df_wims_w_trim0, #unscaled
#                  X=df_wims_w_trim0_scale, #scaled
#                  formula = ~ nh4 + sal_ppt + chla + din + depth + po4 + tempC,
#                  studyDesign = sDsn, row.eff = ~(1|Region),
#                  family="tweedie"
#                  )
# # # saveRDS(m_lvm_4, file="figs/gllvm_traits_nh4SalChlaDinDepPo4Reg_tweed.Rdat") #unscaled #11.623mins
# saveRDS(m_lvm_4, file="figs/gllvm_traits_nh4SalChlaDinDepPo4Reg_tweed_scaled.Rdat") #scaled #6.862054 mins
# Sys.time() - ptm;rm(ptm)
# m_lvm_4 <- readRDS("figs/gllvm_traits_nh4SalChlaDinDepPo4Reg_tweed.Rdat")#unscaled
m_lvm_4 <- readRDS("figs/gllvm_traits_nh4SalChlaDinDepPo4Reg_tweed_scaled.Rdat")#scaled

######
##### Gaussian ####
# ptm <- Sys.time()
# sDsn <- data.frame(Region = df_wims_w_trim0$Region)
# m_lvm_4 <- gllvm(y=df_tx_w_trm, # model with environmental parameters
#                  X=df_wims_w_trim0,
#                  formula = ~ nh4 + sal_ppt + chla + din + depth + po4 + tempC,
#                  studyDesign = sDsn, row.eff = ~(1|Region),
#                  family=gaussian()
#                  )
# saveRDS(m_lvm_4, file="figs/gllvm_traits_nh4SalChlaDinDepPo4Reg_gauss.Rdat")#11.623mins
# Sys.time() - ptm;rm(ptm)
# m_lvm_4 <- readRDS("figs/gllvm_traits_nh4SalChlaDinDepPo4Reg_gauss.Rdat")

tic("Model plots")
cr <- getResidualCor(m_lvm_4)
pdf(file = "figs/m_lvm_4_trt_corrplot.pdf",width=14,height=14)
corrplot::corrplot(cr, diag = FALSE, type = "lower", method = "square",
                   tl.srt = 25)
dev.off()

AIC(m_lvm_0,m_lvm_4)
anova(m_lvm_0,m_lvm_4)

# # Constrained + Random
# sDsn <- data.frame(Region = df_wims_w_trim0$Region)
# ptm <- Sys.time()
# m_lvm_4 <- gllvm(y=df_tx_w_trm, # model with environmental parameters & random FX
#                  X=df_wims_w_trim0,
#                  formula = ~ nh4 + sal_ppt + chla + din + depth + po4 + tempC,
#                  # family = gaussian(),
#                  family = "tweedie",
#                  studyDesign = sDsn, row.eff = ~(1|Region)
# )
# # saveRDS(m_lvm_4, file="figs/gllvm_traits_nh4SalChlaDinDepPo4Tmp_norm.Rdat")#4.1096 mins
# saveRDS(m_lvm_4, file="figs/gllvm_traits_nh4SalChlaDinDepPo4Tmp_tweed.Rdat")#4.1096 mins
# Sys.time() - ptm;rm(ptm) 
# m_lvm_4 <- readRDS("figs/gllvm_traits_nh4SalChlaDinDepPo4Tmp_tweed.Rdat")
# 
# AIC(m_lvm_0,m_lvm_3,m_lvm_4)
# anova(m_lvm_0,m_lvm_3,m_lvm_4)

#### Nested by WB|Region ####
##### Tweedie ####
# ptm <- Sys.time()
# sDsn2 <- data.frame(WB = df_wims_w_trim0$WB)
# m_lvm_5 <- gllvm(y=df_tx_w_trm, # model with environmental parameters
#                  # X=df_wims_w_trim0, #unscaled
#                  X=df_wims_w_trim0_scale, #scaled
#                  formula = ~ nh4 + sal_ppt + chla + din + depth + po4 + tempC,
#                  studyDesign = sDsn2, row.eff = ~(1|WB),
#                  family="tweedie"
#                  )
# # saveRDS(m_lvm_5, file="figs/gllvm_traits_nh4SalChlaDinDepPo4Reg_tweed.Rdat") #unscaled #11.623mins
# saveRDS(m_lvm_5, file="figs/gllvm_traits_nh4SalChlaDinDepPo4WB_tweed_scaled.Rdat") #scaled #6.862054 mins
# Sys.time() - ptm;rm(ptm)
# m_lvm_5 <- readRDS("figs/gllvm_traits_nh4SalChlaDinDepPo4WB_tweed_scaled.Rdat")#scaled

##### GLLVM plots ####
# pdf(file = "figs/m_lvm_3_trt_all_ordered.pdf",width=16,height=8)
# coefplot(m_lvm_3,cex.ylab = 0.3,
#          order=TRUE)
# dev.off()
# 
pdf(file = "figs/coef_trt_all_unordered.pdf",width=16,height=8)
coefplot(m_lvm_4,cex.ylab = 0.7,
         order=FALSE)
dev.off()

pdf(file = "figs/coef_trt_1.pdf",width=7,height=14)
coefplot(m_lvm_4,mfrow = c(1,1),which.Xcoef = 1, cex.ylab = 0.6,
         main="NH4")
dev.off()

pdf(file = "figs/coef_trt_2.pdf",width=7,height=14)
coefplot(m_lvm_4,mfrow = c(1,1),which.Xcoef = 2, cex.ylab = 0.6,
         main="Salinity")
dev.off()

pdf(file = "figs/coef_trt_3.pdf",width=7,height=14)
coefplot(m_lvm_4,mfrow = c(1,1),which.Xcoef = 3, cex.ylab = 0.6,
         main="Chlorophyll")
dev.off()

pdf(file = "figs/coef_trt_4.pdf",width=7,height=14)
coefplot(m_lvm_4,mfrow = c(1,1),which.Xcoef = 4, cex.ylab = 0.6,
         main="DIN")
dev.off()

pdf(file = "figs/coef_trt_5.pdf",width=7,height=14)
coefplot(m_lvm_4,mfrow = c(1,1),which.Xcoef = 5, cex.ylab = 0.6,
         main="Depth")
dev.off()

pdf(file = "figs/coef_trt_6.pdf",width=7,height=14)
coefplot(m_lvm_4,mfrow = c(1,1),which.Xcoef = 6, cex.ylab = 0.6,
         order=TRUE, main="PO4")
dev.off()

pdf(file = "figs/coef_trt_7.pdf",width=7,height=14)
coefplot(m_lvm_4,mfrow = c(1,1),which.Xcoef = 7, cex.ylab = 0.6,
         order=TRUE, main="Temperature")
dev.off()

#### GLLVM model explore ####
# ordiplot.gllvm(m_lvm_0)
# # ordiplot.gllvm(m_lvm_3)
# ordiplot.gllvm(m_lvm_4, biplot = TRUE)

# tail(confint.gllvm(m_lvm_3))
# tail(confint.gllvm(m_lvm_4))

## extract 'significant' model/species terms from model
ci_mod_all <- as.data.frame(confint(m_lvm_4))
ci_mod_var <- ci_mod_all[grep("^X", rownames(ci_mod_all)), ]
rownames(ci_mod_var) <- substring(rownames(ci_mod_var), 7)
ci_mod_var$varTrt <- rownames(ci_mod_var)

sigterms_all <- summary(m_lvm_4)
sigterms_all <- as.data.frame(sigterms_all$Coef.tableX)
sigterms_all$variable <- sub(":.*","",row.names(sigterms_all))
sigterms_all$trt <- sub(".*:","",row.names(sigterms_all))
sigterms_all$varTrt <- rownames(sigterms_all)
sigterms_all <- left_join(sigterms_all, ci_mod_var, by = "varTrt")
sigterms_all$sig <- sigterms_all$`2.5 %`*sigterms_all$`97.5 %`>0

sigterms_sig <- sigterms_all[sigterms_all$`Pr(>|z|)`>0.05,]

# View(sigterms_sig)
# View(sigterms_all)

### plot! ####
## recreate coefplot
# ggplot(sigterms_all[sigterms_all$variable=="nh4",],
#        aes(x=Estimate, y=trt,
#            xmin=`2.5 %`,
#            xmax=`97.5 %`,
#            colour=sig))+
#   geom_vline(xintercept = 0)+
#   geom_errorbar()+
#   geom_point()+
#   scale_y_discrete(limits = rev(levels(as.factor(sigterms_all$trt))))+
#   scale_colour_manual(values = c("grey","black"))+
#   guides(colour="none")

#############
plot_list <- list()
sigterms_all$variable <- as.factor(sigterms_all$variable)
ntrt <- length(unique(sigterms_all$trt))-.5
sigterms_all %>% 
  mutate(flag = case_when(
    !sig ~ "NONE",
    sig & Estimate >0 ~ "pos",
    sig & Estimate <0 ~ "neg"
  )) -> sigterms_all

# Iterate over each level of the factor 'trt'
for (level in levels(sigterms_all$variable)) {
  # Subset the data for the current level
  subset_data <- sigterms_all[sigterms_all$variable == level, ]
  
  # Calculate the mean value of 'Estimate' for the current level
  mean_estimate <- mean(subset_data$Estimate)
  
  # Determine the color of the vertical line based on the mean estimate
  line_color <- ifelse(mean_estimate > 0, "blue", ifelse(mean_estimate < 0, "red", "grey"))
  
  # Create a plot for the current level
  current_plot <- ggplot(subset_data,
                         aes(x=Estimate, y=trt,
                             xmin=`2.5 %`,
                             xmax=`97.5 %`,
                             colour=flag,
                             fill=flag)) +
    geom_hline(yintercept = seq(1.5,ntrt,by=1),col="lightgrey",lty=3)+
    geom_vline(xintercept = 0)+
    # Add vertical line for mean estimate
    geom_vline(xintercept = mean_estimate, color = line_color,linetype="dashed") +
    geom_linerange()+
    labs(title = paste0(level))+
    geom_point(shape=21) +
    scale_y_discrete(limits = rev(levels(as.factor(sigterms_all$trt))))+
    scale_colour_manual(values = c("red",#negative
                                         "grey",#null
                                         "blue"#positive
                                         ))+
    scale_fill_manual(values = c("red",#negative
                                        "white",#null
                                        "blue"#positive
                                        ))+
    guides(colour="none",
           fill="none")+
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust=0.5))
  
  # Add the current plot to the list
  plot_list[[as.character(level)]] <- current_plot
  
  # Iterate over each plot in the list
  for (i in seq_along(plot_list)) {
    # If it's not the first plot, hide y-axis labels
    if (i > 1) {
      plot_list[[i]] <- plot_list[[i]] + theme(axis.text.y = element_blank())
    }
  }
}

# Combine all the individual plots into a single plot
(final_plot <- wrap_plots(plotlist = plot_list,
                          ncol = nlevels(sigterms_all$variable))+  # Adjust the number of columns as needed
    plot_annotation(title="Caterpillar plot of generalised linear latent variable model outputs",
                    # subtitle = "Based on zooplankton taxon lifeforms and water quality parameters",#unscaled
                    subtitle = bquote("Point estimates & 95% confidence intervals of lifeform-specific model coefficients "~italic(hat(beta)[j])~". Based on zooplankton taxon lifeform densities and scaled water quality parameters"), #scaled
                    caption = paste0("Colours indicate lifeform 95% confidence intervals which do (grey) or do not (red/blue) include zero","\n",
                                     "Dashed vertical lines indicate mean point estimate values\n","Lifeforms recorded in fewer than ",n+1," samples removed from data prior to model estimations","\n",
                                     "Model call: ~",as.character(m_lvm_4$formula)[2],
                                     "\nDistribution family: ",as.character(m_lvm_4$family),". ",
                                     "Random row effects: ",as.character(m_lvm_4$call)[7]),
                    theme = theme(plot.title = element_text(size = 16, face="bold"))))

# pdf(file = "figs/coef_trt_all_unordered_v2_unscaled.pdf",width=16,height=8)#unscaled
pdf(file = "figs/coef_trt_all_unordered_v2_scaled.pdf",width=16,height=8) #scaled
print(final_plot)
dev.off()

### compare models w/ w/out env parameters####
# based on:
#https://jenniniku.github.io/gllvm/articles/vignette1.html#studying-co-occurrence-patterns
# getResidualCov function can be used to quantify the amount of variation in the
# data that can be explained by environmental variables.
# Specifically, if we use the trace of the residual covariance matrix Σ as a
# measure of unexplained variation, then we can compare this quantity before
# and after environmental variables are included in the model. The ratio of
# traces suggests that environmental variables explain approximately 40% of
# the (co)variation in ant species abundances.

(rcov0 <- getResidualCov(m_lvm_0, adjust = 0)) # 'null' model
(rcov1 <- getResidualCov(m_lvm_4, adjust = 0)) # model with env variables #REGION
# (rcov1 <- getResidualCov(m_lvm_5, adjust = 0)) # model with env variables #WB
rcov0$trace; rcov1$trace
100 - (rcov1$trace / rcov0$trace*100)
AIC(m_lvm_0,m_lvm_4)

toc(log=TRUE)

#######
tic("generate Pairs plot")
### generate Pairs plot
n <- 200 #number of non zero values

# Function to return points and geom_smooth
# allow for the method to be changed
# my_fn <- function(data, mapping, method="loess", ...){
#   p <- ggplot(data = data, mapping = mapping) + 
#     geom_point() + 
#     geom_smooth(method=method, ...)
#   p
# }
mth <- "gam"
# Custom function to modify the lower panels
my_fn <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) + 
    geom_point(color = "black", alpha = 0.3, size = 1) +
    geom_smooth(method = mth, color = "blue")
}

dfw %>% 
  # remove non-zoop columns
  dplyr::select(.,-c(Pot.Number:WB,
                     WIMS.Code.y:"Zinc, Dissolved_ug/l")) %>% 
  #remove taxa with <n occurences
  dplyr::select(where(~ sum(. != 0) >= n)) %>% 
  # GGally::ggpairs(.,lower = list(continuous = wrap(my_fn, method="gam"))) +
  # GGally::ggpairs(.,lower = list(continuous = wrap("smooth",
  #                                                  method="gam",
  #                                                  colour="blue",
  #                                                  alpha=0.3))) +
  GGally::ggpairs(.,lower = list(continuous = wrap(my_fn))) +
  labs(title = "Pairwise plot of zooplankton lifeform abundances",
       subtitle = paste0("Only lifeforms recorded in ",n,
                         " or more samples included (representing ",
                         round(n/nrow(dfw)*100,1),"% of samples). Blue line represents ",mth," smoother.")) -> pl
  
pdf(file = "figs/pairs_plot_LF.pdf",width=24,height=15) #scaled
print(pl)
dev.off()
rm(mth, my_fn,pl)
toc(log=TRUE)

log_out <- unlist(tictoc::tic.log())

# PRIORITY : REPRODUCE CODE ####
## Currently untidy and seems to produce 'issues'
## Errors alluding to "non-numeric argument to binary operator"

# TIDY UP ####
tic("TIDY UP")
# remove data
rm(list = ls(pattern = "^df"))
rm(list = ls(pattern = "^m_"))
rm(list = ls(pattern = "^mod"))
rm(list = ls(pattern = "^cbPalette"))
rm(list = ls(pattern = "^scores"))
rm(list = ls(pattern = "^sigt"))
rm(list = ls(pattern = "^sp"))
rm(list = ls(pattern = "^logs"))
rm(list = ls(pattern = "^ci_"))
rm(list = ls(pattern = "^max_"))
rm(list = ls(pattern = "^min_"))
rm(list = ls(pattern = "*plot"))
rm(list = ls(pattern = "*consiste"))
rm(datfol,nit,perms, ppi,replace_values,pl_ts_N,sDsn,subset_data,i,level,
   ntrt,sbtt,srt,ttl,cr,n, rcov0,rcov1,mean_estimate,line_color)

# unload packages
detach("package:seas", unload=TRUE)
detach("package:patchwork", unload=TRUE)
detach("package:performance", unload=TRUE)
detach("package:corrplot", unload=TRUE)
detach("package:gllvm", unload=TRUE)
detach("package:ecoCopula", unload=TRUE)
detach("package:gclus", unload=TRUE)
detach("package:vegan", unload=TRUE)
detach("package:mvabund", unload=TRUE)
detach("package:lubridate", unload=TRUE)
detach("package:tidyverse", unload=TRUE)
detach("package:MASS", unload=TRUE)
toc(log=TRUE)
detach("package:tictoc", unload=TRUE)
log_out
rm(log_out)
