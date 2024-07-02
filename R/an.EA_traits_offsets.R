### an.EA_traits_offsets.R ####
#### Initial explorations of data imported in imp.EA_traits.R
## adding offsets term

## set up ####
#### load packages ####
ld_pkgs <- c("tidyverse","MASS","lubridate","vegan","mvabund","seas","patchwork",
             "ecoCopula","performance","gclus","corrplot","gllvm","tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)
tictoc::tic.clearlog();tic("set universals");print("set universals")
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
toc(log=TRUE)

#### load data ####
tic("Load and format data")
dfw0 <- readRDS(paste0(datfol,"processedData/","zoopWideTraitAbund_m3_taxOnly_USE.RDat"))

##coalesce to replace NA net volume values with volume = 1
dfw <- dfw0 %>% 
  mutate(NetUse = coalesce(`Net.volume.sampled.(m3)`, 1)) %>% 
  relocate(NetUse,.after = `Net.volume.sampled.(m3)`)

df_tx_w <- dfw %>% # remove WIMS data from object
  dplyr::select(-c(WIMS.Code.y:"Zinc, Dissolved_ug/l"))
toc(log=TRUE)

# Quick ordination ####
tic("NMDS ordination prep")
### weight abundances by net volume
df_tx_w %>% 
  pivot_longer(cols = -c(Pot.Number:WB),names_to = "taxon",values_to = "abund") %>% 
  filter(abund >0) %>% 
  # mutate(NetUse = coalesce(`Net.volume.sampled.(m3)`, 1)) %>% 
  mutate(abund_m3=NetUse*abund) %>% ## weight abundances by net volume
  dplyr::select(.,-c(abund)) %>% 
  # names(.) #%>% 
  pivot_wider(.,names_from = taxon, values_from = abund_m3,values_fill = 0) %>% 
  dplyr::select(-c(Pot.Number:WB)) %>% ###remove metadata info
  dplyr::select_if(~ !is.numeric(.) || sum(.) !=0) %>% ### drop cols summing to 0
  filter(rowSums(across(where(is.numeric)))!=0) -> dftmp ###remove 'empty' rows
toc(log=TRUE)
### NMDS ####
# tic("Run NMDS")
# set.seed(pi);ord <-   vegan::metaMDS(dftmp,trymax = 500)
# saveRDS(ord, file="figs/nmds_trt.Rdat")
# toc(log=TRUE)
ord <- readRDS("figs/nmds_trt.Rdat")
plot(ord)

### extract site info
tic("Export ggplots for NMDS")
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
  
  # Force 'Region' to be a factor with specific levels
  subset_data$Region <- factor(subset_data$Region, levels = names(consistent_palette))
  
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
(final_plot <- wrap_plots(plotlist = plot_list, ncol = 6)+  # Adjust the number of columns as needed
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

rm(ord, scores_site, scores_species, subset_data,current_plot, df_tx_w,dftmp)
toc(log = TRUE)

#############
# GLLVMs ####
#############
tic("gllvm model setup")
### remove samples without net volume measurements
#dfw %>% filter(.,!is.na(`Net.volume.sampled.(m3)`)) -> dfw

### set number of iterations
runs <- 5

Y <- dfw %>% dplyr::select(.,-c(Pot.Number:WB,WIMS.Code.y:"Zinc, Dissolved_ug/l"))
X <- dfw %>% dplyr::select(.,c(Pot.Number:WB,WIMS.Code.y:"Zinc, Dissolved_ug/l"))

# choose interesting environmental variables
keep <- c(
          "Net.volume.sampled.(m3)",
          "Ammoniacal Nitrogen, Filtered as N_mg/l",
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
kp <- names(X) %in% keep
X <- X[,kp]
rm(kp, keep)

### replace LESS THAN values with value*0.5
# define function to replace "less than" values by  half their value
replace_values <- function(x) {
  if_else(str_detect(x, "^<"), as.numeric(sub("^<", "", x))/2, as.numeric(x))
}

# amend < values in wims data
X %>% 
  mutate_all(.,replace_values) -> X

##append region and wb names
X$WB <- dfw$WB
X$Region <- dfw$Region

#OPTIONAL: remove lifeforms which only appear n>=4 times ####
n <- 10 # 7 terms in 'full' X-model
Y <- Y[,colSums(ifelse(Y==0,0,1))>n]

## replace NA values with mean values for respective column.
## Leave non-numeric cols unchanged
X %>% 
  mutate(across(where(is.numeric),
                ~replace(.,
                         is.na(.),
                         mean(.,
                              na.rm = TRUE)
                )
  )
  ) -> X

## rename colums
X <- X %>% 
  rename(
    wb="WB",
    region="Region",
    net_vol_m3 = "Net.volume.sampled.(m3)",
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

X$gn <- ifelse(X$region == "Southern", "Sth",
                                 ifelse(X$region == "NEast", "NE",
                                        ifelse(X$region == "NWest", "NW",
                                               ifelse(X$region == "Anglian", "An",
                                                      ifelse(X$region == "Thames", "Th",
                                                             ifelse(X$region == "SWest", "SW",NA)
                                                      )))))

### create scaled version for comparison of effects on model
X %>% 
  mutate_if(is.numeric,scale) -> X_scaled
X_scaled$net_vol_m3 <- X$net_vol_m3
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
#### Negbin ####
# tic("Fit Unconstrained Negative binomial model")
# sDsn <- data.frame(Region = X$region)
# m_lvm_0 <- gllvm(Y, # unconstrained model
#                  offset = log(X$net_vol_m3),#set model offset
#                  studyDesign = sDsn,
#                  row.eff = ~(1|Region),
#                  family = "negative.binomial",
#                  starting.val="res",
#                  n.init=runs,#re-run model to get best fit
#                  trace=TRUE,
#                  seed = pi,
#                  num.lv = 2
#                  )
# saveRDS(m_lvm_0, file="figs/gllvm_traits_uncon_negbin.Rdat") #3.265326 mins
# toc(log=TRUE)
m_lvm_0 <- readRDS("figs/gllvm_traits_uncon_negbin.Rdat")
par(mfrow=c(2,2));plot(m_lvm_0, which=1:4);par(mfrow=c(1,1))
ordiplot(m_lvm_0,biplot = TRUE,symbols=TRUE)
cr <- getResidualCor(m_lvm_0)

pdf(file = "figs/m_lvm_0_trt_corrplot.pdf",width=14,height=14)
corrplot::corrplot(cr, diag = FALSE, type = "lower", method = "square",
                   tl.srt = 25)
dev.off();rm(cr)

##########################
# Constrained w/ Random ####
#### Nested by Region ####
##### Negbin ####
# tic("Fit Constrained Negative binomial model")
# sDsn <- data.frame(region = X$region)
# m_lvm_4 <- gllvm(y=Y, # model with environmental parameters
#                  X=X_scaled, #scaled
#                  offset = log(X$net_vol_m3),#set model offset
#                  formula = ~ nh4 + sal_ppt + chla + din + depth + po4 + tempC,
#                  studyDesign = sDsn, row.eff = ~(1|region),
#                  family="negative.binomial",
#                  starting.val="res",
#                  n.init=runs,#re-run model to get best fit
#                  trace=TRUE,
#                  num.lv = 2
#                  )
# 
# saveRDS(m_lvm_4, file="figs/gllvm_traits_nh4SalChlaDinDepPo4Reg_negbin_scaled.Rdat") #scaled #6.862054 mins
# toc(log=TRUE)

tic("Model summaries & comparisons")
m_lvm_4 <- readRDS("figs/gllvm_traits_nh4SalChlaDinDepPo4Reg_negbin_scaled.Rdat")#scaled

cr <- getResidualCor(m_lvm_4)
pdf(file = "figs/m_lvm_4_trt_corrplot.pdf",width=14,height=14)
corrplot::corrplot(cr, diag = FALSE, type = "lower", method = "square",
                   tl.srt = 25)
dev.off()

AIC(m_lvm_0,m_lvm_4)
anova(m_lvm_0,m_lvm_4)

##### GLLVM plots ####
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

# source("R/function_scores.gllvm.R")
# scrs <- scores.gllvm(m_lvm_4)

#### GLLVM model explore ####

#### extract model terms for plotting ####
mod_coefs <- as.data.frame(m_lvm_4$params$Xcoef)
mod_coefs$LF <- row.names(mod_coefs)
mod_coefs <- mod_coefs %>% 
  pivot_longer(.,cols=!LF,names_to = "coefficient", values_to = "Estimate")
sdXcoef <- as.data.frame(m_lvm_4$sd$Xcoef[,  drop = FALSE])
sdXcoef$LF <- row.names(sdXcoef)
sdXcoef <- sdXcoef %>% 
  pivot_longer(.,cols=!LF,names_to = "coefficient", values_to = "sd")

mod_coefs <- dplyr::left_join(mod_coefs,sdXcoef, by=c("LF","coefficient"))
mod_coefs <- mod_coefs %>% 
  mutate(.,lower = Estimate-1.96*sd,
         upper = Estimate+1.96*sd) %>% 
  mutate(.,varTrt=paste0(coefficient,"_",LF))

## extract P values
## extract 'significant' model/species terms from model
# ci_mod_all <- as.data.frame(confint(m_lvm_4))
# ci_mod_var <- ci_mod_all[grep("^X", rownames(ci_mod_all)), ]
# rownames(ci_mod_var) <- substring(rownames(ci_mod_var), 7)
# ci_mod_var$varTrt <- rownames(ci_mod_var)
# ci_mod_var$varTrt <- gsub("\\.","_",ci_mod_var$varTrt)

sigterms_all <- summary(m_lvm_4)
sigterms_all <- as.data.frame(sigterms_all$Coef.tableX)
sigterms_all$variable <- sub(":.*","",row.names(sigterms_all))
sigterms_all$trt <- sub(".*:","",row.names(sigterms_all))
sigterms_all$varTrt <- rownames(sigterms_all)
sigterms_all$varTrt <- gsub(":","_",sigterms_all$varTrt)
# sigterms_all <- left_join(sigterms_all, ci_mod_var, by = "varTrt")
# sigterms_all$sig <- sigterms_all$`2.5 %`*sigterms_all$`97.5 %`>0

# sigterms_sig <- sigterms_all[sigterms_all$`Pr(>|z|)`>0.05,]

mod_coefs <- dplyr::left_join(mod_coefs,
                              sigterms_all[,c("varTrt","Std. Error","z value","Pr(>|z|)")],
                              by="varTrt")
rm(sigterms_all,sdXcoef)
## create flag if CI doesn't cross 0
mod_coefs$sig <- mod_coefs$lower*mod_coefs$upper>0

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
### updated plot ####
plot_list <- list()
mod_coefs$coefficient <- as.factor(mod_coefs$coefficient)
# sigterms_all$variable <- as.factor(sigterms_all$variable)
ntrt <- length(unique(mod_coefs$LF))-.5
mod_coefs %>% 
  mutate(flag = case_when(
    !sig ~ "NONE",
    sig & Estimate >0 ~ "pos",
    sig & Estimate <0 ~ "neg"
  )) -> mod_coefs

# Iterate over each level of the factor 'LF'
for (level in levels(mod_coefs$coefficient)) {
  # Subset the data for the current level
  subset_data <- mod_coefs[mod_coefs$coefficient == level, ]
  #subset_data <- mod_coefs[mod_coefs$coefficient == "nh4", ]
  
  # Calculate the mean value of 'Estimate' for the current level
  mean_estimate <- mean(subset_data$Estimate)
  
  # Determine the color of the vertical line based on the mean estimate
  line_color <- ifelse(mean_estimate > 0, "blue", ifelse(mean_estimate < 0, "red", "grey"))
  
  # Create a plot for the current level
  current_plot <- ggplot(subset_data,
                         aes(x=Estimate, y=LF,
                             xmin=lower,
                             xmax=upper,
                             colour=flag,
                             fill=flag)) +
    geom_hline(yintercept = seq(1.5,ntrt,by=1),col="lightgrey",lty=3)+
    geom_vline(xintercept = 0)+
    # Add vertical line for mean estimate
    geom_vline(xintercept = mean_estimate, color = line_color,linetype="dashed") +
    geom_linerange()+
    labs(title = paste0(level))+
    geom_point(shape=21) +
    scale_y_discrete(limits = rev(levels(as.factor(mod_coefs$LF))))+
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

rcov0 <- getResidualCov(m_lvm_0, adjust = 1) # 'null' model
rcov1 <- getResidualCov(m_lvm_4, adjust = 1) # model with env variables #REGION
btwn <- 100 - (rcov1$trace / rcov0$trace*100)
print(paste0("Including environmental parameters in the model explains ",round(btwn,2),"% of the (co)variation in zooplankton abundances"))
AIC(m_lvm_0,m_lvm_4)

# Combine all the individual plots into a single plot
(final_plot <- wrap_plots(plotlist = plot_list,
                          ncol = nlevels(mod_coefs$coefficient))+  # Adjust the number of columns as needed
    plot_annotation(title="Caterpillar plot of generalised linear latent variable model outputs",
                    subtitle = paste0("Model including environmental variables explains ",round(btwn,2),"% of the (co)variation in lifeform abundances compared to the null (lifeforms-only) model"),
                    caption = paste0("Colours indicate lifeform 95% confidence intervals which do (grey) or do not (red/blue) include zero","\n",
                                     "Dashed vertical lines indicate mean point estimate values\n","Lifeforms recorded in fewer than ",n+1," samples removed from data prior to model estimations","\n",
                                     "Model call: ~",as.character(m_lvm_4$formula)[2],
                                     "\nDistribution family: ",as.character(m_lvm_4$family),". ",
                                     "Random row effects: ",as.character(m_lvm_4$call)[8]),
                    theme = theme(plot.title = element_text(size = 16, face="bold"))))

pdf(file = "figs/coef_trt_all_unordered_v2_scaled.pdf",width=16,height=8) #scaled
print(final_plot)
dev.off()

toc(log=TRUE)

# fun & games with unimodal QUADRATIC model ####
# model is processor-hungry and needs lots of data.  Trim Y
YQuad <- Y[,colSums(ifelse(Y==0,0,1))>100]

#### Nested by Region ####
##### Negbin ####
tic("Fit QUADRATIC unconstrained Negative binomial model")
sDsn <- data.frame(region = X$region)
m_lvm_quad_0 <- gllvm(y = YQuad, 
                      # X = X_scaled, #scaled
                      offset = log(X$net_vol_m3),#set model offset
                      # formula = ~ nh4 + sal_ppt + chla + din + depth + po4 + tempC,
                      studyDesign = sDsn, row.eff = ~(1 | region), 
                      family = "negative.binomial", 
                      starting.val = "zero", 
                      n.init = 100, n.init.max = 10,
                      trace = TRUE, 
                      num.lv = 2,
                      quadratic = TRUE,
                      seed = pi,
                      #sd.errors = FALSE#quicker fitting for model testing only
                      )
saveRDS(m_lvm_quad_0, file="figs/gllvm_traits_unconstrainedQuadraticNegBin.Rdat")
m_lvm_quad_0 <- readRDS("figs/gllvm_traits_unconstrainedQuadraticNegBin.Rdat")
toc(log=TRUE)

ordiplot(m_lvm_quad_0,biplot=TRUE,symbols = TRUE)#species optima are far away from the gradient (hence arrows)
optima(m_lvm_quad_0,sd.errors = FALSE)## inspect optima
tolerances(m_lvm_quad_0,sd.errors = FALSE)## inspect tolerances

# some very large values in the optima and tolerances, indicating linear responses
## Use model to predict some values for visualisation
#LV1
LVs = getLV(m_lvm_quad_0)
newLV <- cbind(LV1=seq(min(LVs[,1]), max(LVs[,1]),length.out=1000),LV2=0)
### causes error:
preds <- gllvm::predict.gllvm(m_lvm_quad_0,
    type = "response", newLV = newLV,
    newX = cbind(Region = rep(0, nrow(newLV))))
plot(NA,
  ylim = range(preds), xlim = c(range(getLV(m_lvm_quad_0))),
  ylab  = "Predicted response", xlab = "LV1")
segments(
  x0 = optima(m_lvm_quad_0, sd.errors = FALSE)[, 1],
  x1 = optima(m_lvm_quad_0, sd.errors = FALSE)[, 1],
  y0 = rep(0, ncol(m_lvm_quad_0$y)),
  y1 = apply(preds, 2, max),
  col = "red", lty = "dashed", lwd = 2)
rug(getLV(m_lvm_quad_0))[,1]
sapply(1:ncol(m_lvm_quad_0$y), function(j)
  lines(sort(newLV[, 1]), preds[order(newLV[, 1]), j], lwd = 2))

## LV2
newLV = cbind(LV1 = 0, LV2 =  seq(min(LVs[, 2]), max(LVs[, 2]),
                                  length.out = 1000))
preds <- gllvm::predict.gllvm(m_lvm_quad_0,
    type = "response", newLV = newLV,
    newX = cbind(Region = rep(0, nrow(newLV))))
plot( NA,
  ylim = range(preds), xlim = c(range(getLV(m_lvm_quad_0))),
  ylab  = "Predicted response", xlab = "LV2")
segments(
  x0 = optima(m_lvm_quad_0, sd.errors = FALSE)[, 2],
  x1 = optima(m_lvm_quad_0, sd.errors = FALSE)[, 2],
  y0 = rep(0, ncol(m_lvm_quad_0$y)),
  y1 = apply(preds, 2, max),
  col = "red", lty = "dashed", lwd = 2)
rug(getLV(m_lvm_quad_0)[, 2])
sapply(1:ncol(m_lvm_quad_0$y), function(j)
  lines(sort(newLV[, 2]), preds[order(newLV[, 2]), j], lwd = 2))

# calculate turnover for these two estimated gradients
# Extract tolerances
tol <- tolerances(m_lvm_quad_0, sd.errors = FALSE)
gradLength <- 4/apply(tol, 2, median)
cat("Gradient length:", gradLength)
turn <- 2*qnorm(.999, sd = apply(tol, 2, median))
cat("Turnover rate:", turn)

paste0("LV1 tunrover rate = ",round(2*qnorm(.999,sd=1/sqrt(-2*m_lvm_quad_0$params$theta[1,3]*m_lvm_quad_0$params$sigma.lv[1]^2)),2))
paste0("LV2 tunrover rate = ",round(2*qnorm(.999,sd=1/sqrt(-2*m_lvm_quad_0$params$theta[1,4]*m_lvm_quad_0$params$sigma.lv[2]^2)),2))

# TIDY UP ####
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
   ntrt,sbtt,srt,ttl,cr,n, rcov0,rcov1)
rm(list = ls(pattern = "^X"))
rm(Y,runs,mean_estimate,line_color,btwn)
rm(LVs, newLV,preds,tol,YQuad,gradLength,turn)

########## TO DO ##############
# * fix caterpillar plot: currently giving off-centre point estimates
# * check data: unconstrained ordination only showing variation on one axis
# * 
