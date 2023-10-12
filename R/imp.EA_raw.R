### imp.EA_raw.R ####
#### Import EA data, attach WIMS data & create a block for later analysis ####
#### analysed Jul 2023
### RH has amended the files returned from the lab to include
### sample volumes

### useful links for plotting MDS with ggplot
# https://chrischizinski.github.io/rstats/vegan-ggplot2/
# https://www.youtube.com/watch?v=Y0GI34S-ZMI

## set up ####
#### install required packages ####
#### set local package library 

#### check and install required packages ####
ptm <- Sys.time()
req_packages <- c("lubridate",# working with dates
                  "patchwork", #arranging multiple plots
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
                  "seas"#converting dates to seasonal factors
)

new_packages <- req_packages[!(req_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages,library=libfolder,type="binary")
Sys.time() - ptm;rm(ptm,req_packages,new_packages)

#### load packages ####
ld_pkgs <- c("tidyverse","MASS","lubridate","vegan","mvabund","seas","patchwork",
             "ecoCopula","performance","gclus","corrplot")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

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

### load data ####
### taxon data
df_tx0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,
                                               "processedData/MBA_Returns_Amalgamated_USE.xlsx"),
                                        # sheet="outR01"))
                                        # sheet="outR02"))
                                        # sheet="outR03"))
                                        sheet="outR04"))
### WIMS chemical data
df_wims0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,
                                                 "/WIMS_Extract_WaterQuality_Zoop_Samples_230809.xlsx"),
                                          sheet="allDat"))

### counts for pots 1&2 have been summed as to have those for pots 3&4
### removed whitespace from AphiaIDs & taxon names

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

## write WIMS wide csv
# write.csv(df_wims_w,file=paste0(datfol,"processedData/WIMS_wide.csv"),row.names = FALSE)

### remove odd data
df_tx <- as_tibble(df_tx0) ### create new data (keep df0 as 'raw')

### convert dates
df_tx$sample.date <- as.Date(df_tx$sample.date, origin = "1899-12-30")

### sites now matched to WBs in Excel doc

# Remove 100 µm data [OPTIONAL] ####
df_tx_100um <- df_tx %>% 
  filter(str_starts(Sample.comments,"100um"))

df_tx %>% 
  filter(!str_starts(Sample.comments,"100um")) -> df_tx

###widen data & fill NAs with 0s ####
df_tx %>% 
  dplyr::select(.,-c("Aphia.ID","AbundanceRaw")) %>% 
  pivot_wider(names_from = "Taxa",values_from = "Abund_m3",
              values_fill = 0) -> df_tx_w

### join & save data ####  
dfw <- left_join(df_tx_w,df_wims_w,by="PRN")

# write.csv(dfw,file=paste0(datfol,"processedData/","zoopWIDEAbund_m3_WIMS_USE.csv"),row.names = FALSE)
# write.csv(df_tx_w,file=paste0(datfol,"processedData/","zoopWIDEAbund_m3_taxOnly_USE.csv"),row.names = FALSE)

###############################
## look at taxon data only ####
###############################

### we now have formatted abundance data standardised by volume of seawater filtered ###

### Initial analyses ####

### tax number
df_tx_w$S <- vegan::specnumber(df_tx_w[, -c(1:21)])
df_tx_w$N <- rowSums(df_tx_w[, -c(1:21, length(df_tx_w))])

### how often do different taxa appear?
xx <- df_tx_w %>% 
  dplyr::select(-c(1:21)) %>% ###remove metadata info
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
  dplyr::select(-c(WB,WB_lb)) %>% ## plot by Region
  # dplyr::select(-c(WB,Region)) %>% ## plot by WB
  group_by(Region) %>% ## plot by Region
  # group_by(WB_lb) %>%  ## plot by WB
  summarise(across(everything(),sum)) %>% # -> df_tx_prev_Rgn
  ### need to pivot to longer for plotting
  pivot_longer(!Region, names_to="taxon",values_to="preval") %>% #by Region
  # pivot_longer(!WB_lb, names_to="taxon",values_to="preval") %>%
  dplyr::filter(preval >0) -> xxprev ## generates prevalence (# observations) per taxon

xxprev %>% 
  ggplot(.,aes(x=preval))+
  geom_histogram(fill = "lightgrey", colour = 1,bins=40)+
  facet_wrap(.~Region)+ #plot by Region
  # facet_wrap(.~WB_lb)+ #plot by WB
  labs(title="Prevalence of zooplankton taxa recorded by Region", #by Region
       # labs(title="Prevalence of zooplankton taxa recorded by Region_Water body", # by WB
       subtitle ="'Rare' taxa are common. Indicates a high degree of heterogeneity in the data", #For Region
       caption=paste0("Prevalence is the number of times that individual taxa have been recorded in a given region (at any abundance)\nZero values excluded\n",#Region
                      # caption=paste0("Prevalence is the number of times that individual taxa have been recorded in a given water body (at any abundance)\nZero values excluded\n",#WB
                      "Samples gathered between ",min(dfw$sample.date)," & ",max(dfw$sample.date)),
       y = "Count",
       x = "Prevalence")
ggsave(filename = "figs/zoopPrevalence_by_Region.pdf",width = 12,height = 7,units = "in")
# ggsave(filename = "figs/zoopPrevalence_by_WB.pdf",width = 12,height = 10,units = "in")

### quick species accumulation curve ####
df_tx_w %>% 
  dplyr::select(-c(1:21)) %>% ###remove metadata info
  dplyr::select(-last_col()) %>% ###remove taxon abundance
  dplyr::select(-last_col()) %>% ###remove taxon richness
  vegan::specaccum(.) -> sp1
plot(sp1, ci.type = "poly",col="blue",lwd=2,ci.lty = 0,ci.col = "lightblue")
# str(sp1)

df_tx_w %>% 
  dplyr::select(-c(1:21)) %>% ###remove metadata info
  dplyr::select(-last_col()) %>% ###remove taxon abundance
  dplyr::select(-last_col()) %>% ###remove taxon richness
  dplyr::select_if(~ !is.numeric(.) || sum(.) !=0) %>% ### drop cols summing to 0
  #vegan::specpool(.) ### estimates of taxon richness
  # vegan::specpool(.,dfw$Region) ### estimates of taxon richness by region
  vegan::poolaccum(.) -> spaccum; plot(spaccum)
# vegan::specnumber(.)

# pdf(file = "figs/spp_accm_all_WBs_200um_Jun22_Apr23.pdf",width=12,height=7)
# plot(spaccum)
# dev.off()

#### quick ordinations ####
df_tx_w %>% 
  dplyr::select(-c(1:21)) %>% ###remove metadata info
  dplyr::select(-last_col()) %>% ###remove taxon abundance
  dplyr::select(-last_col()) %>% ###remove taxon richness
  dplyr::select_if(~ !is.numeric(.) || sum(.) !=0) %>% ### drop cols summing to 0
  filter(rowSums(across(where(is.numeric)))!=0) -> dftmp ###remove 'empty' rows

### NMDS ####
ptm <- Sys.time()###
set.seed(pi+5);ord <-   vegan::metaMDS(dftmp,trymax = 1000)
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
scores_species$lb <-  make.cepnames(row.names(scores_species))#shorten names

# scores_site$loc <- substr(scores_site$`station.name/WIMS.code`,1,3)

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
  labs(title="Non-metric Multidimensional Scaling of zooplankton taxon abundances",
       subtitle="Colours & shapes indicate region",
       caption=paste0("Stress = ",round(ord$stress,3),"\nSamples gathered between ",
                      min(dfw$sample.date)," & ",max(dfw$sample.date)))+
  theme(legend.title = element_blank(),
        axis.title = element_text(face="bold"));pl

ggsave(filename = "figs/nmds_by_Region.pdf",width = 12,height = 12,units = "in",
       plot=pl);rm(pl)

# pl <- ggplot(data=scores_site,aes(x=NMDS1,y=NMDS2,
#                                   fill=Region,
#                                   colour=Region))+
#   geom_hline(colour="grey",yintercept = 0, lty=2)+
#   geom_vline(colour="grey",xintercept = 0, lty=2)+
#   geom_text(data=scores_species, aes(x = NMDS1, y=NMDS2, label=lb),
#             size=3,alpha=0.2,inherit.aes = FALSE)+
#   # stat_ellipse(geom="polygon",data=scores_site,aes(x=NMDS1,y=NMDS2,
#   #                                                  colour=Region, fill=Region),
#   #              level = .8,alpha=0.5)+
#   geom_segment(data=scores_site,aes(x=NMDS1,y=NMDS2,
#                                     colour=Region,
#                                     xend=mn_ax1_Rgn,yend=mn_ax2_Rgn),
#                show.legend = FALSE)+
#   geom_point(size=3)+
#   scale_shape_manual(values = c(24:21))+
#   coord_equal()+
#   guides(colour=guide_legend(override.aes = list(shape=22)))+
#   scale_fill_manual(values=c(cbPalette))+
#   scale_colour_manual(values=c(cbPalette2))+
#   labs(caption=paste0("Stress = ",
#                       round(ord$stress,3)))+
#   theme(legend.title = element_blank(),
#         axis.title = element_text(face="bold"));pl
# 
# ggsave(filename = "figs/nmds_by_Region.pdf",width = 12,height = 12,units = "in",
#        plot=pl);rm(pl)
# rm(pl)

# pl <- ggplot(data=scores_site,aes(x=NMDS1,y=NMDS2,
#                                   fill=Region,
#                                   colour=Region))+
#   geom_hline(colour="grey",yintercept = 0, lty=2)+
#   geom_vline(colour="grey",xintercept = 0, lty=2)+
#   geom_text(data=scores_species, aes(x = NMDS1, y=NMDS2, label=lb),
#             size=3,alpha=0.2,inherit.aes = FALSE)+
#   geom_segment(data=scores_site,aes(x=NMDS1,y=NMDS2,
#                                     colour=Region,
#                                     xend=mn_ax1_Rgn,yend=mn_ax2_Rgn),
#                show.legend = FALSE)+
#   geom_point(size=3)+
#   scale_shape_manual(values = c(24:21))+
#   coord_equal()+
#   guides(colour=guide_legend(override.aes = list(shape=22)))+
#   scale_fill_manual(values=c(cbPalette))+
#   scale_colour_manual(values=c(cbPalette2))+
#   labs(caption=paste0("Stress = ",round(ord$stress,3)))+
#   theme(legend.title = element_blank(),
#         axis.title = element_text(face="bold"));pl
# 
# ggsave(filename = "figs/nmdsJan_23_02.pdf",width = 12,height = 12,units = "in",
#        plot=pl);rm(pl)
# rm(pl)

### variability of ntax
##add 'mesh' to data
dfw %>% 
  dplyr::mutate(mesh=ifelse(grepl("200um",Sample.comments),"200um",
                            ifelse(grepl("100um",Sample.comments),"100um",NA))) -> dfw
# dplyr::mutate(mesh=substr(Sample.comments,1,5)) -> dfstd

### add 'label' variable for Water Bodies
df_tx_w$WB_lb1 <- ifelse(df_tx_w$Region == "Southern","Sth",
                         ifelse(df_tx_w$Region == "Thames","Thm",
                                ifelse(df_tx_w$Region == "Anglian","Ang",
                                       ifelse(df_tx_w$Region == "NWest","NW",
                                              ifelse(df_tx_w$Region == "NEast","NE",
                                                     ifelse(df_tx_w$Region == "SWest","SW",NA)
                                              )))))

df_tx_w$WB_lb2 <- ifelse(df_tx_w$WB == "Solent","Solent",
                         ifelse(df_tx_w$WB == "SOUTHAMPTON WATER","Soton Wtr",
                                ifelse(df_tx_w$WB == "THAMES LOWER","Thm Low",
                                       ifelse(df_tx_w$WB == "Blackwater Outer","Blckw Out",
                                              ifelse(df_tx_w$WB == "Cornwall North","Cornw Nth",
                                                     ifelse(df_tx_w$WB == "Barnstaple Bay","Brnstp B",
                                                            ifelse(df_tx_w$WB == "Kent South","Kent Sth",
                                                                   ifelse(df_tx_w$WB == "Mersey Mouth","Mersey Mth",
                                                                          ifelse(df_tx_w$WB == "Wash Outer","Wash Out",
                                                                                 ifelse(df_tx_w$WB == "Lincolnshire","Lincs",
                                                                                        ifelse(df_tx_w$WB == "Yorkshire South","Yorks Sth",
                                                                                               ifelse(df_tx_w$WB == "TEES","Tees",
                                                                                                      ifelse(df_tx_w$WB == "Northumberland North","Nrthmb Nth",
                                                                                                             ifelse(df_tx_w$WB == "Farne Islands to Newton Haven","Farne Is",
                                                                                                                    ifelse(df_tx_w$WB == "Bristol Channel Inner South","Brist Ch In Sth",
                                                                                                                           NA))))))))))
                                              )))))
df_tx_w$WB_lb <- paste0(df_tx_w$WB_lb1,"_",df_tx_w$WB_lb2)
df_tx_w$WB_lb1 <- NULL; df_tx_w$WB_lb2 <- NULL

meanS <- mean(df_tx_w$S); sdS <- sd(df_tx_w$S)
sdSmin <- meanS-sdS; sdSmax <- meanS+sdS

pl_S <- df_tx_w[df_tx_w$S>0,] %>% 
  # ggplot(.,aes(x=Region,y=S))+
  ggplot(.,aes(x=WB_lb,y=S))+
  geom_hline(yintercept = 0,col="lightgrey",lty=2)+
  geom_ribbon(aes(ymin=sdSmin,ymax=sdSmax), group=1, fill="lightgrey",alpha=50)+
  geom_hline(yintercept = meanS, col = "black", lty = 3)+
  geom_boxplot(aes(fill=Region), outlier.shape = NA,varwidth = TRUE, show.legend = FALSE)+
  # geom_point(aes(fill=WB),position = position_jitterdodge(), show.legend = FALSE)+
  geom_jitter(width=0.2)+
  scale_fill_manual(values=cbPalette) +
  labs(title="Zooplankton taxon richness by water body",
       subtitle ="Colours indicate water body Region",
       caption=paste0("Dotted line & shading indicate global mean ±  standard deviation\n",
                      "Samples gathered between ",min(df_tx_w$sample.date),
                      " & ", max(df_tx_w$sample.date))) +
  xlab("Region_Water Body")+
  ylab("Taxon richness")+
  theme(axis.text.x = element_text(angle=270,hjust=0,vjust=.5));pl_S

ggsave(filename = "figs/taxRich.pdf",width = 12,height = 12,units = "in",
       plot=pl_S)
rm(pl_S,meanS,sdS,sdSmax,sdSmin)

pl_ts_S <- df_tx_w[df_tx_w$S>0,] %>% 
  ggplot(.,aes(#group = WB_lb,
    y= S, x=sample.date)) +
  # geom_line()+
  geom_point(data = transform(df_tx_w[df_tx_w$S>0,],
                              WB_lb=NULL),# create 'grey' points
             col="lightgrey") +
  geom_smooth(se=FALSE, span=0.9, method="loess") + ## Loess smoother
  # geom_smooth(se=FALSE, method = "gam",formula = y ~ s(x, k = 3)) + ## GAM
  geom_point(pch=21,aes(fill=Region), size=3, show.legend=FALSE) +
  scale_fill_manual(values=cbPalette) +
  coord_cartesian(ylim=c(0,NA)) +
  labs(title="Temporal trends in taxon richness of zooplankton assemblages by water body",
       caption = "Grey points display ALL data\nColours indicate Region\nBlue line indicates loess smooth (span=0.9)",
       x="Date")+
  facet_wrap(.~WB_lb);pl_ts_S

ggsave(filename = "figs/taxRichTSbyWB.pdf",width = 12,height = 12,units = "in",
       plot=pl_ts_S)
rm(pl_ts_S)

meanN <- mean(df_tx_w$N); sdN <- sd(df_tx_w$N)
logmeanN <- mean(log(df_tx_w$N)); logsdN <- sd(log(df_tx_w$N))
logsdNmin <- logmeanN-logsdN; logsdNmax <- logmeanN+logsdN

### variability of abund
pl_N <- df_tx_w[df_tx_w$N>0,] %>% 
  ggplot(.,aes(x=WB_lb,y=log(N)))+
  geom_hline(yintercept = 0,col="lightgrey",lty=2)+
  geom_ribbon(aes(ymin=logsdNmin,ymax=logsdNmax), group=1, fill="lightgrey",alpha=50)+
  geom_hline(yintercept = logmeanN, col = "black", lty = 3)+
  geom_boxplot(aes(fill=Region),outlier.shape = NA,varwidth = TRUE, show.legend=FALSE)+
  geom_jitter(width=0.2)+
  # geom_point(aes(fill=WB),position = position_jitterdodge(), show.legend = FALSE)+
  # labs(title="Total abundance (individuals m^3)")+
  scale_fill_manual(values=cbPalette)+
  labs(title = bquote('Total zooplankton abundance'~(log~individuals~m^-3)),
       subtitle ="Colours indicate region",
       caption=paste0("Dotted line & shading indicate global mean ±  standard deviation\n",
                      "Samples gathered between ",min(df_tx_w$sample.date),
                      " & ", max(df_tx_w$sample.date))) +
  xlab("Region_Water Body")+
  ylab("log(abundance)")+
  theme(axis.text.x = element_text(angle=270,hjust=0,vjust=.5)); pl_N

ggsave(filename = "figs/taxAbund.pdf",width = 12,height = 12,units = "in",
       plot=pl_N)

rm(pl_N, meanN, logmeanN, logsdNmin, sdN, meanN)

pl_ts_N <- df_tx_w[df_tx_w$N>0,] %>% 
  ggplot(.,aes(#group = WB_lb,
    y= log(N), x=sample.date))+
  geom_point(data = transform(df_tx_w[df_tx_w$N>0,],
                              WB_lb=NULL),
             col="lightgrey")+
  geom_smooth(se=FALSE, span=0.9)+
  geom_point(pch=21,aes(fill=Region), size=3, show.legend=FALSE) +
  scale_fill_manual(values=cbPalette) +
  coord_cartesian(ylim=c(0,NA)) +
  labs(title="Log(taxon abundance) of zooplankters over time by water body",
       caption = "Grey points display ALL data\nColours indicate Region\nBlue line indicates loess smooth (span=0.9)",
       x="Date")+
  facet_wrap(.~WB_lb);pl_ts_N

ggsave(filename = "figs/taxAbundTSbyWB.pdf",width = 12,height = 12,units = "in",
       plot=pl_ts_N)

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
    # geom_segment(data = subset_data, aes(x = NMDS1, y = NMDS2, colour = Region, xend = mn_ax1_Rgn, yend = mn_ax2_Rgn), show.legend = FALSE) +
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
                    subtitle = "Based on zooplankton taxon abundance data",
                    caption = paste0("Stress = ",round(ord$stress,3),
                                     "\nSamples gathered between ",
                                     min(dfw$sample.date)," & ",
                                     max(dfw$sample.date)),
                    theme = theme(plot.title = element_text(size = 16, face="bold"))))

ggsave(filename = "figs/nmds_by_Region&YYMM_Taxa.pdf",width = 12,height = 12,units = "in",
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
    labs(title = level)+  # Set the title to the 'yymm' level
    theme(axis.title=element_blank(),
          axis.text=element_blank())
  
  # Add the current plot to the list
  plot_list[[as.character(level)]] <- current_plot
}

# Combine all the individual plots into a single plot
(final_plot <- wrap_plots(plotlist = plot_list, ncol = 2)+  # Adjust the number of columns as needed
    plot_annotation(title="Non-metric Multidimensional Scaling Plots by sampling season",
                    subtitle = "Based on zooplankton taxon abundance data",
                    caption = paste0("Stress = ",round(ord$stress,3),
                                     "\nSamples gathered between ",
                                     min(dfw$sample.date)," & ",
                                     max(dfw$sample.date),
                                     "\nDJF = December-February; MAM = March-May,\nJJA=June-August,SON=September-November"),
                    theme = theme(plot.title = element_text(size = 16, face="bold"))))

ggsave(filename = "figs/nmds_by_Region&season_Taxa.pdf",width = 12,height = 12,units = "in",
       plot=final_plot);rm(final_plot, plot_list)

########## STATS ###############

### MVABUND ####
#### unconstrained ordination ####
mv_dftmp <- mvabund::mvabund(dftmp)#create mvabund object
ttl <- "Very strong mean-variance relationship in zooplankton abundance data"
sbtt <- "Variance within the dataset covers *12 orders of magnitude*. Many multivariate analyses (e.g. ANOSIM, PERMANOVA) assume *no mean-variance relationship*\nThis makes interpretation of such analyses potentially erroneous. Model-based approaches offer an alternative, allowing the mean-variance relationship to be incorporated into the model predictions"

png(file = "figs/zoopMeanVar.png",
    width=12*ppi, height=6*ppi, res=ppi)
mvabund::meanvar.plot(mv_dftmp,# mean-variance plot
                      # main = "Mean-variance relationship of zooplankton taxon abundances",
                      # # sub="Many multivariate analyses assume no mean-variance relationship",
                      add.trendline=TRUE,
                      xlab="Mean",
                      ylab="Variance")

mtext(side=3, line = 3, at =-0.07, adj=0, cex = 1, ttl, font=2)
mtext(side=3, line = 0.75, at =-0.07, adj=0, cex = 0.7, sbtt)
dev.off()

## poisson: Intercept only (unconstrained):
system.time(mod1 <- manyglm(mv_dftmp~1,family="poisson"))
# summary(mod1)
plot(mod1)
system.time(ord_glmInt1 <- cord(mod1))
plot(ord_glmInt1, biplot=TRUE)
srt <- order.single(ord_glmInt1$sigma)

pdf(file = "figs/Unconstr_poisson_corrplot_Jun22_Jul23.pdf",width=12,height=12)
# corrplot(ord_glmInt1$sigma[srt,srt],type = "lower",diag = FALSE,method="square",
corrplot(ord_glmInt1$sigma,type = "lower",diag = FALSE,method="square",
         tl.col = 1,tl.cex=0.35)
dev.off()

## negative binomial
system.time(mod2 <- manyglm(mv_dftmp~1,family="negative_binomial"))
# summary(mod2)
plot(mod2)
system.time(ord_glmInt2 <- cord(mod2))
plot(ord_glmInt2, biplot=TRUE)
srt <- order.single(ord_glmInt2$sigma)
# pl_cor <- corrplot(ord_glmInt2$sigma[srt,srt],type = "lower",diag = FALSE,method="square")

pdf(file = "figs/Unconstr_negbin_corrplot_Jun22_Jul23.pdf",width=12,height=12)
# corrplot(ord_glmInt2$sigma[srt,srt],type = "lower",diag = FALSE,method="square",
corrplot(ord_glmInt2$sigma,type = "lower",diag = FALSE,method="square",
         tl.col = 1,tl.cex=0.35)
dev.off()
###

library(mgcv)

tw <- manyany(mv_dftmp~1,"glm",family=Tweedie(var.power=1.2), var.power=1.2)

#######################
#######################
#######################

### initial statistical comparisons ####
### taxon richness ####

### subset standardised data
# dfsub <- dfw[dfw$S>0,] ### removes 'empty' samples

#### check distribution ####
plot(performance::check_distribution(df_tx_w$S)) ### >60% beta-binomial

### taxon richness ####
summary(m_S_0 <- lm(S~WB_lb, data=df_tx_w));visreg::visreg(m_S_0)
performance::check_posterior_predictions(m_S_0)

summary(m_S_pois <- glm(S~WB_lb,
                        family = poisson(),
                        data=df_tx_w)
);visreg::visreg(m_S_pois)
performance::check_posterior_predictions(m_S_pois)

summary(m_S_nb <- MASS::glm.nb(S~WB_lb,
                               data=df_tx_w));visreg::visreg(m_S_nb)
performance::check_posterior_predictions(m_S_nb)

AIC(m_S_0,m_S_pois,m_S_nb) ## NBGLM > LM > Poisson

#################################################################
performance::check_model(m_S_0)
performance::check_model(m_S_nb)## 'best' model
performance::check_model(m_S_pois)

### taxon abundance ####
#### check distribution ####
# plot(performance::check_distribution(dfsub$N)) ### 
plot(performance::check_distribution(log(df_tx_w$N))) # log ~Weibull distribution 
plot(performance::check_distribution(df_tx_w$N)) # lognormal

### run models ####
summary(m_N_0 <- lm(N~WB_lb, data=df_tx_w))
visreg::visreg(m_N_0)
summary(m_N_pois <- glm(N~WB_lb,data=df_tx_w, family = poisson()))
summary(m_N_nb <- MASS::glm.nb(N~WB_lb,data=df_tx_w))
summary(m_N_lgnm <- glm(log(N)~WB_lb, family=gaussian(link="log"),data=df_tx_w))
visreg::visreg(m_N_lgnm)

performance::check_model(m_N_lgnm)
AIC(m_N_0,m_N_pois,m_N_nb,m_N_lgnm)#log-normal is 'best'

### multivariate analysis ####
(df_adonis <- vegan::adonis2(dftmp ~ scores_site$WB, by="onedf"))
#plot(df_adonis)
# dfdist <- vegdist(dftmp)
# summary(df_anosim <- anosim(dfdist,scores_site$WB))
# plot(df_anosim)

### to do:
### look at functional groups(lifeforms)
# master list saved in:
# \\prodds.ntnl\Shared\AN\KFH\Groups\N_Marine\07 Training & Reference Documents\A&R Technical Guidance\Traits, Lifeforms etc\Plankton Lifeform Extraction Tool
# consider reproducing ordination in 3 dimensions using rgl (see: https://riffomonas.org/code_club/2021-03-24-rgl)

