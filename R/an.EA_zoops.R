### an.EA_zoops.R ####
#### Inital explorations of data imported in imp.EA_raw.R

## set up ####
#### load packages ####
ld_pkgs <- c("tidyverse","MASS","lubridate","vegan","mvabund","seas","patchwork",
             "ecoCopula","performance","gclus","corrplot","gllvm")
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

#### load data ####
dfw <- readRDS(paste0(datfol,"processedData/","zoopWIDEAbund_m3_WIMS_USE.RDat"))

df_tx_w <- dfw %>% 
  dplyr::select(-c(WIMS.Code.y:"Zinc, Dissolved_ug/l"))

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
                    ifelse(xx$WB == "Solway Outer South","Solway O",
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
                                                                                                                             NA)))))))))))
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
                      "Samples gathered between ",format(min(dfw$sample.date), "%d/%m/%Y")," & ",format(max(dfw$sample.date), "%d/%m/%Y")),
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

# pdf(file = "figs/spp_accm_all_WBs_200um.pdf",width=12,height=7)
# plot(spaccum)
# dev.off()

### quick ordinations ####
df_tx_w %>% 
  dplyr::select(-c(1:21)) %>% ###remove metadata info
  dplyr::select(-last_col()) %>% ###remove taxon abundance
  dplyr::select(-last_col()) %>% ###remove taxon richness
  #dplyr::select(where(~sum(. !=0)>1)) %>%  ### [OPTIONAL!!!] remove singleton taxa
  dplyr::select_if(~ !is.numeric(.) || sum(.) !=0) %>% ### drop cols summing to 0
  filter(rowSums(across(where(is.numeric)))!=0) -> dftmp ###remove 'empty' rows

### NMDS ####
# ptm <- Sys.time()###
# set.seed(pi+15);ord <-   vegan::metaMDS(dftmp,trymax = 500)
# ord <- vegan::metaMDS(dftmp,trymax = 500, previous.best = ord)
# saveRDS(ord, file="figs/nmds_taxa.Rdat")
# Sys.time() - ptm;rm(ptm)
ord <- readRDS("figs/nmds_taxa.Rdat")
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
                      format(min(dfw$sample.date),"%d/%m/%Y")," & ",format(max(dfw$sample.date),"%d/%m/%Y")))+
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

# Variability Plots ####
### Prep data ####
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
                         ifelse(df_tx_w$WB == "Solway Outer South","Solway O",
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
                                                                                                                                  NA)))))))))))
                                              )))))

df_tx_w$WB_lb <- paste0(df_tx_w$WB_lb1,"_",df_tx_w$WB_lb2)
df_tx_w$WB_lb1 <- NULL; df_tx_w$WB_lb2 <- NULL

### order factors
df_tx_w$WB_lb <- factor(df_tx_w$WB_lb,
                        levels = c(
                          "NE_Nrthmb Nth", "NE_Farne Is", "NE_Tees",
                          "Ang_Yorks Sth", "Ang_Lincs", "Ang_Wash Out",
                          "Ang_Blckw Out", "Thm_Thm Low", "Sth_Kent Sth",
                          "Sth_Solent", "Sth_Soton Wtr", "SW_Cornw Nth",
                          "SW_Brnstp B", "SW_Brist Ch In Sth", "NW_Mersey Mth",
                          "NW_Solway O"
                        ))
#### Ntax (S) ####
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
                      "Samples gathered between ",format(min(df_tx_w$sample.date),"%d/%m/%Y"),
                      " & ", format(max(df_tx_w$sample.date),"%d/%m/%Y"))) +
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
                      "Samples gathered between ",format(min(df_tx_w$sample.date),"%d/%m/%Y"),
                      " & ", format(max(df_tx_w$sample.date),"%d/%m/%Y"))) +
  xlab("Region_Water Body")+
  ylab("log(abundance)")+
  theme(axis.text.x = element_text(angle=270,hjust=0,vjust=.5)); pl_N

ggsave(filename = "figs/taxAbund.pdf",width = 12,height = 12,units = "in",
       plot=pl_N)

rm(pl_N, meanN, logmeanN, logsdNmin, sdN)

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
                                     format(min(dfw$sample.date),"%d/%m/%Y")," & ",
                                     format(max(dfw$sample.date),"%d/%m/%Y")),
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
                                     format(min(dfw$sample.date),"%d/%m/%Y")," & ",
                                     format(max(dfw$sample.date),"%d/%m/%Y"),
                                     "\nDJF = December-February; MAM = March-May,\nJJA=June-August,SON=September-November"),
                    theme = theme(plot.title = element_text(size = 16, face="bold"))))

ggsave(filename = "figs/nmds_by_Region&season_Taxa.pdf",width = 12,height = 12,units = "in",
       plot=final_plot);rm(final_plot, plot_list)

### STATS ###############

### MVABUND ####
#### unconstrained ordination ####
mv_dftmp <- mvabund::mvabund(dftmp)#create mvabund object

png(file = "figs/zoopMeanVar.png",
    width=12*ppi, height=6*ppi, res=ppi)
mvpl <- mvabund::meanvar.plot(mv_dftmp,# mean-variance plot
                      # main = "Mean-variance relationship of zooplankton taxon abundances",
                      # # sub="Many multivariate analyses assume no mean-variance relationship",
                      add.trendline=TRUE,
                      xlab="Mean",
                      ylab="Variance",
                      table=TRUE
                      )

# Find the minimum and maximum values
min_value <- min(mvpl[,2])
max_value <- max(mvpl[,2])

min_order <- floor(log10(min_value))
max_order <- floor(log10(max_value))
orders_of_magnitude_covered <- max_order - min_order

ttl <- "Very strong mean-variance relationship in zooplankton abundance data"
sbtt <- paste0("Variance within the dataset covers ",orders_of_magnitude_covered," orders of magnitude*.\nMany multivariate analyses (e.g. ANOSIM, PERMANOVA) assume *no mean-variance relationship*\nThis makes interpretation of such analyses potentially erroneous. Model-based approaches offer an alternative, allowing the mean-variance relationship to be incorporated into the model predictions")

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

pdf(file = "figs/Unconstr_poisson_corrplot_Jun22_Nov23.pdf",width=12,height=12)
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

pdf(file = "figs/Unconstr_negbin_corrplot_Jun22_Nov23.pdf",width=12,height=12)
# corrplot(ord_glmInt2$sigma[srt,srt],type = "lower",diag = FALSE,method="square",
corrplot(ord_glmInt2$sigma,type = "lower",diag = FALSE,method="square",
         tl.col = 1,tl.cex=0.35)
dev.off()
###

#######################

### initial statistical comparisons ####
### taxon richness ####

### subset standardised data
# dfsub <- dfw[dfw$S>0,] ### removes 'empty' samples

#### check distribution ####
plot(performance::check_distribution(df_tx_w$S)) ### beta-binomial

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

AIC(m_S_0,m_S_pois,m_S_nb) ## LM > NBGLM > Poisson

#################################################################
performance::check_model(m_S_0)## 'best' model
performance::check_model(m_S_nb)
performance::check_model(m_S_pois)

### taxon abundance ####
#### check distribution ####
# plot(performance::check_distribution(dfsub$N)) ### 
plot(performance::check_distribution(log(df_tx_w$N))) # tweedie distribution 
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
# (df_adonis <- vegan::adonis2(dftmp ~ scores_site$WB, by="onedf"))
#plot(df_adonis)
# dfdist <- vegdist(dftmp)
# summary(df_anosim <- anosim(dfdist,scores_site$WB))
# plot(df_anosim)

#################################################################
# GLLVMs ####
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

# df_tx_w %>% 
#   dplyr::select(-c(Pot.Number:Category,S,N,WB_lb)) -> df_tx_w_trm

# extract taxon density data
dfw %>% 
  dplyr::select(-c(Pot.Number:Category,
                   "WIMS.Code.y":"Zinc, Dissolved_ug/l","mesh")
  ) -> df_tx_w_trm

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

## rename columns
df_wims_w_trim0 <- df_wims_w_trim0 %>% 
  rename(
    nh4="Ammoniacal Nitrogen, Filtered as N_mg/l",
    chla ="Chlorophyll : Acetone Extract_ug/l",
    ngr_e="NGR : Easting_NGR",
    ngr_n="NGR : Northing_NGR",
    no3="Nitrate, Filtered as N_mg/l",
    no2="Nitrite, Filtered as N_mg/l",
    din="Nitrogen, Dissolved Inorganic : as N_mg/l",
    ton="Nitrogen, Total Oxidised, Filtered as N_mg/l",
    po4="Orthophosphate, Filtered as P_mg/l",
    o2_dis_mgl="Oxygen, Dissolved as O2_mg/l",
    o2_dis_sat="Oxygen, Dissolved, % Saturation_%",
    sal_ppt="Salinity : In Situ_ppt",
    si="Silicate, Filtered as SiO2_mg/l",
    tempC="Temperature of Water_CEL",
    turb="Turbidity : In Situ_FTU",
    depth="Water Depth_m"
  )

### fit models ####
### unconstrained model ####
#### Tweedie distribution ####
# ptm <- Sys.time()
sDsn <- data.frame(Region = df_wims_w_trim0$Region)
# m_lvm_0 <- gllvm(df_tx_w_trm,
#                  family="tweedie",
#                  studyDesign = sDsn, row.eff = ~(1|Region)
#                  )
# saveRDS(m_lvm_0, file="figs/gllvm_uncon_tweed.Rdat")
# Sys.time() - ptm;rm(ptm)
m_lvm_0 <- readRDS("figs/gllvm_uncon_tweed.Rdat")

#### gaussian distribution ####
# ptm <- Sys.time()
# sDsn <- data.frame(Region = df_wims_w_trim0$Region)
# m_lvm_0 <- gllvm(df_tx_w_trm,
#                  family=gaussian(),
#                  studyDesign = sDsn, row.eff = ~(1|Region),
#                  starting.val="random"
#                  )
# saveRDS(m_lvm_0, file="figs/gllvm_uncon_gauss.Rdat")
# Sys.time() - ptm;rm(ptm)
# m_lvm_0 <- readRDS("figs/gllvm_uncon_gauss.Rdat")

#### gamma distribution ####
# ptm <- Sys.time()
# sDsn <- data.frame(Region = df_wims_w_trim0$Region)
# m_lvm_0 <- gllvm(df_tx_w_trm,
#                  family="gamma",
#                  studyDesign = sDsn, row.eff = ~(1|Region)
# )
# saveRDS(m_lvm_0, file="figs/gllvm_uncon_gamma.Rdat")
# Sys.time() - ptm;rm(ptm)
# m_lvm_0 <- readRDS("figs/gllvm_uncon_gamma.Rdat")

### constrained 1 ####
# ptm <- Sys.time()
# m_lvm_1 <- gllvm(df_tx_w_trm, # model with environmental parameters
#                  df_wims_w_trim0,
#                  formula = ~ nh4 + sal_ppt + chla + din + depth + po4 + tempC + Region,
#                  # family="negative.binomial"
#                  family="tweedie"
#                  )
# saveRDS(m_lvm_1, file="figs/gllvm_env_tweed.Rdat")
# Sys.time() - ptm;rm(ptm)
m_lvm_1 <- readRDS("figs/gllvm_env_tweed.Rdat")

### constrained 2 ####
# ptm <- Sys.time()
# m_lvm_2 <- gllvm(df_tx_w_trm, # model with environmental parameters
#                  df_wims_w_trim0,
#                  formula = ~ nh4 + sal_ppt + chla + din + depth + po4 + tempC + Region,
#                  # family="negative.binomial"
#                  family="tweedie"
#                  )
# saveRDS(m_lvm_2, file="figs/gllvm_envReg_tweed.Rdat")
# Sys.time() - ptm;rm(ptm)
m_lvm_2 <- readRDS("figs/gllvm_envReg_tweed.Rdat")

# pdf(file = "figs/coef_1.pdf",width=7,height=14)
# coefplot(m_lvm_2,mfrow = c(1,1),which.Xcoef = 1, cex.ylab = 0.3,
#          main="NH4")
# dev.off()
# 
# pdf(file = "figs/coef_2.pdf",width=7,height=14)
# coefplot(m_lvm_2,mfrow = c(1,1),which.Xcoef = 2, cex.ylab = 0.3,
#          main="Salinity")
# dev.off()
# 
# pdf(file = "figs/coef_3.pdf",width=7,height=14)
# coefplot(m_lvm_2,mfrow = c(1,1),which.Xcoef = 3, cex.ylab = 0.3,
#          main="Chlorophyll")
# dev.off()
# 
# pdf(file = "figs/coef_4.pdf",width=7,height=14)
# coefplot(m_lvm_2,mfrow = c(1,2),which.Xcoef = 4:5, cex.ylab = 0.3,
#          order=FALSE)
# dev.off()
# 
# pdf(file = "figs/coef_5.pdf",width=7,height=14)
# coefplot(m_lvm_2,mfrow = c(1,2),which.Xcoef = 6:7, cex.ylab = 0.3,
#          order=FALSE)
# dev.off()
# 
# pdf(file = "figs/coef_6.pdf",width=7,height=14)
# coefplot(m_lvm_2,mfrow = c(1,2),which.Xcoef = 8, cex.ylab = 0.3,
#          order=FALSE)
# dev.off()

### constrained 3 ####
#### Tweedie #####
# ptm <- Sys.time()
# sDsn <- data.frame(Region = df_wims_w_trim0$Region)
# m_lvm_3 <- gllvm(y=df_tx_w_trm, # model with environmental parameters
#                  X=df_wims_w_trim0,
#                  formula = ~ nh4 + sal_ppt + chla + din + depth + po4 + tempC,
#                  family="tweedie",
#                  studyDesign = sDsn, row.eff = ~(1|Region)
# )
# saveRDS(m_lvm_3, file="figs/gllvm_nh4SalChlaDinDepPo4Reg_tweed.Rdat")
# Sys.time() - ptm;rm(ptm) #37.7332 mins
m_lvm_3 <- readRDS("figs/gllvm_nh4SalChlaDinDepPo4Reg_tweed.Rdat")

#### Gaussian #####
# ptm <- Sys.time()
# sDsn <- data.frame(Region = df_wims_w_trim0$Region)
# m_lvm_3 <- gllvm(y=df_tx_w_trm, # model with environmental parameters
#                  X=df_wims_w_trim0,
#                  formula = ~ nh4 + sal_ppt + chla + din + depth + po4 + tempC,
#                  family=gaussian(),
#                  studyDesign = sDsn, row.eff = ~(1|Region),
#                  starting.val="random"
#                  )
# saveRDS(m_lvm_3, file="figs/gllvm_nh4SalChlaDinDepPo4Reg_gauss.Rdat")
# Sys.time() - ptm;rm(ptm) #37.7332 mins
# m_lvm_3 <- readRDS("figs/gllvm_nh4SalChlaDinDepPo4Reg_gauss.Rdat")

#### Gamma #####
# FAILS!!!! ###
# ptm <- Sys.time()
# sDsn <- data.frame(Region = df_wims_w_trim0$Region)
# m_lvm_3 <- gllvm(y=df_tx_w_trm, # model with environmental parameters
#                  X=df_wims_w_trim0,
#                  formula = ~ nh4 + sal_ppt + chla + din + depth + po4 + tempC,
#                  family="gamma",
#                  studyDesign = sDsn, row.eff = ~(1|Region)
# )
# saveRDS(m_lvm_3, file="figs/gllvm_nh4SalChlaDinDepPo4Reg_gamma.Rdat")
# Sys.time() - ptm;rm(ptm) #37.7332 mins
# m_lvm_3 <- readRDS("figs/gllvm_nh4SalChlaDinDepPo4Reg_gamma.Rdat")
##########

pdf(file = "figs/m_lvm_3_tx_all.pdf",width=16,height=8)
coefplot(m_lvm_3,cex.ylab = 0.3,
         order=FALSE)
dev.off()

### compare models
AIC(m_lvm_0,m_lvm_1,m_lvm_2,m_lvm_3)
anova(m_lvm_0,m_lvm_1,m_lvm_2,m_lvm_3)

## extract 'significant' model/species terms from model
ci_mod_all <- as.data.frame(confint(m_lvm_3))
ci_mod_var <- ci_mod_all[grep("^X", rownames(ci_mod_all)), ]
rownames(ci_mod_var) <- substring(rownames(ci_mod_var), 7)
ci_mod_var$varTrt <- rownames(ci_mod_var)

sigterms_all <- summary(m_lvm_3)
sigterms_all <- as.data.frame(sigterms_all$Coef.tableX)
sigterms_all$variable <- sub(":.*","",row.names(sigterms_all))
sigterms_all$trt <- sub(".*:","",row.names(sigterms_all))
sigterms_all$varTrt <- rownames(sigterms_all)
sigterms_all <- left_join(sigterms_all, ci_mod_var, by = "varTrt")
sigterms_all$sig <- sigterms_all$`2.5 %`*sigterms_all$`97.5 %`>0

sigterms_sig <- sigterms_all[sigterms_all$`Pr(>|z|)`>0.05,]

### plot! ####
## recreate coefplot
ggplot(sigterms_all[sigterms_all$variable=="nh4",],
       aes(x=Estimate, y=trt,
           xmin=`2.5 %`,
           xmax=`97.5 %`,
           colour=sig))+
  geom_vline(xintercept = 0)+
  geom_errorbar()+
  geom_point()+
  scale_y_discrete(limits = rev(levels(as.factor(sigterms_all$trt))))+
  scale_colour_manual(values = c("grey","black"))+
  guides(colour="none")

#############
plot_list <- list()
sigterms_all$variable <- as.factor(sigterms_all$variable)
ntrt <- length(unique(sigterms_all$trt))-.5

# Iterate over each level of the factor 'trt'
for (level in levels(sigterms_all$variable)) {
  # Subset the data for the current level
  subset_data <- sigterms_all[sigterms_all$variable == level, ]
  
  # Create a plot for the current level
  current_plot <- ggplot(subset_data,
                         aes(x=Estimate, y=trt,
                             xmin=`2.5 %`,
                             xmax=`97.5 %`,
                             colour=sig,
                             fill=sig)) +
    geom_hline(yintercept = seq(1.5,ntrt,by=1),col="lightgrey",lty=3)+
    geom_vline(xintercept = 0)+
    # geom_errorbar()+
    geom_linerange()+
    labs(title = paste0(level))+
    geom_point(shape=21) +
    scale_y_discrete(limits = rev(levels(as.factor(sigterms_all$trt))))+
    scale_colour_manual(values = c("grey","black"))+
    scale_fill_manual(values = c("white","black"))+
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
    if (i == 1){
      plot_list[[i]] <- plot_list[[i]] + theme(axis.text.y = element_text(size=2.5))
    }
  }
}

# Combine all the individual plots into a single plot
final_plot <- wrap_plots(plotlist = plot_list, ncol = nlevels(sigterms_all$variable))+  # Adjust the number of columns as needed
    plot_annotation(title="Generalised linear latent variable model outputs",
                    subtitle = "Based on zooplankton taxon abundance data",
                    caption = paste0("Colours indicate lifeform 95% confidence intervals which do (grey) or do not (black) include zero",
                                     "\nModel call: ~",as.character(m_lvm_3$formula)[2],
                                     "\nFamily: ",as.character(m_lvm_3$family),". ",
                                     "Random row effects: ",as.character(m_lvm_3$call)[7]),
                    theme = theme(plot.title = element_text(size = 16, face="bold")))

pdf(file = "figs/coef_tax_all_unordered_v2.pdf",width=16,height=8)
print(final_plot)
dev.off()

### to do:
### look at functional groups(lifeforms)
# master list saved in:
# \\prodds.ntnl\Shared\AN\KFH\Groups\N_Marine\07 Training & Reference Documents\A&R Technical Guidance\Traits, Lifeforms etc\Plankton Lifeform Extraction Tool
# consider reproducing ordination in 3 dimensions using rgl (see: https://riffomonas.org/code_club/2021-03-24-rgl)

# PRIORITY : REPRODUCE CODE ####
## Currently untidy and seems to produce 'issues'
## Errors alluding to "non-numeric argument to binary operator"

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
rm(list = ls(pattern = "^ord"))
rm(list = ls(pattern = "^mv"))
rm(list = ls(pattern = "^xx"))
rm(list = ls(pattern = "^logs"))
rm(list = ls(pattern = "^ci_"))
rm(list = ls(pattern = "^max_"))
rm(list = ls(pattern = "^min_"))
rm(list = ls(pattern = "*plot"))
rm(list = ls(pattern = "*consiste"))
rm(datfol,nit,perms, ppi,replace_values,pl_ts_N,subset_data,i,level,
   ntrt,sbtt,srt,ttl)
