# an.EA_zoops_trends.R ####
# Analysis of temporal and spatial trends in zooplankton taxa #

# load packages ####
ld_pkgs <- c("tidyverse", "tictoc","gllvm","patchwork", "lubridate")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

#### set universals ####
tictoc::tic.clearlog();tic("Set universals");print("set universals")
source("R/set_meta.R")
toc(log=TRUE)

### load data ####
tic("load data sets")
source("R/imp.load_data_taxa.R")
source("R/imp.load_data_wims.R")
toc(log=TRUE)

# join & save data ####  
tic("Join taxon & WIMS data. Generate LONG version")
dfw <- left_join(df_tx_w,df_wims_w,by="PRN")

### generate LONG version WITH zero values (for calculation of means)
rm(df_tx_w, df_tx, df_tx_100um, df_wims, df_wims_w)
toc(log=TRUE)

tic("calculate indices")
# Calculate indices ####
## taxon richness
S <- vegan::specnumber(dfw %>% dplyr::select(.,-c(Pot.Number:Category,
                                                  WIMS.Code.y:last_col())))
## taxon density (per m3)
dfw %>% 
  dplyr::select(.,-c(Pot.Number:Category,
                     WIMS.Code.y:last_col())) %>%
  rowSums(.) -> Nraw
N <- Nraw/dfw$`Net.volume.sampled.(m3)`

## Shannon diversity
dfw %>% 
  dplyr::select(.,-c(Pot.Number:Category,
                     WIMS.Code.y:last_col())) %>%
  vegan::diversity(x=.,index = "shannon") -> H

## Hill's N1
HillsN1 <- exp(H)

dfw$S <- S
dfw$N <- N
dfw$H <- H
dfw$HillsN1 <- HillsN1
rm(S,N,H,HillsN1)
dfw %>% relocate(., S,N,H,HillsN1, .before = WIMS.Code.y) -> dfw
toc(log=TRUE)

# Initial plots ####
tic("Generate initial plots")
month_levels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# S
dfw %>% 
  ggplot(., aes(x=DJF, y = S)) +
  geom_violin(adjust = 1)+
  geom_point(position=position_jitter(h=0,w=0.3,seed = pi))+
  ylim(0,NA)+
  facet_wrap(.~Region, drop = FALSE)+
  theme(axis.title.x = element_blank(),
        strip.text = element_text(face=2),
        axis.title.y = element_text(face=2),
        axis.text.x = element_text(face=2))+
  labs(title="Seasonal taxon richness of zooplankton monitored in different regions",
       caption=paste0("Samples gathered between ",
                      format(min(dfw$sample.date),
                      "%d/%m/%Y")," & ",format(max(dfw$sample.date),"%d/%m/%Y"))) -> pl
ggsave(plot=pl,filename = "figs/zoopViolinSppRich_by_Region.pdf",width = 12,height = 7,units = "in")
rm(pl)

dfw %>% 
  mutate(month_factor = factor(month_levels[month],
                               levels = month_levels,
                               ordered = TRUE)) %>% 
  ggplot(., aes(x=month_factor, y=S))+  
  coord_polar(start=0,direction=1) +  
  geom_point(aes(colour=DJF),show.legend = FALSE, size=3) +  
  scale_colour_manual(values = c("deepskyblue2","chartreuse2",
                                               "darkorange","sienna"))+
                                                 facet_wrap(.~Region)+
  theme_light()+
  xlab(NULL)+
  labs(title = "Polar plot of zooplankton taxon richness by month across EA survey regions",
       caption=paste0("Samples gathered between ",
                      format(min(dfw$sample.date),
                             "%d/%m/%Y")," & ",format(max(dfw$sample.date),"%d/%m/%Y")))+
  theme(strip.text = element_text(face=2, colour=1),
        axis.title.y = element_text(face=2),
        strip.background = element_rect(fill="white")) -> pl
ggsave(plot=pl,filename = "figs/zoopPolarSppRich_by_Region.pdf",width = 12,height = 7,units = "in")
rm(pl)

#N
dfw %>% 
  ggplot(., aes(x=DJF, y = log(N))) +
  geom_violin(adjust = 1)+
  geom_point(position=position_jitter(h=0,w=0.3,seed = pi))+
  ylim(0,NA)+
  facet_wrap(.~Region, drop = FALSE)+
  theme(axis.title.x = element_blank(),
        strip.text = element_text(face=2),
        axis.title.y = element_text(face=2),
        axis.text.x = element_text(face=2))+
  xlab(NULL)+
  labs(title="Seasonal taxon densities of zooplankton monitored in different regions",
       caption=paste0("Samples gathered between ",
                      format(min(dfw$sample.date),
                             "%d/%m/%Y")," & ",format(max(dfw$sample.date),"%d/%m/%Y")))->pl
ggsave(plot=pl,filename = "figs/zoopViolinSppDens_by_Region.pdf",width = 12,height = 7,units = "in")
rm(pl)

dfw %>% 
  mutate(month_factor = factor(month_levels[month],
                               levels = month_levels,
                               ordered = TRUE)) %>% 
  ggplot(., aes(x=month_factor, y=log(N)))+  
  coord_polar(start=0,direction=1) +  
  geom_point(aes(colour=DJF),show.legend = FALSE, size=3) +  
  scale_colour_manual(values = c("deepskyblue2","chartreuse2",
                                               "darkorange","sienna"))+
                                                 theme_light()+
  facet_wrap(.~Region)+
  labs(title = "Polar plot of (log) zooplankton density by month across EA survey regions",
       caption=paste0("Samples gathered between ",
                      format(min(dfw$sample.date),
                             "%d/%m/%Y")," & ",format(max(dfw$sample.date),"%d/%m/%Y")))+
  xlab(NULL)+
  theme(strip.text = element_text(face=2, colour=1),
        axis.title.y = element_text(face=2),
        strip.background = element_rect(fill="white")) -> pl
ggsave(plot=pl,filename = "figs/zoopPolarSppDens_by_Region.pdf",width = 12,height = 7,units = "in")
rm(pl)

#HillsN1
dfw %>% 
  ggplot(., aes(x=DJF, y = HillsN1)) +
  geom_violin(adjust = 1)+
  geom_point(position=position_jitter(h=0,w=0.3,seed = pi))+
  ylim(0,NA)+
  facet_wrap(.~Region, drop = FALSE)+
  theme(axis.title.x = element_blank(),
        strip.text = element_text(face=2),
        axis.title.y = element_text(face=2),
        axis.text.x = element_text(face=2))+
  ylab("Hill's N1")+
  xlab(NULL)+
  labs(title="Seasonal Hill's N1 values of zooplankton monitored in different regions",
       caption=paste0("Samples gathered between ",
                      format(min(dfw$sample.date),
                             "%d/%m/%Y")," & ",format(max(dfw$sample.date),"%d/%m/%Y"))) -> pl
ggsave(plot=pl,filename = "figs/zoopViolinHillsN1_by_Region.pdf",width = 12,height = 7,units = "in")
rm(pl)

dfw %>% 
  mutate(month_factor = factor(month_levels[month],
                               levels = month_levels,
                               ordered = TRUE)) %>% 
  ggplot(., aes(x=month_factor, y=HillsN1))+  
  coord_polar(start=0,direction=1) +  
  geom_point(aes(colour=DJF),show.legend = FALSE, size=3) +  
  # geom_point(aes(y=pnt), colour="red") +
  scale_colour_manual(values = c("deepskyblue2","chartreuse2",
                                               "darkorange","sienna"))+
  theme_light()+
  facet_wrap(.~Region)+
  ylab("Hill's N1")+
  xlab(NULL)+
  labs(title = "Polar plot of Hill's N1 values of zooplankton by month across EA survey regions",
       caption=paste0("Samples gathered between ",
                      format(min(dfw$sample.date),
                             "%d/%m/%Y")," & ",format(max(dfw$sample.date),"%d/%m/%Y")))+
  theme(strip.text = element_text(face=2, colour=1),
        axis.title.y = element_text(face=2),
        strip.background = element_rect(fill="white")) -> pl
ggsave(plot=pl,filename = "figs/zoopPolarHillsN1_by_Region.pdf",width = 12,height = 7,units = "in")
rm(pl)
toc(log=TRUE)
