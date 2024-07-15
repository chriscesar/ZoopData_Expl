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
source("R/imp.load_data_taxa.R")
source("R/imp.load_data_wims.R")

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
# S
dfw %>% 
  ggplot(., aes(x=DJF, y = S)) +
  geom_violin(adjust = 1)+
  geom_point(position=position_jitter(h=0,w=0.3,seed = pi))+
  ylim(0,NA)+
  facet_wrap(.~Region, drop = FALSE)+
  theme(axis.title.x = element_blank())

#N
dfw %>% 
  ggplot(., aes(x=DJF, y = log(N))) +
  geom_violin(adjust = 1)+
  geom_point(position=position_jitter(h=0,w=0.3,seed = pi))+
  ylim(0,NA)+
  facet_wrap(.~Region, drop = FALSE)+
  theme(axis.title.x = element_blank())

#HillsN1
dfw %>% 
  ggplot(., aes(x=DJF, y = HillsN1)) +
  geom_violin(adjust = 1)+
  geom_point(position=position_jitter(h=0,w=0.3,seed = pi))+
  ylim(0,NA)+
  facet_wrap(.~Region, drop = FALSE)+
  theme(axis.title.x = element_blank())

dfw %>% 
  ggplot(., aes(x=month, y=HillsN1))+  
  coord_polar(start=0,direction=1) +  
  geom_point(aes(color=DJF),show.legend = FALSE) +  
  # geom_point(aes(y=pnt), colour="red") +
  theme_light()+
  facet_wrap(.~Region)+
  theme(axis.title.x = element_blank())
