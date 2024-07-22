# an.EA_lifeforms_trends.R ####
# Analysis of temporal and spatial trends in zooplankton lifeforms #

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
source("R/imp.load_data_lifeforms.R")
source("R/imp.load_data_wims.R")
toc(log=TRUE)

# join & save data ####  
tic("Join taxon & WIMS data. Generate LONG version")
dfw <- left_join(df_tx_w,df_wims_w,by="PRN")

### generate LONG version WITH zero values (for calculation of means)
rm(df_tx_w, df_tx, df_tx_100um, df_wims, df_wims_w)
toc(log=TRUE)

# comparing large and small copepod abundances
summary(dfw$Cop_Lg)
summary(dfw$Cop_Sm)
summary(dfw$Cop_Ambi)
summary(dfw$Cop_NYA)

# Large:Small copepods ratio
hist(log(dfw$Cop_Lg)+1,breaks = 30)
hist(log(dfw$Cop_Sm)+1,breaks = 30)

LgSm <- (dfw$Cop_Lg+.01)/(dfw$Cop_Sm+.01)
SmLg <- (dfw$Cop_Sm+.01)/(dfw$Cop_Lg+.01)

dfw %>% 
  mutate(LgSm = LgSm, SmLg=SmLg) %>% 
  ggplot(., aes(x=DJF,
                y=SmLg))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0,seed = pi))+
  facet_wrap(.~Region)+
  labs(title = "Proportion of small/large copepods recorded in zooplankton data",
       x=NULL,
       y="Small/Large copepods")+
  theme(strip.text = element_text(face=2),
        axis.title.y = element_text(face=2),
        axis.text.x = element_text(face=2))
