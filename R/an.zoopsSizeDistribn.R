# an.zoopsSizeDistribn.R ####
# Initial analyses of size distributions of zooplankters

# 00 set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "tictoc","gllvm","patchwork", "lubridate","ggpubr",
             "cowplot","ggridges", "seas")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

## set universals ####
tictoc::tic.clearlog();tic("Set universals");print("set universals")
source("R/set_meta.R")
toc(log=TRUE)

tic("Load zoop data")
source("R/imp.load_data_all_zoops.R")
toc(log=TRUE)
rm(dfl0)


## initial look at a single sample
dfl %>% 
  filter(.,Pot.Number==5) -> dftmp

dftmp %>% mutate(wt_len = Abund_m3*mnlongMaxAxis_mm) ## weight length by abundance

ggplot(dftmp, aes(x = mnlongMaxAxis_mm)) +
  # geom_histogram(fill = "skyblue", color = "black") +  # Adjust binwidth as needed
  geom_density(size=1)+
  labs(title = "Histogram of Species Lengths (unweighted)",
       subtitle = "Unweighted by abundance (effectvely presence-absence data)",
       x = "Species Length", 
       y = "Total Abundance") +
  ggthemes::theme_few()

dftmp %>% 
  ggplot(., aes(x=mnlongMaxAxis_mm, weight = Abund_m3))+
  # geom_histogram(fill="skyblue", colour = 1)+
  geom_density(size=1, aes(weight=Abund_m3))+
  labs(title = "Histogram of Species Lengths Weighted by Abundance",
       subtitle = "Weighted by taxon abundances",
       x = "Species Length (mm)", 
       y = "Total Abundance") +
  ggthemes::theme_few()

dftmp %>% 
  ggplot(., aes(x=mnlongMaxAxis_mm, weight = mn_carbTot_m3))+
  # geom_histogram(fill="skyblue", colour = 1)+
  geom_density(size=1,aes(weight = mn_carbTot_m3))+
  labs(title = "Histogram of Species Lengths Weighted by total carbon content",
       subtitle = "Weighted by carbon content",
       x = "Species Length (mm)", 
       y = "Total Abundance") +
  ggthemes::theme_few()

dftmp2 <- dftmp  %>% 
  filter(., !is.na(mnlongMaxAxis_mm)) %>% 
  filter(., !is.na(mn_carbTot_m3))

total_weight_carb <- sum(dftmp2$mn_carbTot_m3); rm(dftmp2)

dftmp %>% #names()
  filter(., !is.na(mnlongMaxAxis_mm)) %>% 
  filter(., !is.na(mn_carbTot_m3)) %>% 
  ggplot(., aes(x=mnlongMaxAxis_mm, weight = mn_carbTot_m3))+
  # geom_histogram(fill="skyblue", colour = 1)+
  geom_histogram(
    aes(y= ..density.. * total_weight_carb),
                 fill="skyblue", colour = 1)+
  geom_density(aes(y = ..density.. * total_weight_carb), color = "red", size = 1) +  # Rescaled density curve
    labs(title = "Histogram of Species Lengths Weighted by total carbon content",
       subtitle = "Weighted by carbon content",
       x = "Species Length (mm)", 
       y = "Total Abundance") +
  ggthemes::theme_few()

#### generate plot for each Biosys code 
### 1. generate new variable: yyy_mm
dfl$yyy_mm <- paste0(year(dfl$sample.date),"_",sprintf("%02d", month(dfl$sample.date)))

dftmp <- dfl %>% 
  filter(.,BIOSYS.Code == unique(dfl$BIOSYS.Code)[11])

dftmp <- dftmp %>% 
  filter(., !is.na(mnlongMaxAxis_mm)) %>% filter(., !is.na(mn_carbTot_m3))

twt <- sum(dftmp$mn_carbTot_m3)

dftmp %>% #names()
  ggplot(., aes(
    x=mnlongMaxAxis_mm,
    y=yyy_mm,
    # y=DJF,
    weight = mn_carbTot_m3))+
  geom_density_ridges(alpha=0.7, aes(fill=DJF),
                      jittered_points=TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape= "|",
                      point_size = 3,
                      point_alpha=1
                      )+
  ggthemes::theme_few()+
  xlim(-10,120)+
  labs(x = "Mean taxon length (mm)",
       subtitle = paste0(unique(dftmp$Region)[1],": ", unique(dftmp$WB)[1],
                         " (",unique(dftmp$BIOSYS.Code)[1],")"))+
  scale_fill_manual(values=c("skyblue","lightgreen","orange","sienna"))+
  scale_y_discrete(limits=rev)+
  theme(axis.title.y = element_blank(),
        legend.position = "none")

## generate pdfs of individual BIOSYS sites
