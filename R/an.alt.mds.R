 # an.alt.mds.R####
## alternative calcs of nmMDS plots.  PRIMER analysis generated stress values
## of approx 11.  vegan(metaMDS) stress values ~28!
 
# load packages ####
ld_pkgs <- c("tidyverse","MASS","lubridate","vegan","mvabund","seas","patchwork",
             "ecoCopula","performance","gclus","corrplot","gllvm","tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)
tictoc::tic.clearlog();tic("set universals");print("set universals")

# set universals ####
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

# Load data ####
## load zoops ####
source("R/imp.load_data_all_zoops.R")

## load wims ####
source("R/imp.load_data_wims.R")

# TRAITS ####
## widen data ####
dfl %>% dplyr::select(.,c("Pot.Number":"WB",LF02,Abund_m3)) %>% 
  dplyr::select(.,-c(Aphia.ID,Taxa,CEA.Notes)) %>% 
  group_by(across(c(-Abund_m3))) %>% #names()# group by everything except abundance
  summarise(Abund_m3=sum(Abund_m3), #sum abundances
            .groups="drop") %>% 
  pivot_wider(.,names_from = LF02, values_from = Abund_m3, values_fill = 0) -> dfw_trt

dfw_trt %>% #names()
  dplyr::select(.,-c(Pot.Number:WB)) -> dfw_trt_tmp

## calculate bray curtis distance based on log-transformed data
dfw_trt_tmp_log_dist <- vegdist(log(dfw_trt_tmp+1), method = "bray")

## 1. MASS::isoMDS ####
tic("Traits nmds: MASS")
nmds_result_MASS <- isoMDS(dfw_trt_tmp_log_dist, k = 2)  # 2D NMDS
toc(log=TRUE)
nmds_result_MASS$stress

# Extract site coordinates (points represent sites)
site_coords <- as.data.frame(nmds_result_MASS$points)

# To get species coordinates, you can use weighted averages of site scores
species_coords <- as.data.frame(t(aggregate(dfw_trt_tmp,
                                            by = list(rowMeans(nmds_result_MASS$points)),
                                            FUN = mean)))

# Plot site coordinates with ggplot2
ggplot(site_coords, aes(x = V1, y = V2)) +
  geom_point()
rm(site_coords,species_coords,nmds_result_MASS)

## 2. labdsv::nmds() ####
tic("Traits nmds: labdsv")
nmds_result_labdsv <- labdsv::nmds(dfw_trt_tmp_log_dist, k=2)
toc(log = TRUE)

# Extract scores (species and site coordinates)
species_scores <- nmds_result_labdsv$species
site_scores <- nmds_result_labdsv$points
colnames(site_scores) <- c("nMDS1","nMDS2")
ggplot(site_scores, aes(x=nMDS1,y=nMDS2))+
  geom_point()
nmds_result_labdsv$stress
rm(nmds_result_labdsv,site_scores,species_scores)

## 3. ecodist::nmds() ####
tic("Traits nmds: ecodist")
nmds_result_ecodist <- ecodist::nmds(dfw_trt_tmp_log_dist, mindim = 2, maxdim = 2,trace = TRUE)
str(nmds_result_ecodist)
min(nmds_result_ecodist$stress)
rm(nmds_result_ecodist)
toc(log=TRUE)

## 4. vegan::metaMDS() ####
tic("Traits nmds: vegan")
nmds_result_vegan <- vegan::metaMDS(dfw_trt_tmp_log_dist,k=2, trymax = 500)
nmds_result_vegan$stress
toc(log=TRUE)
