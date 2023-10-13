# mds_zoops_3d.R ####
### generate 3 dimensional NMDS plot

#### load packages ####
ld_pkgs <- c("tidyverse","vegan","rgl")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

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

### load data ####
### taxon data
df0 <- as_tibble(read.csv(file=paste0(datfol,
                                      "processedData/",
                                      "zoopWIDEAbund_m3_taxOnly_USE.csv")))

#### quick ordinations ####
df0 %>% 
  dplyr::select(-c(1:21)) %>% ###remove metadata info
  dplyr::select_if(~ !is.numeric(.) || sum(.) !=0) %>% ### drop cols summing to 0
  filter(rowSums(across(where(is.numeric)))!=0) -> dftmp ###remove 'empty' rows

### 3D NMDS ####
ptm <- Sys.time()###
set.seed(pi+5);ord <-   vegan::metaMDS(dftmp, trymax = 1000, k = 3)
Sys.time() - ptm;rm(ptm)
saveRDS(ord, file = "data/out/ord3d.Rdata")

plot(ord)

#### extract ordination axes ####
scores_site <- df0 %>% 
  dplyr::select(c(1:21))
tmp_sites <- as_tibble(as.data.frame(scores(ord,display = "site")))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
scores_site$NMDS1 <- tmp_sites$NMDS1 #location of individual samples in NMDS space
scores_site$NMDS2 <- tmp_sites$NMDS2 #location of individual samples in NMDS space
scores_site$NMDS3 <- tmp_sites$NMDS3 #location of individual samples in NMDS space
#assign colours to variable 'Region'
scores_site$RegionCol <- case_when(scores_site$Region == "Anglian" ~ cbPalette[1],
                                   scores_site$Region == "NEast" ~ cbPalette[2],
                                   scores_site$Region == "NWest" ~ cbPalette[3],
                                   scores_site$Region == "Southern" ~ cbPalette[4],
                                   scores_site$Region == "SWest" ~ cbPalette[5],
                                   scores_site$Region == "Thames" ~ cbPalette[6],
                                   TRUE ~ NA_character_)
saveRDS(scores_site, file = "data/out/scores_site3d.Rdata")
rm(tmp_sites)

# Using the scores function from vegan to extract the species scores and
# convert to a data.frame
scores_species <- as.data.frame(scores(ord,display = "species"))
scores_species$lbfull <-  row.names(scores_species)
scores_species$lb <-  make.cepnames(row.names(scores_species))#shorten names
saveRDS(scores_species, file = "data/out/scores_species3d.Rdata")

plot3d(x = scores_site$NMDS1,
       y = scores_site$NMDS2,
       z = scores_site$NMDS3,
       col=scores_site$RegionCol,
       type="s", size = 1,
       xlab = "NMDS1",ylab = "NMDS2",zlab = "NMDS3")

# movie3d(spin3d(axis = c(0,0,1),rpm=4), duration = 15, dir="./")
# plot3d(x = scores_species$NMDS1,
#        y = scores_species$NMDS2,
#        z = scores_species$NMDS3,
#        type="s",size = 1,
#        xlab = "NMDS1",ylab = "NMDS2",zlab = "NMDS3")

text3d(x = scores_species$NMDS1,
       y = scores_species$NMDS2,
       z = scores_species$NMDS3,
       scores_species$lb, size=1,
       cex = .75, col="grey")

axes3d();title3d(xlab="NMDS1",
                 ylab="NMDS2",
                 zlab="NMDS3",
                 font=2)
