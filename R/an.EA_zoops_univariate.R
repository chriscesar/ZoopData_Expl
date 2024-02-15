### 111111 ####
#### Univariate explorations of data imported in imp.EA_raw.R

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
S <- vegan::specnumber(df_tx_w[, -c(1:21)])
N <- rowSums(df_tx_w[, -c(1:21, length(df_tx_w))])

### remove species data
dfw %>% 
  dplyr::select(.,Pot.Number:Category,WIMS.Code.y:`Zinc, Dissolved_ug/l`) -> df
df$S <- S;df$N <- N; rm(S,N)


summary(m0 <- lmerTest::lmer(S ~ 1 + (1|Region), data=df))
