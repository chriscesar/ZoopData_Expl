### an.EA_taxa_RDA.R ####
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

#################################################################
# RDA ####
#see:
# https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
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

# extract taxon density data
dfw %>% 
  dplyr::select(-c(Pot.Number:Category,
                   "WIMS.Code.y":"Zinc, Dissolved_ug/l")
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

### create scaled version for comparison of effects on model
df_wims_w_trim0 %>% 
  mutate_if(is.numeric,scale) -> df_wims_w_trim0_scale

colnames(df_tx_w_trm) <- make.cepnames(colnames(df_tx_w_trm))
rda1 <- vegan::rda(df_tx_w_trm ~., data=df_wims_w_trim0_scale[,c(1:16)])

summary(rda1)
summary(rda1)[1]
plot(rda1)

vegan::RsquareAdj(rda1)
vegan::anova.cca(rda1) #'full' model is significant
vegan::anova.cca(rda1,by="term") #variables: chla, NGR_n, DIN, O2, Si & Region are significant
vegan::anova.cca(rda1,by="axis") #[WARNING: Takes ages!] No axes are significant (RDA1 P = 0.052)

# Custom triplot code!

## extract % explained by the first 2 axes
perc <- round(100*(summary(rda1)$cont$importance[2, 1:2]), 2)

## extract scores - these are coordinates in the RDA space
sc_si <- scores(rda1, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(rda1, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(rda1, display="bp", choices=c(1, 2), scaling=1)

plot(rda1)

## Custom triplot, step by step
scores(rda1)
sp_scores <- as.data.frame(scores(rda1)$species)
st_scores <- as.data.frame(scores(rda1)$sites)

sp_scores$nm <- vegan::make.cepnames(row.names(sp_scores))
st_scores$WB <- df_tx_w$WB
st_scores$Region <- df_tx_w$Region
rda_terms <- as.data.frame(scores(rda1)$biplot)


# Set up a blank plot with scaling, axes, and labels

png(file = "figs/zoopRDA_taxa.png",
    width=12*ppi, height=10*ppi, res=ppi)
pl <- plot(rda1,
     scaling = 1, # set scaling type 
     type = "text", # this excludes the plotting of any points from the results
     frame = FALSE,
     # set axis limits
     # xlim = c(-1,1), 
     # ylim = c(-300,100),
     # label the plot (title, and axes)
     # main = "Triplot RDA - scaling 1",
     xlab = paste0("RDA1"), 
     ylab = paste0("RDA2") 
)
dev.off()

# add points for site scores
points(sc_si, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black", # outline colour
       bg = "steelblue", # fill colour
       cex = 1.2) # size
# add points for species scores
points(sc_sp, 
       pch = 22, # set shape (here, square with a fill colour)
       col = "black",
       bg = "#f2bd33", 
       cex = 1.2)
# add text labels for species abbreviations
text(sc_sp + c(0.03, 0.09), # adjust text coordinates to avoid overlap with points 
     labels = vegan::make.cepnames(rownames(sc_sp)), 
     col = "grey40", 
     font = 2, # bold
     cex = 0.6)
# add arrows for effects of the expanatory variables
arrows(0,0, # start them from (0,0)
       sc_bp[,1], sc_bp[,2], # end them at the score value
       col = "red", 
       lwd = 3)
# add text labels for arrows
text(x = sc_bp[,1] -0.1, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bp[,2] - 0.03, 
     labels = rownames(sc_bp), 
     col = "red", 
     cex = 1, 
     font = 2)

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
   ntrt,sbtt,srt,ttl,rcov0,rcov1,sDsn)
