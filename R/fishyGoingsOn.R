# fishyGoingsOn.R ####
## analysis of fish eggs and larvae

#### load packages ####
ld_pkgs <- c("tidyverse")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

### set up folders & import functions ###
source("R/folder.links.R")

theme_set(ggthemes::theme_few())###set theme for all ggplot objects
nit <- 9999 #number of iterations
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

### import data ####
## data are in LONG format with *tentative* lifeform values attached
df0 <- as_tibble(readxl::read_excel(paste0(datfol,
                                      "processedData/",
                                      "MBA_Returns_Amalgamated_USE.xlsx"),
                                    sheet = "outR04_LF"))

### to do ####
# flag 'fish' type rows
# assign variable "type" to indicate whether it's an EGG or LARVAE
# (think about what to do in cases where just the taxon name is provided,
# e.g. "Gobiidae")
# Summarise counts by sample (i.e. Total larvae and Total egg values per sample)
# calculate mean values ± SD by WB
# present, e.g., stacked larvae-egg barchart by WB arranged by region?