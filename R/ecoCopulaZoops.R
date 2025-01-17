# library(ecoCopula)
# library(tidyr)
# library(tidygraph)
# library(ggraph)

# load packages ####
ld_pkgs <- c("tidyverse", "tictoc","ggraph","ecoCopula")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)


### load data ####
tic("load data sets")
source("R/set_meta.R")
source("R/imp.load_data_lifeforms.R")
source("R/imp.load_data_wims.R")
toc(log=TRUE)

# join data ####  
tic("Join taxon & WIMS data. Generate LONG version")
dfw <- left_join(df_tx_w,df_wims_w,by="PRN")

rm(df_tx_w, df_tx, df_tx_100um, df_wims, df_wims_w)
toc(log=TRUE)

# prep data ####
tic("prep data")
## abundance data
dfabd <- dfw[,c(25:61)]
n <- 20
dfabd_trm <- dfabd %>% 
  dplyr::select(which(apply(.,2,function(x) sum(x !=0) >n)))


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
X <- dfw[,kp]
rm(kp, keep)

### replace LESS THAN values with value*0.5
# define function to replace "less than" values by  half their value
replace_values <- function(x) {
  if_else(str_detect(x, "^<"), as.numeric(sub("^<", "", x))/2, as.numeric(x))
}

# amend < values in wims data
X %>% 
  mutate_all(.,replace_values) -> X

## replace NA values with mean values for respective column.
## Leave non-numeric cols unchanged
X %>% 
  mutate(across(where(is.numeric),
                ~replace(.,
                         is.na(.),
                         mean(.,
                              na.rm = TRUE) ) ) ) -> X

## rename columns
X <- X %>% 
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
X %>% 
  mutate_if(is.numeric,scale) -> X_scale

X_use <- X_scale[,c(1:2,7,9,12,14,15)]
toc(log=TRUE)

## graphical model
tic("fit marginal model")
# fit marginal model
fit_nb <- stackedsdm(dfabd_trm,~., data = X_use, family="negative.binomial", ncores = 2) #eqiv. manyglm()
# fit copula ordination 
fit_gr <- cgr(fit_nb, seed=3)
# biplot
plot(fit_gr, pad=0)
# ggsave(filename = "figs/zoop_all_copula.pdf",width = 12,height = 12,units = "in",plot = pl)

igraph_out <- fit_gr$best_graph$igraph_out
igraph_out %>% ggraph('fr') + # see ?layout_tbl_graph_igraph
  geom_edge_fan0(aes( colour = partcor, width=partcor)) +
  scale_edge_width(range = c(0.5, 3))+
  scale_edge_color_gradient2(low="#b2182b",
                             # mid="white",
                             mid="darkgrey",
                             high="#2166ac")+
  geom_node_text(aes(label=name), repel = TRUE)+
  geom_node_point(aes(size=1.3),show.legend = FALSE)+
  labs(caption = "Taxa pairs that have no direct edge between them have no direct association. Any co-occurrence patterns they have are due to associations they both have with other taxa in the data (i.e., mediator species).
       Taxa with direct edges between them have positive (blue), or negative (pink) associations even after controlling for mediation effects of all other taxa.")+
  theme_void() +
  theme(
    # legend.position = 'none',
    panel.background = element_rect(fill="azure2",
                                    colour="azure2")
    ) +
  guides(width = "none") -> pl
ggsave(filename = "figs/zoop_all_copula.pdf",width = 20,height = 12,units = "in",plot = pl)
toc(log=TRUE)
