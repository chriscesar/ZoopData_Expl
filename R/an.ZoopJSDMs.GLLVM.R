# an.ZoopJSDMs.GLLVM.R (file name) #### 
# (Experimental) Fun with GLLVMs
# Initial tinkerings with a Joint Species Distribution Model
# Code/idea stolen from https://jenniniku.github.io/gllvm/articles/vignette1.html

ld_pkgs <- c("stringr", "tidyverse", "ggplot2", "patchwork",
             "gllvm", "corrplot", "gclus") # what packages do we need to load?
vapply(ld_pkgs, library, logical(1L), # load them and display TRUE/FALSE if loaded
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

## set metadata ####
cbPalette <- c("#0072B2","#e79f00","#009E73", "#9ad0f3", "#000000", 
               "#D55E00", "#CC79A7", "#F0E442") # set up colour palette
ppi <- 300 # resolution of exported images
theme_set(ggthemes::theme_few())###set theme for all ggplot objects

# set data folder:
source("R/folder.links.R")

# Import data ####
#keep '0' version unedited
dfw0 <- as_tibble(read.csv(paste0(datfol,
                                  "processedData/zoopWIDEAbund_m3_WIMS_USE.csv")
                           )
                  )


# tidy up data
replace_values <- function(x) { #function to replace "less than" values by  half their value
  if_else(str_detect(x, "^<"), as.numeric(sub("^<", "", x))/2, as.numeric(x))
}

dfw <- dfw0 %>% 
  filter(.,!is.na(Nitrite..Filtered.as.N_mg.l)) %>% # remove chunk of NA's at start of data
  dplyr::select(-c(Benzo.a.Pyrene_ug.l:X1.1.1.Trichloroethane_ug.l)) %>%
  ### if a cell contains "<", then replace that value with (x/2)
  ### E.g. "<0.2" becomes 0.1
  mutate_at(vars(Ammoniacal.Nitrogen..Filtered.as.N_mg.l:Tributyl.Tin.as.Cation_ug.l), replace_values)

## remove columns that sum to zero
dfw_numeric_filtered <- as_tibble(dfw) %>% 
  select_if(is.numeric) %>% 
  select_if(function(x) sum(x, na.rm=TRUE) !=0)

dfw_char <- as_tibble(dfw) %>% 
  dplyr::select_if(is.character)# %>% 

dfw <- as_tibble(cbind(dfw_char,dfw_numeric_filtered))
rm(dfw_char,dfw_numeric_filtered)

# build gllvm ####
# Y - abundances ####
# isolate taxon abundance data & retain taxa which occur in 5 or mode samples
Y <- as_tibble(dfw) %>% 
  dplyr::select(Acartia.clausi:"Daphnia..Freshwater.")

tmp <- Y %>% 
  mutate_if(is.numeric,~1 * (.>0)) %>% ## convert to presence-absence
  colSums()
tmp <- tmp>4 #ID which columns have >4 observations

Y <- Y[,tmp]### keep only those with >4
rm(tmp)

X_all <- as_tibble(dfw) %>%  #all non-biological variables
  dplyr::select(-c(27:224))#### NEED TO RECHECK THIS EACH TIME!

X <- X_all %>% ### env and sea
  dplyr::select(Region,
                Nitrite..Filtered.as.N_mg.l,
                Salinity...In.Situ_ppt,
                Silicate..Filtered.as.SiO2_mg.l,
                Turbidity...In.Situ_FTU,
                Ammoniacal.Nitrogen..Filtered.as.N_mg.l) %>% 
  rename(.,
         SiO3 = Silicate..Filtered.as.SiO2_mg.l,
         Salinity = Salinity...In.Situ_ppt,
         NO3 = Nitrite..Filtered.as.N_mg.l,
         Turb = Turbidity...In.Situ_FTU,
         NH3N = Ammoniacal.Nitrogen..Filtered.as.N_mg.l) %>% 
  mutate(across(c(Region), factor))

X_2 <- X %>% #exclude 'Region' term
  dplyr::select(SiO3,
                Salinity,
                NO3,
                Turb,
                NH3N)
X_2 <- scale(as.matrix(X_2))

# simple gllvm model ####
##identify NA data in environmental data
na_counts <-  X %>%
  summarise_all(~ sum(is.na(.)))

na1 <- which(is.na(X$Turb))
na2 <- which(is.na(X$Salinity))
na_out <- unique(c(na1,na2))
rm(na1,na2)

## Fit null model with no env predictors ####
tmp_longspp <- names(Y) ## create long version of taxon names
names(Y) <- vegan::make.cepnames(names(Y))

ptm <- Sys.time()
gllvm_null <- gllvm(Y, family = "gaussian", seed = pi, num.lv = 2)
Sys.time()-ptm;rm(ptm)
saveRDS(gllvm_null, file = "data/out/gllvm_null.Rdat")
gllvm_null <- readRDS("data/out/gllvm_null.Rdat")
ordiplot(gllvm_null, biplot=TRUE)

## by Region ####
ptm <- Sys.time()
gllvm_rgn <- gllvm(y = Y[-na_out,],
                   X = X[-na_out,],
                   set.seed = pi,
                   formula = ~ Region,
                   num.lv = 2,
                   family="gaussian")
Sys.time()-ptm;rm(ptm)
saveRDS(gllvm_rgn, file = "data/out/gllvm_rgn.Rdata")
gllvm_rgn <- readRDS("data/out/gllvm_rgn.Rdata")
ordiplot(gllvm_rgn, biplot=TRUE)


## with env predictors ####
ptm <- Sys.time()
gllvm_env <- gllvm(y = Y[-na_out,],
                   X = X[-na_out,],
                   set.seed = pi,
                   formula = ~ Region + NO3 + NH3N + SiO3 + Salinity,
                   num.lv = 2,
                   family="gaussian")
Sys.time()-ptm;rm(ptm)
saveRDS(gllvm_env, file = "data/out/gllvm_env.Rdata")
gllvm_env <- readRDS("data/out/gllvm_env.Rdata")
ordiplot(gllvm_env, biplot=TRUE)

# par(mfrow = c(3,2), mar=c(4,4,2,2))
# for(i in 1:ncol(X)){
#   Col <- cbPalette[as.numeric(cut(X[,i], breaks = 20))]
#   ordiplot(gllvm_env, symbols = T, s.colors = Col, main = colnames(X)[i], 
#            biplot = TRUE)
# }
