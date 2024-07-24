# an.EA_lifeforms_JDSM.R ####
# Comparison of Joint Species Distribution Model approaches for zooplankton
# lifeforms data

# load packages ####
ld_pkgs <- c("tidyverse", "tictoc","gllvm","patchwork", "lubridate", "Hmsc",
             "ggpubr")
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

# join data ####  
tic("Join taxon & WIMS data. Generate LONG version")
dfw <- left_join(df_tx_w,df_wims_w,by="PRN")

rm(df_tx_w, df_tx, df_tx_100um, df_wims, df_wims_w)
toc(log=TRUE)

#### 
tic("Split data into respective objects")
## ecological data (Y) ####
dfw %>% 
  dplyr::select(.,-c('Pot.Number':'WB','WIMS.Code.y':last_col())) -> Y

## metadata  ####
mt <- dfw %>% 
  dplyr::select(.,c('Pot.Number':'WB'))

## Environmental parameters (X)  ####
dfw %>% 
  dplyr::select(.,-c("1,1,1-Trichloroethane_ug/l":last_col())) -> X

toc(log=TRUE)

# Model set up: initial parameters ####
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

### append region & WB
X$Region <- dfw$Region
X$WB <- dfw$WB

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
toc(log=TRUE)

################################
# Unconstrained GLLVM model ####
################################
tic("gllvm_uncon_offset_nb: Unconstrained NB with offset")
runs <- 5
sDsn <- data.frame(Region = X$Region)

m_lvm_0nb <- gllvm(y = Y,
                   family = "negative.binomial",
                   offset = log(mt$`Net.volume.sampled.(m3)`),#logging offset to match poisson link
                   studyDesign = sDsn, row.eff = ~(1|Region),
                   num.lv = 2,seed = pi,
                   starting.val = "random",
                   n.init = runs,#re-run model to get best fit
                   trace = TRUE)

saveRDS(m_lvm_0nb, file = "figs/gllvm_logOffset_uncon_nb.Rdat")
toc(log = TRUE)

m_lvm_0nb <- readRDS("figs/gllvm_logOffset_uncon_nb.Rdat")
par(mfrow=c(2,2));plot(m_lvm_0nb,which = 1:4);par(mfrow=c(1,1))
ordiplot.gllvm(m_lvm_0nb,biplot=TRUE)
gllvm::ordiplot(m_lvm_0nb,predict.region=TRUE, symbols=TRUE)

################################
# Unconstrained HMSC model #####
################################
## Hmsc is primarily geared to modelling abundances against
## predictor variables.
## However, if we can model abundances against a constant.  This model essentially
## becomes a latent variable model for the community data.
## Note: HMSC does not allow negative binomial, instead uses
## lognormal poisson

## WARNING: Hmsc does not currently (Jul-24) and hasn't since it's release in ~2019
## incorporate an offset term for (e.g.) incorporating sampling effort as an error term

## In this case therefore, zooplankton count data are not controlled by Net Volume.
## One solution would be to recalc zoop counts to convert to standardised
## abundances per m3 and model as a Tweedie or similar.

tic("hmsc_uncon_offset_nb: Unconstrained NB with offset")
# create model
Yuse <- as.matrix(Y)
XData <- data.frame(intercept = rep(1, nrow(Yuse)))

### create Hmsc object
studyDesign <- mt %>%
  dplyr::select(Region,Pot.Number) %>%
  mutate(Region = as.factor(Region)) %>%
  mutate(Pot.Number = as.factor(Pot.Number)) %>%
  data.frame()
ranlevels <- HmscRandomLevel(units = studyDesign$Region)

simul <- Hmsc(Y = Yuse, XData=XData,distr = "lognormal poisson",
              studyDesign=studyDesign,
              ranLevels = list("Region" = ranlevels))

# Run model
thin <- 10
samples <- 10000
nChains <- 3
transient <- ceiling(thin*samples*.5)

tic("hmsc_uncon_offset_nb: Unconstrained NB with offset")
# mod_HMSC <- sampleMcmc(simul,
#                        samples = samples,
#                        thin = thin,
#                        transient = transient,
#                        nChains = nChains#,
#                        # nParallel = nChains
#                        )
# saveRDS(mod_HMSC, file = paste0("outputs/modUnconHMSC","_smp",samples,
#                                 "_thn",thin,"_trns",transient,"_chn",nChains,
#                                 ".Rdata"))
toc(log = TRUE)
mod_HMSC <- readRDS( paste0("outputs/modUnconHMSC","_smp",samples,
                                "_thn",thin,"_trns",transient,"_chn",nChains,
                                ".Rdata"))

## investigate model outputs ####
mpost <- convertToCodaObject(mod_HMSC) # model diagnostics/convergence
preds <- computePredictedValues(mod_HMSC) # model performance
MF <- evaluateModelFit(hM=mod_HMSC, predY = preds) # r2, etc
VP <- computeVariancePartitioning(mod_HMSC) # variance partitioning

ess.beta <- effectiveSize(mpost$Beta) %>% 
  as_tibble() %>% dplyr::rename(ess_beta=value)
psrf.beta <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf %>% 
  as_tibble() %>% dplyr::rename(psrf_beta = "Point est.")


(diag_all <- ggarrange(ggplot(ess.beta,aes(x=ess_beta))+
                        geom_histogram()+
                        xlab("Effective sample size"),
                      ggplot(psrf.beta,aes(x=psrf_beta))+
                        geom_histogram()+
                        geom_vline(xintercept = 1.001, col=2)+
                        xlab("Gelman diagnostic"), align = "v")+
  ggtitle("All plots"))

# species niches
postBeta <- getPostEstimate(mod_HMSC, parName = "Beta")

plotVariancePartitioning(mod_HMSC, VP=VP, las=2, horiz=TRUE)
plotBeta(mod_HMSC,post=postBeta, param = "Support", #supportLevel=0.95,
         split = .4, spNamesNumbers = c(T,F))
biPlot(mod_HMSC,etaPost = getPostEstimate(mod_HMSC,"Eta"),
       lambdaPost = getPostEstimate(mod_HMSC, "Lambda"))
