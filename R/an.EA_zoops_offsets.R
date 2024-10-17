#### load packages ####
ld_pkgs <- c("tidyverse","MASS","lubridate", "tictoc","gllvm","purrr","patchwork")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

#### set universals ####
tictoc::tic.clearlog();tic("set universals");print("set universals")
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

### load data ####
tic("Load taxon data and correct taxon names")
print("Load taxon data and correct taxon names")
### taxon data
df_tx0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,
                                               "processedData/MBA_Returns_Amalgamated_USE.xlsx"),
                                        sheet="outR04"))

### append updated taxon names
tx_chk0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,
                                                "processedData/MBA_Returns_Amalgamated_USE.xlsx"),
                                         sheet="TaxonomicRaw"))

tx_chk <- tx_chk0 %>% 
  rename(Taxa=ScientificName_accepted,
         Aphia.ID=AphiaID_accepted) %>% 
  #dplyr::select(.,-ScientificName) %>% 
  distinct()

tx_chk %>% dplyr::select(., Taxa, Aphia.ID) %>% 
  distinct() -> tx_chktrm

df_tx0 <- left_join(df_tx0, tx_chktrm, by="Aphia.ID")
df_tx0$Taxa.x <- df_tx0$Taxa.y;df_tx0$Taxa.y <- NULL
df_tx0 %>%
  rename(Taxa=Taxa.x) %>% 
  dplyr::select(.,-Abund_m3) %>% 
  group_by(across(c(!AbundanceRaw))) %>% 
  summarise(.,AbundanceRaw=sum(AbundanceRaw),.groups = "drop") %>%
  ungroup() %>% 
  as_tibble(.) -> df_tx
rm(tx_chk,tx_chk0,tx_chktrm)
toc(log=TRUE)

### WIMS chemical data
tic("Load WIMS data")
print("Load WIMS data")
# WIMS Extract based on:
# Materials = 2HZZ & 2IZZ; SMPT_TYPE = CD, CC, CE; Dates from 01/06/2022-present
# SMP_Code:
# 42100171, 42100174, 42100179, 45400826, 60510027, 73015085, 82510555,
# 82615055, 82615255, 88002837, 88007163, 88007172, 88025879, BE061099,
# E0001449, E0001450, E0004730, G0003532, G0003572, LC544405, LC560357,
# PTTR0026, WA560349, Y0004367, Y0017477, YC536426

df_wims0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,
                                                 # "/WIMS_Extract_WaterQuality_Zoop_Samples_240618.xlsx"),
                                                 "/WIMS_Extract_WaterQuality_Zoop_Samples_240930.xlsx"),
                                          sheet="allDat"))

### counts for pots 1&2 have been summed as to have those for pots 3&4
### removed whitespace from AphiaIDs & taxon names

### prep data ####
### format & widen WIMS data ####
df_wims <- df_wims0

df_wims$PRN <- df_wims$SAMP_SCHEDULE_SAMPLE_ID
df_wims %>%
  dplyr::mutate(det=paste0(DETE_DESC,"_",UNIT_SHORT_DESC)) %>% ##create new variable label
  dplyr::mutate(Result=ifelse(is.na(df_wims$MEAS_SIGN == "<"), MEAS_RESULT,
                              paste0(MEAS_SIGN,MEAS_RESULT))) %>% #add "<" to results
  dplyr::select(.,c(WIMS.Code,REGION,Biosys.ID,
                    SMPT_LONG_NAME, SAMP_SAMPLE_DATE, SAMP_SAMPLE_TIME,
                    SAMP_Notes,
                    PRN,det,Result)) %>% ##only keep variables of interest
  tidyr::pivot_wider(.,names_from=det, values_from = Result) -> df_wims_w###widen data

## write WIMS wide csv
# write.csv(df_wims_w,file=paste0(datfol,"processedData/WIMS_wide.csv"),row.names = FALSE)
toc(log=TRUE)

### remove odd data
#df_tx <- as_tibble(df_tx0) ### create new data (keep df0 as 'raw')
tic("Prep for export & export data");print("Prep for export & export data")

### convert dates
df_tx$sample.date <- as.Date(df_tx$sample.date, origin = "1899-12-30")

# Remove 100 µm data [OPTIONAL] ####
df_tx_100um <- df_tx %>% 
  filter(str_starts(Sample.comments,"100um"))

df_tx %>% 
  filter(!str_starts(Sample.comments,"100um")) -> df_tx

###widen data & fill NAs with 0s ####
df_tx %>% 
  dplyr::select(.,-c("Aphia.ID")) %>% 
  pivot_wider(names_from = "Taxa",values_from = "AbundanceRaw",#Abund_m3
              values_fill = 0) %>% 
  filter(.,!is.na(`Net.volume.sampled.(m3)`)) %>% ##OPTIONAL: remove 'empty' net volumes
  ungroup() -> df_tx_w

### join & save data ####  
dfw <- left_join(df_tx_w,df_wims_w,by="PRN")

### generate LONG version WITH zero values (for calculation of means)
df_tx_w %>% 
  pivot_longer(.,
               cols = !(Pot.Number:Category)
  ) -> df_tx_l
toc(log=TRUE)

##### look at mean-variance relationship ####
## throw away wims stuff
mn <- apply(df_tx_w %>% dplyr::select(.,-c(Pot.Number:Category)),2,mean)
v <- apply(df_tx_w %>% dplyr::select(.,-c(Pot.Number:Category)),2,var)
x <- data.frame(mn,v);rm(mn,v)

x %>% 
  ggplot(.,aes(x=log(mn),y=log(v)))+
  # geom_abline(slope = 1, linetype="dashed")+
  xlab("Mean")+ylab("Variance")+
  geom_point()
  
rm(df_tx_w)
###########
### GLLVMs ####
## create taxon data
tic("Model set up")
dfw %>% 
  dplyr::select(-c(1:21)) %>% ###remove metadata info
  dplyr::select_if(~ !is.numeric(.) || sum(.) !=0) %>% ### drop cols summing to 0
  dplyr::select_if(~ !is.numeric(.) || sum(. != 0) >= 100) %>%  # Drop numeric columns with <10 non-zero values
  filter(rowSums(across(where(is.numeric)))!=0) -> dftmp ###remove 'empty' rows



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
toc(log=TRUE)

### fit models ####
Y <- dftmp %>% dplyr::select(.,-c(WIMS.Code.y:"Zinc, Dissolved_ug/l"))
names(Y) <- vegan::make.cepnames(names(Y))
X <- df_wims_w_trim0_scale
runs <- 3 #set number of reruns for model fitting

### unconstrained model ####
#### Unconstrained NB distribution ####
tic("gllvm_uncon_offset_pois: Unconstrained NB with offset")
sDsn <- data.frame(Region = df_wims_w_trim0$Region)

set.seed(pi);m_lvm_0nb <- gllvm(y=Y,
                   family="negative.binomial",
                   offset = log(dfw$`Net.volume.sampled.(m3)`),#logging offset to match poisson link
                   studyDesign = sDsn, row.eff = ~(1|Region),
                   num.lv = 2,
                   starting.val="random",
                   n.init=runs,#re-run model to get best fit
                   trace=TRUE
                   )
saveRDS(m_lvm_0nb, file="figs/gllvm_logOffset_uncon_nb.Rdat") #12.647 mins
toc(log=TRUE)

m_lvm_0nb <- readRDS("figs/gllvm_logOffset_uncon_nb.Rdat")
par(mfrow=c(2,2));plot(m_lvm_0nb,which = 1:4);par(mfrow=c(1,1))
ordiplot.gllvm(m_lvm_0nb,biplot=TRUE)
gllvm::ordiplot(m_lvm_0nb,predict.region=TRUE)

### extract variables ###
LVs <- getLV(m_lvm_0nb)

unscaledLoadings <- coef(m_lvm_0nb, "loadings")
scaleLoadings <- coef(m_lvm_0nb, "sigma.lv")
sppLoadings <- unscaledLoadings%*%diag(scaleLoadings)
colnames(sppLoadings) <- c("LV1","LV2")
# sppLoadings <- as.data.frame(sppLoadings)
# sppLoadings$spp <- row.names(sppLoadings)

range(LVs[,1]);range(sppLoadings[,1])
range(LVs[,2]);range(sppLoadings[,2])

# ggplot(LVs,aes(x=LVs[,1],y=LVs[,2],col=df_wims_w_trim0_scale$Region))+
#   geom_point()+
#   geom_hline(yintercept = 0,linetype=2,col="grey")+
#   geom_vline(xintercept = 0,linetype=2,col="grey")+
#   xlim(-3.1,3.0)+
#   ylim(-2.5,2.73)+
#   geom_text(data=sppLoadings,aes(x=LV1,y=LV2,label=spp),col=4)+
#   xlab("LV1")+ylab("LV2")+
#   theme(legend.title = element_blank(),
#         axis.title = element_text(face=2))

# Calculate covariances
covariances <- getResidualCov(m_lvm_0nb)
correlations <- getResidualCor(m_lvm_0nb)
## visualise covariance
p1 <- ggplot(reshape2::melt(correlations))+
  geom_tile(aes(x=Var1, y=Var2, fill = value))+
  scale_fill_gradient2(low="red",high="blue",mid = "white")+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 1, size = 12, hjust = 1),
    axis.text.y = element_text(size = 12)) +
  ggplot2::coord_fixed()+xlab(NULL)+ylab(NULL)

do_svd <- svd(LVs)
rotation <- do_svd$v
scales <- sapply(1:ncol(LVs), function(q)sqrt(sum(LVs[,q]^2))*sqrt(sum(sppLoadings[,q]^2)))
newLVs <- apply(LVs,2,function(x)x/sqrt(sum(x^2))*scales^0.5)
newRotatedLVs <- newLVs%*%do_svd$v
newLoadings <- apply(sppLoadings,2,function(x)x/sqrt(sum(x^2))*scales^0.5)
newRotatedLoadings <- sppLoadings%*%do_svd$v

p2 <- ggplot()+
  geom_hline(yintercept = 0,linetype=2,col="grey")+
  geom_vline(xintercept = 0,linetype=2,col="grey")+
  geom_point(data=LVs, aes(y=LV2,x=LV1, col=df_wims_w_trim0_scale$Region))+
  geom_text(data=sppLoadings, aes(y = LV2, x = LV1, label = colnames(Y)),col=4)+
  ggthemes::theme_few()+coord_fixed()+
    theme(legend.title = element_blank(),
          axis.title = element_text(face=2))
  guides(col="none")

p1|p2

#### Unconstrained unimodal NB distribution ####
# following guidance in:
# https://github.com/BertvanderVeen/GLLVM-workshop/blob/main/Practicals/5Practical.html
tic("gllvm_uncon_offset_pois: Unconstrained unimodal NB with offset")
sDsn <- data.frame(Region = df_wims_w_trim0$Region)

set.seed(pi);m_lvm_0nbUnim <- gllvm(y=Y,
                                family="negative.binomial",
                                offset = log(dfw$`Net.volume.sampled.(m3)`),#logging offset to match poisson link
                                num.lv = 2,
                                n.init=runs,#re-run model to get best fit
                                trace=TRUE,
                                quadratic = TRUE
                                )
saveRDS(m_lvm_0nbUnim, file="figs/gllvm_logOffset_uncon_nb_Unim.Rdat")
toc(log=TRUE)

ordiplot(m_lvm_0nbUnim,biplot=TRUE)

## make turnover predictions (NOT WORKING!!!!!)
# LVs = getLV(m_lvm_0nbUnim)
# newLV = cbind(LV1 = seq(min(LVs[,1]), max(LVs[,1]), length.out=1000), LV2 = 0)
# preds <- predict(m_lvm_0nbUnim, type = "response", newLV = newLV)
# plot(NA, ylim = range(preds), xlim = c(range(getLV(m_lvm_0nbUnim))), ylab  = "Predicted response", xlab = "LV1")
# segments(x0=optima(m_lvm_0nbUnim, sd.errors = FALSE)[,1],x1 = optima(m_lvm_0nbUnim, sd.errors = FALSE)[,1], y0 = rep(0, ncol(model1$y)), y1 = apply(preds,2,max), col = "red", lty = "dashed", lwd = 2)
# rug(getLV(m_lvm_0nbUnim)[,1])
# sapply(1:ncol(m_lvm_0nbUnim$y), function(j)lines(sort(newLV[,1]), preds[order(newLV[,1]),j], lwd = 2))
# 
# newLV = cbind(LV1 = 0, LV2 =  seq(min(LVs[,2]), max(LVs[,2]), length.out=1000))
# preds <- predict(m_lvm_0nbUnim, type = "response", newLV = newLV)
# plot(NA, ylim = range(preds), xlim = c(range(getLV(model1))), ylab  = "Predicted response", xlab = "LV2")
# segments(x0=optima(model1, sd.errors = FALSE)[,2],x1 = optima(model1, sd.errors = FALSE)[,2], y0 = rep(0, ncol(model1$y)), y1 = apply(preds,2,max), col = "red", lty = "dashed", lwd = 2)
# rug(getLV(model1)[,2])
# sapply(1:ncol(model1$y), function(j)lines(sort(newLV[,2]), preds[order(newLV[,2]),j], lwd = 2))

##### create DHARMa object ####
# simulate the model x1000
sim.uncon <- do.call("cbind", replicate(1000,
                                        c(as.matrix(
                                          gllvm::simulate(m_lvm_0nb,conditional=TRUE))),
                                        simplify = FALSE))
preds <- c(predict(m_lvm_0nb))
obs <- c(as.matrix(Y))
dharma <- DHARMa::createDHARMa(
  simulatedResponse = sim.uncon,
  observedResponse = obs,
  integerResponse = TRUE,
  fittedPredictedResponse = preds)
plot(dharma)
toc(log=TRUE)

### constrained model ####
#### NB distribution #####
tic("gllvm_offset_nh4SalChlaDinDepPo4Reg_NB_Scaled: Constrained & scaled with offset")
sDsn <- data.frame(Region = df_wims_w_trim0$Region)

# m_lvm_1nb <- gllvm(y=Y,
#                     X=X,#scaled
#                     formula = ~ nh4 + sal_ppt + chla + din + depth + po4 + tempC,
#                     family="negative.binomial",
#                     offset = log(dfw$`Net.volume.sampled.(m3)`), #logging offset to match NB link
#                     studyDesign = sDsn, row.eff = ~(1|Region),
#                     num.lv = 2,
#                    n.init=5, trace=TRUE
#                    )

# updated model with lv.formula
m_lvm_1nb <- gllvm(y=Y,
                   X=X,
                   formula = ~ (tempC+sal_ppt+din+chla),
                   #lv.formula = ~ tempC+sal_ppt+din+chla,#latent variables
                   family="negative.binomial",
                   offset = log(dfw$`Net.volume.sampled.(m3)`),
                   studyDesign = sDsn,
                   row.eff = ~(1|Region),
                   # num.lv.c = 2,#latent variables in current ordination
                   # num.RR = 2,#latent variables in constrained ordination
                   num.lv=2,
                   n.init=runs,#re-run model to get best fit
                   trace=TRUE,
                   randomB = "LV"
                   )
toc(log=TRUE)
## NO PREDICTOR EFFECTS IN CURRENT MODEL ##


saveRDS(m_lvm_1nb, file="figs/gllvm_logOffset_nh4SalChlaDinDepPo4Reg_nb_Scaled_constrained.Rdat") # scaled
m_lvm_1nb <- readRDS("figs/gllvm_logOffset_nh4SalChlaDinDepPo4Reg_nb_Scaled_constrained.Rdat") #scaled
hist(m_lvm_1nb$TMBfn$gr())
summary(m_lvm_1nb)
par(mfrow=c(2,2));plot(m_lvm_1nb,which = 1:4);par(mfrow=c(1,1))
ordiplot.gllvm(m_lvm_1nb,biplot=TRUE,symbols = TRUE)
#summary(m_lvm_1nb, by="all")
toc(log=TRUE)

pdf(file = "figs/m_lvm_1nb_tx_all.pdf",width=16,height=8)
coefplot(m_lvm_1nb,cex.ylab = 0.3,
         order=FALSE)
dev.off()

tic("Extract & plot model estimated")
# extract 'significant' model/species terms from model for plotting ####
ci_mod_all <- as.data.frame(confint(m_lvm_1nb))
ci_mod_var <- ci_mod_all[grep("^X", rownames(ci_mod_all)), ]
rownames(ci_mod_var) <- substring(rownames(ci_mod_var), 7)
ci_mod_var$varTrt <- rownames(ci_mod_var)

sigterms_all <- summary(m_lvm_1nb)
sigterms_all <- as.data.frame(sigterms_all$Coef.tableX)
sigterms_all$variable <- sub(":.*","",row.names(sigterms_all))
sigterms_all$trt <- sub(".*:","",row.names(sigterms_all))
sigterms_all$varTrt <- rownames(sigterms_all)
sigterms_all <- left_join(sigterms_all, ci_mod_var, by = "varTrt")
sigterms_all$sig <- sigterms_all$`2.5 %`*sigterms_all$`97.5 %`>0

sigterms_sig <- sigterms_all[sigterms_all$`Pr(>|z|)`>0.05,]

### plot! ####
## recreate coefplot
ggplot(sigterms_all[sigterms_all$variable=="nh4",],
       aes(x=Estimate, y=trt,
           xmin=`2.5 %`,
           xmax=`97.5 %`,
           colour=sig))+
  geom_vline(xintercept = 0)+
  geom_errorbar()+
  geom_point()+
  scale_y_discrete(limits = rev(levels(as.factor(sigterms_all$trt))))+
  scale_colour_manual(values = c("grey","black"))+
  guides(colour="none")

#############
plot_list <- list()
sigterms_all$variable <- as.factor(sigterms_all$variable)
ntrt <- length(unique(sigterms_all$trt))-.5
sigterms_all %>% 
  mutate(flag = case_when(
    !sig ~ "NONE",
    sig & Estimate >0 ~ "pos",
    sig & Estimate <0 ~ "neg"
  )) -> sigterms_all

# Iterate over each level of the factor 'trt'
for (level in levels(sigterms_all$variable)) {
  # Subset the data for the current level
  subset_data <- sigterms_all[sigterms_all$variable == level, ]
  
  # Calculate the mean value of 'Estimate' for the current level
  mean_estimate <- mean(subset_data$Estimate)
  
  # Determine the color of the vertical line based on the mean estimate
  line_color <- ifelse(mean_estimate > 0, "blue", ifelse(mean_estimate < 0, "red", "grey"))
  
  # Create a plot for the current level
  current_plot <- ggplot(subset_data,
                         aes(x=Estimate, y=trt,
                             xmin=`2.5 %`,
                             xmax=`97.5 %`,
                             colour=flag,
                             fill=flag)) +
    geom_hline(yintercept = seq(1.5,ntrt,by=1),col="lightgrey",lty=3)+
    geom_vline(xintercept = 0)+
    # Add vertical line for mean estimate
    geom_vline(xintercept = mean_estimate, color = line_color,linetype="dashed") +
    geom_linerange()+
    labs(title = paste0(level))+
    geom_point(shape=21) +
    scale_y_discrete(limits = rev(levels(as.factor(sigterms_all$trt))))+
    scale_colour_manual(values = c("red",#negative
                                        "grey",#null
                                        "blue"#positive
    ))+
    scale_fill_manual(values = c("red",#negative
                                      "white",#null
                                      "blue"#positive
    ))+
    guides(colour="none",
           fill="none")+
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust=0.5),
          axis.text.y = element_text(size=3))
  
  # Add the current plot to the list
  plot_list[[as.character(level)]] <- current_plot
  
  # Iterate over each plot in the list
  for (i in seq_along(plot_list)) {
    # If it's not the first plot, hide y-axis labels
    if (i > 1) {
      plot_list[[i]] <- plot_list[[i]] + theme(axis.text.y = element_blank())
    }
  }
}

# Combine all the individual plots into a single plot
final_plot <- wrap_plots(plotlist = plot_list, ncol = nlevels(sigterms_all$variable))+  # Adjust the number of columns as needed
  plot_annotation(title="Caterpillar plot of generalised linear latent variable model outputs",
                  subtitle = bquote("Point estimates & 95% confidence intervals of species-specific coefficients "~italic(hat(beta)[j])~". Based on zooplankton taxon abundance data and scaled water quality parameters"), #scaled
                  caption = paste0("Colours indicate taxon 95% confidence intervals which do (grey) or do not (red/blue) include zero",
                                   "\nMissing covariate values replaced with covariate means prior to model run",
                                   "\nModel call: ~",as.character(m_lvm_1nb$formula)[2],
                                   "\nFamily: ",as.character(m_lvm_1nb$family),". ",
                                   "\nModel offset: sampling volume (m3)",
                                   "Random row effects: ",as.character(m_lvm_1nb$call)[7]),
                  theme = theme(plot.title = element_text(size = 16, face="bold")))

pdf(file = "figs/coef_tax_all_unordered_v2_scaled_nb_offset.pdf",width=16,height=8)
print(final_plot)
dev.off()

### compare models w/ w/out env parameters####
# based on:
#https://jenniniku.github.io/gllvm/articles/vignette1.html#studying-co-occurrence-patterns
# getResidualCov function can be used to quantify the amount of variation in the
# data that can be explained by environmental variables.
# Specifically, if we use the trace of the residual covariance matrix Σ as a
# measure of unexplained variation, then we can compare this quantity before
# and after environmental variables are included in the model. The ratio of
# traces suggests that environmental variables explain approximately 40% of
# the (co)variation in ant species abundances.

## no adjust (for poisson?)
rcov0 <- getResidualCov(m_lvm_0nb, adjust = 0) # 'null' model
rcov1 <- getResidualCov(m_lvm_1nb, adjust = 0) # model with env variables

rcov0 <- getResidualCov(m_lvm_0nb, adjust = 1) # 'null' model
rcov1 <- getResidualCov(m_lvm_1nb, adjust = 1) # model with env variables
rcov0$trace; rcov1$trace
100 - (rcov1$trace / rcov0$trace*100)
AIC(m_lvm_0nb,m_lvm_1nb)

ordiplot(m_lvm_1nb, biplot=TRUE, symbols = TRUE)

# Environmental correlation
# modelmat <- model.matrix(m_lvm_3$formula, data = df_wims_w_trim0) %>% 
#   as.matrix %>% 
#   {.[,-1]} # Remove intercept
# linpred <- tcrossprod(modelmat, m_lvm_3$params$Xcoef)
# envircor <- cov2cor(cov(linpred))


### to do:
### look at functional groups(lifeforms)
# master list saved in:
# \\prodds.ntnl\Shared\AN\KFH\Groups\N_Marine\07 Training & Reference Documents\A&R Technical Guidance\Traits, Lifeforms etc\Plankton Lifeform Extraction Tool
# consider reproducing ordination in 3 dimensions using rgl (see:
# https://riffomonas.org/code_club/2021-03-24-rgl)

# PRIORITY : REPRODUCE CODE ####
## Currently untidy and seems to produce 'issues'
## Errors alluding to "non-numeric argument to binary operator"
toc()
