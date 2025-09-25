# an.GLLVM_revisit_Indiv_WB.R ####
# an.GLLVM_revisit.R ####
# Conducting GLLVM modelling on carbon content values

# 00 set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "tictoc","gllvm", "lubridate")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

## set universals ####
tictoc::tic.clearlog();tic("Set universals");print("set universals")
source("R/set_meta.R")
toc(log=TRUE)

## load data ####
# ### Biological data ####

tic("Load mean carbon content values")
source("R/an_EA_taxa_gllvm_no_offset.R")
toc(log=TRUE)

#############
# GLLVMs ####
#############
## data wrangling ####
tic("gllvm: data wrangling")

## pull out only Solent values
wbdy <- "Solent"
df_carb0 <- df_carb;rm(df_carb)
use <- df_carb0$df_carb_meta$WB == wbdy
meta <- df_carb0$df_carb_meta[use,]
lf <- df_carb0$df_carb_lf[use,] %>% 
  dplyr::select(where(~ any(.x != 0)))
wims <- df_carb0$df_carb_wims[use,]

df_carb <- list(
  df_carb_meta = meta,
  df_carb_lf = lf,
  df_carb_wims = wims
  )
rm(use,meta,lf,wims)

### set number of iterations
runs <- 5

## define X (environmental parameters) & Y (biological responses)
Y <- df_carb[[2]] %>% dplyr::select(.,-c(PRN))
X <- df_carb[[3]] %>% dplyr::select(.,-c(PRN))

X0 <- X
Y0 <- Y

# choose interesting environmental variables
keep <- c(
  # "Net.volume.sampled.(m3)",
  "Ammoniacal Nitrogen, Filtered as N_mg/l",
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
kp <- names(X) %in% keep
X <- X[,kp]
rm(kp, keep)

### replace LESS THAN values with value*0.5
# define function to replace "less than" values by  half their value
replace_values <- function(x) {
  if_else(str_detect(x, "^<"), as.numeric(sub("^<", "", x))/2, as.numeric(x))
}

# amend < values in wims data
X %>% 
  mutate_all(.,replace_values) -> X

##append region and wb names
X$WB <- df_carb[[1]]$WB_lb
X$Region <- df_carb[[1]]$Region
X$Site.code <- df_carb[[1]]$BIOSYS.Code

#OPTIONAL: remove lifeforms which only appear n>=4 times ####
n <- 5 # 7 terms in 'full' X-model
Y <- Y[,colSums(ifelse(Y==0,0,1))>n]

## replace NA values with mean values for respective column.
## Leave non-numeric cols unchanged
X %>% 
  mutate(across(where(is.numeric),
                ~replace(.,
                         is.na(.),
                         mean(.,
                              na.rm = TRUE)
                )
  )
  ) -> X

# # OPTIONAL: subsample data for model building ####
# num_rows <- nrow(X)
# num_true <- 5000
# 
# sample_vec <- rep(FALSE, num_rows)
# sample_vec[sample(num_rows, num_true)] <- TRUE
# rm(num_rows, num_true)
# Y <- Y[sample_vec,]
# X <- X[sample_vec,]
# ## identify and remove 'empty' abundance rows after subsampling
# emptY <- rowSums(Y) != 0
# X <- X[emptY,]
# Y <- Y[emptY,]

## rename colums
X <- X %>% 
  rename(
    wb="WB",
    region="Region",
    # net_vol_m3 = "Net.volume.sampled.(m3)",
    nh4 = "Ammoniacal Nitrogen, Filtered as N_mg/l",
    chla = "Chlorophyll : Acetone Extract_ug/l",
    ngr_e = "NGR : Easting_NGR",
    ngr_n = "NGR : Northing_NGR",
    no3 = "Nitrate, Filtered as N_mg/l",
    no2 = "Nitrite, Filtered as N_mg/l",
    din = "Nitrogen, Dissolved Inorganic : as N_mg/l",
    ton = "Nitrogen, Total Oxidised, Filtered as N_mg/l",
    po4 = "Orthophosphate, Filtered as P_mg/l",
    o2_dis_mgl = "Oxygen, Dissolved as O2_mg/l",
    o2_dis_sat = "Oxygen, Dissolved, % Saturation_%",
    sal_ppt = "Salinity : In Situ_ppt",
    si = "Silicate, Filtered as SiO2_mg/l",
    tempC = "Temperature of Water_CEL",
    turb = "Turbidity : In Situ_FTU",
    depth = "Water Depth_m"
  )

## shorten Region labels
X$rgn <- ifelse(X$region == "Southern", "Sth",
                ifelse(X$region == "NEast", "NE",
                       ifelse(X$region == "NWest", "NW",
                              ifelse(X$region == "Anglian", "An",
                                     ifelse(X$region == "Thames", "Th",
                                            ifelse(X$region == "SWest", "SW",NA)
                                     )))))

### create scaled version for comparison of effects on model
X %>% 
  mutate_if(is.numeric,scale) -> X_scaled
offset_m3 <- df_carb[[1]]$netVol_use
toc(log=TRUE)

tic("Fit Unconstrained model")
## Fit GLLVM ####
## unconstrained ####
### Tweedie ####
sDsn <- data.frame(Site.code = X$Site.code)
m_lvm_0 <- gllvm(as.matrix(Y), # unconstrained model
                 studyDesign = sDsn,
                 row.eff = ~(1|Site.code),
                 family = "tweedie",
                 starting.val="res",
                 n.init = runs, #re-run model to get best fit
                 trace=TRUE,
                 seed = 123,
                 num.lv = 2
)

saveRDS(m_lvm_0, file="figs/2412dd/gllvm_SOLENT_carbon_uncon_WB_tweedie.Rdat")
toc(log=TRUE)

tic("Generate summary plots for unconstrained model")
m_lvm_0 <- readRDS("figs/2412dd/gllvm_SOLENT_carbon_uncon_WB_tweedie.Rdat")
par(mfrow=c(2,2));plot(m_lvm_0, which=1:4);par(mfrow=c(1,1))
gllvm::ordiplot.gllvm(m_lvm_0,biplot = TRUE,symbols=TRUE)
cr <- getResidualCor(m_lvm_0)

pdf(file = "figs/2412dd/m_lvm_0_carb_LF_corrplot_SOLENT.pdf",width=14,height=14)
corrplot::corrplot(cr, diag = FALSE, type = "lower", method = "square",
                   tl.srt = 25)
dev.off();rm(cr)
toc(log=TRUE)

## Constrained ####
### Negbin ####
tic("Fit Constrained Negative binomial model")
# sDsn <- data.frame(region = X$region)
sDsn <- data.frame(Site.code = factor(X$Site.code))
m_lvm_4 <- gllvm::gllvm(y=as.matrix(Y), # model with environmental parameters
                 X=X_scaled, #scaled
                 formula = ~ nh4 + sal_ppt + chla + din + depth + po4 + tempC,
                 studyDesign = sDsn,
                 row.eff = ~(1|Site.code),
                 family="tweedie",
                 starting.val="res",
                 n.init=runs,#re-run model to get best fit
                 trace=TRUE,
                 num.lv = 2
                 )

saveRDS(m_lvm_4, file="figs/2412dd/gllvm_SOLENT_carbon_Con_WB_tweedie.Rdat")
toc(log=TRUE)

tic("Model summaries & comparisons")
m_lvm_4 <- readRDS("figs/2412dd/gllvm_SOLENT_carbon_Con_WB_tweedie.Rdat")#scaled
gllvm::ordiplot.gllvm(m_lvm_4,biplot = TRUE,symbols=TRUE)
cr <- getResidualCor(m_lvm_4)

pdf(file = "figs/2412dd/m_lvm_4_carb_corrplot_SOLENT.pdf",width=14,height=14)
corrplot::corrplot(cr, diag = FALSE, type = "lower", method = "square",
                   tl.srt = 25)
dev.off()

AIC(m_lvm_0,m_lvm_4)
anova(m_lvm_0,m_lvm_4)

## GLLVM plots ####
pdf(file = "figs/2412dd/coef_trt_all_unordered_SOLENT.pdf",width=16,height=8)
coefplot(m_lvm_4,cex.ylab = 0.7,
         order=FALSE)
dev.off()

for(i in 1:ncol(m_lvm_4$X.design)){
  pdf(file = paste0("figs/2412dd/coef_nb_trt_",i,"_SOLENT.pdf"),width = 7, height = 14)
  coefplot(m_lvm_4,mfrow = c(1,1),which.Xcoef = i, cex.ylab = 0.6,
           main=colnames(m_lvm_4$X.design)[i])
  dev.off()
}

### extract model terms for plotting ####
mod_coefs <- as.data.frame(m_lvm_4$params$Xcoef)
mod_coefs$LF <- row.names(mod_coefs)
mod_coefs <- mod_coefs %>% 
  pivot_longer(.,cols=!LF,names_to = "coefficient", values_to = "Estimate")
sdXcoef <- as.data.frame(m_lvm_4$sd$Xcoef[,  drop = FALSE])
sdXcoef$LF <- row.names(sdXcoef)
sdXcoef <- sdXcoef %>% 
  pivot_longer(.,cols=!LF,names_to = "coefficient", values_to = "sd")

mod_coefs <- dplyr::left_join(mod_coefs,sdXcoef, by=c("LF","coefficient"))
mod_coefs <- mod_coefs %>% 
  mutate(.,lower = Estimate-1.96*sd,
         upper = Estimate+1.96*sd) %>% 
  mutate(.,varTrt=paste0(coefficient,"_",LF))

sigterms_all <- summary(m_lvm_4)
sigterms_all <- as.data.frame(sigterms_all$Coef.tableX)
sigterms_all$variable <- sub(":.*","",row.names(sigterms_all))
sigterms_all$trt <- sub(".*:","",row.names(sigterms_all))
sigterms_all$varTrt <- rownames(sigterms_all)
sigterms_all$varTrt <- gsub(":","_",sigterms_all$varTrt)

mod_coefs <- dplyr::left_join(mod_coefs,
                              sigterms_all[,c("varTrt","Std. Error","z value","Pr(>|z|)")],
                              by="varTrt")
rm(sigterms_all,sdXcoef)
## create flag if CI doesn't cross 0
mod_coefs$sig <- mod_coefs$lower*mod_coefs$upper>0

### plot! ####
plot_list <- list()
mod_coefs$coefficient <- as.factor(mod_coefs$coefficient)
ntrt <- length(unique(mod_coefs$LF))-.5
mod_coefs %>% 
  mutate(flag = case_when(
    !sig ~ "NONE",
    sig & Estimate >0 ~ "pos",
    sig & Estimate <0 ~ "neg"
  )) %>% 
  mutate(clr = case_when(
    !sig ~ "grey",
    sig & Estimate >0 ~ "blue",
    sig & Estimate <0 ~ "red"
  )) -> mod_coefs


# Iterate over each level of the factor 'LF'
for (level in levels(mod_coefs$coefficient)) {
  # Subset the data for the current level
  subset_data <- mod_coefs[mod_coefs$coefficient == level, ]
  #subset_data <- mod_coefs[mod_coefs$coefficient == "sal_ppt", ]
  
  # Calculate the mean value of 'Estimate' for the current level
  mean_estimate <- mean(subset_data$Estimate)
  
  # Determine the color of the vertical line based on the mean estimate
  line_color <- ifelse(mean_estimate > 0, "blue", ifelse(mean_estimate < 0, "red", "grey"))
  
  # Create a plot for the current level
  current_plot <- ggplot(subset_data,
                         aes(x=Estimate, y=LF,
                             xmin=lower,
                             xmax=upper,
                             colour=clr,
                             fill=clr)) +
    geom_hline(yintercept = seq(1.5,ntrt,by=1),col="lightgrey",lty=3)+
    geom_vline(xintercept = 0)+
    # Add vertical line for mean estimate
    geom_vline(xintercept = mean_estimate, color = line_color,linetype="dashed") +
    geom_linerange()+
    labs(title = paste0(level))+
    geom_point() +
    scale_y_discrete(limits = rev(levels(as.factor(mod_coefs$LF))))+
    scale_colour_identity()+
    guides(colour="none",
           fill="none")+
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust=0.5))
  
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

rcov0 <- getResidualCov(m_lvm_0, adjust = 1) # 'null' model
rcov1 <- getResidualCov(m_lvm_4, adjust = 1) # model with env variables #REGION
btwn <- 100 - (rcov1$trace / rcov0$trace*100)
print(paste0("Including environmental parameters in the model explains ",round(btwn,2),"% of the (co)variation in zooplankton abundances"))
AIC(m_lvm_0,m_lvm_4)

######################################
# Combine all the individual plots into a single plot
(final_plot <- patchwork::wrap_plots(plotlist = plot_list,
                                     ncol = nlevels(mod_coefs$coefficient))+  # Adjust the number of columns as needed
   patchwork::plot_annotation(title=paste0("Caterpillar plot of generalised linear latent variable model outputs\n",wbdy," water body"),
                              subtitle = paste0("Including environmental variables explains ",round(btwn,2),"% of the (co)variation in lifeform abundances compared to the null (lifeforms-only) model"),
                              caption = paste0("Colours indicate lifeform 95% confidence intervals which do (grey) or do not (red/blue) include zero. ",
                                               "Dashed vertical lines indicate mean point estimate values\n","Lifeforms recorded in ",n," or fewer samples removed from data prior to model estimations","\n",
                                               "Model call: ~",as.character(m_lvm_4$formula)[2],". ",
                                               "Distribution family: ",as.character(m_lvm_4$family),".",
                                               "\nRandom row effects: ",as.character(m_lvm_4$call)[8],"; number of model iterations = ",m_lvm_4$n.init,".",
                                               "\nSamples gathered between ",format(min(dfl$sample.date), "%d/%m/%Y")," & ",format(max(dfl$sample.date), "%d/%m/%Y")),
                              theme = theme(plot.title = element_text(size = 16, face="bold"))))

pdf(file = "figs/2412dd/coef_trt_all_unordered_v2_scaled_SOLENT.pdf",width=16,height=8) #scaled
print(final_plot)
dev.off()
toc(log=TRUE)
