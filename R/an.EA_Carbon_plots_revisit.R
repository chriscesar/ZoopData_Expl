# an.EA_Carbon_plots_revisit.R ####
# Re-examiniation of carbon contents lifeforms data


# 00 set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "tictoc","gllvm","patchwork", "lubridate","ggpubr")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

## set universals ####
tictoc::tic.clearlog();tic("Set universals");print("set universals")
source("R/set_meta.R")
toc(log=TRUE)

tic("Load zoop data")
source("R/imp.load_data_all_zoops.R")
toc(log=TRUE)
rm(dfl0)

# total carbon trend ####
## Medians ####
dfl %>% #names()
  dplyr::select(.,Pot.Number,WB_lb,sample.date,
                DisplayName:md_carbTot_m3) %>% #names()#View()
  ## remove unneeded
  dplyr::select(.,-c(
    DisplayName,
    AbundanceRaw, Abund_m3,
    copNonCop,mnlongMaxAxis_mm, mdlongMaxAxis_mm,
    mnCPerIndiv_ug, mdCPerIndiv_ug, mn_carbTot_raw,
    md_carbTot_raw,mn_carbTot_m3)) %>% #View()
  group_by(across(c(!md_carbTot_m3))) %>% 
  summarise(sum=sum(md_carbTot_m3, na.rm = TRUE), .groups = "drop") %>% 
    ggplot(., aes(x=sample.date, y=log(sum)))+
    geom_point()+
    geom_smooth(se=FALSE)+
    facet_wrap(.~WB_lb)+#, scale="free_y")+
    ylim(0,NA)+
    labs(title = "Trend in total carbon content within zooplankton assemblages by date in EA water bodies",
         y = expression(bold("Log total carbon content (ug m"^-3~")")),
         x = "Date",
         caption=paste0("Samples gathered between ",format(min(dfl$sample.date), "%d/%m/%Y")," & ",format(max(dfl$sample.date), "%d/%m/%Y"),
                        "\nBlue lines inidate loess smooths.",
                        "\nCarbon content values are based on *median* estimates of carbon contents per taxon")) +
    theme(legend.position = "none",
          axis.title = element_text(face=2),
          strip.text = element_text(face=2)) -> pl

ggsave(plot = pl, filename = "figs/2412dd/logtotCByDateByWB_Fixed_Y_median.pdf",
       width = 12,height = 8,units = "in"); rm(pl)

## Means ####
dfl %>% #names()
  dplyr::select(.,Pot.Number,WB_lb,sample.date,
                DisplayName:md_carbTot_m3) %>% #names()#View()
  ## remove unneeded
  dplyr::select(.,-c(
    DisplayName,
    AbundanceRaw, Abund_m3,
    copNonCop,mnlongMaxAxis_mm, mdlongMaxAxis_mm,
    mnCPerIndiv_ug, mdCPerIndiv_ug, mn_carbTot_raw,
    md_carbTot_raw,md_carbTot_m3)) %>% #names()
  group_by(across(c(!mn_carbTot_m3))) %>% 
  summarise(sum=sum(mn_carbTot_m3, na.rm = TRUE), .groups = "drop") %>% 
  ggplot(., aes(x=sample.date, y=log(sum)))+
  geom_point()+
  geom_smooth(se=FALSE)+
  facet_wrap(.~WB_lb)+#, scale="free_y")+
  ylim(0,NA)+
  labs(title = "Trend in total carbon content within zooplankton assemblages by date in EA water bodies",
       y = expression(bold("Log total carbon content (ug m"^-3~")")),
       x = "Date",
       caption=paste0("Samples gathered between ",format(min(dfl$sample.date), "%d/%m/%Y")," & ",format(max(dfl$sample.date), "%d/%m/%Y"),
                      "\nBlue lines inidate loess smooths.",
                      "\nCarbon content values are based on *mean* estimates of carbon contents per taxon")) +
  theme(legend.position = "none",
        axis.title = element_text(face=2),
        strip.text = element_text(face=2)) -> pl

ggsave(plot = pl, filename = "figs/2412dd/logtotCByDateByWB_Fixed_Y_mean.pdf",
       width = 12,height = 8,units = "in"); rm(pl)

# Carbon by WB ####
## Medians ####
dfl %>% #names()
  dplyr::select(.,Pot.Number,WB_lb,sample.date,Region,
                DisplayName:md_carbTot_m3) %>% #names()#View()
  ## remove unneeded
  dplyr::select(.,-c(
    DisplayName,
    AbundanceRaw, Abund_m3,
    copNonCop,mnlongMaxAxis_mm, mdlongMaxAxis_mm,
    mnCPerIndiv_ug, mdCPerIndiv_ug, mn_carbTot_raw,
    md_carbTot_raw,mn_carbTot_m3)) %>% #View()
  group_by(across(-md_carbTot_m3)) %>% 
  summarise(sum=sum(md_carbTot_m3, na.rm = TRUE), .groups = "drop") %>% 
  mutate(mean_log_sum = mean(log(sum), na.rm = TRUE)) %>% 
  ggplot(.,aes(x=WB_lb, y=log(sum), colour=Region))+
  geom_hline(aes(yintercept = mean_log_sum),col=2,linetype="dashed")+
  geom_boxplot(varwidth = TRUE,outlier.shape = NA)+
  geom_point(aes(group=Region), position=position_jitterdodge(),alpha=0.3)+
  scale_colour_manual(values = cbPalette2)+
  labs(title = "Total carbon content in zooplankton assemblages by EA water body",
       y="log(total carbon)",
       caption = paste0("Samples gathered between ",format(min(dfl$sample.date), "%d/%m/%Y")," & ",format(max(dfl$sample.date), "%d/%m/%Y"),
         "\nDashed line indicates global mean log carbon content across all water bodies",
                        "\nBox widths are proportional to the number of observations",
                        "\nCarbon content values are based on *median* estimates of carbon contents per taxon"))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face=2),
        axis.text.x = element_text(face=2)) -> pl
ggsave(plot = pl, filename = "figs/2412dd/carbonByWB_median.pdf",
       width = 20,height = 12,units = "in");rm(pl)

## Means ####
dfl %>% #names()
  dplyr::select(.,Pot.Number,WB_lb,sample.date,Region,
                DisplayName:md_carbTot_m3) %>% #names()#View()
  ## remove unneeded
  dplyr::select(.,-c(
    DisplayName,
    AbundanceRaw, Abund_m3,
    copNonCop,mnlongMaxAxis_mm, mdlongMaxAxis_mm,
    mnCPerIndiv_ug, mdCPerIndiv_ug, mn_carbTot_raw,
    md_carbTot_raw,md_carbTot_m3)) %>% #View()
  group_by(across(-mn_carbTot_m3)) %>% 
  summarise(sum=sum(mn_carbTot_m3, na.rm = TRUE), .groups = "drop") %>% 
  mutate(mean_log_sum = mean(log(sum), na.rm = TRUE)) %>% 
  ggplot(.,aes(x=WB_lb, y=log(sum), colour=Region))+
  geom_hline(aes(yintercept = mean_log_sum),col=2,linetype="dashed")+
  geom_boxplot(varwidth = TRUE,outlier.shape = NA)+
  geom_point(aes(group=Region), position=position_jitterdodge(),alpha=0.3)+
  scale_colour_manual(values = cbPalette2)+
  labs(title = "Total carbon content in zooplankton assemblages by EA water body",
       y="log(total carbon)",
       caption = paste0("Samples gathered between ",format(min(dfl$sample.date), "%d/%m/%Y")," & ",format(max(dfl$sample.date), "%d/%m/%Y"),
                        "\nDashed line indicates global mean log carbon content across all water bodies",
                        "\nBox widths are proportional to the number of observations",
                        "\nCarbon content values are based on *mean* estimates of carbon contents per taxon"))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face=2),
        axis.text.x = element_text(face=2)) -> pl
ggsave(plot = pl, filename = "figs/2412dd/carbonByWB_mean.pdf",
       width = 20,height = 12,units = "in");rm(pl)

# Carbon by zoop type by WB ####
## Median ####
dfl %>% #names(.)
  mutate(.,LF03 = if_else(Copepod == "Y", paste0("Cop_",CopSize),ZooType)) %>%
  mutate(.,LF03 = if_else(is.na(LF03) | LF03 == "NYA", LF02, LF03)) %>%
  dplyr::select(.,c(Pot.Number, sample.date, yday, BIOSYS.Code, WIMS.Code,
                    LF03, PRN, 
                    md_carbTot_m3, Region,WB_lb)) %>% 
  mutate(md_carbTot_m3 = replace_na(md_carbTot_m3, 0)) %>% 
  filter(.,md_carbTot_m3 != 0) %>% 
  group_by(across(c(!md_carbTot_m3))) %>% 
  summarise(md_carbTot_m3 = sum(md_carbTot_m3),.groups = "drop") %>%
  ungroup(.) %>%
  ggplot(., aes(x=LF03, y=log(md_carbTot_m3), col=Region))+
  geom_boxplot(varwidth = FALSE,outlier.shape = NA)+
  coord_flip()+
  facet_wrap(.~WB_lb)+
  # facet_wrap(.~WB_lb)+
  geom_point(aes(group=LF03), position=position_jitter(width = 0.2),alpha=0.2)+
  scale_colour_manual(values = cbPalette2)+
  scale_x_discrete(limits=rev)+
  labs(
    title = "Total carbon content by zooplankton type by EA water body",
    y="log(Total carbon)",
    caption = paste0("Samples gathered between ",format(min(dfl$sample.date), "%d/%m/%Y")," & ",format(max(dfl$sample.date), "%d/%m/%Y"),
                     "\nCarbon content values are based on *median* estimates of carbon contents per taxon"))+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        strip.text = element_text(face=2),
        axis.title.x = element_text(face=2),
        axis.text.x = element_text(face=2),
        axis.text.y = element_text(face=2,
                                   size = 7)) -> pl

ggsave(plot = pl, filename = "figs/2412dd/carbonByZooTypeWB_Median.pdf",
       width = 20,height = 12,units = "in"); rm(pl)

## Mean ####
dfl %>% #names(.)
  mutate(.,LF03 = if_else(Copepod == "Y", paste0("Cop_",CopSize),ZooType)) %>%
  mutate(.,LF03 = if_else(is.na(LF03) | LF03 == "NYA", LF02, LF03)) %>%
  dplyr::select(.,c(Pot.Number, sample.date, yday, BIOSYS.Code, WIMS.Code,
                    LF03, PRN, 
                    mn_carbTot_m3, Region,WB_lb)) %>% 
  mutate(mn_carbTot_m3 = replace_na(mn_carbTot_m3, 0)) %>% 
  filter(.,mn_carbTot_m3 != 0) %>% 
  group_by(across(c(!mn_carbTot_m3))) %>% 
  summarise(mn_carbTot_m3 = sum(mn_carbTot_m3),.groups = "drop") %>%
  ungroup(.) %>%
  ggplot(., aes(x=LF03, y=log(mn_carbTot_m3), col=Region))+
  geom_boxplot(varwidth = FALSE,outlier.shape = NA)+
  coord_flip()+
  facet_wrap(.~WB_lb)+
  # facet_wrap(.~WB_lb)+
  geom_point(aes(group=LF03), position=position_jitter(width = 0.2),alpha=0.2)+
  scale_colour_manual(values = cbPalette2)+
  scale_x_discrete(limits=rev)+
  labs(
    title = "Total carbon content by zooplankton type by EA water body",
    y="log(Total carbon)",
    caption = paste0("Samples gathered between ",format(min(dfl$sample.date), "%d/%m/%Y")," & ",format(max(dfl$sample.date), "%d/%m/%Y"),
                     "\nCarbon content values are based on *mean* estimates of carbon contents per taxon"))+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        strip.text = element_text(face=2),
        axis.title.x = element_text(face=2),
        axis.text.x = element_text(face=2),
        axis.text.y = element_text(face=2,
                                   size = 7)) -> pl

ggsave(plot = pl, filename = "figs/2412dd/carbonByZooTypeWB_Mean.pdf",
       width = 20,height = 12,units = "in"); rm(pl)

# Seahorse food: total C ####
## Median ####
dfl %>% #names()
  dplyr::select(.,BIOSYS.Code,WIMS.Code,Pot.Number,WB_lb,sample.date,
                DisplayName:md_carbTot_m3) %>% #names()#View()
  ## remove unneeded
  dplyr::select(.,-c(
    DisplayName,
    AbundanceRaw, Abund_m3,
    copNonCop,mnlongMaxAxis_mm, mdlongMaxAxis_mm,
    mnCPerIndiv_ug, mdCPerIndiv_ug, mn_carbTot_raw,
    md_carbTot_raw,mn_carbTot_m3)) %>% #View()
  group_by(across(c(!md_carbTot_m3))) %>% 
  summarise(sum=sum(md_carbTot_m3, na.rm = TRUE), .groups = "drop") %>%
  filter(.,WB_lb %in% c("Sth_Solent","Sth_SotonWtr")) %>% 
  filter(.,WIMS.Code != "Y00017477") %>% 
  mutate(.,WIMS.Code = factor(WIMS.Code, levels = c("G0003572",
                                                    "Y0017477",
                                                    "Y0004367",
                                                    "G0003532"))) %>% #View()
  ggplot(.,aes(x=BIOSYS.Code,y=sum, colour=WB_lb))+
  geom_hline(yintercept = seq(from=0, to=150000, by=10000), colour = "grey",
             linetype=2)+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(width = 0.25)+
  labs(title = "Total carbon content of zooplankters recorded in the Solent and Southampton Water",
       x="BIOSYS site code",
       y= "Total carbon (ug C/m3)",
       caption = paste0("Samples gathered between ",format(min(dfl$sample.date), "%d/%m/%Y")," & ",format(max(dfl$sample.date), "%d/%m/%Y"),
                        "\nCarbon content values are based on *median* estimates of carbon contents per taxon"))+
  scale_y_continuous(breaks = c(0,50000,100000,150000),labels = scales::comma_format())+
  coord_flip()+
  theme(legend.title = element_blank(),
        axis.title = element_text(face=2)) -> pl

ggsave(plot = pl, filename = "figs/2412dd/carbonTot_Soton_Median.pdf",
       width = 20,height = 12,units = "in");rm(pl)

## Mean ####
dfl %>% #names()
  filter(.,WB_lb %in% c("Sth_Solent","Sth_SotonWtr")) %>% 
  filter(.,WIMS.Code != "Y00017477") %>% 
  mutate(.,WIMS.Code = factor(WIMS.Code, levels = c("G0003572",
                                                    "Y0017477",
                                                    "Y0004367",
                                                    "G0003532"))) %>% 
  dplyr::select(.,PRN,BIOSYS.Code,WIMS.Code,Pot.Number,WB_lb,sample.date,
                DisplayName:md_carbTot_m3) %>%# View()#names()
  ## remove unneeded
  dplyr::select(.,-c(
    DisplayName,
    AbundanceRaw, Abund_m3,
    copNonCop,mnlongMaxAxis_mm, mdlongMaxAxis_mm,
    mnCPerIndiv_ug, mdCPerIndiv_ug, mn_carbTot_raw,
    md_carbTot_raw,md_carbTot_m3)) %>% #names()
  group_by(across(c(!mn_carbTot_m3))) %>% 
  summarise(sum=sum(mn_carbTot_m3, na.rm = TRUE), .groups = "drop") %>%
  ggplot(.,aes(x=BIOSYS.Code,y=sum, colour=WB_lb))+
  geom_hline(yintercept = seq(from=0, to=330000, by=10000), colour = "grey",
             linetype=2)+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(width = 0.25)+
  labs(title = "Total carbon content of zooplankters recorded in the Solent and Southampton Water",
       x="BIOSYS site code",
       y= "Total carbon (ug C/m3)",
       caption = paste0("Samples gathered between ",format(min(dfl$sample.date), "%d/%m/%Y")," & ",format(max(dfl$sample.date), "%d/%m/%Y"),
                        "\nCarbon content values are based on *mean* estimates of carbon contents per taxon"))+
  scale_y_continuous(breaks = seq(from=0,to=330000,by=50000),labels = scales::comma_format())+
  coord_flip()+
  theme(legend.title = element_blank(),
        axis.title = element_text(face=2)) -> pl

ggsave(plot = pl, filename = "figs/2412dd/carbonTot_Soton_Mean.pdf",
       width = 20,height = 12,units = "in");rm(pl)

# Seahorse food: By taxon ####
## Median ####
dfl %>%
  dplyr::select(.,c("BIOSYS.Code","WIMS.Code","sample.date","PRN",
                    "WB_lb",LF0:DisplayName,"Abund_m3","md_carbTot_m3")) %>% 
  filter(.,WB_lb %in% c("Sth_Solent","Sth_SotonWtr")) %>%
  mutate(.,label = ifelse(LF02 == "Cop_Sm","Cop_Sm",
                          ifelse(LF02 == "Cop_Lg","Cop_Lg",
                                 ifelse(LF02 == "Cop_Ambi","Cop_Ambi",
                                        ifelse(LF02 == "Cop_NYA","Cop_NYA",
                                               ifelse(Order == "Mysida","Mysida","X")))))) %>%#View()
  mutate(.,label = ifelse(DisplayName=="Cirripedia","Cirripedia",label)) %>% 
  mutate(.,label = ifelse(is.na(label),"X",label)) %>% 
  filter(.,label != "X") %>% 
  group_by(BIOSYS.Code, WIMS.Code, sample.date, PRN, WB_lb, label) %>%
  summarise(
    md_carbTot_m3 = sum(md_carbTot_m3, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  ggplot(.,aes(x= BIOSYS.Code,y=log(md_carbTot_m3+1), colour=WB_lb))+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(.~label, scales = "free_y")+
  geom_jitter(width = 0.25)+
  labs(title = "Total carbon content of selected zooplankters recorded in the Solent and Southampton Water",
       x="BIOSYS site code",
       y= "Log (n+1) carbon content (ug C/m3)",
       caption = paste0("Samples gathered between ",format(min(dfl$sample.date), "%d/%m/%Y")," & ",format(max(dfl$sample.date), "%d/%m/%Y"),
                        "\nCarbon content values are based on *median* estimates of carbon contents per taxon"))+
  theme(legend.title = element_blank(),
        axis.title = element_text(face=2))+
  coord_flip() -> pl

ggsave(plot = pl, filename = "figs/2412dd/carbonTax_Soton_Median.pdf",
       width = 20,height = 12,units = "in");rm(pl)

## Mean ####
dfl %>%
  dplyr::select(.,c("BIOSYS.Code","WIMS.Code","sample.date","PRN",
                    "WB_lb",LF0:DisplayName,"Abund_m3","mn_carbTot_m3")) %>% 
  filter(.,WB_lb %in% c("Sth_Solent","Sth_SotonWtr")) %>%
  mutate(.,label = ifelse(LF02 == "Cop_Sm","Cop_Sm",
                          ifelse(LF02 == "Cop_Lg","Cop_Lg",
                                 ifelse(LF02 == "Cop_Ambi","Cop_Ambi",
                                        ifelse(LF02 == "Cop_NYA","Cop_NYA",
                                               ifelse(Order == "Mysida","Mysida","X")))))) %>%#View()
  mutate(.,label = ifelse(DisplayName=="Cirripedia","Cirripedia",label)) %>% 
  mutate(.,label = ifelse(is.na(label),"X",label)) %>% 
  filter(.,label != "X") %>% 
  group_by(BIOSYS.Code, WIMS.Code, sample.date, PRN, WB_lb, label) %>%
  summarise(
    mn_carbTot_m3 = sum(mn_carbTot_m3, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  ggplot(.,aes(x= BIOSYS.Code,y=log(mn_carbTot_m3+1), colour=WB_lb))+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(.~label, scales = "free_y")+
  geom_jitter(width = 0.25)+
  labs(title = "Total carbon content of selected zooplankters recorded in the Solent and Southampton Water",
       x="BIOSYS site code",
       y= "Log (n+1) carbon content (ug C/m3)",
       caption = paste0("Samples gathered between ",format(min(dfl$sample.date), "%d/%m/%Y")," & ",format(max(dfl$sample.date), "%d/%m/%Y"),
                        "\nCarbon content values are based on *mean* estimates of carbon contents per taxon"))+
  theme(legend.title = element_blank(),
        axis.title = element_text(face=2))+
  coord_flip() -> pl

ggsave(plot = pl, filename = "figs/2412dd/carbonTax_Soton_Mean.pdf",
       width = 20,height = 12,units = "in");rm(pl)

############
## export carbon contents per sample
dfl %>% 
  dplyr::select(.,c(
    Pot.Number,
    sample.date:WIMS.Code,
    `Net.volume.sampled.(m3)`,
    PRN,
    Eastings:Category,
    AbundanceRaw:md_carbTot_m3)) %>% 
  dplyr::select(., -copNonCop) %>% 
  group_by(across(-c(
    AbundanceRaw,
    Abund_m3,
    mnlongMaxAxis_mm,
    mdlongMaxAxis_mm,
    mnCPerIndiv_ug,
    mdCPerIndiv_ug,
    mn_carbTot_raw,
    md_carbTot_raw,
    mn_carbTot_m3,
    md_carbTot_m3
  ))) %>% 
  summarise(.,
            AbundanceRaw=sum(AbundanceRaw,na.rm = TRUE),
            Abund_m3=sum(Abund_m3,na.rm = TRUE),
            mnlongMaxAxis_mm=mean(mnlongMaxAxis_mm,na.rm = TRUE),
            mdlongMaxAxis_mm=mean(mdlongMaxAxis_mm,na.rm = TRUE),
            mnCPerIndiv_ug=mean(mnCPerIndiv_ug,na.rm = TRUE),
            mdCPerIndiv_ug=mean(mdCPerIndiv_ug,na.rm = TRUE),
            mn_carbTot_raw=sum(mn_carbTot_raw,na.rm = TRUE),
            md_carbTot_raw=sum(md_carbTot_raw,na.rm = TRUE),
            mn_carbTot_m3=sum(mn_carbTot_m3,na.rm = TRUE),
            md_carbTot_m3=sum(md_carbTot_m3,na.rm = TRUE),
  .groups = "drop") ->df_tot

## write outs
write.csv(dfl,file="outputs/abundances_by_taxon.csv",row.names = FALSE)
write.csv(df_tot,file="outputs/abundances_by_sample.csv",row.names = FALSE)
