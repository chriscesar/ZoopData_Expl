## RichnessAbundancePlots.R ####
# 00 set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

## set universals ####
tictoc::tic.clearlog();tic("Set universals");print("set universals")
source("R/set_meta.R")
toc(log=TRUE)

## load data ####
source("R/imp.load_data_all_zoops.R")
rm(dfl0)

## widen data and calculate taxon richnesses and abundances

names(dfl)
dfl %>% 
  dplyr::select(.,c(Pot.Number:WB_lb,DisplayName,Abund_m3)) %>% #names()
  dplyr::select(.,-c(Aphia.ID,Taxa)) %>% 
  dplyr::select(.,-c(
  max_axis_length_mm, ugC_per_individ_mn,ugC_in_sample_mn,
  ugC_per.m3_mn,mgC_per_individ_mn,mgC_in_sample_mn,
  mgC_per.m3_mn,ugC_per_individ_md,ugC_in_sample_md,
  ugC_per.m3_md, mgC_per_individ_md,mgC_in_sample_md,
  mgC_per.m3_md
  )
  ) %>% 
  group_by(., across(!c(Abund_m3))) %>% 
  summarise(.,Abund_m3 = sum(Abund_m3),.groups = "drop") %>% #names()
  ungroup(.) %>% 
  pivot_wider(.,
              names_from = DisplayName,
              values_from = Abund_m3,
              values_fill = 0
              ) -> dfw

dftmp <- dfw %>% dplyr::select(.,-c(Pot.Number:WB_lb))

S <- vegan::specnumber(dftmp)
N <- rowSums(dftmp)

dfw$S <- S
dfw$N <- N

mnS <- mean(S);mnN <- mean(N)

## plot taxon richness
dfw %>% 
  ggplot(., aes(x= WB_lb, y = S, colour = Region))+
  # geom_hline(yintercept = mnS, lty=2)+
  geom_hline(yintercept = mean(S), lty=2)+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(width=0.25, alpha=0.4)+
  labs(
    title = "Taxon richness in zooplankton samples",
    caption=paste0("\nSamples gathered between ",
                   format(min(dfl$sample.date),
                          "%d/%m/%Y")," & ",format(max(dfl$sample.date), "%d/%m/%Y")),
    y="Taxon richness")+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = 2),
        axis.text = element_text(face = 2)) -> pl
pdf(file = paste0("figs/TaxRich_",
                  format(min(dfl$sample.date), "%y%m%d"),"_",
                  format(max(dfl$sample.date), "%y%m%d"),".pdf"),
    width=18,height=10) #scaled
print(pl)
dev.off();rm(pl)

## plot taxon abundance
dfw %>% 
  ggplot(., aes(x= WB_lb, y = log(N), colour = Region))+
  # geom_hline(yintercept = log(mean(N)), col= "black",lty=2)+
  geom_hline(yintercept = mean(log(N)), col= "black",lty=2)+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(width=0.25, alpha=0.4)+
  labs(
    title = "Log taxon abundances in zooplankton samples",
    caption=paste0("\nSamples gathered between ",
                   format(min(dfl$sample.date),
                          "%d/%m/%Y")," & ",format(max(dfl$sample.date), "%d/%m/%Y")),
    y="Log taxon abundance")+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = 2),
        axis.text = element_text(face = 2)) -> pl
pdf(file = paste0("figs/TaxAbund_",
                  format(min(dfl$sample.date), "%y%m%d"),"_",
                  format(max(dfl$sample.date), "%y%m%d"), ".pdf"),
    width=18,height=10) #scaled
print(pl)
dev.off();rm(pl)

# time series plots ####
# Taxon Richness (S)
dfw %>% mutate(year = lubridate::year(sample.date)) %>% 
  ggplot(.,aes(x=sample.date, y=S))+
  geom_point(aes(colour=Region),show.legend = FALSE)+
  facet_wrap(~WB_lb)+
  geom_smooth()+
  labs(
    title = "Taxon richness in zooplankton samples",
    y="Taxon richness",
    caption=paste0("Samples gathered between ",format(min(dfw$sample.date), "%d/%m/%Y")," & ",format(max(dfw$sample.date), "%d/%m/%Y")))+
theme(axis.title.x = element_blank()) -> pl
pdf(file = paste0("figs/TaxRich_ts_",
                  format(min(dfl$sample.date), "%y%m%d"),"_",
                  format(max(dfl$sample.date), "%y%m%d"),
                  ".pdf"),
    width=18,height=10) #scaled
print(pl)
dev.off();rm(pl)

dfw %>% mutate(year = lubridate::year(sample.date)) %>% 
  ggplot(.,aes(x=year, y=S,group = year))+
  geom_boxplot(outliers = FALSE,varwidth = TRUE)+
  geom_jitter(width=0.25,aes(colour=Region),show.legend = FALSE)+
  facet_wrap(~WB_lb)+
  #geom_smooth()+
  labs(
    title = "Taxon richness in zooplankton samples",
    y="Taxon richness",
    caption=paste0("Samples gathered between ",format(min(dfw$sample.date), "%d/%m/%Y")," & ",format(max(dfw$sample.date), "%d/%m/%Y")))+
  scale_x_continuous(breaks = scales::breaks_pretty(n = 3))+
  theme(axis.title.x = element_blank()) -> pl
pdf(file = paste0("figs/TaxRich_ts_box_",
                  format(min(dfl$sample.date), "%y%m%d"),"_",
                  format(max(dfl$sample.date), "%y%m%d"),
                  ".pdf"),width=18,height=10) #scaled
print(pl)
dev.off();rm(pl)

# Taxon Abundance (N)
dfw %>% mutate(year = lubridate::year(sample.date)) %>% 
  ggplot(.,aes(x=sample.date, y=log(N)))+
  geom_point(aes(colour=Region),show.legend = FALSE)+
  facet_wrap(~WB_lb)+
  geom_smooth()+
  labs(
    title = "Taxon abundance in zooplankton samples",
    y="Log taxon abundance",
    caption=paste0("Samples gathered between ",format(min(dfw$sample.date), "%d/%m/%Y")," & ",format(max(dfw$sample.date), "%d/%m/%Y")))+
  theme(axis.title.x = element_blank()) -> pl
pdf(file = paste0("figs/TaxAbund_ts_",
                  format(min(dfl$sample.date), "%y%m%d"),"_",
                  format(max(dfl$sample.date), "%y%m%d"),
                  ".pdf"),width=18,height=10) #scaled
print(pl)
dev.off();rm(pl)

dfw %>% mutate(year = lubridate::year(sample.date)) %>% 
  ggplot(.,aes(x=year, y=log(N),group = year))+
  geom_boxplot(outliers = FALSE,varwidth = TRUE)+
  geom_jitter(width=0.25,aes(colour=Region),show.legend = FALSE)+
  facet_wrap(~WB_lb)+
  labs(
    title = "Taxon abundance in zooplankton samples",
    y="log taxon abundance",
    caption=paste0("Samples gathered between ",format(min(dfw$sample.date), "%d/%m/%Y")," & ",format(max(dfw$sample.date), "%d/%m/%Y")))+
  scale_x_continuous(breaks = scales::breaks_pretty(n = 3))+
  theme(axis.title.x = element_blank()) -> pl
pdf(file = paste0("figs/TaxAbund_ts_box_",
    format(min(dfl$sample.date), "%y%m%d"),"_",
    format(max(dfl$sample.date), "%y%m%d"),".pdf"),
    width=18,height=10) #scaled
print(pl)
dev.off();rm(pl)


## tidy
rm(list=ls(pattern = "^df"))
rm(list=ls(pattern = "^cb"))
rm(list=ls(pattern = "^mn"))
rm(datfol,N,nit,perms,ppi,S)
