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
  geom_hline(yintercept = mnS, lty=2)+
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
pdf(file = "figs/2412dd/TaxRich.pdf",width=18,height=10) #scaled
print(pl)
dev.off();rm(pl)

## plot taxon abundance
dfw %>% 
  ggplot(., aes(x= WB_lb, y = log(N), colour = Region))+
  geom_hline(yintercept = log(mnN), lty=2)+
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
pdf(file = "figs/2412dd/TaxAbund.pdf",width=18,height=10) #scaled
print(pl)
dev.off();rm(pl)

