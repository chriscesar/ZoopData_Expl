# 2026_08_zoop_trends_expl.01.R ----

# 00 set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "tictoc","dplyr","ggthemes","vegan")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

## set universals ####
tictoc::tic.clearlog();tic("Set universals");print("set universals")
source("R/set_meta.R")
toc(log=TRUE)

# load and format data ####
tic("Load zoop data")
source("R/imp.load_data_all_zoops.R");rm(dfl0)
toc(log=TRUE)

# Retain only fish eggs ----
df_egg <- dfl %>% 
  dplyr::filter(
    STAGE == "EG" &
      phylum == "Chordata"
    )

df_egg %>% dplyr::rename(Taxon = DisplayName) %>% 
  dplyr::filter(Taxon == "Engraulis #eggs") %>% 
  ggplot(.,
         aes(
           x = sample.date,
           y = Abund_m3
         )) +
  geom_point(aes(
    fill = Taxon,
    shape = Taxon
    ),
    size = 3
    )+
  scale_shape_manual(values = c(21,22,23,24)) +
  # guides(colour = guide_legend(title = "Taxon")) +
  labs(
    x="",
    y = "Abundance per m3",
    ) +
  facet_wrap(.~Region)+
  # geom_smooth()+
  theme(
    axis.title = element_text(face = 2),
    axis.text = element_text(face = 2),
    strip.text = element_text(face = 2),
  )

# identify samples where undetermined and anchovy eggs occur
df_egg %>% 
  dplyr::select(Pot.Number, sample.date,WB,Region,DisplayName,Abund_m3) %>% 
  group_by(Pot.Number) %>%
  filter(
    all(c("Fish unident #eggs", "Engraulis #eggs") %in% DisplayName)
  ) %>%
  ungroup()
