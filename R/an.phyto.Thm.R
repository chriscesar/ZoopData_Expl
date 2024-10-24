# an.phyto.Thm.R ####
# analysis and quick plot generation of Phyto data for the Thames

# load packages ####
ld_pkgs <- c("tidyverse", "tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

# set universals ####
tictoc::tic.clearlog();tic("Set universals");print("set universals")
source("R/set_meta.R")
toc(log=TRUE)

# Load data ####
tic("LOAD DATA");print("LOAD DATA")
df0_biom <- readxl::read_xlsx("data/Copy of PhytoChl_2000-2020 (WB+RegionalSeas).xlsx",
                              sheet = "Biomass2000_2020",guess_max = 10000)
df0_abnd <- readxl::read_xlsx("data/Copy of PhytoChl_2000-2020 (WB+RegionalSeas).xlsx",
                              sheet = "Abundance2000_2020",guess_max = 10000)
toc(log=TRUE)

### filter out SNS data
df0_abnd_sns <- df0_abnd %>% 
  filter(
    str_starts(str_to_lower(Waterbody),"thames")|
      str_starts(str_to_lower(Waterbody),"essex")|
      str_starts(str_to_lower(Waterbody),"whitstable")|
      str_starts(str_to_lower(Waterbody),"medway")|
      str_starts(str_to_lower(Waterbody),"swale")|
      str_starts(str_to_lower(Waterbody),"kent n"))

df0_biom_sns <- df0_biom %>% 
  filter(
    str_starts(str_to_lower(Waterbody),"thames")|
      str_starts(str_to_lower(Waterbody),"essex")|
      str_starts(str_to_lower(Waterbody),"whitstable")|
      str_starts(str_to_lower(Waterbody),"medway")|
      str_starts(str_to_lower(Waterbody),"swale")|
      str_starts(str_to_lower(Waterbody),"kent n"))

df0_abnd_sns$Waterbody <- factor(df0_abnd_sns$Waterbody, levels = c(
  "THAMES UPPER", "THAMES MIDDLE","THAMES LOWER",
  "MEDWAY","SWALE",
  "Thames Coastal North","Thames Coastal South", "Whitstable Bay",
  "Essex","Kent North"
  ))

df0_biom_sns$Waterbody <- factor(df0_biom_sns$Waterbody, levels = c(
  "THAMES UPPER", "THAMES MIDDLE","THAMES LOWER",
  "MEDWAY","SWALE",
  "Thames Coastal North","Thames Coastal South", "Whitstable Bay",
  "Essex","Kent North"
))

# plots ####
##BIOMASS ####
df0_biom_sns %>% 
  filter(., Year > 2002) %>% 
  mutate(Biomass_USE = if_else(
    str_detect(Biomass,"^<"), # Detect "<" at the beginning
    as.numeric(str_remove(Biomass, "^<"))/2, # Remove "<" and divide by 2
    as.numeric(Biomass) # Keep numeric values as they are
  )) %>% 
  ggplot(., aes(x = Year, y = Biomass_USE))+
  geom_boxplot(aes(group = Year),outliers = FALSE)+
  geom_jitter(alpha=0.2,height = 0, width = 0.3)+
  facet_wrap(.~Waterbody,ncol = 3)+
  labs(title = "Biomass values over time in water bodies within the Greater Thames embayment",
       y="Biomass (ug/l)")+
  theme(
    axis.title.x = element_blank(),
    strip.text = element_text(face="bold")
    )+
  #geom_smooth(method = "gam")#+
  #geom_smooth(method = "loess", col=3)+
  geom_smooth(method = "lm", col=2)

## ABUNDANCE ####
df0_abnd_sns %>% 
  dplyr::select(., -c(Index,`Size Class`, Taxon, aphiaID, TorN)) %>% 
  group_by(across(c(!Abundance))) %>% 
  summarise(Abundance=sum(Abundance, na.rm = TRUE), .groups = "drop") %>% 
  #ggplot(., aes(x = Year, y = Abundance))+
  ggplot(., aes(x = Year, y = log(Abundance+1)))+
  geom_boxplot(aes(group = Year),outliers = FALSE)+
  geom_jitter(alpha=0.2,height = 0, width = 0.3)+
    #ylim(0,60000)+
  facet_wrap(.~Waterbody,ncol = 3)+
  labs(title = "Phytoplankton abundances over time in water bodies within the Greater Thames embayment",
        y="log abundance (n+1)"
       #y="Abundance"
       )+
    theme(
      axis.title.x = element_blank(),
      strip.text = element_text(face="bold")
    )#+
  geom_smooth(method = "gam")
