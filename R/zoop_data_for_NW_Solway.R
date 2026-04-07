# zoop_data_for_NW_Solway.R ####
## generate figures for 2026 Inshore Plankton Report

# 00 set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "tictoc","dplyr","ggthemes")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

## set universals ####
tictoc::tic.clearlog();tic("Set universals");print("set universals")
source("R/set_meta.R")
outfol <- "//prodds.ntnl/Shared/AN/KFH/Groups/N_Marine/02 Projects_Tasks/04 E&B/Ecology and Ecosystems/Planktonic_Food_Webs/202602_FigsForReport/"
toc(log=TRUE)

# load and format data ####
tic("Load zoop data")
source("R/imp.load_data_all_zoops.R");rm(dfl0)
toc(log=TRUE)

# retain only NW & SOLWAY ####
wbs <- c("Cumbria","Solway Outer South","Mersey Mouth")
dfnw <- dfl %>% 
  filter(WB %in% wbs);rm(wbs)

# convert Excel dates to date
dfnw$sample.date <- (as.Date(dfnw$sample.date, origin = "1899-12-30"))

# table(dfl0$Region)
write.csv(dfnw, file = "outputs/zoops_nw.csv",row.names = FALSE)

# plot ####
# Create a sequence of January 1st for each year spanned by your data
year_lines <- tibble(
  year_start = seq(
    floor_date(min(dfnw$sample.date, na.rm = TRUE), unit = "year"),
    floor_date(max(dfnw$sample.date, na.rm = TRUE), unit = "year"),
    by = "1 year"
  )
)

dfnw %>% 
  dplyr::select(
    sample.date,
    WB,Abund_m3,
    Pot.Number,
    #LF02
    ) %>% 
  group_by(across(-Abund_m3)) %>% 
  summarise(Abund_m3 = sum(Abund_m3, na.rm = TRUE),
            .groups = "drop") %>% ungroup() %>% 
  ggplot(.,aes(
    x = sample.date,
    y = log10(Abund_m3+1)
  )) +
  geom_vline(data = year_lines,
             aes(xintercept = as.numeric(year_start)),
             color = "grey70", linewidth = 0.3,lty=3)+
  geom_point()+
  facet_wrap(.~WB)+
  ggthemes::theme_few()+
  geom_smooth(method="gam")+
  labs(
    title = "Zooplankton abundances recorded in selected English water bodies",
    subtitle = "Data faceted by Environment Agency water body",
    # caption = "EA data filtered to retain taxa belonging to Infraphylum Dinoflagellata",
    y="Log10 individuals per litre (n+1)"
  )+
  theme(
    strip.text = element_text(face=2,size = 14),
    axis.text = element_text(face=2,size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2,size = 14),
    plot.title = element_text(face=2,size = 16),
    plot.subtitle = element_text(face=2,size = 14),
    plot.caption = element_text(face=2,size = 14),
  )

ggsave(plot = get_last_plot(),
       filename = paste0("outputs/SolwayNW_zoops.png"),
       width = 18, height = 8, units = "in"
       )
