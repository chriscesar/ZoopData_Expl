# 202602_zoop_figs_for_rep.R ####
## generate figures for 2026 Inshore Plankton Report

# 00 set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "tictoc","dplyr")
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

tic("Format data")
## keep only 'interesting' variables
df <- dfl %>% 
  dplyr::select(Pot.Number,Region,WB,sample.date,DisplayName,LF02,
         Abund_m3,md_carbTot_m3_ug) %>% 
  # create yyymm variable & move it
  mutate(yyyymm = sprintf("%d%02d",
                          year(sample.date),
                          month(sample.date))) %>% 
  relocate(yyyymm, .after = sample.date) %>% 
  mutate(date_plot = make_date(year = year(sample.date),
                               month = month(sample.date),
                               day = 1)) %>% 
  relocate(date_plot, .after=yyyymm)

toc(log=TRUE)

# summarise counts by sampling event & plot ####
# png(file = paste0(outfol,"zoop_abund_log10.png"),
png(file = paste0(outfol,
                  # "zoop_abund_raw.png"
                  "zoop_abund_log10.png"
                  ),
    width=18*ppi, height=12*ppi, res=ppi)
df %>%
  ## total abundance: individuals per m3
  dplyr::select(Pot.Number,date_plot,Abund_m3) %>% #View()
  group_by(across(-Abund_m3)) %>% 
  # get sum per sample
  summarise(Abund_m3 = sum(Abund_m3,na.rm = TRUE),
            .groups = "drop") %>% 
  # get mean sum by date
  ## drop sample ID
  dplyr::select(-Pot.Number) %>% group_by(date_plot) %>% 
  summarise(mn_Abund_m3 = mean(Abund_m3,na.rm = TRUE),
            sd_Abund_m3 = sd(Abund_m3,na.rm=TRUE),
            n = n(),
            .groups = "drop"
            ) %>% 
  # filter(date_plot <"2025-06-01") %>% 
  ggplot(.,
         aes(
           x = date_plot,
           y = log10(mn_Abund_m3+1)
           # y=mn_Abund_m3
         )
  )+
  geom_rug(sides = "b")+
  geom_vline(xintercept = as.Date(c("2023-01-01","2024-01-01","2025-01-01")),
  lty = 2, col="grey")+
  geom_vline(xintercept = seq(
    from = as.Date("2022-06-01"),
    to   = as.Date("2025-06-01"),
    by   = "1 month"),
    lty = 3, col="grey")+
  geom_smooth(method = "gam")+
  geom_point(aes(
    fill = n
    ),
    size = 6,pch=21)+
  labs(
    title = "Log10 mean zooplankton abundances by year_month",
    # title = "Mean zooplankton abundances by year_month",
    y = "Log10(n+1) mean abundances (individuals per m3) per sample",
    # y = "Mean abundances (individuals per m3) per sample",
    caption=paste0("Values represent mean total counts across all samples gathered in a particular month","<br>",
                   "Point shading reflects number of samples contributing to that mean",
                   "<br>","Blue line represents generalised additive model trend"
                   ),
       fill = "Num.<br>samples"
  )+
  # scale_colour_binned() +
  ylim(0,NA)+
  scale_fill_stepsn(
    # colours = c("red", "yellow", "green", "yellow", "red"),
    # breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
    colours = c("#DEEBF7", "#9ECAE1", "#3182BD"),
    breaks = c(20, 30,40)
  )+
  theme(
    # palette.color.continuous = c("#DEEBF7", "#9ECAE1", "#3182BD"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2,size=12),
    axis.text = element_text(face=2),
    axis.text.x = element_text(size = 12),
    plot.caption = ggtext::element_markdown(face=2,size=12),
    plot.title = element_text(face=2,size = 14),
    legend.title = ggtext::element_markdown(face=2),
  )
rm(pad,ylims,df_sum_tmp)
dev.off()

# summarise carbon by sampling event & plot ####
# png(file = paste0(outfol,"zoop_carbon_log10.png"),
png(file = paste0(outfol,"zoop_carbon_raw.png"),
    width=18*ppi, height=12*ppi, res=ppi)
df %>%
  ## total abundance: cells per litre
  dplyr::select(Pot.Number,date_plot,md_carbTot_m3_ug) %>% #View()
  group_by(across(-md_carbTot_m3_ug)) %>% 
  # get sum per sample
  summarise(md_carbTot_m3_ug = sum(md_carbTot_m3_ug,na.rm = TRUE),
            .groups = "drop") %>% 
  # get mean sum by date
  ## drop sample ID
  dplyr::select(-Pot.Number) %>% group_by(date_plot) %>% 
  summarise(mn_Abund_m3 = mean(md_carbTot_m3_ug,na.rm = TRUE),
            sd_Abund_m3 = sd(md_carbTot_m3_ug,na.rm=TRUE),
            n = n(),
            .groups = "drop"
  ) %>% 
  ggplot(.,
         aes(
           x = date_plot,
           # y = log10(mn_Abund_m3+1)
           y=mn_Abund_m3
         )
  )+
  geom_rug(sides = "b")+
  geom_vline(xintercept = as.Date(c("2023-01-01","2024-01-01","2025-01-01")),
             lty = 2, col="grey")+
  geom_vline(xintercept = seq(
    from = as.Date("2022-06-01"),
    to   = as.Date("2025-06-01"),
    by   = "1 month"),
    lty = 3, col="grey")+
  geom_smooth(method = "gam")+
  geom_point(aes(
    fill = n
    ),
    size=6,pch=21)+
  labs(
    # title = "Log10 mean zooplankton carbon content by year_month",
    title = "Mean zooplankton carbon content by year_month",
    # y = "Log10(n+1) mean total estimated carbon content (ug per m3) per sample",
    y = "Mean total estimated carbon content (ug per m3) per sample",
    caption=paste0("Values represent mean total carbon content within all samples gathered in a particular month","<br>",
                   "Point shading reflects number of samples contributing to that mean",
                   "<br>","Blue line represents generalised additive model trend"
    ),
    fill = "Num.<br>samples"
  )+
  # scale_colour_binned() +
  scale_fill_stepsn(
    # colours = c("red", "yellow", "green", "yellow", "red"),
    # breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
    colours = c("#FEE0D2", "#FC9272", "#DE2D26"),
    breaks = c(20, 30,40)
  )+
  theme(
    # palette.color.continuous = c("#DEEBF7", "#9ECAE1", "#3182BD"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2,size=12),
    axis.text = element_text(face=2),
    axis.text.x = element_text(size = 12),
    plot.caption = ggtext::element_markdown(face=2,size=12),
    plot.title = element_text(face=2,size = 14),
    legend.title = ggtext::element_markdown(face=2),
  )
rm(pad,ylims,df_sum_tmp)
dev.off()

