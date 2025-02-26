# phytoCarbonComparison.R ####
# compare PML and EA estimates of carbon content

# 00 set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "tictoc","ggtext","glue")

vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

## set universals ####
tictoc::tic.clearlog();tic("Set universals");print("set universals")
source("R/set_meta.R")
toc(log=TRUE)

tic("Load phyto data")
df0 <- readxl::read_xlsx(paste0(datfol,"Phytodata/Phyto carbon PML and EA_v7 Feb 2025_CC_edits.xlsx"),
                        sheet = "WORKING")
toc(log=TRUE)

df0 %>% #View()
  dplyr::select(.,-c(`ABS_Diff_EA-PML_Carb`,`APHIA ID`,EA_n,PML_n,
                     EA_Mean_Vol_per_cell_um3, PML_Mean_Vol_per_cell_um3,
                     `ORDERS_OF_MAGNITUDE DIFFERENCE`,
                     `ABS_Diff_EA-PML_Vol`)) %>% #names()
  mutate(., Diff = (EA_Mean_C_per_cell_pgC-PML_Mean_C_per_cell_pgC)) %>% #View()
  mutate(., AbsLogDiff = log10(abs(Diff))) %>% #View()
  mutate(.,LogDiffSign = ifelse(Diff<0,-1*AbsLogDiff,AbsLogDiff)) %>% 
  ## orders of magnitude diffs
  mutate(.,OrdMagDiff = log10(EA_Mean_C_per_cell_pgC)-log10(PML_Mean_C_per_cell_pgC)) %>% 
  mutate(.,OrdMagDiffSign = ifelse(Diff<0,-1*OrdMagDiff,OrdMagDiff)) -> df

# plot ####
# png(file = "figs/PhytoCarbonComparison_v1.png",
#     width=9*ppi, height=12*ppi, res=ppi)
png(file = "figs/PhytoCarbonComparison_v2_wide.png",
    width=12*ppi, height=8*ppi, res=ppi)
df %>% 
  arrange(LogDiffSign) %>% #View() %>% 
  mutate(Name = factor(Name, levels = Name)) %>%  # Reorder factor levels
  ggplot(., aes(y = Name, x = LogDiffSign))+
  geom_rect(aes(xmin = -1,xmax = 1, ymin = -Inf,ymax=Inf),fill="lightgrey")+
  geom_rect(aes(xmin = -Inf,xmax = -1, ymin = -Inf,ymax=Inf), fill = "lightblue")+
  geom_rect(aes(xmin = 1,xmax = Inf, ymin = -Inf,ymax=Inf), fill = "lightgreen")+
  geom_vline(xintercept = seq(-5,5, by = 1),col="grey", lty=3)+
  geom_vline(xintercept = 0, lty=2, col=2)+
  geom_point()+
  scale_x_continuous(breaks = seq(-5, 5, by = 1))+
  labs(title = "Differences in phytoplankton carbon content estimates for Environment Agency\nand Plymouth Marine Laboratory",
    x=expression(bold(Log[10]~"difference in carbon content values")),
       caption = paste0("Values represent the log10 of the absolute differences in estimates between institutions.  Where EA
       Values <0 indicate estimates for which PML values are larger; values >0 indicate estimates for which EA values are larger.
                        Values are based on estimates of carbon content (pg per cell) by each institution.
                        Differences in estimates are provided for the ",nrow(df),
                        " taxa which have carbon content estimates from both institutions.
                        The grey area indicates estimates which are within 1 order of magnitude of each other"))+
  theme(
    axis.text.x = element_text(face=2),
    axis.title.x = element_text(face=2),
    axis.title.y= element_blank(),
    axis.text.y = element_blank(),
    plot.caption = element_text(size=11)
  )
dev.off()

# png(file = "figs/PhytoCarbonComparison_v2_wide.png",
    # width=9*ppi, height=12*ppi, res=ppi)
png(file = "figs/PhytoCarbonComparison_v2_wide.png",
    width=12*ppi, height=8*ppi, res=ppi)
df %>% 
  arrange(OrdMagDiff) %>% 
  mutate(Name = factor(Name, levels = Name)) %>%  # Reorder factor levels
  ggplot(aes(y = Name, x = OrdMagDiff)) +
  geom_rect(aes(xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf), fill = "lightgrey") +
  geom_rect(aes(xmin = -Inf, xmax = -1, ymin = -Inf, ymax = Inf), fill = "lightblue") +
  geom_rect(aes(xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "lightgreen") +
  geom_vline(xintercept = seq(-3, 3, by = 1), col = "grey", lty = 3) +
  geom_vline(xintercept = 0, lty = 2, col = 2) +
  geom_point() +
  scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
  labs(
    title = "Differences in phytoplankton carbon content estimates for Environment Agency<br>and Plymouth Marine Laboratory",
    x = expression(bold(Log[10]~"difference in carbon content values")),
    caption = glue(
      "Orders of magnitude difference calculated using log<sub>10</sub>(EA) - log<sub>10</sub>(PML).<br>",
      "Values < 0 indicate estimates for which PML values are larger; values > 0 indicate estimates for which EA values are larger.<br>",
      "Values are based on estimates of carbon content (pg per cell) by each institution.<br>",
      "Differences in estimates are provided for the {nrow(df)} taxa which have carbon content estimates from both institutions.<br>",
      "The grey area indicates estimates which are within 1 order of magnitude of each other."
    )
  ) +
  theme(
    plot.caption = element_markdown(size=11),  # Enables Markdown rendering for the caption
    axis.text.x = element_text(face = 2),
    axis.title.x = element_text(face = 2),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(face=2)
  )

dev.off()
