# an.EA_wordcloud.R ####
#### Generate wordckouds of life forms and taxon data

#### load packages ####
ld_pkgs <- c("tidyverse","ggwordcloud", "tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

#### set universals ####
tic()
source("R/folder.links.R") ## data folders
theme_set(ggthemes::theme_few())###set theme for all ggplot objects
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

#### load LIFEFORMS data ####
df_LF0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,
                                               "processedData/MBA_Returns_Amalgamated_USE.xlsx"),
                                        sheet="outR04_LF"))

#### load TAXON data ####
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

tx_chk %>%
  dplyr::select(., Taxa, Aphia.ID) %>% 
  distinct() -> tx_chktrm

tx_chk %>%
  dplyr::select(., -c(ReturnedTaxName)) %>%
  distinct() -> tx_chktrm2

df_tx0 <- left_join(df_tx0, tx_chktrm, by="Aphia.ID")
df_tx0$Taxa.x <- df_tx0$Taxa.y;df_tx0$Taxa.y <- NULL
df_tx0 %>%
  rename(Taxa=Taxa.x) %>% 
  dplyr::select(.,-AbundanceRaw) %>% 
  group_by(across(c(!Abund_m3))) %>% 
  summarise(.,Abund_m3=sum(Abund_m3),.groups = "drop") %>%
  ungroup() %>% 
  as_tibble(.) -> df_tx
rm(tx_chk,tx_chk0,tx_chktrm)
toc()

###############
## summarise by counts & plot cloud
# LF version ####
tic();set.seed(271);df_LF0 %>% 
  # group_by(Region) %>% 
  count(.,LF02) %>% 
  mutate(tot=sum(n)) %>% 
  mutate(prop=n/tot) %>% ungroup() %>% 
  ggplot(.,aes(label=LF02, size=prop))+
  # geom_text_wordcloud()+
  geom_text_wordcloud_area()+
  scale_size_area(max_size = 150,
                  trans = power_trans(1/.7)) + #trans = power_trans(1/.7) "better fit human area perception"
  # facet_wrap(.~Region, scales = "free")+
  # scale_size_area(max_size = 20) +
  theme_minimal()+
  theme(strip.text = element_text(face="bold",
                                  size=14,
                                  colour="red")) -> pl_LF_wc

png(file = "figs/wordcloud_LF_v3.5.png",
    width=16*ppi, height=10*ppi, res=ppi)
print(pl_LF_wc)
dev.off();toc()

# taxon version ####
# append taxon info
tic();set.seed(21);df_tx %>% 
  # group_by(Region) %>% 
  count(.,Taxa) %>% 
  mutate(tot=sum(n)) %>% 
  mutate(prop=n/tot) %>% ungroup() %>% 
  ggplot(.,aes(label=Taxa, size=prop))+
  # geom_text_wordcloud()+
  geom_text_wordcloud_area()+
  # scale_size_area(max_size = 50, trans = power_trans(1/.7)) +
  # facet_wrap(.~Region, scales = "free")+
  # scale_size_area(max_size = 20) +
  scale_size_area(max_size = 70,
                  trans = power_trans(1/.7)) + #trans = power_trans(1/.7) "better fit human area perception"
  theme_minimal()+
  theme(strip.text = element_text(face="bold",
                                  size=14,
                                  colour="red")) -> pl_tx_wc

png(file = "figs/wordcloud_tx_v3.5.png",
    width=16*ppi,
    height=10*ppi,
    res=ppi)
print(pl_tx_wc)
dev.off();toc()

################################################################
# WIP ####
# new version ###
#append colours
### append Group variable for colours:
df_LF0$GROUP <- ifelse(df_LF0$Kingdom=="Chromista","Chromista",
                      ifelse(grepl("^Fish",df_LF0$LF02),"Fish",
                             ifelse(grepl("^Cop",df_LF0$LF02),"Copepod",
                                    ifelse(grepl("^Bryo",df_LF0$LF02),"Bryozoan",
                                           ifelse(grepl("^Crus",df_LF0$LF02),"Arthropoda",
                                                  ifelse(grepl("^Clad",df_LF0$LF02),"Arthropoda",
                                                         ifelse(grepl("^Gast",df_LF0$LF02),"Gastropod",
                                                                ifelse(grepl("^Poly",df_LF0$LF02),"Annelid",
                                                                       ifelse(grepl("^Gela",df_LF0$LF02),"Gelatinous",
                                                                              ifelse(grepl("^Acari",df_LF0$LF02),"Arthropoda","Other")
                                                                       )))))))))
df_LF0$cols <- as.factor(ifelse(df_LF0$GROUP=="Annelid", "burlywood4",
                     ifelse(df_LF0$GROUP=="Arthropoda", "blueviolet",
                            ifelse(df_LF0$GROUP=="Bryozoan", "sienna",
                                   ifelse(df_LF0$GROUP=="Chromista", "dodgerblue",
                                          ifelse(df_LF0$GROUP=="Copepod", "chartreuse",
                                                 ifelse(df_LF0$GROUP=="Fish", "cyan3",
                                                        ifelse(df_LF0$GROUP=="Gastropod", "aquamarine2",
                                                               ifelse(df_LF0$GROUP=="Gelatinous", "deeppink2",
                                                                      ifelse(df_LF0$GROUP=="Other", "darkgrey","")
                                                               )))))))))
pal <- c("burlywood4","blueviolet","sienna","dodgerblue","chartreuse","cyan3",
                     "aquamarine2","deeppink2","darkgrey")


df_LF0 %>%
  dplyr::select(
    #Taxa,
    DisplayName,
    GROUP,
    cols) %>%
  distinct()->tmp_col

tic();set.seed(21);df_LF0 %>% 
  # group_by(Region) %>% 
  # count(.,Taxa) %>%
  count(.,DisplayName) %>% 
  mutate(tot=sum(n)) %>% 
  mutate(prop=n/tot) %>% ungroup() -> tmp

tmp <- left_join(tmp,tmp_col, by="DisplayName")
tmp %>% 
  ggplot(.,aes(
    label=DisplayName,
    size=prop,
    col=GROUP
    )
    )+
  # geom_text_wordcloud()+
  geom_text_wordcloud_area()+
  scale_color_discrete(pal)+
  # scale_size_area(max_size = 50, trans = power_trans(1/.7)) +
  # facet_wrap(.~Region, scales = "free")+
  # scale_size_area(max_size = 20) +
  scale_size_area(max_size = 70,
                  trans = power_trans(1/.7)) + #trans = power_trans(1/.7) "better fit human area perception"
  theme_minimal()+
  theme(strip.text = element_text(face="bold",
                                  size=14,
                                  colour="red")) -> pl_tx_wc

png(file = "figs/wordcloud_tx_v4.png",
    width=16*ppi,
    height=10*ppi,
    res=ppi)
print(pl_tx_wc)
dev.off();toc()

#########################
set.seed(21)
words <- c(rep("apple", 50), rep("banana", 30), rep("orange", 20))
groups <- c(rep("Group1", 50), rep("Group2", 30), rep("Group3", 20))
xdf <- data.frame(word = words, group = groups)

# Plot the word cloud
ggplot(xdf, aes(label = word, color = group, size = 1)) +
  geom_text_wordcloud_area() +
  scale_color_manual(values = c("Group1" = "blue", "Group2" = "green", "Group3" = "red")) +
  theme_minimal() -> pltest
png(file = "figs/wordcloud_TEST.png",
    width=16*ppi,
    height=10*ppi,
    res=ppi)
print(pltest)
dev.off()

### tidy up
rm(list = ls(pattern = "^df"))
rm(list = ls(pattern = "^pl"))
rm(list = ls(pattern = c("^cb")))
rm(list = ls(pattern = c("^tmp")))
rm(tx_chktrm2,xdf,datfol,groups,pal,pl,ppi,words)

detach("package:tidyverse", unload=TRUE)
detach("package:tictoc", unload=TRUE)
detach("package:ggwordcloud", unload=TRUE)
