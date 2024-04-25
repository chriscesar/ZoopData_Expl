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

tx_chk %>% dplyr::select(., Taxa, Aphia.ID) %>% 
  distinct() -> tx_chktrm

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

###############
## summarise by counts & plot cloud
## LF
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

png(file = "figs/wordcloud_LF_v3.png",
    width=16*ppi, height=10*ppi, res=ppi)
print(pl_LF_wc)
dev.off();toc()

### taxon version
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

png(file = "figs/wordcloud_tx_v3.png",
    width=16*ppi,
    height=10*ppi,
    res=ppi)
print(pl_tx_wc)
dev.off();toc()

### prep taxon data ####
### remove odd data
df_tx <- as_tibble(df_tx0) ### create new data (keep df0 as 'raw')

### convert dates
df_tx$sample.date <- as.Date(df_tx$sample.date, origin = "1899-12-30")

# Remove 100 µm data [OPTIONAL] ####
df_tx_100um <- df_tx %>% 
  filter(str_starts(Sample.comments,"100um"))

df_tx %>% 
  filter(!str_starts(Sample.comments,"100um")) -> df_tx

###widen data & fill NAs with 0s ####
df_tx %>% 
  dplyr::select(.,-c("Aphia.ID","AbundanceRaw","Taxa","Category":"Unallocated",
                     "LF0","Kingdom":"Subspecies")) %>% #drop unneeded cols
  group_by(across(c(-Abund_m3))) %>% # group by everything except abundance
  summarise(Abund_m3=sum(Abund_m3), #sum abundances
            .groups="drop") %>% 
  pivot_wider(names_from = "LF02",values_from = "Abund_m3", #widen
              values_fill = 0) -> df_tx_w

### join & save data ####  
dfw <- left_join(df_tx_w,df_wims_w,by="PRN")

dftaxa0 <- df_tx %>% 
  dplyr::select(.,Kingdom:DisplayName,Taxa) %>% 
  distinct() %>% 
  mutate_all(as.character) %>%
  mutate_all(~ replace_na(., "NA"))
toc()