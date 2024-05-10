### imp.taxTree.R ####
#### Plot taxonomic tree for fun! ####

#### load packages ####
ld_pkgs <- c("tidyverse","MASS","lubridate","vegan","mvabund","seas","ape","ggtree",
             "ecoCopula","performance","gclus","corrplot", "patchwork","gllvm")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

#### set universals ####
source("R/folder.links.R") ## data folders
perms <- 999 ### number of permutations to run for multivariate analyses
theme_set(ggthemes::theme_few())###set theme for all ggplot objects
nit <- 200 #number of iterations
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
### taxon data
df_tx0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,
                                               "processedData/MBA_Returns_Amalgamated_USE.xlsx"),
                                        sheet="outR04_LF"))
### WIMS chemical data
# WIMS Extract based on:
# Materials = 2HZZ & 2IZZ; SMPT_TYPE = CD, CC, CE; Dates from 01/06/2022-present
# SMP_Code:
# 42100171, 42100174, 42100179, 45400826, 60510027, 73015085, 82510555,
# 82615055, 82615255, 88002837, 88007163, 88007172, 88025879, BE061099,
# E0001449, E0001450, E0004730, G0003532, G0003572, LC544405, LC560357,
# PTTR0026, WA560349, Y0004367, Y0017477, YC536426

df_wims0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,
                                                 # "/WIMS_Extract_WaterQuality_Zoop_Samples_230809.xlsx"),
                                                 # "/WIMS_Extract_WaterQuality_Zoop_Samples_231218.xlsx"),
                                                 # "/WIMS_Extract_WaterQuality_Zoop_Samples_240108.xlsx"),
                                                 "/WIMS_Extract_WaterQuality_Zoop_Samples_240405.xlsx"),
                                          sheet="allDat"))

### prep WIMS data ####
### format & widen WIMS data ###
df_wims <- df_wims0

df_wims$PRN <- df_wims$SAMP_SCHEDULE_SAMPLE_ID
df_wims %>%
  dplyr::mutate(det=paste0(DETE_DESC,"_",UNIT_SHORT_DESC)) %>% ##create new variable label
  dplyr::mutate(Result=ifelse(is.na(df_wims$MEAS_SIGN == "<"), MEAS_RESULT,
                              paste0(MEAS_SIGN,MEAS_RESULT))) %>% #add "<" to results
  dplyr::select(.,c(WIMS.Code,REGION,Biosys.ID,
                    SMPT_LONG_NAME, SAMP_SAMPLE_DATE, SAMP_SAMPLE_TIME,
                    SAMP_Notes,
                    PRN,det,Result)) %>% ##only keep variables of interest
  tidyr::pivot_wider(.,names_from=det, values_from = Result) -> df_wims_w###widen data

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

# View(dftaxa)
# dftaxa <- dftaxa0
# 
# dftaxa$Kingdom <- factor(dftaxa$Kingdom)
# dftaxa$Phylum <- factor(dftaxa$Phylum)
# dftaxa$Class <- factor(dftaxa$Class)
# dftaxa$Order <- factor(dftaxa$Order)
# dftaxa$Family <- factor(dftaxa$Family)
# dftaxa$Genus <- factor(dftaxa$Genus)
# dftaxa$Species <- factor(dftaxa$Species)
# dftaxa$Taxa <- factor(dftaxa$Taxa)

# dftaxa_spp <- dftaxa[!is.na(dftaxa$Species),]
# dftaxa_spp %>% dplyr::select(.,-Subgenus) -> dftaxa_spp
# dftaxa_spp0 <- dftaxa[is.na(dftaxa$Species),]
# dftaxa_spp0 %>% 
#   dplyr::select(.,-Subgenus) %>% 
#   rowwise() %>% 
#   mutate(Species = last(na.omit(c_across()))) %>%
#   as.data.frame(.) %>% distinct() -> dftaxa_sppx
# 
# names(dftaxa_spp);names(dftaxa_sppx)
# 
# dftaxa_spp <- rbind(dftaxa_spp,dftaxa_sppx)
# dftaxa_spp <- dftaxa_spp %>% 
#   distinct()

## assign GROUPS for colouring
# #dftaxa$GROUP <- 
# 
# ### plot!
# taxonomicTree = as.phylo(~Kingdom/Phylum/Class/Order/Family/Genus/Species/Taxa,
#                          # data = dftaxa_spp,
#                          data = dftaxa,
#                          collapse = FALSE)
# taxonomicTree$edge.length = rep(1,length(taxonomicTree$edge))
# plot(taxonomicTree, cex=0.5)


### ver 2 ####
dftaxa0 <- df_tx %>% 
  dplyr::select(.,Kingdom:Species,DisplayName, LF02) %>% 
  distinct() %>% 
  mutate_all(as.character) %>%
  mutate_all(~ replace_na(., "NA"))

dftaxa <- dftaxa0

# assign groups
dftaxa$GROUP <- ifelse(dftaxa$Kingdom=="Chromista","Chromista",
                       ifelse(grepl("^Fish",dftaxa$LF02),"Fish",
                              ifelse(grepl("^Cop",dftaxa$LF02),"Copepod",
                                     ifelse(grepl("^Bryo",dftaxa$LF02),"Bryozoan",
                                            ifelse(grepl("^Crus",dftaxa$LF02),"Arthropoda",
                                                   ifelse(grepl("^Clad",dftaxa$LF02),"Arthropoda",
                                                   ifelse(grepl("^Gast",dftaxa$LF02),"Gastropod",
                                                          ifelse(grepl("^Poly",dftaxa$LF02),"Annelid",
                                                                 ifelse(grepl("^Gela",dftaxa$LF02),"Gelatinous",
                                                                        ifelse(grepl("^Acari",dftaxa$LF02),"Arthropoda","Other")
                                                                 )))))))))
dftaxa$cols <- ifelse(dftaxa$GROUP=="Annelid", "burlywood4",
                      ifelse(dftaxa$GROUP=="Arthropoda", "blueviolet",
                             ifelse(dftaxa$GROUP=="Bryozoan", "sienna",
                                    ifelse(dftaxa$GROUP=="Chromista", "dodgerblue",
                                           ifelse(dftaxa$GROUP=="Copepod", "chartreuse",
                                                  ifelse(dftaxa$GROUP=="Fish", "cyan3",
                                                         ifelse(dftaxa$GROUP=="Gastropod", "aquamarine2",
                                                                ifelse(dftaxa$GROUP=="Gelatinous", "deeppink2",
                                                                       ifelse(dftaxa$GROUP=="Other", "darkgrey","")
                                                                ))))))
                      ))
dftaxa <- dftaxa %>% 
  dplyr::select(.,-LF02) %>% 
  distinct()

dftaxa$Kingdom <- factor(dftaxa$Kingdom)
dftaxa$Phylum <- factor(dftaxa$Phylum)
dftaxa$Class <- factor(dftaxa$Class)
dftaxa$Order <- factor(dftaxa$Order)
dftaxa$Family <- factor(dftaxa$Family)
dftaxa$Genus <- factor(dftaxa$Genus)
dftaxa$Species <- factor(dftaxa$Species)
dftaxa$DisplayName <- factor(dftaxa$DisplayName)

taxonomicTree = as.phylo(~Kingdom/Phylum/Class/Order/Family/Genus/Species/DisplayName,
                         data = dftaxa,
                         collapse = FALSE)
tbl_taxonomicTree <- as_tibble(taxonomicTree)
taxonomicTree$edge.length = rep(1,length(taxonomicTree$edge))

xx <- taxonomicTree$tip.label
xz <- dftaxa[,c(9,11)]
dfcol <- data.frame(DisplayName=xx)
xy <- left_join(dfcol,xz,by="DisplayName")

plot(taxonomicTree, cex=0.5, tip.color = xy$cols)

# col.grp <- dftaxa[,c("Taxa", "GROUP")]
# cols <- ifelse(col.grp$GROUP == "Annelid", "burlywood4",
#                ifelse(col.grp$GROUP == "Arthropoda", "blueviolet",
#                       ifelse(col.grp$GROUP == "Bryozoan", "sienna",
#                              ifelse(col.grp$GROUP == "Chromista", "dodgerblue",
#                                     ifelse(col.grp$GROUP == "Copepod","chartreuse",
#                                            ifelse(col.grp$GROUP == "Fish","cyan3",
#                                                   ifelse(col.grp$GROUP == "Gastropod","darkolivegreen2",
#                                                          ifelse(col.grp$GROUP == "Gelatinous","darkorchid",
#                                                                 ifelse(col.grp$GROUP == "Other","darkgrey","")
#                                                          ))))))))

png(file = "figs/zooptaxa.png",
    width=8*ppi, height=16*ppi, res=ppi)
plot(taxonomicTree, cex=0.5, tip.color = xy$cols)
dev.off()


# un-needed ####

# dftaxa %>% 
#   rowwise() %>% 
#   dplyr::select(.,Kingdom:Genus) %>% 
#   mutate(Species = last(na.omit(c_across()))) %>%
#   ungroup() %>% 
#   mutate(Species = factor(Species)) %>% 
#   as.data.frame(.) %>% 
#   distinct() -> dftaxa2
# 
# write.csv(dftaxa2,file = "./dftaxa2.csv")
# 
# dftaxaTrm <- dftaxa[complete.cases(dftaxa),]
# 
# # CONSTRUCT TAXONOMICAL TREE TO BE USED AS PROXY FOR PHYLOGENETIC TREE
# # taxonomicTree = as.phylo(~Kingdom/Phylum/Class/Order/Family/Genus,
# #                          data = dftaxaTrm, collapse = FALSE)
# # taxonomicTree$edge.length = rep(1,length(taxonomicTree$edge))
# # plot(taxonomicTree, cex=0.5)
# # ggtree((as.phylo(~Kingdom/Phylum/Class/Order/Family/Genus,
# #                               data = dftaxaTrm, collapse = FALSE)))
# 
# taxonomicTree = as.phylo(~Kingdom/Phylum/Class/Order/Family/Genus/Species,
#                          data = dftaxa2, collapse = TRUE)
# taxonomicTree$edge.length = rep(1,length(taxonomicTree$edge))
# plot(taxonomicTree, cex=0.5)
# 
# taxa <- as.phylo(~Kingdom/Phylum/Class/Order/Species, data = dat)
# 
# col.grp <- merge(data.frame(Species = taxa$tip.label), dat[c("Species", "Group")], by = "Species", sort = F)
# 
# cols <- ifelse(col.grp$Group == "Benthos", "burlywood4", ifelse(col.grp$Group == "Zooplankton", "blueviolet", ifelse(col.grp$Group == "Fish", "dodgerblue", ifelse(col.grp$Group == "Phytoplankton", "darkolivegreen2", ""))))
# 
# plot(taxa,, tip.col = cols)
# 
# # 
# # # Reshape the data into long format
# # taxadf_long <- gather(taxadf, Taxonomic_Level, Name, -Species)
# # 
# # # Replace NA values in Species column with empty strings
# # taxadf_long$Species[is.na(taxadf_long$Species)] <- ""
# # 
# # # Create a new column "Parent" which indicates the parent node for each taxonomic level
# # taxadf_long$Parent <- NA
# # taxadf_long$Parent[!is.na(taxadf_long$Species)] <- taxadf_long$Genus[!is.na(taxadf_long$Species)]
# # taxadf_long$Parent[is.na(taxadf_long$Parent) & !is.na(taxadf_long$Genus)] <- taxadf_long$Family[is.na(taxadf_long$Parent) & !is.na(taxadf_long$Genus)]
# # taxadf_long$Parent[is.na(taxadf_long$Parent) & !is.na(taxadf_long$Family)] <- taxadf_long$Order[is.na(taxadf_long$Parent) & !is.na(taxadf_long$Family)]
# # taxadf_long$Parent[is.na(taxadf_long$Parent) & !is.na(taxadf_long$Order)] <- taxadf_long$Class[is.na(taxadf_long$Parent) & !is.na(taxadf_long$Order)]
# # taxadf_long$Parent[is.na(taxadf_long$Parent) & !is.na(taxadf_long$Class)] <- taxadf_long$Phylum[is.na(taxadf_long$Parent) & !is.na(taxadf_long$Class)]
# # taxadf_long$Parent[is.na(taxadf_long$Parent) & !is.na(taxadf_long$Phylum)] <- taxadf_long$Kingdom[is.na(taxadf_long$Parent) & !is.na(taxadf_long$Phylum)]
# # 
# # taxadf_long$Taxonomic_Level <- factor(taxadf_long$Taxonomic_Level,
# #                                       levels = c("Kingdom", "Phylum", "Class",
# #                                                  "Order", "Family","Genus"))
# # 
# # ggplot(taxadf_long, aes(x = Name, y = Taxonomic_Level, group = Name)) +
# #   geom_point() +
# #   geom_segment(aes(xend = Parent, yend = Taxonomic_Level)) +
# #   theme_minimal() +
# #   theme(axis.text.y = element_text(hjust = 1))+
# #   coord_flip()
# 
# ###
# # see:
#   # https://stackoverflow.com/questions/9904361/making-simple-phylogenetic-dendrogram-tree-from-a-list-of-species
# 
# # Group <- c("Benthos","Benthos","Benthos","Benthos","Benthos","Benthos","Zooplankton","Zooplankton","Zooplankton","Zooplankton",
# #            "Zooplankton","Zooplankton","Fish","Fish","Fish","Fish","Fish","Fish","Phytoplankton","Phytoplankton","Phytoplankton","Phytoplankton")
# # Domain <- rep("Eukaryota", length(Group))
# # Kingdom <- c(rep("Animalia", 18), rep("Chromalveolata", 4))
# # Phylum <- c("Annelida","Annelida","Arthropoda","Arthropoda","Porifera","Sipunculida","Arthropoda","Arthropoda","Arthropoda",
# #             "Arthropoda","Echinoidermata","Chorfata","Chordata","Chordata","Chordata","Chordata","Chordata","Chordata","Heterokontophyta",
# #             "Heterokontophyta","Heterokontophyta","Dinoflagellata")
# # Class <- c("Polychaeta","Polychaeta","Malacostraca","Malacostraca","Demospongiae","NA","Malacostraca","Malacostraca",
# #            "Malacostraca","Maxillopoda","Ophiuroidea","Actinopterygii","Chondrichthyes","Chondrichthyes","Chondrichthyes","Actinopterygii",
# #            "Actinopterygii","Actinopterygii","Bacillariophyceae","Bacillariophyceae","Prymnesiophyceae","NA")
# # Order <- c("NA","NA","Amphipoda","Cumacea","NA","NA","Amphipoda","Decapoda","Euphausiacea","Calanioda","NA","Gadiformes",
# #            "NA","NA","NA","NA","Gadiformes","Gadiformes","NA","NA","NA","NA")                     
# # Species <- c("Nephtys sp.","Nereis sp.","Gammarus sp.","Diastylis sp.","Axinella sp.","Ph. Sipunculida","Themisto abyssorum","Decapod larvae (Zoea)",
# #              "Thysanoessa sp.","Centropages typicus","Ophiuroidea larvae","Gadus morhua eggs / larvae","Etmopterus spinax","Amblyraja radiata",
# #              "Chimaera monstrosa","Clupea harengus","Melanogrammus aeglefinus","Gadus morhua","Thalassiosira sp.","Cylindrotheca closterium",
# #              "Phaeocystis pouchetii","Ph. Dinoflagellata")   
# # dat <- data.frame(Group, Domain, Kingdom, Phylum, Class, Order, Species)
# # 
# # dat
# # 
# # dat$Group <- as.factor(dat$Group)
# # dat$Domain <- as.factor(dat$Domain)
# # dat$Kingdom <- as.factor(dat$Kingdom)
# # dat$Phylum <- as.factor(dat$Phylum)
# # dat$Class <- as.factor(dat$Class)
# # dat$Order <- as.factor(dat$Order)
# # dat$Species <- as.factor(dat$Species)
# # 
# # 
# # taxa <- as.phylo(~Kingdom/Phylum/Class/Order/Species, data = dat)
# # 
# # col.grp <- merge(data.frame(Species = taxa$tip.label), dat[c("Species", "Group")], by = "Species", sort = F)
# # 
# # cols <- ifelse(col.grp$Group == "Benthos", "burlywood4", ifelse(col.grp$Group == "Zooplankton", "blueviolet", ifelse(col.grp$Group == "Fish", "dodgerblue", ifelse(col.grp$Group == "Phytoplankton", "darkolivegreen2", ""))))
# # 
# # plot(taxa,, tip.col = cols)

# tidy up ####
### unload packages
detach("package:tidyverse", unload=TRUE)
detach("package:vegan", unload=TRUE)
detach("package:gllvm", unload=TRUE)
detach("package:ecoCopula", unload=TRUE)
detach("package:mvabund", unload=TRUE)
detach("package:performance", unload=TRUE)
detach("package:patchwork", unload=TRUE)
detach("package:ggtree", unload=TRUE)
detach("package:corrplot", unload=TRUE)
detach("package:ape", unload=TRUE)
detach("package:gclus", unload=TRUE)
detach("package:seas", unload=TRUE)
detach("package:lubridate", unload=TRUE)
detach("package:MASS", unload=TRUE)

### unload objects
rm(list = ls(pattern = "^df"))
rm(list = ls(pattern = "^cbP"))
rm(list = ls(pattern = "^x"))
rm(list = ls(pattern = "taxon"))
rm(datfol,nit,perms,ppi)
