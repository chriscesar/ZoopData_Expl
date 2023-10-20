# fishyGoingsOn_ver2.R ####
## analysis of fish eggs and larvae
## incorporating zero fish counts

#### load packages ####
ld_pkgs <- c("tidyverse")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

### set up folders & import functions ###
source("R/folder.links.R")

theme_set(ggthemes::theme_few())###set theme for all ggplot objects
nit <- 9999 #number of iterations
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

### import data ####
## data are in LONG format with *tentative* lifeform values attached
df0 <- as_tibble(readxl::read_excel(paste0(datfol,
                                           "processedData/",
                                           "MBA_Returns_Amalgamated_USE.xlsx"),
                                    sheet = "outR04_LF"))

# convert dates
df0$`sample date` <- as.Date(df0$`sample date`, origin="1899-12-30")

## flag 'fishy' rows
# any life form flag in LF0 that contains the string "fish"
# this outputs a TRUE/FALSE against each row in the data
df0$LF0fish <- grepl("fish", df0$LF0, ignore.case = TRUE)

#### if fish == TRUE, assign "Eggs" to eggy lads and "Larvae" to non-eggsters
df0$fish_type <- ifelse(df0$LF0fish & grepl("egg",
                                            df0$Taxa,
                                            ignore.case = TRUE),
                        "Fish eggs",
                        ifelse(df0$LF0fish, "Fish larvae", "FALSE"))

### incorporate fish info into tax list
df0$taxUSE <- ifelse(df0$fish_type == "FALSE", df0$Taxa, df0$fish_type)

### remove unnesseccary cols
dfl <- df0 %>% 
  dplyr::select(., c(`Pot Number`, `sample date`, `BIOSYS Code`,
                     `WIMS Code`, PRN, `Sample comments`,
                     taxUSE, Abund_m3:Category))
## widen
dfl %>% 
  group_by(across(!Abund_m3)) %>% 
  summarise(Abund_m3 = sum(Abund_m3,
                           na.rm = TRUE),
            .groups = "drop") %>%  # sum fishy things
  pivot_wider(names_from = taxUSE, values_from = Abund_m3,
              values_fill = 0) %>%  ##widen 
  dplyr::select(.,c(`Pot Number`:Category,
                    `Fish eggs`,
                    `Fish larvae`)) %>% #get rid of unneccessary taxa
  pivot_longer(cols = c(`Fish eggs`, `Fish larvae`),
               names_to = "fish_type",
               values_to = "Abund_m3") -> dfl
  

### create shortened WB label
dfl$WBlb <- ifelse(
  dfl$WB == "Solent", "Solent",
  ifelse(
    dfl$WB == "SOUTHAMPTON WATER", "Soton Wtr",
    ifelse(
      dfl$WB == "THAMES LOWER", "Thm Low",
      ifelse(
        dfl$WB == "Blackwater Outer", "Blckw Out",
        ifelse(
          dfl$WB == "Cornwall North", "Cornw Nth",
          ifelse(
            dfl$WB == "Barnstaple Bay", "Brnst B",
            ifelse(
              dfl$WB == "Kent South", "Kent Sth",
              ifelse(
                dfl$WB == "Mersey Mouth", "Mersey Mth",
                ifelse(
                  dfl$WB == "Wash Outer", "Wash Out",
                  ifelse(
                    dfl$WB == "Lincolnshire", "Lincs",
                    ifelse(
                      dfl$WB == "Yorkshire South", "Yorks Sth",
                      ifelse(
                        dfl$WB == "TEES", "Tees",
                        ifelse(
                          dfl$WB == "Northumberland North", "Nrthmb Nth",
                          ifelse(
                            dfl$WB == "Farne Islands to Newton Haven",
                            "Farne Is",
                            ifelse(
                              dfl$WB == "Bristol Channel Inner South",
                              "Brist Ch In Sth",NA
                            )))))))))))))))
dfl$RegSh <- ifelse(dfl$Region == "Southern", "Sth",
                    ifelse(dfl$Region == "Anglian", "Ang",
                           ifelse(dfl$Region == "SWest", "SW",
                                  ifelse(dfl$Region == "NWest", "NW",
                                         ifelse(dfl$Region == "NEast", "NE",
                                                ifelse(dfl$Region == "Thames", "Thm",NA
                                                       ))))))

dfl$RgWB <- paste0(dfl$RegSh,"_",dfl$WBlb)


# plot boxplot
png(file = "figs/fishBoxplot.png",
    width=18*ppi, height=9*ppi, res=ppi)
dfl %>% 
ggplot(., aes(x = as.factor(RgWB), y=Abund_m3, fill=as.factor(fish_type)))+
  geom_hline(yintercept = 0, col="grey")+
  geom_boxplot(outlier.shape = NA, varwidth = TRUE)+
  geom_jitter(aes(shape = fish_type),
              position = position_jitterdodge(),
              alpha=0.3,
              show.legend = FALSE)+
  labs(title = "Fish larvae and egg abundances by WFD water body",
       subtitle = bquote("Values indicate mean recorded abundances " ~m^-3),
       y = bquote("Abundance "~(m^-3)),
       caption = paste0("Samples gathered between ",min(df0$`sample date`),
                        " & ",max(df0$`sample date`),".","\nBox widths are proportional to the number of samples.
                        Individual sample points overlain (Circles = fish eggs; Triangles = fish larvae)."))+
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position="bottom")
dev.off()