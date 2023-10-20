# fishyGoingsOn.R ####
## analysis of fish eggs and larvae

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

### remove unneccessary col
df0 <- df0 %>% 
  dplyr::select(., !date_site)

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
# 
# ### create smaller data for summarising
# wbs <- unique(dfsumm$WB) ## water body names

dfsummary <- df0 %>%
  dplyr::select(.,
                c(`Pot Number`:Category,
                  fish_type)) %>% # retain only pertinent variables
  filter(., fish_type != FALSE) %>% #retain only data flagged as fishy
  dplyr::select(.,!c(Taxa,
                     AbundanceRaw,
                     `Aphia ID`)) %>% ## drop unneccessary
  group_by(across(c(!fish_type,!Abund_m3))) %>% # group by variables of interest
  summarise(Abund_m3 = sum(Abund_m3,
                           na.rm = TRUE),
            .groups = "drop") # sum fishy things

### create shortened WB label
dfsummary$WBlb <- ifelse(
  dfsummary$WB == "Solent", "Solent",
  ifelse(
    dfsummary$WB == "SOUTHAMPTON WATER", "Soton Wtr",
    ifelse(
      dfsummary$WB == "THAMES LOWER", "Thm Low",
      ifelse(
        dfsummary$WB == "Blackwater Outer", "Blckw Out",
        ifelse(
          dfsummary$WB == "Cornwall North", "Cornw Nth",
          ifelse(
            dfsummary$WB == "Barnstaple Bay", "Brnst B",
            ifelse(
              dfsummary$WB == "Kent South", "Kent Sth",
              ifelse(
                dfsummary$WB == "Mersey Mouth", "Mersey Mth",
                ifelse(
                  dfsummary$WB == "Wash Outer", "Wash Out",
                  ifelse(
                    dfsummary$WB == "Lincolnshire", "Lincs",
                    ifelse(
                      dfsummary$WB == "Yorkshire South", "Yorks Sth",
                      ifelse(
                        dfsummary$WB == "TEES", "Tees",
                        ifelse(
                          dfsummary$WB == "Northumberland North", "Nrthmb Nth",
                          ifelse(
                            dfsummary$WB == "Farne Islands to Newton Haven",
                            "Farne Is",
                            ifelse(
                              dfsummary$WB == "Bristol Channel Inner South",
                              "Brist Ch In Sth",NA
                            )))))))))))))))
dfsummary$RegSh <- ifelse(dfsummary$Region == "Southern", "Sth",
                          ifelse(dfsummary$Region == "Anglian", "Ang",
                                 ifelse(dfsummary$Region == "SWest", "SW",
                                        ifelse(dfsummary$Region == "NWest", "NW",
                                               ifelse(dfsummary$Region == "NEast", "NE",
                                                      ifelse(dfsummary$Region == "Thames", "Thm",NA
                                                      ))))))

dfsummary$RgWB <- paste0(dfsummary$RegSh,"_",dfsummary$WBlb)

### calculate mean and SD by WB
dfsummary %>% 
  group_by(RgWB,
           fish_type) %>% 
  summarise(n = n(),
            mean = mean(Abund_m3,
                        na.rm = TRUE),
            sd = sd(Abund_m3,
                    na.rm = TRUE),.groups = "drop") -> dfplot

dfnu <- data.frame("RgWB" = "SW_Brist Ch In Sth",
                   "fish_type" = "Fish larvae",
                   n = 0, mean = 0, sd=NA)
dfplot <- rbind(dfplot, dfnu);rm(dfnu)

png(file = "figs/fishMeans.png",
    width=18*ppi, height=9*ppi, res=ppi)
ggplot(dfplot, aes(x=as.factor(RgWB), y=mean, fill=as.factor(fish_type))) +
  # geom_vline(xintercept = 5)+
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  coord_cartesian(ylim = c(0,NA))+
  labs(title = "Fish larvae and egg abundances by WFD water body",
       subtitle = bquote("Values indicate mean recorded abundances " ~m^-3~" ± standard deviation"),
       y = bquote("Abundance "~(m^-3)),
       caption = paste0("Samples gathered between ",min(df0$`sample date`),
                        " & ",max(df0$`sample date`)))+
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position="bottom",
        axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0))
dev.off()

## as boxplot:
## append 'empty value'
x <- dfsummary[dfsummary$RgWB=="SW_Brist Ch In Sth",]
x$Abund_m3 <- 0;x$fish_type <- "Fish larvae"

dfsummary <- rbind(dfsummary,x)

dfsummary %>% 
  ggplot(., aes(x = as.factor(RgWB), y=Abund_m3, fill=as.factor(fish_type)))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(shape = fish_type),
              position = position_jitterdodge(),
              alpha=0.3,
              show.legend = FALSE)+
  labs(title = "Fish larvae and egg abundances by WFD water body",
       subtitle = bquote("Values indicate mean recorded abundances " ~m^-3~" ± standard deviation"),
       y = bquote("Abundance "~(m^-3)),
       caption = paste0("Samples gathered between ",min(df0$`sample date`),
                        " & ",max(df0$`sample date`)))+
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position="bottom")

# dfsum <- dfsummary %>% 
#   dplyr::select(., !c(Taxa, AbundanceRaw, `Aphia ID`)) %>% 
#   group_by(across(c(!fish_type, !Abund_m3))) %>% 
#   summarise(Abund_m3 = sum(Abund_m3, na.rm = TRUE), .groups = "drop")

# create total Fish Egg and total Fish larvae values per sample
# calculate mean egg and larvae values by WB ± sd

### to do ####
# flag 'fish' type rows - DONE
# assign variable "type" to indicate whether it's an EGG or LARVAE - DONE
# (think about what to do in cases where just the taxon name is provided,
# e.g. "Gobiidae") - DONE
# Summarise counts by sample (i.e. Total larvae and Total egg values per sample)
# calculate mean values ± SD by WB
# present, e.g., stacked larvae-egg barchart by WB arranged by region?