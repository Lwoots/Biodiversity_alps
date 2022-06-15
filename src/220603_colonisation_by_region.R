#Script to connect species to their global affinities

#Packages
library(tidyverse)
library(here)
library(ape)

#Data###

data <- read.csv(here("Data", "220530_distribution_data_androsace.csv"), sep = ",")

affinities <- data %>% 
  mutate(Species = str_extract(libelle, "[A-Za-z]+ [a-z\\-]+"))

fa_data <- read.csv(here("Data", "220608_full_bedrock_data.csv"), sep = ";")

tree <- read.tree("/Users/larawootton/Documents/Doctorate/Tree_dating/Tree_files/220301_treePL_optimisation_run2_final.tre")
final_taxonset <- read.csv(here("Data", "220419_taxonset_with_tips.csv"))
sample_classification <- read.csv(here("Data", "220421_coding_endemics_completed.csv"), sep = ";")
ana <- read.csv(here("Data", "220517_anagenetic_species.csv"))

#Wrangling ####
island_species <- sample_classification %>% 
  filter(!Status %in% c("Mainland", "Endemic"))

island_species <- sample_classification %>% 
  filter(!Status == "Mainland")

island_affinities <- left_join(island_species, affinities, by = c("EuroMedAcceptedName_formatted" = "Species"))

no_data <- island_affinities %>% 
  filter(is.na(code))

#Assign each 
colonists <- island_affinities %>% 
  filter(Status == "Non_endemic_MaxAge" | EuroMedAcceptedName_formatted %in% ana$EuroMedAcceptedName_formatted ) %>% 
  #filter(!is.na(code)) %>% 
  mutate(Area = case_when(
    str_detect(nom, "A1") ~ "South America",
    str_detect(nom, "A2") ~ "North America",
    str_detect(nom, "A3") ~ "Africa",
    str_detect(nom, "A4") ~ "Asia",
    str_detect(nom, "A5") ~ "Australia",
    str_detect(nom, "B1") ~ "Warm zones",
    str_detect(nom, "B2") ~ "Cold zones",
    str_detect(nom, "B3") ~ "Eurasia", #Eurasia, N. am
    str_detect(nom, "B4") ~ "Eurasia", 
    str_detect(nom, "B5") ~ "Arctic",
    str_detect(nom, "B6") ~ "Eurasia",
    str_detect(nom, "C1") ~ "Mediterranean",
    str_detect(nom, "C2") ~ "Mediterranean",
    str_detect(nom, "D1d") ~ "Mediterranean",
    str_detect(nom, "D1e") ~ "Mediterranean",
    str_detect(nom, "D1f") ~ "Mediterranean",
    str_detect(nom, "D3") ~ "Europe", #check for D1d = canary islands and D1f = Sardinia
    str_detect(nom, "D5") ~ "Arctic",
    str_detect(nom, "D6") ~ "Arctic",
    str_detect(nom, "E") ~ "EAS", #Europe mountains
    str_detect(nom, "F") ~ "Alps",
  )) %>% 
  distinct(EuroMedAcceptedName_formatted, .keep_all = T)

ch <- colonists %>% 
  group_by(Area) %>% 
  count()

#Wrangling for flora alpina data

table(fa_data$Distribution)

fa_data <- fa_data %>% 
  mutate(Area = case_when(
    str_detect(Distribution, "sia") ~ "Eurasia",
    str_detect(Distribution, "edit") ~ "Mediterranean",
    str_detect(Distribution, "Eurosib") ~ "Arctic",
    str_detect(Distribution, "N-Am") ~ "NAmerican/Cosmopolitan",
    str_detect(Distribution, "osmop") ~ "NAmerican/Cosmopolitan",
    str_detect(Distribution, "Arct") ~ "Arctic",
    str_detect(Distribution, "Eur-Mont") ~ "EAS",
    str_detect(Distribution, "Euro-Mont") ~ "EAS",
    str_detect(Distribution, "Eur") ~ "Europe",
    str_detect(Distribution, "Alp") ~ "Alps",
  ))


all_dat <- left_join(colonists, fa_data, by = "Tip")

all_dat <- all_dat %>% 
  mutate(Area = case_when(
    is.na(Area.x) ~ Area.y,
    is.na(Area.y) ~ Area.x
  ))

#Get colonisation times
tip_ages <- setNames(tree$edge.length[sapply(1:length(tree$tip.label),function(x,y) which (y==x),y=tree$edge[,2])],tree$tip.label)

tip_branch_ages <- data.frame(Tip = names(tip_ages),
                              Branching_times = tip_ages)

colonists <- left_join(all_dat, tip_branch_ages, by = "Tip")

boxplot(colonists$Branching_times ~ colonists$Area, ylim = c(0,35))

col_filtered <- colonists %>% 
  filter(!Area %in% c("Warm zones", "Cold zones", "NAmerican/Cosmopolitan"),
         !is.na(Area))

ggplot(col_filtered, aes(x = Branching_times, colour = Area)) +
  geom_density() +
  xlim(0,35) +
  scale_color_brewer(palette = "Set1") +
  xlab("Time of arrival") +
  theme_bw()

ggplot(col_filtered, aes(y = Branching_times, colour = Area)) +
  geom_boxplot() +
  theme_bw()+
  ylim(0,35) +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set1")
  
  ggplot(col_filtered, aes(y = Branching_times, x = factor(Area), colour = Area)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    theme_bw()+
    theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set1") +
  ylim(0,35) 
