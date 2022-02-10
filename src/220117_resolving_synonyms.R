#17 Jan 2022
#Script to update synonymy between datasets
# Objectives:
# 1) Update PhyloAlps data with correct names
# 2) Update Flora alpina data with correct names
# 3) Combine Flora alpina and PhyloAlps datasets and connect to tree tip names
# 4) Trim tree/datasets to include only Alps species and only one tip per species

#Set up ####

#Packages

library(tidyverse)
library(here)

#Data

phyloalps <- read.csv(here("Data", "PHYLOALPS_TAXO_62_20nov2021.csv"), sep = ";")
phyloalps <- phyloalps[,1:10]

euro_med <- read.csv(here("Data/PHYLOALPS", "PhyloAlpsDirectCorrespondances.csv"))

flora_alpina <- read.csv(here("Data", "Cristina_Etages_Provinces.csv"))

synonyms_euro <- read.csv(here("Data/PHYLOALPS", "SynonymyEuromed.csv"))
synonyms_taxref <- read.csv(here("Data/PHYLOALPS", "SynonymyTaxref.csv"))
synonyms_plantlist <- read.csv(here("Data/PHYLOALPS", "SynonymyPlantList.csv"))

#Subsp wrangling

#Next get subspecies and vars into same format across datasets

flora_alpina$Full_names <- str_c(flora_alpina$Species, flora_alpina$Sous_espece, sep = " ")
flora_alpina$Full_names <- str_trim(flora_alpina$Full_names, side = "right") #Gets rid of white space at end of text

synonyms_euro$EuroMedAccepted_formatted <- str_replace_all(synonyms_euro$EuroMedAccepted,
                                                           pattern = c(" var. " = " ",
                                                                       " subsp. " = " ")
)

synonyms_euro$EuroMedNonAccepted_formatted <- str_replace_all(synonyms_euro$EuroMedNonAccepted,
                                                              pattern = c(" var. " = " ",
                                                                          " subsp. " = " ")
)

synonyms_taxref$EuroMedAccepted_formatted <- str_replace_all(synonyms_taxref$EuroMedAcceptedName,
                                                             pattern = c(" var. " = " ",
                                                                         " subsp. " = " ")
)

synonyms_taxref$scientificName_formatted <- str_replace_all(synonyms_taxref$scientificName,
                                                            pattern = c(" var. " = " ",
                                                                        " subsp. " = " ")
)

synonyms_plantlist$EuroMedAccepted_formatted <- str_replace_all(synonyms_plantlist$EuroMedAcceptedName,
                                                                pattern = c(" var. " = " ",
                                                                            " subsp. " = " ")
)

synonyms_plantlist$PlantName_formatted <- str_replace_all(synonyms_plantlist$PlantListName,
                                                          pattern = c(" var. " = " ",
                                                                      " subsp. " = " ")
)

euro_med$EuroMedAcceptedName_formatted <- str_replace_all(euro_med$EuroMedAcceptedName,
                                                          pattern = c(" var. " = " ",
                                                                      " subsp. " = " ")
)


#Obj. 1) Update PhloAlps ####

#Check whether sequencing code is being read in the same way

setdiff(euro_med$Sequencing_ID, phyloalps$Sequencing_ID)
setdiff(phyloalps$Sequencing_ID,euro_med$Sequencing_ID) #Seems good



#Reduce euromed to accepted name and sequncing ID to prevent column.x etc when joining
euro_med_reduced <- euro_med %>% 
  select("EuroMedAcceptedName_formatted", "Sequencing_ID")

#Combine the with phyloalps
combined_em_pa <- left_join(euro_med_reduced, phyloalps, by = "Sequencing_ID")


#Obj. 2) Update Flora alpina ####

#Clean typos

flora_alpina[1817, 5] <- "Hippophae" 
flora_alpina[1817, 8] <- "Hippophae rhamnoides"

flora_alpina[4124:4126, 5] <- "Hierochloe"
flora_alpina[4124, 8] <- "Hierochloe australis" 
flora_alpina[4125, 8] <- "Hierochloe odorata"
flora_alpina[4126, 8] <- "Hierochloe arctica"

flora_alpina[4363, 8] <- "Narcissus poeticus"
flora_alpina[4363, 7] <- "poeticus"

flora_alpina[15:17, 5] <- "Isoetes"


flora_alpina[15, 8] <- "Isoetes lacustris"
flora_alpina[16, 8] <- "Isoetes echinospora"
flora_alpina[17, 8] <- "Isoetes malinverniana"


#Make loop that compares each flora alpina species to all the euromed synonyms, plus pulls out correct name
flora_alpina$Current_names <- rep(NA, length(flora_alpina$Full_names))

for (i in 1:length(flora_alpina$Full_names)) {
  for (j in 1:length(synonyms_euro$EuroMedNonAccepted)) {
       if (flora_alpina$Full_names[i] == synonyms_euro$EuroMedNonAccepted_formatted[j]) {
          flora_alpina$Current_names[i] <- synonyms_euro$EuroMedAccepted_formatted[j]
          print(i/length(flora_alpina$Full_names)*100)
       }
  }
}

#Try same loop but in the TaxRef database

flora_alpina$Current_names_taxref <- rep(NA, length(flora_alpina$Full_names))

for (i in 1:length(flora_alpina$Full_names)) {
  for (j in 1:length(synonyms_taxref$EuroMedAccepted_formatted)) {
    if (flora_alpina$Full_names[i] == synonyms_taxref$scientificName_formatted[j]) {
      flora_alpina$Current_names_taxref[i] <- synonyms_taxref$EuroMedAccepted_formatted[j]
      print(i/length(flora_alpina$Full_names)*100)
    }
  }
}

#Same loop but plant list

flora_alpina$Current_names_plantlist <- rep(NA, length(flora_alpina$Full_names))

for (i in 1:length(flora_alpina$Full_names)) {
  for (j in 1:length(synonyms_plantlist$EuroMedAccepted_formatted)) {
    if (flora_alpina$Full_names[i] == synonyms_plantlist$PlantName_formatted[j]) {
      flora_alpina$Current_names_plantlist[i] <- synonyms_plantlist$EuroMedAccepted_formatted[j]
      print(i/length(flora_alpina$Full_names)*100)
    }
  }
}

#Create column with up to date names

flora_alpina <- flora_alpina %>% 
  mutate(Combined_names = case_when(
    !is.na(Current_names) ~ Current_names,
    is.na(Current_names) & !is.na(Current_names_taxref) ~ Current_names_taxref,
    is.na(Current_names) & is.na(Current_names_taxref) & !is.na(Current_names_plantlist) ~ Current_names_plantlist,
    is.na(Current_names) & is.na(Current_names_taxref) & is.na(Current_names_plantlist) ~ Full_names
  ))

#In case of inconsistencies between datasets, use euromed name

#However, there are lots of flora alpina species that aren't in either synonyms or accepted names

#not_in_accepted <- setdiff(flora_alpina$Full_names, synonyms_euro$EuroMedAccepted_formatted)
#not_in_synonyms <- setdiff(flora_alpina$Full_names, synonyms_euro$EuroMedNonAccepted_formatted)
#
#not_in_either <- intersect(not_in_accepted, not_in_synonyms)
#
#not_in_accepted_tax <- setdiff(flora_alpina$Full_names, synonyms_taxref$EuroMedAccepted_formatted)
#not_in_synonyms_tax <- setdiff(flora_alpina$Full_names, synonyms_taxref$scientificName_formatted)
#
#not_in_either_tax <- intersect(not_in_accepted, not_in_synonyms)
#
#intersect(not_in_either_tax, not_in_either)
#
#setdiff(flora_alpina$Combined_names, synonyms_euro$EuroMedAccepted_formatted)
#setdiff(flora_alpina$Combined_names, synonyms_plantlist$EuroMedAccepted_formatted)

#Let's join the corrected names in flora alpina to corrected names in phyloalps to see what we still need

test <- left_join(combined_em_pa, flora_alpina, by = c("EuroMedAcceptedName_formatted" = "Combined_names"))
#test <- test[,c(1:6, 12, 82:86)] #For ease of visualisation

#Now get rid of samples from non-phylo alps projects, in order to have mainly Alpine species

alps_only <- test %>% 
  filter(Project_name == "PhyloAlps")

#Pull out names that don't have flora alpina data

missing_data <- alps_only %>% 
  filter(is.na(Full_names)) #Still 764 entries that aren't connected to flora alpina

missing_data <- missing_data %>% 
  distinct(EuroMedAcceptedName_formatted, .keep_all = T) #582 species that are missing data

odds <- intersect(missing_data$EuroMedAcceptedName_formatted, flora_alpina$Full_names) #Strangely, 300 of the names can be found in the FA data set

#Whats going on?

fa_inconsistencies <- flora_alpina %>% 
                  filter(Full_names %in% odds) 
#fa_inconsistencies <- fa_inconsistencies[,70:75] #It seems that the same synonym in different files has the different reference files

#Okay, lets go with the direct correspondances file names then

flora_alpina <- flora_alpina %>% 
  mutate(Inconsistancies = 
           case_when(
             Full_names %in% odds ~ Full_names
           ))
flora_alpina <- flora_alpina %>% 
  mutate(DirectCorr = 
           case_when(
             !is.na(Inconsistancies) ~ Inconsistancies,
             is.na(Inconsistancies) ~ Combined_names
           ))

#Now recombine with phyloalps

all_dat <- left_join(combined_em_pa, flora_alpina, by = c("EuroMedAcceptedName_formatted" = "DirectCorr"))

alps_only <- all_dat %>% 
  filter(Project_name == "PhyloAlps")

missing_data <- alps_only %>% 
  filter(is.na(Full_names)) #Still 482 entries that aren't connected to flora alpina

missing_data <- missing_data %>% 
  distinct(EuroMedAcceptedName_formatted, .keep_all = T) #331 species that are missing data


more_species <- intersect(missing_data[,1], flora_alpina$Species) #But there's still 100 species that are in FA, mainly due to subspecies inconsistencies I think
#Ie it is a full subspecies in flora alpina but just a species in phyloalps

ms <- missing_data %>% 
  filter(!EuroMedAcceptedName_formatted %in% more_species) %>% 
  select(EuroMedAcceptedName_formatted)
write_csv(ms, file = "Missing_species.csv")

#So join on to the data set using species rather than all subspecies

more_missing <- flora_alpina %>% 
  filter(Species %in% more_species)

#get 

all_dat <- left_join(all_dat, flora_alpina, by = c("EuroMedAccepted_formatted" = "Species"))

alps_only <- all_dat %>% 
  filter(Project_name == "PhyloAlps")

missing_data <- alps_only %>% 
  filter(is.na(Full_names))


# Obj. 3) Connect to tree tip names ####

#Join phyloalps, flora alpina
initial_join <- left_join(combined_em_pa, flora_alpina, by = c("EuroMedAcceptedName_formatted" = "DirectCorr"))

#Get rid of columns that aren't super useful

selected_join <- initial_join %>% 
  select(!c(DB_ID, Sample_ID, Provider_name, Scientific_Name, Nom_Taxon.vérifié.ThePlantList,
            No_Flora_alpina, No_Famille, No_Genre, No_Espece, No_Sous_espece,
            Current_names, Current_names_taxref, Current_names_plantlist, Inconsistancies))

#Read in tree

library(ape)


#Obj. 4) Combine subspecies

#make small test dataset

subsp <- all_dat %>% 
  filter(Species %in% c("Aconitum variegatum", 
                        "Asplenium trichomanes", 
                        "Pulsatilla alpina",
                        "Minuartia cherlerioides")) %>% 
  select(Sequencing_ID,EuroMedAcceptedName_formatted, Species, Collineen, Montagnard, Alpin, Nival )

check <- subsp %>% 

  group_by(Species) %>% 

  summarise(
            c = max(Collineen),
            m = max(Montagnard),
            a = max(Alpin),
            n = max(Nival))

  subsp %>% 
    group_by(Species) %>% 
    summarise_at(vars(
                 Collineen,
                 Montagnard,
                 Alpin,
                 Nival), max)
ch <-  left_join(subsp, check, by = "Species")
ch <- ch %>% 
  select(-c(Collineen,
            Montagnard,
            Alpin,
            Nival)) 

distinct(ch, Species, .keep_all = T)  

#Now to apply it to the full dataset

check <- all_dat %>%  
  group_by(Species) %>% 
  summarise_at(vars(
    Collineen:SLO,
    Nival), max)
che <- left_join(all_dat, check, by = "Species")
che <- distinct(che, Species, .keep_all = T)
