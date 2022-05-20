#14 Feb 2022
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
library(ape)

#Data

phyloalps <- read.csv(here("Data", "PHYLOALPS_HERBARIUM_66_16feb2022.csv"), sep = ";")
phyloalps <- phyloalps %>% 
  select(DB_ID, Sample_ID, Project_name, Provider_Name, Scientific_Name, Nom_Taxon_verifie_ThePlantList,
         Taxon_family, NCBI_ID, code_CBNA, Sequencing_ID) #Picking out columns that correspond to old dataset

euro_med <- read.csv(here("Data/UPDATE2022", "PhyloAlpsDirectCorrespondances.csv"))

flora_alpina <- read.csv(here("Data", "Cristina_Etages_Provinces.csv"))

synonyms_euro <- read.csv(here("Data/UPDATE2022", "SynonymyEuromed.csv"))

incorrect_names <- read.csv(here("Data/list.error.filters.csv")) #Libraries that are incorrect

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
flora_alpina[4278, 7] <- "tubiformis"
flora_alpina[4278, 8] <- "Fritillaria tubiformis"

euro_med[2642, 19] <- "Hierochloe australis" 
euro_med[2643, 19] <- "Hierochloe odorata"
euro_med[1831, 19] <- "Hippophae rhamnoides"

#Subsp wrangling

#Next get subspecies and vars into same format across datasets

#Flora alpina
flora_alpina$FA_Full_names <- str_c(flora_alpina$Species, flora_alpina$Sous_espece, sep = " ") #~Combines species and subsp into trinomial where relevant
flora_alpina$FA_Full_names <- str_trim(flora_alpina$FA_Full_names, side = "right") #Gets rid of white space at end of text

#Direct correspondences file
euro_med$EuroMedAcceptedName_formatted <- str_replace_all(euro_med$EuroMedAcceptedName,
                                                          pattern = c(" var. " = " ",
                                                                      " subsp. " = " "))

#euromed synonyms
synonyms_euro$EuroMedAccepted_formatted <- str_replace_all(synonyms_euro$EuroMedAccepted,
                                                           pattern = c(" var. " = " ",
                                                                       " subsp. " = " "))

synonyms_euro$EuroMedNonAccepted_formatted <- str_replace_all(synonyms_euro$EuroMedNonAccepted,
                                                              pattern = c(" var. " = " ",
                                                                          " subsp. " = " "))
#Correct a few names that still contain the author's name
synonyms_euro[synonyms_euro == "Carex filiformis L."] <- "Carex filiformis" 
synonyms_euro[synonyms_euro == "Veratrum album L."] <- "Veratrum album"
synonyms_euro[synonyms_euro == "Ranunculus auricomus coll."] <- "Ranunculus auricomus"

#Replace NAs in Nonaccepted names for use in loop later

synonyms_euro$EuroMedNonAccepted_formatted <- ifelse(is.na(synonyms_euro$EuroMedNonAccepted_formatted), 
                                                     synonyms_euro$PlantListName,
                                                     synonyms_euro$EuroMedNonAccepted_formatted)
synonyms_euro[synonyms_euro == "Carex atrata subsp. atrata"] <- "Carex atrata atrata"

#Now all names are in trinomial format without var. subsp. across all datasets

#

#Obj. 1) Update PhloAlps ####
#Connect phyloalps data with Seb's names to corrected names provided by Julien

#Check whether sequencing code is being read in the same way

setdiff(euro_med$Sequencing_ID, phyloalps$Sequencing_ID)
setdiff(phyloalps$Sequencing_ID,euro_med$Sequencing_ID) #Seems good


#Reduce euromed to accepted name and sequncing ID to prevent column.x etc when joining
euro_med_reduced <- euro_med %>% 
  select("EuroMedAcceptedName_formatted", "Sequencing_ID")

#Combine with phyloalps
combined_em_pa <- left_join(euro_med_reduced, phyloalps, by = "Sequencing_ID")

#EuroMedAcceptedName_formatted is the accepted name according to Julien for each phyloAlps sample where possible

#Reduce phyloalps down to remove outgroups, norway etc

combined_em_pa_alps <- combined_em_pa %>% 
  filter(Project_name == "PhyloAlps")

#Join to flora_alpina

#pa_fa_combined <- left_join(combined_em_pa_alps, flora_alpina, by = c("EuroMedAcceptedName_formatted" = "FA_Full_names"))

#pa_without_fa <- pa_fa_combined %>% 
#  filter(is.na(Espece)) #1198 entries without flora alpina data


#Make loop that compares each flora alpina species to all the euromed synonyms, plus pulls out correct name
flora_alpina$FA_Updated_by_euro <- rep(NA, length(flora_alpina$FA_Full_names))

for (i in 1:length(flora_alpina$FA_Full_names)) {
  for (j in 1:length(synonyms_euro$EuroMedNonAccepted)) {
    if (flora_alpina$FA_Full_names[i] == synonyms_euro$EuroMedNonAccepted_formatted[j]) {
      flora_alpina$FA_Updated_by_euro[i] <- synonyms_euro$EuroMedAccepted_formatted[j]
      print(i/length(flora_alpina$FA_Full_names)*100)
    }
  }
}

#Create column updating FA with the accepted names
flora_alpina <- flora_alpina %>% 
  mutate(Complete_FA_Updated_by_euro = case_when(
    is.na(FA_Updated_by_euro) ~ FA_Full_names,
    !is.na(FA_Updated_by_euro) ~ FA_Updated_by_euro
  ))

#Now try to join with updated names

pa_fa_updated_combined <- left_join(combined_em_pa_alps, flora_alpina, by = c("EuroMedAcceptedName_formatted" = "Complete_FA_Updated_by_euro"))

pa_without_updated_fa <- pa_fa_updated_combined %>% 
  filter(is.na(Espece)) #now only 796 entries with no flora alpina data

#Work out how many of them are because of subsp issues

#Extract binomials from trinomials
pa_without_updated_fa <- pa_without_updated_fa %>% 
  mutate(PhyloAlps_species_name = str_extract(EuroMedAcceptedName_formatted, "[A-Za-z]+ [a-z\\-]+") )

flora_alpina <-  flora_alpina %>% 
  mutate(Accepted_binomials = str_extract(Complete_FA_Updated_by_euro, "[A-Za-z]+ [a-z\\-]+"))

ch <- intersect(pa_without_updated_fa$PhyloAlps_species_name, flora_alpina$Accepted_binomials) #261 entries 
#ch2 <- intersect(pa_without_updated_fa$PhyloAlps_species_name, flora_alpina$Species) #349


recoverable_names <- pa_without_updated_fa %>% 
  filter(PhyloAlps_species_name %in% ch) #All names that don't have FA because of subsp issues but could
recoverable_names <- recoverable_names[,c(1:11,84)] #Drop FA columns that are artefact of earlier joins

#Attach to fa data
recoverable_names_with_fa <- left_join(recoverable_names, 
                                       flora_alpina, 
                                       by = (c("PhyloAlps_species_name" = "Accepted_binomials")))

#Drop rows repeated due to joining (bse subspecies)
recoverable_names_with_fa <- recoverable_names_with_fa %>% 
  distinct(Sequencing_ID, .keep_all = T)


setdiff(names(recoverable_names_with_fa), names(pa_fa_updated_combined))
recoverable_names_with_fa <- recoverable_names_with_fa %>%
  select(!c(Complete_FA_Updated_by_euro, PhyloAlps_species_name)) #Now column names are the same for rejoining

#First remove recoverable names from the original dataset, to be rejoined later
pa_fa_updated_combined_less_recover <- pa_fa_updated_combined %>% 
  filter(!Sequencing_ID %in% recoverable_names_with_fa$Sequencing_ID)
  

pa_fa_updated_combined2 <- rbind(pa_fa_updated_combined_less_recover, recoverable_names_with_fa)

#What is still missing?
t <- pa_fa_updated_combined2 %>% filter(is.na(Espece)) 
#t <- t %>% 
#  distinct(EuroMedAcceptedName_formatted, .keep_all = T) #286 no alpina data

#Check whether, of these, any of them match between the provider names and fa. This checks for names
#that were updated in the direct correspondences file but fa wasn't updated with the loop

int_sub <- intersect(t$Provider_Name, flora_alpina$FA_Full_names) #64
int_bi <- intersect(t$Provider_Name, flora_alpina$Accepted_binomials)
int_sp <- intersect(t$Provider_Name, flora_alpina$Species)

all_match <- union(int_bi, int_sp)
all_match <- union(all_match, int_sub)

#Do the same routine of removing them from the dataset then readding 

recoverable_names <- pa_fa_updated_combined2 %>% 
  filter(Provider_Name %in% all_match) #All names that don't have FA but could
setdiff(all_match, recoverable_names$Provider_Name)
recoverable_names <- recoverable_names[,c(1:11)] #Drop useless FA data

recoverable_names_with_fa <- left_join(recoverable_names, 
                                       flora_alpina, 
                                       by = (c("Provider_Name" = "FA_Full_names")))
#Make sure all datasets have the same columns
setdiff(names(recoverable_names_with_fa), names(pa_fa_updated_combined2))
setdiff(names(pa_fa_updated_combined2), names(recoverable_names_with_fa))

recoverable_names_with_fa <- recoverable_names_with_fa %>% 
  select(-c("Complete_FA_Updated_by_euro", "Accepted_binomials"))


#Remove newly matched names from original

pa_fa_updated_combined_less_recover <- pa_fa_updated_combined2 %>% #The full dataset, removing rows to be rejoined later
  filter(!Sequencing_ID %in% recoverable_names_with_fa$Sequencing_ID)

pa_fa_updated_combined_less_recover <- pa_fa_updated_combined_less_recover %>% 
  select(!FA_Full_names)

pa_fa_updated_combined3 <- rbind(pa_fa_updated_combined_less_recover, recoverable_names_with_fa)

#What is still missing?
d <- pa_fa_updated_combined3 %>% filter(is.na(Espece)) 

d <- d %>% 
  distinct(EuroMedAcceptedName_formatted, .keep_all = T) #173 no alpina data


taxonset <- pa_fa_updated_combined3 %>% 
  filter(!is.na(Espece))

#Pull out species from the missing list that have been manually checked
corrected_species <- read.csv(here("Data", "220228_manually_coded_species2.csv"), sep = ";")
corrected_species <- corrected_species[,1:13]

setdiff(corrected_species$EuroMedAcceptedName_formatted, d$EuroMedAcceptedName_formatted) #Slight discrepancy bse manually checked on a older iteration of wranglng

added_species <- pa_fa_updated_combined3 %>% 
  filter(Sequencing_ID %in% corrected_species$Sequencing_ID) %>% 
  filter(!Sequencing_ID %in% c("BGN_IHB", "BGN_KCF", "RSZ_RSZAXPI000688-80")) #Ignore species that since matched

taxonset <- rbind(taxonset, added_species)

#Remove libraries that are know to be wrong

taxonset <- taxonset %>% 
  filter(!Sequencing_ID %in% incorrect_names$Sequencing_ID)

#Get rid of invasives ####

data <- read.csv(here("Data", "TRAITS_ALL_2022-03-22.csv"), sep = ";")

invasives <- data %>% 
  select(libelle, IMMIG_INVASIV) %>% 
  filter(str_detect(IMMIG_INVASIV, "neophyte")) %>% 
  mutate(Species = str_extract(libelle, "[A-Za-z]+ [a-z\\-]+"))

aliens <- data.frame(Species = intersect(taxonset$EuroMedAcceptedName_formatted, invasives$Species))


taxonset <- taxonset %>% 
  filter(!EuroMedAcceptedName_formatted %in% aliens$Species) %>% 
  filter(!Provider_Name %in% aliens$Species) %>% 
  filter(!EuroMedAcceptedName_formatted %in% c("Cylindropuntia imbricata"))

#Combine subspecies ####

#Make species column based on accepted names

taxonset <- taxonset %>% 
  mutate(Species_name_nas = str_extract(EuroMedAcceptedName_formatted, "[A-Za-z]+ [a-z\\-]+") ) %>% 
  mutate(Species_name = case_when(!is.na(Species_name_nas) ~ Species_name_nas,
                                  is.na(Species_name_nas) ~ Provider_Name)) #For species that don't have euromed accepted names



check <- taxonset %>%  
  group_by(Species_name) %>% 
  summarise_at(vars(
    Collineen:SLO), max) #Summarise ecological and location data

#Shuffle dataset to prevent biases

set.seed(1333)
taxonset <- taxonset[sample(1:nrow(taxonset)),]

final_taxonset <- distinct(taxonset, Species_name, .keep_all = T)

temp <- final_taxonset %>% 
  select(!Collineen:SLO)

final_taxonset <- left_join(temp, check, by = "Species_name")

#Get rid of not relevant columns
final_taxonset <- final_taxonset %>% 
  select(!c(DB_ID, Sample_ID, Provider_Name, Scientific_Name, Nom_Taxon_verifie_ThePlantList,
            No_Flora_alpina, No_Famille, No_Genre, No_Espece, No_Sous_espece))

##Add in a few species that occur in the alps but were sequenced under other project names
#Plus moehringiodes as a missed synonym

extras <- combined_em_pa %>% 
  filter(EuroMedAcceptedName_formatted %in% c("Ranunculus crenatus",
                                              "Arenaria moehringioides",
                                              "Carex bigelowii dacica",
                                              "Carex fuliginosa",
                                              "Poa nemoralis"),
         !Sequencing_ID == "BGN_CBK")

extras[extras$EuroMedAcceptedName_formatted == "Arenaria moehringioides", 1] <- "Arenaria ciliata"
extras[extras$EuroMedAcceptedName_formatted == "Carex bigelowii dacica", 1] <- "Carex bigelowii"

extras <- left_join(extras, flora_alpina, by = (c("EuroMedAcceptedName_formatted" = "Species")))
extras <- extras %>% 
  mutate(Species = rep(NA, 5),
         Species_name_nas = rep(NA, 5),
         Species_name = c(
                          "Arenaria moehringioides",
                          "Carex bigelowii dacica",
                          "Carex fuliginosa",
                          "Poa nemoralis",
                          "Ranunculus crenatus"))
extras <- extras %>% select(names(final_taxonset))

final_taxonset <- bind_rows(final_taxonset, extras)


write.csv(final_taxonset, file = "220520_taxonset_resolved.csv", row.names = F)
rm(list=ls(pattern="recover"))
rm(list=ls(pattern="pa"))
rm(list=ls(pattern="che"))
rm(list=ls(pattern="euro"))
rm(list=ls(pattern="int"))
rm(temp)
rm(t)
rm(d)
rm(taxonset)
rm(miss)
remove(more)
rm(incorrect_names)
rm(added_species)
rm(corrected_species)

# Join to tree ####

final_taxonset <- read.csv(here("Data", "220419_taxonset_resolved.csv"))
tree <- read.tree("/Users/larawootton/Documents/Doctorate/Tree_dating/Tree_files/220301_treePL_optimisation_run2_final.tre")
tips <- tree$tip.label

seq_ids <- str_match(tips, ".*_(.*?_.*)")[,2] #Pull out seq ids from tree tips

tips <- data.frame(Tip = tips, SeqIDS = seq_ids)

final_taxonset <- final_taxonset %>% 
  mutate(SeqIDS = str_replace_all(final_taxonset$Sequencing_ID, ":", "_")) #Format seq id names

setdiff(final_taxonset$Sequencing_ID, final_taxonset$SeqIDS)
setdiff(final_taxonset$SeqIDS, final_taxonset$Sequencing_ID)
setdiff(final_taxonset$SeqIDS, seq_ids)
setdiff(seq_ids,final_taxonset$SeqIDS)

final_taxonset <- left_join(final_taxonset, tips, by = "SeqIDS") #Join tips to full dataset
#There are about 200 species that are in the data set but not in the tree - ask Seb

taxon_reduced_to_tree <- final_taxonset %>% 
  filter(!is.na(Tip)) #Retain only data connected to a tip

reduced_tree <- keep.tip(tree, taxon_reduced_to_tree$Tip)
length(reduced_tree$tip.label)

#write.tree(reduced_tree, file = "220419_for_Daisie_datedrun2.tre")
#write.csv(taxon_reduced_to_tree, file = "220419_taxonset_with_tips.csv", row.names = F)
#The end



