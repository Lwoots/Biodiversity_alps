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

synonyms_euro <- read.csv(here("Data/PHYLOALPS", "SynonymyEuromed.csv"))

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

flora_alpina$FA_Full_names <- str_c(flora_alpina$Species, flora_alpina$Sous_espece, sep = " ")
flora_alpina$FA_Full_names <- str_trim(flora_alpina$FA_Full_names, side = "right") #Gets rid of white space at end of text

synonyms_euro$EuroMedAccepted_formatted <- str_replace_all(synonyms_euro$EuroMedAccepted,
                                                           pattern = c(" var. " = " ",
                                                                       " subsp. " = " "))

synonyms_euro$EuroMedNonAccepted_formatted <- str_replace_all(synonyms_euro$EuroMedNonAccepted,
                                                              pattern = c(" var. " = " ",
                                                                          " subsp. " = " "))

euro_med$EuroMedAcceptedName_formatted <- str_replace_all(euro_med$EuroMedAcceptedName,
                                                          pattern = c(" var. " = " ",
                                                                      " subsp. " = " "))

#

#Obj. 1) Update PhloAlps ####

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

pa_fa_combined <- left_join(combined_em_pa_alps, flora_alpina, by = c("EuroMedAcceptedName_formatted" = "FA_Full_names"))

pa_without_fa <- pa_fa_combined %>% 
  filter(is.na(Espece)) #1198 entries without flora alpina data


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


pa_without_updated_fa <- pa_without_updated_fa %>% 
  mutate(PhyloAlps_species_name = str_extract(EuroMedAcceptedName_formatted, "[A-Za-z]+ [a-z\\-]+") )

ch <- intersect(pa_without_updated_fa$PhyloAlps_species_name, flora_alpina$Species) #349 entries 

recoverable_names <- pa_without_updated_fa %>% 
  filter(PhyloAlps_species_name %in% ch) #All names that don't have FA but could
recoverable_names <- recoverable_names[,c(1:11,84)] #Drop useless FA data

recoverable_names_with_fa <- left_join(recoverable_names, 
                                       flora_alpina, 
                                       by = (c("PhyloAlps_species_name" = "Species")))

names(recoverable_names_with_fa)[12] <- "Species"

setdiff(names(recoverable_names_with_fa), names(pa_fa_updated_combined))
recoverable_names_with_fa <- recoverable_names_with_fa %>%
  select(!Complete_FA_Updated_by_euro) #Now column names are the same for rejoining

pa_fa_updated_combined_less_recover <- pa_fa_updated_combined %>% #The full dataset, removing rows to be rejoined later
  filter(!Sequencing_ID %in% recoverable_names_with_fa$Sequencing_ID)
  

pa_fa_updated_combined2 <- rbind(pa_fa_updated_combined_less_recover, recoverable_names_with_fa)

#What is still missing?
t <- pa_fa_updated_combined2 %>% filter(is.na(Espece)) 
#t <- t %>% 
  filter(!is.na(EuroMedAcceptedName_formatted))
t <- t %>% 
  distinct(EuroMedAcceptedName_formatted, .keep_all = T) #286 no alpina data

#Check whether, of these, any of them match between the provider names and fa

int_sub <- intersect(t$Provider_Name, flora_alpina$FA_Full_names) #64

#Do the same routine of removing them from the dataset then readding 

t <- pa_fa_updated_combined2 %>% filter(is.na(Espece)) 
int_sub <- intersect(t$Provider_Name, flora_alpina$FA_Full_names)

recoverable_names <- pa_fa_updated_combined2 %>% 
  filter(Provider_Name %in% int_sub) #All names that don't have FA but could
setdiff(int_sub, recoverable_names$Provider_Name)
recoverable_names <- recoverable_names[,c(1:11,83)] #Drop useless FA data

recoverable_names_with_fa <- left_join(recoverable_names, 
                                       flora_alpina, 
                                       by = (c("Provider_Name" = "FA_Full_names")))
#Make sure all datasets have the same columns
setdiff(names(recoverable_names_with_fa), names(pa_fa_updated_combined2))
setdiff(names(pa_fa_updated_combined2), names(recoverable_names_with_fa))

recoverable_names_with_fa <- recoverable_names_with_fa %>% 
  select(-c("FA_Updated_by_euro.x",  "FA_Updated_by_euro.y", "Complete_FA_Updated_by_euro" ))

pa_fa_updated_combined2 <- pa_fa_updated_combined2 %>% 
  select(-c("FA_Full_names", "FA_Updated_by_euro"))

#Remove newly matched names from original

pa_fa_updated_combined_less_recover <- pa_fa_updated_combined2 %>% #The full dataset, removing rows to be rejoined later
  filter(!Sequencing_ID %in% recoverable_names_with_fa$Sequencing_ID)

pa_fa_updated_combined3 <- rbind(pa_fa_updated_combined_less_recover, recoverable_names_with_fa)

#What is still missing?
d <- pa_fa_updated_combined3 %>% filter(is.na(Espece)) 

d <- d %>% 
  distinct(EuroMedAcceptedName_formatted, .keep_all = T) #173 no alpina data



#Look at species without a name from Julien
phyloalps_no_euroname <- pa_fa_updated_combined3 %>% 
  filter(is.na(EuroMedAcceptedName_formatted) & is.na(Espece))

miss <- rbind(d[,1:11], phyloalps_no_euroname[1:11])

#write.csv(miss, 
#          file = "220215Species_without_fa_data.csv", row.names = F)



taxonset <- pa_fa_updated_combined3 %>% 
  filter(!is.na(Espece))

#Pull out species from the missing list that have been manually checked
corrected_species <- read.csv(here("Data", "220228_manually_coded_species2.csv"), sep = ";")
corrected_species <- corrected_species[,1:13]

added_species <- pa_fa_updated_combined3 %>% 
  filter(Sequencing_ID %in% corrected_species$Sequencing_ID)

#d_corrected <- d %>% 
#  filter(Sequencing_ID %in% c("BGN_IBH", "BGN_LFB", "BGN_DII", "BGN_DTT", "BGN_KIE",
#                              "BGN_KHP", "BGN_CDD", "BGN_GPM", "BGN_FEM", "BGN_KGA",
#                              "BGN_APT", "BGN_IGI", "BGN_GQR", "BGN_PQF", "GWM_1539",
#                              "GWM_1231"))
#more <- phyloalps_no_euroname %>% 
#  filter(Sequencing_ID %in% c("BGN_IBH", "BGN_LFB", "BGN_DII", "BGN_DTT", "BGN_KIE",
#                              "BGN_KHP", "BGN_CDD", "BGN_GPM", "BGN_FEM", "BGN_KGA",
#                              "BGN_APT", "BGN_IGI", "BGN_GQR", "BGN_PQF", "GWM_1539",
#                              "GWM_1231"))
#added_species <- rbind(d_corrected, more)

taxonset <- rbind(taxonset, added_species)

#Remove libraries that are know to be wrong

taxonset <- taxonset %>% 
  filter(!Sequencing_ID %in% incorrect_names$Sequencing_ID)


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

set.seed(133)
taxonset <- taxonset[sample(1:nrow(taxonset)),]

final_taxonset <- distinct(taxonset, Species_name, .keep_all = T)

temp <- final_taxonset %>% 
  select(!Collineen:SLO)

final_taxonset <- left_join(temp, check, by = "Species_name")

#Get rid of not relevant columns
final_taxonset <- final_taxonset %>% 
  select(!c(DB_ID, Sample_ID, Provider_Name, Scientific_Name, Nom_Taxon_verifie_ThePlantList,
            No_Flora_alpina, No_Famille, No_Genre, No_Espece, No_Sous_espece))

#write.csv(final_taxonset, file = "220302_taxonset_resolved.csv", row.names = F)
rm(list=ls(pattern="recover"))
rm(list=ls(pattern="pa"))
rm(list=ls(pattern="che"))
rm(list=ls(pattern="euro"))
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

final_taxonset <- read.csv(here("Data", "220302_taxonset_resolved.csv"))
tree <- read.tree("/Users/larawootton/Documents/Doctorate/Tree_dating/Tree_files/211201_treePL_smooth0001.tre")
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

#write.tree(reduced_tree, file = "220302_for_Daisie_dated.tre")
#write.csv(taxon_reduced_to_tree, file = "220302_taxonset_with_tips.csv", row.names = F)
#The end



