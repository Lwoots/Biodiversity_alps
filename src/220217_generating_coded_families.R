#15 Feb 2022
#Script to get clades set up for DAISIE

###Set up ####

library(tidyverse)
library(here)
library(ape)
library(phytools)
library(MonoPhy)
library(ggtree)

#tree <- read.tree(here("Data", "220328_for_Daisie_datedrun2.tre"))
tree <- read.tree("/Users/larawootton/Documents/Doctorate/Tree_dating/Tree_files/220301_treePL_optimisation_run2_final.tre")
final_taxonset <- read.csv(here("Data", "220419_taxonset_with_tips.csv"))
species_manual <- read.csv(here("Data", "220228_manually_coded_species2.csv"), sep = ";")
flora_provences <- read.csv(here("Data", "220304_Provinces_sorted_to_regions.csv"), sep = ";")
incorrect_names <- read.csv(here("Data/list.error.filters.csv"))
classification <- read.csv("/Users/larawootton/Documents/Doctorate/Tree_dating/Data/210512_tip_classificationV1.csv")
phyloalps <- read.csv(here("Data", "PHYLOALPS_HERBARIUM_66_16feb2022.csv"), sep = ";")

#Data wrangling ####

final_taxonset[final_taxonset == "?"] <- 0 #Question marks in the dataset and book
setdiff(flora_provences$Department, names(final_taxonset)[c(19:74)])
setdiff(names(final_taxonset)[c(19:74)], flora_provences$Department)

flora_provences$Department[flora_provences$Department == "0_OB"] <- "O_OB"
flora_provences$Department[flora_provences$Department == "0_SW"] <- "O_SW"

names(final_taxonset)[41] <- "O_OB"
names(final_taxonset)[42] <- "O_SW"

#Drop incorrect tips from tree

tips <- tree$tip.label


seq_ids <- str_match(tips, ".*_(.*?_.*)")[,2] #Pull out seq ids from tree tips
taxon_nos <- str_match(tips, "_[0-9]+_")

tips <- data.frame(Tip = tips, SeqIDS = seq_ids, Taxon = taxon_nos)

incorrect_names <- incorrect_names %>% 
  mutate(SeqIDS = str_replace_all(incorrect_names$Sequencing_ID, ":", "_"))

wrong_tips <- tips %>% 
  filter(SeqIDS %in% incorrect_names$SeqIDS)

tips_single <- tips %>% 
  distinct(Taxon, .keep_all = T)

alps <- final_taxonset %>% 
  select(Tip, Sequencing_ID)
names(alps)[2] <- "SeqIDS"

combined_tips <- rbind(tips_single[,1:2], alps)
combined_tips <- combined_tips %>% 
  distinct(SeqIDS, .keep_all = T)

taxon_nos <- str_match(combined_tips$Tip, "_[0-9]+_")
combined_tips <- cbind(combined_tips, taxon_nos)

dups <- combined_tips %>% 
  group_by(taxon_nos) %>% 
  filter(n() > 1) %>% 
  filter(!SeqIDS %in% alps$SeqIDS)

combined_tips <- combined_tips %>% 
  filter(!SeqIDS %in% dups$SeqIDS) %>% 
  select(!taxon_nos)

tree <- keep.tip(tree, combined_tips$Tip)
tree <- drop.tip(tree, wrong_tips$Tip)

#Connect tips to families


seq_ids <- str_match(classification$Tips, ".*_(.*?_.*)")[,2]
classification$Sequencing_ID <- seq_ids

fam <- classification %>% 
  select(Sequencing_ID, family) 

tip_fams <- tree$tip.label
seq_ids <- str_match(tip_fams, ".*_(.*?_.*)")[,2]
tip_fams <- data.frame(Tip = tip_fams, Sequencing_ID = seq_ids)

tip_fams <- left_join(tip_fams, fam, by = "Sequencing_ID")

phylo_fam <- phyloalps %>% 
  select(Taxon_family, Sequencing_ID)

tip_fams <- left_join(tip_fams, phylo_fam, by = "Sequencing_ID")

tip_fams <- tip_fams %>% 
  mutate(Family = case_when(!is.na(family) ~ family,
                            is.na(family) ~ Taxon_family))

miss <- tip_fams %>% 
  filter(is.na(Family))

tree <- drop.tip(tree, miss$Tip)

tip_fams <- tip_fams %>% 
  filter(!is.na(Family))

#Code tips as alpine non-alpine

mainland <-  species_manual %>% filter(Status %in% c("Mainland")) #These are the manually added taxa
non_endem <- species_manual %>% filter(Status %in% c("Non_endemic_MaxAge"))
alpine <- species_manual %>% filter(Status %in% c("Endemic"))

final_taxonset <- final_taxonset %>% 
  mutate(Origins = case_when(Alpin %in% c(1,2) | Nival %in% c(1,2) ~ "High_elevation",
                     Alpin %in% c(0) & Nival %in% c(0) ~ "Low_elevation"
                     ),
         Status = case_when(
                    (Alpin %in% c(1,2) | Nival %in% c(1,2)) & Subalpin %in% c(0, 1) & Montagnard == 0 ~ "Endemic",
                    Alpin %in% c(1,2) & ( Subalpin > 1 | Montagnard > 0 | Collineen > 0)  ~ "Non_endemic_MaxAge",
                    Alpin == 0 & ( Subalpin > 0 | Montagnard > 0 | Collineen > 0)  ~ "Mainland",
                    Sequencing_ID %in% mainland$Sequencing_ID ~ "Mainland",
                    Sequencing_ID %in% c(non_endem$Sequencing_ID, "BGN_ECD") ~ "Non_endemic_MaxAge",
                    Sequencing_ID %in% alpine$Sequencing_ID ~ "Endemic")
  ) 

#Join to full tree tips

tree_info <- data.frame(Tip = tree$tip.label)



alpine_tips <- final_taxonset %>% 
  select(Tip, Status, Sequencing_ID, EuroMedAcceptedName_formatted, Taxon_family, Genre, Collineen, Montagnard, Subalpin, Alpin, Nival)

final_taxonset <- left_join(tree_info, alpine_tips, by = "Tip")
final_taxonset <- final_taxonset %>% 
  mutate(Status = case_when(!is.na(Status) ~ Status,
                            is.na(Status) ~ "Mainland"))

final_taxonset <- left_join(final_taxonset, tip_fams, by = "Tip")
names(final_taxonset)[15] <- "Taxon_family"
names(final_taxonset)[12] <- "Sequencing_ID"

#How many orders are there 

length(unique(final_taxonset$Famille)) 
setdiff(unique(final_taxonset$Taxon_family), unique(final_taxonset$Famille))

#tips that are misplaced

#It seems that many of them just have been given the wrong family name, rather than being phylogenetically incorrect

final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_DVC"] <- "Asteraceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_EEE"] <- "Asteraceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_CIR"] <- "Asparagaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_CAV"] <- "Fabaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_CBQ"] <- "Apiaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_CHD"] <- "Fabaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_CAN"] <- "Poaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_DRI"] <- "Caryophyllaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_DVK"] <- "Asteraceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_CKA"] <- "Papaveraceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_CIM"] <- "Brassicales"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_DRR"] <- "Polygonaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_VQ"] <- "Caprifoliaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_CKK"] <- "Apiaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_DQS"] <- "Poaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_VR"] <- "Amaryllidaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_DCF"] <- "Ranunculaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_CFT"] <- "Poaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_CRB"] <- "Poaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_DPT"] <- "Poaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_VM"] <- "Oleaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_CIM"] <- "Brassicaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_ETM"] <- "Brassicaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "RSZ_RSZAXPI000732-23"] <- "Brassicaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_ETM"] <- "Brassicaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_CRH"] <- "Alismataceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_EDQ"] <- "Ranunculaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_CHI"] <- "Ranunculaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_CIQ"] <- "Fabaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_DLL"] <- "Fabaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_DVL"] <- "Poaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_CPV"] <- "Poaceae"

final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_CGC"] <- "Asteraceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_DNC"] <- "Salicaceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_DAL"] <- "Alismataceae"
final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_VV"] <- "Lamiaceae"

final_taxonset$Taxon_family[final_taxonset$Sequencing_ID == "BGN_NFD"] <- "Primulaceae"

final_taxonset$Taxon_family[final_taxonset$Taxon_family == "Dipsacaceae"] <- "Caprifoliaceae"
final_taxonset$Taxon_family[final_taxonset$Taxon_family == "Hyacinthaceae"] <- "Asparagaceae"
final_taxonset$Taxon_family[final_taxonset$Taxon_family == "Chenopodiaceae"] <- "Amaranthaceae"
final_taxonset$Taxon_family[final_taxonset$Taxon_family == "Thesiaceae"] <- "Santalaceae"
final_taxonset$Taxon_family[final_taxonset$Taxon_family == "Valerianaceae"] <- "Caprifoliaceae"
#Wrong in tree: Linnaea_borealis BGN_KIS Hedera helix RSZ_RSZAXPI001130-87 Polemonium_caeruleum RSZ_RSZAXPI000710-107
#Osyris_alba, luckily all mainland, so can drop them


tree <- drop.tip(tree, c("Osyris_alba_350585_PHA006326_BGN_PXV", "Polemonium_caeruleum_174663_PHA006879_RSZ_RSZAXPI000710-107",
                         "Linnaea_borealis_77623_PHA005374_BGN_KIS", "Hedera_helix_4052_PHA004286_RSZ_RSZAXPI001130-87",
                         "Polycarpon_tetraphyllum_115622_PHA006882_BGN_BAM"))
final_taxonset <- final_taxonset %>% filter(!Tip %in% c("Osyris_alba_350585_PHA006326_BGN_PXV", "Polemonium_caeruleum_174663_PHA006879_RSZ_RSZAXPI000710-107",
                                                        "Linnaea_borealis_77623_PHA005374_BGN_KIS",
                                                        "Hedera_helix_4052_PHA004286_RSZ_RSZAXPI001130-87",
                                                        "Polycarpon_tetraphyllum_115622_PHA006882_BGN_BAM"))

tree <- drop.tip(tree, c("Kirkia_sp._43702_PHA010604_BGN_PWV",
                         "Festuca_breistrofferi_1532822_PHA003646_BGN_CFT",
                         "Osyris_alba_350585_PHA006325_RSZ_RSZAXPI000676-52",
                         "Adonis_flammea_1532246_PHA000148_BGN_CHI", "Papaver_dubium_subsp_lecoqii_1533249_PHA006418_BGN_CKA",
                         "Valeriana_repens_1534711_PHA009526_BGN_IIB",
                         "Salsola_kali_tragus_1671102_PHA007888_BGN_HCM",
                         "Streptopus_amplexifolius_134860_PHA008978_BGN_BRP"
                         ))
final_taxonset <- final_taxonset %>% filter(!Tip %in% c("Kirkia_sp._43702_PHA010604_BGN_PWV",
                                                        "Festuca_breistrofferi_1532822_PHA003646_BGN_CFT",
                                                        "Osyris_alba_350585_PHA006325_RSZ_RSZAXPI000676-52",
                                                        "Adonis_flammea_1532246_PHA000148_BGN_CHI", "Papaver_dubium_subsp_lecoqii_1533249_PHA006418_BGN_CKA",
                                                        "Valeriana_repens_1534711_PHA009526_BGN_IIB",
                                                        "Salsola_kali_tragus_1671102_PHA007888_BGN_HCM",
                                                        "Streptopus_amplexifolius_134860_PHA008978_BGN_BRP"))


#Check monophyly ####

taxonomy <- data.frame(Tip = final_taxonset$Tip, Family = final_taxonset$Taxon_family)

monophyly <- AssessMonophyly(tree, taxonomy, verbosity = 25)
GetSummaryMonophyly(monophyly)
results <- GetResultMonophyly(monophyly)
results$Family

write.csv(results, file = "220421_Checking_monophyly_trimmed_tree.csv", row.names = T)


#Getting data into daisie format ####

tree_order <- data.frame(Tip=tree$tip.label)

ordered_by_tree <- left_join(tree_order, final_taxonset, by = "Tip") #Organising so is in the same order as the tree

ordered_by_tree$Tip[1:20] 

#Code non-endemics

#First name clade of non-endemics

ordered_by_tree <- ordered_by_tree %>% 
  mutate(Clade_name =
           case_when(
             Status == "Non_endemic_MaxAge" ~ paste(Genre, row.names(.), sep = "_")))
ordered_by_tree <- ordered_by_tree %>% 
  select(Tip, EuroMedAcceptedName_formatted, Taxon_family, Collineen, Montagnard, Subalpin, Alpin, Nival,
         #Origins, 
         Status, Clade_name)

ordered_by_tree %>% filter(is.na(Status))

ordered_by_tree$Status[ordered_by_tree$Tip == "Coincya_cheiranthos_montana_400013_PHA002477_BGN_ECD"] <- "Non_endemic_MaxAge"

write.csv(ordered_by_tree, file = "220421_coding_endemics.csv", row.names = F)


#Make plots of all the families ####
#

node <-7183
clade <- extract.clade(tree, node)
#clade <- extract.clade(clade, 655) #655
plot(clade)

clade_data <- final_taxonset %>% 
  filter(Tip %in% clade$tip.label) %>% 
  select(Tip, Status)


stat <- as.vector(clade_data$Status)
names(stat) <- clade_data$Tip
stat

mtrees<-make.simmap(clade,
                    stat,
                    model="ER",
                    nsim=500)
pd<-summary(mtrees,plot=FALSE)

cols <- setNames(c("#D72000","#1BB6AF","#FFAD0A"),
                 c("Endemic", "Mainland", "Non_endemic_MaxAge"))
plot(pd,
     fsize=0.6,
     ftype="i", 
     colors = cols, 
     xlim = c(0,230), mar=c(3.1,0.1,0.1,0.1), cex = 0.5)
axisPhylo(side = 1)
#abline(v = 42)
add.simmap.legend(colors=cols, fsize = 0.5, x = 1, y = 0.1)


#ggplot for clarity

g <- ggtree(clade, ladderize = F) %<+% clade_data


p <- g + geom_tiplab(hjust = 0, nudge_x = 1) +
  geom_tippoint(aes(color = Status), size = 3) +
  scale_colour_manual(values = c("#D72000","#1BB6AF","#FFAD0A")) +
  theme(legend.position = "n") +
  xlim(0,90)
p


#Assign species to regions ####

west <- flora_provences %>% 
  filter(Geology_regions == "Western Alps") %>% 
  select(Department) %>% 
  pull()
east <- flora_provences %>% 
  filter(Geology_regions == "Eastern Alps") %>% 
  select(Department) %>% 
  pull()
central <- flora_provences %>% 
  filter(Geology_regions == "Central Alps") %>% 
  select(Department) %>% 
  pull()

periph <- flora_provences %>% 
  filter(Edge_regions == "Peripheral") %>% 
  select(Department) %>% 
  pull()
inter <- flora_provences %>% 
  filter(Edge_regions == "Internal") %>% 
  select(Department) %>% 
  pull()

#Sum columns in each region, then code as 0,1
final_taxonset_regions <- final_taxonset %>% 
  mutate(Sum_W = rowSums(across(all_of(west))),
         Sum_C = rowSums(across(all_of(central))),
         Sum_E = rowSums(across(all_of(east))),
         Sum_P = rowSums(across(all_of(periph))),
         Sum_I = rowSums(across(all_of(inter))),
         across(Sum_W:Sum_I, ~ if_else(.x >0, 1, 0))
         )

