#Creating daisie like data from bedrock info

#### Packages ####

library(tidyverse)
library(here)
library(ape)
library(MonoPhy)
library(ggtree)
library(paletteer)

#### Data ####

tree <- read.tree("/Users/larawootton/Documents/Doctorate/Tree_dating/Tree_files/220301_treePL_optimisation_run2_final.tre")
raw_data <- read.csv(here("Data", "220613_bedrock_data.csv"), sep = ";")
incorrect_names <- read.csv(here("Data/list.error.filters.csv"))
classification <- read.csv("/Users/larawootton/Documents/Doctorate/Tree_dating/Data/210512_tip_classificationV1.csv")
phyloalps <- read.csv(here("Data", "PHYLOALPS_HERBARIUM_66_16feb2022.csv"), sep = ";")
final_taxonset <- read.csv(here("Data", "220609_taxonset_with_tips.csv"))

#### Data wrangling ####

# Combine data by accepted name
#Temporary step until I double check data matching/synonyms

sum_data <- raw_data %>%  
  group_by(EuroMedAcceptedName_formatted) %>% 
  summarise_at(vars(
    Ca:Si), max)
temp <- raw_data %>% 
  select(!Ca:Si)

bedrock <- left_join(temp,sum_data, by = "EuroMedAcceptedName_formatted")
bedrock <- distinct(bedrock, EuroMedAcceptedName_formatted, .keep_all = T)
bedrock <- bedrock %>% 
  select(Tip, EuroMedAcceptedName_formatted, Ca:Si)

#Filter out incorrect tips in the tree
tips <- tree$tip.label

seq_ids <- str_match(tips, ".*_(.*?_.*)")[,2] #Pull out seq ids from tree tips
taxon_nos <- str_match(tips, "_[0-9]+_")

tips <- data.frame(Tip = tips, SeqIDS = seq_ids, Taxon = taxon_nos)

incorrect_names <- incorrect_names %>% 
  mutate(SeqIDS = str_replace_all(incorrect_names$Sequencing_ID, ":", "_"))

wrong_tips <- tips %>% 
  filter(SeqIDS %in% incorrect_names$SeqIDS) #Remove incorrect libraries


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
  filter(!is.na(Family)) %>% 
  select(!c(family, Taxon_family))

## Code tips by bedrock ####

combined <- left_join(tip_fams, bedrock, by = "Tip")


combined <- combined %>% 
  mutate(Rock_type = case_when(
    Ca == 2 & Si %in% c(0, 1) ~ "Calcareous",
    Si == 2 & Ca %in% c(0, 1) ~ "Siliceous",
    Si == 2 & Ca == 2 ~ "Generalist",
    Si == 1 & Ca == 1 ~ "Generalist",
    Si  %in% c(0, 1) & Ca  %in% c(0, 1) ~ "Generalist",
    is.na(Si) ~ "Generalist"
  ))

#Making dataset
#Order by tree
tree_order <- data.frame(Tip=tree$tip.label)

ordered_by_tree <- left_join(tree_order, combined, by = "Tip") #Organising so is in the same order as the tree

ordered_by_tree <- ordered_by_tree %>% 
  mutate(Clade_name =
           case_when(
             Rock_type == "Calcareous" ~ paste(EuroMedAcceptedName_formatted, "Calc", row.names(.), sep = "_"),
             Rock_type == "Siliceous" ~ paste(EuroMedAcceptedName_formatted, "Sil", row.names(.), sep = "_"),
             Si > 0 & Rock_type == "Generalist" ~ "TBA"
         ))

write.csv(ordered_by_tree, file = "220613_coding_bedrock.csv", row.names = F)

#Check monophyly ####

taxonomy <- data.frame(Tip = combined$Tip, Family = combined$Family)

monophyly <- AssessMonophyly(tree, taxonomy, verbosity = 25)
GetSummaryMonophyly(monophyly)
results <- GetResultMonophyly(monophyly)
results$Family

#write.csv(results, file = "220613_monophyly_bedrock_tree.csv", row.names = T)

#More wrangling ####

tree <- drop.tip(tree, c("Osyris_alba_350585_PHA006326_BGN_PXV", "Polemonium_caeruleum_174663_PHA006879_RSZ_RSZAXPI000710-107",
                         "Linnaea_borealis_77623_PHA005374_BGN_KIS", "Hedera_helix_4052_PHA004286_RSZ_RSZAXPI001130-87",
                         "Polycarpon_tetraphyllum_115622_PHA006882_BGN_BAM"))
combined <- combined %>% filter(!Tip %in% c("Osyris_alba_350585_PHA006326_BGN_PXV", "Polemonium_caeruleum_174663_PHA006879_RSZ_RSZAXPI000710-107",
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
combined <- combined %>% filter(!Tip %in% c("Kirkia_sp._43702_PHA010604_BGN_PWV",
                                                        "Festuca_breistrofferi_1532822_PHA003646_BGN_CFT",
                                                        "Osyris_alba_350585_PHA006325_RSZ_RSZAXPI000676-52",
                                                        "Adonis_flammea_1532246_PHA000148_BGN_CHI", "Papaver_dubium_subsp_lecoqii_1533249_PHA006418_BGN_CKA",
                                                        "Valeriana_repens_1534711_PHA009526_BGN_IIB",
                                                        "Salsola_kali_tragus_1671102_PHA007888_BGN_HCM",
                                                        "Streptopus_amplexifolius_134860_PHA008978_BGN_BRP"))

combined$Family[combined$Sequencing_ID == "BGN_DVC"] <- "Asteraceae"
combined$Family[combined$Sequencing_ID == "BGN_EEE"] <- "Asteraceae"
combined$Family[combined$Sequencing_ID == "BGN_CIR"] <- "Asparagaceae"
combined$Family[combined$Sequencing_ID == "BGN_CAV"] <- "Fabaceae"
combined$Family[combined$Sequencing_ID == "BGN_CBQ"] <- "Apiaceae"
combined$Family[combined$Sequencing_ID == "BGN_CHD"] <- "Fabaceae"
combined$Family[combined$Sequencing_ID == "BGN_CAN"] <- "Poaceae"
combined$Family[combined$Sequencing_ID == "BGN_DRI"] <- "Caryophyllaceae"
combined$Family[combined$Sequencing_ID == "BGN_DVK"] <- "Asteraceae"
combined$Family[combined$Sequencing_ID == "BGN_CKA"] <- "Papaveraceae"
combined$Family[combined$Sequencing_ID == "BGN_CIM"] <- "Brassicales"
combined$Family[combined$Sequencing_ID == "BGN_DRR"] <- "Polygonaceae"
combined$Family[combined$Sequencing_ID == "BGN_VQ"] <- "Caprifoliaceae"
combined$Family[combined$Sequencing_ID == "BGN_CKK"] <- "Apiaceae"
combined$Family[combined$Sequencing_ID == "BGN_DQS"] <- "Poaceae"
combined$Family[combined$Sequencing_ID == "BGN_VR"] <- "Amaryllidaceae"
combined$Family[combined$Sequencing_ID == "BGN_DCF"] <- "Ranunculaceae"
combined$Family[combined$Sequencing_ID == "BGN_CFT"] <- "Poaceae"
combined$Family[combined$Sequencing_ID == "BGN_CRB"] <- "Poaceae"
combined$Family[combined$Sequencing_ID == "BGN_DPT"] <- "Poaceae"
combined$Family[combined$Sequencing_ID == "BGN_VM"] <- "Oleaceae"
combined$Family[combined$Sequencing_ID == "BGN_CIM"] <- "Brassicaceae"
combined$Family[combined$Sequencing_ID == "BGN_ETM"] <- "Brassicaceae"
combined$Family[combined$Sequencing_ID == "RSZ_RSZAXPI000732-23"] <- "Brassicaceae"
combined$Family[combined$Sequencing_ID == "BGN_ETM"] <- "Brassicaceae"
combined$Family[combined$Sequencing_ID == "BGN_CRH"] <- "Alismataceae"
combined$Family[combined$Sequencing_ID == "BGN_EDQ"] <- "Ranunculaceae"
combined$Family[combined$Sequencing_ID == "BGN_CHI"] <- "Ranunculaceae"
combined$Family[combined$Sequencing_ID == "BGN_CIQ"] <- "Fabaceae"
combined$Family[combined$Sequencing_ID == "BGN_DLL"] <- "Fabaceae"
combined$Family[combined$Sequencing_ID == "BGN_DVL"] <- "Poaceae"
combined$Family[combined$Sequencing_ID == "BGN_CPV"] <- "Poaceae"
combined$Family[combined$Sequencing_ID == "BGN_CGC"] <- "Asteraceae"
combined$Family[combined$Sequencing_ID == "BGN_DNC"] <- "Salicaceae"
combined$Family[combined$Sequencing_ID == "BGN_DAL"] <- "Alismataceae"
combined$Family[combined$Sequencing_ID == "BGN_VV"] <- "Lamiaceae"
combined$Family[combined$Sequencing_ID == "BGN_NFD"] <- "Primulaceae"
combined$Family[combined$Family == "Dipsacaceae"] <- "Caprifoliaceae"
combined$Family[combined$Family == "Hyacinthaceae"] <- "Asparagaceae"
combined$Family[combined$Family == "Chenopodiaceae"] <- "Amaranthaceae"
combined$Family[combined$Family == "Thesiaceae"] <- "Santalaceae"
combined$Family[combined$Family == "Valerianaceae"] <- "Caprifoliaceae"




#Make plots of all the families ####
#

node <- 9074
clade <- extract.clade(tree, node)
#clade <- extract.clade(clade, 655) #655
plot(clade)
#nodelabels()

clade_data <- combined %>% 
  filter(Tip %in% clade$tip.label) %>% 
  select(Tip, Rock_type)
#nodelabels()

stat <- as.vector(clade_data$Rock_type)
names(stat) <- clade_data$Tip
stat

mtrees<-make.simmap(clade,
                    stat,
                    model="ER",
                    nsim=500)
pd<-summary(mtrees,plot=FALSE)

cols <- setNames(c("#548F01FF","#F6F18FFF","#CFA3EEFF"),
                 c("Calcareous", "Generalist", "Siliceous"))
plot(pd,
     fsize=0.6,
     ftype="i", 
     colors = cols, 
     xlim = c(0,75), mar=c(3.1,0.1,0.1,0.1), cex = 0.2)
axisPhylo(side = 1)
#abline(v = 42)
add.simmap.legend(colors=cols, fsize = 0.5, x = 1, y = 0.1)


#Ggtree
paletteer_d("nationalparkcolors::Voyageurs")
g <- ggtree(clade, ladderize = F) %<+% clade_data


p <- g + geom_tiplab(hjust = 0, nudge_x = 1, size = 1.5) +
  geom_tippoint(aes(color = Rock_type), size = 1) +
  scale_colour_manual(values = c("#548F01FF","#F6F18FFF","#CFA3EEFF")) +
  theme(legend.position = "n") +
  xlim(20,70)
p
