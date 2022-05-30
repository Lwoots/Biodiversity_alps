#Plotting alpine history - normalised rates of colonisation and speciation

#Packages ####
library(tidyverse)
library(here)
library(castor)
library(ape)
library(ggtree)
library(phangorn)
library(glue)
library(dispRity)
library(DAISIE)

#Data ####

tree <- read.tree("/Users/larawootton/Documents/Doctorate/Tree_dating/Tree_files/220301_treePL_optimisation_run2_final.tre")
sample_classification <- read.csv(here("Data", "220421_coding_endemics_completed.csv"), sep = ";")
ana <- read.csv(here("Data", "220517_anagenetic_species.csv"))

#Wrangling####

#Reduce tree to relevant tips

length(tree$tip.label)
tree <- keep.tip(tree, sample_classification$Tip)

#Add filler to clade_name NAs

sample_classification$Clade_name[is.na(sample_classification$Clade_name)] <- "main"

#Correct some species that were filtered out incorrectly

sample_classification$Status[sample_classification$Tip == "Ranunculus_crenatus_286892_CAR007457_BGN_MEG"] <- "Endemic"
sample_classification$Clade_name[sample_classification$Tip == "Ranunculus_crenatus_286892_CAR007457_BGN_MEG"] <- "Ranunculus_4011"

sample_classification$Status[sample_classification$Tip == "Aquilegia_alpina_218850_PHA009929_BGN_FQ"] <- "Non_endemic_MaxAge"
sample_classification$Clade_name[sample_classification$Tip == "Aquilegia_alpina_218850_PHA009929_BGN_FQ"] <- "Aquilalp"

sample_classification$Status[sample_classification$Tip == "Arenaria_ciliata_multicaulis_1403275_PHA000797_BGN_FLQ"] <- "Endemic"
sample_classification$Clade_name[sample_classification$Tip == "Arenaria_ciliata_multicaulis_1403275_PHA000797_BGN_FLQ"] <- "Arenmulti"

sample_classification$Status[sample_classification$Tip == "Trifolium_alpinum_97008_PHA010063_BGN_KSK"] <- "Non_endemic_MaxAge"
sample_classification$Clade_name[sample_classification$Tip == "Trifolium_alpinum_97008_PHA010063_BGN_KSK"] <- "Trifalp"

sample_classification$Status[sample_classification$Tip == "Euphrasia_salisburgensis_475020_PHA003604_BGN_BIG"] <- "Endemic"
sample_classification$Status[sample_classification$Tip == "Euphrasia_micrantha_475390_TROM_V_131542_BXA_ALA"] <- "Endemic"
sample_classification$Clade_name[sample_classification$Tip == "Euphrasia_micrantha_475390_TROM_V_131542_BXA_ALA"] <- "Euphrasia_1752"

sample_classification$Status[sample_classification$Tip == "Valeriana_montana_243118_PHA010071_AXZ_T"] <- "Endemic"
sample_classification$Clade_name[sample_classification$Tip == "Valeriana_montana_243118_PHA010071_AXZ_T"] <- "Valeriana_3013"

sample_classification$Status[sample_classification$Tip == "Carex_bigelowii_dacica_1134499_CAR001779_BGN_LMT"] <- "Non_endemic_MaxAge"
sample_classification$Clade_name[sample_classification$Tip == "Carex_bigelowii_dacica_1134499_CAR001779_BGN_LMT"] <- "Carbig"

sample_classification$Status[sample_classification$Tip == "Carex_fuliginosa_241210_CAR001871_BGN_NAG"] <- "Non_endemic_MaxAge"
sample_classification$Clade_name[sample_classification$Tip == "Carex_fuliginosa_241210_CAR001871_BGN_NAG"] <- "Carful"

sample_classification$Status[sample_classification$Tip == "Poa_nemoralis_carpatica_1907909_CAR006858_BGN_NCI"] <- "Non_endemic_MaxAge"
sample_classification$Clade_name[sample_classification$Tip == "Poa_nemoralis_carpatica_1907909_CAR006858_BGN_NCI"] <- "Poanem"

sample_classification$Status[sample_classification$Tip == "Artemisia_umbelliformis_72355_PHA000880_BGN_CTR"] <- "Endemic"
sample_classification$Clade_name[sample_classification$Tip == "Artemisia_umbelliformis_72355_PHA000880_BGN_CTR"] <- "Artemisia_2292"

#Get node ages

node_ages <- tree.age(tree)

#Get list of endemic radiation nodes

multi <- sample_classification %>% 
  group_by(Clade_name) %>% 
  filter(n() > 1) %>% #To be a radiation need more than one tip
  select(Clade_name) %>% 
  pull()

uni_multi <- unique(multi)
uni_multi <- uni_multi[-1]

## Colonisation through time ####

#first get branch lengths of colonist and anagenetic tips

tip_ages <- setNames(tree$edge.length[sapply(1:length(tree$tip.label),function(x,y) which (y==x),y=tree$edge[,2])],tree$tip.label)

tip_branch_ages <- data.frame(Tip = names(tip_ages),
                              Branching_times = tip_ages)
tip_branch_ages <- left_join(sample_classification, tip_branch_ages, by = "Tip")

colonist_ages <- tip_branch_ages %>% 
  filter(Status == "Non_endemic_MaxAge")

ana_ages <- tip_branch_ages %>% 
  filter(Tip %in% ana$Tip)

colonist_ages <- bind_rows(colonist_ages, ana_ages)

#Get stem age of clades






all_stem_ages <- vector()

for (i in 1:length(uni_multi)) {
  
  radiation_node <- getMRCA(tree, sample_classification$Tip[sample_classification$Clade == uni_multi[i] ])
  
  stem <- Ancestors(tree, radiation_node, "parent")
  all_stem_ages[i] <- node_ages[c(stem), 1]

}

all_stem_ages <- data.frame(Branching_times = all_stem_ages)

all_colonist_ages <- bind_rows(colonist_ages, all_stem_ages)
all_colonist_ages[all_colonist_ages$Branching_times > 35, 11] <- 35
events <- hist(all_colonist_ages$Branching_times, breaks = seq(0,35,1), plot = FALSE)$counts
events_per_mya <- events/2300

plot(0:34, events_per_mya)

colon_df <- data.frame(Events = events_per_mya, Time = 0:-34)

ggplot(colon_df, aes(x = Time, y = Events)) +
  geom_point() +
  geom_smooth()


#Cladogenesis

island <- sample_classification %>% 
  filter(Status == "Endemic"
  )

tree_island <- keep.tip(tree, island$Tip)
tree_island <- drop.tip(tree_island, ana$Tip)


is.ultrametric(tree_island)
tree_island <- force.ultrametric(tree_island, method = "extend")

max(nodeHeights(tree_island)) - 35
length(tree_island$tip.label)

no_lineages <- count_lineages_through_time(tree_island,
                                           times = seq(124.4999, 159.4999, 1),
                                           ultrametric = T)

plot(no_lineages$times, no_lineages$lineages)


max(no_lineages$lineages)

rates <- vector()
for (i in 2:length(no_lineages$lineages)) {
  difference <- no_lineages$lineages[i] - no_lineages$lineages[i-1]
  rates[i] <- difference/no_lineages$lineages[i-1]
}


sp_df <- data.frame(Events = rates, Time = -35:0)

clad <- ggplot(sp_df, aes(x = Time, y = Events)) +
  theme_bw() +
  geom_point(color = "#D72000", alpha = 0.4) +
  geom_smooth(colour = "#D72000")+
  ylab("Speciation (normalised by cladogenetic sp)") +
  theme(axis.title.y = element_text(size = 10))

#Nomalising by total lineages

colonist_ages[colonist_ages$Branching_times > 35, 11] <- 35
events <- hist(colonist_ages$Branching_times, breaks = seq(0,35,1), plot = FALSE)$counts

colon_event <- cumsum(rev(events))

all_rates <- vector()
for (i in 2:length(no_lineages$lineages)) {
  difference <- no_lineages$lineages[i] - no_lineages$lineages[i-1]
  all_rates[i] <- difference/(no_lineages$lineages[i-1] + colon_event[i-1])
}

nor_sp_df <- data.frame(Events = all_rates, Time = -35:0)

nor_clado <- ggplot(nor_sp_df, aes(x = Time, y = Events)) +
  theme_bw() +
  geom_point(color = "#D72000", alpha = 0.4) +
  geom_smooth(colour = "#D72000")+
  ylab("Speciation (normalised by all lineages)") +
  theme(axis.title.y = element_text(size = 10))






full_island <- sample_classification %>% 
  filter(Status %in% c("Endemic", "Non_endemic_MaxAge"))
  
full_island <- force.ultrametric(full_island, method = "extend")
  
  max(nodeHeights(tree_island)) - 35
  length(tree_island$tip.label)
  
  no_lineages <- count_lineages_through_time(tree_island,
                                             times = seq(124.4999, 159.4999, 1),
                                             ultrametric = T)
  
  plot(no_lineages$times, no_lineages$lineages)
  
  
  max(no_lineages$lineages)
  
  rates <- vector()
  for (i in 2:length(no_lineages$lineages)) {
    difference <- no_lineages$lineages[i] - no_lineages$lineages[i-1]
    rates[i] <- difference/no_lineages$lineages[i-1]
  }
  
  






colo <- ggplot(colon_df, aes(x = Time, y = Events)) +
  geom_point(color = "#FFAD0A", alpha = 0.4) +
  geom_smooth(color = "#FFAD0A") +
  theme_bw() +
  ylab("Colonisation")

ggpubr::ggarrange(nor_clado, clad, colo, ncol = 1)
