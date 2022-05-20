#Workflow to get DAISIE format 16/03/2022

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


#Wrangling####

#Reduce tree to relevant tips

length(tree$tip.label)
tree <- keep.tip(tree, sample_classification$Tip)

#Add filler to clade_name NAs

sample_classification$Clade_name[is.na(sample_classification$Clade_name)] <- "main"

#Correct some species that were filtered out incorrectly

sample_classification$Status[sample_classification$Tip == "Ranunculus_breyninus_286868_PHA007444_BGN_HAR"] <- "Non_endemic_MaxAge"
sample_classification$Clade_name[sample_classification$Tip == "Ranunculus_breyninus_286868_PHA007444_BGN_HAR"] <- "Ranunbrey"

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

## Extract all branch lengths from endemic clades

#Get list of endemic radiation nodes

multi <- sample_classification %>% 
  group_by(Clade_name) %>% 
  filter(n() > 1) %>% #To be a radiation need more than one tip
  select(Clade_name) %>% 
  pull()

uni_multi <- unique(multi)
uni_multi <- uni_multi[-1]


#Test on a single radiation
radiation_node <- getMRCA(tree, sample_classification$Tip[sample_classification$Clade_name =="AndrosaceA"])

radiation <- extract.clade(tree, radiation_node)
plot(radiation)

#Get branching times

stem <- Ancestors(tree, radiation_node, "parent")
rad_times <- Descendants(tree, radiation_node, type = "all")
rad_tips <- unlist(Descendants(tree, radiation_node, type = "tips"))
rad_nodes <- rad_times[!rad_times %in% rad_tips]
all_daisie_nodes <- as.vector(c(stem, radiation_node, rad_nodes))
node_ages[all_daisie_nodes,]


#now make it a loop

all_daisie_nodes <- vector()

for (i in 1:length(uni_multi)) {
  
  radiation_node <- getMRCA(tree, sample_classification$Tip[sample_classification$Clade == uni_multi[i] ])
  
  stem <- Ancestors(tree, radiation_node, "parent")
  rad_times <- Descendants(tree, radiation_node, type = "all")
  rad_tips <- unlist(Descendants(tree, radiation_node, type = "tips"))
  rad_nodes <- rad_times[!rad_times %in% rad_tips]
  actual_ages <- node_ages[c(stem, radiation_node, rad_nodes), 1]
  
  all_daisie_nodes[i] <- glue::glue_collapse(actual_ages, sep = ",")
  
}


test <- data.frame(Clade_name = uni_multi, 
                   Status = rep("Endemic", length(uni_multi)),
                   Missing_species = rep(0, length(uni_multi)),
                   Branching_times = all_daisie_nodes)

#Singletons

tip_ages <- setNames(tree$edge.length[sapply(1:length(tree$tip.label),function(x,y) which (y==x),y=tree$edge[,2])],tree$tip.label)

tip_branch_ages <- data.frame(Tip = names(tip_ages),
                              Branching_times = tip_ages)

tip_branch_ages <- left_join(sample_classification, tip_branch_ages, by = "Tip")

tip_branch_ages <- tip_branch_ages %>% 
  filter(!Clade_name %in% c("main", uni_multi))

test2 <- data.frame(Clade_name = tip_branch_ages$Clade_name,
                    Status = tip_branch_ages$Status,
                    Missing_species = rep(0, length(tip_branch_ages$Tip)),
                    Branching_times = tip_branch_ages$Branching_times)

daisie_df <- rbind(test, test2)

daisie_df$Clade_name <- as.factor(daisie_df$Clade_name)
daisie_df$Status <- as.factor(daisie_df$Status)
daisie_df$Branching_times <- as.factor(daisie_df$Branching_times)

#Add in missing species to relevent clades

daisie_df$Status[daisie_df$Tip == "Euphrasia_minima_475014_PHA003596_BGN_ADL"] <- "Endemic"
daisie_df$Missing_species[daisie_df$Clade_name == "Euphrasia_1749"] <- 1 #Euphrasia inopinata

daisie_df$Status[daisie_df$Tip == "Juncus_triglumis_223682_PHA004931_BGN_BLE"] <- "Endemic"
daisie_df$Missing_species[daisie_df$Clade_name == "Juncus_4995"] <- 1 #Juncus biglumis

daisie_df$Status[daisie_df$Tip == "Lolium_perenne_4522_PHA005443_BGN_CGE"] <- "Endemic"
daisie_df$Missing_species[daisie_df$Clade_name == "Lolium_4585"] <- 1 #Juncus biglumis

daisie_df$Missing_species[daisie_df$Clade_name == "Ranunculus_4035"] <- 1 #Ranunculus venetus
daisie_df$Missing_species[daisie_df$Clade_name == "Androsace_3359"] <- 1 #Androsace helvetica
daisie_df$Missing_species[daisie_df$Clade_name == "Gentiana_2090"] <- 1 #Gentiana pilosa
daisie_df$Missing_species[daisie_df$Clade_name == "SaxifragaA"] <- 1 #Saxifraga aspera
daisie_df$Missing_species[daisie_df$Clade_name == "Poa_4668"] <- 1 #Poa pumila
daisie_df$Missing_species[daisie_df$Clade_name == "Gymnadenia_4202"] <- 1 #Gymnadenia buschmanniae

#write.csv(daisie_df, file = "220428_DAISE_df.csv", row.names = F)

#Choosing reasonable parameter values ####

#Alp age = 35 Myr

continent <- sample_classification %>% 
  filter(Status == "Mainland") %>% 
  filter(!is.na(Nival))

#Cladogenesis
sp <- sample_classification %>% 
  group_by(Clade_name) %>% 
  filter(n() > 1) %>% 
  filter(!Status %in% "Mainland") #552 species by cladogenesis, 116 lineages (length(uni_multi))
552/116 # 4.758621 events per lineage
log(4.758621)/35 #0.04457023 events per lineage per Myr

clado <- c(0.001, 0.045, 1.01)

#Anagenesis

ana <- tip_branch_ages %>% 
  filter(Status %in% c("Endemic"))
length(ana$Status)/35 #1 event per Myr


#Per lineage, but what lineages? All island? All mainland? Is a lineage a species or event?

35/length(tip_branch_ages$Status == "Non_endemic_MaxAge") #Number anagenetic divided by island lineages 0.07526882
0.1/35 #No. anagenetic events per lineage per year = 0.002

anagen <- c(0.0001, 0.002, 1.01)

#Extinction

0.5*0.04
0.9*0.04

extinct <- c(0.0001, 0.02, 0.036)

#Carrying capacity

total_sp <- sample_classification %>% 
  filter(!Status %in% c("Mainland"))

cc <- c(910, 1500)

#Colonisation

colon <- sample_classification %>% 
  filter(Status %in% c("Non_endemic_MaxAge"))

length(uni_multi)+length(colon$Status) + length(ana$Status)

470/2243
0.2095408/35 #0.006

colon <- c(0.001, 0.2, 1.01)
#


#Running Daisie####



DAISIE_plot_island(daisie_df, island_age = 35)
DAISIE_plot_age_diversity(daisie_df, island_age = 35)

daisie_datalist <- DAISIE_dataprep( 
  datatable = daisie_df, 
  island_age = 35, 
  M = length(sample_classification$Status[sample_classification$Status == "Mainland"]) 
)


DAISIE_ML( 
  datalist = daisie_datalist, 
  initparsopt = c(0.14,0.07,9100,0.006,0.5), 
  ddmodel = 0, 
  idparsopt = c(1,3,2,4,5), 
  parsfix = NULL, 
  idparsfix = NULL
) 


#Maximum likelihood parameter estimates: lambda_c: 0.667480, mu: 0.919318, K: Inf,
#gamma: 0.119875, lambda_a: 0.055941 
#Maximum loglikelihood: -3460.398762 
#lambda_c        mu   K     gamma   lambda_a    loglik df conv
#0.6674796 0.9193183 Inf 0.1198749 0.05594082 -3460.399  5    0

#Maximum likelihood parameter estimates: lambda_c: 0.662345, mu: 0.807555, K: 0.607670,
#gamma: 0.082731, lambda_a: 0.000001 
#Maximum loglikelihood: -3731.004979 
#lambda_c        mu         K      gamma     lambda_a    loglik df conv
#1 0.6623449 0.8075551 0.6076704 0.08273111 1.422221e-06 -3731.005  5    0

DAISIE_tut


#Make matrix of all combinations for running on cluster

expand.grid(clado, cc,colon, anagen)



#Simulate

clado_rate <- 1.127254 # cladogenesis rate
ext_rate <- 1.486281 # extinction rate
clade_carr_cap <- 9100  # clade-level carrying capacity
imm_rate <- 0.2338993 # immigration rate
ana_rate <- 0.07476353 # anagenesis rate

island_replicates_K <- DAISIE_sim_constant_rate( 
  time = 35, 
  M = 2300, 
  pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
  replicates = 10,
  plot_sims = T,
  divdepmodel = 'IW',
  verbose = TRUE
) 
DAISIE_plot_sims(island_replicates_K10)


#Summary stats ####


#Colonisation times

actual_ages <- vector()

for (i in 1:length(uni_multi)) {
  
  radiation_node <- getMRCA(tree, sample_classification$Tip[sample_classification$Clade == uni_multi[i] ])
  
  stem <- Ancestors(tree, radiation_node, "parent")
  actual_ages[i] <- node_ages[c(stem), 1]
  
}

hist(actual_ages, xlim = c(0,35), breaks = 100)

all_colon <- c(actual_ages, tip_branch_ages$Branching_times)
hist(all_colon, xlim = c(0,35), breaks = 200)

cdf <- data.frame(Colonisation = all_colon)

ggplot(cdf, aes(x = Colonisation)) +
  geom_histogram(colour = "white", binwidth = 0.5) +
  xlim(c(0,35)) +
  ylim(c(0,150)) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 18)) +
  xlab("Time of colonisation (MYA)") +
  ylab("Frequency")

#Proportion of lineages that colonised

tip_branch_ages %>% 
  filter(Status == "Endemic") #62 anagenetic 0.096875 of lineages
62/904 #0.068 of species

532 - 62  #470 colonisations 0.73 of lineages
470/904 #0.519

108 #radiations 0.169
904-532
372/904 #0.41

ldf <- data.frame(type = c("Colonisation", "Cladogenesis", "Anagenesis"),
                  value = c(73,16.9, 9.68))
sdf <- data.frame(type = c("Colonisation", "Cladogenesis", "Anagenesis"),
                  value = c(51.9, 41, 6.8))

library(khroma)

p<-ggplot(ldf, aes(x=type, y=value, fill=type)) +
  geom_bar(stat="identity",colour = "black", width = 0.9) +
  theme_classic() +
  ylab("Percentage of lineages (%)") +
  xlab("") +
  ylim(c(0,100)) +
  scale_fill_muted() +
  theme(axis.text = element_text(colour = "black", size = 10),
        legend.position = "none")
p


d<-ggplot(sdf, aes(x=type, y=value, fill=type)) +
  geom_bar(stat="identity", colour = "black", width = 0.9) +
  theme_classic() +
  ylab("Percentage of species (%)") +
  xlab("") +
  ylim(c(0,100)) +
  scale_fill_muted() +
  theme(axis.text = element_text(colour = "black", size = 10),
        legend.position = "none") 
d



#Speciation timing





rad_ages <- list()

for (i in 1:length(uni_multi)) {
  
  radiation_node <- getMRCA(tree, sample_classification$Tip[sample_classification$Clade == uni_multi[i] ])
  
  rad_times <- Descendants(tree, radiation_node, type = "all")
  rad_tips <- unlist(Descendants(tree, radiation_node, type = "tips"))
  rad_nodes <- rad_times[!rad_times %in% rad_tips]
  rad_ages[[i]] <- node_ages[c(radiation_node, rad_nodes), 1]
  
  
}

all_sp <- unlist(rad_ages)
hist(all_sp)

csdf <- data.frame(Sp = all_sp)

ggplot(csdf, aes(x = Sp)) +
  geom_histogram(colour = "white", binwidth = 1) +
  xlim(c(0,35)) +
  ylim(c(0,150)) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 18)) +
  xlab("Time of cladogenetic speciation (MYA)") +
  ylab("Frequency")


#LTT

w <- sample_classification %>% filter(!Status == "Mainland")
west_tree <- keep.tip(tree, w$Tip)
phytools::ltt(west_tree, log.lineages = T)
phytools::ltt(tree, xlim = c(310,350), log.lineages = T)

plot(west_tree)
