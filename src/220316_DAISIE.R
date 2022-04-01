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

tree <- read.tree(here("Data", "220328_for_Daisie_datedrun2.tre"))
sample_classification <- read.csv(here("Data", "220331_coding_endemics.csv"), sep = ";")


#Wrangling####

#Reduce tree to relevant tips

length(tree$tip.label)
tree <- keep.tip(tree, sample_classification$Tip)

#Add filler to clade_name NAs

sample_classification$Clade_name[is.na(sample_classification$Clade_name)] <- "main"

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

#write.csv(daisie_df, file = "220321_DAISE_df_test.csv", row.names = F)

#Choosing reasonable parameter values ####

#Alp age = 35 Myr

#Cladogenesis
sp <- sample_classification %>% 
  group_by(Clade_name) %>% 
  filter(n() > 1) %>% 
  filter(!Status %in% "Mainland") #552 species by cladogenesis, 116 lineages (length(uni_multi))
552/116 # 4.758621 events per lineage
log(4.758621)/35 #0.1359 events per lineage per Myr

[0.001; 0.07; 0.14; 0.28; 1; 2]

#Anagenesis

ana <- tip_branch_ages %>% 
  filter(Status %in% c("Endemic"))
length(ana$Status)/35 #1 event per Myr


#Per lineage, but what lineages? All island? All mainland? Is a lineage a species or event?

35/(length(uni_multi)+length(tip_branch_ages$Status)) #Number anagenetic divided by island lineages 0.07526882
0.07526882/35 #No. anagenetic events per lineage per year = 0.002

[0.001, 0.002, 0.01, 0.1, 0.5, 1]

#Extinction

0.5*0.14
0.9*0.14

[0; 0.07; 0.13]

#Carrying capacity

total_sp <- sample_classification %>% 
  filter(!Status %in% c("Mainland"))

[910; 1500, inf]

#Colonisation

colon <- sample_classification %>% 
  filter(Status %in% c("Non_endemic_MaxAge"))

length(uni_multi)+length(colon$Status) + length(ana$Status)

470/2243
0.2095408/35 #0.006

#


#Running Daisie####



DAISIE_plot_island(daisie_df, island_age = 35)

daisie_datalist <- DAISIE_dataprep( 
  datatable = daisie_df, 
  island_age = 35, 
  M = length(sample_classification$Status[sample_classification$Status == "Mainland"]) 
)


DAISIE_ML( 
  datalist = daisie_datalist, 
  initparsopt = c(0.14,0.07,9100,0.006,0.5), 
  ddmodel = 0, 
  idparsopt = 1:5, 
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