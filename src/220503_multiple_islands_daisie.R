#Data wrangling for assigning species to regions and running archipelago daisie
#03 May 2022

#Packages ####

library(tidyverse)
library(here)
library(DAISIE)
library(ape)
library(castor)
library(dispRity)
library(phangorn)

#Data ####

tree <- read.tree("/Users/larawootton/Documents/Doctorate/Tree_dating/Tree_files/220301_treePL_optimisation_run2_final.tre")
final_taxonset <- read.csv(here("Data", "220419_taxonset_with_tips.csv"))
sample_classification <- read.csv(here("Data", "220421_coding_endemics_completed.csv"), sep = ";")
flora_provences <- read.csv(here("Data", "220304_Provinces_sorted_to_regions.csv"), sep = ";")


#Data wrangling ####

flora_provences[flora_provences$Department == "0_OB", 1] <- "X0_OB"
flora_provences[flora_provences$Department == "0_SW", 1] <- "X0_SW"
names(flora_provences)[3] <- "Regions2"

setdiff(flora_provences$Department, names(final_taxonset))

sample_classification$Clade_name[is.na(sample_classification$Clade_name)] <- "main"

#Subset data by region

west <- flora_provences %>% 
  filter(Geology_regions == "Western Alps") %>% 
  select(Department)
east <- flora_provences %>% 
  filter(Geology_regions == "Eastern Alps") %>% 
  select(Department)
central <- flora_provences %>% 
  filter(Geology_regions == "Central Alps") %>% 
  select(Department)

final_taxonset <- final_taxonset %>% 
  mutate(Western = case_when((CH_FR |
                             CH_V0 |
                             CH_VS |
                              F_04 |
                              F_05 |
                              F_06 |
                              F_26 |
                              F_38 |
                              F_73 |
                              F_74 |
                              F_83 |
                              F_84 |
                              I_AO |
                              I_CN |
                              I_IM |
                              I_SV |
                              I_TO |
                              I_VC) == 1 ~ 1),
         Eastern = case_when((A_B  |
                              A_K |
                              A_N |
                              A_NT |
                              A_O |
                              A_OT |
                              A_S |
                              A_ST |
                             X0_OB |
                              I_BL |
                              I_BZ |
                              I_PN |
                              I_TN |
                              I_TV |
                              I_U0 |
                              I_VI |
                               SLO) == 1 ~ 1 ),
         Central = case_when((A_V  |
                              CH_AP |
                              CH_BE |
                              CH_GL |
                              CH_GR |
                              CH_LU |
                              CH_SG |
                              CH_SZ |
                              CH_TI |
                              CH_UR |
                              CH_UW |
                              X0_SW |
                                FL |
                              I_BG |
                              I_BS |
                              I_CO |
                              I_NO |
                              I_SO |
                              I_VA |
                              I_VR) == 1 ~ 1)
  )


western_species <- final_taxonset %>% 
  filter(Western %in% c(1))

eastern_species <- final_taxonset %>% 
  filter(Eastern %in% c(1))

central_species <- final_taxonset %>% 
  filter(Central %in% c(1))

#

class_w <- sample_classification %>% 
  filter(Tip %in% western_species$Tip)

setdiff(western_species$Tip, class_w$Tip)

class_e <- sample_classification %>% 
  filter(Tip %in% eastern_species$Tip)

class_c <- sample_classification %>% 
  filter(Tip %in% central_species$Tip)

#Get them into daisie format

node_ages <- tree.age(tree)

#Western Alps ####

multi <- class_w %>% 
  group_by(Clade_name) %>% 
  filter(n() > 1) %>% #To be a radiation need more than one tip
  select(Clade_name) %>% 
  pull()

uni_multi <- unique(multi)
uni_multi <- uni_multi[-1]

#Get branching times

west_daisie_nodes <- vector()

for (i in 1:length(uni_multi)) {
  
  radiation_node <- getMRCA(tree, class_w$Tip[class_w$Clade == uni_multi[i]])
  
  stem <- Ancestors(tree, radiation_node, "parent")
  rad_times <- Descendants(tree, radiation_node, type = "all")
  rad_tips <- unlist(Descendants(tree, radiation_node, type = "tips"))
  rad_nodes <- rad_times[!rad_times %in% rad_tips]
  actual_ages <- node_ages[c(stem, radiation_node, rad_nodes), 1]
  
  west_daisie_nodes[i] <- glue::glue_collapse(actual_ages, sep = ",")
  
}

test <- data.frame(Clade_name = uni_multi, 
                   Status = rep("Endemic", length(uni_multi)),
                   Missing_species = rep(0, length(uni_multi)),
                   Branching_times = west_daisie_nodes)

#Singletons
tip_ages <- setNames(tree$edge.length[sapply(1:length(tree$tip.label),function(x,y) which (y==x),y=tree$edge[,2])],tree$tip.label)

tip_branch_ages <- data.frame(Tip = names(tip_ages),
                              Branching_times = tip_ages)

tip_branch_ages <- left_join(class_w, tip_branch_ages, by = "Tip")

tip_branch_ages <- tip_branch_ages %>% 
  filter(!Clade_name %in% c("main", uni_multi))

test2 <- data.frame(Clade_name = tip_branch_ages$Clade_name,
                    Status = tip_branch_ages$Status,
                    Missing_species = rep(0, length(tip_branch_ages$Tip)),
                    Branching_times = tip_branch_ages$Branching_times)

west_daisie_df <- rbind(test, test2)

west_daisie_df$Clade_name <- as.factor(west_daisie_df$Clade_name)
west_daisie_df$Status <- as.factor(west_daisie_df$Status)
west_daisie_df$Branching_times <- as.factor(west_daisie_df$Branching_times)


#Eastern Alps ####

multi <- class_e %>% 
  group_by(Clade_name) %>% 
  filter(n() > 1) %>% #To be a radiation need more than one tip
  select(Clade_name) %>% 
  pull()

uni_multi <- unique(multi)
uni_multi <- uni_multi[-1]

#Get branching times

east_daisie_nodes <- vector()

for (i in 1:length(uni_multi)) {
  
  radiation_node <- getMRCA(tree, class_e$Tip[class_e$Clade == uni_multi[i]])
  
  stem <- Ancestors(tree, radiation_node, "parent")
  rad_times <- Descendants(tree, radiation_node, type = "all")
  rad_tips <- unlist(Descendants(tree, radiation_node, type = "tips"))
  rad_nodes <- rad_times[!rad_times %in% rad_tips]
  actual_ages <- node_ages[c(stem, radiation_node, rad_nodes), 1]
  
  east_daisie_nodes[i] <- glue::glue_collapse(actual_ages, sep = ",")
  
}

test <- data.frame(Clade_name = uni_multi, 
                   Status = rep("Endemic", length(uni_multi)),
                   Missing_species = rep(0, length(uni_multi)),
                   Branching_times = east_daisie_nodes)

#Singletons
tip_ages <- setNames(tree$edge.length[sapply(1:length(tree$tip.label),function(x,y) which (y==x),y=tree$edge[,2])],tree$tip.label)

tip_branch_ages <- data.frame(Tip = names(tip_ages),
                              Branching_times = tip_ages)

tip_branch_ages <- left_join(class_e, tip_branch_ages, by = "Tip")

tip_branch_ages <- tip_branch_ages %>% 
  filter(!Clade_name %in% c("main", uni_multi))

test2 <- data.frame(Clade_name = tip_branch_ages$Clade_name,
                    Status = tip_branch_ages$Status,
                    Missing_species = rep(0, length(tip_branch_ages$Tip)),
                    Branching_times = tip_branch_ages$Branching_times)

east_daisie_df <- rbind(test, test2)

east_daisie_df$Clade_name <- as.factor(east_daisie_df$Clade_name)
east_daisie_df$Status <- as.factor(east_daisie_df$Status)
east_daisie_df$Branching_times <- as.factor(east_daisie_df$Branching_times)

# Central Alps ####

multi <- class_c %>% 
  group_by(Clade_name) %>% 
  filter(n() > 1) %>% #To be a radiation need more than one tip
  select(Clade_name) %>% 
  pull()

uni_multi <- unique(multi)
uni_multi <- uni_multi[-1]

#Get branching times

central_daisie_nodes <- vector()

for (i in 1:length(uni_multi)) {
  
  radiation_node <- getMRCA(tree, class_c$Tip[class_c$Clade == uni_multi[i]])
  
  stem <- Ancestors(tree, radiation_node, "parent")
  rad_times <- Descendants(tree, radiation_node, type = "all")
  rad_tips <- unlist(Descendants(tree, radiation_node, type = "tips"))
  rad_nodes <- rad_times[!rad_times %in% rad_tips]
  actual_ages <- node_ages[c(stem, radiation_node, rad_nodes), 1]
  
  central_daisie_nodes[i] <- glue::glue_collapse(actual_ages, sep = ",")
  
}

test <- data.frame(Clade_name = uni_multi, 
                   Status = rep("Endemic", length(uni_multi)),
                   Missing_species = rep(0, length(uni_multi)),
                   Branching_times = central_daisie_nodes)

#Singletons
tip_ages <- setNames(tree$edge.length[sapply(1:length(tree$tip.label),function(x,y) which (y==x),y=tree$edge[,2])],tree$tip.label)

tip_branch_ages <- data.frame(Tip = names(tip_ages),
                              Branching_times = tip_ages)

tip_branch_ages <- left_join(class_c, tip_branch_ages, by = "Tip")

tip_branch_ages <- tip_branch_ages %>% 
  filter(!Clade_name %in% c("main", uni_multi))

test2 <- data.frame(Clade_name = tip_branch_ages$Clade_name,
                    Status = tip_branch_ages$Status,
                    Missing_species = rep(0, length(tip_branch_ages$Tip)),
                    Branching_times = tip_branch_ages$Branching_times)

central_daisie_df <- rbind(test, test2)

central_daisie_df$Clade_name <- as.factor(central_daisie_df$Clade_name)
central_daisie_df$Status <- as.factor(central_daisie_df$Status)
central_daisie_df$Branching_times <- as.factor(central_daisie_df$Branching_times)

write.csv(west_daisie_df, file = "220504_western_daise_df.csv", row.names = F)
write.csv(east_daisie_df, file = "220504_eastern_daise_df.csv", row.names = F)
write.csv(central_daisie_df, file = "220504_central_daise_df.csv", row.names = F)

#All together

DAISIE_plot_island(central_daisie_df, island_age = 35)

central_datalist <- DAISIE_dataprep( 
  datatable = central_daisie_df, 
  island_age = 35, 
  2300) 
west_datalist <- DAISIE_dataprep( 
  datatable = west_daisie_df, 
  island_age = 35, 
  2300) 

east_datalist <- DAISIE_dataprep( 
  datatable = east_daisie_df, 
  island_age = 35, 
  2300) 


archipelago <- list(west_datalist, east_datalist, central_datalist)

names(archipelago) <- c("WesternAlps", "EasternAlps", "CentralAlps")

DAISIE_ML(
  datalist = archipelago,
  datatype = "multiple",
  initparsopt = c(0.14,0.07,9100,0.006,0.5,0.14,0.14),
  idparsmat = rbind(c(1,2,3,4,5),c(6,2,3,4,5),c(7,2,3,4,5)),
  idparsopt = c(1:7),
  parsfix = NULL,
  idparsfix = NULL
)

rbind(1:5,c(6,2,3,7,5),1:5,1:5)
DAISIE_ML( 
  datalist = daisie_datalist, 
  initparsopt = c(0.14,0.07,9100,0.006,0.5), 
  idparsmat = rbind(1:5,c(6,2,3,7,5),1:5),
  idparsopt = c(2,4,5,6,7),
  parsfix = c(1.1,Inf),
  idparsfix = c(1,3)
) 



utils::data(Galapagos_datalist)
DAISIE_ML(
  datalist = list(west_datalist,east_datalist),
  datatype = 'multiple',
  initparsopt = c(2.5,2.7,20,0.009,1.01,0.009,1.01),
  idparsmat = rbind(1:5,c(1,2,3,6,7)),
  idparsopt = 1:7,
  parsfix = NULL,
  idparsfix = NULL
)


#LTT
phytools::ltt(tree)

w <- class_w %>% filter(!Status == "Mainland")
west_tree <- keep.tip(tree, w$Tip)
phytools::ltt(tree, xlim = c(310,350), log.lineages = F)



west_col <- west_daisie_df %>% 
  filter(Status == "Non_endemic_MaxAge") %>% 
  mutate(Region = rep("West", length(west_col$Status)))
east_col <- east_daisie_df %>% 
  filter(Status == "Non_endemic_MaxAge")%>% 
  mutate(Region = rep("East", length(east_col$Status)))
central_col <- central_daisie_df %>% 
  filter(Status == "Non_endemic_MaxAge")%>% 
  mutate(Region = rep("Central", length(central_col$Status)))

all_dat <- rbind(west_col, east_col, central_col)
all_dat$Branching_times <- as.numeric(as.character(all_dat$Branching_times))

all <- ggplot(all_dat, aes(x = Branching_times, colour = Region)) +
  geom_density() +
  xlim(0, 30) +
  xlab("Time of colonisation")


disti <- final_taxonset %>% 
  mutate(Distribution = case_when(
    Central + Western + Eastern == 3 ~ "WCE",
    Eastern == 1 & is.na(Central) & is.na(Western) ~ "E",
    is.na(Eastern) & is.na(Central) & Western == 1 ~ "W",
    is.na(Eastern) & Central == 1 & is.na(Western) ~ "C",
    Eastern == 1 & Central == 1 & is.na(Western) ~ "EC",
    is.na(Eastern) & Central == 1 & Western == 1 ~ "WC",
    Eastern == 1 & is.na(Central) & Western == 1 ~ "EW"
    
  )) %>% 
  select(Tip, Western, Eastern, Central, Distribution) 

disti %>% 
  group_by(Distribution) %>% 
  count()

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

ana <- tip_branch_ages %>% 
  filter(!Status %in% c("Non_endemic_MaxAge"))


tip_branch_ages <- rbind(tip_branch_ages, ana)


all_col_dist <- left_join(tip_branch_ages, disti, by = "Tip")
all_col_dist <- all_col_dist %>% 
  filter(!is.na(Distribution))

all_col_dist %>% 
  group_by(Distribution) %>% 
  count()

all_reg <- ggplot(all_col_dist, aes(x = Branching_times, colour = Distribution)) +
  geom_density() +
  xlim(0, 30) +
  scale_color_brewer(palette = "Set1") +
  xlab("Time of colonisation")
  
w <- 
ggplot(all_col_dist) +
  geom_density(aes(x = Branching_times, colour = "Western"))+
  geom_density(aes(x = Branching_times, colour = "red"))

ggpubr::ggarrange(all, all_reg, ncol = 1)
