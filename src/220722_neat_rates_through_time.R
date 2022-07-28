#Plotting alpine history - normalised rates of colonisation and speciation take 2
#Neat version of 220530_rates through_time.R 

#Packages ####
library(tidyverse)
library(here)
library(castor)
library(ape)
library(phangorn)
library(glue)
library(dispRity)
library(phytools)


#Data ####

tree <- read.tree("/Users/larawootton/Documents/Doctorate/Tree_dating/Tree_files/220301_treePL_optimisation_run2_final.tre")
sample_classification <- read.csv(here("Data", "220711_coding_endemics_completed.csv"), sep = ";")
final_taxonset <- read.csv(here("Data", "220609_taxonset_with_tips.csv"))

#Wrangling####

#Get rid of blank column
sample_classification <- sample_classification %>% 
  select(-X)

#Tree wrangling ###

tree <- keep.tip(tree, sample_classification$Tip)

## Make list of anagenetic species ##

groups <- sample_classification %>% 
  filter(Status == "Endemic") %>% 
  group_by(Clade_name) %>% 
  summarise(No = n())

ana_names <- groups %>% 
  filter(No == 1)
ana <- sample_classification %>% 
  filter(Clade_name %in% ana_names$Clade_name)

#Checking no sp in each category

all_island <- sample_classification %>% 
  filter(Status %in% c("Endemic", "Non_endemic_MaxAge"))

groups <- all_island %>% 
  filter(Status == "Endemic") %>% 
  group_by(Clade_name) %>% 
  summarise(No = n())


all_island <- all_island %>% 
  mutate(Category = case_when(
    Status == "Non_endemic_MaxAge" ~ "Colonist",
    Clade_name %in% c(ana$Clade_name) ~ "Anagenesis",
    Status == "Endemic" ~"Cladogenesis"
  ))

total_sp <- length(all_island$Tip)

colonists <- all_island %>% 
  filter(Category == "Colonist")
total_col <- length(colonists$Tip)
anagen <- all_island %>% 
  filter(Category == "Anagenesis")
total_ana <- length(anagen$Tip)

clado <- all_island %>% 
  filter(Category == "Cladogenesis")
total_clado <- length(clado$Tip)

## drop unneeded subspecies which would be biasing speciation time ##

#Break tip names into binomials
sample_classification <- sample_classification %>% 
  mutate(Names = str_match(sample_classification$Tip,  "(.+?)(?=(_[0-9]))")[,1]) #Extract names

sample_classification$Names <- str_replace_all(sample_classification$Names, "_", " ") #Remove underscore

sample_classification <- sample_classification %>% 
  mutate(Binomials = str_extract(Names, "[A-Za-z]+ [a-z\\-]+") ) #Convert trinomials to binomials

#Group by binomial

species_groups <- sample_classification %>% 
  group_by(Binomials) %>% 
  summarise(no = n()) #Count number of times each name occurs

species_groups <- species_groups %>% 
  filter(no > 1) #Find which binomials have more than one occurrence

complexes <- sample_classification %>% 
  filter(Binomials %in% species_groups$Binomials)

names_tb_removed <- complexes %>% 
  filter(is.na(EuroMedAcceptedName_formatted) & Status == "Mainland") #Pull out subspecies that aren't in cleaned/FA data set

sample_classification <- sample_classification %>% 
  filter(!Tip %in% names_tb_removed$Tip) #Remove them from data

tree <- drop.tip(tree, names_tb_removed$Tip)

rm(species_groups, complexes, names_tb_removed, groups)

## Make list of radiations ##

#Add filler to clade_name NAs
sample_classification$Clade_name[is.na(sample_classification$Clade_name)] <- "main"

multi <- sample_classification %>% 
  group_by(Clade_name) %>% 
  filter(n() > 1) %>% #To be a radiation need more than one tip
  select(Clade_name) %>% 
  pull()

uni_multi <- unique(multi)
uni_multi <- uni_multi[-1] #gets rid of filler

## Analysis ####

#Get node ages
node_ages <- tree.age(tree)


## Colonisation rate ####

#Arrival of non endemics, anagens, stem age of radiations

#first get branch lengths of colonist and anagenetic tips

tip_ages <- setNames(tree$edge.length[sapply(1:length(tree$tip.label),function(x,y) which (y==x),y=tree$edge[,2])],tree$tip.label)

tip_branch_ages <- data.frame(Tip = names(tip_ages),
                              Branching_times = tip_ages)
tip_branch_ages <- left_join(sample_classification, tip_branch_ages, by = "Tip") #Join to full dataset

#Filter colonists
colonist_ages <- tip_branch_ages %>% 
  filter(Status == "Non_endemic_MaxAge")
#Filter anagenesists
ana_ages <- tip_branch_ages %>% 
  filter(Tip %in% ana$Tip)

colonist_ages <- bind_rows(colonist_ages, ana_ages)

#Then get stem ages of radiations

all_stem_ages <- vector()

for (i in 1:length(uni_multi)) {
  
  radiation_node <- getMRCA(tree, sample_classification$Tip[sample_classification$Clade == uni_multi[i]])
  stem <- Ancestors(tree, radiation_node, "parent")
  all_stem_ages[i] <- node_ages[c(stem), 1]
  
}

all_stem_ages <- data.frame(Branching_times = all_stem_ages) #Convert to df 
hist(all_stem_ages$Branching_times)

all_arrivals <- rbind(data.frame(Branching_times = colonist_ages$Branching_times), all_stem_ages) #Join

#Make all colonisation events that are older than the alpls equivalent to the age of the alps
all_arrivals[all_arrivals$Branching_times > 35, 1] <- 35
max(all_arrivals$Branching_times)
hist(all_arrivals$Branching_times, breaks = seq(0,35,1))

#Count number of arrivals in 1My time intervals
events <- hist(all_arrivals$Branching_times, breaks = seq(0,35,1), plot = FALSE)$counts

#Standardise by regional pool

events_per_mya <- events/2300

colon_df <- data.frame(Events = events_per_mya, Time = 0:-34)

#Plot

colo <- ggplot(colon_df, aes(x = Time, y = Events)) +
  geom_smooth(alpha=0.3, size=0, span=0.5) +
  stat_smooth(geom="line", color = "#FFAD0A",alpha=0.4, size=1, span=0.5) +
  geom_point(color = "#FFAD0A", alpha = 1) +
  theme_bw() +
  ylab("Colonisation rate per Myr")
colo

## Anagenesis rate ####

ana_events <- hist(ana_ages$Branching_times, breaks = seq(0,35,1), plot = FALSE)$counts
ana_events_per_mya <- ana_events/2300
ana_df <- data.frame(Events = ana_events_per_mya, Time = 0:-34)

#Plot

ana_plot <- ggplot(ana_df, aes(x = Time, y = Events)) +
  geom_smooth(alpha=0.3, size=0, span=0.5) +
  stat_smooth(geom="line", color = "dark blue",alpha=0.4, size=1, span=0.5) +
  geom_point(color = "dark blue", alpha = 1) +
  theme_bw() +
  ylab("Anagenesis rate per Myr")
ana_plot

## Cladogenetic rates ####

#In order to not count deep ghost lineages joining independent radiations, need to cut each radiation into mini tree

#Extract crown node for all clades
radiation_node <- vector()
actual_ages <- vector()

for (i in 1:length(uni_multi)) {
  radiation_node[i] <- getMRCA(tree, sample_classification$Tip[sample_classification$Clade == uni_multi[i]])
  actual_ages[i] <- node_ages[radiation_node[i], 1]
}

ggplot(data.frame(Age = actual_ages), aes(x = Age)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  xlab("Age of crown node (Myr)") +
  ylab("No. clades")

#Use crown nodes to extract clades
#But also don't want to count non-island sp within those clades, so drop them
non_island <- sample_classification %>% 
  filter(Status %in% c("Mainland","Non_endemic_MaxAge") | Tip %in% ana$Tip)

clade_list <- list()

for (i in 1:length(radiation_node)) {
  clade <- extract.clade(tree, radiation_node[i], root.edge = T)
  clade_list[[i]] <- drop.tip(clade, tip = non_island$Tip)
}

tips <- list()

for (i in 1:length(clade_list)) {
  
  tips[[i]] <- clade_list[[i]]$tip.label
}

tips <- unlist(tips)

setdiff(tips, clado$Tip)


#Get branching times

all_branches <- list()
all_sp_no <- list()

for (i in 1:length(clade_list)) {
  bt <- branching.times(clade_list[[i]])
  all_branches[[i]] <- c(300, sort(bt, decreasing = T), 0) 
  
  no_sp <- c(1:length(bt)) +1 #+one to get sp no from branch split
  all_sp_no[[i]] <- c(0, no_sp, max(no_sp))
}

#Count species in time slices

lins <- matrix(nrow = 106, ncol = 36)

for (i in 1:length(clade_list)) {
  for(j in 0:35) {
    
    raw_lineages <- data.frame(branch_times = all_branches[[i]],
                               species_no = all_sp_no[[i]])
    lins[i, j+1] <- max(raw_lineages$species_no[raw_lineages$branch_times >= j]) #How many species are there after a given time point
    
  }
}

lins <- as.data.frame(lins)

#Find cumulative number of species over time

cumsp <- colSums(lins)
cumsp <- rev(cumsp) #Reversed for plotting purposes

rev_events <- cumsum(rev(events))
#
all_rates <- vector()
for (i in 2:length(lins)) {
  difference <- cumsp[i] - cumsp[i-1]
  all_rates[i] <- difference/(cumsp[i-1] + rev_events[i-1])
}

clade_rates <- data.frame(Events = all_rates, Time = -35:0)

clado <- ggplot(clade_rates, aes(x = Time, y = Events)) +
  theme_bw() +
  geom_point(color = "#D72000")  +
  geom_smooth(alpha=0.3, size=0, span=0.5) +
  stat_smooth(geom="line", color = "#D72000",alpha=0.4, size=1, span=0.5) +
  ylab("Cladogenesis rate per Myr") +
  theme(axis.title.y = element_text(size = 10))
clado


#Combined plots ####

combo <- ggplot() +
  theme_bw() +
  geom_point(data = clade_rates, aes(x = Time, y = Events), color = "#D72000", alpha = 0.9) +
  #geom_smooth(data = clade_rates, aes(x = Time, y = Events),alpha=0.1, fill = "#D72000",size=0, span=0.5) +
  stat_smooth(data = clade_rates, aes(x = Time, y = Events),geom="line", color = "#D72000",alpha=0.4, size=1, span=0.5) +
  geom_point(data = colon_df, aes(x = Time, y = Events), color = "#FFAD0A", alpha = 0.9) +
  #geom_smooth(data = colon_df, aes(x = Time, y = Events*10),alpha=0.1, fill = "#FFAD0A", size=0, span=0.5) +
  stat_smooth(data = colon_df, aes(x = Time, y = Events),geom="line", color = "#FFAD0A",alpha=0.4, size=1, span=0.5) +
  geom_point(data = ana_df, aes(x = Time, y = Events), color = "dark blue", alpha = 0.9) +
  #geom_smooth(data = ana_df, aes(x = Time, y = Events*100),alpha=0.1, fill = "dark blue", size=0, span=0.5) +
  stat_smooth(data = ana_df, aes(x = Time, y = Events),geom="line", color = "dark blue",alpha=0.4, size=1, span=0.5) +
  ylab("All rates")


ggpubr::ggarrange(ana_plot, colo, clado, combo, ncol = 1)


#Rates by family ####


fam_tree <- tree

fam_names <- c("Asteraceae", "Poaceae", "Brassicaceae", "Caryophyllaceae", "Cyperaceae", "Rosaceae", "Primulaceae")
fam_names <- c("Asteraceae")
colon_df <- list()
clado_df <- list()


for (i in 1:length(fam_names)) {
  
  family_classification <- sample_classification %>% 
    filter(Taxon_family == fam_names[i])
  
  tree <- keep.tip(fam_tree, family_classification$Tip)
  node_ages <- tree.age(tree)
  
  family_classification$Clade_name[is.na(family_classification$Clade_name)] <- "main"
  
  multi <- family_classification %>% 
    group_by(Clade_name) %>% 
    filter(n() > 1) %>% #To be a radiation need more than one tip
    select(Clade_name) %>% 
    pull()
  
  uni_multi <- unique(multi)
  uni_multi <- uni_multi[-1]
  
  #first get branch lengths of colonist and anagenetic tips
  
  tip_ages <- setNames(tree$edge.length[sapply(1:length(tree$tip.label),function(x,y) which (y==x),y=tree$edge[,2])],tree$tip.label)
  
  tip_branch_ages <- data.frame(Tip = names(tip_ages),
                                Branching_times = tip_ages)
  tip_branch_ages <- left_join(family_classification, tip_branch_ages, by = "Tip") #Join to full dataset
  
  #Filter colonists
  colonist_ages <- tip_branch_ages %>% 
    filter(Status == "Non_endemic_MaxAge")
  #Filter anagenesists
  ana_ages <- tip_branch_ages %>% 
    filter(Tip %in% ana$Tip)
  
  colonist_ages <- bind_rows(colonist_ages, ana_ages)
  
  #Then get stem ages of radiations
  
  all_stem_ages <- vector()
  
  for (j in 1:length(uni_multi)) {
    
    radiation_node <- getMRCA(tree, family_classification$Tip[family_classification$Clade == uni_multi[j]])
    stem <- Ancestors(tree, radiation_node, "parent")
    all_stem_ages[j] <- node_ages[c(stem), 1]
    
  }
  
  all_stem_ages <- data.frame(Branching_times = all_stem_ages) #Convert to df 
  #hist(all_stem_ages$Branching_times)
  
  all_arrivals <- rbind(data.frame(Branching_times = colonist_ages$Branching_times), all_stem_ages) #Join
  
  #Make all colonisation events that are older than the alpls equivalent to the age of the alps
  all_arrivals[all_arrivals$Branching_times > 35, 1] <- 35
  
  #Count number of arrivals in 1My time intervals
  events <- hist(all_arrivals$Branching_times, breaks = seq(0,35,1), plot = FALSE)$counts
  
  #Standardise by regional pool
  
  events_per_mya <- events/2300
  
  colon_df[[i]] <- data.frame(Events = events_per_mya, Time = 0:-34, Family = rep(fam_names[i], 35))
  
  #Cladogenesis
  
  #Extract crown node for all clades
  radiation_node <- vector()
  actual_ages <- vector()
  
  for (c in 1:length(uni_multi)) {
    radiation_node[c] <- getMRCA(tree, family_classification$Tip[family_classification$Clade == uni_multi[c]])
    actual_ages[c] <- node_ages[radiation_node[c], 1]
  }
  
  #But also don't want to count non-island sp within those clades, so drop them
  non_island <- family_classification %>% 
    filter(Status %in% c("Mainland","Non_endemic_MaxAge"))
  
  clade_list <- list()
  
  for (b in 1:length(radiation_node)) {
    clade <- extract.clade(tree, radiation_node[b], root.edge = T)
    clade_list[[b]] <- drop.tip(clade, tip = non_island$Tip)
  }
  
  #Get branching times
  
  all_branches <- list()
  all_sp_no <- list()
  
  for (p in 1:length(clade_list)) {
    bt <- branching.times(clade_list[[p]])
    all_branches[[p]] <- c(300, sort(bt, decreasing = T), 0) 
    
    no_sp <- c(1:length(bt)) +1 #+one to get sp no from branch split
    all_sp_no[[p]] <- c(0, no_sp, max(no_sp))
  }
  
  lins <- matrix(nrow = length(uni_multi), ncol = 36)
  
  for (q in 1:length(clade_list)) {
    for(r in 0:35) {
      
      raw_lineages <- data.frame(branch_times = all_branches[[q]],
                                 species_no = all_sp_no[[q]])
      lins[q, r+1] <- max(raw_lineages$species_no[raw_lineages$branch_times >= r]) #How many species are there after a given time point
      
    }
  }
  
  lins <- as.data.frame(lins)
  
  #Find cumulative number of species over time
  
  cumsp <- colSums(lins)
  cumsp <- rev(cumsp) #Reversed for plotting purposes
  
  rev_events <- cumsum(rev(events))
  #
  all_rates <- vector()
  for (s in 2:length(lins)) {
    difference <- cumsp[s] - cumsp[s-1]
    all_rates[s] <- difference/(cumsp[s-1] + rev_events[s-1])
  }
  
  clado_df[[i]] <- data.frame(Events = all_rates, Time = -35:0, Family = rep(fam_names[i], 36))
  
}


full_col_df <- rbind(colon_df[[1]], colon_df[[2]], colon_df[[3]], 
                     colon_df[[4]], colon_df[[5]], colon_df[[6]], 
                     colon_df[[7]])

full_clad_df <- rbind(clado_df[[1]], clado_df[[2]], clado_df[[3]],
                      clado_df[[4]], clado_df[[5]], clado_df[[6]],
                      clado_df[[7]])

#Plot

colo <- ggplot(full_col_df, aes(x = Time, y = Events, fill = Family, color = Family)) +
  #geom_smooth(alpha=0.3, size=0, span=0.5) +
  stat_smooth(geom="line",alpha=0.4, size=1, span=0.5) +
  #geom_point(alpha = 1) +
  theme_bw() +
  ylab("Colonisation rate per Myr")
colo


clad <- ggplot(full_clad_df, aes(x = Time, y = Events, fill = Family, color = Family)) +
  #geom_smooth(alpha=0.3, size=0, span=0.5) +
  stat_smooth(geom="line",alpha=0.4, size=1, span=2) +
  #geom_point(alpha = 1) +
  theme_bw() +
  ylab("Cladogenesis rate per Myr")


clad2 <- ggplot(full_clad_df, aes(x = Time, y = Events, fill = Family, color = Family)) +
  #geom_smooth(alpha=0.3, size=0, span=0.5) +
  stat_smooth(geom="line",alpha=0.4, size=1, span=2) +
  geom_point(alpha = 1) +
  theme_bw() +
  ylab("Cladogenesis rate per Myr")
clad2

ggpubr::ggarrange(colo, clad, clad2, ncol = 1, common.legend = T)


ggplot(clado_df[[7]], aes(x = Time, y = Events, fill = Family, color = Family)) +
  #geom_smooth(alpha=0.3, size=0, span=0.5) +
  stat_smooth(geom="line",alpha=0.4, size=1, span=1) +
  geom_point(alpha = 1) +
  theme_bw() +
  ylab("Cladogenesis rate per Myr")


