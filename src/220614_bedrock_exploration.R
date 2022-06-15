#Initial bedrock exploration

# Packages ####
library(tidyverse)
library(here)
library(ape)
library(castor)
library(phytools)
library(dispRity)

# Data ####
data <- read.csv(here("Data", "220614_coding_bedrock.csv"), sep = ";")
tree <- read.tree("/Users/larawootton/Documents/Doctorate/Tree_dating/Tree_files/220301_treePL_optimisation_run2_final.tre")

#Wrangling ####
sample_classification <- data %>% 
  filter(Family %in% c("Poaceae", "Brassicaceae", "Fabaceae", "Caryophyllaceae", 
                       "Cyperaceae", "Apiaceae", "Ranunculaceae", "Lamiaceae",
                       "Primulaceae", "Saxifragaceae", "Campanulaceae", "Orobanchaceae",
                       "Plantaginaceae", "Caprifoliaceae", "Rubiaceae", "Gentianaceae",
                       "Juncaceae", "Violaceae", "Asparagaceae"))

#Comparison of siliceous and calcareous timings

#Get node ages

node_ages <- tree.age(tree)

#Get list of endemic radiation nodes

multi <- sample_classification %>% 
  group_by(Clade_name) %>% 
  filter(n() > 1) %>% #To be a radiation need more than one tip
  select(Clade_name) %>% 
  pull()

uni_multi <- unique(multi)
uni_multi <- uni_multi[-c(1:2)]

#Get stem ages

all_stem_ages <- vector()


for (i in 1:length(uni_multi)) {
  
  tips <- sample_classification %>% 
    filter(Clade_name == uni_multi[i]) %>% 
    pull(Tip)
  
  radiation_node <- getMRCA(tree, tips)
  
  stem <- Ancestors(tree, radiation_node, "parent")
  all_stem_ages[i] <- node_ages[c(stem), 1]
  
}

all_stem_ages <- data.frame(Branching_times = all_stem_ages,
                            Clade_name = uni_multi)
all_stem_ages <- all_stem_ages %>% 
  mutate(Type = case_when(
    str_detect(Clade_name, "Calc") ~ "Calcareous",
    str_detect(Clade_name, "Sil") ~ "Siliceous"
  ))

boxplot(all_stem_ages$Branching_times ~ all_stem_ages$Type)
t.test(all_stem_ages$Branching_times[all_stem_ages$Type == "Calcareous"], all_stem_ages$Branching_times[all_stem_ages$Type == "Siliceous"])

rad <- ggplot(all_stem_ages, aes(x = Type, y = Branching_times, fill = Type)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Radiations") +
  ylab("Stem age (Ma)") +
  scale_fill_manual(values = c("#548F01FF","#CFA3EEFF")) +
 theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = seq(0, 30, 2.5), limits = c(0,30))

#Singletons

tip_ages <- setNames(tree$edge.length[sapply(1:length(tree$tip.label),function(x,y) which (y==x),y=tree$edge[,2])],tree$tip.label)

tip_branch_ages <- data.frame(Tip = names(tip_ages),
                              Branching_times = tip_ages)
tip_branch_ages <- left_join(sample_classification, tip_branch_ages, by = "Tip")


#pull out singleton specialists

singletons <- tip_branch_ages %>% 
  group_by(Clade_name) %>% 
  filter(n() == 1) %>% 
  select(Clade_name) %>% 
  pull()

singles <- tip_branch_ages %>% 
  filter(Clade_name %in% c(singletons, "check"))

boxplot(singles$Branching_times ~ singles$Rock_type)
t.test(singles$Branching_times[singles$Rock_type == "Calcareous"], singles$Branching_times[singles$Rock_type == "Siliceous"])

sing <- ggplot(singles, aes(x = Rock_type, y = Branching_times, fill = Rock_type)) +
  geom_boxplot() +
  theme_bw() +
  ylab("Stem age (Ma)") +
  xlab("Arrivals") + 
  theme(panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#548F01FF","#F6F18FFF","#CFA3EEFF")) +
  scale_y_continuous(breaks = seq(0, 30, 2.5), limits = c(0,30))


t.test(singles$Branching_times[singles$Rock_type == "Siliceous"], all_stem_ages$Branching_times[all_stem_ages$Type == "Siliceous"])
t.test(singles$Branching_times[singles$Rock_type == "Calcareous"], all_stem_ages$Branching_times[all_stem_ages$Type == "Calcareous"])

ggpubr::ggarrange(sing, rad, ncol = 2, common.legend = T, legend = "right")

#By family ####

multi <- sample_classification %>% 
  filter(Family == "Brassicaceae",
         !Clade_name == "check", 
         !is.na(Clade_name)) %>% 
  group_by(Clade_name) %>% 
  filter(n() > 1) %>% #To be a radiation need more than one tip
  select(Clade_name) %>% 
  pull()

uni_multi <- unique(multi)


poa <- all_stem_ages %>% 
  filter(Clade_name %in% uni_multi) 
ggplot(poa, aes(x = Branching_times, y = Type, colour = Type)) +
  geom_violin() +
  geom_point() +
  scale_colour_manual(values = c("#548F01FF","#CFA3EEFF")) 


clade_dat <- sample_classification %>% 
  filter(!Clade_name == "check", 
         !is.na(Clade_name))

clade_dat <- left_join(clade_dat, all_stem_ages, by = "Clade_name")
clade_dat <- clade_dat %>% 
  filter(!is.na(Branching_times)) %>% 
  distinct(Clade_name, .keep_all = T)

ggplot(clade_dat, aes(x = Branching_times, y = Type, colour = Type)) +
  geom_boxplot(aes(x = Branching_times, y = Type, fill = Type), alpha = 0.4) +
  geom_point(size = 4) +
  scale_colour_manual(values = c("#548F01FF","#CFA3EEFF")) +
  scale_fill_manual(values = c("#548F01FF","#CFA3EEFF")) +
  facet_grid(Family ~ .) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlab("Stem clade age")


#Through time plots ####

multi <- sample_classification %>% 
  group_by(Clade_name) %>% 
  filter(n() > 1) %>% #To be a radiation need more than one tip
  select(Clade_name) %>% 
  filter(str_detect(Clade_name, "Sil")) %>% 
  pull()

uni_multi <- unique(multi)
#uni_multi <- uni_multi[-c(1:2)]

radiation_node <- vector()
actual_ages <- vector()

for (i in 1:length(uni_multi)) {
  tips <- sample_classification %>% 
    filter(Clade_name == uni_multi[i]) %>% 
    pull(Tip)
  
  radiation_node[i] <- getMRCA(tree, tips)
  actual_ages[i] <- node_ages[radiation_node[i], 1]
}

#Use crown nodes to extract clades

non_spe <- sample_classification %>% 
  filter(Clade_name %in% c(NA, "check")) %>% 
  pull(Tip)

clade_list <- list()

for (i in 1:length(radiation_node)) {
  clade <- extract.clade(tree, radiation_node[i], root.edge = T)
  clade_list[[i]] <- drop.tip(clade, non_spe)
}

#Make into species number and and branching times dfs

#Branching times

all_branches <- list()
all_sp_no <- list()

for (i in 1:length(clade_list)) {
  bt <- branching.times(clade_list[[i]])
  #stem <- clade_list[[i]]$root.edge + max(bt)
  bt <- c(#stem, 
    bt)
  all_branches[[i]] <- c(300, sort(bt, decreasing = T), 0) 
  
  no_sp <- c(1:length(bt))
  all_sp_no[[i]] <- c(0, no_sp, max(no_sp))
}


#Count species in time slices

lins <- matrix(nrow = length(uni_multi), ncol = 36)

for (i in 1:length(clade_list)) {
  for(j in 0:35) {
    
    raw_lineages <- data.frame(branch_times = all_branches[[i]],
                               species_no = all_sp_no[[i]])
    lins[i, j+1] <- max(raw_lineages$species_no[raw_lineages$branch_times >= j])
    
  }
  
}

lins <- as.data.frame(lins)

cumsp <- colSums(lins)
cumsp <- rev(cumsp)

all_rates <- vector()
for (i in 2:length(lins)) {
  difference <- cumsp[i] - cumsp[i-1]
  all_rates[i] <- difference/900
}

clade_rates_sil <- data.frame(Events = all_rates, Time = -35:0)
#clade_rates_calc <- data.frame(Events = all_rates, Time = -35:0)

ggplot(clade_rates_calc[20:36,], aes(x = Time, y = Events)) +
  theme_bw() +
  geom_point(color = "#548F01FF", alpha = 0.5) +
  geom_smooth(colour = "#548F01FF", fill = "#548F01FF", alpha = 0.1) +
  geom_point(data = clade_rates_sil[20:36,], aes(x = Time, y = Events),
             color = "#CFA3EEFF", alpha = 0.5) +
  geom_smooth(data = clade_rates_sil[20:36,], aes(x = Time, y = Events),
             color = "#CFA3EEFF", alpha = 0.1) +
  ylab("Speciation by bedrock type") +
  theme(axis.title.y = element_text(size = 10))
