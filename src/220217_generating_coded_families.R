#15 Feb 2022
#Script to get clades set up for DAISIE

###Set up ####

library(tidyverse)
library(here)
library(ape)
library(phytools)


tree <- read.tree(here("Data", "220215_provisional_Daisie_dated.tre"))
final_taxonset <- read.csv(here("Data", "220216_taxonset_with_tips.csv"))

final_taxonset[2574,16] <- 0 #Question marks in the dataset and book
final_taxonset[1322,16] <- 0

#Code tips as alpine non-alpine

final_taxonset <- final_taxonset %>% 
  mutate(Origins = case_when(Alpin %in% c(1,2) | Nival %in% c(1,2) ~ "High_elevation",
                     Alpin %in% c(0) & Nival %in% c(0) ~ "Low_elevation"
                     )
         ) %>% 
  mutate(Status = case_when(
                    (Alpin %in% c(1,2) | Nival %in% c(1,2)) & Subalpin %in% c(0, 1) & Montagnard == 0 ~ "Endemic",
                    Alpin %in% c(1,2) & ( Subalpin > 0 | Montagnard > 0 | Collineen > 0)  ~ "Non_endemic_MaxAge",
                    Alpin == 0 & ( Subalpin > 0 | Montagnard > 0 | Collineen > 0)  ~ "Mainland")
  )


#Set up for ancestral trait reconstruction
tip_classifications <- as.vector(final_taxonset$Status)
names(tip_classifications) <- final_taxonset$Tip

test_anc <- ace(tip_classifications, tree, type = "discrete", model = "ER")

plot(tree)
nodelabels(node=1:tree$Nnode+Ntip(tree),
           pie=test_anc$lik.anc, cex=0.5)

simtree <- make.simmap(tree, tip_classifications, model = "ER")

plot(simtree, type="fan",fsize=0.08,ftype="i")

#Test on smaller segment ##

small_tree <- extract.clade(tree, 3245)
plot(small_tree)

small_data <- final_taxonset %>% 
  filter(Tip %in% small_tree$tip.label)
status <- as.vector(small_data$Status)
names(status) <- small_data$Tip

test_anc <- ace(status, small_tree, type = "discrete", model = "ER")
plot(test_anc)

nodelabels(node=1:small_tree$Nnode+Ntip(small_tree),
           pie=test_anc$lik.anc, cex=0.5)
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=cols,cex=0.3)

#How many orders are there 

length(unique(final_taxonset$Famille)) 
setdiff(unique(final_taxonset$Taxon_family), unique(final_taxonset$Famille))

library("MonoPhy")

taxonomy <- data.frame(Tip = final_taxonset$Tip, Family = final_taxonset$Taxon_family)

monophyly <- AssessMonophyly(tree, taxonomy, verbosity = 15)
GetSummaryMonophyly(monophyly)
results <- GetResultMonophyly(monophyly)
results$Family

write.csv(results, file = "220217_Checking_monophyly_trimmed_tree.csv", row.names = T)


#Families
library(ggtree)
campanulaceae <- extract.clade(tree, 5056)
length(campanulaceae$tip.label)

plot(campanulaceae)

campanulaceae_data <- final_taxonset %>% 
  filter(Tip %in% campanulaceae$tip.label) %>% 
  select(Tip, Status)

g <- ggtree(campanulaceae, ladderize = F) 
g <- g %<+% campanulaceae_data

p <- g + geom_tiplab(hjust = 0, nudge_x = 1) +
  geom_tippoint(aes(colour = Status), size = 3) +
  scale_colour_manual(values = c("#D72000","#1BB6AF","#FFAD0A")) +
  theme(legend.position = "left") +
  xlim(0, 150) 
p

cols <- c("#D72000","#1BB6AF","#FFAD0A")
stat <- as.vector(campanulaceae_data$Status)
names(stat) <- campanulaceae_data$Tip
camp_ace <- ace(stat, campanulaceae, type = "discrete", model = "ER")

plot(campanulaceae)
nodelabels(node=1:campanulaceae$Nnode+Ntip(campanulaceae),
           pie=camp_ace$lik.anc,piecol=cols,cex=0.3)
tiplabels(pie=to.matrix(stat,sort(unique(stat))),piecol=cols,cex=0.3)


mtrees<-make.simmap(campanulaceae,
                    stat,
                    model="ER",
                    nsim=500)
pd<-summary(mtrees,plot=FALSE)


cols <- setNames(c("#D72000","#1BB6AF","#FFAD0A"),
                 c("Endemic", "Mainland", "Non_endemic_MaxAge"))
plot(pd,
     fsize=0.5,
     ftype="i", 
     colors = cols, 
     xlim = c(0,50), mar=c(3.1,0.1,0.1,0.1))

abline(v = 10)
axisPhylo(side = 1)
add.simmap.legend(colors=cols, fsize = 0.5, x = 1, y = 0.1)
