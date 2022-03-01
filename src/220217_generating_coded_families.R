#15 Feb 2022
#Script to get clades set up for DAISIE

###Set up ####

library(tidyverse)
library(here)
library(ape)
library(phytools)


tree <- read.tree(here("Data", "220221_for_Daisie_dated.tre"))
final_taxonset <- read.csv(here("Data", "220221_taxonset_with_tips.csv"))

final_taxonset[2574,16] <- 0 #Question marks in the dataset and book
final_taxonset[1322,16] <- 0

#Code tips as alpine non-alpine

mainland <- c("BGN_IBH", "BGN_LFB", "BGN_DTT", "BGN_KHP", "BGN_CDD",
              "BGN_GPM", "BGN_FEM", "BGN_KGA", "BGN_APT", "BGN_IGI") #These are the manually added taxa
non_endem <- c("BGN_DII", "BGN_KIE", "BGN_GQR")
alpine <- c("BGN_PQF","GWM_1539", "GWM_1231")

final_taxonset <- final_taxonset %>% 
  mutate(Origins = case_when(Alpin %in% c(1,2) | Nival %in% c(1,2) ~ "High_elevation",
                     Alpin %in% c(0) & Nival %in% c(0) ~ "Low_elevation"
                     )
         ) %>% 
  mutate(Status = case_when(
                    (Alpin %in% c(1,2) | Nival %in% c(1,2)) & Subalpin %in% c(0, 1) & Montagnard == 0 ~ "Endemic",
                    Alpin %in% c(1,2) & ( Subalpin > 1 | Montagnard > 0 | Collineen > 0)  ~ "Non_endemic_MaxAge",
                    Alpin == 0 & ( Subalpin > 0 | Montagnard > 0 | Collineen > 0)  ~ "Mainland",
                    Sequencing_ID %in% mainland ~ "Mainland",
                    Sequencing_ID %in% non_endem ~ "Non_endemic_MaxAge",
                    Sequencing_ID %in% alpine ~ "Endemic")
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

#

#How many orders are there 

length(unique(final_taxonset$Famille)) 
setdiff(unique(final_taxonset$Taxon_family), unique(final_taxonset$Famille))

library("MonoPhy")

taxonomy <- data.frame(Tip = final_taxonset$Tip, Family = final_taxonset$Taxon_family)

monophyly <- AssessMonophyly(tree, taxonomy, verbosity = 15)
GetSummaryMonophyly(monophyly)
results <- GetResultMonophyly(monophyly)
results$Family

write.csv(results, file = "220221_Checking_monophyly_trimmed_tree.csv", row.names = T)

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
#Wrong in tree: Linnaea_borealis BGN_KIS Hedera helix RSZ_RSZAXPI001130-87 Polemonium_caeruleum RSZ_RSZAXPI000710-107
#Osyris_alba

final_taxonset$Status[final_taxonset$Tip == "Dianthus_barbatus_278075_PHA002888_BGN_EDR"] <- "Mainland"
final_taxonset$Status[final_taxonset$Tip == "Pimpinella_alpina_1533290_PHA006713_BGN_HHH"] <- "Mainland"

tree <- drop.tip(tree, c("Oxalis_debilis_519149_PHA006334_BGN_VT"))
final_taxonset <- final_taxonset %>% filter(!Tip == "Oxalis_debilis_519149_PHA006334_BGN_VT")

#Make plots of all the families
#

node <-6483
clade <- extract.clade(tree, node)
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
     xlim = c(0,20), mar=c(3.1,0.1,0.1,0.1), cex = 0.25)
axisPhylo(side = 1)
abline(v = 42)
add.simmap.legend(colors=cols, fsize = 0.5, x = 1, y = 0.1)



clade_ace <- ace(stat, clade, type = "discrete", model = "ER")

plotTree(clade, cex = 0.5, xlim = c(0,90))

nodelabels(node=1:clade$Nnode+Ntip(clade),
           pie=clade_ace$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(stat,sort(unique(stat))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)

library(ggtree)

g <- ggtree(clade, ladderize = F) %<+% clade_data


p <- g + geom_tiplab(hjust = 0, nudge_x = 1) +
  geom_tippoint(aes(color = Status), size = 3) +
  scale_colour_manual(values = c("#D72000","#1BB6AF","#FFAD0A")) +
  theme(legend.position = "n") +
  xlim(0,25)
p


#Getting data into daisie format

tree_order <- data.frame(Tip=tree$tip.label)

rev(tree_order)
final_taxonset$Tip[1:20]

ordered_by_tree <- left_join(rev(tree_order), final_taxonset, by = "Tip") #Organising so is in the same order as the tree

ordered_by_tree$Tip[1:20] 

#Code non-endemics

#First name clade of non-endemics

ordered_by_tree <- ordered_by_tree %>% 
  mutate(Clade_name =
           case_when(
             Status == "Non_endemic_MaxAge" ~ paste(Genre, row.names(.), sep = "_")))
ordered_by_tree <- ordered_by_tree %>% 
  select(Tip, EuroMedAcceptedName_formatted, Taxon_family, Collineen, Montagnard, Alpin, Nival,
         Origins, Status, Clade_name)
write.csv(ordered_by_tree, file = "220222_coding_endemics.csv", row.names = F)

