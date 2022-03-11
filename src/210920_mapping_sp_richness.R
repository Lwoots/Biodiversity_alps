##Sept 2021
#Using OriginsAlpes data to make a map of species richness based on the grid cells of Thiel Eganter 2011

require(rgdal)
require(sp)
require(here)
require(tidyverse)

raw_data <- read.csv(here("Data", "Thiel2011_coords.csv")) #Extracted from supplementary data

#Creating the grid
grid <- SpatialPoints(raw_data[, 2:3])
plot(grid)

squares <- SpatialPixels(grid, tolerance = 0.5)
plot(squares)

square_df <- SpatialPixelsDataFrame(
  points = grid,
  tolerance = 0.5,
  data = data.frame(raw_data[, "Code"]),
  proj4string = CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs")
) 

polygons <- as(square_df, "SpatialPolygonsDataFrame")

polygons@proj4string <- CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs")

load(here("Data", "CompleteOccurences.RData"))
occ <- Complete
rm(Complete)
occ <- occ[seq(1, 9749000, 100),]
species <- data.frame(occ[, "EUROMED"])
summary(species)
occ_sp <- SpatialPointsDataFrame(
  coords      = occ[, c("lon", "lat")],
  data        = species,
  proj4string = CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs")
)

#points(occ_sp)


occ_sp$QDS_label <- over(occ_sp, polygons)


QDS_labels <- polygons@data$raw_data....Code..  # or whatever the QDS label column is called
species <- sort(unique(occ_sp$occ....EUROMED..))  # or whatever the species column is called


#Create presence absence matrix

presence <- matrix(nrow = length(QDS_labels), ncol = length(species))
rownames(presence) <- QDS_labels
colnames(presence) <- species

for (i in 1:length(QDS_labels)) {
  for (j in 1:length(species)) {
    presence[i, j] <-
      QDS_labels[i] %in% unique(occ_sp$QDS_label[occ_sp$occ....EUROMED.. == species[j],])
    # i.e., "Is the ith QDS label in the set of QDS where the jth species is found?"
    # (and then storing that resulting logical in the matrix `presence`)
  }
  print((i/length(QDS_labels))*100)
} 

#Convert to binary

presence_binary <- presence

presence_binary[,] <- ifelse(presence_binary %in% c(TRUE), 1, 0)
presence_binary <- as.data.frame(presence_binary)
#write.csv(presence_binary, file = "220307_non_endemics_grid_distribution.csv")

presence_binary$sp_rich <- rowSums(presence_binary)
presence_binary$sp_rich

#Add to polygons

richness <- as.data.frame(presence_binary$sp_rich)
names(richness)[1] <- "sp_rich" 


pid <- sapply(slot(polygons, "polygons"), function(x) slot(x, "ID"))
rownames(richness) <- pid


rich_map <- SpatialPolygonsDataFrame(polygons, richness)

non_endemic_plot <- spplot(rich_map, col = "transparent",
       xlab = "Species richness of non-endemics")

#Need to redo with just high elevation species

#Join euromed numbers to final taxon set

euro_med <- read.csv(here("Data/UPDATE2022", "PhyloAlpsDirectCorrespondances.csv"))
euro_med[2642, 19] <- "Hierochloe australis" 
euro_med[2643, 19] <- "Hierochloe odorata"
euro_med[1831, 19] <- "Hippophae rhamnoides"
euro_med <- euro_med %>% 
  select(PTNameFk, EuroMedAcceptedName) %>% 
  filter(!is.na(EuroMedAcceptedName))
euro_med$EuroMedAcceptedName_formatted <- str_replace_all(euro_med$EuroMedAcceptedName,
                                                          pattern = c(" var. " = " ",
                                                                      " subsp. " = " "))

phylo <- read.csv(here("Data", "220302_taxonset_resolved.csv"))
phylo <- phylo %>% 
  filter(!is.na(EuroMedAcceptedName_formatted))

phylo <- left_join(phylo, euro_med, by = "EuroMedAcceptedName_formatted")

final_taxonset <- phylo %>% 
  mutate(Origins = case_when(Alpin %in% c(1,2) | Nival %in% c(1,2) ~ "High_elevation",
                             Alpin %in% c(0) & Nival %in% c(0) ~ "Low_elevation"
  )
  ) %>% 
  mutate(Status = case_when(
    (Alpin %in% c(1,2) | Nival %in% c(1,2)) & Subalpin %in% c(0, 1) & Montagnard == 0 ~ "Endemic",
    Alpin %in% c(1,2) & ( Subalpin > 1 | Montagnard > 0 | Collineen > 0)  ~ "Non_endemic_MaxAge",
    Alpin == 0 & ( Subalpin > 0 | Montagnard > 0 | Collineen > 0)  ~ "Mainland"
  )
  )

endemic <- final_taxonset %>% 
  filter(Status == "Endemic") %>% 
  select(PTNameFk)

non_endemic <- final_taxonset %>% 
  filter(Status == "Non_endemic_MaxAge") %>% 
  select(PTNameFk)
both <- final_taxonset %>% 
  filter(Status %in% c("Non_endemic_MaxAge", "Endemic")) %>% 
  select(PTNameFk)


occ <- Complete %>% 
  filter(EUROMED %in% c(non_endemic$PTNameFk, endemic$PTNameFk))

occ2 <- occ


#Little pcoa analysis
require(vegan)
require(ape)

presence_binary <- presence_binary[,1:842]
no_empty <- presence_binary[rowSums(presence_binary[])>0,]

com_mat <- as.matrix(no_empty)

jacc <- vegdist(com_mat, method = "jaccard")

pco_jacc <- ape::pcoa(jacc)
plot(pco_jacc)
ape::biplot.pcoa(pco_jacc)
