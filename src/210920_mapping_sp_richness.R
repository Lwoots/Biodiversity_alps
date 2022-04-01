##Sept 2021
#Using OriginsAlpes data to make a map of species richness based on the grid cells of Thiel Eganter 2011

require(rgdal)
require(sp)
require(here)
require(tidyverse)
require(vegan)

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
#rm(Complete)
#occ <- occ[seq(1, 9749000, 100),]
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
write.csv(presence_binary, file = "220315_non_endemics_and_endemics_grid_distribution.csv")

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
euro_med <- euro_med %>% 
  distinct(PTNameFk, .keep_all = T)

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
  filter(EUROMED %in% c(both$PTNameFk))

occ2 <- occ


#Little pcoa analysis ####
require(vegan)
require(ape)
require(viridis)


presence_binary <- read.csv(here("Data", "220315_non_endemics_and_endemics_grid_distribution.csv"))
presence_binary$lon <- grid@coords[1:385]
presence_binary$lat <- grid@coords[386:770]
row.names(presence_binary) <- presence_binary$X


#presence_binary <- presence_binary[,1:842]
just_sp <- presence_binary[,2:844]
no_empty <- just_sp[rowSums(just_sp[])>0,]

com_mat <- as.matrix(no_empty)

jacc <- vegdist(com_mat, method = "jaccard")

pco_jacc <- ape::pcoa(jacc)
plot(pco_jacc)
ape::biplot.pcoa(pco_jacc)

empty_grids <- setdiff(presence_binary$X, row.names(no_empty))

pco_jacc$values
final <- presence_binary %>% 
  filter(!X %in% empty_grids) %>% 
  mutate(Eigenvec = pco_jacc$values$Eigenvalues,
         Rel_eigen = pco_jacc$values$Relative_eig,
         axis1 = pco_jacc$vectors[,1],
         axis2 = pco_jacc$vectors[,2],
         axis3 = pco_jacc$vectors[,3]
         )

a <- ggplot(final, aes(lon, lat, color = axis1)) +
  geom_point(size = 4, shape = "square") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(title = "Axis 1, 45.5%") +
  scale_color_viridis(option = "inferno") +
  theme_bw() 
b <- ggplot(final, aes(lon, lat, color = axis2)) +
  geom_point(size = 4, shape = "square")+
  scale_color_viridis(option = "inferno")+
  xlab("Longitude") +
  ylab("Latitude") +
  labs(title = "Axis 2, 8.2%") +
  theme_bw()
c <- ggplot(final, aes(lon, lat, color = axis3)) +
  geom_point(size = 4, shape = "square") +
  scale_color_viridis(option = "inferno")+
  xlab("Longitude") +
  ylab("Latitude") +
  labs(title = "Axis 3, 6.6%") +
  theme_bw()


abc <- ggpubr::ggarrange(a,b,c, ncol = 1)
ggpubr::annotate_figure(abc, 
                        bottom = text_grob("PCOA of endemic and non-endemic species, \n based on OriginAlps data", 
                                           color = "black",  size = 12))


#Redo with Flora alpina data ####

flora_provences <- read.csv(here("Data", "220304_Provinces_sorted_to_regions.csv"), sep = ";")
flora_provences[flora_provences$Department == "0_OB", 1] <- "X0_OB"
flora_provences[flora_provences$Department == "0_SW", 1] <- "X0_SW"
names(flora_provences)[3] <- "Regions2"

final_taxonset <- read.csv(here("Data", "220302_taxonset_resolved.csv"))

final_taxonset <- final_taxonset %>% 
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


dist_island <- final_taxonset %>% 
  filter(Status %in% c("Endemic","Non_endemic_MaxAge")) %>% 
  drop_na() %>% 
  select(A_B:SLO) 

dist_island <- t(dist_island)
jacc_isl <- vegdist(dist_island, method = "jaccard")

pco_isl <- ape::pcoa(jacc_isl)
biplot.pcoa(pco_isl)



dist_island_df <- as.data.frame(dist_island)
dist_island_df <- dist_island_df %>% 
  mutate(Department = row.names(dist_island_df),
         axis1 = pco_isl$vectors[,1],
         axis2 = pco_isl$vectors[,2],
         axis3 = pco_isl$vectors[,3],
         axis4 = pco_isl$vectors[,4])
dist_island_df$Department

combined_isl <- left_join(dist_island_df, flora_provences, by = "Department")

require(khroma)
require(ggpubr)
my_theme <- theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        legend.key = element_rect(fill = NA, color = "grey"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 40, vjust = 0.9, hjust = 1, colour = "black"),
        axis.title = element_text(size = 12)
  )

a <- ggplot(combined_isl, aes(axis1, axis2, colour = Region, label = Department)) +
  geom_point(size = 3) +
  geom_text(size = 3, hjust = -0.18) +
  my_theme +
  scale_color_highcontrast(reverse = T) +
  xlab("Axis 1, 31.9%") +
  ylab("Axis 2, 16.5%") +
  scale_x_continuous(limits = c(-0.2, 0.6))
b <- ggplot(combined_isl, aes(axis1, axis3, colour = Region, label = Department)) +
  geom_point(size = 3)+
  my_theme +
  scale_color_highcontrast(reverse = T) +
  geom_text(size = 3, hjust = -0.18) +
  xlab("Axis 1, 31.9%") +
  ylab("Axis 3, 9.9%")+
  scale_x_continuous(limits = c(-0.2, 0.6))

ab <- ggpubr::ggarrange(a, b, nrow = 1, common.legend = T )
ggpubr::annotate_figure(ab, bottom = text_grob("PCOA of endemic and non-endemic species", 
                                      color = "black",  size = 14))#6X10


c <- ggplot(combined_isl, aes(axis1, axis2, colour = Regions2, label = Department)) +
  geom_point(size = 3)+
  my_theme +
  geom_text(size = 3, hjust = -0.18) +
  scale_color_highcontrast(reverse = T)+
  xlab("Axis 1, 31.9%") +
  ylab("Axis 2, 16.5%")+
  scale_x_continuous(limits = c(-0.2, 0.6))
d <- ggplot(combined_isl, aes(axis1, axis3, colour = Regions2, label = Department)) +
  geom_point(size = 3)+
  my_theme +
  geom_text(size = 3, hjust = -0.18) +
  scale_color_highcontrast(reverse = T)+
  xlab("Axis 1, 31.9%") +
  ylab("Axis 3, 9.9%")+
  scale_x_continuous(limits = c(-0.2, 0.6))


ab <- ggpubr::ggarrange(a, b, c, d, nrow = 2, ncol = 2, common.legend = T, labels = "AUTO" )
ggpubr::annotate_figure(ab, 
                        bottom = text_grob("PCOA of endemic and non-endemic species, based on FA.\n Rows show different department classifications", 
                                               color = "black",  size = 10))#7x7


#Remove small outlier departments

dist_island <- final_taxonset %>% 
  filter(Status %in% c("Endemic","Non_endemic_MaxAge")) %>% 
  drop_na() %>% 
  select(A_B:SLO) %>% 
  select(!c(I_IM, F_84, F_83, I_SV, A_B, I_VA))

dist_island <- t(dist_island)

jacc_isl <- vegdist(dist_island, method = "jaccard")

pco_isl <- ape::pcoa(jacc_isl)
biplot.pcoa(pco_isl)
pco_isl$values
dist_island_df <- as.data.frame(dist_island)
dist_island_df <- dist_island_df %>% 
  mutate(Department = row.names(dist_island_df),
         axis1 = pco_isl$vectors[,1],
         axis2 = pco_isl$vectors[,2],
         axis3 = pco_isl$vectors[,3],
         axis4 = pco_isl$vectors[,4])
dist_island_df$Department

combined_isl <- left_join(dist_island_df, flora_provences, by = "Department")

a <- ggplot(combined_isl, aes(axis1, axis2, colour = Region, label = Department)) +
  geom_point(size = 3) +
  geom_text(size = 3, hjust = -0.18) +
  my_theme +
  scale_color_highcontrast(reverse = T) +
  xlab("Axis 1, 30.3%") +
  ylab("Axis 2, 17.2%") +
  scale_x_continuous(limits = c(-0.3, 0.5))
b <- ggplot(combined_isl, aes(axis1, axis3, colour = Region, label = Department)) +
  geom_point(size = 3)+
  my_theme +
  scale_color_highcontrast(reverse = T) +
  geom_text(size = 3, hjust = -0.18) +
  xlab("Axis 1, 30.3%") +
  ylab("Axis 3, 14.4%")+
  scale_x_continuous(limits = c(-0.3, 0.5))

c <- ggplot(combined_isl, aes(axis1, axis2, colour = Regions2, label = Department)) +
  geom_point(size = 3)+
  my_theme +
  geom_text(size = 3, hjust = -0.18) +
  scale_color_highcontrast(reverse = T)+
  xlab("Axis 1, 30.3%") +
  ylab("Axis 2, 17.2%")+
  scale_x_continuous(limits = c(-0.3, 0.5))
d <- ggplot(combined_isl, aes(axis1, axis3, colour = Regions2, label = Department)) +
  geom_point(size = 3)+
  my_theme +
  geom_text(size = 3, hjust = -0.18) +
  scale_color_highcontrast(reverse = T)+
  xlab("Axis 1, 30.3%") +
  ylab("Axis 3, 14.4%")+
  scale_x_continuous(limits = c(-0.3, 0.5))


ab <- ggpubr::ggarrange(a, b, c, d, nrow = 2, ncol = 2, common.legend = T, labels = "AUTO" )
ggpubr::annotate_figure(ab, 
                        bottom = text_grob("PCOA of endemic and non-endemic species, based on FA.\n Small outlier departments removed. Rows show different department classifications", 
                                           color = "black",  size = 10))#7x7

#For all species

distributions <- final_taxonset %>% 
  select(A_B:SLO) %>% 
  drop_na()
distributions <- t(distributions)

jacc_flo <- vegdist(distributions, method = "jaccard")

pco_flo <- ape::pcoa(jacc_flo)

biplot.pcoa(pco_flo)
pco_flo$values

dist_df <- as.data.frame(distributions)
dist_df <- dist_df %>% 
  mutate(Department = row.names(dist_island_df),
         axis1 = pco_isl$vectors[,1],
         axis2 = pco_isl$vectors[,2],
         axis3 = pco_isl$vectors[,3],
         axis4 = pco_isl$vectors[,4])
dist_df$Department

combined_flo <- left_join(dist_df, flora_provences, by = "Department")

a <- ggplot(combined_isl, aes(axis1, axis2, colour = Region, label = Department)) +
  geom_point(size = 3) +
  geom_text(size = 3, hjust = -0.18) +
  my_theme +
  scale_color_highcontrast(reverse = T) +
  xlab("Axis 1, 23.9%") +
  ylab("Axis 2, 15.0%") +
  scale_x_continuous(limits = c(-0.25, 0.6))
b <- ggplot(combined_isl, aes(axis1, axis3, colour = Region, label = Department)) +
  geom_point(size = 3)+
  my_theme +
  scale_color_highcontrast(reverse = T) +
  geom_text(size = 3, hjust = -0.18) +
  xlab("Axis 1, 23.9%") +
  ylab("Axis 3, 11.0%")+
  scale_x_continuous(limits = c(-0.25, 0.6))

c <- ggplot(combined_isl, aes(axis1, axis2, colour = Regions2, label = Department)) +
  geom_point(size = 3)+
  my_theme +
  geom_text(size = 3, hjust = -0.18) +
  scale_color_highcontrast(reverse = T)+
  xlab("Axis 1, 23.9%") +
  ylab("Axis 2, 15.0%")+
  scale_x_continuous(limits = c(-0.25, 0.6))
d <- ggplot(combined_isl, aes(axis1, axis3, colour = Regions2, label = Department)) +
  geom_point(size = 3)+
  my_theme +
  geom_text(size = 3, hjust = -0.18) +
  scale_color_highcontrast(reverse = T)+
  xlab("Axis 1, 23.9%") +
  ylab("Axis 3, 11.0%")+
  scale_x_continuous(limits = c(-0.25, 0.6))


ab <- ggpubr::ggarrange(a, b, c, d, nrow = 2, ncol = 2, common.legend = T, labels = "AUTO" )
ggpubr::annotate_figure(ab, 
                        bottom = text_grob("PCOA of full data set, based on FA.\n Rows show different department classifications", 
                                           color = "black",  size = 10))#7x7


#Hclust ####

require(dendextend)

colorCodes <- c(Central="red", West="blue", East="yellow")

hclust_isl <- hclust(jacc_isl, method = "average") 
plot(hclust_isl)

labs <- left_join(data.frame(Department = hclust_isl$labels), flora_provences, by = "Department")
labs$code <- paste(labs$Region, labs$Department)
labs <- labs %>% 
  mutate(colours = case_when(
    Region == "Western Alps" ~ "blue",
    Region == "Central Alps" ~ "red",
    Region == "Eastern Alps" ~ "orange"
  ))

hclust_isl$labels <- labs$code
labs$colours
plot(hclust_isl, cex = 0.5)

dend <- as.dendrogram(hclust_isl)
dend <- dend %>% 
  set("labels_colors", c(7,7,7,2,4,2,7,7,7,7,7,7,4,4,4,4,
                         2,4,4,4,4,4,4,4,2,2,2,2,2,2,2,2,2,
                         2,2,2,7,7,7,7,7,7,7,7,2,2,2,2,2,7,
                         4,4,4,4,4)) %>% 
  assign_values_to_leaves_nodePar(0.5,"lab.cex")
  
plot(dend)

nP <- list(col = 3:2, cex = c(2.0, 0.75), pch =  21:22,
           lab.cex = 0.75, lab.col = c(2,2))

