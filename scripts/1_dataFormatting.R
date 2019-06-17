### Format data for analyses ####

# Packages

library(dplyr)

library(sf)
library(units)

library(survival)
library(survminer)


### DATA ####

# Mortality data

data_muir <- read.csv("data/Muir_data.csv", header = TRUE, row.names = 1)

# Tree_ID: Unique identifier number taken in 1987
# Species: Species code
# Plot_id: New and corrected unique identifier number for the quadrat
# Old_plot: Old unique identifier number for the quadrat
# DBH_1987: Diameter at Breast Height measured in 1987
# X,Y: Spatial coordinates in EPSG:2145
# State: 0 - alive, 1 - dead
# BBD: Beech Barch Scale insect; 0 - not infested, 1 - Infested
# Year_death: Year of death 
# State_xxxx: State at each sampling time
# Ice_storm: Ice storm canopy loss damage 

# Spatial data

# Create sf object with coordinates
xy <- st_as_sf(data_muir, coords = c("X", "Y"), crs = 2145) 

zone_muir <- st_read(dsn = 'data/zone_etude.shp')

plot(st_geometry(xy))
plot(st_geometry(zone_muir), add=T, border="blue")
dev.off()

### SPATIAL MANIPULATION FOR COMPETITION INDEX ####

# compute a 10m buffer around each tree
tree_buffer <- st_buffer(xy, 10)

# intersect 10m buffer with muir study zone 
tree_inters <- st_intersection(zone_muir, tree_buffer)

# compute area of intersection
tree_inters$area_inzone <- st_area(tree_inters)

# compute proportion of area of intersection
tree_inters$prop_inzone <- tree_inters$area_inzone/max(tree_inters$area_inzone)



### Compute competition index ####

# Total tree density
cat("Compute tree density within a 10m radius...")
data_muir$Tree_density <- (st_is_within_distance(xy[,2], dist = 10, sparse = TRUE))
data_muir$Tree_density <- unlist(lapply(data_muir$Tree_density, function(x) length(x)))
cat("Done!\n")

# Hegyi index sum(DBH_comp / (DBH_subj * Dist+1))

cat("Compute Hegyi index within a 10m radius...")

distarb <- st_distance(st_geometry(xy))
distarb[which(distarb == set_units(0, m))] <- 0.1

distdf <- as.data.frame(distarb)
colnames(distdf) <- xy$Tree_ID
row.names(distdf) <- xy$Tree_ID

data_muir$Hegyi <- NA
data_muir$Hegyi_Acer <- NA
data_muir$Hegyi_Fagus <- NA
data_muir$Hegyi_Ostrya <- NA
data_muir$Hegyi_Tilia <- NA
data_muir$Hegyi_Tsuga <- NA

study_sp <- c('ASA', 'FGR', 'OVI', 'TAM','TCA')


for(id in xy$Tree_ID) { 
  
  subject <- which(tree_inters$Tree_ID==id)
  
  tmp <- tree_inters[subject,]

  dhp_1 <- tmp$DBH_1987
  dhp_comp <- tree_inters$DBH_1987[-subject]
  distij <- distdf[-subject, subject]
  
  # Only distances <= 10m
  tree10m <- distij <= set_units(10, m)
  
  distij10 <- distij[tree10m]
  dhp_comp10 <- dhp_comp[tree10m]
  
  # hegyi inside 10m buffer with correction for forest edge
  data_muir$Hegyi[subject] <- sum(dhp_comp10/(dhp_1*(distij10))) / tmp$prop_inzone
  
  hegyi_sp <- list()
  
  for(sp in study_sp){
    dhp_comp_sp <- tree_inters$DBH_1987[which(tree_inters$Species==sp & !(tree_inters$Tree_ID==id))]
    distij_sp <- distdf[which(tree_inters$Species==sp & !(tree_inters$Tree_ID==id)), subject]
    
    # Only distances <= 10m
    tree_sp10m <- distij_sp <= set_units(10, m)
    
    distij_sp10 <- distij_sp[tree_sp10m]
    dhp_comp_sp10 <- dhp_comp_sp[tree_sp10m]
    
    hegyi_tmp <- sum(dhp_comp_sp10/(dhp_1*(distij_sp10))) / tmp$prop_inzone
    
    hegyi_sp[[sp]] <- hegyi_tmp
    
  }
  data_muir[subject, c("Hegyi_Acer", "Hegyi_Fagus", "Hegyi_Ostrya", "Hegyi_Tilia",  "Hegyi_Tsuga")] <- unlist(hegyi_sp)
  
}

cat("Done!\n")

### DATA CLEANING ####

# For trees that are still alive at the end of the study, year of death should be the last year of sampling (to create the survival object)
data_muir$Year_death[is.na(data_muir$Year_death)] <- 2018

# replace NA by 0 in BBD
data_muir$BBD[is.na(data_muir$BBD)] <- 0

# Remove trees that all still alive after 1998 but were not measured for ice storm damage (in the hydric portion of the forest)
data_muir <- filter(data_muir, !(Year_death > 1998 & is.na(Ice_storm))) 

# Relevel Ice storm 

data_muir$Ice_storm[data_muir$Year_death < 1998] <- "0"
data_muir$Ice_storm <- as.factor(data_muir$Ice_storm)
levels(data_muir$Ice_storm) <- c("0", "1", "2", "3", "4", "4")

# Remove trees that died before 1987

data_muir <- subset(data_muir, Year_death > 1987)




data_muir <- data_muir %>% 
  select(Tree_ID, Species, State, Year_death, Initial_DBH = DBH_1987,
         Ice_storm, Tree_density, Hegyi,
         Hegyi_Acer, Hegyi_Fagus, Hegyi_Ostrya , 
         Hegyi_Tilia,  Hegyi_Tsuga,
         BBD, 
         X, Y) %>%
  filter(Species %in% study_sp)

saveRDS(data_muir, "data/data_muir_clean.RDS")


### SURVIVAL DF ####

surv_sp <- list()

for(sp in study_sp) {

  sp_tmp <- filter(data_muir, Species == sp)

  sp_tmp$surv_object <- Surv(time = sp_tmp$Year_death, event = sp_tmp$State)
  
  surv_sp[[sp]] <- sp_tmp
  
}

saveRDS(surv_sp, "data/surv_sp.RDS")

