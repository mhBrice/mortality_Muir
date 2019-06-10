
# Library
library(dplyr)

library(sf)
library(units)

library(survival)
library(survminer)


# Mortality data

data_muir <- read.csv("data/Muir_data.csv", header = TRUE, dec = ".", sep = ",")
data_muir$Annee_mort[is.na(data_muir$Annee_mort)] <- 2018

# Spatial data

xy <- st_as_sf(data_muir, coords = c("coordo_x", "coordo_y"), crs = 2145) #produire un fichier de coordonnees 

zone_muir <- st_read(dsn = 'data/zone_etude.shp')

plot(st_geometry(xy))
plot(st_geometry(zone_muir), add=T, border="blue")

### Spatial manipulation

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
data_muir$index_densite <- (st_is_within_distance(xy[,2], dist = 10, sparse = TRUE))
data_muir$Tree_density <- lapply(data_muir$index_densite, function(x) length(x))

data_muir$index_densite <- NULL

# Hegyi index sum(DBH_comp / (DBH_subj * Dist+1))

distarb <- st_distance(st_geometry(xy))
distarb[which(distarb == set_units(0, m))] <- 0.1

distdf <- as.data.frame(distarb)
colnames(distdf) <- xy$N_OID_19XX
row.names(distdf) <- xy$N_OID_19XX

data_muir$Hegyi = NA
data_muir$Hegyi_Acer = NA
data_muir$Hegyi_Fagus = NA
data_muir$Hegyi_Ostrya = NA
data_muir$Hegyi_Tilia = NA
data_muir$Hegyi_Tsuga = NA

study_sp = c('ASA', 'FGR', 'OVI', 'TAM','TCA')


for(id in xy$N_OID_19XX) { 
  
  subject <- which(tree_inters$N_OID_19XX==id)
  
  tmp <- tree_inters[subject,]

  dhp_1 = tmp$DHP1987
  dhp_comp = tree_inters$DHP1987[-subject]
  distij = distdf[-subject, subject]
  
  # prendre juste les distance en <= 10m
  tree10m <- distij <= set_units(10, m)
  
  distij10 <- distij[tree10m]
  dhp_comp10 <- dhp_comp[tree10m]
  
  # hegyi inside 10m buffer with correction for forest edge
  data_muir$Hegyi[subject] = sum(dhp_comp10/(dhp_1*(distij10))) / tmp$prop_inzone
  
  hegyi_sp = list()
  
  for(sp in study_sp){
    dhp_comp_sp = tree_inters$DHP1987[which(tree_inters$ESPECE==sp & !(tree_inters$N_OID_19XX==id))]
    distij_sp = distdf[which(tree_inters$ESPECE==sp & !(tree_inters$N_OID_19XX==id)), subject]
    
    # prendre juste les distance en <= 10m
    tree_sp10m <- distij_sp <= set_units(10, m)
    
    distij_sp10 <- distij_sp[tree_sp10m]
    dhp_comp_sp10 <- dhp_comp_sp[tree_sp10m]
    
    hegyi_tmp = sum(dhp_comp_sp10/(dhp_1*(distij_sp10))) / tmp$prop_inzone
    
    hegyi_sp[[sp]] = hegyi_tmp
    
  }
  data_muir[subject, c("Hegyi_Acer", "Hegyi_Fagus", "Hegyi_Ostrya", "Hegyi_Tilia",  "Hegyi_Tsuga")] = unlist(hegyi_sp)
  
}

data_muir$Tree_density <- as.numeric(data_muir$Tree_density)
data_muir$ESPECE <- as.factor(data_muir$ESPECE)
data_muir$Annee_mort <- as.numeric(data_muir$Annee_mort)
data_muir <- mutate(data_muir, Time2 = (Annee_mort - 1987))


#seulement 1968 arbres de nos especes d'interets

data_muir$MCH96[is.na(data_muir$MCH96)] <- 0


#Eliminer les arbres qui n'ont pas ete mesure pour le verglas (donc, annee_mort >1998,Verglas cote = NA)

data_muir <- filter(data_muir, !(Annee_mort > 1998 & is.na(Verglas.cote1998))) %>% 
  rename(Ice_storm = Verglas.cote1998)

# Relevel verglas

data_muir$Ice_storm[data_muir$Annee_mort<1998] <- "0"
data_muir$Ice_storm <- as.factor(data_muir$Ice_storm)
levels(data_muir$Ice_storm) <- c("0", "1", "2", "3", "4", "4")

# Remove trees that died before 1987

data_muir <- subset(data_muir, Annee_mort > 1987)





data_muir <- select(data_muir, N_OID_19XX, ESPECE, ETAT, Annee_mort, DHP1987,
                                  Ice_storm, Tree_density, Hegyi,
                                  Hegyi_Acer, Hegyi_Fagus, Hegyi_Ostrya , 
                                  Hegyi_Tilia,  Hegyi_Tsuga,
                                  MCH96, 
                                  coordo_x, coordo_y, Time2) %>%
  rename(ID = N_OID_19XX, Initial_DBH = DHP1987, BBD = MCH96) %>% 
  filter(ESPECE %in% study_sp)

saveRDS(data_muir, "data/data_muir_survie.RDS")


### SURVIVAL DF ####

surv_sp <- list()

for(sp in study_sp) {
  # species
  sp_tmp <- filter(data_muir, ESPECE == sp)
  
  # coordinates
  xy <- sp_tmp[ ,c('coordo_x','coordo_y')]
  
  sp_tmp$surv_object <- Surv(time = sp_tmp$Annee_mort, event = sp_tmp$ETAT)
  
  sp_tmp <- subset(sp_tmp, Time2 > 0)
  
  surv_sp[[sp]] <- sp_tmp
  
}

saveRDS(surv_sp, "data/surv_sp.RDS")

