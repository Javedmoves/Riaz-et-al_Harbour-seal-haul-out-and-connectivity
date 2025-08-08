##############################################################################################################################################
## Shapefile creation and manipulaion
##############################################################################################################################################

## Read in relevant packages and set data path



library(tidyverse)
library(aniMotum)
library(hms)
library(viridis)
library(scales)
library(readr)
library(grid)
library(raster)
library(stringr)
library(lubridate)
library(sp)
library(diveMove)
library(sf)
library(ggpubr)
library(ggOceanMaps)
library(ggspatial)
library(data.table)
library(future)
library(future.apply)
library(adehabitatLT)
library(rnaturalearthhires)
library(geosphere)
library(stars)

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
# Overlap with management areas


##############################################################################################################################################
# Norway
##############################################################################################################################################

setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Shapefiles/Norway/")


NorwayInternalWaters_sf <- read_sf("eez_internal_waters.shp")

NorwayInternalWaters_Map <- ggplot() +
  # geom_sf(data = NorwayCommune_sf, aes(fill = kommunenavn) ) +
  geom_sf(data = NorwayInternalWaters_sf, aes(fill = geoname) ) +
  theme(legend.position = "none") 
NorwayInternalWaters_Map


NorwayCoastalWaters_sf <- read_sf("eez_12nm.shp")

NorwayCoastalWater_Map <- ggplot() +
  # geom_sf(data = NorwayCommune_sf, aes(fill = kommunenavn) ) +
  geom_sf(data = NorwayInternalWaters_sf, aes(fill = geoname) ) +
  
  geom_sf(data = NorwayCoastalWaters_sf, aes(fill = geoname) ) +
  theme(legend.position = "none") 
NorwayCoastalWater_Map

Norway_Skagerrak <- st_union(NorwayCoastalWaters_sf, NorwayInternalWaters_sf) %>% 
  st_make_valid() 


ggplot() +
  geom_sf(data = Norway_Skagerrak, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Cleaned Merged Boundaries")


write_sf(Norway_Skagerrak , "Norway_Skagerrak.shp")




##############################################################################################################################################
# Sweden
##############################################################################################################################################


setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Shapefiles/Sweden/")


SwedenEEZWaters_sf <- read_sf("eez.shp")


SwedenEEZWaters_Map <- ggplot() +
  geom_sf(data = SwedenEEZWaters_sf, aes(fill = geoname) ) +
  theme(legend.position = "none") 
SwedenEEZWaters_Map


SwedenInternalWaters_sf <- read_sf("eez_internal_waters.shp")

SwedenInternalWaters_Map <- ggplot() +
  geom_sf(data = SwedenInternalWaters_sf, aes(fill = geoname) ) +
  theme(legend.position = "none") 
SwedenInternalWaters_Map

SwedenCoastalWaters_sf <- read_sf("eez_12nm.shp")


SwedenCoastalWater_Map <- ggplot() +
  geom_sf(data = SwedenInternalWaters_sf, aes(fill = geoname) ) +
  geom_sf(data = SwedenCoastalWaters_sf, aes(fill = geoname) ) +
  theme(legend.position = "none") 
SwedenCoastalWater_Map


Sweden_Skagerrak <- st_union(SwedenCoastalWaters_sf, SwedenInternalWaters_sf) %>% 
  st_make_valid()



ggplot() +
  geom_sf(data = SwedenEEZWaters_sf, aes(fill = geoname) ) +
  
  geom_sf(data = Sweden_Skagerrak, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Cleaned Merged Boundaries") 
# coord_sf(xlim = c(10,13), ylim = c(56, 60), expand = FALSE)



write_sf(Sweden_Skagerrak , "Sweden_Skagerrak.shp")



##############################################################################################################################################
# Denmark
##############################################################################################################################################



setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Shapefiles/Denmark/")

DenmarkInternalWaters_sf <- read_sf("eez_internal_waters.shp")

DenmarkInternalWaters_Map <- ggplot() +
  geom_sf(data = DenmarkInternalWaters_sf, aes(fill = geoname) ) +
  theme(legend.position = "none") 
DenmarkInternalWaters_Map

DenmarkCoastalWaters_sf <- read_sf("eez_12nm.shp")

DenmarkCoastalWater_Map <- ggplot() +
  geom_sf(data = DenmarkInternalWaters_sf, aes(fill = geoname) ) +
  geom_sf(data = DenmarkCoastalWaters_sf, aes(fill = geoname) ) +
  theme(legend.position = "none") 
DenmarkCoastalWater_Map


Denmark_Skagerrak <- st_union(DenmarkCoastalWaters_sf, DenmarkInternalWaters_sf) %>% 
  st_make_valid() 

ggplot() +
  geom_sf(data = Denmark_Skagerrak, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Cleaned Merged Boundaries")


write_sf(Denmark_Skagerrak , "Denmark_Skagerrak.shp")



##############################################################################################################################################
# Combine
##############################################################################################################################################

setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Shapefiles/Norway/")
Norway_Skagerrak <- read_sf("Norway_Skagerrak.shp")
Norway_Skagerrak$Area <- "Norway"
ggplot() +
  geom_sf(data = Norway_Skagerrak, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Cleaned Merged Boundaries")


setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Shapefiles/Sweden/")
Sweden_Skagerrak <- read_sf("Sweden_Skagerrak.shp")
Sweden_Skagerrak$Area <- "Sweden"


setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Shapefiles/Denmark/")
Denmark_Skagerrak <- read_sf("Denmark_Skagerrak.shp")
Denmark_Skagerrak$Area <- "Denmark"
ggplot() +
  geom_sf(data = Denmark_Skagerrak, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Cleaned Merged Boundaries")



setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Shapefiles/")

Combined_Maritime_Boundaries<- bind_rows(Norway_Skagerrak, Sweden_Skagerrak,Denmark_Skagerrak)
ggplot() +
  geom_sf(data = Combined_Maritime_Boundaries, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Cleaned Merged Boundaries")





Swedend_Denmark_Combined_Maritime_Boundaries<- bind_rows(Sweden_Skagerrak,Denmark_Skagerrak)
ggplot() +
  geom_sf(data = Swedend_Denmark_Combined_Maritime_Boundaries, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Cleaned Merged Boundaries")




write_sf(Combined_Maritime_Boundaries, "Combined_Maritime_Boundaries.shp")

write_sf(Swedend_Denmark_Combined_Maritime_Boundaries, "Swedend_Denmark_Combined_Maritime_Boundaries.shp")



##############################################################################################################################################
# Bathymetry
##############################################################################################################################################

setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Shapefiles/")


bbox <- st_bbox(c(xmin = 6.9, xmax = 12.31, ymin = 56.19, ymax = 59.8), crs = 4326)

bathy_data <- read_stars("EMODnet_bathymetry_2022.nc", proxy = TRUE)

st_crs(bathy_data)

bathy_data <- bathy_data %>%
  st_set_crs(st_crs(4326))


bathy_crop <- st_crop(bathy_data, bbox)

write_sf(bathy_df, "Skagerrak_Bathymetry_EMODnet_115m.shp")


bathy_df <- as.data.frame(bathy_crop, xy = TRUE)

colnames(bathy_df)[3] <- "depth"

bathy_df$Depth <- as.numeric(bathy_df$depth)

bathy_df <- bathy_df  %>%
  dplyr::select(x, y, Depth)


write_rds(bathy_df, "Skagerrak_Bathymetry_EMODnet_115m.rds")


bathy_sf <- st_as_sf(bathy_df, coords = c("x", "y"), crs = 4326)
st_write(bathy_sf, "Skagerrak_Bathymetry_EMODnet_115m.shp")

write_sf(bathy_df, "Skagerrak_Bathymetry_EMODnet_115m.shp")































