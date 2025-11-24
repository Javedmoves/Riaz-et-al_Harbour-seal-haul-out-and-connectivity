##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

# MANUSCRIPT TITLE: Haul-out site use and connectivity of harbour seals between management units in southern Scandinavia

# AUTHOR OF CODE: Javed Riaz

# VERSION DATE: 25/11/2024 

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################


library(tidyverse)
library(hms)
library(viridis)
library(scales)
library(readr)
library(grid)
library(raster)
library(stringr)
library(lubridate)
library(sp)
library(sf)
library(ggpubr)
library(ggOceanMaps)
library(ggspatial)
library(data.table)
library(rnaturalearthhires)
library(ggnewscale)
library(geosphere)
library(glmmTMB)
library(broom.mixed)
library(sjPlot)
library(DHARMa)
library(survival)
library(survminer)
library(spatsoc)
library(mgcv)
library(modelsummary)
library(gt)
library(cowplot)
library(magick)
library(flextable)
library(webshot2)
library(grid)
library(broom)
library(patchwork)
library(DescTools)
library(ggsci)
library(ggsurvfit)
library(ggarchery)
library(leaflet)
library(htmlwidgets)
library(htmltools)
library(asnipe)
library(igraph)
library(RColorBrewer)
library(ggraph)
library(tidygraph)
library(plotrix)
library(RColorBrewer)
library(spatsoc)
library(datawizard)
library(gratia)
library(performance)
library(coxme)
library(kableExtra)

Sys.setenv(LANG = "en")


# Some helpful functions for plotting and calculations that we'll need later. Folded for brevity

annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{layer(data = data, stat = StatIdentity, position = PositionIdentity,geom = ggplot2:::GeomCustomAnn,
       inherit.aes = TRUE, params = list(grob = grob,xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))}


centroid_group <- function(
    DT = NULL,
    coords = NULL,
    group = 'group',
    na.rm = FALSE) {
  
  if (is.null(DT)) {
    stop('input DT required')
  }
  
  if (length(coords) != 2) {
    stop('coords requires a vector of column names for coordinates X and Y')
  }
  
  if (is.null(group)) {
    stop('group column name required')
  }
  
  if (any(!(
    c(coords, group) %in% colnames(DT)
  ))) {
    stop(paste0(
      as.character(paste(setdiff(
        c(coords, group),
        colnames(DT)
      ), collapse = ', ')),
      ' field(s) provided are not present in input DT'
    ))
  }
  
  if (any(!(DT[, vapply(.SD, is.numeric, TRUE), .SDcols = coords]))) {
    stop('coords must be numeric')
  }
  
  if (is.null(na.rm)) {
    stop('na.rm is required')
  }
  
  if (!is.logical(na.rm)) {
    stop('na.rm should be a boolean (TRUE/FALSE), see ?mean')
  }
  
  xcol <- data.table::first(coords)
  ycol <- data.table::last(coords)
  
  out_xcol <- paste0('centroid_', gsub(' ', '', xcol))
  out_ycol <- paste0('centroid_', gsub(' ', '', ycol))
  
  if (out_xcol %in% colnames(DT)) {
    message(paste(out_xcol, 'column will be overwritten by this function'))
    data.table::set(DT, j = out_xcol, value = NULL)
  }
  
  if (out_ycol %in% colnames(DT)) {
    message(paste(out_ycol, 'column will be overwritten by this function'))
    data.table::set(DT, j = out_ycol, value = NULL)
  }
  
  DT[, c(out_xcol) := mean(.SD[[xcol]], na.rm = na.rm), by = c(group)]
  DT[, c(out_ycol) := mean(.SD[[ycol]], na.rm = na.rm), by = c(group)]
  
  return(DT[])
}


arith_mean_ci <- function(x, conf.level = 0.95) {
  x <- x[!is.na(x)]
  mean_x <- mean(x)
  se_x <- sd(x) / sqrt(length(x))
  z <- qnorm(1 - (1 - conf.level) / 2)
  lower <- mean_x - z * se_x
  upper <- mean_x + z * se_x
  list(mean = mean_x, lower_CI = lower, upper_CI = upper)
}


geo_mean_ci <- function(x, conf.level = 0.95) {
  x <- x[!is.na(x) & x > 0]  
  log_x <- log(x)            
  mean_log <- mean(log_x)
  se_log <- sd(log_x) / sqrt(length(log_x))
  
  z <- qnorm(1 - (1 - conf.level) / 2)
  lower_log <- mean_log - z * se_log
  upper_log <- mean_log + z * se_log
  list(
    geometric_mean = exp(mean_log),
    lower_CI = exp(lower_log),
    upper_CI = exp(upper_log)
  )
}


calculate_distance <- function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000  # Convert meters to kilometers
}

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
# Read in relevant data


# setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Harbour seal movement/GPS processing/")
setwd("C:/Users/AADDT/Documents/Javed/Haulouts/Haulouts/Data files/Seal Data/")




QCedRAWDATA <- read_rds("QCed Data.rds")


Seal_GPS_SSM_Clean <- read_rds("Seal_GPS_SSM_Clean.rds")


SupMeta <- readxl::read_xlsx("Seal summary data 2017-2022.xlsx")


SupMeta <- SupMeta %>% 
  dplyr::select(2, 3, PTT, REF, `Weight (kg)`, 14, 'N days', `Tagging year`, Sex, `Capture site`, 6, 7, 8) %>% 
  rename("Capture_site" = "Capture site",
         'Col_lat' = "Col lat",
         "Col_lon" = "Col lon") %>% 
  rename("id" = "REF") %>% 
  rename("TagDays" = "N days") %>% 
  rename("TaggingYear" = `Tagging year`)


QCedRAWDATA <- QCedRAWDATA %>% 
  dplyr::select(- lc) %>% 
  distinct() %>% 
  left_join(SupMeta)



MetaData <- QCedRAWDATA %>% 
  dplyr::select(1,5, 7:15) %>% 
  distinct() 
# left_join(SupMeta)


# write_rds(SupMeta ,"SupMeta.rds")


MetaData1 <- MetaData %>% 
  dplyr::select(Capture_site, Col_lon, Col_lat) %>% 
  distinct() %>% 
  group_by(Capture_site) %>%
  filter(row_number()==1) %>% 
  mutate(colony = Capture_site) %>% 
  mutate(Level = "All tracking data")


# setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Harbour seal movement/Haulouts/")

pb35 <- readxl::read_xlsx("haulout1.xlsx")
pb68 <- readxl::read_xlsx("haulout2.xlsx")
# pb76 <- readxl::read_xlsx("haulout3.xlsx")
pb77 <- readxl::read_xlsx("haulout4.xlsx")
pb78 <- readxl::read_xlsx("haulout5.xlsx")
pb2020 <- read_rds("pb2020_REV.rds")


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
# Compile and clean data

AllHauloutData <- pb35 %>% 
  full_join(pb68) %>% 
  # full_join(pb76) %>% 
  full_join(pb77) %>% 
  full_join(pb2020) %>%
  full_join(pb78) %>% 
  filter(HAULOUT_NUMBER > 0.1) %>% 
  group_by(REF) %>% 
  mutate(Year = first(year(S_DATE))) %>% 
  mutate(id = REF) %>% 
  dplyr::select(id, REF, PTT, S_DATE, E_DATE, HAULOUT_NUMBER, lat, lon, Year) %>% 
  drop_na(lat) %>% 
  mutate(HaulTime = E_DATE - S_DATE) %>% 
  mutate(HaulMins = as.numeric(HaulTime)/60) %>% 
  
  mutate(HaulHours = as.numeric(HaulTime)/60/60) %>% 
  filter(HaulMins > 10) %>% ### Add this in
  filter(HaulHours < 32) %>% ### Add this in
  ungroup() %>% 
  group_by(id) %>%
  arrange(id, S_DATE) %>% 
  mutate(InterHaulOutTime = difftime(S_DATE, lag(E_DATE), units = "hours")) %>%
  mutate(InterHaulOutTime_Hours = as.numeric(InterHaulOutTime)) %>% 
  left_join(MetaData) %>% 
  drop_na(TagDays) %>% 
  ungroup() %>% 
  group_by(id) %>%
  arrange(id, S_DATE) %>% 
  mutate(InterHaulOutTime = difftime(S_DATE, lag(E_DATE), units = "hours")) %>%
  mutate(InterHaulOutTime_Hours = as.numeric(InterHaulOutTime)) 

hist(AllHauloutData$InterHaulOutTime_Hours)



unique(AllHauloutData$Year)


HaulForMerge <-  AllHauloutData %>% 
  dplyr::select(id, REF, PTT, S_DATE, E_DATE, HAULOUT_NUMBER, Year, 13:21) 

# Merge haul-out events within 5 minutes
AllHauloutData_UPDATE <- AllHauloutData %>%
  ungroup() %>% 
  group_by(id) %>%
  arrange(id, S_DATE) %>%
  mutate(group = cumsum(coalesce(InterHaulOutTime_Hours, Inf) > 0.1666667)) %>% 
  group_by(id, group) %>%
  summarise(S_DATE = min(S_DATE),
            E_DATE = max(E_DATE),
            lon = mean(lon),
            lat = mean(lat),
            Weight = unique(`Weight (kg)`),
            Length = unique(`Length (cm)`),
            TagDays = unique(TagDays),
            Sex = unique(Sex),
            Col_lat = unique(Col_lat),
            Col_lon = unique(Col_lon),
            Capture_site = unique(Capture_site)) %>%
  arrange(id, S_DATE) %>%
  group_by(id) %>%
  mutate(Year = first(year(S_DATE))) %>% 
  mutate(InterHaulOutTime = difftime(S_DATE, lag(E_DATE), units = "hours"),
         InterHaulOutTime_Hours = as.numeric(InterHaulOutTime)) %>% 
  mutate(HaulTime = E_DATE - S_DATE) %>% 
  mutate(HaulMins = as.numeric(HaulTime)/60) %>% 
  mutate(HaulHours = as.numeric(HaulTime)/60/60) %>% 
  dplyr:: select(-group) %>% 
  filter(HaulHours < 32) ### Add this in
# left_join(HaulForMerge)




AllHauloutData <- AllHauloutData_UPDATE



ggplot() + geom_point(data = AllHauloutData, aes(x = S_DATE, y = HaulHours)) +
  facet_wrap(~ id, scales = "free_x")


duplicates_df <- AllHauloutData %>%
  group_by(id, S_DATE) %>%
  filter(n() > 1) %>%
  ungroup()


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
# Some visualising of the raw data

# setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Harbour seal movement/Haulouts/ID Checks/")


sf_raw <- QCedRAWDATA %>%                       
  # filter(Capture_site == "Jomfruland") %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

sf_frame <- AllHauloutData %>%                  
  # filter(Capture_site == "Jomfruland") %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  group_by(id) %>% 
  mutate(Num = row_number()) %>% 
  ungroup() %>% 
  mutate(
    tooltip = paste0(
      "<strong>ID:</strong> ", id,
      "<br><strong>S_DATE:</strong> ", format(S_DATE, "%Y‑%m‑%d %H:%M")
    )
  )

removals <- list(
  "pv74-F46_Olivia-20" = c(
    34, 6, 56, 29, 36, 5, 150, 136, 32, 81, 4, 44, 21, 20, 47, 37, 38, 22, 48,
    58, 39, 31, 3, 148, 149, 70, 83, 146, 145, 2, 189, 190, 116, 195, 90, 91, 92,
    196, 197, 208, 184, 188, 186, 185, 152, 84, 102, 183, 40, 19, 41, 18, 49, 42,
    17, 27, 140, 153, 212, 213, 206, 214, 200, 175, 194, 192, 215, 202, 159, 158,
    165, 167, 160, 157, 156, 168, 128:139, 95, 96, 11, 113, 114, 30, 122, 155, 170,
    171, 161, 205, 193, 28, 23, 141, 142, 93, 86, 15, 177:179, 110, 85, 78, 181
  ),
  "pv74-M70_Osito-20" = c(
    141, 142, 143, 134, 32, 31, 28, 27, 5, 6, 145, 133, 24, 23, 38, 39, 40, 41,
    42, 43, 44, 36, 50, 51, 35, 45, 46, 47, 48, 49, 52, 53, 54, 55, 56, 57, 58,
    59, 64, 65, 66, 67, 75, 76, 77, 79, 80, 81, 82, 86, 130, 136, 33, 29, 30, 2,3,4, 135, 131, 60:63
  ),
  "pv74-M86_Bjorn-20" = c(
    82, 140, 141, 102, 105, 106, 144, 137, 146, 216, 217, 195, 196, 193, 242,
    234, 233, 239, 238,
    51:53, 59:80, 87:100, 148:168, 190:192, 208:211, 173:179, 244:252,
    257, 255,254,253, 107,55, 54, 147, 22, 50, 136, 108, 49, 22, 200:202, 212:214, 42, 134, 194, 219:221, 222, 34, 171, 172,
    27, 28, 45, 43,135, 19, 41, 188, 186, 181, 240, 241,
    
    43, 77, 78, 106, 26, 105, 42, 49, 50, 41, 40, 76, 80, 105, 26, 24, 23, 22, 29, 21, 39, 36, 19, 20,
    91, 92, 98, 114
    
  ),
  "pv74-M62_Gamle-Erik-20" = c(
    13, 14,40, 41, 56, 2, 4, 29, 54, 41, 12
  ),
  "pv74-M88_Diego-20" = c(
    65, 70, 71, 4, 5, 50, 80, 79, 170, 172, 85, 77, 76, 169,
    17:20, 99:102, 172:180, 181:187, 22:36, 2, 3, 72,189, 78, 194, 197, 167, 78, 86:88, 166, 199, 163, 161, 90, 91,
    171, 21, 196, 7, 190, 188
  )
)


sf_frame_cleaned <- sf_frame %>%
  filter(!purrr::map2_lgl(id, Num, function(x, y) y %in% removals[[x]] %||% FALSE))



# unique IDs 
ids <- intersect(sf_raw$id, sf_frame_cleaned$id)

gps_col    <- "#6A0DAD"   
haulout_col<- "#FF7F00"   

sf_frame_cleaned <- sf_frame_cleaned %>%
  mutate(
    tooltip = paste0(
      "<strong>id: </strong>", id, "<br>",
      "<strong>Num: </strong>", Num, "<br>",
      "<strong>S_DATE: </strong>", format(S_DATE, "%Y-%m-%d %H:%M")
    ),
    popup = tooltip  
  )


# # loop over each id
# for (this_id in ids) {
#   
#   gps_pts   <- sf_raw   %>% filter(id == this_id)
#   haul_pts  <- sf_frame_cleaned %>% filter(id == this_id)
#   
#   m <- leaflet() %>%
#     addTiles() %>%
#     
#     ## GPS layer  (purple)
#     addCircleMarkers(
#       data        = gps_pts,
#       lng         = ~lon, lat = ~lat,
#       radius      = 4,
#       color       = gps_col, fillColor = gps_col,
#       stroke      = TRUE, fillOpacity = 0.8,
#       group       = "GPS"
#     ) %>%
#     
#     ## Haul‑out layer (orange) 
# addCircleMarkers(
#   data        = haul_pts,
#   lng         = ~lon, lat = ~lat,
#   radius      = 5,
#   color       = haulout_col, fillColor = haulout_col,
#   stroke      = TRUE, fillOpacity = 0.8,
#   label       = ~lapply(tooltip, htmltools::HTML),   
#   popup       = ~htmltools::HTML(tooltip),
#   group       = "Haul‑out"
# )
#     
#     # ## 
#     # addLegend(
#     #   position = "bottomright",
#     #   # colors   = c(gps_col, haulout_col),
#     #   # labels   = c("GPS", "Haul‑out"),
#     #   title    = NULL
#     # )
#   
#     
#   ## 
#   file_name <- paste0("map_id_", this_id, ".html")
#   saveWidget(m, file = file_name, selfcontained = TRUE)
#   message("saved ", file_name)
# }


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
# Let's fix the 2020 data to capture the points that are clearly haul-outs



sf_raw_2020 <- QCedRAWDATA %>%                       
  filter(TaggingYear == "2020") %>% 
  mutate(Type = "GPS",
         Year= TaggingYear) %>% 
  mutate(`Weight (kg)` = "Weight") %>% 
  mutate(S_DATE = date)

sf_frame_2020 <- sf_frame_cleaned %>% 
  filter(Year == "2020") %>% 
  # filter(Capture_site == "Jomfruland")  %>%
  mutate(Type = "Haul") %>% 
  full_join(sf_raw_2020) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  group_by(id) %>%
  mutate(Num = row_number()) %>% 
  ungroup() %>% 
  mutate(
    tooltip = paste0(
      "<strong>id: </strong>", id, "<br>",
      "<strong>Num: </strong>", Num, "<br>",
      "<strong>S_DATE: </strong>", format(S_DATE, "%Y-%m-%d %H:%M")))


ids <- unique(sf_frame_2020$id)

pal_type <- colorFactor(palette = c("blue", "darkorange"), domain = c("GPS", "Haul"))

changes <- list(
  "pv74-F46_Olivia-20" = c(411, 407, 395, 386, 397, 406, 394, 410, 391, 450:453, 455:470),
  
  "pv74-M70_Osito-20" = c(416, 420, 314, 319, 322, 219, 220, 204, 218, 235, 232, 234,
                          212, 237, 265, 266, 272, 269, 270, 267, 268, 271, 295, 277, 298,
                          301, 284, 291, 287, 278, 293, 300, 354, 353, 346, 350, 355, 348, 
                          289, 280, 107, 99, 112, 108, 104, 105, 109, 103, 279, 352, 221, 
                          190, 233, 215, 207, 213, 229, 214, 230, 199, 201, 196, 209, 190, 230),
  
  "pv74-M62_Gamle-Erik-20" = c(306, 305, 302, 301, 303, 291, 292, 304, 204, 298, 277,
                               297, 318, 274, 280, 276, 281, 322, 271,
                               242, 241, 238, 237, 235, 231, 232, 233, 
                               236, 243, 230, 176, 234, 229, 177, 230, 
                               196, 290, 240, 432, 196, 436, 420, 345, 
                               439, 438, 338, 250, 188, 190, 429, 356, 357, 
                               206, 295, 299,294, 297, 317, 278,
                               200, 217, 218, 205, 208, 210, 209, 210, 201, 361)
)


changes1 <- list(
   "pv74-F46_Olivia-20" = c(14, 24, 16, 17, 100, 77, 44, 54, 35, 76, 78, 91, 1, 75, 99, 65, 93, 47, 11, 90, 89, 49, 48, 2, 46, 15, 6, 34, 88, 79, 10, 13, 45, 8, 69 ,66,
                            67, 68, 85:87, 95:97, 43, 36, 9, 12, 7, 62, 74, 101),
  "pv74-M86_Bjorn-20" = c(35:39, 68, 69, 34, 22, 21, 69, 95, 96, 20, 97, 70, 66, 40, 41, 33, 32, 31, 103, 104, 30, 23, 19, 94, 27:30, 89, 24, 
                          67, 80:88, 1:15, 98, 99, 100, 101, 102))




# sf_frame_2020 <- sf_frame_2020 %>%
#   rowwise() %>%
#   mutate(Type = if (id %in% names(changes) && Num %in% changes[[id]]) "Haul" else Type) %>%
#   ungroup()

sf_frame_2020 <- sf_frame_2020 %>%
  rowwise() %>%
  mutate(Type = case_when(
    id %in% names(changes) && Num %in% changes[[id]] ~ "Haul",
     id %in% names(changes1) && Num %in% changes1[[id]] ~ "GPS",
    
    TRUE ~ Type
  )) %>%
  ungroup()


# Loop through each ID 
for (this_id in ids) {
  
  sf_id <- sf_frame_2020 %>% filter(id == this_id)
  
  m <- leaflet(data = sf_id) %>%
    addTiles() %>%
    addCircleMarkers(
      lng         = ~lon, lat = ~lat,
      radius      = 4,
      color       = ~pal_type(Type),
      stroke      = TRUE, fillOpacity = 0.8,
      label       = ~lapply(tooltip, htmltools::HTML),  
      popup       = ~htmltools::HTML(tooltip)           
    ) %>%
    addLegend(
      position = "bottomright",
      pal = pal_type,
      values = ~Type,
      title = "Location Type",
      opacity = 1
    )
  
  # Save HTML 
  file_name <- paste0("FIXED_", this_id, ".html")
  saveWidget(m, file = file_name, selfcontained = TRUE)
  message("Saved: ", file_name)
}



HAUL_sf_frame_2020 <- sf_frame_2020 %>% 
  st_drop_geometry() %>% 
  filter(Type == "Haul")



##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

# MAPPING

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

# Base <- basemap(limits = c(6.9, 12.31, 56.19, 59.8), bathymetry = FALSE, projection.grid = FALSE) +
#   theme_minimal() + 
#   guides(fill = NULL) +
#   theme(
#     legend.position = "none",
#     panel.grid       = element_blank(),    
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border     = element_blank()
#   ) 
# 
# 
# Base


##############################################################################################################################################
# Norway
##############################################################################################################################################

# setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Shapefiles/Norway/")

setwd("C:/Users/AADDT/Documents/Javed/Haulouts/Haulouts/Data files/Terrestrial and Combined Maritime Shapefiles/")


NorwayLandMap_sf <- read_sf("Fylker19.geojson") ## Seal management units

NorwayLandMap_sf <- NorwayLandMap_sf %>%
  st_transform(st_crs(4326))  

NorwayLandMap <- ggplot() +
  geom_sf(data = NorwayLandMap_sf, aes(fill = fylkesnavn) ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    # legend.justification = c(0.1,0.1),
    #     legend.position = c(0.45,0.1),
    # legend.key.size = unit(1.5, "lines"),
    legend.background = element_rect(fill = "lightgrey", colour = "black"),
    legend.key = element_rect(fill = "lightgrey")) +
  labs(fill = "Harbour seal management units\nin Norway") +
  guides(fill=guide_legend(ncol=2))
# NorwayLandMap

# tiff("NorwayLandMap.tiff", width = 6, height = 6, units = 'in', res = 300)
# NorwayLandMap
# dev.off()
# 

NorwayLandMap_crop <- st_crop(NorwayLandMap_sf
                              , xmin = 6.9, xmax = 12.31, ymin = 56.19, ymax = 59.8)

region_labels <- c(
  "Rogaland" = "", "Telemark" = "TL", "Vestfold" = "VF", "Østfold" = "OF", 
  "Oslo" = "", "Vest-Agder" = "VA", "Hordaland" = "", "Buskerud" = "BS", 
  "Aust-Agder" = "AA", "Akershus" = "AK"
)

NorwayLandMap_crop <- NorwayLandMap_crop %>%
  mutate(abbreviation = recode(fylkesnavn, !!!region_labels))  


##############################################################################################################################################
# Norway Martime Extension
##############################################################################################################################################

# setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Shapefiles/Norway/")
setwd("C:/Users/AADDT/Documents/Javed/Haulouts/Haulouts/Data files/Terrestrial and Combined Maritime Shapefiles/")


st_layers("Basisdata_0000_Norge_25833_Fylker_GeoJSON.geojson")

NorwayMaritimeMap_sf <- read_sf("Basisdata_0000_Norge_25833_Fylker_GeoJSON.geojson", layer = "administrative_enheter.fylke")


NorwayMaritimeMap_sf <- NorwayMaritimeMap_sf %>%
  st_transform(st_crs(4326))  

NorwayMartimeMap_crop <- st_crop(NorwayMaritimeMap_sf
                                 , xmin = 6.9, xmax = 12.31, ymin = 56.19, ymax = 59.8)

NorwayMaritimeMap <- ggplot() +
  geom_sf(data = NorwayMartimeMap_crop, aes(fill = objtype), fill = "lightblue", colour= "black" ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    # legend.justification = c(0.1,0.1),
    #     legend.position = c(0.45,0.1),
    # legend.key.size = unit(1.5, "lines"),
    legend.background = element_rect(fill = "lightgrey", colour = "black"),
    legend.key = element_rect(fill = "lightgrey")) +
  labs(fill = "Harbour seal management units\nin Norway") +
  guides(fill=guide_legend(ncol=2))
# NorwayMaritimeMap



##############################################################################################################################################
# Norway Inset
##############################################################################################################################################



InsetMap <- basemap(limits = c(5.5, 23, 53, 72), bathymetry = FALSE , rotate = TRUE, grid.col = NA, land.col = "grey90") +
  theme(legend.position = "right")  +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  # ggspatial::geom_spatial_rect(aes(xmin = 5.4, ymin = 56, xmax = 12.5, ymax = 60),colour = "red", fill = "transparent", linewidth = 0.75) +
  theme_bw() +
  theme_void() +  # Remove all extra spacing
  theme(
    panel.border = element_rect(fill = NA, colour = "black", size = 1.3),  
    plot.margin = margin(0, 0, 0, 0)  # Trim excess space
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white")  
  ) 

# InsetMap



# colors_npg <- pal_npg("nrc")(10)
# colors_npg
# "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF" "#7E6148FF" "#B09C85FF"
# colors_jco <- pal_jco("default")(8)
# colors_jco
# "#0073C2FF" "#EFC000FF" "#868686FF" "#CD534CFF" "#7AA6DCFF" "#003C67FF" "#8F7700FF" "#3B3B3BFF"
# 
# 
# scale_color_manual(values = c( "#5C126EFF","#C43C4EFF", "#000004FF", "orange")) +



# 
# "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "red4" "#7E6148FF",  "orange", "#0073C2FF" "#EFC000FF"  "#5C126EFF", "#C43C4EFF" ,"#7AA6DCFF" ,"#8F7700FF" ,"#3B3B3BFF", "black"




fylkes_levels <- c(
  "Rogaland", "Finnmark", "Troms", "Møre og Romsdal", "Telemark",
  "Vestfold", "Østfold", "Oslo", "Oppland", "Hedmark", "Nordland",
  "Vest-Agder", "Trøndelag", "Sogn og Fjordane", "Hordaland",
  "Buskerud", "Aust-Agder", "Akershus"
)

fylkes_colors_raw <- c(
  "#F39B7FFF", "#7E6148FF", "#00A087FF", "#3C5488FF", "orange",     
  "#C43C4EFF", "black", "#8491B4FF", "#91D1C2FF", "red4", 
  "#4DBBD5FF", "#3B3B3BFF", "#EFC000FF", "#E64B35FF", "#7AA6DCFF",
  "#8F7700FF", "#5C126EFF", "#0073C2FF"  
)


# Assign names
fylkes_colors_Norway <- setNames(fylkes_colors_raw, fylkes_levels)



InsetMap1 <- InsetMap +  
  geom_sf(data = NorwayLandMap_sf, aes(fill = fylkesnavn)) +
  scale_fill_manual(values = fylkes_colors_Norway) +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.background = element_rect(fill = "lightgrey", colour = "black"),
    legend.key = element_rect(fill = "lightgrey"),
    
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(t = 0, b = -0.8 , l = -0.8, 0)) +
  labs(fill = "Harbour seal management units\nin Norway") +
  guides(fill = guide_legend(ncol = 2)) +
  ggspatial::geom_spatial_rect(
    aes(xmin = 5, ymin = 56, xmax = 13, ymax = 60),
    colour = "red", fill = "transparent", linewidth = 0.6)

# InsetMap1

# pdf("InsetMap.pdf", width = 5, height = 5.5)
# InsetMap
# dev.off()




##############################################################################################################################################
# Read in Combined Maritime
##############################################################################################################################################

# setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Shapefiles/")
setwd("C:/Users/AADDT/Documents/Javed/Haulouts/Haulouts/Data files/Terrestrial and Combined Maritime Shapefiles/")

Combined_Maritime_Boundaries <- read_sf("Combined_Maritime_Boundaries.shp")

Combined_Maritime_Boundaries_crop <- st_crop(Combined_Maritime_Boundaries
                                             , xmin = 6.9, xmax = 12.31, ymin = 56.19, ymax = 59.8)



# setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Shapefiles/")

Swedend_Denmark_Combined_Maritime_Boundaries <- read_sf("Swedend_Denmark_Combined_Maritime_Boundaries.shp")

Swedend_Denmark_Combined_Maritime_Boundaries_crop <- st_crop(Swedend_Denmark_Combined_Maritime_Boundaries
                                                             , xmin = 6.9, xmax = 12.31, ymin = 56.19, ymax = 59.8)
# ggplot() +
#   geom_sf(data = Swedend_Denmark_Combined_Maritime_Boundaries_crop, fill = "lightblue", color = "black") +
#   theme_minimal() +
#   labs(title = "Cleaned Merged Boundaries")
# 




##############################################################################################################################################
# HIGH RES BATHYMETRY
##############################################################################################################################################

# setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Shapefiles/")

SkagerrakBathymetry <- read_rds("Skagerrak_Bathymetry_EMODnet_115m.rds")


SkagerrakBathymetry_NA <- SkagerrakBathymetry %>% 
  filter(is.na(Depth)|
           Depth >= 0) %>% 
  mutate(Land = "Land")

# ggplot() +
#   geom_tile(data = SkagerrakBathymetry_NA, aes(x = x, y = y), fill = "grey")

unique(NorwayLandMap_crop$fylkesnavn)

##############################################################################################################################################
# CREATE BASEMAP
##############################################################################################################################################


label_coords <- st_centroid(NorwayLandMap_crop) %>%
  st_coordinates() %>%
  as.data.frame()

NorwayLandMap_crop$X <- label_coords$X
NorwayLandMap_crop$Y <- label_coords$Y

NorwayLandMap_crop <- NorwayLandMap_crop %>%
  mutate(
    Y = case_when(
      abbreviation %in% c("AA") ~ Y - 0.4,
      TRUE ~ Y
    ),
    X = case_when(
      abbreviation %in% c("VA", "AA") ~ X + 0.3,
      TRUE ~ X
    )
  ) %>% 
  mutate(
    Y = case_when(
      abbreviation %in% c("VA") ~ Y - 0.1,
      TRUE ~ Y
    )) %>% 
  mutate(
    Y = case_when(
      abbreviation %in% c("AK") ~ Y + 0.4,
      TRUE ~ Y
    ))



# setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Harbour seal movement/Haulouts/")
setwd("C:/Users/AADDT/Documents/Javed/Haulouts/Haulouts/Data files/Clean Figures/")


# Updated BASEMAP plot
BASEMAP <- ggplot() +
  geom_tile(data = SkagerrakBathymetry_NA, aes(x = x, y = y), fill = "#eeeac4") +
  
  ggnewscale::new_scale_fill() +
  
  geom_sf(data = Combined_Maritime_Boundaries_crop, aes(fill = Area), fill = "lightblue", colour = "black", alpha = 0.3) +
  # geom_sf(data = NorwayMartimeMap_crop, aes(fill = objtype), fill = "lightblue", colour= "black", alpha = 0.3) +
  
  guides(fill = FALSE) +
  
  ggnewscale::new_scale_fill() +
  
  geom_sf(data = NorwayLandMap_crop, fill = "#eeeac4", color = "black", size = 0.3) +
  
  # Adjusted text labels
  geom_text(data = NorwayLandMap_crop, aes(x = X, y = Y, label = abbreviation),
            size = 3, fontface = "bold") +
  
  labs(x = "Longitude", y = "Latitude", fill = "Harbour seal\nmanagement\nunits") +
  coord_sf(
    ylim = c(st_bbox(Combined_Maritime_Boundaries_crop)["ymin"] + 0.1, 
             st_bbox(Combined_Maritime_Boundaries_crop)["ymax"] - 0.1),
    expand = FALSE
  ) +
  theme_minimal()

 # BASEMAP



# tiff("BASEMAP.tiff", width = 5, height = 5.5, units = 'in', res = 300)
# BASEMAP
# dev.off()

# ggsave("myBASEMAP.png", width = 5, height = 5.5, units = "in")

# write_rds(BASEMAP, "BASEMAP.rds")


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
# Plot the haul out data

# setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Harbour seal movement/Haulouts/")

# MetaData1 <- MetaData1 %>% 
#   group_by(Capture_site) %>% 
#   mutate(Plot_col_lon = if_else(Capture_site == "Hvaler", Col_lon + 0.3, Col_lon - 0.3)) 
#   
levels(AllHauloutData$Capture_site)


sf_frame_cleaned <- sf_frame_cleaned %>% 
  filter(!Year == "2020") %>% 
  st_drop_geometry() 


sf_frame_cleaned <- sf_frame_cleaned %>% 
  full_join(HAUL_sf_frame_2020)


AllHauloutData <- sf_frame_cleaned 

MetaData1 <- MetaData1 %>%
  mutate(Capture_site = factor(Capture_site, levels = c("Askerøy", "Jomfruland", "Bolærne", "Hvaler")))

AllHauloutData <- AllHauloutData %>%
  mutate(Capture_site = factor(Capture_site, levels = c("Askerøy", "Jomfruland", "Bolærne", "Hvaler")))

QCedRAWDATA <- QCedRAWDATA %>%
  mutate(Capture_site = factor(Capture_site, levels = c("Askerøy", "Jomfruland", "Bolærne", "Hvaler")))

site_colors <- c("Askerøy" = "#5C126EFF",
                 "Jomfruland" = "orange",
                 "Bolærne" = "#C43C4EFF",
                 "Hvaler" = "#000004FF")

# arrow_dx <- 0.3  
# arrow_dy <- 0.3  
# 
# MetaData1_arrows <- MetaData1 %>%
#   mutate(
#     xend = case_when(
#       Capture_site == "Hvaler"     ~ Col_lon - arrow_dx,  # left
#       TRUE                         ~ Col_lon + arrow_dx   # right
#     ),
#     yend = Col_lat - arrow_dy     
#   )

Fig1 <- BASEMAP +
  geom_spatial_path(data = QCedRAWDATA,
                    aes(x = lon, y = lat, colour = Capture_site, group = id),
                    alpha = 0.5, size = 0.3) +
  guides(fill = FALSE, color = guide_legend(override.aes = list(fill = NA))) +
  # geom_spatial_segment(
  #   data = MetaData1_arrows,
  #   aes(x = Col_lon, y = Col_lat, xend = xend, yend = yend, colour = Capture_site),
  #   arrow = arrow(type = "closed", length = unit(0.2, "cm")),
  #   size = 1.2
  # ) +
  scale_color_manual(values = site_colors) +
  guides(colour = FALSE) +
  ggnewscale::new_scale_color() +
  
  geom_spatial_point(data = filter(AllHauloutData, Capture_site == "Bolærne"),
                     aes(x = lon, y = lat, colour = Capture_site),
                     alpha = 0.5, size = 3.5) +
  geom_spatial_point(data = filter(AllHauloutData, Capture_site == "Askerøy"),
                     aes(x = lon, y = lat, colour = Capture_site),
                     alpha = 0.5, size = 3.5) +
  geom_spatial_point(data = filter(AllHauloutData, Capture_site == "Hvaler"),
                     aes(x = lon, y = lat, colour = Capture_site),
                     alpha = 0.5, size = 3.5) +
  geom_spatial_point(data = filter(AllHauloutData, Capture_site == "Jomfruland"),
                     aes(x = lon, y = lat, colour = Capture_site),
                     alpha = 0.5, size = 3) +
  scale_color_manual(name = "Tagging and\nhaul-out\nlocations",
                     values = site_colors,
                     breaks = c("Askerøy", "Jomfruland", "Bolærne", "Hvaler")) +
  ggspatial::annotation_scale(location = "tr", text_cex = 0.5) +
    labs(x = "Longitude", y = "Latitude") +
  theme_bw(10) +
  theme(
    strip.background = element_rect(fill = "gray85"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black"),
    # legend.key = element_rect(fill = "lightgrey", colour = NA),
    legend.justification = c(0.05, 0.05),
    legend.position = c(0.02, 0.69),
    legend.background = element_rect(fill = "lightgrey", colour = "black"),
    legend.key = element_rect(fill = "lightgrey")
  )

# Fig1

InsetMap_grob <- ggplotGrob(InsetMap1)


Fig1_Inset <- Fig1 + annotation_custom2(
  InsetMap_grob, 
  data = MetaData1,
  
  xmin = -Inf, xmax = 8.15, 
  ymin = -Inf, ymax = 57.4)


## Save the figure and add configure arrows and scalebar text manually  in Inkscape



 tiff("Fig1.tiff", width = 5, height = 5.5, units = 'in', res = 300)
 Fig1_Inset
 dev.off()

# pdf("Fig1_Unedited.pdf", width = 5, height = 5.5)
# Fig1_Inset
# dev.off()




##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
# Track cleaning and pre-processing

AllHauloutData <- AllHauloutData %>% 
  # mutate(Sex=replace_na(Sex, "M")) %>% 
  # mutate(Col_lon=replace_na(Col_lon, 9.571337)) %>% 
  # mutate(Col_lat=replace_na(Col_lat, 58.86441)) %>% 
  mutate(dstep=distHaversine(cbind(Col_lon, Col_lat),cbind(lon,lat))/1000)

unique(AllHauloutData$id)
range(AllHauloutData$TagDays)
range(AllHauloutData$dstep)

mean(AllHauloutData$dstep)
sd(AllHauloutData$dstep)


##############################################################################################################################################
### Number of haulouts summary
##############################################################################################################################################


IDSummary_HauloutsEvents <- AllHauloutData %>% 
  mutate(count = 1) %>% 
  group_by(id) %>% 
  dplyr::select(id, count) %>% 
  summarise(sumevents = sum(count)) %>% 
  distinct() 

HauloutsEvents_summary <- arith_mean_ci(IDSummary_HauloutsEvents$sumevents)
HauloutsEvents_summary
cat("Population-level mean of max number of haulouts:", round(HauloutsEvents_summary$mean, 2), "\n")
cat("95% CI:", round(HauloutsEvents_summary$lower_CI, 2), "-", round(HauloutsEvents_summary$upper_CI, 2), "\n")


##############################################################################################################################################
### Tag duration (days) summary
##############################################################################################################################################

IDSummary_TagDays <- AllHauloutData %>% 
  group_by(id) %>% 
  dplyr::select(id, TagDays) %>% 
  summarise(mean_tagdays = mean(TagDays, na.rm = TRUE)) %>%
  distinct()

tagdays_summary <- arith_mean_ci(IDSummary_TagDays$mean_tagdays)
cat("Population-level mean of individual tays days:", round(tagdays_summary$mean, 2), "\n")
cat("95% CI:", round(tagdays_summary$lower_CI, 2), "-", round(tagdays_summary$upper_CI, 2), "\n")


##############################################################################################################################################
### Distance from colony summary
##############################################################################################################################################

IDSummary_Dstep <- AllHauloutData %>%
  group_by(id) %>%
  summarise(mean_dstep = mean(dstep, na.rm = TRUE)) %>%
  ungroup()


dstep_summary <- arith_mean_ci(IDSummary_Dstep$mean_dstep)
dstep_summary
# cat("Population-level mean of individual dstep means:", round(dstep_summary$mean, 2), "\n")
# cat("95% CI:", round(dstep_summary$lower_CI, 2), "-", round(dstep_summary$upper_CI, 2), "\n")


##############################################################################################################################################
### Duration of haulouts
##############################################################################################################################################

# mean(AllHauloutData$HaulHours)

IDSummary_duration <- AllHauloutData %>%
  mutate(dummy = 1) %>% 
  filter(!Year == "2020") %>% 
  group_by(id) %>%
  # summarise(mean_duration = mean(HaulHours, na.rm = TRUE)) %>%
  ungroup()

range(IDSummary_duration$HaulHours)

FigS1A <- ggplot() + geom_histogram(data = AllHauloutData %>% filter(!Year == "2020"), aes(x = HaulHours, fill = Capture_site), position = "dodge", colour = "black", bins = 10)  +
  scale_fill_manual(values = site_colors) +
  theme_bw(10) +
  theme(
    strip.background = element_rect(fill="gray85"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour="black"),
    legend.key = element_rect(fill="lightgrey", colour=NA)) +
  theme(legend.justification = c(0.4,0.9),
        legend.position = c(0.8,0.95),
        # legend.key.size = unit(1.5, "lines"),
        legend.background = element_rect(fill = "lightgrey", colour = "black"),
        legend.key = element_rect(fill = "lightgrey")) +
  scale_x_continuous(breaks = c(0,5,10,15,20,25, 30)) +
  scale_y_continuous(expand = c(0,10), breaks = seq(0,600, by= 100)) +
  theme(legend.position = "none") +
  labs(fill = "Tagging locations", x = "Haul-out duration (hours)", y = "Frequency") 
FigS1A



mean_duration_summary <- arith_mean_ci(IDSummary_duration$HaulHours)
mean_duration_summary
# cat("Population-level mean of individual duration means:", round(mean_duration_summary$mean), "\n")
# cat("95% CI:", round(mean_duration_summary$lower_CI, 2), "-", round(mean_duration_summary$upper_CI, 2), "\n")
# 



##############################################################################################################################################
### Fig S1
##############################################################################################################################################


FigS1B <- ggplot() + geom_histogram(data = AllHauloutData, aes(x = dstep, fill = Capture_site), position = "dodge", colour = "black", bins = 10)  +
  scale_fill_manual(values = site_colors) +
  theme_bw(10) +
  theme(
    strip.background = element_rect(fill="gray85"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour="black"),
    legend.key = element_rect(fill="lightgrey", colour=NA)) +
  theme(legend.justification = c(0.4,0.9),
        legend.position = c(0.8,0.95),
        # legend.key.size = unit(1.5, "lines"),
        legend.background = element_rect(fill = "lightgrey", colour = "black"),
        legend.key = element_rect(fill = "lightgrey")) +
  # geom_hline(yintercept = 0, colour = "black") +
  scale_x_continuous(breaks = c(0,50,100,150,200,250)) +
  scale_y_continuous(expand = c(0,10), breaks = seq(0,1300, by= 100)) +
  # theme(legend.position = "none") +
  # scale_y_continuous(expand = c(0,10), breaks = c(0,100,200, 300, 400, 500, 600, 700, 800, 900)) +
  labs(fill = "Tagging locations", x = "Haul-out distance from initial tagging location (km)", y = "Frequency") 
FigS1B


FigS1_MS <- FigS1A + FigS1B + plot_layout(ncol = 1)
FigS1_MS


tiff("FigS1_clean.tiff", width = 6, height = 7, units = 'in', res = 300)
FigS1_MS
dev.off()






##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
# Overlap with management areas - data wrangling

management_areas <- NorwayLandMap_crop

management_areas_buffered <- st_buffer(management_areas, dist = 1000)  # 1000 meters = 1 km


# seals_sf <- st_transform(seals_sf, st_crs(management_areas))

seals_sf <- st_as_sf(AllHauloutData, coords = c("lon", "lat"), crs = 4326)

# Assign haul-out locations to management areas
seals_sf <- st_join(seals_sf, management_areas_buffered, join = st_within)

EEZwaters <- Combined_Maritime_Boundaries_crop %>% 
  dplyr::select(Area, geometry)


seals_sf <- st_join(seals_sf, EEZwaters, join = st_within)

seals_sf <- seals_sf %>%
  mutate(Area = ifelse(!is.na(fylkesnavn) & is.na(Area), "Norway", Area)) %>% 
  mutate(fylkesnavn = ifelse(is.na(fylkesnavn) & !is.na(Area), Area, fylkesnavn))
# mutate(Area = ifelse(!is.na(fylkesnavn) & is.na(Area), "Norway", Area))

AssignNA <- seals_sf %>% 
  filter(is.na(fylkesnavn))

seals_sf <- seals_sf %>%
  mutate(
    fylkesnavn = ifelse(is.na(fylkesnavn) & st_coordinates(seals_sf)[,2] > 58.0, "Sweden",
                        ifelse(is.na(fylkesnavn) & st_coordinates(seals_sf)[,2] <= 58.0, "Denmark", fylkesnavn)),
    Area = ifelse(is.na(fylkesnavn) & st_coordinates(seals_sf)[,2] > 58.0, "Sweden",
                  ifelse(is.na(fylkesnavn) & st_coordinates(seals_sf)[,2] <= 58.0, "Denmark", Area))
  ) %>% 
  mutate(
    Area = ifelse(!is.na(fylkesnavn) & is.na(Area), fylkesnavn, Area))
# 
# # Manually correct fylkesnavn for Olivia's specific haulout
# fylkesnavn = ifelse(id == "pv74-F46_Olivia-20" & 
#                       format(S_DATE, "%Y-%m-%d %H:%M:%S") == "2020-11-23 16:36:37",
#                     "Telemark", fylkesnavn)



Test <- seals_sf %>% 
  filter(id == "pv74-F46_Olivia-20")



unique(seals_sf$id)

# Quick summary
Summaryseals_sf <- seals_sf %>%
  group_by(id, Area) %>%
  summarise(count = n(), .groups = "drop")



# Extract the first recorded management area and location for each seal
tagging_locations <- seals_sf %>%
  group_by(id) %>%
  slice_min(order_by = S_DATE) %>% 
  mutate(deployment_management_region = Capture_site) %>%
  mutate(deployment_management_region = case_when(
    trimws(Capture_site) == "Bolærne"  ~ "Vestfold",
    trimws(Capture_site) == "Jomfruland" ~ "Telemark",
    trimws(Capture_site) == "Hvaler"    ~ "Østfold",
    trimws(Capture_site) == "Askerøy"    ~ "Aust-Agder",
    TRUE ~ NA_character_)) %>% 
  #   TRUE ~ NA_character_  # Set NA for any other sites
  # rename(tagging_long = longitude, tagging_lat = latitude) %>%
  ungroup() %>% 
  dplyr::mutate(lat = sf::st_coordinates(.)[,2],
                lon = sf::st_coordinates(.)[,1]) %>% 
  as.data.frame() %>% 
  dplyr::select(id, deployment_management_region, Col_lon, Col_lat) 


HauloutModelFrame <- left_join(seals_sf, tagging_locations) %>% 
  dplyr::mutate(lat = sf::st_coordinates(.)[,2],
                lon = sf::st_coordinates(.)[,1]) 



# Compute distance for each haul-out event
HauloutModelFrame <- HauloutModelFrame %>%
  mutate(distance_km = mapply(calculate_distance, Col_lon, Col_lat, lon, lat))



HauloutModelFrame <- HauloutModelFrame %>%
  mutate(Cross_Boundary_Management = ifelse(fylkesnavn != deployment_management_region, 1, 0)) %>% 
  as.data.frame() 




## Summary Table

SummaryTable <- HauloutModelFrame %>% 
  dplyr::select(id, Year, Weight, Sex, TagDays, deployment_management_region, S_DATE) %>%
  mutate(Tagging_date = as.Date(S_DATE)) %>%
  group_by(id) %>% 
  filter(Tagging_date <= min(Tagging_date)) %>% 
  dplyr::distinct() %>% 
  ungroup() %>% 
  group_by(deployment_management_region) %>%
  mutate(Seal_ID = paste(deployment_management_region, dense_rank(id), sep = "_"),
         "Tag_duration (days)" = TagDays) %>% 
  dplyr::select(Seal_ID, Tagging_date, Weight, Sex, TagDays, `Tag_duration (days)`) %>% 
  arrange(Seal_ID) %>% 
  ungroup(deployment_management_region) %>% 
  dplyr::select(Seal_ID, Weight, Sex, Tagging_date, `Tag_duration (days)`) %>% 
  distinct()

# write_csv2(SummaryTable, "SummaryTable.csv")




##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
# Overlap with management areas - plotting and calculations


# Quick summary of cross boundary ventures
Summaryseals_CB <- HauloutModelFrame %>%
  group_by(id, deployment_management_region, fylkesnavn, Cross_Boundary_Management) %>%
  summarise(count = n(), .groups = "drop") 


OnlyResident <- Summaryseals_CB %>% 
  filter(Cross_Boundary_Management == 0)


sum(OnlyResident$count)
sum(OnlyResident$count)/sum(Summaryseals_CB$count)
unique(OnlyResident$id)


OnlyCrossBoundary <- Summaryseals_CB %>% 
  filter(Cross_Boundary_Management == 1)
sum(OnlyCrossBoundary$count)
sum(OnlyCrossBoundary$count)/sum(Summaryseals_CB$count)
unique(OnlyCrossBoundary$id)

OnlyCrossBoundaryDistance <- HauloutModelFrame %>% 
  filter(Cross_Boundary_Management == 1) ## 12 individuals that crossed boundary

HauloutModelFrame$Month <- lubridate::month(HauloutModelFrame$S_DATE)

# Reorder the Month factor 
HauloutModelFrame$Month1 <- factor(HauloutModelFrame$Month, 
                                   levels = c(8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7), 
                                   labels = c("August", "September", "October", "November", "December", 
                                              "January", "February", "March", "April", "May", "June", "July"))


SummaryMonth <- HauloutModelFrame %>%
  group_by(id, Month) %>%
  summarise(count = n(), .groups = "drop") 


# Create the NewID column 
HauloutModelFrame <- HauloutModelFrame %>%
  group_by(deployment_management_region) %>%
  mutate(NewID = paste(deployment_management_region, dense_rank(id), sep = "_")) %>%
  ungroup()  


NewIDs <- HauloutModelFrame %>% 
  dplyr::select(id, NewID) %>% 
  distinct()


# Prop summary
Summaryseals_CB_proportion <- HauloutModelFrame %>%
  group_by(deployment_management_region, NewID, Cross_Boundary_Management) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) 


NumberPerArea <- Summaryseals_CB %>% 
  mutate(tally = 1) %>% 
  # dplyr::select(- Cross_Boundary_Management, - count) %>% 
  group_by(id) %>% 
  mutate(NZones = sum(tally)) %>% 
  dplyr::select(id, NZones) %>% 
  distinct() %>% 
  left_join(NewIDs)


NumberPerArea_summary <- arith_mean_ci(NumberPerArea$NZones)

cat("Population-level mean of N zones:", round(NumberPerArea_summary$mean, 2), "\n")
cat("95% CI:", round(NumberPerArea_summary$lower_CI, 2), "-", round(NumberPerArea_summary$upper_CI, 2), "\n")



OnlyCrossProps <- Summaryseals_CB_proportion %>% 
  filter(Cross_Boundary_Management == 0) %>% 
  mutate(proportion = (1- proportion) * 100) %>% 
  mutate(
    label = paste0(ceiling(proportion), "%")
  )


# OnlyCrossProps <- bind_rows(
#   OnlyCrossProps,
#   tibble(
#     deployment_management_region = "Østfold",
#     NewID = "Østfold_1",
#     Cross_Boundary_Management = 0,
#     count = 64,
#     proportion = 1,
#     label = "100%"
#   )
# )




OnlyCrossProps$NewID <- factor(OnlyCrossProps$NewID, levels = rev(c(
  "Aust-Agder_1", "Aust-Agder_2", "Aust-Agder_3",  
  "Telemark_1", "Telemark_2", "Telemark_3", "Telemark_4", "Telemark_5", 
  "Telemark_6", "Telemark_7", "Telemark_8", 
  "Vestfold_1", "Vestfold_2", "Vestfold_3", "Vestfold_4", "Vestfold_5", 
  "Vestfold_6", "Vestfold_7", "Vestfold_8", "Vestfold_9", "Vestfold_10", 
  "Vestfold_11", "Vestfold_12", "Vestfold_13", "Vestfold_14", "Østfold_1"
)))



OnlyCrossProps_summary <- arith_mean_ci(OnlyCrossProps$proportion)

cat("Population-level mean of cross-boundary means:", round(OnlyCrossProps_summary$mean, 2), "\n")
cat("95% CI:", round(OnlyCrossProps_summary$lower_CI, 2), "-", round(OnlyCrossProps_summary$upper_CI, 2), "\n")



unique(Summaryseals_CB_proportion$NewID)


TemporalCheck <- HauloutModelFrame %>% 
  group_by(NewID) %>% 
  mutate(JDAY = lubridate::yday(S_DATE)) %>% 
  dplyr::select(NewID, S_DATE, Cross_Boundary_Management, deployment_management_region, fylkesnavn, Month, Month1, JDAY) %>% 
  distinct()


TemporalCheck <- TemporalCheck %>%
  group_by(NewID) %>%
  mutate(date1 = update(S_DATE, year = 2023))

Y1 <- TemporalCheck %>%
  filter(date1 > "2023-08-01") %>%
  mutate(date1 = update(date1, year = 2022))

Y2 <- TemporalCheck %>%
  filter(date1 <= "2023-08-01")

TemporalCheck <- Y1 %>%
  full_join(Y2) %>%
  ungroup() %>% 
  mutate(New_JDAY = lubridate::yday(date1)) %>% 
  mutate(ShortDate = as.Date(date1)) %>% 
  mutate(NumDat = as.numeric(ShortDate))
# group_by(NumDat)


month_labels <- TemporalCheck %>%
  group_by(Month1) %>%
  summarize(NumDat = min(NumDat)) 

unique(month_labels$Month1)

display_months <- c("September", "October", "November", "December", "January", "February", "March")


filtered_month_labels <- month_labels %>%
  filter(Month1 %in% display_months)



fylke_colors <- c(
  "Aust-Agder" = "#5C126EFF",
  "Vestfold"   = "#C45C4EFF",
  "Østfold"    = "#000004FF",
  "Telemark"   = "orange",
  "Vest-Agder"      = "#1E90FF",   
  "Buskerud"   = "#009E73",   
  "Sweden"       = "yellow2",  
  "Denmark"      = "red1"   
)


# Legend order 
TemporalCheck$fylkesnavn <- factor(
  TemporalCheck$fylkesnavn,
  levels = c("Vest-Agder", "Aust-Agder", "Telemark", "Vestfold", "Buskerud", "Østfold", "Sweden", "Denmark")
)


x_max <- max(TemporalCheck$NumDat, na.rm = TRUE)

x_min <- min(TemporalCheck$NumDat, na.rm = TRUE)

y_levels <- rev(c(
  "Aust-Agder_1", "Aust-Agder_2", "Aust-Agder_3",  
  "Telemark_1", "Telemark_2", "Telemark_3", "Telemark_4", "Telemark_5", 
  "Telemark_6", "Telemark_7", "Telemark_8", 
  "Vestfold_1", "Vestfold_2", "Vestfold_3", "Vestfold_4", "Vestfold_5", 
  "Vestfold_6", "Vestfold_7", "Vestfold_8", "Vestfold_9", "Vestfold_10", 
  "Vestfold_11", "Vestfold_12", "Vestfold_13", "Vestfold_14", "Østfold_1"
))

band_df <- data.frame(
  y = y_levels[seq(2, length(y_levels), by = 2)],  # every second level
  ymin = seq(1.5, length(y_levels) - 0.5, by = 2),
  ymax = seq(2.5, length(y_levels) + 0.5, by = 2)
)


x_min <- min(TemporalCheck$NumDat - 1)
x_max_with_padding <- x_max + 11


# Plot the fig
Fig2B <- ggplot() +
  geom_rect(
    data = band_df,
    aes(xmin = x_min, xmax = x_max_with_padding, ymin = ymin, ymax = ymax),
    fill = "grey70",
    alpha = 0.4,
    inherit.aes = FALSE
  ) +
  # Jitter points
  geom_jitter(data = TemporalCheck,
              aes(x = NumDat, y = NewID, group = Month1, color = as.factor(fylkesnavn)),
              width = 0.1, height = 0.1, size = 1.5, alpha = 0.7) +
  geom_text(data = OnlyCrossProps,
            aes(x = x_max + 3, y = NewID, label = label),
            color = "black", hjust = 0, size = 3.2) +  # adjust size/hjust for spacing
  
  scale_x_continuous(
    expand = c(0.01, 0),
    breaks = filtered_month_labels$NumDat,
    labels = filtered_month_labels$Month1,
    limits = c(min(TemporalCheck$NumDat- 1), x_max + 15)  # extra room for labels
  ) +
  
  scale_y_discrete(limits = levels(OnlyCrossProps$NewID)) +
  
  labs(
    x = "Temporal distribution of haul-outs",
    y = "Seal ID",
    color = NULL
  ) +
  annotate("text", x = x_max + 15, y = length(levels(OnlyCrossProps$NewID)) / 2,
           label = "Haul-outs (%) in a different management jurisdiction from tag deployment",
           angle = 270, size = 3.5, color = "black", hjust = 0.5) +
  
  scale_colour_manual(values = fylke_colors) +
  # ggsci::scale_color_flatui()
  # scale_colour_viridis_d(option = "B") +
  theme_bw() +
  theme(
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black"),
    # axis.title.y = element_blank(),
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = unit(c(0, 0, 0, 0), 'lines'),
    legend.position = "top",
    legend.background = element_rect(fill = "grey80", colour = "black"),
    legend.key = element_rect(fill = "grey80"),
    legend.key.size = unit(0.5, "cm"),  # Size of legend key box
    legend.text = element_text(size = 10),
    legend.spacing.x = unit(2, 'cm')
  ) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 3))) +  # Bigger legend dots
  scale_y_discrete(limits = rev(c(
    "Aust-Agder_1", "Aust-Agder_2", "Aust-Agder_3",  
    "Telemark_1", "Telemark_2", "Telemark_3", "Telemark_4", "Telemark_5", 
    "Telemark_6", "Telemark_7", "Telemark_8", 
    "Vestfold_1", "Vestfold_2", "Vestfold_3", "Vestfold_4", "Vestfold_5", 
    "Vestfold_6", "Vestfold_7", "Vestfold_8", "Vestfold_9", "Vestfold_10", 
    "Vestfold_11", "Vestfold_12", "Vestfold_13", "Vestfold_14", "Østfold_1"
  )))  

Fig2B

tiff("Fig2_clean.tiff", width = 9, height = 6, units = 'in', res = 300)
Fig2B
dev.off()



##############################################################################################################################################
## Do a manual check of the delineations
##############################################################################################################################################


HauloutModelFrame_sf <- HauloutModelFrame %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)


HauloutModelFrame_sf <- HauloutModelFrame_sf %>%
  mutate(
    tooltip = paste0(
      "<strong>id: </strong>", id, "<br>",
      "<strong>NewID: </strong>", NewID, "<br>",
      "<strong>S_DATE: </strong>", format(S_DATE, "%Y-%m-%d %H:%M")
    ),
    popup = tooltip  
  )

# Color palette 
manual_colors <- c(
  "#5C126EFF",  # Aust-Agder
  "#009E73",    # Buskerud
  "red2",       # Denmark
  "yellow2",    # Sweden
  "orange",     # Telemark
  "#1E90FF",    # Vest-Agder
  "#C43C4EFF",  # Vestfold
  "#000004FF"   # Østfold
)

unique_fylker <- unique(HauloutModelFrame_sf$fylkesnavn)
pal <- colorFactor(
  palette = manual_colors[seq_along(unique_fylker)],
  domain = unique_fylker
)

# Interactive plot for checking
MapPlot <- leaflet(data = HauloutModelFrame_sf) %>%
  addTiles() %>%
  addCircleMarkers(
    radius = 5,
    color = ~pal(fylkesnavn),
    stroke = TRUE,
    fillOpacity = 0.7,
    label = ~lapply(tooltip, htmltools::HTML), 
    popup = ~tooltip  
  ) %>%
  addLegend(
    position = "bottomright",
    pal = pal,
    values = ~fylkesnavn,
    title = "fylkesnavn"
  )
MapPlot

# HTML
saveWidget(MapPlot, file = "Interactive_HauloutMap.html", selfcontained = TRUE)



##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
# Create Fig.S2 to check the cross-boundary haulouts around January

start_date <- as.Date("2022-12-24")
end_date <- as.Date("2023-01-07")

JanuaryHunt <- TemporalCheck %>%
  filter(ShortDate >= start_date & ShortDate <= end_date) %>%
  distinct(NewID)

JanuaryHunt_filtered <- TemporalCheck %>%
  filter(NewID %in% JanuaryHunt$NewID) %>% 
  filter(ShortDate <= end_date)

Maxi <- JanuaryHunt_filtered %>% 
  group_by(NewID) %>% 
  summarise(MaxDate = max(ShortDate))


January_Summaryseals_CB_proportion <- JanuaryHunt_filtered %>%
  group_by(deployment_management_region, NewID, Cross_Boundary_Management) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) 


Jnuary_OnlyCrossProps <- January_Summaryseals_CB_proportion %>% 
  filter(Cross_Boundary_Management == 0) %>% 
  mutate(proportion = (1- proportion) * 100) %>% 
  mutate(
    label = paste0(ceiling(proportion), "%")
  )

# 
# Jnuary_OnlyCrossProps <- bind_rows(
#   Jnuary_OnlyCrossProps,
#   tibble(
#     deployment_management_region = "Sweden",
#     NewID = "Østfold_1",
#     Cross_Boundary_Management = 0,
#     count = 64,
#     proportion = 1,
#     label = "100%"
#   )
# )



Jnuary_OnlyCrossProps$NewID <- factor(Jnuary_OnlyCrossProps$NewID, levels = rev(c(
  "Aust-Agder_1", "Aust-Agder_2", "Aust-Agder_3",  
  "Telemark_1", "Telemark_5", 
  "Telemark_6", "Telemark_7", "Telemark_8", 
  "Vestfold_1", "Vestfold_2", "Vestfold_3", "Vestfold_4", "Vestfold_5", 
  "Vestfold_6", "Vestfold_7", "Vestfold_8","Vestfold_10", 
  "Vestfold_11", "Vestfold_12", "Vestfold_14", "Østfold_1"
)))


unique(Jnuary_OnlyCrossProps$NewID)


y_levels <- rev(c(
  "Aust-Agder_1", "Aust-Agder_2", "Aust-Agder_3",  
  "Telemark_1", "Telemark_5", 
  "Telemark_6", "Telemark_7", "Telemark_8", 
  "Vestfold_1", "Vestfold_2", "Vestfold_3", "Vestfold_4", "Vestfold_5", 
  "Vestfold_6", "Vestfold_7", "Vestfold_8","Vestfold_10", 
  "Vestfold_11", "Vestfold_12", "Vestfold_14", "Østfold_1"
))

Jnuary_OnlyCrossProps$NewID <- factor(Jnuary_OnlyCrossProps$NewID, levels = y_levels)

background_df <- data.frame(
  y = y_levels,
  index = seq_along(y_levels)
) %>%
  filter(index %% 2 == 0) %>%  # Every second row
  mutate(
    ymin = as.numeric(index) - 0.5,
    ymax = as.numeric(index) + 0.5,
    xmin = -Inf,
    xmax = Inf
  )


vline1 <- as.numeric(as.Date("2022-12-23"))
vline2 <- as.numeric(as.Date("2023-01-08"))


xmax_shade <- as.numeric(end_date + 8)


background_df <- data.frame(
  y = y_levels,
  index = seq_along(y_levels)
) %>%
  filter(index %% 2 == 0) %>%
  mutate(
    ymin = as.numeric(index) - 0.5,
    ymax = as.numeric(index) + 0.5,
    xmin = -Inf,
    xmax = xmax_shade
  )



Jan_HuntingPlot <- ggplot() +
  geom_rect(data = background_df, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey70", inherit.aes = FALSE, alpha = 0.4) +
  
  geom_vline(xintercept = vline1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = vline2, linetype = "dashed", color = "black") +
  
  geom_jitter(data = JanuaryHunt_filtered,
              aes(x = NumDat, y = NewID, group = Month1, color = as.factor(fylkesnavn)),
              width = 0.1, height = 0.1, size = 1.5, alpha = 0.7) +
  
  geom_text(data = Jnuary_OnlyCrossProps,
            aes(x = end_date + 3, y = NewID, label = label),
            color = "black", hjust = 0, size = 3.2) +
  
  scale_x_continuous(
    expand = c(0.01, 0),
    breaks = filtered_month_labels$NumDat,
    labels = filtered_month_labels$Month1,
    limits = c(min(JanuaryHunt_filtered$NumDat - 1), x_max - 68)
  ) +
  scale_y_discrete(limits = y_levels) +
  
  scale_colour_manual(values = fylke_colors) +
  
  labs(
    x = "Temporal distribution of haul-outs",
    y = "Seal ID",
    title = "Haul-out activity during the start of the Norwegian hunting period",
    subtitle = "Data displayed only for individuals transmitting data on 01 January (± 1 week)",
    color = NULL
  ) +
  
  annotate("text", x = x_max - 69, y = length(y_levels) / 1.93,
           label = "Haul-outs (%) in a different management jurisdiction from tag deployment",
           angle = 270, size = 3.5, color = "black", hjust = 0.5) +
  
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black"),
    plot.margin = unit(c(0, 0, 0, 0), 'lines'),
    legend.position = "top",
    legend.background = element_rect(fill = "grey80", colour = "black"),
    legend.key = element_rect(fill = "grey80"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 10),
    legend.spacing.x = unit(2, 'cm')
  ) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 3)))

Jan_HuntingPlot


tiff("FigS3_clean.tiff", width = 9, height = 6, units = 'in', res = 300)
Jan_HuntingPlot
dev.off()





##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

# COX MODELLING

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################


options(scipen=999)

# setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/Data/Harbour seal movement/Haulouts/")

# write_rds(HauloutModelFrame, "ProcessedHauloutData.rds")

HauloutModelFrame$Species <- "HarbourSeal"
HauloutModelFrame$id <- as.factor(HauloutModelFrame$id)
HauloutModelFrame$deployment_management_region <- as.factor(HauloutModelFrame$deployment_management_region)


HauloutModelFrame <- HauloutModelFrame %>%
  group_by(NewID) %>%
  mutate(count = 1) %>% 
  mutate(MaxHaulNuber = sum(count)) %>%
  # # filter(MaxHaulNuber >= 20) %>% 
  # ungroup() %>% 
  # arrange(S_DATE) %>% 
  # mutate(timeLub = lubridate::decimal_date(S_DATE)) %>%
  # mutate(times = glmmTMB::numFactor(timeLub)) %>% 
  # mutate(HNUM = glmmTMB::numFactor(HAULOUT_NUMBER)) %>% 
  # filter(!deployment_management_region == "Østfold") %>% 
  ungroup()


HauloutModelFrame$NewID <- as.character(HauloutModelFrame$NewID)
HauloutModelFrame$deployment_management_region <- as.character(HauloutModelFrame$deployment_management_region)

HauloutModelFrame$NewID <- as.factor(HauloutModelFrame$NewID)
HauloutModelFrame$deployment_management_region <- as.factor(HauloutModelFrame$deployment_management_region)

unique(HauloutModelFrame$deployment_management_region)

HauloutModelFrame$DistCorr <- round(HauloutModelFrame$distance_km)

HauloutModelFrame$DistCorr <- numFactor(HauloutModelFrame$distance_km)

HauloutModelFrame$pos <- numFactor(HauloutModelFrame$lon, HauloutModelFrame$lat)

HauloutModelFrame$pos_scaled <- datawizard::standardise(HauloutModelFrame$pos)

HauloutModelFrame$TagDays_St <- datawizard::normalize(HauloutModelFrame$TagDays)

HauloutModelFrame$Sex <- as.factor(HauloutModelFrame$Sex)


unique(HauloutModelFrame$NewID)

Checks <- HauloutModelFrame %>% 
  dplyr::filter(NewID == "Østfold_1")


ggplot() + geom_point(data = Checks, aes(x = lon, y = lat, colour = fylkesnavn))



CoxModelFrame <- HauloutModelFrame %>% 
  dplyr::select(NewID, deployment_management_region, fylkesnavn, lat, lon, Cross_Boundary_Management, S_DATE, E_DATE, Weight, TagDays) %>% 
  mutate(Type = "Haul")



## Some data wrangling to massage it into a cox format

CoxModelFrame <- CoxModelFrame %>% 
  # full_join(FirstDummy) %>%
  group_by(NewID) %>% 
  arrange(NewID, S_DATE) %>% 
  mutate(Event_Number = row_number()) %>% 
  mutate(TripID = paste(NewID,Event_Number))



ggplot(data = CoxModelFrame, aes(x = Event_Number, y = NewID, color = as.factor(Type))) +
  geom_jitter(width = 0.1, height = 0.1, size = 1.5, alpha = 0.7) + 
  scale_color_brewer(palette = "Set1", direction = -1) +
  labs(
    title = "Temporal Distribution of New Management Area Haul-outs",
    x = "Haul out events through time",
    # y = "Seal ID",
    color = "Type"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) 



CoxModelFrame <- CoxModelFrame %>%
  arrange(NewID, S_DATE) %>%  
  group_by(NewID) %>%  
  mutate(
    prev_region = lag(fylkesnavn, default = deployment_management_region[1]),  # Previous haul-out region
    UniqueManagementEvent = ifelse(fylkesnavn != prev_region & !is.na(prev_region), 1, 0)  # Mark first transition
  ) %>%
  ungroup() %>% 
  mutate(JDAY = lubridate::yday(S_DATE),
         Survival = UniqueManagementEvent) %>% 
  # mutate(Type = "Haulout") %>% 
  group_by(NewID) %>% 
  arrange(NewID, S_DATE) 


ggplot(data = CoxModelFrame, aes(x = Event_Number, y = NewID, color = as.factor(Survival))) +
  geom_jitter(width = 0.1, height = 0.1, size = 1.5, alpha = 0.7) + 
  scale_color_brewer(palette = "Set1", direction = -1) +
  labs(
    # title = "Temporal Distribution of New Management Area Haul-outs",
    x = "Haul out events through time",
    # y = "Seal ID",
    color = "Cross Boundary"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) 


CoxModelFrame_Time <- CoxModelFrame %>%
  # full_join(TagStart) %>% 
  ungroup() %>% 
  group_by(NewID) %>%
  arrange(NewID, S_DATE) %>% 
  mutate(Previous_Haul_Time = lag(E_DATE)) %>% 
  mutate(InterHaulOutTime = difftime(S_DATE, lag(E_DATE), units = "hours")) %>%
  mutate(InterHaulOutTime_Hours = as.numeric(InterHaulOutTime)) %>% 
  mutate(Overall_Start_Time = first(E_DATE),
         Time_Since_Tagging = (as.numeric(S_DATE) - as.numeric(Overall_Start_Time))/60/60/24)



CoxModelFrame_Time <- CoxModelFrame_Time %>%
  mutate(Time_Since_Tagging = ifelse(Time_Since_Tagging < 0, 0, Time_Since_Tagging))

IDs_with_event <- CoxModelFrame_Time %>%
  filter(UniqueManagementEvent == 1) %>%
  pull(NewID) %>% unique()

FilteredData1 <- CoxModelFrame_Time %>%
  filter(NewID %in% IDs_with_event)


SubsetData1 <- FilteredData1 %>%
  group_by(NewID) %>%
  filter(row_number() == 1 | UniqueManagementEvent == 1) %>%
  arrange(NewID, S_DATE) %>%
  mutate(SurvivalTime = as.numeric(difftime(S_DATE, lag(E_DATE), units = "days"))) %>%
  ungroup()


IDs_without_event <- CoxModelFrame_Time %>%
  group_by(NewID) %>%
  filter(all(UniqueManagementEvent != 1)) %>%
  pull(NewID) %>% unique()

FilteredData2 <- CoxModelFrame_Time %>%
  filter(NewID %in% IDs_without_event)

SubsetData2 <- FilteredData2 %>%
  group_by(NewID) %>%
  slice_tail(n = 1) %>%
  ungroup() %>% 
  mutate(SurvivalTime = Time_Since_Tagging)



SealSurvival <- SubsetData1 %>% 
  full_join(SubsetData2)



SealSurvival$Jurisdiction <- SealSurvival$deployment_management_region
SealSurvival$prev_region <- as.factor(SealSurvival$prev_region)


SexMeta <- HauloutModelFrame %>% 
  dplyr::select(id, NewID, Sex) %>% 
  distinct()

SealSurvival <- SealSurvival %>% 
  left_join(SexMeta) %>% 
  dplyr::select(NewID, Survival, SurvivalTime, Time_Since_Tagging, prev_region, Weight, TagDays, Sex,) %>% 
  drop_na() 
# group_by(NewID) %>% 
# filter(Time_Since_Tagging  <= min(Time_Since_Tagging ))

SealSurvival$prev_region <- as.character(SealSurvival$prev_region)
SealSurvival$prev_region <- as.factor(SealSurvival$prev_region)

table(SealSurvival$prev_region)

table(SexMeta$Sex)


print(SealSurvival, n = nrow(SealSurvival))

##############################################################################################################################################
##############################################################################################################################################
## Run the model
##############################################################################################################################################
##############################################################################################################################################

# # In this DF, Survival represents a presence/absense for a cross-boundary haul out events and survival time is time taken to perform cross-boundary events-
# # Only 0's for the individuals that never crossed boundaries (censored individuals)


SealSurvival <- within(SealSurvival, prev_region <- relevel(prev_region, ref = "Sweden")) ## Make Sweden the reference level as the highest coeff


# library(writexl)
# write_xlsx(SealSurvival, "CPH_ModelFrame.xlsx")


cox_mod <- coxme(Surv(SurvivalTime, Survival) ~ prev_region  + Weight + (1 |NewID),
                 data = SealSurvival)

summary(cox_mod)

# ph_test <- cox.zph(cox_mod)
# print(ph_test)
# plot(ph_test)
# martingale_resid <- residuals(cox_mod, type = "martingale")
# plot(SealSurvival$Weight, martingale_resid)


###############################################################################################################################################
# Cross validation exercise of the model
##############################################################################################################################################


regions <- levels(droplevels(SealSurvival$prev_region))


# Loop over each model where there are difference reference levels
coef_results <- map_dfr(regions, function(ref_level) {
  SealSurvival$prev_region_relevel <- relevel(SealSurvival$prev_region, ref = ref_level) 
  
  mod <- coxme(Surv(SurvivalTime, Survival) ~ prev_region_relevel + Weight + (1 | NewID),
               data = SealSurvival)
  
  coefs <- fixef(mod)
  ses <- sqrt(diag(vcov(mod)))
  
  term_idx <- grep("prev_region_relevel", names(coefs))
  
  df <- tibble(
    term = names(coefs)[term_idx],
    estimate = coefs[term_idx],
    std.error = ses[term_idx],
    reference = ref_level
  ) %>%
    mutate(
      term = gsub("prev_region_relevel", "", term),
      comparison = paste(term, "vs", reference),
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      p.value = 2 * (1 - pnorm(abs(estimate / std.error))),
      signif = p.value < 0.05
    )
  
  return(df)
})


# Create the neat cross validation plot

NewSup <- ggplot(coef_results, aes(x = term, y = estimate, color = signif)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), fatten = 2, size = 0.8) +
  facet_wrap(~ reference, scales = "free_x", nrow = 2) +
  scale_color_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"),
                     labels = c("FALSE", "TRUE"),
                     name = "Statistical\nsignificance") +
  labs(
    title = "Cross-validation of CPH models",
    subtitle = "Configured with different jurisdictions as the reference levels",
    x = "Jurisdiction",
    y = "Hazard Ratio"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right")
NewSup


# tiff("Fig_S3.tiff", width = 9, height = 6, units = 'in', res = 300)
# NewSup
# dev.off()


print(coef_results, n = nrow(coef_results))



##############################################################################################################################################
## Survival Curve Plotting
##############################################################################################################################################


fylke_colors <- c(
  "Aust-Agder" = "#5C126EFF",
  "Vestfold"   = "#C45C4EFF",
  "Østfold"    = "#000004FF",
  "Telemark"   = "orange",
  "Vest-Agder"      = "#1E90FF",   
  "Buskerud"   = "#009E73",   
  "Sweden"       = "yellow2",   
  "Denmark"      = "red1"    
)

fylke_colors1 <- c(
  "prev_region=Aust-Agder" = "#5C126EFF",
  "prev_region=Vestfold"   = "#C45C4EFF",
  "prev_region=Østfold"    = "#000004FF",
  "prev_region=Telemark"   = "orange",
  "prev_region=Vest-Agder"      = "#1E90FF",  
  "prev_region=Buskerud"   = "#009E73",   
  "prev_region=Sweden"       = "yellow2",   
  "prev_region=Denmark"      = "red1"    
)


legend_labels <- c(
  "Aust-Agder",
  "Vestfold",
  "Østfold",
  "Telemark",
  "Vest-Agder",
  "Buskerud",
  "Sweden",
  "Denmark"
)


Fig3A <- ggsurvplot(survfit(Surv(SurvivalTime, Survival) ~ prev_region, data = SealSurvival),
                    # Change legends: title & labels
                    legend.title = "Jurisdiction",
                    # legend.labs = c("Vest-Agder", "Aust-Agder", "Buskerud","Sweden", "Telemark","Vestfold","Østfold"),
                    # pval = TRUE,
                    legend = "right",
                    conf.int = FALSE,
                    # Add risk table
                    risk.table = FALSE,
                    xlab = "Time (days)",
                    tables.height = 0.2,
                    tables.theme = theme_cleantable(),
                    censor = FALSE,
                    # surv.median.line = "hv",
                    palette = alpha(c(fylke_colors1), 0.8),
                    # legend.labs = c("Vest-Agder","Telemark", "Aust-Agder", "Buskerud","Sweden","Vestfold","Østfold"),
                    ggtheme = theme_bw())


Fig3A <- Fig3A$plot


Fig3A


##############################################################################################################################################
## Hazards Ratio Plotting
##############################################################################################################################################


multiple_cox_survival <- emmeans::emmeans(cox_mod, list(pairwise~prev_region), adj="Tukey")

summary(multiple_cox_survival)


coeff_cox_survival <- as.data.frame(multiple_cox_survival$`emmeans of prev_region`)%>%
  group_by(prev_region)%>%
  summarise(DCV.Coef=emmean, DCV.Coef.SE=SE)%>%
  mutate(DCV.Coef=DCV.Coef-DCV.Coef[prev_region=="Sweden"],
         prev_region=factor(prev_region, levels=c("Vest-Agder","Telemark", "Aust-Agder", "Buskerud","Sweden","Vestfold","Østfold")))

coeff_cox_survival


Fig3B <- ggplot(coeff_cox_survival, aes(x=reorder(prev_region, desc(prev_region)), y=DCV.Coef, fill=prev_region))+
  geom_col( alpha = 0.8, colour = "black")+
  geom_errorbar(aes(ymin=DCV.Coef-DCV.Coef.SE, ymax=DCV.Coef+DCV.Coef.SE), width=0.2, size=0.7)+
  scale_fill_manual(values= c("#1E90FF", "orange", "#5C126EFF", "#009E73", "yellow2" , "#C45C4EFF", "#000004FF" ))+
  scale_colour_manual(values= c( "orange", "black", "black" ,  "black" , "black" ,  "black"  ,"black"   , "black" ))+
  
  theme_bw() + labs(x = "Jurisdiction", y = "Coefficients ± SE",fill = "Jurisdiction", colour = "Jurisdiction") +
  theme(legend.position="none") + coord_flip() +
  theme(    axis.title.y = element_blank(),
            # axis.text.y = element_blank(),
            axis.ticks.y = element_blank())

Fig3B


Fig3<- ggarrange(Fig3B, Fig3A, nrow = 2, common.legend = TRUE, labels = "AUTO", legend = "right")
Fig3


tiff("Fig3_clean.tiff", width = 8, height = 6, units = 'in', res = 300)
Fig3
dev.off()




##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

# SPATIAL NETWORK ANALYSIS

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################


ggplot() + 
  geom_point(data = HauloutModelFrame, aes(x= lon, y= lat))+
  geom_point(data = MetaData1, aes(x= Col_lon, y= Col_lat), colour = "red")+
  facet_wrap(~Capture_site, scales = "free")


Spattest <- HauloutModelFrame %>% 
  ungroup() %>% 
  left_join(TemporalCheck) %>% 
  mutate(
    date1 = if_else(
      is.na(date1),
      lead(date1) - seconds(1),
      date1
    ),
    ShortDate = if_else(
      is.na(ShortDate),
      lead(ShortDate),
      ShortDate
    )) %>% 
  mutate(ID = as.character(NewID),
         X = lon,
         Y = lat,
         datetime = date1,
         datetime1 = first(datetime),
         population = as.character(deployment_management_region)) %>% 
  dplyr::select(ID, X, Y, datetime,datetime1, population) %>% 
  st_as_sf(coords = c("X", "Y"), crs = 4326) %>%   
  st_transform(crs = 32632) %>% 
  mutate(
    X = st_coordinates(.)[,1],  
    Y = st_coordinates(.)[,2]  
  )


Check <- Spattest %>% 
  filter(population == "Aust-Agder") 
# filter(ID == "Aust-Agder_2")

ggplot() + 
  geom_point(data = Check, aes(x= X, y= Y, colour = ID), size = 5) +
  geom_path(data = Check, aes(x= X, y= Y, colour = ID), size = 0.1) 




unique(Spattest$population)


##############################################################################################################################################
## Spatially grouping into nodes
##############################################################################################################################################


DT <- as.data.table(Spattest)

DT[, datetime1 := as.POSIXct(datetime1)]

DT[, yr := year(datetime)]


## Create dummy temporal group for our purposes 
group_times(DT, datetime = 'datetime1', threshold = '365 days')

utm <- 32632

group_pts(DT, threshold = 10000, id = 'ID',
          coords = c('X', 'Y'), timegroup = "timegroup") 


N_IDs <- DT[, uniqueN(ID), by = group]


DT_Centroid <- DT


centroid_group(DT_Centroid, coords = c('X', 'Y'), group = 'group', na.rm = TRUE)

# Spatial group now
DT_Centroid <- DT_Centroid %>%
  as.data.frame() %>% 
  dplyr::select(ID, group, centroid_X, centroid_Y) %>% 
  st_as_sf(coords = c("centroid_X", "centroid_Y"), crs = 32632) %>% 
  st_transform(crs = 4326) %>% 
  mutate(
    centroid_lon = st_coordinates(.)[,1],
    centroid_lat = st_coordinates(.)[,2]) %>% 
  dplyr::select(ID, group, centroid_lon, centroid_lat) %>% 
  distinct()



plotrr <- DT %>% 
  as.data.frame() %>% 
  st_as_sf(coords = c("X", "Y"), crs = 32632) %>%  
  st_transform(crs = 4326) %>% 
  mutate(SpatialGrouping = as.factor(group)) %>% 
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2])



##############################################################################################################################################
## Calculate some basic spatial network measures
##############################################################################################################################################

## Number of individuals within each node
node_sizes <- DT %>%
  group_by(group) %>%
  summarise(num_individuals = n_distinct(ID))  

## Weighted edges per ID
movements <- DT %>%
  arrange(ID) %>%  
  group_by(ID) %>%
  mutate(next_group = lead(group)) %>%  
  filter(!is.na(next_group) & group != next_group) %>%  
  count(group, next_group, name = "weight")  

individual_transitions <- DT %>%
  arrange(ID, datetime1) %>%  
  group_by(ID) %>%
  mutate(next_group = lead(group)) %>%
  filter(!is.na(next_group) & group != next_group) %>%
  rowwise() %>%
  mutate(edge = paste(sort(c(group, next_group)), collapse = "-")) %>%
  ungroup()

# Count unique individuals per edge
individuals_per_edge <- individual_transitions %>%
  distinct(edge, ID) %>%
  count(edge, name = "unique_individuals")

# Overall weight of edges
edges <- movements %>% 
  rowwise() %>%
  mutate(edge = paste(sort(c(group, next_group)), collapse = "-")) %>% 
  group_by(edge) %>%
  summarise(from = first(group), to = first(next_group), weight = sum(weight))

edges <- edges %>%
  left_join(individuals_per_edge, by = "edge") %>%
  mutate(unique_individuals = replace_na(unique_individuals, 0))


nodes <- DT_Centroid %>%
  left_join(node_sizes, by = "group") %>%
  replace_na(list(num_individuals = 1))  




nodes <- nodes %>% distinct(group, .keep_all = TRUE) %>% 
  st_drop_geometry() %>% 
  as.data.frame() %>% 
  ungroup() %>%
  dplyr::select(-ID)  %>% 
  left_join(DT_Centroid) %>% 
  distinct(group, .keep_all = TRUE) 



edges <- edges %>%
  filter(from %in% nodes$group & to %in% nodes$group)



nodes$group <- as.character(nodes$group)
edges$from <- as.character(edges$from)
edges$to <- as.character(edges$to)

nodes$group <- as.factor(nodes$group)
edges$from <- as.factor(edges$from)
edges$to <- as.factor(edges$to)




# Node Degree 
node_degrees <- edges %>%
  dplyr::select(from, to) %>%
  pivot_longer(cols = everything(), values_to = "group") %>%
  count(group, name = "Degree") 


# Weighted degree 
node_w_degrees <- edges %>%
  pivot_longer(cols = c(from, to), values_to = "group") %>%
  group_by(group) %>%
  summarise(WeightedDegree = sum(weight), .groups = "drop")


# Join both to nodes
nodes <- nodes %>%
  left_join(node_degrees, by = "group") %>%
  left_join(node_w_degrees, by = "group") %>%
  mutate(Degree = replace_na(Degree, 0),
         WeightedDegree = replace_na(WeightedDegree, 0))




distinct_palette <- brewer.pal(12, "Paired")
distinct_palette
# "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"


nodes$group <- factor(nodes$group, levels = as.character(sort(as.numeric(unique(nodes$group)))))



NorwayLandMap_crop1 <- NorwayLandMap_crop %>%
  mutate(
    Y = case_when(
      abbreviation %in% c("VA") ~ Y - 0.8,
      TRUE ~ Y
    ),
    X = case_when(
      abbreviation %in% c("VA") ~ X + 0.3,
      TRUE ~ X
    )
  ) 


BASEMAP1 <- ggplot() +
  geom_tile(data = SkagerrakBathymetry_NA, aes(x = x, y = y), fill = "#eeeac4") +
  
  ggnewscale::new_scale_fill() +
  
  geom_sf(data = Combined_Maritime_Boundaries_crop, aes(fill = Area), fill = "lightblue", colour = "black", alpha = 0.3) +
  guides(fill = FALSE) +
  
  ggnewscale::new_scale_fill() +
  
  geom_sf(data = NorwayLandMap_crop1, fill = "#eeeac4", color = "black", size = 0.3) +
  
  # Adjusted text labels
  geom_text(data = NorwayLandMap_crop, aes(x = X, y = Y, label = abbreviation),
            size = 3, fontface = "bold") +
  
  labs(x = "Longitude", y = "Latitude", fill = "Harbour seal\nmanagement\nunits") +
  coord_sf(
    ylim = c(st_bbox(Combined_Maritime_Boundaries_crop)["ymin"] + 0.1, 
             st_bbox(Combined_Maritime_Boundaries_crop)["ymax"] - 0.1),
    expand = FALSE
  ) +
  theme_minimal()




# ggsave("myBASEMAP1.png", width = 5, height = 5.5, units = "in")





Nodes_Check <- BASEMAP1 + 
  # Plot the nodes
  geom_spatial_point(data = nodes, 
                     aes(x = centroid_lon, y = centroid_lat), 
                     shape = 21, stroke = 1, fill = "white", colour = "white", size = 4, alpha = 0.5) +
  labs(fill = "Spatial\nclustering", x = "Longitude", y = "Latitude") +
  
  theme_bw(10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black"),
    legend.key = element_rect(fill = "lightgrey", colour = NA),
    # legend.justification = c(0.1, 0.1),
    # legend.position = c(0.1, 0.1),
    legend.position = "none",
    legend.background = element_rect(fill = "lightgrey", colour = "black")
  ) + 
  geom_spatial_text(data = nodes, aes(x = centroid_lon, y = centroid_lat, label = group), colour = "red2" , size = 5, fontface = "bold")


# # Display the plot
# Nodes_Check
# 
# tiff("NodesCheck.tiff", width = 5, height = 7, units = 'in', res = 300)
# Nodes_Check
# dev.off()
# # 


##############################################################################################################################################
## Create the NODES figure for manuscript
##############################################################################################################################################


node_sizes$group <- as.factor(node_sizes$group)


nodes <- nodes %>%
  left_join(node_sizes, by = "group")  

edges_plot <- edges %>%
  left_join(nodes, by = c("from" = "group")) %>%
  rename(from_lon = centroid_lon, from_lat = centroid_lat) %>%
  left_join(nodes, by = c("to" = "group")) %>%
  rename(to_lon = centroid_lon, to_lat = centroid_lat)


nodes$group <- factor(nodes$group, levels = as.character(sort(as.numeric(unique(nodes$group)))))
group_levels <- levels(nodes$group)


nodes$group <- factor(nodes$group, levels = sort(as.numeric(as.character(nodes$group))))




####################################################################################################################################
##############################################################################
# some additional tweaks based on reviewer comments
##############################################################################
########################################################################################################################################################################################################################################################################



# Unique numeric to letter codes
node_mapping <- data.frame(
  old_group = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"),
  new_group = c("D", "F", "H", "J", "E", "B", "C", "A", "K", "M", "L", "G", "I"),
  stringsAsFactors = FALSE
)


DT <- DT %>%
  mutate(group = as.character(group)) %>%
  left_join(node_mapping, by = c("group" = "old_group")) %>%
  mutate(group = new_group) %>%
  select(-new_group)


DT_Centroid <- DT_Centroid %>%
  mutate(group = as.character(group)) %>%
  left_join(node_mapping, by = c("group" = "old_group")) %>%
  mutate(group = new_group) %>%
  select(-new_group)


##############################################################################
# Calculate distances between nodes and tracking effort
##############################################################################


# Calculate number of individuals tracked per node
sampling_effort <- DT %>%
  group_by(group) %>%
  summarise(
    n_individuals = n_distinct(ID),
    n_observations = n(),
    .groups = "drop"
  )

# Calculate distances between all node pairs
node_coords <- DT_Centroid %>%
  st_drop_geometry() %>%
  distinct(group, centroid_lon, centroid_lat)

from_nodes <- node_coords %>%
  rename(from = group, from_lon = centroid_lon, from_lat = centroid_lat)

to_nodes <- node_coords %>%
  rename(to = group, to_lon = centroid_lon, to_lat = centroid_lat)

node_distances <- expand.grid(
  from = from_nodes$from,
  to = to_nodes$to,
  stringsAsFactors = FALSE
) %>%
  left_join(from_nodes, by = "from") %>%
  left_join(to_nodes, by = "to") %>%
  filter(from != to) %>%
  rowwise() %>%
  mutate(
    distance_km = distHaversine(c(from_lon, from_lat), c(to_lon, to_lat)) / 1000
  ) %>%
  ungroup()

##############################################################################
## Calculate edge weights (i.e. average trips per ID)
##############################################################################


# Get raw movements
individual_transitions <- DT %>%
  arrange(ID, datetime1) %>%  
  group_by(ID) %>%
  mutate(next_group = lead(group)) %>%
  filter(!is.na(next_group) & group != next_group) %>%
  rowwise() %>%
  mutate(edge = paste(sort(c(group, next_group)), collapse = "-")) %>%
  ungroup()

# Count unique individuals per edge
individuals_per_edge <- individual_transitions %>%
  distinct(edge, ID) %>%
  count(edge, name = "unique_individuals")

# Overall weight of edges
movements <- DT %>%
  arrange(ID) %>%  
  group_by(ID) %>%
  mutate(next_group = lead(group)) %>%  
  filter(!is.na(next_group) & group != next_group) %>%  
  count(group, next_group, name = "weight")

edges <- movements %>% 
  rowwise() %>%
  mutate(edge = paste(sort(c(group, next_group)), collapse = "-")) %>% 
  group_by(edge) %>%
  summarise(
    from = first(group), 
    to = first(next_group), 
    raw_weight = sum(weight),
    .groups = "drop"
  )

edges <- edges %>%
  left_join(individuals_per_edge, by = "edge") %>%
  mutate(unique_individuals = replace_na(unique_individuals, 0))

# Now add distance information
edges <- edges %>%
  left_join(node_distances %>% 
              select(from, to, distance_km), 
            by = c("from", "to"))

# Stamdarise the edge weights
edges <- edges %>%
  mutate(
    trips_per_individual = raw_weight / unique_individuals,
    
    trips_per_individual_per_km = trips_per_individual / distance_km
  )


##############################################################################
## Create network and calculate centrality metrics
##############################################################################


nodes <- DT_Centroid %>%
  left_join(sampling_effort, by = "group") %>%
  replace_na(list(n_individuals = 1, n_observations = 0)) %>%
  distinct(group, .keep_all = TRUE) %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  ungroup() %>%
  dplyr::select(-ID)

nodes$group <- as.character(nodes$group)
edges$from <- as.character(edges$from)
edges$to <- as.character(edges$to)

edges <- edges %>%
  filter(from %in% nodes$group & to %in% nodes$group)

# Distance-corrected weights for betweenness calculation
graph_network <- graph_from_data_frame(
  d = edges %>% select(from, to, weight = trips_per_individual_per_km),
  vertices = nodes %>% 
    select(group), directed = FALSE)

# Node network metrics
nodes_metrics <- data.frame(
  group = V(graph_network)$name,
  stringsAsFactors = FALSE
)

# Degree 
nodes_metrics$degree <- degree(graph_network)

# Betweenness centrality
nodes_metrics$betweenness <- betweenness(
  graph_network, 
  weights = 1/E(graph_network)$weight,
  normalized = TRUE
)

# Closeness centrality
nodes_metrics$closeness <- closeness(
  graph_network,
  weights = 1/E(graph_network)$weight,
  normalized = TRUE
)

# Merge with node data
nodes_complete <- nodes %>%
  left_join(nodes_metrics, by = "group")



##############################################################################
##  Create Fig. 4
##############################################################################


# Plotting frame
edges_plot <- edges %>%
  left_join(nodes_complete %>% select(group, centroid_lon, centroid_lat), 
            by = c("from" = "group")) %>%
  rename(from_lon = centroid_lon, from_lat = centroid_lat) %>%
  left_join(nodes_complete %>% select(group, centroid_lon, centroid_lat), 
            by = c("to" = "group")) %>%
  rename(to_lon = centroid_lon, to_lat = centroid_lat)


distinct_pub_palette1 <- c(
  "#999933",  
  "#B15928",  
  "#B2DF8A",  
  "#4477AA",  
  "#88CCEE",  
  
  "#117733",  
  "#F0E442",  
  "#E31A1C",  
  "#FFFFFF" ,  
  "pink",     
  "#AA4499",  
  "#E69F00",  
  "#332288")



nodes_complete$group <- factor(nodes_complete$group, 
                               levels = sort(unique(nodes_complete$group)))

# Some bounding box things
nodes_sf <- st_as_sf(nodes_complete, coords = c("centroid_lon", "centroid_lat"), crs = 4326)
nodes_sf <- st_transform(nodes_sf, crs = 32632)
convex_hull <- st_convex_hull(st_union(nodes_sf))
buffered_hull <- st_buffer(convex_hull, dist = 10000) 
buffered_hull <- st_transform(buffered_hull, crs = 4326)
bbox <- st_bbox(buffered_hull)


Fig4 <- BASEMAP1 +
  # Edges sized by edge weight
  geom_curve(
    data = edges_plot,
    aes(x = from_lon, y = from_lat, xend = to_lon, yend = to_lat, 
        linewidth = unique_individuals),
    color = "grey40",
    curvature = 0.2,
    alpha = 0.6
  ) +
  scale_linewidth_continuous(
    name = "Edge\noccupancy",
    range = c(0.3, 2.6),
    breaks = c(1, 2, 3, 6, 7),
    guide = guide_legend(override.aes = list(
      color = "grey40",
      linewidth = c(0.8, 1.2, 1.7, 2.3, 2.8)
    ))
  ) +
  
  ggnewscale::new_scale("size") +
  
  # Nodes sized by betweenness centrality
  geom_spatial_point(
    data = nodes_complete,
    aes(x = centroid_lon, y = centroid_lat, 
        size = betweenness, 
        fill = group),
    shape = 21,
    colour = "black",
    stroke = 1,
    alpha = 0.9
  ) +
  
  scale_fill_manual(
    values = distinct_pub_palette1, 
    guide = "none"
  ) +
  
  scale_size_continuous(
    name = "Node\nbetweenness\ncentrality",
    range = c(3, 10),
    # breaks = c(0, 1, 2, 3, 4, 5, 6),
    guide = guide_legend(override.aes = list(shape = 21, fill = "white", color = "black"))
  ) +
  
  labs(x = "Longitude", y = "Latitude") +
  
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black"),
    legend.key = element_rect(fill = "lightgrey", colour = NA),
    legend.background = element_rect(fill = "lightgrey", colour = "black"),
    legend.position = "right"
  ) +
  
  coord_sf(
    xlim = c(bbox["xmin"] - 0.1, bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"]),
    expand = FALSE
  )

tiff("Fig4.tiff", width = 6, height = 6, units = 'in', res = 300)
Fig4
dev.off()



##############################################################################
## Quicke summary tables
##############################################################################


# Table 1: Node metrics
table1_nodes <- nodes_complete %>%
  arrange(desc(betweenness)) %>%
  select(
    Node = group,
    `N individuals` = n_individuals,
    # `N observations` = n_observations,
    # Degree = degree,
    `Betweenness centrality` = betweenness
  ) 


write.csv(table1_nodes, "Table1_Node_Metrics.csv", row.names = FALSE)

# Table 2: Edge metrics
table2_edges <- edges %>%
  arrange(desc(trips_per_individual)) %>%
  select(
    From = from,
    To = to,
    `Raw trips` = raw_weight,
    `N individuals` = unique_individuals,
    `Distance (km)` = distance_km,
    `Trips per individual` = trips_per_individual,
    `Trips per ind per km` = trips_per_individual_per_km
  ) %>%
  mutate(
    across(where(is.numeric) & !matches("Raw trips|N individuals"), 
           ~round(.x, 3)),
    `Distance (km)` = round(`Distance (km)`, 1),
    `Raw trips` = as.integer(`Raw trips`),
    `N individuals` = as.integer(`N individuals`)
  )


write.csv(table2_edges, "Table2_Edge_Metrics.csv", row.names = FALSE)

# Some Network-level statistics
table3_network <- data.frame(
  Metric = c(
    "Number of nodes",
    "Number of edges",
    "Network density",
    "Mean degree",
    "Mean betweenness centrality",
    "Clustering coefficient",
    "Assortativity"
  ),
  Value = c(
    vcount(graph_network),
    ecount(graph_network),
    round(edge_density(graph_network), 3),
    round(mean(degree(graph_network)), 2),
    round(mean(betweenness(graph_network, weights = 1/E(graph_network)$weight, 
                           normalized = TRUE)), 3),
    round(transitivity(graph_network, type = "global"), 3),
    round(assortativity_degree(graph_network), 3)
  )
)


write.csv(table3_network, "Table3_Network_Statistics.csv", row.names = FALSE)



####################################################################################################################################
# Conceptual Figure

## Code to create the basic idea in R. Then export to Inkscape and do manual edits

########################################################################################################################################################################################################################################################################
# 
# library(ggforce)
# 
# nodes <- data.frame(
#   name = LETTERS[1:6],
#   x = c(1, 2, 3, 2, 4.5, 5),
#   y = c(1, 2.5, 1, 0, 2.5, 0.5)
# )
# 
# edges <- data.frame(
#   from = c("A", "A", "B", "C", "D", "E"),
#   to   = c("B", "C", "C", "D", "C", "F")
# )
# 
# graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
# 
# ConceptualPlot <- ggraph(graph, layout = "manual", x = nodes$x, y = nodes$y) +
#   
#   geom_circle(aes(x0 = x, y0 = y, r = 0.8), 
#               data = nodes, linetype = "dashed", color = "grey70") +
#   
#   geom_edge_link(color = "grey40", width = 1.2) +
#   
#   geom_node_point(size = 10, fill = "skyblue", shape = 21, color = "black") +
#   
#   geom_node_text(aes(label = name), size = 4, fontface = "bold", vjust = -1.5) +
#   
#   geom_text(data = data.frame(x = 1.9, y = 3.2, label = "Node degree = number of connections"), 
#             aes(x, y, label = label), size = 4, hjust = 0, fontface = "italic") +
#   
#   geom_text(data = data.frame(x = 1.9, y = 2.8, label = "Edge density = total number of edges"), 
#             aes(x, y, label = label), size = 4, hjust = 0, fontface = "italic") +
#   
#   coord_equal() +
#   theme_void() +
#   labs(title = "Conceptual spatial network analysis of seal haul-out events") +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )
# 
# ConceptualPlot
# 
# 
# # pdf("ConceptualFigure.pdf", width = 8, height = 8)
# # ConceptualPlot
# # dev.off()
# 

