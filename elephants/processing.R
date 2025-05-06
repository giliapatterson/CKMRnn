library(tidyverse)
library(sf)
library(lubridate)

# Read in protected areas in Uganda
pas <- st_read("empirical_data/Ug_Protected-areas2007/Ug_Protected-areas2007.shp")

# Remove all but Kibale and Kisangi (DJM)
kibale <- filter(pas, NAME == "Kibale", GAZTYPE == "NP")
kisangi <- filter(pas, NAME == "Kisangi", GAZTYPE == "DJM")
kisangi2 <- filter(pas, NAME == "Kisangi", GAZTYPE == "CFR")
kibale_new <- st_simplify(st_union(rbind(kibale, kisangi)), preserveTopology = TRUE, dTolerance = 0.01)

# Limits of the map
min_long = st_bbox(kibale_new)$xmin
max_long = st_bbox(kibale_new)$xmax
min_lat = st_bbox(kibale_new)$ymin
max_lat = st_bbox(kibale_new)$ymax
bb <- st_as_sfc(st_bbox(kibale_new))

## Calculate dimensions of the image and of the SLiM simulation map ##
# Dimensions of the kibale map in meters
corners <- tibble(long = c(min_long, max_long, min_long),
                  lat = c(min_lat, min_lat, max_lat)) %>% 
  st_as_sf(coords = c("long", "lat"), crs = st_crs(kibale))
distances <- st_distance(corners)
width_meters <- as.numeric(distances[1, 2])
height_meters <- as.numeric(distances[1, 3])
# Dimensions of the kibale map image in pixels
# Number of pixels for the saved image
w = round(width_meters/10)
h = round(height_meters/10)
# Dimensions to use in slim so that one unit is one kilometer
MAX_WIDTH <- width_meters/1000
MAX_HEIGHT <- height_meters/1000

# Plot black and white image and save
ggplot(kibale_new) +
  geom_sf(data = bb, fill = "white", linewidth = 0) +
  geom_sf(fill = "black") +
  coord_sf(expand = FALSE) +
  theme_void()
ggsave("simulations/kibale.png", width = w, height = h, units = "px", limitsize = FALSE)

## Take survival probabilities for five year increments from:
## Armbruster P, Lande R. A Population Viability Analysis for African Elephant (Loxodonta africana): How Big Should Reserves Be? Conservation Biology. 1993;7:602–10.
## Convert to yearly survival probabilities and write to a file
five_year_s <- tibble(Age = 0:60,
                      F5 = c(rep(0.5,5), rep(0.887,5), rep(0.884,5), rep(0.898,5),
                             rep(0.905,5), rep(0.883,5), rep(0.881,5), rep(0.875,5),
                             rep(0.857,5), rep(0.625,5), rep(0.4,5), rep(0,6)),
                      M5 = c(rep(0.5,5), rep(0.887,5), rep(0.884,5), rep(0.898,5),
                             rep(0.923,5), rep(0.694,5), rep(0.674,5), rep(0.695,5),
                             rep(0.622,5), rep(0.118,5), rep(0.333,5), rep(0,6)))
one_year_s <- mutate(five_year_s, 'F' = F5^0.2, M = M5^0.2)
write.csv(one_year_s, 'simulations/one_year_survival.csv', row.names = F)

## Convert latitude and longitude to slim x and y
lat_to_y <- function(lat, MAX_HEIGHT, h, min_lat, max_lat){
  px_lat = (max_lat-min_lat)/h
  slim_y = MAX_HEIGHT*(lat - min_lat - (1/2)*px_lat)/(max_lat-min_lat - px_lat)
  return(slim_y)
}
long_to_x <- function(long, MAX_WIDTH, w, min_long, max_long){
  px_long = (max_long-min_long)/w
  slim_x = MAX_WIDTH*(long - min_long - (1/2)*px_long)/(max_long-min_long - px_long)
  return(slim_x)
}

## Find potential sampling locations and save ##
# Create a grid of points approximately 1 km on each edge, intersect with kibale
potential_xs = seq(min_long, max_long, length = round(width_meters/1000))
potential_ys = seq(min_lat, max_lat, length = round(height_meters/1000))
grid <- expand_grid(lat = potential_ys, long = potential_xs) %>% 
  st_as_sf(coords = c("long", "lat"), crs = st_crs(kibale))
kibale_grid <- st_intersection(grid, kibale_new)
# Convert lat long grid to slim x y
kibale_grid_slim <- as_tibble(st_coordinates(kibale_grid)) %>% 
  rename(lat = Y, long = X) %>% 
  mutate(x = round(long_to_x(long, MAX_WIDTH, w, min_long, max_long), 3),
         y = round(lat_to_y(lat, MAX_HEIGHT, h, min_lat, max_lat),3)) %>%
  select(x, y)
write.csv(kibale_grid_slim, "simulations/sampling_locations.txt", row.names = FALSE)
# Width and height of sampling grid cells
DW = diff(sort(unique(kibale_grid_slim$x))[c(1,2)])
DH = diff(sort(unique(kibale_grid_slim$y))[c(1,2)])

## Extract sampling times and sample sizes from empirical data ##
samples <- read_csv("empirical_data/Kibale_Species_Classification_Final.csv")
# Create a column that is a dates object and sort by date
samples <- mutate(samples, date = mdy(samples$`Date collected (MM/DD/YY)`)) %>%
  arrange(date)
# Compute number of days since the first sampling date for each sampling date
samples <- mutate(samples, slim_days = c(0, cumsum(int_diff(date)/days(1)))+1)
# SLiM sampling days AND YEARS
if(max(samples$slim_days) <= 365){
  SAMPLE_DAYS = unique(samples$slim_days)
  SAMPLE_YEARS = 200
  samples$slim_years = 200
}
# Write sampling times to file
sample_times <- tibble(SAMPLE_YEARS, SAMPLE_DAYS)
write.csv(sample_times, "simulations/sample_times.csv", row.names = FALSE)

# Compute sample sizes of males and females
n_males <- sum(samples$sex == "male")
n_females <- sum(samples$sex == "female")
SAMPLE_SIZE = n_males + n_females
## Write other parameters to be used in the SLiM simulation to a file ##
slim_parameters <- tibble(MAX_WIDTH, MAX_HEIGHT, DW, DH, w, h, min_long, max_long, min_lat, max_lat,SAMPLE_SIZE)
#write.csv(slim_parameters, "slim_parameters.csv", row.names = FALSE)

## Convert relationships and samples to the format needed to create input images for the neural network ##

relationships <- read.csv("empirical_data/MLRelate_KNP.csv")

# Convert sample locations to slim (dates are already in slim format)
samples_sf = samples %>% st_as_sf(coords = c("GPS_lat_DD", "GPS_long_DD"), crs = st_crs(kibale))
sample_coord = st_coordinates(samples_sf)
samples <- mutate(samples, slim_x = long_to_x(sample_coord[,1], MAX_WIDTH, w, min_long, max_long),
                  slim_y = lat_to_y(sample_coord[,2], MAX_HEIGHT, h, min_lat, max_lat))

# Check that the transformation looks correct
ggplot() +
  geom_point(data = kibale_grid_slim, aes(x=x, y= y)) +
  geom_point(data = samples, aes(x=slim_x, y=slim_y), color = "red") +
  geom_sf()

# Shorten
short_rel <- select(relationships, "Ind1", "Ind2", "R") %>% filter(R != "U")
short_samples <- select(samples, "individual", "slim_days", "slim_x", "slim_y") %>%
  distinct(individual, slim_days, .keep_all = TRUE) %>%
  rename(x = slim_x, y = slim_y)

# Sampling times and locations of sibling and parent-offspring pairs
# Add info about both individuals
ind1_info <- left_join(short_rel, short_samples, by = join_by(Ind1 == individual), multiple = "all")
ind2_info <- left_join(short_rel, short_samples, by = join_by(Ind2 == individual), multiple = "all")
pair_info <- full_join(ind1_info, ind2_info, by = c("Ind1", "Ind2", "R"), multiple = "all", suffix = c("1", "2"))
sibs <- filter(pair_info, R == "HS" | R == "FS")
pops <- filter(pair_info, R == "PO")
unique(paste(pops$Ind1, pops$Ind2))
filter(pops, Ind1 == "M3_KNP" | Ind2 == "M3_KNP")
expand(pops, nesting(Ind1, Ind2))

# Sampling times and locations of recaptures
recap_ind_df <- group_by(short_samples, individual) %>% summarise(n = n()) %>% filter(n > 1)
recap_inds <- recap_ind_df$individual
rows <- list()
for(ind in recap_inds){
  ind_rows = filter(short_samples, individual == ind) %>% arrange(slim_days)
  # Add all pair of recaptures
  for(i in 1:(nrow(ind_rows)-1)){
    for(j in (i+1):nrow(ind_rows)){
      new_rows <- full_join(ind_rows[i,], ind_rows[j,], by = "individual", suffix = c("1", "2"))
      rows <- append(rows, list(new_rows))
    }
  }
}
recaps <- bind_rows(rows)

# Write siblings, pops, and recaptures to file
#write_csv(sibs, "empirical_data/processed/empirical_sibs_df.csv")
write_csv(pops, "empirical_data/processed/empirical_pops_df.csv")
write_csv(recaps, "empirical_data/processed/empirical_recaps_df.csv")

# Numbers
(nrecaps <- nrow(recaps))
(nunique <- length(unique(samples$individual)))
(nunique_pops <- length(unique(paste(pops$Ind1, pops$Ind2))))
(nrecap_pops <- nrow(pops))
(nunique_sibs <- length(unique(paste(sibs$Ind1, sibs$Ind2))))

# Average distance between recaptures
recap_dist <- mutate(recaps, distance = sqrt((x2-x1)^2 + (y2-y1)^2),
                     days = slim_days2 - slim_days1,
                     distperday = distance/days)
summary(recap_dist)
# Average distance between pops
pop_dist <- mutate(pops, distance = sqrt((x2-x1)^2 + (y2-y1)^2))
summary(pop_dist)
pop_dist_unique <- mutate(pops, distance = sqrt((x2-x1)^2 + (y2-y1)^2)) %>%
  distinct(Ind1, Ind2, .keep_all = TRUE)
summary(pop_dist_unique)

# Sampling intensity

# For each sample, find the closest point in the sampling grid
nearest <- kibale_grid[st_nearest_feature(samples_sf, kibale_grid),]
nearest_coord <- st_coordinates(nearest)
samples <- mutate(samples, grid_x = long_to_x(nearest_coord[,1], MAX_WIDTH, w, min_long, max_long),
                  grid_y = lat_to_y(nearest_coord[,2], MAX_HEIGHT, h, min_lat, max_lat))
ggplot() +
  geom_sf(data = samples_sf, color = "red") +
  geom_sf(data = nearest)

# Count the number of samples at each grid cell
sample_intensity <- select(samples, grid_x, grid_y) %>% group_by(grid_x, grid_y) %>% summarise(nsampled = n()) %>%
  rename(x=grid_x, y = grid_y)
write_csv(sample_intensity, "empirical_data/processed/empirical_intensity_df.csv")
ggplot() +
  geom_point(data = kibale_grid_slim, aes(x=x,y=y)) +
  geom_point(data = sample_intensity, aes(x=x, y = y), color = "red") +
  geom_point(data = samples, aes(x=slim_x, y = slim_y), color = "blue", alpha = 0.5) +
  geom_sf()

# Write empirical grid sampling locations to file
empirical_locations <- select(samples, grid_x, grid_y) %>% rename(x = grid_x, y = grid_y) %>% distinct()
write.csv(empirical_locations, "simulations/empirical_sampling_locations.txt", row.names = FALSE)

# Write empirical grid sampling locations, times, and number of samples to file
sample_intensity_time <- select(samples, grid_x, grid_y, slim_days, slim_years) %>% group_by(grid_x, grid_y, slim_days, slim_years) %>% summarise(nsampled = n()) %>%
  rename(x=grid_x, y = grid_y, SAMPLE_DAYS = slim_days, SAMPLE_YEARS = slim_years)
#write.csv(sample_intensity_time, "simulations/empirical_sampling_locations_times.txt", row.names = FALSE)

# Plot locations of samples and save
ggplot() +
  geom_sf(data = kibale_new, fill = NA, linewidth = 0.5) +
  geom_sf(data = samples_sf, alpha = 0.3, size = 5) +
  coord_sf(expand = FALSE) +
  theme_bw(base_size = 30)

ggsave("paper_figures/kibale_sample_locations.png", width = w, height = h, units = "px", limitsize = FALSE)
