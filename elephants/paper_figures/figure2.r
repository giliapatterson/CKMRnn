library(tidyverse)
library(sf)
library(lubridate)
library(magick)

# Read in protected areas in Uganda
pas <- st_read("../empirical_data/Ug_Protected-areas2007/Ug_Protected-areas2007.shp")

# Remove all but Kibale and Kisangi (DJM)
kibale <- filter(pas, NAME == "Kibale", GAZTYPE == "NP")
kisangi <- filter(pas, NAME == "Kisangi", GAZTYPE == "DJM")
kisangi2 <- filter(pas, NAME == "Kisangi", GAZTYPE == "CFR")
kibale_new <- st_simplify(st_union(rbind(kibale, kisangi)), preserveTopology = TRUE, dTolerance = 0.01)
# Limits of the map
min_long = st_bbox(kibale_new)$xmin
max_long = st_bbox(kibale_new)$xmax

# Empirical data
samples <- read_csv("../empirical_data/Kibale_Species_Classification_Final.csv")
samples_sf = samples %>% st_as_sf(coords = c("GPS_lat_DD", "GPS_long_DD"), crs = st_crs(kibale))

# Plot locations of samples
samples <- ggplot() +
  geom_sf(data = kibale_new, fill = NA, linewidth = 0.5) +
  geom_sf(data = samples_sf, alpha = 0.3, size = 1) +
  coord_sf(expand = TRUE) +
  theme_bw(base_size = 14) +
  scale_x_continuous(breaks = c(30.20, 30.37, 30.54))

# Empirical input images
intensity <-  ggdraw() +
  draw_image("../network/empirical_input_images/empirical_intensity.png", scale = 0.8)
recaps <- ggdraw() +
  draw_image("../network/empirical_input_images/empirical_recaps.png", scale = 0.8)
pops <- ggdraw() +
  draw_image("../network/empirical_input_images/empirical_pops.png", scale = 0.8)

samples_row <- plot_grid(NULL, samples, NULL, ncol = 3, labels = c("", "(a)", ""),
                         rel_widths = c(1, 1.5, 1), align = "hv",
                         label_x = 0, label_y = 0,
                         hjust = -0.5, vjust = 0.5)
images_row <- plot_grid(intensity, recaps, pops, ncol = 3,
                        labels = c("(b)", "(c)", "(d)"),
                        label_x = 0, label_y = 0,
                        hjust = 0, vjust = -0.9)
input <- plot_grid(samples_row, images_row, ncol = 1, rel_heights = c(1.2, 1))
ggsave("../paper_figures/figure2.pdf", input, width = 6, height = 6)
