# *******************************************
# This file contains code for making maps of
# the WHO super-regions and regions based on
# clustering by PM2.5.
# *******************************************

# Setup -------------------------------------------------------------------
library(sp)
library(spdep)
library(mapmisc)
library(maptools)
library(rworldmap)
library(scales)
library(RColorBrewer)
library(dplyr)

# load 'GM' SpatialPointsDataFrame
load("bayes-vis.RData")


# Map of WHO super-regions ----------------------------------------------
spdf <- joinCountryData2Map(GM@data, nameJoinColumn="iso3")

# code missings
spdf$super_region[is.na(spdf$super_region) & spdf$continent == "Eurasia"] <- 4
spdf$super_region[is.na(spdf$super_region) & spdf$continent == "South America"] <- 5
spdf$super_region[is.na(spdf$super_region) & spdf$continent == "Africa"] <- 7
spdf$super_region[spdf$NAME %in% c("Algeria", "Libya", "Yemen", "Syria", "Iraq")] <- 2
spdf$super_region[spdf$NAME %in% c("Taiwan", "Vietnam", "Cambodia", "Laos", "Papua New Guinea", "N. Korea", "Sri Lanka")] <- 6
spdf$super_region[spdf$NAME %in% "Nepal"] <- 3
spdf$super_region[spdf$NAME %in% c("Greenland", "W. Sahara")] <- NA

opar <- par()
png(filename = "plots/map-who.png", res = 500, width = 1000, height = 600)
par(mar=c(0,0,0,0))
mapCountryData(
  spdf,
  nameColumnToPlot = "super_region",
  catMethod = "categorical",
  colourPalette = rev(
    scales::hue_pal(
      h = c(0, 360) + 15,
      c = 100,
      l = 65,
      h.start = 0,
      direction = 1
    )(7)
  ),
  addLegend = FALSE,
  mapTitle = "",
  lwd = 0.25
)
dev.off()
par(opar)



# Map of regions obtained from clustering ---------------------------------
average <-
  GM@data %>%
  group_by(iso3) %>%
  summarise(pm25 = mean(pm25))
d <- dist(average)
hh <- hclust(d)
clust <- cutree(hh,k = 6)
GM@data$cluster_region <-
  sapply(GM@data$iso3, function(x) clust[which(average$iso3 == x)])

spdf2 <- joinCountryData2Map(GM@data, nameJoinColumn="iso3")
spdf2$cluster_region[is.na(spdf2$cluster_region) & !is.na(spdf2$ISO3) & spdf$ISO3 != "ATA" & spdf2$ISO3 != "GRL"] <- 7

png(filename = "plots/map-clusters.png", res = 500, width = 1000, height = 600)
par(mar=c(0,0,0,0))
maplist <- mapCountryData(
  spdf2,
  nameColumnToPlot = "cluster_region",
  catMethod = "fixedWidth",
  colourPalette = rev(
    scales::hue_pal(
      h = c(0, 360) + 15,
      c = 100,
      l = 65,
      h.start = 0,
      direction = 1
    )(7)
  ),
  addLegend = FALSE,
  mapTitle = "",
  lwd = 0.25
)
dev.off()
graphics.off()
