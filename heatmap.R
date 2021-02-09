library(maps)
library(mapdata)
library(maptools)
library(rgdal)
library(ggmap)
library(ggplot2)
library(rgeos)
library(broom)
library(plyr)
library(pals)


shapefile = readOGR(dsn = "~/Documents/uk_covid_data_sensitive/Local_Authority_Districts_(December_2019)_Boundaries_UK_BGC-shp/", 
    layer = "Local_Authority_Districts_(December_2019)_Boundaries_UK_BGC")

mapdata = tidy(shapefile, region = "lad19cd")
setDT(mapdata)
mapdata = mapdata[!id %like% "^N"]

missingness = cd[!is.na(specimen_date.x) & specimen_date.x > "2020-11-01" & LTLA_code != "" & pillar == "Pillar 2", 
    .(missingness = mean(is.na(sgtf_under30CT))), by = LTLA_code]
missingness[missingness == 1] # NOTE: exclude Isles of Scilly

mapdata = merge(mapdata, missingness, by.x = "id", by.y = "LTLA_code", all = TRUE)

llabs = fread(
"lab	lat	long	E	N	sgtf_capable
Milton Keynes	52.0406	-0.7594	485189	238748	TRUE
Alderley Park	53.271	-2.233	384558	374916	TRUE
Glasgow	55.8642	-4.2518	259176	665740	TRUE
Cambridge	52.2053	0.1218	545088	258463	FALSE
Newport	51.5842	-2.9977	330970	187732	FALSE
Charnwood	52.7407	-1.1451	457814	316239	FALSE")


# https://github.com/wch/ggplot2/wiki/New-theme-system
new_theme_empty <- theme_bw()
new_theme_empty$line <- element_blank()
new_theme_empty$rect <- element_blank()
new_theme_empty$strip.text <- element_blank()
new_theme_empty$axis.text <- element_blank()
new_theme_empty$plot.title <- element_blank()
new_theme_empty$axis.title <- element_blank()
new_theme_empty$plot.margin <- structure(c(0, 0, -1, -1), unit = "lines", valid.unit = 3L, class = "unit")

g = ggplot(mapdata) + 
    geom_polygon(aes(x = long, y = lat, group = group, fill = missingness), colour = NA) +
    geom_point(data = llabs, aes(x = E, y = N, colour = sgtf_capable), size = 1) +
    geom_text(data = llabs, aes(x = E + 10000, y = N - 10000, label = lab), colour = "#88cc88", size = 2.5, hjust = 0) +
    coord_fixed(1) +
    scale_fill_gradientn(colours = rev(ocean.solar(20)), limits = c(0, 1)) +
    new_theme_empty

ggsave("./output/missing_plot.png", g, width = 15, height = 15, units = "cm")

