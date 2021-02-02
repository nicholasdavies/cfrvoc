library(data.table)
library(ggplot2)

source("./phe_data.R")
source("./hazard_data.R")

rates = cd[pillar == "Pillar 2" & !is.na(sgtf_under30CT), .(sgtf = sum(sgtf_under30CT == 1), non = sum(sgtf_under30CT == 0)), keyby = .(specimen_date.x, NHSER_name)]

ggplot(rates[specimen_date.x >= "2020-12-01"]) + 
    geom_line(aes(x = specimen_date.x, y = sgtf / (sgtf + non), colour = "sgtf")) +
    facet_wrap(~NHSER_name)

pos = phe_positives("20210201")
pos[, specimen_date := dmy(specimen_date)]
pos = pos[NHSER_name != "", .(p1 = sum(pillar == "Pillar 1"), p2 = sum(pillar == "Pillar 2")), keyby = .(specimen_date, NHSER_name)]
pos[, f := zoo::rollmean(p1 / (p1 + p2), 7, fill = "extend"), by = NHSER_name]

ggplot(pos[specimen_date >= "2020-11-01"]) +
    geom_line(aes(x = specimen_date, y = f, colour = NHSER_name))

death = phe_deaths("20210201")

death[, table(sign(finalid), pillars, useNA = "ifany")]

w = death[, .(p1 = sum(pillars == "PILLAR1" | pillars == "PILLAR1 AND PILLAR2", na.rm = T), 
              p2 = sum(pillars == "PILLAR2" | pillars == "PILLAR1 AND PILLAR2", na.rm = T)), keyby = .(nhser_name, dod)]

w[, f := zoo::rollmean(p1 / (p1 + p2), 7, fill = "extend"), by = nhser_name]

ggplot(w[dod >= "2020-11-01"]) +
    geom_line(aes(x = dod, y = f, colour = nhser_name))
