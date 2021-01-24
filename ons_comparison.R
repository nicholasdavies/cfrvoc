# Analysis of ONS SGTF data -- is there a difference in testing behaviour between SGTF and non-SGTF?
library(readxl)
library(data.table)
library(ggplot2)
library(lubridate)
library(stringr)
library(ogwrangler)
CreateCache()
source("hazard_data.R")

# Load ONS data
ons_r = as.data.table(read_excel("~/Dropbox/uk_covid_data/ons-cis/covid19infectionsurveydatasets20210122.xlsx", "6b", skip = 3))[, 1:8]
names(ons_r) = letters[1:8]
ons_r = ons_r[1:(ons_r[, which(a == "Notes:") - 1])]

regions = ons_r[a %like% "^Percentage", str_remove_all(a, "Percentage and CT Values of COVID-19 cases, ")]
starts = ons_r[, which(a %like% "^Percentage")] + 3;
ends = c(tail(starts, -1) - 6, nrow(ons_r) - 1);

ons_data = NULL;
for (i in seq_along(regions))
{
    ons_sub_data = ons_r[starts[i]:ends[i]];
    ons_sub_data[, region := regions[i]];
    ons_data = rbind(ons_data, ons_sub_data);
}

ons_data[, date := dmy(a)];
ons_data[, N := as.numeric(b)];
ons_data[, O := as.numeric(c)];
ons_data[, S := as.numeric(d)];
ons_data[, NO := as.numeric(e)];
ons_data[, OS := as.numeric(f)];
ons_data[, NS := as.numeric(g)];
ons_data[, NOS := as.numeric(h)];

ons_data = ons_data[, .(date, region, N, O, S, NO, OS, NS, NOS)];

# Load Pillar 2 data
pillar2 = sgtf_counts("20210122", c("LTLA_code"))
pillar2 = pillar2[group != "" & specimen_date > "2020-09-01"]
pillar2[group %in% c("E07000004", "E07000005", "E07000006", "E07000007"), group := "E06000060"];
pillar2[, rgn := ogwhat(group, "rgn")]
pillar2 = pillar2[, .(sgtf = sum(sgtf), other = sum(other)), keyby = .(specimen_date, rgn)]
pillar2[, region := ogwhat(rgn)]
pillar2[, sgtf_ct := sgtf / (sgtf + other)]
pillar2[, sgtf_lo := qbeta(0.025, shape1 = sgtf, shape2 = other)]
pillar2[, sgtf_hi := qbeta(0.975, shape1 = sgtf, shape2 = other)]

ggplot() +
    geom_ribbon(data = pillar2, aes(x = specimen_date, ymin = sgtf_lo, ymax = sgtf_hi, fill = "Pillar 2"), alpha = 0.4) +
    geom_line(data = pillar2, aes(x = specimen_date, y = sgtf_ct, colour = "Pillar 2")) +
    geom_line(data = ons_data, aes(x = date + 3.5, y = NO / (NO + NOS), colour = "ONS CIS 1")) +
    geom_line(data = ons_data, aes(x = date + 3.5, y = NO / (NO + NOS + N + O + S + OS + NS), colour = "ONS CIS 2")) +
    scale_colour_manual(aesthetics = c("fill", "colour"), values = c("Pillar 2" = "darkorchid", "ONS CIS 1" = "darkorange", "ONS CIS 2" = "darkgreen")) +
    facet_wrap(~region) +
    labs(x = "Date", y = "Proportion with S-gene dropout", fill = "Source", colour = "Source")


cd = complete_data("20210122")
ct = function(x) as.numeric(ifelse(x == 0, 40, x))
cts = cd[!is.na(specimen_date.x), .(sgtf, sgtf_under30CT, ctORF1ab = ct(P2CH1CQ), ctN = ct(P2CH2CQ), ctS = ct(P2CH3CQ), ctControl = ct(P2CH4CQ)), 
    keyby = .(specimen_date = specimen_date.x, NHSER_name, UTLA_name, LTLA_name)]

ggplot(cts) +
    geom_density2d(aes(ctORF1ab, ctS)) +
    geom_abline(slope = 1, intercept = 0)

ggplot(cts) +
    geom_density2d(aes(ctN, ctS)) +
    geom_abline(slope = 1, intercept = 0)

cts = melt(cts, id.vars = 1:6)

ggplot(cts[!is.na(value) & value < 40]) +
    geom_density(aes(x = value, colour = variable))



ons_data
pillar2

hist(ons_data$NO)
hist(ons_data$OS)
hist(ons_data$NS)

ons_data
starts
ons_r
View(ons_r)

historic = as.data.table(read_excel("~/Dropbox/uk_covid_data/ons-cis/covid19infectionsurveydatasets20210122.xlsx", "8b"));
historic = historic[1:1159, 1:4]
names(historic) = letters[1:4]

date_rows = historic[, which(a == "Date")]
start_rows = date_rows + 2
end_rows = c(tail(date_rows - 5, -1), 1159)
date_ids = historic[date_rows - 3, a]

hist = NULL
for (i in seq_along(start_rows))
{
    hist_subset = historic[start_rows[i]:end_rows[i]]
    hist_subset[, publication := date_ids[i]]
    hist = rbind(hist, hist_subset)
}

names(hist) = c("date", "mid", "lo", "hi", "publication_date")
hist[, date2 := ymd("1899-12-30") + as.numeric(date)]
hist[, mid := as.numeric(mid)]
hist[, lo := as.numeric(lo)]
hist[, hi := as.numeric(hi)]
hist[, set := dmy(str_remove_all(publication_date, "Publication date "))]
hist[, offset := set - min(set)]
hist[, offset := as.numeric(offset) / max(as.numeric(offset))]

ggplot(hist) +
    geom_ribbon(aes(x = date2, ymin = lo * 100, ymax = hi * 100, fill = as.factor(set)), alpha = 0.4, colour = "black", size = 0.125) +
    scale_fill_hue(aesthetics = c("colour", "fill")) +
    scale_x_date(date_breaks = "1 month", date_labels = "%B") +
    annotate("text", x = ymd("2021-01-10"), y = 0.01, label = "Data source: Office for National Statistics", hjust = 1, size = 3.5, fontface = "bold") +
    labs(y = "SARS-CoV-2 prevalence in England, %", x = NULL)

hist
end_rows
View(historic)

# Load new variant compatible modelled data for countries of the UK and for regions of England
ons_c = as.data.table(read_excel("~/Dropbox/uk_covid_data/ons-cis/covid19infectionsurveydatasets202101082.xlsx", "6c", "A5:Y49", na = "N/A"));
ons_r = as.data.table(read_excel("~/Dropbox/uk_covid_data/ons-cis/covid19infectionsurveydatasets202101082.xlsx", "6d", "A5:BC49", na = "N/A"));

# Reshape ONS data into long format
reformat_ons = function(ons_c)
{
    regions = names(ons_c);
    regions = regions[!regions %like% "\\.\\.\\."];
    names(ons_c)[1] = "date";
    
    dates = dmy(ons_c[3:.N, date]);
    output = NULL;
    
    for (r in seq_along(regions))
    {
        col_start = r * 6 - 4;
        col_end = r * 6 + 1;
        ons_sub = ons_c[3:.N, lapply(.SD, as.numeric), .SDcols = col_start:col_end];
        names(ons_sub) = c("mean.sgtf", "lo.sgtf", "hi.sgtf", "mean.other", "lo.other", "hi.other");
        ons_sub[, geo := regions[r]];
        
        output = rbind(output,
            cbind(date = dates, ons_sub))
    }
    
    output[, geo := factor(geo, levels = regions)];
    
    return (output[])
}

approx_quotient = function(mu.a, lo.a, hi.a, mu.b, lo.b, hi.b)
{
    sd.a = (hi.a - lo.a) / (1.96 * 2)
    sd.b = (hi.b - lo.b) / (1.96 * 2)

    var.c = (mu.a^2 / mu.b^2) * (sd.a^2 / mu.a^2 + sd.b^2 / mu.b^2)
    sd.c = sqrt(var.c)
    mu.c = mu.a / mu.b
    
    return (list(mean.q = mu.c, lo.q = mu.c - 1.96 * sd.c, hi.q = mu.c + 1.96 * sd.c))
}

approx_frac = function(mu.a, lo.a, hi.a, mu.b, lo.b, hi.b)
{
    mu.ab = mu.a + mu.b
    
    sd.a = (hi.a - lo.a) / (1.96 * 2)
    sd.b = (hi.b - lo.b) / (1.96 * 2)
    
    var.ab = sd.a^2 + sd.b^2
    sd.ab = sqrt(var.ab)
    
    approx_quotient(mu.a, lo.a, hi.a, mu.ab, mu.ab - 1.96 * sd.ab, mu.ab + 1.96 * sd.ab)
}

mu.c
sd.c

#w = reformat_ons(ons_c)
w = reformat_ons(ons_r)

w[, c("mean.frac", "lo.frac", "hi.frac") := approx_frac(mean.sgtf, lo.sgtf, hi.sgtf, mean.other, lo.other, hi.other)]
w

ggplot(w) +
    geom_pointrange(aes(x = date, y = mean.ratio, ymin = lo.ratio, ymax = hi.ratio, colour = "ratio"), fatten = 0.1) +
    geom_pointrange(aes(x = date, y = mean.sgtf,  ymin = lo.sgtf,  ymax = hi.sgtf,  colour = "sgtf"), fatten = 0.1) +
    geom_pointrange(aes(x = date, y = mean.other, ymin = lo.other, ymax = hi.other, colour = "other"), fatten = 0.1) +
    facet_wrap(~geo)

# p from phe_data.R
p = p[specimen_date >= "2020-11-22"]
p[, sgtf_frac := sgtf / (other + sgtf)]
p0 = p0[specimen_date >= "2020-11-22"]
p0[, sgtf_frac := sgtf / (other + sgtf)]

w[, point.sgtf := mean.sgtf / (mean.other + mean.sgtf)]

ggplot(w) +
    geom_pointrange(aes(x = date, y = mean.frac, ymin = lo.frac, ymax = hi.frac, colour = "ONS"), fatten = 0.1) +
    geom_line(data = p, aes(x = specimen_date, y = sgtf_frac, colour = "Pillar 2")) +
    facet_wrap(~geo) +
    labs(y = "Proportion SGTF") +
    theme(panel.background = element_rect(fill = "#f4f4f4"))


lapply(w, class)
ons_c
