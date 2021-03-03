# Analysis of ONS SGTF data -- is there a difference in testing behaviour between SGTF and non-SGTF?
library(readxl)
library(data.table)
library(ggplot2)
library(lubridate)
library(stringr)
library(ogwrangler)
library(writexl)
library(cowplot)
CreateCache()

source("./phe_data.R")
source("./hazard_data.R")

# Get UK population
ukpopf = data.table(read_excel("./data/ukmidyearestimates20192020ladcodes.xlsx", "MYE2 - Females", "A5:CQ431"))
ukpopm = data.table(read_excel("./data/ukmidyearestimates20192020ladcodes.xlsx", "MYE2 - Males", "A5:CQ431"))
ukpopf = melt(ukpopf[Code %like% "^E06|^E07|^E08|^E09"], id.vars = 1:3, variable.name = "age", value.name = "population")[age != "All ages"]
ukpopm = melt(ukpopm[Code %like% "^E06|^E07|^E08|^E09"], id.vars = 1:3, variable.name = "age", value.name = "population")[age != "All ages"]
ukpopf[age == "90+", age := "90"]
ukpopm[age == "90+", age := "90"]
ukpopf[, age := as.character((as.numeric(as.character(age)) %/% 5) * 5)]
ukpopm[, age := as.character((as.numeric(as.character(age)) %/% 5) * 5)]
ukpopf = ukpopf[, .(pop = sum(population)), by = .(Code, age)]
ukpopm = ukpopm[, .(pop = sum(population)), by = .(Code, age)]
ukpopf[, lad := ogwhat(Code)]
ukpopm[, lad := ogwhat(Code)]
ukpopf[, sex := "Female"]
ukpopm[, sex := "Male"]
ukpop = rbind(ukpopf, ukpopm)
ukpop = ukpop[order(lad, age, sex)]

# Load Pillar 2 data
pillar2 = sgtf_counts("20210225", c("age5", "sex", "LTLA_code", "cat"), criterion = "ONS")
pillar2[, c("age", "sex", "LTLA_code", "cat") := tstrsplit(group, "\\|")]

# Keep only residential accommodation, samples from Nov 30 on, samples with information
pillar2 = pillar2[cat %in% c("Residential dwelling (including houses, flats, sheltered accommodation)", "House in multiple occupancy (HMO)")]
pillar2 = pillar2[!is.na(age)]
pillar2 = pillar2[sex != "Unknown"]
pillar2 = pillar2[LTLA_code != ""]
pillar2 = pillar2[specimen_date >= "2020-11-30"]
pillar2[, specimen_week := floor_date(specimen_date, "week", week_start = 1)]

# Amalgamate by English administrative region
pillar2[LTLA_code %in% c("E07000004", "E07000005", "E07000006", "E07000007"), LTLA_code := "E06000060"];
pillar2 = pillar2[, .(sgtf = sum(sgtf), other = sum(other)), keyby = .(specimen_date, LTLA_code, age, sex)]
pillar2 = merge(pillar2, ukpop, by.x = c("LTLA_code", "age", "sex"), by.y = c("Code", "age", "sex"), all.x = TRUE)
pillar2[, rgn := ogwhat(LTLA_code, "rgn")]
pillar2[, region := ogwhat(rgn)]
pillar2 = pillar2[!is.na(pop)]

pillar2 = pillar2[, .(sgtf_ct = weighted.mean(sgtf / (sgtf + other), pop, na.rm = T)), keyby = .(geo = region, specimen_date)]

# Load new variant compatible modelled data for regions of England
ons_r = as.data.table(read_excel("./data/covid19infectionsurveydatasets2021012928012021212458.xlsx", "6d", "A5:CD49", na = "N/A"));

# Reshape ONS data into long format
reformat_ons = function(ons_c)
{
    regions = names(ons_c)[-1];
    regions = regions[!regions %like% "\\.\\.\\."];
    names(ons_c)[1] = "date";
    
    if ("POSIXct" %in% class(ons_c[3, date])) {
        dates = ymd(ons_c[3:.N, date]);
    } else {
        dates = ymd("1899-12-30") + ons_c[3:.N, as.numeric(date)];
    }
    output = NULL;
    
    for (r in seq_along(regions))
    {
        col_start = r * 9 - 7;
        col_end = r * 9 - 2;
        ons_sub = ons_c[3:.N, lapply(.SD, as.numeric), .SDcols = col_start:col_end];
        names(ons_sub) = c("mean.sgtf", "lo.sgtf", "hi.sgtf", "mean.other", "lo.other", "hi.other");
        ons_sub[, geo := regions[r]];
        
        output = rbind(output,
            cbind(date = dates, ons_sub))
    }
    
    output[, geo := factor(geo, levels = regions)];
    
    return (output[])
}

approx_frac = function(mu.a, lo.a, hi.a, mu.b, lo.b, hi.b)
{
    sd.a = (hi.a - lo.a) / (1.96 * 2)
    sd.b = (hi.b - lo.b) / (1.96 * 2)
    
    a = pmax(0.001, rnorm(1000, mu.a, sd.a));
    b = pmax(0.001, rnorm(1000, mu.b, sd.b));
    f = a / (a + b);
    
    return (list(mean.f = mean(f), lo.f = quantile(f, 0.025), hi.f = quantile(f, 0.975)))
}


w = reformat_ons(ons_r)
w[, c("mean.frac", "lo.frac", "hi.frac") := approx_frac(mean.sgtf, lo.sgtf, hi.sgtf, mean.other, lo.other, hi.other), by = .(date, geo)]

theme_set(theme_cowplot(font_size = 11) + theme(strip.background = element_blank()))

ggplot(w) +
    geom_pointrange(aes(x = date, y = mean.frac,  ymin = lo.frac, ymax = hi.frac, colour = "ONS"), fatten = 0.1) +
    geom_line(data = pillar2, aes(x = specimen_date, y = sgtf_ct, colour = "Pillar 2")) +
    labs(colour = NULL, x = "Date", y = "Fraction ORF1ab + N : ORF1ab + N + S") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~geo)

ww = merge(w[, .(specimen_date = date, mean.ons = mean.frac, lo.ons = lo.frac, hi.ons = hi.frac, geo)], 
    pillar2[, .(geo, specimen_date, mean.pillar2 = sgtf_ct)], by = c("geo", "specimen_date"), all = TRUE)

ggsave("./output/ons_comparison.pdf", width = 20, height = 15, units = "cm", useDingbats = FALSE)
write_xlsx(ww, "./manuscript/sdE_ons_comparison.xlsx")
