library(data.table)
library(ggplot2)
library(lubridate)
library(cowplot)
source("./hazard_data.R")


# Load complete data set
cd = complete_data("20210122")

# Assemble data set
dataS = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0,
    date_min = "2020-11-01", keep_missing = TRUE)
dataS = prep_data(dataS)

# Set age X
dataS[, age_x := as.numeric(Age)]
ages = levels(dataS$Age)

dataS[, Age2 := as.character(Age)]
dataS[age < 55, Age2 := "1–54"]
ages2 = unique(dataS[order(age)]$Age2)
dataS[, Age2 := factor(Age2, levels = ages2)]
dataS[, age2_x := as.numeric(Age2)]

# Plot deaths rates by age and by "what"
death_rates = function(what, region = "")
{
    # Deaths by age
    if (region == "") {
        d_standard = dataS[, .(deaths = sum(died), days = sum(time), rate = 10000 * sum(died) / sum(time)), keyby = .(age_x)]
    } else {
        d_standard = dataS[NHSER_name == region, .(deaths = sum(died), days = sum(time), rate = 10000 * sum(died) / sum(time)), keyby = .(age_x)]
    }
    d_standard[, rate_lo_bayesian := qgamma(0.025, shape = deaths, rate = days / 10000)]
    d_standard[, rate_hi_bayesian := qgamma(0.975, shape = deaths, rate = days / 10000)]
    d_standard[, rate_lo_frequentist := qchisq(0.025, 2 * deaths) / (2 * days / 10000)]
    d_standard[, rate_hi_frequentist := qchisq(0.975, 2 * deaths + 2) / (2 * days / 10000)]

    # Deaths by "what"
    if (region == "") {
        dc = dataS[, .(deaths = sum(died), days = sum(time), rate = 10000 * sum(died) / sum(time)), keyby = .(age_x, what = get(what))]
    } else {
        dc = dataS[NHSER_name == region]
        dc[, what := get(..what)]
        dc = dc[, .(deaths = sum(died), days = sum(time), rate = 10000 * sum(died) / sum(time)), keyby = .(age_x, what)]
    }
    dc[, rate_lo_bayesian := qgamma(0.025, shape = deaths, rate = days / 10000)]
    dc[, rate_hi_bayesian := qgamma(0.975, shape = deaths, rate = days / 10000)]
    dc[, rate_lo_frequentist := qchisq(0.025, 2 * deaths) / (2 * days / 10000)]
    dc[, rate_hi_frequentist := qchisq(0.975, 2 * deaths + 2) / (2 * days / 10000)]

    ggplot() +
        geom_rect(data = d_standard, aes(xmin = age_x - 0.5, xmax = age_x + 0.5, 
            ymin = rate_lo_frequentist, ymax = rate_hi_frequentist), fill = "#cccccc") +
        geom_segment(data = d_standard, aes(x = age_x - 0.5, xend = age_x + 0.5, 
            y = rate, yend = rate), size = 0.25, linetype = "22") +
        geom_pointrange(data = dc, aes(x = age_x, ymin = rate_lo_frequentist, y = rate, ymax = rate_hi_frequentist, 
            colour = what), shape = "–", fatten = 8, position = position_dodge(width = 0.8)) +
        scale_x_continuous(breaks = 1:5, labels = ages) +
        scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = c(0.01, 0.1, 1, 10, 100)) +
        labs(x = "Age", y = "Deaths per 10,000 days of followup", colour = what) +
        theme(legend.position = c(0.025, 0.975), legend.justification = c(0, 1))
}

# Plot death rates by age, by "what_category", and by SGTF status
death_rates_sgtf = function(what_category, region = "")
{
    # Deaths by age
    if (region == "") {
        d_standard = dataS[!is.na(sgtf), .(deaths = sum(died), days = sum(time), rate = 10000 * sum(died) / sum(time)), keyby = .(age2_x)]
    } else {
        d_standard = dataS[!is.na(sgtf) & NHSER_name == region, .(deaths = sum(died), days = sum(time), rate = 10000 * sum(died) / sum(time)), keyby = .(age2_x)]
    }
    d_standard[, rate_lo_bayesian := qgamma(0.025, shape = deaths, rate = days / 10000)]
    d_standard[, rate_hi_bayesian := qgamma(0.975, shape = deaths, rate = days / 10000)]
    d_standard[, rate_lo_frequentist := qchisq(0.025, 2 * deaths) / (2 * days / 10000)]
    d_standard[, rate_hi_frequentist := qchisq(0.975, 2 * deaths + 2) / (2 * days / 10000)]

    # Deaths by "what"
    if (region == "") {
        dc = dataS[!is.na(sgtf)]
        dc[, what := get(..what_category)]
        dc = dc[, .(deaths = sum(died), days = sum(time), rate = 10000 * sum(died) / sum(time)), keyby = .(age2_x, what, sgtf)]
    } else {
        dc = dataS[!is.na(sgtf) & NHSER_name == region]
        dc[, what := get(..what_category)]
        dc = dc[, .(deaths = sum(died), days = sum(time), rate = 10000 * sum(died) / sum(time)), keyby = .(age2_x, what, sgtf)]
    }
    dc[, rate_lo_bayesian := qgamma(0.025, shape = deaths, rate = days / 10000)]
    dc[, rate_hi_bayesian := qgamma(0.975, shape = deaths, rate = days / 10000)]
    dc[, rate_lo_frequentist := qchisq(0.025, 2 * deaths) / (2 * days / 10000)]
    dc[, rate_hi_frequentist := qchisq(0.975, 2 * deaths + 2) / (2 * days / 10000)]
    dc[, z_sgtf_label := ifelse(sgtf == 1, "SGTF", "Non-SGTF")]

    ggplot() +
        geom_rect(data = d_standard, aes(xmin = age2_x - 0.5, xmax = age2_x + 0.5, 
            ymin = rate_lo_frequentist, ymax = rate_hi_frequentist), fill = "#cccccc") +
        geom_segment(data = d_standard, aes(x = age2_x - 0.5, xend = age2_x + 0.5, 
            y = rate, yend = rate), size = 0.25, linetype = "22") +
        geom_pointrange(data = dc, aes(x = age2_x, ymin = rate_lo_frequentist, y = rate, ymax = rate_hi_frequentist, 
            colour = what, shape = z_sgtf_label), fatten = 4, position = position_dodge(width = 0.8)) +
        scale_x_continuous(breaks = 1:4, labels = ages2) +
        scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = c(0.01, 0.1, 1, 10, 100)) +
        labs(x = "Age", y = "Deaths per 10,000 days of followup", colour = what_category, shape = NULL) +
        theme(legend.position = c(0.025, 0.975), legend.justification = c(0, 1)) + 
        guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))
}

pla = death_rates("Sex")
plb = death_rates("Place of residence")
plc = death_rates("Ethnicity")
pld = death_rates("Index of Multiple Deprivation decile")
ple = death_rates("NHS England region")
plf = death_rates("Specimen date")

theme_set(theme_cowplot(font_size = 10))
pl = plot_grid(pla, plb, plc, pld, ple, plf, nrow = 3, labels = letters, label_size = 10)
ggsave("./output/death_rates.png", width = 30, height = 34, units = "cm")


pla = death_rates_sgtf("Sex")
plb = death_rates_sgtf("Place of residence")
plc = death_rates_sgtf("Ethnicity")
pld = death_rates_sgtf("Index of Multiple Deprivation decile")
ple = death_rates_sgtf("NHS England region")
plf = death_rates_sgtf("Specimen date")

theme_set(theme_cowplot(font_size = 10))
pl = plot_grid(pla, plb, plc, pld, ple, plf, nrow = 3, labels = letters, label_size = 10)
ggsave("./output/death_rates_sgtf.png", width = 30, height = 34, units = "cm")

