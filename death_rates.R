library(data.table)
library(ggplot2)
library(lubridate)
library(cowplot)
source("./hazard_data.R")


# Load complete data set
cd = complete_data("20210225")

# Assemble data set
dataD = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0,
    date_min = "2020-11-01", keep_missing = TRUE)
dataD = prep_data(dataD)

# Set age X
dataD[, age_x := as.numeric(Age)]
ages = levels(dataD$Age)

dataD[, Age2 := as.character(Age)]
dataD[age < 55, Age2 := "1-54"]
ages2 = unique(dataD[order(age)]$Age2)
dataD[, Age2 := factor(Age2, levels = ages2)]
dataD[, age2_x := as.numeric(Age2)]

# Plot death rates by age, by "what_category", and by SGTF status
death_rates_sgtf = function(what_category, clrs, show_shape_guide, xlab, ylab)
{
    # Deaths by age
    d_standard = dataD[!is.na(sgtf), .(deaths = sum(died), days = sum(time), rate = 10000 * sum(died) / sum(time)), keyby = .(age2_x)]
    d_standard[, rate_lo_frequentist := qchisq(0.025, 2 * deaths) / (2 * days / 10000)]
    d_standard[, rate_hi_frequentist := qchisq(0.975, 2 * deaths + 2) / (2 * days / 10000)]

    # Deaths by "what"
    dc = dataD[!is.na(sgtf)]
    dc[, what := get(..what_category)]
    dc = dc[, .(deaths = sum(died), days = sum(time), rate = 10000 * sum(died) / sum(time)), keyby = .(age2_x, what, sgtf)]
    dc[, rate_lo_frequentist := qchisq(0.025, 2 * deaths) / (2 * days / 10000)]
    dc[, rate_hi_frequentist := qchisq(0.975, 2 * deaths + 2) / (2 * days / 10000)]
    dc[, z_sgtf_label := factor(ifelse(sgtf == 1, "SGTF", "Non-SGTF"), levels = c("Non-SGTF", "SGTF"))]
    
    # Count comparisons for which SGTF is higher than non-SGTF
    higher = sum(dc[sgtf == 1, rate] > dc[sgtf == 0, rate])
    nothigher = sum(dc[sgtf == 1, rate] <= dc[sgtf == 0, rate])
    cat("SGTF crude death rate higher for ", higher, " / ", higher + nothigher, " comparisons.\n", sep = "")

    ggplot() +
        geom_rect(data = d_standard, aes(xmin = age2_x - 0.5, xmax = age2_x + 0.5, 
            ymin = rate_lo_frequentist, ymax = rate_hi_frequentist), fill = "#cccccc") +
        geom_segment(data = d_standard, aes(x = age2_x - 0.5, xend = age2_x + 0.5, 
            y = rate, yend = rate), size = 0.25, linetype = "22") +
        geom_pointrange(data = dc, aes(x = age2_x, ymin = rate_lo_frequentist, y = rate, ymax = rate_hi_frequentist, 
            colour = what, shape = z_sgtf_label), size = 0.4, fatten = 2.4, position = position_dodge(width = 1), stroke = 0.5) +
        scale_x_continuous(breaks = 1:4, labels = ages2, expand = expansion(0.01)) +
        scale_y_log10(breaks = c(0.1, 1, 10, 100), labels = c(0.1, 1, 10, 100), limits = c(0.03, 150), oob = scales::oob_squish) +
        scale_colour_manual(values = clrs) +
        scale_shape_manual(values = c("SGTF" = 5, "Non-SGTF" = 20)) +
        labs(x = if (xlab) "Age" else NULL, 
            y = if (ylab) "Deaths per 10,000\ndays of followup" else NULL, 
            colour = what_category, shape = NULL) +
        theme(legend.position = c(0.015, 1.025), legend.title = element_text(size = 9, margin = margin(0, 0, -0.05, 0, "cm")), 
            legend.text = element_text(size = 9, margin = margin()), legend.justification = c(0, 1),
            legend.key.size = unit(0.225, "cm")) + 
        guides(colour = guide_legend(order = 1, override.aes = list(size = 0.4)), 
            shape = if (show_shape_guide) guide_legend(order = 2, reverse = TRUE, override.aes = list(size = 0.4)) else "none")
}

colours = c("#6388b4", "#fd57ca", "#47a74e", "#eb1e2c", "#9c5142", "#8175aa", "#ccb22b")

pla = death_rates_sgtf("Sex", colours, TRUE, FALSE, TRUE)
plb = death_rates_sgtf("Place of residence", colours, FALSE, FALSE, FALSE)
plc = death_rates_sgtf("Ethnicity", colours, FALSE, FALSE, TRUE)
pld = death_rates_sgtf("IMD decile", colours, FALSE, FALSE, FALSE)
ple = death_rates_sgtf("NHS England region", colours, FALSE, TRUE, TRUE)
plf = death_rates_sgtf("Specimen date", colours, FALSE, TRUE, FALSE)

theme_set(theme_cowplot(font_size = 11))
pl_deathrates_sgtf = plot_grid(pla, plb, plc, pld, ple, plf, 
    nrow = 3, labels = letters[4:9], label_size = 11, rel_widths = c(1.08, 1), rel_heights = c(1, 1, 1.07))
