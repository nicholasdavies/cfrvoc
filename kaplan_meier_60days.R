## km plots paper
# Kaplan-Meier plots

library(ggplot2)
library(data.table)
library(KMunicate)
library(survival)
library(cowplot)

source("./kmunicate2.R")
source("./phe_data.R")
source("./hazard_data.R")

theme_set(theme_cowplot(font_size = 10))

# Load complete data set
cd = complete_data("20210205")

# 60 DAY VIEW
# Assemble data set
dataS60 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 60, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS60[, sgtf_label := ifelse(sgtf == 0, "Non-SGTF", "SGTF")]
dataS60[, sgtf_label := factor(sgtf_label, c("SGTF", "Non-SGTF"))]

# Summary view
tsc = seq(0, 60, by = 10)
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60)
plAA60 = KMunicate2(fit = km, time_scale = tsc, .ylim = c(0.995, 1), .margin = c(0, 0.1, 0, 0.1), .title = "Overall", .legend_position = c(0.05, 0.1), .risk_table_base_size = 9)

plAA60_inset = KMunicate2(fit = km, time_scale = c(0, 60), .ylim = c(0, 1), .margin = 0.3, .title = "", .risk_table = NULL, .legend_position = "none",.ylab = "", .xlab = "") + scale_y_continuous(breaks = c(0,1), expand = expansion(0)) 

plAA60_both <- ggdraw() + draw_plot(plAA60) + draw_plot(plAA60_inset, x = .5, y = .6, width = .5, height = .4)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[age < 70])
plBB60 = KMunicate2(fit = km, time_scale = tsc, .ylim = c(0.9975, 1), .margin = c(0, 0.1, 0, 0.1), .title = "Under 70", .legend_position = c(0.05, 0.1), .risk_table_base_size = 9)

plBB60_inset = KMunicate2(fit = km, time_scale = c(0, 60), .ylim = c(0, 1), .margin = 0.3, .title = "", .risk_table = NULL, .legend_position = "none",.ylab = "", .xlab = "") + scale_y_continuous(breaks = c(0,1), expand = expansion(0)) 

plBB60_both <- ggdraw() + draw_plot(plBB60) + draw_plot(plBB60_inset, x = 0.6, y = .7, width = .4, height = .3)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[age >= 70])
plCC60 = KMunicate2(fit = km, time_scale = tsc, .ylim = c(0.92, 1), .margin = c(0, 0.1, 0, 0.1), .title = "70 or older", .legend_position = c(0.05, 0.1), .risk_table_base_size = 9)

plCC60_inset = KMunicate2(fit = km, time_scale = c(0, 60), .ylim = c(0, 1), .margin = 0.3, .title = "", .risk_table = NULL, .legend_position = "none",.ylab = "", .xlab = "") + scale_y_continuous(breaks = c(0,1), expand = expansion(0)) 

plCC60_both <- ggdraw() + draw_plot(plCC60) + draw_plot(plCC60_inset, x = 0.6, y = .7, width = .4, height = .3)

pl = cowplot::plot_grid(plAA60, plBB60, plCC60, nrow = 1, labels = letters, label_size = 10)
ggsave("./output/kmcurves_60_main.pdf", pl, width = 45, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves_60_main.png", pl, width = 45, height = 15, units = "cm")

pl60_inset = cowplot::plot_grid(plAA60_both, plBB60_both, plCC60_both, nrow = 1, labels = letters, label_size = 10)
ggsave("./output/kmcurves_60_main_inset.pdf", pl_60_inset, width = 45, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves_60_main_inset.png", pl60_inset, width = 45, height = 15, units = "cm")



# By Sex ------------------------------------------------------------------


km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[sex == "Female"])
plA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Female", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[sex == "Male"])
plB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Male", .risk_table = NULL, .legend_position = c(0.05, 0.1))


pl = cowplot::plot_grid(plA, plB, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves_60_sex.pdf", pl, width = 10, height = 5, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves_60_sex.png", pl, width = 10, height = 5, units = "cm")


# By Age ------------------------------------------------------------------


km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[age_group == "(0,35]"])
plA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.999, 1), .margin = 0.3, .title = "0-34", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[age_group  == "(35,55]"])
plB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.996, 1), .margin = 0.3, .title = "35-54", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[age_group == "(55,70]"])
plC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.98, 1), .margin = 0.3, .title = "55-69", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[age_group  %in% c("(70,85]", "(85,120]")])
plD = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.89, 1), .margin = 0.3, .title = "70+", .risk_table = NULL, .legend_position = c(0.05, 0.1))


pl = cowplot::plot_grid(plA, plB, plC, plD, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves_60_age.pdf", pl, width = 10, height = 10, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves_60_age.png", pl, width = 10, height = 10, units = "cm")


# Place of residence ------------------------------------------------------

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[res_cat == "Residential"])
plA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.995, 1), .margin = 0.3, .title = "Residential", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[res_cat == "Care/Nursing home"])
plB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.75, 1), .margin = 0.3, .title = "Care/Nursing home", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[res_cat == "Other/Unknown"])
plC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.992, 1), .margin = 0.3, .title = "Other/Unknown", .risk_table = NULL, .legend_position = c(0.05, 0.1))

dataS60[, table(res_cat, died)]
pl = cowplot::plot_grid(plA, plB, plC, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves_60_res.pdf", pl, width = 15, height = 10, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves_60_res.png", pl, width = 15, height = 10, units = "cm")


# Ethnicity ---------------------------------------------------------------

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[eth_cat == "W"])
plA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.992, 1), .margin = 0.3, .title = "White", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[eth_cat == "A"])
plB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.992, 1), .margin = 0.3, .title = "Asian", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[eth_cat == "B"])
plC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.992, 1), .margin = 0.3, .title = "Black", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[eth_cat == "O"])
plD = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.992, 1), .margin = 0.3, .title = "Other/Mixed/Unknown", .risk_table = NULL, .legend_position = c(0.05, 0.1))


pl = cowplot::plot_grid(plA, plB, plC, plD, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves_60_eth.pdf", pl, width = 15, height = 10, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves_60_eth.png", pl, width = 15, height = 10, units = "cm")


# IMD ---------------------------------------------------------------------

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd1"])
plA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "IMD 1", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd2"])
plB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "IMD 2", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd3"])
plC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "IMD 3", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd4"])
plD = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "IMD 4", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd5"])
plE = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "IMD 5", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd6"])
plF = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "IMD 6", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd7"])
plG = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "IMD 7", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd8"])
plH = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "IMD8", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd9"])
plI = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "IMD 9", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd10"])
plJ = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "IMD 10", .risk_table = NULL, .legend_position = c(0.05, 0.1))


pl = cowplot::plot_grid(plA, plB, plC, plD, plE, plF, plG, plH, plI, plJ, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves_60_imd.pdf", pl, width = 20, height = 30, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves_60_imd.png", pl, width = 20, height = 30, units = "cm")




# NHS region --------------------------------------------------------------

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[NHSER_name == "London"])
plA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "London", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[NHSER_name == "South East"])
plB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "South East", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[NHSER_name == "East of England"])
plF = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "East of England", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[NHSER_name == "South West"])
plC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "South West", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[NHSER_name == "Midlands"])
plD = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Midlands", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[NHSER_name == "North East and Yorkshire"])
plE = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "North East and Yorkshire", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[NHSER_name == "North West"])
plF = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "North West", .risk_table = NULL, .legend_position = c(0.05, 0.1))


pl = cowplot::plot_grid(plA, plB, plC, plD, plE, plF, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves_60_nhs.pdf", pl, width = 15, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves_60_nhs.png", pl, width = 15, height = 15, units = "cm")

# Speciment date ----------------------------------------------------------

## Group the speciment date into 3 week categories.
start_date <- as.Date("2020-11-01")
dt_fmt <- function(x) format(x, "%d-%m-%Y")
dataS60[between(specimen_date, start_date, start_date + 13), spec_2week := paste0(dt_fmt(start_date), " to ", dt_fmt(start_date + 13))]
start_date <- start_date + 14
dataS60[between(specimen_date, start_date, start_date + 13), spec_2week := paste0(dt_fmt(start_date), " to ", dt_fmt(start_date + 13))]
start_date <- start_date + 14
dataS60[between(specimen_date, start_date, start_date + 13), spec_2week := paste0(dt_fmt(start_date), " to ", dt_fmt(start_date + 13))]
start_date <- start_date + 14
dataS60[between(specimen_date, start_date, start_date + 13), spec_2week := paste0(dt_fmt(start_date), " to ", dt_fmt(start_date + 13))]
start_date <- start_date + 14
dataS60[between(specimen_date, start_date, max(specimen_date)), spec_2week := paste0(dt_fmt(start_date), " to ", dt_fmt(max(specimen_date)))]

spec_dates <- unique(dataS60[order(specimen_date)]$spec_2week)
spec_dates

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[1]])
plA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[1], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[2]])
plB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[2], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[3]])
plC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[3], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[4]])
plD = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[4], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[5]])
plE = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[5], .risk_table = NULL, .legend_position = c(0.05, 0.1))

dataS60[,.N, by =spec_2week][order(spec_2week)]

pl = cowplot::plot_grid(plA, plB, plC, plD, plE, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves_60_spec_week.pdf", pl, width = 15, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves_60_spec_week.png", pl, width = 15, height = 15, units = "cm")






# Specimen date age 35 ----------------------------------------------------------

## Group the speciment date into 3 week categories.

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[1] & !age_group %in% c("(0,35]")])
plA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[1], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[2] & !age_group %in% c("(0,35]")])
plB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[2], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[3] & !age_group %in% c("(0,35]")])
plC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[3], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[4] & !age_group %in% c("(0,35]")])
plD = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[4], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[5] & !age_group %in% c("(0,35]")])
plE = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[5], .risk_table = NULL, .legend_position = c(0.05, 0.1))

dataS60[,.N, by =spec_2week][order(spec_2week)]

pl = cowplot::plot_grid(plA, plB, plC, plD, plE, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves_60_spec_week_age35plus.pdf", pl, width = 15, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves_60_spec_week_age35plus.png", pl, width = 15, height = 15, units = "cm")

# Specimen date age 55 ----------------------------------------------------------

## Group the speciment date into 3 week categories.

age_groups_ = c("(0,35]", "(35,55]")

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[1] & !age_group %in% age_groups_])
plA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.96, 1), .margin = 0.3, .title = spec_dates[1], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[2] & !age_group %in% age_groups_])
plB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.96, 1), .margin = 0.3, .title = spec_dates[2], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[3] & !age_group %in% age_groups_])
plC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.96, 1), .margin = 0.3, .title = spec_dates[3], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[4] & !age_group %in% age_groups_])
plD = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.96, 1), .margin = 0.3, .title = spec_dates[4], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[5] & !age_group %in% age_groups_])
plE = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.96, 1), .margin = 0.3, .title = spec_dates[5], .risk_table = NULL, .legend_position = c(0.05, 0.1))

dataS60[,.N, by =spec_2week][order(spec_2week)]

pl = cowplot::plot_grid(plA, plB, plC, plD, plE, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves_60_spec_week_age55plus.pdf", pl, width = 15, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves_60_spec_week_age55plus.png", pl, width = 15, height = 15, units = "cm")

# Specimen date age 70 ----------------------------------------------------------

## Group the speciment date into 3 week categories.

age_groups_ = c("(0,35]", "(35,55]", "(55,70]")

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[1] & !age_group %in% age_groups_])
plA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.86, 1), .margin = 0.3, .title = spec_dates[1], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[2] & !age_group %in% age_groups_])
plB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.86, 1), .margin = 0.3, .title = spec_dates[2], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[3] & !age_group %in% age_groups_])
plC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.86, 1), .margin = 0.3, .title = spec_dates[3], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[4] & !age_group %in% age_groups_])
plD = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.86, 1), .margin = 0.3, .title = spec_dates[4], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[5] & !age_group %in% age_groups_])
plE = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.86, 1), .margin = 0.3, .title = spec_dates[5], .risk_table = NULL, .legend_position = c(0.05, 0.1))

dataS60[,.N, by =spec_2week][order(spec_2week)]

pl = cowplot::plot_grid(plA, plB, plC, plD, plE, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves_60_spec_week_age70plus.pdf", pl, width = 15, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves_60_spec_week_age70plus.png", pl, width = 15, height = 15, units = "cm")



















