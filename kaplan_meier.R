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

# Assemble data set
dataS = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS[, sgtf_label := ifelse(sgtf == 0, "Non-SGTF", "SGTF")]
dataS[, sgtf_label := factor(sgtf_label, c("SGTF", "Non-SGTF"))]


# SURVIVAL

# Full view

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS)
plA = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "Overall", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[age < 70])
plB = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "Under 70", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[age >= 70])
plC = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.92, 1), .margin = 0.3, .title = "70 or older", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "London"])
plD = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "London", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "South East"])
plE = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "South East", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "East of England"])
plF = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "East of England", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "South West"])
plG = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "South West", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "Midlands"])
plH = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "Midlands", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "North East and Yorkshire"])
plI = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "North East and Yorkshire", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "North West"])
plJ = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "North West", .risk_table = NULL, .legend_position = c(0.05, 0.1))


pl = cowplot::plot_grid(plA, plB, plC, plD, plE, plF, plG, plH, plI, plJ, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves.pdf", pl, width = 20, height = 30, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves.png", pl, width = 20, height = 30, units = "cm")


km = survfit(Surv(time, status) ~ sgtf_label, data = dataS)
plA = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "Overall", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[age < 70])
plB = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "Under 70", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[age >= 70])
plC = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.92, 1), .margin = 0.3, .title = "70 or older", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "London"])
plD = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "London", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "South East"])
plE = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "South East", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "East of England"])
plF = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "East of England", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "South West"])
plG = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "South West", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "Midlands"])
plH = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "Midlands", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "North East and Yorkshire"])
plI = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "North East and Yorkshire", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "North West"])
plJ = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.99, 1), .margin = 0.3, .title = "North West", .risk_table = NULL, .legend_position = c(0.05, 0.1))


pl = cowplot::plot_grid(plA, plB, plC, plD, plE, plF, plG, plH, plI, plJ, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves.pdf", pl, width = 20, height = 30, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves.png", pl, width = 20, height = 30, units = "cm")



# Summary view
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS)
plAA = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.995, 1), .margin = 0.3, .title = "Overall", .legend_position = c(0.05, 0.1), .risk_table_base_size = 9)

plAA_inset = KMunicate2(fit = km, time_scale = seq(0, 28, by = 28), .ylim = c(0, 1), .margin = 0.3, .title = "", .risk_table = NULL, .legend_position = "none",.ylab = "", .xlab = "") + scale_y_continuous(breaks = c(0,1), expand = expansion(0)) 

plAA_both <- ggdraw() + draw_plot(plAA) + draw_plot(plAA_inset, x = 0.7, y = .75, width = .25, height = .2)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[age < 70])
plBB = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.998, 1), .margin = 0.3, .title = "Under 70", .legend_position = c(0.05, 0.1), .risk_table_base_size = 9)

plBB_inset = KMunicate2(fit = km, time_scale = seq(0, 28, by = 28), .ylim = c(0, 1), .margin = 0.3, .title = "", .risk_table = NULL, .legend_position = "none",.ylab = "", .xlab = "") + scale_y_continuous(breaks = c(0,1), expand = expansion(0)) 

plBB_both <- ggdraw() + draw_plot(plBB) + draw_plot(plBB_inset, x = 0.7, y = .75, width = .25, height = .2)


km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[age >= 70])
plCC = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.92, 1), .margin = 0.3, .title = "70 or older", .legend_position = c(0.05, 0.1), .risk_table_base_size = 9)

plCC_inset = KMunicate2(fit = km, time_scale = seq(0, 28, by = 28), .ylim = c(0, 1), .margin = 0.3, .title = "", .risk_table = NULL, .legend_position = "none",.ylab = "", .xlab = "", .risk_table_base_size = 8) + scale_y_continuous(breaks = c(0,1), expand = expansion(0)) 

plCC_both <- ggdraw() + draw_plot(plCC) +draw_plot(plCC_inset, x = 0.7, y = .75, width = .25, height = .2)

pl = cowplot::plot_grid(plAA, plBB, plCC, nrow = 1, labels = letters, label_size = 10)
ggsave("./output/kmcurves2.pdf", pl, width = 45, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves2.png", pl, width = 45, height = 15, units = "cm")

pl_inset = cowplot::plot_grid(plAA_both, plBB_both, plCC_both, nrow = 1, labels = letters, label_size = 10)
ggsave("./output/kmcurves2_inset.pdf", pl_inset, width = 45, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves2_inset.png", pl_inset, width = 45, height = 15, units = "cm")



# HAZARD

# Full view

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS)
plA = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0, 0.01), .margin = 0.3, .title = "Overall", .risk_table = NULL, .legend_position = c(0.05, 0.9), .reverse = TRUE)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[age < 70])
plB = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0, 0.01), .margin = 0.3, .title = "Under 70", .risk_table = NULL, .legend_position = c(0.05, 0.9), .reverse = TRUE)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[age >= 70])
plC = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0, 0.1), .margin = 0.3, .title = "70 or older", .risk_table = NULL, .legend_position = c(0.05, 0.9), .reverse = TRUE)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "London"])
plD = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0, 0.01), .margin = 0.3, .title = "London", .risk_table = NULL, .legend_position = c(0.05, 0.9), .reverse = TRUE)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "South East"])
plE = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0, 0.01), .margin = 0.3, .title = "South East", .risk_table = NULL, .legend_position = c(0.05, 0.9), .reverse = TRUE)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "East of England"])
plF = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0, 0.01), .margin = 0.3, .title = "East of England", .risk_table = NULL, .legend_position = c(0.05, 0.9), .reverse = TRUE)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "South West"])
plG = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0, 0.01), .margin = 0.3, .title = "South West", .risk_table = NULL, .legend_position = c(0.05, 0.9), .reverse = TRUE)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "Midlands"])
plH = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0, 0.01), .margin = 0.3, .title = "Midlands", .risk_table = NULL, .legend_position = c(0.05, 0.9), .reverse = TRUE)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "North East and Yorkshire"])
plI = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0, 0.01), .margin = 0.3, .title = "North East and Yorkshire", .risk_table = NULL, .legend_position = c(0.05, 0.9), .reverse = TRUE)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[NHSER_name == "North West"])
plJ = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0, 0.01), .margin = 0.3, .title = "North West", .risk_table = NULL, .legend_position = c(0.05, 0.9), .reverse = TRUE)


pl = cowplot::plot_grid(plA, plB, plC, plD, plE, plF, plG, plH, plI, plJ, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves3.pdf", pl, width = 20, height = 30, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves3.png", pl, width = 20, height = 30, units = "cm")



# Summary view
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS)
plAA = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0, 0.005), .margin = 0.3, .title = "Overall", .legend_position = c(0.05, 0.9), .reverse = TRUE)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[age < 70])
plBB = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0, 0.0015), .margin = 0.3, .title = "Under 70", .legend_position = c(0.05, 0.9), .reverse = TRUE)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[age >= 70])
plCC = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0, 0.1), .margin = 0.3, .title = "70 or older", .legend_position = c(0.05, 0.9), .reverse = TRUE)


pl = cowplot::plot_grid(plAA, plBB, plCC, nrow = 1, labels = letters, label_size = 10)
ggsave("./output/kmcurves4.pdf", pl, width = 45, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves4.png", pl, width = 45, height = 15, units = "cm")




# 60 DAY VIEW
# Assemble data set
dataS60 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 60, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS60[, sgtf_label := ifelse(sgtf == 0, "Non-SGTF", "SGTF")]
dataS60[, sgtf_label := factor(sgtf_label, c("SGTF", "Non-SGTF"))]

# Summary view
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60)
plAA60 = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.995, 1), .margin = 0.3, .title = "Overall", .legend_position = c(0.05, 0.1), .risk_table_base_size = 8, .rel_heights = c(1, 0.2, 0.2), .align = "v")
plAA60_inset = KMunicate2(fit = km, time_scale = c(0, 60), .ylim = c(0, 1), .margin = 0.3, .title = "", .risk_table = NULL, .legend_position = "none",.ylab = "", .xlab = "") + scale_y_continuous(breaks = c(0,1), expand = expansion(0)) 
plAA60_both <- ggdraw() + draw_plot(plAA60) + draw_plot(plAA60_inset, x = 0.5, y = .65, width = .5, height = .35)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[age < 70])
plBB60 = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.998, 1), .margin = 0.3, .title = "Under 70", .legend_position = c(0.05, 0.1), .risk_table_base_size = 8, .rel_heights = c(1, 0.2, 0.2), .align = "v")
plBB60_inset = KMunicate2(fit = km, time_scale = c(0, 60), .ylim = c(0, 1), .margin = 0.3, .title = "", .risk_table = NULL, .legend_position = "none",.ylab = "", .xlab = "") + scale_y_continuous(breaks = c(0,1), expand = expansion(0)) 
plBB60_both <- ggdraw() + draw_plot(plBB60) + draw_plot(plBB60_inset, x = 0.5, y = .65, width = .5, height = .35)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[age >= 70])
plCC60 = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.92, 1), .margin = 0.3, .title = "70 or older", .legend_position = c(0.05, 0.1), .risk_table_base_size = 8, .rel_heights = c(1, 0.2, 0.2), .align = "v")
plCC60_inset = KMunicate2(fit = km, time_scale = c(0, 60), .ylim = c(0, 1), .margin = 0.3, .title = "", .risk_table = NULL, .legend_position = "none",.ylab = "", .xlab = "") + scale_y_continuous(breaks = c(0,1), expand = expansion(0)) 
plCC60_both <- ggdraw() + draw_plot(plCC60) + draw_plot(plCC60_inset, x = 0.5, y = .65, width = .5, height = .35)

pl = cowplot::plot_grid(plAA60, plBB60, plCC60, nrow = 1, labels = letters, label_size = 10)
ggsave("./output/kmcurves5.pdf", pl, width = 45, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves5.png", pl, width = 45, height = 15, units = "cm")

pl60_inset = cowplot::plot_grid(plAA60_both, plBB60_both, plCC60_both, nrow = 1, labels = letters, label_size = 10)
ggsave("./output/kmcurves5_inset.pdf", pl_60_inset, width = 45, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves5_inset.png", pl60_inset, width = 45, height = 15, units = "cm")




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

# Specimen date ----------------------------------------------------------

## Group the specimen date into 3 week categories.
start_date <- as.Date("2020-10-05")
dataS60[between(specimen_date, start_date, start_date + 27), spec_3week := paste0(start_date, " to ", start_date + 27)]
start_date <- start_date + 28
dataS60[between(specimen_date, start_date, start_date + 21), spec_3week := paste0(start_date, " to ", start_date + 20)]
start_date <- start_date + 21
dataS60[between(specimen_date, start_date, start_date + 20), spec_3week := paste0(start_date, " to ", start_date + 20)]
start_date <- start_date + 21
dataS60[between(specimen_date, start_date, max(specimen_date)), spec_3week := paste0(start_date, " to ", max(specimen_date))]

spec_dates <- sort(unique(dataS60$spec_3week))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_3week == spec_dates[1]])
plA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[1], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_3week == spec_dates[2]])
plB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[2], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_3week == spec_dates[3]])
plC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[3], .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_3week == spec_dates[4]])
plD = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[4], .risk_table = NULL, .legend_position = c(0.05, 0.1))

dataS60[,.N, by =spec_3week][order(spec_3week)]

pl = cowplot::plot_grid(plA, plB, plC, plD, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves_60_spec_week.pdf", pl, width = 15, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves_60_spec_week.png", pl, width = 15, height = 15, units = "cm")




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



# IMD ---------------------------------------------------------------------

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd1"])
plA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "imd1", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd2"])
plB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "imd2", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd3"])
plC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "imd3", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd4"])
plD = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "imd4", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd5"])
plE = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "imd5", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd6"])
plF = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "imd6", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd7"])
plG = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "imd7", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd8"])
plH = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "imd8", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd9"])
plI = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "imd9", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd_group == "imd10"])
plJ = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "imd10", .risk_table = NULL, .legend_position = c(0.05, 0.1))


pl = cowplot::plot_grid(plA, plB, plC, plD, plE, plF, plG, plH, plI, plJ, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves_60_imd.pdf", pl, width = 20, height = 30, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves_60_imd.png", pl, width = 20, height = 30, units = "cm")



# Place of residence ------------------------------------------------------

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[res_cat == "Residential"])
plA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Residential", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[res_cat == "Care/Nursing home"])
plB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.98, 1), .margin = 0.3, .title = "Care/Nursing home", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[res_cat == "Other/Unknown"])
plC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Other/Unknown", .risk_table = NULL, .legend_position = c(0.05, 0.1))


pl = cowplot::plot_grid(plA, plB, plC, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves_60_res.pdf", pl, width = 15, height = 10, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves_60_res.png", pl, width = 15, height = 10, units = "cm")


# Ethnicity ---------------------------------------------------------------

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[eth_cat == "W"])
plA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "White", .risk_table = NULL, .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[eth_cat == "A"])
plB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Asian", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[eth_cat == "B"])
plC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Black", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[eth_cat == "O"])
plD = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Other/Mixed/Unknown", .risk_table = NULL, .legend_position = c(0.05, 0.1))


pl = cowplot::plot_grid(plA, plB, plC, plD, ncol = 2, labels = letters, label_size = 10)
ggsave("./output/kmcurves_60_eth.pdf", pl, width = 15, height = 10, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves_60_eth.png", pl, width = 15, height = 10, units = "cm")

