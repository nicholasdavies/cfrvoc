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
cd = complete_data("20210122")

# Assemble data set
dataS = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, prevalence_cutoff = TRUE)
dataS[, sgtf_label := ifelse(sgtf == 0, "Other", "SGTF")]
dataS[, sgtf_label := factor(sgtf_label, c("SGTF", "Other"))]


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



# Summary view
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS)
plAA = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.995, 1), .margin = 0.3, .title = "Overall", .legend_position = c(0.05, 0.1))

plAA_inset = KMunicate2(fit = km, time_scale = seq(0, 28, by = 28), .ylim = c(0, 1), .margin = 0.3, .title = "", .risk_table = NULL, .legend_position = "none",.ylab = "", .xlab = "") + scale_y_continuous(breaks = c(0,1), expand = expansion(0)) 

plAA_both <- 
    ggdraw() +
    draw_plot(plAA) +
    draw_plot(plAA_inset, x = 0.7, y = .75, width = .25, height = .2)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[age < 70])
plBB = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.9985, 1), .margin = 0.3, .title = "Under 70", .legend_position = c(0.05, 0.1))

plBB_inset = KMunicate2(fit = km, time_scale = seq(0, 28, by = 28), .ylim = c(0, 1), .margin = 0.3, .title = "", .risk_table = NULL, .legend_position = "none",.ylab = "", .xlab = "") + scale_y_continuous(breaks = c(0,1), expand = expansion(0)) 


plBB_both <- 
    ggdraw() +
    draw_plot(plBB) +
    draw_plot(plBB_inset, x = 0.7, y = .75, width = .25, height = .2)


km = survfit(Surv(time, status) ~ sgtf_label, data = dataS[age >= 70])
plCC = KMunicate2(fit = km, time_scale = seq(0, 28, by = 4), .ylim = c(0.92, 1), .margin = 0.3, .title = "70 or older", .legend_position = c(0.05, 0.1))

plCC_inset = KMunicate2(fit = km, time_scale = seq(0, 28, by = 28), .ylim = c(0, 1), .margin = 0.3, .title = "", .risk_table = NULL, .legend_position = "none",.ylab = "", .xlab = "") + scale_y_continuous(breaks = c(0,1), expand = expansion(0)) 

plCC_both <- 
    ggdraw() +
    draw_plot(plCC) +
    draw_plot(plCC_inset, x = 0.7, y = .75, width = .25, height = .2)

pl = cowplot::plot_grid(plAA, plBB, plCC, nrow = 1, labels = letters, label_size = 10)
ggsave("./output/kmcurves2.pdf", pl, width = 45, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves2.png", pl, width = 45, height = 15, units = "cm")

pl_inset = cowplot::plot_grid(plAA_both, plBB_both, plCC_both, nrow = 1, labels = letters, label_size = 10)
ggsave("./output/kmcurves2_inset.pdf", pl_both, width = 45, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves2_inset.png", pl_both, width = 45, height = 15, units = "cm")




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
dataS60 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 60, reg_cutoff = 10, P_voc = 0, prevalence_cutoff = TRUE)
dataS60[, sgtf_label := ifelse(sgtf == 0, "Other", "SGTF")]
dataS60[, sgtf_label := factor(sgtf_label, c("SGTF", "Other"))]

# Summary view
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60)
plAA60 = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.995, 1), .margin = 0.3, .title = "Overall", .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[age < 70])
plBB60 = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.998, 1), .margin = 0.3, .title = "Under 70", .legend_position = c(0.05, 0.1))

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[age >= 70])
plCC60 = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.92, 1), .margin = 0.3, .title = "70 or older", .legend_position = c(0.05, 0.1))


pl = cowplot::plot_grid(plAA60, plBB60, plCC60, nrow = 1, labels = letters, label_size = 10)
ggsave("./output/kmcurves5.pdf", pl, width = 45, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/kmcurves5.png", pl, width = 45, height = 15, units = "cm")
