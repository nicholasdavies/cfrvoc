## km plots paper
# Kaplan-Meier plots

library(ggplot2)
library(data.table)
library(KMunicate)
library(survival)
library(cowplot)
library(writexl)
library(stringr)

source("./kmunicate2.R")
source("./phe_data.R")
source("./hazard_data.R")

theme_set(theme_cowplot(font_size = 11) + theme(plot.title = element_text(size = 10)))

# Load complete data set
cd = complete_data("20210225")

# 60 DAY VIEW
# Assemble data set
dataS60 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 60, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS60[, sgtf_label := ifelse(sgtf == 0, "Non-SGTF", "SGTF")]
dataS60[, sgtf_label := factor(sgtf_label, c("SGTF", "Non-SGTF"))]

# Summary view for Fig. 1
tsc = seq(0, 60, by = 10)
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60)
plAA60 = KMunicate2(fit = km, time_scale = tsc, .xlim = c(-3, 62), .ylim = c(0.994, 1), .margin = c(0, 0.1, 0, 0.15), .title = NULL, 
    .legend_position = c(0.05, 0.1), .risk_table_base_size = 9, .rel_heights = c(3.5, 1, 1), .align = "v")
plAA60_inset = KMunicate2(fit = km, time_scale = c(0, 60), .ylim = c(0, 1), .margin = 0.3, .title = NULL, .risk_table = NULL, .legend_position = "none",.ylab = "", .xlab = "") + scale_y_continuous(breaks = c(0,1), expand = expansion(0)) 
plAA60_both = ggdraw() + draw_plot(plAA60) + draw_plot(plAA60_inset, x = .47, y = .65, width = .5, height = .35)



# By Sex ------------------------------------------------------------------

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[sex == "Female"])
plSA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Sex: Female", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[sex == "Male"])
plSB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Sex: Male", .risk_table = NULL, .legend_position = c(0.05, 0.1))


# By Age ------------------------------------------------------------------


km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[age_group == "[1,35)"])
plAA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.9998, 1), .margin = 0.3, .title = "Age: 1-34", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[age_group  == "[35,55)"])
plAB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.998, 1), .margin = 0.3, .title = "Age: 35-54", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[age_group == "[55,70)"])
plAC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.985, 1), .margin = 0.3, .title = "Age: 55-69", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[age_group == "[70,85)"])
plAD = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.92, 1), .margin = 0.3, .title = "Age: 70-84", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[age_group == "[85,120)"])
plAE = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.75, 1), .margin = 0.3, .title = "Age: 85+", .risk_table = NULL, .legend_position = c(0.05, 0.1))


# Place of residence ------------------------------------------------------

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[res_cat == "Residential"])
plRA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.995, 1), .margin = 0.3, .title = "Residence: Residential", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[res_cat == "Care/Nursing home"])
plRB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.75, 1), .margin = 0.3, .title = "Residence: Care/Nursing home", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[res_cat == "Other/Unknown"])
plRC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.992, 1), .margin = 0.3, .title = "Residence: Other/Unknown", .risk_table = NULL, .legend_position = c(0.05, 0.1))


# Ethnicity ---------------------------------------------------------------

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[eth_cat == "W"])
plEA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.992, 1), .margin = 0.3, .title = "Ethnicity: White", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[eth_cat == "A"])
plEB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.992, 1), .margin = 0.3, .title = "Ethnicity: Asian", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[eth_cat == "B"])
plEC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.992, 1), .margin = 0.3, .title = "Ethnicity: Black", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[eth_cat == "O"])
plED = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.992, 1), .margin = 0.3, .title = "Ethnicity: Other/Mixed/Unknown", .risk_table = NULL, .legend_position = c(0.05, 0.1))


# IMD ---------------------------------------------------------------------

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd %in% c(1, 2)])
plIA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "IMD: 1-2", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd %in% c(3, 4)])
plIB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "IMD: 3-4", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd %in% c(5, 6)])
plIC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "IMD: 5-6", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd %in% c(7, 8)])
plID = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "IMD: 7-8", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[imd %in% c(9, 10)])
plIE = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "IMD: 9-10", .risk_table = NULL, .legend_position = c(0.05, 0.1))



# NHS region --------------------------------------------------------------

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[NHSER_name == "London"])
plNA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Region: London", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[NHSER_name == "South East"])
plNB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Region: South East", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[NHSER_name == "East of England"])
plNC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Region: East of England", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[NHSER_name == "South West"])
plND = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Region: South West", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[NHSER_name == "Midlands"])
plNE = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Region: Midlands", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[NHSER_name == "North East and Yorkshire"])
plNF = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Region: North East and Yorkshire", .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[NHSER_name == "North West"])
plNG = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = "Region: North West", .risk_table = NULL, .legend_position = c(0.05, 0.1))


# Specimen date ----------------------------------------------------------

## Group the specimen date into 3 week categories.
start_date <- as.Date("2020-11-01")
dt_fmt <- function(x) format(x, "%Y-%m-%d")
dataS60[between(specimen_date, start_date, start_date + 20), spec_2week := paste0(dt_fmt(start_date), " to ", dt_fmt(start_date + 20))]
start_date <- start_date + 21
dataS60[between(specimen_date, start_date, start_date + 20), spec_2week := paste0(dt_fmt(start_date), " to ", dt_fmt(start_date + 20))]
start_date <- start_date + 21
dataS60[between(specimen_date, start_date, start_date + 20), spec_2week := paste0(dt_fmt(start_date), " to ", dt_fmt(start_date + 20))]
start_date <- start_date + 21
dataS60[between(specimen_date, start_date, start_date + 20), spec_2week := paste0(dt_fmt(start_date), " to ", dt_fmt(start_date + 20))]
start_date <- start_date + 21
dataS60[between(specimen_date, start_date, max(specimen_date)), spec_2week := paste0(dt_fmt(start_date), " to ", dt_fmt(max(specimen_date)))]

spec_dates <- unique(dataS60[order(specimen_date)]$spec_2week)

km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[1]])
plDA = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[1], .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[2]])
plDB = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[2], .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[3]])
plDC = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[3], .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[4]])
plDD = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[4], .risk_table = NULL, .legend_position = c(0.05, 0.1))
km = survfit(Surv(time, status) ~ sgtf_label, data = dataS60[spec_2week == spec_dates[5]])
plDE = KMunicate2(fit = km, time_scale = seq(0, 60, by = 10), .ylim = c(0.99, 1), .margin = 0.3, .title = spec_dates[5], .risk_table = NULL, .legend_position = c(0.05, 0.1))


# sex SA-SB
# age AA-AE
# res RA-RC
# eth EA-ED
# nhs region NA-NG
# imd IA-IE
# date DA-DE

big_km = plot_grid(
    plSA, plSB, plAA, plAB, plAC, plAD, plAE,
    plRA, plRB, plRC, plEA, plEB, plEC, plED,
    plNA, plNB, plNC, plND, plNE, plNF, plNG,
    plIA, plIB, plIC, plID, plIE, ggdraw(), ggdraw(),
    plDA, plDB, plDC, plDD, plDE, ggdraw(), ggdraw(),
    byrow = FALSE, nrow = 7, align = "hv", label_size = 11, labels = letters)

ggsave("./output/sfigKM.pdf", big_km, width = 40, height = 40, units = "cm", device = cairo_pdf)
ggsave("./output/sfigKM.png", big_km, width = 40, height = 40, units = "cm")

kmplot_to_data = function(kmlist)
{
    dat = list()
    for (km in kmlist)
    {
        dat[[str_replace_all(km$labels$title, "[: /]", "")]] = km$data
    }
    dat
}

km_sdata = kmplot_to_data(list(
    plSA, plSB, plAA, plAB, plAC, plAD, plAE,
    plRA, plRB, plRC, plEA, plEB, plEC, plED,
    plNA, plNB, plNC, plND, plNE, plNF, plNG,
    plIA, plIB, plIC, plID, plIE, 
    plDA, plDB, plDC, plDD, plDE))

write_xlsx(km_sdata, "./manuscript/sdE_km.xlsx")
