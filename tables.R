# Tables
library(data.table)
library(stringr)
library(lubridate)

source("./phe_data.R")
source("./hazard_data.R")

do_table_prevalence = function(dataS, what)
{
    column = function(tbl)
    {
        paste0(prettyNum(tbl$sgtf, big.mark = ","), " / ", prettyNum(tbl$sgtf + tbl$other, big.mark = ","), " (", round(tbl$pct, 1), "%)")
    }

    tbl_prev = dataS[, .(sgtf = sum(sgtf == 1, na.rm = T), other = sum(sgtf == 0, na.rm = T)), keyby = .(characteristic = get(what))][, .(characteristic, sgtf, other, pct = 100 * sgtf / (sgtf + other))]

    result = data.table(
        what = c(if (what == " ") NULL else what, paste0("    ", as.character(tbl_prev[, characteristic]))),
        prev = c(if (what == " ") NULL else NA, column(tbl_prev))
    )
        
    names(result) = c(" ", "SGTF prevalence")

    return (result)
}

do_table_missingness = function(dataS, what)
{
    column = function(tbl)
    {
        paste0(prettyNum(tbl$absent, big.mark = ","), " / ", prettyNum(tbl$absent + tbl$present, big.mark = ","), " (", round(tbl$pct, 1), "%)")
    }

    tbl_prev = dataS[, .(present = sum(!is.na(sgtf)), absent = sum(is.na(sgtf))), keyby = .(characteristic = get(what))][, .(characteristic, present, absent, pct = 100 * absent / (absent + present))]

    result = data.table(
        what = c(if (what == " ") NULL else what, paste0("    ", as.character(tbl_prev[, characteristic]))),
        prev = c(if (what == " ") NULL else NA, column(tbl_prev))
    )
        
    names(result) = c(" ", "Missingness")

    return (result)
}

do_table_specimens = function(dataS, what, keep_missing = FALSE)
{
    column = function(tbl)
    {
        paste0(prettyNum(tbl$N, big.mark = ","), " (", round(tbl$pct, 1), "%)")
    }

    if (keep_missing) {
        tbl_all =   dataS[,            .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
        tbl_miss =  dataS[is.na(sgtf), .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
        tbl_sgtf =  dataS[sgtf == 1,   .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
        tbl_other = dataS[sgtf == 0,   .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
        
        result = data.table(
            what =  c(if (what == " ") NULL else what, paste0("    ", as.character(tbl_all[, characteristic]))),
            all =   c(if (what == " ") NULL else NA, column(tbl_all)),
            miss =  c(if (what == " ") NULL else NA, column(tbl_miss)),
            sgtf =  c(if (what == " ") NULL else NA, column(tbl_sgtf)),
            other = c(if (what == " ") NULL else NA, column(tbl_other))
        )
        
        names(result) = c(" ", "All", "Missing", "SGTF", "Non-SGTF")
    } else {
        tbl_all =   dataS[!is.na(sgtf), .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
        tbl_sgtf =  dataS[sgtf == 1,    .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
        tbl_other = dataS[sgtf == 0,    .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
        
        result = data.table(
            what =  c(if (what == " ") NULL else what, paste0("    ", as.character(tbl_all[, characteristic]))),
            all =   c(if (what == " ") NULL else NA, column(tbl_all)),
            sgtf =  c(if (what == " ") NULL else NA, column(tbl_sgtf)),
            other = c(if (what == " ") NULL else NA, column(tbl_other))
        )
        
        names(result) = c(" ", "All", "SGTF", "Non-SGTF")
    }
        
    return (result)
}

do_table_deaths = function(dataS, what, keep_missing = FALSE)
{
    column = function(tbl)
    {
        paste0(prettyNum(tbl$N, big.mark = ","), " (", round(tbl$pct, 1), "%)")
    }
    
    if (keep_missing) {
        tbl_all =   dataS[              died == TRUE, .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
        tbl_miss =  dataS[is.na(sgtf) & died == TRUE, .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
        tbl_sgtf =  dataS[sgtf == 1   & died == TRUE, .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
        tbl_other = dataS[sgtf == 0   & died == TRUE, .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
        
        result = data.table(
            what =  c(if (what == " ") NULL else what, paste0("    ", as.character(tbl_all[, characteristic]))),
            all =   c(if (what == " ") NULL else NA, column(tbl_all)),
            miss =  c(if (what == " ") NULL else NA, column(tbl_miss)),
            sgtf =  c(if (what == " ") NULL else NA, column(tbl_sgtf)),
            other = c(if (what == " ") NULL else NA, column(tbl_other))
        )
        
        names(result) = c(" ", "All", "Missing", "SGTF", "Non-SGTF")
    } else {
        tbl_all =   dataS[!is.na(sgtf) & died == TRUE, .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
        tbl_sgtf =  dataS[sgtf == 1    & died == TRUE, .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
        tbl_other = dataS[sgtf == 0    & died == TRUE, .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
        
        result = data.table(
            what =  c(if (what == " ") NULL else what, paste0("    ", as.character(tbl_all[, characteristic]))),
            all =   c(if (what == " ") NULL else NA, column(tbl_all)),
            sgtf =  c(if (what == " ") NULL else NA, column(tbl_sgtf)),
            other = c(if (what == " ") NULL else NA, column(tbl_other))
        )
        
        names(result) = c(" ", "All", "SGTF", "Non-SGTF")
    }

    return (result)
}

do_table_deathrate = function(dataS, what, keep_missing = FALSE)
{
    column = function(tbl)
    {
        paste0(prettyNum(tbl$deaths, big.mark = ","), " / ", prettyNum(round(tbl$days), big.mark = ","), " (", round(tbl$rate, 2), ")")
    }
    
    if (keep_missing) {
        tbl_all =   dataS[,            .(deaths = sum(died == TRUE), days = sum(time)), keyby = .(characteristic = get(what))][, .(characteristic, deaths, days, rate_per10K = 10000 * deaths / days)]
        tbl_miss =  dataS[is.na(sgtf), .(deaths = sum(died == TRUE), days = sum(time)), keyby = .(characteristic = get(what))][, .(characteristic, deaths, days, rate_per10K = 10000 * deaths / days)]
        tbl_sgtf =  dataS[sgtf == 1,   .(deaths = sum(died == TRUE), days = sum(time)), keyby = .(characteristic = get(what))][, .(characteristic, deaths, days, rate_per10K = 10000 * deaths / days)]
        tbl_other = dataS[sgtf == 0,   .(deaths = sum(died == TRUE), days = sum(time)), keyby = .(characteristic = get(what))][, .(characteristic, deaths, days, rate_per10K = 10000 * deaths / days)]
        
        result = data.table(
            what =  c(if (what == " ") NULL else what, paste0("    ", as.character(tbl_all[, characteristic]))),
            all =   c(if (what == " ") NULL else NA, column(tbl_all)),
            miss =  c(if (what == " ") NULL else NA, column(tbl_miss)),
            sgtf =  c(if (what == " ") NULL else NA, column(tbl_sgtf)),
            other = c(if (what == " ") NULL else NA, column(tbl_other))
        )
        
        names(result) = c(" ", "All", "Missing", "SGTF", "Non-SGTF")
    } else {
        tbl_all =   dataS[!is.na(sgtf), .(deaths = sum(died == TRUE), days = sum(time)), keyby = .(characteristic = get(what))][, .(characteristic, deaths, days, rate_per10K = 10000 * deaths / days)]
        tbl_sgtf =  dataS[sgtf == 1,    .(deaths = sum(died == TRUE), days = sum(time)), keyby = .(characteristic = get(what))][, .(characteristic, deaths, days, rate_per10K = 10000 * deaths / days)]
        tbl_other = dataS[sgtf == 0,    .(deaths = sum(died == TRUE), days = sum(time)), keyby = .(characteristic = get(what))][, .(characteristic, deaths, days, rate_per10K = 10000 * deaths / days)]
        
        result = data.table(
            what =  c(if (what == " ") NULL else what, paste0("    ", as.character(tbl_all[, characteristic]))),
            all =   c(if (what == " ") NULL else NA, column(tbl_all)),
            sgtf =  c(if (what == " ") NULL else NA, column(tbl_sgtf)),
            other = c(if (what == " ") NULL else NA, column(tbl_other))
        )
        
        names(result) = c(" ", "All", "SGTF", "Non-SGTF")
    }
    
    return (result)
}

# Load complete data set
cd = complete_data("20210225")
cd[, table(sgtf, useNA = "ifany")]

# Basic missingness stats
cd2 = cd[pillar == "Pillar 2" & specimen_date.x >= "2020-11-01"]
cd2[is.na(age), .N]
cd2[sex == "Unknown", .N]
cd2[is.na(imd_decile), .N]
cd2[ethnicity_final.x == "Unknown", .N]
100 * cd2[ethnicity_final.x == "Unknown", .N] / cd2[, .N]
cd2[cat == "Undetermined" | cat == "", .N]
100 * cd2[cat == "Undetermined" | cat == "", .N] / cd2[, .N]

# Types of missingness
dataS_all = model_data(cd, criterion = "all", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0,
    date_min = "2020-11-01", keep_missing = TRUE)
dataS_u30 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0,
    date_min = "2020-11-01", keep_missing = TRUE)

dataS_all[, table(sgtf, useNA = "ifany")]
dataS_u30[, table(sgtf, useNA = "ifany")]
# 20210225:
# dataS_all[, table(sgtf, useNA = "ifany")]
# sgtf
#          0      1    <NA> 
# all 586707 674539  984017 
# u30 471995 674539 1098729 

# 2,245,263 total tests
#   984,017 missing due to not assessed (43.8%) }
#   114,712 missing due to low CT        (5.1%) } 48.9%
#   471,995 confirmed non-SGTF          (21.0%) }
#   674,539 confirmed SGTF              (30.0%) } 51.1%

# Assemble data set. 
dataS = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0,
    date_min = "2020-11-01", keep_missing = TRUE)
dataS = prep_data(dataS)

# SGTF prevalence
prevalence = rbind(
    do_table_prevalence(dataS, " "),
    do_table_prevalence(dataS, "Sex"),
    do_table_prevalence(dataS, "Age"),
    do_table_prevalence(dataS, "Place of residence"),
    do_table_prevalence(dataS, "Ethnicity"),
    do_table_prevalence(dataS, "IMD decile"),
    do_table_prevalence(dataS, "NHS England region"),
    do_table_prevalence(dataS, "Specimen date")
)

prevalence
fwrite(prevalence, "./output/table_prevalence.csv")

# Missingness
missingness_tbl = rbind(
    do_table_missingness(dataS, " "),
    do_table_missingness(dataS, "Sex"),
    do_table_missingness(dataS, "Age"),
    do_table_missingness(dataS, "Place of residence"),
    do_table_missingness(dataS, "Ethnicity"),
    do_table_missingness(dataS, "IMD decile"),
    do_table_missingness(dataS, "NHS England region"),
    do_table_missingness(dataS, "Specimen date")
)

missingness_tbl
fwrite(missingness_tbl, "./output/table_missingness.csv")

# All tested individuals
specimens = rbind(
    do_table_specimens(dataS, " "),
    do_table_specimens(dataS, "Sex"),
    do_table_specimens(dataS, "Age"),
    do_table_specimens(dataS, "Place of residence"),
    do_table_specimens(dataS, "Ethnicity"),
    do_table_specimens(dataS, "IMD decile"),
    do_table_specimens(dataS, "NHS England region"),
    do_table_specimens(dataS, "Specimen date")
)

specimens_missing = rbind(
    do_table_specimens(dataS, " ", keep_missing = TRUE),
    do_table_specimens(dataS, "Sex", keep_missing = TRUE),
    do_table_specimens(dataS, "Age", keep_missing = TRUE),
    do_table_specimens(dataS, "Place of residence", keep_missing = TRUE),
    do_table_specimens(dataS, "Ethnicity", keep_missing = TRUE),
    do_table_specimens(dataS, "IMD decile", keep_missing = TRUE),
    do_table_specimens(dataS, "NHS England region", keep_missing = TRUE),
    do_table_specimens(dataS, "Specimen date", keep_missing = TRUE)
)

fwrite(specimens, "./output/old_table1.csv")
fwrite(specimens_missing, "./output/old_table1_missing.csv")

table1 = cbind(
    specimens_missing,
    prevalence[, 2],
    missingness_tbl[, 2]
)

fwrite(table1, "./output/extended_table_1.csv")


# All deaths
deaths = rbind(
    do_table_deaths(dataS, " "),
    do_table_deaths(dataS, "Sex"),
    do_table_deaths(dataS, "Age"),
    do_table_deaths(dataS, "Place of residence"),
    do_table_deaths(dataS, "Ethnicity"),
    do_table_deaths(dataS, "IMD decile"),
    do_table_deaths(dataS, "NHS England region"),
    do_table_deaths(dataS, "Specimen date")
)

deaths_missing = rbind(
    do_table_deaths(dataS, " ", keep_missing = TRUE),
    do_table_deaths(dataS, "Sex", keep_missing = TRUE),
    do_table_deaths(dataS, "Age", keep_missing = TRUE),
    do_table_deaths(dataS, "Place of residence", keep_missing = TRUE),
    do_table_deaths(dataS, "Ethnicity", keep_missing = TRUE),
    do_table_deaths(dataS, "IMD decile", keep_missing = TRUE),
    do_table_deaths(dataS, "NHS England region", keep_missing = TRUE),
    do_table_deaths(dataS, "Specimen date", keep_missing = TRUE)
)

fwrite(deaths, "./output/old_table2.csv")
fwrite(deaths_missing, "./output/old_table2_missing.csv")

# Deaths within 28 days, per day of followup
deathrate_28 = rbind(
    do_table_deathrate(dataS, " "),
    do_table_deathrate(dataS, "Sex"),
    do_table_deathrate(dataS, "Age"),
    do_table_deathrate(dataS, "Place of residence"),
    do_table_deathrate(dataS, "Ethnicity"),
    do_table_deathrate(dataS, "IMD decile"),
    do_table_deathrate(dataS, "NHS England region"),
    do_table_deathrate(dataS, "Specimen date")
)
deathrate_28_missing = rbind(
    do_table_deathrate(dataS, " ", keep_missing = TRUE),
    do_table_deathrate(dataS, "Sex", keep_missing = TRUE),
    do_table_deathrate(dataS, "Age", keep_missing = TRUE),
    do_table_deathrate(dataS, "Place of residence", keep_missing = TRUE),
    do_table_deathrate(dataS, "Ethnicity", keep_missing = TRUE),
    do_table_deathrate(dataS, "IMD decile", keep_missing = TRUE),
    do_table_deathrate(dataS, "NHS England region", keep_missing = TRUE),
    do_table_deathrate(dataS, "Specimen date", keep_missing = TRUE)
)

fwrite(deathrate_28, "./output/old_table3.csv")
fwrite(deathrate_28_missing, "./output/extended_table_2.csv")

# With full followup
dataS999 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 999, reg_cutoff = 10, P_voc = 0,
    date_min = "2020-11-01", keep_missing = TRUE)
dataS999 = prep_data(dataS999)

deathrate_full = rbind(
    do_table_deathrate(dataS999, " "),
    do_table_deathrate(dataS999, "Sex"),
    do_table_deathrate(dataS999, "Age"),
    do_table_deathrate(dataS999, "Place of residence"),
    do_table_deathrate(dataS999, "Ethnicity"),
    do_table_deathrate(dataS999, "IMD decile"),
    do_table_deathrate(dataS999, "NHS England region"),
    do_table_deathrate(dataS999, "Specimen date")
)

deathrate_full_missing = rbind(
    do_table_deathrate(dataS999, " ", keep_missing = TRUE),
    do_table_deathrate(dataS999, "Sex", keep_missing = TRUE),
    do_table_deathrate(dataS999, "Age", keep_missing = TRUE),
    do_table_deathrate(dataS999, "Place of residence", keep_missing = TRUE),
    do_table_deathrate(dataS999, "Ethnicity", keep_missing = TRUE),
    do_table_deathrate(dataS999, "IMD decile", keep_missing = TRUE),
    do_table_deathrate(dataS999, "NHS England region", keep_missing = TRUE),
    do_table_deathrate(dataS999, "Specimen date", keep_missing = TRUE)
)

fwrite(deathrate_full_missing, "./output/SI_table_1.csv")


# Absolute risk -----------------------------------------------------------

# effects
summaries = fread("./output/model_summary.csv")
cc28_b   = summaries[model_id == "d28.SGTF.ss.2020-11-01..r10.LTLA:date..cc"  & parameter == "sgtf", coef]
cc28_se  = summaries[model_id == "d28.SGTF.ss.2020-11-01..r10.LTLA:date..cc"  & parameter == "sgtf", se_coef]
ipw28_b  = summaries[model_id == "d28.pVOC.ss.2020-11-01..r10.LTLA:date..ipw" & parameter == "p_voc", coef]
ipw28_se = summaries[model_id == "d28.pVOC.ss.2020-11-01..r10.LTLA:date..ipw" & parameter == "p_voc", se_coef]

#
# CC
#

# Assemble data set
dataS = model_data(cd, criterion = "all", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-08-01", date_max = "2020-10-31", keep_missing = TRUE)

cfr_28 <- dataS[, .(.N, died = sum(died)), by = .(age_group, sex)]
cfr_28[, y := cc28_b]
cfr_28[, y_se := cc28_se]
cfr_28[, cfr := died/N]
cfr_28[, cfr_se := cfr*(1-cfr)/N]
cfr_28[, cfr_lci := cfr - 1.96*cfr_se]
cfr_28[, cfr_uci := cfr + 1.96*cfr_se]
cfr_28[, x := 1-cfr]
cfr_28[, x_se := x*(1-x)/N]

## absolute risk
cfr_28[, abs_risk := 1 - x^exp(y)]
cfr_28[, var_abs_risk := ((exp(y)*x^(exp(y)-1))^2)*x_se^2 + (( exp(y)* (x^(exp(y)) *log(x)))^2)*y_se^2]
cfr_28[, var_abs_risk_se := sqrt(var_abs_risk)]
cfr_28[, uci := abs_risk + 1.96*var_abs_risk_se]
cfr_28[, lci := abs_risk - 1.96*var_abs_risk_se]


fmt_pct <- function(a,b,c){
    paste0(signif(a*100,2), "% (", signif(b*100,2), "–", signif(c*100,2),"%)")
}

abs_risk_tab_28 <- cfr_28[, 
    .("Sex" = sex, "Age" = age_group,
      `Baseline CFR` = fmt_pct(cfr, cfr_lci, cfr_uci),
      `Variant CFR`  = fmt_pct(abs_risk, lci, uci))
]

age_lev <- c("[1,35)", "[35,55)", "[55,70)", "[70,85)", "[85,120)")
age_lab <- c("0-34", "35-54", "55-69", "70-84", "85+")

abs_risk_tab_28[, Age := factor(Age, levels = age_lev, labels = age_lab)]

abs_risk_tab1 <- abs_risk_tab_28
abs_risk_tab1


#
# IPW
#

cfr_28 <- dataS[, .(.N, died = sum(died)), by = .(age_group, sex)]
cfr_28[, y := ipw28_b]
cfr_28[, y_se := ipw28_se]
cfr_28[, cfr := died/N]
cfr_28[, cfr_se := cfr*(1-cfr)/N]
cfr_28[, cfr_lci := cfr - 1.96*cfr_se]
cfr_28[, cfr_uci := cfr + 1.96*cfr_se]
cfr_28[, x:= 1-cfr]
cfr_28[, x_se := x*(1-x)/N]

## absolute risk
cfr_28[, abs_risk := 1 - x^exp(y)]
cfr_28[, var_abs_risk := ((exp(y)*x^(exp(y)-1))^2)*x_se^2 + (( exp(y)* (x^(exp(y)) *log(x)))^2)*y_se^2]
cfr_28[, var_abs_risk_se := sqrt(var_abs_risk)]
cfr_28[, uci := abs_risk + 1.96*var_abs_risk_se]
cfr_28[, lci := abs_risk - 1.96*var_abs_risk_se]


abs_risk_tab_28 <- cfr_28[, 
    .("Sex" = sex, "Age" = age_group,
      `Baseline CFR` = fmt_pct(cfr, cfr_lci, cfr_uci),
      `Variant CFR`  = fmt_pct(abs_risk, lci, uci))
]

age_lev <- c("[1,35)", "[35,55)", "[55,70)", "[70,85)", "[85,120)")
age_lab <- c("0-34", "35-54", "55-69", "70-84", "85+")

abs_risk_tab_28[, Age := factor(Age, levels = age_lev, labels = age_lab)]


abs_risk_tab2 <- abs_risk_tab_28
abs_risk_tab2


abs_risk_tab1
abs_risk_tab2

abs_risk_tab = cbind(abs_risk_tab1[, 1:4], abs_risk_tab2[, 4])
names(abs_risk_tab) = c("Sex", "Age", "Baseline CFR", "SGTF CFR (complete cases)", "pVOC CFR (IPW)")
abs_risk_tab = abs_risk_tab[order(Sex, Age)]
abs_risk_tab[Age == "85+", Age := "85 and older"]
abs_risk_tab
fwrite(abs_risk_tab, "./output/table1.csv")



# Cross tabulated table

dataX = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0,
    date_min = "2020-11-01", keep_missing = FALSE)

dates = list(
    ymd(c("2020-11-01", "2020-11-20")),
    ymd(c("2020-11-01", "2020-11-20")) + 21,
    ymd(c("2020-11-01", "2020-11-20")) + 42,
    ymd(c("2020-11-01", "2020-11-20")) + 63,
    ymd(c("2020-11-01", "2020-11-22")) + 84
)

ages = list(
    c(1, 54),
    c(55, 69),
    c(70, 84),
    c(85, 120)
)

IMDs = list(
    c(1, 2),
    c(3, 8),
    c(9, 10)
)

NHSERs = list(
    c("East of England", "London", "South East"),
    c("Midlands", "North East and Yorkshire", "North West", "South West")
)

# provide # SGTF cases, deaths, follow-up days and death-rates and similarly for non-SGTF cases.

xt_cell = function(dataS, dates, ages, IMDs, NHSERs)
{
    dataS[specimen_date %between% dates &
            age %between% ages &
            imd %between% IMDs &
            NHSER_name %in% NHSERs, 
        paste0("SGTF: ", sum(sgtf == 1), " cases, ", sum(sgtf == 1 & died == TRUE), " deaths / ", sum(time[sgtf == 1]),  " days (", 
            round(10000 * sum(sgtf == 1 & died == TRUE)  / sum(time[sgtf == 1]), 2),  ")\n",
           "Non-SGTF: ", sum(sgtf == 0), " cases, ", sum(sgtf == 0 & died == TRUE), " deaths / ", sum(time[sgtf == 0]), " days (", 
            round(10000 * sum(sgtf == 0 & died == TRUE) / sum(time[sgtf == 0]), 2), ")")]
}

xt = matrix("", nrow = 14*5, ncol = 4)

for (di in 1:5) {
    for (ai in 1:4) {
        for (ii in 1:3) {
            for (ni in 1:2) {
                cat(".")
                cell = xt_cell(dataX, dates[[di]], ages[[ai]], IMDs[[ii]], NHSERs[[ni]])
                xt[ai * 3 + ii + di * 14 - 15, ni + 2] = cell;

                xt[di * 14 - 13, ni + 2] = paste0(dates[[di]][1], " - ", dates[[di]][2]);
                xt[di * 14 - 12, ni + 2] = if (ni == 1) "EE, LD, SE" else "ML, NEY, NW, SW";
                xt[ai * 3 + ii + di * 14 - 15, 1] = paste0("Age ", ages[[ai]][1], " - ", ages[[ai]][2]);
                xt[ai * 3 + ii + di * 14 - 15, 2] = paste0("IMD ", IMDs[[ii]][1], " - ", IMDs[[ii]][2]);
            }
        }
    }
}

fwrite(xt, "./output/SI_table_2.csv")
