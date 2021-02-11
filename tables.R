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
cd = complete_data("20210205")
cd[, table(sgtf, useNA = "ifany")]

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
    do_table_prevalence(dataS, "Index of Multiple Deprivation decile"),
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
    do_table_missingness(dataS, "Index of Multiple Deprivation decile"),
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
    do_table_specimens(dataS, "Index of Multiple Deprivation decile"),
    do_table_specimens(dataS, "NHS England region"),
    do_table_specimens(dataS, "Specimen date")
)

specimens_missing = rbind(
    do_table_specimens(dataS, " ", keep_missing = TRUE),
    do_table_specimens(dataS, "Sex", keep_missing = TRUE),
    do_table_specimens(dataS, "Age", keep_missing = TRUE),
    do_table_specimens(dataS, "Place of residence", keep_missing = TRUE),
    do_table_specimens(dataS, "Ethnicity", keep_missing = TRUE),
    do_table_specimens(dataS, "Index of Multiple Deprivation decile", keep_missing = TRUE),
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

fwrite(table1, "./output/table1.csv")


# All deaths
deaths = rbind(
    do_table_deaths(dataS, " "),
    do_table_deaths(dataS, "Sex"),
    do_table_deaths(dataS, "Age"),
    do_table_deaths(dataS, "Place of residence"),
    do_table_deaths(dataS, "Ethnicity"),
    do_table_deaths(dataS, "Index of Multiple Deprivation decile"),
    do_table_deaths(dataS, "NHS England region"),
    do_table_deaths(dataS, "Specimen date")
)

deaths_missing = rbind(
    do_table_deaths(dataS, " ", keep_missing = TRUE),
    do_table_deaths(dataS, "Sex", keep_missing = TRUE),
    do_table_deaths(dataS, "Age", keep_missing = TRUE),
    do_table_deaths(dataS, "Place of residence", keep_missing = TRUE),
    do_table_deaths(dataS, "Ethnicity", keep_missing = TRUE),
    do_table_deaths(dataS, "Index of Multiple Deprivation decile", keep_missing = TRUE),
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
    do_table_deathrate(dataS, "Index of Multiple Deprivation decile"),
    do_table_deathrate(dataS, "NHS England region"),
    do_table_deathrate(dataS, "Specimen date")
)
deathrate_28_missing = rbind(
    do_table_deathrate(dataS, " ", keep_missing = TRUE),
    do_table_deathrate(dataS, "Sex", keep_missing = TRUE),
    do_table_deathrate(dataS, "Age", keep_missing = TRUE),
    do_table_deathrate(dataS, "Place of residence", keep_missing = TRUE),
    do_table_deathrate(dataS, "Ethnicity", keep_missing = TRUE),
    do_table_deathrate(dataS, "Index of Multiple Deprivation decile", keep_missing = TRUE),
    do_table_deathrate(dataS, "NHS England region", keep_missing = TRUE),
    do_table_deathrate(dataS, "Specimen date", keep_missing = TRUE)
)

fwrite(deathrate_28, "./output/old_table3.csv")
fwrite(deathrate_28_missing, "./output/table2.csv")

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
    do_table_deathrate(dataS999, "Index of Multiple Deprivation decile"),
    do_table_deathrate(dataS999, "NHS England region"),
    do_table_deathrate(dataS999, "Specimen date")
)

deathrate_full_missing = rbind(
    do_table_deathrate(dataS999, " ", keep_missing = TRUE),
    do_table_deathrate(dataS999, "Sex", keep_missing = TRUE),
    do_table_deathrate(dataS999, "Age", keep_missing = TRUE),
    do_table_deathrate(dataS999, "Place of residence", keep_missing = TRUE),
    do_table_deathrate(dataS999, "Ethnicity", keep_missing = TRUE),
    do_table_deathrate(dataS999, "Index of Multiple Deprivation decile", keep_missing = TRUE),
    do_table_deathrate(dataS999, "NHS England region", keep_missing = TRUE),
    do_table_deathrate(dataS999, "Specimen date", keep_missing = TRUE)
)

fwrite(deathrate_full_missing, "./output/table_S1.csv")


# Absolute risk -----------------------------------------------------------

# effects
summaries = fread("./output/model_summary.csv")
cc28_b   = summaries[model_id == "d28.SGTF.ss.2020-11-01..r10.LTLA:date..cc"  & parameter == "sgtf", coef]
cc28_se  = summaries[model_id == "d28.SGTF.ss.2020-11-01..r10.LTLA:date..cc"  & parameter == "sgtf", se_coef]
ipw28_b  = summaries[model_id == "d28.SGTF.ss.2020-11-01..r10.LTLA:date..ipw" & parameter == "sgtf", coef]
ipw28_se = summaries[model_id == "d28.SGTF.ss.2020-11-01..r10.LTLA:date..ipw" & parameter == "sgtf", se_coef]

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
cfr_28[, x:= 1-cfr]
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
names(abs_risk_tab) = c("Sex", "Age", "Baseline CFR", "Variant CFR (complete cases)", "Variant CFR (IPW)")
abs_risk_tab = abs_risk_tab[order(Sex, Age)]
abs_risk_tab[Age == "85+", Age := "85 and older"]
abs_risk_tab
fwrite(abs_risk_tab, "./output/table5.csv")
