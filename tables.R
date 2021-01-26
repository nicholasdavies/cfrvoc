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
cd = complete_data("20210122")
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

fwrite(specimens, "./output/table1.csv")
fwrite(specimens_missing, "./output/table1_missing.csv")

# All deaths
deaths = rbind(
    do_table_deaths(dataS, " "),
    do_table_deaths(dataS, "Sex"),
    do_table_deaths(dataS, "Age"),
    do_table_deaths(dataS, "Place of residence"),
    do_table_deaths(dataS, "Ethnicity"),
    do_table_deaths(dataS, "Index of Multiple Deprivation decile"),
    do_table_deaths(dataS, "NHS England region"),
    do_table_specimens(dataS, "Specimen date")
)

deaths_missing = rbind(
    do_table_deaths(dataS, " ", keep_missing = TRUE),
    do_table_deaths(dataS, "Sex", keep_missing = TRUE),
    do_table_deaths(dataS, "Age", keep_missing = TRUE),
    do_table_deaths(dataS, "Place of residence", keep_missing = TRUE),
    do_table_deaths(dataS, "Ethnicity", keep_missing = TRUE),
    do_table_deaths(dataS, "Index of Multiple Deprivation decile", keep_missing = TRUE),
    do_table_deaths(dataS, "NHS England region", keep_missing = TRUE),
    do_table_specimens(dataS, "Specimen date", keep_missing = TRUE)
)

fwrite(deaths, "./output/table2.csv")
fwrite(deaths_missing, "./output/table2_missing.csv")

# Deaths within 28 days, per day of followup
deathrate_28 = rbind(
    do_table_deathrate(dataS, " "),
    do_table_deathrate(dataS, "Sex"),
    do_table_deathrate(dataS, "Age"),
    do_table_deathrate(dataS, "Place of residence"),
    do_table_deathrate(dataS, "Ethnicity"),
    do_table_deathrate(dataS, "Index of Multiple Deprivation decile"),
    do_table_deathrate(dataS, "NHS England region"),
    do_table_specimens(dataS, "Specimen date")
)
deathrate_28_missing = rbind(
    do_table_deathrate(dataS, " ", keep_missing = TRUE),
    do_table_deathrate(dataS, "Sex", keep_missing = TRUE),
    do_table_deathrate(dataS, "Age", keep_missing = TRUE),
    do_table_deathrate(dataS, "Place of residence", keep_missing = TRUE),
    do_table_deathrate(dataS, "Ethnicity", keep_missing = TRUE),
    do_table_deathrate(dataS, "Index of Multiple Deprivation decile", keep_missing = TRUE),
    do_table_deathrate(dataS, "NHS England region", keep_missing = TRUE),
    do_table_specimens(dataS, "Specimen date", keep_missing = TRUE)
)

fwrite(deathrate_28, "./output/table3.csv")
fwrite(deathrate_28_missing, "./output/table3_missing.csv")

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
    do_table_deathrate(dataS999, "NHS England region")
)

deathrate_full_missing = rbind(
    do_table_deathrate(dataS999, " ", keep_missing = TRUE),
    do_table_deathrate(dataS999, "Sex", keep_missing = TRUE),
    do_table_deathrate(dataS999, "Age", keep_missing = TRUE),
    do_table_deathrate(dataS999, "Place of residence", keep_missing = TRUE),
    do_table_deathrate(dataS999, "Ethnicity", keep_missing = TRUE),
    do_table_deathrate(dataS999, "Index of Multiple Deprivation decile", keep_missing = TRUE),
    do_table_deathrate(dataS999, "NHS England region", keep_missing = TRUE)
)

fwrite(deathrate_full, "./output/table4.csv")
fwrite(deathrate_full_missing, "./output/table4_missing.csv")


# Absolute risk -----------------------------------------------------------

# 60 DAY VIEW
# Assemble data set
dataS = model_data(cd, criterion = "all", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-09-01", date_max = "2020-10-01", keep_missing = TRUE)
dataS60 = model_data(cd, criterion = "all", remove_duplicates = TRUE, death_cutoff = 60, reg_cutoff = 10, P_voc = 0, date_min = "2020-09-01", date_max = "2020-10-01", keep_missing = TRUE)


cfr_28 <- dataS[, .(.N, died = sum(died)), by = .(age_group, sex)]
cfr_28[, y := log(1.30)]
cfr_28[, y_se := ((log(1.30) - log(1.09))/1.96)]
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


# repeat for days at 60 ---------------------------------------------------
cfr_60 <- dataS60[, .(.N, died = sum(died)), by = .(age_group, sex)]
cfr_60[, y := log(1.27)]
cfr_60[, y_se := ((log(1.27) - log(1.069))/1.96)]
cfr_60[, cfr := died/N]
cfr_60[, cfr_se := cfr*(1-cfr)/N]
cfr_60[, cfr_lci := cfr - 1.96*cfr_se]
cfr_60[, cfr_uci := cfr + 1.96*cfr_se]
cfr_60[, x:= 1-cfr]
cfr_60[, x_se := x*(1-x)/N]

## absolute risk
cfr_60[, abs_risk := 1 - x^exp(y)]
cfr_60[, var_abs_risk := ((exp(y)*x^(exp(y)-1))^2)*x_se^2 + (( exp(y)* (x^(exp(y)) *log(x)))^2)*y_se^2]
cfr_60[, var_abs_risk_se := sqrt(var_abs_risk)]
cfr_60[, uci := abs_risk + 1.96*var_abs_risk_se]
cfr_60[, lci := abs_risk - 1.96*var_abs_risk_se]


fmt_pct <- function(a,b,c){
    paste0(signif(round(a*100,3),3), " (", signif(round(b*100,3),3), "-", signif(round(c*100,3),3),")")}


abs_risk_tab_28 <- cfr_28[, .("Sex" = sex, "Age" = age_group,
                              `Base CFR 28 days` = fmt_pct(cfr, cfr_lci, cfr_uci),
                              `New CFR 28 days` = fmt_pct(abs_risk, lci, uci))
                          
]
abs_risk_tab_60 <- cfr_60[, .("Sex" = sex, "Age" = age_group,
                              `Base CFR 60 days` = fmt_pct(cfr, cfr_lci, cfr_uci),
                              `New CFR 60 days` = fmt_pct(abs_risk, lci, uci))
                          
]

age_lev <- c("[1,35)", "[35,55)", "[55,70)", "[70,85)", "[85,120)")
age_lab <- c("0-34", "35-54", "55-69", "70-84", "85+")

abs_risk_tab_28[, Age := factor(Age, levels = age_lev, labels = age_lab)]
abs_risk_tab_60[, Age := factor(Age, levels = age_lev, labels = age_lab)]


abs_risk_tab <- merge(abs_risk_tab_28, abs_risk_tab_60, by = c("Sex", "Age"))[order(Sex, Age)]
fwrite(abs_risk_tab, "./output/table5_abs_risk.csv")
