# Tables
library(data.table)

source("./phe_data.R")
source("./hazard_data.R")

# Load complete data set
cd = complete_data("20210120")

# Assemble data set
dataS = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0,
    date_min = "2020-11-04")

# # For full followup
# dataS = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 999, reg_cutoff = 10, P_voc = 0,
#     date_min = "2020-11-04")

dataS[, sgtf_label := ifelse(sgtf == 0, "Other", "SGTF")]
dataS[, sgtf_label := factor(sgtf_label, c("SGTF", "Other"))]

# Create more descriptive categories for table outputs
dataS[, `Sex` := factor(sex)]
dataS[, `Age` := factor(revalue(age_group,
    c(
        "(0,35]" = "0–34",
        "(35,55]" = "35–54",
        "(55,70]" = "55–69",
        "(70,85]" = "70–84",
        "(85,120]" = "85 and older"
    )), levels = c("0–34", "35–54", "55–69", "70–84", "85 and older"))]
dataS[, `Place of residence` := factor(res_cat)]
dataS[, `Index of Multiple Deprivation decile` := factor(imd, levels = 1:10)]
dataS[, `Ethnicity` := factor(revalue(eth_cat,
    c(
        "W" = "White",
        "A" = "Asian",
        "B" = "Black",
        "O" = "Other/Mixed/Unknown"
    )), levels = c("White", "Asian", "Black", "Other/Mixed/Unknown"))]
dataS[, `NHS England region` := factor(NHSER_name)]

do_table = function(dataS, what)
{
    tbl_all = dataS[!is.na(sgtf), .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
    tbl_sgtf = dataS[sgtf == 1, .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
    tbl_other = dataS[sgtf == 0, .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
    
    column = function(tbl)
    {
        paste0(prettyNum(tbl$N, big.mark = ","), " (", round(tbl$pct, 1), "%)")
    }
    
    result = data.table(
        what = c(what, paste0("    ", as.character(tbl_all[, characteristic]))),
        all = c(NA, column(tbl_all)),
        sgtf = c(NA, column(tbl_sgtf)),
        other = c(NA, column(tbl_other))
    )
    
    names(result) = c("Characteristic", 
        paste0("All (N = ",   prettyNum(sum(tbl_all$N), big.mark = ","), ")"), 
        paste0("SGTF (N = ",  prettyNum(sum(tbl_sgtf$N), big.mark = ","), ")"), 
        paste0("Other (N = ", prettyNum(sum(tbl_other$N), big.mark = ","), ")"))
    
    return (result)
}

do_table_deaths = function(dataS, what)
{
    tbl_all = dataS[!is.na(sgtf) & died == TRUE, .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
    tbl_sgtf = dataS[sgtf == 1 & died == TRUE, .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
    tbl_other = dataS[sgtf == 0 & died == TRUE, .N, keyby = .(characteristic = get(what))][, .(characteristic, N, pct = 100 * N / sum(N))]
    
    column = function(tbl)
    {
        paste0(prettyNum(tbl$N, big.mark = ","), " (", round(tbl$pct, 1), "%)")
    }
    
    result = data.table(
        what = c(what, paste0("    ", as.character(tbl_all[, characteristic]))),
        all = c(NA, column(tbl_all)),
        sgtf = c(NA, column(tbl_sgtf)),
        other = c(NA, column(tbl_other))
    )
    
    names(result) = c("Characteristic", 
        paste0("All (N = ",   prettyNum(sum(tbl_all$N), big.mark = ","), ")"), 
        paste0("SGTF (N = ",  prettyNum(sum(tbl_sgtf$N), big.mark = ","), ")"), 
        paste0("Other (N = ", prettyNum(sum(tbl_other$N), big.mark = ","), ")"))
    
    return (result)
}

do_table_deathrate = function(dataS, what)
{
    tbl_all = dataS[!is.na(sgtf), .(deaths = sum(died == TRUE), days = sum(time)), keyby = .(characteristic = get(what))][, .(characteristic, deaths, days, rate_per10K = 10000 * deaths / days)]
    tbl_sgtf = dataS[sgtf == 1,   .(deaths = sum(died == TRUE), days = sum(time)), keyby = .(characteristic = get(what))][, .(characteristic, deaths, days, rate_per10K = 10000 * deaths / days)]
    tbl_other = dataS[sgtf == 0,  .(deaths = sum(died == TRUE), days = sum(time)), keyby = .(characteristic = get(what))][, .(characteristic, deaths, days, rate_per10K = 10000 * deaths / days)]
    
    column = function(tbl)
    {
        paste0(prettyNum(tbl$deaths, big.mark = ","), " / ", prettyNum(round(tbl$days), big.mark = ","), " (", round(tbl$rate, 2), ")")
    }
    
    result = data.table(
        what = c(what, paste0("    ", as.character(tbl_all[, characteristic]))),
        all = c(NA, column(tbl_all)),
        sgtf = c(NA, column(tbl_sgtf)),
        other = c(NA, column(tbl_other))
    )
    
    names(result) = c("Characteristic", 
        paste0("All (",   prettyNum(sum(tbl_all$deaths),   big.mark = ","), " deaths, ", prettyNum(round(sum(tbl_all$days)),   big.mark = ","), " days of followup)"), 
        paste0("SGTF (",  prettyNum(sum(tbl_sgtf$deaths),  big.mark = ","), " deaths, ", prettyNum(round(sum(tbl_sgtf$days)),  big.mark = ","), " days of followup)"), 
        paste0("Other (", prettyNum(sum(tbl_other$deaths), big.mark = ","), " deaths, ", prettyNum(round(sum(tbl_other$days)), big.mark = ","), " days of followup)"))
    
    return (result)
}

final = rbind(
    do_table(dataS, "Sex"),
    do_table(dataS, "Age"),
    do_table(dataS, "Place of residence"),
    do_table(dataS, "Ethnicity"),
    do_table(dataS, "Index of Multiple Deprivation decile"),
    do_table(dataS, "NHS England region")
)

fwrite(final, "./output/table1.csv")

final_deaths = rbind(
    do_table_deaths(dataS, "Sex"),
    do_table_deaths(dataS, "Age"),
    do_table_deaths(dataS, "Place of residence"),
    do_table_deaths(dataS, "Ethnicity"),
    do_table_deaths(dataS, "Index of Multiple Deprivation decile"),
    do_table_deaths(dataS, "NHS England region")
)

fwrite(final_deaths, "./output/table2.csv")

final_deathrate = rbind(
    do_table_deathrate(dataS, "Sex"),
    do_table_deathrate(dataS, "Age"),
    do_table_deathrate(dataS, "Place of residence"),
    do_table_deathrate(dataS, "Ethnicity"),
    do_table_deathrate(dataS, "Index of Multiple Deprivation decile"),
    do_table_deathrate(dataS, "NHS England region")
)

final_deathrate

fwrite(final_deathrate, "./output/table3.csv")


# With full followup
final_deathrate = rbind(
    do_table_deathrate(dataS, "Sex"),
    do_table_deathrate(dataS, "Age"),
    do_table_deathrate(dataS, "Place of residence"),
    do_table_deathrate(dataS, "Ethnicity"),
    do_table_deathrate(dataS, "Index of Multiple Deprivation decile"),
    do_table_deathrate(dataS, "NHS England region")
)

final_deathrate

fwrite(final_deathrate, "./output/table4.csv")

