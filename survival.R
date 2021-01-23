library(survival)
library(rms)
library(KMunicate)
library(qs)
library(ggplot2)
library(cowplot)

source("./hazard_data.R")

# Check total number of strata and total number of useful strata
check_stratification = function(data, strata_cols)
{
    if (length(strata_cols) > 0) {
        data[, stratum := do.call(paste, c(.SD, sep = "|")), .SDcols = strata_cols]
        data[, stratum := factor(stratum)]
    }

    strata_names = paste(strata_cols, collapse = " * ")
    n_strata = data[, uniqueN(stratum)]
    sd = data[, .(informative = sum(died) > 0 & sum(sgtf == 0) > 0 & sum(sgtf == 1) > 0, n_deaths = sum(died)), by = stratum]
    n_informative_strata = sd[, sum(informative == TRUE)];
    n_deaths = sd[informative == TRUE, sum(n_deaths)]
    n_with2deaths = sd[informative == TRUE, sum(n_deaths > 1)]
    cat(strata_names, ": ", n_strata, " strata, ", n_informative_strata, " informative strata containing ", n_deaths, "; ", n_with2deaths, " informative strata with 2+ deaths.\n", sep = "")
}


# Workhorse for Cox PH survival analyses
do_cox = function(mdl_id, formula, data, strata_cols, ...)
{
    if (length(strata_cols) > 0) {
        data[, stratum := do.call(paste, c(.SD, sep = "|")), .SDcols = strata_cols]
        data[, stratum := factor(stratum)]
    }
    model = coxph(formula, data = data, ...)
    
    print(summary(model))

    # Extract effect sizes, SEs, CIs and P values
    summ = as.data.table(keep.rownames = TRUE, summary(model)$coefficients)
    names(summ) = c("parameter", "coef", "HR", "se_coef", "z_score", "P")
    summ[, HR.lo95 := exp(coef - 1.96 * se_coef)]
    summ[, HR.hi95 := exp(coef + 1.96 * se_coef)]
    summ[, model_id := mdl_id]
    summ[, formula := paste0(deparse(formula, width.cutoff = 500L), collapse = "")]
    summ[, strata := paste(strata_cols, collapse = " | ")]
    summ[, data_id := data[1, data_id]]
    summ[, last_run := date()]
    
    # Write to summary file
    if (!mdl_id %like% "^test") {
        if (file.exists("./output/model_summary.csv")) {
            file.copy("./output/model_summary.csv", paste0("./backup/model_summary_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv", collapse = ""));
            ms = fread("./output/model_summary.csv");
            ms = ms[model_id != mdl_id];
            ms = rbind(ms, summ);
            fwrite(ms, "./output/model_summary.csv");
        } else {
            fwrite(summ, "./output/model_summary.csv");
        }
        
        # Save model
        qsave(model, paste0("./output/", make.names(mdl_id), ".qs"))
    }
    
    
    model
}

load_model = function(mdl_id) qread(paste0("./output/", make.names(mdl_id), ".qs"))

# Workhorse for Cox PH survival analyses with interaction terms
# marker: e.g. p_voc or sgtf
# group: e.g. age_group, must be categorical
# terms: other terms for the RHS of the model
# data: as in do_cox
# strata_cols: as in do_cox
do_cox_interaction = function(mdl_id, marker, group, terms, data, strata_cols)
{
    data2 = copy(data);
    if (length(strata_cols) > 0) {
        data2[, stratum := do.call(paste, c(.SD, sep = "|")), .SDcols = strata_cols]
        data2[, stratum := factor(stratum)]
    }
    
    group_levels = data[, unique(get(group))];
    marker_names = make.names(paste0(marker, group_levels));
    if (length(group_levels) > 25) {
        stop("More than 25 group levels.")
    }

    for (g in seq_along(group_levels)) {
        data2[, (marker_names[g]) := ifelse(get(group) == group_levels[g], get(marker), 0)]
    }
    
    # Run both formulas and check for interaction LRT
    formula_int = paste("Surv(time, status) ~", paste(c(marker_names, terms, "strata(stratum)"), collapse = " + "));
    formula_noint = paste("Surv(time, status) ~", paste(c(marker, terms, "strata(stratum)"), collapse = " + "));
    model_int = coxph(as.formula(formula_int), data = data2)
    model_noint = coxph(as.formula(formula_noint), data = data2)
    
    # Write to summary file
    if (!mdl_id %like% "^test") {
        # Send output to file
        sink(paste0("./output/interaction_", mdl_id, ".txt"));
        cat("With interaction: ", formula_int, "\n");
        print(summary(model_int))
        cat("\n\nWithout interaction: ", formula_int, "\n");
        print(summary(model_noint))
        cat("\n\nLikelihood ratio test:\n")
        print(anova(model_int, model_noint))
        sink(NULL)
        
        # Save models
        qsave(list(model_int = model_int, model_noint = model_noint), paste0("./output/", make.names(mdl_id), ".qs"))
    }
}




##########################
# SGTF SURVIVAL ANALYSES #
##########################

# Load complete data set
cd = complete_data("20210122")

# 1. MAIN ANALYSES

# Assemble data set to be used for SGTF-based analyses
dataS = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")

# Check stratification
{
    sink("./output/check_stratification.txt")
    check_stratification(dataS, c("NHSER_name", "specimen_week"))
    check_stratification(dataS, c("UTLA_name", "specimen_week"))
    check_stratification(dataS, c("LTLA_name", "specimen_week"))
    check_stratification(dataS, c("NHSER_name", "specimen_date"))
    check_stratification(dataS, c("UTLA_name", "specimen_date"))
    check_stratification(dataS, c("LTLA_name", "specimen_date"))
    sink(NULL)
}

# Linear analyses
do_cox("SGTF + lin age + lin IMD | NHSE + spec week", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_week"))
do_cox("SGTF + lin age + lin IMD | UTLA + spec week", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("UTLA_name", "specimen_week"))
do_cox("SGTF + lin age + lin IMD | LTLA + spec week", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_week"))
do_cox("SGTF + lin age + lin IMD | NHSE + spec date", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_date"))
do_cox("SGTF + lin age + lin IMD | UTLA + spec date", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("UTLA_name", "specimen_date"))
do_cox("SGTF + lin age + lin IMD | LTLA + spec date", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))

# Spline analyses
do_cox("SGTF + spl age + spl IMD | NHSE + spec week", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_week"))
do_cox("SGTF + spl age + spl IMD | UTLA + spec week", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("UTLA_name", "specimen_week"))
do_cox("SGTF + spl age + spl IMD | LTLA + spec week", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_week"))
do_cox("SGTF + spl age + spl IMD | NHSE + spec date", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_date"))
do_cox("SGTF + spl age + spl IMD | UTLA + spec date", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("UTLA_name", "specimen_date"))
do_cox("SGTF + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))

# 2. INTERACTION TESTS
do_cox_interaction("age", "sgtf", "age_group", c("age_group", "sex", "imd",       "eth_cat", "res_cat"), dataS, c("NHSER_name", "specimen_week"))
do_cox_interaction("sex", "sgtf", "sex",       c("age",       "sex", "imd",       "eth_cat", "res_cat"), dataS, c("NHSER_name", "specimen_week"))
do_cox_interaction("imd", "sgtf", "imd_group", c("age",       "sex", "imd_group", "eth_cat", "res_cat"), dataS, c("NHSER_name", "specimen_week"))
do_cox_interaction("eth", "sgtf", "eth_cat",   c("age",       "sex", "imd",       "eth_cat", "res_cat"), dataS, c("NHSER_name", "specimen_week"))
do_cox_interaction("res", "sgtf", "res_cat",   c("age",       "sex", "imd",       "eth_cat", "res_cat"), dataS, c("NHSER_name", "specimen_week"))

# TODO interaction tests in the model with SGTF:time interaction.

# 3. CENSORING LENGTHS
dataS07 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 7,  reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS14 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 14, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS21 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 21, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS28 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS60 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 60, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS999= model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff =999, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
do_cox("Censoring: SGTF_07 + lin age + lin IMD | LTLA + spec week", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS07, c("LTLA_name", "specimen_week"))
do_cox("Censoring: SGTF_14 + lin age + lin IMD | LTLA + spec week", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS14, c("LTLA_name", "specimen_week"))
do_cox("Censoring: SGTF_21 + lin age + lin IMD | LTLA + spec week", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS21, c("LTLA_name", "specimen_week"))
do_cox("Censoring: SGTF_28 + lin age + lin IMD | LTLA + spec week", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS28, c("LTLA_name", "specimen_week"))
do_cox("Censoring: SGTF_60 + lin age + lin IMD | LTLA + spec week", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS60, c("LTLA_name", "specimen_week"))
do_cox("Censoring: SGTF999 + lin age + lin IMD | LTLA + spec week", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS999,c("LTLA_name", "specimen_week"))

# 4. NON-OVERLAPPING PERIODS
splitS = survSplit(Surv(time, status) ~ ., data = dataS, cut = c(0, 7, 14, 21, 28), start = "tstart", end = "tstop")
setDT(splitS)
splitS[, sgtf_w1 := ifelse(sgtf == 1 & tstart ==  0, 1, 0)]
splitS[, sgtf_w2 := ifelse(sgtf == 1 & tstart ==  7, 1, 0)]
splitS[, sgtf_w3 := ifelse(sgtf == 1 & tstart == 14, 1, 0)]
splitS[, sgtf_w4 := ifelse(sgtf == 1 & tstart == 21, 1, 0)]
do_cox("Periods: SGTF + lin age + lin IMD | NHSE + spec week", Surv(tstart, tstop, status) ~ sgtf_w1 + sgtf_w2 + sgtf_w3 + sgtf_w4 + age + sex + imd + eth_cat + res_cat + strata(stratum), splitS, c("NHSER_name", "specimen_week"))
do_cox("Periods: SGTF + lin age + lin IMD | UTLA + spec week", Surv(tstart, tstop, status) ~ sgtf_w1 + sgtf_w2 + sgtf_w3 + sgtf_w4 + age + sex + imd + eth_cat + res_cat + strata(stratum), splitS, c("UTLA_name", "specimen_week"))
do_cox("Periods: SGTF + lin age + lin IMD | LTLA + spec week", Surv(tstart, tstop, status) ~ sgtf_w1 + sgtf_w2 + sgtf_w3 + sgtf_w4 + age + sex + imd + eth_cat + res_cat + strata(stratum), splitS, c("LTLA_name", "specimen_week"))
do_cox("Periods: SGTF + lin age + lin IMD | NHSE + spec date", Surv(tstart, tstop, status) ~ sgtf_w1 + sgtf_w2 + sgtf_w3 + sgtf_w4 + age + sex + imd + eth_cat + res_cat + strata(stratum), splitS, c("NHSER_name", "specimen_date"))
do_cox("Periods: SGTF + lin age + lin IMD | UTLA + spec date", Surv(tstart, tstop, status) ~ sgtf_w1 + sgtf_w2 + sgtf_w3 + sgtf_w4 + age + sex + imd + eth_cat + res_cat + strata(stratum), splitS, c("UTLA_name", "specimen_date"))
do_cox("Periods: SGTF + lin age + lin IMD | LTLA + spec date", Surv(tstart, tstop, status) ~ sgtf_w1 + sgtf_w2 + sgtf_w3 + sgtf_w4 + age + sex + imd + eth_cat + res_cat + strata(stratum), splitS, c("LTLA_name", "specimen_date"))

# 5. TIME-X INTERACTION TERM
event_times = dataS[!is.na(death_date), sort(unique(as.numeric(death_date - specimen_date)))]
dataST = survSplit(Surv(time, status) ~ ., data = dataS, cut = event_times, start = "tstart", end = "tstop")
setDT(dataST)
dataST[, tstop2 := tstop^2]
dataST[, sexM := ifelse(sex == "Male", 1, 0)]
dataST[, resC := ifelse(res_cat == "Care/Nursing home", 1, 0)]
dataST[, resO := ifelse(res_cat == "Other/Unknown", 1, 0)]
dataST[, ethA := ifelse(eth_cat == "A", 1, 0)]
dataST[, ethB := ifelse(eth_cat == "B", 1, 0)]
dataST[, ethO := ifelse(eth_cat == "O", 1, 0)]


# 5a. TIME-SGTF

# With linear term only
# do_cox("Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | NHSE + spec week", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_week"))
# do_cox("Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | UTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_week"))
# do_cox("Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | LTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_week"))
# do_cox("Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | NHSE + spec date", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_date"))
# do_cox("Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | UTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_date"))
do_cox("Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | LTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_date"))

# With quadratic term
# do_cox("Time^2-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | NHSE + spec week", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + sgtf:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_week"))
# do_cox("Time^2-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | UTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + sgtf:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_week"))
# do_cox("Time^2-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | LTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + sgtf:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_week"))
# do_cox("Time^2-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | NHSE + spec date", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + sgtf:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_date"))
# do_cox("Time^2-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | UTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + sgtf:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_date"))
do_cox("Time^2-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | LTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + sgtf:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_date"))

# 5b. TIME-AGE

# With linear term only
# do_cox("Time-age interaction term | NHSE + spec week", Surv(tstart, tstop, status) ~ sgtf + age + age:tstop + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_week"))
# do_cox("Time-age interaction term | UTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + age:tstop + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_week"))
# do_cox("Time-age interaction term | LTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + age:tstop + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_week"))
# do_cox("Time-age interaction term | NHSE + spec date", Surv(tstart, tstop, status) ~ sgtf + age + age:tstop + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_date"))
# do_cox("Time-age interaction term | UTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + age:tstop + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_date"))
do_cox("Time-age interaction term | LTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + age:tstop + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_date"))

# With quadratic term
# do_cox("Time^2-age interaction term | NHSE + spec week", Surv(tstart, tstop, status) ~ sgtf + age + age:tstop + age:tstop2 + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_week"))
# do_cox("Time^2-age interaction term | UTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + age:tstop + age:tstop2 + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_week"))
# do_cox("Time^2-age interaction term | LTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + age:tstop + age:tstop2 + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_week"))
# do_cox("Time^2-age interaction term | NHSE + spec date", Surv(tstart, tstop, status) ~ sgtf + age + age:tstop + age:tstop2 + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_date"))
# do_cox("Time^2-age interaction term | UTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + age:tstop + age:tstop2 + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_date"))
do_cox("Time^2-age interaction term | LTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + age:tstop + age:tstop2 + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_date"))

# 5c. TIME-SEX

# With linear term only
# do_cox("Time-sex interaction term | NHSE + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sexM + sexM:tstop + imd + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_week"))
# do_cox("Time-sex interaction term | UTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sexM + sexM:tstop + imd + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_week"))
# do_cox("Time-sex interaction term | LTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sexM + sexM:tstop + imd + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_week"))
# do_cox("Time-sex interaction term | NHSE + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sexM + sexM:tstop + imd + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_date"))
# do_cox("Time-sex interaction term | UTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sexM + sexM:tstop + imd + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_date"))
do_cox("Time-sex interaction term | LTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sexM + sexM:tstop + imd + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_date"))

# With quadratic term
# do_cox("Time^2-sex interaction term | NHSE + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sexM + sexM:tstop + sexM:tstop2 + imd + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_week"))
# do_cox("Time^2-sex interaction term | UTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sexM + sexM:tstop + sexM:tstop2 + imd + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_week"))
# do_cox("Time^2-sex interaction term | LTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sexM + sexM:tstop + sexM:tstop2 + imd + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_week"))
# do_cox("Time^2-sex interaction term | NHSE + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sexM + sexM:tstop + sexM:tstop2 + imd + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_date"))
# do_cox("Time^2-sex interaction term | UTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sexM + sexM:tstop + sexM:tstop2 + imd + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_date"))
do_cox("Time^2-sex interaction term | LTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sexM + sexM:tstop + sexM:tstop2 + imd + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_date"))

# 5d. TIME-IMD

# With linear term only
# do_cox("Time-IMD interaction term | NHSE + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + imd:tstop + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_week"))
# do_cox("Time-IMD interaction term | UTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + imd:tstop + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_week"))
# do_cox("Time-IMD interaction term | LTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + imd:tstop + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_week"))
# do_cox("Time-IMD interaction term | NHSE + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + imd:tstop + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_date"))
# do_cox("Time-IMD interaction term | UTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + imd:tstop + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_date"))
do_cox("Time-IMD interaction term | LTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + imd:tstop + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_date"))

# With quadratic term
# do_cox("Time^2-IMD interaction term | NHSE + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + imd:tstop + imd:tstop2 + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_week"))
# do_cox("Time^2-IMD interaction term | UTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + imd:tstop + imd:tstop2 + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_week"))
# do_cox("Time^2-IMD interaction term | LTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + imd:tstop + imd:tstop2 + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_week"))
# do_cox("Time^2-IMD interaction term | NHSE + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + imd:tstop + imd:tstop2 + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_date"))
# do_cox("Time^2-IMD interaction term | UTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + imd:tstop + imd:tstop2 + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_date"))
do_cox("Time^2-IMD interaction term | LTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + imd:tstop + imd:tstop2 + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_date"))

# 5e. TIME-ETHNICITY

# With linear term only
# do_cox("Time-ethnicity interaction term | NHSE + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + ethA + ethA:tstop + ethB + ethB:tstop + ethO + ethO:tstop + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_week"))
# do_cox("Time-ethnicity interaction term | UTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + ethA + ethA:tstop + ethB + ethB:tstop + ethO + ethO:tstop + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_week"))
# do_cox("Time-ethnicity interaction term | LTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + ethA + ethA:tstop + ethB + ethB:tstop + ethO + ethO:tstop + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_week"))
# do_cox("Time-ethnicity interaction term | NHSE + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + ethA + ethA:tstop + ethB + ethB:tstop + ethO + ethO:tstop + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_date"))
# do_cox("Time-ethnicity interaction term | UTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + ethA + ethA:tstop + ethB + ethB:tstop + ethO + ethO:tstop + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_date"))
do_cox("Time-ethnicity interaction term | LTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + ethA + ethA:tstop + ethB + ethB:tstop + ethO + ethO:tstop + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_date"))

# With quadratic term
# do_cox("Time^2-ethnicity interaction term | NHSE + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + ethA + ethA:tstop + ethA:tstop2 + ethB + ethB:tstop + ethB:tstop2 + ethO + ethO:tstop + ethO:tstop2 + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_week"))
# do_cox("Time^2-ethnicity interaction term | UTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + ethA + ethA:tstop + ethA:tstop2 + ethB + ethB:tstop + ethB:tstop2 + ethO + ethO:tstop + ethO:tstop2 + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_week"))
# do_cox("Time^2-ethnicity interaction term | LTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + ethA + ethA:tstop + ethA:tstop2 + ethB + ethB:tstop + ethB:tstop2 + ethO + ethO:tstop + ethO:tstop2 + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_week"))
# do_cox("Time^2-ethnicity interaction term | NHSE + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + ethA + ethA:tstop + ethA:tstop2 + ethB + ethB:tstop + ethB:tstop2 + ethO + ethO:tstop + ethO:tstop2 + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_date"))
# do_cox("Time^2-ethnicity interaction term | UTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + ethA + ethA:tstop + ethA:tstop2 + ethB + ethB:tstop + ethB:tstop2 + ethO + ethO:tstop + ethO:tstop2 + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_date"))
do_cox("Time^2-ethnicity interaction term | LTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + ethA + ethA:tstop + ethA:tstop2 + ethB + ethB:tstop + ethB:tstop2 + ethO + ethO:tstop + ethO:tstop2 + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_date"))

# 5f. TIME-RESIDENCE TYPE

# With linear term only
# do_cox("Time-residence interaction term | NHSE + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + eth_cat + resC + resC:tstop + resO + resO:tstop + strata(stratum), dataST, c("NHSER_name", "specimen_week"))
# do_cox("Time-residence interaction term | UTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + eth_cat + resC + resC:tstop + resO + resO:tstop + strata(stratum), dataST, c("UTLA_name", "specimen_week"))
# do_cox("Time-residence interaction term | LTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + eth_cat + resC + resC:tstop + resO + resO:tstop + strata(stratum), dataST, c("LTLA_name", "specimen_week"))
# do_cox("Time-residence interaction term | NHSE + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + eth_cat + resC + resC:tstop + resO + resO:tstop + strata(stratum), dataST, c("NHSER_name", "specimen_date"))
# do_cox("Time-residence interaction term | UTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + eth_cat + resC + resC:tstop + resO + resO:tstop + strata(stratum), dataST, c("UTLA_name", "specimen_date"))
do_cox("Time-residence interaction term | LTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + eth_cat + resC + resC:tstop + resO + resO:tstop + strata(stratum), dataST, c("LTLA_name", "specimen_date"))

# With quadratic term
# do_cox("Time^2-residence interaction term | NHSE + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + eth_cat + resC + resC:tstop + resC:tstop2 + resO + resO:tstop + resO:tstop2 + strata(stratum), dataST, c("NHSER_name", "specimen_week"))
# do_cox("Time^2-residence interaction term | UTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + eth_cat + resC + resC:tstop + resC:tstop2 + resO + resO:tstop + resO:tstop2 + strata(stratum), dataST, c("UTLA_name", "specimen_week"))
# do_cox("Time^2-residence interaction term | LTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + eth_cat + resC + resC:tstop + resC:tstop2 + resO + resO:tstop + resO:tstop2 + strata(stratum), dataST, c("LTLA_name", "specimen_week"))
# do_cox("Time^2-residence interaction term | NHSE + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + eth_cat + resC + resC:tstop + resC:tstop2 + resO + resO:tstop + resO:tstop2 + strata(stratum), dataST, c("NHSER_name", "specimen_date"))
# do_cox("Time^2-residence interaction term | UTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + eth_cat + resC + resC:tstop + resC:tstop2 + resO + resO:tstop + resO:tstop2 + strata(stratum), dataST, c("UTLA_name", "specimen_date"))
do_cox("Time^2-residence interaction term | LTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + age + sex + imd + eth_cat + resC + resC:tstop + resC:tstop2 + resO + resO:tstop + resO:tstop2 + strata(stratum), dataST, c("LTLA_name", "specimen_date"))

# 6. SENSITIVITY ANALYSES

# 6a. Cutoff by date range
dataS_sep = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-09-01")
do_cox("Sensitivity: Sep onwards (SGTF + lin age + lin IMD | LTLA + spec date)", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS_sep, c("LTLA_name", "specimen_date"))
dataS_oct = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-10-01")
do_cox("Sensitivity: Oct onwards (SGTF + lin age + lin IMD | LTLA + spec date)", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS_oct, c("LTLA_name", "specimen_date"))
dataS_nov = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
do_cox("Sensitivity: Nov onwards (SGTF + lin age + lin IMD | LTLA + spec date)", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS_nov, c("LTLA_name", "specimen_date"))
dataS_dec = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-12-01")
do_cox("Sensitivity: Dec onwards (SGTF + lin age + lin IMD | LTLA + spec date)", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS_dec, c("LTLA_name", "specimen_date"))
dataS_LTLAprev = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, prevalence_cutoff = TRUE)
do_cox("Sensitivity: Prevalence cutoff (SGTF + lin age + lin IMD | LTLA + spec date)", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS_dec, c("LTLA_name", "specimen_date"))

# 6b. No registration cutoff
dataSc = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 0, P_voc = 0, date_min = "2020-11-01")
do_cox("Sensitivity: no registration cutoff (SGTF + lin age + lin IMD | LTLA + spec date)", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataSc, c("LTLA_name", "specimen_date"))

# 6c. Adjusting for, rather than stratifying by, specimen week and NHSE region
dataS2 = copy(dataS)
# Seems we need to exclude weeks with fewer data points -- check when updating data.
dataS2 = dataS2[specimen_week >= "2020-11-02" & specimen_week < "2021-01-11"]
dataS2[, specimen_week_f := factor(specimen_week)]
do_cox("Sensitivity: SGTF + lin age + lin IMD + NHSE * spec week", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + NHSER_name * specimen_week_f, dataS2, c())


###########################
# P_VOC SURVIVAL ANALYSES #
###########################

# 1. MAIN ANALYSES

# Assemble data set to be used for p_voc-based analyses
dataP = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0)

# Linear
do_cox("p_voc + lin age + lin IMD | NHSE + spec week", Surv(time, status) ~ p_voc + age + sex + imd + eth_cat + res_cat + strata(stratum), dataP, c("NHSER_name", "specimen_week"))
do_cox("p_voc + lin age + lin IMD | UTLA + spec week", Surv(time, status) ~ p_voc + age + sex + imd + eth_cat + res_cat + strata(stratum), dataP, c("UTLA_name", "specimen_week"))
do_cox("p_voc + lin age + lin IMD | LTLA + spec week", Surv(time, status) ~ p_voc + age + sex + imd + eth_cat + res_cat + strata(stratum), dataP, c("LTLA_name", "specimen_week"))
do_cox("p_voc + lin age + lin IMD | NHSE + spec date", Surv(time, status) ~ p_voc + age + sex + imd + eth_cat + res_cat + strata(stratum), dataP, c("NHSER_name", "specimen_date"))
do_cox("p_voc + lin age + lin IMD | UTLA + spec date", Surv(time, status) ~ p_voc + age + sex + imd + eth_cat + res_cat + strata(stratum), dataP, c("UTLA_name", "specimen_date"))
do_cox("p_voc + lin age + lin IMD | LTLA + spec date", Surv(time, status) ~ p_voc + age + sex + imd + eth_cat + res_cat + strata(stratum), dataP, c("LTLA_name", "specimen_date"))

# Splines
do_cox("p_voc + spl age + spl IMD | NHSE + spec week", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataP, c("NHSER_name", "specimen_week"))
do_cox("p_voc + spl age + spl IMD | UTLA + spec week", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataP, c("UTLA_name", "specimen_week"))
do_cox("p_voc + spl age + spl IMD | LTLA + spec week", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataP, c("LTLA_name", "specimen_week"))
do_cox("p_voc + spl age + spl IMD | NHSE + spec date", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataP, c("NHSER_name", "specimen_date"))
do_cox("p_voc + spl age + spl IMD | UTLA + spec date", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataP, c("UTLA_name", "specimen_date"))
do_cox("p_voc + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataP, c("LTLA_name", "specimen_date"))

# 2. TIME-P_VOC
event_times = dataP[!is.na(death_date), sort(unique(as.numeric(death_date - specimen_date)))]
dataPT = survSplit(Surv(time, status) ~ ., data = dataP, cut = event_times, start = "tstart", end = "tstop")
setDT(dataPT)
dataPT[, tstop2 := tstop^2]

# With linear term only
# do_cox("Time-p_voc interaction term: p_voc + p_voc:tstop + lin age + lin IMD | NHSE + spec week", Surv(tstart, tstop, status) ~ p_voc + p_voc:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataPT, c("NHSER_name", "specimen_week"))
# do_cox("Time-p_voc interaction term: p_voc + p_voc:tstop + lin age + lin IMD | UTLA + spec week", Surv(tstart, tstop, status) ~ p_voc + p_voc:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataPT, c("UTLA_name", "specimen_week"))
# do_cox("Time-p_voc interaction term: p_voc + p_voc:tstop + lin age + lin IMD | LTLA + spec week", Surv(tstart, tstop, status) ~ p_voc + p_voc:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataPT, c("LTLA_name", "specimen_week"))
# do_cox("Time-p_voc interaction term: p_voc + p_voc:tstop + lin age + lin IMD | NHSE + spec date", Surv(tstart, tstop, status) ~ p_voc + p_voc:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataPT, c("NHSER_name", "specimen_date"))
# do_cox("Time-p_voc interaction term: p_voc + p_voc:tstop + lin age + lin IMD | UTLA + spec date", Surv(tstart, tstop, status) ~ p_voc + p_voc:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataPT, c("UTLA_name", "specimen_date"))
do_cox("Time-p_voc interaction term: p_voc + p_voc:tstop + lin age + lin IMD | LTLA + spec date", Surv(tstart, tstop, status) ~ p_voc + p_voc:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataPT, c("LTLA_name", "specimen_date"))

# With quadratic term
# do_cox("Time^2-p_voc interaction term: p_voc + p_voc:tstop + lin age + lin IMD | NHSE + spec week", Surv(tstart, tstop, status) ~ p_voc + p_voc:tstop + p_voc:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataPT, c("NHSER_name", "specimen_week"))
# do_cox("Time^2-p_voc interaction term: p_voc + p_voc:tstop + lin age + lin IMD | UTLA + spec week", Surv(tstart, tstop, status) ~ p_voc + p_voc:tstop + p_voc:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataPT, c("UTLA_name", "specimen_week"))
# do_cox("Time^2-p_voc interaction term: p_voc + p_voc:tstop + lin age + lin IMD | LTLA + spec week", Surv(tstart, tstop, status) ~ p_voc + p_voc:tstop + p_voc:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataPT, c("LTLA_name", "specimen_week"))
# do_cox("Time^2-p_voc interaction term: p_voc + p_voc:tstop + lin age + lin IMD | NHSE + spec date", Surv(tstart, tstop, status) ~ p_voc + p_voc:tstop + p_voc:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataPT, c("NHSER_name", "specimen_date"))
# do_cox("Time^2-p_voc interaction term: p_voc + p_voc:tstop + lin age + lin IMD | UTLA + spec date", Surv(tstart, tstop, status) ~ p_voc + p_voc:tstop + p_voc:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataPT, c("UTLA_name", "specimen_date"))
do_cox("Time^2-p_voc interaction term: p_voc + p_voc:tstop + lin age + lin IMD | LTLA + spec date", Surv(tstart, tstop, status) ~ p_voc + p_voc:tstop + p_voc:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataPT, c("LTLA_name", "specimen_date"))



# # Try and get some output......
# 
# summaries = fread("./output/model_summary.csv")
# pl = ggplot(summaries) +
#     geom_pointrange(aes(x = parameter, y = HR, ymin = HR.lo95, ymax = HR.hi95), fatten = 0.1) +
#     facet_wrap(~model_id, scales = "free", ncol = 3) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4))
# ggsave("./output/summary.pdf", pl, width = 40, height = 100, units = "cm", limitsize = FALSE)
# 
# 
# fwrite(summaries[parameter == "sgtf"], "./output/temp_output.csv")

# from Ruth: 
# It would be nice to see a plot of the hazard ratio over time, 
# including confidence intervals. The log HR is linear with time 
# according to your model: A+Bt. I would obtain the HR for a series 
# of times (e.g. on a grid from 0 to 28 in increments of 1), plus 
# the 95% CI. The 95% CI at time t is
# 
# (A+Bt) +/- 1.96*SE where
# SE=Sqrt(Var(A)+t^2*var(B)+2*t*cov(A,B))
#  
# Var(A) is the square of the SE for A, and similar for var(B). 
# You can get the covariance between A and B from the covariance 
# matrix for the parameters – i.e. using the appropriate term in 
# vcov(your model). You can also get var(A) and var(B) from the 
# diagonal on the covariance matrix.

plot_hazard = function(marker, constant_model, linear_varying_model, ylimits)
{
    mv = load_model(linear_varying_model)
    
    marker_int = paste0(marker, ":tstop")
    
    A = coefficients(mv)[[marker]]
    B = coefficients(mv)[[marker_int]]
    varA = vcov(mv)[marker, marker]
    varB = vcov(mv)[marker_int, marker_int]
    covAB = vcov(mv)[marker, marker_int]
    
    haz = data.table(t = 0:28)
    haz[, x.ct := A + B*t]
    haz[, x.se := sqrt(varA + t^2 * varB + 2*t*covAB)]
    haz[, x.hi := x.ct + 1.96 * x.se]
    haz[, x.lo := x.ct - 1.96 * x.se]
    
    mc = load_model(constant_model)
    
    haz[, c.ct := coefficients(mc)[[marker]]]
    haz[, c.hi := c.ct + 1.96 * summary(mc)$coefficients[marker, "se(coef)"]]
    haz[, c.lo := c.ct - 1.96 * summary(mc)$coefficients[marker, "se(coef)"]]
    
    ggplot(haz) +
        geom_ribbon(aes(x = t, ymin = exp(c.lo), ymax = exp(c.hi), fill = "Constant HR (proportional hazards)"), alpha = 0.4) +
        geom_ribbon(aes(x = t, ymin = exp(x.lo), ymax = exp(x.hi), fill = "Time-varying HR"), alpha = 0.4) +
        geom_line(aes(x = t, y = exp(c.ct), colour = "Constant HR (proportional hazards)")) +
        geom_line(aes(x = t, y = exp(x.ct), colour = "Time-varying HR")) +
        geom_hline(aes(yintercept = 1), linetype = "33", size = 0.25) +
        labs(x = "Days post specimen", y = "Hazard ratio (HR)", fill = "Model", colour = "Model") +
        scale_x_continuous(breaks = c(0, 7, 14, 21, 28)) +
        theme(legend.position = c(0.1, 0.9)) +
        ylim(ylimits)
}

theme_set(theme_cowplot(font_size = 10))


hp_sgtf1 = plot_hazard("sgtf", "SGTF + lin age + lin IMD | LTLA + spec date", "Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | LTLA + spec date", c(0, NA))
ggsave("./output/time_varying_sgtf.pdf", hp_sgtf1, width = 15, height = 10, units = "cm", useDingbats = FALSE)
ggsave("./output/time_varying_sgtf.png", hp_sgtf1, width = 15, height = 10, units = "cm")

hp_sgtf = plot_hazard("sgtf", "SGTF + lin age + lin IMD | LTLA + spec date", "Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | LTLA + spec date", c(0, 4.5))
hp_voc = plot_hazard("p_voc", "p_voc + lin age + lin IMD | LTLA + spec date", "Time-p_voc interaction term: p_voc + p_voc:tstop + lin age + lin IMD | LTLA + spec date", c(0, 4.5))
hp_compare = plot_grid(hp_sgtf, hp_voc, nrow = 1, labels = letters, label_size = 10)
ggsave("./output/time_varying_compare.pdf", hp_compare, width = 20, height = 10, units = "cm", useDingbats = FALSE)
ggsave("./output/time_varying_compare.png", hp_compare, width = 20, height = 10, units = "cm")


# Effects of SGTF

haz_summ = function(days)
{
    summaries = fread("./output/model_summary.csv");
    models = summaries[, unique(model_id)];
    
    effects = NULL
    for (m in models) {
        subset = summaries[model_id == m];
        subset
    }
}

summaries

pl = ggplot(summaries[parameter %in% c("sgtf", "p_voc") & !model_id %like% "Time.{0,2}-(SGTF|p_voc)"]) +
    geom_pointrange(aes(y = model_id, x = HR, xmin = HR.lo95, xmax = HR.hi95), fatten = 0.1) +
    geom_vline(aes(xintercept = 1), linetype = "33", size = 0.25) +
    labs(x = "Hazard ratio", y = NULL)
pl
ggsave("./output/summary.pdf", pl, width = 20, height = 20, units = "cm")


fwrite(summaries[parameter == "sgtf"], "./output/temp_output.csv")
