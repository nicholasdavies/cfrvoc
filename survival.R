library(survival)
library(rms)
library(KMunicate)
library(qs)

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
    cat(strata_names, ": ", n_strata, " strata, ", n_informative_strata, " informative strata containing ", n_deaths, " deaths.\n", sep = "")
}


# Workhorse for Cox PH survival analyses
do_cox = function(mdl_id, formula, data, strata_cols)
{
    if (length(strata_cols) > 0) {
        data[, stratum := do.call(paste, c(.SD, sep = "|")), .SDcols = strata_cols]
        data[, stratum := factor(stratum)]
    }
    model = coxph(formula, data = data)
    
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

# Load complete data set
cd = complete_data("20210121")




##########################
# SGTF SURVIVAL ANALYSES #
##########################

# 1. MAIN ANALYSES

# # Determining a cutoff point for SGTF data based upon VOC growth modelling
# data = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0)
# data[sgtf == 1, mean(p_voc), keyby = specimen_date][V1 >= 0.5] # Earliest specimen_date: Nov 4 2020

# Assemble data set to be used for SGTF-based analyses
dataS = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, prevalence_cutoff = TRUE)

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

# 3. CENSORING LENGTHS
dataS07 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 7,  reg_cutoff = 10, P_voc = 0, prevalence_cutoff = TRUE)
dataS14 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 14, reg_cutoff = 10, P_voc = 0, prevalence_cutoff = TRUE)
dataS21 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 21, reg_cutoff = 10, P_voc = 0, prevalence_cutoff = TRUE)
dataS28 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, prevalence_cutoff = TRUE)
dataS60 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 60, reg_cutoff = 10, P_voc = 0, prevalence_cutoff = TRUE)
dataS999= model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff =999, reg_cutoff = 10, P_voc = 0, prevalence_cutoff = TRUE)
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

# 5. LINEAR TIME-VARIANT INTERACTION TERM
event_times = dataS[!is.na(death_date), sort(unique(as.numeric(death_date - specimen_date)))]
event_times
dataST = survSplit(Surv(time, status) ~ ., data = dataS, cut = event_times, start = "tstart", end = "tstop")
nrow(dataS)
nrow(dataST)
setDT(dataST)

do_cox("Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | NHSE + spec week", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_week"))
do_cox("Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | UTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_week"))
do_cox("Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | LTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_week"))
do_cox("Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | NHSE + spec date", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_date"))
do_cox("Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | UTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_date"))
do_cox("Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | LTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_date"))

# With quadratic term
dataST[, tstop2 := tstop^2]
do_cox("Time^2-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | NHSE + spec week", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + sgtf:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_week"))
do_cox("Time^2-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | UTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + sgtf:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_week"))
do_cox("Time^2-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | LTLA + spec week", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + sgtf:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_week"))
do_cox("Time^2-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | NHSE + spec date", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + sgtf:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("NHSER_name", "specimen_date"))
do_cox("Time^2-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | UTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + sgtf:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("UTLA_name", "specimen_date"))
do_cox("Time^2-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | LTLA + spec date", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + sgtf:tstop2 + age + sex + imd + eth_cat + res_cat + strata(stratum), dataST, c("LTLA_name", "specimen_date"))

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

# 6b. No registration cutoff
dataSc = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 0, P_voc = 0, prevalence_cutoff = TRUE)
do_cox("Sensitivity: no registration cutoff (SGTF + lin age + lin IMD | LTLA + spec date)", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataSc, c("LTLA_name", "specimen_date"))

# 6c. Adjusting for, rather than stratifying by, specimen week and NHSE region
dataS2 = dataS[specimen_week < "2021-01-04"] # Seems to be not enough data points here for regression -- check when updating data.
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

A = coefficients(m)[["sgtf"]]
B = coefficients(m)[["sgtf:tstop"]]
varA = vcov(m)["sgtf", "sgtf"]
varB = vcov(m)["sgtf:tstop", "sgtf:tstop"]
covAB = vcov(m)["sgtf", "sgtf:tstop"]

haz = data.table(t = 0:28)
haz[, x.ct := A + B*t]
haz[, x.se := sqrt(varA + t^2 * varB + 2*t*covAB)]
haz[, x.hi := x.ct + 1.96 * x.se]
haz[, x.lo := x.ct - 1.96 * x.se]

haz[, c.ct := 0.255365]
haz[, c.hi := c.ct + 1.96 * 0.097451]
haz[, c.lo := c.ct - 1.96 * 0.097451]

ggplot(haz) +
    geom_ribbon(aes(x = t, ymin = exp(c.lo), ymax = exp(c.hi), fill = "Constant"), alpha = 0.4) +
    geom_ribbon(aes(x = t, ymin = exp(x.lo), ymax = exp(x.hi), fill = "Time-varying"), alpha = 0.4) +
    geom_line(aes(x = t, y = exp(c.ct), colour = "Constant")) +
    geom_line(aes(x = t, y = exp(x.ct), colour = "Time-varying")) +
    labs(x = "Days post specimen", y = "Relative hazard of death", fill = "Model", colour = "Model") +
    scale_x_continuous(breaks = c(0, 7, 14, 21, 28)) +
    theme(legend.position = c(0.05, 0.9)) +
    ylim(0, NA)
ggsave("./output/time_varying.pdf", width = 15, height = 10, units = "cm", useDingbats = FALSE)
ggsave("./output/time_varying.png", width = 15, height = 10, units = "cm")


