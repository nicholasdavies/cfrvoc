library(survival)
library(rms)
library(KMunicate)
library(qs)
library(ggplot2)
library(cowplot)
library(survminer)
library(survey)
library(aod)
library(ggstance)
library(stringr)

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
    cat(strata_names, ": ", n_strata, " strata, ", n_informative_strata, " informative strata containing ", n_deaths, " deaths; ", n_with2deaths, " informative strata with 2+ deaths.\n", sep = "")
}

# Perform both unweighted and weighted Cox analyses
do_cox = function(mdl_id, formula, data, strata_cols, cc_only = FALSE)
{
    if (cc_only) {
        model_cc  = do_cox_work(paste0(mdl_id, ".cc"), formula, data, strata_cols, use_weights = FALSE)
        invisible(list(model_cc = model_cc))
    } else {
        model_cc  = do_cox_work(paste0(mdl_id, ".cc"), formula, data, strata_cols, use_weights = FALSE)
        model_ipw = do_cox_work(paste0(mdl_id, ".ipw"), formula, data, strata_cols, use_weights = TRUE)
        invisible(list(model_cc = model_cc, model_ipw = model_ipw))
    }
}

# Workhorse for Cox PH survival analyses
do_cox_work = function(mdl_id, formula, data, strata_cols, use_weights = FALSE)
{
    if (length(strata_cols) > 0) {
        data[, stratum := do.call(paste, c(.SD, sep = "|")), .SDcols = strata_cols]
        data[, stratum := factor(stratum)]
    }
    
    if (use_weights == FALSE) {
        model = coxph(formula, data = data)
    } else {
        fstr = paste0(deparse(formula, width.cutoff = 500L), collapse = "")
        vecstr = function(x) if (length(x) <= 1) as.character(x) else paste0("c(", paste0(x, collapse = ","), ")")
        fstr = str_replace_all(fstr, "AK", vecstr(AK));
        fstr = str_replace_all(fstr, "IK", vecstr(IK));
        model = svycoxph(str2lang(fstr), data = data, svydesign(ids = ~1, weights = ~wt, data = data))
        # model = svycoxph(formula, data = data, svydesign(ids = ~1, weights = ~wt, data = data))
    }
    
    print(summary(model))

    # Extract effect sizes, SEs, CIs and P values
    summ = as.data.table(keep.rownames = TRUE, summary(model)$coefficients)
    if (use_weights == FALSE) {
        names(summ) = c("parameter", "coef", "HR", "se_coef", "z_score", "P")
        summ[, se_nonrobust := se_coef];
        setcolorder(summ, c("parameter", "coef", "HR", "se_nonrobust", "se_coef", "z_score", "P"));
    } else {
        names(summ) = c("parameter", "coef", "HR", "se_nonrobust", "se_coef", "z_score", "P")
    }
    summ[, HR.lo95 := exp(coef - 1.96 * se_coef)]
    summ[, HR.hi95 := exp(coef + 1.96 * se_coef)]
    summ[, model_id := mdl_id]
    summ[, weighted := use_weights]
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
            ms = rbind(ms, summ, fill = TRUE);
            fwrite(ms, "./output/model_summary.csv");
        } else {
            fwrite(summ, "./output/model_summary.csv");
        }
        
        # Save model
        qsave(model, paste0("./output/models/", make.names(mdl_id), ".qs"))
    }
    
    model
}

load_model = function(mdl_id)
{
    if ("coxph" %in% class(mdl_id)) {
        return (mdl_id)
    } else {
        return (qread(paste0("./output/models/", make.names(mdl_id), ".qs")))
    }
}

# Workhorse for Cox PH survival analyses with interaction terms
# marker: e.g. p_voc or sgtf
# group: e.g. age_group, must be categorical
# terms: other terms for the RHS of the model
# data: as in do_cox
# strata_cols: as in do_cox
# use_weights: as in do_cox
do_cox_interaction = function(mdl_id, marker, group, terms, data, strata_cols)
{
    group_levels = data[, unique(get(group))];
    marker_names = make.names(paste0(marker, group_levels));
    if (length(group_levels) > 25) {
        stop("More than 25 group levels.")
    }

    data2 = copy(data);
    for (g in seq_along(group_levels)) {
        data2[, (marker_names[g]) := ifelse(get(group) == group_levels[g], get(marker), 0)]
    }
    
    # Run both formulas and check for interaction by LRT / Wald
    formula_int = paste("Surv(time, status) ~", paste(c(marker_names, terms, "strata(stratum)"), collapse = " + "));
    formula_noint = paste("Surv(time, status) ~", paste(c(marker, terms, "strata(stratum)"), collapse = " + "));
    formula_checkipw = paste("Surv(time, status) ~", paste(c(marker, terms, paste(marker, group, sep = ":")), collapse = " + "));
    
    model_int = list()
    model_noint = list()
    
    model_int[[".cc"]]    = do_cox_work(paste0("int",   mdl_id, ".cc"),  as.formula(formula_int),      data2, strata_cols, use_weights = FALSE)
    model_int[[".ipw"]]   = do_cox_work(paste0("int",   mdl_id, ".ipw"), as.formula(formula_int),      data2, strata_cols, use_weights = TRUE)
    model_noint[[".cc"]]  = do_cox_work(paste0("noint", mdl_id, ".cc"),  as.formula(formula_noint),    data2, strata_cols, use_weights = FALSE)
    model_noint[[".ipw"]] = do_cox_work(paste0("noint", mdl_id, ".ipw"), as.formula(formula_noint),    data2, strata_cols, use_weights = TRUE)
    model_checkipw        = do_cox_work(paste0("test checkipw", mdl_id), as.formula(formula_checkipw), data2, strata_cols, use_weights = TRUE)
    
    # Save model to file
    if (!mdl_id %like% "^test") {
        for (suffix in c(".cc", ".ipw")) {
            mdl_id_s = paste0(mdl_id, suffix);
            
            # Send output to file
            sink(paste0("./output/interaction_", mdl_id_s, ".txt"));
            cat("With interaction: ", formula_int, "\n");
            print(summary(model_int[[suffix]]))
            cat("\n\nWithout interaction: ", formula_int, "\n");
            print(summary(model_noint[[suffix]]))
            
            if (suffix == ".cc") {
                cat("\n\nLikelihood ratio test:\n")
                print(anova(model_int[[suffix]], model_noint[[suffix]]))
            } else if (suffix == ".ipw") {
                cat("\n\nWald test:\n")
                print(waldt(model_checkipw, like = paste0(marker, ":")))
            }
            sink(NULL)
            
            # Save models
            qsave(list(model_int = model_int[[suffix]], model_noint = model_noint[[suffix]]), paste0("./output/models/", make.names(mdl_id_s), ".qs"))
        }
    }
}

# Do likelihood ratio test for two saved models
lrt = function(id_base, id_nested)
{
    model_base = load_model(id_base);
    model_nested = load_model(id_nested);
    anova(model_nested, model_base);
}

# Do joint Wald test for a saved model, specific terms
waldt = function(model_id, terms, like = NULL)
{
    model = load_model(model_id)
    cv = vcov(model)
    cf = coef(model)
    if (is.null(like)) {
        cx = which(names(cf) %in% terms)
    } else {
        cx = which(names(cf) %like% like)
    }
    cat("Joint Wald test for terms:\n")
    print(names(cf)[cx])
    wald.test(cv, cf, cx)
}

# Build missingness weights for a data set and save to outfile
build_weights = function(data, wtid)
{
    gc();
    
    outfile = paste0("./output/weights_", wtid, ".qs")

    data[, missing := ifelse(is.na(sgtf), 1, 0)]
    data[, specimen_week_f := factor(specimen_week)]
    data[, table(specimen_week)]
    
    # Estimate missingness model (takes a while to run)
    # Knots for age here are different than those used for the mortality model, as missingness is more naturally
    # spread across ages.
    missing_model = glm(missing ~ rcs(age, 5) + sex + rcs(imd, 3) + eth_cat + res_cat * asymptomatic + NHSER_name * specimen_week_f, 
        family = binomial(link = "cauchit"), data = data, control = list(trace = TRUE));
    data[, wt := 1 / (1 - predict(missing_model, newdata = .SD, type = "response"))];
    
    # Save missingness model weights
    model_weights = data[, .(person_id, wt)];
    qsave(model_weights, outfile);
    
    gc();
    
    return (invisible())
}

# Add weights from missingness model to data from model_data
weight_data = function(data, wtid)
{
    wtfile = paste0("./output/weights_", wtid, ".qs");
    model_weights = qread(wtfile);
    
    setkey(data, person_id);
    setkey(model_weights, person_id);
    data = merge(data, model_weights, by = "person_id", all.x = TRUE);
    if (data[, any(is.na(wt))]) {
        stop("Incomplete set of weights.")
    }
    return (data)
}








####################
# ESTIMATE WEIGHTS #
####################

# Load data
cd = complete_data("20210225")
#cd = reduced_data("20210225") # Use this line to load a reduced data set supplied with the repo.


# September 1st
dataW = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-09-01", keep_missing = TRUE)
build_weights(dataW, "sep")

# October 1st
dataW = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-10-01", keep_missing = TRUE)
build_weights(dataW, "oct")

# November 1st
dataW = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01", keep_missing = TRUE)
build_weights(dataW, "nov")

# December 1st
dataW = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-12-01", keep_missing = TRUE)
build_weights(dataW, "dec")

# January 1st
dataW = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2021-01-01", keep_missing = TRUE)
build_weights(dataW, "jan")

# LTLA specific start point
dataW = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, prevalence_cutoff = TRUE, keep_missing = TRUE)
build_weights(dataW, "ltla")

# No registration cutoff
dataW = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 0, P_voc = 0, date_min = "2020-11-01", keep_missing = TRUE)
build_weights(dataW, "nocut")

# Full followup only
dataW = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 0, P_voc = 0, date_min = "2020-11-01", date_max = as.character(ymd("2021-02-25") - 38), keep_missing = TRUE)
build_weights(dataW, "ff")



##########################
# SGTF SURVIVAL ANALYSES #
##########################

# 1. MAIN ANALYSES

# Knots for spline models
# Space age knots evenly from 20 to 100, as the default quantile-based knot locations are heavily clustered
# in younger age groups, who are more highly represented among individuals who get tested. We want more detail
# in the upper range of the age distribution, as mortality risk is higher here and small differences in the
# risk make a bigger difference on outcomes.
AK = c(20, 40, 60, 80, 100)
# Evenly spaced knots by quantile are appropriate for IMD.
IK = 3

# Assemble data set to be used for SGTF-based analyses
dataS = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS = weight_data(dataS, "nov");

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

# Linear age, linear IMD
do_cox("d28.SGTF.ll.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))

# Spline age, linear IMD
do_cox("d28.SGTF.sl.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))

# Linear age, spline IMD
do_cox("d28.SGTF.ls.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + age + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))

# Spline age, spline IMD
# This is the main SGTF analysis
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))


# LRT linear age, linear IMD -> spline age, linear IMD
lrt("d28.SGTF.ll.2020-11-01..r10.LTLA:date..cc", "d28.SGTF.sl.2020-11-01..r10.LTLA:date..cc")

# LRT linear age, linear IMD -> linear age, spline IMD
lrt("d28.SGTF.ll.2020-11-01..r10.LTLA:date..cc", "d28.SGTF.ls.2020-11-01..r10.LTLA:date..cc")

# LRT linear age, spline IMD -> spline age, spline IMD
lrt("d28.SGTF.ls.2020-11-01..r10.LTLA:date..cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date..cc")

# LRT spline age, linear IMD -> spline age, spline IMD
lrt("d28.SGTF.sl.2020-11-01..r10.LTLA:date..cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date..cc")

# Alternative stratifications
do_cox("d28.SGTF.ss.2020-11-01..r10.NHSE:week.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_week"))
do_cox("d28.SGTF.ss.2020-11-01..r10.UTLA:week.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS, c("UTLA_name", "specimen_week"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:week.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_week"))
do_cox("d28.SGTF.ss.2020-11-01..r10.NHSE:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.UTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS, c("UTLA_name", "specimen_date"))


# 2. INTERACTION TESTS
do_cox_interaction("age", "sgtf", "age_group",    c("rcs(age, AK)", "sex", "rcs(imd, IK)", "eth_cat", "res_cat"), dataS, c("LTLA_name", "specimen_date"))
do_cox_interaction("sex", "sgtf", "sex",          c("rcs(age, AK)", "sex", "rcs(imd, IK)", "eth_cat", "res_cat"), dataS, c("LTLA_name", "specimen_date"))
do_cox_interaction("imd", "sgtf", "imd_group",    c("rcs(age, AK)", "sex", "rcs(imd, IK)", "eth_cat", "res_cat"), dataS, c("LTLA_name", "specimen_date"))
do_cox_interaction("eth", "sgtf", "eth_cat",      c("rcs(age, AK)", "sex", "rcs(imd, IK)", "eth_cat", "res_cat"), dataS, c("LTLA_name", "specimen_date"))
do_cox_interaction("res", "sgtf", "res_cat",      c("rcs(age, AK)", "sex", "rcs(imd, IK)", "eth_cat", "res_cat"), dataS, c("LTLA_name", "specimen_date"))


# 3. CENSORING LENGTHS / DEATH TYPES
dataS_d07 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 7,  reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS_d14 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 14, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS_d21 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 21, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS_d28 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS_d60 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 60, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS_dNA = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = NA, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS_c28 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01", death_type = "cod")
dataS_e60 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = NA, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01", death_type = "60cod")

dataS_d07 = weight_data(dataS_d07, "nov")
dataS_d14 = weight_data(dataS_d14, "nov")
dataS_d21 = weight_data(dataS_d21, "nov")
dataS_d28 = weight_data(dataS_d28, "nov")
dataS_d60 = weight_data(dataS_d60, "nov")
dataS_dNA = weight_data(dataS_dNA, "nov")
dataS_c28 = weight_data(dataS_c28, "nov")
dataS_e60 = weight_data(dataS_e60, "nov")

do_cox("d07.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS_d07, c("LTLA_name", "specimen_date"))
do_cox("d14.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS_d14, c("LTLA_name", "specimen_date"))
do_cox("d21.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS_d21, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS_d28, c("LTLA_name", "specimen_date"))
do_cox("d60.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS_d60, c("LTLA_name", "specimen_date"))
do_cox("dNA.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS_dNA, c("LTLA_name", "specimen_date"))
do_cox("c28.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS_c28, c("LTLA_name", "specimen_date"))
do_cox("e60.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS_e60, c("LTLA_name", "specimen_date"))


# 4. NON-OVERLAPPING PERIODS
splitS = survSplit(Surv(time, status) ~ ., data = dataS, cut = c(0, 7, 14, 21, 28), start = "tstart", end = "tstop")
setDT(splitS)
splitS[, sgtf_w1 := ifelse(sgtf == 1 & tstart ==  0, 1, 0)]
splitS[, sgtf_w2 := ifelse(sgtf == 1 & tstart ==  7, 1, 0)]
splitS[, sgtf_w3 := ifelse(sgtf == 1 & tstart == 14, 1, 0)]
splitS[, sgtf_w4 := ifelse(sgtf == 1 & tstart == 21, 1, 0)]
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf_by_week", Surv(tstart, tstop, status) ~ sgtf_w1 + sgtf_w2 + sgtf_w3 + sgtf_w4 + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), splitS, c("LTLA_name", "specimen_date"))


# 5. TIME-X INTERACTION TERM
event_times = dataS[!is.na(death_date), sort(unique(as.numeric(death_date - specimen_date)))]
event_times = event_times[event_times <= 28]
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
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0sgtf:tstop", 
    Surv(tstart, tstop, status) ~ sgtf                            + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop               + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop+sgtf:tstop2", 
    Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + sgtf:tstop2 + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))

# LRT time interaction
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0sgtf:tstop.cc",  "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop.cc")
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop.cc",  "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop+sgtf:tstop2.cc")
waldt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop.ipw", "sgtf:tstop")

# 5b. TIME-AGE
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0age:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, AK)                          + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.age:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, AK) + age:tstop              + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.age:tstop+age:tstop2", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, AK) + age:tstop + age:tstop2 + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))

# LRT time interaction
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0age:tstop.cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.age:tstop.cc")
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.age:tstop.cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.age:tstop+age:tstop2.cc")
waldt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.age:tstop.ipw", "age:tstop")

# 5c. TIME-SEX
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0sex:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, AK) + sexM                            + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sex:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, AK) + sexM + sexM:tstop               + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sex:tstop+sex:tstop2", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, AK) + sexM + sexM:tstop + sexM:tstop2 + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))

# LRT time interaction
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0sex:tstop.cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sex:tstop.cc")
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sex:tstop.cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sex:tstop+sex:tstop2.cc")
waldt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sex:tstop.ipw", "sexM:tstop")

# 5d. TIME-IMD
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0imd:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK)                          + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.imd:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + imd:tstop              + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.imd:tstop+imd:tstop2", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + imd:tstop + imd:tstop2 + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))

# LRT time interaction
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0imd:tstop.cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.imd:tstop.cc")
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.imd:tstop.cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.imd:tstop+imd:tstop2.cc")
waldt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.imd:tstop.ipw", "imd:tstop")

# 5e. TIME-ETHNICITY
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0eth:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + ethA + ethB + ethO + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.eth:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + ethA + ethB + ethO + ethA:tstop + ethB:tstop + ethO:tstop + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.eth:tstop+eth:tstop2", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + ethA + ethB + ethO + ethA:tstop + ethB:tstop + ethO:tstop + ethA:tstop2 + ethB:tstop2 + ethO:tstop2 + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))

# LRT time interaction
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0eth:tstop.cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.eth:tstop.cc")
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.eth:tstop.cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.eth:tstop+eth:tstop2.cc")
waldt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.eth:tstop.ipw", like = "eth[A-Z]:tstop")

# 5f. TIME-RESIDENCE TYPE
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0res:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + resC + resO + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.res:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + resC + resO + resC:tstop + resO:tstop + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.res:tstop+res:tstop2", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + resC + resO + resC:tstop + resO:tstop + resC:tstop2 + resO:tstop2 + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))

# LRT time interaction
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0res:tstop.cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.res:tstop.cc")
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.res:tstop.cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.res:tstop+res:tstop2.cc")
waldt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.res:tstop.ipw", like = "res[A-Z]:tstop")


# 6. SENSITIVITY ANALYSES

# 6a. Cutoff by date range
dataS_sep = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-09-01")
dataS_oct = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-10-01")
dataS_nov = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS_dec = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-12-01")
dataS_jan = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2021-01-01")
dataS_sep = weight_data(dataS_sep, "sep")
dataS_oct = weight_data(dataS_oct, "oct")
dataS_nov = weight_data(dataS_nov, "nov")
dataS_dec = weight_data(dataS_dec, "dec")
dataS_jan = weight_data(dataS_jan, "jan")
do_cox("d28.SGTF.ss.2020-09-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS_sep, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-10-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS_oct, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS_nov, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-12-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS_dec, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2021-01-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS_jan, c("LTLA_name", "specimen_date"))

dataS_LTLAprev = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, prevalence_cutoff = TRUE)
dataS_LTLAprev = weight_data(dataS_LTLAprev, "ltla")
do_cox("d28.SGTF.ss.LTLA..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataS_LTLAprev, c("LTLA_name", "specimen_date"))

# 6b. No registration cutoff
dataSc = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 0, P_voc = 0, date_min = "2020-11-01")
dataSc = weight_data(dataSc, "nocut")
do_cox("d28.SGTF.ss.2020-11-01..rNA.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataSc, c("LTLA_name", "specimen_date"))

# 6c. Include individuals with full followup only
dataSd = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 0, P_voc = 0, 
    date_min = "2020-11-01", date_max = as.character(ymd("2021-02-25") - 38))
dataSd = weight_data(dataSd, "ff")
do_cox("d28.SGTF.ss.2020-11-01.tminus38.rNA.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataSd, c("LTLA_name", "specimen_date"))

# 6d. Include asymptomatic indicator
dataS[, asymptomatic := factor(asymptomatic, levels = c("N", "Y"))]
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.asymptomatic", Surv(time, status) ~ sgtf + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + asymptomatic + strata(stratum), dataS, c("LTLA_name", "specimen_date"))



###########################
# P_VOC SURVIVAL ANALYSES #
###########################

# 1. MAIN ANALYSES

# 1a. WHOLE DATA RANGE
dataP = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0)
dataP = weight_data(dataP, "sep")
do_cox("d28.pVOC.ss.2020-09-01..r10.LTLA:date.", Surv(time, status) ~ p_voc + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataP, c("LTLA_name", "specimen_date"))

# 1b. POST-NOV 1 ONLY
dataPn = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataPn = weight_data(dataPn, "nov")
do_cox("d28.pVOC.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ p_voc + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataPn, c("LTLA_name", "specimen_date"))

# 1c. 60 DAY FOLLOWUP P_VOC
dataP60 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 60, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataP60 = weight_data(dataP60, "nov")
do_cox("d60.pVOC.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ p_voc + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataP60, c("LTLA_name", "specimen_date"))

# 1d. P_VOC SENSITIVITY - SEQUENCING
cds = complete_data("20210225", sgtfv_file = "./sgtf_voc_sequencing.csv")
dataPnseq = model_data(cds, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataPnseq = weight_data(dataPnseq, "nov")
do_cox("d28.pVOC2.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ p_voc + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), dataPnseq, c("LTLA_name", "specimen_date"))

# 2. TIME-P_VOC
event_times = dataPn[!is.na(death_date), sort(unique(as.numeric(death_date - specimen_date)))]
event_times = event_times[event_times <= 28]
dataPT = survSplit(Surv(time, status) ~ ., data = dataPn, cut = event_times, start = "tstart", end = "tstop")
setDT(dataPT)
dataPT[, tstop2 := tstop^2]

do_cox("d28.pVOC.ss.2020-11-01..r10.LTLA:date.0p_voc:tstop", 
    Surv(tstart, tstop, status) ~ p_voc                              + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), 
    dataPT, c("LTLA_name", "specimen_date"))
do_cox("d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc:tstop", 
    Surv(tstart, tstop, status) ~ p_voc + p_voc:tstop                + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), 
    dataPT, c("LTLA_name", "specimen_date"))
do_cox("d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc:tstop+p_voc:tstop2", 
    Surv(tstart, tstop, status) ~ p_voc + p_voc:tstop + p_voc:tstop2 + rcs(age, AK) + sex + rcs(imd, IK) + eth_cat + res_cat + strata(stratum), 
    dataPT, c("LTLA_name", "specimen_date"), cc_only = TRUE)

# LRT time interaction
lrt("d28.pVOC.ss.2020-11-01..r10.LTLA:date.0p_voc:tstop.cc", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc:tstop.cc")
lrt("d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc:tstop.cc", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc:tstop+p_voc:tstop2.cc")
waldt("d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc:tstop.ipw", "p_voc:tstop")




#########
# PLOTS #
#########

# Hazard ratio over time
plot_hazard = function(marker, constant_model, linear_varying_model, ylimits, y_label)
{
    cat("Loading varying model...")
    mv = load_model(linear_varying_model)
    cat("OK\n")
    
    marker_int = paste0(marker, ":tstop")
    
    A = coefficients(mv)[[marker]]
    B = coefficients(mv)[[marker_int]]
    varA = vcov(mv)[marker, marker]
    varB = vcov(mv)[marker_int, marker_int]
    covAB = vcov(mv)[marker, marker_int]
    rm(mv)
    gc()
    
    haz = data.table(t = 0:28)
    haz[, x.ct := A + B*t]
    haz[, x.se := sqrt(varA + t^2 * varB + 2*t*covAB)]
    haz[, x.hi := x.ct + 1.96 * x.se]
    haz[, x.lo := x.ct - 1.96 * x.se]

    cat("Loading constant model...")
    mc = load_model(constant_model)
    cat("OK\n")
    
    haz[, c.ct := coefficients(mc)[[marker]]]
    haz[, c.hi := c.ct + 1.96 * summary(mc)$coefficients[marker, "se(coef)"]]
    haz[, c.lo := c.ct - 1.96 * summary(mc)$coefficients[marker, "se(coef)"]]
    rm(mc)
    gc()
    
    ggplot(haz) +
        geom_ribbon(aes(x = t, ymin = exp(c.lo), ymax = exp(c.hi), fill = "Constant"), alpha = 0.4) +
        geom_ribbon(aes(x = t, ymin = exp(x.lo), ymax = exp(x.hi), fill = "Time-varying"), alpha = 0.4) +
        geom_line(aes(x = t, y = exp(c.ct), colour = "Constant")) +
        geom_line(aes(x = t, y = exp(x.ct), colour = "Time-varying")) +
        geom_hline(aes(yintercept = 1), linetype = "33", size = 0.25) +
        labs(x = "Days since positive test", y = y_label, fill = NULL, colour = NULL) +
        scale_x_continuous(breaks = c(0, 7, 14, 21, 28)) +
        scale_colour_manual(aesthetics = c("colour", "fill"), values = c("darkorchid", "#44aa88")) +
        theme(legend.position = c(0.05, 0.9)) +
        ylim(ylimits)
    
    # Uncomment if want to extract time-varying hazards
    # return (haz)
}

# haz = plot_hazard("sgtf", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.0sgtf:tstop.cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop.cc", c(0, 4), "Hazard ratio for SGTF")
# haz[t == 1, .(exp(x.ct), exp(x.lo), exp(x.hi))]
# haz[t == 14, .(exp(x.ct), exp(x.lo), exp(x.hi))]
# haz[t == 28, .(exp(x.ct), exp(x.lo), exp(x.hi))]

theme_set(theme_cowplot(font_size = 11))

hp_sgtf_cc = plot_hazard("sgtf", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.0sgtf:tstop.cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop.cc", 
    c(0, 3), "Hazard ratio")
hp_voc_cc = plot_hazard("p_voc", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.0p_voc:tstop.cc", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc:tstop.cc", 
    c(0, 3), NULL) + theme(legend.position = "none")
hp_sgtf_ipw = plot_hazard("sgtf", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.0sgtf:tstop.ipw", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop.ipw", 
    c(0, 3), NULL) + theme(legend.position = "none")
hp_voc_ipw = plot_hazard("p_voc", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.0p_voc:tstop.ipw", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc:tstop.ipw", 
    c(0, 3), NULL) + theme(legend.position = "none")

#qsave(list(hp_sgtf_cc, hp_voc_cc, hp_sgtf_ipw), "plots_temp.qs")
#qsave(hp_voc_ipw, "plots_temp2.qs")
pt = qread("plots_temp.qs")
hp_sgtf_cc = pt[[1]]
hp_voc_cc = pt[[2]]
hp_sgtf_ipw = pt[[3]]
hp_voc_ipw = qread("plots_temp2.qs")

# Effects of SGTF
summaries = fread("./output/model_summary.csv")
summaries_int = summaries[str_count(model_id, "\\.") != 8]
fwrite(
    summaries_int[parameter %like% "^sgtf.", .(model_id, parameter, HR, HR.lo95, HR.hi95, P)],
    "./output/table_subgroup_effects.csv"
)

summaries = summaries[str_count(model_id, "\\.") == 8] 
summaries[, model_id := str_replace_all(model_id, "\\.ll\\.", ".L.L.")]
summaries[, model_id := str_replace_all(model_id, "\\.sl\\.", ".S.L.")]
summaries[, model_id := str_replace_all(model_id, "\\.ls\\.", ".L.S.")]
summaries[, model_id := str_replace_all(model_id, "\\.ss\\.", ".S.S.")]
summaries[, model_id := str_replace_all(model_id, "\\.r10\\.", ".10.")]
summaries[, model_id := str_replace_all(model_id, "\\.rNA\\.", ".0.")]

summaries[model_id == "d28.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.cc" & parameter == "sgtf_w1",
    model_id := "d00-07.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.cc"]
summaries[model_id == "d28.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.cc" & parameter == "sgtf_w2",
    model_id := "d08-14.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.cc"]
summaries[model_id == "d28.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.cc" & parameter == "sgtf_w3",
    model_id := "d15-21.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.cc"]
summaries[model_id == "d28.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.cc" & parameter == "sgtf_w4",
    model_id := "d22-28.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.cc"]

summaries[model_id == "d28.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.ipw" & parameter == "sgtf_w1",
    model_id := "d00-07.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.ipw"]
summaries[model_id == "d28.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.ipw" & parameter == "sgtf_w2",
    model_id := "d08-14.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.ipw"]
summaries[model_id == "d28.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.ipw" & parameter == "sgtf_w3",
    model_id := "d15-21.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.ipw"]
summaries[model_id == "d28.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.ipw" & parameter == "sgtf_w4",
    model_id := "d22-28.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.ipw"]

summaries = summaries[model_id != "d28.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.cc" &
        model_id != "d28.SGTF.S.S.2020-11-01..10.LTLA:date.sgtf_by_week.ipw"]

summaries[parameter %like% "^sgtf_w[1-4]$", parameter := "sgtf"]

summaries[, c("m_deaths", "m_marker", "m_age", "m_imd", "m_startdate", "m_enddate", "m_reg_cutoff", "m_strata", "m_xvars", "m_weighting") := tstrsplit(model_id, "\\.")]
summaries[m_enddate == "tminus38", m_enddate := as.character(ymd("2021-02-25") - 38)]
summaries = summaries[order(model_id)]

# Sensitivity analyses
sensitivity_groups = list(
    `Death type` = list(
        pattern = "d*/SGTF/S/S/2020-11-01//10/LTLA:date//*",
        highlight = 1
    ),
    `.d2` = list(
        pattern = "c28/SGTF/S/S/2020-11-01//10/LTLA:date//*",
        highlight = 1
    ),
    `.d3` = list(
        pattern = "e60/SGTF/S/S/2020-11-01//10/LTLA:date//*",
        highlight = 1
    ),
    
    `Misclassification adjustment` = list(
        pattern = "*/pVOC/*/*/*/*/*/*//*",
        highlight = c(1, 2, 5)
    ),
    `.pvoc2` = list(
        pattern = "*/pVOC2/*/*/*/*/*/*//*",
        highlight = c(1, 2, 5)
    ),
    
    `Geographical and temporal stratification` = list(
        pattern = "d28/SGTF/S/S/2020-11-01//10/NHSE:week//*",
        highlight = 8
    ),
    `.gt2` = list(
        pattern = "d28/SGTF/S/S/2020-11-01//10/UTLA:week//*",
        highlight = 8
    ),
    `.gt3` = list(
        pattern = "d28/SGTF/S/S/2020-11-01//10/LTLA:week//*",
        highlight = 8
    ),
    `.gt4` = list(
        pattern = "d28/SGTF/S/S/2020-11-01//10/NHSE:date//*",
        highlight = 8
    ),
    `.gt5` = list(
        pattern = "d28/SGTF/S/S/2020-11-01//10/UTLA:date//*",
        highlight = 8
    ),
    `.gt6` = list(
        pattern = "d28/SGTF/S/S/2020-11-01//10/LTLA:date//*",
        highlight = 8
    ),
    
    `Age and IMD terms (linear vs. spline)` = list(
        pattern = "d28/SGTF/*/*/2020-11-01//10/LTLA:date//*",
        highlight = c(3, 4)
    ),
    
    `By week since positive test` = list(
        pattern = "*/SGTF/S/S/2020-11-01//10/LTLA:date/sgtf_by_week/*",
        highlight = 1
    ),
    
    `Analysis start date` = list(
        pattern = "d28/SGTF/S/S/*//10/LTLA:date//*",
        highlight = 5
    ),
    
    `Covariate interactions with time since positive test` = list(
        pattern = "d28/SGTF/S/S/2020-11-01//10/LTLA:date/age:tstop/*",
        highlight = 9
    ),
    `.sex` = list(
        pattern = "d28/SGTF/S/S/2020-11-01//10/LTLA:date/sex:tstop/*",
        highlight = 9
    ),
    `.imd` = list(
        pattern = "d28/SGTF/S/S/2020-11-01//10/LTLA:date/imd:tstop/*",
        highlight = 9
    ),
    `.eth` = list(
        pattern = "d28/SGTF/S/S/2020-11-01//10/LTLA:date/eth:tstop/*",
        highlight = 9
    ),
    `.res` = list(
        pattern = "d28/SGTF/S/S/2020-11-01//10/LTLA:date/res:tstop/*",
        highlight = 9
    ),

    `No registration cutoff` = list(
        pattern = "d28/SGTF/S/S/2020-11-01//0/LTLA:date//*",
        highlight = 7
    ),
    
    `Subjects with full 28-day follow-up only` = list(
        pattern = "d28/SGTF/S/S/2020-11-01/tminus38/0/LTLA:date//*",
        highlight = c(6, 7)
    ),
    
    `Asymptomatic screening indicator included as covariate` = list(
        pattern = "d28/SGTF/S/S/2020-11-01//10/LTLA:date/asymptomatic/*",
        highlight = 9
    )
)


theme_set(theme_cowplot(font_size = 11))

pl_sensitivity = plot_forest(summaries, sensitivity_groups, 13, -2.3, 0.55, x_widths = c(1, 1.6, 1.1, 0.5, 1.6, 2.6, 1.8, 1.8, 2.8))
ggsave("./output/S1-hazard-sensitivity.png", pl_sensitivity, width = 20, height = 25, units = "cm")


plot_forest = function(summaries, groups, left_margin, tab_from, tab_to, x_widths = c(1, 1.5, 1, 1, 2, 2, 1, 2, 3))
{
    interpret = function(str) {
        str = str_replace_all(str, "\\*", "[^.]*");
        str = str_replace_all(str, "\\/", "[.]");
        str
    }
    
    # Assemble data for plotting
    table_cols = c("m_deaths", "m_marker", "m_age", "m_imd", "m_startdate", "m_enddate", "m_reg_cutoff", "m_strata", "m_xvars", "m_weighting");
    pd = data.table()
    headings = data.table()
    summary_table = data.table()
    y_index = 1;
    for (g in seq_along(groups)) {
        if (!names(groups)[g] %like% "^\\.") {
            headings = rbind(headings, 
                data.table(label = names(groups)[g], y_index)
            )
            y_index = y_index + 1;
            
            summary_table = rbind(summary_table, data.table(analysis = names(groups)[g]), fill = TRUE)
        }
        model_ids = summaries[model_id %like% interpret(groups[[g]]$pattern), unique(model_id)];
        model_ids_genus = unique(str_remove_all(model_ids, "\\.cc$|\\.ipw$"))
        
        for (model in model_ids_genus) {
            entry = summaries[(model_id == model | model_id == paste0(model, ".cc") | model_id == paste0(model, ".ipw")) & parameter %like% "^sgtf$|^p_voc$", 
                    .(parameter, h = HR, h0 = HR.lo95, h1 = HR.hi95, 
                      blank = "", highlight = paste(table_cols[groups[[g]]$highlight], collapse = "|"), y_index,
                      m_deaths, m_marker, m_age, m_imd, m_startdate, m_enddate, m_reg_cutoff, m_strata, m_xvars, m_weighting)]
            
            pd = rbind(pd, entry)
            summary_table = rbind(summary_table, entry, use.names = TRUE, fill = TRUE)

            y_index = y_index + 1;
        }
    }
    
    # Make labels
    blank_col = which(names(pd) == "blank");
    labels = melt(pd[m_weighting == "cc", .SD, .SDcols = (blank_col + 1):(ncol(pd) - 1)], id.vars = c("y_index", "highlight"))
    labels[value == "", value := "—"]
    labels[, x_index := as.numeric(variable)]
    x_shifts = cumsum(x_widths) - x_widths[1] * 0.9
    labels[, x_pos := tab_from + (tab_to - tab_from) * x_shifts[x_index] / max(x_shifts)]
    labels[, colour := ifelse(variable %like% highlight, "#000000", "#888888"), by = .(y_index, x_index)]
    
    fwrite(
        summary_table,
        "./output/table_sensitivity_effects.csv"
    )
    
    # Make table structure
    tabr = data.table(xmin = tab_from - 0.14, xmax = 2.6, 
        ymin = labels[, unique(y_index) - 0.5], ymax = labels[, unique(y_index) + 0.5])
    tabr[, fill := rep_len(c("#f0f0f0", "#f8f8f8"), .N)]
    tabt = data.table(label = c("Death\ntype", "Marker", "Age", "IMD", "Start\ndate", "End\ndate", 
        "Registr.\ncutoff", "Strati-\nfication", "Extra\nterms"), x = labels[, unique(x_pos)], y = y_index)
    
    pd[, weighting := ifelse(m_weighting == "cc", "Complete\ncases", "IPW")]
    pd[, weighting := factor(weighting, c("IPW", "Complete\ncases"))]

    ggplot(pd) +
        annotate("rect", xmin = tabr$xmin, xmax = tabr$xmax, ymin = tabr$ymin, ymax = tabr$ymax, fill = tabr$fill) +
        geom_point(aes(x = h, y = y_index, colour = weighting), size = 0.5, position = position_dodgev(height = 0.8)) +
        geom_linerangeh(aes(xmin = h0, xmax = h1, y = y_index, colour = weighting), size = 0.5, position = position_dodgev(height = 0.8)) +
        scale_y_reverse(breaks = pd$y_index, labels = pd$blank) +
        geom_vline(aes(xintercept = 1.0), size = 0.125) +
        geom_vline(aes(xintercept = 1.5), size = 0.125, colour = "#aaaaaa") +
        geom_vline(aes(xintercept = 2.0), size = 0.125, colour = "#aaaaaa") +
        geom_vline(aes(xintercept = 2.5), size = 0.125, colour = "#aaaaaa") +
        scale_colour_manual(values = c("Complete\ncases" = "#bbaa88", "IPW" = "#9999bb"), 
            guide = guide_legend(reverse = TRUE, override.aes = list(fatten = 0.1))) +
        labs(x = "Hazard ratio", y = NULL, colour = NULL) +
        theme(plot.margin = unit(c(0.3, 0.2, 0.5, left_margin), "cm"),
            legend.position = c(0.625, 0.95), 
            legend.background = element_rect(fill = "white", colour = "black", size = 0.25),
            legend.margin = margin(-0.1, 0.3, 0.1, 0.1, "cm")) +
        annotate("text", x = tabt$x, y = tabt$y, label = tabt$label, hjust = 1, vjust = 0.5, angle = 90, lineheight = 0.75, size = 11/ggplot2:::.pt) +
        annotate("text", x = labels$x_pos, y = labels$y_index, label = labels$value, colour = labels$colour, size = 11/ggplot2:::.pt) +
        annotate("text", x = mean(c(tab_from, tab_to)), y = headings$y_index, label = headings$label, size = 11/ggplot2:::.pt) +
        coord_cartesian(xlim = c(0.9, 2.6), clip = 'off', expand = FALSE)
}




# Exploratory data plots
dataS[, SGTF_label := factor(ifelse(sgtf == 1, "SGTF", "Non-SGTF"), levels = c("SGTF", "Non-SGTF"))]
plot_samples = ggplot(dataS[, .N, keyby = .(specimen_date, SGTF_label)]) +
    geom_col(aes(specimen_date, N, fill = SGTF_label), position = "stack") +
    labs(x = "Specimen date", y = "Samples", fill = NULL) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = ymd(c("2020-10-31", "2021-02-15")), expand = expansion(0)) +
    scale_y_continuous(expand = expansion(0)) +
    theme(legend.position = c(0.03, 0.9))

plot_deaths = ggplot(dataS[died == TRUE, .N, keyby = .(specimen_date, SGTF_label)]) +
    geom_col(aes(specimen_date, N, fill = SGTF_label), position = "stack") +
    labs(x = "Specimen date", y = "Deaths", fill = NULL) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = ymd(c("2020-10-31", "2021-02-15")), expand = expansion(0)) +
    scale_y_continuous(expand = expansion(0)) +
    theme(legend.position = "none")


tbl_out = summaries[parameter %in% c("sgtf", "p_voc") & !model_id %like% "Time.{0,2}-(SGTF|p_voc)", 
    .(`Model` = model_id, `Parameter` = parameter, `Hazard ratio` = paste0(round(HR, 2), " (", round(HR.lo95, 2), "–", round(HR.hi95, 2), ")"),
        `P value` = ifelse(P < 0.001, "< 0.001", round(P, 3)))]

fwrite(tbl_out, "./output/table_effects.csv")

# Test proportional hazards assumption
dataS = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS[, stratum := do.call(paste, c(.SD, sep = "|")), .SDcols = c("LTLA_name", "specimen_date")]
dataS[, stratum := factor(stratum)]
model = coxph(Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), data = dataS)
resid = residuals(model, type = "schoenfeld")
resid_names = colnames(resid)
resid = residuals(model, type = "scaledsch")
times = as.numeric(rownames(resid))
resid = as.data.table(resid)
names(resid) = make.names(resid_names)
resid[, time := times]
cox.zph(model)

res_plot2 = function(resid, var, varname)
{
    ggplot(resid) +
        geom_jitter(aes_string(x = "time", y = var), colour = "#ff0000", alpha = 0.3, size = 0.5, width = 0.25) +
        geom_smooth(aes_string(x = "time", y = var), method = "loess", alpha = 0.7, colour = "blue", size = 0.5, fill = "#aaaaaa") +
        geom_hline(aes(yintercept = 0), size = 0.25) +
        labs(x = "Time", y = paste0("Residuals for ", varname))
}

theme_set(theme_cowplot(font_size = 11) + theme(strip.background = element_blank()))

pl1 = res_plot2(resid, "sgtf", "SGTF")
pl2 = res_plot2(resid, "age", "age")
pl3 = res_plot2(resid, "sexMale", "sex\n(Male)")
pl4 = res_plot2(resid, "imd", "IMD")
pl5 = res_plot2(resid, "eth_catA", "ethnicity\n(Asian)")
pl6 = res_plot2(resid, "eth_catB", "ethnicity\n(Black)")
pl7 = res_plot2(resid, "eth_catO", "ethnicity\n(Other/Mixed/Unknown)")
pl8 = res_plot2(resid, "res_catCare.Nursing.home", "residence type\n(Care/Nursing home)")
pl9 = res_plot2(resid, "res_catOther.Unknown", "residence type\n(Other/Unknown)")

pl = cowplot::plot_grid(pl1, pl2, pl3, pl4, pl5, pl6, pl7, pl8, pl9, labels = letters, label_size = 10, nrow = 3)
ggsave("./output/schoenfeld.pdf", width = 20, height = 20, units = "cm", useDingbats = FALSE)
ggsave("./output/schoenfeld.png", width = 20, height = 20, units = "cm")
write_xlsx(list(Residuals = resid), "./manuscript/sdE_schoenfeld.xlsx")
