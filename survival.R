library(survival)
library(rms)
library(KMunicate)
library(qs)
library(ggplot2)
library(cowplot)
library(survminer)

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
do_cox = function(mdl_id, formula, data, strata_cols)
{
    model_cc  = do_cox_work(paste0(mdl_id, ".cc"), formula, data, strata_cols, use_weights = FALSE)
    model_ipw = do_cox_work(paste0(mdl_id, ".ipw"), formula, data, strata_cols, use_weights = TRUE)
    invisible(list(model_cc = model_cc, model_ipw = model_ipw))
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
        wt = data$wt
        model = coxph(formula, data = data, weights = wt)
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
    if (class(mdl_id) == "coxph") {
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
    
    # Run both formulas and check for interaction LRT
    formula_int = paste("Surv(time, status) ~", paste(c(marker_names, terms, "strata(stratum)"), collapse = " + "));
    formula_noint = paste("Surv(time, status) ~", paste(c(marker, terms, "strata(stratum)"), collapse = " + "));
    
    model_int = list()
    model_noint = list()
    
    model_int[[".cc"]]    = do_cox_work(paste0(mdl_id, ".cc"),  as.formula(formula_int),   data2, strata_cols, use_weights = FALSE)
    model_int[[".ipw"]]   = do_cox_work(paste0(mdl_id, ".ipw"), as.formula(formula_int),   data2, strata_cols, use_weights = TRUE)
    model_noint[[".cc"]]  = do_cox_work(paste0(mdl_id, ".cc"),  as.formula(formula_noint), data2, strata_cols, use_weights = FALSE)
    model_noint[[".ipw"]] = do_cox_work(paste0(mdl_id, ".ipw"), as.formula(formula_noint), data2, strata_cols, use_weights = TRUE)
    
    # Write to summary file
    if (!mdl_id %like% "^test") {
        for (suffix in c(".cc", ".ipw")) {
            mdl_id_s = paste0(mdl_id, suffix);
            
            # Send output to file
            sink(paste0("./output/interaction_", mdl_id_s, ".txt"));
            cat("With interaction: ", formula_int, "\n");
            print(summary(model_int[[suffix]]))
            cat("\n\nWithout interaction: ", formula_int, "\n");
            print(summary(model_noint[[suffix]]))
            cat("\n\nLikelihood ratio test:\n")
            print(anova(model_int[[suffix]], model_noint[[suffix]]))
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




####################
# ESTIMATE WEIGHTS #
####################

cd = complete_data("20210205")
dataW = model_data(cd, "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 0, P_voc = 0, 
    date_min = "2020-09-01", date_max = "2100-01-01", prevalence_cutoff = FALSE, sgtfv_cutoff = 0, keep_missing = TRUE, death_type = "all")

# Load data
dataW[, missing := ifelse(is.na(sgtf), 1, 0)]
dataW[, specimen_week_f := factor(specimen_week)]

missing_model = glm(missing ~ rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat * asymptomatic + NHSER_name * specimen_week_f, 
    family = binomial(link = "cauchit"), data = dataW, control = list(trace = TRUE));

weight_data = function(data, missing_model)
{
    data[, specimen_week_f := factor(specimen_week)]
    data[, wt := 1 / (1 - predict(missing_model, newdata = .SD, type = "response"))];
}

# qsave(missing_model, "./output/missing_model.qs")
# missing_model = qread("./output/missing_model.qs")

##########################
# SGTF SURVIVAL ANALYSES #
##########################

# 1. MAIN ANALYSES

# Assemble data set to be used for SGTF-based analyses
dataS = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
weight_data(dataS, missing_model);

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
# 
# dataREDUCED = model_data(rd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
# do_cox("testSGTF + lin age + lin IMD | LTLA + spec date", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataREDUCED, c("LTLA_name", "specimen_date"))
# do_cox("testSGTF + spl age + spl IMD | NHSE + spec week", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_week"))
# do_cox("testSGTF + spl age + spl IMD | NHSE + spec week", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataREDUCED, c("NHSER_name", "specimen_week"))
# 
# dataP = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0)
# dataPR = model_data(rd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0)
# 
# # Linear
# do_cox("test p_voc + lin age + lin IMD | LTLA + spec date", Surv(time, status) ~ p_voc + age + sex + imd + eth_cat + res_cat + strata(stratum), dataPR, c("LTLA_name", "specimen_date"))
# load_model("p_voc + lin age + lin IMD | LTLA + spec date")

# Linear age, linear IMD
do_cox("d28.SGTF.ll.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))

# Spline age, linear IMD
do_cox("d28.SGTF.sl.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))

# Linear age, spline IMD
do_cox("d28.SGTF.ls.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + age + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))

# Spline age, spline IMD
# This is the main SGTF analysis
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))

# LRT linear age, linear IMD -> spline age, linear IMD
lrt("d28.SGTF.ll.2020-11-01..r10.LTLA:date..cc", "d28.SGTF.sl.2020-11-01..r10.LTLA:date..cc")
lrt("d28.SGTF.ll.2020-11-01..r10.LTLA:date..ipw", "d28.SGTF.sl.2020-11-01..r10.LTLA:date..ipw")

# LRT linear age, linear IMD -> linear age, spline IMD
lrt("d28.SGTF.ll.2020-11-01..r10.LTLA:date..cc", "d28.SGTF.ls.2020-11-01..r10.LTLA:date..cc")
lrt("d28.SGTF.ll.2020-11-01..r10.LTLA:date..ipw", "d28.SGTF.ls.2020-11-01..r10.LTLA:date..ipw")

# LRT linear age, spline IMD -> spline age, spline IMD
lrt("d28.SGTF.ls.2020-11-01..r10.LTLA:date..cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date..cc")
lrt("d28.SGTF.ls.2020-11-01..r10.LTLA:date..ipw", "d28.SGTF.ss.2020-11-01..r10.LTLA:date..ipw")

# LRT spline age, linear IMD -> spline age, spline IMD
lrt("d28.SGTF.sl.2020-11-01..r10.LTLA:date..cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date..cc")
lrt("d28.SGTF.sl.2020-11-01..r10.LTLA:date..ipw", "d28.SGTF.ss.2020-11-01..r10.LTLA:date..ipw")

# Alternative stratifications
do_cox("d28.SGTF.ss.2020-11-01..r10.NHSE:week.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_week"))
do_cox("d28.SGTF.ss.2020-11-01..r10.UTLA:week.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("UTLA_name", "specimen_week"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:week.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_week"))
do_cox("d28.SGTF.ss.2020-11-01..r10.NHSE:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.UTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("UTLA_name", "specimen_date"))



# 2. INTERACTION TESTS
do_cox_interaction("age", "sgtf", "age_group",    c("age_group",        "sex", "rcs(imd, nk = 3)", "eth_cat", "res_cat"), dataS, c("LTLA_name", "specimen_date"))
do_cox_interaction("sex", "sgtf", "sex",          c("rcs(age, nk = 3)", "sex", "rcs(imd, nk = 3)", "eth_cat", "res_cat"), dataS, c("LTLA_name", "specimen_date"))
do_cox_interaction("imd", "sgtf", "imd_group",    c("rcs(age, nk = 3)", "sex", "imd_group",        "eth_cat", "res_cat"), dataS, c("LTLA_name", "specimen_date"))
do_cox_interaction("eth", "sgtf", "eth_cat",      c("rcs(age, nk = 3)", "sex", "rcs(imd, nk = 3)", "eth_cat", "res_cat"), dataS, c("LTLA_name", "specimen_date"))
do_cox_interaction("res", "sgtf", "res_cat",      c("rcs(age, nk = 3)", "sex", "rcs(imd, nk = 3)", "eth_cat", "res_cat"), dataS, c("LTLA_name", "specimen_date"))
do_cox_interaction("asy", "sgtf", "asymptomatic", c("rcs(age, nk = 3)", "sex", "rcs(imd, nk = 3)", "eth_cat", "res_cat", "asymptomatic"), dataS, c("LTLA_name", "specimen_date"))


# age - alternative
int_0age = do_cox("test", Surv(time, status) ~ sgtf +
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date")) 
int_1age = do_cox("test", Surv(time, status) ~ sgtf * rcs(age, nk = 3) +
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date")) 

lrt(int_0age$model_cc,  int_1age$model_cc)  # **
lrt(int_0age$model_ipw, int_1age$model_ipw) # ***

# Investigating significant SGTF x residence type term -- does this explain other significant interactions for IPW model?

# age
int_res_0age = do_cox("test", Surv(time, status) ~ sgtf * res_cat +
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date")) 
int_res_1age = do_cox("test", Surv(time, status) ~ sgtf * (res_cat + rcs(age, nk = 3)) +
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date")) 

lrt(int_res_0age$model_cc,  int_res_1age$model_cc)  # **
lrt(int_res_0age$model_ipw, int_res_1age$model_ipw) # ***

# sex
int_res_0sex = do_cox("test", Surv(time, status) ~ sgtf * res_cat +
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date")) 
int_res_1sex = do_cox("test", Surv(time, status) ~ sgtf * (res_cat + sex) +
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date")) 

lrt(int_res_0sex$model_cc,  int_res_1sex$model_cc)
lrt(int_res_0sex$model_ipw, int_res_1sex$model_ipw)

# imd
int_res_0imd = do_cox("test", Surv(time, status) ~ sgtf * res_cat +
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date")) 
int_res_1imd = do_cox("test", Surv(time, status) ~ sgtf * (res_cat + rcs(imd, nk = 3)) +
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date")) 

lrt(int_res_0imd$model_cc,  int_res_1imd$model_cc)
lrt(int_res_0imd$model_ipw, int_res_1imd$model_ipw)

# eth
int_res_0eth = do_cox("test", Surv(time, status) ~ sgtf * res_cat +
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date")) 
int_res_1eth = do_cox("test", Surv(time, status) ~ sgtf * (res_cat + eth_cat) +
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date")) 

lrt(int_res_0eth$model_cc,  int_res_1eth$model_cc)
lrt(int_res_0eth$model_ipw, int_res_1eth$model_ipw)

# TODO interaction tests in the model with SGTF:time interaction.

# 3. CENSORING LENGTHS / DEATH TYPES
dataS_d07 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 7,  reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS_d14 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 14, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS_d21 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 21, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS_d28 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS_d60 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 60, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS_dNA = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = NA, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS_c28 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01", death_type = "cod")
dataS_e60 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = NA, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01", death_type = "60cod")

weight_data(dataS_d07, missing_model)
weight_data(dataS_d14, missing_model)
weight_data(dataS_d21, missing_model)
weight_data(dataS_d28, missing_model)
weight_data(dataS_d60, missing_model)
weight_data(dataS_dNA, missing_model)
weight_data(dataS_c28, missing_model)
weight_data(dataS_e60, missing_model)

do_cox("d07.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_d07, c("LTLA_name", "specimen_date"))
do_cox("d14.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_d14, c("LTLA_name", "specimen_date"))
do_cox("d21.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_d21, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_d28, c("LTLA_name", "specimen_date"))
do_cox("d60.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_d60, c("LTLA_name", "specimen_date"))
do_cox("dNA.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_dNA, c("LTLA_name", "specimen_date"))
do_cox("c28.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_c28, c("LTLA_name", "specimen_date"))
do_cox("e60.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_e60, c("LTLA_name", "specimen_date"))


# 4. NON-OVERLAPPING PERIODS
splitS = survSplit(Surv(time, status) ~ ., data = dataS, cut = c(0, 7, 14, 21, 28), start = "tstart", end = "tstop")
setDT(splitS)
splitS[, sgtf_w1 := ifelse(sgtf == 1 & tstart ==  0, 1, 0)]
splitS[, sgtf_w2 := ifelse(sgtf == 1 & tstart ==  7, 1, 0)]
splitS[, sgtf_w3 := ifelse(sgtf == 1 & tstart == 14, 1, 0)]
splitS[, sgtf_w4 := ifelse(sgtf == 1 & tstart == 21, 1, 0)]
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf_by_week", Surv(tstart, tstop, status) ~ sgtf_w1 + sgtf_w2 + sgtf_w3 + sgtf_w4 + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), splitS, c("LTLA_name", "specimen_date"))


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
    Surv(tstart, tstop, status) ~ sgtf                            + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop               + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop+sgtf:tstop2", 
    Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + sgtf:tstop2 + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))

# LRT time interaction
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0sgtf:tstop.cc",  "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop.cc")
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0sgtf:tstop.ipw", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop.ipw")
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop.cc",   "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop+sgtf:tstop2.cc")
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop.ipw",  "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop+sgtf:tstop2.ipw")

# 5b. TIME-AGE
do_cox("d28.SGTF.ls.2020-11-01..r10.LTLA:date.0age:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + age                          + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ls.2020-11-01..r10.LTLA:date.age:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + age + age:tstop              + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ls.2020-11-01..r10.LTLA:date.age:tstop+age:tstop2", 
    Surv(tstart, tstop, status) ~ sgtf + age + age:tstop + age:tstop2 + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))

# LRT time interaction
lrt("d28.SGTF.ls.2020-11-01..r10.LTLA:date.0age:tstop", "d28.SGTF.ls.2020-11-01..r10.LTLA:date.age:tstop")
lrt("d28.SGTF.ls.2020-11-01..r10.LTLA:date.age:tstop", "d28.SGTF.ls.2020-11-01..r10.LTLA:date.age:tstop+age:tstop2")

# 5c. TIME-SEX
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0sex:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, nk = 3) + sexM                            + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sex:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, nk = 3) + sexM + sexM:tstop               + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sex:tstop+sex:tstop2", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, nk = 3) + sexM + sexM:tstop + sexM:tstop2 + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))

# LRT time interaction
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0sex:tstop", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sex:tstop")
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sex:tstop", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sex:tstop+sex:tstop2")

# 5d. TIME-IMD
do_cox("d28.SGTF.sl.2020-11-01..r10.LTLA:date.0imd:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, nk = 3) + sex + imd                          + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.sl.2020-11-01..r10.LTLA:date.imd:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, nk = 3) + sex + imd + imd:tstop              + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.sl.2020-11-01..r10.LTLA:date.imd:tstop+imd:tstop2", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, nk = 3) + sex + imd + imd:tstop + imd:tstop2 + eth_cat + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))

# LRT time interaction
lrt("d28.SGTF.sl.2020-11-01..r10.LTLA:date.0imd:tstop", "d28.SGTF.sl.2020-11-01..r10.LTLA:date.imd:tstop")
lrt("d28.SGTF.sl.2020-11-01..r10.LTLA:date.imd:tstop", "d28.SGTF.sl.2020-11-01..r10.LTLA:date.imd:tstop+imd:tstop2")

# 5e. TIME-ETHNICITY
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0eth:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + ethA + ethB + ethO + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.eth:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + ethA + ethB + ethO + ethA:tstop + ethB:tstop + ethO:tstop + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.eth:tstop+eth:tstop2", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + ethA + ethB + ethO + ethA:tstop + ethB:tstop + ethO:tstop + ethA:tstop2 + ethB:tstop2 + ethO:tstop2 + res_cat + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))

# LRT time interaction
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0eth:tstop", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.eth:tstop")
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.eth:tstop", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.eth:tstop+eth:tstop2")

# 5f. TIME-RESIDENCE TYPE
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0res:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + resC + resO + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.res:tstop", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + resC + resO + resC:tstop + resO:tstop + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.res:tstop+res:tstop2", 
    Surv(tstart, tstop, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + resC + resO + resC:tstop + resO:tstop + resC:tstop2 + resO:tstop2 + strata(stratum), 
    dataST, c("LTLA_name", "specimen_date"))

# LRT time interaction
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0res:tstop", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.res:tstop")
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.res:tstop", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.res:tstop+res:tstop2")


# 6. SENSITIVITY ANALYSES

# 6a. Cutoff by date range
dataS_sep = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-09-01")
dataS_oct = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-10-01")
dataS_nov = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataS_dec = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-12-01")
dataS_jan = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2021-01-01")
weight_data(dataS_sep, missing_model)
weight_data(dataS_oct, missing_model)
weight_data(dataS_nov, missing_model)
weight_data(dataS_dec, missing_model)
weight_data(dataS_jan, missing_model)
do_cox("d28.SGTF.ss.2020-09-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_sep, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-10-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_oct, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_nov, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-12-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_dec, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2021-01-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_jan, c("LTLA_name", "specimen_date"))

dataS_LTLAprev = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, prevalence_cutoff = TRUE)
weight_data(dataS_LTLAprev, missing_model)
do_cox("d28.SGTF.ss.LTLA..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_LTLAprev, c("LTLA_name", "specimen_date"))

# 6b. No registration cutoff
dataSc = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 0, P_voc = 0, date_min = "2020-11-01")
weight_data(dataSc, missing_model)
do_cox("d28.SGTF.ss.2020-11-01..rNA.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataSc, c("LTLA_name", "specimen_date"))

# 6c. Adjusting for, rather than stratifying by, specimen week and NHSE region
# Seems we need to exclude weeks with fewer data points in order for this to converge; check when updating data.
dataS2 = copy(dataS)
dataS2 = dataS2[specimen_week >= "2020-11-02" & specimen_week < "2021-01-25"]
weight_data(dataS2, missing_model)
do_cox("d28.SGTF.ss.2020-11-02.2020-01-24.r10..NHSE:week", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + NHSER_name:specimen_week_f, dataS2, c())

# 6d. Include individuals with full followup only
dataSd = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 0, P_voc = 0, 
    date_min = "2020-11-01", date_max = as.character(ymd("2021-02-05") - 38))
weight_data(dataSd, missing_model)
do_cox("d28.SGTF.ss.2020-11-01.tminus38.rNA.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataSd, c("LTLA_name", "specimen_date"))

# 6e. Include asymptomatic indicator
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.asymptomatic", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + asymptomatic + strata(stratum), dataS, c("LTLA_name", "specimen_date"))



###########################
# P_VOC SURVIVAL ANALYSES #
###########################

# 1. MAIN ANALYSES

# 1a. WHOLE DATA RANGE
dataP = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0)
weight_data(dataP, missing_model)
do_cox("d28.pVOC.ss.2020-09-01..r10.LTLA:date.", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataP, c("LTLA_name", "specimen_date"))

# 1b. POST-NOV 1 ONLY
dataPn = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
weight_data(dataPn, missing_model)
do_cox("d28.pVOC.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataPn, c("LTLA_name", "specimen_date"))

# 1c. 60 DAY FOLLOWUP P_VOC
dataP60 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 60, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
weight_data(dataP60, missing_model)
do_cox("d60.pVOC.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataP60, c("LTLA_name", "specimen_date"))


# 2. TIME-P_VOC
event_times = dataPn[!is.na(death_date), sort(unique(as.numeric(death_date - specimen_date)))]
event_times = event_times[event_times <= 28]
dataPT = survSplit(Surv(time, status) ~ ., data = dataPn, cut = event_times, start = "tstart", end = "tstop")
setDT(dataPT)
dataPT[, tstop2 := tstop^2]

do_cox("d28.pVOC.ss.2020-11-01..r10.LTLA:date.0p_voc:tstop", 
    Surv(tstart, tstop, status) ~ p_voc                              + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataPT, c("LTLA_name", "specimen_date"))
do_cox("d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc:tstop", 
    Surv(tstart, tstop, status) ~ p_voc + p_voc:tstop                + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataPT, c("LTLA_name", "specimen_date"))
do_cox("d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc:tstop+p_voc:tstop2", 
    Surv(tstart, tstop, status) ~ p_voc + p_voc:tstop + p_voc:tstop2 + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataPT, c("LTLA_name", "specimen_date"))

# LRT time interaction
lrt("d28.pVOC.ss.2020-11-01..r10.LTLA:date.0p_voc:tstop", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc:tstop")
lrt("d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc:tstop", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc:tstop+p_voc:tstop2")



### Additional investigations ###

# HR over time 

dataSt = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataPt = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
weight_data(dataSt, missing_model)
weight_data(dataPt, missing_model)

dataSt[, sgtf_t1 := ifelse(specimen_date >= "2020-11-01" & specimen_date <= "2020-11-30", sgtf, 0)]
dataSt[, sgtf_t2 := ifelse(specimen_date >= "2020-12-01" & specimen_date <= "2020-12-31", sgtf, 0)]
dataSt[, sgtf_t3 := ifelse(specimen_date >= "2021-01-01"                                , sgtf, 0)]
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf_t", Surv(time, status) ~ 
        sgtf_t1 + sgtf_t2 + sgtf_t3 +
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataSt, c("LTLA_name", "specimen_date"))
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date..cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf_t.cc")
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date..ipw", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf_t.ipw")

dataPt[, p_voc_t1 := ifelse(specimen_date >= "2020-11-01" & specimen_date <= "2020-11-30", p_voc, 0)]
dataPt[, p_voc_t2 := ifelse(specimen_date >= "2020-12-01" & specimen_date <= "2020-12-31", p_voc, 0)]
dataPt[, p_voc_t3 := ifelse(specimen_date >= "2021-01-01"                                , p_voc, 0)]
do_cox("d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc_t", Surv(time, status) ~ 
        p_voc_t1 + p_voc_t2 + p_voc_t3 +
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataPt, c("LTLA_name", "specimen_date"))
lrt("d28.pVOC.ss.2020-11-01..r10.LTLA:date..cc", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc_t.cc")
lrt("d28.pVOC.ss.2020-11-01..r10.LTLA:date..ipw", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc_t.ipw")

# HR by region

dataSt[, sgtf_ee := ifelse(NHSER_name == "East of England",          sgtf, 0)]
dataSt[, sgtf_ld := ifelse(NHSER_name == "London",                   sgtf, 0)]
dataSt[, sgtf_ml := ifelse(NHSER_name == "Midlands",                 sgtf, 0)]
dataSt[, sgtf_ne := ifelse(NHSER_name == "North East and Yorkshire", sgtf, 0)]
dataSt[, sgtf_nw := ifelse(NHSER_name == "North West",               sgtf, 0)]
dataSt[, sgtf_se := ifelse(NHSER_name == "South East",               sgtf, 0)]
dataSt[, sgtf_sw := ifelse(NHSER_name == "South West",               sgtf, 0)]

do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf_nhse", Surv(time, status) ~ 
        sgtf_ee + sgtf_ld + sgtf_ml + sgtf_ne + sgtf_nw + sgtf_se + sgtf_sw + 
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataSt, c("LTLA_name", "specimen_date"))
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date..cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf_nhse.cc")
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date..ipw", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf_nhse.ipw")


dataPt[, p_voc_ee := ifelse(NHSER_name == "East of England",          p_voc, 0)]
dataPt[, p_voc_ld := ifelse(NHSER_name == "London",                   p_voc, 0)]
dataPt[, p_voc_ml := ifelse(NHSER_name == "Midlands",                 p_voc, 0)]
dataPt[, p_voc_ne := ifelse(NHSER_name == "North East and Yorkshire", p_voc, 0)]
dataPt[, p_voc_nw := ifelse(NHSER_name == "North West",               p_voc, 0)]
dataPt[, p_voc_se := ifelse(NHSER_name == "South East",               p_voc, 0)]
dataPt[, p_voc_sw := ifelse(NHSER_name == "South West",               p_voc, 0)]

do_cox("d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc_nhse", Surv(time, status) ~ 
        p_voc_ee + p_voc_ld + p_voc_ml + p_voc_ne + p_voc_nw + p_voc_se + p_voc_sw + 
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataPt, c("LTLA_name", "specimen_date"))
lrt("d28.pVOC.ss.2020-11-01..r10.LTLA:date..cc", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc_nhse.cc")
lrt("d28.pVOC.ss.2020-11-01..r10.LTLA:date..ipw", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc_nhse.ipw")


# HR over time by residence type

dataSt[, sgtf_R := ifelse(res_cat == "Residential"      , sgtf, 0)]
dataSt[, sgtf_C := ifelse(res_cat == "Care/Nursing home", sgtf, 0)]
dataSt[, sgtf_O := ifelse(res_cat == "Other/Unknown"    , sgtf, 0)]

dataSt[, sgtf_R_t1 := ifelse(res_cat == "Residential"       & specimen_date >= "2020-11-01" & specimen_date <= "2020-11-30", sgtf, 0)]
dataSt[, sgtf_R_t2 := ifelse(res_cat == "Residential"       & specimen_date >= "2020-12-01" & specimen_date <= "2020-12-31", sgtf, 0)]
dataSt[, sgtf_R_t3 := ifelse(res_cat == "Residential"       & specimen_date >= "2021-01-01"                                , sgtf, 0)]
dataSt[, sgtf_C_t1 := ifelse(res_cat == "Care/Nursing home" & specimen_date >= "2020-11-01" & specimen_date <= "2020-11-30", sgtf, 0)]
dataSt[, sgtf_C_t2 := ifelse(res_cat == "Care/Nursing home" & specimen_date >= "2020-12-01" & specimen_date <= "2020-12-31", sgtf, 0)]
dataSt[, sgtf_C_t3 := ifelse(res_cat == "Care/Nursing home" & specimen_date >= "2021-01-01"                                , sgtf, 0)]
dataSt[, sgtf_O_t1 := ifelse(res_cat == "Other/Unknown"     & specimen_date >= "2020-11-01" & specimen_date <= "2020-11-30", sgtf, 0)]
dataSt[, sgtf_O_t2 := ifelse(res_cat == "Other/Unknown"     & specimen_date >= "2020-12-01" & specimen_date <= "2020-12-31", sgtf, 0)]
dataSt[, sgtf_O_t3 := ifelse(res_cat == "Other/Unknown"     & specimen_date >= "2021-01-01"                                , sgtf, 0)]

do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf_res", Surv(time, status) ~ 
        sgtf_R + 
        sgtf_C + 
        sgtf_O + 
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataSt, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf_resxt", Surv(time, status) ~ 
        sgtf_R_t1 + sgtf_R_t2 + sgtf_R_t3 + 
        sgtf_C_t1 + sgtf_C_t2 + sgtf_C_t3 + 
        sgtf_O_t1 + sgtf_O_t2 + sgtf_O_t3 + 
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataSt, c("LTLA_name", "specimen_date"))
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf_res.cc",  "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf_resxt.cc")
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf_res.ipw", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf_resxt.ipw")

dataPt[, p_voc_R := ifelse(res_cat == "Residential"      , p_voc, 0)]
dataPt[, p_voc_C := ifelse(res_cat == "Care/Nursing home", p_voc, 0)]
dataPt[, p_voc_O := ifelse(res_cat == "Other/Unknown"    , p_voc, 0)]

dataPt[, p_voc_R_t1 := ifelse(res_cat == "Residential"       & specimen_date >= "2020-11-01" & specimen_date <= "2020-11-30", p_voc, 0)]
dataPt[, p_voc_R_t2 := ifelse(res_cat == "Residential"       & specimen_date >= "2020-12-01" & specimen_date <= "2020-12-31", p_voc, 0)]
dataPt[, p_voc_R_t3 := ifelse(res_cat == "Residential"       & specimen_date >= "2021-01-01"                                , p_voc, 0)]
dataPt[, p_voc_C_t1 := ifelse(res_cat == "Care/Nursing home" & specimen_date >= "2020-11-01" & specimen_date <= "2020-11-30", p_voc, 0)]
dataPt[, p_voc_C_t2 := ifelse(res_cat == "Care/Nursing home" & specimen_date >= "2020-12-01" & specimen_date <= "2020-12-31", p_voc, 0)]
dataPt[, p_voc_C_t3 := ifelse(res_cat == "Care/Nursing home" & specimen_date >= "2021-01-01"                                , p_voc, 0)]
dataPt[, p_voc_O_t1 := ifelse(res_cat == "Other/Unknown"     & specimen_date >= "2020-11-01" & specimen_date <= "2020-11-30", p_voc, 0)]
dataPt[, p_voc_O_t2 := ifelse(res_cat == "Other/Unknown"     & specimen_date >= "2020-12-01" & specimen_date <= "2020-12-31", p_voc, 0)]
dataPt[, p_voc_O_t3 := ifelse(res_cat == "Other/Unknown"     & specimen_date >= "2021-01-01"                                , p_voc, 0)]

do_cox("d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc_res", Surv(time, status) ~ 
        p_voc_R + 
        p_voc_C + 
        p_voc_O + 
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataPt, c("LTLA_name", "specimen_date"))
do_cox("d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc_resxt", Surv(time, status) ~ 
        p_voc_R_t1 + p_voc_R_t2 + p_voc_R_t3 + 
        p_voc_C_t1 + p_voc_C_t2 + p_voc_C_t3 + 
        p_voc_O_t1 + p_voc_O_t2 + p_voc_O_t3 + 
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataPt, c("LTLA_name", "specimen_date"))
lrt("d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc_res.cc",  "d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc_resxt.cc")
lrt("d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc_res.ipw", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc_resxt.ipw")





#########
# PLOTS #
#########

# Hazard ratio over time
plot_hazard = function(marker, constant_model, linear_varying_model, ylimits, y_label)
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
        geom_ribbon(aes(x = t, ymin = exp(c.lo), ymax = exp(c.hi), fill = "Constant hazard ratio (proportional hazards)"), alpha = 0.4) +
        geom_ribbon(aes(x = t, ymin = exp(x.lo), ymax = exp(x.hi), fill = "Time-varying hazard ratio"), alpha = 0.4) +
        geom_line(aes(x = t, y = exp(c.ct), colour = "Constant hazard ratio (proportional hazards)")) +
        geom_line(aes(x = t, y = exp(x.ct), colour = "Time-varying hazard ratio")) +
        geom_hline(aes(yintercept = 1), linetype = "33", size = 0.25) +
        labs(x = "Days post specimen", y = y_label, fill = "Model", colour = "Model") +
        scale_x_continuous(breaks = c(0, 7, 14, 21, 28)) +
        scale_colour_manual(aesthetics = c("colour", "fill"), values = c("darkorchid", "#44aa88")) +
        theme(legend.position = c(0.1, 0.8)) +
        ylim(ylimits)
    
    # Uncomment if want to extract time-varying hazards
    # return (haz)
}

# haz = plot_hazard("sgtf", "SGTF + lin age + lin IMD | LTLA + spec date", "Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | LTLA + spec date", c(0, NA))
# haz[t == 1, .(exp(x.ct), exp(x.lo), exp(x.hi))]

theme_set(theme_cowplot(font_size = 10))

hp_sgtf_cc = plot_hazard("sgtf", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.0sgtf:tstop.cc", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop.cc", 
    c(0, 4), "Hazard ratio for SGTF")
hp_voc_cc = plot_hazard("p_voc", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.0p_voc:tstop.cc", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc:tstop.cc", 
    c(0, 4), expression("Hazard ratio for"~p[VOC])) + theme(legend.position = "none")
hp_sgtf_ipw = plot_hazard("sgtf", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.0sgtf:tstop.ipw", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop.ipw", 
    c(0, 4), "Hazard ratio for SGTF") + theme(legend.position = "none")
hp_voc_ipw = plot_hazard("p_voc", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.0p_voc:tstop.ipw", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc:tstop.ipw", 
    c(0, 4), expression("Hazard ratio for"~p[VOC])) + theme(legend.position = "none")


# Effects of SGTF

summaries = fread("./output/model_summary.csv")
summaries = summaries[str_count(model_id, "\\.") == 8] # TODO note this gets rid of some interaction tests
summaries[, model_id := str_replace_all(model_id, "\\.ll\\.", ".L.L.")]
summaries[, model_id := str_replace_all(model_id, "\\.sl\\.", ".S.L.")]
summaries[, model_id := str_replace_all(model_id, "\\.ls\\.", ".L.S.")]
summaries[, model_id := str_replace_all(model_id, "\\.ss\\.", ".S.S.")]
summaries[, model_id := str_replace_all(model_id, "\\.r10\\.", ".10.")]
summaries[, model_id := str_replace_all(model_id, "\\.rNA\\.", ".0.")]
summaries[, c("m_deaths", "m_marker", "m_age", "m_imd", "m_startdate", "m_enddate", "m_reg_cutoff", "m_strata", "m_xvars", "m_weighting") := tstrsplit(model_id, "\\.")]
summaries[m_enddate == "tminus38", m_enddate := "T - 38"]
summaries = summaries[order(model_id)]

groups = list(
    `Geographical and temporal stratification` = list(
        pattern = "d28/SGTF/S/S/2020-11-01//10/*//*",
        highlight = 8
    ),
    `Age and IMD terms (linear vs. spline)` = list(
        pattern = "d28/SGTF/*/*/2020-11-01//10/LTLA:date//*",
        highlight = c(3, 4)
    ),
    # Interactions
    `Death type` = list(
        pattern = "*/SGTF/S/S/2020-11-01//10/LTLA:date//*",
        highlight = 1
    ),
    # SGTF by week / nonoverlapping periods into episode
    # SGTF by calendar period
    # Time-X interaction
    `Analysis start date` = list(
        pattern = "d28/SGTF/S/S/*//10/LTLA:date//*",
        highlight = 5
    ),
    
    `Misclassification adjustment` = list(
        pattern = "*/pVOC/*/*/*/*/*/*//*",
        highlight = c(1, 2, 5)
    ),
    
    `Miscellaneous` = list(
        pattern = "d28/SGTF/S/S/2020-11-01//0/LTLA:date//*",
        highlight = c(7)
    ),
    
    `.1` = list(
        pattern = "d28/SGTF/S/S/2020-11-02/2020-01-24/10//NHSE:week/*",
        highlight = c(9)
    ),
    
    `.2` = list(
        pattern = "d28/SGTF/S/S/2020-11-01/tminus38/0/LTLA:date//*",
        highlight = c(6, 7)
    ),
    
    `.3` = list(
        pattern = "d28/SGTF/S/S/2020-11-01//10/LTLA:date/asymptomatic/*",
        highlight = c(9)
    )
)

library(ggstance)
theme_set(theme_cowplot(font_size = 10))

pl_effects = plot_forest(summaries, groups, 13, -2.8, 0.3, x_widths = c(0.7, 1.3, 1.1, 0.7, 1.6, 2.6, 1.8, 1.8, 2.8))
ggsave("./output/summary.png", pl_effects, width = 20, height = 20, units = "cm")


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
    y_index = 1;
    for (g in seq_along(groups)) {
        if (!names(groups)[g] %like% "^\\.") {
            headings = rbind(headings, 
                data.table(label = names(groups)[g], y_index)
            )
            y_index = y_index + 1;
        }
        model_ids = summaries[model_id %like% interpret(groups[[g]]$pattern), unique(model_id)];
        model_ids_genus = unique(str_remove_all(model_ids, "\\.cc$|\\.ipw$"))
        
        for (model in model_ids_genus) {
            pd = rbind(pd,
                summaries[(model_id == model | model_id == paste0(model, ".cc") | model_id == paste0(model, ".ipw")) & parameter %like% "^sgtf$|^p_voc$", 
                    .(h = HR, h0 = HR.lo95, h1 = HR.hi95, 
                      blank = "", highlight = paste(table_cols[groups[[g]]$highlight], collapse = "|"), y_index,
                      m_deaths, m_marker, m_age, m_imd, m_startdate, m_enddate, m_reg_cutoff, m_strata, m_xvars, m_weighting)]
            )
            y_index = y_index + 1;
        }
    }
    
    # Make labels
    blank_col = which(names(pd) == "blank");
    labels = melt(pd[m_weighting == "cc", .SD, .SDcols = (blank_col + 1):(ncol(pd) - 1)], id.vars = c("y_index", "highlight"))
    labels[value == "", value := "—"]
    labels[, x_index := as.numeric(variable)]
    x_shifts = cumsum(x_widths) - x_widths[1]
    labels[, x_pos := tab_from + (tab_to - tab_from) * x_shifts[x_index] / max(x_shifts)]
    labels[, colour := ifelse(variable %like% highlight, "#000000", "#888888"), by = .(y_index, x_index)]
    
    # Make table structure
    tabr = data.table(xmin = tab_from - 0.2, xmax = 2.6, 
        ymin = labels[, unique(y_index) - 0.5], ymax = labels[, unique(y_index) + 0.5])
    tabr[, fill := rep_len(c("#f0f0f0", "#f8f8f8"), .N)]
    tabt = data.table(label = c("Death\ntype", "Marker", "Age", "IMD", "Start\ndate", "End\ndate", 
        "Registr.\ncutoff", "Strati-\nfication", "Extra\nterms"), x = labels[, unique(x_pos)], y = y_index)
    
    pd[, weighting := ifelse(m_weighting == "cc", "Complete-case", "IPW")]
    pd[, weighting := factor(weighting, c("IPW", "Complete-case"))]

    ggplot(pd) +
        annotate("rect", xmin = tabr$xmin, xmax = tabr$xmax, ymin = tabr$ymin, ymax = tabr$ymax, fill = tabr$fill) +
        geom_point(aes(x = h, y = y_index, colour = weighting), size = 0.5, position = ggstance::position_dodgev(height = 0.8)) +
        geom_linerange(aes(x = h, xmin = h0, xmax = h1, y = y_index, colour = weighting), size = 0.5, position = position_dodge(width = 0.8)) +
        scale_y_reverse(breaks = pd$y_index, labels = pd$blank) +
        geom_vline(aes(xintercept = 1.0), size = 0.125) +
        geom_vline(aes(xintercept = 1.5), size = 0.125, colour = "#aaaaaa") +
        geom_vline(aes(xintercept = 2.0), size = 0.125, colour = "#aaaaaa") +
        geom_vline(aes(xintercept = 2.5), size = 0.125, colour = "#aaaaaa") +
        scale_colour_manual(values = c("Complete-case" = "#bbaa88", "IPW" = "#9999bb"), 
            guide = guide_legend(reverse = TRUE, override.aes = list(fatten = 0.1))) +
        labs(x = "Hazard ratio", y = NULL, colour = NULL) +
        theme(plot.margin = unit(c(0.1, 0.3, 0.1, left_margin), "cm"),
            legend.position = c(0.1, 0.975)) +
        annotate("text", x = tabt$x, y = tabt$y, label = tabt$label, hjust = 1, vjust = 0.5, angle = 90, lineheight = 0.75, size = 3.5) +
        annotate("text", x = labels$x_pos, y = labels$y_index, label = labels$value, colour = labels$colour, size = 3.5) +
        annotate("text", x = mean(c(tab_from, tab_to)), y = headings$y_index, label = headings$label, size = 3.5) +
        coord_cartesian(xlim = c(0.8, 2.5), clip = 'off')
}




# Exploratory data plots
dataS[, SGTF_label := factor(ifelse(sgtf == 1, "SGTF", "Non-SGTF"), levels = c("SGTF", "Non-SGTF"))]
plot_samples = ggplot(dataS[, .N, keyby = .(specimen_date, SGTF_label)]) +
    geom_col(aes(specimen_date, N, fill = SGTF_label), position = "stack") +
    labs(x = "Specimen date", y = "Samples", fill = NULL) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = ymd(c("2020-10-31", "2021-02-01"))) +
    theme(legend.position = c(0.1, 0.9))

plot_deaths = ggplot(dataS[died == TRUE, .N, keyby = .(specimen_date, SGTF_label)]) +
    geom_col(aes(specimen_date, N, fill = SGTF_label), position = "stack") +
    labs(x = "Specimen date", y = "Deaths", fill = NULL) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = ymd(c("2020-10-31", "2021-02-01"))) +
    theme(legend.position = "none")

plot_censoring = ggplot(dataS[!is.na(death_date), .(specimen_date, sgtf, died, raw_death_delay = as.numeric(death_date - specimen_date))]) +
    geom_point(aes(specimen_date, raw_death_delay, shape = ifelse(sgtf == 1, "Failure", "Present"), colour = ifelse(died, "Died", "Censored")), alpha = 0.5) +
    labs(x = "Specimen date", y = "Time to death", shape = "S gene", colour = "Status") +
    scale_colour_manual(values = c("Died" = "#000000", "Censored" = "#aa88cc"))

tbl_out = summaries[parameter %in% c("sgtf", "p_voc") & !model_id %like% "Time.{0,2}-(SGTF|p_voc)", 
    .(`Model` = model_id, `Parameter` = parameter, `Hazard ratio` = paste0(round(HR, 2), " (", round(HR.lo95, 2), "–", round(HR.hi95, 2), ")"),
        `P value` = ifelse(P < 0.001, "< 0.001", round(P, 3)))]

fwrite(tbl_out, "./output/table_effects.csv")

# Test proportional hazards assumption
cd = complete_data("20210122")
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





