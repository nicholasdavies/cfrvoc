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

load_model = function(mdl_id)
{
    if (class(mdl_id) == "coxph") {
        return (mdl_id)
    } else {
        return (qread(paste0("./output/", make.names(mdl_id), ".qs")))
    }
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

# Do likelihood ratio test for two saved models
lrt = function(id_base, id_nested)
{
    model_base = load_model(id_base);
    model_nested = load_model(id_nested);
    anova(model_nested, model_base);
}



# Just residential etc
cd1 = complete_data("20210122")
dataS1 = model_data(cd1, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")[res_cat == "Care/Nursing home"]
do_cox("test 000 SGTF + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + strata(stratum), dataS1, c("LTLA_name", "specimen_date"))

cd2 = complete_data("20210129")
dataS2 = model_data(cd2, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")[res_cat == "Care/Nursing home"]
do_cox("test 000 SGTF + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + strata(stratum), dataS2, c("LTLA_name", "specimen_date"))



cd = complete_data("20210122")


# Individual months
dataS_sep_only = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-09-01", date_max = "2020-09-30")
do_cox("Sensitivity: Sep only (SGTF + lin age + lin IMD | LTLA + spec date)", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS_sep_only, c("LTLA_name", "specimen_date"))
dataS_oct_only = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-10-01", date_max = "2020-10-31")
do_cox("Sensitivity: Oct only (SGTF + lin age + lin IMD | LTLA + spec date)", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS_oct_only, c("LTLA_name", "specimen_date"))
dataS_nov_only = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01", date_max = "2020-11-30")
do_cox("Sensitivity: Nov only (SGTF + lin age + lin IMD | LTLA + spec date)", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS_nov_only, c("LTLA_name", "specimen_date"))
dataS_dec_only = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-12-01", date_max = "2020-12-31")
do_cox("Sensitivity: Dec only (SGTF + lin age + lin IMD | LTLA + spec date)", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS_dec_only, c("LTLA_name", "specimen_date"))
dataS_jan_only = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2021-01-01", date_max = "2021-01-31")
do_cox("Sensitivity: Jan only (SGTF + lin age + lin IMD | LTLA + spec date)", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS_jan_only, c("LTLA_name", "specimen_date"))

# Test: SGTF:pressure interaction
cd = complete_data("20210129", pressure_shift = 9)
dataS = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")

mbase = do_cox("test base", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))
mp_mabs = do_cox("test mabs", Surv(time, status) ~ sgtf + sgtf:medstaff_abs_per_bed + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))
mp_nabs = do_cox("test nabs", Surv(time, status) ~ sgtf + sgtf:nursing_abs_per_bed + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))
mp_mvpr = do_cox("test mvpr", Surv(time, status) ~ sgtf + sgtf:mv_pressure + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))
mp_nipr = do_cox("test nipr", Surv(time, status) ~ sgtf + sgtf:ni_pressure + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))
mp_ospr = do_cox("test ospr", Surv(time, status) ~ sgtf + sgtf:os_pressure + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))
mp_aopr = do_cox("test aopr", Surv(time, status) ~ sgtf + sgtf:ao_pressure + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))

lrt(mbase, mp_mabs) # 0.55
lrt(mbase, mp_nabs) # 0.60
lrt(mbase, mp_mvpr) # 0.35
lrt(mbase, mp_nipr) # 0.31
lrt(mbase, mp_ospr) # 0.39
lrt(mbase, mp_aopr) # 0.76


m2base = do_cox("test base", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))
m2p_mabs = do_cox("test mabs", Surv(time, status) ~ sgtf + medstaff_abs_per_bed + sgtf:medstaff_abs_per_bed + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))
m2p_nabs = do_cox("test nabs", Surv(time, status) ~ sgtf + nursing_abs_per_bed + sgtf:nursing_abs_per_bed + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))
m2p_mvpr = do_cox("test mvpr", Surv(time, status) ~ sgtf + mv_pressure + sgtf:mv_pressure + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))

lrt(m2base, m2p_mabs)
lrt(m2base, m2p_nabs)
lrt(m2base, m2p_mvpr)



##########################
# SGTF SURVIVAL ANALYSES #
##########################

# Load complete data set
cd = complete_data("20210129")
##cd = complete_data("20210201")

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
do_cox("d28.SGTF.ll.2020-11-01..r10.NHSE:week.", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_week"))
do_cox("d28.SGTF.ll.2020-11-01..r10.UTLA:week.", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("UTLA_name", "specimen_week"))
do_cox("d28.SGTF.ll.2020-11-01..r10.LTLA:week.", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_week"))
do_cox("d28.SGTF.ll.2020-11-01..r10.NHSE:date.", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_date"))
do_cox("d28.SGTF.ll.2020-11-01..r10.UTLA:date.", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("UTLA_name", "specimen_date"))
do_cox("d28.SGTF.ll.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))

# Spline age, linear IMD
do_cox("d28.SGTF.sl.2020-11-01..r10.NHSE:week.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_week"))
do_cox("d28.SGTF.sl.2020-11-01..r10.UTLA:week.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("UTLA_name", "specimen_week"))
do_cox("d28.SGTF.sl.2020-11-01..r10.LTLA:week.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_week"))
do_cox("d28.SGTF.sl.2020-11-01..r10.NHSE:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_date"))
do_cox("d28.SGTF.sl.2020-11-01..r10.UTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("UTLA_name", "specimen_date"))
do_cox("d28.SGTF.sl.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + imd + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))

# Linear age, spline IMD
do_cox("d28.SGTF.ls.2020-11-01..r10.NHSE:week.", Surv(time, status) ~ sgtf + age + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_week"))
do_cox("d28.SGTF.ls.2020-11-01..r10.UTLA:week.", Surv(time, status) ~ sgtf + age + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("UTLA_name", "specimen_week"))
do_cox("d28.SGTF.ls.2020-11-01..r10.LTLA:week.", Surv(time, status) ~ sgtf + age + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_week"))
do_cox("d28.SGTF.ls.2020-11-01..r10.NHSE:date.", Surv(time, status) ~ sgtf + age + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_date"))
do_cox("d28.SGTF.ls.2020-11-01..r10.UTLA:date.", Surv(time, status) ~ sgtf + age + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("UTLA_name", "specimen_date"))
do_cox("d28.SGTF.ls.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + age + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))

# Spline age, spline IMD
do_cox("d28.SGTF.ss.2020-11-01..r10.NHSE:week.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_week"))
do_cox("d28.SGTF.ss.2020-11-01..r10.UTLA:week.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("UTLA_name", "specimen_week"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:week.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_week"))
do_cox("d28.SGTF.ss.2020-11-01..r10.NHSE:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("NHSER_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.UTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("UTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS, c("LTLA_name", "specimen_date"))

# LRT linear age, linear IMD -> spline age, linear IMD
lrt("d28.SGTF.ll.2020-11-01..r10.NHSE:week.", "d28.SGTF.sl.2020-11-01..r10.NHSE:week.")
lrt("d28.SGTF.ll.2020-11-01..r10.UTLA:week.", "d28.SGTF.sl.2020-11-01..r10.UTLA:week.")
lrt("d28.SGTF.ll.2020-11-01..r10.LTLA:week.", "d28.SGTF.sl.2020-11-01..r10.LTLA:week.")
lrt("d28.SGTF.ll.2020-11-01..r10.NHSE:date.", "d28.SGTF.sl.2020-11-01..r10.NHSE:date.")
lrt("d28.SGTF.ll.2020-11-01..r10.UTLA:date.", "d28.SGTF.sl.2020-11-01..r10.UTLA:date.")
lrt("d28.SGTF.ll.2020-11-01..r10.LTLA:date.", "d28.SGTF.sl.2020-11-01..r10.LTLA:date.")

# LRT linear age, linear IMD -> linear age, spline IMD
lrt("d28.SGTF.ll.2020-11-01..r10.NHSE:week.", "d28.SGTF.ls.2020-11-01..r10.NHSE:week.")
lrt("d28.SGTF.ll.2020-11-01..r10.UTLA:week.", "d28.SGTF.ls.2020-11-01..r10.UTLA:week.")
lrt("d28.SGTF.ll.2020-11-01..r10.LTLA:week.", "d28.SGTF.ls.2020-11-01..r10.LTLA:week.")
lrt("d28.SGTF.ll.2020-11-01..r10.NHSE:date.", "d28.SGTF.ls.2020-11-01..r10.NHSE:date.")
lrt("d28.SGTF.ll.2020-11-01..r10.UTLA:date.", "d28.SGTF.ls.2020-11-01..r10.UTLA:date.")
lrt("d28.SGTF.ll.2020-11-01..r10.LTLA:date.", "d28.SGTF.ls.2020-11-01..r10.LTLA:date.")

# LRT linear age, spline IMD -> spline age, spline IMD
lrt("d28.SGTF.ls.2020-11-01..r10.NHSE:week.", "d28.SGTF.ss.2020-11-01..r10.NHSE:week.")
lrt("d28.SGTF.ls.2020-11-01..r10.UTLA:week.", "d28.SGTF.ss.2020-11-01..r10.UTLA:week.")
lrt("d28.SGTF.ls.2020-11-01..r10.LTLA:week.", "d28.SGTF.ss.2020-11-01..r10.LTLA:week.")
lrt("d28.SGTF.ls.2020-11-01..r10.NHSE:date.", "d28.SGTF.ss.2020-11-01..r10.NHSE:date.")
lrt("d28.SGTF.ls.2020-11-01..r10.UTLA:date.", "d28.SGTF.ss.2020-11-01..r10.UTLA:date.")
lrt("d28.SGTF.ls.2020-11-01..r10.LTLA:date.", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.")

# LRT spline age, linear IMD -> spline age, spline IMD
lrt("d28.SGTF.sl.2020-11-01..r10.NHSE:week.", "d28.SGTF.ss.2020-11-01..r10.NHSE:week.")
lrt("d28.SGTF.sl.2020-11-01..r10.UTLA:week.", "d28.SGTF.ss.2020-11-01..r10.UTLA:week.")
lrt("d28.SGTF.sl.2020-11-01..r10.LTLA:week.", "d28.SGTF.ss.2020-11-01..r10.LTLA:week.")
lrt("d28.SGTF.sl.2020-11-01..r10.NHSE:date.", "d28.SGTF.ss.2020-11-01..r10.NHSE:date.")
lrt("d28.SGTF.sl.2020-11-01..r10.UTLA:date.", "d28.SGTF.ss.2020-11-01..r10.UTLA:date.")
lrt("d28.SGTF.sl.2020-11-01..r10.LTLA:date.", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.")


# 2. INTERACTION TESTS
do_cox_interaction("age", "sgtf", "age_group", c("age_group",        "sex", "rcs(imd, nk = 3)", "eth_cat", "res_cat"), dataS, c("LTLA_name", "specimen_date"))
do_cox_interaction("sex", "sgtf", "sex",       c("rcs(age, nk = 3)", "sex", "rcs(imd, nk = 3)", "eth_cat", "res_cat"), dataS, c("LTLA_name", "specimen_date"))
do_cox_interaction("imd", "sgtf", "imd_group", c("rcs(age, nk = 3)", "sex", "imd_group",        "eth_cat", "res_cat"), dataS, c("LTLA_name", "specimen_date"))
do_cox_interaction("eth", "sgtf", "eth_cat",   c("rcs(age, nk = 3)", "sex", "rcs(imd, nk = 3)", "eth_cat", "res_cat"), dataS, c("LTLA_name", "specimen_date"))
do_cox_interaction("res", "sgtf", "res_cat",   c("rcs(age, nk = 3)", "sex", "rcs(imd, nk = 3)", "eth_cat", "res_cat"), dataS, c("LTLA_name", "specimen_date"))

# # Investigating marginally significant SGTF x residence type term
# dataR = copy(dataS)
# dataR[, table(res_cat, useNA = "ifany")]
# dataR[, sgtf_resR := ifelse(res_cat == "Residential", sgtf, 0)]
# dataR[, sgtf_resC := ifelse(res_cat == "Care/Nursing home", sgtf, 0)]
# dataR[, sgtf_resO := ifelse(res_cat == "Other/Unknown", sgtf, 0)]
# do_cox("SGTF:res_cat + lin age + lin IMD | LTLA + spec date", Surv(time, status) ~ age + sex + imd + eth_cat + sgtf_resR + sgtf_resC + sgtf_resO + res_cat + strata(stratum), dataR, c("LTLA_name", "specimen_date"))

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
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.0sgtf:tstop", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop")
lrt("d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop+sgtf:tstop2")

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
do_cox("d28.SGTF.ss.2020-09-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_sep, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-10-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_oct, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_nov, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2020-12-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_dec, c("LTLA_name", "specimen_date"))
do_cox("d28.SGTF.ss.2021-01-01..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_jan, c("LTLA_name", "specimen_date"))

dataS_LTLAprev = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, prevalence_cutoff = TRUE)
do_cox("d28.SGTF.ss.LTLA..r10.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataS_LTLAprev, c("LTLA_name", "specimen_date"))

# 6b. No registration cutoff
dataSc = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 0, P_voc = 0, date_min = "2020-11-01")
do_cox("d28.SGTF.ss.2020-11-01..rNA.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataSc, c("LTLA_name", "specimen_date"))

# 6c. Adjusting for, rather than stratifying by, specimen week and NHSE region
# Seems we need to exclude weeks with fewer data points in order for this to converge; check when updating data.
dataS2 = copy(dataS)
dataS2 = dataS2[specimen_week >= "2020-11-02" & specimen_week < "2021-01-18"]
dataS2[, specimen_week_f := factor(specimen_week)]
do_cox("d28.SGTF.ss.2020-11-02.2020-01-17.r10..NHSE:week", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + NHSER_name:specimen_week_f, dataS2, c())

# 6d. Include individuals with full followup only
dataSd = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 0, P_voc = 0, 
    date_min = "2020-11-01", date_max = as.character(ymd("2021-01-29") - 38))
do_cox("d28.SGTF.ss.2020-11-01.tminus38.rNA.LTLA:date.", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataSd, c("LTLA_name", "specimen_date"))

# 6e. Include asymptomatic indicator
do_cox("d28.SGTF.ss.2020-11-01..r10.LTLA:date.asymptomatic", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + asymptomatic + strata(stratum), dataS, c("LTLA_name", "specimen_date"))


### INVESTIGATIONS OF INCREASED HAZARD

# Suggestions for further investigation are:
# 1. Fit the model below (i.e. without time-varying HR) separately by some specimen date categories  e.g. 1 or 2 week periods.
dataSX = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
dataPX = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0)

do_cox("x Sens 1-15 Nov: SGTF + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataSX[specimen_date >= "2020-11-01" & specimen_date <= "2020-11-15"] , c("LTLA_name", "specimen_date"))
do_cox("x Sens 16-30 Nov: SGTF + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataSX[specimen_date >= "2020-11-16" & specimen_date <= "2020-11-30"] , c("LTLA_name", "specimen_date"))
do_cox("x Sens 1-15 Dec: SGTF + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataSX[specimen_date >= "2020-12-01" & specimen_date <= "2020-12-15"] , c("LTLA_name", "specimen_date"))
do_cox("x Sens 16-31 Dec: SGTF + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataSX[specimen_date >= "2020-12-16" & specimen_date <= "2020-12-31"] , c("LTLA_name", "specimen_date"))
do_cox("x Sens 1-15 Jan: SGTF + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ sgtf + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataSX[specimen_date >= "2021-01-01" & specimen_date <= "2021-01-15"] , c("LTLA_name", "specimen_date"))

do_cox("x Sens 1-15 Nov: p_voc + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataPX[specimen_date >= "2020-11-01" & specimen_date <= "2020-11-15"] , c("LTLA_name", "specimen_date"))
do_cox("x Sens 16-30 Nov: p_voc + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataPX[specimen_date >= "2020-11-16" & specimen_date <= "2020-11-30"] , c("LTLA_name", "specimen_date"))
do_cox("x Sens 1-15 Dec: p_voc + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataPX[specimen_date >= "2020-12-01" & specimen_date <= "2020-12-15"] , c("LTLA_name", "specimen_date"))
do_cox("x Sens 16-31 Dec: p_voc + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataPX[specimen_date >= "2020-12-16" & specimen_date <= "2020-12-31"] , c("LTLA_name", "specimen_date"))
do_cox("x Sens 1-15 Jan: p_voc + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataPX[specimen_date >= "2021-01-01" & specimen_date <= "2021-01-15"] , c("LTLA_name", "specimen_date"))

# 2. Fit the model below but adding an interaction between SGTF and specimen date (categorized) and perform a test of the 
# interaction using a likelihood ratio test.
dataSX[, sgtf_t1 := ifelse(specimen_date >= "2020-11-01" & specimen_date <= "2020-11-15", sgtf, 0)]
dataSX[, sgtf_t2 := ifelse(specimen_date >= "2020-11-16" & specimen_date <= "2020-11-30", sgtf, 0)]
dataSX[, sgtf_t3 := ifelse(specimen_date >= "2020-12-01" & specimen_date <= "2020-12-15", sgtf, 0)]
dataSX[, sgtf_t4 := ifelse(specimen_date >= "2020-12-16" & specimen_date <= "2020-12-31", sgtf, 0)]
dataSX[, sgtf_t5 := ifelse(specimen_date >= "2021-01-01" & specimen_date <= "2021-01-18", sgtf, 0)]
do_cox("x Sens SGTF-specimen date interaction: SGTF + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ sgtf_t1 + sgtf_t2 + sgtf_t3 + sgtf_t4 + sgtf_t5 + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataSX, c("LTLA_name", "specimen_date"))
lrt("SGTF + spl age + spl IMD | LTLA + spec date", "x Sens SGTF-specimen date interaction: SGTF + spl age + spl IMD | LTLA + spec date")

do_cox("x Sens SGTF-specimen date-IMD interaction | LTLA + spec date", Surv(time, status) ~ sgtf_t1*imd + sgtf_t2*imd + sgtf_t3*imd + sgtf_t4*imd + sgtf_t5*imd + rcs(age, nk = 3) + sex + eth_cat + res_cat + strata(stratum), dataSX, c("LTLA_name", "specimen_date"))
lrt("SGTF + spl age + spl IMD | LTLA + spec date", "x Sens SGTF-specimen date-IMD interaction | LTLA + spec date")

dataPX[, p_voc_t0 := ifelse(specimen_date <  "2020-11-01",                                 p_voc, 0)]
dataPX[, p_voc_t1 := ifelse(specimen_date >= "2020-11-01" & specimen_date <= "2020-11-15", p_voc, 0)]
dataPX[, p_voc_t2 := ifelse(specimen_date >= "2020-11-16" & specimen_date <= "2020-11-30", p_voc, 0)]
dataPX[, p_voc_t3 := ifelse(specimen_date >= "2020-12-01" & specimen_date <= "2020-12-15", p_voc, 0)]
dataPX[, p_voc_t4 := ifelse(specimen_date >= "2020-12-16" & specimen_date <= "2020-12-31", p_voc, 0)]
dataPX[, p_voc_t5 := ifelse(specimen_date >= "2021-01-01" & specimen_date <= "2021-01-18", p_voc, 0)]
do_cox("x Sens p_voc-specimen date interaction: p_voc + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ p_voc_t0 + p_voc_t1 + p_voc_t2 + p_voc_t3 + p_voc_t4 + p_voc_t5 + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataPX, c("LTLA_name", "specimen_date"))
lrt("p_voc + spl age + spl IMD | LTLA + spec date", "x Sens p_voc-specimen date interaction: p_voc + spl age + spl IMD | LTLA + spec date")

# by epi week . . .
dataSX[, specimen_epiweek := epiweek(specimen_date)]

dataSX[, sgtf_w46 := ifelse(specimen_epiweek == 46, sgtf, 0)]
dataSX[, sgtf_w47 := ifelse(specimen_epiweek == 47, sgtf, 0)]
dataSX[, sgtf_w48 := ifelse(specimen_epiweek == 48, sgtf, 0)]
dataSX[, sgtf_w49 := ifelse(specimen_epiweek == 49, sgtf, 0)]
dataSX[, sgtf_w50 := ifelse(specimen_epiweek == 50, sgtf, 0)]
dataSX[, sgtf_w51 := ifelse(specimen_epiweek == 51, sgtf, 0)]
dataSX[, sgtf_w52 := ifelse(specimen_epiweek == 52, sgtf, 0)]
dataSX[, sgtf_w53 := ifelse(specimen_epiweek == 53, sgtf, 0)]

do_cox("xx Sens SGTF-specimen epiweek interaction: SGTF + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ 
        sgtf +
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataSX[specimen_epiweek %between% c(46, 53)], c("LTLA_name", "specimen_date"))

do_cox("x Sens SGTF-specimen epiweek interaction: SGTF + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ 
        sgtf_w46 + sgtf_w47 + sgtf_w48 + sgtf_w49 + sgtf_w50 + sgtf_w51 + sgtf_w52 + sgtf_w53 +
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataSX[specimen_epiweek %between% c(46, 53)], c("LTLA_name", "specimen_date"))

lrt("xx Sens SGTF-specimen epiweek interaction: SGTF + spl age + spl IMD | LTLA + spec date", "x Sens SGTF-specimen epiweek interaction: SGTF + spl age + spl IMD | LTLA + spec date")


# 3. Fit the model with the time-varying HR separately by some specimen date categories and plot the results, i.e. plot HR 
# overtime for each specimen date category. Some of the more recently specimen date categories will have shorter followup.
event_times = dataSX[!is.na(death_date), sort(unique(as.numeric(death_date - specimen_date)))]
event_times = event_times[event_times <= 28]
dataSTX = survSplit(Surv(time, status) ~ ., data = dataSX, cut = event_times, start = "tstart", end = "tstop")
setDT(dataSTX)

do_cox("d28.S.ss.2020-11-01.2020-11-15.r10.LTLA:date.sgtf:tstop", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataSTX[specimen_date >= "2020-11-01" & specimen_date <= "2020-11-15"] , c("LTLA_name", "specimen_date"))
do_cox("d28.S.ss.2020-11-16.2020-11-30.r10.LTLA:date.sgtf:tstop", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataSTX[specimen_date >= "2020-11-16" & specimen_date <= "2020-11-30"] , c("LTLA_name", "specimen_date"))
do_cox("d28.S.ss.2020-12-01.2020-12-15.r10.LTLA:date.sgtf:tstop", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataSTX[specimen_date >= "2020-12-01" & specimen_date <= "2020-12-15"] , c("LTLA_name", "specimen_date"))
do_cox("d28.S.ss.2020-12-16.2020-12-31.r10.LTLA:date.sgtf:tstop", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataSTX[specimen_date >= "2020-12-16" & specimen_date <= "2020-12-31"] , c("LTLA_name", "specimen_date"))
do_cox("d28.S.ss.2021-01-01.2021-01-15.r10.LTLA:date.sgtf:tstop", Surv(tstart, tstop, status) ~ sgtf + sgtf:tstop + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), 
    dataSTX[specimen_date >= "2021-01-01" & specimen_date <= "2021-01-15"] , c("LTLA_name", "specimen_date"))


pl1 = plot_hazard("sgtf", 
    "x Sens 1-15 Nov: SGTF + spl age + spl IMD | LTLA + spec date", 
    "d28.S.ss.2020-11-01.2020-11-15.r10.LTLA:date.sgtf:tstop", c(0, 4), "Hazard ratio for SGTF, 1-15 Nov")

pl2 = plot_hazard("sgtf", 
    "x Sens 16-30 Nov: SGTF + spl age + spl IMD | LTLA + spec date", 
    "d28.S.ss.2020-11-16.2020-11-30.r10.LTLA:date.sgtf:tstop", c(0, 4), "Hazard ratio for SGTF, 16-30 Nov")

pl3 = plot_hazard("sgtf", 
    "x Sens 1-15 Dec: SGTF + spl age + spl IMD | LTLA + spec date", 
    "d28.S.ss.2020-12-01.2020-12-15.r10.LTLA:date.sgtf:tstop", c(0, 4), "Hazard ratio for SGTF, 1-15 Dec")

pl4 = plot_hazard("sgtf", 
    "x Sens 16-31 Dec: SGTF + spl age + spl IMD | LTLA + spec date", 
    "d28.S.ss.2020-12-16.2020-12-31.r10.LTLA:date.sgtf:tstop", c(0, 4), "Hazard ratio for SGTF, 15-31 Dec")

pl5 = plot_hazard("sgtf", 
    "x Sens 1-15 Jan: SGTF + spl age + spl IMD | LTLA + spec date", 
    "d28.S.ss.2021-01-01.2021-01-15.r10.LTLA:date.sgtf:tstop", c(0, 4), "Hazard ratio for SGTF, 1-15 Jan")

plot_grid(pl1, pl2, pl3, pl4, pl5, nrow = 1)

# 4. Depending on the results of the above you might then investigate a model including both a time-varying HR and an interaction 
# between SGTV and specimen date (and possibly even a 3-way interaction between SGTF, time and specimen date)


# SGTF:region interaction....
dataSX[, sgtf_ee := ifelse(NHSER_name == "East of England",          sgtf, 0)]
dataSX[, sgtf_ld := ifelse(NHSER_name == "London",                   sgtf, 0)]
dataSX[, sgtf_ml := ifelse(NHSER_name == "Midlands",                 sgtf, 0)]
dataSX[, sgtf_ne := ifelse(NHSER_name == "North East and Yorkshire", sgtf, 0)]
dataSX[, sgtf_nw := ifelse(NHSER_name == "North West",               sgtf, 0)]
dataSX[, sgtf_se := ifelse(NHSER_name == "South East",               sgtf, 0)]
dataSX[, sgtf_sw := ifelse(NHSER_name == "South West",               sgtf, 0)]

do_cox("SGTF:region + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ 
        sgtf_ee + sgtf_ld + sgtf_ml + sgtf_ne + sgtf_nw + sgtf_se + sgtf_sw + 
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataSX, c("LTLA_name", "specimen_date"))

lrt("SGTF + spl age + spl IMD | LTLA + spec date", "SGTF:region + spl age + spl IMD | LTLA + spec date")

model = load_model("SGTF:region + spl age + spl IMD | LTLA + spec date")
coeffs = as.data.table(summary(model)$coefficients[1:7, c(1,3)], keep.rownames = "id")
coeffs[, region := str_remove_all(id, "sgtf_")]
coeffs[, hi := coef + 1.96 * `se(coef)`]
coeffs[, lo := coef - 1.96 * `se(coef)`]
ggplot(coeffs) +
    geom_pointrange(aes(x = region, y = coef, ymin = lo, ymax = hi, colour = region)) +
    geom_hline(aes(yintercept = 0), linetype = "33") +
    labs(x = "Region", y = "Log hazard ratio")


# SGTF:region:time interaction . . . .

dataSX[, sgtf_ee1 := ifelse(NHSER_name == "East of England"          & month(specimen_date) == 11, sgtf, 0)]
dataSX[, sgtf_ee2 := ifelse(NHSER_name == "East of England"          & month(specimen_date) == 12, sgtf, 0)]
dataSX[, sgtf_ee3 := ifelse(NHSER_name == "East of England"          & month(specimen_date) ==  1, sgtf, 0)]
dataSX[, sgtf_ld1 := ifelse(NHSER_name == "London"                   & month(specimen_date) == 11, sgtf, 0)]
dataSX[, sgtf_ld2 := ifelse(NHSER_name == "London"                   & month(specimen_date) == 12, sgtf, 0)]
dataSX[, sgtf_ld3 := ifelse(NHSER_name == "London"                   & month(specimen_date) ==  1, sgtf, 0)]
dataSX[, sgtf_ml1 := ifelse(NHSER_name == "Midlands"                 & month(specimen_date) == 11, sgtf, 0)]
dataSX[, sgtf_ml2 := ifelse(NHSER_name == "Midlands"                 & month(specimen_date) == 12, sgtf, 0)]
dataSX[, sgtf_ml3 := ifelse(NHSER_name == "Midlands"                 & month(specimen_date) ==  1, sgtf, 0)]
dataSX[, sgtf_ne1 := ifelse(NHSER_name == "North East and Yorkshire" & month(specimen_date) == 11, sgtf, 0)]
dataSX[, sgtf_ne2 := ifelse(NHSER_name == "North East and Yorkshire" & month(specimen_date) == 12, sgtf, 0)]
dataSX[, sgtf_ne3 := ifelse(NHSER_name == "North East and Yorkshire" & month(specimen_date) ==  1, sgtf, 0)]
dataSX[, sgtf_nw1 := ifelse(NHSER_name == "North West"               & month(specimen_date) == 11, sgtf, 0)]
dataSX[, sgtf_nw2 := ifelse(NHSER_name == "North West"               & month(specimen_date) == 12, sgtf, 0)]
dataSX[, sgtf_nw3 := ifelse(NHSER_name == "North West"               & month(specimen_date) ==  1, sgtf, 0)]
dataSX[, sgtf_se1 := ifelse(NHSER_name == "South East"               & month(specimen_date) == 11, sgtf, 0)]
dataSX[, sgtf_se2 := ifelse(NHSER_name == "South East"               & month(specimen_date) == 12, sgtf, 0)]
dataSX[, sgtf_se3 := ifelse(NHSER_name == "South East"               & month(specimen_date) ==  1, sgtf, 0)]
dataSX[, sgtf_sw1 := ifelse(NHSER_name == "South West"               & month(specimen_date) == 11, sgtf, 0)]
dataSX[, sgtf_sw2 := ifelse(NHSER_name == "South West"               & month(specimen_date) == 12, sgtf, 0)]
dataSX[, sgtf_sw3 := ifelse(NHSER_name == "South West"               & month(specimen_date) ==  1, sgtf, 0)]

do_cox("SGTF:region:month + spl age + spl IMD | LTLA + spec date", Surv(time, status) ~ 
        sgtf_ee1 + sgtf_ld1 + sgtf_ml1 + sgtf_ne1 + sgtf_nw1 + sgtf_se1 + sgtf_sw1 + 
        sgtf_ee2 + sgtf_ld2 + sgtf_ml2 + sgtf_ne2 + sgtf_nw2 + sgtf_se2 + sgtf_sw2 + 
        sgtf_ee3 + sgtf_ld3 + sgtf_ml3 + sgtf_ne3 + sgtf_nw3 + sgtf_se3 + sgtf_sw3 + 
        rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataSX, c("LTLA_name", "specimen_date"))

lrt("SGTF + spl age + spl IMD | LTLA + spec date", "SGTF:region:month + spl age + spl IMD | LTLA + spec date")

model = load_model("SGTF:region:month + spl age + spl IMD | LTLA + spec date")
coeffs = as.data.table(summary(model)$coefficients[1:21, c(1,3)], keep.rownames = "id")
coeffs[, region := str_remove_all(id, "sgtf_|[0-9]")]
coeffs[, month := as.numeric(str_extract(id, "[0-9]"))]
coeffs[, hi := coef + 1.96 * `se(coef)`]
coeffs[, lo := coef - 1.96 * `se(coef)`]
ggplot(coeffs) +
    geom_pointrange(aes(x = month, y = coef, ymin = lo, ymax = hi, colour = region), position = position_dodge(width = 0.8)) +
    geom_hline(aes(yintercept = 0), linetype = "33") +
    labs(x = "Month (1 = Nov, 2 = Dec, 3 = Jan)", y = "Log hazard ratio")



###########################
# P_VOC SURVIVAL ANALYSES #
###########################

# 1. MAIN ANALYSES

# 1a. WHOLE DATA RANGE
dataP = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0)
do_cox("d28.pVOC.ss.2020-09-01..r10.LTLA:date.", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataP, c("LTLA_name", "specimen_date"))

# 1b. POST-NOV 1 ONLY
dataPn = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
do_cox("d28.pVOC.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataPn, c("LTLA_name", "specimen_date"))

# 1c. 60 DAY FOLLOWUP P_VOC
dataP60 = model_data(cd, criterion = "under30CT", remove_duplicates = TRUE, death_cutoff = 60, reg_cutoff = 10, P_voc = 0, date_min = "2020-11-01")
do_cox("d60.pVOC.ss.2020-11-01..r10.LTLA:date.", Surv(time, status) ~ p_voc + rcs(age, nk = 3) + sex + rcs(imd, nk = 3) + eth_cat + res_cat + strata(stratum), dataP60, c("LTLA_name", "specimen_date"))


# 2. TIME-P_VOC
event_times = dataPn[!is.na(death_date), sort(unique(as.numeric(death_date - specimen_date)))]
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
        theme(legend.position = c(0.1, 0.9)) +
        ylim(ylimits)
    
    # Uncomment if want to extract time-varying hazards
    # return (haz)
}

# haz = plot_hazard("sgtf", "SGTF + lin age + lin IMD | LTLA + spec date", "Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | LTLA + spec date", c(0, NA))
# haz[t == 1, .(exp(x.ct), exp(x.lo), exp(x.hi))]

theme_set(theme_cowplot(font_size = 10))

hp_sgtf1 = plot_hazard("sgtf", "SGTF + lin age + lin IMD | LTLA + spec date", "Time-SGTF interaction term: SGTF + SGTF:tstop + lin age + lin IMD | LTLA + spec date", c(0, NA), "Hazard ratio for SGTF")
ggsave("./output/time_varying_sgtf_updated.pdf", hp_sgtf1, width = 15, height = 10, units = "cm", useDingbats = FALSE)
ggsave("./output/time_varying_sgtf_updated.png", hp_sgtf1, width = 15, height = 10, units = "cm")

hp_sgtf = plot_hazard("sgtf", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.0sgtf:tstop", "d28.SGTF.ss.2020-11-01..r10.LTLA:date.sgtf:tstop", c(0, 5), "Hazard ratio for SGTF")
hp_voc = plot_hazard("p_voc", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.0p_voc:tstop", "d28.pVOC.ss.2020-11-01..r10.LTLA:date.p_voc:tstop", c(0, 5), expression("Hazard ratio for"~p[VOC]))
hp_compare = plot_grid(hp_sgtf, hp_voc, nrow = 1, labels = letters, label_size = 10)
hp_compare
ggsave("./output/time_varying_compare.pdf", hp_compare, width = 20, height = 10, units = "cm", useDingbats = FALSE)
ggsave("./output/time_varying_compare.png", hp_compare, width = 20, height = 10, units = "cm")


# Effects of SGTF

summaries = fread("./output/model_summary.csv")
summaries = summaries[str_count(model_id, "\\.") == 7]
summaries[, c("m_deaths", "m_marker", "m_ageimd", "m_startdate", "m_enddate", "m_reg_cutoff", "m_strata", "m_xvars") := tstrsplit(model_id, "\\.")]

summaries
summaries
?str_split

summaries = summaries[!model_id %like% "\\(no interaction\\)"]

model_id_replacements = c(
    "spec date" = "date",
    "spec week" = "week",
    "interaction term" = "interaction",
    "\\^2" = "²",
    "\\(" = "",
    "\\)" = "",
    "\\+ spl age \\+ spl IMD" = "(spline)",
    "\\+ lin age \\+ lin IMD " = "",
    "Sensitivity: SGTF \\+ NHSE \\* week" = "SGTF (adjust for NHSE × week)",
    "^p_voc (.*)$" = "p_voc \\1 [1 Sep]",
    "^Sens 1 Nov: (.*)$" = "\\1 [1 Nov]",
    "^Sensitivity: Sep onwards (.*)$" = "\\1 [1 Sep]",
    "^Sensitivity: Oct onwards (.*)$" = "\\1 [1 Oct]",
    "^Sensitivity: Nov onwards (.*)$" = "\\1 [1 Nov]",
    "^Sensitivity: Dec onwards (.*)$" = "\\1 [1 Dec]",
    "^Sensitivity: Prevalence cutoff (.*)$" = "\\1 [by LTLA]",
    "^Censoring: (SGTF|p_voc)_?0?([0-9]+) (.*)$" = "\\1 (\\2-day followup) \\3",
    "999-day" = "-unlimited",
    "Sensitivity: no registration cutoff (.*)" = "\\1 [1 Nov+]",
    "^Sensitivity: max specimen date 20 Dec (.*)$" = "\\1 [1 Nov-20 Dec]",
    "Sensitivity asymptomatic: SGTF" = "SGTF (symptom status)",
    "^(.* interaction)" = "SGTF (\\1)"
)

for (repl in seq_along(model_id_replacements)) {
    summaries[, model_id := str_replace_all(model_id, names(model_id_replacements)[repl], model_id_replacements[repl])]
}
summaries[!model_id %like% "\\[", model_id := paste0(model_id, " [1 Nov]")]
summaries = summaries[!(model_id %like% "followup" & model_id %like% "week")]


pl_effects = ggplot(summaries[parameter %in% c("sgtf", "p_voc") & !model_id %like% "Time.{0,2}-(SGTF|p_voc)"]) +
    geom_pointrange(aes(y = model_id, x = HR, xmin = HR.lo95, xmax = HR.hi95), fatten = 0.1) +
    geom_vline(aes(xintercept = 1), linetype = "33", size = 0.25) +
    labs(x = "Hazard ratio", y = NULL) +
    theme(panel.background = element_rect(fill = "#f0f0f0"),
        panel.grid.major = element_line(colour = "#ffffff", size = 0.5), 
        plot.margin = unit(c(0.1, 0.3, 0.1, 0.1), "cm"))
pl_effects

ggsave("./output/summary.pdf", pl, width = 20, height = 20, units = "cm")

# Exploratory data plots
dataS[, SGTF_label := factor(ifelse(sgtf == 1, "SGTF", "Other"), levels = c("SGTF", "Other"))]
plot_samples = ggplot(dataS[, .N, keyby = .(specimen_date, SGTF_label)]) +
    geom_col(aes(specimen_date, N, fill = SGTF_label), position = "stack") +
    labs(x = "Specimen date", y = "Samples", fill = NULL) +
    theme(legend.position = c(0.1, 0.9))

plot_deaths = ggplot(dataS[died == TRUE, .N, keyby = .(specimen_date, SGTF_label)]) +
    geom_col(aes(specimen_date, N, fill = SGTF_label), position = "stack") +
    labs(x = "Specimen date", y = "Deaths", fill = NULL) +
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

