# Compare potential missingness models

library(data.table)
library(lubridate)
library(survival)
library(rms)
library(binom)
library(qs)
library(ggplot2)
library(ResourceSelection)
library(glmx)
library(writexl)
library(cowplot)

source("./hazard_data.R")

cd = complete_data("20210225")
dataW = model_data(cd, "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, 
    date_min = "2020-11-01", date_max = "2100-01-01", prevalence_cutoff = FALSE, sgtfv_cutoff = 0, keep_missing = TRUE, death_type = "all")

# Load data
dataW[, missing := ifelse(is.na(sgtf), 1, 0)]
dataW[, specimen_t := as.numeric(specimen_date - ymd("2020-11-01"))]
dataW[, LTLA_name2 := LTLA_name]
dataW[LTLA_name2 == "Isles of Scilly", LTLA_name2 := "Cornwall"]
dataW[, specimen_week_f := factor(specimen_week)]

missingness = function(name, formula, binomial_link, pfunc, data)
{
    cat("Fitting model...\n");
    missing_model = glm(formula, family = binomial(link = binomial_link), data = data, control = list(trace = TRUE));
    
    cat("Calculating output...\n");
    summ = summary(missing_model);
    aic = extractAIC(missing_model);
    data[, missing_lin := predict(missing_model, newdata = .SD)];
    data[, missing_p := pfunc(missing_lin)];
    data[, missing_hinkley := missing_lin^2];
    
    # cat("Hinkley test...\n");
    # formula_hinkley = as.formula(paste(paste(deparse(formula, width.cutoff = 500), collapse = " "), "+ missing_hinkley"));
    # missing_hinkley = glm(formula_hinkley, family = binomial(link = binomial_link), data = data, control = list(trace = TRUE));
    # hinkley = summary(missing_hinkley)$coefficients["missing_hinkley", ];
    # 
    # cat("Hosmer-Lemeshow test...\n");
    # hoslem = hoslem.test(data$missing, data$missing_p, 10);
    
    cat("Hosmer-Lemeshow plot...\n");
    step = 0.025;
    data[, p_bin := cut(missing_p, seq(0.0, 1, 0.025))];
    missing_fit = data[, .(m_yes = sum(missing == 1), m_no = sum(missing == 0)), keyby = p_bin];
    missing_fit[, c("mean", "lower", "upper") := binom.confint(m_yes, m_yes + m_no, method = "exact")[, 4:6]];
    missing_fit[, bin_x := as.numeric(str_remove_all(p_bin, "\\(|,[0-9.]+\\]")) + step / 2];

    cat("Weights plot...\n");
    dataN = data[, .N, by = .(miss = as.factor(missing), wt = floor(1 / (1 - missing_p)))]

    summ$deviance.resid = NULL
    
    results = list(summ = summ, aic = aic, 
        #hinkley = hinkley, hoslem = hoslem, 
        hoslem_data = missing_fit, weights_hist = dataN);
    
    qsave(results, paste0("./output/", name, "_missing.qs"))

    return (results)
}

wC1 = missingness("wC1", missing ~ rcs(age, 5) + sex + rcs(imd, 3) + eth_cat + res_cat * asymptomatic + NHSER_name * specimen_week_f,
    "cauchit", pcauchy, dataW)

wL1 = missingness("wL1", missing ~ rcs(age, 5) + sex + rcs(imd, 3) + eth_cat + res_cat * asymptomatic + NHSER_name * specimen_week_f,
    "logit", plogis, dataW)

wR1 = missingness("wR1", missing ~ rcs(age, 5) + sex + rcs(imd, 3) + eth_cat + res_cat * asymptomatic + NHSER_name * specimen_week_f,
    gosset(nu = 4), function(x) pt(x, df = 4), dataW)

wC2 = missingness("wC2", missing ~ rcs(age, 5) + sex + rcs(imd, 3) + eth_cat + res_cat * asymptomatic + LTLA_name + rcs(specimen_t, 5),
    "cauchit", pcauchy, dataW)

wL2 = missingness("wL2", missing ~ rcs(age, 5) + sex + rcs(imd, 3) + eth_cat + res_cat * asymptomatic + LTLA_name + rcs(specimen_t, 5),
    "logit", plogis, dataW)

wR2 = missingness("wR2", missing ~ rcs(age, 5) + sex + rcs(imd, 3) + eth_cat + res_cat * asymptomatic + LTLA_name + rcs(specimen_t, 5),
    gosset(nu = 4), function(x) pt(x, df = 4), dataW)

diagnostic = function(w)
{
    ph = ggplot(w$hoslem_data) +
        geom_pointrange(aes(x = bin_x, y = mean, ymin = lower, ymax = upper), fatten = 1) +
        geom_abline(aes(intercept = 0, slope = 1), linetype = "22") +
        labs(x = "Modelled missingness", y = "Observed missingness") +
        xlim(0, 1) + ylim(0, 1)
    
    pw = ggplot(w$weights_hist) + 
        geom_col(aes(wt, N, fill = miss)) +
        scale_colour_manual(values = c("0" = "blue", "1" = "red")) +
        scale_y_log10() +
        labs(x = "Weight", y = "Count", fill = "Missing")
    
    return (list(ph = ph, pw = pw))
}

dwC1 = diagnostic(wC1)
dwL1 = diagnostic(wL1)
dwR1 = diagnostic(wR1)
theme_set(theme_cowplot(font_size = 11))
pmiss = plot_grid(
    dwC1$ph, dwC1$pw, 
    dwR1$ph, dwR1$pw,
    dwL1$ph, dwL1$pw,
    labels = c("a", "", "b", "", "c", ""), label_size = 10, nrow = 3, rel_widths = c(1, 1.2))
ggsave("./output/missingness.pdf", pmiss, width = 18, height = 18, units = "cm")

diagnostic(wC2)
diagnostic(wR2)
diagnostic(wL2)

write_xlsx(list(
    A_qq_plot = wC1$hoslem_data, A_weights_plot = wC1$weights_hist[order(miss, wt)],
    B_qq_plot = wR1$hoslem_data, B_weights_plot = wR1$weights_hist[order(miss, wt)],
    C_qq_plot = wL1$hoslem_data, C_weights_plot = wL1$weights_hist[order(miss, wt)]
    ), "./manuscript/sdE_missingness_model.xlsx")

