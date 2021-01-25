# Bayesian misclassification model
library(data.table)
library(RCppMCMC) # remotes::install_github("nicholasdavies/RCppMCMC")
library(ggplot2)
library(lubridate)
library(cowplot)
library(stringr)
library(binom)

source("./phe_data.R")

# Helper functions
nameval = function(names, values)
{
    x = values;
    names(x) = names;
    return (x)
}

logistic = function(x, a, b)
{
    xx = a * (x - b);
    ifelse(xx < -200, 0, ifelse(xx > 200, 1, exp(xx) / (1 + exp(xx))))
}
    
cpp_vec = function(x) paste("{", paste(x, collapse = ", "), "}")
cpp_bbinom = 
'auto bbinom = [](double k, double n, double mode, double conc)
{
    auto lgamma = [](double x) { return gsl_sf_lngamma(x); };
    double a = mode * (conc - 2) + 1;
    double b = (1 - mode) * (conc - 2) + 1;
    
    return (lgamma(n + 1) + lgamma(k + a) + lgamma(n - k + b) + lgamma(a + b))
        - (lgamma(k + 1) + lgamma(n - k + 1) + lgamma(n + a + b) + lgamma(a) + lgamma(b));
};';
cpp_logistic = 
'auto logistic = [](double x, double a, double b)
{
    double xx = a * (x - b);
    if (xx < -200) return 0.0;
    if (xx > 200)  return 1.0;
    return exp(xx) / (1 + exp(xx));
};';

build_model = function(data, falsepos)
{
    groups = data[, unique(group)];
    n_groups = length(groups);
    
    # Set up priors
    prior_intercept_name = paste0("intercept", 1:n_groups);
    prior_intercept_dist = rep("N 0 1000", n_groups);
    prior_intercept_code = paste0("double intercept[] = {", paste(prior_intercept_name, collapse = ", "), "};");
    prior_intercept = nameval(prior_intercept_name, prior_intercept_dist);
    
    prior_conc_name = paste0("conc", 1:n_groups);
    prior_conc_dist = rep("N 0 500 T 2 2000", n_groups);
    prior_conc_code = paste0("double conc[] = {", paste(prior_conc_name, collapse = ", "), "};");
    prior_conc = nameval(prior_conc_name, prior_conc_dist);
    
    if (is.character(falsepos)) {
        prior_falsepos_name = paste0("falsepos", 1:n_groups);
        prior_falsepos_code = paste0("double falsepos[] = {", paste(prior_falsepos_name, collapse = ", "), "};");
        prior_falsepos_dist = rep(falsepos, n_groups);
        prior_falsepos = nameval(prior_falsepos_name, prior_falsepos_dist);
    } else {
        prior_falsepos_code = paste0("double falsepos[] = ", cpp_vec(rep_len(falsepos, n_groups)), ";");
        prior_falsepos = NULL;
    }
    
    # Make priors
    priors = c(
        slope = "N 0 1 T 0 100",
        prior_intercept,
        prior_conc,
        prior_falsepos
    );
    
    # Make code
    code = glue::glue(
        cpp_bbinom,
        cpp_logistic,
        prior_intercept_code,
        prior_conc_code,
        prior_falsepos_code,
        'std::vector<double> t = ${cpp_vec(as.numeric(data$specimen_date - ymd("2020-01-01")))};',
        'std::vector<double> s = ${cpp_vec(data$sgtf)};',
        'std::vector<double> f = ${cpp_vec(data$other)};',
        'std::vector<unsigned int> r = ${cpp_vec(match(data$group, groups) - 1)};',
        'for (unsigned int i = 0; i < t.size(); ++i) {',
        '    double frequency = logistic(t[i], slope, intercept[r[i]]);',
        '    double predicted = frequency + (1 - frequency) * falsepos[r[i]];',
        '    ll += bbinom(s[i], s[i] + f[i], predicted, conc[r[i]]);',
        '}',
    .sep = "\n", .open = "${", .close = "}")
    
    make_model("bbinom", priors, code)
}

# Load data
sgtf = sgtf_counts("20210121", "NHSER_name")[group != ""]

# Shared slope model
sgtf_model = build_model(sgtf, "B 1.5 15")
results = RCppMCMC(sgtf_model, 30000, 1000, threads = 6, verbose = TRUE)
setDT(results)
results

# Independent slopes
groups = sgtf[, unique(group)]
results_1 = NULL;
ll = NULL;
for (i in seq_along(groups))
{
    sgtf_model_1 = build_model(sgtf[group == groups[i]], "B 1.5 15")
    results_1_0 = RCppMCMC(sgtf_model_1, 10000, 5000, threads = 6, verbose = TRUE)
    setnames(results_1_0, "slope", paste0("slope", i))
    setnames(results_1_0, "intercept1", paste0("intercept", i))
    setnames(results_1_0, "conc1", paste0("conc", i))
    setnames(results_1_0, "falsepos1", paste0("falsepos", i))
    ll = cbind(ll, results_1_0$ll)
    if (is.null(results_1)) {
        results_1 = results_1_0
    } else {
        results_1 = cbind(results_1, results_1_0[, 5:ncol(results_1_0)])
    }
}
setDT(results_1)

# Extraction of results
extract_results = function(results, date_min = "2020-09-01", date_max = "2021-01-31")
{
    res_data_all = NULL;
    
    res = melt(results, id.vars = numeric(), measure.vars = 5:ncol(results));
    res[, group := groups[as.numeric(str_remove_all(variable, "[a-z]*"))]];
    res[, var := str_remove_all(variable, "[0-9]*")];
    
    extract = function(res, gr, varname)
    {
        ex = res[group == gr & var == varname, value];
        if (length(ex) > 0)
            return (ex);
        return (res[is.na(group) & var == varname, value]);
    }
    
    set.seed(12345);
    for (i in seq_along(groups))
    {
        # Extract posterior for this group
        res_data = data.table(group = groups[i], date = ymd(date_min) + 0:as.numeric(ymd(date_max) - ymd(date_min)));
        slope = extract(res, groups[i], "slope");
        intercept = extract(res, groups[i], "intercept");
        conc = extract(res, groups[i], "conc");
        falsepos = extract(res, groups[i], "falsepos");

        nsamp = 10000;
        rows = sample(length(slope), nsamp, replace = TRUE);
        slope = slope[rows];
        intercept = intercept[rows];
        conc = conc[rows];
        falsepos = falsepos[rows];
        
        pred0 = matrix(0, nrow = nrow(res_data), ncol = nsamp);
        pred  = matrix(0, nrow = nrow(res_data), ncol = nsamp);
        predr = matrix(0, nrow = nrow(res_data), ncol = nsamp);
        sgtfv = matrix(0, nrow = nrow(res_data), ncol = nsamp);
        fpr   = matrix(0, nrow = nrow(res_data), ncol = nsamp);
        
        for (j in 1:nsamp)
        {
            pred0[, j] = logistic(as.numeric(res_data$date - ymd("2020-01-01")), slope[j], intercept[j]);
            sgtfv[, j] = pred0[, j] / (pred0[, j] + (1 - pred0[, j]) * falsepos[j]);
            pred[, j] = pred0[, j] + (1 - pred0[, j]) * falsepos[j];
            predr[, j] = rbeta(length(pred[, j]), pred[, j] * (conc[j] - 2) + 1, (1 - pred[, j]) * (conc[j] - 2) + 1);
            fpr[, j] = rbeta(length(falsepos[j]), falsepos[j] * (conc[j] - 2) + 1, (1 - falsepos[j]) * (conc[j] - 2) + 1);
        }
        
        res_data[, vlo := apply(pred0, 1, function(x) quantile(x, 0.025))];
        res_data[, vmd := apply(pred0, 1, function(x) quantile(x, 0.500))];
        res_data[, vhi := apply(pred0, 1, function(x) quantile(x, 0.975))];
        res_data[, mlo := apply(pred,  1, function(x) quantile(x, 0.025))];
        res_data[, mmd := apply(pred,  1, function(x) quantile(x, 0.500))];
        res_data[, mhi := apply(pred,  1, function(x) quantile(x, 0.975))];
        res_data[, rlo := apply(predr, 1, function(x) quantile(x, 0.025))];
        res_data[, rmd := apply(predr, 1, function(x) quantile(x, 0.500))];
        res_data[, rhi := apply(predr, 1, function(x) quantile(x, 0.975))];
        res_data[, sgtfv := apply(sgtfv, 1, mean)];
        res_data[, fplo := quantile(falsepos, 0.025)];
        res_data[, fpmd := quantile(falsepos, 0.500)];
        res_data[, fphi := quantile(falsepos, 0.975)];
        res_data[, rfplo := apply(fpr, 1, function(x) quantile(x, 0.025))];
        res_data[, rfpmd := apply(fpr, 1, function(x) quantile(x, 0.500))];
        res_data[, rfphi := apply(fpr, 1, function(x) quantile(x, 0.975))];
        
        res_data_all = rbind(res_data_all, res_data);
    }
    
    res_data_all
}

# Results by NHS region
w = extract_results(results_1)
fwrite(w, "./sgtf_voc.csv")

w2 = copy(sgtf)
w2[, c("mean", "lower", "upper") := binom.confint(w2$sgtf, w2$sgtf + w2$other, method = "exact")[4:6]]
w2 = w2[specimen_date >= "2020-09-01"]

ww = merge(w, w2[, .(date = specimen_date, group, dmean = mean, dlo = lower, dhi = upper)], by = c("date", "group"), all = TRUE)

theme_set(theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

cross_dates = ww[, date[which.min(abs(vmd - fpmd))], by = group]
cross_dates

ggplot(ww) + 
    geom_ribbon(aes(date, ymin = rfplo, ymax = rfphi, fill = "Modelled non-VOC SGTF"), alpha = 0.4) + 
    geom_line(aes(date, y = rfpmd, colour = "Modelled non-VOC SGTF")) + 
    geom_ribbon(aes(date, ymin = rlo, ymax = rhi, fill = "Modelled SGTF"), alpha = 0.4) + 
    geom_line(aes(date, y = rmd, colour = "Modelled SGTF")) + 
    geom_ribbon(aes(date, ymin = vlo, ymax = vhi, fill = "Modelled VOC"), alpha = 0.4) +
    geom_line(aes(date, vmd, colour = "Modelled VOC")) +
    geom_ribbon(aes(date, ymin = dlo, ymax = dhi, fill = "Observed SGTF"), alpha = 0.4) +
    geom_line(aes(date, dmean, colour = "Observed SGTF")) +
    geom_line(aes(date, sgtfv, colour = "P(VOC|SGTF)")) +
    scale_colour_manual(aesthetics = c("fill", "colour"), values = c("Modelled non-VOC SGTF" = "orange", "Modelled SGTF" = "darkorchid", 
        "Modelled VOC" = "blue", "Observed SGTF" = "black", "P(VOC|SGTF)" = "#008888")) +
    geom_vline(data = cross_dates, aes(xintercept = V1), linetype = "33", size = 0.25) +
    labs(x = "Specimen date", y = NULL, colour = NULL, fill = NULL) +
    theme(legend.position = c(0.35, 0.2)) +
    facet_wrap(~group)

ggsave("./output/misclassification.pdf", width = 20, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/misclassification.png", width = 20, height = 15, units = "cm")

# Avoid chopping off parts of the plot
ww2 = copy(ww)
ww2[, rlo := pmin(0.999, rlo)]
ww2[, rmd := pmin(0.999, rmd)]
ww2[, rhi := pmin(0.999, rhi)]
ww2[, vlo := pmin(0.999, vlo)]
ww2[, vmd := pmin(0.999, vmd)]
ww2[, vhi := pmin(0.999, vhi)]
ww2[, vlo := pmax(0.001, vlo)]
ww2[, vhi := pmax(0.001, vhi)]

ggplot(ww2) + 
    geom_ribbon(aes(date, ymin = rfplo, ymax = rfphi, fill = "Modelled non-VOC SGTF"), alpha = 0.4) + 
    geom_line(aes(date, y = rfpmd, colour = "Modelled non-VOC SGTF")) + 
    geom_ribbon(aes(date, ymin = rlo, ymax = rhi, fill = "Modelled SGTF"), alpha = 0.4) + 
    geom_line(aes(date, y = rmd, colour = "Modelled SGTF")) + 
    geom_ribbon(aes(date, ymin = vlo, ymax = vhi, fill = "Modelled VOC"), alpha = 0.4) +
    geom_line(aes(date, vmd, colour = "Modelled VOC")) +
    geom_ribbon(aes(date, ymin = dlo, ymax = dhi, fill = "Observed SGTF"), alpha = 0.4) +
    geom_line(aes(date, dmean, colour = "Observed SGTF")) +
    geom_line(aes(date, sgtfv, colour = "P(VOC|SGTF)")) +
    scale_colour_manual(aesthetics = c("fill", "colour"), values = c("Modelled non-VOC SGTF" = "orange", "Modelled SGTF" = "darkorchid", 
        "Modelled VOC" = "blue", "Observed SGTF" = "black", "P(VOC|SGTF)" = "#008888")) +
    geom_vline(data = cross_dates, aes(xintercept = V1), linetype = "33", size = 0.25) +
    labs(x = "Specimen date", y = NULL, colour = NULL, fill = NULL) +
    theme(legend.position = c(0.5, 0.15)) +
    facet_wrap(~group) +
    scale_y_continuous(trans = scales::logit_trans(), limits = c(0.001, 0.999), breaks = c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99))

ggsave("./output/misclassification_logit.pdf", width = 20, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/misclassification_logit.png", width = 20, height = 15, units = "cm")


# ggplot(w) + geom_line(aes(x = date, y = sgtfv)) + facet_wrap(~group)
