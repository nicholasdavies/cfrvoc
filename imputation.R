library(data.table)
library(mice)

source("./hazard_data.R")

cd = complete_data("20210201")
dataM = model_data(cd, "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, 
    date_min = "2020-11-01", date_max = "2100-01-01", prevalence_cutoff = FALSE, sgtfv_cutoff = 0, keep_missing = TRUE, death_type = "all")

rd

rd2 = rd[!is.na(sgtf), .(sgtf_under30CT, P2CH1CQ, P2CH2CQ, P2CH3CQ, P2CH4CQ, LTLA_name, UTLA_name)]
m = mice(rd2)
w = as.data.table(complete(m, 1))

ggplot(w[sample(.N, 10000)]) +
    geom_point(aes(x = P2CH1CQ, y = P2CH2CQ, colour = sgtf_under30CT))
hist(w$sgtf_under30CT)

cd



library(smcfcs)

cd = complete_data("20210201")
dataM = model_data(cd, "under30CT", remove_duplicates = TRUE, death_cutoff = 28, reg_cutoff = 10, P_voc = 0, 
    date_min = "2020-11-01", date_max = "2100-01-01", prevalence_cutoff = FALSE, sgtfv_cutoff = 0, keep_missing = TRUE, death_type = "all")
w = smcfcs(dataMr, "coxph", "Surv(time,status)~sgtf+age+sex+imd+eth_cat+res_cat",
    c("", "", "logreg", "", "", "", "", ""))


# 1. split model_data into two bits, one which arranges the data and one which filters the data.
#


# Back to mice

dataMr = dataM[, .(time, status,
    sgtf = as.logical(sgtf), age, sex = factor(sex),
    imd, eth_cat = factor(eth_cat), res_cat = factor(res_cat),
    LTLA_name, specimen_date)]
dataMr[, ch_nelson_aalen := nelsonaalen(.SD, "time", "status")]
#dataMr[, stratum := factor(paste0("LTLA_name", "specimen_date", sep = "|"))]


plot(w$smCoefIter[1,4,])

wm = mice(dataMr)

coxph(Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat, data = complete(wm, 1)) # 1.47
coxph(Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat, data = complete(wm, 2)) # 1.50
coxph(Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat, data = complete(wm, 3)) # 1.42
coxph(Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat, data = complete(wm, 4)) # 1.41
coxph(Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat, data = complete(wm, 5)) # 1.36
coxph(Surv(time, status) ~ sgtf + age + sex + imd + eth_cat + res_cat, data = dataMr[!is.na(sgtf)]) # 1.57
