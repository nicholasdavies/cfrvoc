library(data.table)
library(lubridate)
library(ogwrangler)
library(readxl)
library(forecast)
library(zoo)

# Available beds, staff absences, etc
dirs = list.files(path = "~/Documents/uk_covid_data_sensitive/folder_20127408/Archive/", include.dirs = TRUE)
alldat = NULL;
for (dir in dirs) {
    cat(".")
    searchdir = paste0("~/Documents/uk_covid_data_sensitive/folder_20127408/Archive/", dir);
    file = list.files(path = searchdir, pattern = "^Covid.sitrep.report.*(xlsx|XLSX|xlsm)$")
    dat = as.data.table(read_excel(paste0(searchdir, "/", file), "R Data", skip = 1))
    dat = dat[, .(date = as.Date(period), region = Region, sitecode = `Site/Org Code`, type = `Organisation Type`,
        mv_occ_cov = SIT032_OccupiedCov, mv_occ_sus = SIT032_OccupiedSuspected, mv_occ_non = SIT032_OccupiedNonCovNS, mv_unocc = SIT032_Unoccupied, # mechanical ventilation beds
        ni_occ_cov = SIT033_OccupiedCov, ni_occ_sus = SIT033_OccupiedSuspected, ni_occ_non = SIT033_OccupiedNonCovNS, ni_unocc = SIT033_Unoccupied, # noninvasive ventilation beds
        os_occ_cov = SIT034_OccupiedCov, os_occ_sus = SIT034_OccupiedSuspected, os_occ_non = SIT034_OccupiedNonCovNS, os_unocc = SIT034_Unoccupied, # oxygenation support beds
        ao_occ_cov = SIT058_OccupiedCov, ao_occ_sus = SIT058_OccupiedSuspected, ao_occ_non = SIT058_OccupiedNonCovNS, ao_unocc = SIT058_Unoccupied, # any other beds
        medstaff_absent = SIT048_TotalAbsent, nursing_absent = SIT049_TotalAbsent)]
    dat[, occ_all :=
        mv_occ_cov + mv_occ_sus + mv_occ_non +
        ni_occ_cov + ni_occ_sus + ni_occ_non +
        os_occ_cov + os_occ_sus + os_occ_non +
        ao_occ_cov + ao_occ_sus + ao_occ_non]
    dat[, occ_cov :=
        mv_occ_cov + mv_occ_sus +
        ni_occ_cov + ni_occ_sus +
        os_occ_cov + os_occ_sus +
        ao_occ_cov + ao_occ_sus]
    
    # Correct dates improperly entered as 2020 when they should be 2021
    dat[date < "2020-05-01", date := date + 366]

    alldat = rbind(alldat, dat)
}

# 1. AMALGAMATE BY NHS REGION

CreateCache()
etrust = fread("~/Dropbox/uk_covid_data/etrust/etrust.csv", header = FALSE)[, .(sitecode = V1, postcode = V10)]

alldat2 = merge(alldat[type == "Acute Trust"], etrust, by = "sitecode", all.x = TRUE)
alldat2[, ladcd := ogpost(postcode, geo = "lad")]
alldat2[, nhscd := ogpost(postcode, geo = "nhser")]

# Remove outliers in absences
sitecodes = alldat2[, unique(sitecode)]
for (sc in sitecodes)
{
    TS = alldat2[sitecode == sc, ts(nursing_absent)];
    cleaned = tsclean(TS);
    alldat2[sitecode == sc, nursing_absent := as.numeric(cleaned)];

    TS = alldat2[sitecode == sc, ts(medstaff_absent)];
    cleaned = tsclean(TS);
    alldat2[sitecode == sc, medstaff_absent := as.numeric(cleaned)];
}

w = alldat2[, .(
    mv_pressure = 1 - sum(mv_unocc) / (sum(mv_occ_cov + mv_occ_sus + mv_occ_non + mv_unocc)),
    ni_pressure = 1 - sum(ni_unocc) / (sum(ni_occ_cov + ni_occ_sus + ni_occ_non + ni_unocc)),
    os_pressure = 1 - sum(os_unocc) / (sum(os_occ_cov + os_occ_sus + os_occ_non + os_unocc)),
    ao_pressure = 1 - sum(ao_unocc) / (sum(ao_occ_cov + ao_occ_sus + ao_occ_non + ao_unocc)),
    n_beds = sum(mv_occ_cov + mv_occ_sus + mv_occ_non + mv_unocc +
                 ni_occ_cov + ni_occ_sus + ni_occ_non + ni_unocc +
                 os_occ_cov + os_occ_sus + os_occ_non + os_unocc +
                 ao_occ_cov + ao_occ_sus + ao_occ_non + ao_unocc),
    medstaff_abs = sum(medstaff_absent),
    nursing_abs = sum(nursing_absent)
    ), keyby = .(date, nhscd)]

w[, medstaff_abs_per_bed := medstaff_abs / mean(n_beds), by = nhscd]
w[, nursing_abs_per_bed := nursing_abs / mean(n_beds), by = nhscd]

date_max = max(w$date)
nhs_codes = w[, unique(nhscd)]
newdata = data.table(date = rep(date_max + 1:28, each = length(nhs_codes)), nhscd = rep(nhs_codes, 28))
w = merge(w, newdata, all = T)

# Plot data
w2 = melt(w, id.vars = 1:2)
ggplot(w2) + geom_line(aes(x = date, y = value, colour = nhscd)) + facet_wrap(~variable, scales = "free")

# Smooth and extend data
w2[, value := na.fill(value, fill = "extend"), by = .(nhscd, variable)]
w2[, value := rollmean(value, 7, fill = "extend"), by = .(nhscd, variable)]
ggplot(w2) + geom_line(aes(x = date, y = value, colour = nhscd)) + facet_wrap(~variable, scales = "free")

# # Standardise absences
# w2[variable %in% c("medstaff_abs", "nursing_abs"), value := (value - mean(value)) / sd(value), by = .(nhscd, variable)]
# w2[, value := zoo::rollmean(value, 7, fill = "extend"), by = .(nhscd, variable)]
# ggplot(w2) + geom_line(aes(x = date, y = value, colour = nhscd)) + facet_wrap(~variable, scales = "free")

# Write
pressure = dcast(w2, nhscd + date ~ variable)
fwrite(pressure, "~/Documents/uk_covid_data_sensitive/pressure-2021-01-31.csv")





# 2. AMALGAMATE BY LTLA
ltla_site = fread("~/Documents/uk_covid_data_sensitive/spells_by_site_and_LAD_raw.csv");
ltla_site = ltla_site[Der_Postcode_Dist_Unitary_Auth %like% "^E", 
    .(sitecode = Der_Provider_Site_Code, spells), by = .(ltla = Der_Postcode_Dist_Unitary_Auth)]
alldat3 = copy(alldat[type == "Acute Trust" & nchar(sitecode) == 5])

# Remove outliers in absences
sitecodes = alldat3[, unique(sitecode)]
for (sc in sitecodes)
{
    TS = alldat3[sitecode == sc, ts(nursing_absent)];
    cleaned = tsclean(TS);
    alldat3[sitecode == sc, nursing_absent := as.numeric(cleaned)];

    TS = alldat3[sitecode == sc, ts(medstaff_absent)];
    cleaned = tsclean(TS);
    alldat3[sitecode == sc, medstaff_absent := as.numeric(cleaned)];
}

w = alldat3[, .(
    mv_pressure = 1 - sum(mv_unocc) / (sum(mv_occ_cov + mv_occ_sus + mv_occ_non + mv_unocc)),
    ni_pressure = 1 - sum(ni_unocc) / (sum(ni_occ_cov + ni_occ_sus + ni_occ_non + ni_unocc)),
    os_pressure = 1 - sum(os_unocc) / (sum(os_occ_cov + os_occ_sus + os_occ_non + os_unocc)),
    ao_pressure = 1 - sum(ao_unocc) / (sum(ao_occ_cov + ao_occ_sus + ao_occ_non + ao_unocc)),
    n_beds = sum(mv_occ_cov + mv_occ_sus + mv_occ_non + mv_unocc +
                 ni_occ_cov + ni_occ_sus + ni_occ_non + ni_unocc +
                 os_occ_cov + os_occ_sus + os_occ_non + os_unocc +
                 ao_occ_cov + ao_occ_sus + ao_occ_non + ao_unocc),
    medstaff_abs = sum(medstaff_absent),
    nursing_abs = sum(nursing_absent)
    ), keyby = .(date, sitecode)]

w[, medstaff_abs_per_bed := medstaff_abs / mean(n_beds), by = sitecode]
w[, nursing_abs_per_bed := nursing_abs / mean(n_beds), by = sitecode]

# "E07000005" "E07000004" "E07000006" "E07000007" "E06000053"
ltlas = ltla_site[, unique(ltla)]
w2 = NULL
for (l in ltlas)
{
    w_ltla = merge(w, ltla_site[ltla == l], by = "sitecode")
    w_ltla = w_ltla[, lapply(.SD, weighted.mean, w = spells, na.rm = T), .SDcols = c(3:6, 10:11), keyby = .(date, ltla)]
    w2 = rbind(w2, w_ltla)
}

w = w2

date_max = max(w$date)
ltla_codes = w[, unique(ltla)]
newdata = data.table(date = rep(date_max + 1:28, each = length(ltla_codes)), ltla = rep(ltla_codes, 28))
w = merge(w, newdata, all = T)

# Plot data
w2 = melt(w, id.vars = 1:2)
ggplot(w2) + geom_line(aes(x = date, y = value, colour = ltla)) + facet_wrap(~variable, scales = "free") + theme(legend.position = "none")

# Smooth and extend data
w2 = w2[variable %like% "abs" | variable == "mv_pressure"]
w2[, value := na.fill(value, fill = "extend"), by = .(ltla, variable)]
w2[, value := rollmean(value, 7, fill = "extend"), by = .(ltla, variable)]
ggplot(w2) + geom_line(aes(x = date, y = value, colour = ltla)) + facet_wrap(~variable, scales = "free") + theme(legend.position = "none")

# Write
pressure_ltla = dcast(w2, ltla + date ~ variable)
fwrite(pressure_ltla, "~/Documents/uk_covid_data_sensitive/pressure_ltla-2021-01-31.csv")
