library(data.table)
library(readxl)
library(lubridate)
library(plyr)

source("./phe_data.R")

prefer = function(a, b)
{
    ifelse (!is.na(a), a, b)
}

# # Look for optimal hospital pressure time: get avg time from pillar 2 test to hospital admission and to death
# delay_admit = d[!is.na(dateadmission_NHSE) & !is.na(specimen_date.x) & !is.na(dod) & pillar == "Pillar 2", 
#     as.numeric(as.Date(dateadmission_NHSE) - as.Date(specimen_date.x))]
# delay_death = d[!is.na(dateadmission_NHSE) & !is.na(specimen_date.x) & !is.na(dod) & pillar == "Pillar 2", 
#     as.numeric(as.Date(dod) - as.Date(specimen_date.x))]
# 
# mean(delay_admit[delay_admit >= 0])
# median(delay_admit[delay_admit >= 0])
# mean(delay_death[delay_death >= 0])
# median(delay_death[delay_death >= 0])

# # Get latitude and longitude by LTLA from ONS Postcode Directory (available from ONS)
# onspd = fread("~/Downloads/ONSPD_AUG_2020_UK/Data/ONSPD_AUG_2020_UK.csv");
# onspd = onspd[is.na(doterm), .(oslaua, lsoa11, msoa11, lat, long)]
# onspd = onspd[lsoa11 %like% "^E"]
# onspd = onspd[, .(lat = mean(lat), long = mean(long)), by = .(oslaua, lsoa11)]
# onspd[, pop := ogwhat(lsoa11, "pop2019")]
# ltla_geo = onspd[, .(lat = weighted.mean(lat, pop), long = weighted.mean(long, pop)), by = oslaua]
# rm(onspd)

# Load and assemble complete data set
complete_data = function(dateid, pressure_delay = 9)
{
    # Note: These files contain personally identifiable information, so they are not included with the repo.
    d_death = phe_deaths(dateid)
    ll = phe_positives(dateid)
    sgtf = phe_sgtf(dateid)
    
    ll[, specimen_date := dmy(specimen_date)]
    sgtf[, specimen_date := ymd(specimen_date)]
    
    d = merge(ll, sgtf, by = c("FINALID", "specimen_date"), all = TRUE)
    d = merge(d, d_death, by.x = "FINALID", by.y = "finalid", all = TRUE)
    
    # sgtfvoc: from misclassification.R
    sgtfvoc = fread("./sgtf_voc.csv")
    d = merge(d, sgtfvoc[, .(specimen_date.x = date, NHSER_name = group, sgtfv)], by = c("specimen_date.x", "NHSER_name"), all.x = TRUE)

    d[, hosp_date := specimen_date.x + pressure_delay]
    
    # pressure: from build_sitrep_data.R
    # pressure = fread("~/Documents/uk_covid_data_sensitive/pressure.csv");
    # d = merge(d, pressure, by.x = c("hosp_date", "NHSER_code"), by.y = c("date", "nhscd"), all.x = TRUE)
    
    d[, data_id := dateid];
    
    # # geography
    # d[, ltla_new := LTLA_code];
    # d[ltla_new %in% c("E07000004", "E07000005", "E07000006", "E07000007"), ltla_new := "E06000060"]; # LADs amalgamated to Buckinghamshire
    # d = merge(d, ltla_geo, by.x = "ltla_new", by.y = "oslaua", all.x = TRUE);
    
    return (d)
}

# sgtf This is the SGTF indicator from OST, SGTF definition (SGTF=1): P2CH3CQ ==0, P2CH2CQ  <= 30, P2CH1CQ <= 30.
# P2CH1CQ = ORF1ab
# P2CH2CQ = N gene
# P2CH3CQ = S gene
# P2CH4CQ = MS2 Control


# Build data set for modelling
# d: data set from complete_data()
# criterion: "under30CT" or "all" -- "under30CT" recommended
# death_cutoff: e.g. 28 for only considering deaths within 28 days of first positive test
# reg_cutoff: censor data at max_date - reg_cutoff
# P_voc: if between 0.5-1.0, classify as probable voc based upon modelled prevalence estimates; if 0, don't
# keep_missing: if TRUE, keep entries with missing sgtf information
model_data = function(d, criterion, remove_duplicates, death_cutoff, reg_cutoff, P_voc, 
    date_min = "2000-01-01", date_max = "2100-01-01", prevalence_cutoff = FALSE, sgtfv_cutoff = 0, keep_missing = FALSE)
{
    ct = function(x) ifelse(x == 0, 40, x)
    
    if (criterion == "under30CT") {
        sgtf_column = "sgtf_under30CT"
    } else if (criterion == "all") {
        sgtf_column = "sgtf"
    } else {
        stop("criterion must be under30CT or sgtf")
    }
    
    # Exclusion of duplicates
    dupes = numeric();
    if (remove_duplicates) {
        dupes = d[duplicated(FINALID), FINALID]
        dupes = setdiff(unique(dupes), NA);
    }
    
    data = d[!(FINALID %in% dupes) & # Exclude duplicates if requested
             !is.na(FINALID) & # Remove any deaths not linked to a Pillar 2 test
             !is.na(prefer(age.y, age.x)) & prefer(age.y, age.x) != 0 & # Seems like some age 0 individuals are miscoded unknowns. 
             !is.na(sex) & sex != "Unknown" & # Exclude unknown sex
             !is.na(LTLA_name) & LTLA_name != "" & # Exclude unknown LTLA
             !is.na(UTLA_name) & UTLA_name != "" & # Exclude unknown UTLA
             !is.na(NHSER_name) & NHSER_name != "" & # Exclude unknown NHS England region
             !is.na(specimen_date.x), # Exclude any unknown specimen dates
        .(sgtf = get(sgtf_column), p_voc = get(sgtf_column) * sgtfv, sgtfv = sgtfv, 
            age = prefer(age.y, age.x), sex = sex, 
            LTLA_name = factor(LTLA_name), UTLA_name = factor(UTLA_name), NHSER_name = factor(NHSER_name), 
            imd = imd_decile, 
            ethnicity_final = ethnicity_final.x, res = cat,
            specimen_date = as.Date(specimen_date.x), specimen_week = floor_date(as.Date(specimen_date.x), "1 week", week_start = 1),
            death_date = as.Date(dod), admission_date = as.Date(dateadmission_NHSE),
            # mv_pressure, ni_pressure, os_pressure, ao_pressure, medstaff_abs, nursing_abs,
            ctORF1ab = ct(P2CH1CQ), ctN = ct(P2CH2CQ), ctS = ct(P2CH3CQ), ctControl = ct(P2CH4CQ),
            asymptomatic = factor(asymptomatic_indicator),
            data_id)];
    
    if (!keep_missing) {
        data = data[!is.na(sgtf)]
    }
    
    # Set age and IMD groups
    data[, age_group := cut(age, c(1, 35, 55, 70, 85, 120), right = FALSE)]
    data[, imd_group := factor(paste0("imd", imd))]
    
    # Cutoff based upon LTLA prevalence of "false positives" from prior to Oct 15
    if (prevalence_cutoff)
    {
        prior_sgtf = data[specimen_date >= "2020-09-01" & specimen_date <= "2020-10-15", .(nhs_sgtf = mean(sgtf, na.rm = T)), by = NHSER_name]
        baseline_ltla = data[specimen_date >= "2020-09-01" & specimen_date <= "2020-10-15", .(Ns = sum(sgtf == 1, na.rm = T), No = sum(sgtf == 0, na.rm = T)), by = .(NHSER_name, LTLA_name)]
        baseline_ltla = merge(baseline_ltla, prior_sgtf, by = "NHSER_name")
        baseline_ltla = baseline_ltla[, .(priorS = sum(Ns) + mean(nhs_sgtf) * 100, priorO = sum(No) + mean(1 - nhs_sgtf) * 100), by = LTLA_name]
        baseline_ltla[, baseline := priorS / (priorS + priorO)]
        
        trace = data[specimen_date >= "2020-09-01", .(sgtf = mean(sgtf, na.rm = T), nspec = .N), keyby = .(LTLA_name, specimen_date)]
        trace[is.nan(sgtf), sgtf := 0]
        trace = merge(trace, baseline_ltla, by = "LTLA_name")
        trace[, indicator := sgtf > 1 - (1 - baseline)^2]
        trace[, n_good_so_far := cumsum(indicator * nspec), by = LTLA_name]
        trace[, n_bad_remaining := rev(cumsum(rev((!indicator) * nspec))), by = LTLA_name]
        co = trace[n_bad_remaining <= n_good_so_far, .(date_cutoff = min(specimen_date)), by = LTLA_name]
        
        data = merge(data, co, by = "LTLA_name")
        data = data[specimen_date >= date_cutoff]
        data[, date_cutoff := NULL]
    }
    
    # Restrict data based upon date and sgtfv
    data = data[specimen_date >= date_min & specimen_date <= date_max];
    if (!keep_missing) {
        data = data[sgtfv >= sgtfv_cutoff];
    }
    
    # Revalue certain factors
    data[, eth_cat := factor(revalue(ethnicity_final,
        c(
            "African (Black or Black British)"     = "B",
            "Any other Asian background"           = "A",
            "Any other Black background"           = "B",
            "Any other ethnic group"               = "O",
            "Any other Mixed background"           = "O",
            "Any other White background"           = "W",
            "Bangladeshi (Asian or Asian British)" = "A",
            "British (White)"                      = "W",
            "Caribbean (Black or Black British)"   = "B",
            "Chinese (other ethnic group)"         = "A",
            "Indian (Asian or Asian British)"      = "A",
            "Irish (White)"                        = "W",
            "Pakistani (Asian or Asian British)"   = "A",
            "Unknown"                              = "O",
            "White and Asian (Mixed)"              = "O",
            "White and Black African (Mixed)"      = "O",
            "White and Black Caribbean (Mixed)"    = "O"
        )), levels = c("W", "A", "B", "O"))];
    data[, ethnicity_final := factor(ethnicity_final)];
    
    data[res == "", res := " "]
    data[, res_cat := factor(revalue(res,
        c(
            " "                                                                        = "Other/Unknown",
            "Care/Nursing home"                                                        = "Care/Nursing home",
            "House in multiple occupancy (HMO)"                                        = "Residential",
            "Medical facilities (including hospitals and hospices, and mental health)" = "Other/Unknown",
            "No fixed abode"                                                           = "Other/Unknown",
            "Other property classifications"                                           = "Other/Unknown",
            "Overseas address"                                                         = "Other/Unknown",
            "Prisons, detention centres, secure units"                                 = "Other/Unknown",
            "Residential dwelling (including houses, flats, sheltered accommodation)"  = "Residential",
            "Residential institution (including residential education)"                = "Other/Unknown",
            "Undetermined"                                                             = "Other/Unknown"
        )), levels = c("Residential", "Care/Nursing home", "Other/Unknown"))];
    data[, res := factor(res)];
    
    # Remove entries with specimen date after death date
    data = data[is.na(death_date) | (death_date >= specimen_date)];
    
    # Probabilistic assessment of VOC
    if (P_voc %between% c(0.5, 1.0)) {
        data[, probable_voc := ifelse(p_voc >= P_voc, TRUE, ifelse((1 - p_voc) >= P_voc, FALSE, NA))];
        data = data[!is.na(probable_voc)];
    } else if (P_voc != 0) {
        stop("P_voc must be either 0 or between 0.5 and 1.0.");
    }

    # Censor data at death cutoff / registration cutoff
    max_date = data[, max(death_date, na.rm = T)] - reg_cutoff;
    data[, followup_date := pmin(death_date, specimen_date + death_cutoff, max_date, na.rm = T)];
    data = data[followup_date >= specimen_date];
    
    data[, died := !is.na(death_date) & death_date <= followup_date];
    data[, status := ifelse(died, 1, 0)];
    data[, time := as.numeric(followup_date - specimen_date)];
    
    # Add 0.5 to time = 0
    data[time == 0, time := 0.5]

    return (data[])
}

# Add more informative data categories
prep_data = function(data)
{
    data[, sgtf_label := ifelse(sgtf == 0, "Non-SGTF", "SGTF")]
    data[, sgtf_label := factor(sgtf_label, c("SGTF", "Non-SGTF"))]

    # Create more descriptive categories for table outputs
    data[, `Sex` := factor(sex)]
    data[, `Age` := factor(revalue(age_group,
        c(
            "[1,35)" = "1–34",
            "[35,55)" = "35–54",
            "[55,70)" = "55–69",
            "[70,85)" = "70–84",
            "[85,120)" = "85 and older"
        )), levels = c("1–34", "35–54", "55–69", "70–84", "85 and older"))]
    data[, `Place of residence` := factor(res_cat)]
    data[, `Index of Multiple Deprivation decile` := factor(imd, levels = 1:10)]
    data[, `Ethnicity` := factor(revalue(eth_cat,
        c(
            "W" = "White",
            "A" = "Asian",
            "B" = "Black",
            "O" = "Other/Mixed/Unknown"
        )), levels = c("White", "Asian", "Black", "Other/Mixed/Unknown"))]
    data[, `NHS England region` := factor(NHSER_name)]
    data[, spec_date_ind := as.numeric(specimen_date - ymd("2020-11-01")) %/% 14]
    # If last 2-week period contains 7 days or fewer, combine with penultimate period
    if (data[spec_date_ind == max(spec_date_ind), uniqueN(specimen_date) <= 7]) {
        data[spec_date_ind == max(spec_date_ind), spec_date_ind := spec_date_ind - 1];
    }
    data[, `Specimen date` := paste0(str_trim(format(min(specimen_date), "%e %b")), "–", str_trim(format(max(specimen_date), "%e %b"))), by = spec_date_ind]
    data[, `Specimen date` := factor(`Specimen date`, levels = data[order(spec_date_ind), unique(`Specimen date`)])]
    data[, ` ` := ""]
    
    return (data)
}
