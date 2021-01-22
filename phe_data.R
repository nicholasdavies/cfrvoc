library(data.table)
library(stringr)
library(lubridate)

# Assemble PHE data filename
phe_file = function(basename, dateid)
{
    paste0("~/Documents/uk_covid_data_sensitive/phe/", dateid, "/", str_replace_all(basename, "\\$", dateid))
}

# Load various PHE data files
phe_deaths = function(dateid)
{
    d = read_excel(phe_file("$ COVID19 Deaths.xlsx", dateid))
    setDT(d)
    return (d)
}

phe_negatives1 = function(dateid)
{
    fread(phe_file("$ Negatives pillar1.csv", dateid))
}

phe_negatives2 = function(dateid)
{
    fread(phe_file("$ Negatives pillar2.csv", dateid))
}

phe_positives = function(dateid)
{
    fread(phe_file("Anonymised Combined Line List $.csv", dateid))
}

phe_sgtf = function(dateid)
{
    fread(phe_file("SGTF_linelist_$.csv", dateid))
}

sgtf_counts = function(dateid, group_cols, criterion = "under30CT")
{
    pos = phe_positives(dateid)
    sgtf = phe_sgtf(dateid)
    data = merge(pos, sgtf, by = "FINALID", all = TRUE)
    data[, group := do.call(paste, c(.SD, sep = "|")), .SDcols = group_cols]
    
    if (criterion == "under30CT") {
        data[!is.na(sgtf_under30CT), .(other = sum(sgtf_under30CT == 0), sgtf = sum(sgtf_under30CT == 1)),
            keyby = .(specimen_date = dmy(specimen_date.x), group)]
    } else {
        data[!is.na(sgtf), .(other = sum(sgtf == 0), sgtf = sum(sgtf == 1)),
            keyby = .(specimen_date = dmy(specimen_date.x), group)]
    }
}


# # Assemble PHE data in different ways
# dateid = "20210118"
# 
# pos = phe_positives(dateid)
# sgtf = phe_sgtf(dateid)
# neg1 = phe_negatives1(dateid)
# neg2 = phe_negatives2(dateid)
# 
# # First, just looking at number of SGTF positive cases.
# ps = merge(pos, sgtf, by = "FINALID", all = TRUE)
# p = ps[!is.na(sgtf_under30CT), .(other = sum(sgtf_under30CT == 0), sgtf = sum(sgtf_under30CT == 1)), 
#     keyby = .(specimen_date = dmy(specimen_date.x), LTLA_code)]
# p = p[LTLA_code != ""]
# p[LTLA_code %in% c("E07000004", "E07000005", "E07000006", "E07000007"), LTLA_code := "E06000060"] # LADs amalgamated to Buckinghamshire
# 
# p[, rgn_code := ogwhat(LTLA_code, "rgn")]
# p[, geo := ogwhat(rgn_code)]
# p = p[, .(other = sum(other), sgtf = sum(sgtf)), keyby = .(specimen_date, geo)]
# 
# p
# 
# p0 = ps[!is.na(sgtf), .(other = sum(sgtf == 0), sgtf = sum(sgtf == 1)), 
#     keyby = .(specimen_date = dmy(specimen_date.x), LTLA_code)]
# p0 = p0[LTLA_code != ""]
# p0[LTLA_code %in% c("E07000004", "E07000005", "E07000006", "E07000007"), LTLA_code := "E06000060"] # LADs amalgamated to Buckinghamshire
# 
# p0[, rgn_code := ogwhat(LTLA_code, "rgn")]
# p0[, geo := ogwhat(rgn_code)]
# p0 = p0[, .(other = sum(other), sgtf = sum(sgtf)), keyby = .(specimen_date, geo)]
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # SCRATCH . . .
# 
# p[, table(specimen_date)]
# sgtf[, .(sum(sgtf == 0), sum(sgtf == 1))]
# sgtf[, .(sum(sgtf_under30CT == 0, na.rm = T), sum(sgtf_under30CT == 1, na.rm = T))]
# 
# w = phe_deaths("20210118")
# w = phe_negatives1("20210118")
# w = phe_negatives2("20210118")
# w = phe_positives("20210118")
# w = phe_sgtf("20210118")
# 
# pos = phe_positives("20210118")
# imd_ltla_pos = unique(pos[, .(imd_rank, LTLA_name)])
# 
# deaths = phe_deaths("20210118")
# pos
# 
# imd_ltla_pos
# 
# im = fread("Index_of_Multiple_Deprivation_(December_2019)_Lookup_in_England.csv")
# 
# im[which(duplicated(IMD19))]
# 
# w = merge(imd_ltla_pos, im, by.x = "imd_rank", by.y = "IMD19", all = TRUE)
# w[LTLA_name != LAD19NM]
# 
# # Assemble test number and test positivity data
# load_phe_set = function(basename, dateid)
# {
#     filename = paste0("~/Documents/uk_covid_data_sensitive/phe/", dateid, "/", str_replace_all(basename, "$", dateid))
# }
