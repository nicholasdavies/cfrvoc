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
