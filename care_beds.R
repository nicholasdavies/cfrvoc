# NOT USED

library(readxl)
library(data.table)
library(ogwrangler)
CreateCache()

beds = function(filename, reference)
{
    b = data.table(read_excel(filename, 2))
    b = b[`Care home?` == "Y", .(postcode = `Location Postal Code`, nbeds = `Care homes beds`)]
    b[, lsoa11 := ogpost(postcode, "lsoa11")]
    
    b = b[, .(nhomes = .N, nbeds = sum(nbeds)), keyby = lsoa11]

    lsoas = ogwrangler:::lookup[lsoa11 %like% "^E", unique(lsoa11)]
    bb = data.table(lsoa11 = lsoas)
    b = merge(b, bb, by = "lsoa11", all = TRUE)
    
    b[is.na(nhomes), nhomes := 0]
    b[is.na(nbeds), nbeds := 0]
    b[, pop2019 := ogwhat(lsoa11, "pop2019")]
    b[, ref := reference]

    return (b[])
}

beds_ref = rbind(
    beds("./care_homes/3 August 2020 HSCA active locations.xlsx", "2020-08-01"),
    beds("./care_homes/1 September 2020 HSCA Active Locations.xlsx", "2020-09-01"),
    beds("./care_homes/1 October 2020 HSCA Active Locations.xlsx", "2020-10-01"),
    beds("./care_homes/02 November 2020 HSCA Active Locations.xlsx", "2020-11-02"),
    beds("./care_homes/1 December 2020 HSCA Active Locations.xlsx", "2020-12-01"),
    beds("./care_homes/4_January_2021_HSCA_Active_Locations.xlsx", "2021-01-04")
)

im = fread("Index_of_Multiple_Deprivation_(December_2019)_Lookup_in_England.csv")

beds_ref = merge(beds_ref, im[, .(lsoa11 = LSOA11CD, imd_rank = IMD19)], all = TRUE)

fwrite(beds_ref, "./beds_ref.csv")
