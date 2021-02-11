# cfrvoc

Estimate of the hazard of death associated with infection by SARS-CoV-2 VOC 202012/01 relative to preexisting variants.

This repo includes an anonymised data set (in ./dataset/) which allows partial replication of the analyses in our preprint. 

Note that the anonymised data set randomises ages within 5-year age bands, and caps the maximum age at 104 for purposes of not 
publishing personally identifiable information. For this reason, results obtained using this data set will not be exactly 
equal to those found in the paper.

The FINALID column in the anonymised data set is also randomised.

# Guide to the code

The code is designed to operate either on non-public data supplied by PHE, or on the anonymised data set packaged with the repo.
Anonymised data sets live in the `./dataset/` folder.

`survival.R` is the main analysis file for the survival analyses. Calls to the R function `complete_data` can be replaced with 
calls to the R function `reduced_data` to load an anonymised data set. All output goes to the `./output/` folder which you will 
need to create in the project directory.
