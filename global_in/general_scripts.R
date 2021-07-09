## Generic functions used throughout the project


#sampids <- stab$SAMPID

GTEx_SAMPID_to_SUBJID <- function(sampids){
    ## A simple function to quickly turn SAMPIDs into SUBJIDs to connect
    ## individuals to their phenotype
  stringr::str_extract(string = sampids, pattern = "GTEX-[[:alnum:]]*")
}
