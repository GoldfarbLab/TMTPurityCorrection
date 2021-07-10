## code to prepare `default_impurities` dataset goes here

default_impurities <- read_tsv("data-raw/VJ309267.csv")
usethis::use_data(default_impurities, overwrite = TRUE)
