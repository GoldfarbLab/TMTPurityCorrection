#' Prepares the default_impurities dataset
default_impurities <- read_csv("data-raw/VJ309267.csv")
usethis::use_data(default_impurities, overwrite = TRUE)
