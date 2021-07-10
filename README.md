# TMTPurityCorrection

## Installation
```R
library(devtools)
install_github("GoldfarbLab/TMTPurityCorrection")
```

## Usage
```R
# Extract uncorrected reporter intensities from MaxQuant results
data <- read_tsv("evidence,txt") %>% select(matches("Reporter intensity \\d+"))
# Read in lot-specific impurity file
impurities <- read_csv("VJ309267.csv")
# Correct impurities
correctedIntensities <- correctTMTproImpurities(data, impurities)
```
## Creating a lot-specific impurity file

1. Find the data sheet for your lot at the bottom of: https://www.thermofisher.com/order/catalog/product/A44520#/A44520
2. Create a csv file in the same format as data-raw/VJ309267.csv, which is created from the product sheet below.
3. Note, the "extra" column's value is the row sum of impurities that don't correspond to a label.

<img src="https://user-images.githubusercontent.com/8868522/125169892-ea191b00-e171-11eb-828b-01b754a103d1.png" width="600">

