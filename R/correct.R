#' Performs purity correction on TMTpro intensity data.
#' The number of intensity columns must equal the number of reporters in the impurity file
#' Output is a data frame of the same dimensions
#'
#' @param data data frame of uncorrected TMT intensities
#' @param impurities data frame of impurities for a specific lot
#'
#' @export
#' @examples
#' \dontrun{
#' data <- read_tsv("evidence,txt") %>% select(matches("Reporter intensity \\d+"))
#' impurities <- read_csv("VJ309267.csv")
#' correctedIntensities <- correctTMTproImpurities(data, impurities)
#' }
#'
#' @import magrittr
#' @import dplyr
#' @import matlib
#' @import stringr
#'
#'
correctTMTproImpurities <- function(data, impurities)
{
  nlabels <- nrow(impurities)

  impurities <- impurities %>%
    replace(is.na(.), 0) %>%
    select(-.data$label) %>%
    as.matrix()

  # normalize by row sum
  correction.matrix <- apply(impurities, 1, function(i) i/sum(i))[1:nlabels,]

  # compute inverse
  AI <- inv(correction.matrix)

  # compute correction per row
  corrected.intensities <- t(apply(data, 1, function(x) AI %*% as.matrix(x)))

  # set negatives to 0
  corrected.intensities[data <= 0] <- 0

  # set colnames
  colnames(corrected.intensities) <- str_c("Reporter.intensity.corrected.", 1:nlabels)

  return(corrected.intensities)
}
