% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correct.R
\name{correctTMTproImpurities}
\alias{correctTMTproImpurities}
\title{Performs purity correction on TMTpro intensity data.
The number of intensity columns must equal the number of reporters in the impurity file
Output is a data frame of the same dimensions}
\usage{
correctTMTproImpurities(data, impurities)
}
\arguments{
\item{data}{data frame of uncorrected TMT intensities}

\item{impurities}{data frame of impurities for a specific lot}
}
\description{
Performs purity correction on TMTpro intensity data.
The number of intensity columns must equal the number of reporters in the impurity file
Output is a data frame of the same dimensions
}
\examples{
\dontrun{
data <- read_tsv("evidence,txt") \%>\% select(matches("Reporter intensity \\\\d+"))
impurities <- read_csv("VJ309267.csv")
correctedIntensities <- correctTMTproImpurities(data, impurities)
}

}
