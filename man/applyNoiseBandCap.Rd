% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correct.R
\name{applyNoiseBandCap}
\alias{applyNoiseBandCap}
\title{Performs noiseband capping on TMTpro intensity data.}
\usage{
applyNoiseBandCap(data, noise)
}
\arguments{
\item{data}{data frame of TMT intensities}

\item{noise}{data frame of TMT noisebands for all scans}
}
\description{
Performs noiseband capping on TMTpro intensity data.
}
\examples{
\dontrun{
data <- read_tsv("evidence,txt") \%>\% select(matches("Reporter intensity \\\\d+"))
impurities <- read_csv("VJ309267.csv")
correctedIntensities <- correctTMTproImpurities(data, impurities)
}

}
