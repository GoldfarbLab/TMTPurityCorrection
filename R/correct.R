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
    select(-.data$label, -.data$mass) %>%
    as.matrix()

  # normalize by row sum
  correction.matrix <- apply(impurities, 1, function(i) i/sum(i))[1:nlabels,]

  # compute inverse
  AI <- inv(correction.matrix)

  # compute correction per row
  corrected.intensities <- t(apply(data, 1, function(x) AI %*% as.matrix(x)))

  # set negatives to 0
  corrected.intensities[data <= 0] <- 0
  corrected.intensities[corrected.intensities < 0] <- 0

  # set colnames
  colnames(corrected.intensities) <- str_c("Reporter.intensity.corrected.", 1:nlabels)

  return(corrected.intensities)
}








#' Performs noiseband capping on TMTpro intensity data.
#'
#' @param data data frame of TMT intensities
#' @param noise data frame of TMT noisebands for all scans
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
#'
#'
applyNoiseBandCap <- function(data, noise)
{
  # join data by scan number
  data.aligned <- data %>% inner_join(noise, by=c("Raw file", "Scan number"))
  data.left <- data %>% full_join(noise, by=c("Raw file", "Scan number"))
  # check that each quantified scan has matching noise
  if (nrow(data.aligned) != nrow(data))
  {
    warning("Data is missing matching noise values ", nrow(data), " ", nrow(data.aligned))
  }
  # check that the number of columns are the same
  if (ncol(data) != ncol(noise))
  {
    warning("Different number of columns between data and noise: ", ncol(data), ", ", ncol(noise))
  }

  num.quant.cols <- ncol(data)-2
  data.quant <- data.aligned %>% select(-`Raw file`, -`Scan number`)
  # split the aligned data
  data.intensity <- data.quant %>% select(1:num.quant.cols)
  data.noise <- data.quant %>% select((num.quant.cols+1):ncol(data.quant))

  data.intensity[data.intensity < data.noise] <- data.noise[data.intensity < data.noise]

  return(data.intensity)
}







#' Performs purity correction on TMT RTS data.
#' The number of intensity columns must equal the number of reporters in the impurity file
#' Output is a data frame of the same dimensions
#'
#' @param msms data frame of uncorrected TMT intensities, raw file, and MS2 scan number
#' @param impurities data frame of impurities for a specific lot
#' @param noise optional data frame of TMT noiseband values, raw file, and MS2 scan number
#' @param method purity correction method. Options = "NNLS" and "OLS"
#' @param noise.replacement.method should the noiseband replace missing values before or after correction. Options = "pre" and "post"
#' @param remain.missing should missing values be set to 0 (missing) after correction or not. Only used if not doing noiseband imputation. True or False
#' @param remove.missing.rows should rows containing all missing rows be removed or not. True or False
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
#'
#'
correctImpurities.RTS <- function(msms, impurities, noise,
                                  method=c("NNLS","OLS"),
                                  noise.replacement.method=c("pre","post"),
                                  remain.missing=T,
                                  remove.missing.rows=T,
                                  use.razor=T)
{
  return(correctedImpurities.msms(msms, impurities, noise, method=method,
                                  noise.replacement.method=noise.replacement.method,
                                  remaing.missing=remaining,
                                  remove.missing.rows=T))
}




#' Performs purity correction on MaxQuant data at specified levels
#' The number of intensity columns must equal the number of reporters in the impurity file
#' Output is a list of data frames for each level of analysis
#'
#' @param impurities data frame of impurities for a specific lot
#' @param txt.path path to the MaxQuant txt output folder. If using noiseband imputation, it should contain an intensity.txt and noise.txt file.
#' @param levels list of which levels to perform purity correction and noiseband imputation. Options are "msms", "evidence", and "proteinGroups"
#' @param method purity correction method. Options = "NNLS" and "OLS"
#' @param noise.replacement.method should the noiseband replace missing values before or after correction. Options = "pre", "post", "none"
#' @param remain.missing should missing values be set to 0 (missing) after correction or not. Only used if not doing noiseband imputation. True or False
#' @param remove.missing.rows should rows containing all missing rows be removed or not. True or False
#' @param use.razor should razor peptides be used for protein quantification. Only used for proteinGroups and if msms is supplied.
#' @param guess_max guess_max used for reading the csv files. It's the number of rows to read to guess column data types.
#'
#' @export
#' @examples
#' \dontrun{
#' impurities <- read_csv(""~/Box/CellBio-GoldfarbLab/Data/Annotations/TMT lots/A44521_VB294909.csv")
#' txt.path <- "~/Box/CellBio-GoldfarbLab/Data/Mass Spec/Search Results/Dennis/McGinty/txt/"
#' corrected.tables <- correctImpurities.MaxQuant(impurities, txt.path)
#' }
#'
#'
correctImpurities.MaxQuant <- function(impurities, txt.path,
                                       levels=c("msms", "evidence", "proteinGroups"),
                                       method=c("NNLS","LS"),
                                       noise.replacement.method=c("pre","post","none"),
                                       remain.missing=T,
                                       remove.missing.rows=T,
                                       use.razor=T,
                                       guess_max=20000)
{
  # match up arguments to valid options
  levels <- match.arg(levels, several.ok = T)
  method <- match.arg(method)
  noise.replacement.method <- match.arg(noise.replacement.method)
  msms <- NA
  intensity <- NA
  noise <- NA

  # read in necessary files
  if ("msms" %in% levels || noise.replacement.method != "none")
  {
    msms <- read_tsv(str_c(txt.path, "/msms.txt"), guess_max=guess_max)
  }
  if ("evidence" %in% levels)
  {
    evidence <- read_tsv(str_c(txt.path, "/evidence.txt"), guess_max=guess_max)
  }
  if ("proteinGroups" %in% levels)
  {
    proteinGroups <- read_tsv(str_c(txt.path, "/proteinGroups.txt"), guess_max=guess_max)
    summary <- read_tsv(str_c(txt.path, "/summary.txt"))
  }
  if (noise.replacement.method != "none")
  {
    noise <- read_tsv(str_c(txt.path, "/noise.txt"))
    intensity <- read_tsv(str_c(txt.path, "/intensity.txt"))
  }

  output <- list()

  # call requested functions
  if ("msms" %in% levels)
  {
    output$msms <- correctImpurities.msms(msms, impurities,
                                          intensities=intensity,
                                          noise=noise,
                                          method=method,
                                          noise.replacement.method=noise.replacement.method,
                                          remain.missing=remain.missing,
                                          remove.missing.rows=remove.missing.rows)
  }
  if ("evidence" %in% levels)
  {
    output$evidence <- correctImpurities.evidence(evidence, impurities,
                                                  msms=msms,
                                                  intensities=intensity,
                                                  noise=noise,
                                                  method=method,
                                                  noise.replacement.method=noise.replacement.method,
                                                  remain.missing=remain.missing,
                                                  remove.missing.rows=remove.missing.rows)
  }
  if ("proteinGroups" %in% levels)
  {
    output$proteinGroups <- correctImpurities.proteinGroups(summary,
                                                            proteinGroups, impurities,
                                                            msms=msms,
                                                            intensities=intensity,
                                                            noise=noise,
                                                            method=method,
                                                            use.razor=use.razor,
                                                            noise.replacement.method=noise.replacement.method,
                                                            remain.missing=remain.missing,
                                                            remove.missing.rows=remove.missing.rows)
  }

  # return list
  return(output)
}












#' Performs purity correction on MaxQuant msms.txt data.
#'
#' @param msms data frame of msms.txt
#' @param impurities data frame of impurities for a specific lot
#' @param intensities data frame of TMT intensity values, raw file, and MS2 scan number
#' @param noise optional data frame of TMT noiseband values, raw file, and MS2 scan number
#' @param method purity correction method. Options = "NNLS" and "OLS"
#' @param noise.replacement.method should the noiseband replace missing values before or after correction. Options = "pre", "post", "none"
#' @param remain.missing should missing values be set to 0 (missing) after correction or not. Only used if not doing noiseband imputation. True or False
#' @param remove.missing.rows should rows containing all missing rows be removed or not. True or False
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
#'
#'
correctImpurities.msms <- function(msms, impurities, intensities, noise,
                                   method=c("NNLS","OLS"),
                                   noise.replacement.method=c("pre","post","none"),
                                   remain.missing=T,
                                   remove.missing.rows=T)
{
  method <- match.arg(method)
  noise.replacement.method <- match.arg(noise.replacement.method)


  # join data by raw file and scan number
  aligned <- msms %>%
    select(-matches("Reporter intensity \\d+")) %>%
    left_join(intensities, by=c("Raw file", "Scan number"))

  if (!missing(noise) && !is.na(noise))
  {
    aligned <- aligned %>%
      left_join(noise, by=c("Raw file", "Scan number"))
  }

  # exclude scans with no reporter ions since we don't want to impute them all
  if (remove.missing.rows)
  {
    aligned <- aligned %>% filter(rowSums(across(matches("Reporter intensity \\d+"))) > 0)
    msms <- msms %>% filter(id %in% aligned$id)
  }

  # extract intensities
  intensities <- aligned %>% select(matches("Reporter intensity \\d+"))
  # extract noise
  if (!missing(noise) && !is.na(noise))
  {
    noise <- aligned %>% select(matches("Reporter noise \\d+"))

    # check that the number of columns are the same
    if (ncol(intensities) != ncol(noise))
    {
      stop("Different number of columns between data and noise: ", ncol(intensities), ", ", ncol(noise))
    }
  }

  # correct impurities
  corrected.intensities <- correctImpurities(intensities, impurities, noise, method=method,
                                noise.replacement.method=noise.replacement.method,
                                remain.missing=remain.missing)

  # overwrite corrected columns
  corrected.msms <- msms %>%
    select(-matches("Reporter intensity corrected \\d+")) %>%
    cbind(corrected.intensities)

  return(corrected.msms)
}






#' Performs purity correction on MaxQuant evidence.txt data
#'
#' @param evidence data frame of evidence.txt
#' @param msms optional data frame of msms.txt to perform correction at spectrum level. Required for noiseband imputation.
#' @param impurities data frame of impurities for a specific lot
#' @param intensities data frame of TMT intensity values, raw file, and MS2 scan number
#' @param noise optional data frame of TMT noiseband values, raw file, and MS2 scan number
#' @param method purity correction method. Options = "NNLS" and "OLS"
#' @param noise.replacement.method should the noiseband replace missing values before or after correction. Options = "pre", "post", "none"
#' @param remain.missing should missing values be set to 0 (missing) after correction or not. Only used if not doing noiseband imputation. True or False
#' @param remove.missing.rows should rows containing all missing rows be removed or not. True or False
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
#'
#'
correctImpurities.evidence <- function(evidence, impurities, msms, intensities, noise,
                                       method=c("NNLS","OLS"),
                                       noise.replacement.method=c("pre","post","none"),
                                       remain.missing=T,
                                       remove.missing.rows=T)
{
  method <- match.arg(method)
  noise.replacement.method <- match.arg(noise.replacement.method)

  # correct evidence without noise
  if (missing(msms) || is.na(msms))
  {
    intensities <- evidence %>% select(matches("Reporter intensity \\d+"))
    corrected.intensities <- correctImpurities(intensities, impurities,
                                               method=method,
                                               noise.replacement.method=noise.replacement.method,
                                               remain.missing=remain.missing)

    corrected.evidence <- evidence %>%
      select(-matches("Reporter intensity corrected \\d+")) %>%
      cbind(corrected.intensities)
  }

  # correct evidence using msms and possibly noise
  else
  {
    # correct MSMS
    corrected.msms <- correctImpurities.msms(msms, impurities,
                                             intensities=intensities,
                                             noise=noise,
                                             method=method,
                                             noise.replacement.method=noise.replacement.method,
                                             remain.missing=remaing.missing)

    # summarize to evidence level
    summarized.msms <- corrected.msms %>%
      group_by(`Evidence ID`) %>%
      summarise(across(matches("Reporter intensity corrected \\d+"), ~sum(., na.rm=T)),
                across(matches("Is missing \\d+"), ~all(as.logical(.), na.rm=T)),
                across(matches("SN \\d+"), ~sum(., na.rm=T)))

    # overwrite corrected intensity values
    corrected.evidence <- evidence %>%
      select(-matches("Reporter intensity corrected \\d+")) %>%
      inner_join(summarized.msms, by=c("id" = "Evidence ID"))
  }

  # exclude scans with no reporter ions since we don't want to impute them all
  if (remove.missing.rows)
  {
    corrected.evidence <- corrected.evidence %>% filter(rowSums(across(matches("Reporter intensity \\d+"))) > 0)
  }

  return(corrected.evidence)
}







#' Performs purity correction on MaxQuant proteinGroups.txt data
#'
#' @param summary data frame of summary.txt
#' @param proteinGroups data frame of evidence.txt
#' @param msms optional data frame of msms.txt to perform correction at spectrum level. Required for noiseband imputation.
#' @param impurities data frame of impurities for a specific lot
#' @param intensities data frame of TMT intensity values, raw file, and MS2 scan number
#' @param noise optional data frame of TMT noiseband values, raw file, and MS2 scan number
#' @param use.razor should razor peptides be used for protein quantification. Only used if msms is supplied.
#' @param method purity correction method. Options = "NNLS" and "OLS"
#' @param noise.replacement.method should the noiseband replace missing values before or after correction. Options = "pre", "post", and "none"
#' @param remain.missing should missing values be set to 0 (missing) after correction or not. Only used if not doing noiseband imputation. True or False
#' @param remove.missing.rows should rows containing all missing rows be removed or not. True or False
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
#'
#'
correctImpurities.proteinGroups <- function(summary, proteinGroups, impurities, msms, intensities, noise,
                                            use.razor=T,
                                            method=c("NNLS","OLS"),
                                            noise.replacement.method=c("pre","post", "none"),
                                            remain.missing=T,
                                            remove.missing.rows=T)
{
  method <- match.arg(method)
  noise.replacement.method <- match.arg(noise.replacement.method)

  # get raw file to experiment mapping
  summary <- summary %>%
    filter(`Raw file` != "Total") %>%
    select(`Raw file`, Experiment)
  experiments <- unique(summary$Experiment)

  corrected.proteinGroups <- proteinGroups %>%
    select(-matches("Reporter intensity corrected \\d+"))

  # correct proteinGroups without noise
  if (missing(msms) || is.na(msms))
  {
    for (exp in experiments)
    {
      if (exp == "")
      {
        intensities <- proteinGroups %>% select(matches("Reporter intensity \\d+$"))
      }
      else
      {
        intensities <- proteinGroups %>% select(matches(str_c("Reporter intensity \\d+ ", exp)))
      }

      corrected.intensities <- correctImpurities(intensities, impurities,
                                                 method=method,
                                                 noise.replacement.method=noise.replacement.method,
                                                 remain.missing=remain.missing)

      if (exp != "")
      {
        corrected.intensities <- corrected.intensities %>% as.data.frame() %>%
          rename_with(.cols = matches("Reporter intensity corrected \\d+"), .fn = ~ str_c(.x, exp, sep=" "))
      }

      corrected.proteinGroups <- corrected.proteinGroups %>% cbind(corrected.intensities)

    }
  }


  # correct proteinGroups using msms and possibly noise
  else
  {
    # correct MSMS
    corrected.msms <- correctImpurities.msms(msms, impurities,
                                             intensities=intensities,
                                             noise=noise,
                                             method=method,
                                             noise.replacement.method=noise.replacement.method,
                                             remain.missing=remaing.missing)

    # get razor peptide IDs and their protein ID from proteinGroups
    razor.peptide.ids <- proteinGroups %>%
      select(`Peptide is razor`, `Peptide IDs`, `id`) %>%
      separate_rows(everything(), sep=";", convert=T) %>%
      filter(`Peptide is razor` == "True") %>%
      select(`Peptide IDs`, `id`)

    if (!use.razor)
    {
      corrected.msms <- corrected.msms %>% filter(!str_detect(`Protein group IDs`, ";"))
    }

    for (exp in experiments)
    {
      # filter corrected msms for raw files with this experiment
      if (exp == "")
      {
        corrected.msms.exp <- corrected.msms
      }
      else
      {
        corrected.msms.exp <- corrected.msms %>%
          inner_join(summary) %>%
          filter(Experiment == exp)
      }



      # summarize to proteinGroup level
      summarized.msms.exp <- corrected.msms.exp %>%
        separate_rows(`Protein group IDs`, sep=";") %>%
        mutate(`Protein group IDs` = as.numeric(`Protein group IDs`)) %>%
        inner_join(razor.peptide.ids, by=c("Peptide ID" = "Peptide IDs", "Protein group IDs" = "id")) %>%
        group_by(`Protein group IDs`) %>%
        summarise(across(matches("Reporter intensity corrected \\d+"), ~sum(., na.rm=T)),
                  across(matches("Is missing \\d+"), ~all(as.logical(.), na.rm=T)),
                  across(matches("SN \\d+"), ~sum(., na.rm=T)))

      # update columns if multiple experiments
      if (exp != "" & length(experiments) > 1)
      {
        summarized.msms.exp <- summarized.msms.exp %>%
          rename_with(.cols = matches("Reporter intensity corrected \\d+"), .fn = ~ str_c(.x, exp, sep=" ")) %>%
          rename_with(.cols = matches("Is missing \\d+"), .fn = ~ str_c(.x, exp, sep=" ")) %>%
          rename_with(.cols = matches("SN \\d+"), .fn = ~ str_c(.x, exp, sep=" "))
      }

      # overwrite corrected intensity values.
      corrected.proteinGroups <- corrected.proteinGroups %>%
        left_join(summarized.msms.exp, by=c("id" = "Protein group IDs"))
    }
  }

  # exclude scans with no reporter ions since we don't want to impute them all
  if (remove.missing.rows)
  {
    corrected.proteinGroups <- corrected.proteinGroups %>% filter(rowSums(across(matches("Reporter intensity \\d+"))) > 0)
  }

  return(corrected.proteinGroups)
}















#' Performs purity correction on generic TMT data.
#'
#' @param data data frame of uncorrected TMT intensities.
#' @param impurities data frame of impurities for a specific lot
#' @param noise optional data frame of TMT noiseband values. Must line up row by row with data
#' @param method purity correction method. Options = "NNLS" and "OLS"
#' @param noise.replacement.method should the noiseband replace missing values before or after correction. Options = "pre", "post", "none"
#' @param remain.missing should missing values be set to 0 (missing) after correction or not. Only used if not doing noiseband imputation. True or False
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
#'
#'
correctImpurities <- function(data, impurities, noise,
                              method=c("NNLS","OLS"),
                              noise.replacement.method=c("pre","post","none"),
                              remain.missing=T)
{
  method <- match.arg(method)
  noise.replacement.method <- match.arg(noise.replacement.method)

  ##########################
  # Parameter sanity checks
  ##########################
  if (!missing(noise) && !is.na(noise) && dim(noise) != dim(data))
  {
    stop("Data and noise are different dimensions")
  }

  # Reformat impurities
  impurities <- impurities %>%
    replace(is.na(.), 0) %>%
    select(-.data$label, -.data$mass) %>%
    as.matrix()

  data <- as.matrix(data)

  # Generate matrix of isMissing
  isMissing <- data == 0
  colnames(isMissing) <- str_replace(colnames(isMissing), "Reporter intensity", "Is missing")

  # Do not perform noise band correction
  if (missing(noise) || is.na(noise))
  {
    corrected.data <- .correct(data, impurities, method)
    if (remain.missing)
    {
      # anything that was missing should still be missing
      corrected.data[data <= 0 | is.na(data)] <- 0
    }
  }

  # Perform noise band correction
  else if (!missing(noise) && !is.na(noise))
  {
    noise <- as.matrix(noise)
    # Generate matrix of S/N
    sn <- data / noise
    colnames(sn) <- str_replace(colnames(sn), "Reporter intensity", "SN")

    # Add noise before purity correction
    if (noise.replacement.method=="pre")
    {
      missing <- is.na(data) | data==0
      data[missing] <- noise[missing]
    }

    # perform purity correction
    corrected.data <- .correct(data, impurities, method)

    # Add noise after purity correction
    if (noise.replacement.method=="post")
    {
      missing.or.low <- is.na(data) | data==0 | corrected.data < noise
      corrected.data[missing.or.low] <- noise[missing.or.low]
    }

    corrected.data <- corrected.data %>% cbind(sn)
  }

  corrected.data <- corrected.data %>% cbind(isMissing)

  return(corrected.data)

}


###
.correct <- function(data, impurities, method=c("NNLS", "OLS"))
{
  method <- match.arg(method)

  if (method == "NNLS") return(.NNLS.correction(data, impurities))
  else if (method == "OLS") return(.OLS.correction(data, impurities))
  else return(.OLS.correction(data, impurities))
}

###
.NNLS.correction <- function(data, impurities)
{
  nlabels <- nrow(impurities)

  # normalize by row sum
  correction.matrix <- apply(impurities, 1, function(i) i/sum(i))[1:nlabels,]

  # perform NNLS on each row
  corrected.intensities <- t(apply(data, 1, function(row) {nnls::nnls(correction.matrix, row)$x}))

  # set colnames
  colnames(corrected.intensities) <- str_c("Reporter intensity corrected ", 1:nlabels)

  return(corrected.intensities)
}

###
.OLS.correction <- function(data, impurities)
{
  nlabels <- nrow(impurities)

  # normalize by row sum
  correction.matrix <- apply(impurities, 1, function(i) i/sum(i))[1:nlabels,]

  # compute inverse
  AI <- inv(correction.matrix)

  # compute correction per row
  corrected.intensities <- t(apply(data, 1, function(x) AI %*% as.matrix(x)))

  # set negatives to 0
  corrected.intensities[corrected.intensities < 0] <- 0

  # set colnames
  colnames(corrected.intensities) <- str_c("Reporter intensity corrected ", 1:nlabels)

  return(corrected.intensities)
}


