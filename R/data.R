#' Result of Aldmck Analysis for France EES 2009 Data
#'
#' This dataset contains the result of the `aldmck` analysis applied to the `franceEES2009` dataset.
#'
#' @details
#' The `result.france` object is generated using the `aldmck` function with the following parameters:
#' \code{result.france <- aldmck(franceEES2009, respondent = 1, polarity = 2, missing = c(77, 88, 89), verbose = FALSE)}
#'
#' @format An object of class \code{"aldmck"}.
#'
#' @usage data(result.france)
#'
#' @examples
#' data(result.france)
#' summary(result.france)
"result.france"


#' Selected Issues Sweden 2010 Dataset
#'
#' This dataset, `issues.sweden`, is a matrix created from the `Sweden2010` dataset, specifically using columns 7 to 56. It contains issue-related data from Sweden's 2010 election study.
#'
#' @format A matrix with rows representing respondents and columns representing different issues or variables.
#' @source Sweden 2010 Election Study
#' @examples
#' data(issues.sweden)
#' head(issues.sweden)
#' 
#' # Example usage in analysis
#' result <- aldmck(issues.sweden, respondent=1, polarity=2, missing=c(77,88,89), verbose=FALSE)
#' summary(result)
"issues.sweden"


#' Rankings Data from France EES 2009
#'
#' This dataset, `rankings`, contains the rankings data extracted from the France 2009 European Election Study (EES). It is a matrix of numeric rankings for various political parties.
#'
#' @format A numeric matrix with rows representing respondents and columns representing the rankings of different political parties.
#' @details
#' The `rankings` matrix was created from the `franceEES2009` dataset by selecting columns 2 to 9, which correspond to the respondents' rankings of different political parties. The data was then converted to numeric mode for further analysis.
#' 
#' The following code was used to create the `rankings` object:
#' 
#' \code{
#' rankings <- as.matrix(franceEES2009[,2:9])
#' mode(rankings) <- "numeric"
#' }
#'
#' This dataset is used as an input for the `blackbox_transpose` function to perform dimensional analysis.
#'
#' @source France 2009 European Election Study (EES)
#' @examples
#' \dontrun{
#' data(rankings)
#' original <- blackbox_transpose(rankings, missing = c(77, 88, 89), dims = 3, minscale = 5, verbose = FALSE)
#' }
"rankings"

#' BAM Data from France EES 2009
#'
#' This dataset, `bamdata`, was prepared using the `bamPrep` function on the France 2009 European Election Study (EES) data. It is used for Bayesian Aldrich-McKelvey (BAM) scaling.
#'
#' @format A list of class `bamPrep` with two components:
#' \describe{
#'   \item{stims}{A matrix of stimuli placements (excluding self-placement) with missing values handled.}
#'   \item{self}{A vector of self-placement values.}
#' }
#' @details
#' The `bamdata` object was created by applying the `bamPrep` function to the `franceEES2009` dataset with specified missing values (`c(77, 88, 89)`), focusing on self-placement (`self = 1`), and requiring a minimum of 5 non-missing values per respondent.
#'
#' @source France 2009 European Election Study (EES)
#' @examples
#' \dontrun{
#' data(bamdata)
#' result <- BAM(bamdata, polarity = 1, zhatSave = TRUE)
#' }
"bamdata"

