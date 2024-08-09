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

"bamdata"

#' BAM Analysis of French Political Data
#'
#' The object `bam.france` is a result of applying the `BAM()` function to a dataset (`bamdata`) with specific parameters 
#' to analyze French political data. This object contains the Bayesian Aldrich-McKelvey scaling results.
#'
#' @format A list of class `BAM` containing the following components:
#' \describe{
#'   \item{polarity}{The polarity of the analysis, set to 2 in this case. This indicates the polarity constraint applied during the scaling.}
#'   \item{n.adapt}{The number of iterations used for adaptation, which was 2500 in this case.}
#'   \item{n.sample}{The number of MCMC samples collected, which was 5000 in this case.}
#'   \item{zhat}{A logical value indicating whether the ideal points should be adjusted for mean-zero scaling (`zhat=TRUE`).}
#'   \item{ab}{A logical value indicating whether the response ideal points are used (`ab=TRUE`).}
#'   \item{resp.idealpts}{A logical value indicating whether to estimate respondent ideal points (`resp.idealpts=TRUE`).}
#'   \item{data}{The original data used in the analysis.}
#'   \item{idealpoints}{The estimated ideal points from the BAM analysis.}
#'   \item{posteriors}{Posterior distributions of the ideal points and other parameters.}
#'   \item{convergence}{Convergence diagnostics for the MCMC chains.}
#'   \item{other_components}{Additional components that store the results and diagnostics of the BAM analysis.}
#' }
#'
#' @details
#' The `bam.france` object was created using the `BAM()` function with the following parameters:
#' \code{bam.france <- BAM(bamdata, polarity=2, n.adapt=2500, n.sample=5000, zhat=TRUE, ab=TRUE, resp.idealpts=TRUE)}
#' 
#' This analysis uses a Bayesian Aldrich-McKelvey scaling model to estimate ideal points for French political data, 
#' capturing the political preferences and scaling them accordingly.
#'
#'
#' @seealso \code{\link[BAM]{BAM}} for more details on the BAM function and its parameters.
#'
#' @source The BAM model was applied to a dataset `bamdata` with specific settings to generate `bam.france`.
#'
#' @keywords datasets
"bam.france"

