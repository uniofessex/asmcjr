#' @encoding UTF-8
#' @title Result of Aldmck Analysis for France EES 2009 Data
#' @description This dataset contains the result of the `aldmck` analysis applied to the `franceEES2009` dataset.
#'
#' @details
#' The `result.france` object is generated using the `aldmck` function with the following parameters:
#' \code{result.france <- aldmck(franceEES2009, respondent = 1, polarity = 2, missing = c(77, 88, 89), verbose = FALSE)}
#'
#' @format An object of class \code{"aldmck"}.
#'
#' @usage data(result.france)
#' @keywords datasets
#' @name result.france
#' @docType data
NULL

#' @encoding UTF-8
#' @title Selected Issues Sweden 2010 Dataset
#' @description This dataset, `issues.sweden`, is a matrix created from the `Sweden2010` dataset, specifically using columns 7 to 56. It contains issue-related data from Sweden's 2010 election study.
#'
#' @format A matrix with rows representing respondents and columns representing different issues or variables.
#' @source Sweden 2010 Election Study
#' @name issues.sweden
#' @keywords datasets
#' @docType data
NULL


#' @encoding UTF-8
#' @title Rankings Data from France EES 2009
#' @description This dataset, `rankings`, contains the rankings data extracted from the France 2009 European Election Study (EES). It is a matrix of numeric rankings for various political parties.
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
#' @keywords datasets
#' @source France 2009 European Election Study (EES)
#' @name rankings
#' @docType data
NULL

#' @encoding UTF-8
#' @title BAM Data from France EES 2009
#'
#' @description This dataset, `bamdata`, was prepared using the `bamPrep` function on the France 2009 European Election Study (EES) data. It is used for Bayesian Aldrich-McKelvey (BAM) scaling.
#'
#' @format A list of class `bamPrep` with two components:
#' \describe{
#'   \item{stims}{A matrix of stimuli placements (excluding self-placement) with missing values handled.}
#'   \item{self}{A vector of self-placement values.}
#' }
#' @details
#' The `bamdata` object was created by applying the `bamPrep` function to the `franceEES2009` dataset with specified missing values (`c(77, 88, 89)`), focusing on self-placement (`self = 1`), and requiring a minimum of 5 non-missing values per respondent.
#' @keywords datasets
#' @source France 2009 European Election Study (EES)
#' @name bamdata
#' @docType data
NULL


#' @encoding UTF-8
#' @title BAM Analysis of French Political Data
#'
#' @description The object `bam.france` is a result of applying the `BAM()` function to a dataset `bamdata`) with specific parameters 
#' to analyze French political data. This object contains the Bayesian Aldrich-McKelvey scaling results.
#'
#' @format A list of class `BAM` containing the following components:
#' \describe{
#'   \item{polarity}{The polarity of the analysis, set to 2 in this case. This indicates the polarity constraint applied during the scaling.}
#'   \item{n.adapt}{The number of iterations used for adaptation, which was 2500 in this case.}
#'   \item{n.sample}{The number of MCMC samples collected, which was 5000 in this case.}
#'   \item{zhat}{A logical value indicating whether the ideal points should be adjusted for mean-zero scaling (zhat=TRUE).}
#'   \item{ab}{A logical value indicating whether the response ideal points are used (ab=TRUE).}
#'   \item{resp.idealpts}{A logical value indicating whether to estimate respondent ideal points (resp.idealpts=TRUE).}
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
#' @seealso \code{\link{BAM}} for more details on the BAM function and its parameters.
#'
#' @source The BAM model was applied to a dataset \code{bamdata} with specific settings to generate \code{bam.france}.
#'
#' @keywords datasets
#' @name bam.france
#' @docType data
NULL

#' @encoding UTF-8
#' @title Issues Matrix from CDS2000 Dataset
#' @description This object, `issues`, is a matrix extracted from the `CDS2000` dataset. It contains selected columns that represent various issues.
#' @format A numeric matrix with rows corresponding to observations and columns representing different issues.
#' @source Extracted from the `CDS2000` dataset.
#' @examples
#' \dontrun{
#' data(issues)
#' }
#' @keywords datasets
#' @name issues
#' @docType data
NULL


#' @encoding UTF-8
#' @title Blackbox Analysis Results for Issues Matrix
#' @description The object `result.repdem` contains the results of applying the `blackbox` function to the `issues` matrix. This analysis was performed to extract dimensions that represent the underlying structure of the issues, specifically for a dataset containing political data.
#' @format An object of class `blackbox` containing the results of the dimensional analysis, including the identified dimensions and other related statistics.
#'
#' @return The `result.repdem` object contains the extracted dimensions, along with other statistics generated by the `blackbox` function.
#' 
#' @examples
#' \dontrun{
#' data(result.repdem)
#' }
#' @keywords datasets
#' @name result.repdem
#' @docType data
NULL

#' @encoding UTF-8
#' @title Prepare Input Matrix from interest1981 Dataset
#' @description The `input2` object is a matrix created from the `interest1981` dataset for further analysis. The process involves extracting relevant columns, filtering rows, transforming values, and handling missing data.
#'
#' @details
#' The `input2` matrix is derived from the `interest1981` dataset through the following steps:
#'
#' 1. **Extract relevant columns**: 
#'    \code{input <- as.matrix(interest1981[, 9:38])}
#'
#' 2. **Filter rows with sufficient data**: 
#'    \code{input <- input[rowSums(!is.na(input)) >= 5, ]}
#'
#' 3. **Transform the matrix**: 
#'    \code{input <- (100 - input) / 50}
#'
#' @format A numeric matrix with rows corresponding to filtered observations and columns representing transformed variables.
#'
#' @examples
#' \dontrun{
#' data(input)
#' }
#'
#' @seealso \code{\link[base]{as.matrix}}, \code{\link[base]{rowSums}}, \code{\link[base]{mean}}
#' @name input
#' @docType data
#' @keywords datasets
NULL


#' @encoding UTF-8
#' @title Prepare Input Matrix from interest1981 Dataset
#' @description The `input2` object is a matrix created from the `interest1981` dataset for further analysis. The process involves extracting relevant columns, filtering rows, transforming values, and handling missing data.
#'
#' @details
#' The `input2` matrix is derived from the `interest1981` dataset through the following steps:
#'
#' 1. **Extract relevant columns**: 
#'    \code{input <- as.matrix(interest1981[, 9:38])}
#'
#' 2. **Filter rows with sufficient data**: 
#'    \code{input <- input[rowSums(!is.na(input)) >= 5, ]}
#'
#' 3. **Transform the matrix**: 
#'    \code{input <- (100 - input) / 50}
#'
#' 4. **Square the matrix to create `input2`**: 
#'    \code{input2 <- input * input}
#'
#' 5. **Handle missing values**: 
#'    \code{input2[is.na(input)] <- (mean(input, na.rm = TRUE))^2}
#'
#' @format A numeric matrix with rows corresponding to filtered observations and columns representing transformed variables.
#'
#' @examples
#' \dontrun{
#' data(input2)
#' input <- as.matrix(interest1981[, 9:38])
#' }
#'
#' @seealso \code{\link[base]{as.matrix}}, \code{\link[base]{rowSums}}, \code{\link[base]{mean}}
#' @name input2
#' @docType data
#' @keywords datasets
NULL


#' @encoding UTF-8
#' @title MDS Solution from mlsmu6
#' @description The `mlsmu6_out` dataset contains the multidimensional scaling (MDS) solution generated by applying the `mlsmu6` function to a subset of the `interest1981` dataset. The MDS solution is computed using two dimensions and a cutoff value of 5, with data grouped by political party labels.
#'
#' @format A matrix or data frame, depending on the structure of the output from `mlsmu6`, typically containing coordinates in the reduced dimensional space.
#'
#' @details
#' The `mlsmu6_out` object is generated by running the following command:
#' \preformatted{
#' mlsmu6_out <- mlsmu6(input = interest1981[, 9:38], ndim = 2, cutoff = 5,
#'                      id = factor(interest1981$party, labels = c("D", "R")))
#' }
#' The `mlsmu6` function applies multidimensional scaling to reduce the dimensionality of the input data, with the number of dimensions set to 2 and a cutoff value of 5. The input data is a subset of the `interest1981` dataset, specifically columns 9 to 38, grouped by the political party (`party`) variable.
#'
#' @examples
#' \dontrun{
#' data(mlsmu6_out)
#' }
#' @seealso \code{\link{mlsmu6}}, \code{\link{interest1981}}
#' @source Generated using the `mlsmu6` function on the `interest1981` dataset.
#' @name mlsmu6_out
#' @docType data
#' @keywords datasets
NULL


#' @encoding UTF-8
#' @title ANES Input Data
#' @description The `anes.input` object is a subset of the American National Election Study (ANES) 1968 dataset. It contains selected variables used as input for further analysis.
#'
#' @details
#' The `anes.input` object is created by loading the `ANES1968` dataset from the `asmcjr` package and selecting the first 12 columns. The data is then converted to a matrix format for analysis.
#'
#' The steps to create `anes.input` are as follows:
#' \preformatted{
#' data(ANES1968)
#' anes.input <- ANES1968[, 1:12]
#' anes.input <- as.matrix(anes.input)
#' }
#'
#' @format A numeric matrix with 12 columns, representing selected variables from the ANES 1968 dataset.
#'
#' @source The data comes from the American National Election Study (ANES) 1968, available in the `asmcjr` package.
#'
#' @examples
#' \dontrun{
#' data(anes.input)
#' }
#' @name anes.input
#' @docType data
#' @keywords datasets
NULL


#' @encoding UTF-8
#' @title ANES 1968 Feeling Thermometers and Voting Data
#' @description The `ANES1968` dataset includes feeling thermometers data and information on whether respondents reported voting in the 1968 elections, as well as their reported presidential vote choice.
#'
#' @format A dataframe with several variables:
#' \describe{
#'   \item{vote.turnout}{Binary variable indicating whether the respondent reported voting in the 1968 election (1 = voted, 0 = did not vote).}
#'   \item{presidential.vote}{Categorical variable representing the respondent's reported presidential vote choice (1 = Humphrey, 2 = Nixon, 3 = Wallace).}
#'   \item{T}{A matrix of feeling thermometers ranging from 0 ("very cold or unfavorable") to 100 ("very warm or favorable").}
#' }
#'
#' @details
#' The feeling thermometers measure respondents' attitudes towards various political figures and entities, with higher values indicating warmer or more favorable feelings. Values of 98, 99, and 100 were recoded in the original ANES file as 97, allowing these values to represent refused or missing responses.
#'
#' The `ANES1968` dataset is a valuable resource for analyzing voter behavior and sentiment during the 1968 U.S. presidential election.
#'
#' @source American National Election Studies (ANES) 1968 dataset.
#'
#' @examples
#' \dontrun{
#' data(ANES1968)
#' summary(ANES1968)
#' }
#' @name ANES1968
#' @keywords datasets
#' @docType data
NULL


#' @encoding UTF-8
#' @title ANES 2004 Issue Scales Data
#' @description The `ANES2004` dataset includes responses to various issue scales from the 2004 American National Election Studies (ANES). The dataset contains respondents' positions on several key political issues, measured on scales with varying ranges.
#'
#' @format A dataframe with several variables representing different political issues:
#' \describe{
#'   \item{libcon}{Liberal–Conservative scale, ranging from 1 (left) to 7 (right).}
#'   \item{diplomacy}{Diplomacy–Military Force scale, ranging from 1 (left) to 7 (right).}
#'   \item{iraqwar}{Bush’s Handling of Iraq War scale, ranging from 1 (right) to 4 (left).}
#'   \item{govtspend}{Government Spending/Services scale, ranging from 1 (right) to 7 (left).}
#'   \item{defense}{Defense Spending scale, ranging from 1 (left) to 7 (right).}
#'   \item{bushtaxcuts}{Bush Tax Cuts scale, ranging from 1 (right) to 4 (left).}
#'   \item{healthinsurance}{Government Health Insurance scale, ranging from 1 (left) to 7 (right).}
#'   \item{govtjobs}{Guaranteed Jobs scale, ranging from 1 (left) to 7 (right).}
#'   \item{aidblacks}{Government Aid to Blacks scale, ranging from 1 (left) to 7 (right).}
#'   \item{govtfundsabortion}{Government Abortion Funding scale, ranging from 1 (left) to 4 (right).}
#'   \item{partialbirthabortion}{Partial-Birth Abortion Ban scale, ranging from 1 (right) to 4 (left).}
#'   \item{environmentjobs}{Environment–Jobs scale, ranging from 1 (left) to 7 (right).}
#'   \item{deathpenalty}{Death Penalty scale, ranging from 1 (right) to 4 (left).}
#'   \item{gunregulations}{Gun Regulations scale, ranging from 1 (left) to 5 (right).}
#'   \item{womenrole}{Women’s Role scale, ranging from 1 (left) to 7 (right).}
#'   \item{gaymarriage}{Gay Marriage scale, ranging from 1 (left) to 3 (right).}
#' }
#'
#' @details
#' The `ANES2004` dataset captures respondents' views on a variety of political and social issues during the 2004 U.S. election period. The scales vary in their directionality, with some scales placing liberal positions on the left and others placing conservative positions on the right.
#'
#' The dataset is useful for analyzing public opinion on key issues and understanding the political landscape during the 2004 election.
#'
#' @source American National Election Studies (ANES) 2004 dataset.
#'
#' @examples
#' \dontrun{
#' data(ANES2004)
#' summary(ANES2004)
#' }
#' @name ANES2004
#' @docType data
#' @keywords datasets
NULL


#' @encoding UTF-8
#' @title Danish Module of the 2009 European Election Study (EES)
#' @description The \code{denmarkEES2009} dataset contains data from the Danish module of the 2009 European Election Study (EES). 
#' This dataset includes responses from 1,000 Danish participants who rated their propensity to vote for each of 
#' eight political parties on a 0–10 point scale. A score of 0 denotes "not at all possible," while a score of 10 denotes "very probable."
#' 
#' The dataset can be used to demonstrate the \code{smacofRect()} function.
#' 
#' @format A matrix with 1,000 rows and 8 columns:
#' \describe{
#'   \item{V1}{The propensity to vote for Party 1.}
#'   \item{V2}{The propensity to vote for Party 2.}
#'   \item{V3}{The propensity to vote for Party 3.}
#'   \item{V4}{The propensity to vote for Party 4.}
#'   \item{V5}{The propensity to vote for Party 5.}
#'   \item{V6}{The propensity to vote for Party 6.}
#'   \item{V7}{The propensity to vote for Party 7.}
#'   \item{V8}{The propensity to vote for Party 8.}
#' }
#' @details
#' There are only 61 missing ratings in the 1000×8 matrix.
#' 
#' The dataset is particularly useful for demonstrating the use of the \code{smacofRect()} function in multidimensional scaling analysis.
#'
#' @usage data(denmarkEES2009)
#' @examples
#' data(denmarkEES2009) 
#' @seealso \code{\link[smacof]{smacofRect}} for more details on the \code{smacofRect} function.
#' 
#' @keywords datasets
#' @name denmarkEES2009
#' @docType data
NULL


#' @encoding UTF-8
#' @title Interest Group Ratings of Members of Congress (1959-1981)
#' @description The `interest1981` dataset is a subset of a much larger dataset compiled by Keith Poole, containing nearly 200,000 
#' interest group ratings of members of Congress between 1959 and 1981. These data have been extensively analyzed 
#' by Poole (1981, 1984, 1990) and Poole and Daniels (1985). The dataset has been used to perform the MLSMU6 
#' unfolding procedure in two dimensions, following the methodologies established in the mentioned studies.
#'
#' @details
#' This dataset contains interest group ratings that were used to perform the MLSMU6 unfolding procedure in two dimensions, 
#' as described in the works of Poole (1981, 1984, 1990) and Poole and Daniels (1985).
#'
#' @source
#' The data were compiled by Keith Poole and have been analyzed in several studies:
#'  Poole, K. T. (1981, 1984, 1990) and Poole, K. T., & Daniels, P. (1985).
#'
#' @usage data(interest1981)
#'
#' @examples
#' data(interest1981)
#'
#' @keywords datasets
#' @name interest1981
#' @docType data
NULL


#' @encoding UTF-8
#' @title 2004 American National Election Study (ANES) Data
#' @description The `ANES2004_OOC` dataset contains data from the 2004 American National Election Study (ANES). The 2004 ANES asked respondents about their policy preferences on issues ranging from diplomacy and defense spending to government spending and abortion.
#'
#' @format A data frame with rows representing respondents and columns representing their policy preferences on various issues.
#'
#' @details
#' This dataset is part of the 2004 American National Election Study (ANES), which surveyed respondents on a wide range of political and social issues. The data includes responses related to diplomacy, defense spending, government spending, and abortion, among others.
#'
#' @source
#' The data was collected as part of the 2004 American National Election Study (ANES).
#'
#' @usage data(ANES2004_OOC)
#'
#' @examples
#' \dontrun{
#' data(ANES2004_OOC)
#' }
#'
#' @keywords datasets
#' @name ANES2004_OOC
#' @docType data
NULL


#' @encoding UTF-8
#' @title Roll Call Data for U.S. Congress
#' @description The `rc_ep` object is a roll call dataset compiled by Poole and Rosenthal. This dataset represents roll call votes in the U.S. Congress and has been processed into an object of class `rollcall()` for analysis.
#'
#' @details
#' Poole and Rosenthal have compiled House and Senate roll call datasets covering the history of the U.S. Congress. These datasets are maintained at \url{http://www.voteview.com}. The `rc_ep` object represents a specific subset of these data and is formatted as a `rollcall` object, suitable for various forms of legislative analysis.
#'
#' @format An object of class \code{"rollcall"} containing roll call vote data.
#'
#' @source Poole, K. T., & Rosenthal, H. (Various Years). House and Senate Roll Call Data. Retrieved from \url{http://www.voteview.com}.
#'
#' @usage data(rc_ep)
#'
#' @examples
#' \dontrun{
#' data(rc_ep)
#' }
#'
#' @keywords datasets
#' @name rc_ep
#' @docType data
NULL


#' @encoding UTF-8
#' @title State of the Union Address Corpus
#' @description The `SOTUcorpus` dataset contains the text of each presidential State of the Union address since 1790. 
#' These data were collected and assembled by The American Presidency Project at the University of California, Santa Barbara.
#' The full dataset can be accessed and downloaded at \url{https://www.presidency.ucsb.edu/sou.php}.
#'
#' @details
#' The `SOTUcorpus` dataset provides a comprehensive collection of State of the Union addresses delivered by U.S. presidents from 1790 to the present.
#' These speeches are essential primary sources for studying American political rhetoric, policy priorities, and historical context across different administrations.
#'
#' @source
#' These data were collected and assembled by The American Presidency Project at the University of California, Santa Barbara.
#' The full dataset is available at \url{https://www.presidency.ucsb.edu/sou.php}.
#'
#' @usage data(SOTUcorpus)
#' @examples
#' \dontrun{
#' data(SOTUcorpus)
#' }
#'
#' @keywords datasets
#' @name SOTUcorpus
#' @docType data
NULL


#' @encoding UTF-8
#' @title 2000 Convention Delegate Study (CDS)
#' @description The `CDS2000` dataset contains data from the 2000 Convention Delegate Study (CDS), which interviewed delegates to the 
#' Republican and Democratic National Conventions. The survey included a battery of issue scales on which delegates were 
#' asked to place their policy preferences and those of major political figures (e.g., Al Gore and George W. Bush).
#'
#' @format A data frame with the following 14 variables:
#' \describe{
#'   \item{Party}{The party affiliation of the delegate (Democratic or Republican).}
#'   \item{Preferred Presidential Nominee}{The delegate's preferred presidential nominee.}
#'   \item{Race}{The race of the delegate.}
#'   \item{Religious Tradition}{The religious tradition of the delegate.}
#'   \item{Lib-Con}{Liberal-Conservative self-placement scale.}
#'   \item{Abortion}{Policy preference on abortion.}
#'   \item{Govt Services}{Policy preference on government services.}
#'   \item{Defense Spending}{Policy preference on defense spending.}
#'   \item{Aid to Blacks}{Policy preference on aid to Black Americans.}
#'   \item{Health Insurance}{Policy preference on health insurance.}
#'   \item{Protect Homosexuals}{Policy preference on protecting homosexual rights.}
#'   \item{Affirmative Action}{Policy preference on affirmative action.}
#'   \item{Surplus for Tax Cuts}{Policy preference on using budget surplus for tax cuts.}
#'   \item{Free Trade}{Policy preference on free trade.}
#' }
#'
#' @details
#' In 2000, the CDS received completed questionnaires from 1,907 delegates to the Democratic National Convention and 985 
#' delegates to the Republican National Convention. This dataset has been used in various studies analyzing delegate behavior 
#' and policy preferences, including works by Stone and Abramowitz (1983), Layman (2001), and Layman et al. (2010).
#'
#' @source
#' Data collected from the 2000 Convention Delegate Study (CDS). Studies that have analyzed CDS data include:
#'  Stone, W. J., & Abramowitz, A. I. (1983), Layman, G. C. (2001) and Layman, G. C., et al. (2010).
#'  
#' @usage data(CDS2000)
#'
#' @examples
#' \dontrun{
#' data(CDS2000)
#' }
#' @keywords datasets
#' @name CDS2000
#' @docType data
NULL


#' @encoding UTF-8
#' @title French Module of the 2009 European Election Study (EES)
#' @description The `franceEES2009` dataset contains data from the French module of the 2009 European Election Study (EES). The EES surveyed 1,000 French citizens, asking them to place themselves and eight major political parties on a 0-10 left-right scale (0 representing the most left-wing position, 10 representing the most right-wing position).
#'
#' @format A data frame with 1,000 rows and 9 columns:
#' \describe{
#'   \item{self}{Numeric, self-placement on the left-right scale.}
#'   \item{Extreme Left}{Numeric, placement of the "Extreme Left" party on the left-right scale.}
#'   \item{Communist}{Numeric, placement of the "Communist" party on the left-right scale.}
#'   \item{Socialist}{Numeric, placement of the "Socialist" party on the left-right scale.}
#'   \item{Greens}{Numeric, placement of the "Greens" party on the left-right scale.}
#'   \item{UDF (Bayrou)}{Numeric, placement of the "UDF (Bayrou)" party on the left-right scale.}
#'   \item{UMP (Sarkozy)}{Numeric, placement of the "UMP (Sarkozy)" party on the left-right scale.}
#'   \item{National Front}{Numeric, placement of the "National Front" party on the left-right scale.}
#'   \item{Left Party}{Numeric, placement of the "Left Party" on the left-right scale.}
#' }
#'
#' @details
#' Responses are coded as 77, 88, or 89 if the respondent refuses to answer or does not know the party or where to place it. The dataset is useful for analyzing the political landscape in France during 2009, especially in terms of how citizens and parties are positioned on the left-right spectrum.
#'
#' @source European Election Study (EES) 2009, French module.
#'
#' @usage data(franceEES2009)
#'
#' @examples
#' \dontrun{
#' data(franceEES2009)
#' }
#'
#' @keywords datasets
#' @name franceEES2009
#' @docType data
NULL


#' @encoding UTF-8
#' @title Mexican Module of the Comparative Study of Electoral Systems (CSES) 2000 and 2006
#' @description The `mexicoCSES2006` dataset contains data from the 2000 and 2006 Mexican modules of the Comparative Study of Electoral Systems (CSES). In these surveys, Mexican citizens were asked to place the major political parties on an 11-point left-right scale.
#'
#' @format A data frame with rows representing respondents and 8 columns representing different political parties:
#' \describe{
#'   \item{PAN}{Numeric, placement of the PAN (National Action Party) on the left-right scale.}
#'   \item{PRD}{Numeric, placement of the PRD (Party of the Democratic Revolution) on the left-right scale.}
#'   \item{PRI}{Numeric, placement of the PRI (Institutional Revolutionary Party) on the left-right scale.}
#'   \item{Greens}{Numeric, placement of the Green Party on the left-right scale.}
#'   \item{PT}{Numeric, placement of the PT (Labor Party) on the left-right scale.}
#'   \item{Convergencia}{Numeric, placement of the Convergencia (Convergence) party on the left-right scale.}
#'   \item{Nueva Alianza}{Numeric, placement of the Nueva Alianza (New Alliance) party on the left-right scale.}
#'   \item{PSD}{Numeric, placement of the PSD (Social Democratic Party) on the left-right scale.}
#' }
#'
#' @details
#' The dataset includes responses from the 2000 and 2006 Mexican modules of the Comparative Study of Electoral Systems (CSES), where respondents were asked to place major political parties on an 11-point left-right scale.
#'
#' @source Comparative Study of Electoral Systems (CSES), Mexican modules for 2000 and 2006.
#'
#' @usage data(mexicoCSES2006)
#'
#' @examples
#' \dontrun{
#' data(mexicoCSES2006)
#' }
#'
#' @keywords datasets
#' @name mexicoCSES2006
#' @docType data
NULL


#' @encoding UTF-8
#' @title 2008 U.S. Presidential Vote Data
#' @description The `presvote2008` dataset contains data on the voting behavior in the 2008 U.S. Presidential election, where a vote 
#' for John McCain is coded as 0 and a vote for Barack Obama is coded as 1.
#'
#' @details
#' This dataset provides a binary coding for the voting choice in the 2008 U.S. Presidential election, which can be used for various 
#' analyses, such as logistic regression, to study voting behavior and preferences.
#'
#' @usage data(presvote2008)
#'
#' @examples
#' \dontrun{
#' data(presvote2008)
#' }
#'
#' @keywords datasets
#' @name presvote2008
#' @docType data
NULL


#' @encoding UTF-8
#' @title Roll Call Data from the 108th US House of Representatives (2003-2005)
#' @description The `hr108` dataset contains data from the 108th US House of Representatives, covering the period from 2003 to 2005. During this session, the House conducted 843 recorded roll call votes, with 440 Representatives serving in the chamber. The roll call matrix omits President George W. Bush.
#'
#' @format A data frame or matrix with 440 rows (representing Representatives) and 843 columns (representing roll call votes). Each entry in the matrix indicates the vote of a Representative on a specific roll call.
#'
#' @details
#' This dataset provides a detailed record of the roll call votes conducted in the 108th US House of Representatives. It is useful for analyzing voting patterns, party alignment, and legislative behavior during this congressional session. Note that the roll call matrix excludes President George W. Bush.
#'
#' @source US House of Representatives, 108th Congress (2003-2005).
#'
#' @usage data(hr108)
#'
#' @examples
#' \dontrun{
#' data(hr108)
#' }
#'
#' @keywords datasets
#' @name hr108
#' @docType data
NULL


#' @encoding UTF-8
#' @title Vietnam War Issue Scales from the 1968 National Election Study (NES)
#' @description The `nes1968_vietnam` dataset contains responses from the 1968 National Election Study (NES) where respondents were asked to place themselves, President Lyndon Johnson, and the three major presidential candidates—Democrat Hubert Humphrey, Republican Richard Nixon, and American Independent George Wallace—on two seven-point issue scales regarding the Vietnam War.
#'
#' @format A data frame with rows representing respondents and the following columns:
#' \describe{
#'   \item{vote.choice}{Categorical, respondent's reported vote choice in the 1968 presidential election. This variable indicates which of the three major presidential candidates (Humphrey, Nixon, or Wallace) the respondent voted for, if any.}
#'   \item{johnson}{Numeric, placement of President Lyndon Johnson on the Vietnam War scale (1 to 7).}
#'   \item{humphrey}{Numeric, placement of Democratic candidate Hubert Humphrey on the Vietnam War scale (1 to 7).}
#'   \item{nixon}{Numeric, placement of Republican candidate Richard Nixon on the Vietnam War scale (1 to 7).}
#'   \item{wallace}{Numeric, placement of American Independent candidate George Wallace on the Vietnam War scale (1 to 7).}
#'   \item{self}{Numeric, respondent's self-placement on the Vietnam War scale (1 to 7).}
#' }
#' @details
#' This dataset is part of the 1968 National Election Study (NES). Respondents were asked to place themselves and key political figures on two seven-point scales relating to the Vietnam War. These scales measure opinions on how the war should be conducted or resolved.
#'
#' @source 1968 American National Election Study (ANES).
#'
#' @usage data(nes1968_vietnam)
#'
#' @examples
#' \dontrun{
#' data(nes1968_vietnam)
#' }
#'
#' @keywords datasets
#' @name nes1968_vietnam
#' @docType data
NULL


#' @encoding UTF-8
#' @title Urban Unrest Issue Scales from the 1968 National Election Study (NES)
#' @description The `nes1968_urbanunrest` dataset contains responses from the 1968 National Election Study (NES) where respondents were asked to place themselves, President Lyndon Johnson, and the three major presidential candidates—Democrat Hubert Humphrey, Republican Richard Nixon, and American Independent George Wallace—on two seven-point issue scales regarding urban unrest.
#'
#' @format A data frame with rows representing respondents and the following columns:
#' \describe{
#'   \item{vote.choice}{Categorical, respondent's reported vote choice in the 1968 presidential election. This variable indicates which of the three major presidential candidates (Humphrey, Nixon, or Wallace) the respondent voted for, if any.}
#'   \item{johnson}{Numeric, placement of President Lyndon Johnson on the Vietnam War scale (1 to 7).}
#'   \item{humphrey}{Numeric, placement of Democratic candidate Hubert Humphrey on the Vietnam War scale (1 to 7).}
#'   \item{nixon}{Numeric, placement of Republican candidate Richard Nixon on the Vietnam War scale (1 to 7).}
#'   \item{wallace}{Numeric, placement of American Independent candidate George Wallace on the Vietnam War scale (1 to 7).}
#'   \item{self}{Numeric, respondent's self-placement on the Vietnam War scale (1 to 7).}
#' }
#'
#' @details
#' This dataset is part of the 1968 National Election Study (NES). Respondents were asked to place themselves and key political figures on two seven-point scales relating to urban unrest. These scales measure opinions on how the issue of urban unrest should be addressed.
#'
#' @source 1968 American National Election Study (ANES).
#'
#' @usage data(nes1968_urbanunrest)
#'
#' @examples
#' \dontrun{
#' data(nes1968_urbanunrest)
#' }
#'
#' @keywords datasets
#' @name nes1968_urbanunrest
#' @docType data
NULL


#' @encoding UTF-8
#' @title 2010 Swedish Parliamentary Candidate Survey
#' @description The `Sweden2010` dataset contains data from the 2010 Swedish Parliamentary Candidate Survey, conducted by the Swedish 
#' public broadcasting network Sveriges Television (SVT). The survey targeted all 5,627 parliamentary candidates, with 
#' completed interviews from 2,830 candidates, including 289 of the 349 candidates who were elected.
#'
#' @format A data frame with the following variables:
#' \describe{
#'   \item{id}{Unique identifier for each candidate.}
#'   \item{elected}{Indicator of whether the candidate was elected (1) or not (0).}
#'   \item{party.name}{Name of the political party the candidate belongs to.}
#'   \item{party.code}{Numeric code representing the political party.}
#'   \item{govt.party}{Indicator of whether the candidate's party is part of the government coalition (1) or not (0).}
#'   \item{left.right.self.fivept}{Left-right self-placement on a 1-5 scale.}
#'   \item{congestion.taxes}{Opinion on congestion taxes.}
#'   \item{highspeed.trains}{Opinion on high-speed trains.}
#'   \item{hunt.wolves}{Opinion on hunting wolves.}
#'   \item{nuclear.power}{Opinion on nuclear power.}
#'   \item{gasoline.taxes}{Opinion on gasoline taxes.}
#'   \item{museum.fees}{Opinion on museum fees.}
#'   \item{online.piracy}{Opinion on online piracy.}
#'   \item{state.TV}{Opinion on state TV.}
#'   \item{refugee.cities}{Opinion on establishing refugee cities.}
#'   \item{asylum.seekers}{Opinion on asylum seekers.}
#'   \item{refugee.healthcare}{Opinion on healthcare for refugees.}
#'   \item{teacher.veils}{Opinion on teachers wearing veils.}
#'   \item{paternal.leave}{Opinion on paternal leave.}
#'   \item{affirmative.action.universities}{Opinion on affirmative action in universities.}
#'   \item{child.raising.allowance}{Opinion on child-raising allowance.}
#'   \item{property.taxes.wealthy}{Opinion on property taxes for the wealthy.}
#'   \item{wealth.tax}{Opinion on wealth tax.}
#'   \item{tax.wealthy}{Opinion on taxing the wealthy.}
#'   \item{tax.pensions}{Opinion on taxing pensions.}
#'   \item{household.services.deduction}{Opinion on household services deduction.}
#'   \item{work.income.tax}{Opinion on work income tax.}
#'   \item{sex.purchase}{Opinion on purchasing sex.}
#'   \item{DUI.penalty}{Opinion on penalties for driving under the influence (DUI).}
#'   \item{criminal.sentences}{Opinion on criminal sentences.}
#'   \item{wiretaps}{Opinion on wiretapping.}
#'   \item{retirement.age}{Opinion on retirement age.}
#'   \item{health.insurance.time}{Opinion on health insurance time limits.}
#'   \item{dental.insurance}{Opinion on dental insurance.}
#'   \item{competition.public.sector}{Opinion on competition in the public sector.}
#'   \item{mandatory.unemployment.insurance}{Opinion on mandatory unemployment insurance.}
#'   \item{municipal.home.care}{Opinion on municipal home care services.}
#'   \item{circumcision}{Opinion on circumcision.}
#'   \item{private.healthcare.profits}{Opinion on profits in private healthcare.}
#'   \item{govt.alcohol.monopoly}{Opinion on the government's alcohol monopoly.}
#'   \item{employment.protection}{Opinion on employment protection.}
#'   \item{sell.public.corporations}{Opinion on selling public corporations.}
#'   \item{Aghanistan.withdrawal}{Opinion on withdrawal from Afghanistan.}
#'   \item{exporting.arms}{Opinion on exporting arms.}
#'   \item{aid.undemocratic.countries}{Opinion on aid to undemocratic countries.}
#'   \item{compulsory.military.service}{Opinion on compulsory military service.}
#'   \item{leave.EU}{Opinion on leaving the European Union.}
#'   \item{transfer.students}{Opinion on transferring students.}
#'   \item{local.control.education}{Opinion on local control of education.}
#'   \item{student.grades}{Opinion on student grading systems.}
#'   \item{number.private.schools}{Opinion on the number of private schools.}
#'   \item{university.eligibility}{Opinion on university eligibility.}
#'   \item{criminalize.racist.organizations}{Opinion on criminalizing racist organizations.}
#'   \item{abolish.monarchy}{Opinion on abolishing the monarchy.}
#'   \item{ballot.order}{Opinion on the order of ballots.}
#'   \item{referenda.elections}{Opinion on referenda during elections.}
#' }
#'
#' @details
#' Candidates were asked 50 Likert-type questions, using a 4-point scale (from "strongly disagree" to "strongly agree") to register 
#' their opinions on a series of policy statements. Most issue scales focus on economic/social welfare issues, but questions related 
#' to foreign policy, social/cultural matters, law and order, immigration, and environmental issues are also included. 
#' Missing responses are coded as 8.
#'
#' @source
#' Data collected by the Swedish public broadcasting network Sveriges Television (SVT).
#'
#' @usage data(Sweden2010)
#'
#' @examples
#' \dontrun{
#' data(Sweden2010)
#' }
#'
#' @keywords datasets
#' @name Sweden2010
#' @docType data
NULL


#' @encoding UTF-8
#' @title Candidate Favorability Ratings from the 2008 ANES
#' @description The `candidatetherms2008` dataset contains favorability ratings from the 2008 American National Election Study (ANES). 
#' Respondents were asked to rate their favorability towards nine political figures and parties on a 0-100 scale.
#'
#' @format A matrix with rows representing respondents and columns representing the following nine political stimuli:
#' \describe{
#'   \item{mccain}{Favorability rating for Sen. John McCain.}
#'   \item{bush}{Favorability rating for Pres. George W. Bush.}
#'   \item{obama}{Favorability rating for Sen. Barack Obama.}
#'   \item{biden}{Favorability rating for Sen. Joe Biden.}
#'   \item{palin}{Favorability rating for Gov. Sarah Palin.}
#'   \item{hclinton}{Favorability rating for Sen. Hillary Clinton.}
#'   \item{bclinton}{Favorability rating for former Pres. Bill Clinton.}
#'   \item{demparty}{Favorability rating for the Democratic Party.}
#'   \item{repparty}{Favorability rating for the Republican Party.}
#' }
#'
#' @details
#' The 2008 ANES asked respondents to rate their favorability towards nine political figures and parties on a scale from 0 to 100, 
#' where 0 represents the least favorable and 100 represents the most favorable. Missing values are coded as `NA`.
#'
#' @source
#' Data from the 2008 American National Election Study (ANES).
#'
#' @usage data(candidatetherms2008)
#'
#' @examples
#' \dontrun{
#' data(candidatetherms2008)
#' }
#'
#' @keywords datasets
#' @name candidatetherms2008
#' @docType data
NULL


#' @encoding UTF-8
#' @title CHES EU Dataset: Party Positions from Chapel Hill Expert Survey
#' @description Party means and standard deviations from the 2010 Chapel Hill Expert Survey 
#' (CHES). Contains expert placements of 118 political parties from 14 EU 
#' member countries, including three anchoring vignette parties for scale 
#' calibration.
#'
#' @format A data frame with 121 rows (118 actual parties + 3 vignettes) and 
#'   5 variables:
#' \describe{
#'   \item{party}{Name or identifier of the political party (character)}
#'   \item{mean}{Mean party placement across all experts who evaluated this 
#'         party (numeric)}
#'   \item{sd}{Standard deviation of party placements across experts (numeric)}
#'   \item{country}{Country where the party is based. ISO country codes or 
#'         full country names (character)}
#'   \item{vignette}{Indicator for vignette parties: 1 = vignette party used 
#'         for anchoring, 0 = actual political party (numeric or logical)}
#' }
#'
#' @details
#' The 2010 Chapel Hill Expert Survey (CHES) included 224 experts evaluating 
#' parties from 14 EU member countries. Each actual party was evaluated by 
#' 8-17 experts, while all three vignette parties were evaluated by over 160 
#' experts. Vignette parties are hypothetical parties used as anchoring points 
#' to calibrate expert placements and adjust for differential item functioning 
#' (DIF).
#' 
#' This dataset is used to demonstrate anchoring vignette methods in spatial 
#' analysis, as discussed in Section 2.5 (pages 58-60) of the accompanying 
#' textbook. The vignettes allow researchers to adjust for systematic 
#' differences in how experts use the rating scales.
#'
#' @source
#' Bakker, R., De Vries, C., Edwards, E., Hooghe, L., Jolly, S., Marks, G., 
#' Polk, J., Rovny, J., Steenbergen, M., & Vachudova, M. A. (2015). 
#' Measuring party positions in Europe: The Chapel Hill expert survey trend 
#' file, 1999-2010. \emph{Party Politics}, 21(1), 143-152. 
#' \doi{10.1177/1354068812462931}
#' 
#' CHES data: \url{https://www.chesdata.eu/}
#'
#' @references
#' Bakker, R., et al. (2015). Measuring party positions in Europe: The Chapel 
#' Hill expert survey trend file, 1999-2010. \emph{Party Politics}, 21(1), 
#' 143-152.
#' 
#' King, G., & Wand, J. (2007). Comparing Incomparable Survey Responses: 
#' Evaluating and Selecting Anchoring Vignettes. \emph{Political Analysis}, 
#' 15(1), 46-66.
#'
#' @usage data(ches_eu)
#'
#' @examples
#' \dontrun{
#' data(ches_eu)
#' 
#' # Dataset structure
#' str(ches_eu)
#' dim(ches_eu)  # 121 parties (118 actual + 3 vignettes)
#' 
#' # Separate vignettes from actual parties
#' vignettes <- ches_eu[ches_eu$vignette == 1, ]
#' actual_parties <- ches_eu[ches_eu$vignette == 0, ]
#' 
#' # Summary statistics
#' summary(ches_eu$mean)
#' summary(ches_eu$sd)
#' 
#' # Parties by country
#' table(actual_parties$country)
#' 
#' # Compare uncertainty: actual parties vs. vignettes
#' boxplot(sd ~ vignette, data = ches_eu,
#'         names = c("Actual Parties", "Vignettes"),
#'         ylab = "Standard Deviation",
#'         main = "Expert Agreement: Actual vs. Vignette Parties")
#' 
#' # Plot party positions
#' plot(ches_eu$mean, ches_eu$sd,
#'      col = ifelse(ches_eu$vignette == 1, "red", "blue"),
#'      pch = ifelse(ches_eu$vignette == 1, 17, 16),
#'      xlab = "Mean Position", ylab = "Standard Deviation",
#'      main = "Party Positions from CHES 2010")
#' legend("topright", legend = c("Actual", "Vignette"),
#'        col = c("blue", "red"), pch = c(16, 17))
#' }
#'
#' @keywords datasets
#' @name ches_eu
#' @docType data
NULL


#' @encoding UTF-8
#' @title Mexican Political Party Positions on Left-Right Scale (2000 & 2006)
#' @description The `mexicoCSES2000` dataset contains data from the 2000 and 2006 Mexican modules of the Comparative Study of Electoral Systems (CSES).
#' In these surveys, Mexican citizens were asked to place the major political parties on an 11-point left-right scale.
#'
#' @format A data frame with the following variables representing major political parties in Mexico:
#' \describe{
#'   \item{PAN}{Left-right placement of the National Action Party (Partido Acción Nacional).}
#'   \item{PRI}{Left-right placement of the Institutional Revolutionary Party (Partido Revolucionario Institucional).}
#'   \item{PRD}{Left-right placement of the Party of the Democratic Revolution (Partido de la Revolución Democrática).}
#'   \item{PT}{Left-right placement of the Labor Party (Partido del Trabajo).}
#'   \item{Greens}{Left-right placement of the Green Party (Partido Verde Ecologista de México).}
#'   \item{PARM}{Left-right placement of the Authentic Party of the Mexican Revolution (Partido Auténtico de la Revolución Mexicana).}
#' }
#'
#' @details
#' The data were collected as part of the 2000 and 2006 Mexican modules of the Comparative Study of Electoral Systems (CSES). 
#' In these surveys, respondents were asked to rate the major political parties in Mexico on an 11-point left-right ideological scale.
#'
#' @source
#' Data from the 2000 and 2006 Mexican modules of the Comparative Study of Electoral Systems (CSES).
#'
#' @usage data(mexicoCSES2000)
#'
#' @examples
#' \dontrun{
#' data(mexicoCSES2000)
#' }
#'
#' @keywords datasets
#' @name mexicoCSES2000
#' @docType data
NULL


#' @encoding UTF-8
#' @title Transposed Rankings Data using Blackbox Method
#' @description  The `original` object is created by applying the `blackbox_transpose` function to a dataset of rankings.
#' This function is used to perform multidimensional scaling on a set of rankings, handling missing values 
#' and using specified dimensions and scaling parameters.
#'
#' @format A list object containing the results of the `blackbox_transpose` function.
#'
#' @details
#' The `original` object was created using the following call:
#' \preformatted{
#' original <- blackbox_transpose(rankings,
#'                                missing = c(77, 88, 89),
#'                                dims = 3,
#'                                minscale = 5,
#'                                verbose = FALSE)
#' }
#'
#' - **rankings**: A dataset containing the rankings to be analyzed.
#' - **missing**: A vector of values that are treated as missing data.
#' - **dims**: The number of dimensions to be used in the multidimensional scaling.
#' - **minscale**: The minimum scale to be applied in the scaling process.
#' - **verbose**: A logical flag indicating whether to print detailed output during the process.
#'
#' The `blackbox_transpose` function is typically used to analyze and visualize multidimensional data, allowing researchers 
#' to understand the underlying structure of ranked data by transforming it into a lower-dimensional space.
#'
#' @usage data(original)
#'
#' @examples
#' \dontrun{
#' data(original)
#' }
#'
#' @keywords datasets
#' @name original
#' @docType data
NULL

#' @encoding UTF-8
#' @title Roll Call Voting Data from the First European Parliament (1979-1984)
#' #' @description Roll call voting data from the first elected European Parliament (1979-1984), 
#' assembled by Hix, Noury, and Roland (2006). Contains voting records of 
#' Members of the European Parliament (MEPs) on 886 roll call votes.
#'
#' @format A data frame with 410 MEPs (rows) and 891 columns:
#' @format A data frame with 548 rows and 891 columns:
#' \describe{
#'   \item{EPG}{European Parliament Group (character)}
#'   \item{MEPID}{Member of European Parliament ID (numeric)}
#'   \item{MEPNAME}{Member of European Parliament name (character)}
#'   \item{MS}{Member State (character)}
#'   \item{NP}{National Party (character)}
#'   \item{V1}{Roll call vote (numeric)}
#'   \item{V2}{Roll call vote (numeric)}
#'   \item{V3}{Roll call vote (numeric)}
#'   \item{V4}{Roll call vote (numeric)}
#'   \item{V5}{Roll call vote (numeric)}
#'   \item{V6}{Roll call vote (numeric)}
#'   \item{V7}{Roll call vote (numeric)}
#'   \item{V8}{Roll call vote (numeric)}
#'   \item{V9}{Roll call vote (numeric)}
#'   \item{V10}{Roll call vote (numeric)}
#'   \item{V11}{Roll call vote (numeric)}
#'   \item{V12}{Roll call vote (numeric)}
#'   \item{V13}{Roll call vote (numeric)}
#'   \item{V14}{Roll call vote (numeric)}
#'   \item{V15}{Roll call vote (numeric)}
#'   \item{V16}{Roll call vote (numeric)}
#'   \item{V17}{Roll call vote (numeric)}
#'   \item{V18}{Roll call vote (numeric)}
#'   \item{V19}{Roll call vote (numeric)}
#'   \item{V20}{Roll call vote (numeric)}
#'   \item{V21}{Roll call vote (numeric)}
#'   \item{V22}{Roll call vote (numeric)}
#'   \item{V23}{Roll call vote (numeric)}
#'   \item{V24}{Roll call vote (numeric)}
#'   \item{V25}{Roll call vote (numeric)}
#'   \item{V26}{Roll call vote (numeric)}
#'   \item{V27}{Roll call vote (numeric)}
#'   \item{V28}{Roll call vote (numeric)}
#'   \item{V29}{Roll call vote (numeric)}
#'   \item{V30}{Roll call vote (numeric)}
#'   \item{V31}{Roll call vote (numeric)}
#'   \item{V32}{Roll call vote (numeric)}
#'   \item{V33}{Roll call vote (numeric)}
#'   \item{V34}{Roll call vote (numeric)}
#'   \item{V35}{Roll call vote (numeric)}
#'   \item{V36}{Roll call vote (numeric)}
#'   \item{V37}{Roll call vote (numeric)}
#'   \item{V38}{Roll call vote (numeric)}
#'   \item{V39}{Roll call vote (numeric)}
#'   \item{V40}{Roll call vote (numeric)}
#'   \item{V41}{Roll call vote (numeric)}
#'   \item{V42}{Roll call vote (numeric)}
#'   \item{V43}{Roll call vote (numeric)}
#'   \item{V44}{Roll call vote (numeric)}
#'   \item{V45}{Roll call vote (numeric)}
#'   \item{V46}{Roll call vote (numeric)}
#'   \item{V47}{Roll call vote (numeric)}
#'   \item{V48}{Roll call vote (numeric)}
#'   \item{V49}{Roll call vote (numeric)}
#'   \item{V50}{Roll call vote (numeric)}
#'   \item{V51}{Roll call vote (numeric)}
#'   \item{V52}{Roll call vote (numeric)}
#'   \item{V53}{Roll call vote (numeric)}
#'   \item{V54}{Roll call vote (numeric)}
#'   \item{V55}{Roll call vote (numeric)}
#'   \item{V56}{Roll call vote (numeric)}
#'   \item{V57}{Roll call vote (numeric)}
#'   \item{V58}{Roll call vote (numeric)}
#'   \item{V59}{Roll call vote (numeric)}
#'   \item{V60}{Roll call vote (numeric)}
#'   \item{V61}{Roll call vote (numeric)}
#'   \item{V62}{Roll call vote (numeric)}
#'   \item{V63}{Roll call vote (numeric)}
#'   \item{V64}{Roll call vote (numeric)}
#'   \item{V65}{Roll call vote (numeric)}
#'   \item{V66}{Roll call vote (numeric)}
#'   \item{V67}{Roll call vote (numeric)}
#'   \item{V68}{Roll call vote (numeric)}
#'   \item{V69}{Roll call vote (numeric)}
#'   \item{V70}{Roll call vote (numeric)}
#'   \item{V71}{Roll call vote (numeric)}
#'   \item{V72}{Roll call vote (numeric)}
#'   \item{V73}{Roll call vote (numeric)}
#'   \item{V74}{Roll call vote (numeric)}
#'   \item{V75}{Roll call vote (numeric)}
#'   \item{V76}{Roll call vote (numeric)}
#'   \item{V77}{Roll call vote (numeric)}
#'   \item{V78}{Roll call vote (numeric)}
#'   \item{V79}{Roll call vote (numeric)}
#'   \item{V80}{Roll call vote (numeric)}
#'   \item{V81}{Roll call vote (numeric)}
#'   \item{V82}{Roll call vote (numeric)}
#'   \item{V83}{Roll call vote (numeric)}
#'   \item{V84}{Roll call vote (numeric)}
#'   \item{V85}{Roll call vote (numeric)}
#'   \item{V86}{Roll call vote (numeric)}
#'   \item{V87}{Roll call vote (numeric)}
#'   \item{V88}{Roll call vote (numeric)}
#'   \item{V89}{Roll call vote (numeric)}
#'   \item{V90}{Roll call vote (numeric)}
#'   \item{V91}{Roll call vote (numeric)}
#'   \item{V92}{Roll call vote (numeric)}
#'   \item{V93}{Roll call vote (numeric)}
#'   \item{V94}{Roll call vote (numeric)}
#'   \item{V95}{Roll call vote (numeric)}
#'   \item{V96}{Roll call vote (numeric)}
#'   \item{V97}{Roll call vote (numeric)}
#'   \item{V98}{Roll call vote (numeric)}
#'   \item{V99}{Roll call vote (numeric)}
#'   \item{V100}{Roll call vote (numeric)}
#'   \item{V101}{Roll call vote (numeric)}
#'   \item{V102}{Roll call vote (numeric)}
#'   \item{V103}{Roll call vote (numeric)}
#'   \item{V104}{Roll call vote (numeric)}
#'   \item{V105}{Roll call vote (numeric)}
#'   \item{V106}{Roll call vote (numeric)}
#'   \item{V107}{Roll call vote (numeric)}
#'   \item{V108}{Roll call vote (numeric)}
#'   \item{V109}{Roll call vote (numeric)}
#'   \item{V110}{Roll call vote (numeric)}
#'   \item{V111}{Roll call vote (numeric)}
#'   \item{V112}{Roll call vote (numeric)}
#'   \item{V113}{Roll call vote (numeric)}
#'   \item{V114}{Roll call vote (numeric)}
#'   \item{V115}{Roll call vote (numeric)}
#'   \item{V116}{Roll call vote (numeric)}
#'   \item{V117}{Roll call vote (numeric)}
#'   \item{V118}{Roll call vote (numeric)}
#'   \item{V119}{Roll call vote (numeric)}
#'   \item{V120}{Roll call vote (numeric)}
#'   \item{V121}{Roll call vote (numeric)}
#'   \item{V122}{Roll call vote (numeric)}
#'   \item{V123}{Roll call vote (numeric)}
#'   \item{V124}{Roll call vote (numeric)}
#'   \item{V125}{Roll call vote (numeric)}
#'   \item{V126}{Roll call vote (numeric)}
#'   \item{V127}{Roll call vote (numeric)}
#'   \item{V128}{Roll call vote (numeric)}
#'   \item{V129}{Roll call vote (numeric)}
#'   \item{V130}{Roll call vote (numeric)}
#'   \item{V131}{Roll call vote (numeric)}
#'   \item{V132}{Roll call vote (numeric)}
#'   \item{V133}{Roll call vote (numeric)}
#'   \item{V134}{Roll call vote (numeric)}
#'   \item{V135}{Roll call vote (numeric)}
#'   \item{V136}{Roll call vote (numeric)}
#'   \item{V137}{Roll call vote (numeric)}
#'   \item{V138}{Roll call vote (numeric)}
#'   \item{V139}{Roll call vote (numeric)}
#'   \item{V140}{Roll call vote (numeric)}
#'   \item{V141}{Roll call vote (numeric)}
#'   \item{V142}{Roll call vote (numeric)}
#'   \item{V143}{Roll call vote (numeric)}
#'   \item{V144}{Roll call vote (numeric)}
#'   \item{V145}{Roll call vote (numeric)}
#'   \item{V146}{Roll call vote (numeric)}
#'   \item{V147}{Roll call vote (numeric)}
#'   \item{V148}{Roll call vote (numeric)}
#'   \item{V149}{Roll call vote (numeric)}
#'   \item{V150}{Roll call vote (numeric)}
#'   \item{V151}{Roll call vote (numeric)}
#'   \item{V152}{Roll call vote (numeric)}
#'   \item{V153}{Roll call vote (numeric)}
#'   \item{V154}{Roll call vote (numeric)}
#'   \item{V155}{Roll call vote (numeric)}
#'   \item{V156}{Roll call vote (numeric)}
#'   \item{V157}{Roll call vote (numeric)}
#'   \item{V158}{Roll call vote (numeric)}
#'   \item{V159}{Roll call vote (numeric)}
#'   \item{V160}{Roll call vote (numeric)}
#'   \item{V161}{Roll call vote (numeric)}
#'   \item{V162}{Roll call vote (numeric)}
#'   \item{V163}{Roll call vote (numeric)}
#'   \item{V164}{Roll call vote (numeric)}
#'   \item{V165}{Roll call vote (numeric)}
#'   \item{V166}{Roll call vote (numeric)}
#'   \item{V167}{Roll call vote (numeric)}
#'   \item{V168}{Roll call vote (numeric)}
#'   \item{V169}{Roll call vote (numeric)}
#'   \item{V170}{Roll call vote (numeric)}
#'   \item{V171}{Roll call vote (numeric)}
#'   \item{V172}{Roll call vote (numeric)}
#'   \item{V173}{Roll call vote (numeric)}
#'   \item{V174}{Roll call vote (numeric)}
#'   \item{V175}{Roll call vote (numeric)}
#'   \item{V176}{Roll call vote (numeric)}
#'   \item{V177}{Roll call vote (numeric)}
#'   \item{V178}{Roll call vote (numeric)}
#'   \item{V179}{Roll call vote (numeric)}
#'   \item{V180}{Roll call vote (numeric)}
#'   \item{V181}{Roll call vote (numeric)}
#'   \item{V182}{Roll call vote (numeric)}
#'   \item{V183}{Roll call vote (numeric)}
#'   \item{V184}{Roll call vote (numeric)}
#'   \item{V185}{Roll call vote (numeric)}
#'   \item{V186}{Roll call vote (numeric)}
#'   \item{V187}{Roll call vote (numeric)}
#'   \item{V188}{Roll call vote (numeric)}
#'   \item{V189}{Roll call vote (numeric)}
#'   \item{V190}{Roll call vote (numeric)}
#'   \item{V191}{Roll call vote (numeric)}
#'   \item{V192}{Roll call vote (numeric)}
#'   \item{V193}{Roll call vote (numeric)}
#'   \item{V194}{Roll call vote (numeric)}
#'   \item{V195}{Roll call vote (numeric)}
#'   \item{V196}{Roll call vote (numeric)}
#'   \item{V197}{Roll call vote (numeric)}
#'   \item{V198}{Roll call vote (numeric)}
#'   \item{V199}{Roll call vote (numeric)}
#'   \item{V200}{Roll call vote (numeric)}
#'   \item{V201}{Roll call vote (numeric)}
#'   \item{V202}{Roll call vote (numeric)}
#'   \item{V203}{Roll call vote (numeric)}
#'   \item{V204}{Roll call vote (numeric)}
#'   \item{V205}{Roll call vote (numeric)}
#'   \item{V206}{Roll call vote (numeric)}
#'   \item{V207}{Roll call vote (numeric)}
#'   \item{V208}{Roll call vote (numeric)}
#'   \item{V209}{Roll call vote (numeric)}
#'   \item{V210}{Roll call vote (numeric)}
#'   \item{V211}{Roll call vote (numeric)}
#'   \item{V212}{Roll call vote (numeric)}
#'   \item{V213}{Roll call vote (numeric)}
#'   \item{V214}{Roll call vote (numeric)}
#'   \item{V215}{Roll call vote (numeric)}
#'   \item{V216}{Roll call vote (numeric)}
#'   \item{V217}{Roll call vote (numeric)}
#'   \item{V218}{Roll call vote (numeric)}
#'   \item{V219}{Roll call vote (numeric)}
#'   \item{V220}{Roll call vote (numeric)}
#'   \item{V221}{Roll call vote (numeric)}
#'   \item{V222}{Roll call vote (numeric)}
#'   \item{V223}{Roll call vote (numeric)}
#'   \item{V224}{Roll call vote (numeric)}
#'   \item{V225}{Roll call vote (numeric)}
#'   \item{V226}{Roll call vote (numeric)}
#'   \item{V227}{Roll call vote (numeric)}
#'   \item{V228}{Roll call vote (numeric)}
#'   \item{V229}{Roll call vote (numeric)}
#'   \item{V230}{Roll call vote (numeric)}
#'   \item{V231}{Roll call vote (numeric)}
#'   \item{V232}{Roll call vote (numeric)}
#'   \item{V233}{Roll call vote (numeric)}
#'   \item{V234}{Roll call vote (numeric)}
#'   \item{V235}{Roll call vote (numeric)}
#'   \item{V236}{Roll call vote (numeric)}
#'   \item{V237}{Roll call vote (numeric)}
#'   \item{V238}{Roll call vote (numeric)}
#'   \item{V239}{Roll call vote (numeric)}
#'   \item{V240}{Roll call vote (numeric)}
#'   \item{V241}{Roll call vote (numeric)}
#'   \item{V242}{Roll call vote (numeric)}
#'   \item{V243}{Roll call vote (numeric)}
#'   \item{V244}{Roll call vote (numeric)}
#'   \item{V245}{Roll call vote (numeric)}
#'   \item{V246}{Roll call vote (numeric)}
#'   \item{V247}{Roll call vote (numeric)}
#'   \item{V248}{Roll call vote (numeric)}
#'   \item{V249}{Roll call vote (numeric)}
#'   \item{V250}{Roll call vote (numeric)}
#'   \item{V251}{Roll call vote (numeric)}
#'   \item{V252}{Roll call vote (numeric)}
#'   \item{V253}{Roll call vote (numeric)}
#'   \item{V254}{Roll call vote (numeric)}
#'   \item{V255}{Roll call vote (numeric)}
#'   \item{V256}{Roll call vote (numeric)}
#'   \item{V257}{Roll call vote (numeric)}
#'   \item{V258}{Roll call vote (numeric)}
#'   \item{V259}{Roll call vote (numeric)}
#'   \item{V260}{Roll call vote (numeric)}
#'   \item{V261}{Roll call vote (numeric)}
#'   \item{V262}{Roll call vote (numeric)}
#'   \item{V263}{Roll call vote (numeric)}
#'   \item{V264}{Roll call vote (numeric)}
#'   \item{V265}{Roll call vote (numeric)}
#'   \item{V266}{Roll call vote (numeric)}
#'   \item{V267}{Roll call vote (numeric)}
#'   \item{V268}{Roll call vote (numeric)}
#'   \item{V269}{Roll call vote (numeric)}
#'   \item{V270}{Roll call vote (numeric)}
#'   \item{V271}{Roll call vote (numeric)}
#'   \item{V272}{Roll call vote (numeric)}
#'   \item{V273}{Roll call vote (numeric)}
#'   \item{V274}{Roll call vote (numeric)}
#'   \item{V275}{Roll call vote (numeric)}
#'   \item{V276}{Roll call vote (numeric)}
#'   \item{V277}{Roll call vote (numeric)}
#'   \item{V278}{Roll call vote (numeric)}
#'   \item{V279}{Roll call vote (numeric)}
#'   \item{V280}{Roll call vote (numeric)}
#'   \item{V281}{Roll call vote (numeric)}
#'   \item{V282}{Roll call vote (numeric)}
#'   \item{V283}{Roll call vote (numeric)}
#'   \item{V284}{Roll call vote (numeric)}
#'   \item{V285}{Roll call vote (numeric)}
#'   \item{V286}{Roll call vote (numeric)}
#'   \item{V287}{Roll call vote (numeric)}
#'   \item{V288}{Roll call vote (numeric)}
#'   \item{V289}{Roll call vote (numeric)}
#'   \item{V290}{Roll call vote (numeric)}
#'   \item{V291}{Roll call vote (numeric)}
#'   \item{V292}{Roll call vote (numeric)}
#'   \item{V293}{Roll call vote (numeric)}
#'   \item{V294}{Roll call vote (numeric)}
#'   \item{V295}{Roll call vote (numeric)}
#'   \item{V296}{Roll call vote (numeric)}
#'   \item{V297}{Roll call vote (numeric)}
#'   \item{V298}{Roll call vote (numeric)}
#'   \item{V299}{Roll call vote (numeric)}
#'   \item{V300}{Roll call vote (numeric)}
#'   \item{V301}{Roll call vote (numeric)}
#'   \item{V302}{Roll call vote (numeric)}
#'   \item{V303}{Roll call vote (numeric)}
#'   \item{V304}{Roll call vote (numeric)}
#'   \item{V305}{Roll call vote (numeric)}
#'   \item{V306}{Roll call vote (numeric)}
#'   \item{V307}{Roll call vote (numeric)}
#'   \item{V308}{Roll call vote (numeric)}
#'   \item{V309}{Roll call vote (numeric)}
#'   \item{V310}{Roll call vote (numeric)}
#'   \item{V311}{Roll call vote (numeric)}
#'   \item{V312}{Roll call vote (numeric)}
#'   \item{V313}{Roll call vote (numeric)}
#'   \item{V314}{Roll call vote (numeric)}
#'   \item{V315}{Roll call vote (numeric)}
#'   \item{V316}{Roll call vote (numeric)}
#'   \item{V317}{Roll call vote (numeric)}
#'   \item{V318}{Roll call vote (numeric)}
#'   \item{V319}{Roll call vote (numeric)}
#'   \item{V320}{Roll call vote (numeric)}
#'   \item{V321}{Roll call vote (numeric)}
#'   \item{V322}{Roll call vote (numeric)}
#'   \item{V323}{Roll call vote (numeric)}
#'   \item{V324}{Roll call vote (numeric)}
#'   \item{V325}{Roll call vote (numeric)}
#'   \item{V326}{Roll call vote (numeric)}
#'   \item{V327}{Roll call vote (numeric)}
#'   \item{V328}{Roll call vote (numeric)}
#'   \item{V329}{Roll call vote (numeric)}
#'   \item{V330}{Roll call vote (numeric)}
#'   \item{V331}{Roll call vote (numeric)}
#'   \item{V332}{Roll call vote (numeric)}
#'   \item{V333}{Roll call vote (numeric)}
#'   \item{V334}{Roll call vote (numeric)}
#'   \item{V335}{Roll call vote (numeric)}
#'   \item{V336}{Roll call vote (numeric)}
#'   \item{V337}{Roll call vote (numeric)}
#'   \item{V338}{Roll call vote (numeric)}
#'   \item{V339}{Roll call vote (numeric)}
#'   \item{V340}{Roll call vote (numeric)}
#'   \item{V341}{Roll call vote (numeric)}
#'   \item{V342}{Roll call vote (numeric)}
#'   \item{V343}{Roll call vote (numeric)}
#'   \item{V344}{Roll call vote (numeric)}
#'   \item{V345}{Roll call vote (numeric)}
#'   \item{V346}{Roll call vote (numeric)}
#'   \item{V347}{Roll call vote (numeric)}
#'   \item{V348}{Roll call vote (numeric)}
#'   \item{V349}{Roll call vote (numeric)}
#'   \item{V350}{Roll call vote (numeric)}
#'   \item{V351}{Roll call vote (numeric)}
#'   \item{V352}{Roll call vote (numeric)}
#'   \item{V353}{Roll call vote (numeric)}
#'   \item{V354}{Roll call vote (numeric)}
#'   \item{V355}{Roll call vote (numeric)}
#'   \item{V356}{Roll call vote (numeric)}
#'   \item{V357}{Roll call vote (numeric)}
#'   \item{V358}{Roll call vote (numeric)}
#'   \item{V359}{Roll call vote (numeric)}
#'   \item{V360}{Roll call vote (numeric)}
#'   \item{V361}{Roll call vote (numeric)}
#'   \item{V362}{Roll call vote (numeric)}
#'   \item{V363}{Roll call vote (numeric)}
#'   \item{V364}{Roll call vote (numeric)}
#'   \item{V365}{Roll call vote (numeric)}
#'   \item{V366}{Roll call vote (numeric)}
#'   \item{V367}{Roll call vote (numeric)}
#'   \item{V368}{Roll call vote (numeric)}
#'   \item{V369}{Roll call vote (numeric)}
#'   \item{V370}{Roll call vote (numeric)}
#'   \item{V371}{Roll call vote (numeric)}
#'   \item{V372}{Roll call vote (numeric)}
#'   \item{V373}{Roll call vote (numeric)}
#'   \item{V374}{Roll call vote (numeric)}
#'   \item{V375}{Roll call vote (numeric)}
#'   \item{V376}{Roll call vote (numeric)}
#'   \item{V377}{Roll call vote (numeric)}
#'   \item{V378}{Roll call vote (numeric)}
#'   \item{V379}{Roll call vote (numeric)}
#'   \item{V380}{Roll call vote (numeric)}
#'   \item{V381}{Roll call vote (numeric)}
#'   \item{V382}{Roll call vote (numeric)}
#'   \item{V383}{Roll call vote (numeric)}
#'   \item{V384}{Roll call vote (numeric)}
#'   \item{V385}{Roll call vote (numeric)}
#'   \item{V386}{Roll call vote (numeric)}
#'   \item{V387}{Roll call vote (numeric)}
#'   \item{V388}{Roll call vote (numeric)}
#'   \item{V389}{Roll call vote (numeric)}
#'   \item{V390}{Roll call vote (numeric)}
#'   \item{V391}{Roll call vote (numeric)}
#'   \item{V392}{Roll call vote (numeric)}
#'   \item{V393}{Roll call vote (numeric)}
#'   \item{V394}{Roll call vote (numeric)}
#'   \item{V395}{Roll call vote (numeric)}
#'   \item{V396}{Roll call vote (numeric)}
#'   \item{V397}{Roll call vote (numeric)}
#'   \item{V398}{Roll call vote (numeric)}
#'   \item{V399}{Roll call vote (numeric)}
#'   \item{V400}{Roll call vote (numeric)}
#'   \item{V401}{Roll call vote (numeric)}
#'   \item{V402}{Roll call vote (numeric)}
#'   \item{V403}{Roll call vote (numeric)}
#'   \item{V404}{Roll call vote (numeric)}
#'   \item{V405}{Roll call vote (numeric)}
#'   \item{V406}{Roll call vote (numeric)}
#'   \item{V407}{Roll call vote (numeric)}
#'   \item{V408}{Roll call vote (numeric)}
#'   \item{V409}{Roll call vote (numeric)}
#'   \item{V410}{Roll call vote (numeric)}
#'   \item{V411}{Roll call vote (numeric)}
#'   \item{V412}{Roll call vote (numeric)}
#'   \item{V413}{Roll call vote (numeric)}
#'   \item{V414}{Roll call vote (numeric)}
#'   \item{V415}{Roll call vote (numeric)}
#'   \item{V416}{Roll call vote (numeric)}
#'   \item{V417}{Roll call vote (numeric)}
#'   \item{V418}{Roll call vote (numeric)}
#'   \item{V419}{Roll call vote (numeric)}
#'   \item{V420}{Roll call vote (numeric)}
#'   \item{V421}{Roll call vote (numeric)}
#'   \item{V422}{Roll call vote (numeric)}
#'   \item{V423}{Roll call vote (numeric)}
#'   \item{V424}{Roll call vote (numeric)}
#'   \item{V425}{Roll call vote (numeric)}
#'   \item{V426}{Roll call vote (numeric)}
#'   \item{V427}{Roll call vote (numeric)}
#'   \item{V428}{Roll call vote (numeric)}
#'   \item{V429}{Roll call vote (numeric)}
#'   \item{V430}{Roll call vote (numeric)}
#'   \item{V431}{Roll call vote (numeric)}
#'   \item{V432}{Roll call vote (numeric)}
#'   \item{V433}{Roll call vote (numeric)}
#'   \item{V434}{Roll call vote (numeric)}
#'   \item{V435}{Roll call vote (numeric)}
#'   \item{V436}{Roll call vote (numeric)}
#'   \item{V437}{Roll call vote (numeric)}
#'   \item{V438}{Roll call vote (numeric)}
#'   \item{V439}{Roll call vote (numeric)}
#'   \item{V440}{Roll call vote (numeric)}
#'   \item{V441}{Roll call vote (numeric)}
#'   \item{V442}{Roll call vote (numeric)}
#'   \item{V443}{Roll call vote (numeric)}
#'   \item{V444}{Roll call vote (numeric)}
#'   \item{V445}{Roll call vote (numeric)}
#'   \item{V446}{Roll call vote (numeric)}
#'   \item{V447}{Roll call vote (numeric)}
#'   \item{V448}{Roll call vote (numeric)}
#'   \item{V449}{Roll call vote (numeric)}
#'   \item{V450}{Roll call vote (numeric)}
#'   \item{V451}{Roll call vote (numeric)}
#'   \item{V452}{Roll call vote (numeric)}
#'   \item{V453}{Roll call vote (numeric)}
#'   \item{V454}{Roll call vote (numeric)}
#'   \item{V455}{Roll call vote (numeric)}
#'   \item{V456}{Roll call vote (numeric)}
#'   \item{V457}{Roll call vote (numeric)}
#'   \item{V458}{Roll call vote (numeric)}
#'   \item{V459}{Roll call vote (numeric)}
#'   \item{V460}{Roll call vote (numeric)}
#'   \item{V461}{Roll call vote (numeric)}
#'   \item{V462}{Roll call vote (numeric)}
#'   \item{V463}{Roll call vote (numeric)}
#'   \item{V464}{Roll call vote (numeric)}
#'   \item{V465}{Roll call vote (numeric)}
#'   \item{V466}{Roll call vote (numeric)}
#'   \item{V467}{Roll call vote (numeric)}
#'   \item{V468}{Roll call vote (numeric)}
#'   \item{V469}{Roll call vote (numeric)}
#'   \item{V470}{Roll call vote (numeric)}
#'   \item{V471}{Roll call vote (numeric)}
#'   \item{V472}{Roll call vote (numeric)}
#'   \item{V473}{Roll call vote (numeric)}
#'   \item{V474}{Roll call vote (numeric)}
#'   \item{V475}{Roll call vote (numeric)}
#'   \item{V476}{Roll call vote (numeric)}
#'   \item{V477}{Roll call vote (numeric)}
#'   \item{V478}{Roll call vote (numeric)}
#'   \item{V479}{Roll call vote (numeric)}
#'   \item{V480}{Roll call vote (numeric)}
#'   \item{V481}{Roll call vote (numeric)}
#'   \item{V482}{Roll call vote (numeric)}
#'   \item{V483}{Roll call vote (numeric)}
#'   \item{V484}{Roll call vote (numeric)}
#'   \item{V485}{Roll call vote (numeric)}
#'   \item{V486}{Roll call vote (numeric)}
#'   \item{V487}{Roll call vote (numeric)}
#'   \item{V488}{Roll call vote (numeric)}
#'   \item{V489}{Roll call vote (numeric)}
#'   \item{V490}{Roll call vote (numeric)}
#'   \item{V491}{Roll call vote (numeric)}
#'   \item{V492}{Roll call vote (numeric)}
#'   \item{V493}{Roll call vote (numeric)}
#'   \item{V494}{Roll call vote (numeric)}
#'   \item{V495}{Roll call vote (numeric)}
#'   \item{V496}{Roll call vote (numeric)}
#'   \item{V497}{Roll call vote (numeric)}
#'   \item{V498}{Roll call vote (numeric)}
#'   \item{V499}{Roll call vote (numeric)}
#'   \item{V500}{Roll call vote (numeric)}
#'   \item{V501}{Roll call vote (numeric)}
#'   \item{V502}{Roll call vote (numeric)}
#'   \item{V503}{Roll call vote (numeric)}
#'   \item{V504}{Roll call vote (numeric)}
#'   \item{V505}{Roll call vote (numeric)}
#'   \item{V506}{Roll call vote (numeric)}
#'   \item{V507}{Roll call vote (numeric)}
#'   \item{V508}{Roll call vote (numeric)}
#'   \item{V509}{Roll call vote (numeric)}
#'   \item{V510}{Roll call vote (numeric)}
#'   \item{V511}{Roll call vote (numeric)}
#'   \item{V512}{Roll call vote (numeric)}
#'   \item{V513}{Roll call vote (numeric)}
#'   \item{V514}{Roll call vote (numeric)}
#'   \item{V515}{Roll call vote (numeric)}
#'   \item{V516}{Roll call vote (numeric)}
#'   \item{V517}{Roll call vote (numeric)}
#'   \item{V518}{Roll call vote (numeric)}
#'   \item{V519}{Roll call vote (numeric)}
#'   \item{V520}{Roll call vote (numeric)}
#'   \item{V521}{Roll call vote (numeric)}
#'   \item{V522}{Roll call vote (numeric)}
#'   \item{V523}{Roll call vote (numeric)}
#'   \item{V524}{Roll call vote (numeric)}
#'   \item{V525}{Roll call vote (numeric)}
#'   \item{V526}{Roll call vote (numeric)}
#'   \item{V527}{Roll call vote (numeric)}
#'   \item{V528}{Roll call vote (numeric)}
#'   \item{V529}{Roll call vote (numeric)}
#'   \item{V530}{Roll call vote (numeric)}
#'   \item{V531}{Roll call vote (numeric)}
#'   \item{V532}{Roll call vote (numeric)}
#'   \item{V533}{Roll call vote (numeric)}
#'   \item{V534}{Roll call vote (numeric)}
#'   \item{V535}{Roll call vote (numeric)}
#'   \item{V536}{Roll call vote (numeric)}
#'   \item{V537}{Roll call vote (numeric)}
#'   \item{V538}{Roll call vote (numeric)}
#'   \item{V539}{Roll call vote (numeric)}
#'   \item{V540}{Roll call vote (numeric)}
#'   \item{V541}{Roll call vote (numeric)}
#'   \item{V542}{Roll call vote (numeric)}
#'   \item{V543}{Roll call vote (numeric)}
#'   \item{V544}{Roll call vote (numeric)}
#'   \item{V545}{Roll call vote (numeric)}
#'   \item{V546}{Roll call vote (numeric)}
#'   \item{V547}{Roll call vote (numeric)}
#'   \item{V548}{Roll call vote (numeric)}
#'   \item{V549}{Roll call vote (numeric)}
#'   \item{V550}{Roll call vote (numeric)}
#'   \item{V551}{Roll call vote (numeric)}
#'   \item{V552}{Roll call vote (numeric)}
#'   \item{V553}{Roll call vote (numeric)}
#'   \item{V554}{Roll call vote (numeric)}
#'   \item{V555}{Roll call vote (numeric)}
#'   \item{V556}{Roll call vote (numeric)}
#'   \item{V557}{Roll call vote (numeric)}
#'   \item{V558}{Roll call vote (numeric)}
#'   \item{V559}{Roll call vote (numeric)}
#'   \item{V560}{Roll call vote (numeric)}
#'   \item{V561}{Roll call vote (numeric)}
#'   \item{V562}{Roll call vote (numeric)}
#'   \item{V563}{Roll call vote (numeric)}
#'   \item{V564}{Roll call vote (numeric)}
#'   \item{V565}{Roll call vote (numeric)}
#'   \item{V566}{Roll call vote (numeric)}
#'   \item{V567}{Roll call vote (numeric)}
#'   \item{V568}{Roll call vote (numeric)}
#'   \item{V569}{Roll call vote (numeric)}
#'   \item{V570}{Roll call vote (numeric)}
#'   \item{V571}{Roll call vote (numeric)}
#'   \item{V572}{Roll call vote (numeric)}
#'   \item{V573}{Roll call vote (numeric)}
#'   \item{V574}{Roll call vote (numeric)}
#'   \item{V575}{Roll call vote (numeric)}
#'   \item{V576}{Roll call vote (numeric)}
#'   \item{V577}{Roll call vote (numeric)}
#'   \item{V578}{Roll call vote (numeric)}
#'   \item{V579}{Roll call vote (numeric)}
#'   \item{V580}{Roll call vote (numeric)}
#'   \item{V581}{Roll call vote (numeric)}
#'   \item{V582}{Roll call vote (numeric)}
#'   \item{V583}{Roll call vote (numeric)}
#'   \item{V584}{Roll call vote (numeric)}
#'   \item{V585}{Roll call vote (numeric)}
#'   \item{V586}{Roll call vote (numeric)}
#'   \item{V587}{Roll call vote (numeric)}
#'   \item{V588}{Roll call vote (numeric)}
#'   \item{V589}{Roll call vote (numeric)}
#'   \item{V590}{Roll call vote (numeric)}
#'   \item{V591}{Roll call vote (numeric)}
#'   \item{V592}{Roll call vote (numeric)}
#'   \item{V593}{Roll call vote (numeric)}
#'   \item{V594}{Roll call vote (numeric)}
#'   \item{V595}{Roll call vote (numeric)}
#'   \item{V596}{Roll call vote (numeric)}
#'   \item{V597}{Roll call vote (numeric)}
#'   \item{V598}{Roll call vote (numeric)}
#'   \item{V599}{Roll call vote (numeric)}
#'   \item{V600}{Roll call vote (numeric)}
#'   \item{V601}{Roll call vote (numeric)}
#'   \item{V602}{Roll call vote (numeric)}
#'   \item{V603}{Roll call vote (numeric)}
#'   \item{V604}{Roll call vote (numeric)}
#'   \item{V605}{Roll call vote (numeric)}
#'   \item{V606}{Roll call vote (numeric)}
#'   \item{V607}{Roll call vote (numeric)}
#'   \item{V608}{Roll call vote (numeric)}
#'   \item{V609}{Roll call vote (numeric)}
#'   \item{V610}{Roll call vote (numeric)}
#'   \item{V611}{Roll call vote (numeric)}
#'   \item{V612}{Roll call vote (numeric)}
#'   \item{V613}{Roll call vote (numeric)}
#'   \item{V614}{Roll call vote (numeric)}
#'   \item{V615}{Roll call vote (numeric)}
#'   \item{V616}{Roll call vote (numeric)}
#'   \item{V617}{Roll call vote (numeric)}
#'   \item{V618}{Roll call vote (numeric)}
#'   \item{V619}{Roll call vote (numeric)}
#'   \item{V620}{Roll call vote (numeric)}
#'   \item{V621}{Roll call vote (numeric)}
#'   \item{V622}{Roll call vote (numeric)}
#'   \item{V623}{Roll call vote (numeric)}
#'   \item{V624}{Roll call vote (numeric)}
#'   \item{V625}{Roll call vote (numeric)}
#'   \item{V626}{Roll call vote (numeric)}
#'   \item{V627}{Roll call vote (numeric)}
#'   \item{V628}{Roll call vote (numeric)}
#'   \item{V629}{Roll call vote (numeric)}
#'   \item{V630}{Roll call vote (numeric)}
#'   \item{V631}{Roll call vote (numeric)}
#'   \item{V632}{Roll call vote (numeric)}
#'   \item{V633}{Roll call vote (numeric)}
#'   \item{V634}{Roll call vote (numeric)}
#'   \item{V635}{Roll call vote (numeric)}
#'   \item{V636}{Roll call vote (numeric)}
#'   \item{V637}{Roll call vote (numeric)}
#'   \item{V638}{Roll call vote (numeric)}
#'   \item{V639}{Roll call vote (numeric)}
#'   \item{V640}{Roll call vote (numeric)}
#'   \item{V641}{Roll call vote (numeric)}
#'   \item{V642}{Roll call vote (numeric)}
#'   \item{V643}{Roll call vote (numeric)}
#'   \item{V644}{Roll call vote (numeric)}
#'   \item{V645}{Roll call vote (numeric)}
#'   \item{V646}{Roll call vote (numeric)}
#'   \item{V647}{Roll call vote (numeric)}
#'   \item{V648}{Roll call vote (numeric)}
#'   \item{V649}{Roll call vote (numeric)}
#'   \item{V650}{Roll call vote (numeric)}
#'   \item{V651}{Roll call vote (numeric)}
#'   \item{V652}{Roll call vote (numeric)}
#'   \item{V653}{Roll call vote (numeric)}
#'   \item{V654}{Roll call vote (numeric)}
#'   \item{V655}{Roll call vote (numeric)}
#'   \item{V656}{Roll call vote (numeric)}
#'   \item{V657}{Roll call vote (numeric)}
#'   \item{V658}{Roll call vote (numeric)}
#'   \item{V659}{Roll call vote (numeric)}
#'   \item{V660}{Roll call vote (numeric)}
#'   \item{V661}{Roll call vote (numeric)}
#'   \item{V662}{Roll call vote (numeric)}
#'   \item{V663}{Roll call vote (numeric)}
#'   \item{V664}{Roll call vote (numeric)}
#'   \item{V665}{Roll call vote (numeric)}
#'   \item{V666}{Roll call vote (numeric)}
#'   \item{V667}{Roll call vote (numeric)}
#'   \item{V668}{Roll call vote (numeric)}
#'   \item{V669}{Roll call vote (numeric)}
#'   \item{V670}{Roll call vote (numeric)}
#'   \item{V671}{Roll call vote (numeric)}
#'   \item{V672}{Roll call vote (numeric)}
#'   \item{V673}{Roll call vote (numeric)}
#'   \item{V674}{Roll call vote (numeric)}
#'   \item{V675}{Roll call vote (numeric)}
#'   \item{V676}{Roll call vote (numeric)}
#'   \item{V677}{Roll call vote (numeric)}
#'   \item{V678}{Roll call vote (numeric)}
#'   \item{V679}{Roll call vote (numeric)}
#'   \item{V680}{Roll call vote (numeric)}
#'   \item{V681}{Roll call vote (numeric)}
#'   \item{V682}{Roll call vote (numeric)}
#'   \item{V683}{Roll call vote (numeric)}
#'   \item{V684}{Roll call vote (numeric)}
#'   \item{V685}{Roll call vote (numeric)}
#'   \item{V686}{Roll call vote (numeric)}
#'   \item{V687}{Roll call vote (numeric)}
#'   \item{V688}{Roll call vote (numeric)}
#'   \item{V689}{Roll call vote (numeric)}
#'   \item{V690}{Roll call vote (numeric)}
#'   \item{V691}{Roll call vote (numeric)}
#'   \item{V692}{Roll call vote (numeric)}
#'   \item{V693}{Roll call vote (numeric)}
#'   \item{V694}{Roll call vote (numeric)}
#'   \item{V695}{Roll call vote (numeric)}
#'   \item{V696}{Roll call vote (numeric)}
#'   \item{V697}{Roll call vote (numeric)}
#'   \item{V698}{Roll call vote (numeric)}
#'   \item{V699}{Roll call vote (numeric)}
#'   \item{V700}{Roll call vote (numeric)}
#'   \item{V701}{Roll call vote (numeric)}
#'   \item{V702}{Roll call vote (numeric)}
#'   \item{V703}{Roll call vote (numeric)}
#'   \item{V704}{Roll call vote (numeric)}
#'   \item{V705}{Roll call vote (numeric)}
#'   \item{V706}{Roll call vote (numeric)}
#'   \item{V707}{Roll call vote (numeric)}
#'   \item{V708}{Roll call vote (numeric)}
#'   \item{V709}{Roll call vote (numeric)}
#'   \item{V710}{Roll call vote (numeric)}
#'   \item{V711}{Roll call vote (numeric)}
#'   \item{V712}{Roll call vote (numeric)}
#'   \item{V713}{Roll call vote (numeric)}
#'   \item{V714}{Roll call vote (numeric)}
#'   \item{V715}{Roll call vote (numeric)}
#'   \item{V716}{Roll call vote (numeric)}
#'   \item{V717}{Roll call vote (numeric)}
#'   \item{V718}{Roll call vote (numeric)}
#'   \item{V719}{Roll call vote (numeric)}
#'   \item{V720}{Roll call vote (numeric)}
#'   \item{V721}{Roll call vote (numeric)}
#'   \item{V722}{Roll call vote (numeric)}
#'   \item{V723}{Roll call vote (numeric)}
#'   \item{V724}{Roll call vote (numeric)}
#'   \item{V725}{Roll call vote (numeric)}
#'   \item{V726}{Roll call vote (numeric)}
#'   \item{V727}{Roll call vote (numeric)}
#'   \item{V728}{Roll call vote (numeric)}
#'   \item{V729}{Roll call vote (numeric)}
#'   \item{V730}{Roll call vote (numeric)}
#'   \item{V731}{Roll call vote (numeric)}
#'   \item{V732}{Roll call vote (numeric)}
#'   \item{V733}{Roll call vote (numeric)}
#'   \item{V734}{Roll call vote (numeric)}
#'   \item{V735}{Roll call vote (numeric)}
#'   \item{V736}{Roll call vote (numeric)}
#'   \item{V737}{Roll call vote (numeric)}
#'   \item{V738}{Roll call vote (numeric)}
#'   \item{V739}{Roll call vote (numeric)}
#'   \item{V740}{Roll call vote (numeric)}
#'   \item{V741}{Roll call vote (numeric)}
#'   \item{V742}{Roll call vote (numeric)}
#'   \item{V743}{Roll call vote (numeric)}
#'   \item{V744}{Roll call vote (numeric)}
#'   \item{V745}{Roll call vote (numeric)}
#'   \item{V746}{Roll call vote (numeric)}
#'   \item{V747}{Roll call vote (numeric)}
#'   \item{V748}{Roll call vote (numeric)}
#'   \item{V749}{Roll call vote (numeric)}
#'   \item{V750}{Roll call vote (numeric)}
#'   \item{V751}{Roll call vote (numeric)}
#'   \item{V752}{Roll call vote (numeric)}
#'   \item{V753}{Roll call vote (numeric)}
#'   \item{V754}{Roll call vote (numeric)}
#'   \item{V755}{Roll call vote (numeric)}
#'   \item{V756}{Roll call vote (numeric)}
#'   \item{V757}{Roll call vote (numeric)}
#'   \item{V758}{Roll call vote (numeric)}
#'   \item{V759}{Roll call vote (numeric)}
#'   \item{V760}{Roll call vote (numeric)}
#'   \item{V761}{Roll call vote (numeric)}
#'   \item{V762}{Roll call vote (numeric)}
#'   \item{V763}{Roll call vote (numeric)}
#'   \item{V764}{Roll call vote (numeric)}
#'   \item{V765}{Roll call vote (numeric)}
#'   \item{V766}{Roll call vote (numeric)}
#'   \item{V767}{Roll call vote (numeric)}
#'   \item{V768}{Roll call vote (numeric)}
#'   \item{V769}{Roll call vote (numeric)}
#'   \item{V770}{Roll call vote (numeric)}
#'   \item{V771}{Roll call vote (numeric)}
#'   \item{V772}{Roll call vote (numeric)}
#'   \item{V773}{Roll call vote (numeric)}
#'   \item{V774}{Roll call vote (numeric)}
#'   \item{V775}{Roll call vote (numeric)}
#'   \item{V776}{Roll call vote (numeric)}
#'   \item{V777}{Roll call vote (numeric)}
#'   \item{V778}{Roll call vote (numeric)}
#'   \item{V779}{Roll call vote (numeric)}
#'   \item{V780}{Roll call vote (numeric)}
#'   \item{V781}{Roll call vote (numeric)}
#'   \item{V782}{Roll call vote (numeric)}
#'   \item{V783}{Roll call vote (numeric)}
#'   \item{V784}{Roll call vote (numeric)}
#'   \item{V785}{Roll call vote (numeric)}
#'   \item{V786}{Roll call vote (numeric)}
#'   \item{V787}{Roll call vote (numeric)}
#'   \item{V788}{Roll call vote (numeric)}
#'   \item{V789}{Roll call vote (numeric)}
#'   \item{V790}{Roll call vote (numeric)}
#'   \item{V791}{Roll call vote (numeric)}
#'   \item{V792}{Roll call vote (numeric)}
#'   \item{V793}{Roll call vote (numeric)}
#'   \item{V794}{Roll call vote (numeric)}
#'   \item{V795}{Roll call vote (numeric)}
#'   \item{V796}{Roll call vote (numeric)}
#'   \item{V797}{Roll call vote (numeric)}
#'   \item{V798}{Roll call vote (numeric)}
#'   \item{V799}{Roll call vote (numeric)}
#'   \item{V800}{Roll call vote (numeric)}
#'   \item{V801}{Roll call vote (numeric)}
#'   \item{V802}{Roll call vote (numeric)}
#'   \item{V803}{Roll call vote (numeric)}
#'   \item{V804}{Roll call vote (numeric)}
#'   \item{V805}{Roll call vote (numeric)}
#'   \item{V806}{Roll call vote (numeric)}
#'   \item{V807}{Roll call vote (numeric)}
#'   \item{V808}{Roll call vote (numeric)}
#'   \item{V809}{Roll call vote (numeric)}
#'   \item{V810}{Roll call vote (numeric)}
#'   \item{V811}{Roll call vote (numeric)}
#'   \item{V812}{Roll call vote (numeric)}
#'   \item{V813}{Roll call vote (numeric)}
#'   \item{V814}{Roll call vote (numeric)}
#'   \item{V815}{Roll call vote (numeric)}
#'   \item{V816}{Roll call vote (numeric)}
#'   \item{V817}{Roll call vote (numeric)}
#'   \item{V818}{Roll call vote (numeric)}
#'   \item{V819}{Roll call vote (numeric)}
#'   \item{V820}{Roll call vote (numeric)}
#'   \item{V821}{Roll call vote (numeric)}
#'   \item{V822}{Roll call vote (numeric)}
#'   \item{V823}{Roll call vote (numeric)}
#'   \item{V824}{Roll call vote (numeric)}
#'   \item{V825}{Roll call vote (numeric)}
#'   \item{V826}{Roll call vote (numeric)}
#'   \item{V827}{Roll call vote (numeric)}
#'   \item{V828}{Roll call vote (numeric)}
#'   \item{V829}{Roll call vote (numeric)}
#'   \item{V830}{Roll call vote (numeric)}
#'   \item{V831}{Roll call vote (numeric)}
#'   \item{V832}{Roll call vote (numeric)}
#'   \item{V833}{Roll call vote (numeric)}
#'   \item{V834}{Roll call vote (numeric)}
#'   \item{V835}{Roll call vote (numeric)}
#'   \item{V836}{Roll call vote (numeric)}
#'   \item{V837}{Roll call vote (numeric)}
#'   \item{V838}{Roll call vote (numeric)}
#'   \item{V839}{Roll call vote (numeric)}
#'   \item{V840}{Roll call vote (numeric)}
#'   \item{V841}{Roll call vote (numeric)}
#'   \item{V842}{Roll call vote (numeric)}
#'   \item{V843}{Roll call vote (numeric)}
#'   \item{V844}{Roll call vote (numeric)}
#'   \item{V845}{Roll call vote (numeric)}
#'   \item{V846}{Roll call vote (numeric)}
#'   \item{V847}{Roll call vote (numeric)}
#'   \item{V848}{Roll call vote (numeric)}
#'   \item{V849}{Roll call vote (numeric)}
#'   \item{V850}{Roll call vote (numeric)}
#'   \item{V851}{Roll call vote (numeric)}
#'   \item{V852}{Roll call vote (numeric)}
#'   \item{V853}{Roll call vote (numeric)}
#'   \item{V854}{Roll call vote (numeric)}
#'   \item{V855}{Roll call vote (numeric)}
#'   \item{V856}{Roll call vote (numeric)}
#'   \item{V857}{Roll call vote (numeric)}
#'   \item{V858}{Roll call vote (numeric)}
#'   \item{V859}{Roll call vote (numeric)}
#'   \item{V860}{Roll call vote (numeric)}
#'   \item{V861}{Roll call vote (numeric)}
#'   \item{V862}{Roll call vote (numeric)}
#'   \item{V863}{Roll call vote (numeric)}
#'   \item{V864}{Roll call vote (numeric)}
#'   \item{V865}{Roll call vote (numeric)}
#'   \item{V866}{Roll call vote (numeric)}
#'   \item{V867}{Roll call vote (numeric)}
#'   \item{V868}{Roll call vote (numeric)}
#'   \item{V869}{Roll call vote (numeric)}
#'   \item{V870}{Roll call vote (numeric)}
#'   \item{V871}{Roll call vote (numeric)}
#'   \item{V872}{Roll call vote (numeric)}
#'   \item{V873}{Roll call vote (numeric)}
#'   \item{V874}{Roll call vote (numeric)}
#'   \item{V875}{Roll call vote (numeric)}
#'   \item{V876}{Roll call vote (numeric)}
#'   \item{V877}{Roll call vote (numeric)}
#'   \item{V878}{Roll call vote (numeric)}
#'   \item{V879}{Roll call vote (numeric)}
#'   \item{V880}{Roll call vote (numeric)}
#'   \item{V881}{Roll call vote (numeric)}
#'   \item{V882}{Roll call vote (numeric)}
#'   \item{V883}{Roll call vote (numeric)}
#'   \item{V884}{Roll call vote (numeric)}
#'   \item{V885}{Roll call vote (numeric)}
#'   \item{V886}{Roll call vote (numeric)}
#' }
#'
#' @details
#' This dataset covers the first directly elected European Parliament 
#' (1979-1984), a formative period for EU legislative politics. The first 
#' five columns contain MEP identification and affiliation variables; the 
#' remaining 886 columns represent roll call votes.
#' 
#' Each row represents one MEP, and each vote column (V1 through V886) 
#' represents a specific legislative vote. The dataset is useful for analyzing 
#' transnational party behavior, national vs. European party loyalty, and 
#' coalition formation in the early European Parliament.
#'
#' @source
#' Hix, S., Noury, A., & Roland, G. (2006). Dimensions of Politics in the 
#' European Parliament. \emph{American Journal of Political Science}, 50(2), 
#' 494-520. \doi{10.1111/j.1540-5907.2006.00198.x}
#' 
#' Original data: \url{http://personal.lse.ac.uk/hix/HixNouryRolandEPdata.htm}
#'
#' @references
#' Hix, S., Noury, A., & Roland, G. (2006). Dimensions of Politics in the 
#' European Parliament. \emph{American Journal of Political Science}, 50(2), 
#' 494-520.
#'
#' @usage data(rcv_ep1)
#'
#' @examples
#' \dontrun{
#' data(rcv_ep1)
#' 
#' # Dataset dimensions
#' dim(rcv_ep1)  # 410 MEPs, 891 columns
#' 
#' # View MEP information
#' head(rcv_ep1[, 1:5])
#' 
#' # Distribution by country
#' table(rcv_ep1$MS)
#' 
#' # Distribution by EP group
#' table(rcv_ep1$EPG)
#' 
#' # View first few votes
#' head(rcv_ep1[, 6:10])
#' 
#' # Check voting patterns
#' summary(rcv_ep1[, 6:891])
#' }
#'
#' @keywords datasets
#' @name rcv_ep1
#' @docType data
NULL


#' @encoding UTF-8
#' @title Roll Call Data from the 111th U.S. Senate
#' @description The `hr111` dataset contains roll call voting data from the 111th U.S. Senate. This dataset is formatted as a `rollcall` object, which is typically used for analyzing voting behavior in legislative bodies.
#'
#' @format An object of class `rollcall` with the following components:
#' \describe{
#'   \item{votes}{A matrix or data frame where rows represent senators and columns represent individual roll call votes. Each cell contains the vote cast by a senator.}
#'   \item{legis.names}{A vector containing the names of the senators.}
#'   \item{party}{A vector indicating the party affiliation of each senator.}
#'   \item{desc}{A brief description of the roll call data (e.g., "111th U.S. Senate").}
#'   \item{codes}{A list of codes used to represent different types of votes (e.g., Yea, Nay, Abstain).}
#'   \item{source}{A description of the source of the data.}
#' }
#'
#' @details
#' The `hr111` dataset is used to analyze the voting patterns and behavior of senators during the 111th U.S. Senate session. 
#' This dataset is stored as a `rollcall` object, which is a common format for legislative voting data analysis in R.
#'
#' The `votes` component contains the actual voting records, where each row corresponds to a senator and each column to a specific roll call vote. 
#' The `legis.names` and `party` components provide additional context, including the names of the senators and their party affiliations.
#'
#' @usage data(hr111)
#'
#' @examples
#' \dontrun{
#' data(hr111)
#' }
#'
#' @keywords datasets
#' @name hr111
#' @docType data
NULL

#' @encoding UTF-8
#' @title Nation Similarity Ratings Dataset
#' @description The `nation` dataset contains similarity ratings between twelve nations, as collected by Wish (1971). 
#' In 1968, Wish asked 18 students in his psychological measurement class to rate the perceived similarity between each pair of twelve nations using a 9-point scale, 
#' where '1' indicates "very different" and '9' indicates "very similar". The data in this dataset represent the average similarity ratings between these nations.
#'
#' @format A 12x12 matrix with the following column and row names representing nations:
#' \describe{
#'   \item{Brazil}{Average similarity ratings involving Brazil.}
#'   \item{Congo}{Average similarity ratings involving Congo.}
#'   \item{Cuba}{Average similarity ratings involving Cuba.}
#'   \item{Egypt}{Average similarity ratings involving Egypt.}
#'   \item{France}{Average similarity ratings involving France.}
#'   \item{India}{Average similarity ratings involving India.}
#'   \item{Israel}{Average similarity ratings involving Israel.}
#'   \item{Japan}{Average similarity ratings involving Japan.}
#'   \item{China}{Average similarity ratings involving China.}
#'   \item{USSR}{Average similarity ratings involving the Soviet Union (USSR).}
#'   \item{USA}{Average similarity ratings involving the United States of America (USA).}
#'   \item{Yugoslavia}{Average similarity ratings involving Yugoslavia.}
#' }
#'
#' @details
#' The `nations` dataset was constructed by averaging the similarity ratings provided by 18 students in Wish's 1968 psychological measurement class. 
#' The dataset is a symmetric matrix where both the rows and columns represent the twelve nations, and each cell contains the average similarity rating 
#' between the corresponding pair of nations.
#'
#' @source
#' Wish, M. (1971). Individual differences in perceptions and preferences among nations. *In Proceedings of the 79th Annual Convention of the American Psychological Association*.
#'
#' @usage data(nations)
#'
#' @examples
#' \dontrun{
#' data(nations)
#' }
#'
#' @keywords datasets
#' @name nations
#' @docType data
NULL

#' @encoding UTF-8
#' @title DW-NOMINATE Scores for the U.S. Congress
#' @description The `rcx` dataset is a matrix containing DW-NOMINATE scores for members of the U.S. Congress. DW-NOMINATE scores are used 
#' to measure the ideological positions of legislators based on their roll-call voting behavior.
#'
#' @format A data frame with the following 16 variables:
#' \describe{
#'   \item{cong}{Congress number.}
#'   \item{id}{Unique identifier for each member of Congress.}
#'   \item{state}{State abbreviation for the member's state.}
#'   \item{dist}{District number for House members, or 0 for Senators.}
#'   \item{lstate}{Full name of the member's state.}
#'   \item{party}{Political party affiliation (e.g., Democrat, Republican).}
#'   \item{name}{Name of the member of Congress.}
#'   \item{dwnom1}{The first dimension DW-NOMINATE score, typically representing the liberal-conservative spectrum.}
#'   \item{dwnom2}{The second dimension DW-NOMINATE score, often representing regional or other secondary factors.}
#'   \item{dwnom1bse}{Standard error of the first dimension DW-NOMINATE score.}
#'   \item{dwnom2bse}{Standard error of the second dimension DW-NOMINATE score.}
#'   \item{corrbse}{Correlation between the first and second dimension standard errors.}
#'   \item{LogL}{Log-likelihood of the legislator's votes under the DW-NOMINATE model.}
#'   \item{nchoice}{Number of roll-call votes the legislator participated in.}
#'   \item{nerror}{Number of classification errors (votes incorrectly predicted by the model).}
#'   \item{gmp}{Geometric mean probability, a measure of fit for the DW-NOMINATE model.}
#' }
#'
#' @details
#' The `rcx` dataset provides DW-NOMINATE scores for members of the U.S. Congress, which are commonly used to analyze 
#' ideological alignment and voting behavior across different legislative sessions.
#'
#' @source
#' Data collected and processed using the DW-NOMINATE algorithm. For more details, see the [DW-NOMINATE project website](https://voteview.com/dwnominate).
#'
#' @usage data(rcx)
#'
#' @examples
#' \dontrun{
#' data(rcx)
#' }
#'
#' @keywords datasets
#' @name rcx
#' @docType data
NULL

#' @encoding UTF-8
#' @title French Party Placement Data from the 2009 European Election Study (EES)
#' @description  The `french.parties.individuals` dataset contains party placement data from the French module of the 2009 European Election Study (EES). 
#' This dataset was used in Chapter 2 to directly scale respondents’ placements of eight major political parties on a ten-point left-right ideological scale.
#' 
#' @format A matrix with 1,000 rows (representing respondents) and 8 columns (representing political parties):
#' \describe{
#'   \item{extremeleft}{Placement of the "Extreme Left" party.}
#'   \item{communist}{Placement of the Communist Party.}
#'   \item{socialist}{Placement of the Socialist Party.}
#'   \item{greens}{Placement of the Greens Party.}
#'   \item{udfbayrou}{Placement of the Union for French Democracy (Bayrou's party).}
#'   \item{umpsarkozy}{Placement of the Union for a Popular Movement (Sarkozy's party).}
#'   \item{nationalfront}{Placement of the National Front party.}
#'   \item{leftparty}{Placement of the Left Party.}
#' }
#'
#' @details
#' The French module of the 2009 EES asked 1,000 respondents to place eight major political parties on a ten-point left-right ideological scale. 
#' This dataset, `french.parties.individuals`, contains those placements. In Chapter 2, these placements were scaled directly to compare the results 
#' obtained from scaling different forms of data from the same set of observations.
#'
#' In addition to directly scaling respondents' placements, each set of placements can be arranged into a similarities matrix for further analysis.
#'
#' @source
#' Data from the French module of the 2009 European Election Study (EES).
#'
#' @usage data(french.parties.individuals)
#'
#' @examples
#' \dontrun{
#' data(french.parties.individuals)
#' }
#'
#' @keywords datasets
#' @name french.parties.individuals
#' @docType data
NULL

#' @encoding UTF-8
#' @title 7th Legislative Yuan Roll Call Data from Taiwan (2008-2012)
#' @description This dataset contains roll call voting data from the 7th Legislative Yuan 
#' (National Congress) of Taiwan. The dataset includes the names of legislators 
#' and their corresponding votes on various bills.
#'
#' @format A data frame with 113 legislators and 1228 variables:
#' \describe{
#'   \item{party}{Character. The political party of each legislator}
#'   \item{7-1}{Roll call vote (numeric)}
#'   \item{7-10}{Roll call vote (numeric)}
#'   \item{7-100}{Roll call vote (numeric)}
#'   \item{7-1000}{Roll call vote (numeric)}
#'   \item{7-1001}{Roll call vote (numeric)}
#'   \item{7-1002}{Roll call vote (numeric)}
#'   \item{7-1003}{Roll call vote (numeric)}
#'   \item{7-1004}{Roll call vote (numeric)}
#'   \item{7-1005}{Roll call vote (numeric)}
#'   \item{7-1006}{Roll call vote (numeric)}
#'   \item{7-1007}{Roll call vote (numeric)}
#'   \item{7-1008}{Roll call vote (numeric)}
#'   \item{7-1009}{Roll call vote (numeric)}
#'   \item{7-101}{Roll call vote (numeric)}
#'   \item{7-1010}{Roll call vote (numeric)}
#'   \item{7-1011}{Roll call vote (numeric)}
#'   \item{7-1012}{Roll call vote (numeric)}
#'   \item{7-1013}{Roll call vote (numeric)}
#'   \item{7-1014}{Roll call vote (numeric)}
#'   \item{7-1015}{Roll call vote (numeric)}
#'   \item{7-1016}{Roll call vote (numeric)}
#'   \item{7-1017}{Roll call vote (numeric)}
#'   \item{7-1018}{Roll call vote (numeric)}
#'   \item{7-1019}{Roll call vote (numeric)}
#'   \item{7-102}{Roll call vote (numeric)}
#'   \item{7-1020}{Roll call vote (numeric)}
#'   \item{7-1021}{Roll call vote (numeric)}
#'   \item{7-1022}{Roll call vote (numeric)}
#'   \item{7-1023}{Roll call vote (numeric)}
#'   \item{7-1024}{Roll call vote (numeric)}
#'   \item{7-1025}{Roll call vote (numeric)}
#'   \item{7-1026}{Roll call vote (numeric)}
#'   \item{7-1027}{Roll call vote (numeric)}
#'   \item{7-1028}{Roll call vote (numeric)}
#'   \item{7-1029}{Roll call vote (numeric)}
#'   \item{7-103}{Roll call vote (numeric)}
#'   \item{7-1030}{Roll call vote (numeric)}
#'   \item{7-1031}{Roll call vote (numeric)}
#'   \item{7-1032}{Roll call vote (numeric)}
#'   \item{7-1033}{Roll call vote (numeric)}
#'   \item{7-1034}{Roll call vote (numeric)}
#'   \item{7-1035}{Roll call vote (numeric)}
#'   \item{7-1036}{Roll call vote (numeric)}
#'   \item{7-1037}{Roll call vote (numeric)}
#'   \item{7-1038}{Roll call vote (numeric)}
#'   \item{7-1039}{Roll call vote (numeric)}
#'   \item{7-104}{Roll call vote (numeric)}
#'   \item{7-1040}{Roll call vote (numeric)}
#'   \item{7-1041}{Roll call vote (numeric)}
#'   \item{7-1042}{Roll call vote (numeric)}
#'   \item{7-1043}{Roll call vote (numeric)}
#'   \item{7-1044}{Roll call vote (numeric)}
#'   \item{7-1045}{Roll call vote (numeric)}
#'   \item{7-1046}{Roll call vote (numeric)}
#'   \item{7-1047}{Roll call vote (numeric)}
#'   \item{7-1048}{Roll call vote (numeric)}
#'   \item{7-1049}{Roll call vote (numeric)}
#'   \item{7-105}{Roll call vote (numeric)}
#'   \item{7-1050}{Roll call vote (numeric)}
#'   \item{7-1051}{Roll call vote (numeric)}
#'   \item{7-1052}{Roll call vote (numeric)}
#'   \item{7-1053}{Roll call vote (numeric)}
#'   \item{7-1054}{Roll call vote (numeric)}
#'   \item{7-1055}{Roll call vote (numeric)}
#'   \item{7-1056}{Roll call vote (numeric)}
#'   \item{7-1057}{Roll call vote (numeric)}
#'   \item{7-1058}{Roll call vote (numeric)}
#'   \item{7-1059}{Roll call vote (numeric)}
#'   \item{7-106}{Roll call vote (numeric)}
#'   \item{7-1060}{Roll call vote (numeric)}
#'   \item{7-1061}{Roll call vote (numeric)}
#'   \item{7-1062}{Roll call vote (numeric)}
#'   \item{7-1063}{Roll call vote (numeric)}
#'   \item{7-1064}{Roll call vote (numeric)}
#'   \item{7-1065}{Roll call vote (numeric)}
#'   \item{7-1066}{Roll call vote (numeric)}
#'   \item{7-1067}{Roll call vote (numeric)}
#'   \item{7-1068}{Roll call vote (numeric)}
#'   \item{7-1069}{Roll call vote (numeric)}
#'   \item{7-107}{Roll call vote (numeric)}
#'   \item{7-1070}{Roll call vote (numeric)}
#'   \item{7-1071}{Roll call vote (numeric)}
#'   \item{7-1072}{Roll call vote (numeric)}
#'   \item{7-1073}{Roll call vote (numeric)}
#'   \item{7-1074}{Roll call vote (numeric)}
#'   \item{7-1075}{Roll call vote (numeric)}
#'   \item{7-1076}{Roll call vote (numeric)}
#'   \item{7-1077}{Roll call vote (numeric)}
#'   \item{7-1078}{Roll call vote (numeric)}
#'   \item{7-1079}{Roll call vote (numeric)}
#'   \item{7-108}{Roll call vote (numeric)}
#'   \item{7-1080}{Roll call vote (numeric)}
#'   \item{7-1081}{Roll call vote (numeric)}
#'   \item{7-1082}{Roll call vote (numeric)}
#'   \item{7-1083}{Roll call vote (numeric)}
#'   \item{7-1084}{Roll call vote (numeric)}
#'   \item{7-1085}{Roll call vote (numeric)}
#'   \item{7-1086}{Roll call vote (numeric)}
#'   \item{7-1087}{Roll call vote (numeric)}
#'   \item{7-1088}{Roll call vote (numeric)}
#'   \item{7-1089}{Roll call vote (numeric)}
#'   \item{7-109}{Roll call vote (numeric)}
#'   \item{7-1090}{Roll call vote (numeric)}
#'   \item{7-1091}{Roll call vote (numeric)}
#'   \item{7-1092}{Roll call vote (numeric)}
#'   \item{7-1093}{Roll call vote (numeric)}
#'   \item{7-1094}{Roll call vote (numeric)}
#'   \item{7-1095}{Roll call vote (numeric)}
#'   \item{7-1096}{Roll call vote (numeric)}
#'   \item{7-1097}{Roll call vote (numeric)}
#'   \item{7-1098}{Roll call vote (numeric)}
#'   \item{7-1099}{Roll call vote (numeric)}
#'   \item{7-11}{Roll call vote (numeric)}
#'   \item{7-110}{Roll call vote (numeric)}
#'   \item{7-1100}{Roll call vote (numeric)}
#'   \item{7-1101}{Roll call vote (numeric)}
#'   \item{7-1102}{Roll call vote (numeric)}
#'   \item{7-1103}{Roll call vote (numeric)}
#'   \item{7-1104}{Roll call vote (numeric)}
#'   \item{7-1105}{Roll call vote (numeric)}
#'   \item{7-1106}{Roll call vote (numeric)}
#'   \item{7-1107}{Roll call vote (numeric)}
#'   \item{7-1108}{Roll call vote (numeric)}
#'   \item{7-1109}{Roll call vote (numeric)}
#'   \item{7-111}{Roll call vote (numeric)}
#'   \item{7-1110}{Roll call vote (numeric)}
#'   \item{7-1111}{Roll call vote (numeric)}
#'   \item{7-1112}{Roll call vote (numeric)}
#'   \item{7-1113}{Roll call vote (numeric)}
#'   \item{7-1114}{Roll call vote (numeric)}
#'   \item{7-1115}{Roll call vote (numeric)}
#'   \item{7-1116}{Roll call vote (numeric)}
#'   \item{7-1117}{Roll call vote (numeric)}
#'   \item{7-1118}{Roll call vote (numeric)}
#'   \item{7-1119}{Roll call vote (numeric)}
#'   \item{7-112}{Roll call vote (numeric)}
#'   \item{7-1120}{Roll call vote (numeric)}
#'   \item{7-1121}{Roll call vote (numeric)}
#'   \item{7-1122}{Roll call vote (numeric)}
#'   \item{7-1123}{Roll call vote (numeric)}
#'   \item{7-1124}{Roll call vote (numeric)}
#'   \item{7-1125}{Roll call vote (numeric)}
#'   \item{7-1126}{Roll call vote (numeric)}
#'   \item{7-1127}{Roll call vote (numeric)}
#'   \item{7-1128}{Roll call vote (numeric)}
#'   \item{7-1129}{Roll call vote (numeric)}
#'   \item{7-113}{Roll call vote (numeric)}
#'   \item{7-1130}{Roll call vote (numeric)}
#'   \item{7-1131}{Roll call vote (numeric)}
#'   \item{7-1132}{Roll call vote (numeric)}
#'   \item{7-1133}{Roll call vote (numeric)}
#'   \item{7-1134}{Roll call vote (numeric)}
#'   \item{7-1135}{Roll call vote (numeric)}
#'   \item{7-1136}{Roll call vote (numeric)}
#'   \item{7-1137}{Roll call vote (numeric)}
#'   \item{7-1138}{Roll call vote (numeric)}
#'   \item{7-1139}{Roll call vote (numeric)}
#'   \item{7-114}{Roll call vote (numeric)}
#'   \item{7-1140}{Roll call vote (numeric)}
#'   \item{7-1141}{Roll call vote (numeric)}
#'   \item{7-1142}{Roll call vote (numeric)}
#'   \item{7-1143}{Roll call vote (numeric)}
#'   \item{7-1144}{Roll call vote (numeric)}
#'   \item{7-1145}{Roll call vote (numeric)}
#'   \item{7-1146}{Roll call vote (numeric)}
#'   \item{7-1147}{Roll call vote (numeric)}
#'   \item{7-1148}{Roll call vote (numeric)}
#'   \item{7-1149}{Roll call vote (numeric)}
#'   \item{7-115}{Roll call vote (numeric)}
#'   \item{7-1150}{Roll call vote (numeric)}
#'   \item{7-1151}{Roll call vote (numeric)}
#'   \item{7-1152}{Roll call vote (numeric)}
#'   \item{7-1153}{Roll call vote (numeric)}
#'   \item{7-1154}{Roll call vote (numeric)}
#'   \item{7-1155}{Roll call vote (numeric)}
#'   \item{7-1156}{Roll call vote (numeric)}
#'   \item{7-1157}{Roll call vote (numeric)}
#'   \item{7-1158}{Roll call vote (numeric)}
#'   \item{7-1159}{Roll call vote (numeric)}
#'   \item{7-116}{Roll call vote (numeric)}
#'   \item{7-1160}{Roll call vote (numeric)}
#'   \item{7-1161}{Roll call vote (numeric)}
#'   \item{7-1162}{Roll call vote (numeric)}
#'   \item{7-1163}{Roll call vote (numeric)}
#'   \item{7-1164}{Roll call vote (numeric)}
#'   \item{7-1165}{Roll call vote (numeric)}
#'   \item{7-1166}{Roll call vote (numeric)}
#'   \item{7-1167}{Roll call vote (numeric)}
#'   \item{7-1168}{Roll call vote (numeric)}
#'   \item{7-1169}{Roll call vote (numeric)}
#'   \item{7-117}{Roll call vote (numeric)}
#'   \item{7-1170}{Roll call vote (numeric)}
#'   \item{7-1171}{Roll call vote (numeric)}
#'   \item{7-1172}{Roll call vote (numeric)}
#'   \item{7-1173}{Roll call vote (numeric)}
#'   \item{7-1174}{Roll call vote (numeric)}
#'   \item{7-1175}{Roll call vote (numeric)}
#'   \item{7-1176}{Roll call vote (numeric)}
#'   \item{7-1177}{Roll call vote (numeric)}
#'   \item{7-1178}{Roll call vote (numeric)}
#'   \item{7-1179}{Roll call vote (numeric)}
#'   \item{7-118}{Roll call vote (numeric)}
#'   \item{7-1180}{Roll call vote (numeric)}
#'   \item{7-1181}{Roll call vote (numeric)}
#'   \item{7-1182}{Roll call vote (numeric)}
#'   \item{7-1183}{Roll call vote (numeric)}
#'   \item{7-1184}{Roll call vote (numeric)}
#'   \item{7-1185}{Roll call vote (numeric)}
#'   \item{7-1186}{Roll call vote (numeric)}
#'   \item{7-1187}{Roll call vote (numeric)}
#'   \item{7-1188}{Roll call vote (numeric)}
#'   \item{7-1189}{Roll call vote (numeric)}
#'   \item{7-119}{Roll call vote (numeric)}
#'   \item{7-1190}{Roll call vote (numeric)}
#'   \item{7-1191}{Roll call vote (numeric)}
#'   \item{7-1192}{Roll call vote (numeric)}
#'   \item{7-1193}{Roll call vote (numeric)}
#'   \item{7-1194}{Roll call vote (numeric)}
#'   \item{7-1195}{Roll call vote (numeric)}
#'   \item{7-1196}{Roll call vote (numeric)}
#'   \item{7-1197}{Roll call vote (numeric)}
#'   \item{7-1198}{Roll call vote (numeric)}
#'   \item{7-1199}{Roll call vote (numeric)}
#'   \item{7-12}{Roll call vote (numeric)}
#'   \item{7-120}{Roll call vote (numeric)}
#'   \item{7-1200}{Roll call vote (numeric)}
#'   \item{7-1201}{Roll call vote (numeric)}
#'   \item{7-1202}{Roll call vote (numeric)}
#'   \item{7-1203}{Roll call vote (numeric)}
#'   \item{7-1204}{Roll call vote (numeric)}
#'   \item{7-1205}{Roll call vote (numeric)}
#'   \item{7-1206}{Roll call vote (numeric)}
#'   \item{7-1207}{Roll call vote (numeric)}
#'   \item{7-1208}{Roll call vote (numeric)}
#'   \item{7-1209}{Roll call vote (numeric)}
#'   \item{7-121}{Roll call vote (numeric)}
#'   \item{7-1210}{Roll call vote (numeric)}
#'   \item{7-1211}{Roll call vote (numeric)}
#'   \item{7-1212}{Roll call vote (numeric)}
#'   \item{7-1213}{Roll call vote (numeric)}
#'   \item{7-1214}{Roll call vote (numeric)}
#'   \item{7-1215}{Roll call vote (numeric)}
#'   \item{7-1216}{Roll call vote (numeric)}
#'   \item{7-1217}{Roll call vote (numeric)}
#'   \item{7-1218}{Roll call vote (numeric)}
#'   \item{7-1219}{Roll call vote (numeric)}
#'   \item{7-122}{Roll call vote (numeric)}
#'   \item{7-1220}{Roll call vote (numeric)}
#'   \item{7-1221}{Roll call vote (numeric)}
#'   \item{7-1223}{Roll call vote (numeric)}
#'   \item{7-1224}{Roll call vote (numeric)}
#'   \item{7-1225}{Roll call vote (numeric)}
#'   \item{7-1226}{Roll call vote (numeric)}
#'   \item{7-123}{Roll call vote (numeric)}
#'   \item{7-124}{Roll call vote (numeric)}
#'   \item{7-125}{Roll call vote (numeric)}
#'   \item{7-126}{Roll call vote (numeric)}
#'   \item{7-127}{Roll call vote (numeric)}
#'   \item{7-128}{Roll call vote (numeric)}
#'   \item{7-129}{Roll call vote (numeric)}
#'   \item{7-13}{Roll call vote (numeric)}
#'   \item{7-130}{Roll call vote (numeric)}
#'   \item{7-131}{Roll call vote (numeric)}
#'   \item{7-132}{Roll call vote (numeric)}
#'   \item{7-133}{Roll call vote (numeric)}
#'   \item{7-134}{Roll call vote (numeric)}
#'   \item{7-135}{Roll call vote (numeric)}
#'   \item{7-136}{Roll call vote (numeric)}
#'   \item{7-137}{Roll call vote (numeric)}
#'   \item{7-138}{Roll call vote (numeric)}
#'   \item{7-139}{Roll call vote (numeric)}
#'   \item{7-14}{Roll call vote (numeric)}
#'   \item{7-140}{Roll call vote (numeric)}
#'   \item{7-141}{Roll call vote (numeric)}
#'   \item{7-142}{Roll call vote (numeric)}
#'   \item{7-143}{Roll call vote (numeric)}
#'   \item{7-144}{Roll call vote (numeric)}
#'   \item{7-145}{Roll call vote (numeric)}
#'   \item{7-146}{Roll call vote (numeric)}
#'   \item{7-147}{Roll call vote (numeric)}
#'   \item{7-148}{Roll call vote (numeric)}
#'   \item{7-149}{Roll call vote (numeric)}
#'   \item{7-15}{Roll call vote (numeric)}
#'   \item{7-150}{Roll call vote (numeric)}
#'   \item{7-151}{Roll call vote (numeric)}
#'   \item{7-152}{Roll call vote (numeric)}
#'   \item{7-153}{Roll call vote (numeric)}
#'   \item{7-154}{Roll call vote (numeric)}
#'   \item{7-155}{Roll call vote (numeric)}
#'   \item{7-156}{Roll call vote (numeric)}
#'   \item{7-157}{Roll call vote (numeric)}
#'   \item{7-158}{Roll call vote (numeric)}
#'   \item{7-159}{Roll call vote (numeric)}
#'   \item{7-16}{Roll call vote (numeric)}
#'   \item{7-160}{Roll call vote (numeric)}
#'   \item{7-161}{Roll call vote (numeric)}
#'   \item{7-162}{Roll call vote (numeric)}
#'   \item{7-163}{Roll call vote (numeric)}
#'   \item{7-164}{Roll call vote (numeric)}
#'   \item{7-165}{Roll call vote (numeric)}
#'   \item{7-166}{Roll call vote (numeric)}
#'   \item{7-167}{Roll call vote (numeric)}
#'   \item{7-168}{Roll call vote (numeric)}
#'   \item{7-169}{Roll call vote (numeric)}
#'   \item{7-17}{Roll call vote (numeric)}
#'   \item{7-170}{Roll call vote (numeric)}
#'   \item{7-171}{Roll call vote (numeric)}
#'   \item{7-172}{Roll call vote (numeric)}
#'   \item{7-173}{Roll call vote (numeric)}
#'   \item{7-174}{Roll call vote (numeric)}
#'   \item{7-175}{Roll call vote (numeric)}
#'   \item{7-176}{Roll call vote (numeric)}
#'   \item{7-177}{Roll call vote (numeric)}
#'   \item{7-178}{Roll call vote (numeric)}
#'   \item{7-179}{Roll call vote (numeric)}
#'   \item{7-18}{Roll call vote (numeric)}
#'   \item{7-180}{Roll call vote (numeric)}
#'   \item{7-181}{Roll call vote (numeric)}
#'   \item{7-182}{Roll call vote (numeric)}
#'   \item{7-183}{Roll call vote (numeric)}
#'   \item{7-184}{Roll call vote (numeric)}
#'   \item{7-185}{Roll call vote (numeric)}
#'   \item{7-186}{Roll call vote (numeric)}
#'   \item{7-187}{Roll call vote (numeric)}
#'   \item{7-188}{Roll call vote (numeric)}
#'   \item{7-189}{Roll call vote (numeric)}
#'   \item{7-19}{Roll call vote (numeric)}
#'   \item{7-190}{Roll call vote (numeric)}
#'   \item{7-191}{Roll call vote (numeric)}
#'   \item{7-192}{Roll call vote (numeric)}
#'   \item{7-193}{Roll call vote (numeric)}
#'   \item{7-194}{Roll call vote (numeric)}
#'   \item{7-195}{Roll call vote (numeric)}
#'   \item{7-196}{Roll call vote (numeric)}
#'   \item{7-197}{Roll call vote (numeric)}
#'   \item{7-198}{Roll call vote (numeric)}
#'   \item{7-199}{Roll call vote (numeric)}
#'   \item{7-2}{Roll call vote (numeric)}
#'   \item{7-20}{Roll call vote (numeric)}
#'   \item{7-200}{Roll call vote (numeric)}
#'   \item{7-201}{Roll call vote (numeric)}
#'   \item{7-202}{Roll call vote (numeric)}
#'   \item{7-203}{Roll call vote (numeric)}
#'   \item{7-204}{Roll call vote (numeric)}
#'   \item{7-205}{Roll call vote (numeric)}
#'   \item{7-206}{Roll call vote (numeric)}
#'   \item{7-207}{Roll call vote (numeric)}
#'   \item{7-208}{Roll call vote (numeric)}
#'   \item{7-209}{Roll call vote (numeric)}
#'   \item{7-21}{Roll call vote (numeric)}
#'   \item{7-210}{Roll call vote (numeric)}
#'   \item{7-211}{Roll call vote (numeric)}
#'   \item{7-212}{Roll call vote (numeric)}
#'   \item{7-213}{Roll call vote (numeric)}
#'   \item{7-214}{Roll call vote (numeric)}
#'   \item{7-215}{Roll call vote (numeric)}
#'   \item{7-216}{Roll call vote (numeric)}
#'   \item{7-217}{Roll call vote (numeric)}
#'   \item{7-218}{Roll call vote (numeric)}
#'   \item{7-219}{Roll call vote (numeric)}
#'   \item{7-22}{Roll call vote (numeric)}
#'   \item{7-220}{Roll call vote (numeric)}
#'   \item{7-221}{Roll call vote (numeric)}
#'   \item{7-222}{Roll call vote (numeric)}
#'   \item{7-223}{Roll call vote (numeric)}
#'   \item{7-224}{Roll call vote (numeric)}
#'   \item{7-225}{Roll call vote (numeric)}
#'   \item{7-226}{Roll call vote (numeric)}
#'   \item{7-227}{Roll call vote (numeric)}
#'   \item{7-228}{Roll call vote (numeric)}
#'   \item{7-229}{Roll call vote (numeric)}
#'   \item{7-23}{Roll call vote (numeric)}
#'   \item{7-230}{Roll call vote (numeric)}
#'   \item{7-231}{Roll call vote (numeric)}
#'   \item{7-232}{Roll call vote (numeric)}
#'   \item{7-233}{Roll call vote (numeric)}
#'   \item{7-234}{Roll call vote (numeric)}
#'   \item{7-235}{Roll call vote (numeric)}
#'   \item{7-236}{Roll call vote (numeric)}
#'   \item{7-237}{Roll call vote (numeric)}
#'   \item{7-238}{Roll call vote (numeric)}
#'   \item{7-239}{Roll call vote (numeric)}
#'   \item{7-24}{Roll call vote (numeric)}
#'   \item{7-240}{Roll call vote (numeric)}
#'   \item{7-241}{Roll call vote (numeric)}
#'   \item{7-242}{Roll call vote (numeric)}
#'   \item{7-243}{Roll call vote (numeric)}
#'   \item{7-244}{Roll call vote (numeric)}
#'   \item{7-245}{Roll call vote (numeric)}
#'   \item{7-246}{Roll call vote (numeric)}
#'   \item{7-247}{Roll call vote (numeric)}
#'   \item{7-248}{Roll call vote (numeric)}
#'   \item{7-249}{Roll call vote (numeric)}
#'   \item{7-25}{Roll call vote (numeric)}
#'   \item{7-250}{Roll call vote (numeric)}
#'   \item{7-251}{Roll call vote (numeric)}
#'   \item{7-252}{Roll call vote (numeric)}
#'   \item{7-253}{Roll call vote (numeric)}
#'   \item{7-254}{Roll call vote (numeric)}
#'   \item{7-255}{Roll call vote (numeric)}
#'   \item{7-256}{Roll call vote (numeric)}
#'   \item{7-257}{Roll call vote (numeric)}
#'   \item{7-258}{Roll call vote (numeric)}
#'   \item{7-259}{Roll call vote (numeric)}
#'   \item{7-26}{Roll call vote (numeric)}
#'   \item{7-260}{Roll call vote (numeric)}
#'   \item{7-261}{Roll call vote (numeric)}
#'   \item{7-262}{Roll call vote (numeric)}
#'   \item{7-263}{Roll call vote (numeric)}
#'   \item{7-264}{Roll call vote (numeric)}
#'   \item{7-265}{Roll call vote (numeric)}
#'   \item{7-266}{Roll call vote (numeric)}
#'   \item{7-267}{Roll call vote (numeric)}
#'   \item{7-268}{Roll call vote (numeric)}
#'   \item{7-269}{Roll call vote (numeric)}
#'   \item{7-27}{Roll call vote (numeric)}
#'   \item{7-270}{Roll call vote (numeric)}
#'   \item{7-271}{Roll call vote (numeric)}
#'   \item{7-272}{Roll call vote (numeric)}
#'   \item{7-273}{Roll call vote (numeric)}
#'   \item{7-274}{Roll call vote (numeric)}
#'   \item{7-275}{Roll call vote (numeric)}
#'   \item{7-276}{Roll call vote (numeric)}
#'   \item{7-277}{Roll call vote (numeric)}
#'   \item{7-278}{Roll call vote (numeric)}
#'   \item{7-279}{Roll call vote (numeric)}
#'   \item{7-28}{Roll call vote (numeric)}
#'   \item{7-280}{Roll call vote (numeric)}
#'   \item{7-281}{Roll call vote (numeric)}
#'   \item{7-282}{Roll call vote (numeric)}
#'   \item{7-283}{Roll call vote (numeric)}
#'   \item{7-284}{Roll call vote (numeric)}
#'   \item{7-285}{Roll call vote (numeric)}
#'   \item{7-286}{Roll call vote (numeric)}
#'   \item{7-287}{Roll call vote (numeric)}
#'   \item{7-288}{Roll call vote (numeric)}
#'   \item{7-289}{Roll call vote (numeric)}
#'   \item{7-29}{Roll call vote (numeric)}
#'   \item{7-290}{Roll call vote (numeric)}
#'   \item{7-291}{Roll call vote (numeric)}
#'   \item{7-292}{Roll call vote (numeric)}
#'   \item{7-293}{Roll call vote (numeric)}
#'   \item{7-294}{Roll call vote (numeric)}
#'   \item{7-295}{Roll call vote (numeric)}
#'   \item{7-296}{Roll call vote (numeric)}
#'   \item{7-297}{Roll call vote (numeric)}
#'   \item{7-298}{Roll call vote (numeric)}
#'   \item{7-299}{Roll call vote (numeric)}
#'   \item{7-3}{Roll call vote (numeric)}
#'   \item{7-30}{Roll call vote (numeric)}
#'   \item{7-300}{Roll call vote (numeric)}
#'   \item{7-301}{Roll call vote (numeric)}
#'   \item{7-302}{Roll call vote (numeric)}
#'   \item{7-303}{Roll call vote (numeric)}
#'   \item{7-304}{Roll call vote (numeric)}
#'   \item{7-305}{Roll call vote (numeric)}
#'   \item{7-306}{Roll call vote (numeric)}
#'   \item{7-307}{Roll call vote (numeric)}
#'   \item{7-308}{Roll call vote (numeric)}
#'   \item{7-309}{Roll call vote (numeric)}
#'   \item{7-31}{Roll call vote (numeric)}
#'   \item{7-310}{Roll call vote (numeric)}
#'   \item{7-311}{Roll call vote (numeric)}
#'   \item{7-312}{Roll call vote (numeric)}
#'   \item{7-313}{Roll call vote (numeric)}
#'   \item{7-314}{Roll call vote (numeric)}
#'   \item{7-315}{Roll call vote (numeric)}
#'   \item{7-316}{Roll call vote (numeric)}
#'   \item{7-317}{Roll call vote (numeric)}
#'   \item{7-318}{Roll call vote (numeric)}
#'   \item{7-319}{Roll call vote (numeric)}
#'   \item{7-32}{Roll call vote (numeric)}
#'   \item{7-320}{Roll call vote (numeric)}
#'   \item{7-321}{Roll call vote (numeric)}
#'   \item{7-322}{Roll call vote (numeric)}
#'   \item{7-323}{Roll call vote (numeric)}
#'   \item{7-324}{Roll call vote (numeric)}
#'   \item{7-325}{Roll call vote (numeric)}
#'   \item{7-326}{Roll call vote (numeric)}
#'   \item{7-327}{Roll call vote (numeric)}
#'   \item{7-328}{Roll call vote (numeric)}
#'   \item{7-329}{Roll call vote (numeric)}
#'   \item{7-33}{Roll call vote (numeric)}
#'   \item{7-330}{Roll call vote (numeric)}
#'   \item{7-331}{Roll call vote (numeric)}
#'   \item{7-332}{Roll call vote (numeric)}
#'   \item{7-333}{Roll call vote (numeric)}
#'   \item{7-334}{Roll call vote (numeric)}
#'   \item{7-335}{Roll call vote (numeric)}
#'   \item{7-336}{Roll call vote (numeric)}
#'   \item{7-337}{Roll call vote (numeric)}
#'   \item{7-338}{Roll call vote (numeric)}
#'   \item{7-339}{Roll call vote (numeric)}
#'   \item{7-34}{Roll call vote (numeric)}
#'   \item{7-340}{Roll call vote (numeric)}
#'   \item{7-341}{Roll call vote (numeric)}
#'   \item{7-342}{Roll call vote (numeric)}
#'   \item{7-343}{Roll call vote (numeric)}
#'   \item{7-344}{Roll call vote (numeric)}
#'   \item{7-345}{Roll call vote (numeric)}
#'   \item{7-346}{Roll call vote (numeric)}
#'   \item{7-347}{Roll call vote (numeric)}
#'   \item{7-348}{Roll call vote (numeric)}
#'   \item{7-349}{Roll call vote (numeric)}
#'   \item{7-35}{Roll call vote (numeric)}
#'   \item{7-350}{Roll call vote (numeric)}
#'   \item{7-351}{Roll call vote (numeric)}
#'   \item{7-352}{Roll call vote (numeric)}
#'   \item{7-353}{Roll call vote (numeric)}
#'   \item{7-354}{Roll call vote (numeric)}
#'   \item{7-355}{Roll call vote (numeric)}
#'   \item{7-356}{Roll call vote (numeric)}
#'   \item{7-357}{Roll call vote (numeric)}
#'   \item{7-358}{Roll call vote (numeric)}
#'   \item{7-359}{Roll call vote (numeric)}
#'   \item{7-36}{Roll call vote (numeric)}
#'   \item{7-360}{Roll call vote (numeric)}
#'   \item{7-361}{Roll call vote (numeric)}
#'   \item{7-362}{Roll call vote (numeric)}
#'   \item{7-363}{Roll call vote (numeric)}
#'   \item{7-364}{Roll call vote (numeric)}
#'   \item{7-365}{Roll call vote (numeric)}
#'   \item{7-366}{Roll call vote (numeric)}
#'   \item{7-367}{Roll call vote (numeric)}
#'   \item{7-368}{Roll call vote (numeric)}
#'   \item{7-369}{Roll call vote (numeric)}
#'   \item{7-37}{Roll call vote (numeric)}
#'   \item{7-370}{Roll call vote (numeric)}
#'   \item{7-371}{Roll call vote (numeric)}
#'   \item{7-372}{Roll call vote (numeric)}
#'   \item{7-373}{Roll call vote (numeric)}
#'   \item{7-374}{Roll call vote (numeric)}
#'   \item{7-375}{Roll call vote (numeric)}
#'   \item{7-376}{Roll call vote (numeric)}
#'   \item{7-377}{Roll call vote (numeric)}
#'   \item{7-378}{Roll call vote (numeric)}
#'   \item{7-379}{Roll call vote (numeric)}
#'   \item{7-38}{Roll call vote (numeric)}
#'   \item{7-380}{Roll call vote (numeric)}
#'   \item{7-381}{Roll call vote (numeric)}
#'   \item{7-382}{Roll call vote (numeric)}
#'   \item{7-383}{Roll call vote (numeric)}
#'   \item{7-384}{Roll call vote (numeric)}
#'   \item{7-385}{Roll call vote (numeric)}
#'   \item{7-386}{Roll call vote (numeric)}
#'   \item{7-387}{Roll call vote (numeric)}
#'   \item{7-388}{Roll call vote (numeric)}
#'   \item{7-389}{Roll call vote (numeric)}
#'   \item{7-39}{Roll call vote (numeric)}
#'   \item{7-390}{Roll call vote (numeric)}
#'   \item{7-391}{Roll call vote (numeric)}
#'   \item{7-392}{Roll call vote (numeric)}
#'   \item{7-393}{Roll call vote (numeric)}
#'   \item{7-394}{Roll call vote (numeric)}
#'   \item{7-395}{Roll call vote (numeric)}
#'   \item{7-396}{Roll call vote (numeric)}
#'   \item{7-397}{Roll call vote (numeric)}
#'   \item{7-398}{Roll call vote (numeric)}
#'   \item{7-399}{Roll call vote (numeric)}
#'   \item{7-4}{Roll call vote (numeric)}
#'   \item{7-40}{Roll call vote (numeric)}
#'   \item{7-400}{Roll call vote (numeric)}
#'   \item{7-401}{Roll call vote (numeric)}
#'   \item{7-402}{Roll call vote (numeric)}
#'   \item{7-403}{Roll call vote (numeric)}
#'   \item{7-404}{Roll call vote (numeric)}
#'   \item{7-405}{Roll call vote (numeric)}
#'   \item{7-406}{Roll call vote (numeric)}
#'   \item{7-407}{Roll call vote (numeric)}
#'   \item{7-408}{Roll call vote (numeric)}
#'   \item{7-409}{Roll call vote (numeric)}
#'   \item{7-41}{Roll call vote (numeric)}
#'   \item{7-410}{Roll call vote (numeric)}
#'   \item{7-411}{Roll call vote (numeric)}
#'   \item{7-412}{Roll call vote (numeric)}
#'   \item{7-413}{Roll call vote (numeric)}
#'   \item{7-414}{Roll call vote (numeric)}
#'   \item{7-415}{Roll call vote (numeric)}
#'   \item{7-416}{Roll call vote (numeric)}
#'   \item{7-417}{Roll call vote (numeric)}
#'   \item{7-418}{Roll call vote (numeric)}
#'   \item{7-419}{Roll call vote (numeric)}
#'   \item{7-42}{Roll call vote (numeric)}
#'   \item{7-420}{Roll call vote (numeric)}
#'   \item{7-421}{Roll call vote (numeric)}
#'   \item{7-422}{Roll call vote (numeric)}
#'   \item{7-423}{Roll call vote (numeric)}
#'   \item{7-424}{Roll call vote (numeric)}
#'   \item{7-425}{Roll call vote (numeric)}
#'   \item{7-426}{Roll call vote (numeric)}
#'   \item{7-427}{Roll call vote (numeric)}
#'   \item{7-428}{Roll call vote (numeric)}
#'   \item{7-429}{Roll call vote (numeric)}
#'   \item{7-43}{Roll call vote (numeric)}
#'   \item{7-430}{Roll call vote (numeric)}
#'   \item{7-431}{Roll call vote (numeric)}
#'   \item{7-432}{Roll call vote (numeric)}
#'   \item{7-433}{Roll call vote (numeric)}
#'   \item{7-434}{Roll call vote (numeric)}
#'   \item{7-435}{Roll call vote (numeric)}
#'   \item{7-436}{Roll call vote (numeric)}
#'   \item{7-437}{Roll call vote (numeric)}
#'   \item{7-438}{Roll call vote (numeric)}
#'   \item{7-439}{Roll call vote (numeric)}
#'   \item{7-44}{Roll call vote (numeric)}
#'   \item{7-440}{Roll call vote (numeric)}
#'   \item{7-441}{Roll call vote (numeric)}
#'   \item{7-442}{Roll call vote (numeric)}
#'   \item{7-443}{Roll call vote (numeric)}
#'   \item{7-444}{Roll call vote (numeric)}
#'   \item{7-445}{Roll call vote (numeric)}
#'   \item{7-446}{Roll call vote (numeric)}
#'   \item{7-447}{Roll call vote (numeric)}
#'   \item{7-448}{Roll call vote (numeric)}
#'   \item{7-449}{Roll call vote (numeric)}
#'   \item{7-45}{Roll call vote (numeric)}
#'   \item{7-450}{Roll call vote (numeric)}
#'   \item{7-451}{Roll call vote (numeric)}
#'   \item{7-452}{Roll call vote (numeric)}
#'   \item{7-453}{Roll call vote (numeric)}
#'   \item{7-454}{Roll call vote (numeric)}
#'   \item{7-455}{Roll call vote (numeric)}
#'   \item{7-456}{Roll call vote (numeric)}
#'   \item{7-457}{Roll call vote (numeric)}
#'   \item{7-458}{Roll call vote (numeric)}
#'   \item{7-459}{Roll call vote (numeric)}
#'   \item{7-46}{Roll call vote (numeric)}
#'   \item{7-460}{Roll call vote (numeric)}
#'   \item{7-461}{Roll call vote (numeric)}
#'   \item{7-462}{Roll call vote (numeric)}
#'   \item{7-463}{Roll call vote (numeric)}
#'   \item{7-464}{Roll call vote (numeric)}
#'   \item{7-465}{Roll call vote (numeric)}
#'   \item{7-466}{Roll call vote (numeric)}
#'   \item{7-467}{Roll call vote (numeric)}
#'   \item{7-468}{Roll call vote (numeric)}
#'   \item{7-469}{Roll call vote (numeric)}
#'   \item{7-47}{Roll call vote (numeric)}
#'   \item{7-470}{Roll call vote (numeric)}
#'   \item{7-471}{Roll call vote (numeric)}
#'   \item{7-472}{Roll call vote (numeric)}
#'   \item{7-473}{Roll call vote (numeric)}
#'   \item{7-474}{Roll call vote (numeric)}
#'   \item{7-475}{Roll call vote (numeric)}
#'   \item{7-476}{Roll call vote (numeric)}
#'   \item{7-477}{Roll call vote (numeric)}
#'   \item{7-478}{Roll call vote (numeric)}
#'   \item{7-479}{Roll call vote (numeric)}
#'   \item{7-48}{Roll call vote (numeric)}
#'   \item{7-480}{Roll call vote (numeric)}
#'   \item{7-481}{Roll call vote (numeric)}
#'   \item{7-482}{Roll call vote (numeric)}
#'   \item{7-483}{Roll call vote (numeric)}
#'   \item{7-484}{Roll call vote (numeric)}
#'   \item{7-485}{Roll call vote (numeric)}
#'   \item{7-486}{Roll call vote (numeric)}
#'   \item{7-487}{Roll call vote (numeric)}
#'   \item{7-488}{Roll call vote (numeric)}
#'   \item{7-489}{Roll call vote (numeric)}
#'   \item{7-49}{Roll call vote (numeric)}
#'   \item{7-490}{Roll call vote (numeric)}
#'   \item{7-491}{Roll call vote (numeric)}
#'   \item{7-492}{Roll call vote (numeric)}
#'   \item{7-493}{Roll call vote (numeric)}
#'   \item{7-494}{Roll call vote (numeric)}
#'   \item{7-495}{Roll call vote (numeric)}
#'   \item{7-496}{Roll call vote (numeric)}
#'   \item{7-497}{Roll call vote (numeric)}
#'   \item{7-498}{Roll call vote (numeric)}
#'   \item{7-499}{Roll call vote (numeric)}
#'   \item{7-5}{Roll call vote (numeric)}
#'   \item{7-50}{Roll call vote (numeric)}
#'   \item{7-500}{Roll call vote (numeric)}
#'   \item{7-501}{Roll call vote (numeric)}
#'   \item{7-502}{Roll call vote (numeric)}
#'   \item{7-503}{Roll call vote (numeric)}
#'   \item{7-504}{Roll call vote (numeric)}
#'   \item{7-505}{Roll call vote (numeric)}
#'   \item{7-506}{Roll call vote (numeric)}
#'   \item{7-507}{Roll call vote (numeric)}
#'   \item{7-508}{Roll call vote (numeric)}
#'   \item{7-509}{Roll call vote (numeric)}
#'   \item{7-51}{Roll call vote (numeric)}
#'   \item{7-510}{Roll call vote (numeric)}
#'   \item{7-511}{Roll call vote (numeric)}
#'   \item{7-512}{Roll call vote (numeric)}
#'   \item{7-513}{Roll call vote (numeric)}
#'   \item{7-514}{Roll call vote (numeric)}
#'   \item{7-515}{Roll call vote (numeric)}
#'   \item{7-516}{Roll call vote (numeric)}
#'   \item{7-517}{Roll call vote (numeric)}
#'   \item{7-518}{Roll call vote (numeric)}
#'   \item{7-519}{Roll call vote (numeric)}
#'   \item{7-52}{Roll call vote (numeric)}
#'   \item{7-520}{Roll call vote (numeric)}
#'   \item{7-521}{Roll call vote (numeric)}
#'   \item{7-522}{Roll call vote (numeric)}
#'   \item{7-523}{Roll call vote (numeric)}
#'   \item{7-524}{Roll call vote (numeric)}
#'   \item{7-525}{Roll call vote (numeric)}
#'   \item{7-526}{Roll call vote (numeric)}
#'   \item{7-527}{Roll call vote (numeric)}
#'   \item{7-528}{Roll call vote (numeric)}
#'   \item{7-529}{Roll call vote (numeric)}
#'   \item{7-53}{Roll call vote (numeric)}
#'   \item{7-530}{Roll call vote (numeric)}
#'   \item{7-531}{Roll call vote (numeric)}
#'   \item{7-532}{Roll call vote (numeric)}
#'   \item{7-533}{Roll call vote (numeric)}
#'   \item{7-534}{Roll call vote (numeric)}
#'   \item{7-535}{Roll call vote (numeric)}
#'   \item{7-536}{Roll call vote (numeric)}
#'   \item{7-537}{Roll call vote (numeric)}
#'   \item{7-538}{Roll call vote (numeric)}
#'   \item{7-539}{Roll call vote (numeric)}
#'   \item{7-54}{Roll call vote (numeric)}
#'   \item{7-540}{Roll call vote (numeric)}
#'   \item{7-541}{Roll call vote (numeric)}
#'   \item{7-542}{Roll call vote (numeric)}
#'   \item{7-543}{Roll call vote (numeric)}
#'   \item{7-544}{Roll call vote (numeric)}
#'   \item{7-545}{Roll call vote (numeric)}
#'   \item{7-546}{Roll call vote (numeric)}
#'   \item{7-547}{Roll call vote (numeric)}
#'   \item{7-548}{Roll call vote (numeric)}
#'   \item{7-549}{Roll call vote (numeric)}
#'   \item{7-55}{Roll call vote (numeric)}
#'   \item{7-550}{Roll call vote (numeric)}
#'   \item{7-551}{Roll call vote (numeric)}
#'   \item{7-552}{Roll call vote (numeric)}
#'   \item{7-553}{Roll call vote (numeric)}
#'   \item{7-554}{Roll call vote (numeric)}
#'   \item{7-555}{Roll call vote (numeric)}
#'   \item{7-556}{Roll call vote (numeric)}
#'   \item{7-557}{Roll call vote (numeric)}
#'   \item{7-558}{Roll call vote (numeric)}
#'   \item{7-559}{Roll call vote (numeric)}
#'   \item{7-56}{Roll call vote (numeric)}
#'   \item{7-560}{Roll call vote (numeric)}
#'   \item{7-561}{Roll call vote (numeric)}
#'   \item{7-562}{Roll call vote (numeric)}
#'   \item{7-563}{Roll call vote (numeric)}
#'   \item{7-564}{Roll call vote (numeric)}
#'   \item{7-565}{Roll call vote (numeric)}
#'   \item{7-566}{Roll call vote (numeric)}
#'   \item{7-567}{Roll call vote (numeric)}
#'   \item{7-568}{Roll call vote (numeric)}
#'   \item{7-569}{Roll call vote (numeric)}
#'   \item{7-57}{Roll call vote (numeric)}
#'   \item{7-570}{Roll call vote (numeric)}
#'   \item{7-571}{Roll call vote (numeric)}
#'   \item{7-572}{Roll call vote (numeric)}
#'   \item{7-573}{Roll call vote (numeric)}
#'   \item{7-574}{Roll call vote (numeric)}
#'   \item{7-575}{Roll call vote (numeric)}
#'   \item{7-576}{Roll call vote (numeric)}
#'   \item{7-577}{Roll call vote (numeric)}
#'   \item{7-578}{Roll call vote (numeric)}
#'   \item{7-579}{Roll call vote (numeric)}
#'   \item{7-58}{Roll call vote (numeric)}
#'   \item{7-580}{Roll call vote (numeric)}
#'   \item{7-581}{Roll call vote (numeric)}
#'   \item{7-582}{Roll call vote (numeric)}
#'   \item{7-583}{Roll call vote (numeric)}
#'   \item{7-584}{Roll call vote (numeric)}
#'   \item{7-585}{Roll call vote (numeric)}
#'   \item{7-586}{Roll call vote (numeric)}
#'   \item{7-587}{Roll call vote (numeric)}
#'   \item{7-588}{Roll call vote (numeric)}
#'   \item{7-589}{Roll call vote (numeric)}
#'   \item{7-59}{Roll call vote (numeric)}
#'   \item{7-590}{Roll call vote (numeric)}
#'   \item{7-591}{Roll call vote (numeric)}
#'   \item{7-592}{Roll call vote (numeric)}
#'   \item{7-593}{Roll call vote (numeric)}
#'   \item{7-594}{Roll call vote (numeric)}
#'   \item{7-595}{Roll call vote (numeric)}
#'   \item{7-596}{Roll call vote (numeric)}
#'   \item{7-597}{Roll call vote (numeric)}
#'   \item{7-598}{Roll call vote (numeric)}
#'   \item{7-599}{Roll call vote (numeric)}
#'   \item{7-6}{Roll call vote (numeric)}
#'   \item{7-60}{Roll call vote (numeric)}
#'   \item{7-600}{Roll call vote (numeric)}
#'   \item{7-601}{Roll call vote (numeric)}
#'   \item{7-602}{Roll call vote (numeric)}
#'   \item{7-603}{Roll call vote (numeric)}
#'   \item{7-604}{Roll call vote (numeric)}
#'   \item{7-605}{Roll call vote (numeric)}
#'   \item{7-606}{Roll call vote (numeric)}
#'   \item{7-607}{Roll call vote (numeric)}
#'   \item{7-608}{Roll call vote (numeric)}
#'   \item{7-609}{Roll call vote (numeric)}
#'   \item{7-61}{Roll call vote (numeric)}
#'   \item{7-610}{Roll call vote (numeric)}
#'   \item{7-611}{Roll call vote (numeric)}
#'   \item{7-612}{Roll call vote (numeric)}
#'   \item{7-613}{Roll call vote (numeric)}
#'   \item{7-614}{Roll call vote (numeric)}
#'   \item{7-615}{Roll call vote (numeric)}
#'   \item{7-616}{Roll call vote (numeric)}
#'   \item{7-617}{Roll call vote (numeric)}
#'   \item{7-618}{Roll call vote (numeric)}
#'   \item{7-619}{Roll call vote (numeric)}
#'   \item{7-62}{Roll call vote (numeric)}
#'   \item{7-620}{Roll call vote (numeric)}
#'   \item{7-621}{Roll call vote (numeric)}
#'   \item{7-622}{Roll call vote (numeric)}
#'   \item{7-623}{Roll call vote (numeric)}
#'   \item{7-624}{Roll call vote (numeric)}
#'   \item{7-625}{Roll call vote (numeric)}
#'   \item{7-626}{Roll call vote (numeric)}
#'   \item{7-627}{Roll call vote (numeric)}
#'   \item{7-628}{Roll call vote (numeric)}
#'   \item{7-629}{Roll call vote (numeric)}
#'   \item{7-63}{Roll call vote (numeric)}
#'   \item{7-630}{Roll call vote (numeric)}
#'   \item{7-631}{Roll call vote (numeric)}
#'   \item{7-632}{Roll call vote (numeric)}
#'   \item{7-633}{Roll call vote (numeric)}
#'   \item{7-634}{Roll call vote (numeric)}
#'   \item{7-635}{Roll call vote (numeric)}
#'   \item{7-636}{Roll call vote (numeric)}
#'   \item{7-637}{Roll call vote (numeric)}
#'   \item{7-638}{Roll call vote (numeric)}
#'   \item{7-639}{Roll call vote (numeric)}
#'   \item{7-64}{Roll call vote (numeric)}
#'   \item{7-640}{Roll call vote (numeric)}
#'   \item{7-641}{Roll call vote (numeric)}
#'   \item{7-642}{Roll call vote (numeric)}
#'   \item{7-643}{Roll call vote (numeric)}
#'   \item{7-644}{Roll call vote (numeric)}
#'   \item{7-645}{Roll call vote (numeric)}
#'   \item{7-646}{Roll call vote (numeric)}
#'   \item{7-647}{Roll call vote (numeric)}
#'   \item{7-648}{Roll call vote (numeric)}
#'   \item{7-649}{Roll call vote (numeric)}
#'   \item{7-65}{Roll call vote (numeric)}
#'   \item{7-650}{Roll call vote (numeric)}
#'   \item{7-651}{Roll call vote (numeric)}
#'   \item{7-652}{Roll call vote (numeric)}
#'   \item{7-653}{Roll call vote (numeric)}
#'   \item{7-654}{Roll call vote (numeric)}
#'   \item{7-655}{Roll call vote (numeric)}
#'   \item{7-656}{Roll call vote (numeric)}
#'   \item{7-657}{Roll call vote (numeric)}
#'   \item{7-658}{Roll call vote (numeric)}
#'   \item{7-659}{Roll call vote (numeric)}
#'   \item{7-66}{Roll call vote (numeric)}
#'   \item{7-660}{Roll call vote (numeric)}
#'   \item{7-661}{Roll call vote (numeric)}
#'   \item{7-662}{Roll call vote (numeric)}
#'   \item{7-663}{Roll call vote (numeric)}
#'   \item{7-664}{Roll call vote (numeric)}
#'   \item{7-665}{Roll call vote (numeric)}
#'   \item{7-666}{Roll call vote (numeric)}
#'   \item{7-667}{Roll call vote (numeric)}
#'   \item{7-668}{Roll call vote (numeric)}
#'   \item{7-669}{Roll call vote (numeric)}
#'   \item{7-67}{Roll call vote (numeric)}
#'   \item{7-670}{Roll call vote (numeric)}
#'   \item{7-671}{Roll call vote (numeric)}
#'   \item{7-672}{Roll call vote (numeric)}
#'   \item{7-673}{Roll call vote (numeric)}
#'   \item{7-674}{Roll call vote (numeric)}
#'   \item{7-675}{Roll call vote (numeric)}
#'   \item{7-676}{Roll call vote (numeric)}
#'   \item{7-677}{Roll call vote (numeric)}
#'   \item{7-678}{Roll call vote (numeric)}
#'   \item{7-679}{Roll call vote (numeric)}
#'   \item{7-68}{Roll call vote (numeric)}
#'   \item{7-680}{Roll call vote (numeric)}
#'   \item{7-681}{Roll call vote (numeric)}
#'   \item{7-682}{Roll call vote (numeric)}
#'   \item{7-683}{Roll call vote (numeric)}
#'   \item{7-684}{Roll call vote (numeric)}
#'   \item{7-685}{Roll call vote (numeric)}
#'   \item{7-686}{Roll call vote (numeric)}
#'   \item{7-687}{Roll call vote (numeric)}
#'   \item{7-688}{Roll call vote (numeric)}
#'   \item{7-689}{Roll call vote (numeric)}
#'   \item{7-69}{Roll call vote (numeric)}
#'   \item{7-690}{Roll call vote (numeric)}
#'   \item{7-691}{Roll call vote (numeric)}
#'   \item{7-692}{Roll call vote (numeric)}
#'   \item{7-693}{Roll call vote (numeric)}
#'   \item{7-694}{Roll call vote (numeric)}
#'   \item{7-695}{Roll call vote (numeric)}
#'   \item{7-696}{Roll call vote (numeric)}
#'   \item{7-697}{Roll call vote (numeric)}
#'   \item{7-698}{Roll call vote (numeric)}
#'   \item{7-699}{Roll call vote (numeric)}
#'   \item{7-7}{Roll call vote (numeric)}
#'   \item{7-70}{Roll call vote (numeric)}
#'   \item{7-700}{Roll call vote (numeric)}
#'   \item{7-701}{Roll call vote (numeric)}
#'   \item{7-702}{Roll call vote (numeric)}
#'   \item{7-703}{Roll call vote (numeric)}
#'   \item{7-704}{Roll call vote (numeric)}
#'   \item{7-705}{Roll call vote (numeric)}
#'   \item{7-706}{Roll call vote (numeric)}
#'   \item{7-707}{Roll call vote (numeric)}
#'   \item{7-708}{Roll call vote (numeric)}
#'   \item{7-709}{Roll call vote (numeric)}
#'   \item{7-71}{Roll call vote (numeric)}
#'   \item{7-710}{Roll call vote (numeric)}
#'   \item{7-711}{Roll call vote (numeric)}
#'   \item{7-712}{Roll call vote (numeric)}
#'   \item{7-713}{Roll call vote (numeric)}
#'   \item{7-714}{Roll call vote (numeric)}
#'   \item{7-715}{Roll call vote (numeric)}
#'   \item{7-716}{Roll call vote (numeric)}
#'   \item{7-717}{Roll call vote (numeric)}
#'   \item{7-718}{Roll call vote (numeric)}
#'   \item{7-719}{Roll call vote (numeric)}
#'   \item{7-72}{Roll call vote (numeric)}
#'   \item{7-720}{Roll call vote (numeric)}
#'   \item{7-721}{Roll call vote (numeric)}
#'   \item{7-722}{Roll call vote (numeric)}
#'   \item{7-723}{Roll call vote (numeric)}
#'   \item{7-724}{Roll call vote (numeric)}
#'   \item{7-725}{Roll call vote (numeric)}
#'   \item{7-726}{Roll call vote (numeric)}
#'   \item{7-727}{Roll call vote (numeric)}
#'   \item{7-728}{Roll call vote (numeric)}
#'   \item{7-729}{Roll call vote (numeric)}
#'   \item{7-73}{Roll call vote (numeric)}
#'   \item{7-730}{Roll call vote (numeric)}
#'   \item{7-731}{Roll call vote (numeric)}
#'   \item{7-732}{Roll call vote (numeric)}
#'   \item{7-733}{Roll call vote (numeric)}
#'   \item{7-734}{Roll call vote (numeric)}
#'   \item{7-735}{Roll call vote (numeric)}
#'   \item{7-736}{Roll call vote (numeric)}
#'   \item{7-737}{Roll call vote (numeric)}
#'   \item{7-738}{Roll call vote (numeric)}
#'   \item{7-739}{Roll call vote (numeric)}
#'   \item{7-74}{Roll call vote (numeric)}
#'   \item{7-740}{Roll call vote (numeric)}
#'   \item{7-741}{Roll call vote (numeric)}
#'   \item{7-742}{Roll call vote (numeric)}
#'   \item{7-743}{Roll call vote (numeric)}
#'   \item{7-744}{Roll call vote (numeric)}
#'   \item{7-745}{Roll call vote (numeric)}
#'   \item{7-746}{Roll call vote (numeric)}
#'   \item{7-747}{Roll call vote (numeric)}
#'   \item{7-748}{Roll call vote (numeric)}
#'   \item{7-749}{Roll call vote (numeric)}
#'   \item{7-75}{Roll call vote (numeric)}
#'   \item{7-750}{Roll call vote (numeric)}
#'   \item{7-751}{Roll call vote (numeric)}
#'   \item{7-752}{Roll call vote (numeric)}
#'   \item{7-753}{Roll call vote (numeric)}
#'   \item{7-754}{Roll call vote (numeric)}
#'   \item{7-755}{Roll call vote (numeric)}
#'   \item{7-756}{Roll call vote (numeric)}
#'   \item{7-757}{Roll call vote (numeric)}
#'   \item{7-758}{Roll call vote (numeric)}
#'   \item{7-759}{Roll call vote (numeric)}
#'   \item{7-76}{Roll call vote (numeric)}
#'   \item{7-760}{Roll call vote (numeric)}
#'   \item{7-761}{Roll call vote (numeric)}
#'   \item{7-762}{Roll call vote (numeric)}
#'   \item{7-763}{Roll call vote (numeric)}
#'   \item{7-764}{Roll call vote (numeric)}
#'   \item{7-765}{Roll call vote (numeric)}
#'   \item{7-766}{Roll call vote (numeric)}
#'   \item{7-767}{Roll call vote (numeric)}
#'   \item{7-768}{Roll call vote (numeric)}
#'   \item{7-769}{Roll call vote (numeric)}
#'   \item{7-77}{Roll call vote (numeric)}
#'   \item{7-770}{Roll call vote (numeric)}
#'   \item{7-771}{Roll call vote (numeric)}
#'   \item{7-772}{Roll call vote (numeric)}
#'   \item{7-773}{Roll call vote (numeric)}
#'   \item{7-774}{Roll call vote (numeric)}
#'   \item{7-775}{Roll call vote (numeric)}
#'   \item{7-776}{Roll call vote (numeric)}
#'   \item{7-777}{Roll call vote (numeric)}
#'   \item{7-778}{Roll call vote (numeric)}
#'   \item{7-779}{Roll call vote (numeric)}
#'   \item{7-78}{Roll call vote (numeric)}
#'   \item{7-780}{Roll call vote (numeric)}
#'   \item{7-781}{Roll call vote (numeric)}
#'   \item{7-782}{Roll call vote (numeric)}
#'   \item{7-783}{Roll call vote (numeric)}
#'   \item{7-784}{Roll call vote (numeric)}
#'   \item{7-785}{Roll call vote (numeric)}
#'   \item{7-786}{Roll call vote (numeric)}
#'   \item{7-787}{Roll call vote (numeric)}
#'   \item{7-788}{Roll call vote (numeric)}
#'   \item{7-789}{Roll call vote (numeric)}
#'   \item{7-79}{Roll call vote (numeric)}
#'   \item{7-790}{Roll call vote (numeric)}
#'   \item{7-791}{Roll call vote (numeric)}
#'   \item{7-792}{Roll call vote (numeric)}
#'   \item{7-793}{Roll call vote (numeric)}
#'   \item{7-794}{Roll call vote (numeric)}
#'   \item{7-795}{Roll call vote (numeric)}
#'   \item{7-796}{Roll call vote (numeric)}
#'   \item{7-797}{Roll call vote (numeric)}
#'   \item{7-798}{Roll call vote (numeric)}
#'   \item{7-799}{Roll call vote (numeric)}
#'   \item{7-8}{Roll call vote (numeric)}
#'   \item{7-80}{Roll call vote (numeric)}
#'   \item{7-800}{Roll call vote (numeric)}
#'   \item{7-801}{Roll call vote (numeric)}
#'   \item{7-802}{Roll call vote (numeric)}
#'   \item{7-803}{Roll call vote (numeric)}
#'   \item{7-804}{Roll call vote (numeric)}
#'   \item{7-805}{Roll call vote (numeric)}
#'   \item{7-806}{Roll call vote (numeric)}
#'   \item{7-807}{Roll call vote (numeric)}
#'   \item{7-808}{Roll call vote (numeric)}
#'   \item{7-809}{Roll call vote (numeric)}
#'   \item{7-81}{Roll call vote (numeric)}
#'   \item{7-810}{Roll call vote (numeric)}
#'   \item{7-811}{Roll call vote (numeric)}
#'   \item{7-812}{Roll call vote (numeric)}
#'   \item{7-813}{Roll call vote (numeric)}
#'   \item{7-814}{Roll call vote (numeric)}
#'   \item{7-815}{Roll call vote (numeric)}
#'   \item{7-816}{Roll call vote (numeric)}
#'   \item{7-817}{Roll call vote (numeric)}
#'   \item{7-818}{Roll call vote (numeric)}
#'   \item{7-819}{Roll call vote (numeric)}
#'   \item{7-82}{Roll call vote (numeric)}
#'   \item{7-820}{Roll call vote (numeric)}
#'   \item{7-821}{Roll call vote (numeric)}
#'   \item{7-822}{Roll call vote (numeric)}
#'   \item{7-823}{Roll call vote (numeric)}
#'   \item{7-824}{Roll call vote (numeric)}
#'   \item{7-825}{Roll call vote (numeric)}
#'   \item{7-826}{Roll call vote (numeric)}
#'   \item{7-827}{Roll call vote (numeric)}
#'   \item{7-828}{Roll call vote (numeric)}
#'   \item{7-829}{Roll call vote (numeric)}
#'   \item{7-83}{Roll call vote (numeric)}
#'   \item{7-830}{Roll call vote (numeric)}
#'   \item{7-831}{Roll call vote (numeric)}
#'   \item{7-832}{Roll call vote (numeric)}
#'   \item{7-833}{Roll call vote (numeric)}
#'   \item{7-834}{Roll call vote (numeric)}
#'   \item{7-835}{Roll call vote (numeric)}
#'   \item{7-836}{Roll call vote (numeric)}
#'   \item{7-837}{Roll call vote (numeric)}
#'   \item{7-838}{Roll call vote (numeric)}
#'   \item{7-839}{Roll call vote (numeric)}
#'   \item{7-84}{Roll call vote (numeric)}
#'   \item{7-840}{Roll call vote (numeric)}
#'   \item{7-841}{Roll call vote (numeric)}
#'   \item{7-842}{Roll call vote (numeric)}
#'   \item{7-843}{Roll call vote (numeric)}
#'   \item{7-844}{Roll call vote (numeric)}
#'   \item{7-845}{Roll call vote (numeric)}
#'   \item{7-846}{Roll call vote (numeric)}
#'   \item{7-847}{Roll call vote (numeric)}
#'   \item{7-848}{Roll call vote (numeric)}
#'   \item{7-849}{Roll call vote (numeric)}
#'   \item{7-85}{Roll call vote (numeric)}
#'   \item{7-850}{Roll call vote (numeric)}
#'   \item{7-851}{Roll call vote (numeric)}
#'   \item{7-852}{Roll call vote (numeric)}
#'   \item{7-853}{Roll call vote (numeric)}
#'   \item{7-854}{Roll call vote (numeric)}
#'   \item{7-855}{Roll call vote (numeric)}
#'   \item{7-856}{Roll call vote (numeric)}
#'   \item{7-857}{Roll call vote (numeric)}
#'   \item{7-858}{Roll call vote (numeric)}
#'   \item{7-859}{Roll call vote (numeric)}
#'   \item{7-86}{Roll call vote (numeric)}
#'   \item{7-860}{Roll call vote (numeric)}
#'   \item{7-861}{Roll call vote (numeric)}
#'   \item{7-862}{Roll call vote (numeric)}
#'   \item{7-863}{Roll call vote (numeric)}
#'   \item{7-864}{Roll call vote (numeric)}
#'   \item{7-865}{Roll call vote (numeric)}
#'   \item{7-866}{Roll call vote (numeric)}
#'   \item{7-867}{Roll call vote (numeric)}
#'   \item{7-868}{Roll call vote (numeric)}
#'   \item{7-869}{Roll call vote (numeric)}
#'   \item{7-87}{Roll call vote (numeric)}
#'   \item{7-870}{Roll call vote (numeric)}
#'   \item{7-871}{Roll call vote (numeric)}
#'   \item{7-872}{Roll call vote (numeric)}
#'   \item{7-873}{Roll call vote (numeric)}
#'   \item{7-874}{Roll call vote (numeric)}
#'   \item{7-875}{Roll call vote (numeric)}
#'   \item{7-876}{Roll call vote (numeric)}
#'   \item{7-877}{Roll call vote (numeric)}
#'   \item{7-878}{Roll call vote (numeric)}
#'   \item{7-879}{Roll call vote (numeric)}
#'   \item{7-88}{Roll call vote (numeric)}
#'   \item{7-880}{Roll call vote (numeric)}
#'   \item{7-881}{Roll call vote (numeric)}
#'   \item{7-882}{Roll call vote (numeric)}
#'   \item{7-883}{Roll call vote (numeric)}
#'   \item{7-884}{Roll call vote (numeric)}
#'   \item{7-885}{Roll call vote (numeric)}
#'   \item{7-886}{Roll call vote (numeric)}
#'   \item{7-887}{Roll call vote (numeric)}
#'   \item{7-888}{Roll call vote (numeric)}
#'   \item{7-889}{Roll call vote (numeric)}
#'   \item{7-89}{Roll call vote (numeric)}
#'   \item{7-890}{Roll call vote (numeric)}
#'   \item{7-891}{Roll call vote (numeric)}
#'   \item{7-892}{Roll call vote (numeric)}
#'   \item{7-893}{Roll call vote (numeric)}
#'   \item{7-894}{Roll call vote (numeric)}
#'   \item{7-895}{Roll call vote (numeric)}
#'   \item{7-896}{Roll call vote (numeric)}
#'   \item{7-897}{Roll call vote (numeric)}
#'   \item{7-898}{Roll call vote (numeric)}
#'   \item{7-899}{Roll call vote (numeric)}
#'   \item{7-9}{Roll call vote (numeric)}
#'   \item{7-90}{Roll call vote (numeric)}
#'   \item{7-900}{Roll call vote (numeric)}
#'   \item{7-901}{Roll call vote (numeric)}
#'   \item{7-902}{Roll call vote (numeric)}
#'   \item{7-903}{Roll call vote (numeric)}
#'   \item{7-904}{Roll call vote (numeric)}
#'   \item{7-905}{Roll call vote (numeric)}
#'   \item{7-906}{Roll call vote (numeric)}
#'   \item{7-907}{Roll call vote (numeric)}
#'   \item{7-908}{Roll call vote (numeric)}
#'   \item{7-909}{Roll call vote (numeric)}
#'   \item{7-91}{Roll call vote (numeric)}
#'   \item{7-910}{Roll call vote (numeric)}
#'   \item{7-911}{Roll call vote (numeric)}
#'   \item{7-912}{Roll call vote (numeric)}
#'   \item{7-913}{Roll call vote (numeric)}
#'   \item{7-914}{Roll call vote (numeric)}
#'   \item{7-915}{Roll call vote (numeric)}
#'   \item{7-916}{Roll call vote (numeric)}
#'   \item{7-917}{Roll call vote (numeric)}
#'   \item{7-918}{Roll call vote (numeric)}
#'   \item{7-919}{Roll call vote (numeric)}
#'   \item{7-92}{Roll call vote (numeric)}
#'   \item{7-920}{Roll call vote (numeric)}
#'   \item{7-921}{Roll call vote (numeric)}
#'   \item{7-922}{Roll call vote (numeric)}
#'   \item{7-923}{Roll call vote (numeric)}
#'   \item{7-924}{Roll call vote (numeric)}
#'   \item{7-925}{Roll call vote (numeric)}
#'   \item{7-926}{Roll call vote (numeric)}
#'   \item{7-927}{Roll call vote (numeric)}
#'   \item{7-928}{Roll call vote (numeric)}
#'   \item{7-929}{Roll call vote (numeric)}
#'   \item{7-93}{Roll call vote (numeric)}
#'   \item{7-930}{Roll call vote (numeric)}
#'   \item{7-931}{Roll call vote (numeric)}
#'   \item{7-932}{Roll call vote (numeric)}
#'   \item{7-933}{Roll call vote (numeric)}
#'   \item{7-934}{Roll call vote (numeric)}
#'   \item{7-935}{Roll call vote (numeric)}
#'   \item{7-936}{Roll call vote (numeric)}
#'   \item{7-937}{Roll call vote (numeric)}
#'   \item{7-938}{Roll call vote (numeric)}
#'   \item{7-939}{Roll call vote (numeric)}
#'   \item{7-94}{Roll call vote (numeric)}
#'   \item{7-940}{Roll call vote (numeric)}
#'   \item{7-941}{Roll call vote (numeric)}
#'   \item{7-942}{Roll call vote (numeric)}
#'   \item{7-943}{Roll call vote (numeric)}
#'   \item{7-944}{Roll call vote (numeric)}
#'   \item{7-945}{Roll call vote (numeric)}
#'   \item{7-946}{Roll call vote (numeric)}
#'   \item{7-947}{Roll call vote (numeric)}
#'   \item{7-948}{Roll call vote (numeric)}
#'   \item{7-949}{Roll call vote (numeric)}
#'   \item{7-95}{Roll call vote (numeric)}
#'   \item{7-950}{Roll call vote (numeric)}
#'   \item{7-951}{Roll call vote (numeric)}
#'   \item{7-952}{Roll call vote (numeric)}
#'   \item{7-953}{Roll call vote (numeric)}
#'   \item{7-954}{Roll call vote (numeric)}
#'   \item{7-955}{Roll call vote (numeric)}
#'   \item{7-956}{Roll call vote (numeric)}
#'   \item{7-957}{Roll call vote (numeric)}
#'   \item{7-958}{Roll call vote (numeric)}
#'   \item{7-959}{Roll call vote (numeric)}
#'   \item{7-96}{Roll call vote (numeric)}
#'   \item{7-960}{Roll call vote (numeric)}
#'   \item{7-961}{Roll call vote (numeric)}
#'   \item{7-962}{Roll call vote (numeric)}
#'   \item{7-963}{Roll call vote (numeric)}
#'   \item{7-964}{Roll call vote (numeric)}
#'   \item{7-965}{Roll call vote (numeric)}
#'   \item{7-966}{Roll call vote (numeric)}
#'   \item{7-967}{Roll call vote (numeric)}
#'   \item{7-968}{Roll call vote (numeric)}
#'   \item{7-969}{Roll call vote (numeric)}
#'   \item{7-97}{Roll call vote (numeric)}
#'   \item{7-970}{Roll call vote (numeric)}
#'   \item{7-971}{Roll call vote (numeric)}
#'   \item{7-972}{Roll call vote (numeric)}
#'   \item{7-973}{Roll call vote (numeric)}
#'   \item{7-974}{Roll call vote (numeric)}
#'   \item{7-975}{Roll call vote (numeric)}
#'   \item{7-976}{Roll call vote (numeric)}
#'   \item{7-977}{Roll call vote (numeric)}
#'   \item{7-978}{Roll call vote (numeric)}
#'   \item{7-979}{Roll call vote (numeric)}
#'   \item{7-98}{Roll call vote (numeric)}
#'   \item{7-980}{Roll call vote (numeric)}
#'   \item{7-981}{Roll call vote (numeric)}
#'   \item{7-982}{Roll call vote (numeric)}
#'   \item{7-983}{Roll call vote (numeric)}
#'   \item{7-984}{Roll call vote (numeric)}
#'   \item{7-985}{Roll call vote (numeric)}
#'   \item{7-986}{Roll call vote (numeric)}
#'   \item{7-987}{Roll call vote (numeric)}
#'   \item{7-988}{Roll call vote (numeric)}
#'   \item{7-989}{Roll call vote (numeric)}
#'   \item{7-99}{Roll call vote (numeric)}
#'   \item{7-990}{Roll call vote (numeric)}
#'   \item{7-991}{Roll call vote (numeric)}
#'   \item{7-992}{Roll call vote (numeric)}
#'   \item{7-993}{Roll call vote (numeric)}
#'   \item{7-994}{Roll call vote (numeric)}
#'   \item{7-995}{Roll call vote (numeric)}
#'   \item{7-996}{Roll call vote (numeric)}
#'   \item{7-997}{Roll call vote (numeric)}
#'   \item{7-998}{Roll call vote (numeric)}
#'   \item{7-999}{Roll call vote (numeric)}
#' }
#'
#' @details
#' The data captures the legislative behavior during the 7th session of the 
#' Legislative Yuan of Taiwan, providing valuable insights into the political 
#' dynamics and decision-making processes. Each row represents one legislator, 
#' and each vote column (7-1 through 7-1226) represents a specific bill or 
#' motion voted upon.
#'
#' @source
#' Yen-Chieh Liao (2024). Electoral Reform and Fragmented Polarization: 
#' New Evidence from Taiwan Legislative Roll Call. Legislative Studies Quarterly. 
#' \doi{10.1111/lsq.12459}
#'
#' @usage data(legis_7th_Taiwan)
#'
#' @examples
#' \dontrun{
#' data(legis_7th_Taiwan)
#' 
#' # View structure
#' str(legis_7th_Taiwan)
#' 
#' # First few legislators and votes
#' head(legis_7th_Taiwan[, c(1:5, 1227:1228)])
#' 
#' # Summary by party
#' table(legis_7th_Taiwan$party)
#' }
#'
#' @keywords datasets
#' @name legis_7th_Taiwan
#' @docType data
NULL


#' @encoding UTF-8
#' @title 90th US Senate Agreement Score Matrix (1967-1968)
#' @description This dataset contains the agreement score matrix of the 90th US Senate, covering the years 1967-1968.
#' The dataset includes 102 legislators: 100 Senators, President Lyndon Johnson (who "voted" on select bills by announcing a position),
#' and Senator Charles Goodell (R-NY), who replaced Senator Robert F. Kennedy after his assassination in June 1968.
#'
#' @format A data frame with the following variables:
#' \describe{
#'   \item{congress}{Congressional session number. Represents the 90th Congress.}
#'   \item{id}{Legislator's unique identifier within the dataset.}
#'   \item{statecode}{Two-letter state abbreviation representing the legislator's state.}
#'   \item{statename}{Full name of the state the legislator represents.}
#'   \item{party}{Political party affiliation of the legislator, typically 1 for Democrat and 2 for Republican.}
#'   \item{election}{Year of the legislator's election to the Senate.}
#'   \item{name}{Name of the legislator. This includes both first and last names.}
#'   \item{johnson}{Agreement score with President Lyndon Johnson, reflecting how often the legislator's votes aligned with Johnson's positions.}
#'   \item{sparkman}{Agreement score with Senator John Sparkman (D-AL). Similar to other agreement scores, it measures the voting alignment with this senator.}
#'   \item{hill}{Agreement score with Senator Lister Hill (D-AL).}
#'   \item{gruening}{Agreement score with Senator Ernest Gruening (D-AK).}
#'   \item{bartlett}{Agreement score with Senator E.L. Bartlett (D-AK).}
#'   \item{hayden}{Agreement score with Senator Carl Hayden (D-AZ).}
#'   \item{fannin}{Agreement score with Senator Paul Fannin (R-AZ).}
#'   \item{fulbright}{Agreement score with Senator J. William Fulbright (D-AR).}
#'   \item{mcclellan}{Agreement score with Senator John L. McClellan (D-AR).}
#'   \item{kuchel}{Agreement score with Senator Thomas Kuchel (R-CA).}
#'   \item{murphy}{Agreement score with Senator George Murphy (R-CA).}
#'   \item{dominick}{Agreement score with Senator Peter Dominick (R-CO).}
#'   \item{allott}{Agreement score with Senator Gordon Allott (R-CO).}
#'   \item{dodd}{Agreement score with Senator Thomas Dodd (D-CT).}
#'   \item{ribicoff}{Agreement score with Senator Abraham Ribicoff (D-CT).}
#'   \item{boggs}{Agreement score with Senator J. Caleb Boggs (R-DE).}
#'   \item{williamsj}{Agreement score with Senator John Williams (R-DE).}
#'   \item{smathers}{Agreement score with Senator George Smathers (D-FL).}
#'   \item{holland}{Agreement score with Senator Spessard Holland (D-FL).}
#'   \item{russell}{Agreement score with Senator Richard Russell (D-GA).}
#'   \item{talmadge}{Agreement score with Senator Herman Talmadge (D-GA).}
#'   \item{fong}{Agreement score with Senator Hiram Fong (R-HI).}
#'   \item{inouye}{Agreement score with Senator Daniel Inouye (D-HI).}
#'   \item{church}{Agreement score with Senator Frank Church (D-ID).}
#'   \item{jordanl}{Agreement score with Senator Len Jordan (R-ID).}
#'   \item{dirksen}{Agreement score with Senator Everett Dirksen (R-IL).}
#'   \item{percy}{Agreement score with Senator Charles Percy (R-IL).}
#'   \item{hartke}{Agreement score with Senator Vance Hartke (D-IN).}
#'   \item{bayh}{Agreement score with Senator Birch Bayh (D-IN).}
#'   \item{miller}{Agreement score with Senator Jack Miller (R-IA).}
#'   \item{hickenloope}{Agreement score with Senator Bourke Hickenlooper (R-IA).}
#'   \item{carlson}{Agreement score with Senator Frank Carlson (R-KS).}
#'   \item{pearson}{Agreement score with Senator James Pearson (R-KS).}
#'   \item{cooper}{Agreement score with Senator John Sherman Cooper (R-KY).}
#'   \item{morton}{Agreement score with Senator Thruston Morton (R-KY).}
#'   \item{ellender}{Agreement score with Senator Allen Ellender (D-LA).}
#'   \item{longr}{Agreement score with Senator Russell Long (D-LA).}
#'   \item{muskie}{Agreement score with Senator Edmund Muskie (D-ME).}
#'   \item{smith}{Agreement score with Senator Margaret Chase Smith (R-ME).}
#'   \item{brewster}{Agreement score with Senator Daniel Brewster (D-MD).}
#'   \item{tydings}{Agreement score with Senator Joseph Tydings (D-MD).}
#'   \item{brooke}{Agreement score with Senator Edward Brooke (R-MA).}
#'   \item{kennedye}{Agreement score with Senator Edward Kennedy (D-MA).}
#'   \item{griffin}{Agreement score with Senator Robert Griffin (R-MI).}
#'   \item{hart}{Agreement score with Senator Philip Hart (D-MI).}
#'   \item{mondale}{Agreement score with Senator Walter Mondale (D-MN).}
#'   \item{mccarthy}{Agreement score with Senator Eugene McCarthy (D-MN).}
#'   \item{stennis}{Agreement score with Senator John Stennis (D-MS).}
#'   \item{eastland}{Agreement score with Senator James Eastland (D-MS).}
#'   \item{symington}{Agreement score with Senator Stuart Symington (D-MO).}
#'   \item{longe}{Agreement score with Senator Edward Long (D-MO).}
#'   \item{metcalf}{Agreement score with Senator Lee Metcalf (D-MT).}
#'   \item{mansfield}{Agreement score with Senator Mike Mansfield (D-MT).}
#'   \item{curtis}{Agreement score with Senator Carl Curtis (R-NE).}
#'   \item{hruska}{Agreement score with Senator Roman Hruska (R-NE).}
#'   \item{bible}{Agreement score with Senator Alan Bible (D-NV).}
#'   \item{cannon}{Agreement score with Senator Howard Cannon (D-NV).}
#'   \item{mcintyre}{Agreement score with Senator Thomas McIntyre (D-NH).}
#'   \item{cotton}{Agreement score with Senator Norris Cotton (R-NH).}
#'   \item{case}{Agreement score with Senator Clifford Case (R-NJ).}
#'   \item{williamsh}{Agreement score with Senator Harrison Williams (D-NJ).}
#'   \item{anderson}{Agreement score with Senator Clinton Anderson (D-NM).}
#'   \item{montoya}{Agreement score with Senator Joseph Montoya (D-NM).}
#'   \item{javits}{Agreement score with Senator Jacob Javits (R-NY).}
#'   \item{goodell}{Agreement score with Senator Charles Goodell (R-NY), who replaced Robert F. Kennedy.}
#'   \item{kennedyr}{Agreement score with Senator Robert F. Kennedy (D-NY).}
#'   \item{jordanb}{Agreement score with Senator B. Everett Jordan (D-NC).}
#'   \item{ervin}{Agreement score with Senator Sam Ervin (D-NC).}
#'   \item{youngm}{Agreement score with Senator Milton Young (R-ND).}
#'   \item{burdick}{Agreement score with Senator Quentin Burdick (D-ND).}
#'   \item{youngs}{Agreement score with Senator Stephen Young (D-OH).}
#'   \item{lausche}{Agreement score with Senator Frank Lausche (D-OH).}
#'   \item{monroney}{Agreement score with Senator A.S. Mike Monroney (D-OK).}
#'   \item{harris}{Agreement score with Senator Fred Harris (D-OK).}
#'   \item{morse}{Agreement score with Senator Wayne Morse (D-OR).}
#'   \item{hatfield}{Agreement score with Senator Mark Hatfield (R-OR).}
#'   \item{scott}{Agreement score with Senator Hugh Scott (R-PA).}
#'   \item{clark}{Agreement score with Senator Joseph Clark (D-PA).}
#'   \item{pastore}{Agreement score with Senator John Pastore (D-RI).}
#'   \item{pell}{Agreement score with Senator Claiborne Pell (D-RI).}
#'   \item{hollings}{Agreement score with Senator Ernest Hollings (D-SC).}
#'   \item{thurmond}{Agreement score with Senator Strom Thurmond (R-SC).}
#'   \item{mcgovern}{Agreement score with Senator George McGovern (D-SD).}
#'   \item{mundt}{Agreement score with Senator Karl Mundt (R-SD).}
#'   \item{baker}{Agreement score with Senator Howard Baker (R-TN).}
#'   \item{gore}{Agreement score with Senator Albert Gore (D-TN).}
#'   \item{yarborough}{Agreement score with Senator Ralph Yarborough (D-TX).}
#'   \item{tower}{Agreement score with Senator John Tower (R-TX).}
#'   \item{bennett}{Agreement score with Senator Wallace Bennett (R-UT).}
#'   \item{moss}{Agreement score with Senator Frank Moss (D-UT).}
#'   \item{prouty}{Agreement score with Senator Winston Prouty (R-VT).}
#'   \item{aiken}{Agreement score with Senator George Aiken (R-VT).}
#'   \item{spong}{Agreement score with Senator William Spong (D-VA).}
#'   \item{byrdh}{Agreement score with Senator Harry F. Byrd Jr. (I-VA).}
#'   \item{magnuson}{Agreement score with Senator Warren Magnuson (D-WA).}
#'   \item{jackson}{Agreement score with Senator Henry Jackson (D-WA).}
#'   \item{randolph}{Agreement score with Senator Jennings Randolph (D-WV).}
#'   \item{byrdr}{Agreement score with Senator Robert Byrd (D-WV).}
#'   \item{proxmire}{Agreement score with Senator William Proxmire (D-WI).}
#'   \item{nelson}{Agreement score with Senator Gaylord Nelson (D-WI).}
#'   \item{hansen}{Agreement score with Senator Clifford Hansen (R-WY).}
#'   \item{mcgee}{Agreement score with Senator Gale McGee (D-WY).}
#' }
#' 
#'  @details
#' The matrix is used to analyze the dimensions of voting behavior in the Senate during this period, with
#' two primary dimensions identified: liberal-conservative and region/civil-rights.
#' The data were used in the analysis by Poole and Rosenthal (1997).
#'
#' @source
#' Poole, Keith T., and Howard Rosenthal. (1997). *Congress: A Political-Economic History of Roll Call Voting*.
#' Oxford University Press.
#'
#' @usage data(senate.90)
#'
#' @examples
#' \dontrun{
#' data(senate.90)
#' }
#'
#' @keywords datasets
#' @name senate.90
#' @docType data
NULL


#' @encoding UTF-8
#' @title Bootstrapped Blackbox Analysis Output
#' @description The `outbb` object contains the results of a bootstrapped Blackbox analysis performed on the `issues.sweden` dataset. The analysis was conducted using the `boot.blackbox` function with specific parameters for handling missing data, dimensionality, scaling, and stimulus positioning. 
#' 
#' The object also includes the `prerun` output, which provides preliminary diagnostics and summary statistics prior to the bootstrapping process.
#'
#' @format A list containing the results of the bootstrapped Blackbox analysis.
#' @details 
#' The `outbb` object captures the outputs from the bootstrapped analysis, including estimated dimensions, scaling factors, and the placement of stimuli. The `prerun` output embedded within `outbb` contains initial diagnostics that help in assessing the quality and stability of the analysis before bootstrapping is applied.
#'
#' @source Generated by the `boot.blackbox` function.
#'
#' @usage outbb
#' @keywords datasets
#' @name outbb
#' @docType data
NULL


#' @encoding UTF-8
#' @title Blackbox Transpose Analysis Result for Mexico CSES 2000
#' @description  The `result_2000` object contains the results of a Blackbox transpose analysis performed on the `mexicoCSES2000` dataset. This analysis was conducted to explore the dimensional structure of the dataset with specific handling for missing data, scaling, and dimensionality reduction.
#'
#' @format A list containing the results of the Blackbox transpose analysis, including estimated dimensions, scaling factors, and other relevant metrics.
#'
#' @details
#' The `result_2000` object captures the outputs from the Blackbox transpose analysis, where the `mexicoCSES2000` dataset was analyzed across 3 dimensions with a minimum scaling factor of 5. Missing data were coded as `99`, and the analysis was performed with verbose output to provide detailed information on the process.
#'
#' @source Generated by the `blackbox_transpose` function applied to the `mexicoCSES2000` dataset.
#'
#' @usage result_2000
#' @keywords datasets
#' @name result_2000
#' @docType data
NULL


#' @encoding UTF-8
#' @title Blackbox Transpose Analysis Result for Mexico CSES 2006
#' @description The `result_2006` object contains the results of a Blackbox transpose analysis performed on the `mexicoCSES2006` dataset. This analysis was carried out to examine the dimensional structure of the dataset with specific handling for missing data, scaling, and dimensionality reduction.
#'
#' @format A list containing the results of the Blackbox transpose analysis, including estimated dimensions, scaling factors, and other relevant metrics.
#'
#' @details
#' The `result_2006` object captures the outputs from the Blackbox transpose analysis, where the `mexicoCSES2006` dataset was analyzed across 3 dimensions with a minimum scaling factor of 5. Missing data were coded as `99`, and the analysis was performed with verbose output to provide detailed information on the process.
#'
#' @source Generated by the `blackbox_transpose` function applied to the `mexicoCSES2006` dataset.
#'
#' @usage result_2006
#' @keywords datasets
#' @name result_2006
#' @docType data
NULL


#' @encoding UTF-8
#' @title Bootstrapped Blackbox Transpose Analysis Output
#' @description The `outbbt` object contains the results of a bootstrapped Blackbox transpose analysis performed on the `rankings` dataset. This analysis was conducted to explore the dimensional structure of the data, handling specified missing values, and applying dimensionality reduction with bootstrapping.
#'
#' @format A list containing several components, typically including:
#' \describe{
#'   \item{boot_samples}{A list of bootstrapped samples, each containing the results of the Blackbox transpose analysis for a specific bootstrap iteration.}
#'   \item{original_fit}{The original fit of the Blackbox transpose analysis without bootstrapping.}
#'   \item{convergence}{Convergence diagnostics for each bootstrapped sample, indicating how well the bootstrapping process converged for the specified dimensions.}
#'   \item{errors}{Error metrics or diagnostics generated during the bootstrapping process.}
#'   \item{summary}{A summary of the bootstrapped analysis, including averaged results and variability across bootstrapped samples.}
#' }
#'
#' @details
#' The `outbbt` object is generated by performing a bootstrapped Blackbox transpose analysis on the `rankings` dataset, with dimensionality reduction to 3 dimensions. Missing data values specified as `77`, `88`, and `89` were handled during the analysis. The process was repeated `R=5` times to generate a robust estimate of the dimensional structure of the data.
#'
#' @source Generated by the `boot.blackbox_transpose` function applied to the `rankings` dataset.
#'
#' @usage outbbt
#'
#' @keywords datasets
#' @docType data
#' @name outbbt
NULL


#' @encoding UTF-8
#' @title Roll Call Data from the French Fourth Republic
#' @description Roll call voting data from the French Fourth Republic, as analyzed by 
#' Rosenthal and Voeten (2004). This dataset was used to estimate a 
#' party-switcher model where a separate ideal point is estimated each time 
#' a legislator changes party affiliation.
#' @format A data frame with 1,416 rows and 2,177 columns. The first 5 columns 
#'   contain legislator information, and columns 6-2177 contain roll call votes:
#' \describe{
#'   \item{CASEID}{Unique identifier for each deputy's party affiliation (integer)}
#'   \item{MID}{Unique ID for each deputy, constant across party switches (integer)}
#'   \item{NAME}{Name of the deputy/legislator (character)}
#'   \item{PAR}{Party affiliation of the deputy (character)}
#'   \item{PARSEQ}{Sequence number of party affiliation for deputies who switched parties (integer)}
#'   \item{V1001}{Roll call vote (numeric)}
#'   \item{V1002}{Roll call vote (numeric)}
#'   \item{V1003}{Roll call vote (numeric)}
#'   \item{V1004}{Roll call vote (numeric)}
#'   \item{V1005}{Roll call vote (numeric)}
#'   \item{V1006}{Roll call vote (numeric)}
#'   \item{V1007}{Roll call vote (numeric)}
#'   \item{V1008}{Roll call vote (numeric)}
#'   \item{V1009}{Roll call vote (numeric)}
#'   \item{V1010}{Roll call vote (numeric)}
#'   \item{V1011}{Roll call vote (numeric)}
#'   \item{V1012}{Roll call vote (numeric)}
#'   \item{V1013}{Roll call vote (numeric)}
#'   \item{V1014}{Roll call vote (numeric)}
#'   \item{V1015}{Roll call vote (numeric)}
#'   \item{V1016}{Roll call vote (numeric)}
#'   \item{V1017}{Roll call vote (numeric)}
#'   \item{V1018}{Roll call vote (numeric)}
#'   \item{V1019}{Roll call vote (numeric)}
#'   \item{V1020}{Roll call vote (numeric)}
#'   \item{V1021}{Roll call vote (numeric)}
#'   \item{V1022}{Roll call vote (numeric)}
#'   \item{V1023}{Roll call vote (numeric)}
#'   \item{V1024}{Roll call vote (numeric)}
#'   \item{V1025}{Roll call vote (numeric)}
#'   \item{V1026}{Roll call vote (numeric)}
#'   \item{V1027}{Roll call vote (numeric)}
#'   \item{V1028}{Roll call vote (numeric)}
#'   \item{V1029}{Roll call vote (numeric)}
#'   \item{V1030}{Roll call vote (numeric)}
#'   \item{V1031}{Roll call vote (numeric)}
#'   \item{V1032}{Roll call vote (numeric)}
#'   \item{V1033}{Roll call vote (numeric)}
#'   \item{V1034}{Roll call vote (numeric)}
#'   \item{V1035}{Roll call vote (numeric)}
#'   \item{V1036}{Roll call vote (numeric)}
#'   \item{V1037}{Roll call vote (numeric)}
#'   \item{V1038}{Roll call vote (numeric)}
#'   \item{V1039}{Roll call vote (numeric)}
#'   \item{V1040}{Roll call vote (numeric)}
#'   \item{V1041}{Roll call vote (numeric)}
#'   \item{V1042}{Roll call vote (numeric)}
#'   \item{V1043}{Roll call vote (numeric)}
#'   \item{V1044}{Roll call vote (numeric)}
#'   \item{V1045}{Roll call vote (numeric)}
#'   \item{V1046}{Roll call vote (numeric)}
#'   \item{V1047}{Roll call vote (numeric)}
#'   \item{V1048}{Roll call vote (numeric)}
#'   \item{V1049}{Roll call vote (numeric)}
#'   \item{V1050}{Roll call vote (numeric)}
#'   \item{V1051}{Roll call vote (numeric)}
#'   \item{V1052}{Roll call vote (numeric)}
#'   \item{V1053}{Roll call vote (numeric)}
#'   \item{V1054}{Roll call vote (numeric)}
#'   \item{V1055}{Roll call vote (numeric)}
#'   \item{V1056}{Roll call vote (numeric)}
#'   \item{V1057}{Roll call vote (numeric)}
#'   \item{V1058}{Roll call vote (numeric)}
#'   \item{V1059}{Roll call vote (numeric)}
#'   \item{V1060}{Roll call vote (numeric)}
#'   \item{V1061}{Roll call vote (numeric)}
#'   \item{V1062}{Roll call vote (numeric)}
#'   \item{V1063}{Roll call vote (numeric)}
#'   \item{V1064}{Roll call vote (numeric)}
#'   \item{V1065}{Roll call vote (numeric)}
#'   \item{V1066}{Roll call vote (numeric)}
#'   \item{V1067}{Roll call vote (numeric)}
#'   \item{V1068}{Roll call vote (numeric)}
#'   \item{V1069}{Roll call vote (numeric)}
#'   \item{V1070}{Roll call vote (numeric)}
#'   \item{V1071}{Roll call vote (numeric)}
#'   \item{V1072}{Roll call vote (numeric)}
#'   \item{V1073}{Roll call vote (numeric)}
#'   \item{V1074}{Roll call vote (numeric)}
#'   \item{V1075}{Roll call vote (numeric)}
#'   \item{V1076}{Roll call vote (numeric)}
#'   \item{V1077}{Roll call vote (numeric)}
#'   \item{V1078}{Roll call vote (numeric)}
#'   \item{V1079}{Roll call vote (numeric)}
#'   \item{V1080}{Roll call vote (numeric)}
#'   \item{V1081}{Roll call vote (numeric)}
#'   \item{V1082}{Roll call vote (numeric)}
#'   \item{V1083}{Roll call vote (numeric)}
#'   \item{V1084}{Roll call vote (numeric)}
#'   \item{V1085}{Roll call vote (numeric)}
#'   \item{V1086}{Roll call vote (numeric)}
#'   \item{V1087}{Roll call vote (numeric)}
#'   \item{V1088}{Roll call vote (numeric)}
#'   \item{V1089}{Roll call vote (numeric)}
#'   \item{V1090}{Roll call vote (numeric)}
#'   \item{V1091}{Roll call vote (numeric)}
#'   \item{V1092}{Roll call vote (numeric)}
#'   \item{V1093}{Roll call vote (numeric)}
#'   \item{V1094}{Roll call vote (numeric)}
#'   \item{V1095}{Roll call vote (numeric)}
#'   \item{V1096}{Roll call vote (numeric)}
#'   \item{V1097}{Roll call vote (numeric)}
#'   \item{V1098}{Roll call vote (numeric)}
#'   \item{V1099}{Roll call vote (numeric)}
#'   \item{V1100}{Roll call vote (numeric)}
#'   \item{V1101}{Roll call vote (numeric)}
#'   \item{V1102}{Roll call vote (numeric)}
#'   \item{V1103}{Roll call vote (numeric)}
#'   \item{V1104}{Roll call vote (numeric)}
#'   \item{V1105}{Roll call vote (numeric)}
#'   \item{V1106}{Roll call vote (numeric)}
#'   \item{V1107}{Roll call vote (numeric)}
#'   \item{V1108}{Roll call vote (numeric)}
#'   \item{V1109}{Roll call vote (numeric)}
#'   \item{V1110}{Roll call vote (numeric)}
#'   \item{V1111}{Roll call vote (numeric)}
#'   \item{V1112}{Roll call vote (numeric)}
#'   \item{V1113}{Roll call vote (numeric)}
#'   \item{V1114}{Roll call vote (numeric)}
#'   \item{V1115}{Roll call vote (numeric)}
#'   \item{V1116}{Roll call vote (numeric)}
#'   \item{V1117}{Roll call vote (numeric)}
#'   \item{V1118}{Roll call vote (numeric)}
#'   \item{V1119}{Roll call vote (numeric)}
#'   \item{V1120}{Roll call vote (numeric)}
#'   \item{V1121}{Roll call vote (numeric)}
#'   \item{V1122}{Roll call vote (numeric)}
#'   \item{V1123}{Roll call vote (numeric)}
#'   \item{V1124}{Roll call vote (numeric)}
#'   \item{V1125}{Roll call vote (numeric)}
#'   \item{V1126}{Roll call vote (numeric)}
#'   \item{V1127}{Roll call vote (numeric)}
#'   \item{V1128}{Roll call vote (numeric)}
#'   \item{V1129}{Roll call vote (numeric)}
#'   \item{V1130}{Roll call vote (numeric)}
#'   \item{V1131}{Roll call vote (numeric)}
#'   \item{V1132}{Roll call vote (numeric)}
#'   \item{V1133}{Roll call vote (numeric)}
#'   \item{V1134}{Roll call vote (numeric)}
#'   \item{V1135}{Roll call vote (numeric)}
#'   \item{V1136}{Roll call vote (numeric)}
#'   \item{V1137}{Roll call vote (numeric)}
#'   \item{V1138}{Roll call vote (numeric)}
#'   \item{V1139}{Roll call vote (numeric)}
#'   \item{V1140}{Roll call vote (numeric)}
#'   \item{V1141}{Roll call vote (numeric)}
#'   \item{V1142}{Roll call vote (numeric)}
#'   \item{V1143}{Roll call vote (numeric)}
#'   \item{V1144}{Roll call vote (numeric)}
#'   \item{V1145}{Roll call vote (numeric)}
#'   \item{V1146}{Roll call vote (numeric)}
#'   \item{V1147}{Roll call vote (numeric)}
#'   \item{V1148}{Roll call vote (numeric)}
#'   \item{V1149}{Roll call vote (numeric)}
#'   \item{V1150}{Roll call vote (numeric)}
#'   \item{V1151}{Roll call vote (numeric)}
#'   \item{V1152}{Roll call vote (numeric)}
#'   \item{V1153}{Roll call vote (numeric)}
#'   \item{V1154}{Roll call vote (numeric)}
#'   \item{V1155}{Roll call vote (numeric)}
#'   \item{V1156}{Roll call vote (numeric)}
#'   \item{V1157}{Roll call vote (numeric)}
#'   \item{V1158}{Roll call vote (numeric)}
#'   \item{V1159}{Roll call vote (numeric)}
#'   \item{V1160}{Roll call vote (numeric)}
#'   \item{V1161}{Roll call vote (numeric)}
#'   \item{V1162}{Roll call vote (numeric)}
#'   \item{V1163}{Roll call vote (numeric)}
#'   \item{V1164}{Roll call vote (numeric)}
#'   \item{V1165}{Roll call vote (numeric)}
#'   \item{V1166}{Roll call vote (numeric)}
#'   \item{V1167}{Roll call vote (numeric)}
#'   \item{V1168}{Roll call vote (numeric)}
#'   \item{V1169}{Roll call vote (numeric)}
#'   \item{V1170}{Roll call vote (numeric)}
#'   \item{V1171}{Roll call vote (numeric)}
#'   \item{V1172}{Roll call vote (numeric)}
#'   \item{V1173}{Roll call vote (numeric)}
#'   \item{V1174}{Roll call vote (numeric)}
#'   \item{V1175}{Roll call vote (numeric)}
#'   \item{V1176}{Roll call vote (numeric)}
#'   \item{V1177}{Roll call vote (numeric)}
#'   \item{V1178}{Roll call vote (numeric)}
#'   \item{V1179}{Roll call vote (numeric)}
#'   \item{V1180}{Roll call vote (numeric)}
#'   \item{V1181}{Roll call vote (numeric)}
#'   \item{V1182}{Roll call vote (numeric)}
#'   \item{V1183}{Roll call vote (numeric)}
#'   \item{V1184}{Roll call vote (numeric)}
#'   \item{V1185}{Roll call vote (numeric)}
#'   \item{V1186}{Roll call vote (numeric)}
#'   \item{V1187}{Roll call vote (numeric)}
#'   \item{V1188}{Roll call vote (numeric)}
#'   \item{V1189}{Roll call vote (numeric)}
#'   \item{V1190}{Roll call vote (numeric)}
#'   \item{V1191}{Roll call vote (numeric)}
#'   \item{V1192}{Roll call vote (numeric)}
#'   \item{V1193}{Roll call vote (numeric)}
#'   \item{V1194}{Roll call vote (numeric)}
#'   \item{V1195}{Roll call vote (numeric)}
#'   \item{V1196}{Roll call vote (numeric)}
#'   \item{V1197}{Roll call vote (numeric)}
#'   \item{V1198}{Roll call vote (numeric)}
#'   \item{V1199}{Roll call vote (numeric)}
#'   \item{V1200}{Roll call vote (numeric)}
#'   \item{V1201}{Roll call vote (numeric)}
#'   \item{V1202}{Roll call vote (numeric)}
#'   \item{V1203}{Roll call vote (numeric)}
#'   \item{V1204}{Roll call vote (numeric)}
#'   \item{V1205}{Roll call vote (numeric)}
#'   \item{V1206}{Roll call vote (numeric)}
#'   \item{V1207}{Roll call vote (numeric)}
#'   \item{V1208}{Roll call vote (numeric)}
#'   \item{V1209}{Roll call vote (numeric)}
#'   \item{V1210}{Roll call vote (numeric)}
#'   \item{V1211}{Roll call vote (numeric)}
#'   \item{V1212}{Roll call vote (numeric)}
#'   \item{V1213}{Roll call vote (numeric)}
#'   \item{V1214}{Roll call vote (numeric)}
#'   \item{V1215}{Roll call vote (numeric)}
#'   \item{V1216}{Roll call vote (numeric)}
#'   \item{V1217}{Roll call vote (numeric)}
#'   \item{V1218}{Roll call vote (numeric)}
#'   \item{V1219}{Roll call vote (numeric)}
#'   \item{V1220}{Roll call vote (numeric)}
#'   \item{V1221}{Roll call vote (numeric)}
#'   \item{V1222}{Roll call vote (numeric)}
#'   \item{V1223}{Roll call vote (numeric)}
#'   \item{V1224}{Roll call vote (numeric)}
#'   \item{V1225}{Roll call vote (numeric)}
#'   \item{V1226}{Roll call vote (numeric)}
#'   \item{V1227}{Roll call vote (numeric)}
#'   \item{V1228}{Roll call vote (numeric)}
#'   \item{V1229}{Roll call vote (numeric)}
#'   \item{V1230}{Roll call vote (numeric)}
#'   \item{V1231}{Roll call vote (numeric)}
#'   \item{V1232}{Roll call vote (numeric)}
#'   \item{V1233}{Roll call vote (numeric)}
#'   \item{V1234}{Roll call vote (numeric)}
#'   \item{V1235}{Roll call vote (numeric)}
#'   \item{V1236}{Roll call vote (numeric)}
#'   \item{V1237}{Roll call vote (numeric)}
#'   \item{V1238}{Roll call vote (numeric)}
#'   \item{V1239}{Roll call vote (numeric)}
#'   \item{V1240}{Roll call vote (numeric)}
#'   \item{V1241}{Roll call vote (numeric)}
#'   \item{V1242}{Roll call vote (numeric)}
#'   \item{V1243}{Roll call vote (numeric)}
#'   \item{V1244}{Roll call vote (numeric)}
#'   \item{V1245}{Roll call vote (numeric)}
#'   \item{V1246}{Roll call vote (numeric)}
#'   \item{V1247}{Roll call vote (numeric)}
#'   \item{V1248}{Roll call vote (numeric)}
#'   \item{V1249}{Roll call vote (numeric)}
#'   \item{V1250}{Roll call vote (numeric)}
#'   \item{V1251}{Roll call vote (numeric)}
#'   \item{V1252}{Roll call vote (numeric)}
#'   \item{V1253}{Roll call vote (numeric)}
#'   \item{V1254}{Roll call vote (numeric)}
#'   \item{V1255}{Roll call vote (numeric)}
#'   \item{V1256}{Roll call vote (numeric)}
#'   \item{V1257}{Roll call vote (numeric)}
#'   \item{V1258}{Roll call vote (numeric)}
#'   \item{V1259}{Roll call vote (numeric)}
#'   \item{V1260}{Roll call vote (numeric)}
#'   \item{V1261}{Roll call vote (numeric)}
#'   \item{V1262}{Roll call vote (numeric)}
#'   \item{V1263}{Roll call vote (numeric)}
#'   \item{V1264}{Roll call vote (numeric)}
#'   \item{V1265}{Roll call vote (numeric)}
#'   \item{V1266}{Roll call vote (numeric)}
#'   \item{V1267}{Roll call vote (numeric)}
#'   \item{V1268}{Roll call vote (numeric)}
#'   \item{V1269}{Roll call vote (numeric)}
#'   \item{V1270}{Roll call vote (numeric)}
#'   \item{V1271}{Roll call vote (numeric)}
#'   \item{V1272}{Roll call vote (numeric)}
#'   \item{V1273}{Roll call vote (numeric)}
#'   \item{V1274}{Roll call vote (numeric)}
#'   \item{V1275}{Roll call vote (numeric)}
#'   \item{V1276}{Roll call vote (numeric)}
#'   \item{V1277}{Roll call vote (numeric)}
#'   \item{V1278}{Roll call vote (numeric)}
#'   \item{V1279}{Roll call vote (numeric)}
#'   \item{V1280}{Roll call vote (numeric)}
#'   \item{V1281}{Roll call vote (numeric)}
#'   \item{V1282}{Roll call vote (numeric)}
#'   \item{V1283}{Roll call vote (numeric)}
#'   \item{V1284}{Roll call vote (numeric)}
#'   \item{V1285}{Roll call vote (numeric)}
#'   \item{V1286}{Roll call vote (numeric)}
#'   \item{V1287}{Roll call vote (numeric)}
#'   \item{V1288}{Roll call vote (numeric)}
#'   \item{V1289}{Roll call vote (numeric)}
#'   \item{V1290}{Roll call vote (numeric)}
#'   \item{V1291}{Roll call vote (numeric)}
#'   \item{V1292}{Roll call vote (numeric)}
#'   \item{V1293}{Roll call vote (numeric)}
#'   \item{V1294}{Roll call vote (numeric)}
#'   \item{V1295}{Roll call vote (numeric)}
#'   \item{V1296}{Roll call vote (numeric)}
#'   \item{V1297}{Roll call vote (numeric)}
#'   \item{V1298}{Roll call vote (numeric)}
#'   \item{V1299}{Roll call vote (numeric)}
#'   \item{V1300}{Roll call vote (numeric)}
#'   \item{V1301}{Roll call vote (numeric)}
#'   \item{V1302}{Roll call vote (numeric)}
#'   \item{V1303}{Roll call vote (numeric)}
#'   \item{V1304}{Roll call vote (numeric)}
#'   \item{V1305}{Roll call vote (numeric)}
#'   \item{V1306}{Roll call vote (numeric)}
#'   \item{V1307}{Roll call vote (numeric)}
#'   \item{V1308}{Roll call vote (numeric)}
#'   \item{V1309}{Roll call vote (numeric)}
#'   \item{V1310}{Roll call vote (numeric)}
#'   \item{V1311}{Roll call vote (numeric)}
#'   \item{V1312}{Roll call vote (numeric)}
#'   \item{V1313}{Roll call vote (numeric)}
#'   \item{V1314}{Roll call vote (numeric)}
#'   \item{V1315}{Roll call vote (numeric)}
#'   \item{V1316}{Roll call vote (numeric)}
#'   \item{V1317}{Roll call vote (numeric)}
#'   \item{V1318}{Roll call vote (numeric)}
#'   \item{V1319}{Roll call vote (numeric)}
#'   \item{V1320}{Roll call vote (numeric)}
#'   \item{V1321}{Roll call vote (numeric)}
#'   \item{V1322}{Roll call vote (numeric)}
#'   \item{V1323}{Roll call vote (numeric)}
#'   \item{V1324}{Roll call vote (numeric)}
#'   \item{V1325}{Roll call vote (numeric)}
#'   \item{V1326}{Roll call vote (numeric)}
#'   \item{V1327}{Roll call vote (numeric)}
#'   \item{V1328}{Roll call vote (numeric)}
#'   \item{V1329}{Roll call vote (numeric)}
#'   \item{V1330}{Roll call vote (numeric)}
#'   \item{V1331}{Roll call vote (numeric)}
#'   \item{V1332}{Roll call vote (numeric)}
#'   \item{V1333}{Roll call vote (numeric)}
#'   \item{V1334}{Roll call vote (numeric)}
#'   \item{V1335}{Roll call vote (numeric)}
#'   \item{V1336}{Roll call vote (numeric)}
#'   \item{V1337}{Roll call vote (numeric)}
#'   \item{V1338}{Roll call vote (numeric)}
#'   \item{V1339}{Roll call vote (numeric)}
#'   \item{V1340}{Roll call vote (numeric)}
#'   \item{V1341}{Roll call vote (numeric)}
#'   \item{V1342}{Roll call vote (numeric)}
#'   \item{V1343}{Roll call vote (numeric)}
#'   \item{V1344}{Roll call vote (numeric)}
#'   \item{V1345}{Roll call vote (numeric)}
#'   \item{V1346}{Roll call vote (numeric)}
#'   \item{V1347}{Roll call vote (numeric)}
#'   \item{V1348}{Roll call vote (numeric)}
#'   \item{V1349}{Roll call vote (numeric)}
#'   \item{V1350}{Roll call vote (numeric)}
#'   \item{V1351}{Roll call vote (numeric)}
#'   \item{V1352}{Roll call vote (numeric)}
#'   \item{V1353}{Roll call vote (numeric)}
#'   \item{V1354}{Roll call vote (numeric)}
#'   \item{V1355}{Roll call vote (numeric)}
#'   \item{V1356}{Roll call vote (numeric)}
#'   \item{V1357}{Roll call vote (numeric)}
#'   \item{V1358}{Roll call vote (numeric)}
#'   \item{V1359}{Roll call vote (numeric)}
#'   \item{V1360}{Roll call vote (numeric)}
#'   \item{V1361}{Roll call vote (numeric)}
#'   \item{V1362}{Roll call vote (numeric)}
#'   \item{V1363}{Roll call vote (numeric)}
#'   \item{V1364}{Roll call vote (numeric)}
#'   \item{V1365}{Roll call vote (numeric)}
#'   \item{V2001}{Roll call vote (numeric)}
#'   \item{V2002}{Roll call vote (numeric)}
#'   \item{V2003}{Roll call vote (numeric)}
#'   \item{V2004}{Roll call vote (numeric)}
#'   \item{V2005}{Roll call vote (numeric)}
#'   \item{V2006}{Roll call vote (numeric)}
#'   \item{V2007}{Roll call vote (numeric)}
#'   \item{V2008}{Roll call vote (numeric)}
#'   \item{V2009}{Roll call vote (numeric)}
#'   \item{V2010}{Roll call vote (numeric)}
#'   \item{V2011}{Roll call vote (numeric)}
#'   \item{V2012}{Roll call vote (numeric)}
#'   \item{V2013}{Roll call vote (numeric)}
#'   \item{V2014}{Roll call vote (numeric)}
#'   \item{V2015}{Roll call vote (numeric)}
#'   \item{V2016}{Roll call vote (numeric)}
#'   \item{V2017}{Roll call vote (numeric)}
#'   \item{V2018}{Roll call vote (numeric)}
#'   \item{V2019}{Roll call vote (numeric)}
#'   \item{V2020}{Roll call vote (numeric)}
#'   \item{V2021}{Roll call vote (numeric)}
#'   \item{V2022}{Roll call vote (numeric)}
#'   \item{V2023}{Roll call vote (numeric)}
#'   \item{V2024}{Roll call vote (numeric)}
#'   \item{V2025}{Roll call vote (numeric)}
#'   \item{V2026}{Roll call vote (numeric)}
#'   \item{V2027}{Roll call vote (numeric)}
#'   \item{V2028}{Roll call vote (numeric)}
#'   \item{V2029}{Roll call vote (numeric)}
#'   \item{V2030}{Roll call vote (numeric)}
#'   \item{V2031}{Roll call vote (numeric)}
#'   \item{V2032}{Roll call vote (numeric)}
#'   \item{V2033}{Roll call vote (numeric)}
#'   \item{V2034}{Roll call vote (numeric)}
#'   \item{V2035}{Roll call vote (numeric)}
#'   \item{V2036}{Roll call vote (numeric)}
#'   \item{V2037}{Roll call vote (numeric)}
#'   \item{V2038}{Roll call vote (numeric)}
#'   \item{V2039}{Roll call vote (numeric)}
#'   \item{V2040}{Roll call vote (numeric)}
#'   \item{V2041}{Roll call vote (numeric)}
#'   \item{V2042}{Roll call vote (numeric)}
#'   \item{V2043}{Roll call vote (numeric)}
#'   \item{V2044}{Roll call vote (numeric)}
#'   \item{V2045}{Roll call vote (numeric)}
#'   \item{V2046}{Roll call vote (numeric)}
#'   \item{V2047}{Roll call vote (numeric)}
#'   \item{V2048}{Roll call vote (numeric)}
#'   \item{V2049}{Roll call vote (numeric)}
#'   \item{V2050}{Roll call vote (numeric)}
#'   \item{V2051}{Roll call vote (numeric)}
#'   \item{V2052}{Roll call vote (numeric)}
#'   \item{V2053}{Roll call vote (numeric)}
#'   \item{V2054}{Roll call vote (numeric)}
#'   \item{V2055}{Roll call vote (numeric)}
#'   \item{V2056}{Roll call vote (numeric)}
#'   \item{V2057}{Roll call vote (numeric)}
#'   \item{V2058}{Roll call vote (numeric)}
#'   \item{V2059}{Roll call vote (numeric)}
#'   \item{V2060}{Roll call vote (numeric)}
#'   \item{V2061}{Roll call vote (numeric)}
#'   \item{V2062}{Roll call vote (numeric)}
#'   \item{V2063}{Roll call vote (numeric)}
#'   \item{V2064}{Roll call vote (numeric)}
#'   \item{V2065}{Roll call vote (numeric)}
#'   \item{V2066}{Roll call vote (numeric)}
#'   \item{V2067}{Roll call vote (numeric)}
#'   \item{V2068}{Roll call vote (numeric)}
#'   \item{V2069}{Roll call vote (numeric)}
#'   \item{V2070}{Roll call vote (numeric)}
#'   \item{V2071}{Roll call vote (numeric)}
#'   \item{V2072}{Roll call vote (numeric)}
#'   \item{V2073}{Roll call vote (numeric)}
#'   \item{V2074}{Roll call vote (numeric)}
#'   \item{V2075}{Roll call vote (numeric)}
#'   \item{V2076}{Roll call vote (numeric)}
#'   \item{V2077}{Roll call vote (numeric)}
#'   \item{V2078}{Roll call vote (numeric)}
#'   \item{V2079}{Roll call vote (numeric)}
#'   \item{V2080}{Roll call vote (numeric)}
#'   \item{V2081}{Roll call vote (numeric)}
#'   \item{V2082}{Roll call vote (numeric)}
#'   \item{V2083}{Roll call vote (numeric)}
#'   \item{V2084}{Roll call vote (numeric)}
#'   \item{V2085}{Roll call vote (numeric)}
#'   \item{V2086}{Roll call vote (numeric)}
#'   \item{V2087}{Roll call vote (numeric)}
#'   \item{V2088}{Roll call vote (numeric)}
#'   \item{V2089}{Roll call vote (numeric)}
#'   \item{V2090}{Roll call vote (numeric)}
#'   \item{V2091}{Roll call vote (numeric)}
#'   \item{V2092}{Roll call vote (numeric)}
#'   \item{V2093}{Roll call vote (numeric)}
#'   \item{V2094}{Roll call vote (numeric)}
#'   \item{V2095}{Roll call vote (numeric)}
#'   \item{V2096}{Roll call vote (numeric)}
#'   \item{V2097}{Roll call vote (numeric)}
#'   \item{V2098}{Roll call vote (numeric)}
#'   \item{V2099}{Roll call vote (numeric)}
#'   \item{V2100}{Roll call vote (numeric)}
#'   \item{V2101}{Roll call vote (numeric)}
#'   \item{V2102}{Roll call vote (numeric)}
#'   \item{V2103}{Roll call vote (numeric)}
#'   \item{V2104}{Roll call vote (numeric)}
#'   \item{V2105}{Roll call vote (numeric)}
#'   \item{V2106}{Roll call vote (numeric)}
#'   \item{V2107}{Roll call vote (numeric)}
#'   \item{V2108}{Roll call vote (numeric)}
#'   \item{V2109}{Roll call vote (numeric)}
#'   \item{V2110}{Roll call vote (numeric)}
#'   \item{V2111}{Roll call vote (numeric)}
#'   \item{V2112}{Roll call vote (numeric)}
#'   \item{V2113}{Roll call vote (numeric)}
#'   \item{V2114}{Roll call vote (numeric)}
#'   \item{V2115}{Roll call vote (numeric)}
#'   \item{V2116}{Roll call vote (numeric)}
#'   \item{V2117}{Roll call vote (numeric)}
#'   \item{V2118}{Roll call vote (numeric)}
#'   \item{V2119}{Roll call vote (numeric)}
#'   \item{V2120}{Roll call vote (numeric)}
#'   \item{V2121}{Roll call vote (numeric)}
#'   \item{V2122}{Roll call vote (numeric)}
#'   \item{V2123}{Roll call vote (numeric)}
#'   \item{V2124}{Roll call vote (numeric)}
#'   \item{V2125}{Roll call vote (numeric)}
#'   \item{V2126}{Roll call vote (numeric)}
#'   \item{V2127}{Roll call vote (numeric)}
#'   \item{V2128}{Roll call vote (numeric)}
#'   \item{V2129}{Roll call vote (numeric)}
#'   \item{V2130}{Roll call vote (numeric)}
#'   \item{V2131}{Roll call vote (numeric)}
#'   \item{V2132}{Roll call vote (numeric)}
#'   \item{V2133}{Roll call vote (numeric)}
#'   \item{V2134}{Roll call vote (numeric)}
#'   \item{V2135}{Roll call vote (numeric)}
#'   \item{V2136}{Roll call vote (numeric)}
#'   \item{V2137}{Roll call vote (numeric)}
#'   \item{V2138}{Roll call vote (numeric)}
#'   \item{V2139}{Roll call vote (numeric)}
#'   \item{V2140}{Roll call vote (numeric)}
#'   \item{V2141}{Roll call vote (numeric)}
#'   \item{V2142}{Roll call vote (numeric)}
#'   \item{V2143}{Roll call vote (numeric)}
#'   \item{V2144}{Roll call vote (numeric)}
#'   \item{V2145}{Roll call vote (numeric)}
#'   \item{V2146}{Roll call vote (numeric)}
#'   \item{V2147}{Roll call vote (numeric)}
#'   \item{V2148}{Roll call vote (numeric)}
#'   \item{V2149}{Roll call vote (numeric)}
#'   \item{V2150}{Roll call vote (numeric)}
#'   \item{V2151}{Roll call vote (numeric)}
#'   \item{V2152}{Roll call vote (numeric)}
#'   \item{V2153}{Roll call vote (numeric)}
#'   \item{V2154}{Roll call vote (numeric)}
#'   \item{V2155}{Roll call vote (numeric)}
#'   \item{V2156}{Roll call vote (numeric)}
#'   \item{V2157}{Roll call vote (numeric)}
#'   \item{V2158}{Roll call vote (numeric)}
#'   \item{V2159}{Roll call vote (numeric)}
#'   \item{V2160}{Roll call vote (numeric)}
#'   \item{V2161}{Roll call vote (numeric)}
#'   \item{V2162}{Roll call vote (numeric)}
#'   \item{V2163}{Roll call vote (numeric)}
#'   \item{V2164}{Roll call vote (numeric)}
#'   \item{V2165}{Roll call vote (numeric)}
#'   \item{V2166}{Roll call vote (numeric)}
#'   \item{V2167}{Roll call vote (numeric)}
#'   \item{V2168}{Roll call vote (numeric)}
#'   \item{V2169}{Roll call vote (numeric)}
#'   \item{V2170}{Roll call vote (numeric)}
#'   \item{V2171}{Roll call vote (numeric)}
#'   \item{V2172}{Roll call vote (numeric)}
#'   \item{V2173}{Roll call vote (numeric)}
#'   \item{V2174}{Roll call vote (numeric)}
#'   \item{V2175}{Roll call vote (numeric)}
#'   \item{V2176}{Roll call vote (numeric)}
#'   \item{V2177}{Roll call vote (numeric)}
#'   \item{V2178}{Roll call vote (numeric)}
#'   \item{V2179}{Roll call vote (numeric)}
#'   \item{V2180}{Roll call vote (numeric)}
#'   \item{V2181}{Roll call vote (numeric)}
#'   \item{V2182}{Roll call vote (numeric)}
#'   \item{V2183}{Roll call vote (numeric)}
#'   \item{V2184}{Roll call vote (numeric)}
#'   \item{V2185}{Roll call vote (numeric)}
#'   \item{V2186}{Roll call vote (numeric)}
#'   \item{V2187}{Roll call vote (numeric)}
#'   \item{V2188}{Roll call vote (numeric)}
#'   \item{V2189}{Roll call vote (numeric)}
#'   \item{V2190}{Roll call vote (numeric)}
#'   \item{V2191}{Roll call vote (numeric)}
#'   \item{V2192}{Roll call vote (numeric)}
#'   \item{V2193}{Roll call vote (numeric)}
#'   \item{V2194}{Roll call vote (numeric)}
#'   \item{V2195}{Roll call vote (numeric)}
#'   \item{V2196}{Roll call vote (numeric)}
#'   \item{V2197}{Roll call vote (numeric)}
#'   \item{V2198}{Roll call vote (numeric)}
#'   \item{V2199}{Roll call vote (numeric)}
#'   \item{V2200}{Roll call vote (numeric)}
#'   \item{V2201}{Roll call vote (numeric)}
#'   \item{V2202}{Roll call vote (numeric)}
#'   \item{V2203}{Roll call vote (numeric)}
#'   \item{V2204}{Roll call vote (numeric)}
#'   \item{V2205}{Roll call vote (numeric)}
#'   \item{V2206}{Roll call vote (numeric)}
#'   \item{V2207}{Roll call vote (numeric)}
#'   \item{V2208}{Roll call vote (numeric)}
#'   \item{V2209}{Roll call vote (numeric)}
#'   \item{V2210}{Roll call vote (numeric)}
#'   \item{V2211}{Roll call vote (numeric)}
#'   \item{V2212}{Roll call vote (numeric)}
#'   \item{V2213}{Roll call vote (numeric)}
#'   \item{V2214}{Roll call vote (numeric)}
#'   \item{V2215}{Roll call vote (numeric)}
#'   \item{V2216}{Roll call vote (numeric)}
#'   \item{V2217}{Roll call vote (numeric)}
#'   \item{V2218}{Roll call vote (numeric)}
#'   \item{V2219}{Roll call vote (numeric)}
#'   \item{V2220}{Roll call vote (numeric)}
#'   \item{V2221}{Roll call vote (numeric)}
#'   \item{V2222}{Roll call vote (numeric)}
#'   \item{V2223}{Roll call vote (numeric)}
#'   \item{V2224}{Roll call vote (numeric)}
#'   \item{V2225}{Roll call vote (numeric)}
#'   \item{V2226}{Roll call vote (numeric)}
#'   \item{V2227}{Roll call vote (numeric)}
#'   \item{V2228}{Roll call vote (numeric)}
#'   \item{V2229}{Roll call vote (numeric)}
#'   \item{V2230}{Roll call vote (numeric)}
#'   \item{V2231}{Roll call vote (numeric)}
#'   \item{V2232}{Roll call vote (numeric)}
#'   \item{V2233}{Roll call vote (numeric)}
#'   \item{V2234}{Roll call vote (numeric)}
#'   \item{V2235}{Roll call vote (numeric)}
#'   \item{V2236}{Roll call vote (numeric)}
#'   \item{V2237}{Roll call vote (numeric)}
#'   \item{V2238}{Roll call vote (numeric)}
#'   \item{V2239}{Roll call vote (numeric)}
#'   \item{V2240}{Roll call vote (numeric)}
#'   \item{V2241}{Roll call vote (numeric)}
#'   \item{V2242}{Roll call vote (numeric)}
#'   \item{V2243}{Roll call vote (numeric)}
#'   \item{V2244}{Roll call vote (numeric)}
#'   \item{V2245}{Roll call vote (numeric)}
#'   \item{V2246}{Roll call vote (numeric)}
#'   \item{V2247}{Roll call vote (numeric)}
#'   \item{V2248}{Roll call vote (numeric)}
#'   \item{V2249}{Roll call vote (numeric)}
#'   \item{V2250}{Roll call vote (numeric)}
#'   \item{V2251}{Roll call vote (numeric)}
#'   \item{V2252}{Roll call vote (numeric)}
#'   \item{V2253}{Roll call vote (numeric)}
#'   \item{V2254}{Roll call vote (numeric)}
#'   \item{V2255}{Roll call vote (numeric)}
#'   \item{V2256}{Roll call vote (numeric)}
#'   \item{V2257}{Roll call vote (numeric)}
#'   \item{V2258}{Roll call vote (numeric)}
#'   \item{V2259}{Roll call vote (numeric)}
#'   \item{V2260}{Roll call vote (numeric)}
#'   \item{V2261}{Roll call vote (numeric)}
#'   \item{V2262}{Roll call vote (numeric)}
#'   \item{V2263}{Roll call vote (numeric)}
#'   \item{V2264}{Roll call vote (numeric)}
#'   \item{V2265}{Roll call vote (numeric)}
#'   \item{V2266}{Roll call vote (numeric)}
#'   \item{V2267}{Roll call vote (numeric)}
#'   \item{V2268}{Roll call vote (numeric)}
#'   \item{V2269}{Roll call vote (numeric)}
#'   \item{V2270}{Roll call vote (numeric)}
#'   \item{V2271}{Roll call vote (numeric)}
#'   \item{V2272}{Roll call vote (numeric)}
#'   \item{V2273}{Roll call vote (numeric)}
#'   \item{V2274}{Roll call vote (numeric)}
#'   \item{V2275}{Roll call vote (numeric)}
#'   \item{V2276}{Roll call vote (numeric)}
#'   \item{V2277}{Roll call vote (numeric)}
#'   \item{V2278}{Roll call vote (numeric)}
#'   \item{V2279}{Roll call vote (numeric)}
#'   \item{V2280}{Roll call vote (numeric)}
#'   \item{V2281}{Roll call vote (numeric)}
#'   \item{V2282}{Roll call vote (numeric)}
#'   \item{V2283}{Roll call vote (numeric)}
#'   \item{V2284}{Roll call vote (numeric)}
#'   \item{V2285}{Roll call vote (numeric)}
#'   \item{V2286}{Roll call vote (numeric)}
#'   \item{V2287}{Roll call vote (numeric)}
#'   \item{V2288}{Roll call vote (numeric)}
#'   \item{V2289}{Roll call vote (numeric)}
#'   \item{V2290}{Roll call vote (numeric)}
#'   \item{V2291}{Roll call vote (numeric)}
#'   \item{V2292}{Roll call vote (numeric)}
#'   \item{V2293}{Roll call vote (numeric)}
#'   \item{V2294}{Roll call vote (numeric)}
#'   \item{V2295}{Roll call vote (numeric)}
#'   \item{V2296}{Roll call vote (numeric)}
#'   \item{V2297}{Roll call vote (numeric)}
#'   \item{V2298}{Roll call vote (numeric)}
#'   \item{V2299}{Roll call vote (numeric)}
#'   \item{V2300}{Roll call vote (numeric)}
#'   \item{V2301}{Roll call vote (numeric)}
#'   \item{V2302}{Roll call vote (numeric)}
#'   \item{V2303}{Roll call vote (numeric)}
#'   \item{V2304}{Roll call vote (numeric)}
#'   \item{V2305}{Roll call vote (numeric)}
#'   \item{V2306}{Roll call vote (numeric)}
#'   \item{V2307}{Roll call vote (numeric)}
#'   \item{V2308}{Roll call vote (numeric)}
#'   \item{V2309}{Roll call vote (numeric)}
#'   \item{V2310}{Roll call vote (numeric)}
#'   \item{V2311}{Roll call vote (numeric)}
#'   \item{V2312}{Roll call vote (numeric)}
#'   \item{V2313}{Roll call vote (numeric)}
#'   \item{V2314}{Roll call vote (numeric)}
#'   \item{V2315}{Roll call vote (numeric)}
#'   \item{V2316}{Roll call vote (numeric)}
#'   \item{V2317}{Roll call vote (numeric)}
#'   \item{V2318}{Roll call vote (numeric)}
#'   \item{V2319}{Roll call vote (numeric)}
#'   \item{V2320}{Roll call vote (numeric)}
#'   \item{V2321}{Roll call vote (numeric)}
#'   \item{V2322}{Roll call vote (numeric)}
#'   \item{V2323}{Roll call vote (numeric)}
#'   \item{V2324}{Roll call vote (numeric)}
#'   \item{V2325}{Roll call vote (numeric)}
#'   \item{V2326}{Roll call vote (numeric)}
#'   \item{V2327}{Roll call vote (numeric)}
#'   \item{V2328}{Roll call vote (numeric)}
#'   \item{V2329}{Roll call vote (numeric)}
#'   \item{V2330}{Roll call vote (numeric)}
#'   \item{V2331}{Roll call vote (numeric)}
#'   \item{V2332}{Roll call vote (numeric)}
#'   \item{V2333}{Roll call vote (numeric)}
#'   \item{V2334}{Roll call vote (numeric)}
#'   \item{V2335}{Roll call vote (numeric)}
#'   \item{V2336}{Roll call vote (numeric)}
#'   \item{V2337}{Roll call vote (numeric)}
#'   \item{V2338}{Roll call vote (numeric)}
#'   \item{V2339}{Roll call vote (numeric)}
#'   \item{V2340}{Roll call vote (numeric)}
#'   \item{V2341}{Roll call vote (numeric)}
#'   \item{V2342}{Roll call vote (numeric)}
#'   \item{V2343}{Roll call vote (numeric)}
#'   \item{V2344}{Roll call vote (numeric)}
#'   \item{V2345}{Roll call vote (numeric)}
#'   \item{V2346}{Roll call vote (numeric)}
#'   \item{V2347}{Roll call vote (numeric)}
#'   \item{V2348}{Roll call vote (numeric)}
#'   \item{V2349}{Roll call vote (numeric)}
#'   \item{V2350}{Roll call vote (numeric)}
#'   \item{V2351}{Roll call vote (numeric)}
#'   \item{V2352}{Roll call vote (numeric)}
#'   \item{V3001}{Roll call vote (numeric)}
#'   \item{V3002}{Roll call vote (numeric)}
#'   \item{V3003}{Roll call vote (numeric)}
#'   \item{V3004}{Roll call vote (numeric)}
#'   \item{V3005}{Roll call vote (numeric)}
#'   \item{V3006}{Roll call vote (numeric)}
#'   \item{V3007}{Roll call vote (numeric)}
#'   \item{V3008}{Roll call vote (numeric)}
#'   \item{V3009}{Roll call vote (numeric)}
#'   \item{V3010}{Roll call vote (numeric)}
#'   \item{V3011}{Roll call vote (numeric)}
#'   \item{V3012}{Roll call vote (numeric)}
#'   \item{V3013}{Roll call vote (numeric)}
#'   \item{V3014}{Roll call vote (numeric)}
#'   \item{V3015}{Roll call vote (numeric)}
#'   \item{V3016}{Roll call vote (numeric)}
#'   \item{V3017}{Roll call vote (numeric)}
#'   \item{V3018}{Roll call vote (numeric)}
#'   \item{V3019}{Roll call vote (numeric)}
#'   \item{V3020}{Roll call vote (numeric)}
#'   \item{V3021}{Roll call vote (numeric)}
#'   \item{V3022}{Roll call vote (numeric)}
#'   \item{V3023}{Roll call vote (numeric)}
#'   \item{V3024}{Roll call vote (numeric)}
#'   \item{V3025}{Roll call vote (numeric)}
#'   \item{V3026}{Roll call vote (numeric)}
#'   \item{V3027}{Roll call vote (numeric)}
#'   \item{V3028}{Roll call vote (numeric)}
#'   \item{V3029}{Roll call vote (numeric)}
#'   \item{V3030}{Roll call vote (numeric)}
#'   \item{V3031}{Roll call vote (numeric)}
#'   \item{V3032}{Roll call vote (numeric)}
#'   \item{V3033}{Roll call vote (numeric)}
#'   \item{V3034}{Roll call vote (numeric)}
#'   \item{V3035}{Roll call vote (numeric)}
#'   \item{V3036}{Roll call vote (numeric)}
#'   \item{V3037}{Roll call vote (numeric)}
#'   \item{V3038}{Roll call vote (numeric)}
#'   \item{V3039}{Roll call vote (numeric)}
#'   \item{V3040}{Roll call vote (numeric)}
#'   \item{V3041}{Roll call vote (numeric)}
#'   \item{V3042}{Roll call vote (numeric)}
#'   \item{V3043}{Roll call vote (numeric)}
#'   \item{V3044}{Roll call vote (numeric)}
#'   \item{V3045}{Roll call vote (numeric)}
#'   \item{V3046}{Roll call vote (numeric)}
#'   \item{V3047}{Roll call vote (numeric)}
#'   \item{V3048}{Roll call vote (numeric)}
#'   \item{V3049}{Roll call vote (numeric)}
#'   \item{V3050}{Roll call vote (numeric)}
#'   \item{V3051}{Roll call vote (numeric)}
#'   \item{V3052}{Roll call vote (numeric)}
#'   \item{V3053}{Roll call vote (numeric)}
#'   \item{V3054}{Roll call vote (numeric)}
#'   \item{V3055}{Roll call vote (numeric)}
#'   \item{V3056}{Roll call vote (numeric)}
#'   \item{V3057}{Roll call vote (numeric)}
#'   \item{V3058}{Roll call vote (numeric)}
#'   \item{V3059}{Roll call vote (numeric)}
#'   \item{V3060}{Roll call vote (numeric)}
#'   \item{V3061}{Roll call vote (numeric)}
#'   \item{V3062}{Roll call vote (numeric)}
#'   \item{V3063}{Roll call vote (numeric)}
#'   \item{V3064}{Roll call vote (numeric)}
#'   \item{V3065}{Roll call vote (numeric)}
#'   \item{V3066}{Roll call vote (numeric)}
#'   \item{V3067}{Roll call vote (numeric)}
#'   \item{V3068}{Roll call vote (numeric)}
#'   \item{V3069}{Roll call vote (numeric)}
#'   \item{V3070}{Roll call vote (numeric)}
#'   \item{V3071}{Roll call vote (numeric)}
#'   \item{V3072}{Roll call vote (numeric)}
#'   \item{V3073}{Roll call vote (numeric)}
#'   \item{V3074}{Roll call vote (numeric)}
#'   \item{V3075}{Roll call vote (numeric)}
#'   \item{V3076}{Roll call vote (numeric)}
#'   \item{V3077}{Roll call vote (numeric)}
#'   \item{V3078}{Roll call vote (numeric)}
#'   \item{V3079}{Roll call vote (numeric)}
#'   \item{V3080}{Roll call vote (numeric)}
#'   \item{V3081}{Roll call vote (numeric)}
#'   \item{V3082}{Roll call vote (numeric)}
#'   \item{V3083}{Roll call vote (numeric)}
#'   \item{V3084}{Roll call vote (numeric)}
#'   \item{V3085}{Roll call vote (numeric)}
#'   \item{V3086}{Roll call vote (numeric)}
#'   \item{V3087}{Roll call vote (numeric)}
#'   \item{V3088}{Roll call vote (numeric)}
#'   \item{V3089}{Roll call vote (numeric)}
#'   \item{V3090}{Roll call vote (numeric)}
#'   \item{V3091}{Roll call vote (numeric)}
#'   \item{V3092}{Roll call vote (numeric)}
#'   \item{V3093}{Roll call vote (numeric)}
#'   \item{V3094}{Roll call vote (numeric)}
#'   \item{V3095}{Roll call vote (numeric)}
#'   \item{V3096}{Roll call vote (numeric)}
#'   \item{V3097}{Roll call vote (numeric)}
#'   \item{V3098}{Roll call vote (numeric)}
#'   \item{V3099}{Roll call vote (numeric)}
#'   \item{V3100}{Roll call vote (numeric)}
#'   \item{V3101}{Roll call vote (numeric)}
#'   \item{V3102}{Roll call vote (numeric)}
#'   \item{V3103}{Roll call vote (numeric)}
#'   \item{V3104}{Roll call vote (numeric)}
#'   \item{V3105}{Roll call vote (numeric)}
#'   \item{V3106}{Roll call vote (numeric)}
#'   \item{V3107}{Roll call vote (numeric)}
#'   \item{V3108}{Roll call vote (numeric)}
#'   \item{V3109}{Roll call vote (numeric)}
#'   \item{V3110}{Roll call vote (numeric)}
#'   \item{V3111}{Roll call vote (numeric)}
#'   \item{V3112}{Roll call vote (numeric)}
#'   \item{V3113}{Roll call vote (numeric)}
#'   \item{V3114}{Roll call vote (numeric)}
#'   \item{V3115}{Roll call vote (numeric)}
#'   \item{V3116}{Roll call vote (numeric)}
#'   \item{V3117}{Roll call vote (numeric)}
#'   \item{V3118}{Roll call vote (numeric)}
#'   \item{V3119}{Roll call vote (numeric)}
#'   \item{V3120}{Roll call vote (numeric)}
#'   \item{V3121}{Roll call vote (numeric)}
#'   \item{V3122}{Roll call vote (numeric)}
#'   \item{V3123}{Roll call vote (numeric)}
#'   \item{V3124}{Roll call vote (numeric)}
#'   \item{V3125}{Roll call vote (numeric)}
#'   \item{V3126}{Roll call vote (numeric)}
#'   \item{V3127}{Roll call vote (numeric)}
#'   \item{V3128}{Roll call vote (numeric)}
#'   \item{V3129}{Roll call vote (numeric)}
#'   \item{V3130}{Roll call vote (numeric)}
#'   \item{V3131}{Roll call vote (numeric)}
#'   \item{V3132}{Roll call vote (numeric)}
#'   \item{V3133}{Roll call vote (numeric)}
#'   \item{V3134}{Roll call vote (numeric)}
#'   \item{V3135}{Roll call vote (numeric)}
#'   \item{V3136}{Roll call vote (numeric)}
#'   \item{V3137}{Roll call vote (numeric)}
#'   \item{V3138}{Roll call vote (numeric)}
#'   \item{V3139}{Roll call vote (numeric)}
#'   \item{V3140}{Roll call vote (numeric)}
#'   \item{V3141}{Roll call vote (numeric)}
#'   \item{V3142}{Roll call vote (numeric)}
#'   \item{V3143}{Roll call vote (numeric)}
#'   \item{V3144}{Roll call vote (numeric)}
#'   \item{V3145}{Roll call vote (numeric)}
#'   \item{V3146}{Roll call vote (numeric)}
#'   \item{V3147}{Roll call vote (numeric)}
#'   \item{V3148}{Roll call vote (numeric)}
#'   \item{V3149}{Roll call vote (numeric)}
#'   \item{V3150}{Roll call vote (numeric)}
#'   \item{V3151}{Roll call vote (numeric)}
#'   \item{V3152}{Roll call vote (numeric)}
#'   \item{V3153}{Roll call vote (numeric)}
#'   \item{V3154}{Roll call vote (numeric)}
#'   \item{V3155}{Roll call vote (numeric)}
#'   \item{V3156}{Roll call vote (numeric)}
#'   \item{V3157}{Roll call vote (numeric)}
#'   \item{V3158}{Roll call vote (numeric)}
#'   \item{V3159}{Roll call vote (numeric)}
#'   \item{V3160}{Roll call vote (numeric)}
#'   \item{V3161}{Roll call vote (numeric)}
#'   \item{V3162}{Roll call vote (numeric)}
#'   \item{V3163}{Roll call vote (numeric)}
#'   \item{V3164}{Roll call vote (numeric)}
#'   \item{V3165}{Roll call vote (numeric)}
#'   \item{V3166}{Roll call vote (numeric)}
#'   \item{V3167}{Roll call vote (numeric)}
#'   \item{V3168}{Roll call vote (numeric)}
#'   \item{V3169}{Roll call vote (numeric)}
#'   \item{V3170}{Roll call vote (numeric)}
#'   \item{V3171}{Roll call vote (numeric)}
#'   \item{V3172}{Roll call vote (numeric)}
#' }
#' @details
#' The dataset includes 1,416 separate ideal points, as each party-switching 
#' deputy has a separate entry for each party affiliation. Rosenthal and 
#' Voeten (2004) found that the latent ideological space remained stable over 
#' the course of the French Fourth Republic, so the roll call data is not 
#' segmented by legislative session.
#' 
#' The first five columns contain legislator-specific variables, and the 
#' remaining 2,172 columns represent roll call votes. Each row represents 
#' either a unique legislator (if they never switched parties) or a unique 
#' legislator-party combination (if they switched parties during their tenure).
#'
#' @source
#' Rosenthal, H., & Voeten, E. (2004). Analyzing Roll Calls with Perfect 
#' Spatial Voting: France 1946-1958. \emph{American Journal of Political Science}, 
#' 48(3), 620-632. \doi{10.1111/j.0092-5853.2004.00094.x}
#' 
#' Original data and documentation: 
#' \url{http://www9.georgetown.edu/faculty/ev42/france.htm}
#'
#' @references
#' Rosenthal, H., & Voeten, E. (2004). Analyzing Roll Calls with Perfect 
#' Spatial Voting: France 1946-1958. \emph{American Journal of Political Science}, 
#' 48(3), 620-632.
#'
#' @usage data(france4)
#'
#' @examples
#' \dontrun{
#' data(france4)
#' 
#' # Dataset dimensions
#' dim(france4)  # 1416 rows, 2177 columns
#' 
#' # View legislator info
#' head(france4[, 1:5])
#' 
#' # Check for party switchers
#' switchers <- france4[france4$PARSEQ > 1, ]
#' nrow(switchers)  # Number of party-switch instances
#' 
#' # Party distribution
#' table(france4$PAR)
#' 
#' # View first few votes
#' head(france4[, 6:10])
#' }
#'
#' @keywords datasets
#' @name france4
#' @docType data
NULL

