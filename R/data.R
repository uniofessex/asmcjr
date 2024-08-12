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
#' @name result.france
#' @docType data
NULL


#' Selected Issues Sweden 2010 Dataset
#'
#' This dataset, `issues.sweden`, is a matrix created from the `Sweden2010` dataset, specifically using columns 7 to 56. It contains issue-related data from Sweden's 2010 election study.
#'
#' @format A matrix with rows representing respondents and columns representing different issues or variables.
#' @source Sweden 2010 Election Study
#' @name issues.sweden
#' @docType data
NULL


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
#' @name rankings
#' @docType data
NULL


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
#' @name bamdata
#' @docType data
NULL


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
#' @seealso \code{\link[BAM]{BAM}} for more details on the BAM function and its parameters.
#'
#' @source The BAM model was applied to a dataset `bamdata` with specific settings to generate `bam.france`.
#'
#' @keywords datasets
#' @name bam.france
#' @docType data
NULL


#' Issues Matrix from CDS2000 Dataset
#'
#' This object, `issues`, is a matrix extracted from the `CDS2000` dataset. It contains selected columns that represent various issues.
#'
#' @format A numeric matrix with rows corresponding to observations and columns representing different issues.
#' @source Extracted from the `CDS2000` dataset.
#' @examples
#' \dontrun{
#' data(issues)
#' }
#' @name issues
#' @docType data
NULL


#' Blackbox Analysis Results for Issues Matrix
#'
#' The object `result.repdem` contains the results of applying the `blackbox` function to the `issues` matrix. This analysis was performed to extract dimensions that represent the underlying structure of the issues, specifically for a dataset containing political data.
#'
#' @format An object of class `blackbox` containing the results of the dimensional analysis, including the identified dimensions and other related statistics.
#'
#' @return The `result.repdem` object contains the extracted dimensions, along with other statistics generated by the `blackbox` function.
#' 
#' @examples
#' \dontrun{
#' data(result.repdem)
#' }
#' 
#' @name result.repdem
#' @docType data
NULL

#' Prepare Input Matrix from interest1981 Dataset
#'
#' The `input2` object is a matrix created from the `interest1981` dataset for further analysis. The process involves extracting relevant columns, filtering rows, transforming values, and handling missing data.
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


#' Prepare Input Matrix from interest1981 Dataset
#'
#' The `input2` object is a matrix created from the `interest1981` dataset for further analysis. The process involves extracting relevant columns, filtering rows, transforming values, and handling missing data.
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


#' MDS Solution from mlsmu6
#'
#' The `mlsmu6_out` dataset contains the multidimensional scaling (MDS) solution generated by applying the `mlsmu6` function to a subset of the `interest1981` dataset. The MDS solution is computed using two dimensions and a cutoff value of 5, with data grouped by political party labels.
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


#' ANES Input Data
#'
#' The `anes.input` object is a subset of the American National Election Study (ANES) 1968 dataset. It contains selected variables used as input for further analysis.
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
NULL



#' ANES 1968 Feeling Thermometers and Voting Data
#'
#' The `ANES1968` dataset includes feeling thermometers data and information on whether respondents reported voting in the 1968 elections, as well as their reported presidential vote choice.
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
#' @docType data
NULL

#' ANES 2004 Issue Scales Data
#'
#' The `ANES2004` dataset includes responses to various issue scales from the 2004 American National Election Studies (ANES). The dataset contains respondents' positions on several key political issues, measured on scales with varying ranges.
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
NULL


#' Danish Module of the 2009 European Election Study (EES)
#'
#' The \code{denmarkEES2009} dataset contains data from the Danish module of the 2009 European Election Study (EES). 
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

#' Interest Group Ratings of Members of Congress (1959-1981)
#'
#' The `interest1981` dataset is a subset of a much larger dataset compiled by Keith Poole, containing nearly 200,000 
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


#' 2004 American National Election Study (ANES) Data
#'
#' The `ANES2004_OOC` dataset contains data from the 2004 American National Election Study (ANES). The 2004 ANES asked respondents about their policy preferences on issues ranging from diplomacy and defense spending to government spending and abortion.
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


#' Roll Call Data for U.S. Congress
#'
#' The `rc_ep` object is a roll call dataset compiled by Poole and Rosenthal. This dataset represents roll call votes in the U.S. Congress and has been processed into an object of class `rollcall()` for analysis.
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

#' State of the Union Address Corpus
#'
#' The `SOTUcorpus` dataset contains the text of each presidential State of the Union address since 1790. 
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

#' 2000 Convention Delegate Study (CDS)
#'
#' The `CDS2000` dataset contains data from the 2000 Convention Delegate Study (CDS), which interviewed delegates to the 
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

#' French Module of the 2009 European Election Study (EES)
#'
#' The `franceEES2009` dataset contains data from the French module of the 2009 European Election Study (EES). The EES surveyed 1,000 French citizens, asking them to place themselves and eight major political parties on a 0-10 left-right scale (0 representing the most left-wing position, 10 representing the most right-wing position).
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

#' Mexican Module of the Comparative Study of Electoral Systems (CSES) 2000 and 2006
#'
#' The `mexicoCSES2006` dataset contains data from the 2000 and 2006 Mexican modules of the Comparative Study of Electoral Systems (CSES). In these surveys, Mexican citizens were asked to place the major political parties on an 11-point left-right scale.
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

#' 2008 U.S. Presidential Vote Data
#'
#' The `presvote2008` dataset contains data on the voting behavior in the 2008 U.S. Presidential election, where a vote 
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

#' Roll Call Data from the 108th US House of Representatives (2003-2005)
#'
#' The `hr108` dataset contains data from the 108th US House of Representatives, covering the period from 2003 to 2005. During this session, the House conducted 843 recorded roll call votes, with 440 Representatives serving in the chamber. The roll call matrix omits President George W. Bush.
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

#' Vietnam War Issue Scales from the 1968 National Election Study (NES)
#'
#' The `nes1968_vietnam` dataset contains responses from the 1968 National Election Study (NES) where respondents were asked to place themselves, President Lyndon Johnson, and the three major presidential candidates—Democrat Hubert Humphrey, Republican Richard Nixon, and American Independent George Wallace—on two seven-point issue scales regarding the Vietnam War.
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

#' Urban Unrest Issue Scales from the 1968 National Election Study (NES)
#'
#' The `nes1968_urbanunrest` dataset contains responses from the 1968 National Election Study (NES) where respondents were asked to place themselves, President Lyndon Johnson, and the three major presidential candidates—Democrat Hubert Humphrey, Republican Richard Nixon, and American Independent George Wallace—on two seven-point issue scales regarding urban unrest.
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

#' 2010 Swedish Parliamentary Candidate Survey
#'
#' The `Sweden2010` dataset contains data from the 2010 Swedish Parliamentary Candidate Survey, conducted by the Swedish 
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


#' Roll Call Data from the French Fourth Republic
#'
#' The `france4` dataset contains roll call data from the French Fourth Republic, as analyzed by Rosenthal and Voeten (2004). 
#' This dataset was used to estimate a party-switcher model where a separate ideal point is estimated each time a legislator 
#' changes party affiliation. The dataset includes legislator-specific variables and roll call votes. Rosenthal and Voeten (2004) 
#' found that the latent ideological space remained stable over the course of the French Fourth Republic, so the roll call data 
#' is not segmented by legislative session.
#'
#' @format A data frame with the following variables:
#' \describe{
#'   \item{NAME}{Name of the deputy (legislator).}
#'   \item{MID}{Deputy ID, which remains constant even if the deputy switches party.}
#'   \item{CASEID}{Deputy ID, which is unique for each party affiliation if the deputy switches party.}
#'   \item{PAR}{Party affiliation of the deputy.}
#'   \item{PARSEQ}{Sequence of party affiliation for deputies who switched parties.}
#'   \item{VOTE1}{Result of the first roll call vote (and so on for subsequent votes).}
#'   \item{...}{Other variables representing additional roll call votes.}
#' }
#'
#' @details
#' The dataset includes 1,416 separate ideal points, as each party-switching deputy has a separate entry for each party affiliation. 
#' The first five columns of the dataset are legislator-specific variables: NAME (deputy name), MID (constant deputy ID), CASEID 
#' (unique ID for each party affiliation), PAR (party affiliation), and PARSEQ (sequence number of party affiliations for party-switching deputies). 
#' The remaining columns represent the roll call votes.
#'
#' The dataset and its extensive documentation were originally made available online by Rosenthal and Voeten (2004). 
#' For more information and access to the original data, visit \url{http://www9.georgetown.edu/faculty/ev42/france.htm}.
#'
#' @source
#' Data from Rosenthal and Voeten (2004). The original dataset and documentation can be found online at 
#' \url{http://www9.georgetown.edu/faculty/ev42/france.htm}.
#'
#' @usage data(france4)
#'
#' @examples
#' \dontrun{
#' data(france4)
#' }
#'
#' @keywords datasets
#' @name france4
#' @docType data
NULL


#' Candidate Favorability Ratings from the 2008 ANES
#'
#' The `candidatetherms2008` dataset contains favorability ratings from the 2008 American National Election Study (ANES). 
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

#' CHES EU Dataset: Party Means and Standard Deviations
#'
#' The `ches_eu` dataset contains data used to calculate party means and standard deviations across all CHES experts. 
#' The dataset is based on the 2010 wave of the Chapel Hill Expert Survey (CHES), which includes 118 parties and 224 experts 
#' from 14 member countries of the European Union. Over 160 experts placed all 3 vignette parties, and between 8 and 17 experts 
#' placed each of the actual parties.
#'
#' @format A data frame with the following variables:
#' \describe{
#'   \item{party}{The name or identifier of the political party.}
#'   \item{mean}{The mean placement of the party as calculated from the responses of the experts.}
#'   \item{sd}{The standard deviation of the party placements across the experts.}
#'   \item{country}{The country where the party is based.}
#'   \item{vignette}{Indicator of whether the party is one of the vignette parties (1) or an actual party (0).}
#' }
#'
#' @details
#' The dataset was used in the analysis presented in the book chapter "2.5 Using Anchoring Vignettes" (pages 58-60). 
#' The data helps in understanding how experts from the Chapel Hill Expert Survey (CHES) placed different political parties 
#' on various dimensions, providing a basis for calculating mean positions and standard deviations.
#'
#' The 2010 wave of the CHES consists of 118 parties and 224 experts across 14 European Union member countries. Over 160 experts 
#' placed all three vignette parties, while between 8 and 17 experts placed each of the actual parties.
#'
#' @source
#' Data from the 2010 wave of the Chapel Hill Expert Survey (CHES).
#'
#' @usage data(ches_eu)
#'
#' @examples
#' \dontrun{
#' data(ches_eu)
#' }
#'
#' @keywords datasets
#' @name sub.europe
#' @docType data
NULL

#' Mexican Political Party Positions on Left-Right Scale (2000 & 2006)
#'
#' The `mexicoCSES2000` dataset contains data from the 2000 and 2006 Mexican modules of the Comparative Study of Electoral Systems (CSES).
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


#' Transposed Rankings Data using Blackbox Method
#'
#' The `original` object is created by applying the `blackbox_transpose` function to a dataset of rankings.
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

#' Roll Call Voting Data from the First European Parliament (1979–1984)
#'
#' The `rcv_ep1` dataset contains roll call voting data from the first European Parliament (1979–1984), as assembled by Hix, Noury, and Roland (2006). 
#' The dataset includes information on Members of the European Parliament (MEPs) and their votes on various issues.
#'
#' @format A data frame (or matrix) with the following variables:
#' \describe{
#'   \item{MEPID}{Unique identifier for each Member of the European Parliament (MEP).}
#'   \item{MEPNAME}{Name of the MEP.}
#'   \item{MS}{Country (Member State) of the MEP.}
#'   \item{NP}{National party affiliation of the MEP.}
#'   \item{EPG}{European Parliament group (party) affiliation of the MEP.}
#'   \item{V1}{Result of the first roll call vote (and so on for subsequent votes).}
#'   \item{V2-V886}{Results of the roll call votes from the second to the 886th vote. Each column represents a specific vote, and the value indicates the MEP's vote on that issue.}
#' }
#'
#' @details
#' The first five columns of the dataset are legislator-specific variables: MEPID (ID number), MEPNAME (name), MS (country), NP (national party affiliation), and EPG (European Parliament group affiliation). 
#' The remaining columns (V1 to V886) represent the roll call votes. Each vote is represented by a column, and the value in each cell indicates how the MEP voted on that particular issue.
#'
#' This dataset was used to analyze voting behavior in the European Parliament during its first term. The raw data can be accessed from \url{http://personal.lse.ac.uk/hix/HixNouryRolandEPdata.htm}.
#'
#' @source
#' Data from Hix, Noury, and Roland (2006), available online at \url{http://personal.lse.ac.uk/hix/HixNouryRolandEPdata.htm}.
#'
#' @usage data(rcv_ep1)
#'
#' @examples
#' \dontrun{
#' data(rcv_ep1)
#' }
#'
#' @keywords datasets
#' @name rcv_ep1
#' @docType data
NULL


#' Roll Call Data from the 111th U.S. Senate
#'
#' The `hr111` dataset contains roll call voting data from the 111th U.S. Senate. This dataset is formatted as a `rollcall` object, which is typically used for analyzing voting behavior in legislative bodies.
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

#' Nation Similarity Ratings Dataset
#'
#' The `nation` dataset contains similarity ratings between twelve nations, as collected by Wish (1971). 
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
#' @name nation
#' @docType data
NULL

#' DW-NOMINATE Scores for the U.S. Congress
#'
#' The `rcx` dataset is a matrix containing DW-NOMINATE scores for members of the U.S. Congress. DW-NOMINATE scores are used 
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

#' French Party Placement Data from the 2009 European Election Study (EES)
#'
#' The `french.parties.individuals` dataset contains party placement data from the French module of the 2009 European Election Study (EES). 
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

#' 7th Legislative Yuan Roll Call Data from Taiwan (2008-2012)
#'
#' This dataset contains roll call voting data from the 7th Legislative Yuan (National Congress) of Taiwan.
#' The dataset includes the names of legislators and their corresponding votes on various bills.
#'
#' @format A data frame with the following variables:
#' \describe{
#'   \item{legis.names}{The names of the legislators.}
#'   \item{party}{The political party of each legislator.}
#'   \item{\code{7-1} to \code{7-999}}{Columns representing votes on various bills, where each column corresponds to a specific bill.}
#' }
#'
#' @details
#' The data captures the legislative behavior during the 7th session of the Legislative Yuan of Taiwan, 
#' providing valuable insights into the political dynamics and decision-making processes.
#'
#' @source
#' Yen-Chihe Liao (2024). *Electoral Reform and Fragmented Polarization: New Evidence from Taiwan Legislative Roll Call*. Legislative Studies Quarterly. 
#' Available at: \url{https://onlinelibrary.wiley.com/doi/full/10.1111/lsq.12459}.
#'
#' @usage data(legis_7th_Taiwan)
#'
#' @examples
#' \dontrun{
#' data(legis_7th_Taiwan)
#' }
#'
#' @keywords datasets
#' @name legis_7th_Taiwan
#' @docType data
NULL


#' 90th US Senate Agreement Score Matrix (1967-1968)
#'
#' This dataset contains the agreement score matrix of the 90th US Senate, covering the years 1967-1968.
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
#' @details
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