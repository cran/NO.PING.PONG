\name{CONVERT_ES}
\alias{CONVERT_ES}
\title{Converts between r, d, g, and OR effect sizes}
\description{Converts effect sizes, including (a) r to d, g, or OR, (b) d or g to r, d, g or OR,
	and (c) OR to r, d, or g.}
\usage{CONVERT_ES(donnes, ES, ES_type_IN='r', ES_var = NULL, 
           totalN = NULL, grp1_n = NULL, grp2_n = NULL,
           gvar_type_OUT = 'd', 
           CI_level_out = 95,
           CI_level_in = 95, CI_in_lb = NULL, CI_in_ub = NULL, 
           verbose = TRUE)}
\arguments{
  \item{donnes}{
  \code{}A dataframe or matrix with effect sizes and corresponding information, such as the  
  effect size variances, sample sizes, confidence intervals.}

  \item{ES}{
  \code{}The name of the column in donnes with the effect sizes.}

  \item{ES_type_IN}{
  \code{}The type of effect sizes in ES. The options are 'r' (the default), 'z' for Fishers z 
  transformation of r, 'd', 'g', and 'OR'.}

  \item{ES_var}{
  \code{}(optional) The name of the column in donnes with the variances of the effect sizes.}

  \item{totalN}{
  \code{}(optional) The name of the column in donnes with the total N for each study.}

  \item{grp1_n}{
  \code{}(optional) The name of the column in donnes with the Ns for group 1.}

  \item{grp2_n}{
  \code{}(optional) The name of the column in donnes with the Ns for group 2.}

  \item{gvar_type_OUT}{
  \code{}(optional) The kind of SMD variance. The options are 'd' or 'g'.}

  \item{CI_level_out}{
  \code{}(optional) The confidence interval for the output (in whole numbers). 
  The default is 95 (for 95 percent).}

  \item{CI_level_in}{
  \code{}(optional) The confidence interval for the input (in whole numbers). 
  The default is 95 (for 95 percent).}

  \item{CI_in_lb}{
  \code{}(optional) The name of the column in donnes with the lower bound confidence intervals 
  for input, if provided (in whole numbers).}

  \item{CI_in_ub}{
  \code{}(optional) The name of the column in donnes with the upper bound confidence intervals 
  for the input, if provided (in whole numbers).}

  \item{verbose}{
  \code{}(optional) Should detailed results be displayed in console? TRUE (default) or FALSE.}
}
\details{
	This function converts r, z (Fishers z transformation of r), d, g, and OR effect sizes to 
	r, z, d, g, and OR effect sizes using
	conventional formulas (Borenstein & Hedges, 2019; Borenstein, Hedges, 
	Higgins, & Rothstein, 2009). The effect size variances and confidence intervals are also 
	computed if sufficient data are provided as input.
	
	When the input effect sizes are d or g values, it is helpful if the group Ns are
	also provided whenever possible.

    The gvar_type_OUT argument provides a choice between d or g effect size variances whenever
    ES_type_OUT is set to g. The reason for this option is that
    authors of published meta-analyses sometimes report d variances when their analyses were
    conducted on g effect sizes. This is presumably not a wise practice, but it also does
    not make much difference in the computed values.
}
\value{A matrix containing the
	following statistics, if they can be computed using the input informaiton:
	the r, z, d, g, and OR effect sizes and their corresponding variances,
	confidence intervals, and the totalNs used in the analyses.
}
\references{
       {Borenstein, M., & Hedges, L. V. (2019). Effect sizes for meta-analysis. 
	    In H. Cooper, L. V., Hedges, & J. C. Valentine (Eds). \emph{The handbook of 
		research synthesis and meta-analysis (pp. 207-244).} (3rd. edition). 
		New York, NY: Russell Sage Foundation.}
\cr\cr {Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). 
	    Converting among effect sizes. In, \emph{Introduction to meta-analysis pp. 45-49.}
	    Chichester, UK: John Wiley & Sons.}
\cr\cr {Valentine, J. C. & Cooper, H. (2003). \emph{Effect size substantive 
		interpretation guidelines: Issues in the interpretation of effect sizes.} 
		Washington, DC: What Works Clearinghouse.}
}
\author{Brian P. O'Connor}
\examples{
# convert d effect sizes to r	
head(data_NPP$Omega3_Depression)
CONVERT_ES(donnes = data_NPP$Omega3_Depression, ES = 'SMD', ES_type_IN='d', ES_var = NULL, 
           grp1_n = 'CN', grp2_n = 'EN')

\donttest{ 
# convert r effect sizes to g	
head(data_NPP$Math_Performance)
CONVERT_ES(donnes = data_NPP$Math_Performance, ES = 'r', ES_type_IN='r', ES_var = NULL,
           totalN = 'N')
}
}
