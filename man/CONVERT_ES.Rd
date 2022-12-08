\name{CONVERT_ES}
\alias{CONVERT_ES}
\title{Converts between r, d, and g effect sizes}
\description{Converts effect sizes, including r to d or g, and d or g to r, d, or g.}
\usage{CONVERT_ES(ES, ES_var = NULL, ES_type_IN='r', ES_type_OUT='g',
           totalN = NULL, grp1_N = NULL, grp2_N = NULL,
           gvar_type_OUT = 'd', verbose = TRUE)}
\arguments{
  \item{ES}{
  \code{}The effect sizes.}

  \item{ES_var}{
  \code{}(optional) The variances of the effect sizes.}

  \item{ES_type_IN}{
  \code{}The type of effect sizes in ES. The options are 'r' (the default), 'd', or 'g'.}

  \item{ES_type_OUT}{
  \code{}The type of effect sizes for the output. The options are 'r' (the default), 'd', and 'g'.}

  \item{totalN}{
  \code{}(optional) The total N for each study in ES.}

  \item{grp1_N}{
  \code{}(optional) The N for group 1.}

  \item{grp2_N}{
  \code{}(optional) The N for group 2.}

  \item{gvar_type_OUT}{
  \code{}(optional) The kind of SMD variance.}

  \item{verbose}{
  \code{}(optional) Should detailed results be displayed in console? TRUE (default) or FALSE}
}
\details{
	This function converts r, d, and g effect sizes to r, d, or g effect sizes using
	conventional formulas (Borenstein & Hedges, 2019; Borenstein, Hedges, 
	Higgins, & Rothstein, 2009). The effect size variances are also computed if 
	sufficient data are provided as input.

    The gvar_type_OUT argument provides a choice between d or g effect size variances whenever
    ES_type_OUT is set to g. The reason for this option is that
    authors of published meta-analyses sometimes report d variances when their analyses were
    conducted on g effect sizes. This is presumably not a wise practice, but it also does
    not make much difference in the computed values.
    }
\value{An object of class "data.frame". The object is a list containing the
	following possible components:
 \item{g}{The g effect sizes.}
 \item{Vg}{The variances of the g effect sizes.}
 \item{d}{The d effect sizes.}
 \item{Vd}{The variances of the d effect sizes.}
 \item{r}{The r effect sizes.}
 \item{Vr}{The variances of the r effect sizes.}
 \item{totalNs}{The totalNs used in the analyses.}
}
\references{
       {Borenstein, M., & Hedges, L. V. (2019). Effect sizes for meta-analysis. 
	    In H. Cooper, L. V., Hedges, & J. C. Valentine (Eds). \emph{The handbook of 
		research synthesis and meta-analysis (pp. 207-244).} (3rd. edition). 
		New York, NY: Russell Sage Foundation.}
\cr\cr {Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). 
	    Converting among effect sizes. In, \emph{Introduction to meta-analysis} pp. 45-49.
	    Chichester, UK: John Wiley & Sons.}
\cr\cr {Valentine, J. C. & Cooper, H. (2003). \emph{Effect size substantive 
		interpretation guidelines: Issues in the interpretation of effect sizes.} 
		Washington, DC: What Works Clearinghouse.}
}
\author{Brian P. O'Connor}
\examples{
# convert d effect sizes to r	
head(data_NPP$Omega3_Depression)
CONVERT_ES(ES = data_NPP$Omega3_Depression$SMD, ES_var = NULL, ES_type_IN='d', ES_type_OUT='r', 
           grp1_N = data_NPP$Omega3_Depression$CN, grp2_N = data_NPP$Omega3_Depression$EN)

\donttest{ 
# convert r effect sizes to g	
head(data_NPP$Math_Performance)
CONVERT_ES(ES = data_NPP$Math_Performance$r, ES_var = NULL, ES_type_IN='r', ES_type_OUT='g', 
           totalN = data_NPP$Math_Performance$N)
}
}
