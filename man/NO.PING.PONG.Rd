\name{NO.PING.PONG}
\alias{NO.PING.PONG}
\title{Incorporating previous findings when evaluating new data}
\description{A function for revealing what happens when previous findings 
	are taken into account when analyzing new data.}
\usage{NO.PING.PONG(donnes, ES_type_IN=NULL, ES_type_OUT='r',
             rma_method='REML',
             Bayes_type = c('Schmidt_Raju', 'generated', 'raw'), 
             prior_type='MA', CI = 95,
             ES = NULL, N = NULL, vi = NULL,
             grp1_mn = NULL, grp1_sd = NULL, grp1_n = NULL, 
             grp2_mn = NULL, grp2_sd = NULL, grp2_n = NULL,
             gvar_type_OUT = 'd', 
             nitt=53000, burnin=3000, thin=10, verbose=TRUE)}
\arguments{
  \item{donnes}{
  \code{}Either (1) a list with raw data for two variables in each of the list elements
  (wherein each list element consists of rows of raw data from an individual study), or,  
  (2) a dataframe or matrix wherein each row has summary data for a single study
  (e.g., has the group means, SDs, & Ns, or d or g or r values, or effect size and variance
  estimates).}

  \item{ES_type_IN}{
  \code{}The type of effect sizes, if donnes consists of study-level data rather than raw data. 
  The options are 'r' (the default), 'd', or 'g'.}

  \item{ES_type_OUT}{
  \code{}The type of effect sizes for the output. The options are 'r' (the default), 'd', and 'g'.}

  \item{rma_method}{
  \code{}(optional) The method option for the rma function from the metafor package: 'A character 
  string specifying whether a fixed- or a random/mixed-effects model should be fitted. 
  A fixed-effects model (with or without moderators) is fitted when using method="FE". 
  Random/mixed-effects models are fitted by setting method equal to one of the following: 
  "DL", "HE", "SJ", "ML", "REML", "EB", "HS", or "GENQ". The default is "REML".'}

  \item{Bayes_type}{
  \code{}(optional) The kind(s) of Bayesian analyses to be conducted. The options are one or all
  of 'Schmidt_Raju', 'generated' and/or 'raw'. The default is all three, i.e., 
  Bayes_type = c('Schmidt_Raju', 'generated', 'raw'). See the Details.}

  \item{prior_type}{
  \code{}(optional) The type of prior data to be used in the updating analyses (both cumulative and
  Bayesian). The options are 'Bayes' or 'MA'. The default is 'MA'. See the Details.}

  \item{CI}{
  \code{}(optional) The confidence interval for the analyses (in whole numbers). The default is 95.}

  \item{ES}{
  \code{}The name of the column in donnes with the effect sizes.}

  \item{N}{
  \code{}(optional) The name of the column in donnes with the total Ns.}

  \item{vi}{
  \code{}(optional) The name of the column in donnes with the variances for the effect sizes.}

  \item{grp1_mn}{
  \code{}(optional) The name of the column in donnes with the mean for group 1.}

  \item{grp1_sd}{
  \code{}(optional) The name of the column in donnes with the standard deviation for group 1.}

  \item{grp1_n}{
  \code{}(optional) The name of the column in donnes with the number of cases in group 1.}

  \item{grp2_mn}{
  \code{}(optional) The name of the column in donnes with the mean for group 2.}

  \item{grp2_sd}{
  \code{}(optional) The name of the column in donnes with the standard deviation for group 2.}

  \item{grp2_n}{
  \code{}(optional) The name of the column in donnes with the number of cases in group 2.}

  \item{gvar_type_OUT}{
  \code{}(optional) The kind of SMD variance.}

  \item{nitt}{
  \code{}(optional) The number of iterations for Bayesian analyses. The default = 53000}

  \item{burnin}{
  \code{}(optional) The burn-in period for Bayesian analyses. The default = 3000}

  \item{thin}{
  \code{}(optional) The thinning interval for the Bayesian analyses. The default = 10}

  \item{verbose}{
  \code{}(optional) Should detailed results be displayed in console? TRUE (default) or FALSE}
}
\details{
	This function reveals what happens when effect size estimates from previous studies 
	are taken into account when evaluating each new dataset in a study sequence. 
	The analyses are conducted for cumulative meta-analyses and for 
	Bayesian data analyses. The findings from these methods can be contrasted with
	those from null hypothesis statistical sigificance testing (NHST), which does
	not take previous findings into account and which often results in a back and forth
	ping-pong of conclusions about a phenomenon.
	
    The analyses are always conducted using correlation coefficient effect sizes. When the
    input data are not correlation coefficients, then the correlation coefficient equivalents
    are computed from the input data for the analyses. 

    When 'd' or 'g' are specified as the type of effect sizes for the output (using the 
    ES_type_OUT argument), then the output is for conversions of correlation coefficient
    effect sizes to 'd' or 'g' values using the conventional formulas for such conversions.
    
    The gvar_type_OUT argument provides a choice between d or g effect size variances whenever
    ES_type_OUT is set to g. The reason for this option is that
    authors of published meta-analyses sometimes report d variances when their analyses were
    conducted on g effect sizes. This is presumably not a wise practice, but it also does
    not make much difference in the computed values.
    
	The function relies on the metafor package for the meta-analyses and on the
	MCMCglmm package for the Bayesian analyses.
	
	If the function is provided with only the effect sizes and corresponding variance
	estimates, then the analyses will be conducted for meta-analyses and NHST. Bayesian analyses
	are not possible. However, any kind of effect size coefficient can be used in this case.
	
	If the function is provided with d or g effect sizes along with the corresponding total 
	sample sizes, but not with the corresponding variances for the effect sizes, then the
	d or g values will be converted to correlation coefficient equivalents and the
	subsequent analyses will be for r effect sizes. Bayesian analyses are possible in
	this case.
	
	For the "Bayes_type" argument: The Schmidt_Raju method is the computationally
	fastest method of conducting the Bayesian analyses, but it is based on correlation
	coefficient effect sizes, or the correlation coefficient equivalents. 
	The Bayesian analyses can also be conducted using raw
	data for the sequence of studies, or using data that are generated to have the same
	effect size as those in a raw dataset. The raw data option is generally the better approach
	but it is only feasible when the raw data are available. By default, the function
	will run the analyses for the Schmidt_Raju method, for generated data, and for 
	raw data whenever possible.
	
	For the "prior_type" argument: The priors used in the updating analyses (both 
	cumulative and Bayesian) can be the effect size and standard error from previous
	studies based on a meta-analysis or based on Bayesian analyses. 
	The options are 'Bayes' or 'MA', and 'MA' is the default.
	
	By default, the function provides output for cumulative meta-analyses,
    for Bayesian data analyses, and for NHST.
    
    Consistency and agreement rates are computed for the NHST analyses, for the 
    updating MAs, and for the Bayesian analyses. The consistency rate is the 
    proportion of times that the most common conclusion is reached for a pool of 
    effect sizes. Three conclusions are possible for each effect size: 
    a positive effect, a negative effect, and no effect. The signs of the effect 
    sizes and the possible inclusion of a zero value in a confidence interval are 
    used to make these categorizations (e.g., a negative effect conclusion is 
    when a negative effect size has a confidence interval that does not include zero). 
    The number of times each of the three possible conclusions occurs for a pool 
    of effect sizes is counted, and the consistency rate is based on the most 
    common conclusion. The agreement rate for a pool of effect sizes is the proportion 
    of times that the conclusions for individual studies are identical to the 
    conclusion (re: the same three categories) of the final, all-studies-combined MA. 
    More detailed descriptions of the analytic methods were provided by 
    O'Connor and Ermacora (2021).
    
    The output from this function can be entered into the PLOT_NO.PING.PONG function in
    order to obtain a graphical display of the findings across a sequence of studies.
}
\value{An object of class "NO.PING.PONG". The object is a list containing the
	following possible components:
\item{ES_MA}{The final effect size from the cumulative meta-analysis.}
 \item{ES_MA_lb}{The lower bound of the confidence interval for the final effect size 
 	from the cumulative meta-analysis}
 \item{ES_MA_ub}{The upper bound of the confidence interval for the final effect size 
 	from the cumulative meta-analysis}
 \item{Q}{The Q heterogeneity statistic.}
 \item{p_Q}{The statistical significance value for the Q heterogeneity statistic.}
 \item{tau2}{tau2 (or tau-squared) is the variation in effect sizes (between-study 
 	variance) in a random-effects meta-analysis. It is the variance in the true effect sizes.}
 \item{tau2LB}{The lower bound of the confidence interval for tau2.}
 \item{tau2UB}{The upper bound of the confidence interval for tau2.}
 \item{tau}{tau is the square root of tau-squared. tau is the standard deviation 
 	of the true effect sizes.}
 \item{tauLB}{The lower bound of the confidence interval for tau.}
 \item{tauUB}{The upper bound of the confidence interval for tau.}
 \item{isq}{isq estimates (in percent) how much of the total variability in the 
 	effect size estimates (which is composed of heterogeneity plus sampling 
 	variability) can be attributed to heterogeneity among the true effects.}
 \item{isqLB}{The lower bound of the confidence interval for isq.}
 \item{isqUB}{The upper bound of the confidence interval for isq.}
 \item{hsq}{hsq estimates the ratio of the total amount of variability in the 
 	effect size estimates to the amount of sampling variability.}
 \item{hsqLB}{The lower bound of the confidence interval for hsq.}
 \item{hsqUB}{The upper bound of the confidence interval for hsq.}
 \item{results_NHST}{The results for the NHST analyses.}
 \item{consistNHST}{The consistency rate for the NHST analyses.}
 \item{agreeNHST}{The agreement rate for the NHST analyses.}
 \item{results_CUM }{The results for the cumulative meta-analyses.}
 \item{consistCUM}{The consistency rate for the cumulative meta-analyses.}
 \item{agreeCUM}{The agreement rate for the cumulative meta-analyses.}
 \item{results_BA_SR }{The results for the Schmidt-Raju Bayesian analyses.}
 \item{consistBA_SR}{The consistency rate for the Schmidt-Raju Bayesian analyses.}
 \item{agreeBA_SR}{The agreement rate for the Schmidt-Raju Bayesian analyses.}
 \item{results_BA_GEN}{The results for the generated data Bayesian analyses.}
 \item{consistBA_GEN}{The consistency rate for the generated data Bayesian analyses.}
 \item{agreeBA_GEN}{The agreement rate for the generated data Bayesian analyses.}
 \item{results_BA_RAW}{The results for the raw data Bayesian analyses.}
 \item{consistBA_RAW}{The consistency rate for the raw data Bayesian analyses.}
 \item{agreeBA_RAW}{The agreement rate for the raw data Bayesian analyses.}
}
\references{
       {O'Connor, B. P., & Ermacora, D. (2021). Illustrations of why 
       	previous findings should be taken into account when evaluating new datasets. 
       	\emph{Canadian Journal of Behavioral Science}. 
       	Advance online publication. https://doi.org/10.1037/cbs0000259}       	
\cr\cr {O'Connor, B. P., & Khattar, N. (2021, in press). Controversies regarding 
	    null hypothesis testing. In W. O'Donohue, A. Masuda, & S. O. Lilienfeld (Eds.). 
	    \emph{Questionable research practices: Designing, conducting, and reporting sound 
	    research in clinical psychology.} Cham, Switzerland: Springer.}    	
}
\author{Brian P. O'Connor}
\examples{
\donttest{ 
# the IAT Discrimination data
head(data_IAT_Discrimination)

NO.PING.PONG(data_IAT_Discrimination, ES_type_IN='r', ES_type_OUT='r', 
             rma_method='REML',
             Bayes_type = c('Schmidt_Raju'), 
             prior_type='MA', CI = 95,
             ES = 'R', N = 'N', vi = NULL,
             nitt=13000, burnin=3000, thin=10, verbose=TRUE) 
}             

# data from SchmidtRaju (2007, p. 303)
data_Schmidt_Raju <- '
1    60   .44 
2    75   .20 
3    85   .60 
4   110   .32 
5    50   .41 
6    90   .25 
7   100   .12 
8    65   .35 
9    80   .35 
10   65   .19 '
data_Schmidt_Raju <- data.frame(read.table(text=data_Schmidt_Raju, fill=TRUE))
colnames(data_Schmidt_Raju) <- c('Study','N','r')
data_Schmidt_Raju <- data_Schmidt_Raju[,2:3]  

NO.PING.PONG(data_Schmidt_Raju, ES_type_IN='r', ES_type_OUT='r', 
             rma_method='REML',
             Bayes_type = c('generated', 'Schmidt_Raju'), 
             prior_type='MA', CI = 95,
             ES = 'r', N = 'N', vi = NULL,
             nitt=13000, burnin=3000, thin=10, verbose=TRUE)      
                   
}
