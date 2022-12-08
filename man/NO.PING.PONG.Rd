\name{NO.PING.PONG}
\alias{NO.PING.PONG}
\title{Incorporating previous findings when evaluating new data}
\description{A function for revealing what happens when previous findings 
	are taken into account when analyzing new data.}
\usage{NO.PING.PONG(donnes, 
             ES_type_IN=NULL, ES_type_OUT='r',
             rawdata_type = 'for_correl',
             rma_method='REML',
             Bayes_type = c('Schmidt_Raju'), 
             prior_type='META', CI = 95,
             ES = NULL, N = NULL, vi = NULL,
             grp1_mn = NULL, grp1_sd = NULL, grp1_n = NULL, 
             grp2_mn = NULL, grp2_sd = NULL, grp2_n = NULL,
             gvar_type_OUT = 'd',
             paired_samples_ES_type = NULL,
             funnel_plot=FALSE, funnel_plot_type='png', funnel_plot_title=NULL, 
             nitt=53000, burnin=3000, thin=10, 
             verbose=TRUE)}
\arguments{
  \item{donnes}{
  \code{}Either (1) a list with raw data for two numeric variables in each of the multiple list elements
  (wherein each list element consists of rows of raw data from an individual study), or,  
  (2) a dataframe or matrix wherein each row has summary data for a single study
  (e.g., has the group means, SDs, & Ns, or d or g or r values, or effect size and variance
  estimates).}

  \item{ES_type_IN}{
  \code{}(optional) The type of effect sizes, if donnes is a matrix of study-level data rather than raw data. 
  The options are 'r', 'd', or 'g'.}

  \item{ES_type_OUT}{
  \code{}(optional) The type of effect sizes for the output. The options are 'r' (the default), 'd', and 'g'.}

  \item{rawdata_type}{
  \code{}(optional) The type of raw data, if donnes is a list with raw datasets. The options are 'for_correl'
   (for when correlations should be computed), 'indep_groups' (for when the raw data are based on 
   independent groups, as in for an independent groups t-test), and 'paired_samples' (for when 
   the raw data are paired samples, as in for a paired samples t-test).}

  \item{rma_method}{
  \code{}(optional) The method option for the rma function from the metafor package: 'A character 
  string specifying whether a fixed- or a random/mixed-effects model should be fitted. 
  A fixed-effects model (with or without moderators) is fitted when using method="FE". 
  Random/mixed-effects models are fitted by setting method equal to one of the following: 
  "DL", "HE", "SJ", "ML", "REML", "EB", "HS", or "GENQ". The default is "REML".'}

  \item{Bayes_type}{
  \code{}(optional) The kind(s) of Bayesian analyses to be conducted. The options are one or all
  of 'Schmidt_Raju', 'generated' and/or 'raw'. The default is 'Schmidt_Raju', which is the fastest,  
  i.e., Bayes_type = c('Schmidt_Raju'). See the Details.}

  \item{prior_type}{
  \code{}(optional) The type of prior data to be used in the updating analyses (both cumulative and
  Bayesian). The options are 'BAYES' or 'META'. The default is 'META'. See the Details.}

  \item{CI}{
  \code{}(optional) The confidence interval for the analyses (in whole numbers). The default is 95.}

  \item{ES}{
  \code{}(optional) The name of the column in donnes with the effect sizes.}

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

  \item{paired_samples_ES_type}{
  \code{}(optional) The kind of effect size for the analyses if rawdata_type = 'paired_samples'.
  Following the escalc function from the metafor package,
  the options are 'MC' (for raw mean change), 'SMCC' (for the standardized mean change using 
  change score standardization; Gibbons et al., 1993), 'SMCR' (for the standardized mean 
  change using raw score standardization; Becker, 1988), 'SMCRH' (the default, for the 
  standardized mean change using raw score standardization with heteroscedastic population 
  variances at the two measurement occasions, Bonett, 2008), or 'ROMC' ( for the log 
  transformed ratio of means; Lajeunesse, 2011).}

  \item{funnel_plot}{
  \code{}(optional) Should a funnel plot be produced and saved? TRUE or FALSE (default)}

  \item{funnel_plot_type}{
  \code{}(optional) The output format if funnel_plot_save = TRUE. The options are 'bitmap', 'tiff', 
         'png' (the default), 'jpeg', and 'bmp'.}

  \item{funnel_plot_title}{
  \code{}(optional) Text that will be used to name the funnel plot file and appear as the plot main title.}

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
	
	The function relies on the metafor package for the meta-analyses and on the
	MCMCglmm package for the Bayesian analyses.
	
	If the function is provided with only non-correlation coefficient effect sizes and corresponding 
	variance estimates, then the analyses will be conducted for meta-analyses and NHST. Bayesian 
	analyses are not possible. However, any kind of effect size coefficient can be used in this case.
	
	If the function is provided with correlation coefficient effect sizes and the corresponding 
	variance estimates, then the analyses can be conducted for NHST, meta-analyses, and Bayesian 
	analyses.
		
	If the function is provided with d or g effect sizes along with the corresponding total 
	sample sizes, but not with the corresponding variances for the effect sizes, then the
	d or g values will be converted to correlation coefficient equivalents and the
	subsequent analyses will be for r effect sizes. Bayesian analyses are possible in
	this case.

    The Bayesian analyses are always conducted using correlation coefficient effect sizes. When the
    input data are not raw data nor correlation coefficients, then the correlation coefficient
    equivalents are computed from the input data for the analyses. 
    
    When the input data consists of raw datasets and the focus is on group
    mean differences in scores (re: rawdata_type = 'indep_groups'), then it is assumed 
    that the grouping variable is in the first column of the input data ('donnes') and that there 
    are only two levels of the grouping variable.    

    When 'd' or 'g' are specified as the type of effect sizes for the output (using the 
    ES_type_OUT argument), then the output is for conversions of correlation coefficient
    effect sizes to 'd' or 'g' values using the conventional formulas for such conversions.
    
    The gvar_type_OUT argument provides a choice between d or g effect size variances whenever
    ES_type_OUT is set to g. The reason for this option is that
    authors of published meta-analyses sometimes report d variances when their analyses were
    conducted on g effect sizes. This is presumably not a wise practice, but it also does
    not make much difference in the computed values.
    	
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
	The options are 'BAYES' or 'META', and 'META' is the default.
	
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
 \item{biasStats}{Publication bias statistics.}
}
\references{
       {OConnor, B. P., & Ermacora, D. (2021). Unnecessary ping-pong: Illustrations of 
       	why previous findings should be taken into account when evaluating new datasets. 
       	\emph{Canadian Journal of Behavioural Science, 53(3),} 328-341. https://doi.org/10.1037/cbs0000259}       	
\cr\cr {O'Connor, B. P., & Khattar, N. (2022, in press). Controversies regarding 
	    null hypothesis testing. In W. O'Donohue, A. Masuda, & S. O. Lilienfeld (Eds.). 
	    \emph{Avoiding Questionable Research Practices in Applied Psychology} (pp. 147-174). 
	    Cham, Switzerland: Springer Nature Switzerland.}    	
}
\author{Brian P. O'Connor}
\examples{
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

NO.PING.PONG(data_Schmidt_Raju, ES_type_IN='r', ES_type_OUT='r', 
             rma_method='REML',
             Bayes_type = c('Schmidt_Raju', 'generated'), 
             prior_type='META', CI = 95,
             ES = 'r', N = 'N', vi = NULL,
             nitt=13000, burnin=3000, thin=10)      
                          	
\donttest{
	
# using only ES & vi (the effect size & the effect size variance for each study)
# Anxiety_Therapy
# input data = Hedges.g & vi for each study
# Kampmanna (2016). Meta-analysis of technology-assisted interventions for social anxiety disorder
NO.PING.PONG(data_NPP$Anxiety_Therapy, 
             ES = 'Hedges.g', N = NULL, vi = 'Var') 


# Alcohol_Intake
# input data = Hedges g & the Ns for the experimental & control groups
# Prestwich (2016). Does Changing Social Influence Engender Changes in Alcohol Intake  Table 2
NO.PING.PONG(data_NPP$Alcohol_Intake, ES_type_IN='g', 
             ES = 'g.Alcohol.Intake', grp1_n = 'Exp.n', grp2_n = 'Ctrl.n') 


# Anxiety_Therapy   
# input data = Hedges g & the N for each study
# Kampmanna (2016). Meta-analysis of technology-assisted interventions for social anxiety disorder
NO.PING.PONG(data_NPP$Anxiety_Therapy, ES_type_IN='g', ES_type_OUT='r',
             ES = 'Hedges.g', N = 'N',  rma_method='FE') 


# Cannabis_Psychosis   
# input data = Cohen's d & the N for each study
# Large (2001). Cannabis Use and Earlier Onset of Psychosis - A Systematic Meta-analysis
# Supplementary Online Content - The data are from the "eFigure".
NO.PING.PONG(data_NPP$Cannabis_Psychosis, ES_type_IN='d',
             ES = 'Std_diff_in_mean', N = 'N') 


# CBT_Autism
# input data = group means, SDs, & Ns
# Weston (2016). Effectiveness of cognitive behavioural therapy with people who have autistic 
# spectrum disorders  A systematic review and meta-analysis
NO.PING.PONG(data_NPP$CBT_Autism,   
             grp1_mn = 'Con_Mean', grp1_sd = 'Con_SD', grp1_n = 'Con_N', 
             grp2_mn = 'CBT_Mean', grp2_sd = 'CBT_SD', grp2_n = 'CBT_N') 


# CBT_Social_Anxiety    
# input data = Hedges g & the Ns for the experimental & control groups
# Kampmanna (2016). Meta-analysis of technology-assisted interventions for social anxiety  fig 4
NO.PING.PONG(data_NPP$CBT_Social_Anxiety, ES_type_IN='g', 
             ES = 'Hedgesg', grp1_n = 'N.exp', grp2_n = 'N.ctrl') 


# Ego_Depletion
# input data = group means, SDs, & Ns
# Hagger (2016). A multilab preregistered replication of the ego-depletion effect
NO.PING.PONG(data_NPP$Ego_Depletion,  
             grp1_mn = 'E.Mean', grp1_sd = 'E.SD', grp1_n = 'E.N', 
             grp2_mn = 'H.Mean', grp2_sd = 'H.SD', grp2_n = 'H.N') 


# Hypomanic_BIS
# input data = group means, SDs, & Ns
# Katz (2021). The Dual-System Theory of Bipolar Spectrum Disorders   r (BIS) = -.04
NO.PING.PONG(data_NPP$Hypomanic_BIS,  
             grp1_mn = 'BIS_B_MClin', grp1_sd = 'BIS_B_SDClin', grp1_n = 'Nclin', 
             grp2_mn = 'BIS_B_MHC',   grp2_sd = 'BIS_B_SDHC',   grp2_n = 'Nhealthy') 


# IAT_Discrimination
# input data = r & N for each study
# Oswald (2013). Predicting Ethnic and Racial Discrimination - A Meta-Analysis of IAT Studies
NO.PING.PONG(data_NPP$IAT_Discrimination, ES_type_IN='r', 
             ES = 'R', N = 'N') 


# Many_Labs$Anchoring
# input data = group means, SDs, & Ns
# Klein et al. (2014). Investigating variation in replicability: A many labs replication project
NO.PING.PONG(data_NPP$Many_Labs$Anchoring,  
             grp1_mn = 'grp1_mn', grp1_sd = 'grp1_sd', grp1_n = 'grp1_N', 
             grp2_mn = 'grp2_mn', grp2_sd = 'grp2_sd', grp2_n = 'grp2_N') 


# Many_Labs$Gamblers_Fallacy
# input data = group means, SDs, & Ns
# Klein et al. (2014). Investigating variation in replicability: A many labs replication project
NO.PING.PONG(data_NPP$Many_Labs$Gamblers_Fallacy,  
             grp1_mn = 'grp1_mn', grp1_sd = 'grp1_sd', grp1_n = 'grp1_N', 
             grp2_mn = 'grp2_mn', grp2_sd = 'grp2_sd', grp2_n = 'grp2_N') 


# Many_Labs$Math_Attitudes
# input data = group means, SDs, & Ns
# Klein et al. (2014). Investigating variation in replicability: A many labs replication project
NO.PING.PONG(data_NPP$Many_Labs$Math_Attitudes,  
             grp1_mn = 'grp1_mn', grp1_sd = 'grp1_sd', grp1_n = 'grp1_N', 
             grp2_mn = 'grp2_mn', grp2_sd = 'grp2_sd', grp2_n = 'grp2_N') 


# Many_Labs$Sunk_Costs
# input data = group means, SDs, & Ns
# Klein et al. (2014). Investigating variation in replicability: A many labs replication project
NO.PING.PONG(data_NPP$Many_Labs$Sunk_Costs,  
             grp1_mn = 'grp1_mn', grp1_sd = 'grp1_sd', grp1_n = 'grp1_N', 
             grp2_mn = 'grp2_mn', grp2_sd = 'grp2_sd', grp2_n = 'grp2_N') 


# Many_Labs$Quote_Attribution
# input data = group means, SDs, & Ns
# Klein et al. (2014). Investigating variation in replicability: A many labs replication project
NO.PING.PONG(data_NPP$Many_Labs$Quote_Attribution,  
             grp1_mn = 'grp1_mn', grp1_sd = 'grp1_sd', grp1_n = 'grp1_N', 
             grp2_mn = 'grp2_mn', grp2_sd = 'grp2_sd', grp2_n = 'grp2_N') 


# Many_Labs$Flag_Priming
# input data = group means, SDs, & Ns
# Klein et al. (2014). Investigating variation in replicability: A many labs replication project
NO.PING.PONG(data_NPP$Many_Labs$Flag_Priming,  
             grp1_mn = 'grp1_mn', grp1_sd = 'grp1_sd', grp1_n = 'grp1_N', 
             grp2_mn = 'grp2_mn', grp2_sd = 'grp2_sd', grp2_n = 'grp2_N') 


# Math_Performance
# input data = r & N for each study
# Chen (2014). Association between individual differences in non-symbolic number acuity and 
# math performance - A meta-analysis Table 1
NO.PING.PONG(data_NPP$Math_Performance, ES_type_IN='r', 
             ES = 'r', N = 'N') 


# Omega3_Depression
# Grosso (2014).  Role of omega-3 fatty acids in the treatment of depressive disorders
# Findings are for Omega 3 predicting DSM major depression improvement - Fig 2, top portion, p. 9
# input data = group means, SDs, & Ns
NO.PING.PONG(data_NPP$Omega3_Depression,  
             grp1_mn = 'Cmn', grp1_sd = 'Csd', grp1_n = 'CN', 
             grp2_mn = 'Emn', grp2_sd = 'Esd', grp2_n = 'EN') 


# paired samples             
# input data = raw data
NO.PING.PONG(donnes=data_NPP$Paired_Samples, 
             rawdata_type = 'paired_samples',
             paired_samples_ES_type = 'SMCRH') 


# PopulationR.02 - generated data for 2 variables with a population r = .02
# input data = raw data
# O'Connor & Ermacora (2012). Unnecessary ping-pong
NO.PING.PONG(data_NPP$PopulationR.02) 


# PopulationR.10 - generated data for 2 variables with a population r = .10
# input data = raw data
# O'Connor & Ermacora (2012). Unnecessary ping-pong
NO.PING.PONG(data_NPP$PopulationR.10) 


# PopulationR.11 - generated data for 2 variables with a population r = .11
# input data = raw data
# O'Connor & Khattar (2022). Controversies regarding null hypothesis testing
NO.PING.PONG(data_NPP$PopulationR.11) 


# PopulationR.20 - generated data for 2 variables with a population r = .20
# input data = raw data
# O'Connor & Ermacora (2012). Unnecessary ping-pong
NO.PING.PONG(data_NPP$PopulationR.20) 


head(data_NPP$PopulationR.32)
# PopulationR.32 - generated data for 2 variables with a population r = .32
# input data = raw data
# O'Connor & Khattar (2022). Controversies regarding null hypothesis testing
NO.PING.PONG(data_NPP$PopulationR.32) 


head(data_NPP$PopulationR.57)
# PopulationR.57 - generated data for 2 variables with a population r = .57
# input data = raw data
# O'Connor & Khattar (2022). Controversies regarding null hypothesis testing
NO.PING.PONG(data_NPP$PopulationR.57) 


head(data_NPP$PopulationR.077)
# PopulationR.077 - generated data for 2 variables with a population r = .077
# input data = raw data
# O'Connor & Khattar (2022). Controversies regarding null hypothesis testing
NO.PING.PONG(data_NPP$PopulationR.077) 


head(data_NPP$PopulationRneg.04)
# PopulationRneg.04 - generated data for 2 variables with a population r = -04
# input data = raw data
# O'Connor & Khattar (2022). Controversies regarding null hypothesis testing
NO.PING.PONG(data_NPP$PopulationRneg.04) 


# SelfEsteem_Depression
# input data = r & N for each study
# 2013 Sowislo - Does Low Self-Esteem Predict Depression and Anxiety  r = .57
NO.PING.PONG(data_NPP$SelfEsteem_Depression, ES_type_IN='r', 
             ES = 'rSED', N = 'N') 
            
}             

}
