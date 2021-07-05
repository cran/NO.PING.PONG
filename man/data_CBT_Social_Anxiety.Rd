\name{data_CBT_Social_Anxiety}
\alias{data_CBT_Social_Anxiety}
\docType{data}
\title{data_CBT_Social_Anxiety}
\description{Kampmanna, Emmelkamp, and Morina (2016, Figure 4) reported findings 
	from a meta-analysis of 24 studies on internet delivered cognitive behavior 
	therapy (vs. control conditions) for social anxiety disorder (median N = 65). 
	The final effect size, in correlation coefficient metric, was .32.
}
\usage{data(data_CBT_Social_Anxiety)}

\source{
   Kampmann, I. L., Emmelkamp, P. M. G., & Morina, N. (2016). Meta-analysis of 
   technology-assisted interventions for social anxiety disorder. 
   \emph{Journal of Anxiety Disorders, 42,} 7184. 
   \cr {https://doi.org/10.1016/j.janxdis.2016.06.007}
}

\examples{
head(data_CBT_Social_Anxiety)

\donttest{ 
data_CBT_Social_Anxiety_r <- CONVERT_ES(
     data_CBT_Social_Anxiety$Hedgesg, ES_var = data_CBT_Social_Anxiety$Std.Error, 
     ES_type_IN='g', ES_type_OUT='r', 
     grp1_N = data_CBT_Social_Anxiety$N.exp, grp2_N = data_CBT_Social_Anxiety$N.ctrl,
     gvar_type_OUT = 'd', verbose = TRUE)        

nppOutput <-NO.PING.PONG(data_CBT_Social_Anxiety_r, ES_type_IN='r', ES_type_OUT='r', 
                         rma_method='REML',
                         Bayes_type = c('Schmidt_Raju'), 
                         prior_type='MA', CI = 95,
                         ES = 'r', N = 'totalN', vi = 'Vr',
                         gvar_type_OUT = 'd',
                         nitt=13000, burnin=3000, thin=10, verbose=TRUE) 

PLOT_NO.PING.PONG(nppOutput, plot_this = c('NHST','MA','Bayes_SR'), Xrange=c(-1,1))

#PLOT_NO.PING.PONG(nppOutput, plot_this = c('NHST','MA','Bayes_SR'),
#                  plot_save=TRUE, plot_save_type='png', 
#                  plot_title="ICBT for Social Anxiety", 
#                  Xrange=c(-1,1))
}
}




