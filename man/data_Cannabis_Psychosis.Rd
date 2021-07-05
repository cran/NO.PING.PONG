\name{data_Cannabis_Psychosis}
\alias{data_Cannabis_Psychosis}
\docType{data}
\title{data_Cannabis_Psychosis}
\description{Data from a meta-analysis on 38 effect sizes (standardized mean 
	differences) from comparisons of the age at onset of psychosis in 
	cannabis-using groups with the age at onset of psychosis in 
	noncannabis-using groups (Large et al., 2001). The data are from the 
	"eFigure" in the Supplementary Material.
}
\usage{data(data_Cannabis_Psychosis)}
\source{
   Large, M., Sharma, S., Compton, M. T., Slade, T., Nielssen, O. (2011). Cannabis 
   use and earlier onset of psychosis: a systematic meta-analysis. 
   \emph{Archives of General Psychiatry, 68(6),} 555-561.
}
\examples{
head(data_Cannabis_Psychosis)

\donttest{ 
nppOutput <- NO.PING.PONG(data_Cannabis_Psychosis, ES_type_IN='d', ES_type_OUT='d', 
                          rma_method='REML',
                          Bayes_type = c('Schmidt_Raju', 'generated'), 
                          prior_type='MA', CI = 95,
                          ES = 'Std_diff_in_mean', N = 'N', vi = 'Variance',
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 

PLOT_NO.PING.PONG(nppOutput)
}
}
