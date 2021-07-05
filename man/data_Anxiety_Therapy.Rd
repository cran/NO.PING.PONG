\name{data_Anxiety_Therapy}
\alias{data_Anxiety_Therapy}
\docType{data}
\title{data_Anxiety_Therapy}
\description{Data from a meta-analysis on 81 effect sizes (Hedges g values) for the  
	pre-post difference between transdiagnostic cognitive behavior therapy and anxiety 
	(Pearl & Norton, 2017, Figure 2, p. 14). 
}
\usage{data(data_Anxiety_Therapy)}

\source{
   Pearl, S. B., & Norton, P. J. (2017). Transdiagnostic versus diagnosis 
   specific cognitive behavioural therapies for anxiety: A meta-analysis. 
   \emph{Journal of Anxiety Disorders, 46,} 11-24.
}

\examples{
head(data_Anxiety_Therapy)

\donttest{ 
nppOutput <- NO.PING.PONG(data_Anxiety_Therapy, ES_type_IN='g', ES_type_OUT='g', 
                          rma_method='REML',
                          Bayes_type = c('Schmidt_Raju','generated'), 
                          prior_type='MA', CI = 95,
                          ES = 'Hedges.g', N = 'N', vi = 'Var',
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 

PLOT_NO.PING.PONG(nppOutput)
}
}
