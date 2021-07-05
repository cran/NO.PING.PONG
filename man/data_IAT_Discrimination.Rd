\name{data_IAT_Discrimination}
\alias{data_IAT_Discrimination}
\docType{data}
\title{data_IAT_Discrimination}
\description{Data from a meta-analysis on 298 effect sizes (r values) for the association  
	between scores on the Implicit Association Test (IAT) and ethnic and racial discrimination 
	(Oswald, Mitchell, Blanton, Jaccard, & Tetlock (2013).
}
\usage{data(data_IAT_Discrimination)}

\source{
   Oswald, F. L., Mitchell, G., Blanton, H., Jaccard, J., & Tetlock, P. E. (2013). Predicting 
   ethnic and racial discrimination: A meta-analysis of IAT criterion studies. 
   \emph{Journal of Personality and Social Psychology, 105,} 171-192.
}

\examples{
head(data_IAT_Discrimination)

\donttest{ 
nppOutput <- NO.PING.PONG(data_IAT_Discrimination, ES_type_IN='r', ES_type_OUT='r', 
                          rma_method='REML',
                          Bayes_type = c('generated', 'Schmidt_Raju'), 
                          prior_type='MA', CI = 95,
                          ES = 'R', N = 'N', vi = NULL,
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 

PLOT_NO.PING.PONG(nppOutput)
}
}
