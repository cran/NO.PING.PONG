\name{data_Math_Performance}
\alias{data_Math_Performance}
\docType{data}
\title{data_Math_Performance}
\description{Data from a meta-analysis on 35 effect sizes (r values) for the association  
	between individual differences in non-symbolic number acuity and math performance 
	(Chen & Li, 2014, Table 1, p. 167).
}
\usage{data(data_Math_Performance)}

\source{
   Chen, Q., & Li, J. (2014). Association between individual differences in 
   non-symbolic number acuity and math performance: A meta-analysis. 
   \emph{Acta Psychologica, 148,} 163-172.
}

\examples{
head(data_Math_Performance)

\donttest{ 
nppOutput <- NO.PING.PONG(data_Math_Performance, ES_type_IN='r', ES_type_OUT='r', 
                          rma_method='REML',
                          Bayes_type = c('generated', 'Schmidt_Raju'), 
                          prior_type='MA', CI = 95,
                          ES = 'r', N = 'N', vi = NULL,
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 

PLOT_NO.PING.PONG(nppOutput)
}
}



             
             
             
             
             
             