\name{data_PopulationR.20}
\alias{data_PopulationR.20}
\docType{data}
\title{data_PopulationR.20}
\description{A list with raw data for two variables from 50 random samples,  
	N = 100 each, that were randomly drawn from a population of 100,000 cases in 
	which the correlation between the two variables was .20.}
\usage{data(data_PopulationR.20)}

\examples{
head(data_PopulationR.20)

\donttest{ 
nppOutput <- NO.PING.PONG(data_PopulationR.20, ES_type_IN='r', ES_type_OUT='r',
                          rma_method='REML',
                          Bayes_type = c('generated', 'Schmidt_Raju'), 
                          prior_type='MA', CI = 95,
                          ES = 'r', N = 'N', vi = NULL,
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 

PLOT_NO.PING.PONG(nppOutput)
}
}

