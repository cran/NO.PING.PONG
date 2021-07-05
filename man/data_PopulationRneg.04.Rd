\name{data_PopulationRneg.04}
\alias{data_PopulationRneg.04}
\docType{data}
\title{data_PopulationRneg.04}
\description{A list with raw data for two variables from 50 random samples,  
	N = 230 each, that were randomly drawn from a population of 100,000 cases in 
	which the correlation between the two variables was -.04.}
\usage{data(data_PopulationRneg.04)}

\examples{
head(data_PopulationRneg.04)

\donttest{ 
nppOutput <- NO.PING.PONG(data_PopulationRneg.04, ES_type_IN='r', ES_type_OUT='r',
                          rma_method='REML',
                          Bayes_type = c('generated', 'Schmidt_Raju'), 
                          prior_type='MA', CI = 95,
                          ES = 'r', N = 'N', vi = NULL,
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 

PLOT_NO.PING.PONG(nppOutput, plot_this = c('NHST','MA','Bayes_SR'), Xrange=c(-1,1))

#PLOT_NO.PING.PONG(nppOutput, plot_this = c('NHST','MA','Bayes_SR'),
#                  plot_save=TRUE, plot_save_type='png', plot_title='Pop r = -.04', Xrange=c(-1,1))

}
}

