\name{data_Omega3_Depression}
\alias{data_Omega3_Depression}
\docType{data}
\title{data_Omega3_Depression}
\description{Data from a meta-analysis on 12 effect sizes for the role of omega-3 
	fatty acids in the treatment of major depression (Grosso et al., 2014, 
	Figure 2, top portion, p. 9).

}
\usage{data(data_Omega3_Depression)}

\source{
   Grosso, G., Pajak, A., Marventano, S., Castellano, S., Galvano, F., Bucolo, C., & 
   Caraci, F. (2014). Role of omega-3 fatty acids in the treatment of depressive 
   disorders: A comprehensive meta-analysis of randomized clinical trials. 
   \emph{PLoS ONE, 9(5),} e96905.
}

\examples{
head(data_Omega3_Depression)

\donttest{ 
nppOutput <- NO.PING.PONG(data_Omega3_Depression, ES_type_IN='d', ES_type_OUT='g', 
                          rma_method='REML',
                          Bayes_type = c('generated', 'Schmidt_Raju'), 
                          prior_type='MA', CI = 95,
                          ES = NULL, N = NULL, vi = NULL,
                          grp1_mn = 'Cmn', grp1_sd = 'Csd', grp1_n = 'CN', 
                          grp2_mn = 'Emn', grp2_sd = 'Esd', grp2_n = 'EN', 
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 

PLOT_NO.PING.PONG(nppOutput)}
}
