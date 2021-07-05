\name{data_Ego_Depletion}
\alias{data_Ego_Depletion}
\docType{data}
\title{data_Ego_Depletion}
\description{Data from the pre-registered, replication studies for 
	the ego-depletion effect (Hagger et al., 2016). These data were from 23 
	laboratories (N = 2141).
}
\usage{data(data_Ego_Depletion)}

\source{
   Hagger, M. S., Chatzisarantis, N. L., Alberts, H. A., Anggono, C. O., Batailler, C. B., 
   Birt, A. R., et al. (2016). A multilab preregistered replication of the ego-depletion 
   effect. \emph{Perspectives on Psychological Science, 11(4),} 546-573.
}

\examples{
head(data_Ego_Depletion)

\donttest{ 
nppOutput <- NO.PING.PONG(data_Ego_Depletion, ES_type_IN=NULL, ES_type_OUT='g',
                          rma_method='REML',
                          Bayes_type = c('Schmidt_Raju', 'generated', 'raw'), 
                          prior_type='MA', CI = 95,
                          ES = NULL, N = NULL, vi = NULL,
                          grp1_mn = 'H.Mean', grp1_sd = 'H.SD', grp1_n = 'H.N', 
                          grp2_mn = 'E.Mean', grp2_sd = 'E.SD', grp2_n = 'E.N', 
                          gvar_type_OUT = 'd',
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 
                                                 
PLOT_NO.PING.PONG(nppOutput)
}
}

