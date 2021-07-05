\name{data_SelfEsteem_Depression}
\alias{data_SelfEsteem_Depression}
\docType{data}
\title{data_SelfEsteem_Depression}
\description{Sowislo and Orth (2013, Table 2) reported findings from a meta-analysis of 
	77 longitudinal studies that provided effect sizes for self-esteem predicting future 
	depressive symptoms (median N = 224). The final effect size, in correlation 
	coefficient metric, was .57.
}
\usage{data(data_SelfEsteem_Depression)}

\source{
   Sowislo, J. F., & Orth, U. (2013). Does low self-esteem predict depression 
   and anxiety? A meta-analysis of longitudinal studies. 
   \emph{Psychological Bulletin, 139(1),} 213240. https://doi.org/10.1037/a0028931
}

\examples{
head(data_SelfEsteem_Depression)

\donttest{ 
nppOutput <- NO.PING.PONG(data_SelfEsteem_Depression, ES_type_IN='r', ES_type_OUT='r', 
                          rma_method='REML',
                          Bayes_type = c('generated', 'Schmidt_Raju'), 
                          prior_type='MA', CI = 95,
                          ES = 'rSED', N = 'N', vi = NULL,
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 

PLOT_NO.PING.PONG(nppOutput)
}
}


             
             
             
             
             
             