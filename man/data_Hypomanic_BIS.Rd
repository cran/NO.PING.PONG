\name{data_Hypomanic_BIS}
\alias{data_Hypomanic_BIS}
\docType{data}
\title{data_Hypomanic_BIS}
\description{Katz, Naftalovich, Matanky, and Yovel (2021) reported findings 
	from a meta-analysis of 19 studies on hypomanic personality tendencies 
	and behavioral inhibition system sensitivity (median N = 230). The final 
	effect size, in correlation coefficient metric, was -.04.
}
\usage{data(data_Hypomanic_BIS)}

\source{
   Katz, B. A., Naftalovich, H., Matanky, K., & Yovel, I. (2021). The 
   dual-system theory of bipolar spectrum disorders: A meta-analysis. 
   \emph{Clinical Psychology Review, 83,} Article 101945. 
   \cr {https://doi.org/10.1016/j.cpr.2020.101945}
}

\examples{
head(data_Hypomanic_BIS)

\donttest{ 
nppOutput <- NO.PING.PONG(data_Hypomanic_BIS, ES_type_IN=NULL, ES_type_OUT='r', 
                          rma_method='REML', CI = 95,
                          Bayes_type = c('Schmidt_Raju'), prior_type='MA',
                          ES = NULL, N = NULL, vi = NULL,
                          grp1_mn = 'BIS_B_MClin', grp1_sd = 'BIS_B_SDClin', grp1_n = 'Nclin', 
                          grp2_mn = 'BIS_B_MHC',   grp2_sd = 'BIS_B_SDHC',   grp2_n = 'Nhealthy', 
                          nitt=50000, burnin=3000, thin=10, verbose=TRUE) 
                         
PLOT_NO.PING.PONG(nppOutput, plot_this = c('NHST','MA','Bayes_SR'), Xrange=c(-1,1))

#PLOT_NO.PING.PONG(nppOutput, plot_this = c('NHST','MA','Bayes_SR'),
#                  plot_save=TRUE, plot_save_type='png', plot_title="Bipolar & BIS", Xrange=c(-1,1))
}
}


