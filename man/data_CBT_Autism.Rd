\name{data_CBT_Autism}
\alias{data_CBT_Autism}
\docType{data}
\title{data_CBT_Autism}
\description{Weston, Hodgekins, and Langdon (2016) reported findings from 
	a meta-analysis of 17 studies on the effectiveness of CBT on affective 
	symptoms for people with autistic spectrum disorders (median N = 36). 
	The final effect size, in correlation coefficient metric, was .11. 
}
\usage{data(data_CBT_Autism)}

\source{
   Weston, L., Hodgekins, J., & Langdon, P. E. (2016). Effectiveness of 
   cognitive behavioural therapy with people who have autistic spectrum 
   disorders: A systematic review and meta-analysis. 
   \emph{Clinical Psychology Review, 49,} 4154. 
   https://doi.org/10.1016/j.cpr.2016.08.001
}

\examples{
head(data_CBT_Autism)

\donttest{ 
nppOutput <- NO.PING.PONG(data_CBT_Autism[2:nrow(data_CBT_Autism),], 
                          ES_type_IN=NULL, ES_type_OUT='r', 
                          rma_method='REML', CI = 95,
                          Bayes_type = c('generated','Schmidt_Raju'), prior_type='MA',
                          ES = NULL, N = NULL, vi = NULL,
                          grp1_mn = 'Con_Mean', grp1_sd = 'Con_SD', grp1_n = 'Con_N', 
                          grp2_mn = 'CBT_Mean', grp2_sd = 'CBT_SD', grp2_n = 'CBT_N', 
                          nitt=50000, burnin=3000, thin=10, verbose=TRUE) 

PLOT_NO.PING.PONG(nppOutput, plot_this = c('NHST','MA','Bayes_SR'), Xrange=c(-1,1))

# PLOT_NO.PING.PONG(nppOutput, plot_this = c('NHST','MA','Bayes_SR'),
                  # plot_save=TRUE, plot_save_type='png', plot_title="CBT for Autism", 
                  # Xrange=c(-1,1))
}
}


