\name{PLOT_NO.PING.PONG}
\alias{PLOT_NO.PING.PONG}
\title{Plots of output from the NO.PING.PONG function}
\description{Provides plots of the output from the NO.PING.PONG function for a 
	sequence of studies}
\usage{PLOT_NO.PING.PONG(nopingpongOutput,
	                  plot_this = c('NHST','MA','Bayes_SR'),
	                  plot_save = FALSE, plot_save_type = 'png', 
	                  plot_title=NULL, Xrange=NULL)}
\arguments{
  \item{nopingpongOutput}{
  \code{}Output from the NO.PING.PONG function}

  \item{plot_this}{
  \code{}The kind of output data to be plotted. The options are one or any
   combination of 'NHST', 'MA', 'Bayes_gen', 'Bayes_raw', and 'Bayes_SR'. 
   'MA' is for meta-analysis results;
   'Bayes_gen' is for results based on Bayesian generated data analyses;
   'Bayes_raw' is for results based on Bayesian raw data analyses;
   'Bayes_SR' is for results based on the Schmidt-Raju (2007) Bayesian method.
   The default is
   plot_this = c('NHST', 'MA', 'Bayes_SR').}
  
  \item{plot_save}{
  \code{}Should a plot be saved to disk? TRUE or FALSE (the default).}

  \item{plot_save_type}{
  \code{}The output format if plot_save = TRUE. The options are 'bitmap', 'tiff', 
         'png' (the default), 'jpeg', and 'bmp'.}

  \item{plot_title}{
  \code{}optional. A title for the plot that will appear in the saved file name.}

  \item{Xrange}{
  \code{}optional. A range for the x axis in the plots.}
}
\details{This function provides plots of the output from the NO.PING.PONG 
	function for a sequence of studies, with options for specifying
	the kind of results to be plotted (via the plot_this argument),
	whether to save the plot to disc, the file type of the saved
	plot (via the plot_save_type argument), the plot title,
	and the x axis range for the plot.}
\value{A plot is produced, but there are no returned values.}
\references{
       {O'Connor, B. P., & Ermacora, D. (2021). Illustrations of why 
       	previous findings should be taken into account when evaluating new datasets. 
       	\emph{Canadian Journal of Behavioral Science}. 
       	Advance online publication. https://doi.org/10.1037/cbs0000259}       	
\cr\cr {O'Connor, B. P., & Khattar, N. (2021, in press). Controversies regarding 
	    null hypothesis testing. In W. O'Donohue, A. Masuda, & S. O. Lilienfeld (Eds.). 
	    \emph{Questionable research practices: Designing, conducting, and reporting sound 
	    research in clinical psychology.} Cham, Switzerland: Springer.}    	
}
\author{Brian P. O'Connor}
\examples{
# data from SchmidtRaju (2007, p. 303)
data_Schmidt_Raju <- '
1    60   .44 
2    75   .20 
3    85   .60 
4   110   .32 
5    50   .41 
6    90   .25 
7   100   .12 
8    65   .35 
9    80   .35 
10   65   .19 '
data_Schmidt_Raju <- data.frame(read.table(text=data_Schmidt_Raju, fill=TRUE))
colnames(data_Schmidt_Raju) <- c('Study','N','r')
data_Schmidt_Raju <- data_Schmidt_Raju[,2:3]  

NO.PING.PONG(data_Schmidt_Raju, ES_type_IN='r', ES_type_OUT='r', 
             rma_method='REML',
             Bayes_type = c('generated', 'Schmidt_Raju'), 
             prior_type='MA', CI = 95,
             ES = 'r', N = 'N', vi = NULL,
             nitt=13000, burnin=3000, thin=10, verbose=TRUE)      

\donttest{ 
# Cannabis Psychosis data
nppOutput <- NO.PING.PONG(data_Cannabis_Psychosis, ES_type_IN='d', rma_method='REML', CI = 95,
                          Bayes_type = c('Schmidt_Raju', 'generated'), prior_type='MA',
                          ES = 'Std_diff_in_mean', N = 'N', vi = 'Variance',
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 

PLOT_NO.PING.PONG(nppOutput, plot_this = c('NHST','MA'))
PLOT_NO.PING.PONG(nppOutput, plot_this = c('NHST','MA','Bayes_gen'))
PLOT_NO.PING.PONG(nppOutput, plot_this = c('NHST','MA','Bayes_gen','Bayes_raw'))


# PopulationR.20 data (has raw data)
nppOutput = NO.PING.PONG(data_PopulationR.20, Bayes_type = c('generated', 'raw'), prior_type='MA', 
                         nitt=13000, burnin=3000, thin=10, verbose=TRUE)
                         
PLOT_NO.PING.PONG(nppOutput, plot_this = c('NHST','MA','Bayes_gen','Bayes_raw'))
}
}
