\name{data_Alcohol_Intake}
\alias{data_Alcohol_Intake}
\docType{data}
\title{data_Alcohol_Intake}

\description{Data from a meta-analysis on 36 effect sizes (Hedges g values) for the effect of  
	manipulations of social influences on perceived norms regarding alcohol intake and on
	reductions in alcohol intake (Prestwich, Kellar, Conner, 
	Lawton, Gardner,& Turgut, 2016, Table 2, p. 850).
}
\usage{data(data_Alcohol_Intake)}

\source{
   Prestwich, A., Kellar, I., Conner, M., Lawton, R., Gardner, P., & Turgut, L. (2016). 
   Does changing social influence engender changes in alcohol intake? A meta-analysis. 
   \emph{Journal of Consulting & Clinical Psychology, 84,} 845-860.
}

\examples{
head(data_Alcohol_Intake)

\donttest{ 
# compute totalN
data_Alcohol_Intake$totalN <- rowSums(data_Alcohol_Intake[,c('Exp.n','Ctrl.n')])


# changes in alcohol intake
nppOutput <- NO.PING.PONG(data_Alcohol_Intake, ES_type_IN='g', ES_type_OUT='g', 
                          rma_method='REML', 
                          Bayes_type = c('generated', 'Schmidt_Raju'), 
                          prior_type='MA', CI = 95,
                          ES = 'g.Alcohol.Intake', N = 'totalN',
                          gvar_type_OUT = 'd',
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 
PLOT_NO.PING.PONG(nppOutput)


# changes in social influence (perceived norms)
nppOutput <- NO.PING.PONG(data_Alcohol_Intake, ES_type_IN='g', ES_type_OUT='g', 
                          rma_method='REML', 
                          Bayes_type = c('generated', 'Schmidt_Raju'), 
                          prior_type='MA', CI = 95,
                          ES = 'g.Social.Influence', N = 'totalN',
                          gvar_type_OUT = 'd',
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 
PLOT_NO.PING.PONG(nppOutput)
}
}
