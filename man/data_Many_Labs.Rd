\name{data_Many_Labs}
\alias{data_Many_Labs}
\docType{data}
\title{data_Many_Labs}
\description{A list with summary data for eight research questions examined
	in the Many Labs Replication Project
}
\usage{data(data_Many_Labs)}
\details{
A list with summary data for eight pre-registered datasets described by 
Klein et al. (2014; \cr {https://osf.io/wx7ck/).} 
The datasets were from the Many Labs Replication Project.
For all eight projects, there were attempts to replicate the original effect across 
36 samples (studies), involving over 6,000 participants.

The data are for: (1) Sunk Costs (Oppenheimer, Meyvis, & Davidenko, 2009), 
wherein participants are supposedly more likely to go see their favorite football 
team play an important game on a freezing cold day if they had paid for the ticket 
than if the ticket had been free; (2) Anchoring (Jacowitz & Kahneman, 1995), wherein 
participants subsequent judgments are supposedly affected by anchoring information 
provided by researchers; (3) Gamblers Fallacy (Oppenheimer & Monin, 2009), wherein 
the rarity of an independent, chance observation supposedly influences beliefs about 
what occurred before the event; (4) Quote Attribution (Lorge & Curtiss, 1936), 
wherein participants opinions are supposedly more influenced by whether they are 
exposed to the opinions of a highly regarded than a lowly regarded other person; 
(5) Flag Priming (Carter, Ferguson, & Hassin, 2011; Study 2), wherein subtle 
exposure to the American flag supposedly increases conservative opinions among 
US participants; (6) Currency Priming (Caruso, Vohs, Baxter, & Waytz, 2013), 
wherein mere exposure to money symbols increases participants endorsements of 
the current social system; (7) Imagined Contact (Husnu & Crisp, 2010; Study 1), 
wherein merely imagining contact with members of ethnic outgroups is supposedly 
sufficient to reduce prejudice toward those groups; and (8) Math Attitudes (Nosek, 
Banaji, & Greenwald, 2002), wherein women supposedly have more negative implicit 
attitudes toward math compared to arts than do men.
}

\references{
	{Carter, T. J., Ferguson, M. J., & Hassin, R. R. (2011). A single exposure to the 
	American flag shifts support toward Republicanism up to 8 months later. 
	\emph{Psychological science, 22(8),} 1011-1018.}
\cr\cr {Caruso, E. M., Vohs, K. D., Baxter, B., & Waytz, A. (2013). Mere exposure to 
	money increases endorsement of free-market systems and social in- equality. 
	\emph{Journal of Experimental Psychology: General, 142,} 301-306.}
\cr\cr {Husnu, S., & Crisp, R. J. (2010). Elaboration enhances the imagined contact 
	effect. \emph{Journal of Experimental Social Psychology, 46(6),} 943-950.}
\cr\cr {Jacowitz, K. E., & Kahneman, D. (1995). Measures of anchoring in estimation tasks. 
	\emph{Personality and Social Psychology Bulletin, 21(11),} 1161-1166.}
\cr\cr {Klein, R. A., et al. (2014). Investigating variation in replicability: A many 
	labs replication project. \emph{Social Psychology, 45,} 142152.}
\cr\cr {Lorge, I., & Curtiss, C. C. (1936). Prestige, suggestion, and attitudes. 
	\emph{Journal of Social Psychology, 7(4),} 386-402.}
\cr\cr {Nosek, B. A., Banaji, M. R., & Greenwald, A. G. (2002). Math = Male, Me = Female, 
	therefore Math (not equal) Me. \emph{Journal of Personality and Social Psychology, 83(1),} 44-59.}
\cr\cr {Oppenheimer, D. M., & Monin, B. (2009). The retrospective gamblers fallacy: 
	Unlikely events, constructing the past, and multiple universes. 
	\emph{Judgment and Decision Making, 4,} 326-334.}
\cr\cr {Oppenheimer, D. M., Meyvis, T., & Davidenko, N. (2009). Instructional manipulation checks: 
	Detecting satisficing to increase statistical power. 
	\emph{Journal of Experimental Social Psychology, 45(4),} 867-872.}
}

\examples{
names(data_Many_Labs)

\donttest{ 

# the Anchoring data
head(data_Many_Labs$Anchoring)
nppOutput <- NO.PING.PONG(data_Many_Labs$Anchoring, ES_type_IN=NULL, ES_type_OUT='r', 
                          rma_method='REML',
                          Bayes_type = c('generated', 'raw', 'Schmidt_Raju'), 
                          prior_type='MA', CI = 95,
                          ES = NULL, N = NULL, vi = NULL,
                          grp1_mn = 'grp2_mn', grp1_sd = 'grp2_sd', grp1_n = 'grp2_N', 
                          grp2_mn = 'grp1_mn', grp2_sd = 'grp1_sd', grp2_n = 'grp1_N', 
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 
PLOT_NO.PING.PONG(nppOutput)
	
	
# the Gamblers_Fallacy data
head(data_Many_Labs$Gamblers_Fallacy)
nppOutput <- NO.PING.PONG(data_Many_Labs$Gamblers_Fallacy, ES_type_IN=NULL, ES_type_OUT='r', 
                          rma_method='REML',
                          Bayes_type = c('generated', 'raw', 'Schmidt_Raju'),
                          prior_type='MA', CI = 95,
                          ES = NULL, N = NULL, vi = NULL,
                          grp1_mn = 'grp2_mn', grp1_sd = 'grp2_sd', grp1_n = 'grp2_N', 
                          grp2_mn = 'grp1_mn', grp2_sd = 'grp1_sd', grp2_n = 'grp1_N', 
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 
PLOT_NO.PING.PONG(nppOutput)
		
		
# the Math_Attitudes data
head(data_Many_Labs$Math_Attitudes)
nppOutput <- NO.PING.PONG(data_Many_Labs$Math_Attitudes, ES_type_IN=NULL, ES_type_OUT='r', 
                          rma_method='REML',
                          Bayes_type = c('generated', 'raw', 'Schmidt_Raju'), 
                          prior_type='MA', CI = 95,
                          ES = NULL, N = NULL, vi = NULL,
                          grp1_mn = 'grp2_mn', grp1_sd = 'grp2_sd', grp1_n = 'grp2_N', 
                          grp2_mn = 'grp1_mn', grp2_sd = 'grp1_sd', grp2_n = 'grp1_N', 
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 
PLOT_NO.PING.PONG(nppOutput)
	
	
# the Sunk_Costs data
head(data_Many_Labs$Sunk_Costs)
nppOutput <- NO.PING.PONG(data_Many_Labs$Sunk_Costs, ES_type_IN=NULL, ES_type_OUT='r', 
                          rma_method='REML', 
                          Bayes_type = c('generated', 'raw', 'Schmidt_Raju'), 
                          prior_type='MA', CI = 95,
                          ES = NULL, N = NULL, vi = NULL,
                          grp1_mn = 'grp1_mn', grp1_sd = 'grp1_sd', grp1_n = 'grp1_N', 
                          grp2_mn = 'grp2_mn', grp2_sd = 'grp2_sd', grp2_n = 'grp2_N', 
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 
PLOT_NO.PING.PONG(nppOutput)
	
		
# the Quote_Attribution data
head(data_Many_Labs$Quote_Attribution)
nppOutput <- NO.PING.PONG(data_Many_Labs$Quote_Attribution, ES_type_IN=NULL, ES_type_OUT='r', 
                          rma_method='REML', 
                          Bayes_type = c('generated', 'raw', 'Schmidt_Raju'), 
                          prior_type='MA', CI = 95,
                          ES = NULL, N = NULL, vi = NULL,
                          grp1_mn = 'grp2_mn', grp1_sd = 'grp2_sd', grp1_n = 'grp2_N', 
                          grp2_mn = 'grp1_mn', grp2_sd = 'grp1_sd', grp2_n = 'grp1_N', 
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 
PLOT_NO.PING.PONG(nppOutput)
	
	
# the Flag_Priming data
head(data_Many_Labs$Flag_Priming)
nppOutput <- NO.PING.PONG(data_Many_Labs$Flag_Priming, ES_type_IN=NULL, ES_type_OUT='r', 
                          rma_method='REML', 
                          Bayes_type = c('generated', 'raw', 'Schmidt_Raju'), 
                          prior_type='MA', CI = 95,
                          ES = NULL, N = NULL, vi = NULL,
                          grp1_mn = 'grp2_mn', grp1_sd = 'grp2_sd', grp1_n = 'grp2_N', 
                          grp2_mn = 'grp1_mn', grp2_sd = 'grp1_sd', grp2_n = 'grp1_N', 
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 
PLOT_NO.PING.PONG(nppOutput)
	
	
# the Currency_Priming data
head(data_Many_Labs$Currency_Priming)
nppOutput <- NO.PING.PONG(data_Many_Labs$Currency_Priming, ES_type_IN=NULL, ES_type_OUT='r', 
                          rma_method='REML', 
                          Bayes_type = c('generated', 'raw', 'Schmidt_Raju'), 
                          prior_type='MA', CI = 95,
                          ES = NULL, N = NULL, vi = NULL,
                          grp1_mn = 'grp2_mn', grp1_sd = 'grp2_sd', grp1_n = 'grp2_N', 
                          grp2_mn = 'grp1_mn', grp2_sd = 'grp1_sd', grp2_n = 'grp1_N', 
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 
PLOT_NO.PING.PONG(nppOutput)

	
# the Imagined_Contact data
head(data_Many_Labs$Imagined_Contact)
nppOutput <- NO.PING.PONG(data_Many_Labs$Imagined_Contact, ES_type_IN=NULL, ES_type_OUT='r', 
                          rma_method='REML', 
                          Bayes_type = c('generated', 'raw', 'Schmidt_Raju'), 
                          prior_type='MA', CI = 95,
                          ES = NULL, N = NULL, vi = NULL,
                          grp1_mn = 'grp2_mn', grp1_sd = 'grp2_sd', grp1_n = 'grp2_N', 
                          grp2_mn = 'grp1_mn', grp2_sd = 'grp1_sd', grp2_n = 'grp1_N', 
                          nitt=13000, burnin=3000, thin=10, verbose=TRUE) 
PLOT_NO.PING.PONG(nppOutput)	
}
}
