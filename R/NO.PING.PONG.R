 

NO.PING.PONG <- function(donnes, 
                         ES_type_IN = NULL, ES_type_OUT = 'r', 
                         rawdata_type = 'for_correl',                          
                         ma_method = 'HS',
                         cor_stat = 'ZCOR',
                         Bayes_type = c('Schmidt_Raju'), 
                         prior_type = 'META', 
                         ES = NULL, N = NULL, ES_var = NULL,
                         CI_level_out = 95, CI_level_in = 95, CI_in_lb = NULL, CI_in_ub = NULL,
                         grp1_mn = NULL, grp1_sd = NULL, grp1_n = NULL, 
                         grp2_mn = NULL, grp2_sd = NULL, grp2_n = NULL, 
                         gvar_type_OUT = 'd',
                         paired_samples_ES_type = NULL,
                         funnel_plot=FALSE, funnel_plot_type='png', funnel_plot_title=NULL,
                         nitt = 53000, burnin = 3000, thin = 10, 
                         verbose = TRUE) {


if (is.null(ES_type_IN))              ES_type_IN <- 'unspecified'
if (is.null(paired_samples_ES_type))  paired_samples_ES_type <- 'SMCRH'


# determine if donnes is a list (= if it is raw data)
dontype <- ifelse (inherits(donnes, "list"), 'raw data', 'matrix') 


# if donnes is a matrix, convert it to a dataframe
if (dontype == 'matrix' & !is.data.frame(donnes))  donnes <- as.data.frame(donnes)


# the critical z value (that corresponds to the specified CI) to be used in the CI computations
zforCI_out <- qnorm((1 + CI_level_out * .01) / 2) 
zforCI_in  <- qnorm((1 + CI_level_in  * .01) / 2) 



# is complete info (M, SD, & N) for 2 groups provided?
if (!is.null(grp1_mn) & !is.null(grp1_sd) & !is.null(grp1_n) &
    !is.null(grp2_mn) & !is.null(grp2_sd) & !is.null(grp2_n)  ) 
    {grpinfoALL = TRUE} else {grpinfoALL = FALSE}

if (is.null(N) & !is.null(grp1_n) & !is.null(grp2_n)) {
	donnes$N <- donnes[,grp1_n] + donnes[,grp2_n]
	N <- 'N'
}


# warning if both N & grp Ns are missing & dontype is not a list
if (dontype == 'matrix' & is.null(N) & is.null(grp1_n) & is.null(grp2_n)) {
	message('\nThe entered data is a matrix (and not raw data points), and')		
	message('both "N", and "grp1_n" & "grp2_n", are missing/NULL. Bayesian analyses cannot')
	message('be conducted without sample size information.\n')
}


# is one of CI_in_lb or CI_in_ub available
ci_dat <- NULL
if (!is.null(CI_in_lb))                      ci_dat <- donnes[,CI_in_lb]
if (is.null(CI_in_lb) & !is.null(CI_in_ub))  ci_dat <- donnes[,CI_in_ub]



######################  compute ESdat for meta-analyses  &  donnes_RN for Bayesian  #####################


ESdat <- NULL   # ESdat will have the yi & vi values for the meta-analyses

ESdat_type <- 'r'   # the ES type of yi - will be r most of the time

donnes_RN <- NULL # has r & N for the Bayesian analyses; won't do Bayesian if this remains null


	
if (dontype == 'matrix' & ES_type_IN == 'r') {

	# donnes_RN & ESdat if ES & N are provided
	if (!is.null(ES) & !is.null(N) ) {

		donnes_RN <- data.frame(N = donnes[,N], ES = donnes[,ES])

		ESdat <- escalc(measure=cor_stat, ni=donnes_RN$N, ri=donnes_RN$ES )
	}

	# donnes_RN & ESdat if ES & ES_var are provided, but not N
	if (!is.null(ES) & !is.null(ES_var) & is.null(N)) {
		
		donnes[,N] <- (1 - donnes[,ES]^2) / donnes[,ES_var]  + 2
				
		donnes_RN <- data.frame(N = donnes[,N], ES = donnes[,ES])

		ESdat <- escalc(measure=cor_stat, ni=donnes_RN$N, ri=donnes_RN$ES )
	}

	# when ES & CIs are available, but not ES_var or N  
	# PROB = CIs for r are not symmetric & should have been computed from z, & can't be sure it was done
	if (!is.null(ES) & !is.null(ci_dat) & is.null(N) & is.null(ES_var)) {
		
		message('\n\nr and confidence intervals were provided, but not the study Ns or r variances.')
		message('The analyses cannot proceed because ES_var cannot be reliably computed from the provided data.\n')	
		
		# # using one of CI_in_lb or CI_in_ub
		# if (!is.null(CI_in_lb))                      ci_dat <- donnes[,CI_in_lb]
		# if (is.null(CI_in_lb) & !is.null(CI_in_ub))  ci_dat <- donnes[,CI_in_ub]
		
		# es_var <- ( abs(donnes[,ES] - ci_dat) / zforCI_in )^2

		# donnes[,N] <- (1 - donnes[,ES]^2) / es_var  + 2
				
		# donnes_RN <- data.frame(N = donnes[,N], ES = donnes[,ES])

		# ESdat <- escalc(measure=cor_stat, ni=donnes_RN$N, ri=donnes_RN$ES )
	}		
}



if (dontype == 'matrix' & ES_type_IN == 'z') {

	# donnes_RN & ESdat if ES & N are provided
	if (!is.null(ES) & !is.null(N) ) {
		
		# convert z to r
		donnes_RN <- data.frame(N = donnes[,N], ES = tanh(donnes[,ES])) 
						
		ESdat <- escalc(measure=cor_stat, ni=donnes_RN$N, ri=donnes_RN$ES )
	}
}



if (dontype == 'matrix' & (ES_type_IN == 'd' | ES_type_IN == 'g')) {

	# compute ES_var if ES & CIs are provided, but not ES_var or N  
	if (is.null(ES_var) & is.null(N) & !is.null(ci_dat)) {
	 	donnes$ES_var <- (abs(donnes[,ES] - ci_dat) / zforCI_in)^2
		ES_var <- 'ES_var'
	}	
	
	rfromdg <- CONVERT_ES(donnes, ES = ES, ES_type_IN = ES_type_IN, ES_var = ES_var, 
	                      totalN = N, grp1_n = grp1_n, grp2_n = grp2_n, 
	                      CI_in_lb = CI_in_lb, CI_in_ub = CI_in_ub, verbose = FALSE)

	# 
	# CONVERT_ES <- function(donnes, ES, ES_type_IN='r', ES_var = NULL, 
	#                        totalN = NULL, grp1_n = NULL, grp2_n = NULL,
	#                        gvar_type_OUT = 'd', 
	#                        CI_level_out = 95, CI_level_in = 95, CI_in_lb = NULL, CI_in_ub = NULL, 
	#                        verbose = TRUE) {
	# 
	#   CI_in_lb = NULL, CI_in_ub = NULL	  
	#   
	  
	# if have totalN
	if (all(!is.na(rfromdg$totalN))) {
		
		donnes_RN <- data.frame(N = rfromdg$totalN, ES = rfromdg$r)
	
		ESdat <- escalc(measure=cor_stat, ni=donnes_RN$N, ri=donnes_RN$ES )
	} 

	# if do not have totalN, but CIs were provided
	if (all(is.na(rfromdg$totalN))) {
				
		# r & r_var can be computed from d & d_var
		if (ES_type_IN == 'd')  ESdat <- escalc(yi=rfromdg$r, vi=rfromdg$r_var)

		# r & r_var canNOT be computed from g & g_var, & ESdat will be computed in a section below
	}
}



# grp1 & grp2 info is available
if (dontype == 'matrix' & (grpinfoALL)) { 
	
	ESdat_grpinfo <- ES_from_GRPINFO(
	                    grp1_MN=donnes[,grp1_mn], grp1_SD=donnes[,grp1_sd], grp1_n=donnes[,grp1_n], 
		                grp2_MN=donnes[,grp2_mn], grp2_SD=donnes[,grp2_sd], grp2_n=donnes[,grp2_n],
		                ES_type_OUT = 'r')

	# from the escalc documentation: The positive bias in the standardized mean difference is 
	# automatically corrected for within the function, yielding Hedges' g for measure="SMD" (Hedges, 1981)	
	# I am therefore using ES_from_GRPINFO, because donnes can sometimes be d and not g values	

	donnes_RN <- data.frame(N = ESdat_grpinfo$totalN, ES = ESdat_grpinfo$r) 

	ESdat <- escalc(measure=cor_stat, ni=donnes_RN$N, ri=donnes_RN$ES )
}		




# converting OR to r or d or g   -- need CIs,  and Ns are required for donnes_RN
if (dontype == 'matrix' & ES_type_IN == 'OR') {

	# need CIs to obtain logOR_var
	if (!is.null(CI_in_lb) | !is.null(CI_in_ub) ) {

		# using one of CI_in_lb or CI_in_ub
		if (!is.null(CI_in_lb))                      donnes$ci_dat <- donnes[,CI_in_lb]
		if (is.null(CI_in_lb) & !is.null(CI_in_ub))  donnes$ci_dat <- donnes[,CI_in_ub]
			
		# have Ns
		if (!is.null(N)) {
							
			outp <- CONVERT_ES(donnes, ES=ES, ES_type_IN = 'OR', CI_in_lb = 'ci_dat', totalN = 'N', 
			                   verbose = FALSE)

			donnes_RN <- data.frame(N = outp$totalN, ES = outp$r)
				
			ESdat <- escalc(measure=cor_stat, ni=donnes_RN$N, ri=donnes_RN$ES )
		}

		# do not have Ns -- can obtain r & d, but not g or donnes_RN
		if (is.null(N)) {
							
			outp <- CONVERT_ES(donnes, ES=ES, ES_type_IN = 'OR', CI_in_lb = 'ci_dat',  
			                   verbose = FALSE)
			                   
			ESdat <- escalc(yi=outp$r, vi=outp$r_var)
		}
	}

	if (is.null(CI_in_lb) | is.null(CI_in_ub) ) {
		message('\n\nThe input effect sizes are odds ratios but one of CI_in_lb or CI_in_ub is')
		message('also required for the effect size variances to be computed. These values are')
		message('missing and so the analyses cannot be conducted. Errors will be produced.\n')
	}
}




if (dontype == 'raw data' & rawdata_type == 'for_correl') {
	donnes_RN <- data.frame(N = rep(NA,length(donnes)), ES = rep(NA,length(donnes)))
	for (lupe in 1:length(donnes)) {
		donnes_RN$N[lupe]  <- nrow(donnes[[lupe]])
		donnes_RN$ES[lupe] <- cor(na.omit(donnes[[lupe]]))[2,1]
	}
	ESdat <- escalc(measure=cor_stat, ni=donnes_RN$N, ri=donnes_RN$ES )	
}



if (dontype == 'raw data' & rawdata_type == 'indep_groups') {
			
	# ESdat_grpinfo  <- data.frame(yi = rep(NA,length(donnes)), vi = rep(NA,length(donnes)))
	donnes_RN <- data.frame(N =  rep(NA,length(donnes)), ES = rep(NA,length(donnes)))
	for (lupe in 1:length(donnes)) {
		
		means <- aggregate(as.matrix(donnes[[lupe]][,2]), list(donnes[[lupe]][,1]), mean)			
		SDs   <- aggregate(as.matrix(donnes[[lupe]][,2]), list(donnes[[lupe]][,1]), sd)
		Ns    <- table(donnes[[lupe]][,1])
			
		# using just the 1st 2 values of the groups variable (col 1 of donnes)
		ESdat_grpinfo[lupe,] <- ES_from_GRPINFO(
		                           grp1_MN=means[1,2], grp1_SD=SDs[1,2], grp1_n=Ns[1], 
			                       grp2_MN=means[2,2], grp2_SD=SDs[2,2], grp2_n=Ns[2],
			                       ES_type_OUT='r')

		# from the escalc Rd:
		# The positive bias in the standardized mean difference is automatically 
		# corrected for within the function, yielding Hedges' g for measure="SMD" 
		# (Hedges, 1981). Similarly, the same bias correction is applied for 
		# measure="SMDH" (Bonett, 2009). For measure="SMD", one can choose between 
		# vtype="LS" (the default) and vtype="UB". The former uses the usual large-sample 
		# approximation to compute the sampling variances. The latter provides unbiased 
		# estimates of the sampling variances.			
							
		donnes_RN[lupe,] <- data.frame(N = ESdat_grpinfo$totalN, ES = ESdat_grpinfo$r) 
	}
	ESdat <- escalc(measure=cor_stat, ni=donnes_RN$N, ri=donnes_RN$ES )
}



if (dontype == 'raw data' & rawdata_type == 'paired_samples') {
		
		studyStats <- c()
		for (lupe in 1:length(donnes)) {

			studyMNs <- colMeans(donnes[[lupe]])						
			studySDs <- apply(donnes[[lupe]], 2, sd)			
			studyN   <- nrow(donnes[[lupe]])
			studyR   <- cor(donnes[[lupe]])[1,2]

			studyStats <- rbind(studyStats, c(studyMNs, studySDs, studyN, studyR))
		}
		colnames(studyStats) <- c('studyMN1','studyMN2','studySD1','studySD2','studyN','studyR')
		studyStats <- data.frame(studyStats)
		
		ESdat <- escalc(measure=paired_samples_ES_type,  
				        m1i=studyStats$studyMN1,  sd1i=studyStats$studySD1,   
					    m2i=studyStats$studyMN2,  sd2i=studyStats$studySD2, 
					    ni=studyStats$studyN, ri=studyStats$studyR )
					    
		ESdat_type <- ES_type_IN <- ES_type_OUT <- 'unspecified'
}



# when only ES, ES_type_IN, & ES_var are provided 
if (dontype == 'matrix' & is.null(ESdat) & is.null(N) & !grpinfoALL & !is.null(ES) & !is.null(ES_type_IN) & !is.null(ES_var)) {

	# if (ES_type_IN == 'r')   # not necessary, covered above

	# if (ES_type_IN == 'OR')  # not necessary, covered above

	if (ES_type_IN == 'z' | ES_type_IN == 'd' | ES_type_IN == 'g')  ESdat_type <- ES_type_IN
	
	if (ES_type_IN == 'g') {
		
		message('\nOnly the g effect size and the effect size variance are available.')
		message('The meta-analysis output will therefore be in the g metric.')

		ES_type_OUT <- ES_type_IN
	}
		
	ESdat <- escalc(yi=donnes[,ES], vi=donnes[,ES_var])
}



# when only ES & ES_var are provided & ES_type_IN = NULL
if (dontype == 'matrix' & is.null(ESdat) & is.null(N) & is.null(ES_type_IN) & !grpinfoALL & !is.null(ES) & !is.null(ES_var)) {
	
	ESdat <- escalc(yi=donnes[,ES], vi=donnes[,ES_var])
	
	ESdat_type <- ES_type_IN
	
	message('\nOnly the effect size and the effect size variance were provided as input.
The meta-analysis output will therefore be in the metric of the input effect size, which is: ', ES_type_IN)

	ES_type_OUT <- ES_type_IN

}

	
# the total Ns, if available
if (!is.null(donnes_RN)) {totalNs <- donnes_RN$N} else {totalNs <- rep(NA,nrow(ESdat))}


# identifying & removing any row with an NA from ESdat,  & the same rows from totalNs, & donnes_RN
rowswithNAs <- unique(which(is.na(cbind(ESdat)), arr.ind=TRUE)[,1])  
if (length(rowswithNAs) >= 1) {
	ESdat <- ESdat[-rowswithNAs,]
	donnes_RN <- donnes_RN[-rowswithNAs,] 
	totalNs <- totalNs[-rowswithNAs]
	message('\n\nMissing values were found and removed from the data matrix.\n\n')
}

                                       
 Nstudies <- nrow(ESdat)
 
               

 
######################  data & argument specifications  #####################
 
if (verbose) {

	if (dontype == 'raw data') {
		
		message('\n\nThe input data are raw data (and not study effect sizes)')
		
		if (rawdata_type == 'for_correl') 
			message('\nThe raw data were specified to be for computations of correlations between two variables.')
	
		if (rawdata_type == 'indep_groups') 
			message('\nThe raw data were specified to be based on independent groups (in the first column of donnes).')
	
		if (rawdata_type == 'paired_samples') {
			message('\nThe raw data were specified to be based on paired samples, with paired_samples_ES_type = ', 
			        paired_samples_ES_type)			
		}	
	}

	if (dontype == 'matrix') {

		message('\n\nThe input data are a matrix of study effect sizes (not raw data)')
		
		message('\nThe specified type of input effect sizes is: ', ES_type_IN)
		
		if (grpinfoALL == TRUE) 
			message('\nComplete group info (i.e., M, SD, & N) for 2 groups was provided for the analyses')
			
		if (ES_type_IN == 'OR' & is.null(donnes_RN)) {
			message('\nThe input effect sizes were odds ratios but N was not provided. The odds ratios were')
			message('converted to correlation coefficients for the meta-analyses. Bayesian analyses are')
			message('not possible without the Ns.')
		}
	}

	message('\nNumber of studies: ', Nstudies)

	# message('\nThe type of output effect size is: ', ES_type_OUT)
	
	if (ES_type_OUT == 'r' & ESdat_type == 'r')  
		message('\nThe requested kind of statistic (cor_stat) for the meta-analyses: ', cor_stat)
	
	message('\nThe specified method option for the meta-analyses: ',  ma_method)

	message('\nThe requested kind of output effect size: ', ES_type_OUT)

	if (!is.null(donnes_RN)) {
		
		message('\nThe requested kind(s) of Bayesian analyses: ', paste(Bayes_type, collapse=", "))
	
		message('\nThe type of prior data used in the updating analyses (both cumulative and Bayesian): ',  prior_type)
	} 

	message('\nThe confidence interval for the output: ', CI_level_out, '%')
}


               
#############################################  NHST  #######################################


summESdat <- summary(ESdat)  # get CIs, only works if escalc was used

# place results in ES_type_OUT metric
if (ESdat_type == 'r') {

	if (ES_type_OUT == 'r' & cor_stat == 'COR')	
		results_NHST <- cbind(1:Nstudies, summESdat$ci.lb, summESdat$yi, summESdat$ci.ub)
	
	if (ES_type_OUT == 'r' & cor_stat == 'ZCOR') {
		zes <- cbind(summESdat$ci.lb, summESdat$yi, summESdat$ci.ub)
		res <- tanh(zes)
		results_NHST <- cbind(1:Nstudies, res)
	}
	
	if (ES_type_OUT == 'z' & cor_stat == 'COR')	{
		res <- cbind(summESdat$ci.lb, summESdat$yi, summESdat$ci.ub)
		zes <- atanh(res)
		results_NHST <- cbind(1:Nstudies, zes)
	}
	
	if (ES_type_OUT == 'z' & cor_stat == 'ZCOR')
		results_NHST <- cbind(1:Nstudies, summESdat$ci.lb, summESdat$yi, summESdat$ci.ub)
	
	if (ES_type_OUT == 'd' & cor_stat == 'COR') {
		res <- cbind(summESdat$ci.lb, summESdat$yi, summESdat$ci.ub)
		des <- (2 * res) / sqrt(1 - res^2)		
		results_NHST <- cbind(1:Nstudies, des)
	}
	
	if (ES_type_OUT == 'd' & cor_stat == 'ZCOR') {
		zes <- cbind(summESdat$ci.lb, summESdat$yi, summESdat$ci.ub)
		res <- tanh(zes)
		des <- (2 * res) / sqrt(1 - res^2)		
		results_NHST <- cbind(1:Nstudies, des)
	}
	
	 if (ES_type_OUT == 'g' & cor_stat == 'COR' & all(!is.na(totalNs))) {
		res <- cbind(summESdat$ci.lb, summESdat$yi, summESdat$ci.ub)
		des <- (2 * res) / sqrt(1 - res^2)		
		ges <- g_from_d(d = des[,1], totalN = totalNs)[1]
		ges[,2] <- g_from_d(d = des[,2], totalN = totalNs)[1]
		ges[,3] <- g_from_d(d = des[,3], totalN = totalNs)[1]
		results_NHST <- cbind(1:Nstudies, ges)
	}
	
	if (ES_type_OUT == 'g' & cor_stat == 'ZCOR' & all(!is.na(totalNs))) {
		zes <- cbind(summESdat$ci.lb, summESdat$yi, summESdat$ci.ub)
		res <- tanh(zes)
		des <- (2 * res) / sqrt(1 - res^2)		
		ges <- g_from_d(d = des[,1], totalN = totalNs)[1]
		ges[,2] <- g_from_d(d = des[,2], totalN = totalNs)[1]
		ges[,3] <- g_from_d(d = des[,3], totalN = totalNs)[1]
		results_NHST <- cbind(1:Nstudies, ges)
	}

	results_NHST <- data.frame(results_NHST)
	colnames(results_NHST) <- c('Study','ES_lb','ES','ES_ub')

	if (verbose) {
		message('\n\nStudy effect sizes and conventional confidence intervals (ES metric = ', ES_type_OUT,'):\n')
		print(round(results_NHST,3), print.gap=4, row.names = FALSE)
	}
} 


if (ESdat_type == 'z') {  # totalNs are not available, so g not possible

	zes <- cbind(summESdat$ci.lb, summESdat$yi, summESdat$ci.ub)

	if (ES_type_OUT == 'z')  results_NHST <- cbind(1:Nstudies, zes)
		
	if (ES_type_OUT == 'r') {
		res <- tanh(zes)
		results_NHST <- cbind(1:Nstudies, res)
	}
	
	if (ES_type_OUT == 'd') {
		res <- tanh(zes)
		des <- (2 * res) / sqrt(1 - res^2)		
		results_NHST <- cbind(1:Nstudies, des)
	}
	
	# dimnames(results_NHST) <-list(rep("", dim(results_NHST)[1]))
	results_NHST <- data.frame(results_NHST)
	colnames(results_NHST) <- c('Study','ES_lb','ES','ES_ub')

	if (verbose) {
		message('\n\nStudy effect sizes and confidence intervals (ES metric = ', ES_type_OUT,'):\n')
		print(round(results_NHST,3), print.gap=4, row.names = FALSE)
	}
} 


if (ESdat_type == 'd') {  # totalNs are not available, so g not possible

	des <- cbind(summESdat$ci.lb, summESdat$yi, summESdat$ci.ub)

	if (ES_type_OUT == 'd')  results_NHST <- cbind(1:Nstudies, des)
		
	if (ES_type_OUT == 'r')  {  
		res <- r_from_d(d = des, d_var = ESdat$vi)[1:3]		
		results_NHST <- cbind(1:Nstudies, res)
	}
	
	if (ES_type_OUT == 'z') {
		res <- r_from_d(d = des, d_var = ESdat$vi)[1:3]		
		zes <- atanh(res)
		results_NHST <- cbind(1:Nstudies, zes)
	}
	
	results_NHST <- data.frame(results_NHST)
	colnames(results_NHST) <- c('Study','ES_lb','ES','ES_ub')

	if (verbose) {
		message('\n\nStudy effect sizes and confidence intervals (ES metric = ', ES_type_OUT,'):\n')
		print(round(results_NHST,3), print.gap=4, row.names = FALSE)
	}
} 


if (ESdat_type != 'r' & ESdat_type != 'z' & ESdat_type != 'd') {
	results_NHST <- data.frame(cbind(1:Nstudies, summESdat$ci.lb, summESdat$yi, summESdat$ci.ub))
	colnames(results_NHST) <- c('Study','ES_lb','ES','ES_ub')

	if (verbose) {
		message('\n\nStudy effect sizes and confidence intervals (ES metric = ', ES_type_IN,'):\n')
		print(round(results_NHST,3), print.gap=4, row.names = FALSE)
	}
}



####################################  regular meta-analysis  #######################################

outp_MA_REG <- rma(yi=ESdat$yi, vi=ESdat$vi, method=ma_method)  # random-effects model


if (ma_method != 'FE') {
	# confidence intervals for tau & Isq
	heteroests <- confint(outp_MA_REG)
	
	tau2    <- heteroests[[1]][1,1]
	tau2LB  <- heteroests[[1]][1,2]
	tau2UB  <- heteroests[[1]][1,3]
	
	tau     <- heteroests[[1]][2,1]
	tauLB   <- heteroests[[1]][2,2]
	tauUB   <- heteroests[[1]][2,3]
	
	isq     <- heteroests[[1]][3,1] * .01   # because rma outputs as percentages
	isqLB   <- heteroests[[1]][3,2] * .01   # because rma outputs as percentages
	isqUB   <- heteroests[[1]][3,3] * .01   # because rma outputs as percentages
	
	hsq     <- heteroests[[1]][4,1]
	hsqLB   <- heteroests[[1]][4,2]
	hsqUB   <- heteroests[[1]][4,3]
}

if (ma_method == 'FE') {   # for the output object
	tau2 = NULL; tau2LB = NULL; tau2UB = NULL;
	tau = NULL; tauLB = NULL; tauUB = NULL;
	isq = NULL; isqLB = NULL; isqUB = NULL;
	hsq = NULL; hsqLB = NULL; hsqUB = NULL;
}


if (verbose) {

	if (ESdat_type == 'r') {
		
		if (cor_stat == 'COR') { 	
			res <- cbind(outp_MA_REG$b, outp_MA_REG$ci.lb, outp_MA_REG$ci.ub)
			esdat <- escalc(yi = outp_MA_REG$b, vi = outp_MA_REG$vb) # getting z from r using escalc
			esdat <- summary(esdat, level = CI_level_out, transf=transf.rtoz)
			zes <- cbind(esdat$yi, esdat$ci.lb, esdat$ci.ub)			
		}
	
		if (cor_stat == 'ZCOR') { 	
			zes <- cbind(outp_MA_REG$b, outp_MA_REG$ci.lb, outp_MA_REG$ci.ub)
			esdat <- escalc(yi = outp_MA_REG$b, vi = outp_MA_REG$vb) # getting r from z using escalc
			esdat <- summary(esdat, level = CI_level_out, transf=transf.ztor)
			res <- cbind(esdat$yi, esdat$ci.lb, esdat$ci.ub)			
		}	
	  
		des <- d_from_r(r = res)[1]    # des <- r_to_d(esmat)
	
		ges <- g_from_d(d=des, totalN = outp_MA_REG$k)[1]   # ges <- d_to_g(d=des, N=outp_MA_REG$k)
	
		esmat <- rbind(res, zes, t(des), t(ges))
	 	colnames(esmat) <- c('ES','ES_lb','ES_ub')
		rownames(esmat) <- c('r','z','d','g')
	
		message('\n\nResults from a regular, all-studies-at-once, random-effects meta-analysis')
		message('in r, Fisher\'s z, d, & g effect size metrics:\n')
		print(round(esmat[,c('ES','ES_lb','ES_ub')],3), print.gap=4)	
	}


	if (ESdat_type == 'z') {  # totalNs are not available, so g not possible
	
		zes <- cbind(outp_MA_REG$b, outp_MA_REG$ci.lb, outp_MA_REG$ci.ub)
	
		res <- tanh(zes)
		
		des <- d_from_r(r = res)[1]
		
		esmat <- rbind(res, zes, t(des))
	 	colnames(esmat) <- c('ES','ES_lb','ES_ub')
		rownames(esmat) <- c('r','z','d')
	
		message('\n\nResults from a regular, all-studies-at-once, random-effects meta-analysis')
		message('in r, Fisher\'s z, & d effect size metrics:\n')
		print(round(esmat[,c('ES','ES_lb','ES_ub')],3), print.gap=4)	
	} 


	if (ESdat_type == 'd') {  # totalNs are not available, so g not possible
	
		des <- cbind(outp_MA_REG$b, outp_MA_REG$ci.lb, outp_MA_REG$ci.ub)
	
		res <- r_from_d(d = c(des))[,1]
		
		zes <- atanh(res)
				
		esmat <- mapply(c, res, zes, des)
	 	colnames(esmat) <- c('ES','ES_lb','ES_ub')
		rownames(esmat) <- c('r','z','d')
	
		message('\n\nResults from a regular, all-studies-at-once, random-effects meta-analysis')
		message('in r, Fisher\'s z, & d effect size metrics:\n')
		print(round(esmat[,c('ES','ES_lb','ES_ub')],3), print.gap=4)	
	} 


	if (ESdat_type != 'r' & ESdat_type != 'z' & ESdat_type != 'd') {	
		
		esmat <- cbind(outp_MA_REG$b, outp_MA_REG$ci.lb, outp_MA_REG$ci.ub)
		colnames(esmat) <- c('ES','ES_lb','ES_ub')		
		
		if (ES_type_IN == 'OR' & is.null(donnes_RN)) { 
			message('\n\nResults from a regular, all-studies-at-once, random-effects meta-analysis (ES metric = r)\n')
		} else {
		message('\n\nResults from a regular, all-studies-at-once, random-effects meta-analysis (ES metric = ', 
		        ES_type_IN,'):\n')
		}
		rownames(esmat) <- c(' ')
		print(round(esmat[,c('ES','ES_lb','ES_ub')],3), print.gap=4)		
	}
}



##################################  publication bias  ###########################################


# cannot run some tests for publication bias when the study Ns are all the same; using jitter in these cases
	

# Egger's Regression Test for Funnel Plot Asymmetry -- metafor
if (var(ESdat$vi) == 0) {
	# if the vi variance is 0, redo the rma while jittering vi
	outp_MA_REG_j <- rma(yi=ESdat$yi, vi=jitter(ESdat$vi), method=ma_method)  # random-effects model
	funreg <- regtest(outp_MA_REG_j)
} else {
	funreg <- regtest(outp_MA_REG)
}
	
# Rank Correlation Test for Funnel Plot Asymmetry  -- metafor
if (var(ESdat$vi) == 0) {
	funrank <- suppressWarnings(ranktest(outp_MA_REG_j) )
} else {
	funrank <- suppressWarnings(ranktest(outp_MA_REG))
}

# The trim-and-fill method is a nonparametric method to assess selection bias/publication bias.
# The method provides an estimate of the number of missing studies  
# The basic idea of the trim-and-fill method is to add 
# studies to the funnel plot until it becomes symmetric. (from the the meta package documentation, not using)
if (var(ESdat$vi) == 0) {
	trimfill_out_R0 <- trimfill(outp_MA_REG_j, estimator = 'R0')   # metafor
	trimfill_out_L0 <- trimfill(outp_MA_REG_j, estimator = 'L0') #, control=list(maxiter=1000))   # metafor
} else {
	trimfill_out_R0 <- trimfill(outp_MA_REG, estimator = 'R0')   # metafor
	trimfill_out_L0 <- trimfill(outp_MA_REG, estimator = 'L0') #, control=list(maxiter=1000))   # metafor
}

# Test of Excess Significance  -- metafor
tes_res <- metafor::tes(x=ESdat$yi, vi=ESdat$vi)  #,  test="chi2", tes.alternative = "two.sided")	
tes_res_obs <- tes_res$O  # Observed Number of Significant Findings
tes_res_exp <- tes_res$E  # Expected Number of Significant Findings
tes_res_OEratio	<- tes_res$OEratio  # Observed Number / Expected Number
tes_res_theta <- tes_res$theta
pq <- quantile(tes_res$power)		
tes_res_power_median <- pq[3]	# Median Power
tes_res_power_min <- pq[1]	# Minimum Power 
tes_res_power_max <- pq[5]	# Maximum Power
tes_res_pval <- tes_res$pval
tes_res_X2 <- tes_res$X2
tes_res_theta_lim <- tes_res$theta.lim  # Limit Estimate


# metafor - selmodel    -- not using

# # the metabias (& metacor) function is from the meta package, so not using
# m1 <- metacor(cor=as.matrix(donnes_RN[,2]), as.matrix(donnes_RN[,1]),  method.bias='linreg') # -- meta
# thomsharp <- metabias(m1, plotit=F, method='mm') # a variant of Eggers's test allowing for between-study heterogeneity	
# thomsharp$statistic; thomsharp$p.value     Thom/Sharp t','   Thom/Sharp p')
# message('\n   Thom/Sharp:  t = ', round(thomsharp$statistic,2), '   p = ', round(thomsharp$p.value))


biasStats <- list(
	reg_test_z = funreg$zval,
	reg_test_p = funreg$pval,
	rank_test_tau = funrank$tau,
	rank_test_p = funrank$pval,
	tf_missing_studies_R0 = trimfill_out_R0$k0,
	tf_missing_studies_L0 = trimfill_out_L0$k0,
	tes_res_obs = tes_res_obs,
	tes_res_exp = tes_res_exp,
	tes_res_OEratio = tes_res_OEratio,
	tes_res_theta = tes_res_theta,
	tes_res_power_median = tes_res_power_median,
	tes_res_power_min = tes_res_power_min, 
	tes_res_power_max = tes_res_power_max,
	tes_res_pval = tes_res_pval,
	tes_res_X2 = tes_res_X2,
	tes_res_theta_lim = tes_res_theta_lim
)		


if (verbose) {

	message('\n\nTests for Publication Bias:')
	
	message('\n     Regression Test for Funnel Plot Asymmetry:  z = ', round(funreg$zval,2), 
	        '   p = ', round(funreg$pval),4)

	message('\n     Rank Correlation Test for Funnel Plot Asymmetry:  tau = ', 
	        round(funrank$tau,2), '   p = ', round(funrank$pval))

	message('\n     Trim and fill method estimate of the # of missing studies (on the ', 
	        trimfill_out_R0$side,' side) using the R0 method = ',  trimfill_out_R0$k0)

	message('\n     Trim and fill method estimate of the # of missing studies (on the ', 
	        trimfill_out_L0$side,' side) using the L0 method = ',  trimfill_out_L0$k0)
	        
	message('\n     Test of Excess Significance:')
	
	message('\n         Observed Number of Significant Findings: ', tes_res_obs, ' (out of ', tes_res$k, ')')
	
	message('\n         Expected Number of Significant Findings: ', round(tes_res_exp,2))
	
	message('\n         Observed Number / Expected Number: ', round(tes_res_OEratio,2))
	
	message('\n         Estimated Power of Tests (based on theta = ', round(tes_res_theta,2), '):')

	message('\n              Median Power =  ', round(tes_res_power_median,2))
	message('\n              Minimum Power = ', round(tes_res_power_min,2))
	message('\n              Maximum Power = ', round(tes_res_power_max,2))

	message('\n         Test of Excess Significance: p = ', round(tes_res_pval,4), '  (X^2 = ', round(tes_res_X2,3), ', df = 1)' )

	message('\n         Limit Estimate (theta_lim): ', round(tes_res_theta_lim,4), '  (where p = 0.1)' )
}

if (funnel_plot) {

	if (is.null(funnel_plot_title))  funnel_plot_title = deparse(substitute(donnes))

	if (is.null(funnel_plot_type))  funnel_plot_type = 'png'
	
	if (funnel_plot_type == 'bitmap')
		bitmap(paste("Figure - ",funnel_plot_title,".bitmap",sep=""), height=7, width=9, units='in', res=1200, pointsize=12)

	if (funnel_plot_type == 'tiff')
		tiff(paste("Figure - ",funnel_plot_title,".tiff",sep=""), height=7, width=9, units='in', res=1200, pointsize=12)
		
	if (funnel_plot_type == 'png')
		png(paste("Figure - ",funnel_plot_title,".png",sep=""), height=7, width=9, units='in', res=1200, pointsize=12)
		
	if (funnel_plot_type == 'jpeg')
		jpeg(paste("Figure - ",funnel_plot_title,".jpeg",sep=""), height=7, width=9, units='in', res=1200, pointsize=12)
		
	if (funnel_plot_type == 'bmp')
		bmp(paste("Figure - ",funnel_plot_title,".bmp",sep=""), height=7, width=9, units='in', res=1200, pointsize=12)
		        
	# png(paste("Funnel Plot - ",funnel_plot_title,".png",sep=""), width=9, height=7, units="in", res = 600, pointsize=12)
	par( pty="m", mar=c(3,2,3,2) + 2.6)  #     1.8

	oldpar <- par(no.readonly = TRUE)
	on.exit(par(oldpar))

	metafor::funnel(outp_MA_REG, main=funnel_plot_title)
	dev.off()
}

# the commands below can be used to obtain bias-adjusted estimates of meta-analysis effect sizes

# # 2015 Schwarzer - Meta-Analysis with R (book) p. 124
# # all use the meta package
# m1 <- metacor(cor=as.matrix(donnes_RN[,2]), as.matrix(donnes_RN[,1]),  method.bias='linreg') # -- meta

# # The trim-and-fill method is a nonparametric method to assess selection bias/publication bias.
# # The method provides an estimate of (1) the number of missing studies and (2) the treatment 
# # effect adjusted for selection bias. The basic idea of the trim-and-fill method is to add 
# # studies to the funnel plot until it becomes symmetric. 
# tf1 <- trimfill(m1)
# class(tf1)
# tf1 > print(tf1, digits=2, comb.fixed=TRUE) 

# # Copas Selection Model 
# # In contrast to the trim-and-fill method, the selection model by Copas
# # explicitly models publication bias 
# # install.packages("metasens")
# # library(metasens)
# c1 <- copas(m1)
# plot(c1)
# print(summary(c1), digits=2)

# # Adjustment by Regression 
# # a regression-based treatment effect estimate adjusting for small-study effects 
# l1 <- limitmeta(m1)
# print(l1, digits=2)




############################  cumulative meta-analysis  ###############################


# rma using REML can produce NAs for studies in cumul
# is so, then redo the rma using ma_method = 'HS' instead


outp_CUM_META_1 <- cumul(outp_MA_REG)  # cumulative meta-analysis (in the order of publication year)

if (any( apply(data.frame(outp_CUM_META_1), 2, function(x) any(is.na(x))) )) {
		
	outp_MA_REG <- rma(yi=ESdat$yi, vi=ESdat$vi, method='HS')  # random-effects model

	outp_CUM_META_1 <- cumul(outp_MA_REG)  # cumulative meta-analysis (in the order of publication year)
}

outp_CUM_META_2 <- data.frame(cbind(ESdat$yi, outp_CUM_META_1$ci.lb, outp_CUM_META_1$estimate, outp_CUM_META_1$ci.ub))


if (ESdat_type == 'r') {

	if (ES_type_OUT == 'r' & cor_stat == 'COR')	
		results_CUM_META <- cbind(1:Nstudies, totalNs, outp_CUM_META_2)
	
	if (ES_type_OUT == 'r' & cor_stat == 'ZCOR') {
		res <- tanh(outp_CUM_META_2)
		results_CUM_META <- cbind(1:Nstudies, totalNs, res)
	}
	
	if (ES_type_OUT == 'z' & cor_stat == 'COR')	{
		res <- outp_CUM_META_2
		zes <- atanh(res)
		results_CUM_META <- cbind(1:Nstudies, totalNs, zes)
	}
	
	if (ES_type_OUT == 'z' & cor_stat == 'ZCOR')
		results_CUM_META <- cbind(1:Nstudies, totalNs, outp_CUM_META_2)
	
	if (ES_type_OUT == 'd' & cor_stat == 'COR') {
		res <- outp_CUM_META_2
		des <- (2 * res) / sqrt(1 - res^2)
		results_CUM_META <- cbind(1:Nstudies, totalNs, des)
	}
	
	if (ES_type_OUT == 'd' & cor_stat == 'ZCOR') {
		zes <- outp_CUM_META_2
		res <- tanh(zes)
		des <- (2 * res) / sqrt(1 - res^2)
		results_CUM_META <- cbind(1:Nstudies, totalNs, des)
	}
	
	 if (ES_type_OUT == 'g' & cor_stat == 'COR' & all(!is.na(totalNs))) {
	 	res <- outp_CUM_META_2
		des <- (2 * res) / sqrt(1 - res^2)		
		ges <- g_from_d(d = des[,1], totalN = totalNs)[1]
		ges[,2] <- g_from_d(d = des[,2], totalN = totalNs)[1]
		ges[,3] <- g_from_d(d = des[,3], totalN = totalNs)[1]
		ges[,4] <- g_from_d(d = des[,4], totalN = totalNs)[1]
		results_CUM_META <- cbind(1:Nstudies, totalNs, ges)
	}
	
	if (ES_type_OUT == 'g' & cor_stat == 'ZCOR' & all(!is.na(totalNs))) {
		zes <- outp_CUM_META_2
		res <- tanh(zes)		
		des <- (2 * res) / sqrt(1 - res^2)		
		ges <- g_from_d(d = des[,1], totalN = totalNs)[1]
		ges[,2] <- g_from_d(d = des[,2], totalN = totalNs)[1]
		ges[,3] <- g_from_d(d = des[,3], totalN = totalNs)[1]
		ges[,4] <- g_from_d(d = des[,4], totalN = totalNs)[1]
		results_CUM_META <- cbind(1:Nstudies, totalNs, ges)
	}

	colnames(results_CUM_META) <- c('Study','Study N','Study ES','ES_lb','ES','ES_ub')

	if (verbose) {
		message('\n\nCumulative meta-analysis estimates (ES metric = ', ES_type_OUT,'):\n')
		print(round(results_CUM_META,3), print.gap=4, row.names = FALSE)
	}
}


if (ESdat_type == 'z') {  # totalNs are not available, so g not possible

	if (ES_type_OUT == 'z')  results_CUM_META <- outp_CUM_META_2
		
	if (ES_type_OUT == 'r') {
		res <- tanh(outp_CUM_META_2)
		results_CUM_META <- cbind(1:Nstudies, res)
	}
	
	if (ES_type_OUT == 'd') {
		res <- tanh(outp_CUM_META_2)
		des <- (2 * res) / sqrt(1 - res^2)		
		results_CUM_META <- cbind(1:Nstudies, des)
	}
	
	colnames(results_CUM_META) <- c('Study','ES_lb','ES','ES_ub')

	if (verbose) {
		message('\n\nCumulative meta-analysis estimates (ES metric = ', ES_type_OUT,'):\n')
		print(round(results_CUM_META,3), print.gap=4, row.names = FALSE)
	}
} 


if (ESdat_type == 'd') {  # totalNs are not available, so g not possible

	if (ES_type_OUT == 'd')  results_CUM_META <- outp_CUM_META_2
		
	if (ES_type_OUT == 'r')  {
		des <- outp_CUM_META_2		
		res     <- r_from_d(d = des[,1], totalN = totalNs)[1]
		res[,2] <- r_from_d(d = des[,2], totalN = totalNs)[1]
		res[,3] <- r_from_d(d = des[,3], totalN = totalNs)[1]
		res[,4] <- r_from_d(d = des[,4], totalN = totalNs)[1]		
		results_CUM_META <- cbind(1:Nstudies, res)
	}
	
	if (ES_type_OUT == 'z') {
		des <- outp_CUM_META_2		
		res     <- r_from_d(d = des[,1], totalN = totalNs)[1]
		res[,2] <- r_from_d(d = des[,2], totalN = totalNs)[1]
		res[,3] <- r_from_d(d = des[,3], totalN = totalNs)[1]
		res[,4] <- r_from_d(d = des[,4], totalN = totalNs)[1]		
		zes <- atanh(res)
		results_CUM_META <- cbind(1:Nstudies, zes)
	}
	
	colnames(results_CUM_META) <- c('Study','Study ES','ES_lb','ES','ES_ub')

	if (verbose) {
		message('\n\nCumulative meta-analysis estimates (ES metric = ', ES_type_OUT,'):\n')
		print(round(results_CUM_META,3), print.gap=4, row.names = FALSE)
	}
} 


if (ESdat_type != 'r' & ESdat_type != 'z' & ESdat_type != 'd') {
	
	results_CUM_META <- cbind(1:Nstudies, outp_CUM_META_2)
	
	colnames(results_CUM_META) <- c('Study','Study ES','ES_lb','ES','ES_ub')
	
	message('\n\nCumulative meta-analysis estimates (ES metric = ', ES_type_IN,'):\n') 

	print(round(results_CUM_META,3), print.gap=4, row.names = FALSE)
}




#############################  Bayes - create a version of ESdat that is in r metric  ##############################

# necessary because need yi & vi for r for the MA priors

if (cor_stat == 'COR')	 ESdat_r <- ESdat

if (cor_stat == 'ZCOR')	 {
	ESdat_r <- r_from_z(z = ESdat$yi, z_var = ESdat$vi, totalN = donnes_RN$N, zforCI_out=zforCI_out) 
	ESdat_r <- escalc(yi = ESdat_r$r, vi = ESdat_r$r_var)
}



##########################################  Bayes - Schmidt-Raju (2007) method  ####################################

# random effects empirical Bayes meta-analysis 
# 2007 Schmidt, Raju - Updating meta-analytic research findings - Bayesian approaches versus the medical model

if (is.element('Schmidt_Raju', Bayes_type) & !is.null(donnes_RN)) {


results_BAYES_SR <- matrix(NA, Nstudies, 5)

# the 1st row of results
results_BAYES_SR[1,3:5] <- cbind(outp_CUM_META_1$ci.lb[1], outp_CUM_META_1$estimate[1], outp_CUM_META_1$ci.ub[1])


# loop through the studies, treating each subsequent study as the Likelihood
for (luper in 2:Nstudies) {

	priordat_SR <- ESdat_r[1:(luper-1),]



# does the SR method assume that ES = r (and not z)?  If so, then must convert ES, in ESdat, to r if cor_stat = ZCOR
# if not, then must convert z to r after the meta-analysis?

	
	
	outp_MA_SR <- rma(yi=priordat_SR$yi, vi=priordat_SR$vi, method=ma_method, 
	                  control=list(stepadj=0.5, maxiter=1000))
	
	BAYES_SR_priorES <- outp_MA_SR$b
	BAYES_SR_priorSE <- outp_MA_SR$se
	BAYES_SR_priorPopV <- outp_MA_SR$tau2
	BAYES_SR_Vprior <- BAYES_SR_priorSE**2
	k <- nrow(priordat_SR)
	
	rk <- donnes_RN[luper,'ES']
	Nk <- donnes_RN[luper,'N']
	
	Vrk <- ( (1 - rk**2)**2 / (Nk - 1) )  + BAYES_SR_priorPopV
	
	BAYES_SR_r <- ( (BAYES_SR_priorES / BAYES_SR_Vprior) + (rk / Vrk) ) / 
	              ( (1 / BAYES_SR_Vprior) + (1 / Vrk))  # Schmidt & Raju, 2007, p 302, formula 6
	
	BAYES_SR_postV <- (BAYES_SR_Vprior * Vrk) / (BAYES_SR_Vprior + Vrk)  # Schmidt & Raju, 2007, p 302, formula 7
	
	BAYES_SR_postPopV <- ( (k - 1) * BAYES_SR_priorPopV + (k**2 * BAYES_SR_postV - (k - 1)**2 * 
	                      BAYES_SR_Vprior) - Vrk) / k # Schmidt & Raju, 2007, p 303, formula 8 -- same as cc$tau2
	
	# credibility intervals -- see Field 2005 p 448, formula 16
	# Hunter and Schmidt recommend correcting this estimate for artifacts 
	# (see Hunter & Schmidt, 2004, or Hall & Brannick, 2002 for details) and then
	# constructing what they call credibility intervals. These intervals are based
	# on taking the average correlation (see Equation 12) and adding to or subtracting
	# from it the square root of the estimated population variance in 
	# Equation 15 multiplied by ... (1.96 for a 95% interval) 
	
	BAYES_SR_rub <- BAYES_SR_r + zforCI_out * sqrt(BAYES_SR_postV)  # sqrt(BAYES_SR_postV) = cc$se
	BAYES_SR_rlb <- BAYES_SR_r - zforCI_out * sqrt(BAYES_SR_postV)
	
	# confidence intervals -- see Field 2005 p 448, formula 17
	# If confidence intervals are required (rather than credibility intervals) these can
	# be obtained by using the standard error of the mean correlation. To obtain this
	# standard error simply divide the variance of sample correlations (given in Equation 13)
	# by the number of studies in the meta-analysis, k, and take the square root: 
	# BAYES_SR_rub <- BAYES_SR_r + zforCI_out * sqrt(BAYES_SR_postV / k)
	# BAYES_SR_rlb <- BAYES_SR_r - zforCI_out * sqrt(BAYES_SR_postV / k)
	
	results_BAYES_SR[luper,] <- cbind(BAYES_SR_priorES, BAYES_SR_priorSE, # BAYES_SR_postV, BAYES_SR_postPopV, 
	                                  BAYES_SR_rlb, BAYES_SR_r, BAYES_SR_rub)

	postN_Vs <- cbind(donnes_RN[luper,1], BAYES_SR_postV, BAYES_SR_postPopV)
}

results_BAYES_SR <- data.frame(cbind( 1:Nstudies, as.matrix(donnes_RN), results_BAYES_SR ))

colnames(results_BAYES_SR) <- c('Study','Study N','Study ES','prior ES','prior SE','ES_lb','ES','ES_ub')  


# convert r to z, d or g, if requested
if (ES_type_OUT == 'z' | ES_type_OUT == 'd' | ES_type_OUT == 'g') {

	# leaving out prior SE
	results_BAYES_SR_zdg <- results_BAYES_SR[,c('Study','Study N','Study ES','prior ES','ES_lb','ES','ES_ub')]

	if (ES_type_OUT == 'z')	{
		results_BAYES_SR_zdg[3:7] <- atanh(results_BAYES_SR_zdg[3:7])
	}
		
	if (ES_type_OUT == 'd') {
		des <- (2 * results_BAYES_SR_zdg[3:7]) / sqrt(1 - results_BAYES_SR_zdg[3:7]^2)				
		results_BAYES_SR_zdg[3:7] <- des
	}
	
	 if (ES_type_OUT == 'g') {
		des <- (2 * results_BAYES_SR_zdg[3:7]) / sqrt(1 - results_BAYES_SR_zdg[3:7]^2)				
		ges <- g_from_d(d = des[,1], totalN = totalNs)[1]
		ges[,2] <- g_from_d(d = des[,2], totalN = totalNs)[1]
		ges[,3] <- g_from_d(d = des[,3], totalN = totalNs)[1]
		ges[,4] <- g_from_d(d = des[,4], totalN = totalNs)[1]
		results_BAYES_SR_zdg[3:7] <- ges
	}

	results_BAYES_SR <- results_BAYES_SR_zdg  
}


if (verbose) {
	message('\n\n\nBayesian estimates based on the Schmidt & Raju (2007) method (ES metric = ', ES_type_OUT,'):\n')
	print(round(results_BAYES_SR,3), print.gap=4, row.names = FALSE)
}

}



##########################################  Bayes generated data  #############################################

 
if (is.element('generated', Bayes_type) & !is.null(donnes_RN)) {

# generate likelihood raw data for 2 variables with a correlation = the current effect size

# run MCMCglmm, using the effect size & sampling error variance from a MA of previous data

# loop through the studies, treating each subsequent study as the Likelihood
results_BAYES_GEN   <- matrix(NA,Nstudies,5)
BAYES_GEN_post_ests <- matrix(NA,Nstudies,2)

for (luper in 1:Nstudies) {
	
	correl <- donnes_RN[luper,'ES']
	N      <- round(donnes_RN[luper,'N'])
	
	# generating variables with an exact correlation
	dataset1 <- data.frame(mvrnorm(n=N,mu=c(0,0),
						   Sigma=matrix(c(1,correl,correl, 1),nrow=2),empirical=TRUE))
	colnames(dataset1) <- c('varIV','varDV')
	
	# # saving the generated data for the Bayes Factor analyses (when there is no raw data)
	# if (dontype == 'matrix') donnes[[luper]] <- dataset1
	
	if (luper == 1) { 
		model1 <- MCMCglmm(varDV ~ varIV, data=dataset1, nitt=nitt, burnin=burnin, thin=thin, verbose=FALSE)
		model1sum <- summary.MCMCglmm(model1)
		BAYES_GEN_r   <- model1sum$solutions[2,1]
		BAYES_GEN_rlb <- model1sum$solutions[2,2]
		BAYES_GEN_rub <- model1sum$solutions[2,3]
		BAYES_GEN_priorES <- NA
		BAYES_GEN_priorSE <- NA
		BAYES_GEN_post_ests[1,] <- c(BAYES_GEN_r, (diag(var(model1$Sol)))[2])  # saving the ES & the posterior variance
	}	
		  
	if (luper > 1) { 
		if (prior_type == 'META') { 		
			priordat <- ESdat_r[1:(luper-1),]
			outp_MA_GEN <- rma(yi=priordat$yi, vi=priordat$vi, method=ma_method, 
	                           control=list(stepadj=0.5, maxiter=1000))  # random-effects model

			BAYES_GEN_priorES <- outp_MA_GEN$b
			BAYES_GEN_priorSE <- outp_MA_GEN$se
#			BAYES_GEN_priorPopV <- outp_MA_GEN$tau2
			BAYES_GEN_Vprior <- BAYES_GEN_priorSE**2			
		}
		if (prior_type == 'BAYES') {	
			# Jarrod Hadfield  https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q4/019436.html
			# v <- var(model$Sol)  the variance-covariance matrix
			# sqrt(diag(v))  the posterior standard deviations (akin to the standard errors) 
			BAYES_GEN_priorES <- BAYES_GEN_post_ests[(luper-1),1]
			BAYES_GEN_Vprior  <- BAYES_GEN_post_ests[(luper-1),2]
			BAYES_GEN_priorSE <- sqrt(BAYES_GEN_Vprior)
		}
		priors <- list(B=list(mu=c(0, BAYES_GEN_priorES), V=diag(c(1, BAYES_GEN_Vprior))))		
		model2 <- MCMCglmm(varDV ~ varIV, data=dataset1, prior=priors, nitt=nitt, burnin=burnin, 
		                   thin=thin, verbose=FALSE)
		model2sum  <- summary.MCMCglmm(model2)
		BAYES_GEN_r   <- model2sum$solutions[2,1]
		BAYES_GEN_rlb <- model2sum$solutions[2,2]
		BAYES_GEN_rub <- model2sum$solutions[2,3]
		BAYES_GEN_post_ests[luper,] <- c(BAYES_GEN_r, (diag(var(model2$Sol)))[2])  # the ES & the posterior variance
	}

results_BAYES_GEN[luper,] <- cbind(BAYES_GEN_priorES, BAYES_GEN_priorSE, BAYES_GEN_rlb, BAYES_GEN_r, BAYES_GEN_rub)

}

results_BAYES_GEN <- data.frame(cbind( 1:Nstudies, as.matrix(donnes_RN), results_BAYES_GEN ))

colnames(results_BAYES_GEN) <- 
    c('Study','Study N','Study ES','prior ES','prior SE','ES_lb','ES','ES_ub')


# convert r to z, d or g, if requested
if (ES_type_OUT == 'z' | ES_type_OUT == 'd' | ES_type_OUT == 'g') {

	# leaving out prior SE
	results_BAYES_GEN_zdg <- results_BAYES_GEN[,c('Study','Study N','Study ES','prior ES','ES_lb','ES','ES_ub')]

	if (ES_type_OUT == 'z')	{
		results_BAYES_GEN_zdg[3:7] <- atanh(results_BAYES_GEN_zdg[3:7])
	}
		
	if (ES_type_OUT == 'd') {
		des <- (2 * results_BAYES_GEN_zdg[3:7]) / sqrt(1 - results_BAYES_GEN_zdg[3:7]^2)				
		results_BAYES_GEN_zdg[3:7] <- des
	}
	
	 if (ES_type_OUT == 'g') {
		des <- (2 * results_BAYES_GEN_zdg[3:7]) / sqrt(1 - results_BAYES_GEN_zdg[3:7]^2)				
		ges <- g_from_d(d = des[,1], totalN = totalNs)[1]
		ges[,2] <- g_from_d(d = des[,2], totalN = totalNs)[1]
		ges[,3] <- g_from_d(d = des[,3], totalN = totalNs)[1]
		ges[,4] <- g_from_d(d = des[,4], totalN = totalNs)[1]
		results_BAYES_GEN_zdg[3:7] <- ges
	}

	results_BAYES_GEN <- results_BAYES_GEN_zdg  
}


if (verbose) {
	message('\n\n\nBayesian estimates based on generated data (ES metric = ', ES_type_OUT,'):\n')
	print(round(results_BAYES_GEN,3), print.gap=4, row.names = FALSE)
}

}




######################################  Bayes raw data  ########################################


if (is.element('raw', Bayes_type) & dontype != 'raw data') {
	message('\n\nBayesian raw data analyses were requested but donnes is not a list with raw data.')
	message('Bayesian raw data analyses cannot be conducted without raw data.\n')
}

results_BAYES_RAW <- -9999

if (is.element('raw', Bayes_type) & dontype == 'raw data' & 
               (rawdata_type == 'for_correl' | (rawdata_type == 'indep_groups'))) {

BAYES_RAW_rlb <- BAYES_RAW_r <- BAYES_RAW_rub <- NA

# run MCMCglmm, using the effect size & sampling error variance from a MA of previous data

# loop through the studies, treating each subsequent study as the Likelihood
results_BAYES_RAW   <- matrix(NA,Nstudies,5)
BAYES_RAW_post_ests <- matrix(NA,Nstudies,2)

for (luper in 1:Nstudies) {

	dataset1 <- data.frame(donnes[[luper]]); colnames(dataset1) <- c('varIV','varDV')
	
	# standardize data so that the regression estimate is a correlation
	dataset1 <- data.frame(scale(dataset1))
		
	# converting the raw coeff & CI to standardized using conventional formula comes close, not close enough
	# sdDV <- sd(dataset1$varDV)
	# sdIV <- sd(dataset1$varIV)
	# beta <- model1sum$solutions[2,1] * sdIV / sdDV

	# rlb <- model1sum$solutions[2,2] * sdIV / sdDV
	# rub <- model1sum$solutions[2,3] * sdIV / sdDV

	if (luper == 1) { 
		model1 <- MCMCglmm(varDV ~ varIV, data=dataset1, nitt=nitt, burnin=burnin, thin=thin, verbose=FALSE)
		model1sum <- summary.MCMCglmm(model1)
		BAYES_RAW_r   <- model1sum$solutions[2,1]
		BAYES_RAW_rlb <- model1sum$solutions[2,2]
		BAYES_RAW_rub <- model1sum$solutions[2,3]
		BAYES_RAW_priorES <- NA
		BAYES_RAW_priorSE <- NA
		BAYES_RAW_post_ests[1,] <- c(BAYES_RAW_r, (diag(var(model1$Sol)))[2])  # the ES & posterior variance
	}	
		  
	if (luper > 1) { 		
		if (prior_type == 'META') {	
			priordat <- ESdat_r[1:(luper-1),]
			outp_MA_RAW <- rma(yi=priordat$yi, vi=priordat$vi, method=ma_method, 
	                         control=list(stepadj=0.5, maxiter=1000))  # random-effects model

			BAYES_RAW_priorES <- outp_MA_RAW$b
			BAYES_RAW_priorSE <- outp_MA_RAW$se
	#		BAYES_RAW_priorPopV <- outp_MA_RAW$tau2						
			BAYES_RAW_Vprior <- BAYES_RAW_priorSE**2	
		}
		if (prior_type == 'BAYES') {	
			# Jarrod Hadfield  https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q4/019436.html
			# v <- var(model$Sol)  the variance-covariance matrix
			# sqrt(diag(v))  the posterior standard deviations (akin to the standard errors) 
			BAYES_RAW_priorES <- BAYES_RAW_post_ests[(luper-1),1]
			BAYES_RAW_Vprior  <- BAYES_RAW_post_ests[(luper-1),2]
			BAYES_RAW_priorSE <- sqrt(BAYES_RAW_Vprior)
		}
		priors <- list(B=list(mu=c(0, BAYES_RAW_priorES), V=diag(c(1, BAYES_RAW_Vprior))))		
		model2 <- MCMCglmm(varDV ~ varIV, data=dataset1, prior=priors, nitt=nitt, burnin=burnin, 
		                   thin=thin, verbose=FALSE)
		model2sum <- summary.MCMCglmm(model2)
		BAYES_RAW_r   <- model2sum$solutions[2,1]
		BAYES_RAW_rlb <- model2sum$solutions[2,2]
		BAYES_RAW_rub <- model2sum$solutions[2,3]
		BAYES_RAW_post_ests[luper,] <- c(BAYES_RAW_r, (diag(var(model2$Sol)))[2]) # ES & the posterior variance
	}	
	results_BAYES_RAW[luper,] <- cbind(BAYES_RAW_priorES, BAYES_RAW_priorSE, BAYES_RAW_rlb, BAYES_RAW_r, BAYES_RAW_rub)
}

results_BAYES_RAW <- data.frame(cbind(1:Nstudies, as.matrix(donnes_RN), results_BAYES_RAW))

colnames(results_BAYES_RAW) <- 
  c('Study','Study N','Study ES','prior ES','prior SE','ES_lb','ES','ES_ub')


# convert r to z, d or g, if requested
if (ES_type_OUT == 'z' | ES_type_OUT == 'd' | ES_type_OUT == 'g') {

	# leaving out prior SE
	results_BAYES_RAW_zdg <- results_BAYES_RAW[,c('Study','Study N','Study ES','prior ES','ES_lb','ES','ES_ub')]

	if (ES_type_OUT == 'z')	{
		results_BAYES_RAW_zdg[3:7] <- atanh(results_BAYES_RAW_zdg[3:7])
	}
		
	if (ES_type_OUT == 'd') {
		des <- (2 * results_BAYES_RAW_zdg[3:7]) / sqrt(1 - results_BAYES_RAW_zdg[3:7]^2)				
		results_BAYES_RAW_zdg[3:7] <- des
	}
	
	 if (ES_type_OUT == 'g') {
		des <- (2 * results_BAYES_RAW_zdg[3:7]) / sqrt(1 - results_BAYES_RAW_zdg[3:7]^2)				
		ges <- g_from_d(d = des[,1], totalN = totalNs)[1]
		ges[,2] <- g_from_d(d = des[,2], totalN = totalNs)[1]
		ges[,3] <- g_from_d(d = des[,3], totalN = totalNs)[1]
		ges[,4] <- g_from_d(d = des[,4], totalN = totalNs)[1]
		results_BAYES_RAW_zdg[3:7] <- ges
	}

	results_BAYES_RAW <- results_BAYES_RAW_zdg  
}


if (verbose) {
	message('\n\n\nBayesian estimates based on the raw data:\n')
	print(round(results_BAYES_RAW,3), print.gap=4, row.names = FALSE)
}

}




###########################  agreement & consistency  ##########################


# the CUM_META MA conclusion
if (results_CUM_META[,'ES_lb'][Nstudies] < 0 & results_CUM_META[,'ES_ub'][Nstudies] < 0)  effect <- 'negeff'
if (results_CUM_META[,'ES_lb'][Nstudies] > 0 & results_CUM_META[,'ES_ub'][Nstudies] > 0)  effect <- 'poseff'
if (results_CUM_META[,'ES_lb'][Nstudies] < 0 & results_CUM_META[,'ES_ub'][Nstudies] > 0)  effect <- 'noeff'

# NHST -- agreement with final, & consistency
sigposNHST <- signegNHST <- nonsigNHST <- 0
for (luper in 1:nrow(results_NHST)) {
	if (results_NHST[luper,'ES_lb'] < 0 & results_NHST[luper,'ES_ub'] < 0) signegNHST <- signegNHST + 1
	if (results_NHST[luper,'ES_lb'] > 0 & results_NHST[luper,'ES_ub'] > 0) sigposNHST <- sigposNHST + 1
	if (results_NHST[luper,'ES_lb'] < 0 & results_NHST[luper,'ES_ub'] > 0) nonsigNHST <- nonsigNHST + 1
}

# percentage of NHST studies that agreed with the CUM_META conclusion
if (effect == 'negeff') agree_NHST <- signegNHST / nrow(results_NHST)
if (effect == 'poseff') agree_NHST <- sigposNHST / nrow(results_NHST)
if (effect == 'noeff')  agree_NHST <- nonsigNHST / nrow(results_NHST)

# the most common conclusion, as a proportion
consist_NHST <- max(c(sigposNHST,signegNHST,nonsigNHST)) / nrow(results_NHST) 



# CUM_META  -- agreement with final, & consistency
sigposCUM_META <- signegCUM_META <- nonsigCUM_META <- 0
for (luper in 1:nrow(results_CUM_META)) {
	if (results_CUM_META[luper,'ES_lb'] < 0 & results_CUM_META[luper,'ES_ub'] < 0) signegCUM_META <- signegCUM_META + 1
	if (results_CUM_META[luper,'ES_lb'] > 0 & results_CUM_META[luper,'ES_ub'] > 0) sigposCUM_META <- sigposCUM_META + 1
	if (results_CUM_META[luper,'ES_lb'] < 0 & results_CUM_META[luper,'ES_ub'] > 0) nonsigCUM_META <- nonsigCUM_META + 1
}

# percentage of CUM_META studies that agreed with the CUM_META conclusion
if (effect == 'negeff') agree_CUM_META <- signegCUM_META / nrow(results_CUM_META)
if (effect == 'poseff') agree_CUM_META <- sigposCUM_META / nrow(results_CUM_META)
if (effect == 'noeff')  agree_CUM_META <- nonsigCUM_META / nrow(results_CUM_META)

# the most common conclusion, as a proportion
consist_CUM_META <- max(c(sigposCUM_META,signegCUM_META,nonsigCUM_META)) / nrow(results_CUM_META)



if (!is.null(donnes_RN)) {

	if (is.element('Schmidt_Raju', Bayes_type)) {
		agree_BAYES_SR <- consist_BAYES_SR <- NA
		# BAYES_SR -- agreement with final, & consistency
		sigposBAYES_SR <- signegBAYES_SR <- nonsigBAYES_SR <- 0
		for (luper in 1:nrow(results_BAYES_SR)) {
			if (results_BAYES_SR[luper,'ES_lb'] < 0 & results_BAYES_SR[luper,'ES_ub'] < 0) signegBAYES_SR <- signegBAYES_SR + 1
			if (results_BAYES_SR[luper,'ES_lb'] > 0 & results_BAYES_SR[luper,'ES_ub'] > 0) sigposBAYES_SR <- sigposBAYES_SR + 1
			if (results_BAYES_SR[luper,'ES_lb'] < 0 & results_BAYES_SR[luper,'ES_ub'] > 0) nonsigBAYES_SR <- nonsigBAYES_SR + 1
		}
		# percentage of BAYES_SR studies that agreed with the CUM_META conclusion
		if (effect == 'negeff') agree_BAYES_SR <- signegBAYES_SR / nrow(results_BAYES_SR)
		if (effect == 'poseff') agree_BAYES_SR <- sigposBAYES_SR / nrow(results_BAYES_SR)
		if (effect == 'noeff')  agree_BAYES_SR <- nonsigBAYES_SR / nrow(results_BAYES_SR)
		
		# the most common conclusion, as a proportion
		consist_BAYES_SR <- max(c(sigposBAYES_SR,signegBAYES_SR,nonsigBAYES_SR)) / nrow(results_BAYES_SR) 
	}
	
	
	if (is.element('generated', Bayes_type)) {
		# generated BAYES -- agreement with final, & consistency
		sigposBAYES_GEN <- signegBAYES_GEN <- nonsigBAYES_GEN <- 0
		for (luper in 1:nrow(results_BAYES_GEN)) {
			if (results_BAYES_GEN[luper,'ES_lb'] < 0 & results_BAYES_GEN[luper,'ES_ub'] < 0) signegBAYES_GEN <- signegBAYES_GEN + 1
			if (results_BAYES_GEN[luper,'ES_lb'] > 0 & results_BAYES_GEN[luper,'ES_ub'] > 0) sigposBAYES_GEN <- sigposBAYES_GEN + 1
			if (results_BAYES_GEN[luper,'ES_lb'] < 0 & results_BAYES_GEN[luper,'ES_ub'] > 0) nonsigBAYES_GEN <- nonsigBAYES_GEN + 1
		}
		# percentage of BAYES_GEN studies that agreed with the CUM_META conclusion
		if (effect == 'negeff') agree_BAYES_GEN <- signegBAYES_GEN / nrow(results_BAYES_GEN)
		if (effect == 'poseff') agree_BAYES_GEN <- sigposBAYES_GEN / nrow(results_BAYES_GEN)
		if (effect == 'noeff')  agree_BAYES_GEN <- nonsigBAYES_GEN / nrow(results_BAYES_GEN)
		
		# the most common conclusion, as a proportion
		consist_BAYES_GEN <- max(c(sigposBAYES_GEN,signegBAYES_GEN,nonsigBAYES_GEN)) / nrow(results_BAYES_GEN) 
	}
	
	
	#if (is.element('raw', Bayes_type) & dontype == 'raw data') {
	if (is.element('raw', Bayes_type) & results_BAYES_RAW != -9999) {
		agree_BAYES_RAW <- consist_BAYES_RAW <- NA
		# raw BAYES -- agreement with final, & consistency
		sigposBAYES_RAW <- signegBAYES_RAW <- nonsigBAYES_RAW <- 0
		for (luper in 1:nrow(results_BAYES_RAW)) {
			if (results_BAYES_RAW[luper,'ES_lb'] < 0 & results_BAYES_RAW[luper,'ES_ub'] < 0) signegBAYES_RAW <- signegBAYES_RAW + 1
			if (results_BAYES_RAW[luper,'ES_lb'] > 0 & results_BAYES_RAW[luper,'ES_ub'] > 0) sigposBAYES_RAW <- sigposBAYES_RAW + 1
			if (results_BAYES_RAW[luper,'ES_lb'] < 0 & results_BAYES_RAW[luper,'ES_ub'] > 0) nonsigBAYES_RAW <- nonsigBAYES_RAW + 1
		}
		# percentage of BAYES_RAW studies that agreed with the CUM_META conclusion
		if (effect == 'negeff') agree_BAYES_RAW <- signegBAYES_RAW / nrow(results_BAYES_RAW)
		if (effect == 'poseff') agree_BAYES_RAW <- sigposBAYES_RAW / nrow(results_BAYES_RAW)
		if (effect == 'noeff')  agree_BAYES_RAW <- nonsigBAYES_RAW / nrow(results_BAYES_RAW)
		
		# the most common conclusion, as a proportion
		consist_BAYES_RAW <- max(c(sigposBAYES_RAW,signegBAYES_RAW,nonsigBAYES_RAW)) / nrow(results_BAYES_RAW) 
	}

}


####################################  output  #######################################

if (verbose) {

	message('\n\nFinal, all-studies-combined results (ES metric = ', ES_type_OUT,'):')
	
	message('\n   Cumulative Meta-Analysis:    ES = ', round(results_CUM_META[,'ES'][Nstudies],3),
	        '   ES_lb = ', round(results_CUM_META[,'ES_lb'][Nstudies],3),
	        '   ES_ub = ', round(results_CUM_META[,'ES_ub'][Nstudies],3))
	
if (!is.null(donnes_RN)) {

	if (is.element('Schmidt_Raju', Bayes_type)) 
		message('\n   Bayesian (Schmidt-Raju):     ES = ', round(results_BAYES_SR$ES[Nstudies],3),
		        '   ES_lb = ', round(results_BAYES_SR$ES_lb[Nstudies],3),
		        '   ES_ub = ', round(results_BAYES_SR$ES_ub[Nstudies],3))
	
	if (is.element('generated', Bayes_type)) 
		message('\n   Bayesian (generated data):   ES = ', round(results_BAYES_GEN$ES[Nstudies],3),
		        '   ES_lb = ', round(results_BAYES_GEN$ES_lb[Nstudies],3),
		        '   ES_ub = ', round(results_BAYES_GEN$ES_ub[Nstudies],3))
	
	if (is.element('raw', Bayes_type) & results_BAYES_RAW != -9999) 
		message('\n   Bayesian (raw data):         ES = ', round(results_BAYES_RAW$ES[Nstudies],3),
		        '   ES_lb = ', round(results_BAYES_RAW$ES_lb[Nstudies],3),
		        '   ES_ub = ', round(results_BAYES_RAW$ES_ub[Nstudies],3))

}	
	
	if (ma_method != 'FE') {
		message('\n\n\nMeta-analysis heterogeneity statistics:')
		message('\n   Total Var. = ',round(tau2,3))
		message('\n   Q = ', round(outp_MA_REG$QE,3), '    p = ', round(outp_MA_REG$QEp,5))
		message('\n   tau2 = ',  round(tau2,3), '    tau2.lb = ', round(tau2LB,3), '    tau2.ub = ', round(tau2UB,3))
		message('\n   tau  = ',   round(tau,3), '    tau.lb = ', round(tauLB,3),  '    tau.ub = ', round(tauUB,3))
		message('\n   Isq. = ',  round(isq,3),  '    Isq.lb = ', round(isqLB,3),  '    Isq.ub = ', round(isqUB,3))
		message('\n   Hsq. = ',  round(hsq,3),  '    Hsq.lb = ', round(hsqLB,3),  '    Hsq.ub = ', round(hsqUB,3))
	
		message('\n\ntau2 (or tau-squared) is the variation in effect sizes')
		message('(between-study variance) in a random-effects meta-analysis.')
		message('It is the variance in the true effect sizes.')
		
		message('\ntau is the square root of tau-squared. tau is the standard')
		message('deviation of the true effect sizes.')
		
		message('\nIsq. estimates (in percent) how much of the total variability')
		message('in the effect size estimates (which is composed of heterogeneity plus')
		message('sampling variability) can be attributed to heterogeneity among the true effects.')
		
		message('\nHsq. estimates the ratio of the total amount of variability in')
		message('the effect size estimates to the amount of sampling variability.\n')
	}	
	
	
	message('\n\nConsistency & agreement rates:')
	
	message('\n   NHST:                               consistency = ', 
		round(consist_NHST,3), '    agreement = ', round(agree_NHST,3))
	
	message('\n   Cumulative Meta-Analysis:           consistency = ', 
		round(consist_CUM_META,3),  '    agreement = ', round(agree_CUM_META,3))
	
if (!is.null(donnes_RN)) {

	if (is.element('Schmidt_Raju', Bayes_type)) 
		message('\n   Bayesian (Schmidt & Raju, 2007):    consistency = ', 
			round(consist_BAYES_SR,3),   '    agreement = ', round(agree_BAYES_SR,3))
	
	if (is.element('generated', Bayes_type)) 
		message('\n   Bayesian (generated data):          consistency = ', 
			round(consist_BAYES_GEN,3),   '    agreement = ', round(agree_BAYES_GEN,3))
	
	if (is.element('raw', Bayes_type) & results_BAYES_RAW != -9999) 
		message('\n   Bayesian (raw data):                consistency = ', 
			round(consist_BAYES_RAW,3),   '    agreement = ', round(agree_BAYES_RAW,3))

}
	
	message('\n\nThe consistency rate is the proportion of times that the most common')
	message('conclusion was reached for an analytic method for a pool of effect sizes.')
	message('Three conclusions are possible for each effect size: a positive effect, a')
	message('negative effect, and no effect. The signs of the effect sizes and the')
	message('possible inclusion of a zero value in a confidence interval were used to')
	message('make these categorizations (e.g., a "negative effect" conclusion was when a')
	message('negative effect size had a confidence interval that did not include zero).')
	
	message('\nThe number of times each of the three possible conclusions occurred for a')
	message('pool of effect sizes was counted, and the consistency rate was based on the')
	message('most common conclusion.')
	
	message('\nThe agreement rate for a pool of effect sizes is the proportion of times')
	message('that the conclusions for individual studies were identical to the')
	message('conclusion (re: the same three categories) of the final,')
	message('all-studies-combined meta-analysis.')
	
	message('\n')

}

###########################################################################################


nopingpongOutput <- list(
	results_NHST = results_NHST, 	
	results_CUM_META = results_CUM_META, 
	ES_MA    = results_CUM_META[,'ES'][Nstudies],  
	ES_MA_lb = results_CUM_META[,'ES_lb'][Nstudies],  
	ES_MA_ub = results_CUM_META[,'ES_ub'][Nstudies],
	Q = outp_MA_REG$QE, p_Q = outp_MA_REG$QEp,
	tau2 = tau2, tau2LB = tau2LB, tau2UB = tau2UB,
	tau = tau, tauLB = tauLB, tauUB = tauUB,
	isq = isq, isqLB = isqLB, isqUB = isqUB,
	hsq = hsq, hsqLB = hsqLB, hsqUB = hsqUB,
	consist_NHST = consist_NHST, agree_NHST = agree_NHST,
	consist_CUM_META = consist_CUM_META, agree_CUM_META = agree_CUM_META,
	biasStats = biasStats,
	results_BAYES_SR  = NULL, consist_BAYES_SR  = NULL, agree_BAYES_SR  = NULL,
	results_BAYES_GEN = NULL, consist_BAYES_GEN = NULL, agree_BAYES_GEN = NULL,
	results_BAYES_RAW = NULL, consist_BAYES_RAW = NULL, agree_BAYES_RAW = NULL	
 )


if (!is.null(donnes_RN)) {

	if (is.element('Schmidt_Raju', Bayes_type)) {
	    nopingpongOutput$results_BAYES_SR = results_BAYES_SR 
		nopingpongOutput$consist_BAYES_SR = consist_BAYES_SR 
		nopingpongOutput$agree_BAYES_SR   = agree_BAYES_SR 
	}
	
	if (is.element('generated', Bayes_type)) { 
	    nopingpongOutput$results_BAYES_GEN = results_BAYES_GEN
		nopingpongOutput$consist_BAYES_GEN = consist_BAYES_GEN 
		nopingpongOutput$agree_BAYES_GEN   = agree_BAYES_GEN 
	}
			
	if (is.element('raw', Bayes_type) & results_BAYES_RAW != -9999) {
	    nopingpongOutput$results_BAYES_RAW = results_BAYES_RAW
		nopingpongOutput$consist_BAYES_RAW = consist_BAYES_RAW
		nopingpongOutput$agree_BAYES_RAW   = agree_BAYES_RAW 
	}

}

class(nopingpongOutput) <- "NO.PING.PONG"

return(invisible(nopingpongOutput))

}




