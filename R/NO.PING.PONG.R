 

NO.PING.PONG <- function(donnes, ES_type_IN=NULL, ES_type_OUT='r', 
                         rma_method='REML',
                         Bayes_type = c('Schmidt_Raju', 'generated', 'raw'), 
                         prior_type='MA', CI = 95,
                         ES = NULL, N = NULL, vi = NULL,
                         grp1_mn = NULL, grp1_sd = NULL, grp1_n = NULL, 
                         grp2_mn = NULL, grp2_sd = NULL, grp2_n = NULL, 
                         gvar_type_OUT = 'd',
                         nitt=53000, burnin=3000, thin=10, verbose=TRUE) {


# get the object name of donnes, to be used in the name of the saved plot file
titre <- deparse(substitute(donnes))


# determine if donnes is a list (if it is raw data)
dtype <- ifelse (inherits(donnes, "list") == TRUE, 'list', 'not list') 

#if (dtype == 'not list') donnes_RN <- donnes
if (dtype == 'list')     donnesRAW <- donnes


zforCI <- qnorm((1 + CI * .01) / 2) # the z value that corresponds to the specified CI


# is complete info (M, SD, & N) for 2 groups provided?
if (!is.null(grp1_mn) & !is.null(grp1_sd) & !is.null(grp1_n) &
    !is.null(grp2_mn) & !is.null(grp2_sd) & !is.null(grp2_n)  ) 
    {grpinfo = TRUE} else {grpinfo = FALSE}

# notice if both totalN & grp Ns are missing & dtype is not a list
if (dtype == 'not list' & grpinfo == FALSE & is.null(N)) {
	message('\nBoth "N", and "grp1_n" & "grp2_n", are missing/NULL. The analyses cannot be')
	message('conducted. Provide either "N", or both "grp1_n" & "grp2_n".\n')
}

if (grpinfo == TRUE & is.null(ES_type_IN))  ES_type_IN = 'g'
	
if (grpinfo == FALSE & is.null(ES_type_IN)) ES_type_IN = 'r'



######################  compute ESdat for meta-analyses  &  donnesRN for Bayesian  #####################

# ESdat has the yi & vi values for meta-analyses; yi (ES) is in correlation coef metric

# donnesRN has the effect size (r) & N for the Bayesian analyses

donnes_RN <- NULL # won't do Bayesian if this remains null

#totalNs <- NULL


if (ES_type_IN == 'r') {	
	if (dtype == 'not list') {

		# ESdat if ES & N are provided
		# Vr may be available, but it is easier for ESdat to be an escalc object (e.g., for summESdat)
		if (!is.null(ES) & !is.null(N) ) {
			donnes_RN <- data.frame(N = donnes[,N], ES = donnes[,ES])
			ESdat <- escalc(measure='COR', ni=donnes_RN$N, ri=donnes_RN$ES )
			# totalNs <- attr(ESdat$yi, 'ni')
		}
	}

	if (dtype == 'list') {
		donnes_RN <- data.frame(N = rep(NA,length(donnes)), ES = rep(NA,length(donnes)))
		for (lupe in 1:length(donnes)) {
			donnes_RN$N[lupe]  <- nrow(donnes[[lupe]])
			donnes_RN$ES[lupe] <- cor(na.omit(donnes[[lupe]]))[2,1]
		}
		ESdat <- escalc(measure='COR', ni=donnes_RN$N, ri=donnes_RN$ES )	
#		totalNs <- attr(ESdat$yi, 'ni')
	}
}
	

if (ES_type_IN == 'd' | ES_type_IN == 'g') { 
	
	if (dtype == 'not list') { # need ES & vi, & also need N1 & N2, or totalN, for Bayesian

		# if ES & vi & totalN are provided, but group info is NOT provided
		if (!is.null(ES) & !is.null(vi) & !is.null(N) & grpinfo == FALSE) { 

			rfromgd <- CONVERT_ES(ES= donnes[,ES], ES_var = donnes[,vi], 
			                     ES_type_IN=ES_type_IN, ES_type_OUT='r', totalN = donnes[,N], verbose = FALSE)
			
			donnes_RN <- data.frame(N = donnes[,N], ES = rfromgd$r) 

			ESdat <- escalc(measure='COR', ni=donnes_RN$N, ri=donnes_RN$ES )	
		}		
			
		# if vi is not provided, but grp1 & grp2 info is available
		if (is.null(vi) & grpinfo) {  

			ESdat_grpinfo <- ES_from_GRPINFO(
			                    grp1_MN=donnes[,grp1_mn], grp1_SD=donnes[,grp1_sd], grp1_N=donnes[,grp1_n], 
				                grp2_MN=donnes[,grp2_mn], grp2_SD=donnes[,grp2_sd], grp2_N=donnes[,grp2_n],
				                ES_type_OUT = 'r')

			# from the escalc documentation: The positive bias in the standardized mean difference is 
			# automatically corrected for within the function, yielding Hedges' g for measure="SMD" (Hedges, 1981)	
			# I am therefore using ES_from_GRPINFO, because donnes can sometimes be d and not g values	

			# ESdat_g <- escalc(measure='SMD', vtype='UB', 
			                  # m1i=donnes[,grp1_mn], sd1i=donnes[,grp1_sd], n1i=donnes[,grp1_n], 
				              # m2i=donnes[,grp2_mn], sd2i=donnes[,grp2_sd], n2i=donnes[,grp2_n])

			# totalNs <- attr(ESdat_grpinfo$yi, 'ni')

			# rfromgd <- CONVERT_ES(ES = ESdat_grpinfo$yi, Vd = ESdat_grpinfo$vi, dg_type = ES_type_OUT, 
			                     # grp1_N = donnes[,grp1_n], grp2_N = donnes[,grp2_n], verbose = FALSE)
			                     			
			# donnes_RN <- data.frame(N = rfromgd$totalN, ES = rfromgd$r) 

			donnes_RN <- data.frame(N = ESdat_grpinfo$totalN, ES = ESdat_grpinfo$r) 

			ESdat <- escalc(measure='COR', ni=donnes_RN$N, ri=donnes_RN$ES )
		}

		# when vi & group info are NOT provided, but ES & totalN are available
		if (!is.null(ES) & is.null(vi) & grpinfo == FALSE & !is.null(N)) {

			# there is no way of obtaining vi for d or g, so converting both ES & ES_type_IN to r

			# message('\n\n\nA d or g effect size was specified, but without vi or sufficient information')
			# message('about the two groups. vi cannot be computed, and meta-analyses cannot be conducted.')
			# message('The total N and ES were provided, and so the ES and all subsequent analyses')
			# message('will be for correlation coefficient equivalents, which can be computed for the provided data.\n\n')

			rfromgd <- CONVERT_ES(ES = donnes[,ES], ES_var = NULL, 
			                     ES_type_IN=ES_type_IN, ES_type_OUT='r', totalN = donnes[,N], verbose = FALSE)
			
			donnes_RN <- data.frame(N = donnes[,N], ES = rfromgd$r) 

			ESdat <- escalc(measure='COR', ni=donnes_RN$N, ri=donnes_RN$ES )

			ES_type_IN = 'r'
		}
	}		

	if (dtype == 'list') {
		# ESdat_grpinfo  <- data.frame(yi = rep(NA,length(donnes)), vi = rep(NA,length(donnes)))
		donnes_RN <- data.frame(N =  rep(NA,length(donnes)), ES = rep(NA,length(donnes)))
		for (lupe in 1:length(donnes)) {
			
			means <- aggregate(as.matrix(donnes[[lupe]][,2]), list(donnes[[lupe]][,1]), mean)			
			SDs   <- aggregate(as.matrix(donnes[[lupe]][,2]), list(donnes[[lupe]][,1]), sd)
			Ns    <- table(donnes[[lupe]][,1])
				
			# using just the 1st 2 values of the groups variable (col 1 of donnes)
			ESdat_grpinfo[lupe,] <- ES_from_GRPINFO(
			                           grp1_MN=means[1,2], grp1_SD=SDs[1,2], grp1_N=Ns[1], 
				                       grp2_MN=means[2,2], grp2_SD=SDs[2,2], grp2_N=Ns[2],
				                       ES_type_OUT='r')

			# from the escalc Rd:
			# The positive bias in the standardized mean difference is automatically 
			# corrected for within the function, yielding Hedges' g for measure="SMD" 
			# (Hedges, 1981). Similarly, the same bias correction is applied for 
			# measure="SMDH" (Bonett, 2009). For measure="SMD", one can choose between 
			# vtype="LS" (the default) and vtype="UB". The former uses the usual large-sample 
			# approximation to compute the sampling variances. The latter provides unbiased 
			# estimates of the sampling variances.			
								
			# using just the 1st 2 values of the groups variable (col 1 of donnes)
			# ESdat_g[lupe,] <- escalc(measure='SMD', vtype='UB', 
			                         # m1i=means[1,2], sd1i=SDs[1,2], n1i=Ns[1], 
			                         # m2i=means[2,2], sd2i=SDs[2,2], n2i=Ns[2])

			# rfromgd = CONVERT_ES(ES = ESdat_grpinfo$yi[lupe], Vd = ESdat_grpinfo$vi[lupe], dg_type = ES_type_OUT, 
			                    # grp1_N = Ns[1], grp2_N = Ns[2], verbose = FALSE)						
			
			# donnes_RN[lupe,] <- data.frame(N = sum(Ns), ES = rfromgd$r) 

			donnes_RN[lupe,] <- data.frame(N = ESdat_grpinfo$totalN, ES = ESdat_grpinfo$r) 
		}
		ESdat <- escalc(measure='COR', ni=donnes_RN$N, ri=donnes_RN$ES )
	}
}


#totalNs <- attr(ESdat$yi, 'ni')

totalNs <- donnes_RN$N

# identifying & removing any row with an NA from ESdat, totalNs, & donnes_RN
rowswithNAs <- unique(which(is.na(cbind(ESdat, totalNs, donnes_RN)), arr.ind=TRUE)[,1])
if (length(rowswithNAs) >= 1) {
	ESdat <- ESdat[-rowswithNAs,]
	totalNs <- totalNs[-rowswithNAs]
	donnes_RN <- donnes_RN[-rowswithNAs,] 
	message('\n\nMissing values were found and removed from the data matrix.\n\n')
}

                   
Nstudies <- nrow(ESdat)
                      
               
                
##############################  NHST  #######################################


summESdat <- summary(ESdat)  # get CIs, only works if escalc was used

results_NHST <- cbind(1:Nstudies, summESdat$ci.lb, summESdat$yi, summESdat$ci.ub)


# if (class(ESdat)[1] == 'escalc')  {
	# summESdat <- summary(ESdat)  # get CIs, only works if escalc was used
	# results_NHST <- cbind(1:Nstudies, summESdat$ci.lb, summESdat$yi, summESdat$ci.ub)
# } else {	
	# ES.lb <- ESdat$yi - 1.96 * sqrt(ESdat$vi)
	# ES.ub <- ESdat$yi + 1.96 * sqrt(ESdat$vi)		
	# results_NHST <- cbind(1:Nstudies, ES.lb, ESdat$yi, ES.ub)	

	# # computing CIs for r, assuming they are not symetrical -- Loftus
	# FisherZ <- as.matrix(.5 * log( (1 + ESdat$yi) / (1 - ESdat$yi))
	# sdr <- sqrt( 1 / ( totalNs - 3))
	# zUBnhst <- FisherZ + 1.96 * sdr
	# zLBnhst <- FisherZ - 1.96 * sdr
	# ES.ub <- tanh(zUBnhst)
	# ES.lb <- tanh(zLBnhst)
# }

dimnames(results_NHST) <-list(rep("", dim(results_NHST)[1]))
colnames(results_NHST) <- c('Study','ES.lb','ES','ES.ub')


if (verbose) {
	message('\n\nStudy effect sizes and conventional confidence intervals in r metric:\n')
	print(round(results_NHST,3), print.gap=4)
}



####################################  regular meta-analysis  #######################################

outp_MA_1 <- rma(yi=ESdat$yi, vi=ESdat$vi, method=rma_method)  # random-effects model

# confidence intervals for tau & Isq
heteroests <- confint(outp_MA_1)

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


if (verbose) {
	
	message('\n\nResults from a regular, all-studies-at-once, random-effects meta-analysis')
	message('in r, d, & g effect size metrics:\n')
	
	# message('\n   r = ',round(outp_MA_1$b,2),'   r.lb = ',round(outp_MA_1$ci.lb,2),
	        # '    r.ub = ',round(outp_MA_1$ci.ub,2))

	rres <- rbind(outp_MA_1$b, outp_MA_1$ci.lb, outp_MA_1$ci.ub)
	esmat <- t(rres)
		
	dres <- CONVERT_ES(ES = rres, ES_var = rep(outp_MA_1$vb,3), ES_type_IN='r', ES_type_OUT='d', 
                       verbose = FALSE)
	esmat <- rbind(esmat, dres$d)	
	
	gres <- CONVERT_ES(ES = rres, ES_var = rep(outp_MA_1$vb,3), ES_type_IN='r', ES_type_OUT='g', 
                       totalN = rep(outp_MA_1$k,3), verbose = FALSE)
	esmat <- rbind(esmat, gres$g)	
		
	colnames(esmat) <- c('ES','ES.lb','ES.ub')
	rownames(esmat) <- c('r','d','g')
	print(round(esmat,3), print.gap=4)	
}



##################################  publication bias  ###########################################

# cannot run tests for publication bias when the study Ns are all the same
#if (var(donnes_RN[,1]) != 0) {
	funreg <- regtest(outp_MA_1)  # Regression Test for Funnel Plot Asymmetry -- metafor
	funrank <- suppressWarnings(ranktest(outp_MA_1))  #  Rank Correlation Test for Funnel Plot Asymmetry  -- metafor

	# # the metabias (& metacor) function is from the meta package, so not using
	# m1 <- metacor(cor=as.matrix(donnes_RN[,2]), as.matrix(donnes_RN[,1]),  method.bias='linreg') # -- meta
	# thomsharp <- metabias(m1, plotit=F, method='mm') # a variant of Eggers's test allowing for between-study heterogeneity	

	# biasStats <- cbind(funreg$zval, funreg$pval, funrank$tau, funrank$pval, thomsharp$statistic, thomsharp$p.value)
	# colnames(biasStats) <- 
	  # c('  reg test z','  reg test p','     rank test tau','  rank test p','     Thom/Sharp t','   Thom/Sharp p')
	# message('\n\nTests for Publication Bias:\n\n')
	# print(round(biasStats,3))	

	if (verbose) {
		message('\n\nTests for Publication Bias:')
		
		message('\n   Regression Test for Funnel Plot Asymmetry:  z = ', round(funreg$zval,2), 
		        '   p = ', round(funreg$pval))
	
		message('\n   Rank Correlation Test for Funnel Plot Asymmetry:  tau = ', 
		        round(funrank$tau,2), '   p = ', round(funrank$pval))
	
		# message('\n   Thom/Sharp:  t = ', round(thomsharp$statistic,2), 
		        # '   p = ', round(thomsharp$p.value))
	}

#	write.table(round(biasStats,3), "biasStats.csv", sep=',', append=TRUE, col.names=F, row.names=titre)
#	bitmap(paste("Funnel Plot - ",titre,".png",sep=""), height=7, width=9, units='in', res=1200, pointsize=12)
#	par( pty="m", mar=c(3,2,3,2) + 2.6)  #     1.8
#	metafor::funnel(outp_MA_1, main=titre)
#	dev.off()
#}

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


outp_CUM <- cumul(outp_MA_1)  # cumulative meta-analysis (in the order of publication year)

# message('\n\n Cum MA in r metric:'); 
# print(outp_CUM)


CUM_ES    <- outp_CUM$estimate
CUM_tau2  <- outp_CUM$tau2
CUM_tau   <- sqrt(outp_CUM$tau2)
CUM_se    <- outp_CUM$se

# CUM_ES_lb <- outp_CUM$ci.lb
# CUM_ES_ub <- outp_CUM$ci.ub

# # computing CIs assuming they are not symetrical -- Loftus
CUM_z <- as.matrix(.5 * log((1 + CUM_ES) / (1 - CUM_ES)))
#CUM_tau <- sqrt(CUM_tau2)
zUB <- CUM_z + zforCI * CUM_se
zLB <- CUM_z - zforCI * CUM_se
CUM_ES_ub <- tanh(zUB)
CUM_ES_lb <- tanh(zLB)


results_CUM <- cbind(1:Nstudies, totalNs, ESdat$yi, CUM_ES_lb, CUM_ES, CUM_ES_ub)

	# # setting values > 1 to NA
	# results_CUM[,4:6] <- ifelse( results_CUM[,4:6] >  1,  NA, results_CUM[,4:6])
	# results_CUM[,4:6] <- ifelse( results_CUM[,4:6] < -1,  NA, results_CUM[,4:6])

if (ES_type_OUT == 'd' | ES_type_OUT == 'g') {
	results_CUM_dg <- results_CUM
	for (lupe in 3:6) {
		cumdres <- CONVERT_ES(ES = results_CUM[,lupe], ES_var = NULL, 
		                      ES_type_IN='r', ES_type_OUT=ES_type_OUT, totalN = totalNs, verbose = FALSE)
		                      
		if (ES_type_OUT == 'd') results_CUM_dg[,lupe] <- cumdres$d
		if (ES_type_OUT == 'g') results_CUM_dg[,lupe] <- cumdres$g
		                      
#		results_CUM_dg[,lupe] <- t(cumdres[[,1]])	
	}	
	# replacing results_CUM with results_CUM_dg
	results_CUM <- results_CUM_dg
}	

if (verbose) {
	message('\n\n\nThe specified type of effect size for the output below = ', ES_type_OUT)
}

dimnames(results_CUM) <-list(rep("", dim(results_CUM)[1]))
colnames(results_CUM) <- 
         c('Study','Study N','Study ES','ES.lb','ES','ES.ub')

if (verbose) {
	message('\n\n\nRandom-effects cumulative meta-analysis estimates:\n')
	print(round(results_CUM,3), print.gap=4)
}



##############################  ESdat, based on r & N, for the Bayesian analyses  ##########################


# # summESdat <- summary(ESdat)   


# if (is.null(donnes_RN) == FALSE) {

	# ESdat_for_BA <- escalc(measure='COR', ni=donnes_RN$N, ri=donnes_RN$ES )

	# summESdat_for_BA <- summary(ESdat_for_BA)   
# }



##########################################  Bayes - Schmidt-Raju (2007) method  ####################################

# random effects empirical Bayes meta-analysis 
# as in 2007 Schmidt, Raju - Updating meta-analytic research findings - Bayesian approaches versus the medical model

if (is.element('Schmidt_Raju', Bayes_type) & !is.null(donnes_RN)) {


# loop through the studies, treating each subsequent study as the Likelihood
results_BA_SR <- matrix(NA,nrow(ESdat),5)

# the 1st row of results = regular CIs
results_BA_SR[1,3:5] <- cbind(summESdat$ci.lb[1], summESdat$yi[1], summESdat$ci.ub[1])


for (luper in 2:nrow(ESdat)) {

	priordatSR <- ESdat[1:(luper-1),]
	
	outp_MA_3 <- rma(yi=priordatSR$yi, vi=priordatSR$vi, method=rma_method, 
	                 control=list(stepadj=0.5, maxiter=1000))  # random-effects model
	
	BA_SR_priorES <- outp_MA_3$b
	BA_SR_priorSE <- outp_MA_3$se
	BA_SR_priorPopV <- outp_MA_3$tau2
	BA_SR_Vprior <- BA_SR_priorSE**2
	k <- nrow(priordatSR)
	
	rk <- donnes_RN[luper,'ES']
	Nk <- donnes_RN[luper,'N']
	
	Vrk <- ( (1 - rk**2)**2 / (Nk - 1) )  + BA_SR_priorPopV
	
	BA_SR_r <- ( (BA_SR_priorES / BA_SR_Vprior) + (rk / Vrk) ) / 
	          ( (1 / BA_SR_Vprior) + (1 / Vrk))  # Schmidt & Raju, 2007, p 302, formula 6
	
	BA_SR_postV <- (BA_SR_Vprior * Vrk) / (BA_SR_Vprior + Vrk)  # Schmidt & Raju, 2007, p 302, formula 7
	
	BA_SR_postPopV <- ( (k - 1) * BA_SR_priorPopV + (k**2 * BA_SR_postV - (k - 1)**2 * 
	                   BA_SR_Vprior) - Vrk) / k # Schmidt & Raju, 2007, p 303, formula 8 -- same as cc$tau2
	
	# credibility intervals -- see Field 2005 p 448, formula 16
	# Hunter and Schmidt recommend correcting this estimate for artifacts 
	# (see Hunter & Schmidt, 2004, or Hall & Brannick, 2002 for details) and then
	# constructing what they call credibility intervals. These intervals are based
	# on taking the average correlation (see Equation 12) and adding to or subtracting
	# from it the square root of the estimated population variance in 
	# Equation 15 multiplied by ... (1.96 for a 95% interval) 
	
	BA_SR_rub <- BA_SR_r + 1.96 * sqrt(BA_SR_postV)  # sqrt(BA_SR_postV) = cc$se
	BA_SR_rlb <- BA_SR_r - 1.96 * sqrt(BA_SR_postV)
	
	# confidence intervals -- see Field 2005 p 448, formula 17
	# If confidence intervals are required (rather than credibility intervals) these can
	# be obtained by using the standard error of the mean correlation. To obtain this
	# standard error simply divide the variance of sample correlations (given in Equation 13)
	# by the number of studies in the meta-analysis, k, and take the square root: 
	# BA_SR_rub <- BA_SR_r + 1.96 * sqrt(BA_SR_postV / k)
	# BA_SR_rlb <- BA_SR_r - 1.96 * sqrt(BA_SR_postV / k)
	
	results_BA_SR[luper,] <- cbind(BA_SR_priorES, BA_SR_priorSE, # BA_SR_postV, BA_SR_postPopV, 
	                               BA_SR_rlb, BA_SR_r, BA_SR_rub)

	postN_Vs <- cbind(donnes_RN[luper,1], BA_SR_postV, BA_SR_postPopV)
	}

results_BA_SR <- data.frame(cbind( 1:Nstudies, as.matrix(donnes_RN), results_BA_SR ))
#dimnames(results_BA_SR) <- list(rep("", dim(results_BA_SR)[1]))
colnames(results_BA_SR) <- 
    c('Study','Study N','Study ES','prior ES','prior SE','ES.lb','ES','ES.ub')  


# convert r to d or g effect sizes, if requested
if (ES_type_OUT == 'd' | ES_type_OUT == 'g') {
	results_BA_SR_dg <- results_BA_SR[,c('Study','Study N','Study ES','prior ES','ES.lb','ES','ES.ub')]
	for (lupe in 3:7) {
		results_BA_SR_dg[,lupe] <-
			CONVERT_ES(ES = results_BA_SR_dg[,lupe], ES_var = NULL, 
		              ES_type_IN='r', ES_type_OUT=ES_type_OUT, 
		              totalN=results_BA_SR_dg[,'Study N'], verbose = FALSE)[paste(ES_type_OUT)]
	}
	results_BA_SR <- results_BA_SR_dg
}


if (verbose) {
	message('\n\n\nBayesian estimates based on the Schmidt & Raju (2007) method:\n')

	print(round(results_BA_SR,3), print.gap=4)

	# if (ES_type_OUT == 'r') print(round(results_BA_SR,3), print.gap=4)
	
	# if (ES_type_OUT == 'd' | ES_type_OUT == 'g') print(round(results_BA_SR_dg,3), print.gap=4)
}

}



##########################################  Bayes generated data  #############################################

 
if (is.element('generated', Bayes_type) & !is.null(donnes_RN)) {

# generate likelihood raw data for 2 variables with a correlation = the current effect size

# run MCMCglmm, using the effect size & sampling error variance from a MA of previous data

# loop through the studies, treating each subsequent study as the Likelihood
results_BA_GEN <- matrix(NA,Nstudies,5)
BA_GEN_post_ests <- matrix(NA,Nstudies,2)
#if (dtype == 'not list') donnesRAW <- list()

for (luper in 1:Nstudies) {
	
	correl <- donnes_RN[luper,'ES']
	N      <- round(donnes_RN[luper,'N'])
	
	# generating variables with an exact correlation
	dataset1 <- data.frame(mvrnorm(n=N,mu=c(0,0),
						   Sigma=matrix(c(1,correl,correl, 1),nrow=2),empirical=TRUE))
	colnames(dataset1) <- c('varIV','varDV')
	
	# # saving the generated data for the Bayes Factor analyses (when there is no raw data)
	# if (dtype == 'not list') donnesRAW[[luper]] <- dataset1
	
	if (luper == 1) { 
		model1 <- MCMCglmm(varDV ~ varIV, data=dataset1, nitt=nitt, burnin=burnin, thin=thin, verbose=FALSE)
		model1sum <- summary.MCMCglmm(model1)
		BA_GEN_r   <- model1sum$solutions[2,1]
		BA_GEN_rlb <- model1sum$solutions[2,2]
		BA_GEN_rub <- model1sum$solutions[2,3]
		BA_GEN_priorES <- NA
		BA_GEN_priorSE <- NA
		BA_GEN_post_ests[1,] <- c(BA_GEN_r, (diag(var(model1$Sol)))[2])  # saving the ES & the posterior variance
	}	
		  
	if (luper > 1) { 
		if (prior_type == 'MA') { 		
			priordat <- ESdat[1:(luper-1),]
			outp_MA_2 <- rma(yi=priordat$yi, vi=priordat$vi, method=rma_method, 
	                         control=list(stepadj=0.5, maxiter=1000))  # random-effects model

			BA_GEN_priorES <- outp_MA_2$b
			BA_GEN_priorSE <- outp_MA_2$se
#			BA_GEN_priorPopV <- outp_MA_2$tau2
			BA_GEN_Vprior <- BA_GEN_priorSE**2			
		}
		if (prior_type == 'Bayes') {	
			# Jarrod Hadfield  https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q4/019436.html
			# v <- var(model$Sol)  the variance-covariance matrix
			# sqrt(diag(v))  the posterior standard deviations (akin to the standard errors) 
			BA_GEN_priorES <- BA_GEN_post_ests[(luper-1),1]
			BA_GEN_Vprior  <- BA_GEN_post_ests[(luper-1),2]
			BA_GEN_priorSE <- sqrt(BA_GEN_Vprior)
		}
		priors <- list(B=list(mu=c(0, BA_GEN_priorES), V=diag(c(1, BA_GEN_Vprior))))		
		model2 <- MCMCglmm(varDV ~ varIV, data=dataset1, prior=priors, nitt=nitt, burnin=burnin, 
		                   thin=thin, verbose=FALSE)
		model2sum  <- summary.MCMCglmm(model2)
		BA_GEN_r   <- model2sum$solutions[2,1]
		BA_GEN_rlb <- model2sum$solutions[2,2]
		BA_GEN_rub <- model2sum$solutions[2,3]
		BA_GEN_post_ests[luper,] <- c(BA_GEN_r, (diag(var(model2$Sol)))[2])  # the ES & the posterior variance
	}

results_BA_GEN[luper,] <- cbind(BA_GEN_priorES, BA_GEN_priorSE, BA_GEN_rlb, BA_GEN_r, BA_GEN_rub)

}

results_BA_GEN <- data.frame(cbind( 1:Nstudies, as.matrix(donnes_RN), results_BA_GEN ))
#dimnames(results_BA_GEN) <-list(rep("", dim(results_BA_GEN)[1]))
colnames(results_BA_GEN) <- 
    c('Study','Study N','Study ES','prior ES','prior SE','ES.lb','ES','ES.ub')


# convert r to d or g effect sizes, if requested
if (ES_type_OUT == 'd' | ES_type_OUT == 'g') {
	results_BA_GEN_dg <- results_BA_GEN[,c('Study','Study N','Study ES','prior ES','ES.lb','ES','ES.ub')]
	for (lupe in 3:7) {
		results_BA_GEN_dg[,lupe] <-
			CONVERT_ES(ES = results_BA_GEN_dg[,lupe], ES_var = NULL, 
		              ES_type_IN='r', ES_type_OUT=ES_type_OUT, 
		              totalN=results_BA_GEN_dg[,'Study N'], verbose = FALSE)[paste(ES_type_OUT)]
	}
	results_BA_GEN <- results_BA_GEN_dg
}


if (verbose) {
	message('\n\n\nBayesian estimates based on generated data:\n')

	print(round(results_BA_GEN,3), print.gap=4)

	# if (ES_type_OUT == 'r') print(round(results_BA_GEN,3), print.gap=4)
	
	# if (ES_type_OUT == 'd' | ES_type_OUT == 'g') print(round(results_BA_GEN_dg,3), print.gap=4)
}

}




###########################  Bayes raw data  #############################


if (is.element('raw', Bayes_type) & dtype != 'list') {
	message('\n\nBayesian raw data analyses were requested but donnes is not a list with raw data.')
	message('Bayesian raw data analyses cannot be conducted without raw data.\n')
}

if (is.element('raw', Bayes_type) & dtype == 'list') {

BA_RAW_rlb <- BA_RAW_r <- BA_RAW_rub <- NA

# run MCMCglmm, using the effect size & sampling error variance from a MA of previous data

# loop through the studies, treating each subsequent study as the Likelihood
results_BA_RAW <- matrix(NA,Nstudies,5)
BA_RAW_post_ests <- matrix(NA,Nstudies,2)


for (luper in 1:Nstudies) {

	dataset1 <- data.frame(donnesRAW[[luper]]); colnames(dataset1) <- c('varIV','varDV')
	
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
		BA_RAW_r   <- model1sum$solutions[2,1]
		BA_RAW_rlb <- model1sum$solutions[2,2]
		BA_RAW_rub <- model1sum$solutions[2,3]
		BA_RAW_priorES <- NA
		BA_RAW_priorSE <- NA
		BA_RAW_post_ests[1,] <- c(BA_RAW_r, (diag(var(model1$Sol)))[2])  # the ES & posterior variance
	}	
		  
	if (luper > 1) { 		
		if (prior_type == 'MA') {	
			priordat <- ESdat[1:(luper-1),]
			outp_MA_3 <- rma(yi=priordat$yi, vi=priordat$vi, method=rma_method, 
	                         control=list(stepadj=0.5, maxiter=1000))  # random-effects model


			BA_RAW_priorES <- outp_MA_3$b
			BA_RAW_priorSE <- outp_MA_3$se
	#		BA_RAW_priorPopV <- outp_MA_3$tau2						
			BA_RAW_Vprior <- BA_RAW_priorSE**2	
		}
		if (prior_type == 'Bayes') {	
			# Jarrod Hadfield  https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q4/019436.html
			# v <- var(model$Sol)  the variance-covariance matrix
			# sqrt(diag(v))  the posterior standard deviations (akin to the standard errors) 
			BA_RAW_priorES <- BA_RAW_post_ests[(luper-1),1]
			BA_RAW_Vprior  <- BA_RAW_post_ests[(luper-1),2]
			BA_RAW_priorSE <- sqrt(BA_RAW_Vprior)
		}
		priors <- list(B=list(mu=c(0, BA_RAW_priorES), V=diag(c(1, BA_RAW_Vprior))))		
		model2 <- MCMCglmm(varDV ~ varIV, data=dataset1, prior=priors, nitt=nitt, burnin=burnin, 
		                   thin=thin, verbose=FALSE)
		model2sum <- summary.MCMCglmm(model2)
		BA_RAW_r   <- model2sum$solutions[2,1]
		BA_RAW_rlb <- model2sum$solutions[2,2]
		BA_RAW_rub <- model2sum$solutions[2,3]
		BA_RAW_post_ests[luper,] <- c(BA_RAW_r, (diag(var(model2$Sol)))[2]) # ES & the posterior variance
	}	
	results_BA_RAW[luper,] <- cbind(BA_RAW_priorES, BA_RAW_priorSE, BA_RAW_rlb, BA_RAW_r, BA_RAW_rub)
}

results_BA_RAW <- data.frame(cbind(1:Nstudies, as.matrix(donnes_RN), results_BA_RAW))
#dimnames(results_BA_RAW) <-list(rep("", dim(results_BA_RAW)[1]))
colnames(results_BA_RAW) <- 
  c('Study','Study N','Study ES','prior ES','prior SE','ES.lb','ES','ES.ub')


# convert r to d or g effect sizes, if requested
if (ES_type_OUT == 'd' | ES_type_OUT == 'g') {
	results_BA_RAW_dg <- results_BA_RAW[,c('Study','Study N','Study ES','prior ES','ES.lb','ES','ES.ub')]
	for (lupe in 3:7) {
		results_BA_RAW_dg[,lupe] <-
			CONVERT_ES(ES = results_BA_RAW_dg[,lupe], ES_var = NULL, 
		              ES_type_IN='r', ES_type_OUT=ES_type_OUT, 
		              totalN=results_BA_RAW_dg[,'Study N'], verbose = FALSE)[paste(ES_type_OUT)]
	}
	results_BA_RAW <- results_BA_RAW_dg
}


if (verbose) {
	message('\n\n\nBayesian estimates based on the raw data:\n')

	print(round(results_BA_RAW,3), print.gap=4)

	# if (ES_type_OUT == 'r') print(round(results_BA_RAW,3), print.gap=4)
	
	# if (ES_type_OUT == 'd' | ES_type_OUT == 'g') print(round(results_BA_RAW_dg,3), print.gap=4)
}

}




###########################  agreement & consistency  ##########################


# the CUM MA conclusion
if (results_CUM[,'ES.lb'][Nstudies] < 0 & results_CUM[,'ES.ub'][Nstudies] < 0)  effect <- 'negeff'
if (results_CUM[,'ES.lb'][Nstudies] > 0 & results_CUM[,'ES.ub'][Nstudies] > 0)  effect <- 'poseff'
if (results_CUM[,'ES.lb'][Nstudies] < 0 & results_CUM[,'ES.ub'][Nstudies] > 0)  effect <- 'noeff'

# NHST -- agreement with final, & consistency
sigposNHST <- signegNHST <- nonsigNHST <- 0
for (luper in 1:nrow(results_NHST)) {
	if (results_NHST[luper,'ES.lb'] < 0 & results_NHST[luper,'ES.ub'] < 0) signegNHST <- signegNHST + 1
	if (results_NHST[luper,'ES.lb'] > 0 & results_NHST[luper,'ES.ub'] > 0) sigposNHST <- sigposNHST + 1
	if (results_NHST[luper,'ES.lb'] < 0 & results_NHST[luper,'ES.ub'] > 0) nonsigNHST <- nonsigNHST + 1
}

# percentage of NHST studies that agreed with the CUM MA conclusion
if (effect == 'negeff') agreeNHST <- signegNHST / nrow(results_NHST)
if (effect == 'poseff') agreeNHST <- sigposNHST / nrow(results_NHST)
if (effect == 'noeff')  agreeNHST <- nonsigNHST / nrow(results_NHST)

# the most common conclusion, as a proportion
consistNHST <- max(c(sigposNHST,signegNHST,nonsigNHST)) / nrow(results_NHST) 



# CUM  -- agreement with final, & consistency
sigposCUM <- signegCUM <- nonsigCUM <- 0
for (luper in 1:nrow(results_CUM)) {
	if (results_CUM[luper,'ES.lb'] < 0 & results_CUM[luper,'ES.ub'] < 0) signegCUM <- signegCUM + 1
	if (results_CUM[luper,'ES.lb'] > 0 & results_CUM[luper,'ES.ub'] > 0) sigposCUM <- sigposCUM + 1
	if (results_CUM[luper,'ES.lb'] < 0 & results_CUM[luper,'ES.ub'] > 0) nonsigCUM <- nonsigCUM + 1
}

# percentage of CUM studies that agreed with the CUM MA conclusion
if (effect == 'negeff') agreeCUM <- signegCUM / nrow(results_CUM)
if (effect == 'poseff') agreeCUM <- sigposCUM / nrow(results_CUM)
if (effect == 'noeff')  agreeCUM <- nonsigCUM / nrow(results_CUM)

# the most common conclusion, as a proportion
consistCUM <- max(c(sigposCUM,signegCUM,nonsigCUM)) / nrow(results_CUM)


#print(results_BA_SR)


if (is.element('Schmidt_Raju', Bayes_type)) {
	agreeBA_SR <- consistBA_SR <- NA
	# BAYES_SR -- agreement with final, & consistency
	sigposBA_SR <- signegBA_SR <- nonsigBA_SR <- 0
	for (luper in 1:nrow(results_BA_SR)) {
		if (results_BA_SR[luper,'ES.lb'] < 0 & results_BA_SR[luper,'ES.ub'] < 0) signegBA_SR <- signegBA_SR + 1
		if (results_BA_SR[luper,'ES.lb'] > 0 & results_BA_SR[luper,'ES.ub'] > 0) sigposBA_SR <- sigposBA_SR + 1
		if (results_BA_SR[luper,'ES.lb'] < 0 & results_BA_SR[luper,'ES.ub'] > 0) nonsigBA_SR <- nonsigBA_SR + 1
	}
	# percentage of BA_SR studies that agreed with the CUM MA conclusion
	if (effect == 'negeff') agreeBA_SR <- signegBA_SR / nrow(results_BA_SR)
	if (effect == 'poseff') agreeBA_SR <- sigposBA_SR / nrow(results_BA_SR)
	if (effect == 'noeff')  agreeBA_SR <- nonsigBA_SR / nrow(results_BA_SR)
	
	# the most common conclusion, as a proportion
	consistBA_SR <- max(c(sigposBA_SR,signegBA_SR,nonsigBA_SR)) / nrow(results_BA_SR) 
}


if (is.element('generated', Bayes_type)) {
	# generated BAYES -- agreement with final, & consistency
	sigposBA_GEN <- signegBA_GEN <- nonsigBA_GEN <- 0
	for (luper in 1:nrow(results_BA_GEN)) {
		if (results_BA_GEN[luper,'ES.lb'] < 0 & results_BA_GEN[luper,'ES.ub'] < 0) signegBA_GEN <- signegBA_GEN + 1
		if (results_BA_GEN[luper,'ES.lb'] > 0 & results_BA_GEN[luper,'ES.ub'] > 0) sigposBA_GEN <- sigposBA_GEN + 1
		if (results_BA_GEN[luper,'ES.lb'] < 0 & results_BA_GEN[luper,'ES.ub'] > 0) nonsigBA_GEN <- nonsigBA_GEN + 1
	}
	# percentage of BA_GEN studies that agreed with the CUM MA conclusion
	if (effect == 'negeff') agreeBA_GEN <- signegBA_GEN / nrow(results_BA_GEN)
	if (effect == 'poseff') agreeBA_GEN <- sigposBA_GEN / nrow(results_BA_GEN)
	if (effect == 'noeff')  agreeBA_GEN <- nonsigBA_GEN / nrow(results_BA_GEN)
	
	# the most common conclusion, as a proportion
	consistBA_GEN <- max(c(sigposBA_GEN,signegBA_GEN,nonsigBA_GEN)) / nrow(results_BA_GEN) 
}


if (is.element('raw', Bayes_type) & dtype == 'list') {
	agreeBA_RAW <- consistBA_RAW <- NA
	# raw BAYES -- agreement with final, & consistency
	sigposBA_RAW <- signegBA_RAW <- nonsigBA_RAW <- 0
	for (luper in 1:nrow(results_BA_RAW)) {
		if (results_BA_RAW[luper,'ES.lb'] < 0 & results_BA_RAW[luper,'ES.ub'] < 0) signegBA_RAW <- signegBA_RAW + 1
		if (results_BA_RAW[luper,'ES.lb'] > 0 & results_BA_RAW[luper,'ES.ub'] > 0) sigposBA_RAW <- sigposBA_RAW + 1
		if (results_BA_RAW[luper,'ES.lb'] < 0 & results_BA_RAW[luper,'ES.ub'] > 0) nonsigBA_RAW <- nonsigBA_RAW + 1
	}
	# percentage of BA_RAW studies that agreed with the CUM MA conclusion
	if (effect == 'negeff') agreeBA_RAW <- signegBA_RAW / nrow(results_BA_RAW)
	if (effect == 'poseff') agreeBA_RAW <- sigposBA_RAW / nrow(results_BA_RAW)
	if (effect == 'noeff')  agreeBA_RAW <- nonsigBA_RAW / nrow(results_BA_RAW)
	
	# the most common conclusion, as a proportion
	consistBA_RAW <- max(c(sigposBA_RAW,signegBA_RAW,nonsigBA_RAW)) / nrow(results_BA_RAW) 
}





####################################  output  #######################################

if (verbose) {

message('\n\nNumber of studies: ', Nstudies)

message('\n\nThe confidence interval for the analyses: ', CI)

message('\n\nFinal, all-studies-combined results:')

message('\n   Cumulative Meta-Analysis:    ES = ', round(results_CUM[,'ES'][Nstudies],3),
        '   ES.lb = ', round(results_CUM[,'ES.lb'][Nstudies],3),
        '   ES.ub = ', round(results_CUM[,'ES.ub'][Nstudies],3))

if (is.element('Schmidt_Raju', Bayes_type)) 
	message('\n   Bayesian (Schmidt-Raju):     ES = ', round(results_BA_SR$ES[Nstudies],3),
	        '   ES.lb = ', round(results_BA_SR$ES.lb[Nstudies],3),
	        '   ES.ub = ', round(results_BA_SR$ES.ub[Nstudies],3))

if (is.element('generated', Bayes_type)) 
	message('\n   Bayesian (generated data):   ES = ', round(results_BA_GEN$ES[Nstudies],3),
	        '   ES.lb = ', round(results_BA_GEN$ES.lb[Nstudies],3),
	        '   ES.ub = ', round(results_BA_GEN$ES.ub[Nstudies],3))

if (is.element('raw', Bayes_type) & dtype == 'list') 
	message('\n   Bayesian (raw data):         ES = ', round(results_BA_RAW$ES[Nstudies],3),
	        '   ES.lb = ', round(results_BA_RAW$ES.lb[Nstudies],3),
	        '   ES.ub = ', round(results_BA_RAW$ES.ub[Nstudies],3))


message('\n\n\nMeta-analysis heterogeneity statistics:')
message('\n   Total Var. = ',round(tau2,3))
message('\n   Q = ', round(outp_MA_1$QE,3), '    p = ', round(outp_MA_1$QEp,5))
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


message('\n\nConsistency & agreement rates:')

message('\n   NHST:                               consistency = ', 
	round(consistNHST,3), '    agreement = ', round(agreeNHST,3))

message('\n   Cumulative Meta-Analysis:           consistency = ', 
	round(consistCUM,3),  '    agreement = ', round(agreeCUM,3))

if (is.element('Schmidt_Raju', Bayes_type)) 
	message('\n   Bayesian (Schmidt & Raju, 2007):    consistency = ', 
		round(consistBA_SR,3),   '    agreement = ', round(agreeBA_SR,3))

if (is.element('generated', Bayes_type)) 
	message('\n   Bayesian (generated data):          consistency = ', 
		round(consistBA_GEN,3),   '    agreement = ', round(agreeBA_GEN,3))

if (is.element('raw', Bayes_type) & dtype == 'list') 
	message('\n   Bayesian (raw data):                consistency = ', 
		round(consistBA_RAW,3),   '    agreement = ', round(agreeBA_RAW,3))

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
	results_CUM = results_CUM, 
	ES_MA    = results_CUM[,'ES'][Nstudies],  
	ES_MA_lb = results_CUM[,'ES.lb'][Nstudies],  
	ES_MA_ub = results_CUM[,'ES.ub'][Nstudies],
	Q = outp_MA_1$QE, p_Q = outp_MA_1$QEp,
	tau2 = tau2, tau2LB = tau2LB, tau2UB = tau2UB,
	tau = tau, tauLB = tauLB, tauUB = tauUB,
	isq = isq, isqLB = isqLB, isqUB = isqUB,
	hsq = hsq, hsqLB = hsqLB, hsqUB = hsqUB,
	consistNHST = consistNHST, agreeNHST = agreeNHST,
	consistCUM = consistCUM, agreeCUM = agreeCUM )

if (is.element('Schmidt_Raju', Bayes_type))
    nopingpongOutput <- c(nopingpongOutput, list(results_BA_SR = results_BA_SR, 
	consistBA_SR = consistBA_SR, agreeBA_SR = agreeBA_SR) )

if (is.element('generated', Bayes_type))  
    nopingpongOutput <- c(nopingpongOutput, list(results_BA_GEN = results_BA_GEN,
	consistBA_GEN = consistBA_GEN, agreeBA_GEN = agreeBA_GEN ) )
		
if (is.element('raw', Bayes_type) & dtype == 'list')
    nopingpongOutput <- c(nopingpongOutput, list(results_BA_RAW = results_BA_RAW,
	consistBA_RAW = consistBA_RAW, agreeBA_RAW = agreeBA_RAW) )

class(nopingpongOutput) <- "NO.PING.PONG"

return(invisible(nopingpongOutput))

}




