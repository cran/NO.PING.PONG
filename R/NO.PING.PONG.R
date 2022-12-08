 

NO.PING.PONG <- function(donnes, 
                         ES_type_IN = NULL, ES_type_OUT = 'r', 
                         rawdata_type = 'for_correl', 
                         rma_method = 'REML',
                         Bayes_type = c('Schmidt_Raju'), 
                         prior_type = 'META', CI = 95,
                         ES = NULL, N = NULL, vi = NULL,
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


zforCI <- qnorm((1 + CI * .01) / 2) # the z value that corresponds to the specified CI


# is complete info (M, SD, & N) for 2 groups provided?
if (!is.null(grp1_mn) & !is.null(grp1_sd) & !is.null(grp1_n) &
    !is.null(grp2_mn) & !is.null(grp2_sd) & !is.null(grp2_n)  ) 
    {grpinfoALL = TRUE} else {grpinfoALL = FALSE}

if (is.null(N) & !is.null(grp1_n) & !is.null(grp2_n)) {
	donnes$N <- donnes[,grp1_n] + donnes[,grp2_n]
	N <- 'N'
}


# warning if both totalN & grp Ns are missing & dontype is not a list
if (dontype == 'matrix' & is.null(N) & is.null(grp1_n) & is.null(grp2_n)) {
	message('\nThe entered data is a matrix (and not raw data points), and')		
	message('both "N", and "grp1_n" & "grp2_n", are missing/NULL. Bayesian analyses cannot')
	message('be conducted without sample size information.\n')
}




######################  compute ESdat for meta-analyses  &  donnesRN for Bayesian  #####################


ESdat <- NULL   # ESdat will have the yi & vi values for the meta-analyses

# donnesRN has the effect size (r) & N for the Bayesian analyses

donnes_RN <- NULL # won't do Bayesian if this remains null



	
if (dontype == 'matrix' & ES_type_IN == 'r') {

	# ESdat & donnes_RN if ES & N are provided
	if (!is.null(ES) & !is.null(N) ) {

		donnes_RN <- data.frame(N = donnes[,N], ES = donnes[,ES])

		ESdat <- escalc(measure='COR', ni=donnes_RN$N, ri=donnes_RN$ES )
		# totalNs <- attr(ESdat$yi, 'ni')
	}

	# ESdat & donnes_RN if ES & vi are provided, but not N
	if (!is.null(ES) & !is.null(vi)  & is.null(N)) {
		
		donnes[,N] <- (1 - donnes[,ES]^2) / donnes[,vi]  + 2
				
		donnes_RN <- data.frame(N = donnes[,N], ES = donnes[,ES])

		ESdat <- escalc(measure='COR', ni=donnes_RN$N, ri=donnes_RN$ES )
	}
}



if (dontype == 'matrix' & (grpinfoALL | (ES_type_IN == 'd' | ES_type_IN == 'g'))) { 
	
	# need ES & vi, & also need N1 & N2, or totalN, for Bayesian

	# if ES & vi & totalN are provided, but group info is NOT provided
	if (!is.null(ES) & !is.null(vi) & !is.null(N) & grpinfoALL == FALSE) { 

		rfromgd <- CONVERT_ES(ES= donnes[,ES], ES_var = donnes[,vi], 
		                     ES_type_IN=ES_type_IN, ES_type_OUT='r', totalN = donnes[,N], verbose = FALSE)
		
		donnes_RN <- data.frame(N = donnes[,N], ES = rfromgd$r) 

		ESdat <- escalc(measure='COR', ni=donnes_RN$N, ri=donnes_RN$ES )	
	}		
		
	# if vi is not provided, but grp1 & grp2 info is available
	if (is.null(vi) & grpinfoALL) {  

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
	if (!is.null(ES) & is.null(vi) & grpinfoALL == FALSE & !is.null(N)) {

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




if (dontype == 'raw data' & rawdata_type == 'for_correl') {
	donnes_RN <- data.frame(N = rep(NA,length(donnes)), ES = rep(NA,length(donnes)))
	for (lupe in 1:length(donnes)) {
		donnes_RN$N[lupe]  <- nrow(donnes[[lupe]])
		donnes_RN$ES[lupe] <- cor(na.omit(donnes[[lupe]]))[2,1]
	}
	ESdat <- escalc(measure='COR', ni=donnes_RN$N, ri=donnes_RN$ES )	
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
					    m2i=studyStats$studyMN2,  sd2i=studyStats$studySD2, ni=studyStats$studyN, ri=studyStats$studyR )
}



# when only ES & vi are provided
if (dontype == 'matrix' & is.null(N) & !grpinfoALL & !is.null(ES) & !is.null(vi))   ESdat <- escalc(yi=donnes[,ES], vi=donnes[,vi])

	

# the total Ns, if available
if (!is.null(donnes_RN)) {totalNs <- donnes_RN$N} else {totalNs <- rep(NA,nrow(ESdat))}



# identifying & removing any row with an NA from ESdat,   & the same rows from totalNs, & donnes_RN
rowswithNAs <- unique(which(is.na(cbind(ESdat)), arr.ind=TRUE)[,1])   # , totalNs, donnes_RN
if (length(rowswithNAs) >= 1) {
	ESdat <- ESdat[-rowswithNAs,]
	totalNs <- totalNs[-rowswithNAs]
	donnes_RN <- donnes_RN[-rowswithNAs,] 
	message('\n\nMissing values were found and removed from the data matrix.\n\n')
}

                   
Nstudies <- nrow(ESdat)
                      
               

 
######################  data & argument specifications  #####################
 
if (verbose) {

	if (dontype == 'raw data') {
		
		message('\n\nThe input data are raw data (and not study effect sizes)')
		
		if (rawdata_type == 'for_correl') message('\nThe raw data were specified to be for computations of correlations between two variables.')
	
		if (rawdata_type == 'indep_groups') message('\nThe raw data were specified to be based on independent groups (in the first column of donnes).')
	
		if (rawdata_type == 'paired_samples') message('\nThe raw data were specified to be based on paired samples.')
	}

	if (dontype == 'matrix') {

		message('\n\nThe input data are a matrix of study effect sizes (not raw data)')
		
		message('\nThe specified type of input effect sizes is: ', ES_type_IN)
		
		if (grpinfoALL == TRUE) message('\nComplete group info (i.e., M, SD, & N) for 2 groups was provided for the analyses')
	}

	message('\nNumber of studies: ', Nstudies)

	message('\nThe requested type of output effect size is: ', ES_type_OUT)

	message('\nThe specified method option for the meta-analyses: ',  rma_method)

	message('\nThe requested kind(s) of Bayesian analyses: ', paste(Bayes_type, collapse=", "))
	
	message('\nThe type of prior data used in the updating analyses (both cumulative and Bayesian): ',  prior_type) 

	message('\nThe confidence interval for the analyses: ', CI, '%')
}


               
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
	if (!is.null(donnes_RN))  message('\n\nStudy effect sizes and conventional confidence intervals in r metric:\n')
	if (is.null(donnes_RN))  message('\n\nStudy effect sizes and conventional confidence intervals:\n')
	print(round(results_NHST,3), print.gap=4)
}



####################################  regular meta-analysis  #######################################

outp_MA_1 <- rma(yi=ESdat$yi, vi=ESdat$vi, method=rma_method)  # random-effects model

if (rma_method != 'FE') {
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
}

if (rma_method == 'FE') {   # for the output object
	tau2 = NULL; tau2LB = NULL; tau2UB = NULL;
	tau = NULL; tauLB = NULL; tauUB = NULL;
	isq = NULL; isqLB = NULL; isqUB = NULL;
	hsq = NULL; hsqLB = NULL; hsqUB = NULL;
}


if (verbose) {

	if (!is.null(donnes_RN)) {	
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
	
	if (is.null(donnes_RN)) {	
		message('\n\nResults from a regular, all-studies-at-once, random-effects meta-analysis\n')
		rres <- rbind(outp_MA_1$b, outp_MA_1$ci.lb, outp_MA_1$ci.ub)
		esmat <- t(rres)		
		colnames(esmat) <- c('ES','ES.lb','ES.ub')
		rownames(esmat) <- c(' ')
		print(round(esmat,3), print.gap=4)		}
}



##################################  publication bias  ###########################################



#if (var(donnes_RN[,1]) != 0) {  # cannot run some tests for publication bias when the study Ns are all the same
	
	funreg <- regtest(outp_MA_1)  # Regression Test for Funnel Plot Asymmetry -- metafor
	
	funrank <- suppressWarnings(ranktest(outp_MA_1))  # Rank Correlation Test for Funnel Plot Asymmetry  -- metafor

	biasStats <- cbind(funreg$zval, funreg$pval, funrank$tau, funrank$pval)
	colnames(biasStats) <- c('  reg test z','  reg test p','     rank test tau','  rank test p')

	# # the metabias (& metacor) function is from the meta package, so not using
	# m1 <- metacor(cor=as.matrix(donnes_RN[,2]), as.matrix(donnes_RN[,1]),  method.bias='linreg') # -- meta
	# thomsharp <- metabias(m1, plotit=F, method='mm') # a variant of Eggers's test allowing for between-study heterogeneity	

	# biasStats <- cbind(funreg$zval, funreg$pval, funrank$tau, funrank$pval, thomsharp$statistic, thomsharp$p.value)
	# colnames(biasStats) <- 
	  # c('  reg test z','  reg test p','     rank test tau','  rank test p','     Thom/Sharp t','   Thom/Sharp p')


	if (verbose) {

		message('\n\nTests for Publication Bias:')
		
		message('\n   Regression Test for Funnel Plot Asymmetry:  z = ', round(funreg$zval,2), 
		        '   p = ', round(funreg$pval))
	
		message('\n   Rank Correlation Test for Funnel Plot Asymmetry:  tau = ', 
		        round(funrank$tau,2), '   p = ', round(funrank$pval))
	
		# message('\n   Thom/Sharp:  t = ', round(thomsharp$statistic,2), 
		        # '   p = ', round(thomsharp$p.value))
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

		metafor::funnel(outp_MA_1, main=funnel_plot_title)
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


outp_CUM_META <- cumul(outp_MA_1)  # cumulative meta-analysis (in the order of publication year)

# message('\n\n Cum MA in r metric:'); print(outp_CUM_META)


CUM_META_ES    <- outp_CUM_META$estimate
# CUM_META_tau2  <- outp_CUM_META$tau2
# CUM_META_tau   <- sqrt(outp_CUM_META$tau2)
# CUM_META_se    <- outp_CUM_META$se

CUM_META_ES_lb <- outp_CUM_META$ci.lb
CUM_META_ES_ub <- outp_CUM_META$ci.ub

# # # computing CIs assuming they are not symetrical -- Loftus
# CUM_META_z <- as.matrix(.5 * log((1 + CUM_META_ES) / (1 - CUM_META_ES)))
# #CUM_META_tau <- sqrt(CUM_META_tau2)
# zUB <- CUM_META_z + zforCI * CUM_META_se
# zLB <- CUM_META_z - zforCI * CUM_META_se
# CUM_META_ES_ub <- tanh(zUB)
# CUM_META_ES_lb <- tanh(zLB)



if (!is.null(donnes_RN)) {
	
	results_CUM_META <- cbind(1:Nstudies, totalNs, ESdat$yi, CUM_META_ES_lb, CUM_META_ES, CUM_META_ES_ub)

	if (ES_type_OUT == 'd' | ES_type_OUT == 'g') {
	
		results_CUM_META_dg <- results_CUM_META
		for (lupe in 3:6) {
			cumdres <- CONVERT_ES(ES = results_CUM_META[,lupe], ES_var = NULL, 
			                      ES_type_IN='r', ES_type_OUT=ES_type_OUT, totalN = totalNs, verbose = FALSE)
			                      
			if (ES_type_OUT == 'd') results_CUM_META_dg[,lupe] <- cumdres$d
			if (ES_type_OUT == 'g') results_CUM_META_dg[,lupe] <- cumdres$g
		}		
		results_CUM_META <- results_CUM_META_dg  # replacing results_CUM_META with results_CUM_META_dg
	}	
	
	dimnames(results_CUM_META) <-list(rep("", dim(results_CUM_META)[1]))
	colnames(results_CUM_META) <- c('Study','Study N','Study ES','ES.lb','ES','ES.ub')
	
	if (verbose)  message('\n\n\nThe specified type of effect size for the output below = ', ES_type_OUT)
}


if (is.null(donnes_RN)) {

	results_CUM_META <- cbind(1:Nstudies, ESdat$yi, CUM_META_ES_lb, CUM_META_ES, CUM_META_ES_ub)

	dimnames(results_CUM_META) <-list(rep("", dim(results_CUM_META)[1]))
	colnames(results_CUM_META) <- c('Study','Study ES','ES.lb','ES','ES.ub')
}

	# # setting values > 1 to NA
	# results_CUM_META[,4:6] <- ifelse( results_CUM_META[,4:6] >  1,  NA, results_CUM_META[,4:6])
	# results_CUM_META[,4:6] <- ifelse( results_CUM_META[,4:6] < -1,  NA, results_CUM_META[,4:6])


if (verbose) {
	message('\n\n\nRandom-effects cumulative meta-analysis estimates:\n')
	print(round(results_CUM_META,3), print.gap=4)
}




##########################################  Bayes - Schmidt-Raju (2007) method  ####################################

# random effects empirical Bayes meta-analysis 
# as in 2007 Schmidt, Raju - Updating meta-analytic research findings - Bayesian approaches versus the medical model

if (is.element('Schmidt_Raju', Bayes_type) & !is.null(donnes_RN)) {


# loop through the studies, treating each subsequent study as the Likelihood
results_BAYES_SR <- matrix(NA,nrow(ESdat),5)

# the 1st row of results = regular CIs
results_BAYES_SR[1,3:5] <- cbind(summESdat$ci.lb[1], summESdat$yi[1], summESdat$ci.ub[1])


for (luper in 2:nrow(ESdat)) {

	priordatSR <- ESdat[1:(luper-1),]
	
	outp_MA_3 <- rma(yi=priordatSR$yi, vi=priordatSR$vi, method=rma_method, 
	                 control=list(stepadj=0.5, maxiter=1000))  # random-effects model
	
	BAYES_SR_priorES <- outp_MA_3$b
	BAYES_SR_priorSE <- outp_MA_3$se
	BAYES_SR_priorPopV <- outp_MA_3$tau2
	BAYES_SR_Vprior <- BAYES_SR_priorSE**2
	k <- nrow(priordatSR)
	
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
	
	BAYES_SR_rub <- BAYES_SR_r + 1.96 * sqrt(BAYES_SR_postV)  # sqrt(BAYES_SR_postV) = cc$se
	BAYES_SR_rlb <- BAYES_SR_r - 1.96 * sqrt(BAYES_SR_postV)
	
	# confidence intervals -- see Field 2005 p 448, formula 17
	# If confidence intervals are required (rather than credibility intervals) these can
	# be obtained by using the standard error of the mean correlation. To obtain this
	# standard error simply divide the variance of sample correlations (given in Equation 13)
	# by the number of studies in the meta-analysis, k, and take the square root: 
	# BAYES_SR_rub <- BAYES_SR_r + 1.96 * sqrt(BAYES_SR_postV / k)
	# BAYES_SR_rlb <- BAYES_SR_r - 1.96 * sqrt(BAYES_SR_postV / k)
	
	results_BAYES_SR[luper,] <- cbind(BAYES_SR_priorES, BAYES_SR_priorSE, # BAYES_SR_postV, BAYES_SR_postPopV, 
	                               BAYES_SR_rlb, BAYES_SR_r, BAYES_SR_rub)

	postN_Vs <- cbind(donnes_RN[luper,1], BAYES_SR_postV, BAYES_SR_postPopV)
}

results_BAYES_SR <- data.frame(cbind( 1:Nstudies, as.matrix(donnes_RN), results_BAYES_SR ))
# dimnames(results_BAYES_SR) <- list(rep("", dim(results_BAYES_SR)[1]))
colnames(results_BAYES_SR) <- 
    c('Study','Study N','Study ES','prior ES','prior SE','ES.lb','ES','ES.ub')  


	# dimnames(results_CUM_META) <-list(rep("", dim(results_CUM_META)[1]))
	# colnames(results_CUM_META) <- c('Study','Study N','Study ES','ES.lb','ES','ES.ub')


# convert r to d or g effect sizes, if requested
if (ES_type_OUT == 'd' | ES_type_OUT == 'g') {
	results_BAYES_SR_dg <- results_BAYES_SR[,c('Study','Study N','Study ES','prior ES','ES.lb','ES','ES.ub')]
	for (lupe in 3:7) {
		results_BAYES_SR_dg[,lupe] <-
			CONVERT_ES(ES = results_BAYES_SR_dg[,lupe], ES_var = NULL, 
		              ES_type_IN='r', ES_type_OUT=ES_type_OUT, 
		              totalN=results_BAYES_SR_dg[,'Study N'], verbose = FALSE)[paste(ES_type_OUT)]
	}
	results_BAYES_SR <- results_BAYES_SR_dg
}


if (verbose) {
	message('\n\n\nBayesian estimates based on the Schmidt & Raju (2007) method:\n')

	print(round(results_BAYES_SR,3), print.gap=4)
}

}



##########################################  Bayes generated data  #############################################

 
if (is.element('generated', Bayes_type) & !is.null(donnes_RN)) {

# generate likelihood raw data for 2 variables with a correlation = the current effect size

# run MCMCglmm, using the effect size & sampling error variance from a MA of previous data

# loop through the studies, treating each subsequent study as the Likelihood
results_BAYES_GEN <- matrix(NA,Nstudies,5)
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
			priordat <- ESdat[1:(luper-1),]
			outp_MA_2 <- rma(yi=priordat$yi, vi=priordat$vi, method=rma_method, 
	                         control=list(stepadj=0.5, maxiter=1000))  # random-effects model

			BAYES_GEN_priorES <- outp_MA_2$b
			BAYES_GEN_priorSE <- outp_MA_2$se
#			BAYES_GEN_priorPopV <- outp_MA_2$tau2
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
#dimnames(results_BAYES_GEN) <-list(rep("", dim(results_BAYES_GEN)[1]))
colnames(results_BAYES_GEN) <- 
    c('Study','Study N','Study ES','prior ES','prior SE','ES.lb','ES','ES.ub')


# convert r to d or g effect sizes, if requested
if (ES_type_OUT == 'd' | ES_type_OUT == 'g') {
	results_BAYES_GEN_dg <- results_BAYES_GEN[,c('Study','Study N','Study ES','prior ES','ES.lb','ES','ES.ub')]
	for (lupe in 3:7) {
		results_BAYES_GEN_dg[,lupe] <-
			CONVERT_ES(ES = results_BAYES_GEN_dg[,lupe], ES_var = NULL, 
		              ES_type_IN='r', ES_type_OUT=ES_type_OUT, 
		              totalN=results_BAYES_GEN_dg[,'Study N'], verbose = FALSE)[paste(ES_type_OUT)]
	}
	results_BAYES_GEN <- results_BAYES_GEN_dg
}


if (verbose) {
	message('\n\n\nBayesian estimates based on generated data:\n')

	print(round(results_BAYES_GEN,3), print.gap=4)

	# if (ES_type_OUT == 'r') print(round(results_BAYES_GEN,3), print.gap=4)
	
	# if (ES_type_OUT == 'd' | ES_type_OUT == 'g') print(round(results_BAYES_GEN_dg,3), print.gap=4)
}

}




###########################  Bayes raw data  #############################


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
results_BAYES_RAW <- matrix(NA,Nstudies,5)
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
			priordat <- ESdat[1:(luper-1),]
			outp_MA_3 <- rma(yi=priordat$yi, vi=priordat$vi, method=rma_method, 
	                         control=list(stepadj=0.5, maxiter=1000))  # random-effects model


			BAYES_RAW_priorES <- outp_MA_3$b
			BAYES_RAW_priorSE <- outp_MA_3$se
	#		BAYES_RAW_priorPopV <- outp_MA_3$tau2						
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
#dimnames(results_BAYES_RAW) <-list(rep("", dim(results_BAYES_RAW)[1]))
colnames(results_BAYES_RAW) <- 
  c('Study','Study N','Study ES','prior ES','prior SE','ES.lb','ES','ES.ub')


# convert r to d or g effect sizes, if requested
if (ES_type_OUT == 'd' | ES_type_OUT == 'g') {
	results_BAYES_RAW_dg <- results_BAYES_RAW[,c('Study','Study N','Study ES','prior ES','ES.lb','ES','ES.ub')]
	for (lupe in 3:7) {
		results_BAYES_RAW_dg[,lupe] <-
			CONVERT_ES(ES = results_BAYES_RAW_dg[,lupe], ES_var = NULL, 
		              ES_type_IN='r', ES_type_OUT=ES_type_OUT, 
		              totalN=results_BAYES_RAW_dg[,'Study N'], verbose = FALSE)[paste(ES_type_OUT)]
	}
	results_BAYES_RAW <- results_BAYES_RAW_dg
}


if (verbose) {
	message('\n\n\nBayesian estimates based on the raw data:\n')

	print(round(results_BAYES_RAW,3), print.gap=4)
}

}




###########################  agreement & consistency  ##########################


# the CUM_META MA conclusion
if (results_CUM_META[,'ES.lb'][Nstudies] < 0 & results_CUM_META[,'ES.ub'][Nstudies] < 0)  effect <- 'negeff'
if (results_CUM_META[,'ES.lb'][Nstudies] > 0 & results_CUM_META[,'ES.ub'][Nstudies] > 0)  effect <- 'poseff'
if (results_CUM_META[,'ES.lb'][Nstudies] < 0 & results_CUM_META[,'ES.ub'][Nstudies] > 0)  effect <- 'noeff'

# NHST -- agreement with final, & consistency
sigposNHST <- signegNHST <- nonsigNHST <- 0
for (luper in 1:nrow(results_NHST)) {
	if (results_NHST[luper,'ES.lb'] < 0 & results_NHST[luper,'ES.ub'] < 0) signegNHST <- signegNHST + 1
	if (results_NHST[luper,'ES.lb'] > 0 & results_NHST[luper,'ES.ub'] > 0) sigposNHST <- sigposNHST + 1
	if (results_NHST[luper,'ES.lb'] < 0 & results_NHST[luper,'ES.ub'] > 0) nonsigNHST <- nonsigNHST + 1
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
	if (results_CUM_META[luper,'ES.lb'] < 0 & results_CUM_META[luper,'ES.ub'] < 0) signegCUM_META <- signegCUM_META + 1
	if (results_CUM_META[luper,'ES.lb'] > 0 & results_CUM_META[luper,'ES.ub'] > 0) sigposCUM_META <- sigposCUM_META + 1
	if (results_CUM_META[luper,'ES.lb'] < 0 & results_CUM_META[luper,'ES.ub'] > 0) nonsigCUM_META <- nonsigCUM_META + 1
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
			if (results_BAYES_SR[luper,'ES.lb'] < 0 & results_BAYES_SR[luper,'ES.ub'] < 0) signegBAYES_SR <- signegBAYES_SR + 1
			if (results_BAYES_SR[luper,'ES.lb'] > 0 & results_BAYES_SR[luper,'ES.ub'] > 0) sigposBAYES_SR <- sigposBAYES_SR + 1
			if (results_BAYES_SR[luper,'ES.lb'] < 0 & results_BAYES_SR[luper,'ES.ub'] > 0) nonsigBAYES_SR <- nonsigBAYES_SR + 1
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
			if (results_BAYES_GEN[luper,'ES.lb'] < 0 & results_BAYES_GEN[luper,'ES.ub'] < 0) signegBAYES_GEN <- signegBAYES_GEN + 1
			if (results_BAYES_GEN[luper,'ES.lb'] > 0 & results_BAYES_GEN[luper,'ES.ub'] > 0) sigposBAYES_GEN <- sigposBAYES_GEN + 1
			if (results_BAYES_GEN[luper,'ES.lb'] < 0 & results_BAYES_GEN[luper,'ES.ub'] > 0) nonsigBAYES_GEN <- nonsigBAYES_GEN + 1
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
			if (results_BAYES_RAW[luper,'ES.lb'] < 0 & results_BAYES_RAW[luper,'ES.ub'] < 0) signegBAYES_RAW <- signegBAYES_RAW + 1
			if (results_BAYES_RAW[luper,'ES.lb'] > 0 & results_BAYES_RAW[luper,'ES.ub'] > 0) sigposBAYES_RAW <- sigposBAYES_RAW + 1
			if (results_BAYES_RAW[luper,'ES.lb'] < 0 & results_BAYES_RAW[luper,'ES.ub'] > 0) nonsigBAYES_RAW <- nonsigBAYES_RAW + 1
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

	message('\n\nFinal, all-studies-combined results:')
	
	message('\n   Cumulative Meta-Analysis:    ES = ', round(results_CUM_META[,'ES'][Nstudies],3),
	        '   ES.lb = ', round(results_CUM_META[,'ES.lb'][Nstudies],3),
	        '   ES.ub = ', round(results_CUM_META[,'ES.ub'][Nstudies],3))
	
if (!is.null(donnes_RN)) {

	if (is.element('Schmidt_Raju', Bayes_type)) 
		message('\n   Bayesian (Schmidt-Raju):     ES = ', round(results_BAYES_SR$ES[Nstudies],3),
		        '   ES.lb = ', round(results_BAYES_SR$ES.lb[Nstudies],3),
		        '   ES.ub = ', round(results_BAYES_SR$ES.ub[Nstudies],3))
	
	if (is.element('generated', Bayes_type)) 
		message('\n   Bayesian (generated data):   ES = ', round(results_BAYES_GEN$ES[Nstudies],3),
		        '   ES.lb = ', round(results_BAYES_GEN$ES.lb[Nstudies],3),
		        '   ES.ub = ', round(results_BAYES_GEN$ES.ub[Nstudies],3))
	
	if (is.element('raw', Bayes_type) & results_BAYES_RAW != -9999) 
		message('\n   Bayesian (raw data):         ES = ', round(results_BAYES_RAW$ES[Nstudies],3),
		        '   ES.lb = ', round(results_BAYES_RAW$ES.lb[Nstudies],3),
		        '   ES.ub = ', round(results_BAYES_RAW$ES.ub[Nstudies],3))

}	
	
	if (rma_method != 'FE') {
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
	ES_MA_lb = results_CUM_META[,'ES.lb'][Nstudies],  
	ES_MA_ub = results_CUM_META[,'ES.ub'][Nstudies],
	Q = outp_MA_1$QE, p_Q = outp_MA_1$QEp,
	tau2 = tau2, tau2LB = tau2LB, tau2UB = tau2UB,
	tau = tau, tauLB = tauLB, tauUB = tauUB,
	isq = isq, isqLB = isqLB, isqUB = isqUB,
	hsq = hsq, hsqLB = hsqLB, hsqUB = hsqUB,
	consist_NHST = consist_NHST, agree_NHST = agree_NHST,
	consist_CUM_META = consist_CUM_META, agree_CUM_META = agree_CUM_META,
	biasStats = biasStats )

if (!is.null(donnes_RN)) {

	if (is.element('Schmidt_Raju', Bayes_type))
	    nopingpongOutput <- c(nopingpongOutput, list(results_BAYES_SR = results_BAYES_SR, 
		consist_BAYES_SR = consist_BAYES_SR, agree_BAYES_SR = agree_BAYES_SR) )
	
	if (is.element('generated', Bayes_type))  
	    nopingpongOutput <- c(nopingpongOutput, list(results_BAYES_GEN = results_BAYES_GEN,
		consist_BAYES_GEN = consist_BAYES_GEN, agree_BAYES_GEN = agree_BAYES_GEN ) )
			
	if (is.element('raw', Bayes_type) & results_BAYES_RAW != -9999)
	    nopingpongOutput <- c(nopingpongOutput, list(results_BAYES_RAW = results_BAYES_RAW,
		consist_BAYES_RAW = consist_BAYES_RAW, agree_BAYES_RAW = agree_BAYES_RAW) )

}

class(nopingpongOutput) <- "NO.PING.PONG"

return(invisible(nopingpongOutput))

}




