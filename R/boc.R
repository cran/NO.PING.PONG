




r_from_z <- function (z, z_var = NULL, totalN = NULL, zforCI_out=1.96) {

	r_var <- r_lb <- r_ub <- as.numeric(rep(NA, length(z))) 

	r <- tanh(z)

	if (!is.null(totalN))  r_var <- (1 - r**2)**2 / (totalN - 1)  # 2009 Borenstein - Introduction to Meta-Analysis p 41
	
	if (is.null(z_var) & !is.null(totalN))  z_var <- 1 / (totalN - 3)
	
	if (!is.null(z_var)) {
		
		z_lb <- z - zforCI_out * sqrt(z_var)
		
		z_ub <- z + zforCI_out * sqrt(z_var)
	
		r_lb <- tanh(z_lb)
	
		r_ub <- tanh(z_ub)
	}

	if (is.null(totalN)) totalN <- as.numeric(rep(NA, length(z)))
	
	Output <- data.frame(r=r, r_var=r_var, r_lb=r_lb, r_ub=r_ub, totalN=totalN)	

	return(invisible(Output)) 
}






z_from_r <- function (r, totalN = NULL, zforCI_out=1.96) {

	z_var <- z_lb <- z_ub <- as.numeric(rep(NA, length(r))) 

	z <- atanh(r)     # z <-  0.5*log((1 + r)/(1-r))  # r to z' 

	if (!is.null(totalN)) {
		
		z_var <- 1 / (totalN - 3)    

		z_lb <- z - zforCI_out * sqrt(z_var)
		
		z_ub <- z + zforCI_out * sqrt(z_var)
	}

	if (is.null(totalN))  totalN <- as.numeric(rep(NA, length(r)))
	
	Output <- data.frame(z=z, z_var=z_var, z_lb=z_lb, z_ub=z_ub, totalN=totalN)	

	return(invisible(Output)) 
}






r_from_d <- function (d, d_var = NULL, totalN = NULL, grp1_n = NULL, grp2_n = NULL, zforCI_out=1.96) {

	r_var <- r_lb <- r_ub <- as.numeric(rep(NA, length(d))) 

	# are group Ns known?
	if (!is.null(grp1_n) & !is.null(grp2_n))  {grpNs = TRUE} else {grpNs = FALSE}
	
	if (!is.null(totalN) & grpNs)  totalN <- grp1_n + grp2_n
	

	# have the grp Ns & d_var
	if (grpNs & !is.null(d_var)) {	
		# r from d -- Meta-analysis Converting among effect sizes.pdf -- also provides r_var
		# aa corrects for unequal Ns
		aa <- (grp1_n + grp2_n)^2 / (grp1_n * grp2_n)
		r <- d / sqrt(d^2 + aa)
		r_var <- (aa^2 * d_var) / (d^2 + aa)^3
	}
	                 	
	# have d_var, but not group Ns, use n1=n2, which will yield a=4 
	# Meta-analysis Converting among effect sizes.pdf
	# also in R compute.es.pdf p 7
	# gets pretty close, even for quite unequal Ns
	if (!is.null(d_var) & !grpNs) {    	
		aa <- 4
		r <- d / sqrt(d^2 + aa)
		r_var <- (aa^2 * d_var) / (d^2 + aa)^3
	}	

	# have grp Ns but not d_var
	if (grpNs & is.null(d_var)) {    	
		aa <- (grp1_n + grp2_n)^2 / (grp1_n * grp2_n)
		r <- d / sqrt(d^2 + aa)
	}	

	# do not have grp Ns or d_var
	if (!grpNs & is.null(d_var)) {    	
		aa <- 4
		# r <- d / sqrt(d^2 + aa)
		r <- sqrt( d^2 / (d^2 + 4) )		 
	}	


	# r CIs - via r-to-z    -- need totalN
	if (!is.null(totalN)) {

		z <- atanh(r)  # r to z	    # z <-  0.5*log((1 + r)/(1-r))  # r to z' 
		
		z_var <- 1 / (totalN - 3)    

		z_lb <- z - zforCI_out * sqrt(z_var)
		
		z_ub <- z + zforCI_out * sqrt(z_var)
	
		r_lb <- tanh(z_lb)  # z to r
	
		r_ub <- tanh(z_ub)  # z to r
	}

	if (is.null(totalN)) totalN <- as.numeric(rep(NA, length(d)))
	
	Output <- data.frame(r=r, r_var=r_var, r_lb=r_lb, r_ub=r_ub, totalN=totalN)	

	return(invisible(Output)) 
}





d_from_r <- function (r, r_var = NULL, grp1_n = NULL, grp2_n = NULL, zforCI_out=1.96) {

	d <- (2 * r) / sqrt(1 - r^2)  # Meta-analysis Converting among effect sizes (2009 Borenstein).pdf

	d_var <- d_lb <- d_ub <- NULL  # NA 

	# for d_var, best to use the group Ns (if available) than the r_var conversion formula
	if (!is.null(grp1_n) & !is.null(grp2_n))  
		# Effect Sizes Based on Means (2009 Borenstein) 4.20
		d_var <- ( (grp1_n + grp2_n) /  (grp1_n * grp2_n) ) + ( d**2 / (2 * (grp1_n + grp2_n)) )

	# have r_var but not the group Ns (use the r_var conversion formula) 
	# Meta-analysis Converting among effect sizes.pdf
	if (!is.null(r_var) & (is.null(grp1_n) & is.null(grp2_n)))	
		d_var <- (4 * r_var) / (1 - r^2)^3  # Meta-analysis Converting among effect sizes.pdf

	if (!is.null(d_var)) {
		d_lb <- d - zforCI_out * sqrt(d_var)
		d_ub <- d + zforCI_out * sqrt(d_var)
	}

	if (is.null(d_var))  d_var <- d_lb <- d_ub <- as.numeric(rep(NA, length(d)))
	
	Output <- data.frame(d=as.numeric(d), d_var=d_var, d_lb=d_lb, d_ub=d_ub)

	return(invisible(Output)) 
}





d_from_g <- function (g, g_var = NULL, totalN = NULL, grp1_n = NULL, grp2_n = NULL, zforCI_out=1.96) {

	# need either totalN, or grp1_n & grp2_n
	
	d <- d_var <- d_lb <- d_ub <- NULL   # NA 

	if (is.null(totalN) & !is.null(grp1_n) & !is.null(grp2_n))  totalN <- grp1_n + grp2_n

	if (!is.null(totalN)) {
		df <- totalN - 2
		bigJ <- 1 - (3 / (4 * df - 1))
		d <- g / bigJ
	
		if (!is.null(g_var))  d_var <- g_var / bigJ^2

		# no g_var, but have group Ns
		if (is.null(g_var) & !is.null(grp1_n) & !is.null(grp2_n))  
			# Effect Sizes Based on Means (2009 Borenstein) 4.20
			d_var <- ( (grp1_n + grp2_n) /  (grp1_n * grp2_n) ) + ( d**2 / (2 * (grp1_n + grp2_n)) )
	}

	if (!is.null(d_var)) {
		d_lb <- d - zforCI_out * sqrt(d_var)
		d_ub <- d + zforCI_out * sqrt(d_var)
	}

	if (is.null(d_var))  d <- d_var <- d_lb <- d_ub <- as.numeric(rep(NA, length(g)))

	Output <- data.frame(d=d, d_var=d_var, d_lb=d_lb, d_ub=d_ub)

	return(invisible(Output)) 
}




r_to_d <- function(r) { d <- (2 * r) / sqrt(1 - (r^2)); return(invisible(d)) }  # Meta-analysis Converting among effect sizes.pdf



d_to_g <- function(d, N) {
		df <- N - 2			
		bigJ <- 1 - (3 / (4 * df - 1))			  
		g <- bigJ * d   
		return(invisible(g)) 
}




g_from_d <- function (d, d_var = NULL, totalN, grp1_n = NULL, grp2_n = NULL, zforCI_out=1.96, 
                      gvar_type_OUT = 'g') {

	# compute.es.pdf -- p 7
	
	g <- g_var <- g_lb <- g_ub <- NULL 

	# no d_var, but have group Ns
	if (is.null(d_var) & !is.null(grp1_n) & !is.null(grp2_n))  
		# Effect Sizes Based on Means (2009 Borenstein) 4.20
		d_var <- ( (grp1_n + grp2_n) /  (grp1_n * grp2_n) ) + ( d**2 / (2 * (grp1_n + grp2_n)) )

	# no totalN, but have group Ns
	if (is.null(totalN) & !is.null(grp1_n) & !is.null(grp2_n))  
		totalN <- grp1_n + grp2_n

	if (is.null(totalN)) {
		# message('\n\ng values cannot be computed because the total N is not available.')
	} else {
		df <- totalN - 2			
		bigJ <- 1 - (3 / (4 * df - 1))			  
		g <- bigJ * d   
		names(g) <- 'g'	
	
		if (!is.null(d_var)) {
			if (gvar_type_OUT == 'g')  g_var <- bigJ^2 * d_var  
			if (gvar_type_OUT == 'd')  g_var <- d_var
		}
	}

	if (!is.null(g_var)) {
		g_lb <- g - zforCI_out * sqrt(g_var)
		g_ub <- g + zforCI_out * sqrt(g_var)
	}

	if (is.null(g_var))  g_var <- g_lb <- g_ub <- as.numeric(rep(NA, length(d)))

	Output <- data.frame(g=g, g_var=g_var, g_lb=g_lb, g_ub=g_ub)

	return(invisible(Output)) 
}





OR_from_d <- function (d, d_var = NULL, grp1_n = NULL, grp2_n = NULL, zforCI_out=1.96) {

	# *** Meta-analysis Converting among effect sizes (2009 Borenstein)  pp 3-4

	OR <- OR_lb <- OR_ub <- logOR_var <- logOR_lb <- logOR_ub   <-  rep(NA, length(d))
	
	logOR <- d * (pi / sqrt(3))
	
	# d_var if NULL & grp Ns are available 
	if (is.null(d_var) & !is.null(grp1_n) & !is.null(grp2_n)) {
		# derived from the (hidden) esc.vd function in the esc package
		d_var <- (grp1_n + grp2_n) / (grp1_n * grp2_n) + (d * d) / (2 * (grp1_n + grp2_n))
	}
 
	if (!is.null(d_var))  logOR_var <- d_var * (pi**2 / 3)

	# CIs for logOR
	if (!all(is.na(logOR_var))) {  
		logOR_lb <- logOR - zforCI_out * sqrt(logOR_var)
		logOR_ub <- logOR + zforCI_out * sqrt(logOR_var)
	}

	# logOR to OR
	OR <- exp(logOR)
	if (!all(is.na(logOR_lb))) {  	
		OR_lb <- exp(logOR_lb)
		OR_ub <- exp(logOR_ub)
	}
	
	Output <- data.frame(OR=OR, OR_lb=OR_lb, OR_ub=OR_ub,
	                     logOR=logOR, logOR_var=logOR_var, logOR_lb=logOR_lb, logOR_ub=logOR_ub)

	return(invisible(Output)) 
}










ES_from_GRPINFO <- function(grp1_MN = NULL, grp1_SD = NULL, grp1_n = NULL, 
                            grp2_MN = NULL, grp2_SD = NULL, grp2_n = NULL, 
                            ES_type_OUT='g', gvar_type_OUT='d') {

	totalN <- grp1_n + grp2_n

	SDpooled <- sqrt((grp1_SD**2 * (grp1_n - 1) + grp2_SD**2 * (grp2_n - 1))/(grp1_n + grp2_n - 2))
	
	dvalues <- (grp1_MN - grp2_MN) / SDpooled

	# compute variance of d -- Effect Sizes Based on Means.pdf 
	Vd <- (grp1_n + grp2_n) / (grp1_n * grp2_n) + ((dvalues**2) / (2 * (grp1_n + grp2_n)))

	if (ES_type_OUT == 'd')  Output <- data.frame(d=dvalues, Vd=Vd, totalN)

	# g & Vg, if requested
	if (ES_type_OUT == 'g') {

		# after lots of time, I learned that:
		# the escalc function (metafor), the esc_mean_sd function (esc), 
		# the mes function (compute.es), and my commands all produce identical g values
		
		# but the variance of g values from escalc & esc_mean_sd are not the same as those
		# from mes. I discovered that mes computes the variance of g using
		# bigJ**2 * Vd, which is what Hedges & the texts/pdfs say to do.
		# but escalc & esc_mean_sd set Vg = Vd

		# I therefore make it an option, via gvar_type_OUT='g'

		# below, gets the exact gvalues but the Vg values are just slightly off
		# g & Vg, if requested -- from R compute.es.pdf -- p 7
		df <- totalN - 2
		bigJ <- 1 - (3 / (4 * df - 1))  
		gvalues <- bigJ * dvalues   

		if (gvar_type_OUT == 'g') Vg <- bigJ**2 * Vd  
		if (gvar_type_OUT == 'd') Vg <- Vd

		Output <- data.frame(g=gvalues, Vg=Vg, totalN)
	}

	if (ES_type_OUT == 'r') {
		# r from d -- Meta-analysis Converting among effect sizes.pdf 
		aa <- (grp1_n + grp2_n)**2 / (grp1_n * grp2_n)
		rfromd <- dvalues / sqrt(dvalues**2 + aa)
		Vr <- (aa**2 * Vd) / (dvalues**2 + aa)^3

		Output <- data.frame(r=rfromd, Vr=Vr, totalN)
	}

return(invisible(Output))	

#print(Output)
}




