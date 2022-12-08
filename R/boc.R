





ES_from_GRPINFO <- function(grp1_MN = NULL, grp1_SD = NULL, grp1_N = NULL, 
                            grp2_MN = NULL, grp2_SD = NULL, grp2_N = NULL, 
                            ES_type_OUT='g', gvar_type_OUT='d') {

	totalN <- grp1_N + grp2_N

	SDpooled <- sqrt((grp1_SD**2 * (grp1_N - 1) + grp2_SD**2 * (grp2_N - 1))/(grp1_N + grp2_N - 2))
	
	dvalues <- (grp1_MN - grp2_MN) / SDpooled

	# compute variance of d -- Effect Sizes Based on Means.pdf 
	Vd <- (grp1_N + grp2_N) / (grp1_N * grp2_N) + ((dvalues**2) / (2 * (grp1_N + grp2_N)))

	if (ES_type_OUT == 'd')  ES_from_GRPINFOOutput <- data.frame(d=dvalues, Vd=Vd, totalN)

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

		ES_from_GRPINFOOutput <- data.frame(g=gvalues, Vg=Vg, totalN)
	}

	if (ES_type_OUT == 'r') {
		# r from d -- Meta-analysis Converting among effect sizes.pdf 
		aa <- (grp1_N + grp2_N)**2 / (grp1_N * grp2_N)
		rfromd <- dvalues / sqrt(dvalues**2 + aa)
		Vr <- (aa**2 * Vd) / (dvalues**2 + aa)^3

		ES_from_GRPINFOOutput <- data.frame(r=rfromd, Vr=Vr, totalN)
	}

return(invisible(ES_from_GRPINFOOutput))	

#print(ES_from_GRPINFOOutput)
}




