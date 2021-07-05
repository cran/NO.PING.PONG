




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






# # d from r - using r, z and Vr/Vd to compute d
# dfromrVrz <- function(r, Vr, z) {

	# Vd <- (4 * Vr) / (1 - r**2)**3  # Meta-analysis Converting among effect sizes.pdf
	
	# SEd <- sqrt(Vd)
		
	# dvalues <- SEd * z
	
	# print(cbind(dvalues, Vd))
	
	# return(invisible(dvalues))	
# }
               
# dfromrVrz(r = outp_CUM$estimate, Vr = outp_CUM$se^2, z = outp_CUM$zval)

# dd <- esc_rpb(r = outp_CUM$estimate, totaln=totalNs, es.type='g')
# dat <- as.matrix(cbind( as.numeric(dd[[1]]),  as.numeric(dd[[3]]),  as.numeric(dd[[2]]) ))
# colnames(dat) <- c('g', 'var', 'se'); round(dat,4)




# # 
            
# # d or g from r
# gdFROMr <- function(r, Vr = NULL, dg_type_OUT = 'g', gvar_type_OUT = 'd', totalN) {

	# # d & Vd
	# dvalues <- (2 * r) / sqrt(1 - (r**2))  # Meta-analysis Converting among effect sizes.pdf

	# Vd <- Vg <- gvalues <- rep(NA, length(r))   # NULL # may remain NULL 

	# if (!is.null(Vr))	
		# Vd <- (4 * Vr) / (1 - r**2)**3  # Meta-analysis Converting among effect sizes.pdf

	# if (dg_type_OUT == 'd')  gdFROMrOutput <- data.frame(d=dvalues, Vd, totalN)

	# # g & Vg, if g are requested -- from R compute.es.pdf -- p 7
	# if (dg_type_OUT == 'g') {
		# df <- totalN - 2
		# bigJ <- 1 - (3 / (4 * df - 1))  
		# gvalues <- bigJ * dvalues   

		# if(!is.null(Vd)) {  
			# if (gvar_type_OUT == 'g') Vg <- bigJ**2 * Vd  
			# if (gvar_type_OUT == 'd') Vg <- Vd
		# }

		# gdFROMrOutput <- data.frame(g=gvalues, Vg, totalN)
	# }

	# return(invisible(gdFROMrOutput))	
# } 
 
 
 



# # r (& Vr) from g or d  
# rVrFROMgd <- function(ES, Vd = NULL, dg_type = 'g', totalN = NULL, grp1_N = NULL, grp2_N = NULL) {
                     	
	# # are group Ns known?
	# if (!is.null(grp1_N) & !is.null(grp2_N)) {grpNs = TRUE} else {grpNs = FALSE}
	                     	
	# # convert g to d -- from R compute.es.pdf
	# if (dg_type == 'g') {
		# if (!is.null(totalN)) {
			# df <- totalN - 2
			# bigJ <- 1 - (3 / (4 * df - 1))
			# dvalues <- ES / bigJ
		# }
		# if (grpNs == TRUE) {
			# df <- (grp1_N + grp2_N) - 2
			# bigJ <- 1 - (3 / (4 * df - 1))
			# dvalues <- ES / bigJ			
		# }
	# } else {dvalues <- ES}

	# # compute variance of d -- derived from the (hidden) esc.vd function in the esc package
	# if (is.null(Vd) & grpNs == TRUE)
		# Vd <- (grp1_N + grp2_N) / (grp1_N * grp2_N) + (dvalues * dvalues) / (2 * (grp1_N + grp2_N))

	# Vr <- rfromd <- NULL  # rep(NA, length(dvalues))   # will remain NULL if it cannot be computed
	
	# if (grpNs == TRUE & !is.null(Vd)) {	
		# # r from d -- Meta-analysis Converting among effect sizes.pdf -- also provides Vr
		# aa <- (grp1_N + grp2_N)**2 / (grp1_N * grp2_N)
		# rfromd <- dvalues / sqrt(dvalues**2 + aa)
		# Vr <- (aa**2 * Vd) / (dvalues**2 + aa)^3
	# }                 	
	# # if totalN is available & group Ns are not known, use n1=n2, which will yield a=4 
	# # Meta-analysis Converting among effect sizes.pdf
	# # also in R compute.es.pdf p 7
	# # gets pretty close, even for quite unequal Ns
	# if (grpNs == FALSE & !is.null(totalN) & !is.null(Vd)) {    	
		# aa <- 4
		# rfromd <- dvalues / sqrt(dvalues**2 + aa)
		# Vr <- (aa**2 * Vd) / (dvalues**2 + aa)^3
	# }	

	# if (grpNs == TRUE & is.null(Vd) == TRUE) {    	
		# aa <- (grp1_N + grp2_N)**2 / (grp1_N * grp2_N)
		# rfromd <- dvalues / sqrt(dvalues**2 + aa)
	# }	

	# if (grpNs == FALSE & !is.null(totalN) & is.null(Vd)) {    	
		# aa <- 4
		# rfromd <- dvalues / sqrt(dvalues**2 + aa)
	# }	

    # # if (is.na(all(totalN)))  totalN <- rep(NA, length(dvalues))             	
    # # if (is.na(all(Vr)))      Vr     <- rep(NA, length(dvalues))             	

    # if (is.null(totalN))   totalN <- rep(NA, length(dvalues))             	
    # if (is.null(Vr))       Vr     <- rep(NA, length(dvalues))             	
	# if (grpNs == TRUE)     totalN <- grp1_N + grp2_N
                       	
	# rVrFROMgdOutput <- data.frame(r=rfromd, Vr, totalN)

	# # rVrFROMgdOutput <- data.frame(rfromgd=rfromd, Vr=Vr, totalN=totalN)

	# # rVrFROMgdOutput <- list(rfromgd=rfromd, Vr=Vr, totalN=totalN)
	
	# return(invisible(rVrFROMgdOutput))
# }








# # 
            
# CONVERT_ES <- function(ES, ES_var = NULL, ES_type_IN='r', ES_type_OUT='g', 
                      # totalN = NULL, grp1_N = NULL, grp2_N = NULL,
                      # gvar_type_OUT = 'd') {

# # are group Ns known?
# if (!is.null(grp1_N) & !is.null(grp2_N)) {grpNs = TRUE} else {grpNs = FALSE}

# if (is.null(totalN) & grpNs) totalN <- grp1_N + grp2_N

# if (ES_type_IN == 'r') {
	
	# r <- ES
	
	# if (ES_type_OUT == 'd' | ES_type_OUT == 'g') {
		
		# d <- (2 * r) / sqrt(1 - (r^2))  # Meta-analysis Converting among effect sizes.pdf
	
		# Vd <- Vg <- rep(NA, length(r))   # may remain NA 
	
		# if (!is.null(ES_var))	
			# Vd <- (4 * ES_var) / (1 - r^2)^3  # Meta-analysis Converting among effect sizes.pdf
	
		# if (ES_type_OUT == 'd')  Output <- data.frame(d=d, Vd=Vd, totalN)				
	
		# if (ES_type_OUT == 'g') {  # from R compute.es.pdf -- p 7
			
			# df <- totalN - 2			
			# bigJ <- 1 - (3 / (4 * df - 1))			  
			# g <- bigJ * d   
	
			# if(!is.null(Vd)) {  
				# if (gvar_type_OUT == 'g') Vg <- bigJ^2 * Vd  
				# if (gvar_type_OUT == 'd') Vg <- Vd
			# }

			# Output <- data.frame(g=g, Vg=Vg, totalN)			
		# }
	# }
# }


# if (ES_type_IN == 'd' | ES_type_IN == 'g') {

	# # compute Vd if NULL & grp Ns are available 
	# if (ES_type_IN == 'd' & is.null(ES_var) & grpNs & ES_type_OUT == 'd') {
		# # derived from the (hidden) esc.vd function in the esc package
		# Vd <- (grp1_N + grp2_N) / (grp1_N * grp2_N) + (ES * ES) / (2 * (grp1_N + grp2_N))
	# }
		
	# # convert g to d -- from R compute.es.pdf
	# if (ES_type_IN == 'g' & ES_type_OUT == 'd') {

		# df <- totalN - 2
		# bigJ <- 1 - (3 / (4 * df - 1))
		# d <- ES / bigJ

		# if (!is.null(ES_var)) {  
			# if (gvar_type_OUT == 'd')  Vd <- ES_var / bigJ^2
			# if (gvar_type_OUT == 'g')  Vd <- ES_var
		# }

		# Output <- data.frame(d=d, Vd=Vd, totalN)
	# }

	# # convert d to g -- from R compute.es.pdf
	# if (ES_type_IN == 'd' & ES_type_OUT == 'g') {

		# df <- totalN - 2			
		# bigJ <- 1 - (3 / (4 * df - 1))			  
		# g <- bigJ * ES   

		# Vg <- rep(NA, length(ES))   # may remain NA 
		# if(!is.null(Vd)) {  
			# if (gvar_type_OUT == 'g')  Vg <- bigJ^2 * Vd  
			# if (gvar_type_OUT == 'd')  Vg <- Vd
		# }

		# Output <- data.frame(g=g, Vg=Vg, totalN)			
	# }



	# if (ES_type_OUT == 'r') {
	
		# if (is.null(ES_var))  Vd = NULL 
	
		
		# # r will be obtained from d, so if ES = g, need to convert g to d -- from R compute.es.pdf
		# if (ES_type_IN == 'g') {
	
			# df <- totalN - 2
			# bigJ <- 1 - (3 / (4 * df - 1))
			# d <- ES / bigJ
	
			# if (!is.null(ES_var)) {  
				# if (gvar_type_OUT == 'd')  Vd <- ES_var / bigJ^2
				# if (gvar_type_OUT == 'g')  Vd <- ES_var
			# }
		# }
	
		# if (ES_type_IN == 'd')  d <- ES
	
		# if (ES_type_IN == 'd' & !is.null(ES_var))  Vd <- ES_var
	
	
		# Vr <- r <- rep(NA, length(ES))   # may remain NA
		
		# # if have the grp Ns & Vd
		# if (grpNs & !is.null(Vd)) {	
			# # r from d -- Meta-analysis Converting among effect sizes.pdf -- also provides Vr
			# aa <- (grp1_N + grp2_N)^2 / (grp1_N * grp2_N)
			# r <- d / sqrt(d^2 + aa)
			# Vr <- (aa^2 * Vd) / (d^2 + aa)^3
		# }
		                 	
		# # if have Vd, but not group Ns, use n1=n2, which will yield a=4 
		# # Meta-analysis Converting among effect sizes.pdf
		# # also in R compute.es.pdf p 7
		# # gets pretty close, even for quite unequal Ns
		# if (!is.null(Vd) & !grpNs) {    	
			# aa <- 4
			# r <- d / sqrt(d^2 + aa)
			# Vr <- (aa^2 * Vd) / (d^2 + aa)^3
		# }	
	
		# # have grp Ns but not Vd
		# if (grpNs & is.null(Vd)) {    	
			# aa <- (grp1_N + grp2_N)^2 / (grp1_N * grp2_N)
			# r <- d / sqrt(d^2 + aa)
		# }	
	
		# # do not have grp Ns or Vd
		# if (!grpNs & is.null(Vd)) {    	
			# aa <- 4
			# r <- d / sqrt(d^2 + aa)
		# }	
	
	
	    # if (is.null(totalN))   totalN <- rep(NA, length(ES))             	
	    # if (is.null(Vr))       Vr     <- rep(NA, length(ES))             	
	                       	
		# Output <- data.frame(r=r, Vr=Vr, totalN=totalN)
	
	# }
# }

# return(invisible(Output)) 
# }


