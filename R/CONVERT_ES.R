

            
CONVERT_ES <- function(ES, ES_var = NULL, ES_type_IN='r', ES_type_OUT='g', 
                      totalN = NULL, grp1_N = NULL, grp2_N = NULL,
                      gvar_type_OUT = 'd', verbose = TRUE) {

# are group Ns known?
if (!is.null(grp1_N) & !is.null(grp2_N)) {grpNs = TRUE} else {grpNs = FALSE}

if (is.null(totalN) & grpNs) totalN <- grp1_N + grp2_N

if (is.null(totalN)) totalN <- as.numeric(rep(NA, length(ES)))


if (ES_type_IN == 'r') {
	
	r <- ES
	
	# setting r values > 1 to NA
	r <- ifelse( r >  1,   .9999, r)
	r <- ifelse( r < -1,  -.9999, r)
	
	if (ES_type_OUT == 'd' | ES_type_OUT == 'g') {
		
		d <- (2 * r) / sqrt(1 - (r^2))  # Meta-analysis Converting among effect sizes.pdf
	
		Vd <- Vg <- as.numeric(rep(NA, length(r)))   # may remain NA 
	
		if (!is.null(ES_var))	
			Vd <- (4 * ES_var) / (1 - r^2)^3  # Meta-analysis Converting among effect sizes.pdf
	
		if (ES_type_OUT == 'd')  Output <- data.frame(d=d, Vd=Vd, totalN)				
	
		if (ES_type_OUT == 'g') {  # from R compute.es.pdf -- p 7

			if (is.null(totalN)) {
				message('\n\ng values cannot be computed because totalN is NULL.')
				} else {
			
				df <- totalN - 2			
				bigJ <- 1 - (3 / (4 * df - 1))			  
				g <- bigJ * d   
		
				if(!is.null(Vd)) {  
					if (gvar_type_OUT == 'g') Vg <- bigJ^2 * Vd  
					if (gvar_type_OUT == 'd') Vg <- Vd
				}
	
				Output <- data.frame(g=g, Vg=Vg, totalN)			
			}
		}
	}
}


if (ES_type_IN == 'd' | ES_type_IN == 'g') {

	if (!is.null(ES_var))  {Vd <- Vg <- ES_var} else {Vd <- Vg <- NULL}

	# compute Vd if NULL & grp Ns are available 
	if (ES_type_IN == 'd' & is.null(Vd) & grpNs) {
		# derived from the (hidden) esc.vd function in the esc package
		Vd <- (grp1_N + grp2_N) / (grp1_N * grp2_N) + (ES * ES) / (2 * (grp1_N + grp2_N))
	}
		
	# convert g to d -- from R compute.es.pdf
	if (ES_type_IN == 'g' & ES_type_OUT == 'd') {

		df <- totalN - 2
		bigJ <- 1 - (3 / (4 * df - 1))
		d <- ES / bigJ

		if (!is.null(ES_var)) {  
			if (gvar_type_OUT == 'd')  Vd <- ES_var / bigJ^2
			if (gvar_type_OUT == 'g')  Vd <- ES_var
		}

		Output <- data.frame(d=d, Vd=Vd, totalN)
	}

	# convert d to g -- from R compute.es.pdf
	if (ES_type_IN == 'd' & ES_type_OUT == 'g') {

		df <- totalN - 2			
		bigJ <- 1 - (3 / (4 * df - 1))			  
		g <- bigJ * ES   

		Vg <- as.numeric(rep(NA, length(ES)))   # may remain NA 
		if(!is.null(Vd)) {  
			if (gvar_type_OUT == 'g')  Vg <- bigJ^2 * Vd  
			if (gvar_type_OUT == 'd')  Vg <- Vd
		}

		Output <- data.frame(g=g, Vg=Vg, totalN)			
	}



	if (ES_type_OUT == 'r') {
	
		# if (is.null(ES_var))  Vd = NULL 
	
		
		# r will be obtained from d, so if ES = g, need to convert g to d -- from R compute.es.pdf
		if (ES_type_IN == 'g') {
	
			df <- totalN - 2
			bigJ <- 1 - (3 / (4 * df - 1))
			d <- ES / bigJ
	
			if (!is.null(ES_var)) {  
				if (gvar_type_OUT == 'd')  Vd <- ES_var / bigJ^2
				if (gvar_type_OUT == 'g')  Vd <- ES_var
			}
		}
	
		if (ES_type_IN == 'd')  d <- ES
	
		if (ES_type_IN == 'd' & !is.null(ES_var))  Vd <- ES_var
	
	
		Vr <- r <- as.numeric(rep(NA, length(ES)))   # may remain NA
		
		# if have the grp Ns & Vd
		if (grpNs & !is.null(Vd)) {	
			# r from d -- Meta-analysis Converting among effect sizes.pdf -- also provides Vr
			aa <- (grp1_N + grp2_N)^2 / (grp1_N * grp2_N)
			r <- d / sqrt(d^2 + aa)
			Vr <- (aa^2 * Vd) / (d^2 + aa)^3
		}
		                 	
		# if have Vd, but not group Ns, use n1=n2, which will yield a=4 
		# Meta-analysis Converting among effect sizes.pdf
		# also in R compute.es.pdf p 7
		# gets pretty close, even for quite unequal Ns
		if (!is.null(Vd) & !grpNs) {    	
			aa <- 4
			r <- d / sqrt(d^2 + aa)
			Vr <- (aa^2 * Vd) / (d^2 + aa)^3
		}	
	
		# have grp Ns but not Vd
		if (grpNs & is.null(Vd)) {    	
			aa <- (grp1_N + grp2_N)^2 / (grp1_N * grp2_N)
			r <- d / sqrt(d^2 + aa)
		}	
	
		# do not have grp Ns or Vd
		if (!grpNs & is.null(Vd)) {    	
			aa <- 4
			r <- d / sqrt(d^2 + aa)
		}	
	
	
	    if (is.null(totalN))   totalN <- as.numeric(rep(NA, length(ES)))            	
	    if (is.null(Vr))       Vr     <- as.numeric(rep(NA, length(ES)))            	
	                       	
		Output <- data.frame(r=r, Vr=Vr, totalN=totalN)
	
	}
}

if (verbose) {
	message('\n\n Effect sizes:\n')
	print(round(Output,3))
}

class(Output) <- c("data.frame")

return(invisible(Output)) 
}



