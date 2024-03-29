

# # donnes = data_nPP$Omega3_Depression; ES = 'SMD'; ES_var = NULL; ES_type_IN='d'; 
           # grp1_n = 'CN'; grp2_n = 'EN';
            
                       # totalN = NULL; 
                       # gvar_type_OUT = 'd'; 
                       # CI_level_out = 95; CI_level_in = 95; CI_in_lb = NULL; CI_in_ub = NULL; 
                       # verbose = TRUE      
                       
                       
# donnes = data_nPP$Math_Performance; ES = 'r'; ES_var = NULL; ES_type_IN='r';
           # totalN = 'N';
           # grp1_n = NULL; grp2_n = NULL;
                       # gvar_type_OUT = 'd'; 
                       # CI_level_out = 95; CI_level_in = 95; CI_in_lb = NULL; CI_in_ub = NULL; 
                       # verbose = TRUE
           
           
           
            
CONVERT_ES <- function(donnes, ES, ES_type_IN='r', ES_var = NULL, 
                       totalN = NULL, grp1_n = NULL, grp2_n = NULL,
                       gvar_type_OUT = 'd', 
                       CI_level_out = 95, CI_level_in = 95, CI_in_lb = NULL, CI_in_ub = NULL, 
                       verbose = TRUE) {

                       ES     <- donnes[,ES]
if (!is.null(ES_var))  ES_var <- donnes[,ES_var]
if (!is.null(totalN))  totalN <- donnes[,totalN]
if (!is.null(grp1_n))  grp1_n <- donnes[,grp1_n]
if (!is.null(grp2_n))  grp2_n <- donnes[,grp2_n]

# using one of CI_in_lb or CI_in_ub
ci_dat <- NULL
if (!is.null(CI_in_lb))                      ci_dat <- donnes[,CI_in_lb]
if (is.null(CI_in_lb) & !is.null(CI_in_ub))  ci_dat <- donnes[,CI_in_ub]


# are group Ns known?
if (!is.null(grp1_n) & !is.null(grp2_n))  {grpNs = TRUE} else {grpNs = FALSE}

if (is.null(totalN) & grpNs)  totalN <- grp1_n + grp2_n


# the critical z value (that corresponds to the specified CI) to be used in the CI computations
zforCI_out <- qnorm((1 + CI_level_out * .01) / 2) 
zforCI_in  <- qnorm((1 + CI_level_in  * .01) / 2) 


# initializing output variables
r <- r_var <- r_lb <- r_ub <-
d <- d_var <- d_lb <- d_ub <-
g <- g_var <- g_lb <- g_ub <-
z <- z_var <- z_lb <- z_ub <-
OR <- OR_lb <- OR_ub <-
logOR <- logOR_var <- logOR_lb <- logOR_ub    <-  rep(NA, nrow(donnes))



if (ES_type_IN == 'r') {
	
	r  <- ES
	
	# # setting r values > 1 to NA
	# r <- ifelse( r >  1,  .9999, r)
	# r <- ifelse( r < -1, -.9999, r)

	# r_var			
	if (!is.null(ES_var))  r_var <- ES_var
	if (is.null(ES_var) & !is.null(totalN))  r_var <- escalc(measure='ZCOR', ni=totalN, ri=r )$vi

	# r CIs - via r-to-z then z-to-r -- need totalN	
	if (!is.null(totalN)) {
	
		z_stats <- z_from_r(r=r, totalN=totalN, zforCI_out=zforCI_out)
		
		r_stats <- r_from_z(z=z_stats$z, z_var=z_stats$z_var, totalN=z_stats$totalN, zforCI_out=zforCI_out)
		
		if (is.null(r_var))  r_var <- r_stats$r_var
		
		r_lb    <- r_stats$r_lb
		r_ub    <- r_stats$r_ub
	}
				
	# z from r
	z_stats <- z_from_r(r=r, totalN=totalN, zforCI_out=zforCI_out)
	z       <- z_stats$z
	z_var   <- z_stats$z_var
	z_lb    <- z_stats$z_lb
	z_ub    <- z_stats$z_ub
		
	# d from r 
	d_stats <- d_from_r(r=r, r_var=r_var, grp1_n=grp1_n, grp2_n=grp2_n, zforCI_out=zforCI_out) 
	d       <- d_stats$d
	d_var   <- d_stats$d_var
	d_lb    <- d_stats$d_lb
	d_ub    <- d_stats$d_ub

	# g from d
	if (!is.null(totalN)) {
		g_stats  <- g_from_d(d=d, d_var=d_var, totalN=totalN, zforCI_out=zforCI_out, gvar_type_OUT=gvar_type_OUT) 
		g        <- g_stats$g
		g_var    <- g_stats$g_var
		g_lb     <- g_stats$g_lb
		g_ub     <- g_stats$g_ub
	}

	# OR & logOR from d
	OR_stats  <- OR_from_d(d=d, d_var=d_var, grp1_n=grp1_n, grp2_n=grp2_n, zforCI_out=zforCI_out)
	OR        <- OR_stats$OR
	OR_lb     <- OR_stats$OR_lb
	OR_ub     <- OR_stats$OR_ub
	logOR     <- OR_stats$logOR
	logOR_var <- OR_stats$logOR_var
	logOR_lb  <- OR_stats$logOR_lb
	logOR_ub  <- OR_stats$logOR_ub
}



if (ES_type_IN == 'd') {
	
	d  <- ES
	
	if (!is.null(ES_var))  d_var <- ES_var
	
	# if d_var is NULL & grp Ns are available 
	if (all(is.na(d_var)) & grpNs) {
		# derived from the (hidden) esc.vd function in the esc package
		d_var <- (grp1_n + grp2_n) / (grp1_n * grp2_n) + (d * d) / (2 * (grp1_n + grp2_n))
	}

	if (!is.null(d_var)) {
		d_lb <- d - zforCI_out * sqrt(d_var)
		d_ub <- d + zforCI_out * sqrt(d_var)
	}
		
	# g from d
	if (!is.null(totalN)) {
		g_stats  <- g_from_d(d=d, d_var=d_var, totalN=totalN, grp1_n=grp1_n, grp2_n=grp2_n, 
		                     zforCI_out=zforCI_out, gvar_type_OUT=gvar_type_OUT) 
		g     <- g_stats$g
		g_var <- g_stats$g_var
		g_lb  <- g_stats$g_lb
		g_ub  <- g_stats$g_ub
	}

	# r from d
	r_stats <- r_from_d(d=d, d_var=d_var, totalN=totalN, zforCI_out=zforCI_out, grp1_n=grp1_n, grp2_n=grp2_n)
	r       <- r_stats$r
	r_var   <- r_stats$r_var
	r_lb    <- r_stats$r_lb
	r_ub    <- r_stats$r_ub
z
	# z from r
	z_stats <- z_from_r(r=r, totalN=totalN, zforCI_out=zforCI_out)
	z       <- z_stats$z
	z_var   <- z_stats$z_var
	z_lb    <- z_stats$z_lb
	z_ub    <- z_stats$z_ub

	# OR from d
	OR_stats  <- OR_from_d(d=d, d_var=d_var, grp1_n, grp2_n, zforCI_out=zforCI_out)
	OR        <- OR_stats$OR
	OR_lb     <- OR_stats$OR_lb
	OR_ub     <- OR_stats$OR_ub
	logOR     <- OR_stats$logOR
	logOR_var <- OR_stats$logOR_var
	logOR_lb  <- OR_stats$logOR_lb
	logOR_ub  <- OR_stats$logOR_ub
}




if (ES_type_IN == 'g') {
	
	g <- ES
	
	if (!is.null(ES_var))  g_var <- ES_var

	# if g_var is NULL & grp Ns are available 
	if (all(is.na(g_var)) & grpNs) {
		# 2009 Borenstein - Introduction to Meta-Analysis - book.pdf   - 27
		d_var <- (grp1_n + grp2_n) / (grp1_n * grp2_n) + (d * d) / (2 * (grp1_n + grp2_n))
		df <- (grp1_n + grp2_n) - 2			
		bigJ <- 1 - (3 / (4 * df - 1))			  
		if (gvar_type_OUT == 'g')  g_var <- bigJ^2 * d_var  
		if (gvar_type_OUT == 'd')  g_var <- d_var
	}

	if (!is.null(g_var)) {
		g_lb <- g - zforCI_out * sqrt(g_var)
		g_ub <- g + zforCI_out * sqrt(g_var)
	}

	# d from g
	d_stats <- d_from_g(g=g, g_var=g_var, totalN=totalN, grp1_n=grp1_n, grp2_n=grp2_n, zforCI_out=zforCI_out)
	d       <- d_stats$d
	d_var   <- d_stats$d_var
	d_lb    <- d_stats$d_lb
	d_ub    <- d_stats$d_ub
	
	# r from d
	r_stats <- r_from_d(d=d, d_var=d_var, totalN=totalN, grp1_n=grp1_n, grp2_n=grp2_n)
	r       <- r_stats$r
	r_var   <- r_stats$r_var
	r_lb    <- r_stats$r_lb
	r_ub    <- r_stats$r_ub

	# z from r
	z_stats <- z_from_r(r=r, totalN=totalN, zforCI_out=zforCI_out)
	z       <- z_stats$z
	z_var   <- z_stats$z_var
	z_lb    <- z_stats$z_lb
	z_ub    <- z_stats$z_ub

	# OR from d
	OR_stats  <- OR_from_d(d=d, d_var=d_var, grp1_n=grp1_n, grp2_n=grp2_n, zforCI_out=zforCI_out)
	OR        <- OR_stats$OR
	OR_lb     <- OR_stats$OR_lb
	OR_ub     <- OR_stats$OR_ub
	logOR     <- OR_stats$logOR
	logOR_var <- OR_stats$logOR_var
	logOR_lb  <- OR_stats$logOR_lb
	logOR_ub  <- OR_stats$logOR_ub
}



if (ES_type_IN == 'OR' | ES_type_IN == 'or') {

	OR <- ES
	if (!is.null(CI_in_lb))  OR_lb <- donnes[,CI_in_lb]
	if (!is.null(CI_in_ub))  OR_ub <- donnes[,CI_in_ub]
	
	# Converting between effect sizes.pdf
	# Meta-analysis Converting among effect sizes (2009 Borenstein).pdf
	# compute.es.pdf
		
	# convert OR to log OR
	logOR <- log(ES)

	if (!is.null(ci_dat)) {  
		# https://stats.stackexchange.com/questions/10375/how-to-calculate-standard-error-of-odds-ratios
		# https://www.ncbi.nlm.nih.gov/books/NBK431098/
		logOR_var <- (abs(logOR - log(ci_dat)) / zforCI_out)^2
	}

	# CIs for logOR
	if (!is.null(logOR_var)) {  
		logOR_lb <- logOR - zforCI_out * sqrt(logOR_var)
		logOR_ub <- logOR + zforCI_out * sqrt(logOR_var)
	}

	# d from logOR
	d <- logOR * (sqrt(3)/pi) 

	# d_var
	if (!is.null(logOR_var))  d_var <- logOR_var * (3/(pi^2))

	# d_var if ci_dat is NULL & grp Ns are available 
	if (is.null(ci_dat) & grpNs) {
		# derived from the (hidden) esc.vd function in the esc package
		d_var <- (grp1_n + grp2_n) / (grp1_n * grp2_n) + (d * d) / (2 * (grp1_n + grp2_n))
	}

	if (!is.null(d_var)) {
		d_lb <- d - zforCI_out * sqrt(d_var)
		d_ub <- d + zforCI_out * sqrt(d_var)
	}

	# no totalN, but have group Ns
	if (is.null(totalN) & !is.null(grp1_n) & !is.null(grp2_n))  
		totalN <- grp1_n + grp2_n

	# g from d
	if (!is.null(totalN)) {
		g_stats <- g_from_d(d=d, d_var=d_var, totalN=totalN, grp1_n=grp1_n, grp2_n=grp2_n, 
		                    zforCI_out=zforCI_out, gvar_type_OUT=gvar_type_OUT) 
		g       <- g_stats$g
		g_var   <- g_stats$g_var
		g_lb    <- g_stats$g_lb
		g_ub    <- g_stats$g_ub
	}

	# r from d
	r_stats <- r_from_d(d=d, d_var=d_var, totalN=totalN, grp1_n=grp1_n, grp2_n=grp2_n)
	r       <- r_stats$r
	r_var   <- r_stats$r_var
	r_lb    <- r_stats$r_lb
	r_ub    <- r_stats$r_ub

	# z from r
	z_stats <- z_from_r(r=r, totalN=totalN, zforCI_out=zforCI_out)
	z       <- z_stats$z
	z_var   <- z_stats$z_var
	z_lb    <- z_stats$z_lb
	z_ub    <- z_stats$z_ub
}	


if (is.null(totalN))  totalN <- as.numeric(rep(NA, length(ES)))


Output <- data.frame(totalN,
  r, r_var, r_lb, r_ub,
  d, d_var, d_lb, d_ub,
  g, g_var, g_lb, g_ub,
  z, z_var, z_lb, z_ub,
  OR, OR_lb, OR_ub,
  logOR, logOR_var, logOR_lb, logOR_ub)


if (verbose) {
	message('\n\n Effect size statistics:\n')
	print(round(Output,2))
}

class(Output) <- c("data.frame")

return(invisible(Output)) 
}





