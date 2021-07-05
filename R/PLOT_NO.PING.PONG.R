

PLOT_NO.PING.PONG <- function(nopingpongOutput, 
                              plot_this = c('NHST','MA','Bayes_SR'),
                              plot_save = FALSE, plot_save_type = 'png', plot_title=NULL, Xrange=NULL) {
	

noms <- names(nopingpongOutput)
noms2 <- -9999
if ( is.element('results_NHST',   noms) ) noms2 <- c(noms2, 'results_NHST') 
if ( is.element('results_CUM',    noms) ) noms2 <- c(noms2, 'results_CUM') 
if ( is.element('results_BA_SR',  noms) ) noms2 <- c(noms2, 'results_BA_SR') 
if ( is.element('results_BA_GEN', noms) ) noms2 <- c(noms2, 'results_BA_GEN') 
if ( is.element('results_BA_RAW', noms) ) noms2 <- c(noms2, 'results_BA_RAW') 
noms <- noms2[-1] 


# get the min & max LBs & UBs
dum <- nopingpongOutput[noms[1]][[1]]
minLB <- min(dum[,'ES.lb'])
maxUB <- max(dum[,'ES.ub'])

if (length(noms) > 1) {
	for (lupeL in 2:length(noms)) {		
		
		dum <- nopingpongOutput[noms[lupeL]][[1]]
				
		if (min(dum[,'ES.lb']) < minLB)  minLB <- min(dum[,'ES.lb'])
					
		if (max(dum[,'ES.ub']) > maxUB)  maxUB <- max(dum[,'ES.ub'])
	}
}



# # get the min & max LBs & UBs
# minLB <- min(nopingpongOutput[[1]][,'ES.lb'])
# maxUB <- max(nopingpongOutput[[1]][,'ES.ub'])

# if (length(nopingpongOutput) > 1) {
	# for (lupeL in 2:length(nopingpongOutput)) {		
		# if (min(nopingpongOutput[[lupeL]][,'ES.lb']) < minLB)  
			# minLB <- min(nopingpongOutput[[lupeL]][,'ES.lb'])		
		# if (max(nopingpongOutput[[lupeL]][,'ES.ub']) > maxUB)  
			# maxUB <- max(nopingpongOutput[[lupeL]][,'ES.ub'])
	# }
# }

if (is.null(Xrange)) {
	# range for xlim -- adding a % to each end
	intsize <- maxUB - minLB
	toadd <- intsize * .1
	Xrange <- c((minLB - toadd), (maxUB + toadd))
}

labES  <- 'Effect Size'
labseq <- 'Study Sequence'


if (plot_save == TRUE) {
	
	if (is.null(plot_save_type))  plot_save_type = 'png'
	
	if (plot_save_type == 'bitmap')
		bitmap(paste("Figure - ",plot_title,".png",sep=""), height=7, width=9, units='in', res=1200, pointsize=12)

	if (plot_save_type == 'tiff')
		tiff(paste("Figure - ",plot_title,".tiff",sep=""), height=7, width=9, units='in', res=1200, pointsize=12)
		
	if (plot_save_type == 'png')
		png(paste("Figure - ",plot_title,".tiff",sep=""), height=7, width=9, units='in', res=1200, pointsize=12)
		
	if (plot_save_type == 'jpeg')
		jpeg(paste("Figure - ",plot_title,".png",sep=""), height=7, width=9, units='in', res=1200, pointsize=12)
		
	if (plot_save_type == 'bmp')
		bmp(paste("Figure - ",plot_title,".tiff",sep=""), height=7, width=9, units='in', res=1200, pointsize=12)
}


if (length(plot_this) == 1)  par(mfrow=c(1,1), pty="m", mar=c(3,2,3,2) + 2.6)

if (length(plot_this) == 2)  par(mfrow=c(1,2), pty="m", mar=c(3,2,3,2) + 2.6)

if (length(plot_this) == 3)  par(mfrow=c(1,3), pty="m", mar=c(3,2,3,2) + 2.6)

if (length(plot_this) == 4)  par(mfrow=c(1,4), pty="m", mar=c(3,2,3,2) + 2.6)

if (length(plot_this) == 5)  par(mfrow=c(1,5), pty="m", mar=c(3,2,3,2) + 2.6)


oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))


# # if (length(plot_this) == 1) {layout(matrix(1,1, nrow=1, byrow=T)); layout.show(n=1)}

# if (length(plot_this) == 2) {layout(matrix(c(1, 2, 1, 2), nrow=2, byrow=T)); layout.show(n=2)}		

# if (length(plot_this) == 3) {layout(matrix(c(1,2,3,1,2,3,1,2,3), nrow=3, byrow=T)); 
	# layout.show(n=3)}	

# if (length(plot_this) == 4) {layout(matrix(c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4), nrow=4, byrow=T)); 
	# layout.show(n=4)}		

# if (length(plot_this) == 5) {layout(matrix(c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,
	# 1,2,3,4,5), nrow=5, byrow=T)); layout.show(n=5)}		



if ( (is.element('NHST', plot_this)) & (is.element('results_NHST', noms)) ) { 
	results_NHST <- nopingpongOutput$results_NHST
	plot(0,0, type="n", main="NHST", xlab=labES, ylab=labseq,
		font.main=1,font.lab=1,cex.main=1.9,cex.lab=1.6,cex.axis=1.3, 
		ylim=c(max(results_NHST[,'Study']),1), xlim=Xrange) 
		lines(c(0,0), c(max((results_NHST[,'Study'])),1), col=2, lty="dashed") 
		for (luper in 1:nrow(results_NHST)) { 
			lines(rbind(results_NHST[luper,'ES.lb'],results_NHST[luper,'ES.ub']), c(luper,luper), lwd=1)
			points(results_NHST[luper,'ES'], results_NHST[luper,'Study'], pch=19, cex=.5) 
		}
}

if ( (is.element('MA', plot_this)) & (is.element('results_CUM', noms)) ) { 
	results_CUM    <- nopingpongOutput$results_CUM
	plot(0,0, type="n", main='Cumul. MA', xlab=labES, ylab=labseq,
		font.main=1,font.lab=1,cex.main=1.9,cex.lab=1.6,cex.axis=1.3, 
		ylim=c(max((results_CUM[,'Study'])),1), xlim=Xrange)
		lines(c(0,0), c(max((results_CUM[,'Study'])),1), col=2, lty="dashed") 
		for (luper in 1:nrow(results_CUM)) { 
			lines(rbind(results_CUM[luper,'ES.lb'],results_CUM[luper,'ES.ub']), c(luper,luper), lwd=1)
			points(results_CUM[luper,'ES'], results_CUM[luper,'Study'], pch=19, cex=.5) 
		}
}

if ( (is.element('Bayes_SR', plot_this)) & (is.element('results_BA_SR', noms)) ) { 
	results_BA_SR <- nopingpongOutput$results_BA_SR
	plot(0,0, type="n", main='Bayesian', xlab=labES, ylab=labseq,
		font.main=1,font.lab=1,cex.main=1.9,cex.lab=1.6,cex.axis=1.3, 
		ylim=c(max((results_BA_SR[,'Study'])),1), xlim=Xrange)
		lines(c(0,0), c(max((results_BA_SR[,'Study'])),1), col=2, lty="dashed") 
		for (luper in 1:nrow(results_BA_SR)) { 
			lines(rbind(results_BA_SR[luper,'ES.lb'],results_BA_SR[luper,'ES.ub']), c(luper,luper), lwd=1)
			points(results_BA_SR[luper,'ES'], results_BA_SR[luper,'Study'], pch=19, cex=.5) 
		}
}

if ( (is.element('Bayes_gen', plot_this)) & (is.element('results_BA_GEN', noms)) ) { 
	results_BA_GEN <- nopingpongOutput$results_BA_GEN
	plot(0,0, type="n", main='Bayes (gen)', xlab=labES, ylab=labseq,
		font.main=1,font.lab=1,cex.main=1.9,cex.lab=1.6,cex.axis=1.3, 
		ylim=c(max((results_BA_GEN[,'Study'])),1), xlim=Xrange)
		lines(c(0,0), c(max((results_BA_GEN[,'Study'])),1), col=2, lty="dashed") 
		for (luper in 1:nrow(results_BA_GEN)) { 
			lines(rbind(results_BA_GEN[luper,'ES.lb'],results_BA_GEN[luper,'ES.ub']), c(luper,luper), lwd=1)
			points(results_BA_GEN[luper,'ES'], results_BA_GEN[luper,'Study'], pch=19, cex=.5) 
		}
}

if ( (is.element('Bayes_raw', plot_this)) & (is.element('results_BA_RAW', noms)) ) { 
	results_BA_RAW <- nopingpongOutput$results_BA_RAW
	plot(0,0, type="n", main='Bayes (raw)', xlab=labES, ylab=labseq,
		font.main=1,font.lab=1,cex.main=1.9,cex.lab=1.6,cex.axis=1.3, 
		ylim=c(max((results_BA_RAW[,'Study'])),1), xlim=Xrange)
		lines(c(0,0), c(max((results_BA_RAW[,'Study'])),1), col=2, lty="dashed") 
		for (luper in 1:nrow(results_BA_RAW)) { 
			lines(rbind(results_BA_RAW[luper,'ES.lb'],results_BA_RAW[luper,'ES.ub']), c(luper,luper), lwd=1)
			points(results_BA_RAW[luper,'ES'], results_BA_RAW[luper,'Study'], pch=19, cex=.5) 
		}
}

if (plot_save == TRUE)  dev.off()

}




