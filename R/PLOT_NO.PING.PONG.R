

PLOT_NO.PING.PONG <- function(nopingpongOutput, 
                              plot_this = c('NHST','CUM_META','BAYES_SR'),
                              plot_save = FALSE, plot_save_type = 'png', plot_title=NULL, Xrange=NULL) {
	

# not all elements in plot_this may exist in nopingpongOutput, so find the plot-able data

# first, identify what data appears in nopingpongOutput
noms <- names(nopingpongOutput)
noms2 <- c()
if ( is.element('results_NHST',      noms) )   noms2 <- c(noms2, 'NHST') 
if ( is.element('results_CUM_META',  noms) )   noms2 <- c(noms2, 'CUM_META') 
if ( is.element('results_BAYES_SR',  noms) )   noms2 <- c(noms2, 'BAYES_SR') 
if ( is.element('results_BAYES_GEN', noms) )   noms2 <- c(noms2, 'BAYES_GEN') 
if ( is.element('results_BAYES_RAW', noms) )   noms2 <- c(noms2, 'BAYES_RAW') 

# then find the common elements
# plot_this = c('NHST','CUM_META')
# noms = c('results_NHST','results_CUM_META','results_BAYES_SR')
plot_this_2 <- Reduce(intersect, list(plot_this,noms2))
 
 
noms3 <- paste("results_", noms2, sep='')


 
# get the min & max LBs & UBs
dum <- nopingpongOutput[noms3[1]][[1]]
minLB <- min(dum[,'ES.lb'])
maxUB <- max(dum[,'ES.ub'])

if (length(noms3) > 1) {
	for (lupeL in 2:length(noms3)) {		
		
		dum <- nopingpongOutput[noms3[lupeL]][[1]]
				
		if (min(dum[,'ES.lb']) < minLB)  minLB <- min(dum[,'ES.lb'])
					
		if (max(dum[,'ES.ub']) > maxUB)  maxUB <- max(dum[,'ES.ub'])
	}
}


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
		bitmap(paste("Figure - ",plot_title,".bitmap",sep=""), height=7, width=9, units='in', res=1200, pointsize=12)

	if (plot_save_type == 'tiff')
		tiff(paste("Figure - ",plot_title,".tiff",sep=""), height=7, width=9, units='in', res=1200, pointsize=12)
		
	if (plot_save_type == 'png')
		png(paste("Figure - ",plot_title,".png",sep=""), height=7, width=9, units='in', res=1200, pointsize=12)
		
	if (plot_save_type == 'jpeg')
		jpeg(paste("Figure - ",plot_title,".jpeg",sep=""), height=7, width=9, units='in', res=1200, pointsize=12)
		
	if (plot_save_type == 'bmp')
		bmp(paste("Figure - ",plot_title,".bmp",sep=""), height=7, width=9, units='in', res=1200, pointsize=12)
}


if (length(plot_this_2) == 1)  par(mfrow=c(1,1), pty="m", mar=c(3,2,3,2) + 2.6)

if (length(plot_this_2) == 2)  par(mfrow=c(1,2), pty="m", mar=c(3,2,3,2) + 2.6)

if (length(plot_this_2) == 3)  par(mfrow=c(1,3), pty="m", mar=c(3,2,3,2) + 2.6)

if (length(plot_this_2) == 4)  par(mfrow=c(1,4), pty="m", mar=c(3,2,3,2) + 2.6)

if (length(plot_this_2) == 5)  par(mfrow=c(1,5), pty="m", mar=c(3,2,3,2) + 2.6)


oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))


if ( is.element('NHST', plot_this_2) )  { 
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

if ( is.element('CUM_META', plot_this_2) ) { 
	results_CUM_META    <- nopingpongOutput$results_CUM_META
	plot(0,0, type="n", main='Cumul. MA', xlab=labES, ylab=labseq,
		font.main=1,font.lab=1,cex.main=1.9,cex.lab=1.6,cex.axis=1.3, 
		ylim=c(max((results_CUM_META[,'Study'])),1), xlim=Xrange)
		lines(c(0,0), c(max((results_CUM_META[,'Study'])),1), col=2, lty="dashed") 
		for (luper in 1:nrow(results_CUM_META)) { 
			lines(rbind(results_CUM_META[luper,'ES.lb'],results_CUM_META[luper,'ES.ub']), c(luper,luper), lwd=1)
			points(results_CUM_META[luper,'ES'], results_CUM_META[luper,'Study'], pch=19, cex=.5) 
		}
}

if ( is.element('BAYES_SR', plot_this_2) ) { 
	results_BAYES_SR <- nopingpongOutput$results_BAYES_SR
	plot(0,0, type="n", main='Bayes (SR)', xlab=labES, ylab=labseq,
		font.main=1,font.lab=1,cex.main=1.9,cex.lab=1.6,cex.axis=1.3, 
		ylim=c(max((results_BAYES_SR[,'Study'])),1), xlim=Xrange)
		lines(c(0,0), c(max((results_BAYES_SR[,'Study'])),1), col=2, lty="dashed") 
		for (luper in 1:nrow(results_BAYES_SR)) { 
			lines(rbind(results_BAYES_SR[luper,'ES.lb'],results_BAYES_SR[luper,'ES.ub']), c(luper,luper), lwd=1)
			points(results_BAYES_SR[luper,'ES'], results_BAYES_SR[luper,'Study'], pch=19, cex=.5) 
		}
}

if ( is.element('BAYES_GEN', plot_this_2) ) { 
	results_BAYES_GEN <- nopingpongOutput$results_BAYES_GEN
	plot(0,0, type="n", main='Bayes (gen)', xlab=labES, ylab=labseq,
		font.main=1,font.lab=1,cex.main=1.9,cex.lab=1.6,cex.axis=1.3, 
		ylim=c(max((results_BAYES_GEN[,'Study'])),1), xlim=Xrange)
		lines(c(0,0), c(max((results_BAYES_GEN[,'Study'])),1), col=2, lty="dashed") 
		for (luper in 1:nrow(results_BAYES_GEN)) { 
			lines(rbind(results_BAYES_GEN[luper,'ES.lb'],results_BAYES_GEN[luper,'ES.ub']), c(luper,luper), lwd=1)
			points(results_BAYES_GEN[luper,'ES'], results_BAYES_GEN[luper,'Study'], pch=19, cex=.5) 
		}
}

if ( is.element('BAYES_RAW', plot_this_2) ) { 
	results_BAYES_RAW <- nopingpongOutput$results_BAYES_RAW
	plot(0,0, type="n", main='Bayes (raw)', xlab=labES, ylab=labseq,
		font.main=1,font.lab=1,cex.main=1.9,cex.lab=1.6,cex.axis=1.3, 
		ylim=c(max((results_BAYES_RAW[,'Study'])),1), xlim=Xrange)
		lines(c(0,0), c(max((results_BAYES_RAW[,'Study'])),1), col=2, lty="dashed") 
		for (luper in 1:nrow(results_BAYES_RAW)) { 
			lines(rbind(results_BAYES_RAW[luper,'ES.lb'],results_BAYES_RAW[luper,'ES.ub']), c(luper,luper), lwd=1)
			points(results_BAYES_RAW[luper,'ES'], results_BAYES_RAW[luper,'Study'], pch=19, cex=.5) 
		}
}

if (plot_save == TRUE)  dev.off()

}




