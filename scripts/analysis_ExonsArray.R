#setwd('~/Documents/Work/_x_MicroArray/code_and_results/')

source('functions.R')


fig.types = c('F', 'S', 'O') # F for main figures, S for supplementary figures, O for other figures
panels = c('A','B','C','D','E','F','G','H','I','J')

main.fig.nb = 1 ;main.fig.panel.nb = 1
supp.fig.nb = 1 ;supp.fig.panel.nb = 1
other.fig.nb = 1 ;other.fig.panel.nb = 1



load = FALSE
do.step = FALSE
plot = FALSE
load.for.plot = TRUE

# process of the original data
if(do.step){	
	cat('PROCESS THE ORIGINAL DATA \n')
	e = read.table('../complete_run_experiments/RMA_2_days_together/rma-sketch.summary.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
	p = read.table('../complete_run_experiments/RMA_2_days_together/dabg.summary.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
	
	annot.affy = read.table('../complementary_data/annotation_tables/MoEx-1_0-st-v1.na32.mm9.probeset.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE)
	
	amit.symbol = read.table('../complementary_data/annotation_tables/symbol_probes.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)
	amit.ex_in = read.table('../complementary_data/annotation_tables/probesetID_exon_intron.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)
	
	colnames.e = c(paste('ZT0',seq(0,8,by = 2), sep = ''), paste('ZT',seq(10,46,by = 2), sep = ''))
	colnames.p = paste(colnames.e, '.pval', sep ='')
	perm = c(1,2,10:13,3:9,14,22:25,15:21); e = e[,perm]; p = p[,perm]; colnames(e)[2:25] = colnames.e;  colnames(p)[2:25] = colnames.p; 
	N = nrow(e)

	m = match(e$probeset_id, amit.symbol[,2])
	gene = amit.symbol[m,1]
	
	m = match(e$probeset_id, amit.ex_in[,1])
	ex_in = amit.ex_in[m,2]
	
	m = match(e$probeset_id, annot.affy$probeset_id)
	
	
	t = data.frame(e[,2:25], p[,2:25], probeset.ID = e$probeset_id, gene = gene, ex_in = ex_in, chr = annot.affy$seqname[m] ,strand = annot.affy$strand[m] , start = annot.affy$start[m], stop = annot.affy$stop[m] , transcript.cluster.ID =  annot.affy$transcript_cluster_id[m], probeset.type = annot.affy$probeset_type[m], cross.hyb.type = annot.affy$crosshyb_type[m], level = annot.affy$level[m]  , kept = TRUE, rejection.code = 'OK' ,stringsAsFactors = FALSE)
	save(t,file = 'Rdata/raw.Rdata')  ##### the raw data of exon and intron and dabg summary with p values
}

#write.table(t, "Table_Raw_data.txt",sep='\t',quote=F,col.names=TRUE, row.names=F)
my.do.step = TRUE
if(my.do.step)
{
	load('Rawdata.Rdata')
	t.raw = t
	save(t.raw,file = 'my_Rawdata.Rdata')
	load('my_Rawdata.Rdata')
	tt = t.raw
	e = read.table('Raw_data_Annotation/rma-sketch.summary.txt',header = TRUE, sep = '\t', stringsAsFactors = FALSE)
}
if(plot){
	if(load.for.plot){load('Rdata/raw.Rdata')}
	probeset.types = c("main","normgene->exon","normgene->intron","control->affx","control->bgp->antigenomic","control->bgp->genomic","FLmRNA->unmapped")
	for(type  in probeset.types){cat(type, '&', sum(t$probeset.type == type), '\\\\ \n')}
	boxplot(t[,2]~t$probeset.type, col = rainbow(n = 7, s = 0.2, v = 0.9), border = rainbow(n = 7, s = 0.5, v = 0.5), pch = 16, lwd = 2, ylab = 'RMA - log2 scale')
}


# #do the quantile normalisation at the beginning

# if(do.step){
	# if(load){load('Rdata/raw.Rdata')}
	# t$ex_in[is.na(t$ex_in)] = 2

	# library(preprocessCore)
	
	# ZT = 1:24

	# i.ex = which(t$kept & (t$ex_in == 1))
	# i.int = which(t$kept & (t$ex_in == -1))

	# #t.norm.ex = normalize.quantiles(as.matrix(t[i.ex, ZT]))
	# #t[i.ex, ZT] = t.norm.ex; rm(t.norm.ex)
	
	# #t.norm.int = normalize.quantiles(as.matrix(t[i.int, ZT]))
	# #t[i.int, ZT] = t.norm.int; rm(t.norm.int)

	# t.norm = normalize.quantiles(as.matrix(t[c(i.ex,i.int), ZT]))
	# t[c(i.ex,i.int), ZT] = t.norm; rm(t.norm)
	
	# t.unlog = 2^t[,ZT]; colnames(t.unlog) = paste(colnames(t)[ZT],'unlog',sep = '.')
	# t = data.frame(t, t.unlog, stringsAsFactors = FALSE); rm(t.unlog)
	# ZT.unlog = grep('.unlog', colnames(t))
	# t$mean = apply(t[, ZT.unlog], 1, mean);
	# t_rel.ampl = t[, ZT.unlog]/t$mean; colnames(t_rel.ampl) = paste(colnames(t)[ZT],'rel.ampl',sep = '.')
	# t = data.frame(t, t_rel.ampl, stringsAsFactors = FALSE); rm(t_rel.ampl)	
	# save(t, file = 'Rdata/raw_quant_norm_together_and_unloged_and_rel.Rdata')
# }


# if(plot){
	# if(load.for.plot){load('Rdata/raw_quant_norm_together_and_unloged_and_rel.Rdata')}
	# source('functions.R')
	# ZT = c(1:24)
	# zt.rel = grep('.rel.ampl', colnames(t))

	# pdf('plots/distib_raw_and_quant_norm_together.pdf',width = 12, height = 16)
	# par(mfcol = c(4,3))
	# plot.distrib.t(t = t, zt = ZT, zt.rel = zt.rel, main = "just out from RMA - quant norm together")
	# dev.off()
# }




# #do the quantile normalisation at the beginning

# if(do.step){
	# if(load){load('Rdata/raw.Rdata')}
	# t$ex_in[is.na(t$ex_in)] = 2

	# library(preprocessCore)
	
	# ZT = 1:24

	# i.ex = which(t$kept & (t$ex_in == 1))
	# i.int = which(t$kept & (t$ex_in == -1))

	# t.norm.ex = normalize.quantiles(as.matrix(t[i.ex, ZT]))
	# t[i.ex, ZT] = t.norm.ex; rm(t.norm.ex)
	
	# t.norm.int = normalize.quantiles(as.matrix(t[i.int, ZT]))
	# t[i.int, ZT] = t.norm.int; rm(t.norm.int)

	# #t.norm = normalize.quantiles(as.matrix(t[c(i.ex,i.int), ZT]))
	# #t[c(i.ex,i.int), ZT] = t.norm; rm(t.norm)
	
	# t.unlog = 2^t[,ZT]; colnames(t.unlog) = paste(colnames(t)[ZT],'unlog',sep = '.')
	# t = data.frame(t, t.unlog, stringsAsFactors = FALSE); rm(t.unlog)
	# ZT.unlog = grep('.unlog', colnames(t))
	# t$mean = apply(t[, ZT.unlog], 1, mean);
	# t_rel.ampl = t[, ZT.unlog]/t$mean; colnames(t_rel.ampl) = paste(colnames(t)[ZT],'rel.ampl',sep = '.')
	# t = data.frame(t, t_rel.ampl, stringsAsFactors = FALSE); rm(t_rel.ampl)	
	# save(t, file = 'Rdata/raw_quant_norm_sep_and_unloged_and_rel.Rdata')
# }


# if(plot){
	# if(load.for.plot){load('Rdata/raw_quant_norm_sep_and_unloged_and_rel.Rdata')}
	# source('functions.R')
	# ZT = c(1:24)
	# zt.rel = grep('.rel.ampl', colnames(t))

	# pdf('plots/distib_raw_and_quant_norm_sep.pdf',width = 12, height = 16)
	# par(mfcol = c(4,3))
	# plot.distrib.t(t = t, zt = ZT, zt.rel = zt.rel, main = "just out from RMA - quant norm sep")
	# dev.off()
# }



if(do.step){
	if(load){load('Rdata/raw.Rdata')}
	t$ex_in[is.na(t$ex_in)] = 2
	ZT = c(1:24)
	t.unlog = 2^t[,ZT]; colnames(t.unlog) = paste(colnames(t)[ZT],'unlog',sep = '.')
	t = data.frame(t, t.unlog, stringsAsFactors = FALSE); rm(t.unlog)
	ZT.unlog = grep('.unlog', colnames(t))
	t$mean = apply(t[, ZT.unlog], 1, mean); #####compute the ratio between linear singals and their means
	t_rel.ampl = t[, ZT.unlog]/t$mean; colnames(t_rel.ampl) = paste(colnames(t)[ZT],'rel.ampl',sep = '.')
	t = data.frame(t, t_rel.ampl, stringsAsFactors = FALSE); rm(t_rel.ampl)	
	save(t, file = 'Rdata/raw_and_unloged_and_rel.Rdata')
	
}

#write.table(t, "Table_Raw_data_Unlog_Rel.txt",sep='\t',quote=F,col.names=TRUE, row.names=F)
if(my.do.step)
{
	t.raw.unlog.rel = t
	save(t.raw.unlog.rel,file = 'Rawdata_Unlog_Rel.Rdata')
}

if(plot){
	if(load.for.plot){load('Rdata/raw_and_unloged_and_rel.Rdata')}
	source('functions.R')
	ZT = c(1:24)
	zt.rel = grep('.rel.ampl', colnames(t))

	pdf('plots/distib_raw.pdf',width = 12, height = 16)
	par(mfcol = c(4,3))
	plot.distrib.t(t = t, zt = ZT, zt.rel = zt.rel, main = "just out from RMA")
	dev.off()
}

# compute rel. ampl. on each day separately
if(do.step){
	if(load){load('Rdata/raw_and_unloged_and_rel.Rdata')}	
	zt.unlog = grep('.unlog', colnames(t))
	zt = 1:24
	mean1 = apply(t[,zt.unlog[1:12]],1,mean)
	mean2 = apply(t[,zt.unlog[13:24]],1,mean)
	t_rel.ampl1 = t[, zt.unlog[1:12]]/mean1; colnames(t_rel.ampl1) = paste(colnames(t)[zt[1:12]],'rel.ampl',sep = '.')
	t_rel.ampl2 = t[, zt.unlog[13:24]]/mean2; colnames(t_rel.ampl2) = paste(colnames(t)[zt[13:24]],'rel.ampl',sep = '.')
	zt.rel = grep('.rel.ampl', colnames(t))
	t[, zt.rel[1:12]] = t_rel.ampl1
	t[, zt.rel[13:24]] = t_rel.ampl2
	save(t,file = 'Rdata/raw_and_unloged_and_rel_sep_day.Rdata')
}



if(plot){
	if(load.for.plot){load('Rdata/raw_and_unloged_and_rel_sep_day.Rdata')}
	source('functions.R')
	ZT = c(1:24)
	zt.rel = grep('.rel.ampl', colnames(t))

	pdf('plots/distib_raw_rel_sep_day.pdf',width = 12, height = 16)
	par(mfcol = c(4,3))
	plot.distrib.t(t = t, zt = ZT, zt.rel = zt.rel, main = "just out from RMA - rel ampl sep")
	dev.off()
}


if(plot){
	if(load.for.plot){load('Rdata/raw_and_unloged_and_rel_sep_day.Rdata')}
	
	pdf('plots/tp_vs_tp_after_sep_quant_norm.pdf')
	
	for(i in 2:24){
		for(j in 1:(i-1)){
			cat(i, '\t', j, '\n')
			#plot(t[,i],t[,j], main = paste(colnames(t)[i],'vs',colnames(t)[j]), cex = 0.5, pch = 16, col = rgb(0,0,0,0.2))
			#abline(a = 0, b = 1, col ='red')
			
			#d = kde2d(x = t[,i]-t[,j], y = t[,j], n =  100)
			#contour(d, nlevels = 25, drawlabel = FALSE, add = TRUE)
			#filled.contour(x = d$x,y = d$y, z = d$z )
			#abline(a = 0, b = 1, col ='red')
			
			# d = kde2d(x = t[,i], y = t[,j], n =  100)
			# image(d, col = paste('gray',seq(100,0,by = -5),sep = ''), main = paste(colnames(t)[i],'vs',colnames(t)[j]), xlab  =colnames(t)[i], ylab =  colnames(t)[j])
			# abline(a = 0, b = 1, col ='red')
			
			par(mfrow = c(2,2))
			
			d = kde2d(x = t[,i]-t[,j], y = t[,j], n =  100)
			image(d, col = paste('gray',round(seq(95,1,length = 20)),sep = ''), main = paste(colnames(t)[i],'vs',colnames(t)[j]), xlab  =colnames(t)[i], ylab =  colnames(t)[j], zlim = c(2*10^(-4),max(d$z)))
			abline(v = 0, col ='red')
					
			d = kde2d(x = t[(t$ex_in == 1),i]-t[(t$ex_in == 1),j], y = t[(t$ex_in == 1),j], n =  100)
			image(d, col = paste('gray',round(seq(95,1,length = 20)),sep = ''), main = paste(colnames(t)[i],'vs',colnames(t)[j],'\nexon'), xlab  =paste(colnames(t)[i],'-',colnames(t)[j]), ylab =  colnames(t)[j], zlim = c(2*10^(-4),max(d$z)))
			abline(v = 0, col ='red')
			text(-2, 14, paste(colnames(t)[i],'<',colnames(t)[j]), pos = 2)
			text(2, 14, paste(colnames(t)[i],'>',colnames(t)[j]), pos = 4)
			
			d = kde2d(x = t[(t$ex_in == -1),i]-t[(t$ex_in == -1),j], y = t[(t$ex_in == -1),j], n =  100)
			image(d, col = paste('gray',round(seq(95,1,length = 20)),sep = ''), main = paste(colnames(t)[i],'vs',colnames(t)[j],'\nintron'), xlab  =paste(colnames(t)[i],'-',colnames(t)[j]), ylab =  colnames(t)[j], zlim = c(2*10^(-4),max(d$z)))
			abline(v = 0, col ='red')
			text(-2, 14, paste(colnames(t)[i],'<',colnames(t)[j]), pos = 2)
			text(2, 14, paste(colnames(t)[i],'>',colnames(t)[j]), pos = 4)
			
			
			d = kde2d(x = t[((t$ex_in%%2) == 0),i]-t[((t$ex_in%%2) == 0),j], y = t[((t$ex_in%%2) == 0),j], n =  100)
			image(d, col = paste('gray',round(seq(95,1,length = 20)),sep = ''), main = paste(colnames(t)[i],'vs',colnames(t)[j],'\nother'), xlab  =paste(colnames(t)[i],'-',colnames(t)[j]), ylab =  colnames(t)[j], zlim = c(2*10^(-4),max(d$z)))
			abline(v = 0, col ='red')
			text(-2, 14, paste(colnames(t)[i],'<',colnames(t)[j]), pos = 2)
			text(2, 14, paste(colnames(t)[i],'>',colnames(t)[j]), pos = 4)
	
	
		}
	}
	dev.off()
}

############################################### start the selection of probesets

# remove unannotated probes
if(do.step){
	cat('REMOVE UNANNOTATED PROBESETS \n')
	if(load){load('Rdata/raw_and_unloged_and_rel.Rdata')}	#### still use the ratio normalized by means of two days
	t = remove.useless.line(t = t)
	save(t,file = 'Rdata/raw_annotated.Rdata')
}




# remove high pvalues
if(do.step){
	cat('REMOVE HIGH PVALUES \n')
	if(load){load('Rdata/raw_annotated.Rdata')}
	ZT = c(1:24)
	ind = passes.fdr(t = t, fdr = 0.1); t$rejection.code[!ind] = 'pval'; t$kept = t$kept&ind; 
	save(t,file = 'Rdata/raw_annotated_low_pvalues.Rdata')
}


if(plot){
	if(load.for.plot){load('Rdata/raw_annotated_low_pvalues.Rdata')}
	source('functions.R')
	ZT = c(1:24)
	zt.rel = grep('.rel.ampl', colnames(t))

	pdf('plots/distib_raw_low_pvalues.pdf',width = 12, height = 16)
	par(mfcol = c(4,3))
	plot.distrib.t(t = t, zt = ZT, zt.rel = zt.rel, main = "unannotated removed - low pval")
	dev.off()
}



# remove low probesets
if(do.step){
	cat('REMOVE LOW PROBESETS \n')
	if(load){load('Rdata/raw_annotated_low_pvalues.Rdata')}
	source('functions.R')
	ZT = c(1:24)
	ind = remove.low.probesets(t = t, th = 3); t$rejection.code[t$kept & !ind] = 'too low'; t$kept = t$kept & ind; 
	save(t,file = 'Rdata/raw_annotated_above_th.Rdata')
}



if(plot){
	if(load.for.plot){load('Rdata/raw_annotated_above_th.Rdata')}
	source('functions.R')
	ZT = c(1:24)
	zt.rel = grep('.rel.ampl', colnames(t))

	pdf('plots/distib_raw_above_th.pdf',width = 12, height = 16)
	par(mfcol = c(4,3))
	plot.distrib.t(t = t, zt = ZT, zt.rel = zt.rel, main = "unannotated removed - above th")
	dev.off()
}



# remove high Xhyb
if(do.step){
	cat('REMOVE HIGH X-HYB \n')
	if(load){load('Rdata/raw_annotated_above_th.Rdata')}
	source('functions.R')
	ind = remove.high.Xhyb(t = t, th = 3); t$rejection.code[t$kept & !ind] = 'Xhyb'; t$kept = t$kept & ind;
	save(t,file = 'Rdata/raw_annotated_no_Xhyb.Rdata')
}


if(plot){
	if(load.for.plot){load('Rdata/raw_annotated_no_Xhyb.Rdata')}
	source('functions.R')
	ZT = c(1:24)
	zt.rel = grep('.rel.ampl', colnames(t))

	pdf('plots/distib_raw_no_Xhyb.pdf',width = 12, height = 16)
	par(mfcol = c(4,3))
	plot.distrib.t(t = t, zt = ZT, zt.rel = zt.rel, main = "unannotated removed - no Xhyb")
	dev.off()
}


if(plot){
	if(load.for.plot){load('Rdata/raw_annotated_no_Xhyb.Rdata')}
	ZT = 1:24
	pdf('plots/Kept_and_rejected_probeset_histogram.pdf', width = 7, height = 5)
	
	ha = hist(as.numeric(unlist(t[t$kept,ZT])), breaks = seq(0,16,by = 0.2), plot = FALSE)
	hr.pval = hist(as.numeric(unlist(t[t$rejection.code == 'pval',ZT])), breaks = seq(0,16,by = 0.2), plot = FALSE)
	hr.bg = hist(as.numeric(unlist(t[t$rejection.code == 'too low',ZT])), breaks = seq(0,16,by = 0.2), plot = FALSE)
	hr.Xhyb = hist(as.numeric(unlist(t[t$rejection.code == 'Xhyb',ZT])), breaks = seq(0,16,by = 0.2), plot = FALSE)
	
	plot(ha$mids, ha$counts, type = 'l', lwd = 2, col = 'green4', ylim = range(ha$counts, hr$counts), xlab = 'signal [log2 scale]', ylab = 'counts',main = 'pool of all time-points', axes = FALSE)
	points(hr.pval$mids, hr.pval$counts, type = 'l', lwd = 2, col = 'coral2')
	points(hr.bg$mids, hr.bg$counts, type = 'l', lwd = 2, col = 'coral3')
	points(hr.Xhyb$mids, hr.Xhyb$counts, type = 'l', lwd = 2, col = 'coral4')
	legend('topright', c('probeset accepted','pval too high', 'below bg', ' Xhyb' ), col = c('green4', 'coral2', 'coral3', 'coral4'), lty = 1, lwd= 2, bty = 'n')
	axis(1, at = seq(0,16,by = 2))
	axis(2)
	dev.off()
}


# # unlog 
# if(do.step){
	# cat('UNLOG THE DATA \n')
	# if(load){load('Rdata/raw_annotated_no_Xhyb.Rdata')}
	# ZT = c(1:24)
	# t[,ZT] = 2^t[,ZT]
	# save(t,file = 'Rdata/raw_unloged.Rdata')
# } 
if(my.do.step)
{
	t.annotated_no_Xhyb = t
	save(t.annotated_no_Xhyb, file="my_raw_annotated_no_Xhyb.Rdata")
}
	
############################################### finish the selection of probesets


############################################### start to remove background #### not really understood
mynewdo = TRUE
if(mynewdo) ## we remove the background with 25% quantile of signals of each time points
{
	load('my_raw_annotated_no_Xhyb.Rdata')
	tt = t.annotated_no_Xhyb
	rm(t.annotated_no_Xhyb)
	data = as.matrix(tt[,c(1:24)])
	
	mu = rep(NA, 24)
	ii = c(1:24)
	for(i in ii)
	{
		mu[i] = quantile(data[,i], 0.25)
		hist(data[,i],breaks=100)
		abline(v=mu,col='red',lwd=2.0)
	}
	
	remove.background = function(x,mu)
	{
		#bb = 2^4
		x = x - 2^mu
		kk = which(x<=1.0)
		x[kk] = 1.0
		#c(x, nb.zero=length(kk))
		x
	}

	ZT.unlog = grep('.unlog', colnames(tt))
	tt.xx = as.matrix(tt[,ZT.unlog])
	for(i in ii)
	{
		tt.xx[,i] = remove.background(t(tt.xx[,i]), mu[i])
	}
	tt.xx = data.frame(tt.xx)
	
	which.zero = function(x)
	{
		kk = which(x<=1.0)
		length(kk)
	}
	nb.zero = apply(tt.xx, 1, which.zero)
	
	ii = which(nb.zero>8)
	tt[, ZT.unlog] = tt.xx
	tt$kept[ii] = FALSE
	tt$rejection.code[ii] = 'lower.background'
	#tt.temp = tt
	
	#load('my_raw_annotated_no_Xhyb.Rdata')
	#tt.s = t.annotated_no_Xhyb
	#rm(t.annotated_no_Xhyb)
	#tt$rejection.code = tt.s$rejection.code
	#tt = tt.temp
	#tt$rejection.code[ii] = 'lower.background'
	#tt = tt[ii,]
	#keep = keep[ii,]
	#tt[, ZT.unlog] = keep[,c(1:24)];	
	#rm(keep)
	#save(tt,file = 'raw_remove_background_quantile_0.20.Rdata')
	
	ZT = 1:24; 
	tt[,ZT] = log2(tt[, ZT.unlog])

	ZT.rel.ampl = grep('.rel.ampl', colnames(tt))
	tt$mean = apply(tt[, ZT.unlog], 1, mean);
	tt[, ZT.rel.ampl] = tt[, ZT.unlog]/tt$mean;
	save(tt,file = 'raw_remove_background_quantile_0.25.Rdata')
}

if(my.do.step) ## we remove the background with a constant of 2^4 each time points
{
	load('my_raw_annotated_no_Xhyb.Rdata')
	tt = t.annotated_no_Xhyb
	fun.counting = function(x)
	{
		bb = 2^4
		x = x - bb
		kk = which(x<=1.0)
		x[kk] = 1.0
		c(x, nb.zero=length(kk))
	}
	
	ZT.unlog = grep('.unlog', colnames(tt))
	tt.xx = as.matrix(tt[,ZT.unlog])
	keep = t(apply(tt.xx, 1, fun.counting))
	keep = data.frame(keep)
	ii = which(keep$nb.zero<12)
	tt = tt[ii,]
	keep = keep[ii,]
	tt[, ZT.unlog] = keep[,c(1:24)];	
	rm(keep)
	
	ZT = 1:24; 
	tt[,ZT] = log2(tt[, ZT.unlog])
	
	ZT.rel.ampl = grep('.rel.ampl', colnames(tt))
	tt$mean = apply(tt[, ZT.unlog], 1, mean);
	tt[, ZT.rel.ampl] = tt[, ZT.unlog]/tt$mean;
	save(tt,file = 'raw_remove_background_16.Rdata')
}


if(do.step){
	cat('REMOVE THE BACKGROUND \n')
	if(load){load('Rdata/raw_annotated_no_Xhyb.Rdata')}
	ZT.unlog = grep('.unlog', colnames(t))
	bg = 2^5;
	t_norm = bg*exp(-t[, ZT.unlog]/bg)+(t[, ZT.unlog]-bg)+1
	t[, ZT.unlog] = t_norm;	
	ZT = 1:24; t[,ZT] = log2(t[, ZT.unlog])
	ZT.rel.ampl = grep('.rel.ampl', colnames(t))
	t$mean = apply(t[, ZT.unlog], 1, mean);
	t[, ZT.rel.ampl] = t[, ZT.unlog]/t$mean;

	save(t,file = 'Rdata/raw_norm.Rdata')
	}

###### test the formula of removing background 
xx = seq(1,1600)/100
xx = 2^xx
bb = 2^4
yy = bb*exp(-xx/bb)+xx-bb+1
yy1 = xx-bb
xx = log2(xx)
plot(xx,yy,type='l',col='red',log='y')
points(xx,yy1,type='l',col='blue')
###### test the formula of removing background 

if(plot){
	if(load.for.plot){load('Rdata/raw_norm.Rdata')}
	source('functions.R')
	ZT = c(1:24)
	zt.rel = grep('.rel.ampl', colnames(t))

	pdf('plots/distib_raw_bg_removed.pdf',width = 12, height = 16)
	par(mfcol = c(4,3))
	plot.distrib.t(t = t, zt = ZT, zt.rel = zt.rel, main = "unannotated removed - bg removed")
	dev.off()
}


# remove background --- DIFFERENT BG per t-p
if(do.step){
	cat('REMOVE THE BACKGROUND \n')
	if(load){load('Rdata/raw_annotated_no_Xhyb.Rdata')}
	ZT.unlog = grep('.unlog', colnames(t))
	ZT = 1:24
	
	st = apply(t[t$kept, ZT],2,sort)
	bgs = st[round(sum(t$kept,na.rm = TRUE)/5),]
		
	t_norm = bgs*exp(-t[, ZT.unlog]/bgs)+(t[, ZT.unlog]-bgs)+1
	
	t[, ZT.unlog] = t_norm;	rm(t_norm)
	ZT = 1:24; t[,ZT] = log2(t[, ZT.unlog])
	ZT.rel.ampl = grep('.rel.ampl', colnames(t))
	t$mean = apply(t[, ZT.unlog], 1, mean);
	t[, ZT.rel.ampl] = t[, ZT.unlog]/t$mean;

	save(t,file = 'Rdata/raw_norm_adj_bg.Rdata')
	}



if(plot){
	if(load.for.plot){load('Rdata/raw_norm_adj_bg.Rdata')}
	source('functions.R')
	ZT = c(1:24)
	zt.rel = grep('.rel.ampl', colnames(t))

	pdf('plots/distib_raw_bg_removed_adj_values.pdf',width = 12, height = 16)
	par(mfcol = c(4,3))
	plot.distrib.t(t = t, zt = ZT, zt.rel = zt.rel, main = "unannotated removed - bg removed - adj bg values")
	dev.off()
}


#############################################################finishing removing the background

# remove genes with less than 2 probesets
if(my.do.step)
{
	load('raw_remove_background_quantile_0.25.Rdata')
	source('functions.R')
	ind = enough.probesets(t = tt); 
	tt$rejection.code[tt$kept & !ind] = 'not enough probesets'; 
	tt$kept = tt$kept & ind;
	save(tt,file = 'my_raw_enough_probesets.Rdata')
}

if(do.step){
	cat('REMOVE GENES WITH LESS THAN 2 PROBESETS \n')
	if(load){load('Rdata/raw_norm.Rdata')}
	source('functions.R')
	ind = enough.probesets(t = t); t$rejection.code[t$kept & !ind] = 'not enough probesets'; t$kept = t$kept & ind;
	save(t,file = 'Rdata/raw_enough_probesets.Rdata')
	}



if(plot){
	if(load.for.plot){load('Rdata/raw_enough_probesets.Rdata')}
	source('functions.R')
	ZT = c(1:24)
	zt.rel = grep('.rel.ampl', colnames(t))

	pdf('plots/distib_raw_enough_probesets.pdf',width = 12, height = 16)
	par(mfcol = c(4,3))
	plot.distrib.t(t = t, zt = ZT, zt.rel = zt.rel, main = "unannotated removed - enough probesets")
	dev.off()
}


# This has to be changed. The quantile normalisation should not be done at that step but before or not at all
###### Quantile Normalisation (Laura old parts of code)

if(do.step){
	cat('QUANTILE NORMALISATION \n')
	if(load){load('Rdata/raw_enough_probesets.Rdata')}
	
	library(preprocessCore)
	
	ZT.unlog = grep('.unlog', colnames(t))

	
	t.norm = normalize.quantiles(as.matrix(t[t$kept, ZT.unlog]))
	t[t$kept, ZT.unlog] = t.norm; rm(t.norm)
	
	ZT = 1:24; t[,ZT] = log2(t[, ZT.unlog])
	ZT.rel.ampl = grep('.rel.ampl', colnames(t))
	t$mean = apply(t[, ZT.unlog], 1, mean);
	t[, ZT.rel.ampl] = t[, ZT.unlog]/t$mean;
	
	save(t,file = 'Rdata/raw_quant_norm.Rdata')
	}


if(plot){
	if(load.for.plot){load('Rdata/raw_quant_norm.Rdata')}
	source('functions.R')
	ZT = c(1:24)
	zt.rel = grep('.rel.ampl', colnames(t))

	pdf('plots/distib_raw_quant_norm_all_together.pdf',width = 12, height = 16)
	par(mfcol = c(4,3))
	plot.distrib.t(t = t, zt = ZT, zt.rel = zt.rel, main = "unannotated removed - quant norm (together)")
	dev.off()
}





# do quantile normalisation separately for introns and exons
if(do.step){
	cat('QUANTILE NORMALISATION \n')
	if(load){load('Rdata/raw_enough_probesets.Rdata')}
	
	library(preprocessCore)
	
	ZT.unlog = grep('.unlog', colnames(t))

	i.ex = which(t$kept & (t$ex_in == 1))
	i.int = which(t$kept & (t$ex_in == -1))

	t.norm.ex = normalize.quantiles(as.matrix(t[i.ex, ZT.unlog]))
	t[i.ex, ZT.unlog] = t.norm.ex; rm(t.norm.ex)
	
	t.norm.int = normalize.quantiles(as.matrix(t[i.int, ZT.unlog]))
	t[i.int, ZT.unlog] = t.norm.int; rm(t.norm.int)
	
	ZT = 1:24; t[,ZT] = log2(t[, ZT.unlog])
	ZT.rel.ampl = grep('.rel.ampl', colnames(t))
	t$mean = apply(t[, ZT.unlog], 1, mean);
	t[, ZT.rel.ampl] = t[, ZT.unlog]/t$mean;
	
	save(t,file = 'Rdata/raw_quant_norm_sep_ex_int.Rdata')
	}


if(plot){
	if(load.for.plot){load('Rdata/raw_quant_norm_sep_ex_int.Rdata')}
	source('functions.R')
	ZT = c(1:24)
	zt.rel = grep('.rel.ampl', colnames(t))

	pdf('plots/distib_raw_quant_norm_sep_ex_int.pdf',width = 12, height = 16)
	par(mfcol = c(4,3))
	plot.distrib.t(t = t, zt = ZT, zt.rel = zt.rel, main = "unannotated removed - quant norm (sep ex int)")
	dev.off()
}




# do quantile normalisation separately for introns and exons on REL. AMPL.
if(do.step){
	cat('QUANTILE NORMALISATION \n')
	if(load){load('Rdata/raw_enough_probesets.Rdata')}
	
	library(preprocessCore)
	
	ZT.unlog = grep('.unlog', colnames(t))
	ZT.rel.ampl = grep('.rel.ampl', colnames(t))


	i.ex = which(t$kept & (t$ex_in == 1))
	i.int = which(t$kept & (t$ex_in == -1))

	t.norm.ex = normalize.quantiles(as.matrix(t[i.ex, ZT.rel.ampl]))
	t[i.ex, ZT.rel.ampl] = t.norm.ex; rm(t.norm.ex)
	
	t.norm.int = normalize.quantiles(as.matrix(t[i.int, ZT.rel.ampl]))
	t[i.int, ZT.rel.ampl] = t.norm.int; rm(t.norm.int)
	
	t[,ZT.unlog] = t[, ZT.rel.ampl]*t$mean
	ZT = 1:24; t[,ZT] = log2(t[, ZT.unlog])
	
	save(t,file = 'Rdata/raw_quant_norm_sep_ex_int_rel_ampl.Rdata')
	}


if(plot){
	if(load.for.plot){load('Rdata/raw_quant_norm_sep_ex_int_rel_ampl.Rdata')}
	source('functions.R')
	ZT = c(1:24)
	zt.rel = grep('.rel.ampl', colnames(t))

	pdf('plots/distib_raw_quant_norm_sep_ex_int_rel_ampl.pdf',width = 12, height = 16)
	par(mfcol = c(4,3))
	plot.distrib.t(t = t, zt = ZT, zt.rel = zt.rel, main = "unannotated removed - quant norm (sep ex int)_rel_ampl")
	dev.off()
}

################################################################################finish quantiles normalization


####################################
########## clean gene names

if(my.do.step)
{
	load('my_raw_enough_probesets.Rdata')
	
	RS = read.table('Raw_data_Annotation/RefSeq_Mouse_mm9.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
	tt$gene.A = tt$gene
	CHR = grep('_chr.', tt$gene)
	genes.chr.list = unique(sort(tt$gene[CHR]))
	
	not.in.RS = rep(FALSE, length(genes.chr.list)); 
	several.chr.in.RS = rep(FALSE, length(genes.chr.list))
	
	for(i  in 1:length(genes.chr.list))
	{
		gene_chr = genes.chr.list[i]
		cat(i, gene_chr, '\t')
		j = which(tt$gene == gene_chr)
		gene = unlist(strsplit(gene_chr,'_'))[1]; 
		chr = unlist(strsplit(gene_chr,'_'))[2]
		
		k = RS$name2 == gene
		if(sum(k)==0){cat('!!!!!!  ALERT : no gene with this name in REFSEQ \n'); not.in.RS[i] = TRUE;}
		else{
			chr.RS = unique(RS$chrom[k])
			if(length(chr.RS)!=1){cat('!!!!!!  ALERT : several chromosomes for the same gene name \n'); several.chr.in.RS[i] = TRUE}
			else{ if(chr == chr.RS){tt$gene[j] = gene}; cat('\n') }
		}
		
	}
	
	save(tt, file= 'my_raw_name_cleaned.Rdata')

	yy = tt[,c(1:61,111)]
	ZT.unlog = grep('.unlog', colnames(tt))
	yy[,c(1:24)] = tt[,ZT.unlog]
	tt = yy
	save(tt, file= 'my_raw_name_cleaned_linear_scale.Rdata')
	
}


if(do.step){
	if(load){load('Rdata/raw_quant_norm.Rdata')}
	RS = read.table('../complementary_data/RefSeq_Mouse_mm9.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
	
	
	
	# ######## T
	
	# T$gene.A = T$gene
	# CHR = grep('_chr.', T$gene)
	
	# not.in.RS = rep(FALSE, nrow(T)); several.chr.in.RS = rep(FALSE, nrow(T))
	# for(i in 1:length(CHR)){
		
		# j = CHR[i]
		# gene_chr = T$gene[j]; gene = unlist(strsplit(gene_chr,'_'))[1]; chr = unlist(strsplit(gene_chr,'_'))[2]
		# cat(i, gene_chr, '\t')
		
		# k = RS$name2 == gene
		# if(sum(k)==0){cat('!!!!!!  ALERT : no gene with this name in REFSEQ \n'); not.in.RS[j] = TRUE}
		# else{
			# chr.RS = unique(RS$chrom[k])
			# if(length(chr.RS)!=1){cat('!!!!!!  ALERT : several chromosomes for the same gene name \n'); several.chr.in.RS[j] = TRUE}
			# else{ if(chr == chr.RS){T$gene[j] = gene}; cat('\n') }
		# }
	# }
	
	######## t
	
	t$gene.A = t$gene
	CHR = grep('_chr.', t$gene)
	genes.chr.list = unique(sort(t$gene[CHR]))
	
	not.in.RS = rep(FALSE, length(genes.chr.list)); several.chr.in.RS = rep(FALSE, length(genes.chr.list))
	for(i  in 1:length(genes.chr.list)){
		gene_chr = genes.chr.list[i]
		cat(i, gene_chr, '\t')
		j = which(t$gene == gene_chr)
		gene = unlist(strsplit(gene_chr,'_'))[1]; chr = unlist(strsplit(gene_chr,'_'))[2]
		
		k = RS$name2 == gene
		if(sum(k)==0){cat('!!!!!!  ALERT : no gene with this name in REFSEQ \n'); not.in.RS[i] = TRUE}
		else{
			chr.RS = unique(RS$chrom[k])
			if(length(chr.RS)!=1){cat('!!!!!!  ALERT : several chromosomes for the same gene name \n'); several.chr.in.RS[i] = TRUE}
			else{ if(chr == chr.RS){t$gene[j] = gene}; cat('\n') }
		}
		
	}
	
	save(t, file=  'Rdata/raw_name_cleaned.Rdata')

}



# compute rel amplitude

if(my.do.step)
{
	cat('COMPUTE RELATIVE AMPLITUDE \n')
	
	load('my_raw_name_cleaned_linear_scale.Rdata')
	
	ZT = c(1:24)
	tt$mean = apply(tt[,ZT], 1, mean);
	tt_rel.ampl = tt[,ZT]/tt$mean; 
	colnames(tt_rel.ampl) = paste(colnames(tt)[ZT],'rel.ampl',sep = '.')
	tt = data.frame(tt, tt_rel.ampl, stringsAsFactors = FALSE)
	
	save(tt,file = 'my_raw_rel_ampl.Rdata')
}


if(do.step){
	cat('COMPUTE RELATIVE AMPLITUDE \n')
	if(load){load('Rdata/raw_name_cleaned.Rdata')}
	ZT = c(1:24)
	t$mean = apply(t[,ZT], 1, mean);
	t_rel.ampl = t[,ZT]/t$mean; colnames(t_rel.ampl) = paste(colnames(t)[ZT],'rel.ampl',sep = '.')
	t = data.frame(t, t_rel.ampl, stringsAsFactors = FALSE)
	save(t,file = 'Rdata/raw_rel_ampl.Rdata')
	}


####################################### signal selection and contruct the final table T

if(my.do.step)
{
	cat('MAKE SIGNAL SELECTION \n')
	
	load('my_raw_rel_ampl.Rdata')
	source('functions.R')

	#load('gene_verify_exon_intron.Rdata')
	
	#load('my_morning_genes_to_check')
	
	ZT = c(1:24); 
	ZTr = grep('.rel.ampl',colnames(tt))
	gene_list = unique(tt$gene[tt$kept]);
	

	#gene.clock = c("Per1","Per2","Per3","Cry1","Cry2","Dbp","Tef","Hlf","Nr1d1","Nr1d2","Rora","Rorb","Rorc","Arntl","Clock","Npas2","Bhlhe40","Bhlhe41","Cirbp","Tfrc","Hamp","Hamp2","Nr4a2", "Cdkn1a")
	#test.gene = c('Nedd4l', 'Cbs','Gsta3','Loxl4','Rcan1','Cdkn1a','Tubb2a','Mfsd2','Ppard', gene.clock, gene.check[1:20])
	
	
	## to test the signal selection function
	#test.gene = gene.very[1:100]
	
	#test.keep = c()
	#for(gene in test.gene)
	#{test.keep = c(test.keep, which(tt$gene==gene))}
	#test.table = tt[test.keep,]
	#gene_list = unique(test.table$gene[test.table$kept]);
	
	
	
	
	# 1st version is to normalize the data to have the same averges for two days
	#Tt = signal.selection.myv(t = tt, gene_list = gene_list); 
	#Tt = Tt$T; 
	
	#2nd version is just to select and then to use quantile normalization
	Tt = signal.selection.myv2(t = tt, gene_list = gene_list); 
	
	save(Tt,tt, file='my_Tt_tt_signal_selection_temp.Rdata')
	
	
	
	load('my_Tt_tt_signal_selection_temp.Rdata')
	
	ind.ex =  Tt$flag.ex; 
	ind.int = Tt$flag.int; 
	
	Tt = Tt$T
	
	library(preprocessCore)
	# quantile normalize the exon and intron after signal selections
	exon = as.matrix(Tt[,c(2:25)])
	exon = normalize.quantiles(exon)
	intron = as.matrix(Tt[,c(26:49)])
	intron = normalize.quantiles(intron)
	Tt[,c(2:25)] = exon
	Tt[,c(26:49)] = intron
	
	xx = as.matrix(Tt[,c(2:25)])
	cor(xx)
	#ss = apply(xx,1,mean)
	#xx = xx/ss
	#cor(xx)

	tt$rejection.code[tt$kept & (tt$ex_in == 1) & !ind.ex] = 'exon';
	tt$rejection.code[tt$kept & (tt$ex_in == -1) & !ind.int] = 'intron';
	tt$kept = tt$kept&(ind.ex|ind.int)
	save(tt,Tt,file= 'my_genes_rel_ampl_new.Rdata')

	load('my_genes_rel_ampl_new.Rdata')
	# remove gene witout signal
	cat('REMOVE GENES WITHOUT SUMMARIZED SIGNAL \n')

	source('functions.R')
	ind = remove.gene.witout.signal(T = Tt); 
	gene.to.remove = Tt$gene[!ind]; 
	ind.t = !is.na(match(tt$gene, gene.to.remove)); 
	tt$rejection.code[tt$kept & ind.t] = 'only 1 good probeset - exon'; 
	tt$kept[tt$kept & ind.t] = FALSE; 
	Tt = Tt[ind,]
	save(tt,Tt,file = 'my_genes_rel_ampl_clean_new.Rdata')
	
	xx = as.matrix(Tt[,c(2:25)])
	ss = apply(xx,1,mean)
	Tt$exon.median = ss
	xx = xx/ss
	Tt[,c(2:25)] = xx
	ZT.int = grep('.int', colnames(Tt))
	yy = as.matrix(Tt[,ZT.int])
	ss = apply(yy, 1, mean)
	Tt$intron.median = ss
	yy = yy/ss
	Tt[,ZT.int] = yy
	ZT.ex = grep('.ex', colnames(Tt))
	var.ex = apply(Tt[,ZT.ex], 1, sd)
	var.int = apply(Tt[,ZT.int], 1, sd)
	Tt$exon.var = var.ex
	Tt$intron.var = var.int
	
	xx = Tt
	ZT.ex = grep('.ex',colnames(xx))
	ZT.int = grep('.int',colnames(xx))
	colnames(xx)[ZT.ex] = colnames(T)[ZT.ex]
	colnames(xx)[ZT.int] = colnames(T)[ZT.int]
	#aa = seq(0,46, by=2)
	#aa[1:4] = c(00, 02, 04, 06, 08)
	#strsplit(colnames(xx)[ZT.ex],'.')
	#colnames(xx)[ZT.ex] = paste( , 'rel.ampl',sep = '.')
	save(tt,Tt,file = 'my_genes_rel_ampl_clean_new2.Rdata') ## relative signals and newly calculated averages and sd.
	
	###############################################
	############################################### Compare the distribution of absolute signals and relative signals for each time points
	###############################################
	
	pdf("gene_distribution_check_with_quantile.pdf",width=10,height=16)
	par(mfcol = c(3,2))
	
	#### exon log2 absolute value
	rainbow = rainbow(12,s = 0.85, v = 0.85)
	dd = log2(Ttt[,c(2:25)])
	d = density(dd[,zt[1]])
	zt = c(1:24)
	type='abs'
	if(type == 'abs'){xlab =  'abs. signal (log2)'; main = 'abs.signal (log2) Exon';xlim = c(0,15)}else{xlab =  'rel. signal'; xlim = c(0,2.2)}
	ylim = c(0,1.2*max(d$y))
	
	
	plot(1,1,type = 'n', xlab = xlab, ylab = 'density', xlim = xlim, ylim = ylim, main = main)
	abline(v = 1)
	
	for(i in zt){
		i.rel = i-zt[1]+1; 
		#cat(i.rel,"\n")
		d = density(dd[,i])
		points(d, type = 'l', col = rainbow[(i.rel-1)%%12+1], lty = (i.rel-1)%/%12+1)
		legend(x = 0.8*xlim[2],y = ylim[2]*((i.rel+1)/length(zt)), legend = paste('ZT',2*(i.rel-1)), col =  rainbow[(i.rel-1)%%12+1], lty = (i.rel-1)%/%12+1, bty = 'n' )
		#rainbow = rainbow(12,s = 0.85, v = 0.85)
	}
	
	### exon relative signal normalized together
	rainbow = rainbow(12,s = 0.85, v = 0.85)
	dd = Ttt[,c(2:25)]
	ss = apply(dd, 1, mean)
	dd = dd/ss
	d = density(dd[,zt[1]])
	zt = c(1:24)
	type='rel.signal'
	if(type == 'abs'){xlab =  'abs. signal (log2)'; main = 'abs.signal (log2)'; xlim = c(0,15)}else{xlab =  'rel. signal'; main = 'rel.signal.normalized.together';xlim = c(0,2.2)}
	ylim = c(0,1.5*max(d$y))
	
	plot(1,1,type = 'n', xlab = xlab, ylab = 'density', xlim = xlim, ylim = ylim, main = main)
	abline(v = 1)
	
	for(i in zt){
		i.rel = i-zt[1]+1; 
		#cat(i.rel,"\n")
		d = density(dd[,i])
		points(d, type = 'l', col = rainbow[(i.rel-1)%%12+1], lty = (i.rel-1)%/%12+1)
		legend(x = 0.8*xlim[2],y = ylim[2]*((i.rel+1)/length(zt)), legend = paste('ZT',2*(i.rel-1)), col =  rainbow[(i.rel-1)%%12+1], lty = (i.rel-1)%/%12+1, bty = 'n' )
	}
	
	### exon relative signal normalized seperately for two days
	rainbow = rainbow(12,s = 0.85, v = 0.85)
	dd = Ttt[,c(2:25)]
	ss1 = apply(dd[,c(1:12)], 1, mean)
	dd[,c(1:12)] = dd[,c(1:12)]/ss1
	ss2 = apply(dd[,c(13:24)], 1, mean)
	dd[,c(13:24)] = dd[,c(13:24)]/ss2

	d = density(dd[,zt[1]])
	zt = c(1:24)
	type='rel.signal'
	if(type == 'abs'){xlab =  'abs. signal (log2)'; main = 'abs.signal (log2)'; xlim = c(0,15)}else{xlab =  'rel. signal'; main = 'rel.signal.exon.normalized.seperately';xlim = c(0,2.2)}
	#ylim = c(0,1.2*max(d$y))
	
	plot(1,1,type = 'n', xlab = xlab, ylab = 'density', xlim = xlim, ylim = ylim, main = main)
	abline(v = 1)
	
	for(i in zt){
		i.rel = i-zt[1]+1; 
	#cat(i.rel,"\n")
		d = density(dd[,i])
		points(d, type = 'l', col = rainbow[(i.rel-1)%%12+1], lty = (i.rel-1)%/%12+1)
		legend(x = 0.8*xlim[2],y = ylim[2]*((i.rel+1)/length(zt)), legend = paste('ZT',2*(i.rel-1)), col =  rainbow[(i.rel-1)%%12+1], lty = (i.rel-1)%/%12+1, bty = 'n' )
	}
	
	
	
	### intron log2 absolute value
	rainbow = rainbow(12,s = 0.85, v = 0.85)
	dd = log2(Ttt[,c(26:49)])
	d = density(dd[,zt[1]])
	zt = c(1:24)
	type='abs'
	if(type == 'abs'){xlab =  'abs. signal (log2)'; main = 'abs.signal (log2) Intron'; xlim = c(0,15)}else{xlab =  'rel. signal'; main = 'rel.signal';xlim = c(0,2.2)}
	ylim = c(0,1.2*max(d$y))
	
	plot(1,1,type = 'n', xlab = xlab, ylab = 'density', xlim = xlim, ylim = ylim, main = main)
	abline(v = 1)
	
	for(i in zt){
		i.rel = i-zt[1]+1; 
		#cat(i.rel,"\n")
		d = density(dd[,i])
		points(d, type = 'l', col = rainbow[(i.rel-1)%%12+1], lty = (i.rel-1)%/%12+1)
		legend(x = 0.8*xlim[2],y = ylim[2]*((i.rel+1)/length(zt)), legend = paste('ZT',2*(i.rel-1)), col =  rainbow[(i.rel-1)%%12+1], lty = (i.rel-1)%/%12+1, bty = 'n' )
	}
	
	### intron relative signal normalized together
	rainbow = rainbow(12,s = 0.85, v = 0.85)
	dd = Ttt[,c(26:49)]
	ss = apply(dd, 1, mean)
	dd = dd/ss
	d = density(dd[,zt[1]])
	zt = c(1:24)
	type='rel.signal'
	if(type == 'abs'){xlab =  'abs. signal (log2)'; main = 'abs.signal (log2)'; xlim = c(0,15)}else{xlab =  'rel. signal'; main = 'rel.signal.normalized.together';xlim = c(0,2.2)}
	ylim = c(0,1.5*max(d$y))
	
	plot(1,1,type = 'n', xlab = xlab, ylab = 'density', xlim = xlim, ylim = ylim, main = main)
	abline(v = 1)
	
	for(i in zt){
		i.rel = i-zt[1]+1; 
		#cat(i.rel,"\n")
		d = density(dd[,i])
		points(d, type = 'l', col = rainbow[(i.rel-1)%%12+1], lty = (i.rel-1)%/%12+1)
		legend(x = 0.8*xlim[2],y = ylim[2]*((i.rel+1)/length(zt)), legend = paste('ZT',2*(i.rel-1)), col =  rainbow[(i.rel-1)%%12+1], lty = (i.rel-1)%/%12+1, bty = 'n' )
	}
	
	### exon relative signal normalized by two days seperately
	rainbow = rainbow(12,s = 0.85, v = 0.85)
	dd = Ttt[,c(26:49)]
	ss1 = apply(dd[,c(1:12)], 1, mean)
	dd[,c(1:12)] = dd[,c(1:12)]/ss1
	ss2 = apply(dd[,c(13:24)], 1, mean)
	dd[,c(13:24)] = dd[,c(13:24)]/ss2
	
	d = density(dd[,zt[1]])
	zt = c(1:24)
	type='rel.signal'
	if(type == 'abs'){xlab =  'abs. signal (log2)'; main = 'abs.signal (log2)'; xlim = c(0,15)}else{xlab =  'rel. signal'; main = 'rel.signal.intron.normalized.seperately';xlim = c(0,2.2)}
	#ylim = c(0,1.2*max(d$y))
	
	plot(1,1,type = 'n', xlab = xlab, ylab = 'density', xlim = xlim, ylim = ylim, main = main)
	abline(v = 1)
	
	for(i in zt){
		i.rel = i-zt[1]+1; 
	#cat(i.rel,"\n")
		d = density(dd[,i])
		points(d, type = 'l', col = rainbow[(i.rel-1)%%12+1], lty = (i.rel-1)%/%12+1)
		legend(x = 0.8*xlim[2],y = ylim[2]*((i.rel+1)/length(zt)), legend = paste('ZT',2*(i.rel-1)), col =  rainbow[(i.rel-1)%%12+1], lty = (i.rel-1)%/%12+1, bty = 'n' )
	}
	
	dev.off()
	
	
	
	# source('functions.R')
	# ZT = c(1:24)
	# zt.rel = grep('.rel.ampl', colnames(t))
	# pdf('plots/distib_raw_and_quant_norm_sep.pdf',width = 12, height = 16)
	# par(mfcol = c(4,3))
	# plot.distrib.t(t = t, zt = ZT, zt.rel = zt.rel, main = "just out from RMA - quant norm sep")
	# dev.off()
	
	
}



if(do.step){
	cat('MAKE SIGNAL SELECTION \n')
	if(load){load('Rdata/raw_rel_ampl.Rdata')}
	
	source('functions.R')
	ZT = c(1:24); ZTr = grep('.rel.ampl',colnames(t))
	
	gene_list = unique(t$gene[t$kept]);
	T = signal.selection(t = t, gene_list = gene_list); ind.ex =  T$flag.ex; ind.int = T$flag.int; T = T$T; 
	t$rejection.code[t$kept & (t$ex_in == 1) & !ind.ex] = 'exon'; t$rejection.code[t$kept & (t$ex_in == -1) & !ind.int] = 'intron'; t$kept = t$kept&(ind.ex|ind.int)
	save(t,T,file= 'Rdata/genes_rel_ampl.Rdata' )
	}

	
# remove gene witout signal
if(do.step){
	cat('REMOVE GENES WITHOUT SUMMARIZED SIGNAL \n')
	if(load){load('Rdata/genes_rel_ampl.Rdata')}
	source('functions.R')
	ind = remove.gene.witout.signal(T = T); gene.to.remove = T$gene[!ind]; ind.t = !is.na(match(t$gene, gene.to.remove)); t$rejection.code[t$kept & ind.t] = 'only 1 good probeset - exon'; t$kept[t$kept & ind.t] = FALSE; T = T[ind,]
	save(t,T,file = 'Rdata/genes_rel_ampl_light.Rdata')
}
	
	
	
if(plot){
	if(load.for.plot){load('Rdata/genes_rel_ampl_light.Rdata')}
	
	ZT.ex = c(2:25); ZT.int = ZT.ex+24
	ZT = c(1:24); ZTr = grep('.rel.ampl',colnames(t))
	zt = seq(0,46,by=2)

	circ = c("Per1","Per2","Per3","Cry1","Cry2","Dbp","Tef","Hlf","Nr1d1","Nr1d2","Rora","Rorb","Rorc","Arntl","Clock","Npas2","Bhlhe40","Bhlhe41","Cirbp","Tfrc","Hamp","Hamp2","Nr4a2", "Cdkn1a")
	gene_list = c(circ, T$gene[1:10])
	
	source('functions.R')
	pdf('plots/selection_selected_genes.pdf', width = 20, height =5)
	par(mfrow = c(1,4))
	for(gene in gene_list){
		
		cat(gene, '\n')
		plot.selection(gene = gene, t = t, T = T, probe.type = 'intron')
		plot.selection(gene = gene, t = t, T = T, probe.type = 'exon')
	}
	dev.off()
}		
		





# do rhythmic analysis of probesets and check if the probeset selection make sense.
if(my.do.step)
{
	cat('DO RHYTHMIC ANALYSIS OF PROBESETS \n')
	load('my_genes_rel_ampl_clean_new2.Rdata')
	source("f24_modified_1.0.r")
	
	#ii = which(tt$kept=="TRUE"&tt$ex_in==1)
	#jj = which(tt$kept=="TRUE"&tt$ex_in==-1)
	
	#load('Rdata/genes_rel_ampl_light.Rdata')
	
	
	#intron = as.matrix(Tt[,c(1:24)])
`	
	#library(preprocessCore)
	#exon = as.matrix(Tt[,c(2:25)])
	#exon = normalize.quantiles(exon)
	
	#exon = as.matrix(Tt[,c(2:25)])
	#exon = log2(exon)
	boxplot(exon)
	abline(h=mean(exon[,1]),col='red',lwd=2.0)
	
	#intron = as.matrix(Tt[,c(26:49)])
	#intron = normalize.quantiles(intron)
	#ss = apply(intron,1,mean)
	intron = log2(intron)
	boxplot(intron)
	abline(h=mean(intron[,1]),col='red',lwd=2.0)
	
	
	exon = as.matrix(Tt[,c(2:25)])
	exon = log2(exon)
	intron = as.matrix(Tt[,c(26:49)])
	intron = log2(intron)
	res.exon = t(apply((exon),1, f24_R2_alt2, t=c(0:23)*2))
	res.intron = t(apply((intron),1, f24_R2_alt2, t=c(0:23)*2))
	res.exon = cbind(res.exon, qv=qvals(res.exon[,6]))
	res.intron = cbind(res.intron, qv=qvals(res.intron[,6]))
	
	
	
	kk = which(res.exon[,7]<0.05)
	length(kk)
	#plot(res.exon[kk,5], res.exon[kk,4])
	#abline(h=0.05,col='red',lwd=2.0)
	
	length(which(res.exon[,7]<0.05))
	hist(res.exon[which(res.exon[,7]<0.05),5],breaks=24)
	
	hist(res.intron[which(res.intron[,7]<0.05),5],breaks=24)
	
	hist(res.exon[which(res.exon[,7]<0.05 & res.exon[,3]>=0.05),5])
	
	## check the difference between this version with the Laura old version
	load('gene_verify_exon_intron.Rdata')
	gene.check = gene.very
	gene.clock = c("Per1","Per2","Per3","Cry1","Cry2","Dbp","Tef","Hlf","Nr1d1","Nr1d2","Rora","Rorb","Rorc","Arntl","Clock","Npas2","Bhlhe40","Bhlhe41","Cirbp","Tfrc","Hamp","Hamp2","Nr4a2", "Cdkn1a")
	test.gene = c('Nedd4l', 'Cbs','Gsta3','Loxl4','Rcan1','Cdkn1a','Tubb2a','Mfsd2','Ppard', gene.clock, gene.check[1:100])
	
	
	
	pdf("gene_compare_Laura.pdf",width=6,height=4)
	for(n in 1:length(test.gene))
	#for(n in 1:10)
	{
		#par(mfrow = c(1,2))
		genename = test.gene[n]
		ii = which(Tt[,1]==genename)
		if(length(ii)>0){
		test = as.numeric(Tt[ii, c(2:25)])
		test2 = as.numeric(Tt[ii, c(26:49)])
		#genename = Tt[check[n],1]
		lims = range(c(test, test2))
		plot(c(0:23)*2, test,type='b',col='darkblue',lwd=2.0,cex=0.7,main=genename,ylim=c(0,2.5))
		points(c(0:23)*2, test2,type='b',col='darkgreen',lwd=2.0,cex=0.7)
		abline(h=1,col='gray',lwd=2.0)
		}
	}		
	dev.off()
	
	### compare the mrna profiles from Exon array and RNA-seq
	
	mrna = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/mRNA_RNA_Seq_all.txt')
	
	length(which(mrna$Q.RNA.Seq<0.05 & mrna$relamp>0.05))
	hist(mrna$phase[which(mrna$Q.RNA.Seq<0.05 & mrna$relamp>0.0)],breaks=24)
	
	kk = which(mrna$Q.RNA.Seq<0.05 & mrna$relamp>=0.02)
	rownames(res.exon) = Tt$gene
	kk.exon = which(res.exon[,7]<0.05)
	
	jj = match(rownames(mrna)[kk],rownames(res.exon))
	kk1 = kk[which(!is.na(jj)==T)]
	kk1.exon = jj[which(!is.na(jj)==T)]
	
	ii = match(rownames(res.exon)[kk.exon], rownames(mrna)[kk])
	kk.exon = kk.exon[which(!is.na(ii)==T)]
	kk = kk[ii[which(!is.na(ii)==T)]]
	#hist(mrna[kk,]$phase, breaks=24)
	plot(mrna[kk1,]$phase, res.exon[kk1.exon,5],col='gray',xlab='mrna phase from RNA-seq', ylab='mrna phase from Exon array')
	points(mrna[kk,]$phase, res.exon[kk.exon,5],col='orange')
	abline(0,1,col='red', lwd=2.0)
	
	gene.overlap = rownames(res.exon)[kk1.exon]
	gene.overlap.ryth = rownames(res.exon)[kk.exon]
	check = gene.overlap[which(is.na(match(gene.overlap, gene.overlap.ryth))==T)]
	check = match(check, rownames(res.exon))
	verify = match(rownames(res.exon)[check], rownames(mrna))
	gene.check = Tt$gene[check]
	
	save(gene.check, file='my_morning_genes_to_check.Rdata')
	cch = match(gene.check, tt$gene.A)
	cch2 = check
	
	######### check the difference between two replicats
	ZT=2
	plot(tt[,(ZT)], tt[,(ZT+12)],log='xy')
	points(tt[cch,(ZT)], tt[cch,(ZT+12)], col='blue')
	abline(0,1,col='red',lwd=2.0)
	
	ZT=3
	plot(Tt[,(ZT)], Tt[,(ZT+1)],log='xy')
	points(Tt[cch2,(ZT)], Tt[cch2,(ZT+12)], col='blue')
	abline(0,1,col='red',lwd=2.0)
	
	par(mfrow = c(1,2))
	hist(res.exon[check,5],breaks=24)
	hist(mrna[verify,]$phase, breaks=24,col='gray')
	
	pdf("gene_check.pdf",width=15,height=8)
	for(n in 1:length(check))
	#for(n in 1:10)
	{
		par(mfrow = c(1,2))
		test = as.numeric(Tt[check[n], c(2:25)])
		genename = Tt[check[n],1]
		plot(c(0:23)*2, test,type='b',col='blue',lwd=2.0,main=paste("Exon Aarry: ",genename, " phase= ", signif(res.exon[check[n],5],d=2), " Q= ", signif(res.exon[check[n],7],d=3)))
		abline(v=24,col='red',lwd=2.0)
		abline(v=22,col='red',lwd=2.0)
		abline(h=mean(test[c(1:12)], col='black', lwd=2.0))
		abline(h=mean(test[c(13:24)], col='black', lwd=2.0))
		#abline(v=20,col='red',lwd=2.0)
		very = as.numeric(mrna[verify[n], c(1:24)])
		plot(c(0:23)*2, very,type='b',col='darkblue',lwd=2.0,main=paste("RNA-seq:",rownames(mrna)[verify[n]], " phase=", signif(mrna$phase[verify[n]],d=2), " Q= ", signif(mrna$Q.RNA.Seq[verify[n]],d=3)))
		abline(v=24,col='red',lwd=2.0)
		abline(v=22,col='red',lwd=2.0)
		abline(h=mean(very[c(1:12)], col='black', lwd=2.0))
		abline(h=mean(very[c(13:24)], col='black', lwd=2.0))
	}
	dev.off()
	
	
}

if(my.do.step)
{
	cat('DO RHYTHMIC ANALYSIS OF PROBESETS \n')
	load = TRUE
	if(load){load('my_genes_rel_ampl_clean_new2.Rdata')}
	
	ZT = c(1:24); ZTr = grep('rel.ampl',colnames(tt))
	source('functions.R')
	
	tt$phase[tt$kept] = apply(tt[tt$kept,ZTr],1,phase.n, n = 2)/2/pi*24
	tt$pval[tt$kept] = apply(tt[tt$kept,ZTr],1,pval.n, n = 2)
	save(tt,Tt, file = 'my_genes_RAnalysis_probesets.Rdata')
}

if(do.step){
	cat('DO RHYTHMIC ANALYSIS OF PROBESETS \n')
	if(load){load('Rdata/genes_rel_ampl_light.Rdata')}
	
	ZT = c(1:24); ZTr = grep('rel.ampl',colnames(t))
	source('functions.R')
	
	t$phase[t$kept] = apply(t[t$kept,ZTr],1,phase.n, n = 2)/2/pi*24
	t$pval[t$kept] = apply(t[t$kept,ZTr],1,pval.n, n = 2)
	save(t,T, file = 'Rdata/genes_RAnalysis_probesets.Rdata')
	}


############# BED FILES FOR THE VISUALISATION ON UCSC GENOME BROWSER

if(do.step){
	cat('MAKE BED FILES \n')
	if(load.for.plots){load('Rdata/gene_RAnalysis_probesets.Rdata')}
	source('functions.R')
	ZT = c(1:24); ZTr = grep('rel.ampl',colnames(t))
	make.beds()
}



############# statistics on accepted and rejected probesets

if(plot){
	if(load.for.plot){load('Rdata/gene_RAnalysis_probesets.Rdata')}
	rejection.list = unique(t$rejection.code); rejection.list = rejection.list[c(1,7,3,6,4,5,8,2)]
	rejection.counts = vector(mode = 'numeric', len = length(rejection.list))
	for(i in 1:length(rejection.list)){rejection.counts[i] = sum(t$rejection.code == rejection.list[i])}

	pdf('plots/statistics_on_accepted_and_rejected_probesets.pdf')
	pie(rejection.counts, labels = rejection.list, col = c('orange','darkorange','orangered','red2','red3','red4','black','chartreuse3'), init.angle = 90, density = c(rep(-1,7),30), angle = 12, border = NA)
	dev.off()
}

############# How many accepted probeset per gene ?

if(do.step){
	cat('STATS : HOW MANY ACCEPTED PROBESETS PER GENE\n')
	if(load){load('Rdata/genes_RAnalysis_probesets.Rdata');}
	numbers = compute.statistics.on.probeset.numbers(); numbers = data.frame(T$gene, numbers, stringsAsFactors = FALSE)
	save(numbers, file = 'Rdata/statistics_kept_and_rejected_probesets_per_gene.Rdata')
}

if(plot){
	if(load.for.plot){load('Rdata/genes_RAnalysis_probesets.Rdata'); load('Rdata/statistics_kept_and_rejected_probesets_per_gene.Rdata')}
	pdf('plots/statistics_on_probeset_numbers_per_genes.pdf', width = 3.3, height = 2.5)
	par(cex = 0.5, mar = c(2,3,1,0.2))

	boxplot(numbers[,2:8]+1, pch = 16, cex = 0.5, log = 'y', axes = FALSE, col = c('seagreen1','yellowgreen','darkolivegreen1','cadetblue1','lavender','lightskyblue1','gray70'), border = c('seagreen4','green4','green4','turquoise4','slateblue','steelblue4','gray40'))
	box()
	axis(1, at = c(1:7),labels = gsub('\\.','\n',colnames(numbers[,2:8])), tick = FALSE)
	labels =c(0,1,2,5,10,20,50,100,200,500); axis(2, at = labels+1, labels = labels, las = 1)
	abline(v = c(3.5,6.5), col =  'black')

	boxplot(numbers[,2:8], pch = 16, cex = 0.5, ylim = c(0,90), axes = FALSE, col = c('seagreen3','yellowgreen','green3','turquoise3','mediumslateblue','steelblue','gray'))
	box()
	axis(1, at = c(1:7),labels = gsub('\\.','\n',colnames(numbers[,2:8])), tick = FALSE)
	axis(2, las = 1, at = seq(0,100,by = 10))
	abline(v = c(3.5,6.5), col =  'black')

	dev.off()


	RS = read.table('../complementary_data/RefSeq_mouse_mm9.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
	RS$length = RS$txEnd- RS$txStart


	pdf('plots/number_of_intronic_probesets_correlates_with_gene_length_number_of_exonic_probesets_correlate_with_number_of_exons.pdf', width = 20, height = 10)
	par(mfrow = c(2,4))
	j = match(numbers$T.gene, RS$name2)
	plot(numbers$total.probeset, RS$length[j], pch = 16, cex = 0.4, log = 'xy', main = 'ALL probesets', xlab = 'number of probeset per gene', ylab = 'gene length')
	plot(numbers$intronic.probeset, RS$length[j], pch = 16, cex = 0.4, log = 'xy', main = 'INTRONS', xlab = 'number of intronic probeset per gene', ylab = 'gene length')
	plot(numbers$exonic.probeset, RS$length[j], pch = 16, cex = 0.4, log = 'xy', main = 'EXONS', xlab = 'number of exonic probeset per gene', ylab = 'gene length')
	plot(numbers$exonic.probeset, RS$exonCount[j], pch = 16, cex = 0.4, log = 'xy', main = 'EXONS', xlab = 'number of exonic probeset per gene', ylab = 'number of exons')
	abline(a = 0, b = 1, col = 'green3')

	j = match(numbers$T.gene, RS$name2)
	plot(numbers$kept.intronic + numbers$kept.exonic, RS$length[j], pch = 16, cex = 0.4, log = 'xy', main = 'ALL probesets', xlab = 'number of KEPT probeset per gene', ylab = 'gene length')
	plot(numbers$kept.intronic, RS$length[j], pch = 16, cex = 0.4, log = 'xy', main = 'INTRONS', xlab = 'number of intronic KEPT probeset per gene', ylab = 'gene length')
	plot(numbers$kept.exonic, RS$length[j], pch = 16, cex = 0.4, log = 'xy', main = 'EXONS', xlab = 'number of exonic KEPT probeset per gene', ylab = 'gene length')
	plot(numbers$kept.exonic, RS$exonCount[j], pch = 16, cex = 0.4, log = 'xy', main = 'EXONS', xlab = 'number of exonic KEPT probeset per gene', ylab = 'number of exons')
	abline(a = 0, b = 1, col = 'green3')

	dev.off()
}




############### Fold change and Lag between Intron and Exons.



if(plot){
	if(load.for.plot){load('Rdata/genes_RAnalysis_probesets.Rdata')}
	ZT.ex = grep('.rel.ampl.ex', colnames(T)); ZT.int = grep('.rel.ampl.int', colnames(T));

	pdf('plots/unbiaised_analysis.pdf', height = 4.4, width = 6.6)
	
	mar1 = c(3,1,2,0)+0.1
	mar2 = c(3,0,2,0)+0.1
	mar3 = c(3,3,2,0)+0.1

	par(mar = mar1, cex = 0.7, mgp = c(2,0.6,0), las = 1, tcl = -0.3) 
	layout(matrix(1:6,2,3,byrow= TRUE),widths = c(1,0.8,1))

	### difference of absolute levels between introns and exons.
	par(pty = 's')
	color.function = colorRampPalette(c('gray80', 'black'))
	colors = densCols(log(cbind(T$exon.median, T$intron.median)), colramp = color.function, nbin = 500)#
	plot(log2(T$intron.median), log2(T$exon.median), cex = 0.5, pch = 16, col = colors,ylab = 'Exons - Expression level (log2)' , xlab = 'Introns - expression level (log2)', asp = 1, xlim = c(0,14.5), ylim = c(0,14.5));
	abline(a= 0, b = 1, col = 'slategray1', lwd = 3)
	mtext('A',side = 3, outer = FALSE, line =0.5, at = -3, font = 2)



	#par(pty = 'm')
	#boxplot(log2(T$intron.median), log2(T$exon.median), cex = 0.5, pch = 16, col = 'gray80', border = 'gray50', names = c('introns','exons'), ylab = 'Intensity (log2)');
	
	par(pty = 'm', mar = mar2)
	
	diff = log2(T$exon.median)-log2(T$intron.median)
	maxi = ceiling(max(diff,-diff))
	len = 23
	hist(diff, breaks = seq(-maxi, maxi, len= len), col = rep(c('gray80','gray70'),each = (len-1)/2), border= NA, main = '', xlab = ' Difference in expression', ylab = '', axes = FALSE, ylim = c(-0.05,0.25), freq = FALSE);  #abline(v = c(-1,1), lwd = 0.5, col = 'gray');
	start = 3*floor(maxi/3)
	axis(1,at = seq(-start, start,by = floor(maxi/3)))
	mtext('B',side = 3, outer = FALSE, line =0.5, at = -9, font = 2)
	boxplot(diff, pch = 16, add = TRUE, horizontal = TRUE, at = -0.025, boxwex = c(0.05), col = 'gray', border = 'gray50', axes = FALSE, cex = 0.5)
	abline(v = 0, lwd = 3, col = 'slategray1')

	par(pty = 's', mar = mar1)
	
	color.function = colorRampPalette(c('gray80', 'black'))
	colors = densCols(cbind(diff,log(T$exon.median)), colramp = color.function, nbin = 500)#
	plot(log2(T$exon.median),diff, cex = 0.5, pch = 16, col = colors,ylab = 'Difference in expression level' , xlab = 'Exons - expression level (log2)', asp = 1, xlim = c(0,14.5), ylim = c(-maxi,maxi));
	abline(h = 0, col = 'slategray1', lwd = 3)
	mtext('C',side = 3, outer = FALSE, line =0.5, at = -7, font = 2)
	
	# color.function = colorRampPalette(c('gray80', 'black'))
	# colors = densCols(cbind(diff,log(T$intron.median)), colramp = color.function, nbin = 500)#
	# plot(log2(T$intron.median),diff, cex = 0.5, pch = 16, col = colors,ylab = 'Difference in expression level' , xlab = 'introns - expression level (log2)', asp = 1, xlim = c(0,14.5), ylim = c(-maxi,maxi));
	# abline(h = 0, col = 'slategray1', lwd = 3)

	### Fold Change and Delay Analysis
	FC.ex = apply(T[,ZT.ex],1, max) - apply(T[,ZT.ex],1, min)
	FC.int = apply(T[,ZT.int],1, max) - apply(T[,ZT.int],1, min)

	lim = 0.5; j = (FC.int> lim)&(FC.ex> lim)&(T$exon.median > 2^5) #& (T$intron.median > 2^4) ; #T = T[j,]; FC.int = FC.int[j]; FC.ex =FC.ex[j]
	
	par(pty = 's')
	color.function = colorRampPalette(c('gray80', 'black'))
	colors = densCols(log(cbind(FC.ex[j], FC.int[j])), colramp = color.function, nbin = 500)#
	plot( FC.ex[j], FC.int[j], cex = 0.7, pch = 16, col = colors,xlab = 'Exons - Fold change' , ylab = 'Introns - Fold change', asp = 1, xlim = c(0.5,4.7), ylim = c(0.5,4.7), log = 'xy');abline(a= 0, b = 1, col = 'slategray1', lwd = 3)
	mtext('D',side = 3, outer = FALSE, line =0.5, at = 0.3, font = 2)

	par(pty = 'm', mar = mar2)

	ratio = log2(FC.int[j]/FC.ex[j])
	maxi = ceiling(max(ratio,-ratio))
	hist(ratio, breaks = seq(-maxi,maxi, len = 19), col = rep(c('gray80','gray70'),each = 9), border= NA, main = '', xlab = 'Fold change ratio', ylab = '', axes = FALSE, ylim = c(-0.1,0.7), freq = FALSE);  #abline(v = c(-1,1), lwd = 0.5, col = 'gray');
	axis(1,at = seq(-maxi,maxi,by = 1))

	boxplot(ratio, pch = 16, add = TRUE, horizontal = TRUE, at = -0.05, boxwex = c(0.1), col = 'gray', border = 'gray50', axes = FALSE, cex = 0.5)
	abline(v = 0, lwd = 3, col = 'slategray1')
	mtext('E',side = 3, outer = FALSE, line =0.5, at = -3, font = 2)

	par(mar = mar3)
	
	t.max.ex = (apply(T[j,ZT.ex], 1, which.max)-1)*2; t.max.int = (apply(T[j,ZT.int], 1, which.max)-1)*2; diff = (t.max.ex-t.max.int)%%24
	hist(diff, breaks = seq(0,24,by = 2), col = c('gray70'), border = NA, axes = FALSE, freq = TRUE, main = '', xlab = 'Delay between introns and exons [h]', ylab = '')
	axis(1,at = seq(0,24,by = 2))
	axis(2, las = 1)
	mtext('F',side = 3, outer = FALSE, line =0.5, at = -3, font = 2)

	dev.off()
}



if(plot){
	if(load.for.plot){load('Rdata/genes_RAnalysis_probesets.Rdata')}
	ZT.ex = grep('.rel.ampl.ex', colnames(T)); ZT.int = grep('.rel.ampl.int', colnames(T));

	pdf('plots/unbiaised_analysis_for_publication.pdf', height = 4.4, width = 6.6)
	
	mar1 = c(3,1,2,0)+0.1
	mar2 = c(3,0,2,0)+0.1
	mar3 = c(3,3,2,0)+0.1

	par(mar = mar1, cex = 0.7, mgp = c(2,0.6,0), las = 1, tcl = -0.3) 
	layout(matrix(1:6,2,3,byrow= TRUE),widths = c(1,0.8,1))


	par(pty = 'm', mar = mar2)
	
	diff = log2(T$exon.median)-log2(T$intron.median)
	maxi = ceiling(max(diff,-diff))
	len = 23
	hist(diff, breaks = seq(-maxi, maxi, len= len), col = rep(c('gray80','gray70'),each = (len-1)/2), border= NA, main = '', xlab = 'log2(Ratio exp. Exon/Intron)', ylab = '', axes = FALSE, ylim = c(-0.05,0.25), freq = FALSE);  #abline(v = c(-1,1), lwd = 0.5, col = 'gray');
	start = 3*floor(maxi/3)
	axis(1,at = seq(-start, start,by = floor(maxi/3)))
	boxplot(diff, pch = 16, add = TRUE, horizontal = TRUE, at = -0.025, boxwex = c(0.05), col = 'gray', border = 'gray50', axes = FALSE, cex = 0.5)
	abline(v = 0, lwd = 3, col = 'slategray1')

	par(pty = 's', mar = mar1)
	
	
	### Fold Change and Delay Analysis
	FC.ex = apply(T[,ZT.ex],1, max) - apply(T[,ZT.ex],1, min)
	FC.int = apply(T[,ZT.int],1, max) - apply(T[,ZT.int],1, min)

	f24.ex = apply(T[,ZT.ex], 1, f24_R2); f24.ex = t(f24.ex)
	f24.int = apply(T[,ZT.int], 1, f24_R2); f24.int = t(f24.int)

	plim.ex = fdr2pval(f24.ex[,6])
	plim.int = fdr2pval(f24.int[,6])
	ind = which((f24.ex[,6]< plim.ex)&(f24.int[,6]<plim.int))
	ind.l = which((f24.ex[,6]< plim.ex)|(f24.int[,6]<plim.int))

	ratio = f24.ex[,4]/f24.int[,4];
	phase.diff = f24.ex[,5]-f24.int[,5]; phase.diff = phase.diff%%24; phase.diff[phase.diff>18] = phase.diff[phase.diff>18]-24
	
	plot(phase.diff[ind] ,ratio[ind],pch = 16, cex = 0.5, col = 'gray', log = 'y', ylab = 'rel.ampl. ratio [exon/intron]', xlab ='phase diff. [exon-intron] (h)', axes = FALSE)
	x = seq(0,6,by = 0.1)
	points(x,cos(x*2*pi/24), type = 'l', col = 'slategray1', lwd = 3)
	axis(1, at = seq(-6,18,by = 6)); axis(2); box()
	points(phase.diff[ind] ,ratio[int],pch = 16, cex = 0.5)
	
	breaks = seq(-1.5,1.5,by = 0.03)
	hist(log10(ratio), breaks = breaks, axes = FALSE, xlim = c(-1.5,1.5), col = 'gray95', border = NA)
	hist(log10(ratio[ind.l]), breaks = breaks, add = TRUE, col = 'gray', border = NA)
	hist(log10(ratio[ind]), breaks = breaks, add = TRUE, col = 'black', border = NA)
	abline(v = 0, col = 'slategray1', lwd = 3 )
	axis(1, at =seq(-1,1,by = 0.5), labels = 10^(seq(-1,1,by = 0.5))); axis(2)
	
	breaks = seq(floor(min(phase.diff)),ceiling(max(phase.diff)),by = 0.5)
	hist(phase.diff, breaks = breaks, axes = FALSE,  col = 'gray95', border = NA)
	hist(phase.diff[ind.l], breaks = breaks, add = TRUE, col = 'gray', border = NA)
	hist(phase.diff[ind], breaks = breaks, add = TRUE, col = 'black', border = NA)
	abline(v = 0, col = 'slategray1', lwd = 3 )
	axis(1, at =seq(-24,24,by = 6)); axis(2)


	dev.off()
}


#################################################################################################################################################################################################################################
######################### finishing line of processing raw data
#################################################################################################################################################################################################################################

# table t is not modified anymore beyond this point.
###################################
###################################
###################################

if(my.do.step)
{
	cat('RHYTHMIC ANALYSIS OF SUMMARIZED SIGNAL\n')
	load()
	if(load){load('my_genes_RAnalysis_probesets.Rdata')}
	
	ZT.ex = grep('.rel.ampl.ex', colnames(T)); 
	ZT.int = grep('.rel.ampl.int', colnames(T));
	source('functions.R')
	
	T.r.ex = do.rhythmic.analysis(T = Tt, ZT = ZT.ex); 
	colnames(T.r.ex) = paste(colnames(T.r.ex), '.ex',sep = '')
	T.r.int = do.rhythmic.analysis(T = Tt, ZT = ZT.int); 
	colnames(T.r.int) = paste(colnames(T.r.int), '.int',sep = '')
	Tt = data.frame(Tt, T.r.ex, T.r.int, stringsAsFactors = FALSE)
	
	save(tt,Tt, file = 'my_genes_RAnalysis.Rdata')
	## check the noise in exon and intron
	ii = which(Tt$pval.ex>0.05)
	hist(Tt$exon.var[ii],breaks=50)
	jj= which(Tt$pval.int>0.05)
	hist(Tt$intron.var[jj],breaks=50)
	
	load('my_genes_RAnalysis.Rdata')
	ii = which(Tt$exon.median>Tt$intron.median)
	Tt = Tt[ii,]
	save(tt,Tt, file = 'my_genes_RAnalysis_larger_exons.Rdata')
}



if(my.do.step)
{
	cat('RHYTHMIC GENES SELECTION \n')
	load = TRUE
	if(load){load('my_genes_RAnalysis_larger_exons_2.Rdata')}
	
	ZT.ex = c(2:25)
	ZT.int = c(26:49)
	ex.var = apply(as.matrix(Tt[,ZT.ex]), 1, sd)
	int.var = apply(as.matrix(Tt[,ZT.int]), 1, sd)
	#Tt$exon.var = ex.var
	#Tt$intron.var = int.var
	sd.ex = mean(ex.var[which(Tt$pval.int>0.5)])/2
	sd.int = mean(int.var[which(Tt$pval.int>0.5)])/2
	
	FDR = 0.1; 
	#ampl.lim = 0.2; 
	fact = 0
	
	T = Tt
	
	N = nrow(T)
	source('functions.R')
	T$rhythmic.ex.24 = select.rhythmic.genes(probe.type = 'exon', period = 24, fdr = FDR, ampl.lim=sd.ex,fact = fact)
	T$rhythmic.int.24 = select.rhythmic.genes(probe.type = 'intron', period = 24, fdr = FDR, ampl.lim=sd.int,fact = fact)
	T$rhythmic.ex.12 = select.rhythmic.genes(probe.type = 'exon', period = 12, fdr = FDR, ampl.lim=sd.ex, fact = fact)
	T$rhythmic.int.12 = select.rhythmic.genes(probe.type = 'intron', period = 12, fdr = FDR, ampl.lim=sd.int, fact = fact)
	
	Tt = T; 
	
	sum(Tt$rhythmic.ex.24);sum(Tt$rhythmic.int.24)
	
	save(tt,Tt, file = 'my_genes_RAnalysis_larger_exons_2_rhythmic_analysis.Rdata')
}



#### PHASE HISTOGRAM
if(my.plot)
{
	load.for.plot =TRUE
	if(load.for.plot){load('my_genes_RAnalysis_larger_exons_2_rhythmic_analysis.Rdata')}
	
	plot.phase.histogram = function(probe.type = 'exon', period = 24)
	{
		if(probe.type == 'exon'){
			if(period == 24)
			{
				phase = Tt$phase.ex; 
				j = Tt$rhythmic.ex.24
			}else{
				phase = Tt$phase.12.ex; 
				j = Tt$rhythmic.ex.12
			}; 
			col = 'steelblue3'
		}else{if(period == 24){phase = Tt$phase.int; j = Tt$rhythmic.int.24}else{phase = Tt$phase.12.int; j = Tt$rhythmic.int.12}; col = 'green3'}
		
		if(period == 24){breaks = c(0:24); at.axis = seq(0,24,by = 4)}else{breaks = c(0:12); at.axis = seq(0,12,by = 3)}
		h = hist(phase[j], breaks = breaks, border = 'white', col = col, xlab = 'time (hours)', ylab = 'number of genes', main = paste('Histogram of phases -', probe.type,'- period :',period,'h'), axes = FALSE)
		axis(1, at = at.axis)
		axis(2, las = 1)
	}
	
	
	pdf('myplots/phases_histograms.pdf', width = 5, height = 7)
	par(mfrow = c(2,1))
	plot.phase.histogram(probe.type = 'intron', period = 24)
	plot.phase.histogram(probe.type = 'exon', period = 24)
	dev.off()
}



if(do.step){
	cat('RHYTHMIC ANALYSIS OF SUMMARIZED SIGNAL\n')
	if(load){load('Rdata/genes_RAnalysis_probesets.Rdata')}
	
	ZT.ex = grep('.rel.ampl.ex', colnames(T)); ZT.int = grep('.rel.ampl.int', colnames(T));
	source('functions.R')
	
	T.r.ex = do.rhythmic.analysis(T = T, ZT = ZT.ex); colnames(T.r.ex) = paste(colnames(T.r.ex), '.ex',sep = '')
	T.r.int = do.rhythmic.analysis(T = T, ZT = ZT.int); colnames(T.r.int) = paste(colnames(T.r.int), '.int',sep = '')
	T = data.frame(T, T.r.ex, T.r.int, stringsAsFactors = FALSE)

	save(t,T, file = 'Rdata/genes_RAnalysis.Rdata')
	}

save(t, file = 'Table_t_allprobset.R')
save(T, file = 'Table_T_alldata.R')
#write.table(t,'Table_t_allprobsets.txt',sep='\t',quote=F,col.names=TRUE, row.names=F)
#write.table(T,'Table_T_alldata.txt',sep='\t',quote=F,col.names=TRUE, row.names=F)

source("f24_modified_1.0.r")
exon = as.matrix(T[,c(2:25)])
intron = as.matrix(T[,c(26:49)])
res.exon = (apply(exon,1, f24_R2_alt2, t=c(0:23)*2))
res.intron = apply(intron,1, f24_R2_alt2, t=c(0:23)*2)
res.exon = t(res.exon)
res.intron = t(res.intron)
res.exon = cbind(res.exon, qv=qvals(res.exon[,6]))
res.intron = cbind(res.intron, qv=qvals(res.intron[,6]))

hist(res.exon[which(res.exon[,7]<0.05),5])
hist(res.intron[which(res.intron[,7]<0.05),5])

hist(res.exon[which(res.exon[,7]<0.05 & res.exon[,3]>=0.05),5])

mrna = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/mRNA_RNA_Seq_all.txt')

kk = which(mrna$Q.RNA.Seq<0.05)

hist(mrna[kk,]$phase, breaks=24)
plot(mrna[kk,]$phase, mrna[kk,]$mean)

############################################
############### PLOT INDIVIDUAL GENE SUMMARY
####### # load('Rdata/gene_RAnalysis.Rdata')



if(plot){
	if(load.for.plot){load('Rdata/genes_RAnalysis.Rdata')}
	RefSeq = read.table('../complementary_data/RefSeq_Mouse_mm9.txt', stringsAsFactors = FALSE, sep = '\t', header = TRUE)
	source('functions.R') ;
	pdf.folder = 'plots/individual_genes_full_plots'
	circ_genes = c("Per1","Per2","Per3","Cry1","Cry2","Dbp","Tef","Hlf","Nr1d1","Nr1d2","Rora","Rorb","Rorc","Arntl","Bmal1","Clock","Npas2","Bhlhe40","Bhlhe41","Cirbp","Tfrc","Hamp","Hamp2","Nr4a2","Wee1", "Por")
	for(gene in circ_genes){cat(gene, '\n');plot.this.gene.full(pdf.folder = pdf.folder,gene = gene, ylim.rel.type = 'adjusted', ylim.rel = c(0,2), log2.rel = FALSE)}
}


###########################################################
############### PLOT INDIVIDUAL GENE PROBESETS ANIMATED GIF

if(plot){
	if(load.for.plot){load('Rdata/gene_RAnalysis.Rdata')}
	RefSeq = read.table('../complementary_data/RefSeq_Mouse_mm9.txt', stringsAsFactors = FALSE, sep = '\t', header = TRUE)

	source('functions.R') ;
	gif.folder = 'plots/gif'
	dir.create(gif.folder)
	make.gif(gif.folder = gif.folder, gene = 'Nr1d1')
	make.gif(gif.folder = gif.folder, gene = 'Rorc')
	make.gif(gif.folder = gif.folder, gene = 'Dbp')
	make.gif(gif.folder = gif.folder, gene = 'Dmd')
}



###########################################################
############### PLOT INDIVIDUAL GENE - COMPARISON WITH GUILLAUME qPCR


if(plot){
	if(load.for.plot){load('Rdata/gene_RAnalysis.Rdata')}
	source('functions.R') ;

	qPCR = read.table('../complementary_data/qPCR_Guillaume/17-1811_cdna.txt',sep='\t', col.names = c('probes','R1','R2','R3'),stringsAsFactors=F)
	#probes = c("Cry2 mRNA","Cry2 pre","Per1 mRNA","Per1 pre","Per2 mRNA","Per2 pre","Nr1d1 mRNA","Nr1d1 pre","Rorc mRNA","Rorc pre","Cry1 mRNA","Cry1 pre","Dbp mRNA","Dbp I2","Gys2 mRNA","Gys2 pre","March8 mRNA","March8 pre","Qdpr mRNA","Qdpr pre","Nfil3 mRNA","Nfil3 pre","Bhlhb2 mRNA","Bhlhb2 pre","Bhlhb3 mRNA","Bhlhb3 pre","Nr1d2 mRNA","Nr1d2 pre","Hlf mRNA","Hlf pre","Tef mRNA","Tef pre")
	zt.qPCR = seq(2,22,by=4) ; zt = seq(0,46,by = 2)
	d.zt.qPCR = c(zt.qPCR, zt.qPCR +24)
	ZT.ex = grep('.rel.ampl.ex', colnames(T)); ZT.int = grep('.rel.ampl.int', colnames(T))
	probes = unique(qPCR$probes)

	pdf('plots/comparison_with_qPCR_Guillaume.pdf', width = 10, height = 5)
	REL.AMPL.INT = matrix(NA, length(probes), 2); colnames(REL.AMPL.INT) = c('ExonArrays','qPCR'); REL.AMPL.EX = REL.AMPL.INT; genes = rep('', length(probes))
	par(mfrow = c(1,2))
	for(i in 1:length(probes)){
		probe = probes[i]
		qpcr = t(qPCR[qPCR$probes == probe,c('R1','R2','R3')])
		qpcr = qpcr/apply(qpcr,1,mean)
		gene = unlist(strsplit(probe, split = ' ')); type = gene[2]; gene = gene[1]
		j = grep(gene, T$gene); genes[i] = gene
		if(length(j)==0){cat('gene not found!    ', gene)}
		if(length(j)>2){cat('many lines for the same gene!    ', gene); j = j[1]}
		
		if(type == 'mRNA'){ZT = ZT.ex; col.qpcr = rgb(col2rgb('steelblue')[1]/255,col2rgb('steelblue')[2]/255,col2rgb('steelblue')[3]/255,0.3); col.exonarray = 'steelblue'}else{ZT = ZT.int; col.qpcr = rgb(col2rgb('green3')[1]/255,col2rgb('green3')[2]/255,col2rgb('green3')[3]/255,0.3); col.exonarray = 'green3'}
		
		
		matplot(d.zt.qPCR, t(cbind(qpcr,qpcr)) , type = 'l', lty = 1, col = col.qpcr, lwd = 5, axes = FALSE, main = probe, xlab = 'time (h)', ylab = 'relative signal', xlim = range(zt), ylim = range(qpcr,T[j,ZT], na.rm = TRUE))
		axis(1, at = d.zt.qPCR); axis(2,las = 1); box()
		points(zt,T[j,ZT],type = 'b', col = col.exonarray, lwd = 2, pch = 16)
		
		
		rel.ampl = 0
		for(k in 1:3){fft = fft(qpcr[k,]); rel.ampl = rel.ampl + 2*abs(fft[2])/6/3 }
		if(ZT == ZT.ex){
			if(length(j)>0){REL.AMPL.EX[i,1] = T$rel.ampl.ex[j]}else{REL.AMPL.EX[i,1] = NA}
			REL.AMPL.EX[i,2] = rel.ampl
		}else{
			if(length(j)>0){REL.AMPL.INT[i,1] = T$rel.ampl.int[j]}else{REL.AMPL.EX[i,1] = NA}
			REL.AMPL.INT[i,2] = rel.ampl};
		}
		
	plot(REL.AMPL.EX[,1],REL.AMPL.EX[,2],pch = 16, col = 'steelblue', cex = 1, xlim = range(REL.AMPL.EX, na.rm = TRUE)+c(0,0.1), ylim = range(REL.AMPL.EX, na.rm = TRUE)+c(0,0.1), xlab = 'ExonArrays', ylab = 'qPCR', main = 'Exons - Relative Amplitude'); abline(a = 0, b = 1)
	text(REL.AMPL.EX[,1],REL.AMPL.EX[,2],genes, pos = 3 , col = 'steelblue', cex = 0.5)
	plot(REL.AMPL.INT[,1], REL.AMPL.INT[,2],pch = 16, col = 'green3', cex = 1, xlim = range(REL.AMPL.INT, na.rm = TRUE)+c(0,0.1), ylim = range(REL.AMPL.INT, na.rm = TRUE)+c(0,0.1), xlab = 'ExonArrays', ylab = 'qPCR', main = 'Introns - Relative Amplitude'); abline(a = 0, b = 1)
	text(REL.AMPL.INT[,1], REL.AMPL.INT[,2],genes, pos = 3 , col = 'green3', cex = 0.5)
	dev.off()
	
	
	
}

### CT Values


dup = read.table('../complementary_data/qPCR_Guillaume/All_CT_values_duplicates.txt', header = FALSE, sep = '\t', fill = TRUE, stringsAsFactors = FALSE)
dup = dup[!is.na(dup[,4]),]
trip = read.table('../complementary_data/qPCR_Guillaume/All_CT_values_triplicates.txt', header = FALSE, sep = '\t', fill = TRUE, stringsAsFactors = FALSE)
trip = trip[!is.na(trip[,4]),]

qPCR = rbind(dup, trip)

qPCR$diff = qPCR[,8]-qPCR[,4] 


palette(rainbow(n = max(rank(qPCR[,2])), s = 0.8, v = 0.85))
plot(qPCR$diff, pch = 16, col = rank(qPCR[,2]))
index = seq(1,nrow(qPCR)-6, by = 6)+2
text(c(1:nrow(qPCR))[index], qPCR$diff[index],qPCR[index,2], pos = 3,col = rank(qPCR[,2])[index] )


par(mfrow = c(2,1))

plot(2^qPCR$diff, pch = 1, col = rank(qPCR[,2]), log = 'y',cex = 1.8)
text(2^qPCR$diff, substr(qPCR[,1],3,5), col = rank(qPCR[,2]), cex = 0.6)
index = seq(1,nrow(qPCR)-6, by = 6)+1
text(c(1:nrow(qPCR))[index], 2^qPCR$diff[index],qPCR[index,2], pos = 3,col = rank(qPCR[,2])[index], offset = 3)

boxplot(2^qPCR$diff ~ qPCR[,2], col = rank(qPCR[,2])[index], log = 'y')


##############
##############
##############
# RHYTHMIC GENES SELECTION

FDR = 0.1; ampl.lim = 0.15; fact = 1

if(do.step){
	cat('RHYTHMIC GENES SELECTION \n')
	if(load){load('Rdata/gene_RAnalysis.Rdata')}
	
	N = nrow(T)
	source('functions.R')
	T$rhythmic.ex.24 = select.rhythmic.genes(probe.type = 'exon', period = 24, fdr = FDR, fact = fact)
	T$rhythmic.int.24 = select.rhythmic.genes(probe.type = 'intron', period = 24, fdr = FDR, fact = fact)
	T$rhythmic.ex.12 = select.rhythmic.genes(probe.type = 'exon', period = 12, fdr = FDR, fact = fact)
	T$rhythmic.int.12 = select.rhythmic.genes(probe.type = 'intron', period = 12, fdr = FDR, fact = fact)

	save(t,T, file = 'Rdata/genes_selection_rhythmic.Rdata')
}



if(plot){
	if(load.for.plot){load('Rdata/gene_selection_rhythmic.Rdata')}
	plot.pval.histogram = function(ymax = 1, xmax = 1, period = 24, index.ex = rep(TRUE, N), index.int = rep(TRUE, N), FDR = 0.1){
		N = nrow(T)
		plot(1,1,type = 'n', xlim =c(0,xmax), ylim =c(0,ymax), xlab = 'genes sorted by pval', ylab = 'sorted pval', main = paste('period :', period,'h'))
		abline(a = 0, b = 1, col = 'gray', lwd = 0.5)
		if(period == 24){pval.ex = T$pval.ex;pval.int = T$pval.int}else{pval.ex = T$pval.12.ex;pval.int = T$pval.12.int}
		o.ex = order(pval.ex); points(c(1:N)/N, pval.ex[o.ex], pch = 16, cex = 0.3, col = c(rgb(0.9,0.9,0.9,0.2),'steelblue')[index.ex[o.ex]+1])
		o.int = order(pval.int); points(c(1:N)/N, pval.int[o.int], pch = 16, cex = 0.3, col = c(rgb(0.9,0.9,0.9,0.2),'green3')[index.int[o.int]+1])
		abline(a = 0, b = FDR)
		text(xmax,xmax * FDR, paste('FDR :', FDR), adj = c(1,-1))
		legend('topleft', legend = c('exon', 'intron'), col= c('steelblue', 'green3'), bty = 'n', lwd =3)
	}
		
	pdf('plots/pval_histogram.pdf', width = 15, height = 10)
	par(mfrow = c(2,3))
	plot.pval.histogram(ymax = 1, index.ex = T$rhythmic.ex.24,index.int = T$rhythmic.int.24, FDR = FDR)
	plot.pval.histogram(ymax = 0.2, xmax = 0.25,index.ex = T$rhythmic.ex.24,index.int = T$rhythmic.int.24, FDR = FDR)
	plot.pval.histogram(ymax = 0.03, xmax = 0.06,index.ex = T$rhythmic.ex.24,index.int = T$rhythmic.int.24, FDR = FDR)
	
	plot.pval.histogram(ymax = 1, index.ex = T$rhythmic.ex.12,index.int = T$rhythmic.int.12, FDR = FDR, period = 12)
	plot.pval.histogram(ymax = 0.2, xmax = 0.25,index.ex = T$rhythmic.ex.12, index.int = T$rhythmic.int.12, FDR = FDR, period = 12)
	plot.pval.histogram(ymax = 0.03, xmax = 0.06,index.ex = T$rhythmic.ex.12, index.int = T$rhythmic.int.12, FDR = FDR, period = 12)
	dev.off()
}

if(plot){
	if(load.for.plot){load('Rdata/gene_selection_rhythmic.Rdata')}
	N = nrow(T);ZT.ex = grep('.rel.ampl.ex', colnames(T));ZT.int = grep('.rel.ampl.int', colnames(T))
	n.comp = 8; PVAL.ex = matrix(1, nrow = N, ncol = n.comp); colnames(PVAL.ex) = 48/c(1: n.comp); PVAL.int = PVAL.ex
	for(i in 1:n.comp){cat(i,'\n'); PVAL.ex[,i] = apply(T[,ZT.ex],1,pval.n, n = i); PVAL.int[,i] = apply(T[,ZT.int],1,pval.n, n = i);}
	
	pdf('plots/pval_histogram_other_periods.pdf', width = 15, height = 10)
	par(mfcol = c(2,3))
	for(i in 2:(n.comp-1)){
		plot(c(1:N)/N, sort(PVAL.ex[,i]), type = 'l', col = 'steelblue', lwd = 2, main = colnames(PVAL.ex)[i], xlab = 'genes sorted by pval', ylab = 'sorted pval', las = 1)
		points(c(1:N)/N, sort(PVAL.int[,i]), type = 'l', col = 'green3', lwd = 2); abline(a = 0, b = 1, col = 'gray', lwd = 0.5)
		legend('topleft', legend = c('exon', 'intron'), col= c('steelblue', 'green3'), bty = 'n', lwd =3)
	}
	dev.off()
}



#### PHASE HISTOGRAM


if(plot)
{
	if(load.for.plot){load('Rdata/gene_selection_rhythmic.Rdata')}
	
	plot.phase.histogram = function(probe.type = 'exon', period = 24){
		if(probe.type == 'exon'){if(period == 24){phase = T$phase.ex; j = T$rhythmic.ex.24}else{phase = T$phase.12.ex; j = T$rhythmic.ex.12}; col = 'steelblue3'}else{if(period == 24){phase = T$phase.int; j = T$rhythmic.int.24}else{phase = T$phase.12.int; j = T$rhythmic.int.12}; col = 'green3'}

		if(period == 24){breaks = c(0:24); at.axis = seq(0,24,by = 4)}else{breaks = c(0:12); at.axis = seq(0,12,by = 3)}
		h = hist(phase[j], breaks = breaks, border = 'white', col = col, xlab = 'time (hours)', ylab = 'number of genes', main = paste('Histogram of phases -', probe.type,'- period :',period,'h'), axes = FALSE)
		axis(1, at = at.axis)
		axis(2, las = 1)
	}
	
	
	pdf('plots/phases_histograms.pdf', width = 5, height = 7)
	par(mfrow = c(2,1))
	plot.phase.histogram(probe.type = 'intron', period = 24)
	plot.phase.histogram(probe.type = 'exon', period = 24)
	dev.off()
}


#numbers of rhythmic genes in Exon and Intron
if(plot){
	if(load.for.plot){load('Rdata/gene_selection_rhythmic.Rdata')}

	pdf('plots/stat_rhythmic_genes_24h.pdf', width = 12, height = 5)
		circ_genes = c("Per1","Per2","Per3","Cry1","Cry2","Dbp","Tef","Hlf","Nr1d1","Nr1d2","Rora","Rorb","Rorc","Arntl","Bmal1","Clock","Npas2","Bhlhe40","Bhlhe41","Cirbp","Tfrc","Hamp","Hamp2","Nr4a2","Wee1", "Por")
		xmax = 12; ymax = 5
		par(mar = c(0,0,0,0)+0.4)
		plot(1,1,type = 'n',xlim = c(0,xmax), ylim = c(0, ymax), xlab = '', ylab = '', axes = FALSE)
		y = c(0,0,ymax,ymax); lwd = 5
		x.int = c(0,2*xmax/3,2*xmax/3,0); polygon(x.int, y, col = rgb(col2rgb('green3')[1]/255,col2rgb('green3')[2]/255,col2rgb('green3')[3]/255,0.5), border= 'white', lwd = lwd)
		x.ex = c(xmax/3,xmax, xmax,xmax/3); polygon(x.ex, y, col = rgb(col2rgb('steelblue')[1]/255,col2rgb('steelblue')[2]/255,col2rgb('steelblue')[3]/255,0.5), border= 'white', lwd = lwd)
		x.both = c(xmax/3,2*xmax/3,2*xmax/3,xmax/3); polygon(x.both, y, col = 'transparent', border= 'white', lwd = lwd)
		text(c(xmax/6,5*xmax/6), rep(ymax,2), c('Intron','Exon'),pos = 1, offset = 1, col = 'white', font = 2, cex = 2)
		text(c(xmax/6,5*xmax/6), rep(ymax,2), c(sum(T$rhythmic.int.24),sum(T$rhythmic.ex.24)),pos = 1, offset = 3, col = 'white', font = 2, cex = 1.5)
		text(c(xmax/6,xmax/2,5*xmax/6), rep(0,2), c(sum(T$rhythmic.int.24&!T$rhythmic.ex.24),sum(T$rhythmic.int.24&T$rhythmic.ex.24),sum(!T$rhythmic.int.24&T$rhythmic.ex.24)),pos = 3, offset = c(1,1,1), col = 'white', font = 2, cex = 1.5)
		gene.list = sort(intersect(circ_genes,T$gene[T$rhythmic.int.24&!T$rhythmic.ex.24])); text(c(xmax/6),(c(length(gene.list):1)-1)+ymax/3, gene.list, col = 'white', pos = 1)
		gene.list = sort(intersect(circ_genes,T$gene[T$rhythmic.int.24&T$rhythmic.ex.24])); text(c(xmax/2),(c(length(gene.list):1)-1)*ymax/20+ymax/3, gene.list, col = 'white', pos = 1)
		gene.list = sort(intersect(circ_genes,T$gene[!T$rhythmic.int.24&T$rhythmic.ex.24])); text(c(5*xmax/6),(c(length(gene.list):1)-1)*ymax/20+ymax/3, gene.list, col = 'white', pos = 1)
	dev.off()
}





############# Comparison of phases with Gachon data

if(plot){
	if(load.for.plot){load('Rdata/gene_selection_rhythmic.Rdata')}
	load('../complementary_data/Gachon/gachon.Rdata') ### gachon = ....

	gachon.RA = do.rhythmic.analysis(T = gachon, ZT = 1:24)
	gachon$pval_24 = gachon.RA$pval	
	gachon$phase24 = gachon.RA$phase
	gachon$rel_ampl = gachon.RA$rel.ampl

	pdf('plots/Gachon_phases_histograms.pdf', width = 5, height = 3.5)
	pval.lim = fdr2pval(P = gachon$pval_24, fdr = FDR)
	j = (gachon$pval_24 <= pval.lim)&(gachon$rel_ampl >=0.1)
	hist(gachon$phase24[j], breaks = c(0:24), main = 'Gachon', axes = FALSE, xlab = 'time (hours)', ylab = 'number of genes', col= 'slategray3', border = 'white')
	axis(1, at = seq(0,24,by = 4)); axis(2)
	dev.off()

	pdf('plots/Gachon_phases_comparison.pdf', width = 5, height = 5)
	par(pty = 's')
	j = (gachon$pval_24 <= pval.lim)&(gachon$rel_ampl >=0.1)
	k = match(gachon$gene[j], T$gene)
	plot(gachon$phase24[j], T$phase.ex[k], pch = 16, cex = 0.5, axes = FALSE, col = 'gray', xlab = 'phase - Gachon', ylab = 'phase - Exon Array', main  = 'selection : rhythmic in Gachon', xlim = c(0,24), ylim = c(0,24))
	axis(1, at = seq(0,24,by = 4));axis(2, at = seq(0,24,by = 4)); box()
	i = T$rhythmic.ex.24[k]
	points(gachon$phase24[j][i], T$phase.ex[k][i], pch = 16, cex = 0.5, col = 'orangered'); abline(a = 0, b = 1, col = 'black')
	text(25,4,'rhythmic only in Gachon', col = 'gray', pos = 2)
	text(25,2.5,'rhythmic in both', col = 'orangered', pos = 2)
	dev.off()
	
	
	
	pdf('plots/Gachon_phases_comparison_more_details.pdf', width = 10, height = 10)
	par(mar = c(2,2,0,0))
	m = match(T$gene, gachon$gene)
	G = gachon[m,]

zones = matrix(c(0,0,0,8,0,7,5,3,1,9,10,11,12,2,0,10,11,12,4,0,10,11,12,6,0), nrow = 5, ncol = 5)
layout(zones, widths = c(0.5,4,1,1,1),heights = c(1,1,1,4,0.5))


xylim = c(0,24)
plot(T$phase.ex, G$phase24, cex = 0.5, pch = 16, col = rgb(0.2,0.5,0.7,0.1), xlim = xylim, ylim = xylim, xlab = 'ExonArray',ylab = 'Gachon', axes = FALSE)
axis(1); axis(2)
pval.lim = fdr2pval(P = gachon$pval_24, fdr = 0.2)
j = (G$pval_24 <= pval.lim)&(G$rel_ampl >=0.1)
k = T$rhythmic.ex.24
sum(j,na.rm = TRUE)
points(T$phase.ex[j], G$phase24[j], cex = 0.5, pch = 16, col = rgb(0.1,0.4,0.6,1))
points(T$phase.ex[k], G$phase24[k], cex = 0.5, pch = 16, col = rgb(0.9,0.4,0.4,1))
points(T$phase.ex[j&k], G$phase24[j&k], cex = 0.7, pch = 18, col = rgb(0,0,0,1))
abline(a = 0, b = 1)

breaks = c(0:24)
hG = hist(G$phase24, breaks = breaks, plot = FALSE)
hT = hist(T$phase.ex, breaks = breaks, plot = FALSE)
hG.RG = hist(G$phase24[j], breaks = breaks, plot = FALSE)
hT.RG = hist(T$phase.ex[j], breaks = breaks, plot = FALSE)
hG.RT = hist(G$phase24[k], breaks = breaks, plot = FALSE)
hT.RT = hist(T$phase.ex[k], breaks = breaks, plot = FALSE)
hG.Rb = hist(G$phase24[j&k], breaks = breaks, plot = FALSE)
hT.Rb = hist(T$phase.ex[j&k], breaks = breaks, plot = FALSE)

barplot(hG$counts, horiz = TRUE, space = 0, col = rgb(0.2,0.5,0.7,0.1), border = NA, axes = FALSE)
barplot(hT$counts, horiz = FALSE, space = 0, col = rgb(0.2,0.5,0.7,0.1), border = NA, axes = FALSE)

barplot(hG.RG$counts, horiz = TRUE, space = 0, col = rgb(0.1,0.4,0.6,0.3), border = NA, axes = FALSE)
barplot(hG.RT$counts, horiz = TRUE, space = 0, col = rgb(0.9,0.4,0.4,0.3), border = NA, axes = FALSE, add = TRUE)

barplot(hT.RG$counts, horiz = FALSE, space = 0, col = rgb(0.1,0.4,0.6,0.3), border = NA, axes = FALSE)
barplot(hT.RT$counts, horiz = FALSE, space = 0, col = rgb(0.9,0.4,0.4,0.3), border = NA, axes = FALSE, add = TRUE)

barplot(hG.Rb$counts, horiz = TRUE, space = 0, col = 'black', border = NA, axes = FALSE)
barplot(hT.Rb$counts, horiz = FALSE, space = 0, col = 'black', border = NA, axes = FALSE)

plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1), axes = FALSE)
text(0,0,paste('Gachon'), cex=1.2, srt=90)

plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1), axes = FALSE)
text(0,0,paste('ExonArray'), cex=1.2, srt=0)


barplot(hT.RG$counts, horiz = FALSE, space = 0, col = rgb(0.1,0.4,0.6,0.3), border = 'tomato4', axes = FALSE)
barplot(hG.RG$counts, horiz = FALSE, space = 0, col = rgb(0.1,0.4,0.6,0.3), border = NA, axes = FALSE, add = TRUE)

barplot(hT.RT$counts, horiz = FALSE, space = 0, col = rgb(0.9,0.4,0.4,0.3), border = 'tomato4', axes = FALSE)
barplot(hG.RT$counts, horiz = FALSE, space = 0, col = rgb(0.9,0.4,0.4,0.3), border = NA, axes = FALSE, add = TRUE)

barplot(hG.Rb$counts, horiz = FALSE, space = 0, col = rgb(0,0,0,0.3), border = 'tomato4', axes = FALSE)
barplot(hT.Rb$counts, horiz = FALSE, space = 0, col = rgb(0,0,0,0.3), border = NA, axes = FALSE, add = TRUE)

	dev.off()
	
	
}





############## comparison with Cyclix Data


if(do.step){
	cat('COMPARISON WITH CYCLIX DATA\n')
	if(load){load('Rdata/gene_selection_rhythmic.Rdata')}
	
	
	zt.cyclix = c('ZT02','ZT06','ZT10','ZT14','ZT18','ZT22')
	
	cyclix.prom = read.table('../complementary_data/Cyclix Hi-Seq/tss_density_MQ.txt', stringsAsFactors = FALSE, header = FALSE, sep ='\t')
	
	
	
	cyclix.prom = matrix(unlist(apply(cyclix.prom, 1, strsplit, split ='   ')),nrow = nrow(cyclix.prom), ncol = 8)
	j = which(!is.na(cyclix.prom[1,])); cyclix.prom = cyclix.prom[,j]; 
	cyclix.prom[,1] = (cyclix.prom[,1]+cyclix.prom[,7])/2 ; cyclix.prom = cyclix.prom[,-7]
	colnames(cyclix.prom) = paste(zt.cyclix,'.prom', sep = ''); ZT.prom = c(1:6)
	
	
	cyclix.body = read.table('../complementary_data/Cyclix Hi-Seq/body_qnormden.txt', stringsAsFactors = FALSE, header = FALSE, sep ='\t')
	j = which(!is.na(cyclix.body[1,])); cyclix.body = cyclix.body[,j]; 
	cyclix.body[,1] = (cyclix.body[,1]+ cyclix.body[,7])/2 ; cyclix.body = cyclix.body[,-7]
	colnames(cyclix.body) = paste(zt.cyclix,'.body', sep = ''); ZT.body = c(7:12)
	
	cyclix.transcripts = read.table('../complementary_data/Cyclix Hi-Seq/trnscid.txt', stringsAsFactors = FALSE, header = FALSE, sep ='\t')
	
	cyclix = data.frame(cyclix.prom, cyclix.body, row.names = cyclix.transcripts[,1], stringsAsFactors = FALSE)
	cyclix$mean.prom = apply(cyclix[,ZT.prom], 1, mean)
	cyclix$mean.body = apply(cyclix[, ZT.body], 1, mean)
	

	
	cyclix.prom.norm = cyclix[,ZT.prom]/cyclix$mean.prom; colnames(cyclix.prom.norm) = paste(colnames(cyclix.prom.norm), '.rel.ampl',sep = '')
	cyclix.body.norm = cyclix[,ZT.body]/cyclix$mean.body; colnames(cyclix.body.norm) = paste(colnames(cyclix.body.norm), '.rel.ampl',sep = '')
	
	cyclix  = data.frame(cyclix, cyclix.prom.norm, cyclix.body.norm, stringsAsFactors = FALSE)
	ZT.prom = grep('prom.rel.ampl', colnames(cyclix)); ZT.body = grep('body.rel.ampl', colnames(cyclix))
	
	cyclix$ampl.prom = (apply(cyclix[,ZT.prom], 1, max) - apply(cyclix[,ZT.prom], 1, min))/2
	cyclix$ampl.body = (apply(cyclix[,ZT.body], 1, max) - apply(cyclix[,ZT.body], 1, min))/2
	
	
	
	source('functions.R')
	
	cyclix$phase.prom = (apply(cyclix[,ZT.prom], 1, phase.n, n = 1)/2/pi*24 + 2)%%24
	cyclix$phase.body = (apply(cyclix[,ZT.body], 1, phase.n, n = 1)/2/pi*24 + 2)%%24
	cyclix$pval.prom = apply(cyclix[,ZT.prom], 1, pval.n, n = 1)
	cyclix$pval.body = apply(cyclix[, ZT.body], 1, pval.n, n = 1)
	cyclix$rel.ampl.prom = 2*abs(t(apply(cyclix[,ZT.prom], 1, fft)))[,2]/length(ZT.prom)
	cyclix$rel.ampl.body = 2*abs(t(apply(cyclix[, ZT.body], 1, fft)))[,2]/length(ZT.body)
	
	sum(cyclix$pval.prom< 0.05, na.rm = TRUE)
	sum(cyclix$pval.body< 0.05, na.rm = TRUE)
	sum((cyclix$pval.prom< 0.05)&(cyclix$pval.body< 0.05), na.rm = TRUE)
	
	load('../complementary_data/Cyclix Hi-Seq/FilterGeneTranscripts.Rdata')
	m = match(rownames(cyclix), filter$transcript)
	cyclix$gene = filter$name[m]
	
	
	
#	library(biomaRt)
#	ensembl = useMart("ensembl", dataset = 'mmusculus_gene_ensembl')
#	annot = getBM(attributes = c('external_gene_id','ensembl_transcript_id'), filters = 'ensembl_transcript_id', values = rownames(cyclix), mart = ensembl) # other possible attributes : 'external_transcript_id','mgi_symbol','wikigene_name',
#	ii = match(rownames(cyclix), annot[,2])
#	cyclix.gene = annot[ii,1]
		
	k = match(T$gene, cyclix$gene)
	C = cyclix[k,]
	
	save(cyclix,C, file = '../complementary_data/Cyclix Hi-Seq/Cyclix_hi_seq.Rdata')
}	
	
if(plot){
	cat('plot : COMPARISON WITH CYCLIX DATA\n')
	if(load.for.plot){load('Rdata/gene_selection_rhythmic.Rdata'); load('../complementary_data/Cyclix Hi-Seq/Cyclix_hi_seq.Rdata')}
	comparison.with.cyclix = function(exon.array = c('intron','exon'), Cyclix = c('body', 'prom'), with.names = FALSE){
		exon.array = exon.array[1]; Cyclix = Cyclix[1];
		if(exon.array == 'intron'){x = log2(T$intron.median);  rel.ampl = T$rel.ampl.int; phase = T$phase.int; exp.level = T$intron.median; pval = T$pval.int}else{x = log2(T$exon.median); rel.ampl = T$rel.ampl.ex; phase = T$phase.ex; exp.level = T$exon.median; pval = T$pval.ex}
		if(Cyclix == 'body'){y = log2(C$mean.body); rel.ampl.c = C$rel.ampl.body; phase.c = C$phase.body; exp.level.c = C$mean.body; pval.c = C$pval.body}else{y = log2(C$mean.prom); rel.ampl.c = C$rel.ampl.prom; phase.c = C$phase.prom; exp.level.c = C$mean.prom; pval.c = C$pval.prom}
		
		
		rbc = rainbow(n = 100,s =0.8, v = 0.85, start = 0.6, alpha = 0.4 )
		col = rbc[round(x/max(x)*99)+1]
	
		# absolute levels
		ok = !is.na(x)&!is.na(y)&(abs(x) != Inf)&(abs(y) != Inf); x = x[ok]; y = y[ok]; cor = cor(x, y)
		
		plot(x,y, pch = 16, col = col[ok], xlab = paste('MicroArray -',exon.array), ylab = paste('Cyclix - Pol II',Cyclix,'loading'),main = paste('median ABSOLUTE level over time\n corr =',round(cor, digits = 2)))
		abline(lm(y~ x))
		if(with.names){text(x,y,T$gene[ok], col = col[ok],pos = 3, cex = 0.7)}
		
		#relative amplitudes
		ok = (pval<0.2)&(pval.c<0.2)

		sup.lim = seq(0,4,by = 0.01)
		lim = 0.12^2/sup.lim
	
		plot(rel.ampl, rel.ampl.c, pch = 16, cex = 1, ylim = c(0,1.2), col = rgb(0,0,0,0.1), xlab = paste('MicroArray -',exon.array), ylab = paste('Cyclix - Pol II',Cyclix,'loading'), main = 'Relative Amplitudes')
		points(rel.ampl[ok], rel.ampl.c[ok], pch = 16, col = col[ok])
		points(sup.lim, lim, type = 'l', col = 'white', lwd = 2)
		abline(a = 0, b = 1 , col = 'black', lwd = 1)
		if(with.names){text(rel.ampl, rel.ampl.c,T$gene, col = col,pos = 3, cex = 0.7)}

		#### phase


		ok = (pval<0.2)&(pval.c<0.2)&(rel.ampl>0.1)&(rel.ampl.c>0.1);
		phase.diff = (phase.c - phase)%%24;  phase.diff[!is.na(phase.diff)&(phase.diff>12)] = phase.diff[!is.na(phase.diff)&(phase.diff>12)]-24
		
		plot(phase, phase.c, pch = 16, cex = 0.5, col=  rgb(0,0,0,0.1), xlab = exon.array, ylab = Cyclix, main = 'Phases')
		points(phase[ok], phase.c[ok], pch = 16, cex = 1, col= col[ok])
		abline(a = 0, b = 1 , col = 'black', lwd = 1)
		if(with.names){text(phase[ok], phase.c[ok],T$gene[ok], col = col[ok],pos = 3, cex = 0.7)}


	}
	
	
	
	pdf('plots/cyclix_hi_seq.pdf', width = 6.6, height = 3.3)
	par(mar = c(3, 3, 2, 0.2) + 0.1, mgp = c(1.6, 0.5, 0), cex = 0.7, tcl = -0.3, pty = 's', mfrow = c(1,2))

	y = cyclix$mean.body; x = cyclix$mean.prom; ok = (x>0) & (y>0);  x = x[ok]; y = y[ok]
	lim = range(cyclix$mean.body, cyclix$mean.prom); lim = range(x,y)
	color.function = colorRampPalette(c('gray80', 'black'))
	cols = densCols(log(cbind(x,y)), colramp = color.function, nbin = 500)
	plot(x,y, pch = 16, cex = 0.5,col =cols, xlim = lim,ylim = lim, log ='xy' , xlab = 'prom', ylab = 'body', main = 'cyclix - F3A with Hi-seq data', axes = FALSE )
	abline(a = 0, b = 1, col = 'black', lty = 2, lwd = 2); box()
	
	y = cyclix$ampl.body; x = cyclix$ampl.prom; ok= (x<2)&(y<2); x = x[ok]; y = y[ok]
	cols = densCols(log(cbind(x,y)), colramp = color.function, nbin = 500)
	plot(x,y, pch = 16, cex = 0.5, col = cols, xlim = range(x,y, na.rm = TRUE), ylim = range(x,y, na.rm = TRUE), xlab ='amplitude - promoter', ylab = 'amplitude - body', log = 'xy')
	abline(a = 0, b = 1, col = 'black', lty = 2, lwd = 2)

	dev.off()
	
	
	
	pdf('plots/comparison_with_cyclix.pdf', width = 15, height = 12)
	par(mfcol = c(3,4), pty = 's')
	comparison.with.cyclix()
	comparison.with.cyclix(exon.array = 'intron', Cyclix = 'prom')
	comparison.with.cyclix(exon.array = 'exon', Cyclix = 'body')
	comparison.with.cyclix(exon.array = 'exon', Cyclix = 'prom')

	comparison.with.cyclix(with.names=TRUE)
	comparison.with.cyclix(exon.array = 'intron', Cyclix = 'prom',with.names=TRUE)
	comparison.with.cyclix(exon.array = 'exon', Cyclix = 'body',with.names=TRUE)
	comparison.with.cyclix(exon.array = 'exon', Cyclix = 'prom',with.names=TRUE)


	dev.off()
	
	
	#### Show some examples : circadian genes + good and bad examples
	
	
	circ_genes = c("Per1","Per2","Per3","Cry1","Cry2","Dbp","Tef","Hlf","Nr1d1","Nr1d2_chr14","Rora","Rorb","Rorc","Arntl","Clock","Npas2","Bhlhe40","Bhlhe41","Cirbp","Tfrc","Nr4a2","Wee1", "Por", 'Alb');
	good_examples = c("Gys2", "Ces3","Prei4","Igfbp2","Ly6e","Alas1")
	bad_examples = c("Fkrp","Car3","Pitx3","Hdc","Clec2h","Pnpla3","Thrsp","Ifi44","Cpxm1" ,"Gm71","Ndufs7")
	
	gene.list = c(circ_genes, good_examples, bad_examples)
	zt.int = grep('.rel.ampl.int', colnames(T))
	zt.body = grep('.body.rel.ampl', colnames(C))
	zt.cyclix = seq(2,22,by = 4); zt.EA = seq(0,22, by = 2)
	
	pdf('plots/comparison_with_cyclix_individual_genes.pdf', width = 7, height = 5)

	for(gene in gene.list){
		cat(gene,'\n')
		j = which(T$gene == gene)
		plot(zt.cyclix,C[j,zt.body], col = 'green4', lwd = 2, type = 'b', pch = 1, xlab = 'time', ylab = 'relative amplitude', main = gene, ylim = c(0,2.5), xlim = c(0,24))
		points(zt.EA, T[j, zt.int[1:12]], col= 'green3', lwd = 2, type = 'b', pch =  16)
		points(zt.EA, T[j, zt.int[13:24]], col= 'green3', lwd = 2, type = 'b', pch =  18)
		legend('topright', legend = c('Cyclix Hi-Seq\nPol II loading in gene body', 'ExonArray Day1','ExonArray Day2'), pch = c(1,16,18), col = c('green4','green3','green3'),lwd = 2)
		abline(h = 1, col = 'gray')
		}
	dev.off()
	
}


# Comparison with Cyclix : mRNA

if(plot){
	if(load.for.plot){load('Rdata/gene_selection_rhythmic.Rdata'); load('../complementary_data/Cyclix Hi-Seq/Cyclix_hi_seq.Rdata'); load('../complementary_data/Cyclix/Cyclix_microarray.Rdata'); load('../complementary_data/Cyclix/CyclixFilterGeneTranscripts.Rdata')}
	source('functions.R')
	
	Cm = CycliX.mRNA
	m = match(Cm$Probe.Set.ID, filter$probeset)
	Cm$gene = filter$name[m]
	m = which(is.na(Cm$gene))
	Cm = Cm[-m,]
	Cm[,2] = (Cm[,2]+Cm[,8])/2
	Cm = Cm[,-8] 
	
	cm = 2^Cm[,2:7]; colnames(cm) = paste('ZT',c('02','06','10','14','18','22'), sep = '')
	mean = apply(cm,1,mean)
	cm.rel = cm/mean;  colnames(cm.rel) = paste(colnames(cm),'.rel', sep = '')
	cm.rel.RA = do.rhythmic.analysis(T = cm.rel, ZT = 1:6, delta.t = 4)
	
	Cm = data.frame(Cm, cm, cm.rel, cm.rel.RA, stringsAsFactors = FALSE)
	colnames(Cm)[2:7] =  paste(colnames(cm),'.log2', sep = '')
	
	m = match(T$gene, Cm$gene)
	CM = Cm[m,]
	
	save(CM, file = '../complementary_data/Cyclix/MicroArray_processed.Rdata')
	
	pdf('plots/comparison_phase_cyclix_mRNA.pdf', width =10, height = 10)
	par(mar = c(2,2,0,0))
	zones = matrix(c(0,0,0,8,0,7,5,3,1,9,10,11,12,2,0,10,11,12,4,0,10,11,12,6,0), nrow = 5, ncol = 5)
	layout(zones, widths = c(0.5,4,1,1,1),heights = c(1,1,1,4,0.5))

	xylim = c(0,24)
	plot(T$phase.ex, CM$phase, cex = 0.5, pch = 16, col = rgb(0,0,0,0.2),xlim = xylim, ylim = xylim, axes = FALSE)
	axis(1); axis(2)
	abline(a = 0, b= 1)
	pval.lim = 0.1
	j = (CM$pval <= pval.lim)&(CM$rel.ampl >=0.1)
	k = T$rhythmic.ex.24
	points(T$phase.ex[j], CM$phase[j], cex = 0.5, pch = 16, col = rgb(0,0,1))
	points(T$phase.ex[k], CM$phase[k], cex = 0.5, pch = 16, col = rgb(1,0,0))
	points(T$phase.ex[j&k], CM$phase[j&k], cex = 0.9, pch = 18, col = rgb(0,0,0))


	breaks = c(0:24)
	hG = hist(CM$phase, breaks = breaks, plot = FALSE)
	hT = hist(T$phase.ex, breaks = breaks, plot = FALSE)
	hG.RG = hist(CM$phase[j], breaks = breaks, plot = FALSE)
	hT.RG = hist(T$phase.ex[j], breaks = breaks, plot = FALSE)
	hG.RT = hist(CM$phase[k], breaks = breaks, plot = FALSE)
	hT.RT = hist(T$phase.ex[k], breaks = breaks, plot = FALSE)
	hG.Rb = hist(CM$phase[j&k], breaks = breaks, plot = FALSE)
	hT.Rb = hist(T$phase.ex[j&k], breaks = breaks, plot = FALSE)

	barplot(hG$counts, horiz = TRUE, space = 0, col = rgb(0.2,0.5,0.7,0.1), border = NA, axes = FALSE)
	barplot(hT$counts, horiz = FALSE, space = 0, col = rgb(0.2,0.5,0.7,0.1), border = NA, axes = FALSE)

	barplot(hG.RG$counts, horiz = TRUE, space = 0, col = rgb(0.1,0.4,0.6,0.3), border = NA, axes = FALSE)
	barplot(hG.RT$counts, horiz = TRUE, space = 0, col = rgb(0.9,0.4,0.4,0.3), border = NA, axes = FALSE, add = TRUE)

	barplot(hT.RG$counts, horiz = FALSE, space = 0, col = rgb(0.1,0.4,0.6,0.3), border = NA, axes = FALSE)
	barplot(hT.RT$counts, horiz = FALSE, space = 0, col = rgb(0.9,0.4,0.4,0.3), border = NA, axes = FALSE, add = TRUE)

	barplot(hG.Rb$counts, horiz = TRUE, space = 0, col = 'black', border = NA, axes = FALSE)
	barplot(hT.Rb$counts, horiz = FALSE, space = 0, col = 'black', border = NA, axes = FALSE)

	plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1), axes = FALSE)
	text(0,0,paste('Cyclix'), cex=1.2, srt=90)

	plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1), axes = FALSE)
	text(0,0,paste('ExonArray'), cex=1.2, srt=0)


	barplot(hT.RG$counts, horiz = FALSE, space = 0, col = rgb(0.1,0.4,0.6,0.3), border = 'tomato4', axes = FALSE)
	barplot(hG.RG$counts, horiz = FALSE, space = 0, col = rgb(0.1,0.4,0.6,0.3), border = NA, axes = FALSE, add = TRUE)
	
	barplot(hT.RT$counts, horiz = FALSE, space = 0, col = rgb(0.9,0.4,0.4,0.3), border = 'tomato4', axes = FALSE)
	barplot(hG.RT$counts, horiz = FALSE, space = 0, col = rgb(0.9,0.4,0.4,0.3), border = NA, axes = FALSE, add = TRUE)

	barplot(hG.Rb$counts, horiz = FALSE, space = 0, col = rgb(0,0,0,0.3), border = 'tomato4', axes = FALSE)
	barplot(hT.Rb$counts, horiz = FALSE, space = 0, col = rgb(0,0,0,0.3), border = NA, axes = FALSE, add = TRUE)

	dev.off()
	
	
	
	
	
	}



############## comparison with RNA-Seq Data

if(do.step){
	cat('COMPARISON WITH RNA-seq DATA\n')
	if(load){load('Rdata/gene_selection_rhythmic.Rdata')}
	source('functions.R')
	
	RNA.ex = read.table(file = '../complementary_data/RNA_seq/Exons_all.txt', sep = '\t', col.names = c('chr','start','end','ID', paste('ZT0',seq(0,8,by=2),'.ex',sep = ''), paste('ZT',seq(10,22,by=2),'.ex',sep = '')))
	RNA.int = read.table(file = '../complementary_data/RNA_seq/Introns_all.txt', sep = '\t', col.names = c('chr','start','end','ID', paste('ZT0',seq(0,8,by=2),'.int',sep = ''), paste('ZT',seq(10,22,by=2),'.int',sep = '')))
	
	length.ex = RNA.ex$end - RNA.ex$start
	RNA.ex.dens = RNA.ex[5:16]/length.ex; colnames(RNA.ex.dens) = paste(colnames(RNA.ex.dens),'dens', sep = '.')
	length.int = RNA.int$end - RNA.int$start
	RNA.int.dens = RNA.int[5:16]/length.int; colnames(RNA.int.dens) = paste(colnames(RNA.int.dens),'dens', sep = '.')
	
	gene.names.ex = matrix(unlist(strsplit(as.character(RNA.ex$ID),'\\|')),3,nrow(RNA.ex))[3,]
	gene.names.int = matrix(unlist(strsplit(as.character(RNA.int$ID),'\\|')),3,nrow(RNA.int))[3,]

	gene.list = intersect(intersect(unique(gene.names.ex), unique(gene.names.int)),T$gene)
	
	RNA.seq = data.frame(gene = gene.list, matrix(NA,nrow = length(gene.list), ncol = 24)); colnames(RNA.seq)[2:13] = colnames(RNA.ex)[5:16]; colnames(RNA.seq)[14:25] = colnames(RNA.int)[5:16]; 
	
	

	
	
	for(i in 1:length(gene.list)){
		gene = as.character(RNA.seq$gene[i]);
		j = (gene == gene.names.ex)
		RNA.seq[i,2:13] = apply(RNA.ex.dens[j,],2,median)
		j = (gene == gene.names.int)
		RNA.seq[i,14:25] = apply(RNA.int.dens[j,],2,median)
		}
		
		
		
	#relative amplitude	
	
	
	RNA.seq$mean.ex = apply(RNA.seq[,2:13],1, mean)
	RNA.seq$mean.int = apply(RNA.seq[,14:25],1, mean)
	RNA.seq.ex.rel = RNA.seq[,2:13]/RNA.seq$mean.ex; colnames(RNA.seq.ex.rel) =  paste(colnames(RNA.seq)[2:13],'rel.ampl',sep = '.'); 
	RNA.seq.int.rel = RNA.seq[,14:25]/RNA.seq$mean.int; colnames(RNA.seq.int.rel) =  paste(colnames(RNA.seq)[14:25],'rel.ampl',sep = '.')
	RNA.seq = cbind(RNA.seq, RNA.seq.ex.rel, RNA.seq.int.rel)
	
	
	# rhythmic analysis
	
	ZT.RNA.ex = grep('.ex.rel.ampl', colnames(RNA.seq))
	ZT.RNA.int = grep('.int.rel.ampl', colnames(RNA.seq))
	
	RNA.seq$phase.ex = (apply(RNA.seq[,ZT.RNA.ex], 1, phase.n, n = 1)/2/pi*24)%%24
	RNA.seq$phase.int = (apply(RNA.seq[,ZT.RNA.int], 1, phase.n, n = 1)/2/pi*24)%%24
	RNA.seq$pval.ex = apply(RNA.seq[, ZT.RNA.ex], 1, pval.n, n = 1)
	RNA.seq$pval.int = apply(RNA.seq[, ZT.RNA.int], 1, pval.n, n = 1)
	RNA.seq$rel.ampl.ex = 2*abs(t(apply(RNA.seq[, ZT.RNA.ex], 1, fft)))[,2]/length(ZT.RNA.ex)
	RNA.seq$rel.ampl.int = 2*abs(t(apply(RNA.seq[, ZT.RNA.int], 1, fft)))[,2]/length(ZT.RNA.int)	
				
	m = match(T$gene, RNA.seq$gene)	
	R = RNA.seq[m,]
	

	save(RNA.seq, R , file = '../complementary_data/RNA_seq/RNA_seq_average_density_per_gene.Rdata');
	}


	
if(plot){
	cat('plot : COMPARISON WITH RNA-Seq DATA\n')
	if(load.for.plot){load('Rdata/gene_selection_rhythmic.Rdata'); load('../complementary_data/RNA_seq/RNA_seq_average_density_per_gene.Rdata')}
	
	
	
	pdf('plots/RNA_seq.pdf', width = 6.6, height = 3.3)
	par(mar = c(3, 3, 2, 0.2) + 0.1, mgp = c(1.6, 0.5, 0), cex = 0.7, tcl = -0.3, pty = 'm', mfrow = c(1,2), las = 1)

	x = log10(R$mean.int); y = log10(R$mean.ex); 
	color.function = colorRampPalette(c('gray90','turquoise4'))
	cols = densCols(cbind(x,y), colramp = color.function, nbin = 500)
	plot(x,y, pch = 16, cex = 0.5,col =cols, xlab = 'log10(median density in introns)', ylab = 'log10(median density in exons)', main = 'RNA-seq' )
	abline(a = 1, b = 1, col = 'black', lty = 3, lwd = 1)
	abline(a = 1.5, b = 1, col = 'black', lty = 1, lwd = 1)
	abline(a = 2, b = 1, col = 'black', lty = 2, lwd = 1)
	legend('topleft',legend = c('10x','~30x','100x'), lty = c(3,1,2) , bty = 'n')


	ZT.RNA.int = grep('.int.rel.ampl',colnames(R)); ZT.RNA.ex = grep('.ex.rel.ampl',colnames(R))

	x = apply(R[, ZT.RNA.int], 1, max) / apply(R[, ZT.RNA.int], 1, min) 
	y = apply(R[, ZT.RNA.ex], 1, max) / apply(R[, ZT.RNA.ex], 1, min)
	breaks = seq(-7,7, by = 0.2)
	hist(log2(x/y), breaks = breaks, axes = FALSE, col = c(rep('gray80',floor(length(breaks)/2)),rep('gray50',floor(length(breaks)/2))), border = NA, xlab = 'log2(FC int / FC ex)', main = '')
	axis(1);box()


	dev.off()
	
	
	comparison.with.RNA.seq = function(exon.array = c('intron','exon'), RNA = c('intron', 'exon'), with.names = FALSE){
		exon.array = exon.array[1]; RNA = RNA[1];
		if(exon.array == 'intron'){x = log2(T$intron.median);  rel.ampl = T$rel.ampl.int; phase = T$phase.int; exp.level = T$intron.median; pval = T$pval.int}else{x = log2(T$exon.median); rel.ampl = T$rel.ampl.ex; phase = T$phase.ex; exp.level = T$exon.median; pval = T$pval.ex}
		if(RNA == 'intron'){y = log10(R$mean.int); rel.ampl.c = R$rel.ampl.int; phase.c = R$phase.int; exp.level.c = R$mean.int; pval.c = R$pval.int}else{y = log10(R$mean.ex); rel.ampl.c = R$rel.ampl.ex; phase.c = R$phase.ex; exp.level.c = R$mean.ex; pval.c = R$pval.ex}
		
		
		rbc = rainbow(n = 100,s =0.8, v = 0.85, start = 0.6, alpha = 0.4 )
		col = rbc[round(x/max(x)*99)+1]
	
		# absolute levels
		ok = !is.na(x)&!is.na(y)&(abs(x) != Inf)&(abs(y) != Inf); x = x[ok]; y = y[ok]; cor = cor(x, y)
		
		plot(x,y, pch = 16, col = col[ok], xlab = paste('MicroArray -',exon.array), ylab = paste('RNA-Seq - ',RNA),main = paste('median ABSOLUTE level over time\n corr =',round(cor, digits = 2)))
		abline(lm(y~ x))
		if(with.names){text(x,y,T$gene[ok], col = col[ok],pos = 3, cex = 0.7)}
		
		#relative amplitudes

		sup.lim = seq(0,4,by = 0.01)
		lim = 0.12^2/sup.lim
	
		ok = (pval<0.2)&(pval.c<0.2)
	
		plot(rel.ampl, rel.ampl.c, pch = 16, cex = 1, ylim = c(0,1.2), col = rgb(0,0,0,0.1), xlab = paste('MicroArray -',exon.array), ylab = paste('RNA-Seq - ',RNA), main = 'Relative Amplitudes')
		points(rel.ampl[ok], rel.ampl.c[ok], pch = 16, col = col[ok])
		points(sup.lim, lim, type = 'l', col = 'white', lwd = 2)
		abline(a = 0, b = 1 , col = 'black', lwd = 1)
		if(with.names){text(rel.ampl[ok], rel.ampl.c[ok],T$gene[ok], col = col[ok],pos = 3, cex = 0.7)}

		#### phase


		ok = (pval<0.2)&(pval.c<0.2)&(rel.ampl>0.1)&(rel.ampl.c>0.1);
		phase.diff = (phase.c - phase)%%24;  phase.diff[!is.na(phase.diff)&(phase.diff>12)] = phase.diff[!is.na(phase.diff)&(phase.diff>12)]-24
		
		plot(phase, phase.c, pch = 16, cex = 0.5, col=  rgb(0,0,0,0.1), xlab = paste('MicroArray -',exon.array), ylab = paste('RNA-Seq - ',RNA), main = 'Phases')
		points(phase[ok], phase.c[ok], pch = 16, cex = 1, col= col[ok])
		abline(a = 0, b = 1 , col = 'black', lwd = 1)
		if(with.names){text(phase[ok], phase.c[ok],T$gene[ok], col = col[ok],pos = 3, cex = 0.7)}


	}
	
	
	

	
	
	pdf('plots/comparison_with_RNA_seq.pdf', width = 15, height = 12)
	par(mfcol = c(3,4), pty = 's')
	comparison.with.RNA.seq(exon.array = 'intron', RNA = 'intron')
	comparison.with.RNA.seq(exon.array = 'exon', RNA = 'exon')

	comparison.with.RNA.seq(exon.array = 'intron', RNA = 'intron', with.names=TRUE)
	comparison.with.RNA.seq(exon.array = 'exon', RNA = 'exon', with.names=TRUE)

	dev.off()
	
	
	#### Show some examples : circadian genes + good and bad examples
	
	
	circ_genes = c("Per1","Per2","Per3","Cry1","Cry2","Dbp","Tef","Hlf","Nr1d1","Nr1d2_chr14","Rora","Rorb","Rorc","Arntl","Clock","Npas2","Bhlhe40","Bhlhe41","Cirbp","Tfrc","Nr4a2","Wee1", "Por", 'Alb');
	good_examples = c("Gys2", "Ces3","Prei4","Igfbp2","Ly6e","Alas1")
	bad_examples = c("Fkrp","Car3","Pitx3","Hdc","Clec2h","Pnpla3","Thrsp","Ifi44","Cpxm1" ,"Gm71","Ndufs7")
	
	gene.list = c(circ_genes, good_examples, bad_examples)
	zt.int = grep('.rel.ampl.int', colnames(T)); zt.ex = grep('.rel.ampl.ex', colnames(T)); 
	ZT.RNA.int = grep('.int.rel.ampl',colnames(R)); ZT.RNA.ex = grep('.ex.rel.ampl',colnames(R))
	zt = seq(0,22, by = 2)
	
	pdf('plots/comparison_with_RNA_seq_individual_genes.pdf', width = 7, height = 10)
	par(mfrow = c(2,1))

	for(gene in gene.list){
		cat(gene,'\n')
		j = which(T$gene == gene)
		plot(zt,R[j, ZT.RNA.int], col = 'green4', lwd = 2, type = 'b', pch = 1, xlab = 'time', ylab = 'relative amplitude', main = paste(gene, '- Introns'), ylim = c(0,2.5), xlim = c(0,24))
		points(zt, T[j, zt.int[1:12]], col= 'green3', lwd = 1, type = 'b', pch =  16, lty = 2)
		points(zt, T[j, zt.int[13:24]], col= 'green3', lwd = 1, type = 'b', pch =  18, lty = 2)
		legend('topright', legend = c('RNA-Seq', 'Array Day1','Array Day2'), pch = c(1,16,18), col = c('green4','green3','green3'),lwd = c(2,1,1), lty = c(1,2,2))
		#legend('topleft', legend = round(c(R$rel.ampl.int[j], T$rel.ampl.int[j]), digits = 2),lty = 0, bty = 'n')
		abline(h = 1, col = 'gray')
		
		plot(zt,R[j, ZT.RNA.ex], col = 'steelblue4', lwd = 2, type = 'b', pch = 1, xlab = 'time', ylab = 'relative amplitude', main = paste(gene, '- Exons'), ylim = c(0,2.5), xlim = c(0,24))
		points(zt, T[j, zt.ex[1:12]], col= 'steelblue', lwd = 1, type = 'b', pch =  16, lty = 2)
		points(zt, T[j, zt.ex[13:24]], col= 'steelblue', lwd = 1, type = 'b', pch =  18, lty = 2)
		legend('topright', legend = c('RNA-Seq', 'Array Day1','Array Day2'), pch = c(1,16,18), col = c('steelblue4','steelblue','steelblue'),lwd = c(2,1,1), lty = c(1,2,2))
		#legend('topleft', legend = round(c(R$rel.ampl.ex[j], T$rel.ampl.ex[j]), digits = 2),lty = 0, bty = 'n')
		abline(h = 1, col = 'gray')
		}
	dev.off()
	
}

###### CHECK FOR OTHER PATTERNS

if(plot){
	cat('plot : CHECK FOR OTHER PATTERNS\n')	
	if(load.for.plot){load('Rdata/gene_selection_rhythmic.Rdata')}
	source('functions.R')
	pdf('plots/other_patterns.pdf', width = 14, height = 16)
	check.for.other.pattern(probe.type = 'exon')
	check.for.other.pattern(probe.type = 'intron')
	dev.off()

}


####### Plot the Ratio vs. Delay scatter

if(plot){
	cat('plot : Ratio vs. Delay scatter \n')	
	if(load.for.plot){load('Rdata/gene_selection_rhythmic.Rdata')}
	source('functions.R'); N = nrow(T)
	
		
	r.or = T$rhythmic.ex.24|T$rhythmic.int.24
	r.and = T$rhythmic.ex.24&T$rhythmic.int.24

	Ratio = T$rel.ampl.ex/T$rel.ampl.int
	Delay = (T$phase.ex-T$phase.int)%%24; Delay[Delay>18] = Delay[Delay>18]-24
	
	pdf('plots/ratio_phase.pdf', height = 9, width = 6)
	
	layout(matrix(c(1,2),2,1), height = c(1,1.5))
	par(mar = c(4,3,1,1)+0.4, mgp = c(2,0.6,0))
	breaks = c(-6:18)
	hist(Delay[r.or], breaks = breaks, axes = FALSE, col = c(rep('gray85',6), rep('gray75',6), rep('gray85',12)), border = NA, xlab = 'Delay [hours]', ylab = 'number of genes', main = '')
	hist(Delay[r.and], breaks = breaks, axes = FALSE, col = c(rep('slategray2',6), rep('slategray3',6), rep('slategray2',12)), border = NA, add = TRUE)
	axis(1,at = seq(-6,18,by = 3)); axis(2)
	legend('topright', legend = c('Genes rhythmic in exon OR intron','idem - outside predicted phase','Genes rhythmic in exon AND intron','idem - outside predicted phase'), fill = c('gray75','gray85','slategray3','slategray2'), bty = 'n', border = NA, text.col = c('black','gray75','black','gray75'), cex = 0.7)
	
	
	
	colors =rep('transparent', N); colors[r.or] = 'gray85'; colors[r.and] ='steelblue4'; #colors = rep(rgb(0,0,0,0.3), N)
	ylim = range(Ratio[r.or])
	plot(1,1,type = 'n',xlim = c(-6,18), ylim = ylim, ylab = 'Ratio [rel. ampl. exon / rel. ampl. intron]', xlab ='Delay [phase exon - phase intron]', log = 'y', axes = FALSE)
	base = seq(0,6,by = 0.1); cos = cos(base/24 * 2 * pi)
	upper.lim = cos + 0.5; lower.lim = cos - 0.5
	polygon(c(base[upper.lim>0],6,6,0),c(upper.lim[upper.lim>0],ylim[1]/100,1000,1000),col = 'mistyrose', border = 'NA')
	polygon(c(base[lower.lim>0],6,0),c(lower.lim[lower.lim>0],ylim[1]/100,ylim[1]/100),col = 'antiquewhite1', border = 'NA')
	B =(length(base)-1);for(i in 1:B){polygon(c(base[i],base[i+1],base[i+1], base[i]),c(max(ylim[1]/100,lower.lim[i]),max(ylim[1]/100,lower.lim[i+1]), upper.lim[i+1], upper.lim[i]), border = 'NA',col = 'azure1')}# col = rgb(0.8*B+0.2*i,i,i,maxColorValue = B))}
	points( Delay[r.or], Ratio[r.or], col = colors[r.or], pch = 16, cex = 0.5 )
	points(Delay[r.and], Ratio[r.and], col = colors[r.and],pch = 16, cex = 0.5)
	points(base,cos, type = 'l', lwd = 2)
	legend('topright', legend = c('Genes rhythmic in exon OR intron','Genes rhythmic in exon AND intron'), col = c('gray85','steelblue4'), pch = 16, cex = 0.7, bg = 'white')
	axis(1,at = seq(-6,18,by = 3)); axis(2, at = c(ylim[1],0.1,0.3,1,3,10, ylim[2]), labels = c('',0.1,0.3,1,3,10,''))

	dev.off()
	
#	plot(1,1,type = 'n',xlim = c(0,6), ylim = c(0.001,2), ylab = 'Ratio [rel. ampl. exon / rel. ampl. intron]', xlab ='Delay [phase exon - phase intron]')
#	bx.and = boxplot(Ratio[r.and]~round(Delay[r.and]), plot = FALSE)
#	bx.or = boxplot(Ratio[r.or]~round(Delay[r.or]), plot = FALSE)
#	polygon(c(bx.or$names, rev(bx.or$names)), c(bx.or$stats[2,],rev(bx.or$stats[4,])), col = 'gray85', border = NA)
#	polygon(c(bx.and$names, rev(bx.and$names)), c(bx.and$stats[2,],rev(bx.and$stats[4,])), col = 'slategray1', border = NA)
#	points(bx.or$names, bx.or$stats[3,], col = 'gray75', lwd = 3, type = 'l')
#	points(bx.and$names, bx.and$stats[3,], col = 'slategray3', lwd = 3, type = 'l')
#	points(seq(0,6,by = 0.01),cos(seq(0,6,by = 0.01)/24 * 2 * pi) , type = 'l', lwd = 2)
#
#	
#
#	colors =rep('transparent', N); colors[r.or] = 'gray85'; colors[r.and] ='slategray4'; #colors = rep(rgb(0,0,0,0.3), N)
#	plot(1,1,type = 'n',xlim = c(-6,18), ylim = c(0.001,2), ylab = 'Ratio [rel. ampl. exon / rel. ampl. intron]', xlab ='Delay [phase exon - phase intron]')
#	points( Delay[r.or], Ratio[r.or], col = colors[r.or], pch = 16, cex = 0.5 )
#	points(Delay[r.and], Ratio[r.and], col = colors[r.and],pch = 16, cex = 0.5)
#	points(seq(0,6,by = 0.01),cos(seq(0,6,by = 0.01)/24 * 2 * pi) , type = 'l', lwd = 2)
#	legend('topright', legend = c('Genes rhythmic in exon OR intron','Genes rhythmic in exon AND intron'), col = c('gray85','slategray4'), pch = 16)	

	}
	
########################################################################################################################################################
########################################################################################################################################################
############################################################################################### finish preliminary analysis of exon and intron profiles and different comparision of databases
########################################################################################################################################################
########################################################################################################################################################


########################################################### 
###########################################################

# FITS  : Validation of the method with simulated data . Noise is log-normal distributed!

################################################simulated data first
####################################

########## ESTIMATE the standard deviation of noise in introns and exons for real data

load = TRUE
my.do.step = TRUE

if(my.do.step)
{
	cat('ESTIMATION OF STANDARD DEVIATION OF NOISE\n')
	if(load){load('my_genes_RAnalysis_larger_exons_2.Rdata')}
	#load('my_genes_RAalysis_all.Rdata')
	source('functions.R')
	source("f24_modified_1.0.r")
	T = Tt
	
	exon = as.matrix(T[,c(2:25)])
	intron = as.matrix(T[,c(26:49)])
	res.ex = (apply(exon,1, f24_R2_alt2, t=c(0:23)*2))
	res.int = apply(intron,1, f24_R2_alt2, t=c(0:23)*2)
	res.ex = t(res.ex)
	res.int = t(res.int)
	res.ex = cbind(res.ex, qv=qvals(res.ex[,6]))
	res.int = cbind(res.int, qv=qvals(res.int[,6]))
	
	## Different criterion for the selection of rhythmic genes
	length(which(res.ex[,7]<0.05))
	length(which(res.int[,7]<0.05))
	
	ii = which(res.ex[,7]<0.10 & res.ex[,4]>=0.1);
	hist(res.ex[ii,5],breaks=c(0:24))
	ii = which(res.int[,7]<0.10 & res.int[,4]>=0.1);
	hist(res.int[ii,5],breaks=c(0:24))
	
	## replace the results calculated by FISH_24 function by f24_R2_alt2 funciton
	Tt$rel.ampl.ex = res.ex[,4]
	Tt$rel.ampl.int = res.int[,4]
	Tt$phase.ex = res.ex[,5]
	Tt$phase.int = res.int[,5]
	Tt$pval.ex = res.ex[,6]
	Tt$pval.int = res.int[,6]
	Tt$Qval.ex = res.ex[,7]
	Tt$Qval.int = res.int[,7]
	
	#save(Tt, tt, file="my_genes_RAalysis_all_final.Rdata")
	save(Tt, tt, file="my_genes_RAnalysis_larger_exons_final.Rdata");
	
	load(file="my_genes_RAnalysis_larger_exons_final.Rdata");
	
	index.outliers = function(data.xx)
	{
		#data.xx = c(2, 3, 6, 9, 13, 18, 21, 106)
		Q1 = quantile(data.xx, 0.25,type=5)
		Q3 = quantile(data.xx, 0.75, type=5)
		IQD = Q3 - Q1
		lower = Q1 - 1.5*IQD
		upper = Q3 + 1.5*IQD
		index = which(data.xx<lower|data.xx>upper)
		#boxplot(data.xx);abline(h=Q1);abline(h=Q3);
	}
	
	ZT.ex = c(2:25)
	ZT.int = c(26:49)
	
	cutoff = c(4:8)/100
	sigma.s = c()
	sigma.m = c()
	sigma.ms = c()
	
	pdf("myplots/Estimation_Sigma_real_data.pdf",width = 15, height = 6)
	par(mfcol = c(1,3), pty = 's')
	
	for(pval.cutoff in cutoff)
	{
		
		cat('Pval is ', pval.cutoff, '\n');
		data.ex = log(as.matrix(Tt[,ZT.ex]))[which(Tt$pval.ex>pval.cutoff),]
		data.int = log(as.matrix(Tt[,ZT.int]))[which(Tt$pval.int>pval.cutoff),]
		print(length(which(Tt$pval.ex>pval.cutoff)));
		print(length(which(Tt$pval.int>pval.cutoff)))
	
		data.ex = as.vector(data.ex)
		data.int = as.vector(data.int)
		data.ex.int = c(data.ex, data.int)
		
		#cutoff.quantile = 0.02
		
		#new.data.ex = data.ex[which(data.ex>quantile(data.ex,cutoff.quantile) & data.ex<quantile(data.ex, 1-cutoff.quantile))];
		#new.data.int = data.int[which(data.int>quantile(data.int,cutoff.quantile) & data.int<quantile(data.int, 1-cutoff.quantile))];
		#new.data.ex.int = data.ex.int[which(data.ex.int>quantile(data.ex.int,cutoff.quantile) & data.ex.int<quantile(data.ex.int, 1-cutoff.quantile))];
		
		new.data.ex = data.ex[-index.outliers(data.ex)]
		new.data.int = data.int[-index.outliers(data.int)]
		new.data.ex.int = data.ex.int[-index.outliers(data.ex.int)]
		
		# not remove outliers with sd
		cat('sigma estimated by sd with outliers...');
		print(c(sd(data.int), sd(data.ex), sd(data.ex.int)));
		
		cat('sigma estimated with median and outliers...');
		
		print(c(sqrt(1/(length(data.int)-1)*median(data.int^2)*length(data.int)), sqrt(1/(length(data.ex)-1)*median(data.ex^2)*length(data.ex)), sqrt(1/(length(data.ex.int)-1)*median(data.ex.int^2)*length(data.ex.int))));
		
		
		# remove outliers with sd
		cat('sigma estimated by sd without outliers...means is ');
		print(c(mean(new.data.int), mean(new.data.ex), mean(new.data.ex.int)));
		cat('sigma estimated by sd without outliers...');
		print(c(sd(new.data.int), sd(new.data.ex), sd(new.data.ex.int)));
		
		new.sigma.s = sqrt(1/(length(new.data.int)-1)*mean(new.data.int^2)*length(new.data.int))
		new.sigma.m = sqrt(1/(length(new.data.ex)-1)*mean(new.data.ex^2)*length(new.data.ex))
		new.sigma.ms = sqrt(1/(length(new.data.ex.int)-1)*mean(new.data.ex.int^2)*length(new.data.ex.int))
		
		# remove outliers with mean to estimate sigma
		cat('sigma estimated by mean without outliers...');
		print(c(new.sigma.s, new.sigma.m, new.sigma.ms));
		
		sigma.s = c(sigma.s, new.sigma.s)
		sigma.m = c(sigma.m, new.sigma.m)
		sigma.ms = c(sigma.ms, new.sigma.ms)
		
		
		plot(density(data.ex), col='black',type='p',cex=0.5,main= paste('Exon, pval>', pval.cutoff, sep=''))
		points(density(rnorm(n=length(data.ex), mean=0, sd=sd(data.ex))), col='orange',type='l')
		points(density(rnorm(n=length(data.ex), mean=0, sd=sd(new.data.ex))), col='blue',type='l')
		points(density(rnorm(n=length(data.ex), mean=0, sd=new.sigma.m)), col='darkred',type='l')
		
		plot(density(data.int), col='black',type='p',cex=0.5,main= paste('Intron, pval>', pval.cutoff, sep=''))
		points(density(rnorm(n=length(data.int), mean=0, sd=sd(data.int))), col='orange',type='l')
		points(density(rnorm(n=length(data.int), mean=0, sd=sd(new.data.int))), col='blue',type='l')
		points(density(rnorm(n=length(data.int), mean=0, sd=new.sigma.s)), col='darkred',type='l')
		
		plot(density(data.ex.int), col='black',type='p',cex=0.5,main= paste('Exon and Intron, pval>', pval.cutoff, sep=''))
		points(density(rnorm(n=length(data.ex.int), mean=0, sd=sd(data.ex.int))), col='orange',type='l')
		points(density(rnorm(n=length(data.ex.int), mean=0, sd=sd(new.data.ex.int))), col='blue',type='l')
		points(density(rnorm(n=length(data.ex.int), mean=0, sd=new.sigma.ms)), col='darkred',type='l')
		

	}
	dev.off()
	
	sigma.s.real = mean(sigma.s)
	sigma.m.real = mean(sigma.m)
	sigma.ms.real = mean(sigma.ms)
	
	save(sigma.s.real, sigma.m.real, sigma.ms.real, file="Sigma_my_genes_real_data.Rdata");
	
	pdf("myplots/Standard_Deviation_Exons.pdf",width = 6, height = 6)
	hist(data.ex, breaks=200, main='Normal distriubtion of log Exon (pval>0.2)')
	abline(v=0,col='red', lwd=2.0)
	abline(v=-sd(new.data.ex), col='red',lwd=1.0)
	abline(v=sd(new.data.ex), col='red',lwd=1.0)
	dev.off()
	
	pdf("myplots/Standard_Deviation_Introns.pdf",width = 6, height = 6)
	hist(data.int, breaks=200, main='Normal distriubtion of log Intron (pval>0.2)')
	abline(v=0,col='red', lwd=2.0)
	abline(v=-sd(new.data.int), col='red',lwd=1.0)
	abline(v=sd(new.data.int), col='red',lwd=1.0)
	dev.off()
	
	
}

#################### generate fake data using the standard deviation of noise in real data 

load = TRUE
my.do.step = TRUE
zt = seq(0,46,by = 2)
sd.noise = c(0.1,0.25,0.5)
#sd.noise = c(0.1)
X = 100
load = TRUE
intense.debug <-- FALSE
my.do.step = TRUE
#intense.debug <-- TRUE

if(my.do.step)
{
	cat('GENERATE THE FAKE DATA\n')
	#if(load){load('Rdata/genes_selection_rhythmic.Rdata')}
	load('my_genes_RAalysis_all_final.Rdata') ## all expressed genes
	## This is the my data complete
	#load('my_genes_RAnalysis.Rdata')
	T = Tt
	source('model_modify_new/functions.R')
	
	#set.global.variable()
	#xx = Tt$exon.median/Tt$intron.median/splicing.rate
	#xx = Tt$exon.median*(log(2)/5)/Tt$intron.median/(log(2)/8*60)
	#xx = log(xx)
	#hist(xx, breaks=100)
	
	
	ZT.int = grep('.rel.ampl.int', colnames(T)); 
	ZT.ex = grep('.rel.ampl.ex', colnames(T)); 
	
	fake.version = '_diff_sigma_signal_ratio_final_'

	for(sd in sd.noise)
	{
		cat('sd noise : ',sd, '\n')
		F = T[1:(4*X),]; 
		F$gamma = NA; 
		F$eps.gamma = NA ; 
		F$phase.gamma = NA;
		for(model in c(1:4))
		{
			cat('model : ',model, '\n')
			Fm = generate.fake.data(T = F[((model-1)*X+1):(model*X),],X = X, sd = sd, model = model, parametrization = 'sigmoid', sum.species = FALSE, random.splicing=FALSE) # parametrization = 'sigmoid'
			if(model == 1){FF = Fm}
			else{FF = rbind(FF,Fm)}
			#F[((model-1)*X+1):(model*X),] = Fm
		}
		F = FF;
				
		F.r.ex = do.rhythmic.analysis(T = F, ZT = ZT.ex); 
		colnames(F.r.ex) = paste(colnames(F.r.ex), '.ex',sep = '');
		F.r.int = do.rhythmic.analysis(T = F, ZT = ZT.int); 
		colnames(F.r.int) = paste(colnames(F.r.int), '.int',sep = '');

		F[,colnames(F.r.ex)] = F.r.ex; 

		F[,c('fscore.int','pval.int','pval.12.int')] = F.r.int[,c('fscore.int','pval.int','pval.12.int')]; 
		colnames(F.r.int) = paste(colnames(F.r.int), '.RA',sep = ''); 

		F = data.frame(F, F.r.int[,2:5], stringsAsFactors = FALSE);
		write.table(F, file = paste('simulated_data/my_simulated_data_diff_sigma_signal_ratio_',100*sd,'.txt',sep = ''), quote = FALSE, sep = '\t', row.names = FALSE)
		
		eval(parse(text = paste('T',100*sd,'= F',sep ='')))
	}
	
	filename = paste(file = 'my_simulated_data', fake.version, '.Rdata', sep='')
	eval(parse(text = paste("save(",paste(paste("T",100*sd.noise,sep = ""),collapse = ","),", file=filename)", sep = "")))
	#eval(parse(text = paste('save(',paste(paste('T',100*sd.noise,sep = ''),collapse = ','),', file = "my_simulated_data_diff_sigma_signal_ratio.Rdata")', sep = '')))

}

################## check standard deviation of noise in the generated fake data
if(my.do.step)
{
	cat('CHECK THE FAKE DATA\n');
	load = TRUE
	
	#fake.version = '_diff_sigma_signal_ratio_2'
	filename = paste(file = 'my_simulated_data', fake.version, '.Rdata', sep='')
	if(load){load(file=filename)}
	source('model_modify_new/functions.R')
	source("f24_modified_1.0.r")
	
	sd.noise = c(0.1,0.25,0.5)
	
	verify.fake.data = FALSE
	if(verify.fake.data)
	{
		sd = 0.1
		
		cat(sd, '\n')
		eval(parse(text = paste('T = T',100*sd,sep ='')))
		ZT.int = grep('.rel.ampl.int', colnames(T))
		ZT.ex = grep('.rel.ampl.ex', colnames(T))
		zt = seq(0,46,by = 2)
		zt.p = seq(0,48, by=0.5)
		
		j = 103
		
		source('model_modify_new/functions.R')
		#source('functions.R')
		intense.debug = FALSE
		param.fits.results = make.fits.with.all.models.for.one.gene(T = T, gene.index = j, debug = TRUE, parametrization = 'sigmoid', method = 'integration', sum.species = FALSE, zt = zt, i.ex = ZT.ex, i.int = ZT.int, absolute.signal = FALSE);
		param.fits.results
		my.BIC.this.gene(param.fits.results)
		
		pdf.name = paste('myplots/Plots_fake_data', fake.version, sd*100, '.pdf', sep='')
		pdf(pdf.name, width = 10, height = 5)
		par(mfrow = c(1,2))
		for(i in 1:nrow(T))
		{
			eval(parse(text = paste('gamma = T$gamma[i];eps.gamma = T$eps.gamma[i];phase.gamma = T$phase.gamma[i];fold.change = T$fold.change.int[i]; phase.int = T$phase.int[i]; up.time = T$up.time.int[i]; down.time = T$up.time.int[i]', sep = '')))
		
			mt = compute.m.sigmoid(t = zt.p, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)
			st = compute.sigmoid(t = zt.p, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)
			
			plot(zt, as.numeric(T[i,ZT.int]), col='blue',type='b',lty=1, ylim=range(c(as.numeric(T[i,ZT.int]), st, c(0.5,1.5))), main=paste(T$gene[i], ', Intron', sep=''))
			points(zt.p, st, col='darkred', type='l',lty=1, lwd=2.0)
			
			plot(zt, as.numeric(T[i,ZT.ex]), col='green', type='b',lty=1, ylim=range(c(as.numeric(T[i,ZT.ex]), mt, c(0.5,1.5))), main=paste(T$gene[i], ', Exon', sep=''))
			points(zt.p, mt, col='darkred', type='l',lty=1,lwd=2.0)
			
			abline(h=1, col='darkgray')
		}
		dev.off()
		
	}
	
	## check the sigma of noise in the fake data
	#T = T10
	sigma.ss = c()
	sigma.mm = c()
	sigma.ms = c()
	
	for(sd in sd.noise)
	{
		cat(sd, '\n')
		eval(parse(text = paste('T = T',100*sd,sep ='')))
		exon = as.matrix(T[,c(2:25)])
		intron = as.matrix(T[,c(26:49)])
		res.ex = (apply(exon,1, f24_R2_alt2, t=c(0:23)*2))
		res.int = apply(intron,1, f24_R2_alt2, t=c(0:23)*2)
		res.ex = t(res.ex)
		res.int = t(res.int)
		res.ex = cbind(res.ex, qv=qvals(res.ex[,6]))
		res.int = cbind(res.int, qv=qvals(res.int[,6]))
		T$exon.var = apply(log(exon), 1, sd)
		#T$exon.median = res.ex[,2]
		T$intron.var = apply(log(intron), 1, sd)
		#T$intron.median = res.int[,2]
		
		set.global.variable()
		#xx = Tt$exon.median/Tt$intron.median/splicing.rate
		xx = T$exon.median*(log(2)/5)/T$intron.median/(log(2)/8*60)
		xx = log(xx)
		hist(xx, breaks=10)
		
	## Different criterion for the selection of rhythmic genes
	length(which(res.ex[,7]<0.05))
	length(which(res.int[,7]<0.05))
	
	ii = which(res.ex[,7]<0.10 & res.ex[,4]>=0.1);
	length(ii)
	#hist(res.ex[ii,5],breaks=c(0:24))
	#ii = which(res.int[,7]<0.10 & res.int[,4]>=0.1);
	#hist(res.int[ii,5],breaks=c(0:24))
	
	## replace the results calculated by FISH_24 function by f24_R2_alt2 funciton
	T$rel.ampl.ex = res.ex[,4]
	T$rel.ampl.int = res.int[,4]
	T$phase.ex = res.ex[,5]
	T$phase.int = res.int[,5]
	T$pval.ex = res.ex[,6]
	T$pval.int = res.int[,6]
	T$Qval.ex = res.ex[,7]
	T$Qval.int = res.int[,7]
	
	ZT.ex = c(2:25)
	ZT.int = c(26:49)
	
	cutoff = c(1:4)/10
	sigma.s = c()
		sigma.m = c();
		sigma.all = c()
	for(pval.cutoff in cutoff)
	{
		#pval.cutoff = 0.5
		data.ex = log(as.matrix(T[,ZT.ex]))[which(T$pval.ex>pval.cutoff),]
		data.int = log(as.matrix(T[,ZT.int]))[which(T$pval.int>pval.cutoff),]
		print(length(which(T$pval.ex>pval.cutoff)));
		print(length(which(T$pval.int>pval.cutoff)))
		
		data.ex = as.vector(data.ex)
		data.int = as.vector(data.int)
		data.ex.int = c(data.ex, data.int)
		
		#new.data.ex = data.ex[which(data.ex>quantile(data.ex,0.001) & data.ex<quantile(data.ex, 0.999))];
		#new.data.int = data.int[which(data.int>quantile(data.int,0.001) & data.int<quantile(data.int, 0.999))];
		#new.data.ex.int = data.ex.int[which(data.ex.int>quantile(data.ex.int,0.001) & data.ex.int<quantile(data.ex.int, 0.999))];
		
		new.data.ex = data.ex;
		new.data.int = data.int;
		new.data.ex.int = data.ex.int
		#cat('sigma estimated without outliers\n');
		#print(c(sd(new.data.ex),sd(new.data.int)));
		#cat('sigma estimated with outliers\n');
		#print(c(sd(data.ex), sd(data.int)));
		sigma.s = c(sigma.s, sd(new.data.int))
		sigma.m = c(sigma.m, sd(new.data.ex))
		sigma.all = c(sigma.all, sd(new.data.ex.int))
		
	}
	
		sigma.s = mean(sigma.s)
		sigma.m = mean(sigma.m)
		sigma.all = mean(sigma.all)
		sigma.ss = c(sigma.ss, sigma.s)
		sigma.mm = c(sigma.mm, sigma.m)
		sigma.ms = c(sigma.ms, sigma.all)
		cat('sigams estimated :\n')
		print(sigma.s);print(sigma.m);print(sigma.all)
		set.global.sigma();
		if(sd==0.1) {sigma.s = sigma.s/5; sigma.m = sigma.m/5;}
		if(sd==0.25) {sigma.s = sigma.s/2; sigma.m = sigma.m/2;}
		if(sd==0.5) {sigma.s = sigma.s; sigma.m = sigma.m;}
		cat('true sigams:\n')
		print(sigma.s);print(sigma.m);
		
		eval(parse(text = paste('T',100*sd,'= T',sep ='')))
	}
	
	save(sigma.ss, sigma.mm, sigma.ms, file="Sigma_fake_data.Rdata")
}


######################## fitting different models for the fake data considering exons and introns have the same noise, thereby we do not need to specify different sigma in the optimization
if(my.do.step)
{
	
	cat('FIT THE FAKE DATA ON VITAL-IT WITH OLD METHODS (SAME SIGMA FOR EXONS AND INTRON)\n')
	load = TRUE
	if(load){load('my_simulated_data_diff_sigma_all.Rdata')}
	
	source('functions.R')
	
	sd.noise = c(0.1,0.25,0.5)
	#sd.noise = c(0.25,0.5)
	#sd = sd.noise[1]
	for(sd in sd.noise)
	{
		eval(parse(text = paste('T = T',100*sd,sep ='')))
		
		save(T, file = 'stuff_for_vital_IT/genes_ready_for_fits.Rdata')
		gene_list = T$gene
		write.table(gene_list, file = 'stuff_for_vital_IT/gene_list.txt', quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE)
		
		## ON VITAL-IT 
		
		#clear the directory
		cmd = "ssh jwang@dev.vital-it.ch 'rm -r /scratch/cluster/monthly/jingkui/microarray/*' " 
		system(cmd)

		#copy the data and the code on vital-IT
		cmd = 'scp -r stuff_for_vital_IT/* jwang@dev.vital-it.ch:/scratch/cluster/monthly/jingkui/microarray/.'
		system(cmd)
		cmd = 'scp -r functions.R jwang@dev.vital-it.ch:/scratch/cluster/monthly/jingkui/microarray/.'
		system(cmd)
	
		# run the fits : 
		cmd = "ssh jwang@dev.vital-it.ch 'bash -s' < stuff_for_vital_IT/to_run_on_vital_IT.sh"
		system(cmd)
		
		# concatenate all fits results together when fits are all done
		cmd = "ssh jwang@dev.vital-it.ch 'perl /scratch/cluster/monthly/jingkui/microarray/cat_all_fits_results.pl'"
		system(cmd)

		# wait for everything to be done and copy the data on the local computer
		cmd = "ssh jwang@dev.vital-it.ch 'ls -l /scratch/cluster/monthly/jingkui/microarray/fits/all_genes_fits_results_with_header.txt'"
		res = system(cmd, intern = TRUE)
	
		while(length(res) == 0)
		{
			Sys.sleep(60*5) # wait for 5 minutes
			cat('... wait for the fits to finish\n')
			cmd = "ssh jwang@dev.vital-it.ch 'ls -l /scratch/cluster/monthly/jingkui/microarray/fits/all_genes_fits_results_with_header.txt'"
			res = system(cmd, intern = TRUE)
		}
		
		cat('........................................ fits finished on vital-IT\n')	
	
		cmd = paste("scp jwang@dev.vital-it.ch:/scratch/cluster/monthly/jingkui/microarray/fits/all_genes_fits_results_with_header.txt simulated_data/fit_simulated_data_diff_sigma_all/Old_Methods_fits_simulated_data_",100*sd,".txt", sep = '')
		system(cmd)
	
		cat('........................................ fits results imported on local computer\n')
		
		
		fits.results = read.table(file = paste("simulated_data/fit_simulated_data_diff_sigma_all/Old_Methods_fits_simulated_data_",100*sd,".txt", sep = ''), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
		if(nrow(fits.results)<0.9*length(gene_list)){cat('........................................ fits results missing for more than 10% of genes\n')}
		m = match(gene_list, fits.results$gene)
		if(sum(is.na(m))>0){cat('........................................ fits results missing for',sum(is.na(m)),'genes\n')}
		fits.results = fits.results[m,]
	
		T = data.frame(T, fits.results[,-1],stringsAsFactors = FALSE)
		
		eval(parse(text = paste('T',100*sd,'= T',sep ='')))
	}
	eval(parse(text = paste("save(",paste(paste("T",100*sd.noise,sep = ""),collapse = ","),", file = 'my_simulated_data_with_fits_results_diff_sigma_all_Old_Method.Rdata')", sep = "")))
		
}

######################## Fitting different models for the fake data using new method in which noise in exons and introns are considered seperately in the optimization
if(my.do.step)
{
	
	cat('FIT THE FAKE DATA ON VITAL-IT WITH NEW METHODS\n')
	load = TRUE
	if(load){load('my_simulated_data_diff_sigma_all.Rdata')}
	
	source('functions.R')
	
	sd.noise = c(0.1,0.25,0.5);
	load("Sigma_fake_data.Rdata");
	#sd.noise = c(0.25,0.5)
	#sd = sd.noise[1]
	#for(sd in sd.noise)
	{	
		ii = 1
		sd = sd.noise[ii]
		sigma.ss[ii];
		sigma.mm[ii];
		
		
		eval(parse(text = paste('T = T',100*sd,sep ='')))
		
		save(T, file = 'stuff_for_vital_IT/genes_ready_for_fits.Rdata')
		gene_list = T$gene
		write.table(gene_list, file = 'stuff_for_vital_IT/gene_list.txt', quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE)
		
		## ON VITAL-IT 
		
		#clear the directory
		cmd = "ssh jwang@dev.vital-it.ch 'rm -r /scratch/cluster/monthly/jingkui/microarray/*' " 
		system(cmd)
		
		#copy the data and the code on vital-IT
		cmd = 'scp -r stuff_for_vital_IT/* jwang@dev.vital-it.ch:/scratch/cluster/monthly/jingkui/microarray/.'
		system(cmd)
		cmd = 'scp -r functions.R jwang@dev.vital-it.ch:/scratch/cluster/monthly/jingkui/microarray/.'
		system(cmd)
		
		# run the fits directly on vital-it
		## bash -s <to_run_on_vital_IT.sh
		#cmd = "ssh jwang@dev.vital-it.ch 'bash -s' < stuff_for_vital_IT/to_run_on_vital_IT.sh"
		#system(cmd)
		
		# concatenate all fits results together when fits are all done
		cmd = "ssh jwang@dev.vital-it.ch 'perl /scratch/cluster/monthly/jingkui/microarray/cat_all_fits_results.pl'"
		system(cmd)
		
		# wait for everything to be done and copy the data on the local computer
		cmd = "ssh jwang@dev.vital-it.ch 'ls -l /scratch/cluster/monthly/jingkui/microarray/fits/all_genes_fits_results_with_header.txt'"
		res = system(cmd, intern = TRUE)
		
		while(length(res) == 0)
		{
			Sys.sleep(60*5) # wait for 5 minutes
			cat('... wait for the fits to finish\n')
			cmd = "ssh jwang@dev.vital-it.ch 'ls -l /scratch/cluster/monthly/jingkui/microarray/fits/all_genes_fits_results_with_header.txt'"
			res = system(cmd, intern = TRUE)
		}
		
		cat('........................................ fits finished on vital-IT\n')	
		
		cmd = paste("scp jwang@dev.vital-it.ch:/scratch/cluster/monthly/jingkui/microarray/fits/all_genes_fits_results_with_header.txt simulated_data/fit_simulated_data_diff_sigma_all/New_Methods_fits_simulated_data_",100*sd,".txt", sep = '')
		system(cmd)
		
		cat('........................................ fits results imported on local computer\n')
		
		
		fits.results = read.table(file = paste("simulated_data/fit_simulated_data_diff_sigma_all/New_Methods_fits_simulated_data_",100*sd,".txt", sep = ''), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
		if(nrow(fits.results)<0.9*length(gene_list)){cat('........................................ fits results missing for more than 10% of genes\n')}
		m = match(gene_list, fits.results$gene)
		if(sum(is.na(m))>0){cat('........................................ fits results missing for',sum(is.na(m)),'genes\n')}
		fits.results = fits.results[m,]
		
		T = data.frame(T, fits.results[,-1],stringsAsFactors = FALSE)
		
		eval(parse(text = paste('T',100*sd,'= T',sep ='')))
		
		#save(T10,T25, file='my_simulated_data_with_fits_results_diff_sigma_all_New_Method_T10_T25_temp.Rdata')
	}
	
	
	eval(parse(text = paste("save(",paste(paste("T",100*sd.noise,sep = ""),collapse = ","),", file = 'my_simulated_data_with_fits_results_diff_sigma_all_New_Method.Rdata')", sep = "")))
	
}

####### OLD METHODS of LRT, BIC, AIC to select the model for the fake data


load = TRUE
my.do.step = TRUE
sd.noise = c(0.1,0.25,0.5)

final = TRUE
fake.version = '_diff_sigma_signal_ratio'
fake.version = '_diff_sigma_signal_ratio_final_model3_not_treated_'

if(my.do.step)
{
	cat('LRT, AIC, BIC WITH FAKE DATA\n')
	#load('my_simulated_data_with_fits_results_diff_sigma_all_Old_Method_T10.Rdata')
	
	if(final)
	{
		filename = paste(file = 'my_simulated_data_with_fits_results', fake.version, '.Rdata', sep='')
		load(file=filename)
		load("Sigma_fake_data.Rdata")
		source('model_modify_new/functions.R')
		X = 100
		pval = 0.05
		
		for(sd in sd.noise)
		{
			cat(sd,'\n')
			
			eval(parse(text = paste('T = T',100*sd,sep ='')))
			LRT = Likelihood.Ratio.Test.LAURA(T = T, pval = pval)
			LRT.FDR = Likelihood.Ratio.Test.LAURA(T = T, pval = pval, FDR=TRUE)
			colnames(LRT.FDR) = paste(colnames(LRT.FDR), '.FDR', sep='');
			
			AIC = AIC.LAURA(T = T, correction = FALSE)
			AICc = AIC.LAURA(T = T, correction = TRUE)
			BIC = BIC.LAURA(T = T)
			
			T = data.frame(T, AIC, AICc, BIC, LRT, LRT.FDR, stringsAsFactors = FALSE);
			eval(parse(text = paste('T',100*sd,'= T',sep ='')))
		}
		
		filename = paste(file = 'my_simulated_data_with_fits_results_BIC_AIC_LRT', fake.version, '.Rdata', sep='')
		eval(parse(text = paste('save(',paste(paste('T',100*sd.noise,sep = ''),collapse = ','),', file = filename)', sep = '')))
		
	}else{
	if(load){
		load('my_simulated_data_with_fits_results_diff_sigma_all_Old_Method_v2.Rdata');
		TT10 = T10;
		TT25 = T25;
		TT50 = T50;
		rm(list=c("T10", "T25", "T50"));
		load('my_simulated_data_with_fits_results_diff_sigma_all_New_Method.Rdata');
	}
	
	
	#eval(parse(text = paste("save(",paste(paste("T",100*sd.noise,sep = ""),collapse = ","),", file = 'my_simulated_data_with_fits_results_diff_sigma_all_Old_Method_v2.Rdata')", sep = "")))
	
	load("Sigma_fake_data.Rdata")
	source('functions.R')
	X = 100
	
	
		
	pdf('myplots/simulated_data_diff_sigma_all_Methods.pdf', width = 21, height = 17.5)
	par(mfrow = c(3,1))
	
	for(sd in sd.noise)
	{
		cat(sd,'\n')
		
		fact = sigma.ms[which(sd.noise==sd)]
		
		eval(parse(text = paste('T = T',100*sd,sep ='')))
		## caculate the rhythmicities with f24_R2_alt2
		source("f24_modified_1.0.r")
		exon = as.matrix(T[,c(2:25)])
		intron = as.matrix(T[,c(26:49)])
		res.ex = (apply(exon,1, f24_R2_alt2, t=c(0:23)*2))
		res.int = apply(intron,1, f24_R2_alt2, t=c(0:23)*2)
		res.ex = t(res.ex)
		res.int = t(res.int)
		res.ex = cbind(res.ex, qv=qvals(res.ex[,6]))
		res.int = cbind(res.int, qv=qvals(res.int[,6]))
			
		T$rel.ampl.ex = res.ex[,4]
		T$rel.ampl.int = res.int[,4]
		T$phase.ex = res.ex[,5]
		T$phase.int = res.int[,5]
		T$pval.ex = res.ex[,6]
		T$pval.int = res.int[,6]
		T$Qval.ex = res.ex[,7]
		T$Qval.int = res.int[,7]
						
		eval(parse(text = paste('TT = TT',100*sd,sep ='')))
		
		
		pval = 0.05
		## function of likelihood ratio test
		LRT.best.model = my.Likelihood.Ratio.Test(T = T, pval = pval)
		LRT.best.model.fdr = my.Likelihood.Ratio.Test(T = T, pval = pval, FDR=TRUE)
		colnames(LRT.best.model.fdr) = paste(colnames(LRT.best.model.fdr), '.FDR', sep='');
		
		LRT.best.model.old = my.Likelihood.Ratio.Test.Old(T = TT, pval = pval, fact=fact)
		LRT.best.model.old.fdr = my.Likelihood.Ratio.Test.Old(T = TT, pval = pval, fact=fact, FDR=TRUE)
		colnames(LRT.best.model.old.fdr) = paste(colnames(LRT.best.model.old.fdr), '.FDR', sep='');
		
		LRT.best.model.laura = Likelihood.Ratio.Test.LAURA(T = TT, pval = pval)
		LRT.best.model.laura.fdr = Likelihood.Ratio.Test.LAURA(T = TT, pval = pval, FDR=TRUE)
		colnames(LRT.best.model.laura.fdr) = paste(colnames(LRT.best.model.laura.fdr), '.FDR', sep='');
		#colnames(LRT.and.best.model.laura) = paste(colnames(LRT.and.best.model.laura), '.laura', sep='')
		#T = data.frame(T, LRT.and.best.model.old, LRT.and.best.model.laura, stringsAsFactors = FALSE)
		
		
		#AIC = my.AIC(T = T, correction = FALSE, diff.sigma=TRUE)
		AIC = my.AIC(T = T, correction = FALSE)
		AIC.old = my.AIC.Old(T=TT, correction=FALSE, fact=fact)
		AIC.laura = AIC.LAURA(T = TT,correction=FALSE)
		
		plot(AIC.laura$AIC.best.model.Laura, main = paste('AIC LAURA, Noise: ', 100*sd, sep=''))
		plot(AIC.old$AIC.best.model.Old, main = paste('AIC OLD, Noise: ', 100*sd, sep=''))
		plot(AIC$AIC.best.model, main = paste('AIC NEW, Noise: ', 100*sd, sep=''))
		
		
		#AICc = my.AIC(T = T,correction = TRUE,diff.sigma=TRUE)
		AICc = my.AIC(T = T, correction = TRUE)
		AICc.old = my.AIC.Old(T=TT, correction=TRUE, fact=fact)
		AICc.laura = AIC.LAURA(T = TT, correction=TRUE)
		
		plot(AICc.laura$AICc.best.model.Laura, main = paste('AICc LAURA, Noise:', 100*sd, sep=''))
		plot(AICc.old$AICc.best.model.Old, main = paste('AICc OLD, Noise:', 100*sd, sep=''))
		plot(AICc$AICc.best.model, main = paste('AICc NEW, Noise:', 100*sd, sep=''))
		
		
		#BIC = my.BIC(T = T,diff.sigma=TRUE)
		BIC = my.BIC(T = T)
		BIC.old = my.BIC.Old(T = TT, fact =fact)
		BIC.laura = BIC.LAURA(T=TT)
		
		plot(BIC.laura$BIC.best.model.Laura, main = paste('BIC LAURA, Noise:', 100*sd, sep=''))
		plot(BIC.old$BIC.best.model.Old, main = paste('BIC OLD, Noise:', 100*sd, sep=''))
		plot(BIC$BIC.best.model, main = paste('BIC NEW, Noise:', 100*sd, sep=''))
				
		plot(LRT.best.model.laura$LRT.best.model.Laura, main = paste('LRT LAURA, Noise:', 100*sd, sep=''))
		plot(LRT.best.model.old$LRT.best.model.Old, main = paste('LRT OLD, Noise:', 100*sd, sep=''))
		plot(LRT.best.model$LRT.best.model, main = paste('LRT NEW, Noise:', 100*sd, sep=''))
		
		colnames(TT)[c(84:102)] = paste(colnames(TT)[c(84:102)], '.old', sep='')
		
		T = data.frame(T, AIC, AICc, BIC, LRT.best.model, LRT.best.model.fdr, TT[,c(84:102)], AIC.old, AICc.old, BIC.old, LRT.best.model.old, LRT.best.model.old.fdr, AIC.laura, AICc.laura, BIC.laura, LRT.best.model.laura, LRT.best.model.laura.fdr, stringsAsFactors = FALSE)
		eval(parse(text = paste('T',100*sd,'= T',sep ='')))
		
	}
	dev.off()
	eval(parse(text = paste('save(',paste(paste('T',100*sd.noise,sep = ''),collapse = ','),', file = "my_simulated_data_with_fits_results_LRT_and_AIC_BIC_all_versions_v2.Rdata")', sep = '')))
		
	}
}


######### summary of comparison between estimation of Likelihood and real parameters and models

## quantify the difference of these model-selecting methods
my.do.step = TRUE
load = TRUE

if(my.do.step)
{
	cat('COMPARISION OF DIFFERENT MODEL SELECTING METHODS BY QUANTIFICAITON\n')
	if(load){load(file = "my_simulated_data_with_fits_results_LRT_and_AIC_BIC_all_versions_v2.Rdata")}
	source('functions.R')
	library(stringr)
	
	sd.noise = c(0.1,0.25,0.5)
	selecting.approach = c('BIC.best.model','AIC.best.model', 'AICc.best.model', 'LRT.best.model',
						   'BIC.best.model.Old','AIC.best.model.Old','AICc.best.model.Old', 'LRT.best.model.Old',
						   'BIC.best.model.Laura', 'AIC.best.model.Laura', 'AICc.best.model.Laura', 'LRT.best.model.Laura', 'LRT.best.model.Laura.FDR')
	
	X = 100
	#measures = c('M1', 'M2', 'M3', 'M4', 'M34','M234')
	measures = c('M3', 'M4', 'M34','M234')
	error.type = c('Fasle.positive','False.negative', 'FDR', 'FNR', 'F.0.5', 'MCC')
	model = c(3:4)
	mmeasures = c()
	for(test in measures) mmeasures = c(mmeasures, paste(test, error.type, sep='.'))
	measures = mmeasures
	
	#sd = sd.noise[2]
	#sd.noise = c(0.25, 0.1, 0.5)
	
	comparison.final = c()
	
	#pdf('myplots/Simulated_data_fits_results_comparision_diff_methods_quantification.pdf', width = 12, height = 4)
	#par(mfcol = c(1,3), cex =0.45, mar = c(4,3,1,1)+0.4, mgp = c(2,0.6,0),las = 1)
	#par(mar = c(4,3,1,1)+0.4, mgp = c(2,0.6,0))
	for(sd in sd.noise)
	{
		
	cat(sd,'\n')
		pdf.name = paste('myplots/Simulated_data_fits_results_comparision_diff_methods_quantification_sd_',sd, '_v2.pdf', sep='')
		pdf(pdf.name, width = 12, height = 3)
		par(mfcol = c(1,4), cex =0.45, mar = c(4,3,1,1)+0.4, mgp = c(2,0.6,0),las = 1)
	
		
		eval(parse(text = paste('T = T',100*sd,sep ='')))
	#T$LRT.best.model = T$best.model
	index = c(1:400)
	T$true.model = ceiling(index/100)
	
	## HERE Hypothesis H0: THIS GENE IS NOT REGULATED BY CIRCADIAN DEGRADATION OR CIRCADIAN TRANSCRIPTION. UNLESS YOU HAVE ENOUGH EVIDENT, WE WEILL HOLD THIS H0
	comparison = c()
	counts = 1
	for(method in selecting.approach)
	{
		eval(parse(text = paste('best.model=T$', method, sep='')))
		true.model=T$true.model
		err = c()
		
		for(m in model)
		{
			tp = length(which(true.model==m & best.model == m))
			fp = length(which(true.model!=m & best.model == m)) # False positive: reject the hypothesis which is actually true
			tn = length(which(true.model!=m & best.model != m))
			fn = length(which(true.model==m & best.model != m)) # False negative: fail to reject the hypothesis which is actually false
			fdr = fp/(fp+tp)
			fnr = fn/(fn+tp)
			acc = (tp+tn)/(tp+tn+fp+fn)
			f1 = 1.25*tp/(1.25*tp+0.25*fn+fp)
			mcc = (tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
			
			err = c(err, c(fp, fn, fdr, fnr, f1, mcc))
			#err = c(err, err.2-err.1)
			
		}
		
		tp = length(which((true.model==3|true.model==4) & (best.model == 3|best.model == 4)))
		fp = length(which((true.model==1|true.model==2) & (best.model == 3|best.model == 4)))
		tn = length(which((true.model==1|true.model==2) & (best.model == 1|best.model == 2)))
		fn = length(which((true.model==3|true.model==4) & (best.model == 1|best.model == 2)))
		fdr = fp/(fp+tp)
		fnr = fn/(fn+tp)
		acc = (tp+tn)/(tp+tn+fp+fn)
		f1 = 1.25*tp/(1.25*tp+0.25*fn+fp)
		mcc = (tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
		
		err = c(err, c(fp, fn, fdr, fnr, f1, mcc))
		
		tp = length(which(true.model!=1 & best.model != 1))
		fp = length(which(true.model==1 & best.model != 1))
		tn = length(which(true.model==1 & best.model == 1))
		fn = length(which(true.model!=1 & best.model == 1))
		fdr = fp/(fp+tp)
		fnr = fn/(fn+tp)
		acc = (tp+tn)/(tp+tn+fp+fn)
		f1 = 1.25*tp/(1.25*tp+0.25*fn+fp)
		mcc = (tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
		#print(c(tp,fp,tn,fn,mcc))
		
		err = c(err, c(fp, fn, fdr, fnr, f1, mcc))
				
		comparison = cbind(comparison, err)
		counts = counts + 1;
	}
	comparison = data.frame(comparison, stringsAsFactors = FALSE)
	colnames(comparison) = selecting.approach
	rownames(comparison) = measures
	comparison = t(comparison)
		
		## plots of FDR, F_0.5, MCC for different model-selecting methods:
		
		ii.fdr = grep('FDR', colnames(comparison))
		ii.fnr = grep('FNR', colnames(comparison))
		ii.f.0.5 = grep('F.0.5', colnames(comparison))
		ii.mcc = grep('MCC', colnames(comparison))
		
		#points(d, type = 'l', col = rainbow[(i.rel-1)%%12+1], lty = (i.rel-1)%/%12+1)
		#axis(1, at = 1:4, labels = c('CS-CD','RS-CD','CS-RD','RS-RD'), lwd = 0)
		#axis(2, at = 1:4, labels = c('CS-CD','RS-CD','CS-RD','RS-RD'), lwd = 0, las = 3)
		#legend('bottomright',legend = c('MCMC prediction','Optim prediction') , pch = c(1,16))
		#legend(x = 0.8*xlim[2],y = ylim[2]*((i.rel+1)/length(zt)), legend = paste('ZT',2*(i.rel-1)), col =  rainbow[(i.rel-1)%%12+1], lty = (i.rel-1)%/%12+1, bty = 'n' )
		
		cexx = c(2.0, 2.0, rep(1.5, 11))
		ltyy = c(1, 1, rep(2, 11))
		lwdd = c(1.5, 1.5, rep(1.2, 11))
		nb.methods = nrow(comparison)
		test = comparison[, ii.fdr]
		methods.legend = c('BIC','AIC', 'AICc', 'LRT','BIC.Old','AIC.Old','AICc.Old', 'LRT.Old','BIC.Laura', 'AIC.Laura', 'AICc.Laura', 'LRT.Laura', 'LRT.Laura.FDR')
		rainbow = rainbow(nb.methods,s = 0.85, v = 0.85)
		
		i.rel = c(1:nrow(comparison))
		ylim = range(test[which(!is.na(test)==TRUE)])
		matplot(c(1:4), t(test), type='b', lwd=lwdd, cex=cexx, pch = c(1:nrow(comparison)), col = rainbow[(i.rel-1)%%nb.methods+1], lty = ltyy, xlab='Classification of Genes', ylab='FDR',axes = FALSE)
		for(ii in i.rel)legend(x = 0.8*4,y = ylim[2]-(ylim[2]-ylim[1])*(ii-1)/nrow(comparison), legend = methods.legend[ii], col =  rainbow[(ii-1)%%nb.methods+1], pch = ii, lty = ltyy[ii], bty = 'n',cex=1.0)
		axis(1, at = c(1:4), labels = c('M3','M4','M34','M234'), lwd = 0)
		axis(2)
		box()
		
		test = comparison[, ii.fnr]
		i.rel = c(1:nrow(comparison))
		ylim = range(test[which(!is.na(test)==TRUE)])
		matplot(c(1:4), t(test), type='b', lwd=lwdd, cex=cexx, pch = c(1:nrow(comparison)), col = rainbow[(i.rel-1)%%nb.methods+1], lty = ltyy, xlab='Classification of Genes', ylab='False Negative Rate (FNR)',axes = FALSE)
		for(ii in i.rel)legend(x = 0.8*4,y = ylim[2]-(ylim[2]-ylim[1])*(ii-1)/nrow(comparison), legend = methods.legend[ii], col =  rainbow[(ii-1)%%nb.methods+1], pch = ii, lty = ltyy[ii], bty = 'n',cex=1.0)
		axis(1, at = c(1:4), labels = c('M3','M4','M34','M234'), lwd = 0)
		axis(2)
		box()
		
		
		test = comparison[, ii.f.0.5]
		ylim = range(test[which(!is.na(test)==TRUE)])
		matplot(c(1:4), t(test), type='b', lwd=lwdd, cex=cexx, pch = c(1:nrow(comparison)), col = rainbow[(i.rel-1)%%nb.methods+1], lty = ltyy, xlab='Classification of Genes', ylab='F_0.5 Score',axes = FALSE)
		for(ii in i.rel)legend(x = 0.8*4,y = ylim[2]-(ylim[2]-ylim[1])*(ii-1)/nrow(comparison), legend = methods.legend[ii], col =  rainbow[(ii-1)%%nb.methods+1], pch = ii, lty = ltyy[ii], bty = 'n',cex=1.0)
		axis(1, at = c(1:4), labels = c('M3','M4','M34','M234'), lwd = 0)
		axis(2)
		box()
		
		test = comparison[, ii.mcc]
		ylim = range(test[which(!is.na(test)==TRUE)])
		matplot(c(1:4), t(test), type='b', lwd=lwdd, cex=cexx, pch = c(1:nrow(comparison)), col = rainbow[(i.rel-1)%%nb.methods+1], lty = ltyy, xlab='Classification of Genes', ylab='MCC',axes = FALSE)
		for(ii in i.rel)legend(x = 0.8*4,y = ylim[2]-(ylim[2]-ylim[1])*(ii-1)/nrow(comparison), legend = methods.legend[ii], col =  rainbow[(ii-1)%%nb.methods+1], pch = ii, lty = ltyy[ii], bty = 'n',cex=1.0)
		axis(1, at = c(1:4), labels = c('M3','M4','M34','M234'), lwd = 0)
		axis(2)
		box()
		
		dev.off()
		
	comparison.final = rbind(comparison.final, comparison)
	}
	
	
	
	write.table(comparison.final, file='Power_comparison_model_selecting_methods.txt', sep='\t',col.names=TRUE, row.names=TRUE, quote=FALSE)
	
	
}

## plot the comparison of different model-selecting methods: BIC, AIC, AICc, LRT
load.for.plot = TRUE
plot=TRUE
palette(c('gray','green3','tomato','black'))
if(plot)
{
	cat('plot : FITS RESULTS FOR SIMULATED DATA\n')
	filename = paste(file = 'my_simulated_data_with_fits_results_BIC_AIC_LRT', fake.version, '.Rdata', sep='')
	if(load.for.plot){load(file = filename)}

	source('model_modify_new/functions.R')
	library(stringr)
	X = 100
	
	#selecting.approach = c('LRT.best.model','BIC.best.model','AICc.best.model','AIC.best.model','LRT.best.model.Old','BIC.best.model.Old','AICc.best.model.Old','AIC.best.model.Old', 'LRT.best.model.Laura', 'LRT.best.model.Laura.FDR', 'BIC.best.model.Laura','AICc.best.model.Laura','AIC.best.model.Laura')
	
	selecting.approach = c('LRT.best.model','LRT.best.model.FDR','BIC.best.model','AICc.best.model','AIC.best.model')
	
	pdf.name = paste(file = 'myplots/my_simulated_data_with_fits_results_BIC_AIC_LRT', fake.version, '.pdf', sep='')
	pdf(pdf.name, width = 10, height = 7.5)
	par(mfcol = c(4,6), cex =0.45, mar = c(3,3,2,0.2)+0.1, mgp = c(1.6,0.5,0),las = 1, tcl = -0.3)
	for(select.method in selecting.approach)
	{
	
		for(sd in sd.noise)
		{
		cat(sd,'\n')
		eval(parse(text = paste('T = T',100*sd,sep ='')))
		#T$LRT.best.model = T$best.model
		index = c(1:400)
		T$true.model = ceiling(index/100)
		
		kk = which(!is.na(T$BIC.best.model)==TRUE)
		T = T[kk,]
		
		eval(parse(text = paste('T$LRT.best.model=T$', select.method, sep='')))
		
		Prediction.Quality = matrix(NA, nrow = 4, ncol =4)
		for(model in 1:4)
		{
			#h = hist(T$LRT.best.model[((model-1)*X+1):(model*X)], breaks = seq(0.5,4.5,by = 1), plot = FALSE)
			h = hist(T$LRT.best.model[which(T$true.model==model)], breaks = seq(0.5,4.5,by = 1), plot = FALSE)
			Prediction.Quality[model,] = h$counts 
		}
		barplot(t(Prediction.Quality), col = 1:4, names.arg = c('CS-CD','RS-CD','CS-RD','RS-RD'), border = NA, main =paste(select.method, ': Noise:',sd*100,'%'))
		compare.parameters(T, quality = FALSE, parametrization = 'sigmoid')
		}
	}
	dev.off()
	
	### check the issue of fake data
	check = TRUE
	if(check)
	{
		load = TRUE
		my.do.step = TRUE
		sd.noise = c(0.1,0.25,0.5)
		
		final = TRUE
		#fake.version = '_diff_sigma_signal_ratio_2'
		filename = paste(file = 'my_simulated_data_with_fits_results', fake.version, '.Rdata', sep='')
		load(file=filename)
		load("Sigma_fake_data.Rdata")
		source('model_modify_new/functions.R')
		X = 100
		pval = 0.05
		
		for(sd in sd.noise)
		{
			cat(sd,'\n')
			
			eval(parse(text = paste('T = T',100*sd,sep ='')))
			LRT = Likelihood.Ratio.Test.LAURA(T = T, pval = pval)
			AIC = AIC.LAURA(T = T, correction = FALSE)
			AICc = AIC.LAURA(T = T, correction = TRUE)
			BIC = BIC.LAURA(T = T)
			index = c(1:400)
			T$true.model = ceiling(index/100)
			
			T = data.frame(T, AIC, AICc, BIC, LRT, stringsAsFactors = FALSE);
			eval(parse(text = paste('T',100*sd,'= T',sep ='')))
		}
		
		sd = 0.5
		eval(parse(text = paste('T = T',100*sd,sep ='')))
		index = c(1:400)
		T$true.model = ceiling(index/100)
		ii = which(T$true.model==4 & T$BIC.best.model!=4)
		kk = which(T$true.model==3 & T$BIC.best.model==3)
		
		j = ii[2]
		ZT.int = grep('.rel.ampl.int', colnames(T))
		ZT.ex = grep('.rel.ampl.ex', colnames(T))
		zt = seq(0,46,by = 2)
		zt.p = seq(0,48, by=0.5)
		
		source('model_modify_new/functions.R')
		#source('functions.R')
		
		intense.debug = FALSE
		
		param.fits.results = make.fits.with.all.models.for.one.gene(T = T, gene.index = j, debug = TRUE, parametrization = 'sigmoid', method = 'integration', underdetermination.model3=TRUE, sum.species = FALSE, zt = zt, i.ex = ZT.ex, i.int = ZT.int, absolute.signal = FALSE);
				
		param.fits.results2 = make.fits.with.all.models.for.one.gene(T = T, gene.index = j, debug = TRUE, parametrization = 'sigmoid', method = 'integration', underdetermination.model3=FALSE, sum.species = FALSE, zt = zt, i.ex = ZT.ex, i.int = ZT.int, absolute.signal = FALSE);
		
		param.fits.results
		param.fits.results2
		
		my.BIC.this.gene(param.fits.results)
		my.BIC.this.gene(param.fits.results2)
		
		
		
				
		best.model = 2
		i = j
		eval(parse(text = paste('gamma = T$gamma.m',best.model,'[i];fold.change = T$fold.change.int.m',best.model,'[i]; phase.int = T$phase.int.m',best.model,'[i]; up.time = T$up.time.int.m', best.model,'[i]; down.time = T$up.time.int.m', best.model,'[i]', sep = '')))
		eps.gamma = 1
		phase.gamma = 0
		#eval(parse(text = paste('fold.change.t = T$fold.change.int[i]; phase.int.t = T$phase.int[i]; up.time.t = T$up.time.int[i]; down.time.t = T$down.time.int[i]',sep= '')))
		
		m2 = compute.m.sigmoid(t = zt.p, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)
		s2 = compute.sigmoid(t = zt.p, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)
		
		best.model = 4
		i = j
		eval(parse(text = paste('gamma = T$gamma.m',best.model,'[i];eps.gamma = T$eps.gamma.m',best.model,'[i];phase.gamma = T$phase.gamma.m',best.model,'[i];fold.change = T$fold.change.int.m',best.model,'[i]; phase.int = T$phase.int.m',best.model,'[i]; up.time = T$up.time.int.m', best.model,'[i]; down.time = T$up.time.int.m', best.model,'[i]', sep = '')))
		
		#eval(parse(text = paste('fold.change.t = T$fold.change.int[i]; phase.int.t = T$phase.int[i]; up.time.t = T$up.time.int[i]; down.time.t = T$down.time.int[i]',sep= '')))
		
		m4 = compute.m.sigmoid(t = zt.p, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)
		s4 = compute.sigmoid(t = zt.p, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)
		
		i = j
		eval(parse(text = paste('gamma = T$gamma[i];eps.gamma = T$eps.gamma[i];phase.gamma = T$phase.gamma[i];fold.change = T$fold.change.int[i]; phase.int = T$phase.int[i]; up.time = T$up.time.int[i]; down.time = T$up.time.int[i]', sep = '')))
		
		mt = compute.m.sigmoid(t = zt.p, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)
		st = compute.sigmoid(t = zt.p, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)

		plot(zt, as.numeric(T[j,ZT.int]), col='blue',type='b',lty=1)
		points(zt, as.numeric(T[j,ZT.ex]), col='blue', type='b',lty=1)
		points(zt.p, mt, col='red', type='l',lty=1,lwd=2.0)
		points(zt.p, st, col='red', type='l',lty=1, lwd=2.0)
		points(zt.p, m2, col='darkgreen', type='l',lty=2)
		points(zt.p, s2, col='darkgreen', type='l',lty=2)
		points(zt.p, m4, col='black', type='l',lty=3)
		points(zt.p, s4, col='black', type='l',lty=3)
		abline(h=1, col='darkgray')
		
				
	}
}


## check the underdeterminations
load.for.plot = TRUE
plot=TRUE
palette(c('gray','green3','tomato','black'))
if(plot)
{
	cat('plot : FITS RESULTS FOR SIMULATED DATA\n')
	if(load.for.plot){load(file = "my_simulated_data_with_fits_results_LRT_and_AIC_BIC_all_versions_v2.Rdata")}
	
	source('functions.R')
	sd.noise = c(0.1,0.25,0.5)
	sd = 0.25
	
	eval(parse(text = paste('T = T',100*sd,sep ='')))
	#T$LRT.best.model = T$best.model
	index = c(1:400)
	T$true.model = ceiling(index/100)
	
	exon = as.matrix(T[,c(2:25)])
	intron = as.matrix(T[,c(26:49)])
	T$exon.var = apply(log(exon), 1, sd)
	T$intron.var = apply(log(intron), 1, sd)
	
	xx = cleaning.model.selection.results(T=T)
	ex.int = match(c("gene", "exon.var", "intron.var", "rel.ampl.ex", "phase.ex", "rel.ampl.int", "phase.int","Qval.ex", "Qval.int", "gamma", "eps.gamma", "phase.gamma","fold.change.int", "up.time.int", "down.time.int", 
					 "error.m1", "error.m2", "gamma.m2", "fold.change.int.m2","phase.int.m2", "up.time.int.m2", "down.time.int.m2",
					 "error.m4", "gamma.m4", "eps.gamma.m4","phase.gamma.m4","fold.change.int.m4", "phase.int.m4","up.time.int.m4","down.time.int.m4",       
					 "BIC.m1", "BIC.m2", "BIC.m3","BIC.m4", "BIC.best.model", "BIC.best.model.eff", "true.model"), colnames(xx))
	
	kk = which(is.na(xx$BIC.best.model.eff)==TRUE & xx$BIC.best.model==2)
	test = xx[kk, ex.int]
	
	plot(c(0:23)*2, intron[1,], type='b',col='blue')
	
	pdf('myplots/underdetermination_model2_v2.pdf', width = 6.6, height = 9)
	par(mfcol = c(4,2), las = 1, cex = 0.6, mar = c(3,3,2.2,0.5), mgp = c(1.7,0.5,0), tcl = -0.3)	
	
	genes = T[kk,1]
	label = LETTERS[1:9]
	zt = seq(0,46,by = 2); 
	zt.p = seq(0,48,by = 0.05)
	attribute.global.variable.colors()
	
	for(ii in 1:length(genes))
	{
		gene = genes[ii]
		cat(gene,'\n')
		i = T$gene == gene; 
		
		best.model = unlist(T$BIC.best.model[i]);
		ex = unlist(T[i,grep('.rel.ampl.ex', colnames(T))]); 
		int = unlist(T[i,grep('.rel.ampl.int', colnames(T))]);
		
		## calculate the profile of mrna
		eval(parse(text =paste('gamma = T$gamma.m',best.model,'[i]',sep = '')));
		#eval(parse(text =paste('eps.gamma = T$eps.gamma.m',best.model,'[i]',sep = '')));
		#eval(parse(text =paste('phase.gamma = T$phase.gamma.m',best.model,'[i]',sep = '')));
		eps.gamma = 0; 
		phase.gamma = 0; 
		gamma.t = T$gamma[i]
		eps.gamma.t = 0
		phase.gamma.t = 0
		
		eval(parse(text = paste('fold.change = T$fold.change.int.m',best.model,'[i]; phase.int = T$phase.int.m',best.model,'[i]; up.time = T$up.time.int.m', best.model,'[i]; down.time = T$down.time.int.m', best.model,'[i]', sep = '')))
		eval(parse(text = paste('fold.change.t = T$fold.change.int[i]; phase.int.t = T$phase.int[i]; up.time.t = T$up.time.int[i]; down.time.t = T$down.time.int[i]',sep= '')))
		
		
		m = compute.m.sigmoid(t = zt, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)
		m.t = compute.m.sigmoid(t = zt, gamma = gamma.t, eps.gamma = eps.gamma.t, phase.gamma = phase.gamma.t, fold.change = fold.change.t, phase = phase.int.t, up.time = up.time.t, down.time = down.time.t)
		s = compute.sigmoid(t = zt, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)
		s.t = compute.sigmoid(t = zt, fold.change = fold.change.t, phase = phase.int.t, up.time = up.time.t, down.time = down.time.t)
		
		intense.debug = FALSE
		S = int
		M = ex
		model = 2
		par = c(gamma, fold.change, phase.int, up.time, down.time)
		par.t = c(gamma.t, fold.change.t, phase.int.t, up.time.t, down.time.t)
		err = errors(S =S ,s = s, M = M, m = m, model = 2, param = par, parametrization =c('sigmoid'), log = TRUE); 
		
		err.t = errors(S =S ,s = s.t, M = M, m = m.t, model = 2, param = par, parametrization =c('sigmoid'), log = TRUE); 
		print('true parameters:');
		print(c(gamma.t, fold.change.t, phase.int.t, up.time.t, down.time.t))
		print('estimated parameters:');
		print(c(gamma, fold.change, phase.int, up.time, down.time))
		print(c(T$error.m2[i], err, err.t))
		
	}
		ylim = range(int)
		plot(1,1, type = 'n', xlab = 'time (hours)', axes = FALSE, xlim = range(zt),ylim=ylim, main=paste(gene,': ', which(T$gene==gene), sep=''));
		box(); 
		#time.axis(dt = dt); 
		axis(2)
		abline(h = 1, col = 'gray')
		#points(zt, ex, col = col.ex, type = 'p', lwd = 2, pch = 16); 
		points(zt, int, col = 'black', type = 'b', lwd = 1, pch = 1)
		#eval(parse(text =paste('gamma = T$gamma.m',best.model,'[i]',sep = '')));
		#eval(parse(text =paste('eps.gamma = T$eps.gamma.m',best.model,'[i]',sep = '')));
		#eval(parse(text =paste('phase.gamma = T$phase.gamma.m',best.model,'[i]',sep = '')));
		if(best.model > 1)
		{
			profile.int = compute.sigmoid(t = zt.p, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)
			profile.int.t = compute.sigmoid(t = zt.p, fold.change = fold.change.t, phase = phase.int.t, up.time = up.time.t, down.time = down.time.t)
		}
		points(zt.p, profile.int, type='l', lwd=1.5, col='darkred')
		points(zt.p, profile.int.t, type='l', lwd=1.5, col='darkblue')

	}
	dev.off()
	
	
	
	
}

#################################################################################### finishing line for the simulated data

###########################################################
###########################################################  
###########################################################

# FITS WITH REAL DATA


########
#### TEST of the fits
if(my.do.step)
{
	cat('TEST THE FIT on THE REAL DATA\n')
	load('my_genes_RAnalysis_larger_exons.Rdata'); ## eliminate genes in which exonic singal is smaller than intronic signal
	#load('my_genes_RAnalysis.Rdata') ## all expressed genes
	source('functions.R')
	#load('my_genes_RAnalysis.Rdata')
	## change the colnames of Tt
	names = colnames(tt)[grep('.rel.ampl',colnames(tt))]
	colnames(Tt)[2:25] = paste(names, 'ex', sep='.')
	colnames(Tt)[26:49] = paste(names, 'int', sep='.')
	
	#save(Tt, tt, file='my_genes_RAalysis_all.Rdata')
	save(Tt, tt, file='my_genes_RAnalysis_larger_exons_2.Rdata')
	
	ZT.int = grep('.rel.ampl.int', colnames(Tt))
	ZT.ex = grep('.rel.ampl.ex', colnames(Tt))
	zt = seq(0,46,by = 2)
	intense.debug = FALSE
	set.global.variable();
	k = splicing.rate;
	
	
	
	s = compute.sigmoid(t = zt)
	m = compute.m.sigmoid(t = zt, gamma = log(2)/0.5,eps.gamma = 10)
	m2 = compute.m.sigmoid(t = zt, gamma = log(2)/0.5,eps.gamma = 10, rescale = FALSE, synthesis.factor = k)
	
	par(mfrow = c(1,2))
	plot(s, type = 'b', ylim = range(0,s,m,m2, s+m2))
	points(m, type = 'b', col = 'red')
	points(m2, type = 'b', col = 'blue')
	points(s+m2, type = 'b', col = 'green')

	plot(s, type = 'b', ylim = range(0,s,m))
	points(m, type = 'l', col = 'red')
	points((s+m2)/mean(s+m2), type = 'l', col = 'green')
	
	mean(m2/s)

	param.fits.results.no.sum = make.fits.with.all.models.for.one.gene(T = T, gene.index = which(T$gene == 'Dbp'), debug = TRUE, parametrization  = 'sigmoid', method = 'integration')
	param.fits.results.with.sum = make.fits.with.all.models.for.one.gene(T = T, gene.index = which(T$gene == 'Dbp'), debug = TRUE, parametrization  = 'sigmoid', method = 'integration', sum.species = TRUE)
	
	source('functions.R')
	param.fits.results.with.sum_15min = make.fits.with.all.models.for.one.gene(T = T, gene.index = which(T$gene == 'Dbp'), debug = TRUE, parametrization  = 'sigmoid', method = 'integration', sum.species = TRUE)
	source('functions.R')
	set.global.variable();
	cat(splicing.rate)
	param.fits.results.with.sum_8min = make.fits.with.all.models.for.one.gene(T = T, gene.index = which(T$gene == 'Dbp'), debug = TRUE, parametrization  = 'sigmoid', method = 'integration', sum.species = TRUE)

	
	T.before = T;
	load('Rdata/genes_fit_results_and_LRT_and_AIC.Rdata')
	T.after = T;
	T = T.before
	
	T.after[which(T$gene == 'Dbp'), 74:92]
}


plot = TRUE; 
load.for.plot = TRUE
do.step = TRUE; 
load = TRUE

######## FITS ON VITAL-IT

version = '_my_5th_same_up_down_time'

if(do.step)
{
	cat('FIT THE REAL DATA ON VITAL-IT\n')
	load(file="my_genes_RAnalysis_larger_exons_final.Rdata"); ## eliminate genes in which exonic singal is smaller than intronic signal and also replace FISH_24 test by f24_R2_alt2 function
	#load('my_genes_RAnalysis.Rdata') ## all expressed genes
	T = Tt
	load(file="Sigma_my_genes_real_data.Rdata");
	
	#prepare the data for the fits
	save(T, file = 'stuff_for_vital_IT/genes_ready_for_fits.Rdata')
	gene_list = T$gene; 
	#gene_list = gene_list[1:10]
	write.table(gene_list, file = 'stuff_for_vital_IT/gene_list.txt', quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE)

	## ON VITAL-IT 

	#clear the directory
	cmd = "ssh jwang@dev.vital-it.ch 'rm -r /scratch/cluster/monthly/jingkui/microarray/*' " 
	system(cmd)

	#copy the data and the code on vital-IT
	cmd = 'scp -r stuff_for_vital_IT/* jwang@dev.vital-it.ch:/scratch/cluster/monthly/jingkui/microarray/.'
	system(cmd)
	cmd = 'scp -r model_modify/functions.R jwang@dev.vital-it.ch:/scratch/cluster/monthly/jingkui/microarray/.'
	system(cmd)

	# Insted of run fitting in R, we run it directly in Vital-it.
	# run the fits : bash -s <to_run_on_vital_IT.sh
	#cmd = "ssh jwang@dev.vital-it.ch 'bash -s' < stuff_for_vital_IT/to_run_on_vital_IT.sh"
	#system(cmd)
	
	# concatenate all fits results together when fits are all done
	cmd = "ssh jwang@dev.vital-it.ch 'perl /scratch/cluster/monthly/jingkui/microarray/cat_all_fits_results.pl'"
	system(cmd)

	# wait for everything to be done and copy the data on the local computer
	cmd = "ssh jwang@dev.vital-it.ch 'ls -l /scratch/cluster/monthly/jingkui/microarray/fits/all_genes_fits_results_with_header.txt'"
	res = system(cmd, intern = TRUE)

	while(length(res) == 0)
	{
		Sys.sleep(60*5) # wait for 5 minutes
		cat('... wait for the fits to finish\n')
		cmd = "ssh jwang@dev.vital-it.ch 'ls -l /scratch/cluster/monthly/jingkui/microarray/fits/all_genes_fits_results_with_header.txt'"
		res = system(cmd, intern = TRUE)
	}
	

	cat('........................................ fits finished on vital-IT\n')	

	cmd = paste("scp jwang@dev.vital-it.ch:/scratch/cluster/monthly/jingkui/microarray/fits/all_genes_fits_results_with_header.txt Vital-IT_fits_optim/all_genes_fits_results_with_header",version,".txt", sep = '')
	system(cmd)

	cat('........................................ fits results imported on local computer\n')
	
	file =  paste("Vital-IT_fits_optim/all_genes_fits_results_with_header",version,".txt", sep = '')	
	fits.results = read.table(file = file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
	
	if(nrow(fits.results)<0.9*length(gene_list)){cat('........................................ fits results missing for more than 10% of genes\n')}
	m = match(gene_list, fits.results$gene)
	if(sum(is.na(m))>0){cat('........................................ fits results missing for',sum(is.na(m)),'genes\n')}
	fits.results = fits.results[m,]

	Tt = data.frame(Tt, fits.results[,-1],stringsAsFactors = FALSE)	
	save(Tt,tt, file = paste('my_genes_fits_optim_vital-IT',version,'.Rdata',sep = ''))	

}


## Different model selecting methods for real data

my.do.step = TRUE
load.two.version = FALSE

version_old = '_my_1st'
version_new = '_my_4nd'

version = '_my_7th_same_up_down_time_all_genes'
version = '_my_8th_same_up_down_time_all_genes_model3_no_treated'
version = '_my_8th_same_up_down_time_all_genes_sum'

if(my.do.step)
{
	cat('LRT, AIC, BIC WITH REAL DATA\n')
	
	load("Sigma_my_genes_real_data.Rdata")
	source('model_modify_new/functions.R')
	
	## The only sigma of real data without distinction of exons and introns	
	fact = sigma.ms.real
	pval = 0.05
	
	#load('my_simulated_data_with_fits_results_diff_sigma_all_Old_Method_T10.Rdata')
	if(!load.two.version)
	{
		load(file = paste('my_genes_fits_optim_vital-IT',version,'.Rdata',sep = ''));
		T = Tt
		kk = which(is.na(T$error.m1)==TRUE)
		if(length(kk)>0) T = T[-kk,]
		
		LRT = Likelihood.Ratio.Test.LAURA(T = T, pval = pval)
		LRT.FDR = Likelihood.Ratio.Test.LAURA(T = T, pval = pval, FDR=TRUE)
		colnames(LRT.FDR) = paste(colnames(LRT.FDR), '.FDR', sep='');
		
		AIC = AIC.LAURA(T = T, correction = FALSE)
		AICc = AIC.LAURA(T = T, correction = TRUE)
		BIC = BIC.LAURA(T = T)
		T = data.frame(T, AIC, AICc, BIC, LRT, LRT.FDR, stringsAsFactors = FALSE);
		Tt = T
		save(Tt, file = paste('my_genes_fits_results', version,'.Rdata', sep=''))
	}
	
	if(load.two.version)
	{
		load(file = paste('my_genes_fits_optim_vital-IT',version_old,'.Rdata',sep = ''));
		TT = Tt
		rm(list=c("Tt", "tt"));
		load(file = paste('my_genes_fits_optim_vital-IT',version_new,'.Rdata',sep = ''));
		T = Tt
		m = match(TT$gene, T$gene)
		TT = TT[m,]
		
		## likelihood ratio test
		LRT.best.model = my.Likelihood.Ratio.Test(T = T, pval = pval)
		
		LRT.best.model.old = my.Likelihood.Ratio.Test.Old(T = TT, pval = pval, fact=fact)
		
		LRT.best.model.laura = Likelihood.Ratio.Test.LAURA(T = TT, pval = pval)
		LRT.best.model.laura.fdr = Likelihood.Ratio.Test.LAURA(T = TT, pval = pval, FDR=TRUE)
		colnames(LRT.best.model.laura.fdr) = paste(colnames(LRT.best.model.laura.fdr), '.FDR', sep='');
		
		## AIC
		AIC = my.AIC(T = T, correction = FALSE)
		AIC.old = my.AIC.Old(T=TT, correction=FALSE, fact=fact)
		AIC.laura = AIC.LAURA(T = TT,correction=FALSE)
		
		## AIC with correction
		AICc = my.AIC(T = T, correction = TRUE)
		AICc.old = my.AIC.Old(T=TT, correction=TRUE, fact=fact)
		AICc.laura = AIC.LAURA(T = TT, correction=TRUE)
		
		## BIC
		BIC = my.BIC(T = T)
		BIC.old = my.BIC.Old(T = TT, fact =fact)
		BIC.laura = BIC.LAURA(T=TT)
		
		## Combining results by methods of different versions
		
		ii = match(c("error.m1", "error.m2","gamma.m2","fold.change.int.m2", "phase.int.m2", "up.time.int.m2", "down.time.int.m2","error.m3", "gamma.m3", "eps.gamma.m3", "phase.gamma.m3",
					 "error.m4", "gamma.m4", "eps.gamma.m4", "phase.gamma.m4", "fold.change.int.m4", "phase.int.m4", "up.time.int.m4","down.time.int.m4"), colnames(TT));
		
		colnames(TT)[ii] = paste(colnames(TT)[ii], '.old', sep='')
		
		T = data.frame(T, AIC, AICc, BIC, LRT.best.model, TT[,ii], AIC.old, AICc.old, BIC.old, LRT.best.model.old,
					   AIC.laura, AICc.laura, BIC.laura, LRT.best.model.laura, LRT.best.model.laura.fdr, stringsAsFactors = FALSE);
		
		Tt = T;
		
		save(Tt, tt, file = "my_genes_fits_results_LRT_and_AIC_BIC_all_Laura_Old_New_v2.Rdata")
	}
	
}

####
#### Check if the optimization and model selection work (in particular to deal with the optimization issues)
####
version = '_LRT_and_AIC_BIC_all_Laura_Old_New_v2'
version = '_my_5th_same_up_down_time'
version = '_my_6th_same_up_down_time_problem_genes'
version = '_my_7th_same_up_down_time_all_genes'

version = '_my_8th_same_up_down_time_all_genes_model3_no_treated'
version = '_my_8th_same_up_down_time_all_genes_sum'

my.do = TRUE
if(my.do)
{
	version.compare = FALSE
	if(version.compare)
	{
		source('model_modify_new/functions.R')
		version.1 = '_my_8th_same_up_down_time_all_genes_model3_no_treated'
		version.2 = '_my_8th_same_up_down_time_all_genes_sum'
		load(paste('my_genes_fits_results', version.1,'.Rdata', sep = ''))
		T1 = Tt
		load(paste('my_genes_fits_results', version.2,'.Rdata', sep = ''))
		T2 = Tt
		
		m = match(T2$gene, T1$gene)
		T1 = T1[m,]
		plot(as.numeric(T2$BIC.best.model+rnorm(nrow(T2),sd = 0.1)), as.numeric(T1$BIC.best.model)+rnorm(nrow(T1),sd = 0.1), xlab='Classification with sum', 
			 ylab='Classification without sum', main='Comprison of results with or without sum', pch = 16, cex = 1, col = rgb(0,0.5,0.8,0.3))
		
		length(which(T1$BIC.best.model==1))
		length(which(T1$BIC.best.model==2))
		length(which(T1$BIC.best.model==3))
		length(which(T1$BIC.best.model==4))
		
		length(which(T2$BIC.best.model==1))
		length(which(T2$BIC.best.model==2))
		length(which(T2$BIC.best.model==3))
		length(which(T2$BIC.best.model==4))
		
	}
	load(paste('my_genes_fits_results', version,'.Rdata', sep = ''))
	source('model_modify_new/functions.R')
	
	T = Tt
	#Tt$BIC.best.model = Tt$BIC.best.model.Laura
	gamma.max = log(2)/10*60 ; 
	gamma.min = log(2)/24;
	
	### uptime issue solved
	length(which(Tt$BIC.best.model==2))
	length(which(Tt$BIC.best.model==3))
	length(which(Tt$BIC.best.model==4))
	
	length(which((Tt$BIC.best.model==2 & Tt$up.time.int.m2==2)))
	length(which((Tt$BIC.best.model==4 & Tt$up.time.int.m4==2)))
		
	#### Half-life estimated issue
	length(which(Tt$BIC.best.model==2 & (abs(Tt$gamma.m2-gamma.max)<0.00001)))
	length(which(Tt$BIC.best.model==3 & (abs(Tt$gamma.m3-gamma.max)<0.00001)))
	length(which(Tt$BIC.best.model==4 & (abs(Tt$gamma.m4-gamma.max)<0.00001)))
	
	hh = c(which(Tt$BIC.best.model==2 & (abs(Tt$gamma.m2-gamma.max)<0.00001)), which(Tt$BIC.best.model==3 & (abs(Tt$gamma.m3-gamma.max)<0.00001)), which(Tt$BIC.best.model==4 & (abs(Tt$gamma.m4-gamma.max)<0.00001)))
	#hh = c(which(Tt$BIC.best.model==2 & (abs(Tt$gamma.m2-gamma.max)<0.00001|abs(Tt$gamma.m2-gamma.min)<0.00001)), which(Tt$BIC.best.model==4 & (abs(Tt$gamma.m4-gamma.max)<0.00001|abs(Tt$gamma.m4-gamma.min)<0.00001)))
	jj = c(which(Tt$BIC.best.model==2 & (abs(Tt$gamma.m2-gamma.max)>0.00001)),  which(Tt$BIC.best.model==3 & (abs(Tt$gamma.m3-gamma.max)<0.00001)), which(Tt$BIC.best.model==4 & (abs(Tt$gamma.m4-gamma.max)>0.00001)))
	
	
	ex.int = match(c("gene", "exon.median", "intron.median", "exon.var", "intron.var", "rel.ampl.ex", "rel.ampl.int","pval.ex", "pval.int","phase.ex", "phase.int",
					 "gamma.m2", "gamma.stderr.m2", "phase.int.m2", "gamma.m4","eps.gamma.m4", "gamma.stderr.m4","eps.gamma.stderr.m4","phase.gamma.m4","BIC.best.model"), colnames(Tt))
	test = Tt[hh, ex.int]
	test1 = Tt[jj, ex.int]
	aa = abs((as.numeric(Tt$phase.ex[hh])-as.numeric(Tt$phase.int[hh]))%%24)
	bb = abs((Tt$phase.ex[jj]-Tt$phase.int[jj])%%24)
	
	par(mfcol = c(3,2))
	hist(aa,breaks=48)
	#hist(Tt$rel.ampl.ex[hh],breaks=20)
	#hist(Tt$rel.ampl.int[hh],breaks=20)
	plot(aa, Tt$rel.ampl.ex[hh]/Tt$rel.ampl.int[hh],cex=0.5)
	abline(h=1, col='red')
	hist(Tt$rel.ampl.ex[hh]/Tt$rel.ampl.int[hh],breaks=50)
	abline(v=1, col='red')
	#plot(aa, Tt$rel.ampl.int[hh],cex=0.5)
	
	
	### Test whether the modified model works
	ZT.int = grep('.rel.ampl.int', colnames(T))
	ZT.ex = grep('.rel.ampl.ex', colnames(T))
	zt = seq(0,46,by = 2)
	zt.p = seq(0,48, by=0.5)
	
	#j = which(Tt$gene=='Nup155')
	mm = c(which(Tt$BIC.best.model==2 & (abs(Tt$gamma.m2-gamma.max)<0.00001)), which(Tt$BIC.best.model==3 & (abs(Tt$gamma.m3-gamma.max)<0.00001)), which(Tt$BIC.best.model==4 & (abs(Tt$gamma.m4-gamma.max)<0.00001)))
	nn = c(which(Tt$BIC.best.model==2 & (abs(Tt$gamma.m2-gamma.max)>0.00001)),  which(Tt$BIC.best.model==3 & (abs(Tt$gamma.m3-gamma.max)<0.00001)), which(Tt$BIC.best.model==4 & (abs(Tt$gamma.m4-gamma.max)>0.00001)))
	
	j = mm[1]	
	source('model_modify_new/functions.R')
	#source('functions.R')
	
	intense.debug = FALSE
	
	param.fits.results = make.fits.with.all.models.for.one.gene(T = T, gene.index = j, debug = TRUE, parametrization = 'sigmoid', method = 'integration', underdetermination.model3=FALSE, sum.species = FALSE, zt = zt, i.ex = ZT.ex, i.int = ZT.int, absolute.signal = FALSE);
	
	param.fits.results2 = make.fits.with.all.models.for.one.gene(T = T, gene.index = j, debug = TRUE, parametrization = 'sigmoid', method = 'integration', underdetermination.model3=FALSE, sum.species = TRUE, zt = zt, i.ex = ZT.ex, i.int = ZT.int, absolute.signal = FALSE);
	
	param.fits.results
	param.fits.results2
	
	my.BIC.this.gene(param.fits.results)
	my.BIC.this.gene(param.fits.results2)
	
	
	## check the half lives estimated by the ratio of exons and introns 
	#load('my_genes_RAnalysis.Rdata') ## all expressed genes
	shar = read.table('/Users/jiwang/Degradation_Laura/_x_MicroArray/complementary_data/half-lives/Sharova.txt', stringsAsFactors = FALSE, sep = '\t', header = TRUE)
	fried = read.table('/Users/jiwang/Degradation_Laura/_x_MicroArray/complementary_data/half-lives/Friedel.txt', stringsAsFactors = FALSE, sep = '\t', header = TRUE)
	colnames(shar) = c('gene', 'hl')
	colnames(fried) = c('gene', 'probeset.ID','hlmin','hl')
	
	par(mfcol = c(1,2))
	ii = match(fried$gene, Tt$gene)
	length(which(!is.na(ii)))
	
	jj = which(!is.na(ii))
	ii = ii[which(!is.na(ii))]
	plot(Tt$exon.median[ii]/Tt$intron.median[ii], fried[jj,4], log='xy',cex=0.5)
	res1 = lm(Tt$exon.median[ii]/Tt$intron.median[ii] ~ fried[jj,4]- 1)
	
	ii = match(shar$gene, Tt$gene)
	length(which(!is.na(ii)))
	
	jj = which(!is.na(ii))
	ii = ii[which(!is.na(ii))]
	plot(Tt$exon.median[ii]/Tt$intron.median[ii], shar[jj,2], log='xy',cex=0.5)
	res2 = lm(Tt$exon.median[ii]/Tt$intron.median[ii] ~ shar[jj,2]- 1)
	
	
	hist(bb,breaks=48)
#hist(Tt$rel.ampl.ex[jj],breaks=20)
#hist(Tt$rel.ampl.int[jj],breaks=20)
	plot(bb, Tt$rel.ampl.ex[jj]/Tt$rel.ampl.int[jj],cex=0.5)
	abline(h=1, col='red')
	hist(Tt$rel.ampl.ex[jj]/Tt$rel.ampl.int[jj],breaks=50)
	abline(v=1, col='red')
#plot(bb, Tt$rel.ampl.ex[jj],cex=0.5)
#plot(bb, Tt$rel.ampl.int[jj],cex=0.5)
#T = Tt
#T$BIC.best.model = T$BIC.best.model
	if(plot)
	{
		pdf('myplots/underdetermination_real_data_model24_improved_fitting_separatly.pdf', width = 6.6, height = 9)
		par(mfcol = c(2,1), las = 1, cex = 0.6, mar = c(3,3,2.2,0.5), mgp = c(1.7,0.5,0), tcl = -0.3)	
		
#kk = which((Tt$BIC.best.model==2|Tt$BIC.best.model==4)  & (abs(Tt$gamma.m2-gamma.max)<0.00001|Tt$up.time.int.m2==22|Tt$up.time.int.m2==2|Tt$down.time.int.m2==22|Tt$down.time.int.m2==2))
		kk = mm
		genes = T[kk,1]
		label = LETTERS[1:9]
		zt = seq(0,46,by = 2); 
		zt.p = seq(0,48,by = 0.5)
		attribute.global.variable.colors()
		
		for(ii in 1:length(genes))
		{
			gene = genes[ii]
			cat(ii,'\n')
			i = T$gene == gene; 
			
			best.model = unlist(T$BIC.best.model[i]);
			ex = unlist(T[i,grep('.rel.ampl.ex', colnames(T))]); 
			int = unlist(T[i,grep('.rel.ampl.int', colnames(T))]);
			
			eval(parse(text = paste('gamma = T$gamma.m',best.model,'[i]; fold.change = T$fold.change.int.m',best.model,'[i]; phase.int = T$phase.int.m',best.model,'[i]; up.time = T$up.time.int.m', best.model,'[i]; down.time = T$up.time.int.m', best.model,'[i]', sep = '')))
#eval(parse(text =paste('gamma = T$gamma.m',best.model,'[i]',sep = '')));
			
			if(best.model==2) {eps.gamma = 1; phase.gamma = 0;}
			if(best.model==4){eval(parse(text =paste('eps.gamma = T$eps.gamma.m',best.model,'[i]',sep = '')));eval(parse(text =paste('phase.gamma = T$phase.gamma.m',best.model,'[i]',sep = '')));}
			
			
#eval(parse(text = paste('fold.change = T$fold.change.int.m',best.model,'[i]; phase.int = T$phase.int.m',best.model,'[i]; up.time = T$up.time.int.m', best.model,'[i]; down.time = T$up.time.int.m', best.model,'[i]', sep = '')))
#eval(parse(text = paste('fold.change.t = T$fold.change.int[i]; phase.int.t = T$phase.int[i]; up.time.t = T$up.time.int[i]; down.time.t = T$down.time.int[i]',sep= '')))
			
			
			m = compute.m.sigmoid(t = zt.p, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)
			s = compute.sigmoid(t = zt.p, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)
			
			ylim = range(c(int,s))
			plot(1,1, type = 'n', xlab = 'time (hours)', axes = FALSE, xlim = range(zt),ylim=ylim, main=paste('Intron: ', gene,': index ', which(T$gene==gene),', model:', best.model, ', half-life: ', round(log(2)/gamma,d=3), ', up and down time: ', round(up.time,d=3), ', ', round(down.time,d=3), sep=''));
			box(); 
			time.axis(dt = 4); 
			axis(2)
			abline(h = 1, col = 'darkgray')
#points(zt, ex, col = col.ex, type = 'p', lwd = 2, pch = 16); 
			points(zt, int, col = 'black', type = 'b', lwd = 1, pch = 1)
			points(zt.p, s, type='l', lwd=1.5, col='darkred')
			
			
			ylim = range(c(ex,m))
#ylim = c(0,2.5)
			plot(1,1, type = 'n', xlab = 'time (hours)', axes = FALSE, xlim = range(zt),ylim=ylim, main=paste('Exon: ',gene,': index ', which(T$gene==gene),', model:', best.model, ', half-life: ', round(log(2)/gamma,d=3), ', up and down time: ', round(up.time,d=3), ', ', round(down.time,d=3), sep=''));
			box(); 
			time.axis(dt = 4); 
			axis(2)
			abline(h = 1, col = 'darkgray')
#points(zt, ex, col = col.ex, type = 'p', lwd = 2, pch = 16); 
			points(zt, ex, col = 'steelblue', type = 'b', lwd = 1, pch = 1)
			points(zt.p, m, type='l', lwd=1.5, col='steelblue2')
			
			
		}
		dev.off()
	}
	
#Ttt = cleaning.model.selection.results(T=Tt)
	
## boundary of parameters
	gamma.max = log(2)/10*60 ; 
	gamma.min = log(2)/24;
	
	ii = which((Tt$BIC.best.model==2| Tt$BIC.best.model==2) & Tt$Qval.int<0.1)
	length(ii)
	length(which((Tt$BIC.best.model==2| Tt$BIC.best.model==4)))
	
#ii = which((Tt$BIC.best.model==2|Tt$BIC.best.model==4) & (Tt$up.time.int.m2==12|Tt$up.time.int.m2==4))
	ii = which((Tt$BIC.best.model==2|Tt$BIC.best.model==4) & (Tt$up.time.int.m2==4))
	
	
	
	test = T[ii, ex.int]
	
	length(which(Tt$BIC.best.model==2))
	length(which(Tt$BIC.best.model==2 & (Tt$up.time.int.m2==12|Tt$up.time.int.m2==2|Tt$down.time.int.m2==22|Tt$down.time.int.m2==2)))
	length(which(Tt$BIC.best.model==2 & (abs(Tt$gamma.m2-gamma.max)<0.00001|Tt$up.time.int.m2==22|Tt$up.time.int.m2==2|Tt$down.time.int.m2==22|Tt$down.time.int.m2==2)))
	
	length(which(Tt$BIC.best.model.Laura==2 & Tt$up.time.int.m2==4))
	length(which(Tt$BIC.best.model.Laura==2))
	length(which(Tt$BIC.best.model.Laura==4 & Tt$up.time.int.m4==4))
	length(which(Tt$BIC.best.model.Laura==4))
	
#mm = which((Tt$BIC.best.model.Laura==2 & Tt$up.time.int.m2==4)|(Tt$BIC.best.model.Laura==4 & Tt$up.time.int.m4==4))
	
	
	gene.prob = Tt$gene[mm]
	
	kk = which((Tt$BIC.best.model.Laura==2 & (Tt$up.time.int.m2.old==22|Tt$up.time.int.m2.old==2|Tt$down.time.int.m2.old==22|Tt$down.time.int.m2.old==2))
			   |(Tt$BIC.best.model.Laura==4 & (Tt$up.time.int.m4.old==22|Tt$up.time.int.m4.old==2|Tt$down.time.int.m4.old==22|Tt$down.time.int.m4.old==2)))
	
	gene.prob2 = Tt$gene[kk]
	gene.list = unique(c(gene.prob, gene.prob2))
	
	save(gene.list, file='genes_to_test.Rdata')
	
	test = T[mm, ex.int]
	
	ii = which((Tt$BIC.best.model==2|Tt$BIC.best.model==4) & (Tt$up.time.int.m2==22|Tt$up.time.int.m2==2|Tt$down.time.int.m2==22|Tt$down.time.int.m2==2))
	
	kk = which((Tt$BIC.best.model.Laura==2 & (Tt$up.time.int.m2.old==22|Tt$up.time.int.m2.old==2|Tt$down.time.int.m2.old==22|Tt$down.time.int.m2.old==2))|(Tt$BIC.best.model.Laura==4 & (Tt$up.time.int.m4.old==22|Tt$up.time.int.m4.old==2|Tt$down.time.int.m4.old==22|Tt$down.time.int.m4.old==2)))
	
	kk1 = which((Tt$BIC.best.model.Laura==2) & (Tt$up.time.int.m2.old==22|Tt$up.time.int.m2.old==2|Tt$down.time.int.m2.old==22|Tt$down.time.int.m2.old==2))
	kk2 = which((Tt$BIC.best.model.Laura==4) & (Tt$up.time.int.m4.old==22|Tt$up.time.int.m4.old==2|Tt$down.time.int.m4.old==22|Tt$down.time.int.m4.old==2))
	
	jj = which(Tt$BIC.best.model==2|Tt$BIC.best.model==4)
	jj = jj[which(is.na(match(jj, ii))==TRUE)]
	
	
	
	
	ex.int = match(c("gene", "exon.var", "intron.var", "rel.ampl.ex", "phase.ex", "rel.ampl.int", "phase.int","pval.ex", "pval.int",
					 "error.m1", "error.m2", "gamma.m2", "fold.change.int.m2","phase.int.m2", "up.time.int.m2", "down.time.int.m2","up.time.int.m2.old", "down.time.int.m2.old",
					 "error.m4", "gamma.m4", "eps.gamma.m4","phase.gamma.m4","fold.change.int.m4", "phase.int.m4","up.time.int.m4","down.time.int.m4", "up.time.int.m4.old","down.time.int.m4.old",      
					 "BIC.m1", "BIC.m2", "BIC.m3","BIC.m4", "BIC.best.model"), colnames(Tt))
	
	test = T[ii, ex.int]
	test2 = T[jj, ex.int]
	
## check if genes with optimization issue have large intron noise than those without
	load(file = "my_simulated_data_with_fits_results_LRT_and_AIC_BIC_all_versions_v2.Rdata")
	sd.noise = c(0.1,0.25,0.5)
	sd = 0.25
	eval(parse(text = paste('T.fake = T',100*sd,sep ='')))
	index = c(1:400)
	T.fake$true.model = ceiling(index/100)
	exon = as.matrix(T.fake[,c(2:25)])
	intron = as.matrix(T.fake[,c(26:49)])
	T.fake$exon.var = apply(log(exon), 1, sd)
	T.fake$intron.var = apply(log(intron), 1, sd)
	
	par(mfcol = c(1,3))
	lims = range(c(test2$intron.var, test$intron.var))
	hist(test2$intron.var, xlim=lims, breaks=100,main='Distribution of intron sd for real data with optimization issue');abline(v=0.215, col='red');abline(v=0.215*2, col='red');
	hist(test$intron.var, xlim=lims, breaks=100, col='darkgray', main='Distribution of intron sd for real data with optimization issue');abline(v=0.215, col='red');abline(v=0.215*2, col='red');
	hist(T.fake$intron.var[which(T.fake$true.model==2)], xlim=lims, breaks=20, col='darkgreen',main='Distribution of intron sd for simulated data');abline(v=0.27, col='red');abline(v=0.27*2, col='red');
	
	
#### check the pvals and variances for fake data
	sd.noise = c(0.1,0.25,0.5)
	for(sd in sd.noise)
	{
		eval(parse(text = paste('T.fake = T',100*sd,sep ='')))
		index = c(1:400)
		T.fake$true.model = ceiling(index/100)
		exon = as.matrix(T.fake[,c(2:25)])
		intron = as.matrix(T.fake[,c(26:49)])
		T.fake$exon.var = apply(log(exon), 1, sd)
		T.fake$intron.var = apply(log(intron), 1, sd)
		eval(parse(text = paste('T',100*sd, '=T.fake', sep ='')))
	}
	
	par(mfcol = c(1,2))
#kk = which((T10$BIC.best.model==2|T10$BIC.best.model==4) & (T10$up.time.int.m2==22|T10$up.time.int.m2==2|T10$down.time.int.m2==22|T10$down.time.int.m2==2))
	xlims = range(c(log10(T10$intron.var), log10(T25$intron.var), log10(T50$intron.var)))
	ylims = range(c(-log10(T10$pval.int),-log10(T25$pval.int),-log10(T50$pval.int)))
	
	T.fake = T10
	plot(log10(T.fake$intron.var[which(T.fake$true.model==1|T.fake$true.model==3)]), -log10(T.fake$pval.int[which(T.fake$true.model==1|T.fake$true.model==3)]),xlim=xlims, ylim=ylims, cex=0.7, col='darkgray',main='Introns fake')
	points(log10(T.fake$intron.var[which(T.fake$true.model==2|T.fake$true.model==4)]), -log10(T.fake$pval.int[which(T.fake$true.model==2|T.fake$true.model==4)]), cex=0.7, col='darkgreen')
	T.fake = T25
	points(log10(T.fake$intron.var[which(T.fake$true.model==1|T.fake$true.model==3)]), -log10(T.fake$pval.int[which(T.fake$true.model==1|T.fake$true.model==3)]), cex=0.7,pch=16, col='darkgray')
	points(log10(T.fake$intron.var[which(T.fake$true.model==2|T.fake$true.model==4)]), -log10(T.fake$pval.int[which(T.fake$true.model==2|T.fake$true.model==4)]), cex=0.7, pch=16, col='darkgreen')
	T.fake = T50
	points(log10(T.fake$intron.var[which(T.fake$true.model==1|T.fake$true.model==3)]), -log10(T.fake$pval.int[which(T.fake$true.model==1|T.fake$true.model==3)]), cex=0.7,pch=17, col='darkgray')
	points(log10(T.fake$intron.var[which(T.fake$true.model==2|T.fake$true.model==4)]), -log10(T.fake$pval.int[which(T.fake$true.model==2|T.fake$true.model==4)]), cex=0.7, pch=17, col='darkgreen')
	
	abline(5,-0.5,col='red', lwd=2.0)
	
	
	xlims = range(c(log10(T10$exon.var), log10(T25$exon.var), log10(T50$exon.var)))
	ylims = range(c(-log10(T10$pval.ex),-log10(T25$pval.ex),-log10(T50$pval.ex)))
	
	T.fake = T10
	plot(log10(T.fake$exon.var[which(T.fake$true.model==1)]), -log10(T.fake$pval.ex[which(T.fake$true.model==1)]),xlim=lims, ylim=ylims, cex=0.5, col=1,main='Exons fake')
	points(log10(T.fake$exon.var[which(T.fake$true.model==2)]), -log10(T.fake$pval.ex[which(T.fake$true.model==2)]), cex=0.5, col=2)
	points(log10(T.fake$exon.var[which(T.fake$true.model==3)]), -log10(T.fake$pval.ex[which(T.fake$true.model==3)]), cex=0.5, col=3)
	points(log10(T.fake$exon.var[which(T.fake$true.model==4)]), -log10(T.fake$pval.ex[which(T.fake$true.model==4)]), cex=0.5, col=4)
	T.fake = T25
	points(log10(T.fake$exon.var[which(T.fake$true.model==1)]), -log10(T.fake$pval.ex[which(T.fake$true.model==1)]), pch=16, cex=0.5, col=1)
	points(log10(T.fake$exon.var[which(T.fake$true.model==2)]), -log10(T.fake$pval.ex[which(T.fake$true.model==2)]), pch=16,cex=0.5, col=2)
	points(log10(T.fake$exon.var[which(T.fake$true.model==3)]), -log10(T.fake$pval.ex[which(T.fake$true.model==3)]), pch=16,cex=0.5, col=3)
	points(log10(T.fake$exon.var[which(T.fake$true.model==4)]), -log10(T.fake$pval.ex[which(T.fake$true.model==4)]), pch=16,cex=0.5, col=4)
	T.fake = T50
	points(log10(T.fake$exon.var[which(T.fake$true.model==1)]), -log10(T.fake$pval.ex[which(T.fake$true.model==1)]), pch=16, cex=0.5, col=1)
	points(log10(T.fake$exon.var[which(T.fake$true.model==2)]), -log10(T.fake$pval.ex[which(T.fake$true.model==2)]), pch=16,cex=0.5, col=2)
	points(log10(T.fake$exon.var[which(T.fake$true.model==3)]), -log10(T.fake$pval.ex[which(T.fake$true.model==3)]), pch=16,cex=0.5, col=3)
	points(log10(T.fake$exon.var[which(T.fake$true.model==4)]), -log10(T.fake$pval.ex[which(T.fake$true.model==4)]), pch=16,cex=0.5, col=4)
	
	
	
	
##
## check the pvals and vairances for introns and exons
##
	T.fake = T25
	index = c(1:400)
	T.fake$true.model = ceiling(index/100)
	exon = as.matrix(T.fake[,c(2:25)])
	intron = as.matrix(T.fake[,c(26:49)])
	T.fake$exon.var = apply(log(exon), 1, sd)
	T.fake$intron.var = apply(log(intron), 1, sd)
	
	ii = which((Tt$BIC.best.model.Laura==2|Tt$BIC.best.model.Laura==4) & (Tt$up.time.int.m2.old==22|Tt$up.time.int.m2.old==2|Tt$down.time.int.m2.old==22|Tt$down.time.int.m2.old==2))
	jj = which(Tt$BIC.best.model.Laura==2|Tt$BIC.best.model.Laura==4)
	jj = jj[which(is.na(match(jj, ii))==TRUE)]
	
	par(mfcol = c(2,2))
## fake data
	kk = which((T.fake$BIC.best.model==2|T.fake$BIC.best.model==4) & (T.fake$up.time.int.m2==22|T.fake$up.time.int.m2==2|T.fake$down.time.int.m2==22|T.fake$down.time.int.m2==2))
	lims = range(log10(T.fake$intron.var))
	plot(log10(T.fake$intron.var[which(T.fake$true.model==1|T.fake$true.model==3)]), -log10(T.fake$pval.int[which(T.fake$true.model==1|T.fake$true.model==3)]),xlim=lims, ylim=range(-log10(T.fake$pval.int)), cex=0.5, col='darkgray',main='Introns fake')
	points(log10(T.fake$intron.var[which(T.fake$true.model==2|T.fake$true.model==4)]), -log10(T.fake$pval.int[which(T.fake$true.model==2|T.fake$true.model==4)]), cex=0.5, col='darkgreen')
	points(log10(T.fake$intron.var[kk]), -log10(T.fake$pval.int)[kk], col='darkblue', cex=1.0,pch=16)
	
	lims = range(T.fake$exon.var)
	plot(T.fake$exon.var[which(T.fake$true.model==1|T.fake$true.model==3)], -log10(T.fake$pval.ex[which(T.fake$true.model==1|T.fake$true.model==3)]),xlim=lims, ylim=range(-log10(T.fake$pval.ex)), cex=0.5, col='darkgray',main='Exons fake')
	points(T.fake$exon.var[which(T.fake$true.model==2|T.fake$true.model==4)], -log10(T.fake$pval.ex[which(T.fake$true.model==2|T.fake$true.model==4)]), cex=0.5, col='darkgreen')
	points(T.fake$exon.var[kk], -log10(T.fake$pval.ex)[kk], col='darkblue', cex=1.0,pch=16)
	
##real data
	kk = which(Tt$BIC.best.model.Laura==1|Tt$BIC.best.model.Laura==3)
	kk1 = which(Tt$BIC.best.model.Laura==1)
	plot(log10(Tt$intron.var[kk]), -log10(Tt$pval.int[kk]), cex=0.5, col='darkgray', ylim=range(-log10(T.fake$pval.int)),main='Introns real')
	points(log10(Tt$intron.var[jj]), -log10(Tt$pval.int[jj]), col='darkgreen',cex=0.5,pch=16)
	points(log10(Tt$intron.var[ii]), -log10(Tt$pval.int[ii]), col='darkred',cex=0.2)
	abline(-1.,-4,col='red', lwd=2.0)
	
	plot(Tt$exon.var[kk], -log10(Tt$pval.ex[kk]), cex=0.5, col='darkgray', ylim=range(-log10(T.fake$pval.ex)),main='Exons real')
	points(Tt$exon.var[jj], -log10(Tt$pval.ex[jj]), col='darkgreen',cex=0.5,pch=16)
	points(Tt$exon.var[ii], -log10(Tt$pval.ex[ii]), col='darkred',cex=0.2)
	
	
### check the fake data whether the new method to estimate the noise sigma works
	noise.estimate = function(data)
	{
		nb = length(data)
		half = nb/2
		xx = data[1:half]-data[(half+1):nb]
		return(sd(xx)/sqrt(2))		
	}
## fake data 
	T.fake = T25
	exon = as.matrix(T.fake[,c(2:25)])
	intron = as.matrix(T.fake[,c(26:49)])
	sigma.ex = apply((exon), 1, noise.estimate)
	sigma.int = apply((intron), 1, noise.estimate)
	dens.ex = density(sigma.ex)
	dens.int = density(sigma.int)
	
## real data
	exon = as.matrix(Tt[,c(2:25)])
	intron = as.matrix(Tt[,c(26:49)])
	sigma.ex = apply(log(exon), 1, noise.estimate)
	sigma.int = apply(log(intron), 1, noise.estimate)
	
	lims1 = range(c(dens.int$y, hist(sigma.int, breaks=200, probability = TRUE)$density))
	lims2 = range(c(dens.ex$y, hist(sigma.ex, breaks=200, probability = TRUE)$density))
	
	par(mfcol = c(1,2))
	hist(sigma.int, breaks=200, probability = TRUE, ylim=lims1)
	points(dens.int, type='l',col='red', lwd=2.0)
	
	hist(sigma.ex, breaks=200, probability = TRUE, ylim=lims2)
	points(dens.ex, type='l',col='red', lwd=2.0)
## plot genes with problems in the boundaries
	
	
### Test whether the modified model works
	
	j = which(Tt$gene=='Nup155')
	
	ii = which((Tt$BIC.best.model.Laura==2|Tt$BIC.best.model.Laura==4) & (Tt$up.time.int.m2.old==22|Tt$up.time.int.m2.old==2|Tt$down.time.int.m2.old==22|Tt$down.time.int.m2.old==2))
	
	j = mm[1]
	source('model_modify_new/functions.R')
	T = Tt
	param.fits.results = make.fits.with.all.models.for.one.gene(T = T, gene.index = j, debug = FALSE, parametrization = 'sigmoid', method = 'integration', sum.species = FALSE, zt = zt, i.ex = ZT.ex, i.int = ZT.int, absolute.signal = FALSE);
	
#kk = which(Tt$BIC.best.model==2 & (abs(Tt$gamma.m2-gamma.max)<0.00001|Tt$up.time.int.m2==12|Tt$up.time.int.m2==6))
#kk = which(Tt$BIC.best.model==2 & (abs(Tt$gamma.m2-gamma.max)<0.00001|Tt$up.time.int.m2==6))
	kk = which(Tt$BIC.best.model==2 & (Tt$up.time.int.m2==12|Tt$up.time.int.m2==6))
	ex.int = match(c("gene", "exon.var", "intron.var", "rel.ampl.ex", "phase.ex", "rel.ampl.int", "phase.int","Qval.ex", "Qval.int",
					 "error.m1", "error.m2", "gamma.m2", "fold.change.int.m2","phase.int.m2", "up.time.int.m2", 
					 "error.m4", "gamma.m4", "eps.gamma.m4","phase.gamma.m4","fold.change.int.m4", "phase.int.m4","up.time.int.m4",       
					 "BIC.m1", "BIC.m2", "BIC.m3","BIC.m4", "BIC.best.model"), colnames(T))
	
	check = T[kk, ex.int]
	
	
## check the variances of introns which have optimization problems
	index.outliers = function(data.xx)
	{
#data.xx = c(2, 3, 6, 9, 13, 18, 21, 106)
		Q1 = quantile(data.xx, 0.25,type=5)
		Q3 = quantile(data.xx, 0.75, type=5)
		IQD = Q3 - Q1
		lower = Q1 - 1.5*IQD
		upper = Q3 + 1.5*IQD
		index = which(data.xx<lower|data.xx>upper)
#boxplot(data.xx);abline(h=Q1);abline(h=Q3);
	}
	
	ZT.ex = c(2:25)
	ZT.int = c(26:49)
	
	data.ex = log(as.matrix(Tt[, ZT.ex]))
	data.int = log(as.matrix(Tt[, ZT.int]))
	sd.ex = apply(data.ex, 1, sd)
	sd.int = apply(data.int, 1, sd)
	
	Tt$exon.var = sd.ex
	Tt$intron.var = sd.int
	
	cutoff = c(10:80)/100
	sigma.s = c()
	sigma.m = c()
	sigma.ms = c()
	
	outliers = c()
	
	for(pval.cutoff in cutoff)
	{
		
		cat('Pval is ', pval.cutoff, '\n');
		data.ex = log(as.matrix(Tt[,ZT.ex]))[which(Tt$pval.ex>pval.cutoff),]
		data.int = log(as.matrix(Tt[,ZT.int]))[which(Tt$pval.int>pval.cutoff),]
		print(length(which(Tt$pval.ex>pval.cutoff)));
		print(length(which(Tt$pval.int>pval.cutoff)))
		
		data.ex = as.vector(data.ex)
		data.int = as.vector(data.int)
		data.ex.int = c(data.ex, data.int)
		
		new.data.ex = data.ex[-index.outliers(data.ex)]
		new.data.int = data.int[-index.outliers(data.int)]
		outliers = 
		new.data.ex.int = data.ex.int[-index.outliers(data.ex.int)]
		
# not remove outliers with sd
		cat('sigma estimated by sd with outliers...');
		print(c(sd(data.int), sd(data.ex), sd(data.ex.int)));
		
		cat('sigma estimated with median and outliers...');
		
		print(c(sqrt(1/(length(data.int)-1)*median(data.int^2)*length(data.int)), sqrt(1/(length(data.ex)-1)*median(data.ex^2)*length(data.ex)), sqrt(1/(length(data.ex.int)-1)*median(data.ex.int^2)*length(data.ex.int))));
		
		
# remove outliers with sd
		cat('sigma estimated by sd without outliers...means is ');
		print(c(mean(new.data.int), mean(new.data.ex), mean(new.data.ex.int)));
		cat('sigma estimated by sd without outliers...');
		print(c(sd(new.data.int), sd(new.data.ex), sd(new.data.ex.int)));
		
		new.sigma.s = sqrt(1/(length(new.data.int)-1)*mean(new.data.int^2)*length(new.data.int))
		new.sigma.m = sqrt(1/(length(new.data.ex)-1)*mean(new.data.ex^2)*length(new.data.ex))
		new.sigma.ms = sqrt(1/(length(new.data.ex.int)-1)*mean(new.data.ex.int^2)*length(new.data.ex.int))
		
# remove outliers with mean to estimate sigma
		cat('sigma estimated by mean without outliers...');
		print(c(new.sigma.s, new.sigma.m, new.sigma.ms));
		
		sigma.s = c(sigma.s, new.sigma.s)
		sigma.m = c(sigma.m, new.sigma.m)
		sigma.ms = c(sigma.ms, new.sigma.ms)
		
	}
	
	
	length(which(Tt$BIC.best.model==2 & Tt$Qval.int<0.1))
	length(which(Tt$BIC.best.model==2 & Tt$Qval.int<0.1 & (Tt$up.time.int.m2==22|Tt$up.time.int.m2==2|Tt$down.time.int.m2==22|Tt$down.time.int.m2==2)))
	length(which(Tt$BIC.best.model==2 & Tt$Qval.int<0.1 & (abs(Tt$gamma.m2-gamma.max)<0.00001|abs(Tt$gamma.m2-gamma.min)<0.00001|Tt$up.time.int.m2==22|Tt$up.time.int.m2==2|Tt$down.time.int.m2==22|Tt$down.time.int.m2==2)))
	
	
	length(which(Tt$BIC.best.model==4))
	length(which(Tt$BIC.best.model==4 & (Tt$up.time.int.m4==22|Tt$up.time.int.m4==2|Tt$down.time.int.m4==22|Tt$down.time.int.m4==2)))
	
	xx = cbind(ii, unlist(T$gamma.m2)[ii], T$fold.change.int.m2[ii], unlist(T$phase.int.m2)[ii], unlist(T$up.time.int.m2)[ii], unlist(T$down.time.int.m2)[ii], 
			   unlist(T$BIC.m1)[ii], unlist(T$BIC.m2)[ii], unlist(T$BIC.m3)[ii], unlist(T$BIC.m4)[ii], unlist(T$Qval.ex)[ii], unlist(T$Qval.int)[ii])
	colnames(xx) = c('index.gene','gamma.m2', 'fold.change.int.m2', 'phase.int.m2', 'up.time.int.m2', 'down.time.int.m2', 'BIC.m1', 'BIC.m2', 'BIC.m3', 'BIC.m4', 'Qval.ex', 'Qval.int')
	
	ii = which(Tt$BIC.best.model==3)
	xx = cbind(unlist(T$gamma.m3)[ii], unlist(T$eps.gamma.m3)[ii], unlist(T$phase.gamma.m3)[ii],
			   unlist(T$BIC.m1)[ii], unlist(T$BIC.m2)[ii], unlist(T$BIC.m3)[ii], unlist(T$BIC.m4)[ii], unlist(T$Qval.ex)[ii], unlist(T$Qval.int)[ii])
	
	ii = which(Tt$BIC.best.model==4)
	xx = cbind(unlist(T$gamma.m4)[ii],unlist(T$eps.gamma.m4)[ii], unlist(T$phase.gamma.m4)[ii], T$fold.change.int.m4[ii], unlist(T$phase.int.m4)[ii], unlist(T$up.time.int.m4)[ii], unlist(T$down.time.int.m4)[ii])
	
	
	### check if model 3 really has underdeterminations
	source('model_modify_new/functions.R')
	
	zt = seq(0,46,by = 2); 
	zt.p = seq(0,48,by = 0.5)
	eps.gamma = 5; phase.gamma = 12;
	fold.change = 1.0;phase.int=0;up.time=12;down.time=6;
	m.th = c()
	m.noise = c()
	gamma.interval = lseq(log(2)/1, log(2)/12, length=10)
	for(gamma in gamma.interval)
	{
		#gamma = log(2)/(n/2);
		m1 = compute.m.sigmoid(t = zt.p, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)
		m1 = log(m1)
		m.th = cbind(m.th, m1)
		m1 = m1 + rnorm(length(m1), mean=0, sd=0.176)
		#gamma = log(2)/3;
		#m2 = compute.m.sigmoid(t = zt.p, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)
		m.noise = cbind(m.noise, m1)
	}
	
	par(mfcol = c(1,2))
	matplot(zt.p, m.th, type='l')
	abline(h=0,col='darkgray')
	matplot(zt.p, m.noise, type='l')
	abline(h=0,col='darkgray')
	
}

###
### Analyse Results of Real Data
###
my.do.step = TRUE
load = TRUE

if(my.do.step)
{
	cat('ANALYSIS OF RESULTS FOR REAL DATA \n')
	load(paste('my_genes_fits_results', version,'.Rdata', sep = ''))
	source('model_modify_new/functions.R')
	
	T = Tt
	
	kk = which(!is.na(Tt$BIC.best.model)==TRUE)
	Tt = Tt[kk,]
	
	## main results in one plot
	selecting.methods = c('BIC.best.model', 'AIC.best.model', 'AICc.best.model', 'LRT.best.model', 'LRT.best.model.FDR')
	
	percent.rhy.degr = c()
	percent.rhy.mrna = c()
	for(method in selecting.methods)
	{
		#Percentage of rhythmic degradation for rhythmic mRNA
		best.model = eval(parse(text = paste('best.model = Tt$',method,sep ='')))
		
		percent.rhy.degr = c(percent.rhy.degr, length(which(best.model==3|best.model==4))/length(which(best.model!=1)))
	
		#Percentage of rhythmic mRNAs
		percent.rhy.mrna = c(percent.rhy.mrna,length(which(best.model==2|best.model==3|best.model==4))/nrow(Tt))
	}
	
	legend0 = c('BIC', 'AIC', 'AICc', 'LRT', 'LRT.FDR')
		
	for(n in 1:length(legend0))
	{
		legend0[n] = paste(legend0[n], '-- M34:', round(percent.rhy.degr[n], d=3)*100, '%, M234:', round(percent.rhy.mrna[n], d=3)*100, '%', sep='')
	}
	
	h1 = hist(as.numeric(Tt$BIC.best.model), breaks = seq(0.5,4.5,by = 1), plot = FALSE)
	h2 = hist(as.numeric(Tt$AIC.best.model), breaks = seq(0.5,4.5,by = 1), plot = FALSE)
	h3 = hist(as.numeric(Tt$AICc.best.model), breaks = seq(0.5,4.5,by = 1), plot = FALSE)
	h4 = hist(as.numeric(Tt$LRT.best.model), breaks = seq(0.5,4.5,by = 1), plot = FALSE)
	h5 = hist(as.numeric(Tt$LRT.best.model.FDR), breaks = seq(0.5,4.5,by = 1), plot = FALSE)
	
	#par(mfrow = c(1,3))
	#rainbow = rainbow(12,s = 0.85, v = 0.85)
	col1 = c('steelblue2','yellowgreen','indianred2','green','darkgray')
	lims = range(0,12000)
	
	pdf.version = paste('myplots/Classification_by_BIC_AIC_LRT_LRT_FDR_Real_Data', version, '.pdf', sep='')
	pdf(pdf.version, width = 5, height = 5)
	H = cbind(h1$counts, h2$counts, h3$counts, h4$counts, h5$counts)
	barplot(t(H), names.arg = paste('Model',1:4), ylim=lims, beside = TRUE, col = col1, border = NA)
	legend('topright', legend = legend0 , fill = col1, border = NA, bty = 'n', cex=0.8)
	abline(h=500,col='gray',lwd =1.5)
	abline(h=1000,col='gray',lwd =1.5)
	abline(h=2000,col='gray',lwd =1.5)
	dev.off()
	
	#### Determine the percentage of rhythmic mRNAa being reuglated by rhythmic degradation
	
	load('Sigma_my_genes_real_data.Rdata')
	selecting.methods = c('BIC.best.model', 'AIC.best.model', 'AICc.best.model', 'LRT.best.model', 'LRT.best.model.FDR')
	for(method in selecting.methods)
	{
		best.model = eval(parse(text = paste('best.model = Tt$',method,sep ='')))
	
		xx = cbind(best.model, Tt$rel.ampl.ex, Tt$pval.ex, Tt$Qval.ex)
		colnames(xx) = c(method, 'rel.ampl.ex', 'pval.ex', 'Qval.ex')
	
		ssigma = sqrt((exp(sigma.m.real^2)-1)*exp(sigma.m.real^2)) 
	
		res = c()
		cutoff.fdr = c(1, 0.1, 0.05)
		cutoff.rel.ampl = c(0, 0.5*ssigma, ssigma)
		for(cut1 in cutoff.fdr)
		{
			if(cut1==1)
			{
				cut2 = 0; 
				rhy.mrna = length(which(xx[,4]<cut1 & xx[,2]>cut2 & xx[,1]!=1))
				m2 = length(which(xx[,4]<cut1 & xx[,2]>cut2 & xx[,1]==2))
				m3 = length(which(xx[,4]<cut1 & xx[,2]>cut2 & xx[,1]==3))
				m4 = length(which(xx[,4]<cut1 & xx[,2]>cut2 & xx[,1]==4))
				res = rbind(res, c(rhy.mrna/nrow(xx), m2/rhy.mrna, m3/rhy.mrna, m4/rhy.mrna, (m3+m4)/rhy.mrna))
			}
			if(cut1<1)
			{
				for(cut2 in cutoff.rel.ampl)
				{
					rhy.mrna = length(which(xx[,4]<cut1 & xx[,2]>cut2 & xx[,1]!=1))
					m2 = length(which(xx[,4]<cut1 & xx[,2]>cut2 & xx[,1]==2))
					m3 = length(which(xx[,4]<cut1 & xx[,2]>cut2 & xx[,1]==3))
					m4 = length(which(xx[,4]<cut1 & xx[,2]>cut2 & xx[,1]==4))
					res = rbind(res, c(rhy.mrna/nrow(xx), m2/rhy.mrna, m3/rhy.mrna, m4/rhy.mrna, (m3+m4)/rhy.mrna))
				}
			}
		}
		res = data.frame(res, stringsAsFactors = FALSE)
		colnames(res) = c('percent.rhy.mRNA', 'M2', 'M3', 'M4', 'M34')
	
		cexx = 1.5
		#ltyy = c(1, rep(2, 14))
		methods.legend = c('No selection','FDR<0.1', 'FDR<0.1 & rel.ampl<0.5*sd', 'FDR<0.1 & rel.ampl<sd', 
					   'FDR<0.05', 'FDR<0.05 & rel.ampl<0.5*sd', 'FDR<0.05 & rel.ampl<sd')
		pdfname = paste('myplots/Main_results_', method, '.pdf', sep='')
		pdf(pdfname, width = 5, height = 5)
		rainbow = rainbow(nrow(res),s = 0.85, v = 0.85)
		i.rel = c(1:nrow(res))
		ylim = range(res)
		matplot(c(1:5), t(res), type='b', lwd=1.0, main=method, cex=cexx, pch = c(1:nrow(res)), col = rainbow[(i.rel-1)%%nrow(res)+1], lty = 1, xlab='Classification', ylab='Percentages',axes = FALSE)
		for(ii in i.rel)legend(x = 0.65*5,y = ylim[2]-(ylim[2]-ylim[1])*(ii-1)/nrow(res), legend = methods.legend[ii], col =  rainbow[(ii-1)%%nrow(res)+1], pch = ii, lty = 1, bty = 'n',cex=0.6)
		axis(1, at = c(1:5), labels = c('Rhythmic.mRNA','M2','M3','M4','M34'), lwd = 1)
		axis(2)
		box()
		abline(h=0.2, col='darkgray',lwd=1.5)
		abline(h=0.5, col='darkgray',lwd=1.5)
	
		dev.off()
	}
	
	###
	### check each model inferred by BIC
	###	
	xx = cbind(Tt$BIC.best.model, Tt$rel.ampl.ex, Tt$pval.ex, Tt$Qval.ex, Tt$rel.ampl.int, Tt$pval.int, Tt$Qval.int)
	colnames(xx) = c('BIC.Model', 'rel.ampl.ex', 'pval.ex', 'Qval.ex', 'rel.ampl.int', 'pval.int', 'Qval.int')
	ssigma.m = sqrt((exp(sigma.m.real^2)-1)*exp(sigma.m.real^2))
	ssigma.s = sqrt((exp(sigma.s.real^2)-1)*exp(sigma.s.real^2))
	
	pdf("myplots/Rhythimicity_Rel_ampl_mrna_premrna_models.pdf", width = 8, height = 12)
	yy = xx
	par(mfrow = c(3,2))
	
	yy[,3] = -log10(as.numeric(xx[,3]))
	boxplot(as.numeric(yy[,3])~as.numeric(yy[,1]),main="mRNA Rhythimicity", ylab="-log10 pval of mRNA", xlab="Model",ylim=c(min(as.numeric(yy[,3])),max(as.numeric(yy[,3]))))
	abline(h=-log10(0.05),lwd=1.5, col='red')
	
	yy[,6] = -log10(as.numeric(xx[,6]))
	boxplot(as.numeric(yy[,6])~as.numeric(yy[,1]),main="Pre-mRNA Rhythimicity", ylab="-log10 pval of premRNA", xlab="Model",ylim=c(min(as.numeric(yy[,6])),max(as.numeric(yy[,6]))))
	abline(h=-log10(0.05),lwd=1.5, col='red')
	
	boxplot(as.numeric(yy[,2])~as.numeric(yy[,1]), main="mRNA Rel.Ampl", ylab="rel.ampl of mRNA", xlab="Model",ylim=c(min(as.numeric(yy[,2])),max(as.numeric(yy[,2]))))
	abline(h=ssigma.m,lwd=1.5, col='red')
	abline(h=ssigma.m/2,lwd=1.5, col='red')
	
	boxplot(as.numeric(yy[,5])~as.numeric(yy[,1]),main="Pre-mRNA Rel.Ampl", ylab="rel.ampl of premRNA", xlab="Model",ylim=c(min(as.numeric(yy[,5])),max(as.numeric(yy[,5]))))
	abline(h=ssigma.s,lwd=1.5, col='red')
	abline(h=ssigma.s/2,lwd=1.5, col='red')
	
	ii.m1 = which(yy[,1]==1 & xx[,3]<0.05)
	ii.m2 = which(yy[,1]==2 & xx[,3]>0.05)
	ii.m3 = which(yy[,1]==3 & xx[,3]>0.05)
	ii.m4 = which(yy[,1]==4 & xx[,3]>0.05)
	ii = c(ii.m1, ii.m2, ii.m3, ii.m4)
	print(c(length(which(yy[,1]==1)), length(which(yy[,1]==2)), length(which(yy[,1]==3)), length(which(yy[,1]==4))))
	print(c(length(ii.m1),length(ii.m2),length(ii.m3),length(ii.m4)))
	
	boxplot(as.numeric(yy[ii,2])~as.numeric(yy[ii,1]), main="Anormal mRNA Rel.Ampl", ylab="rel.ampl of mRNA", xlab="Model",ylim=c(min(as.numeric(yy[,2])),max(as.numeric(yy[,2]))))
	abline(h=ssigma.m,lwd=1.5, col='red')
	abline(h=ssigma.m/2,lwd=1.5, col='red')
	
	ii.m1 = which(yy[,1]==1 & xx[,6]<0.05)
	ii.m2 = which(yy[,1]==2 & xx[,6]>0.05)
	ii.m3 = which(yy[,1]==3 & xx[,6]<0.05)
	ii.m4 = which(yy[,1]==4 & xx[,6]>0.05)
	ii = c(ii.m1, ii.m2, ii.m3, ii.m4)
	print(c(length(which(yy[,1]==1)), length(which(yy[,1]==2)), length(which(yy[,1]==3)), length(which(yy[,1]==4))))
	print(c(length(ii.m1),length(ii.m2),length(ii.m3),length(ii.m4)))
	
	boxplot(as.numeric(yy[ii,5])~as.numeric(yy[ii,1]),main="Anormal Pre-mRNA Rel.Ampl", ylab="rel.ampl of premRNA", xlab="Model",ylim=c(min(as.numeric(yy[,5])),max(as.numeric(yy[,5]))))
	abline(h=ssigma.s,lwd=1.5, col='red')
	abline(h=ssigma.s/2,lwd=1.5, col='red')
	
	dev.off()
		
	ii = which(xx[,1]==3 & xx[,3]<0.05)
	jj = which(xx[,1]==3 & xx[,3]>0.05)
	
	lims = range(xx[c(ii,jj),2])
	par(mfrow = c(1,2))
	hist(as.numeric(xx[ii,2]),breaks=20,xlim=lims,main='BIC.M1 exons for pval<0.05',xlab='rel.amp');abline(v=ssigma.m/2,col='red');abline(v=ssigma.m,col='red')
	hist(as.numeric(xx[jj,5]),breaks=20, xlim=lims,main='BIC.M1 rel.ampl of introns for pval<0.05',xlab='rel.amp');abline(v=ssigma.s/2,col='red');abline(v=ssigma.s,col='red')
	
		
	ii = which(xx[,1]==2)
	hist(as.numeric(xx[ii,2]),breaks=20,xlim=lims,main='BIC.M2 rel.ampl of exons for FDR<0.05');abline(v=ssigma.m/2,col='red');abline(v=ssigma.m,col='red')
	ii = which(xx[,1]==2 & xx[,3]<0.05)
	jj = which(xx[,1]==2 & xx[,6]<0.05)
	lims = range(xx[,2])
	par(mfrow = c(1,2))
	hist(as.numeric(xx[ii,2]),breaks=20,xlim=lims,main='BIC.M2 rel.ampl of exons for FDR<0.05');abline(v=ssigma.m/2,col='red');abline(v=ssigma.m,col='red')
	hist(as.numeric(xx[jj,5]),breaks=20, xlim=lims,main='BIC.M2 rel.ampl of introns for FDR<0.05');abline(v=ssigma.s/2,col='red');abline(v=ssigma.s,col='red')
	
	
	
	o1 = order(as.numeric(yy[,3]))
	ii = ii[o1]
	yy = yy[o1,]
	
	
	pdf("myplots/BIC_model_1_but_rhythmic.pdf", width = 5, height = 5)
	ZT.ex = grep('.rel.ampl.ex', colnames(Tt))
	zt = seq(0,46,by = 2)
	for(kk in ii)
	{
		test = as.numeric(Tt[kk, ZT.ex])
		plot(zt, test, col='darkblue', type='b', lwd=0.7, main=paste(Tt[kk,1], ', rel.ampl=', round(Tt$rel.ampl.ex[kk],d=2), ', qv=', round(Tt$Qval.ex[kk],d=5), sep=''))
		abline(h=1, lwd=2.0, col='gray')
		
	}
	dev.off()
	
	cutoff = 0.1
	ii = which(Tt$Qval.ex<cutoff & Tt$rel.ampl.ex>sigma.m.real)
	length(ii)
	
	xx = xx[ii,]
	xx[,2] = -log10(as.numeric(xx[,2]))
	
	boxplot(as.numeric(xx[,2])~as.numeric(xx[,1]), ylab="log10 pval", xlab="Model",ylim=c(min(as.numeric(xx[,2])),max(as.numeric(xx[,2]))))
	abline(h=-log10(0.05),lwd=1.5, col='red')
	abline(h=-log10(0.01),lwd=1.5, col='red')
	
	length(which(xx[,1]==1))
	length(which(xx[,1]==2))
	length(which(xx[,1]==3))
	length(which(xx[,1]==4))
	length(ii)/nrow(Tt)
	(length(which(xx[,1]==2)))/(length(which(xx[,1]==2))+length(which(xx[,1]==3))+ length(which(xx[,1]==4)))
	(length(which(xx[,1]==3)))/(length(which(xx[,1]==2))+length(which(xx[,1]==3))+ length(which(xx[,1]==4)))
	(length(which(xx[,1]==4)))/(length(which(xx[,1]==2))+length(which(xx[,1]==3))+ length(which(xx[,1]==4)))
	(length(which(xx[,1]==3))+ length(which(xx[,1]==4)))/(length(which(xx[,1]==2))+length(which(xx[,1]==3))+ length(which(xx[,1]==4)))
	
	## select a special subset of genes in which premrna is rhythmic but mrna is flat
	yy = cbind(Tt$BIC.best.model, Tt$Qval.ex, Tt$rel.ampl.ex, Tt$Qval.int, Tt$rel.ampl.int)
	colnames(yy) = c('Model', 'Qval.ex', 'rel.ampl.ex', 'Qval.int', 'rel.ampl.int')
	
	cutoff.ex = 0.5
	cutoff.int = 0.1
	ii = which(Tt$BIC.best.model==4 & Tt$Qval.int<cutoff.int & Tt$rel.ampl.int>sigma.s.real & Tt$Qval.ex>cutoff.ex)
	#ii = which(Tt$Qval.int<cutoff.int & Tt$rel.ampl.int>sigma.s.real)
	length(ii)
	
	par(mfrow = c(1,3))
	boxplot(Tt$LRT.quality ~ Tt$LRT.best.model, col = c(1:4), pch = 16, cex = 0.6); abline(h= 1)
	h = hist(Tt$LRT.best.model, breaks = seq(0.5,4.5,by=1), col = 1:4, border = 'NA')
	barplot(h$counts[2:4],col=2:4)
	
	CIRC = c("Per1","Per2","Per3","Cry1","Cry2","Dbp","Tef","Hlf","Nr1d1","Nr1d2","Rora","Rorb","Rorc","Arntl","Bmal1","Clock","Npas2","Bhlhe40","Bhlhe41","Cirbp","Hamp","Hamp2","Nr4a2","Tfrc","Wee1", "Por", 
			 "Abcb11","Tfrc", "Per2","Nedd4l", "Cbs", "Gsta3", "Loxl4","Rcan1", "Cdkn1a", "Tubb2a", "Mfsd2", "Ppard")
	k = match(CIRC,Tt$gene)
	cbind(Tt$gene[k], Tt$BIC.best.model[k], Tt$AIC.best.model[k], Tt$AICc.best.model[k], Tt$LRT.best.model[k])

	## compare the nb of models selected by LRT, AIC and AICc
	par(mfrow = c(2,3))
		
	plot(as.numeric(Tt$LRT.best.model)+rnorm(nrow(Tt),sd = 0.1),Tt$BIC.best.model+rnorm(nrow(Tt),sd = 0.1),pch = 16, cex = 1, col = rgb(0,0.5,0.8,0.3))
	plot(as.numeric(Tt$AICc.best.model)+rnorm(nrow(Tt),sd = 0.1),Tt$BIC.best.model+rnorm(nrow(Tt),sd = 0.1),pch = 16, cex = 1, col = rgb(0.3,0.8,0,0.3))
	plot(as.numeric(Tt$AIC.best.model)+rnorm(nrow(Tt),sd = 0.1),Tt$BIC.best.model+rnorm(nrow(Tt),sd = 0.1),pch = 16, cex = 1, col = rgb(0.3,0,0.8,0.3))
	plot(as.numeric(Tt$LRT.best.model)+rnorm(nrow(Tt),sd = 0.1),Tt$BIC.best.model+rnorm(nrow(Tt),sd = 0.1),pch = 16, cex = 1, col = rgb(0,0.5,0,0.3))
	plot(as.numeric(Tt$AICc.best.model)+rnorm(nrow(Tt),sd = 0.1),Tt$LRT.best.model+rnorm(nrow(Tt),sd = 0.1),pch = 16, cex = 1, col = rgb(0.8,0.3,0,0.3))

}


########
######## Analysis of phases, amplitudes, half-lifes for different classes.
########

#Precise exon time-course for results used BIC

do.step = TRUE
if(do.step)
{
	cat('PRECISE EXON TIME-COURSE BIC\n')
	load(paste('my_genes_fits_results', version,'.Rdata', sep = ''))
	source('model_modify_new/functions.R')
	
	kk = which(!is.na(Tt$BIC.best.model)==TRUE)
	Tt = Tt[kk,]
	T = Tt
	
	tp = seq(0,23.9, by = 0.1)	
	
	P = data.frame(matrix(1,nrow = nrow(T), ncol = length(tp)), T$gene, stringsAsFactors = FALSE)	
	ind = c(1:length(tp))
	
	best.model = T$BIC.best.model
	for(j in which(best.model  == 2))
	{
		ww = which(best.model == 2)
		cat(which(ww==j)/length(ww),'\n')
		ex = compute.m.sigmoid(t = tp, gamma = T$gamma.m2[j], eps.gamma = 1, phase.gamma = 0, fold.change = T$fold.change.int.m2[j], phase = T$phase.int.m2[j], up.time = T$up.time.int.m2[j], down.time = T$up.time.int.m2[j])
		P[j,ind] = ex
	}	
	for(j in which(best.model == 3))
	{
		ww = which(best.model  == 3);
		if((which(ww==j)%%50)==0){cat(round(100*which(ww==j)/length(ww), digits = 1),'%\n')}
		ex = compute.m.sigmoid(t = tp, gamma = T$gamma.m3[j], eps.gamma = T$eps.gamma.m3[j], phase.gamma = T$phase.gamma.m3[j], fold.change = 1, phase = 0, up.time = 12, down.time = 12)
		P[j,ind] = ex
	}	
	for(j in which( best.model  == 4))
	{
		ww = which(best.model  == 4);
		i = which(ww==j)
		if((i%%50)==0){cat(round(100*which(ww==j)/length(ww), digits = 1),'%\n')}
		ex = compute.m.sigmoid(t = tp, gamma = T$gamma.m4[j], eps.gamma = T$eps.gamma.m4[j], phase.gamma = T$phase.gamma.m4[j], fold.change = T$fold.change.int.m4[j], phase = T$phase.int.m4[j], up.time = T$up.time.int.m4[j], down.time = T$up.time.int.m4[j])
		P[j,ind] = ex
	}	
	
	save(P,file = paste('my_presice_exon_time_course_BIC_fits_results',version,'.Rdata', sep = ''))
}

## summary the half-lives, amplitude, phase of different classes in on plot
my.plot = TRUE
load.for.plot = TRUE

if(my.plot)
{
	version.1 = '_my_8th_same_up_down_time_all_genes_model3_no_treated'
	version.2 = '_my_8th_same_up_down_time_all_genes_sum'

	version = version.2
	cat('plot : FITS RESULTS ON REAL DATA -- FIGURE 5\n')
	
	load(paste('my_genes_fits_results', version,'.Rdata', sep = ''))
	source('model_modify_new/functions.R')
	
	kk = which(!is.na(Tt$LRT.best.model)==TRUE)
	Tt = Tt[kk,]
	T = Tt
	
	best.model = as.numeric(Tt$BIC.best.model)
	
	### Extract parameters inferred
	TEST = FALSE
	if(TEST)
	{
		
		ex.int = match(c("gene", 
						 "gamma.m2", "gamma.stderr.m2",  "up.time.int.m2", "up.time.int.stderr.m2",  
						 "gamma.m3","gamma.stderr.m3", "eps.gamma.m3" ,"eps.gamma.stderr.m3",
						 "gamma.m4","gamma.stderr.m4", "eps.gamma.m4" ,"eps.gamma.stderr.m4","up.time.int.m4", "up.time.int.stderr.m4", "BIC.best.model"), colnames(Tt))
		
		jj = which(T$BIC.best.model==4)
		test = Tt[jj,ex.int]
		coeffs <- test$gamma.m4
		stderr <- test$gamma.stderr.m4
		zscore <- coeffs/stderr
		length(jj)
		length(which(zscore<1.036))/length(jj)
		
	}
	
	## up time of pre-mrna for model 2 and 4
	up.time.int = rep(1, nrow(Tt));
	up.time.int.zscore = rep(1, nrow(Tt));
	up.time.int[which(best.model == 2)] = Tt$up.time.int.m2[which(best.model == 2)]
	up.time.int[which(best.model == 4)] = Tt$up.time.int.m4[which(best.model == 4)]
	up.time.int.zscore[which(best.model == 2)] = Tt$up.time.int.m2[which(best.model == 2)]/Tt$up.time.int.stderr.m2[which(best.model==2)]
	up.time.int.zscore[which(best.model == 4)] = Tt$up.time.int.m4[which(best.model == 4)]/Tt$up.time.int.stderr.m4[which(best.model==4)]
	
	## phase of pre-mrna for model 2 and model 4
	phase.int = rep(1, nrow(Tt));
	phase.int.zscore = rep(1, nrow(Tt));
	phase.int[which(best.model == 2)] = Tt$phase.int.m2[which(best.model == 2)]
	phase.int[which(best.model == 4)] = Tt$phase.int.m4[which(best.model == 4)]
	phase.int.zscore[which(best.model == 2)] = Tt$phase.int.m2[which(best.model == 2)]/Tt$phase.int.stderr.m2[which(best.model==2)]
	phase.int.zscore[which(best.model == 4)] = Tt$phase.int.m4[which(best.model == 4)]/Tt$phase.int.stderr.m4[which(best.model==4)]
	
	## gamma for model 2, 3 , 4
	gamma = rep(1, nrow(Tt));
	gamma.zscore = rep(1, nrow(Tt));
	gamma[which(best.model == 2)] = Tt$gamma.m2[which(best.model == 2)]
	gamma[which(best.model == 4)] = Tt$gamma.m4[which(best.model == 4)]
	gamma[which(best.model == 3)] = Tt$gamma.m3[which(best.model == 3)]
	
	gamma.zscore[which(best.model == 2)] = Tt$gamma.m2[which(best.model == 2)]/Tt$gamma.stderr.m2[which(best.model==2)]
	gamma.zscore[which(best.model == 4)] = Tt$gamma.m4[which(best.model == 4)]/Tt$gamma.stderr.m4[which(best.model==4)]
	gamma.zscore[which(best.model == 3)] = Tt$gamma.m3[which(best.model == 3)]/Tt$gamma.stderr.m3[which(best.model==3)]
	
	## fold change of gamma in model 3 and 4
	gamma.fc = rep(1, nrow(Tt));
	gamma.fc.zscore = rep(1, nrow(Tt));
	gamma.fc[which(best.model == 3)] = Tt$eps.gamma.m3[which(best.model == 3)]
	gamma.fc[which(best.model == 4)] = Tt$eps.gamma.m4[which(best.model == 4)]
	gamma.fc.zscore[which(best.model == 3)] = Tt$eps.gamma.m3[which(best.model == 3)]/Tt$eps.gamma.stderr.m3[which(best.model==3)]
	gamma.fc.zscore[which(best.model == 4)] = Tt$eps.gamma.m4[which(best.model == 4)]/Tt$eps.gamma.stderr.m4[which(best.model==4)]
	
	## phase of gamma for model 3 and 4 
	phase.gamma = rep(1, nrow(Tt));
	phase.gamma.zscore = rep(1, nrow(Tt));
	phase.gamma[which(best.model == 3)] = Tt$phase.gamma.m3[which(best.model == 3)]
	phase.gamma[which(best.model == 4)] = Tt$phase.gamma.m4[which(best.model == 4)]
	phase.gamma.zscore[which(best.model == 3)] = Tt$phase.gamma.m3[which(best.model == 3)]/Tt$phase.gamma.stderr.m3[which(best.model==3)]
	phase.gamma.zscore[which(best.model == 4)] = Tt$phase.gamma.m4[which(best.model == 4)]/Tt$phase.gamma.stderr.m4[which(best.model==4)]
	
	## Figure for the main results 
	#version.old = '_my_8th_same_up_down_time_all_genes_model3_no_treated'
	load(paste('my_presice_exon_time_course_BIC_fits_results',version,'.Rdata', sep = ''))
	
	source('model_modify_new/functions.R')
	
	pdf(paste('myplots/fits_results_BIC',version,'.pdf',sep = ''), width = 6.6, height = 4)
	plot.figure.fit.results(T =T, best.model = best.model , gamma = gamma, gamma.fc = gamma.fc, phase.gamma = phase.gamma, phase.int = phase.int, up.time.int = up.time.int, P = P)
	dev.off()
	
	source('model_modify_new/functions.R')
	#pdf(paste('myplots/fits_results_BIC',version,'_Zscore.pdf',sep = ''), width = 6.6, height = 4)
	plot.figure.fit.results.zscore(T =T, best.model = best.model, gamma = gamma, gamma.zscore = gamma.zscore, gamma.fc = gamma.fc, gamma.fc.zscore = gamma.fc.zscore, 
								   phase.gamma = phase.gamma,phase.gamma.score = phase.gamma.score, phase.int = phase.int, phase.int.score = phase.int.score, 
								   up.time.int = up.time.int, up.time.int.score = up.time.int.score, P = P)
	#dev.off()
	
	circ.genes = c("Per1","Per2","Per3","Cry1","Cry2","Dbp","Tef","Hlf","Nr1d1","Nr1d2","Rora","Rorb","Rorc","Arntl","Clock","Npas2","Bhlhe40","Bhlhe41","Cirbp","Tfrc","Nr4a2","Wee1", "Por")	
	
	m = match(circ.genes, T$gene)
	CIRC = cbind(T$gene, T$model.MCMC,T$LRT.best.model, T$best.model)[m,]
	colnames(CIRC) = c('gene','MCMC','LRT','AIC')
	CIRC
	write.table(CIRC, file = paste('list_and_tables/fits_results_for_circadian_genes',version,'.txt',sep = ''), quote = FALSE, sep = '\t', row.names = FALSE)
	
	
	## save list of genes of different classes
	write.table(sort(T$gene[which(T$LRT.best.model == 1)]), file = paste('list_and_tables/model1_LRT',version,'.txt',sep = ''), quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE)
	write.table(sort(T$gene[which(T$LRT.best.model == 2)]), file = paste('list_and_tables/model2_LRT',version,'.txt',sep = ''), quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE)
	write.table(sort(T$gene[which(T$LRT.best.model == 3)]), file = paste('list_and_tables/model3_LRT',version,'.txt',sep = ''), quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE)
	write.table(sort(T$gene[which(T$LRT.best.model == 4)]), file = paste('list_and_tables/model4_LRT',version,'.txt',sep = ''), quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE)
	
	write.table(sort(T$gene[which(T$best.model == 1)]), file = paste('list_and_tables/model1_AIC',version,'.txt',sep = ''), quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE)
	write.table(sort(T$gene[which(T$best.model == 2)]), file = paste('list_and_tables/model2_AIC',version,'.txt',sep = ''), quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE)
	write.table(sort(T$gene[which(T$best.model == 3)]), file = paste('list_and_tables/model3_AIC',version,'.txt',sep = ''), quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE)
	write.table(sort(T$gene[which(T$best.model == 4)]), file = paste('list_and_tables/model4_AIC',version,'.txt',sep = ''), quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE)
	
	
	source('functions.R') ;
	ok.model2 = which(T$BIC.best.model == 2)
	#sort(T$gene[ok.model2])
	ok.model3 = which(T$BIC.best.model == 3)
	#sort(T$gene[ok.model3])
	ok.model4 = which(T$BIC.best.model == 4)
	#sort(T$gene[ok.model4])
	RefSeq = read.table('/Users/jiwang/Degradation_Laura/_x_MicroArray/complementary_data/RefSeq_Mouse_mm9.txt', stringsAsFactors = FALSE, sep = '\t', header = TRUE)
	T = Tt
	t = tt
	
	###
	### Plots of individual genes in different classes showing how to summarize the exon and intron signals
	###
	if(!file.exists('myplots/gene_examples/')){dir.create('myplots/gene_examples/')}
	
	gene_list = T$gene[ok.model2]
	pdf.folder = paste('myplots/gene_examples/model2',version,sep = '')
	if(!file.exists(pdf.folder)){dir.create(pdf.folder)}	
	for(gene in gene_list){cat(gene, '\n');plot.this.gene.full(pdf.folder = pdf.folder,gene = gene, ylim.rel.type = 'adjusted', ylim.rel = c(0,2), log2.rel = FALSE, with.fit = TRUE, fitting.method = 'Optim', parametrization = 'sigmoid')}
	
	gene_list = T$gene[ok.model3]
	pdf.folder = paste('plots/gene_examples/model3',version,sep = '')
	if(!file.exists(pdf.folder)){dir.create(pdf.folder)}	
	for(gene in gene_list){cat(gene, '\n');plot.this.gene.full(pdf.folder = pdf.folder,gene = gene, ylim.rel.type = 'adjusted', ylim.rel = c(0,2), log2.rel = FALSE, with.fit = TRUE, fitting.method = 'Optim', parametrization = 'sigmoid')}
	
	gene_list = T$gene[ok.model4]
	pdf.folder = paste('plots/gene_examples/model4',version,sep = '')
	if(!file.exists(pdf.folder)){dir.create(pdf.folder)}	
	for(gene in gene_list){cat(gene, '\n');plot.this.gene.full(pdf.folder = pdf.folder,gene = gene, ylim.rel.type = 'adjusted', ylim.rel = c(0,2), log2.rel = FALSE, with.fit = TRUE, fitting.method = 'Optim', parametrization = 'sigmoid')}
	
	
	## Main Figure 6 fitting examples
	source('functions.R')
	zt = seq(0,46,by = 2); 
	attribute.global.variable.colors()
	
	pdf(paste('myplots/figures_main/microarray_fig6_examples',version,'.pdf',sep =''), width = 6.6, height = 4)
	par(mfcol = c(2,3), las = 1, cex = 0.6, mar = c(3,3,2.2,0.5), mgp = c(1.7,0.5,0), tcl = -0.3)	
	genes = c('Abcb11', 'Tfrc','Per2'); label = c('A','B','C')
	for(i in 1:length(genes))
	{
		gene = genes[i]
		show.summarized.signal(gene = gene, ylim.rel.type = 'personalized', ylim.rel = c(0,2.5), with.fit = TRUE, fitting.method = 'Optim', parametrization = 'sigmoid', abline = TRUE, dt = 6, title = c(gene,'model'))
		mtext(label[i],side =3, line = 17.2, at = -6, font = 2)
	}
	dev.off()
	
	## supplementary examples for fitting and model selection
	pdf(paste('myplots/figures_supplement/microarray_S_fit_examples',version,'.pdf', sep = ''), width = 6.6, height = 9)
	par(mfcol = c(6,3), las = 1, cex = 0.6, mar = c(3,3,2.2,0.5), mgp = c(1.7,0.5,0), tcl = -0.3)	
	genes = c('Nedd4l', 'Cbs','Gsta3','Loxl4','Rcan1','Cdkn1a','Tubb2a','Mfsd2','Ppard'); label = LETTERS[1:9]
	for(i in 1:length(genes)){
		cat(gene,'\n')
		gene = genes[i]
		show.summarized.signal(gene = gene, ylim.rel.type = 'personalized', ylim.rel = c(0,2.5), with.fit = TRUE, fitting.method = 'Optim', parametrization = 'sigmoid', abline = TRUE, dt = 6, title = c(gene,'model'))
		mtext(label[i],side =3, line = 13.2, at = -6, font = 2)
	}
	dev.off()
	
	## supplementary examples for circadian genes
	pdf(paste('myplots/figures_supplement/microarray_S_circadian_examples',version,'.pdf',sep = ''), width = 6.6, height = 9)
	par(mfcol = c(6,4), las = 1, cex = 0.6, mar = c(3,3,2.2,0.5), mgp = c(1.7,0.5,0), tcl = -0.3)	
	genes = c('Per1','Per3', 'Cry1','Nr1d1','Nr1d2','Rorc','Arntl','Tef','Hlf','Npas2','Clock','Bhlhe40'); label = LETTERS[c(1,5,9,2,6,10,3,7,11,4,8,12)]
	for(i in 1:length(genes))
	{
		gene = genes[i]
		show.summarized.signal(gene = gene, ylim.rel.type = 'personalized', ylim.rel = c(0,2.5), with.fit = TRUE, fitting.method = 'Optim', parametrization = 'sigmoid', abline = TRUE, dt = 12, title = c(gene,'model'))
		mtext(label[i],side =3, line = 13.2, at = -6, font = 2)
	}
	dev.off()
	
	## supplementary examples for RNA binding proteins
	genes.RBP = c('Khdrbs3','Snrpa_chr7','Pabpc1','Sfrs1','Sfrs2_chr17','Sfrs2_chr15', 'Sfrs2_chr11','Sfrs9','Nova2','Elavl1','Elavl2','Ptbp1_chr10','Pum2','Ncl','Khsrp_chr17','Ybx1_chr4','Hnrnpa1_chr15','Nono', 'A2bp1')
	
	pdf(paste('myplots/figures_supplement/microarray_S_regulators_examples',version,'.pdf',sep = ''), width = 6.6, height = 9)
	par(mfcol = c(6,3), las = 1, cex = 0.6, mar = c(3,3,2.2,0.5), mgp = c(1.7,0.5,0), tcl = -0.3)	
	genes = c(T$gene[grep('Cpeb', T$gene)],'Ccrn4l_chr3', genes.RBP); label = LETTERS
	for(i in 1:length(genes)){
		gene = genes[i]
		show.summarized.signal(gene = gene, ylim.rel.type = 'personalized', ylim.rel = c(0,2.5), with.fit = TRUE, fitting.method = 'Optim', parametrization = 'sigmoid', abline = TRUE, dt = 12, title = c(gene,'model'))
		mtext(label[i],side =3, line = 13.2, at = -6, font = 2)
	}
	dev.off()
	
	
	gene_list = T$gene[ok.model3]
	
	pdf(paste('myplots/figures_supplement/summary_model3_genes/microarray_model3_examples',version,'.pdf',sep =''), width = 6.6, height = 4)
	par(mfcol = c(2,3), las = 1, cex = 0.6, mar = c(3,3,2.2,0.5), mgp = c(1.7,0.5,0), tcl = -0.3)	
	genes = gene_list
#label = c('A','B','C')
	for(i in 1:length(genes))
	{
		cat(gene, '\n');
		gene = genes[i]
		show.summarized.signal(gene = gene, ylim.rel.type = 'personalized', ylim.rel = c(0,2.5), with.fit = TRUE, fitting.method = 'Optim', parametrization = 'sigmoid', abline = TRUE, dt = 6, title = c(gene,'model'))
#mtext(label[i],side =3, line = 17.2, at = -6, font = 2)
	}
	dev.off()
	
	gene_list = T$gene[ok.model4]
	
	pdf(paste('myplots/figures_supplement/summary_model4_genes/microarray_model4_examples',version,'.pdf',sep =''), width = 6.6, height = 4)
	par(mfcol = c(2,3), las = 1, cex = 0.6, mar = c(3,3,2.2,0.5), mgp = c(1.7,0.5,0), tcl = -0.3)	
	genes = gene_list
	for(i in 1:length(genes))
	{
		cat(gene, '\n');
		gene = genes[i]
		show.summarized.signal(gene = gene, ylim.rel.type = 'personalized', ylim.rel = c(0,2.5), with.fit = TRUE, fitting.method = 'Optim', parametrization = 'sigmoid', abline = TRUE, dt = 6, title = c(gene,'model'))
	}
	dev.off()
	
	gene_list = T$gene[ok.model2]
	
	pdf(paste('myplots/figures_supplement/summary_model2_genes/microarray_model2_examples',version,'.pdf',sep =''), width = 6.6, height = 4)
	par(mfcol = c(2,3), las = 1, cex = 0.6, mar = c(3,3,2.2,0.5), mgp = c(1.7,0.5,0), tcl = -0.3)	
	genes = gene_list
	for(i in 1:length(genes))
	{
		cat(gene, '\n');
		gene = genes[i]
		show.summarized.signal(gene = gene, ylim.rel.type = 'personalized', ylim.rel = c(0,2.5), with.fit = TRUE, fitting.method = 'Optim', parametrization = 'sigmoid', abline = TRUE, dt = 6, title = c(gene,'model'))
	}
	dev.off()
	
}


########################################
########################################  Data validation by qPCR from David Gatfield
if(my.do.step)
{
	
# november2013 measures
qPCR = read.table('/Users/jiwang/Degradation_Laura/_x_MicroArray/complementary_data/qPCR_DGatfield/november2013/final qPCR results.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
qPCR = t(qPCR[,-1])
colnames(qPCR) = c(paste('ZT0',seq(0,8,by = 2), sep = ''), paste('ZT',seq(10,46,by = 2), sep = ''), 'RT')
rownames(qPCR) = substr(rownames(qPCR),6,max(nchar(rownames(qPCR))))

qPCR1 = qPCR[,1:24]
j = which(rownames(qPCR1) == '')
if(length(j)>0){qPCR1 = qPCR1[-j,]}

write.table(qPCR1,'qPCR1.txt', quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

#december2013 measures

qPCR = read.table('/Users/jiwang/Degradation_Laura/_x_MicroArray/complementary_data/qPCR_DGatfield/december2013/20131203-result table.csv', header = TRUE, sep = ',', dec = '.', stringsAsFactors = FALSE)
qPCR = t(qPCR[,-1])
colnames(qPCR) = c(paste('ZT0',seq(0,8,by = 2), sep = ''), paste('ZT',seq(10,46,by = 2), sep = ''))

rownames(qPCR) = substr(rownames(qPCR),6,max(nchar(rownames(qPCR))))

qPCR2 = qPCR

write.table(qPCR2,'qPCR2.txt', quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

#######

rownames(qPCR1) = paste(rownames(qPCR1), '.epx1', sep = '')
rownames(qPCR2) = paste(rownames(qPCR2), '.epx2', sep = '')

qPCR = rbind(qPCR1, qPCR2)
o = order(rownames(qPCR))
qPCR  = qPCR[o,]


full.gene.list  = sapply(strsplit(rownames(qPCR),split ='\\.'), "[",1)
gene.list = unique(full.gene.list)
exon = rep(FALSE, nrow(qPCR))
exon[grep('.E', rownames(qPCR))] = TRUE
intron = rep(FALSE, nrow(qPCR))
intron[grep('.I', rownames(qPCR))] = TRUE

exp = c(rep(1,nrow(qPCR1)), rep(2,nrow(qPCR1))); exp = exp[o]

qpcr = data.frame(gene = full.gene.list, exon = exon, intron = intron, exp = exp  , qPCR, stringsAsFactors = FALSE)
ZT = grep('ZT', colnames(qpcr))
qpcr$mean = apply(qpcr[,ZT],1,mean, na.rm = TRUE)
qpcr[,ZT] = qpcr[,ZT]/qpcr$mean


load(paste('my_genes_fits_results', version,'.Rdata', sep = ''))
source('model_modify_new/functions.R')

kk = which(!is.na(Tt$LRT.best.model)==TRUE)
Tt = Tt[kk,]
T = Tt
	
T.ex = grep(".rel.ampl.ex",colnames(T))
T.int = grep(".rel.ampl.int",colnames(T))
zt = seq(0,46,by = 2)
attribute.global.variable.colors()
zt.p = seq(0,48,by = 0.1)


#profile comparison
pdf('myplots/qPCR_validation_plots_profiles_only.pdf', width = 6.6, height = 6.6)

for(i in 1:length(gene.list))
	{	
	cat(i,'\n')
	
	gene = gene.list[i]
	
	j  = which(T$gene== gene)
	k.int = which((qpcr$gene == gene)& qpcr$intron )
	k.ex = which((qpcr$gene == gene)& qpcr$exon )
	
	if(length(j)>0)
	{
		
		par(mfrow = c(2,2),cex = 0.5, mar = c(3,3,2,0.2)+0.1, mgp = c(1.9,0.5,0), las = 1)
		lwd = 1.2
		lty = '31'
		
		plot(1,1, type = 'n' , xlim = c(0,46), ylim = range(0,qpcr[k.int, ZT], T[j, T.int],2, na.rm = TRUE), main = paste(gene,' - intron'), xlab = 'zt time (h)', ylab = 'rel. signal')
		for(k  in k.int){points(zt,qpcr[k, ZT], col = 'green2', type = 'l' , lwd = lwd, lty = qpcr$exp[k])}
		points(zt,T[j, T.int], col = col.int, type = 'l' , lwd = lwd, pch = 16, lty = lty)
		abline(h = 1, col = 'gray', lty = 3)
		
		
		plot(1,1, type = 'n' , xlim = c(0,46), ylim = range(0,qpcr[k.ex, ZT], T[j, T.ex],2, na.rm = TRUE), main = paste(gene,' - exon'), xlab = 'zt time (h)', ylab = 'rel. signal')
		for(k  in k.ex){ points(zt,qpcr[k, ZT], col = 'steelblue3', type = 'l' , lwd = lwd, lty = qpcr$exp[k])}
		points(zt,T[j, T.ex], col = col.ex, type = 'l' , lwd = lwd, pch = 16, lty = lty)
		abline(h = 1, col = 'gray', lty = 3)
		
		
		plot(1,1, type = 'n' , xlim = c(0,46), ylim = range(0,qpcr[k.int, ZT], qpcr[k.ex, ZT], T[j, T.ex], T[j, T.int],2, na.rm = TRUE), main = paste(gene,' - qPCR'), xlab = 'zt time (h)', ylab = 'rel. signal')
		for(k  in k.int){ points(zt,qpcr[k, ZT], col = 'green2', type = 'l' , lwd = lwd, lty = qpcr$exp[k])}
		for(k  in k.ex){ points(zt,qpcr[k, ZT], col = 'steelblue3', type = 'l' , lwd = lwd, lty = qpcr$exp[k])}
		abline(h = 1, col = 'gray', lty = 3)
		
		plot(zt,T[j, T.ex], col = col.ex, type = 'l' , lwd = lwd, pch = 16, lty = lty, ylim = range(0,qpcr[k.int, ZT], qpcr[k.ex, ZT], T[j, T.ex], T[j, T.int],2, na.rm = TRUE), main = paste(gene,' - exonarray'), xlab = 'zt time (h)', ylab = 'rel. signal')
		points(zt,T[j, T.int], col = col.int, type = 'l' , lwd = lwd, pch = 16, lty = lty)
		abline(h = 1, col = 'gray', lty = 3)
		
		cat('here\n')
	}
	
	}
	dev.off()

}

### Results validation by comparing with mrnas regulated by miRNAs from David Gatfield elife
if(my.do.step)
{
	cat('COMPARSION BETWEEN OUR RESULTS AND DAVID GATFIELD \n')
	
	load(paste('my_genes_fits_results', version,'.Rdata', sep = ''))
	source('model_modify_new/functions.R')
	
	kk = which(!is.na(Tt$LRT.best.model)==TRUE)
	Tt = Tt[kk,]
	T = Tt
	
	load("Sigma_my_genes_real_data.Rdata")
	
	selecting.methods = c('BIC.best.model', 'AIC.best.model', 'AICc.best.model', 'LRT.best.model',
						  'BIC.best.model.Old', 'AIC.best.model.Old','AICc.best.model.Old', 'LRT.best.model.Old', 
						  'BIC.best.model.Laura','AIC.best.model.Laura', 'AICc.best.model.Laura', 'LRT.best.model.Laura','LRT.best.model.Laura.FDR')
	
	## check the fractions of interest and try to understand the discrepencies between Koike, Menet and Cyclix
	xx = cbind(Tt$gene, Tt$BIC.best.model, Tt$rel.ampl.ex, Tt$pval.ex, Tt$Qval.ex, Tt$rel.ampl.int, Tt$pval.int, Tt$Qval.int)
	xx = data.frame(matrix(unlist(xx), nrow=nrow(Tt), byrow=FALSE))
	colnames(xx) = c('Gene.names','BIC.Model','rel.ampl.ex', 'pval.ex', 'Qval.ex', 'rel.ampl.int', 'pval.int', 'Qval.int')
	
	print(length(which(as.numeric(Tt$Qval.ex)<0.05)))
	print(length(which(as.numeric(Tt$Qval.ex)<0.05 & as.numeric(Tt$BIC.best.model)==1)))
	
	## compare our results with David Gatfield results
	xx = cbind(Tt$gene, Tt$BIC.best.model, Tt$AIC.best.model, Tt$AICc.best.model, Tt$LRT.best.model, 
			   Tt$BIC.best.model.Old, Tt$AIC.best.model.Old, Tt$AICc.best.model.Old, Tt$LRT.best.model.Old, 
			   Tt$BIC.best.model.Laura, Tt$AIC.best.model.Laura, Tt$AICc.best.model.Laura, Tt$LRT.best.model.Laura, Tt$LRT.best.model.Laura.FDR,
			   Tt$rel.ampl.ex, Tt$pval.ex, Tt$Qval.ex, Tt$rel.ampl.int, Tt$pval.int, Tt$Qval.int)
	
	xx = data.frame(matrix(unlist(xx), nrow=nrow(Tt), byrow=FALSE))
	colnames(xx) = c('Gene.names','BIC.Model', 'AIC.Model','AICc.Model','LRT.Model','BIC.Model.Old', 'AIC.Model.Old','AICc.Model.Old',
					 'LRT.Model.Old','BIC.Model.Laura', 'AIC.Model.Laura','AICc.Model.Laura','LRT.Model.Laura','LRT.Model.Laura.FDR', 'rel.ampl.ex', 'pval.ex', 'Qval.ex', 'rel.ampl.int', 'pval.int', 'Qval.int')
	ssigma.m = sqrt((exp(sigma.m.real^2)-1)*exp(sigma.m.real^2))
	ssigma.s = sqrt((exp(sigma.s.real^2)-1)*exp(sigma.s.real^2))
	
	## DATA OF David
	david.0 = c('Ddx17', 'Slc1a5', 'Stx2', 'Uba6', 'Zfp697', 'Crat', 'Pcsk4', 'Tra2a', 'Hnrnpc', 'Lpin3', 'Exog', 'Srf', 'Clec4b1', 'Rprd1a', 'Grk4', 'Pxmp4', 'Tbc1d24')
	david.1 = read.table('/Users/jiwang/Degradation_Laura/Code_Laura/Data_DGatfield/setI.txt',header=FALSE,as.is=c(1))
	mapping.ens = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/TF_inference/Elastic-net_analysis/genes_mapping.txt', header=FALSE, as.is=c(1,2,7))
	david.2 = read.csv('/Users/jiwang/Degradation_Laura/Code_Laura/Data_DGatfield/sig_stable_transcripts.csv',as.is=c(1:3))
	
	kk = match(david.1[,1], mapping.ens[,1])
	david.1$Gene.names = mapping.ens[kk,2]
	colnames(david.1)[1] = 'Ensemble'
	which(is.na(david.1[,2])==TRUE)
	
	## mRNAs rhythmically regulated by miRNAs (17 genes with high confidence)
	jj = match(david.0, xx[,1])
	length(which(!is.na(jj)==TRUE))
	jj = jj[which(!is.na(jj)==TRUE)]
	
	xx[jj,c(1:14)]
	
	for(index in c(2:14)) {print(length(which(xx[jj,index]==3|xx[jj,index]==4)));}
	
	counts = matrix(0, nrow=length(jj), ncol=2)
	counts[which(xx[jj,2]==3|xx[jj,2]==4),] = 3
	counts[-which(xx[jj,2]==3|xx[jj,2]==4),1] = 6
	counts[-which(xx[jj,2]==3|xx[jj,2]==4),2] = 5
	venn(list(GrpA=counts[,1],GrpB=counts[,2]))
	
	## mRNAs rhythmically regulated by miRNAs (190 genes with rhythmic mrna and non-rhythmic premrna in WT and nonrhythmic mrna and premrna in Dicer)
	ii = match(david.1[,2], xx[,1])
	length(which(!is.na(ii)==TRUE))
	ii = ii[which(!is.na(ii)==TRUE)]
	
	for(index in c(2:14))
	{
	#test = length(which(xx[ii,index]==3|xx[ii,2]==4));
		test = length(which(xx[ii,index]==3));
		
		print(test)
	}
	
	## mRNAs affected by miRNAs (167 genes showing difference in WT and Dicer)
	kk = match(david.2[,2], xx[,1])
	length(which(!is.na(kk)==TRUE))
	kk = kk[which(!is.na(kk)==TRUE)]
	
	for(index in c(2:14))
	{	
		#test = length(which(xx[kk,index]==3|xx[kk,index]==4));
		test = length(which(xx[kk,index]==3));
		print(test)
		
	}	
}



#### Analysis of David Gatfield Total RNA-seq raw data
version = "_my_8th_same_up_down_time_all_genes_sum"
version = '_my_8th_same_up_down_time_all_genes_model3_no_treated'
if(do.step)
{
	load(paste('my_genes_fits_results', version,'.Rdata', sep = ''))
	source('model_modify_new/functions.R')

	kk = which(!is.na(Tt$LRT.best.model)==TRUE)
	Tt = Tt[kk,]
	T = Tt


	#RNAseq.mRNA = read.table('../complementary_data/RNA-seq_DG/ctrl_mRNA.txt', header = TRUE)
	#RNAseq.premRNA = read.table('../complementary_data/RNA-seq_DG/ctrl_premRNA.txt', header = TRUE, fill = TRUE)
	#load('../complementary_data/RNA-seq_DG/eff.sizes.RData') # eff.size

	load('/Users/jiwang/Degradation_Laura/_x_MicroArray/complementary_data/RNA-seq_DG/Rdata_for_Laura.Rdata') #ctrl.premRNA
	R = ctrl.premRNA

	density.premRNA = R[,2:13]/R$premRNA.effsize
	density.mRNA = R[,14:25]/R$mRNA.effsize
	colnames(density.premRNA) = paste('dens',colnames(density.premRNA), sep = '.')
	colnames(density.mRNA) = paste('dens',colnames(density.mRNA), sep = '.')
	R = data.frame(R, density.premRNA, density.mRNA, stringsAsFactors = FALSE)
	original.R = R

	m = match(T$gene, R$gene_name)
	R = R[m,]


	j.dens.premRNA = grep('dens.ctrl.premRNA.ZT', colnames(R))
	j.dens.mRNA = grep('dens.ctrl.mRNA.ZT', colnames(R))

	R$mean.mRNA.dens  = apply(R[, j.dens.mRNA], 1, mean)
	R$mean.premRNA.dens  = apply(R[, j.dens.premRNA], 1, mean)

	rel.dens.mRNA = R[, j.dens.mRNA]/R$mean.mRNA.dens
	colnames(rel.dens.mRNA) = paste('rel',colnames(rel.dens.mRNA),sep = '.')
	rel.dens.premRNA = R[, j.dens.premRNA]/R$mean.premRNA.dens
	colnames(rel.dens.premRNA) = paste('rel',colnames(rel.dens.premRNA),sep = '.')

	R = data.frame(R, rel.dens.premRNA, rel.dens.mRNA, stringsAsFactors = FALSE)

	j.rel.premRNA = grep('rel.dens.ctrl.premRNA.ZT', colnames(R))
	j.rel.mRNA = grep('rel.dens.ctrl.mRNA.ZT', colnames(R))


	j.int = grep('.rel.ampl.int', colnames(T))
	j.ex = grep('.rel.ampl.ex', colnames(T))

	R$corr.with.microarray.premRNA = rep(0,nrow(R));
	R$corr.with.microarray.mRNA = rep(0,nrow(R));
	for(i in 1:nrow(R))
	{
		if((i%%100) == 0){cat(100*i/nrow(R),'\n')}
		R$corr.with.microarray.premRNA[i] = cor(unlist(R[i, j.rel.premRNA]),unlist(T[i, j.int[seq(1,23,by=2)]]))
		R$corr.with.microarray.mRNA[i] = cor(unlist(R[i, j.rel.mRNA]),unlist(T[i, j.ex[seq(1,23,by=2)]]))
	}



R$premRNA.FC = apply(R[,j.dens.premRNA],1,max)/ apply(R[,j.dens.premRNA],1,min)
R$mRNA.FC = apply(R[, j.dens.mRNA],1,max)/ apply(R[, j.dens.mRNA],1,min)


T$int.FC = apply(T[, j.int],1,max)/apply(T[, j.int],1,min)
T$ex.FC = apply(T[, j.ex],1,max)/apply(T[, j.ex],1,min)

j2 = which(T$LRT.best.model == 2)
j3 = which(T$LRT.best.model == 3)
j4 = which(T$LRT.best.model == 4)



zt.m = seq(0,46,by=2); #time-points for the microarray experiment
zt.r = seq(0,44,by=4); #time-points for the RNA-seq experiment




adjusted.dens.mRNA = R[,j.dens.mRNA]/R$mean.premRNA.dens
colnames(adjusted.dens.mRNA) = paste('adj',colnames(adjusted.dens.mRNA) , sep = '.')

R = data.frame(R, adjusted.dens.mRNA, stringsAsFactors = FALSE)

j.mRNA = grep('adj.dens.ctrl.mRNA.', colnames(R) )


save(R, file = 'Rdata/RNA-seq_DGatfield.Rdata')

	par(mfcol = c(1,2))
	hist(log2(R$mean.mRNA.dens))
	hist(log2(R$mean.premRNA.dens))



pdf('plots/RNA_seq_DGatfield_general_comparison.pdf', width = 6.6, height = 9)
par(cex = 0.7)
par(mfrow = c(3,2))


plot(R$mean.mRNA.dens, R$mean.premRNA.dens, pch = 16, cex = 0.5, col = rgb(0,0,0,0.3), log = 'xy')
abline(a = 0, b = 1, col = 'green3')

plot(R$mean.mRNA.dens-R$mean.premRNA.dens, R$mean.premRNA.dens, pch = 16, cex = 0.5, col = rgb(0,0,0,0.3), log = 'xy')
abline(a = 0, b = 1, col = 'green3')

plot(R$mean.mRNA.dens, T$exon.median, pch = 16, cex = 0.5, col = rgb(0,0,0,0.3), log = 'xy')
plot(R$mean.premRNA.dens, T$intron.median, pch = 16, cex = 0.5, col = rgb(0,0,0,0.3), log = 'xy')

plot( T$exon.median, T$intron.median, pch = 16, cex = 0.5, col = rgb(0,0,0,0.3), log = 'xy')

par(mfcol = c(2,2))


plot(R$premRNA.FC,R$mRNA.FC, pch = 16, cex = 0.5,col = rgb(0,0,0,0.3) , log = 'xy')
abline(a = 0, b = 1, col = 'green3')


plot(T$int.FC,T$ex.FC, pch = 16, cex = 0.5,col = rgb(0,0,0,0.3) , log = 'xy')
abline(a = 0, b = 1, col = 'green3')

plot(T$int.FC,R$premRNA.FC, pch = 16, cex = 0.5,col = rgb(0,0,0,0.3) , log = 'xy')
abline(a = 0, b = 1, col = 'green3')

plot(T$ex.FC,R$mRNA.FC, pch = 16, cex = 0.5,col = rgb(0,0,0,0.3) , log = 'xy')
abline(a = 0, b = 1, col = 'green3')


par(mfrow = c(3,2))


plot(T$int.FC[j2],R$premRNA.FC[j2], pch = 16, cex = 1,col = rgb(0,0.7,0,0.3) , log = 'xy')
abline(a = 0, b = 1, col = 'green3')
plot(T$ex.FC[j2],R$mRNA.FC[j2], pch = 16, cex = 1,col = rgb(0,0.7,0,0.3) , log = 'xy')
abline(a = 0, b = 1, col = 'green3')


plot(T$int.FC[j3],R$premRNA.FC[j3], pch = 16, cex = 1,col = rgb(0.8,0,0,0.3) , log = 'xy')
abline(a = 0, b = 1, col = 'green3')
plot(T$ex.FC[j3],R$mRNA.FC[j3], pch = 16, cex = 1,col = rgb(0.8,0,0,0.3) , log = 'xy')
abline(a = 0, b = 1, col = 'green3')

plot(T$int.FC[j4],R$premRNA.FC[j4], pch = 16, cex = 1,col = rgb(0,0,0,0.3) , log = 'xy')
abline(a = 0, b = 1, col = 'green3')
plot(T$ex.FC[j4],R$mRNA.FC[j4], pch = 16, cex = 1,col = rgb(0,0,0,0.3) , log = 'xy')
abline(a = 0, b = 1, col = 'green3')


par(mfrow = c(4,2))
layout(matrix(c(1:12),nrow = 4, ncol = 3, byrow = TRUE), width = c(1,1,0.2))
par(mar = c(3,3,2.5,0.2)+0.1, mgp = c(2,0.5,0), las = 1, tcl = -0.3, cex.main = 1,cex = 0.6)
n.col = 5


o = order(T$phase.int[j2])
image(t(as.matrix(T[j2[o],j.int])), zlim = c(0,2), main = 'micro-array intron - phase order intron')
image(t(as.matrix(R[j2[o],j.rel.premRNA])), zlim = c(0,2), main = 'RNA-seq intron')
image(t(R$corr.with.microarray.premRNA[j2[o]]), col = cm.colors(n.col), axes = FALSE, main = 'corr', zlim = c(-1,1))


image(t(as.matrix(T[j2[o],j.ex])), zlim = c(0,2), main = 'micro-array exon - phase order intron')
image(t(as.matrix(R[j2[o],j.rel.mRNA])), zlim = c(0,2), main = 'RNA-seq exon')
image(t(R$corr.with.microarray.mRNA[j2[o]]), col = cm.colors(n.col), axes = FALSE, main = 'corr', zlim = c(-1,1))


o = order(T$phase.ex[j2])
image(t(as.matrix(T[j2[o],j.ex])), zlim = c(0,2), main = 'micro-array exon - phase order exon')
image(t(as.matrix(R[j2[o],j.rel.mRNA])), zlim = c(0,2), main = 'RNA-seq exon')
image(t(R$corr.with.microarray.mRNA[j2[o]]), col = cm.colors(n.col), axes = FALSE, main = 'corr', zlim = c(-1,1))


image(t(as.matrix(T[j2[o],j.int])), zlim = c(0,2), main = 'micro-array intron - phase order exon')
image(t(as.matrix(R[j2[o],j.rel.premRNA])), zlim = c(0,2), main = 'RNA-seq intron')
image(t(R$corr.with.microarray.premRNA[j2[o]]), col = cm.colors(n.col), axes = FALSE, main = 'corr', zlim = c(-1,1))



o = order(T$phase.int[j3])
image(t(as.matrix(T[j3[o],j.int])), zlim = c(0,2), main = 'micro-array intron - phase order intron')
image(t(as.matrix(R[j3[o],j.rel.premRNA])), zlim = c(0,2), main = 'RNA-seq intron')
image(t(R$corr.with.microarray.premRNA[j3[o]]), col = cm.colors(n.col), axes = FALSE, main = 'corr', zlim = c(-1,1))


image(t(as.matrix(T[j3[o],j.ex])), zlim = c(0,2), main = 'micro-array exon - phase order intron')
image(t(as.matrix(R[j3[o],j.rel.mRNA])), zlim = c(0,2), main = 'RNA-seq exon')
image(t(R$corr.with.microarray.mRNA[j3[o]]), col = cm.colors(n.col), axes = FALSE, main = 'corr', zlim = c(-1,1))

o = order(T$phase.ex[j3])
image(t(as.matrix(T[j3[o],j.ex])), zlim = c(0,2), main = 'micro-array exon - phase order exon')
image(t(as.matrix(R[j3[o],j.rel.mRNA])), zlim = c(0,2), main = 'RNA-seq exon')
image(t(R$corr.with.microarray.mRNA[j3[o]]), col = cm.colors(n.col), axes = FALSE, main = 'corr', zlim = c(-1,1))

image(t(as.matrix(T[j3[o],j.int])), zlim = c(0,2), main = 'micro-array intron - phase order exon')
image(t(as.matrix(R[j3[o],j.rel.premRNA])), zlim = c(0,2), main = 'RNA-seq intron')
image(t(R$corr.with.microarray.premRNA[j3[o]]), col = cm.colors(n.col), axes = FALSE, main = 'corr', zlim = c(-1,1))



o = order(T$phase.int[j4])
image(t(as.matrix(T[j4[o],j.int])), zlim = c(0,2), main = 'micro-array intron - phase order intron')
image(t(as.matrix(R[j4[o],j.rel.premRNA])), zlim = c(0,2), main = 'RNA-seq intron')
image(t(R$corr.with.microarray.premRNA[j4[o]]), col = cm.colors(n.col), axes = FALSE, main = 'corr', zlim = c(-1,1))

image(t(as.matrix(T[j4[o],j.ex])), zlim = c(0,2), main = 'micro-array exon - phase order intron')
image(t(as.matrix(R[j4[o],j.rel.mRNA])), zlim = c(0,2), main = 'RNA-seq exon')
image(t(R$corr.with.microarray.mRNA[j4[o]]), col = cm.colors(n.col), axes = FALSE, main = 'corr', zlim = c(-1,1))

o = order(T$phase.ex[j4])
image(t(as.matrix(T[j4[o],j.ex])), zlim = c(0,2), main = 'micro-array exon - phase order exon')
image(t(as.matrix(R[j4[o],j.rel.mRNA])), zlim = c(0,2), main = 'RNA-seq exon')
image(t(R$corr.with.microarray.mRNA[j4[o]]), col = cm.colors(n.col), axes = FALSE, main = 'corr', zlim = c(-1,1))

image(t(as.matrix(T[j4[o],j.int])), zlim = c(0,2), main = 'micro-array intron - phase order exon')
image(t(as.matrix(R[j4[o],j.rel.premRNA])), zlim = c(0,2), main = 'RNA-seq intron')
image(t(R$corr.with.microarray.premRNA[j4[o]]), col = cm.colors(n.col), axes = FALSE, main = 'corr', zlim = c(-1,1))

dev.off()


plot(T$phase.int[j2], R$corr.with.microarray.premRNA[j2], pch =16, col = 2, ylim = c(-1,1))
points(T$phase.int[j3], R$corr.with.microarray.premRNA[j3], pch =16, col = 3)
points(T$phase.int[j4], R$corr.with.microarray.premRNA[j4], pch =16, col = 4)

plot(T$phase.ex[j2], R$corr.with.microarray.mRNA[j2], pch =16, col = 2, ylim = c(-1,1))
points(T$phase.ex[j3], R$corr.with.microarray.mRNA[j3], pch =16, col = 3)
points(T$phase.ex[j4], R$corr.with.microarray.mRNA[j4], pch =16, col = 4)


pdf('plots/RNA_seq_DGatfield_gene_comparison.pdf', width = 8, height = 12)
par(mar = c(3,3,2.5,0.2)+0.1, mgp = c(2,0.5,0), las = 1, tcl = -0.3, cex.main = 1,cex = 0.6)

for(m in 2:4)
{
	eval(parse(text = paste('j = j',m,sep = '')))
	cat(length(j),'\n')
	
	par(mfrow = c(6,4))
	
	for(i in j){
		plot(zt.r, zt.r,type = 'n', ylim = c(0,2.5), xlab = 'time', ylab = 'rel. signal', axes = FALSE, main = paste(T$gene[i],'\n corr : ',round(cor(unlist(R[i, j.rel.premRNA]),unlist(T[i, j.int[seq(1,23,by=2)]])),digits = 2)))
		axis(1,at = zt.r); axis(2); box()
		points(zt.r, R[i, j.rel.premRNA], col = 'green3', type = 'b', pch = 17, lwd = 2)
		points(zt.m, T[i, j.int], col = 'green3', type = 'b', pch = 1, lty = 2, cex = 0.6)
		
		plot(zt.r, zt.r,type = 'n', ylim = c(0,2.5), xlab = 'time', ylab = 'rel. signal', axes = FALSE, main = paste(T$gene[i],'\n corr : ',round(cor(unlist(R[i, j.rel.mRNA]),unlist(T[i, j.ex[seq(1,23,by=2)]])),digits = 2)))
		axis(1,at = zt.r); axis(2); box()
		points(zt.r, R[i, j.rel.mRNA], col = 'steelblue', type = 'b', pch = 17, lwd = 2)
		points(zt.m, T[i, j.ex], col = 'steelblue', type = 'b', pch = 1, lty = 2, cex = 0.6)
		
		plot(zt.r, zt.r,type = 'n', ylim = c(0,2.5), xlab = 'time', ylab = 'rel. signal', axes = FALSE)
		axis(1,at = zt.r); axis(2); box()
		points(zt.r, R[i, j.rel.premRNA], col = 'green3', type = 'b', pch = 17, lwd = 2)
		points(zt.r, R[i, j.rel.mRNA], col = 'steelblue', type = 'b', pch = 17, lwd = 2)
		
		plot(zt.r, zt.r,type = 'n', ylim = c(0,2.5), xlab = 'time', ylab = 'rel. signal', axes = FALSE,  main = paste('LRT best model : ',m))
		axis(1,at = zt.r); axis(2); box()
		points(zt.m, T[i, j.int], col = 'green3', type = 'b', pch = 1, lty = 2, cex = 0.6)
		points(zt.m, T[i, j.ex], col = 'steelblue', type = 'b', pch = 1, lty = 2, cex = 0.6)
	}
	
	
}

dev.off()




m = match(c('Tmem33','P4ha2','Sec14l1','Ndrg1','Exog','Tfrc','Pcsk4','Crip2','Loxl4', 'Erbb4','Acot11','Tspan4','Ppard','Upp2','Per2','Acacb','Oas1a','Ak2','Slc24a6','Cdcp1','Tubb2a','Dhrs9'), R$gene)


pdf('plots/RNA_seq_DGatfield_gene_comparison_candidates.pdf', width = 8, height = 12)
par(mar = c(3,3,2.5,0.2)+0.1, mgp = c(2,0.5,0), las = 1, tcl = -0.3, cex.main = 1,cex = 0.6)

for(i in m){
	par(mfrow = c(3,2))
	
	
	plot(zt.r, zt.r,type = 'n', ylim = c(0,2.5), xlab = 'time', ylab = 'rel. signal', axes = FALSE, main = paste(T$gene[i],'\n corr : ',round(cor(unlist(R[i, j.rel.premRNA]),unlist(T[i, j.int[seq(1,23,by=2)]])),digits = 2)))
	axis(1,at = zt.r); axis(2); box()
	points(zt.r, R[i, j.rel.premRNA], col = 'green3', type = 'b', pch = 17, lwd = 2)
	points(zt.m, T[i, j.int], col = 'green3', type = 'b', pch = 1, lty = 2, cex = 0.6)
	
	plot(zt.r, zt.r,type = 'n', ylim = c(0,2.5), xlab = 'time', ylab = 'rel. signal', axes = FALSE, main = paste(T$gene[i],'\n corr : ',round(cor(unlist(R[i, j.rel.mRNA]),unlist(T[i, j.ex[seq(1,23,by=2)]])),digits = 2)))
	axis(1,at = zt.r); axis(2); box()
	points(zt.r, R[i, j.rel.mRNA], col = 'steelblue', type = 'b', pch = 17, lwd = 2)
	points(zt.m, T[i, j.ex], col = 'steelblue', type = 'b', pch = 1, lty = 2, cex = 0.6)
	
	plot(zt.r, zt.r,type = 'n', ylim = c(0,2.5), xlab = 'time', ylab = 'rel. signal', axes = FALSE)
	axis(1,at = zt.r); axis(2); box()
	points(zt.r, R[i, j.rel.premRNA], col = 'green3', type = 'b', pch = 17, lwd = 2)
	points(zt.r, R[i, j.rel.mRNA], col = 'steelblue', type = 'b', pch = 17, lwd = 2)
	
	plot(zt.r, zt.r,type = 'n', ylim = c(0,2.5), xlab = 'time', ylab = 'rel. signal', axes = FALSE,  main = paste('LRT best model : ',T$LRT.best.model[i]))
	axis(1,at = zt.r); axis(2); box()
	points(zt.m, T[i, j.int], col = 'green3', type = 'b', pch = 1, lty = 2, cex = 0.6)
	points(zt.m, T[i, j.ex], col = 'steelblue', type = 'b', pch = 1, lty = 2, cex = 0.6)
	
	plot(zt.r, zt.r,type = 'n', ylim = range(0, R[i, j.dens.premRNA], R[i, j.dens.mRNA]), xlab = 'time', ylab = 'abs. signal', axes = FALSE)
	axis(1,at = zt.r); axis(2); box()
	points(zt.r, R[i, j.dens.premRNA], col = 'green3', type = 'b', pch = 17, lwd = 2)
	points(zt.r, R[i, j.dens.mRNA], col = 'steelblue', type = 'b', pch = 17, lwd = 2)
}

dev.off()




########

#test to fit 4h time resolution data

s.r = compute.sigmoid(t = zt.r)
s.m = compute.sigmoid(t = zt.m)

m.r = compute.m.sigmoid(t = zt.r, eps.gamma = 1)
m.m = compute.m.sigmoid(t = zt.m, eps.gamma = 1)


par(mfrow = c(3,1))

plot(zt.r, zt.r,type = 'n', ylim = c(0,2), xlab = 'time', ylab = 'rel. signal')

points(zt.r, s.r, type = 'b', col = 'green3', pch = 15)
points(zt.m, s.m, type = 'b', col = 'green3', pch = 16, lty = 2)
points(zt.r, m.r, type = 'b', col = 'steelblue', pch = 15)
points(zt.m, m.m, type = 'b', col = 'steelblue', pch = 16,lty = 2)

plot(zt.r, zt.r,type = 'n', ylim = c(0,2), xlab = 'time', ylab = 'rel. signal')
points(zt.r, s.r, type = 'b', col = 'green3', pch = 15)
points(zt.r, m.r, type = 'b', col = 'steelblue', pch = 15)

plot(zt.r, zt.r,type = 'n', ylim = c(0,2), xlab = 'time', ylab = 'rel. signal')
points(zt.m, s.m, type = 'b', col = 'green3', pch = 16, lty = 2)
points(zt.m, m.m, type = 'b', col = 'steelblue', pch = 16,lty = 2)

true.par = c(2,2,12,8,8)
par.init = c(2,3,18,8,8)

zt = zt.r
f2min(par = true.par, model = 2, M = m.r, S = s.r, parametrization = 'sigmoid', method.intern = "integration", sum.species = FALSE)
f2min(par = par.init, model = 2, M = m.r, S = s.r, parametrization = 'sigmoid', method.intern = "integration", sum.species = FALSE)
bounds = set.bounds(model = 2, parametrization = 'sigmoid'); upper = bounds$upper; lower = bounds$lower
opt = optim(par.init, f2min, M = m.r, S = s.r, model = 2, parametrization = 'sigmoid', method.intern = "integration", debug = TRUE, sum.species = FALSE, method = 'L-BFGS-B', lower = lower, upper = upper)
f2min(par = opt$par, model = 2, M = m.r, S = s.r, parametrization = 'sigmoid', method.intern = "integration", sum.species = FALSE)


zt = zt.m
f2min(par = true.par, model = 2, M = m.m, S = s.m, parametrization = 'sigmoid', method.intern = "integration", sum.species = FALSE)
f2min(par = par.init, model = 2, M = m.m, S = s.m, parametrization = 'sigmoid', method.intern = "integration", sum.species = FALSE)
opt = optim(par.init, f2min, M = m.m, S = s.m, model = 2, parametrization = 'sigmoid', method.intern = "integration", debug = TRUE, sum.species = FALSE, method = 'L-BFGS-B', lower = lower, upper = upper)
f2min(par = opt$par, model = 2, M = m.m, S = s.m, parametrization = 'sigmoid', method.intern = "integration", sum.species = FALSE)


m = match(c('Ppard','Upp2','Acot11','Tspan4','Tfrc','Loxl4', 'Erbb4','Pcsk4','Crip2'), R$gene)

source('functions.R')

intense.debug = FALSE

make.optimization(T = R, i = m[1], model = 3, parametrization = 'sigmoid', debug = TRUE, method = 'integration', sum.species = FALSE, zt = zt.r, i.ex = j.mRNA, i.int = j.rel.premRNA, absolute.signal = TRUE)

make.fits.with.all.models.for.one.gene(T = R, gene.index = m[1], debug = TRUE, parametrization = 'sigmoid',  method = 'integration', sum.species = FALSE, zt = zt.r, i.ex = j.mRNA, i.int = j.rel.premRNA, absolute.signal = TRUE)



######## FITS ON VITAL-IT


plot = TRUE; load.for.plot = TRUE
do.step = TRUE; load = TRUE


version = '_RNAseq_DG'

if(do.step){
	cat('FIT RNA-seq DATA ON VITAL-IT\n')
	load('Rdata/RNA-seq_DGatfield.Rdata'); cat('data loaded \n')
	
#prepare the data for the fits
	save(R, file = 'stuff_for_vital_IT_RNA-seq_FIT/genes_ready_for_fits.Rdata')
	gene_list = R$gene[!is.na(R$gene)]; #gene_list = gene_list[1:10]
	write.table(gene_list, file = 'stuff_for_vital_IT_RNA-seq_FIT/gene_list.txt', quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE)
	
	
## ON VITAL-IT 
	
#clear the directory
	cmd = "ssh lsymul@rserv.vital-it.ch 'rm -r /home/lsymul/microarray/*' " 
	system(cmd)
	
#copy the data and the code on vital-IT
	cmd = 'scp stuff_for_vital_IT_RNA-seq_FIT/* lsymul@rserv.vital-it.ch:/home/lsymul/microarray/.'
	system(cmd)
	cmd = 'scp -r functions.R lsymul@rserv.vital-it.ch:/home/lsymul/microarray/.'
	system(cmd)
	
# run the fits : 
	cmd = "ssh lsymul@rserv.vital-it.ch 'bash -s' < stuff_for_vital_IT_RNA-seq_FIT/to_run_on_vital_IT.sh"
	system(cmd)
	
# concatenate all fits results together when fits are all done
	cmd = "ssh lsymul@rserv.vital-it.ch 'perl /home/lsymul/microarray/cat_all_fits_results.pl'"
	system(cmd)
	
# wait for everything to be done and copy the data on the local computer
	cmd = "ssh lsymul@rserv.vital-it.ch 'ls -l /home/lsymul/microarray/fits/all_genes_fits_results_with_header.txt'"
	res = system(cmd, intern = TRUE)
	
	while(length(res) == 0){
		Sys.sleep(60*5) # wait for 10 minutes
		cat('... wait for the fits to finish\n') 
		cmd = "ssh lsymul@rserv.vital-it.ch 'ls -l /home/lsymul/microarray/fits/all_genes_fits_results_with_header.txt'"
		res = system(cmd, intern = TRUE)
	}
	
	cat('........................................ fits finished on vital-IT\n')	
	
	cmd = paste("scp lsymul@rserv.vital-it.ch:/home/lsymul/microarray/fits/all_genes_fits_results_with_header.txt ../complementary_data/Vital-IT_fits_optim/all_genes_fits_results_with_header",version,".txt", sep = '')
	system(cmd)
	
	cat('........................................ fits results imported on local computer\n')
	
	file =  paste("../complementary_data/Vital-IT_fits_optim/all_genes_fits_results_with_header",version,".txt", sep = '')	
	fits.results = read.table(file = file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
	if(nrow(fits.results)<0.9*length(gene_list)){cat('........................................ fits results missing for more than 10% of genes\n')}
	m = match(gene_list, fits.results$gene)
	if(sum(is.na(m))>0){cat('........................................ fits results missing for',sum(is.na(m)),'genes\n')}
	m = match(R$gene, fits.results$gene)
	
	fits.results = fits.results[m,]
	
	R = data.frame(R, fits.results[,-1],stringsAsFactors = FALSE)	
	save(R, file = paste('Rdata/RNA-seq_DGatfield_genes_fits_optim_vital-IT',version,'.Rdata',sep = ''))
	
	
}



if(do.step){
	cat('LIKELIHOOD RATIO TEST\n')
	if(load){load(paste('Rdata/RNA-seq_DGatfield_genes_fits_optim_vital-IT',version,'.Rdata',sep = ''))}
	source('functions.R')
	
	pval = 0.05
	LRT.and.best.model = Likelihood.Ratio.Test(T = R, pval = pval, direct = FALSE, FDR = TRUE, zt = seq(0,44,by=4))
	R = data.frame(R, LRT.and.best.model, stringsAsFactors = FALSE)
	
	par(mfrow = c(1,3))
	boxplot(R$LRT.quality ~ R$LRT.best.model, col = c(1:4), pch = 16, cex = 0.6); abline(h= 1)
	h = hist( R$LRT.best.model, breaks = seq(0.5,4.5,by=1), col = 1:4, border = 'NA')
	barplot(h$counts[2:4])
	
	CIRC = c("Per1","Per2","Per3","Cry1","Cry2","Dbp","Tef","Hlf","Nr1d1","Nr1d2","Rora","Rorb","Rorc","Arntl","Bmal1","Clock","Npas2","Bhlhe40","Bhlhe41","Cirbp","Tfrc","Hamp","Hamp2","Nr4a2","Wee1", "Por")
	k = match(CIRC,R$gene)
	cbind(as.character(R$gene[k]), R$LRT.best.model[k])
	
	save(R, file = paste("Rdata/gene_fit_results_and_LRT",version,".Rdata",sep = ''))
}


noise.T = abs(T$ZT00.rel.ampl.ex-T$ZT24.rel.ampl.ex)/(T$ZT00.rel.ampl.ex+T$ZT24.rel.ampl.ex)
noise.R = abs(R$rel.dens.ctrl.mRNA.ZT0 - R$rel.dens.ctrl.mRNA.ZT24)/(R$rel.dens.ctrl.mRNA.ZT0+R$rel.dens.ctrl.mRNA.ZT24)

hist(noise.R, breaks = 1000, col = 5, border= 5)
hist(noise.T, breaks = 1000, add = TRUE, col = 7, border= 7)


if(do.step){
	cat('AIC\n')
	if(load){load(paste("Rdata/gene_fit_results_and_LRT",version,".Rdata",sep = ''))}
	source('functions.R')
	
	noise = sqrt(0.1)
	AIC = my.AIC(T = R, number.of.data.point= 24, fact = 1/noise^2)
	AIC.c = my.AIC(T = R, number.of.data.point= 24, fact = 1/noise^2, correction = TRUE)
	
	h1 = hist(as.numeric(AIC$best.model), breaks = seq(0.5,4.5,by = 1), plot = FALSE)
	h2 = hist(as.numeric(AIC.c$best.model), breaks = seq(0.5,4.5,by = 1), plot = FALSE)
	h3 = hist(R$LRT.best.model, breaks = seq(0.5,4.5,by = 1), plot = FALSE)
	
	par(mfrow = c(1,4))
	
	H = cbind(h1$counts, h2$counts, h3$counts)
	b = barplot(t(H), names.arg = paste('model',1:4), beside = TRUE, col = c('steelblue2','yellowgreen','indianred2'), border = NA, ylim = c(0,max(H)*1.05))
	text(b, t(H),t(H), col = c('steelblue2','yellowgreen','indianred2'), pos = 3)
	legend('topright', legend = c('AIC','corr AIC','LRT') , fill = c('steelblue2','yellowgreen','indianred2'), border = NA, bty = 'n')
	
	
	plot(as.numeric(AIC$best.model)+rnorm(nrow(R),sd = 0.1), R$LRT.best.model+rnorm(nrow(R),sd = 0.1), pch = 16, cex = 1, col = rgb(0,0.5,0.8,0.3))
	
	plot(as.numeric(AIC.c$best.model)+rnorm(nrow(R),sd = 0.1), R$LRT.best.model+rnorm(nrow(R),sd = 0.1), pch = 16, cex = 1, col = rgb(0.3,0.8,0,0.3))
	
	plot(as.numeric(AIC.c$best.model)+rnorm(nrow(R),sd = 0.1),as.numeric(AIC$best.model)+rnorm(nrow(R),sd = 0.1), pch = 16, cex = 1, col = rgb(0.3,0.8,0,0.3))
	
	
	R = data.frame(R, AIC, stringsAsFactors = FALSE)
	
	save(R, file = paste("Rdata/gene_fit_results_and_LRT_and_AIC",version,".Rdata",sep = ''))
}


load('Rdata/gene_fit_results_and_LRT_and_AIC_RNAseq_DG.Rdata')

load('Rdata/genes_fits_results_quant_norm.Rdata')
source('functions.R')


m = match( T$gene, R$gene)

R = R[m, ]


R$best.model = as.numeric(as.character(as.matrix(R$best.model)))
T$best.model = as.numeric(as.character(as.matrix(T$best.model)))

R.bm = R$LRT.best.model


j1 = which((R.bm ==1)&(R$gamma.m1!=Inf))
j2 = which(R.bm ==2)
j3 = which(R.bm ==3)
j4 = which(R.bm ==4)

j1.T = which(T$LRT.best.model==1)
j2.T = which(T$LRT.best.model==2)
j3.T = which(T$LRT.best.model==3)
j4.T = which(T$LRT.best.model==4)


par(mfrow = c(4,4))

plot(R$corr.with.microarray.premRNA,R$corr.with.microarray.mRNA, pch = 21, cex = 1, bg = R.bm, col = T$LRT.best.model )
abline(a = 0, b = 1)


breaks = seq(-1,1,len =51)
h1 = hist(R$corr.with.microarray.premRNA[j1], breaks = breaks, plot = FALSE)
h2 = hist(R$corr.with.microarray.premRNA[j2], breaks = breaks, plot = FALSE)
h3 = hist(R$corr.with.microarray.premRNA[j3], breaks = breaks, plot = FALSE)
h4 = hist(R$corr.with.microarray.premRNA[j4], breaks = breaks, plot = FALSE)


plot(c(h1$breaks[1], h1$breaks, h1$breaks[length(h1$breaks)]),c(0, h1$counts, h1$counts[length(h1$counts)],0), type = 'n', main = 'Correlation with exon-array - pre-mRNA',  xlab = 'correlation', ylab = '# of genes', ylim = range(h1$counts, h2$counts, h3$counts, h4$counts), axes = TRUE)
polygon(c(rep(h1$breaks, each = 2)), c(0, rep(h1$counts, each=2),0), col = 1+4, border = 1, lwd = 1)
polygon(c(rep(h2$breaks, each = 2)), c(0, rep(h2$counts, each=2),0), col = 2+4, border = 2, lwd = 1)
polygon(c(rep(h3$breaks, each = 2)), c(0, rep(h3$counts, each=2),0), col = 3+4, border = 3, lwd = 1)
polygon(c(rep(h4$breaks, each = 2)), c(0, rep(h4$counts, each=2),0), col = 4+4, border = 4, lwd = 1)
abline(v = 0)


breaks = seq(-1,1,len =51)
h1 = hist(R$corr.with.microarray.mRNA[j1], breaks = breaks, plot = FALSE)
h2 = hist(R$corr.with.microarray.mRNA[j2], breaks = breaks, plot = FALSE)
h3 = hist(R$corr.with.microarray.mRNA[j3], breaks = breaks, plot = FALSE)
h4 = hist(R$corr.with.microarray.mRNA[j4], breaks = breaks, plot = FALSE)


plot(c(h1$breaks[1], h1$breaks, h1$breaks[length(h1$breaks)]),c(0, h1$counts, h1$counts[length(h1$counts)],0), type = 'n', main = 'Correlation with exon-array - mRNA',  xlab = 'correlation', ylab = '# of genes', ylim = range(h1$counts, h2$counts, h3$counts, h4$counts), axes = TRUE)
polygon(c(rep(h1$breaks, each = 2)), c(0, rep(h1$counts, each=2),0), col = 1+4, border = 1, lwd = 1)
polygon(c(rep(h2$breaks, each = 2)), c(0, rep(h2$counts, each=2),0), col = 2+4, border = 2, lwd = 1)
polygon(c(rep(h3$breaks, each = 2)), c(0, rep(h3$counts, each=2),0), col = 3+4, border = 3, lwd = 1)
polygon(c(rep(h4$breaks, each = 2)), c(0, rep(h4$counts, each=2),0), col = 4+4, border = 4, lwd = 1)
abline(v = 0)





corr = rep(NA, nrow(T));
for(i in 1:nrow(T)){corr[i] = cor(c(R$AIC.m1[i],R$AIC.m2[i],R$AIC.m3[i],R$AIC.m4[i]),c(T$AIC.m1[i],T$AIC.m2[i],T$AIC.m3[i],T$AIC.m4[i]))}
j = c(j2,j3,j4)
#plot(R$corr.with.microarray.mRNA,corr,  bg = R$best.model+4, col = T$LRT.best.model+4, pch = 21, cex = 1, lwd = 2)
plot(R$corr.with.microarray.mRNA[j],corr[j],  bg = R.bm[j], col = T$LRT.best.model[j], pch = 21, cex = 1, lwd = 2)



plot(R.bm+rnorm(nrow(T),sd = 0.1), T$LRT.best.model+rnorm(nrow(T),sd = 0.1), pch = 16, cex = 1, col = rgb(0,0.5,0.8,0.3))
for(i in 1:4){for(j in 1:4){text(i,j,sum((R.bm ==i)&(T$LRT.best.model==j),na.rm = TRUE), col = 'black',font = 2)}}

plot(R.bm +rnorm(nrow(T),sd = 0.1), T$best.model+rnorm(nrow(T),sd = 0.1), pch = 16, cex = 1, col = rgb(0,0.5,0.8,0.3))
for(i in 1:4){for(j in 1:4){text(i,j,sum((R$best.model==i)&(T$best.model ==j),na.rm = TRUE), col = 'black',font = 2)}}


R.test = sample(R.bm)
plot(R.test +rnorm(nrow(T),sd = 0.1), T$LRT.best.model +rnorm(nrow(T),sd = 0.1), pch = 16, cex = 1, col = rgb(0,0.5,0.8,0.3))
for(i in 1:4){for(j in 1:4){text(i,j,sum((R.test ==i)&(T$LRT.best.model ==j),na.rm = TRUE), col = 'black',font = 2)}}

R.test = sample(R.bm)
plot(R.test +rnorm(nrow(T),sd = 0.1), T$LRT.best.model +rnorm(nrow(T),sd = 0.1), pch = 16, cex = 1, col = rgb(0,0.5,0.8,0.3))
for(i in 1:4){for(j in 1:4){text(i,j,sum((R.test ==i)&(T$LRT.best.model ==j),na.rm = TRUE), col = 'black',font = 2)}}

h = hist(R.bm, breaks = seq(0.5,4.5), plot = FALSE)
barplot(cbind(h$counts/sum(h$counts), c(0,h$counts[2:4]/sum(h$counts[2:4]))), beside = FALSE, col=1:4, border= NA, main = 'RNA-seq')

h = hist(T$LRT.best.model, breaks = seq(0.5,4.5), plot = FALSE)
barplot(cbind(h$counts/sum(h$counts), c(0,h$counts[2:4]/sum(h$counts[2:4]))), beside = FALSE, col=1:4, border= NA, main = 'Exonarray')


R.gamma = rep(NA, nrow(T))
R.gamma[j1] = R$gamma.m1[j1]
R.gamma[j2] = R$gamma.m2[j2]
R.gamma[j3] = R$gamma.m3[j3]
R.gamma[j4] = R$gamma.m4[j4]

R.hl = log(2)/R.gamma

breaks = lseq(0.001,300, len = 71)

h1 = hist(R.hl[j1], breaks = breaks, plot = FALSE)
h2 = hist(R.hl[j2], breaks = breaks, plot = FALSE)
h3 = hist(R.hl[j3], breaks = breaks, plot = FALSE)
h4 = hist(R.hl[j4], breaks = breaks, plot = FALSE)

plot(c(h1$breaks[1], h1$breaks, h1$breaks[length(h1$breaks)]),c(0, h1$counts, h1$counts[length(h1$counts)],0), log = 'x', type = 'n', main = 'Maximal half-lives',  xlab = 'half-lives [h]', ylab = '# of genes', ylim = range(h1$counts, h2$counts, h3$counts, h4$counts), axes = FALSE)
polygon(c(rep(h1$breaks, each = 2)), c(0, rep(h1$counts, each=2),0), col = 1+4, border = 1, lwd = 1)
polygon(c(rep(h2$breaks, each = 2)), c(0, rep(h2$counts, each=2),0), col = 2+4, border = 2, lwd = 1)
polygon(c(rep(h3$breaks, each = 2)), c(0, rep(h3$counts, each=2),0), col = 3+4, border = 3, lwd = 1)
polygon(c(rep(h4$breaks, each = 2)), c(0, rep(h4$counts, each=2),0), col = 4+4, border = 4, lwd = 1)
axis(1, at = c(0.05,0.1,0.2,0.5,1,2,5,10,20,50))
axis(2); box()

plot(c(h1$breaks[1], h1$breaks, h1$breaks[length(h1$breaks)]),c(0, h1$counts/sum(h1$counts), h1$counts[length(h1$counts)]/sum(h1$counts),0), log = 'x', type = 'n', main = 'Maximal half-lives',  xlab = 'half-lives [h]', ylab = '# of genes', ylim = range(h1$counts/sum(h1$counts), h2$counts/sum(h2$counts), h3$counts/sum(h3$counts), h4$counts/sum(h4$counts)), axes = FALSE)
polygon(c(rep(h1$breaks, each = 2)), c(0, rep(h1$counts/sum(h1$counts), each=2),0), col = 1+4, border = 1, lwd = 1)
polygon(c(rep(h2$breaks, each = 2)), c(0, rep(h2$counts/sum(h2$counts), each=2),0), col = 2+4, border = 2, lwd = 1)
polygon(c(rep(h3$breaks, each = 2)), c(0, rep(h3$counts/sum(h3$counts), each=2),0), col = 3+4, border = 3, lwd = 1)
polygon(c(rep(h4$breaks, each = 2)), c(0, rep(h4$counts/sum(h4$counts), each=2),0), col = 4+4, border = 4, lwd = 1)
axis(1, at = c(0.05,0.1,0.2,0.5,1,2,5,10,20,50))
axis(2); box()

plot(R.hl[intersect(j2,j2.T)], log(2)/T$gamma.m2[intersect(j2,j2.T)], pch = 16, col = 2, log = 'xy', xlab = 'RNA-seq', ylab = 'Exonarray', main = 'half-lives')
points(R.hl[intersect(j3,j3.T)], log(2)/T$gamma.m3[intersect(j3,j3.T)], pch = 16, col = 3)
points(R.hl[intersect(j4,j4.T)], log(2)/T$gamma.m4[intersect(j4,j4.T)], pch = 16, col = 4)
abline(a = 0, b = 1, col = 'gray')

# plot(R.hl[j2], log(2)/T$gamma.m2[j2], pch = 16, col = 2+4, log = 'xy', xlab = 'RNA-seq', ylab = 'Exonarray', main = 'half-lives')
# points(R.hl[intersect(j2,j2.T)], log(2)/T$gamma.m2[intersect(j2,j2.T)], pch = 16, col = 2)
# points(R.hl[j3], log(2)/T$gamma.m3[j3], pch = 16, col = 3+4)
# points(R.hl[intersect(j3,j3.T)], log(2)/T$gamma.m3[intersect(j3,j3.T)], pch = 16, col = 3)
# points(R.hl[j4], log(2)/T$gamma.m4[j4], pch = 16, col = 4+4)
# points(R.hl[intersect(j4,j4.T)], log(2)/T$gamma.m4[intersect(j4,j4.T)], pch = 16, col = 4)
# abline(a = 0, b = 1, col = 'gray')


breaks = seq(0,24,by = 1)

hist(R$phase.gamma.m4[j4], breaks = breaks, col = 4+4, border = 4)
hist(R$phase.gamma.m3[j3], breaks = breaks, col = 3+4, border = 3, add = TRUE)

hist(T$phase.gamma.m4[j4.T], breaks = breaks, col = 4+4, border = 4)
hist(T$phase.gamma.m3[j3.T], breaks = breaks, col = 3+4, border = 3, add = TRUE)
plot(R$phase.gamma.m3[intersect(j3,j3.T)],T$phase.gamma.m3[intersect(j3,j3.T)], col = 3, pch = 16, xlim = c(0,24), ylim = c(0,24))
points(R$phase.gamma.m4[intersect(j4,j4.T)],T$phase.gamma.m4[intersect(j4,j4.T)], col = 4, pch = 16)
text(R$phase.gamma.m3[intersect(j3,j3.T)],T$phase.gamma.m3[intersect(j3,j3.T)],R$gene[intersect(j3,j3.T)], col = 3, pos = 3)
text(R$phase.gamma.m4[intersect(j4,j4.T)],T$phase.gamma.m4[intersect(j4,j4.T)],R$gene[intersect(j4,j4.T)], col = 4, pos = 3)

j34 = intersect(j3,j4.T)
j43 = intersect(j4,j3.T)
points(R$phase.gamma.m3[j34],T$phase.gamma.m4[j34], col = 3+4, pch = 16)
text(R$phase.gamma.m3[j34],T$phase.gamma.m4[j34],R$gene[j34], col = 3+4, pos = 3)
points(R$phase.gamma.m3[j43],T$phase.gamma.m4[j43], col = 4+4, pch = 16)
text(R$phase.gamma.m3[j43],T$phase.gamma.m4[j43],R$gene[j43], col = 4+4, pos = 3)

abline(a = 0, b = 1, col = 'gray')



plot(R$AIC.m2/apply(cbind(R$AIC.m1,R$AIC.m2,R$AIC.m3,R$AIC.m4),1,sum), R$AIC.m3/apply(cbind(R$AIC.m1,R$AIC.m2,R$AIC.m3,R$AIC.m4),1,sum), col = R$best.model, pch = 16)

xaxis = -R$AIC.m2/apply(cbind(R$AIC.m1,R$AIC.m2,R$AIC.m3,R$AIC.m4),1,sum) + R$AIC.m3/apply(cbind(R$AIC.m1,R$AIC.m2,R$AIC.m3,R$AIC.m4),1,sum)
yaxis = R$AIC.m1/apply(cbind(R$AIC.m1,R$AIC.m2,R$AIC.m3,R$AIC.m4),1,sum) - R$AIC.m4/apply(cbind(R$AIC.m1,R$AIC.m2,R$AIC.m3,R$AIC.m4),1,sum)


plot(xaxis, yaxis, col = R$best.model, pch = 16, cex = 0.5, xlim = c(-1,1), ylim = c(-1,1))
abline(a = 0, b = 1, col = 'gray')
abline(a = 0, b = -1, col = 'gray')


AIC.test = rbind(c(1,0,0,0),c(0,1,0,0), c(0,0,1,0), c(0,0,0,1), c(0.25,0.25,0.25,0.25),c(0.1,0.3,0.3,0.3),c(0.3,0.1,0.3,0.3),c(0.3,0.3,0.1,0.3),c(0.3,0.3,0.3,0.1), c(0.3,0.3,0.2,0.2), c(0.25,0.35,0.2,0.2))

xaxis = -AIC.test[,2] + AIC.test[,3]
yaxis = -AIC.test[,4] + AIC.test[,1]

text(xaxis, yaxis, 1:nrow(AIC.test))


T$gene[which((T$LRT.best.model == 3) & (R.bm==3))]

T$gene[which((T$best.model == 3) & (R$best.model==3))]


T$gene[which((T$LRT.best.model == 4) & (R$best.model==4))]

T$gene[which((T$LRT.best.model == 2) & (R$best.model==2))]


plot.RNA.seq.comparison.for.this.gene(gene = 'Tfrc',model = NA)

plot.RNA.seq.comparison.for.this.gene(gene = 'Rbm3',model = NA)

plot.RNA.seq.comparison.for.this.gene(gene = 'Cry1',model = NA)
plot.RNA.seq.comparison.for.this.gene()


load('../complementary_data/Cyclix Low-Seq/Cyclix_Q_masked_pol2rb.Rdata')
load('../complementary_data/Cyclix/CyclixFilterGeneTranscripts.Rdata')





j = which(R.bm==2)
o = order(T$LRT.best.model[j])
pdf('plots/RNA_seq_individual_genes_after_fit_model_2.pdf', width = 12, height = 15)
for(i in j[o]){
	plot.RNA.seq.comparison.for.this.gene(gene = R$gene[i])
}
dev.off()


j = which(R.bm==3)
o = order(T$LRT.best.model[j])
pdf('plots/RNA_seq_individual_genes_after_fit_model_3.pdf', width = 12, height = 15)
for(i in j[o]){
	plot.RNA.seq.comparison.for.this.gene(gene = R$gene[i])
}
dev.off()


j = which(R.bm==4)
o = order(T$LRT.best.model[j])
pdf('plots/RNA_seq_individual_genes_after_fit_model_4.pdf', width = 12, height = 15)
for(i in j[o]){
	plot.RNA.seq.comparison.for.this.gene(gene = R$gene[i])
}
dev.off()


plot.RNA.seq.comparison.for.this.gene = function(gene = 'Loxl4', model = NA){
	
	
	cat(gene, '\n')
	j = which(T$gene == gene)
	
	T$best.model[j]
	T$LRT.best.model[j]
	R$best.model[j]
	R$LRT.best.model[j]
	
	if(is.na(model)){model = R$LRT.best.model[j]}
	
	
	int = grep('.rel.ampl.int', colnames(T))
	ex = grep('.rel.ampl.ex', colnames(T))
	
	i.int = grep('rel.dens.ctrl.premRNA.ZT', colnames(R))
	i.ex = grep('rel.dens.ctrl.mRNA.ZT', colnames(R))
	i.adj.ex = grep('adj.dens.ctrl.mRNA.', colnames(R) )
	
	zt.R = seq(0,44,by = 4)
	zt.T = seq(0,46,by = 2)
	zt.p = seq(0,48, by = 0.1)
	zt.C = seq(2,48,by = 4)
#cyclix data
	
	i.pol2.prom =grep('TxS2K',colnames(all)); i.pol2.prom = i.pol2.prom[1:6]
	i.pol2.body =grep('CBtot.CBlen',colnames(all)); i.pol2.body = i.pol2.body[1:6]
	k1 = which(filter$name == gene)
	k2 = which(rownames(all) == filter$transcript[k1])
	
	lwd = 2
	attribute.global.variable.colors()
	par(mfrow = c(4,2))
	
	
	plot(1,1,type = 'n', xlab = 'time (h)', ylab = 'rel. signal', main = paste(T$gene[j],'intron'), xlim = c(0,48), ylim = range(0,2,R[j,i.int],R[j, i.ex],T[j,int],T[j,ex]))
	abline(h = 1, col = 'gray')
	points(zt.R, R[j,i.int], col = col.int, type = 'b', lwd = lwd)
	points(zt.T, T[j,int], col = col.int, type = 'l', lty = 3, lwd = lwd)
	if(length(k2)>0){k2 = k2[1];
		points(zt.C, rep(all[k2,i.pol2.prom],2)/mean(all[k2,i.pol2.prom]), col = 'red', type = 'l', lty = 2, lwd = lwd)
		points(zt.C, rep(all[k2,i.pol2.body],2)/mean(all[k2,i.pol2.body]), col = 'red3', type = 'l', lty = 2, lwd = lwd/2)
	}
	
	plot(1,1,type = 'n', xlab = 'time (h)', ylab = 'rel. signal', main = paste(T$gene[j],'exon'), xlim = c(0,48), ylim = range(0,2,R[j,i.int],R[j, i.ex],T[j,int],T[j,ex]))
	abline(h = 1, col = 'gray')
	points(zt.R, R[j, i.ex], col = col.ex, type = 'b', lwd = lwd)
	points(zt.T, T[j,ex], col = col.ex, type = 'l', lty = 3, lwd = lwd)
	
	plot(1,1,type = 'n', xlab = 'time (h)', ylab = 'rel. signal', main = paste(T$gene[j],'RNA-seq - model',R.bm[j]), xlim = c(0,48), ylim = range(0,2,R[j,i.int],R[j, i.ex],T[j,int],T[j,ex]))
	abline(h = 1, col = 'gray')
	points(zt.R, R[j,i.int], col = col.int, type = 'b', lwd = lwd)
	points(zt.R, R[j, i.ex], col = col.ex, type = 'b', lwd = lwd)
	
	plot(1,1,type = 'n', xlab = 'time (h)', ylab = 'rel. signal', main = paste(T$gene[j],'Exonarray - model',T$LRT.best.model[j]), xlim = c(0,48), ylim = range(0,2,R[j,i.int],R[j, i.ex],T[j,int],T[j,ex]))
	abline(h = 1, col = 'gray')
	points(zt.T, T[j,int], col = col.int, type = 'l', lty = 3, lwd = lwd)
	points(zt.T, T[j,ex], col = col.ex, type = 'l', lty = 3, lwd = lwd)
	
	
	plot(1,1,type = 'n', xlab = 'time (h)', ylab = 'adj. signal', main = paste(T$gene[j],'RNA-seq'), xlim = c(0,48), ylim = range(0,2,R[j,i.int],R[j, i.ex],T[j,int],T[j,ex],R[j, i.adj.ex] ))
	abline(h = 1, col = 'gray')
	points(zt.R, R[j,i.int], col = col.int, type = 'b', lwd = lwd)
	points(zt.R, R[j, i.adj.ex], col = col.ex, type = 'b', lwd = lwd)
	
	
	
	fold.change = 1
	phase = 0
	up.time = 4
	down.time = 4
	gamma = 2
	eps.gamma = 1
	phase.gamma = 0
	if(model == 1){ gamma = R$gamma.m1[j]}
	if(model == 2){ gamma = R$gamma.m2[j]; fold.change = R$fold.change.int.m2[j]; phase = R$phase.int.m2[j]; up.time = R$up.time.int.m2[j] ; down.time = R$down.time.int.m2[j]}
	if(model == 3){ gamma = R$gamma.m3[j]; eps.gamma = R$eps.gamma.m3[j]; phase.gamma = R$phase.gamma.m3[j] }
	if(model == 4){ gamma = R$gamma.m4[j]; eps.gamma = R$eps.gamma.m4[j]; phase.gamma = R$phase.gamma.m4[j]; fold.change = R$fold.change.int.m4[j]; phase = R$phase.int.m4[j]; up.time = R$up.time.int.m4[j] ; down.time = R$down.time.int.m4[j]}
	
	set.global.variable();syn.f = splicing.rate
	
	
	int.fit = compute.sigmoid(t = zt.p, fold.change = fold.change, phase = phase, up.time = up.time, down.time = down.time)
	ex.fit = compute.m.sigmoid(t = zt.p, gamma = gamma, eps.gamma = eps.gamma , phase.gamma = phase.gamma, fold.change = fold.change, phase = phase, up.time = up.time, down.time =  down.time, rescale = FALSE, synthesis.factor = syn.f )
	
	plot(1,1,type = 'n', xlab = 'time (h)', ylab = 'adj. signal', main = paste(T$gene[j],'RNA-seq + fit - t1/2 = ', round(log(2)/gamma, digits = 2)), xlim = c(0,48), ylim = range(0,2,R[j,i.int],R[j, i.ex],T[j,int],T[j,ex],R[j, i.adj.ex], int.fit, ex.fit ))
	abline(h = 1, col = 'gray')
	points(zt.R, R[j,i.int], col = col.int, type = 'b', lwd = lwd)
	points(zt.R, R[j, i.adj.ex], col = col.ex, type = 'b', lwd = lwd)
	points(zt.p, int.fit, col = col.int, type = 'l', lwd = lwd)
	points(zt.p, ex.fit+int.fit, col = col.ex, type = 'l', lwd = lwd)
	
	plot(1,1,type = 'n', xlab = 'time (h)', ylab = 'adj. signal (log2)', main = paste(T$gene[j],'RNA-seq + fit'), xlim = c(0,48), ylim = log2(range(2,R[j,i.int],R[j, i.ex],T[j,int],T[j,ex],R[j, i.adj.ex], int.fit, ex.fit )))
	abline(h = 0, col = 'gray')
	points(zt.R, log2(R[j,i.int]), col = col.int, type = 'b', lwd = lwd)
	points(zt.R, log2(R[j, i.adj.ex]), col = col.ex, type = 'b', lwd = lwd)
	points(zt.p, log2(int.fit), col = col.int, type = 'l', lwd = lwd)
	points(zt.p, log2(ex.fit+int.fit), col = col.ex, type = 'l', lwd = lwd)
	
	sum.ex = ex.fit+int.fit
	
	plot(1,1,type = 'n', xlab = 'time (h)', ylab = 'rel. signal', main = paste(T$gene[j],'RNA-seq + fit'), xlim = c(0,48), ylim = range(0,2,R[j,i.int],R[j, i.ex],T[j,int],T[j,ex], int.fit/mean(int.fit), sum.ex/mean(sum.ex) ))
	abline(h = 1, col = 'gray')
	points(zt.R, R[j,i.int], col = col.int, type = 'b', lwd = lwd)
	points(zt.R, R[j, i.ex], col = col.ex, type = 'b', lwd = lwd)
	points(zt.p, int.fit/mean(int.fit), col = col.int, type = 'l', lwd = 1)
	points(zt.p, sum.ex/mean(sum.ex), col = col.ex, type = 'l', lwd = 1)
	
	
# param = make.fits.with.all.models.for.one.gene(T = R, gene.index = j, debug = TRUE, parametrization = 'sigmoid', sum.species = TRUE, zt = zt.R, i.ex = i.adj.ex, i.int = i.int, absolute.signal = TRUE)
# #source('functions.R')
# if(model == 1){ gamma = param['gamma.m1']}
# if(model == 2){ gamma = param['gamma.m2']; fold.change = param['fold.change.int.m2']; phase = param['phase.int.m2']; up.time = param['up.time.int.m2 ']; down.time = param['down.time.int.m2']}
# if(model == 3){ gamma = param['gamma.m3']; eps.gamma = param['eps.gamma.m3']; phase.gamma = param['phase.gamma.m3'] }
# if(model == 4){ gamma = param['gamma.m4']; eps.gamma = param['eps.gamma.m4']; phase.gamma = param['phase.gamma.m4']; fold.change = param['fold.change.int.m4']; phase = param['phase.int.m4']; up.time = param['up.time.int.m4']; down.time = param['down.time.int.m4']}
	
# int.fit = compute.sigmoid(t = zt.p, fold.change = fold.change, phase = phase, up.time = up.time, down.time = down.time)
# ex.fit = compute.m.sigmoid(t = zt.p, gamma = gamma, eps.gamma = eps.gamma , phase.gamma = phase.gamma, fold.change = fold.change, phase = phase, up.time = up.time, down.time =  down.time, rescale = FALSE, synthesis.factor = syn.f )
# points(zt.p, int.fit, col = 'red', type = 'l', lwd = lwd)
# points(zt.p, ex.fit+int.fit, col = 'red', type = 'l', lwd = lwd)
	
# fold.change = 2.32424651
# phase = 20.50228770
# up.time = 4.01420724
# down.time = 10.90872393
# gamma = 0.3
# eps.gamma = 5
# phase.gamma = 10
	
# int.fit = compute.sigmoid(t = zt.p, fold.change = fold.change, phase = phase, up.time = up.time, down.time = down.time)
# ex.fit = compute.m.sigmoid(t = zt.p, gamma = gamma, eps.gamma = eps.gamma , phase.gamma = phase.gamma, fold.change = fold.change, phase = phase, up.time = up.time, down.time =  down.time, rescale = FALSE, synthesis.factor = syn.f )
# points(zt.p, int.fit, col = 'cyan', type = 'l', lwd = lwd)
# points(zt.p, ex.fit+int.fit, col = 'cyan', type = 'l', lwd = lwd)
	
# plot(1,1,type = 'n', xlab = 'time (h)', ylab = 'adj. signal', main = paste(T$gene[j],'RNA-seq + fit'), xlim = c(0,48), ylim = log2(range(2,R[j,i.int],R[j, i.ex],T[j,int],T[j,ex],R[j, i.adj.ex], int.fit, ex.fit )))
# abline(h = 0, col = 'gray')
# points(zt.R, log2(R[j,i.int]), col = col.int, type = 'b', lwd = lwd)
# points(zt.R, log2(R[j, i.adj.ex]), col = col.ex, type = 'b', lwd = lwd)
# points(zt.p, log2(int.fit), col = col.int, type = 'l', lwd = lwd)
# points(zt.p, log2(ex.fit+int.fit), col = col.ex, type = 'l', lwd = lwd)
	
	
# s = compute.sigmoid(t = zt.R, fold.change = fold.change, phase = phase, up.time = up.time, down.time = down.time)
# m = compute.m.sigmoid(t = zt.R, gamma = gamma, eps.gamma = eps.gamma , phase.gamma = phase.gamma, fold.change = fold.change, phase = phase, up.time = up.time, down.time =  down.time, rescale = FALSE, synthesis.factor = syn.f )
# points(zt.R, s, col = 'cyan', type = 'l', lwd = lwd)
# points(zt.R, s+m, col = 'cyan', type = 'l', lwd = lwd)
	
	
# par = c(gamma, eps.gamma, phase.gamma, fold.change, phase, up.time, down.time)
# f2min(par = par, M =R[j, i.adj.ex], S = R[j,i.int], model = model, parametrization = 'sigmoid', method.intern = 'integration', debug = TRUE, sum.species = TRUE, zt = zt.R, absolute.signal = TRUE )
# errors(S = R[j,i.int],s = s,M = R[j, i.adj.ex],m =  s+m ,param = par, model = model,parametrization = 'sigmoid', log = TRUE)
# error(S =R[j,i.int], s = s, M =R[j, i.adj.ex], m =s+m, log  = TRUE ) 
	
}


hist(R$phase.int.m2[j2], breaks = breaks, col = 2+4, border = 2)
hist(T$phase.int.m2[j2.T], breaks = breaks, col = 2+4, border = 2)
plot(R$phase.int.m2[intersect(j2,j2.T)],T$phase.int.m2[intersect(j2,j2.T)], col = 2, pch = 16)
abline(a = 0, b = 1, col = 'gray')

hist(R$phase.int.m4[j4], breaks = breaks, col = 4+4, border = 4)
hist(T$phase.int.m4[j4.T], breaks = breaks, col = 4+4, border = 4)
plot(R$phase.int.m4[intersect(j4,j4.T)],T$phase.int.m4[intersect(j4,j4.T)], col = 4, pch = 16)
text(R$phase.int.m4[intersect(j4,j4.T)],T$phase.int.m4[intersect(j4,j4.T)],R$gene[intersect(j4,j4.T)], col = 4, pos = 3)
abline(a = 0, b = 1, col = 'gray')



plot(T$phase.int.m2[j2.T],T$fold.change.int.m2[j2.T],pch = 16, col = 2, log = 'y')
plot(T$phase.int.m4[j4.T],T$fold.change.int.m4[j4.T],pch = 16, col = 4, log = 'y')

plot(T$phase.gamma.m3[j3.T],T$eps.gamma.m3[j3.T],pch = 16, col = 3, log = 'y')






hist(T$phase.int.m2[j2.T], breaks = seq(0,24,by = 0.1), col = 2+4, border = 2)
j = which((T$phase.int.m2[j2.T] == 0)|(T$phase.int.m2[j2.T] == 24))
j = j2.T[j]
plot.RNA.seq.comparison.for.this.gene(gene = R$gene[j[1]])


T[j,]

plot(T$phase.int.m2[j], T$phase.int.m4[j], pch = 16, col = 2)



i = 1
#param = make.fits.with.all.models.for.one.gene(T = T, gene.index = j[i], debug = TRUE, parametrization = 'sigmoid', sum.species = TRUE, zt = zt.T, i.ex = grep('.rel.ampl.ex', colnames(T)), i.int = grep('.rel.ampl.int', colnames(T)), absolute.signal = FALSE)
#param = make.fit.spec.model(T = T, gene.index = j[i],model = 2, debug = TRUE, parametrization = 'sigmoid', method = 'integration', sum.species = TRUE, zt = zt.T, i.ex = grep('.rel.ampl.ex', colnames(T)), i.int = grep('.rel.ampl.int', colnames(T)), absolute.signal = FALSE)
param = make.optimization(T = T,i = j[i], model = 2, Nfit = 8, debug = TRUE, parametrization = 'sigmoid', method = 'integration', sum.species = TRUE, zt = zt.T, i.ex = grep('.rel.ampl.ex', colnames(T)), i.int = grep('.rel.ampl.int', colnames(T)), absolute.signal = FALSE)

param = make.optimization(T = T,i = 4, model = 2, Nfit = 10, debug = TRUE, parametrization = 'sigmoid', method = 'integration', sum.species = TRUE, zt = zt.T, i.ex = grep('.rel.ampl.ex', colnames(T)), i.int = grep('.rel.ampl.int', colnames(T)), absolute.signal = FALSE)


i = 36
param = make.optimization(T = T,i = i, model = 2, Nfit = 5, debug = TRUE, parametrization = 'sigmoid', method = 'integration', sum.species = TRUE, zt = zt.T, i.ex = grep('.rel.ampl.ex', colnames(T)), i.int = grep('.rel.ampl.int', colnames(T)), absolute.signal = FALSE)
param = make.optimization(T = T,i = i, model = 2, Nfit = NA, debug = TRUE, parametrization = 'sigmoid', method = 'integration', sum.species = TRUE, zt = zt.T, i.ex = grep('.rel.ampl.ex', colnames(T)), i.int = grep('.rel.ampl.int', colnames(T)), absolute.signal = FALSE)


T[i,74:92]
}



################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis
################################################################################################################################################################## finishing line of analysis


########## COMPARISON WITH H-L from Rosie Pilot Exp

if(plot){
	load(paste('Rdata/genes_fits_results',version,'.Rdata',sep = ''))
	source('functions.R')
	R = read.table('../complementary_data/rosiepilot/decay.txt', col.names = c('probeset',paste('inter',c('est','std.error','t.val','pr'),sep = '.'),paste('slope',c('est','std.error','t.val','pr'),sep = '.'), 'gene'))
	m = match(T$gene, R$gene)
	R = R[m,]
	degr.rates = -log(2)*R$slope.est*60
	degr.rates = pmax(degr.rates,log(2)/100)
	half.lives = log(2)/degr.rates

	best.model = T$LRT.best.model
	gamma = rep(1,nrow(T));
	gamma[which(best.model == 2)] = T$gamma.m2[which(best.model == 2)]
	gamma[which(best.model == 3)] = T$gamma.m3[which(best.model == 3)]
	gamma[which(best.model == 4)] = T$gamma.m4[which(best.model == 4)]
	
	plot(gamma[best.model>1], degr.rates[best.model>1], pch = 16, col = best.model[best.model>1], log = 'xy')
	abline(a = 0, b = 1)
	
	plot(log(2)/gamma[best.model>1], half.lives[best.model>1], pch = 16, col = best.model[best.model>1], log = 'xy')
	abline(a = 0, b = 1)
}


if(plot){
	cat('plot : FITS RESULTS ON REAL DATA : COMPARISON OPTIM AND MCMC\n')
	if(load.for.plot){load(file = "Rdata/genes_fits_results.Rdata")}

	source('functions.R')
	
	pdf('plots/fits_results_comparison_optim_and_MCMC.pdf',width = 6.6, height = 12)
	par(mfcol = c(4,2), cex = 0.6, mar = c(3,3,2,0.2)+0.1, mgp = c(1.6,0.5,0),las = 1, tcl = -0.3)

	cats = paste(T$LRT.best.model,T$model.MCMC)
	rle = rle(sort(cats))
	m = match( cats, rle$values)
	plot(T$LRT.best.model, T$model.MCMC, pch = 21, bg = T$LRT.best.model, col =  T$model.MCMC, cex =rle$length[m]/20, axes = FALSE, lwd = 2, xlim = c(0.8,4.2), ylim = c(0.8,4.2), xlab = 'Simple Optim method', ylab = 'MCMC method')
	axis(1, at = 1:4, labels = c('CS-CD','RS-CD','CS-RD','RS-RD'), lwd = 0)
	axis(2, at = 1:4, labels = c('CS-CD','RS-CD','CS-RD','RS-RD'), lwd = 0, las = 3)
	legend('bottomright',legend = c('MCMC prediction','Optim prediction') , pch = c(1,16))
	
	
	par(mfrow = c(1,2))
	sd = 0.12
	plot(as.numeric(T$best.model)+rnorm(nrow(T), sd = sd) , T$model.MCMC+rnorm(nrow(T), sd = sd), pch = 18, col = rgb(0.5,0.5,0.7,0.5), cex = 0.5, xlab = 'AIC', ylab = 'MCMC')
	abline(h = seq(1.5,3.5), col = 'gray'); abline(v = seq(1.5,3.5), col = 'gray')
	plot(as.numeric(T$LRT.best.model)+rnorm(nrow(T), sd = sd) , T$model.MCMC+rnorm(nrow(T), sd = sd), pch = 18, col = rgb(0.5,0.5,0.7,0.5), cex = 0.5, xlab = 'LRT', ylab = 'MCMC')
	abline(h = seq(1.5,3.5), col = 'gray'); abline(v = seq(1.5,3.5), col = 'gray')

#	ind = (T$LRT.best.model == 2)|(T$model.MCMC == 2)
#	plot(T$gamma.m2[ind], T$gamma.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind], ylim = range(T$gamma.m2[ind], na.rm = TRUE), main = 'degradation rate', xlab = 'Simple Optim method', ylab = 'MCMC method')
#	ind = (T$LRT.best.model == 3)|(T$model.MCMC == 3)
#	points(T$gamma.m3[ind], T$gamma.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind])
#	ind = (T$LRT.best.model == 4)|(T$model.MCMC == 4)
#	points(T$gamma.m4[ind], T$gamma.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind])
#	abline(a = 0, b = 1)
	
	ind = (T$LRT.best.model == 2)|(T$model.MCMC == 2)
	plot(log(2)/T$gamma.m2[ind], log(2)/T$gamma.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind], ylim = range(log(2)/T$gamma.m2[ind], na.rm = TRUE), main = 'Half-lives', xlab = 'Simple Optim method', ylab = 'MCMC method')
	ind = (T$LRT.best.model == 3)|(T$model.MCMC == 3)
	points(log(2)/T$gamma.m3[ind], log(2)/T$gamma.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind], ylim = range(T$gamma.m2[ind], na.rm = TRUE))
	ind = (T$LRT.best.model == 4)|(T$model.MCMC == 4)
	points(log(2)/T$gamma.m4[ind], log(2)/T$gamma.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind], ylim = range(T$gamma.m2[ind], na.rm = TRUE))
	abline(a = 0, b = 1)
	
	ind = (T$LRT.best.model == 3)|(T$model.MCMC == 3)
	plot(T$eps.gamma.m3[ind], T$eps.gamma.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind], ylim = c(0,1), xlim = c(0,1), main = 'degradation amplitude', xlab = 'Simple Optim method', ylab = 'MCMC method')
	ind = (T$LRT.best.model == 4)|(T$model.MCMC == 4)
	points(T$eps.gamma.m4[ind], T$eps.gamma.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind])
	abline(a = 0, b = 1)
	
	ind = (T$LRT.best.model == 3)|(T$model.MCMC == 3)
	plot(T$phase.gamma.m3[ind], T$phase.gamma.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind], ylim = c(0,24), xlim = c(0,24), main = 'degradation phase', xlab = 'Simple Optim method', ylab = 'MCMC method')
	ind = (T$LRT.best.model == 4)|(T$model.MCMC == 4)
	points(T$phase.gamma.m4[ind], T$phase.gamma.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind])
	abline(a = 0, b = 1)
	
	
		
	ind = (T$LRT.best.model == 2)|(T$model.MCMC == 2)
	plot(T$rel.ampl.int.m2[ind], T$rel.ampl.int.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind], ylim = c(0,1),xlim = c(0,1), main = 'Relative Amplitude Intron', xlab = 'Simple Optim method', ylab = 'MCMC method')
	ind = (T$LRT.best.model == 4)|(T$model.MCMC == 4)
	points(T$rel.ampl.int.m4[ind], T$rel.ampl.int.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind], ylim = range(T$gamma.m2[ind], na.rm = TRUE))
	abline(a = 0, b = 1)

		
		
		
	ind = (T$LRT.best.model == 2)|(T$model.MCMC == 2)
	plot(T$phase.int.m2[ind], T$phase.int.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind], ylim = c(0,24), xlim = c(0,24), main = 'Phase Intron', xlab = 'Simple Optim method', ylab = 'MCMC method')
	ind = (T$LRT.best.model == 4)|(T$model.MCMC == 4)
	points(T$phase.int.m4[ind], T$phase.int.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind], ylim = range(T$gamma.m2[ind], na.rm = TRUE))
	abline(a = 0, b = 1)	
			
	ind = (T$LRT.best.model == 2)|(T$model.MCMC == 2)
	plot(T$rel.ampl.12.int.m2[ind], T$rel.ampl.12.int.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind], ylim = c(0,1),xlim = c(0,1), main = 'Relative Amplitude Intron (12h)', xlab = 'Simple Optim method', ylab = 'MCMC method')
	ind = (T$LRT.best.model == 4)|(T$model.MCMC == 4)
	points(T$rel.ampl.12.int.m4[ind], T$rel.ampl.12.int.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind], ylim = range(T$gamma.m2[ind], na.rm = TRUE))
	abline(a = 0, b = 1)

	ind = (T$LRT.best.model == 2)|(T$model.MCMC == 2)
	plot(T$phase.12.int.m2[ind], T$phase.12.int.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind], ylim = c(0,12), xlim = c(0,12), main = 'Phase Intron (12h)', xlab = 'Simple Optim method', ylab = 'MCMC method')
	ind = (T$LRT.best.model == 4)|(T$model.MCMC == 4)
	points(T$phase.12.int.m4[ind], T$phase.12.int.MCMC[ind], pch = 21,bg = T$LRT.best.model[ind],  col = T$model.MCMC[ind], ylim = range(T$gamma.m2[ind], na.rm = TRUE))
	abline(a = 0, b = 1)	
	
	
	
	dev.off()
	
	
}



if(plot){
	cat('plot : FITS RESULTS ON REAL DATA\n')
	if(load.for.plot){
		#if(gene.set == gene.sets[1]){load(file = "Rdata/gene_fit_results_example_genes_and_LRT.Rdata")}
		#if(gene.set == gene.sets[2]){load(file = "Rdata/gene_fit_results_for_rhythmic_genes_and_LRT.Rdata")}
		load(file = "Rdata/genes_fit_results_and_LRT.Rdata")
	}
	ZT.int = grep('.rel.ampl.int', colnames(T)); ZT.ex = grep('.rel.ampl.ex', colnames(T)); 
	source('functions.R')
	best.model = as.numeric(T$best.model)
		
	#if(gene.set == gene.sets[1]){pdf('plots/fits_results_example_genes.pdf', width = 6.6, height = 4)}
	#if(gene.set == gene.sets[2]){pdf('plots/fits_results_rhythmic_genes.pdf', width = 6.6, height = 4)}
	
	pdf('plots/fits_results_sum_species.pdf', width = 6.6, height = 4)
	par(mfrow = c(2,3), cex = 0.6, mar = c(3,3,2,0.2)+0.1, mgp = c(1.6,0.5,0),las = 1, tcl = -0.3)
	
#	var = apply(T[,c(ZT.int, ZT.ex)],1,var)
#	errors =  as.matrix(T[,grep('error.m',colnames(T))]); error =errors[(best.model-1)*nrow(T)+(1:nrow(T))]; error = error/length(c(ZT.int, ZT.ex))
#	perc.of.variance.explained = 100*(var-error)/var; perc.of.variance.explained  = apply(cbind(perc.of.variance.explained,-5),1,max)
#	o = order(perc.of.variance.explained)
#	plot(perc.of.variance.explained[o], ylim =c(-5,100), pch = '_', col = best.model[o], axes = FALSE, main = '% of var. explained by the fits', ylab = '', xlab = '')
#	abline(h = 0); axis(2); box();

	h = hist(best.model,breaks = seq(0.5,4.5,by = 1), plot = FALSE)
	pie(h$counts, clockwise = FALSE, col = 1:4, border = FALSE, labels = c('CS-CD','RS-CD','CS-RD','RS-RD'))


	fractions.ex = matrix(NA, nrow = 14, ncol = 4); 
	numbers.ex = rep(0,14)
	fractions.int = matrix(NA, nrow = 14, ncol = 4); 
	numbers.int = rep(0,14)
	for(r in 1:14){
		ok = (log2(T$exon.median)>(r-1))&(log2(T$exon.median)<=r)
		fractions.ex[r,] =  hist( T$LRT.best.model[ok], breaks = seq(0.5,4.5,by=1),plot = FALSE)$density
		numbers.ex[r] = sum(ok, na.rm = TRUE)
		ok = (log2(T$intron.median)>(r-1))&(log2(T$intron.median)<=r)
		fractions.int[r,] =  hist( T$LRT.best.model[ok], breaks = seq(0.5,4.5,by=1),plot = FALSE)$density
		numbers.int[r] = sum(ok, na.rm = TRUE)
		}
		
	par(mgp = c(0,0.5,0), mar = c(1,3,2,0.2)+0.1)
	
	b = barplot(t(fractions.ex), col = c(1:4), border = FALSE, ylim = c(-0.15,1.15), space = 0.05, xlab = 'exonic probes expression level')
	text(b,rep(1.005,14), numbers.ex, pos = 4, cex = 0.7, srt = 90, offset = 0.1)
	text(c(b[1]-diff(b)[1],b)+diff(b)[1]/2,rep(0,15), 0:14, pos = 1)
	
	b = barplot(t(fractions.int), col = c(1:4), border = FALSE,  ylim = c(-0.15,1.15), space = 0.05, xlab = 'intronic probes expression level')
	text(b,rep(1.005,14), numbers.int, pos = 4, cex = 0.7, srt = 90, offset = 0.1)
	text(c(b[1]-diff(b)[1],b)+diff(b)[1]/2,rep(0,15), 0:14, pos = 1)
	

	par(mar = c(3,3,2,0.2)+0.1, mgp = c(1.6,0.5,0))


	bounds = set.bounds(model = 4, parametrization = 'sigmoid')

	# half-lives 
	breaks = lseq(log(2)/bounds$upper[1],log(2)/bounds$lower[1], len = 31)
	
	
	h = hist(log(2)/T$gamma.m2[best.model ==2], breaks = breaks, col = 2+4, border = 2, xlab = '[h]', ylab = '', main = 'Half-lives', plot = FALSE)
	plot(c(h$breaks[1],h$breaks,h$breaks[length(h$breaks)]),c(0,h$counts,h$counts[length(h$counts)],0), log = 'x', type = 'n', main = 'Half-lives',  xlab = 'half-lives [h]', ylab = 'counts')
	polygon(c(rep(h$breaks, each = 2)), c(0, rep(h$counts, each=2),0), col = 2+4, border = 2, lwd = 1)
	
	h = hist(log(2)/T$gamma.m3[best.model ==3], breaks = breaks, col = 3+4, border = 3, add = TRUE, freq = FALSE, plot = FALSE)
	polygon(c(rep(h$breaks, each = 2)), c(0, rep(h$counts, each=2),0), col = 3+4, border = 3, lwd = 1)
	
	
	h = hist(log(2)/T$gamma.m4[best.model ==4], breaks = breaks, col = 4+4, border = 4, add = TRUE, freq = FALSE , plot = FALSE)
	polygon(c(rep(h$breaks, each = 2)), c(0, rep(h$counts, each=2),0), col = 4+4, border = 4, lwd = 1)

	
	# amplitude of degradation 
	breaks = lseq(bounds$lower[2],bounds$upper[2], len = 51)
	
	h = hist(T$eps.gamma.m3[best.model ==3], breaks = breaks,  plot = FALSE)
	h4 = hist(T$eps.gamma.m4[best.model ==4], breaks = breaks,  plot = FALSE)

	plot(c(h$breaks[1],h$breaks,h$breaks[length(h$breaks)]),c(0,h$counts,h$counts[length(h$counts)],0), log = 'x', type = 'n', main = 'Fold changes of degradation',  xlab = 'Fold change', ylab = 'counts', ylim = range(h$counts, h4$counts))
	polygon(c(rep(h$breaks, each = 2)), c(0, rep(h$counts, each=2),0), col = 3+4, border = 3, lwd = 1)

	#plot(c(h$breaks[1],h$breaks,h$breaks[length(h$breaks)]),c(0,h$counts,h$counts[length(h$counts)],0), log = 'x', type = 'n', main = 'Fold changes of degradation',  xlab = 'Fold change', ylab = 'counts')
	polygon(c(rep(h4$breaks, each = 2)), c(0, rep(h4$counts, each=2),0), col = 4+4, border = 4, lwd = 1)


	# phase of maximal degradation
	breaks = seq(bounds$upper[3],bounds$lower[3], len = 25)
	h = hist(T$phase.gamma.m3[best.model ==3], breaks = breaks, col = 3+4, border = 3, xlab = '[h]', ylab = '', main = 'Phase of maximal degradation', axes = FALSE)
	h = hist(T$phase.gamma.m4[best.model ==4], breaks = breaks, col = 4+4, border = 4, add = TRUE)
	axis(1, at = seq(0,24,by = 4)); axis(2)
	
	
	
	
	dev.off()
	
	
	# p-values? how to select model 4? with combined Pvalues?

	x =  sign(T$LRT.pval.p-T$LRT.pval.d)*(log10(apply(cbind(T$LRT.pval.d,T$LRT.pval.p),1,min)))
	y = -log10(apply(cbind(T$LRT.pval.P,T$LRT.pval.D),1,max))
	LRT.combined = pchisq(-2*log(T$LRT.pval.P * T$LRT.pval.D), df = 4, lower.tail = FALSE)
	y = -log10(LRT.combined)
	y = -log10(T$LRT.pval.A)

	plot(x,y, pch = 16, col = T$LRT.best.model, cex=  0.5)#, ylim = c(0,10))
	polygon(x = c(log10(0.05),-log10(0.05),-log10(0.05),log10(0.05)), y = c(0,0,-log10(0.05),-log10(0.05)), lty = 3)
	j = (x<= -5)|(x>7)|(y>2)
	text(x[j], y[j], T$gene[j], col = T$LRT.best.model[j], pos = 3, cex = 0.6, offset = 0.1)
	
	
	library(fdrtool)
	fdr = fdrtool(T$LRT.pval.p[!is.na(T$LRT.pval.p)], statistic = 'pvalue')
	
	
	plot(T$LRT.pval.P, LRT.combined, log = 'xy')
	plot(T$LRT.pval.D, LRT.combined, log = 'xy')

	
	# rejection?
	
	hist(T$exon.var, breaks = 100)
	plot(T$exon.var, T$LRT.best.model, pch = 16, cex = 0.5)
	boxplot(T$exon.var~T$LRT.best.model, pch = 16, cex = 0.5)
	
	plot(T$exon.var[T$LRT.best.model == 2], (T$fold.change.int.m2[T$LRT.best.model == 2]-1)/(T$fold.change.int.m2[T$LRT.best.model == 2]+1), pch = 16)
	abline(a = 0, b = 1)
	text(T$exon.var[T$LRT.best.model == 2], (T$fold.change.int.m2[T$LRT.best.model == 2]-1)/(T$fold.change.int.m2[T$LRT.best.model == 2]+1), T$gene[T$LRT.best.model == 2])
	
	plot(T$exon.var[T$LRT.best.model == 4], (T$fold.change.int.m4[T$LRT.best.model == 4]-1)/(T$fold.change.int.m4[T$LRT.best.model == 4]+1), pch = 16)
	abline(a = 0, b = 1)
	text(T$exon.var[T$LRT.best.model == 4], (T$fold.change.int.m4[T$LRT.best.model == 4]-1)/(T$fold.change.int.m4[T$LRT.best.model == 4]+1), T$gene[T$LRT.best.model == 4], pos = 3)
	
	
	# percentage of variance explained
	
	var.ex = apply(T[,ZT.ex], 1, var)
	var.ex = apply(log(T[,ZT.ex]),1,var)
	var.int = apply(T[,ZT.int], 1, var)
	var.int = apply(log(T[, ZT.int]),1,var)
	VAR = var.ex+var.int
	
	VAR = T$error.m1/24
	
	plot(VAR, T$error.m1/24, pch = 21, col = 'lightblue3', bg = 'lightblue', cex = 0.5)
	abline(a = 0, b = 1, col = 'gray')
	
	best.model = T$LRT.best.model

	j2 = (best.model == 2)
	j3 = (best.model == 3)
	j4 = (best.model == 4)
	
	plot(VAR[j2], T$error.m2[j2]/24, pch = 21, col = 'lightblue3', bg = 'lightblue', cex = 0.5)
	abline(a = 0, b = 1, col = 'gray')
	
	par(mfcol = c(3,2))
	
	perc.m2 = (VAR-T$error.m2/24)/VAR; th = 0.33; n = sum(perc.m2[j2]<th, na.rm = TRUE); N = sum(j2,na.rm = TRUE)
	plot(T$LRT.pval.p, perc.m2, pch = 21, col = 'gray80', bg = 'gray90', cex = 0.5, log = 'x', ylim = range(perc.m2,1, na.rm = TRUE), main = paste(n,'/',N,'=',100*round(n/N,digits = 4),'%'))
	points(T$LRT.pval.p[j2], perc.m2[j2], pch = 16, col = 2, cex = 0.5 )
	points(T$LRT.pval.p[j3], perc.m2[j3], pch = 16, col = 3, cex = 0.5 )
	points(T$LRT.pval.p[j4], perc.m2[j4], pch = 16, col = 4, cex = 0.5 )
	abline(v = 0.05); abline(h = th, col = 'black', lty = 3)
	axis(2, th, las = 1)
	
	
	
	perc.m3 = (VAR-T$error.m3/24)/VAR; th = th; n = sum(perc.m3[j3]<th, na.rm = TRUE); N = sum(j3,na.rm = TRUE)
	plot(T$LRT.pval.d, perc.m3, pch = 21, col = 'gray80', bg = 'gray90', cex = 0.5, log = 'x', ylim = range(perc.m3,1, na.rm = TRUE), main = paste(n,'/',N,'=',100*round(n/N,digits = 4),'%'))
	points(T$LRT.pval.d[j2], perc.m3[j2], pch = 16, col = 2, cex = 0.5 )
	points(T$LRT.pval.d[j3], perc.m3[j3], pch = 16, col = 3, cex = 0.5 )
	points(T$LRT.pval.d[j4], perc.m3[j4], pch = 16, col = 4, cex = 0.5 )
	abline(v = 0.05); abline(h = th, col = 'black', lty = 3)
	axis(2, th, las = 1)

	perc.m4 = (VAR-T$error.m4/24)/VAR; th = th; n = sum(perc.m4[j4]<th, na.rm = TRUE); N = sum(j4,na.rm = TRUE)
	#plot(apply(cbind(T$LRT.pval.P, T$LRT.pval.D),1,max), perc.m4, pch = 21, col = 'gray80', bg = 'gray90', cex = 0.5, log = 'x', ylim = range(perc.m4,1, na.rm = TRUE), main = paste(n,'/',N,'=',100*round(n/N,digits = 4),'%'))
	plot(T$LRT.pval.A, perc.m4, pch = 21, col = 'gray80', bg = 'gray90', cex = 0.5, log = 'x', ylim = range(perc.m4,1, na.rm = TRUE), main = paste(n,'/',N,'=',100*round(n/N,digits = 4),'%'))
	points(T$LRT.pval.A[j2], perc.m4[j2], pch = 16, col = 2, cex = 0.5 )
	points(T$LRT.pval.A[j3], perc.m4[j3], pch = 16, col = 3, cex = 0.5 )
	points(T$LRT.pval.A[j4], perc.m4[j4], pch = 16, col = 4, cex = 0.5 )
	abline(v = 0.05); abline(h = th, col = 'black', lty = 3)
	axis(2, th, las = 1)
	
	
	# transition between the 2 days?
	
	par(mfcol = c(2,2))
	breaks = seq(0,24,by = 0.5)
	lim = 4
	lim.sharp = 2.5

	sharp = (T$up.time.int.m2<3)|(T$down.time.int.m2<3)
	peak = (T$up.time.int.m2<4)&(T$down.time.int.m2<4)
	ratio =  T$up.time.int.m2/T$down.time.int.m2
	
	plot(T$phase.int.m2, ratio, cex = 0.1, pch = 21, col = 'gray80', bg = 'gray90', log = 'y')
	points(T$phase.int.m2[j2], ratio[j2], cex = 0.5, pch = 16, col = 2)
	points(T$phase.int.m2[j2&sharp], ratio[j2&sharp], cex = 0.5, pch = 16, col = 'red')
	points(T$phase.int.m2[j2& peak], ratio[j2& peak], cex = 0.5, pch = 16, col = 'blue')
	text(T$phase.int.m2[j2& peak], ratio[j2& peak],T$gene[j2&peak], cex = 0.5,pos = 3, col = 'blue')

	abline(h = 1); abline(h = c(lim,1/lim), lty = 3)
	abline(v = 22)
	hist(T$phase.int.m2[j2], breaks = breaks)
	hist(T$phase.int.m2[j2][(ratio[j2]> lim)|(ratio[j2]<(1/lim))], breaks = breaks, add = TRUE, col = 'gray90')
	hist(T$phase.int.m2[j2&!sharp], breaks = breaks, add = TRUE, col = rgb(1,0,0,0.3))



	plot(T$phase.int, ratio, cex = 0.3, pch = 21, col = 'gray80', bg = 'gray90', log = 'y')
	abline(h = 1); abline(h = c(lim,1/lim), lty = 3)
	points(T$phase.int[j2], ratio[j2], cex = 0.5, pch = 16, col = 2)
	points(T$phase.int[j2&sharp], ratio[j2&sharp], cex = 0.5, pch = 16, col = 'red')

	hist(T$phase.int[j2], breaks = breaks)
	hist(T$phase.int[j2][(ratio[j2]> lim)|(ratio[j2]<(1/lim))], breaks = breaks, add = TRUE, col = 'gray90')
	hist(T$phase.int[j2&!sharp], breaks = breaks, add = TRUE, col = rgb(1,0,0,0.3))







	plot(perc.m2[j2],ratio[j2],cex = 0.5, pch = 21, col = 'gray80', bg = 'gray90',log = 'y')
	abline(v = 0.4)
	points(perc.m2[j2&sharp],ratio[j2&sharp],cex = 0.5, pch = 16, col = 'red')


	plot(T$phase.int.m2[j2], perc.m2[j2],cex = 0.5, pch = 21, col = 'gray80', bg = 'gray90')
	hist(T$phase.int.m2[j2], breaks = breaks)
	hist(T$phase.int.m2[j2&(perc.m2>0.5)], breaks = breaks, add = TRUE, col = rgb(0.9,0.5,0,0.2))
	
	cols = rainbow(n = 220,s = 0.8, v = 0.8, start = 0.5)
	plot(T$phase.int.m2[j2], T$up.time.int.m2[j2],cex = 1, pch = 16, col = cols[round(10*T$down.time.int.m2[j2])])
	
	
	
	
	
	
	
	
	par(mfcol = c(1,2))
	
	T.tot = T$up.time.int.m2+T$down.time.int.m2
	ratio = T$up.time.int.m2/(T.tot)
	lim.perc = 0.4

	Ncol = 10
	colorlow = 'steelblue';  colormed = 'gray95' ; colorhigh  = 'tomato'
	col.ramp <- colorRampPalette(c(colorlow, colormed, colorhigh), space = "Lab")
	cols = col.ramp(Ncol) 
	
	plot(T$phase.int.m2[j2], T.tot[j2],cex = 1, pch = 16, col = cols[round(10*ratio[j2])])
	points(T$phase.int.m2[j2&(perc.m2> lim.perc)], T.tot[j2&(perc.m2> lim.perc)],cex = 1, pch = 21,col = 'black', bg = cols[round(10*ratio[j2])])
	
	Ncol = 22
	colorlow = 'tomato';  colormed = 'gray95' ; colorhigh  = 'green3'
	col.ramp <- colorRampPalette(c(colorlow, colormed, colorhigh), space = "Lab")
	cols = col.ramp(Ncol) 
	
	plot(T$phase.int.m2[j2], ratio[j2],cex = 1, pch = 16, col = cols[round(T.tot[j2])])
	points(T$phase.int.m2[j2&(perc.m2> lim.perc)], ratio[j2&(perc.m2> lim.perc)],cex = 1, pch = 21,col = 'black', bg = cols[round(T.tot[j2])])
	

	boxplot(T.tot[j2]~round(T$phase.int.m2[j2]/2))
	boxplot(T.tot[j2&(perc.m2> lim.perc)]~round(T$phase.int.m2[j2&(perc.m2> lim.perc)]/2))




	perc.var = rep(1,nrow(T))
	perc.var[which(T$LRT.best.model == 2)] = perc.m2[which(T$LRT.best.model == 2)]
	perc.var[which(T$LRT.best.model == 3)] = perc.m3[which(T$LRT.best.model == 3)]
	perc.var[which(T$LRT.best.model == 4)] = perc.m4[which(T$LRT.best.model == 4)]
	
	
	h = hist(T$LRT.best.model, breaks = seq(0.5,4.5,by = 1), plot = TRUE, border = FALSE, col = 5:8)
	h = hist(T$LRT.best.model[perc.var>0.4], breaks = seq(0.5,4.5,by = 1), plot = TRUE, add = TRUE, col = 1:4, border = FALSE)

	100*sum(h$counts[2:4])/sum(h$counts[1:4])
	100*h$counts[2]/sum(h$counts[2:4])
	100*h$counts[3]/sum(h$counts[2:4])
	100*h$counts[4]/sum(h$counts[2:4])
	
	
	par(mfcol = c(2,2))
	
	breaks = seq(0,24,by = 1)
	hist(T$phase.gamma.m3[j3], breaks = breaks, col = 3+4, border = NA, freq = TRUE)
	hist(T$phase.gamma.m3[j3&(perc.var>0.4)], breaks = breaks, add = TRUE, col = 'transparent', border = 3, freq = TRUE)
	
	hist(T$phase.gamma.m3[j3], breaks = breaks, col = 3+4, border = NA, freq = FALSE, ylim = c(0,0.3))
	hist(T$phase.gamma.m3[j3&(perc.var>0.4)], breaks = breaks, add = TRUE, col = 'transparent', border = 3, freq = FALSE)
	

	hist(T$phase.gamma.m4[j4], breaks = breaks, col = 4+4, border = NA, freq = TRUE)
	hist(T$phase.gamma.m4[j4&(perc.var>0.4)], breaks = breaks, add = TRUE, col = 'transparent', border = 4, freq = TRUE)
	
	hist(T$phase.gamma.m4[j4], breaks = breaks, col = 4+4, border = NA, freq = FALSE, ylim = c(0,0.3))
	hist(T$phase.gamma.m4[j4&(perc.var>0.4)], breaks = breaks, add = TRUE, col = 'transparent', border = 4, freq = FALSE)
	
	
	
	#### cross-cor between the two days
	
	
	cross.cor.two.days = function(x){x = unlist(x); return(cor(x[1:12],x[13:24]))}
	cross.cor.int = apply(T[,ZT.int], 1, cross.cor.two.days)
	cross.cor.ex = apply(T[,ZT.ex], 1, cross.cor.two.days)
	
	par(mfcol = c(3,2))
	plot(cross.cor.int, cross.cor.ex, pch = 16, col = T$LRT.best.model+(-sign(T$LRT.best.model-1.5)+1)*2, cex = 0.7)
	
	boxplot(cross.cor.int~T$LRT.best.model, col = 5:8, border = 1:4, main = 'cross-correlation between the two days\nintronic signal', pch = 16)
	boxplot(cross.cor.ex~T$LRT.best.model, col = 5:8, border = 1:4, main = 'cross-correlation between the two days\nexonic signal', pch = 16)

	plot(cross.cor.int, T$LRT.quality, col = T$LRT.best.model+(-sign(T$LRT.best.model-1.5)+1)*2, cex = 0.7, pch = 16)
	plot(cross.cor.ex, T$LRT.quality, col = T$LRT.best.model+(-sign(T$LRT.best.model-1.5)+1)*2, cex = 0.7, pch = 16)
	
	sum(T$LRT.best.model == 3, na.rm = TRUE)
	sum((T$LRT.best.model == 3)&(cross.cor.int>0), na.rm = TRUE)
	sum((T$LRT.best.model == 3)&(cross.cor.ex>0), na.rm = TRUE)
	sum((T$LRT.best.model == 3)&(cross.cor.int>0)&(cross.cor.ex>0), na.rm = TRUE)

	sum(T$LRT.best.model == 2, na.rm = TRUE)
	sum((T$LRT.best.model == 2)&(cross.cor.int>0), na.rm = TRUE)
	sum((T$LRT.best.model == 2)&(cross.cor.ex>0), na.rm = TRUE)
	sum((T$LRT.best.model == 2)&(cross.cor.int>0)&(cross.cor.ex>0), na.rm = TRUE)
	
	
		

	var.int = apply(T[,ZT.int], 1, var)
	sd.int = sqrt(var.int)
	var.ex = apply(T[,ZT.ex], 1, var)
	sd.ex = sqrt(var.ex)
	
	par(mfcol = c(4,3))
	plot(sd.int, cross.cor.int,pch = 16, col = T$LRT.best.model+(-sign(T$LRT.best.model-1.5)+1)*2, cex = 0.7, log = 'x')
	plot(sd.int[j2], cross.cor.int[j2],pch = 16, col = 2, cex = 0.7, log = 'x', xlim = range(sd.int),ylim = range(cross.cor.int))
	abline(a = -0.5, b = -1)
	plot(sd.int[j3], cross.cor.int[j3],pch = 16, col = 3, cex = 0.7, log = 'x', xlim = range(sd.int),ylim = range(cross.cor.int))
	#plot(sd.int[j3], cross.cor.int[j3],pch = 16, col = sign(T$LRT.pval.d[j3]-T$LRT.pval.A[j3])+2, cex = 0.7, log = 'x', xlim = range(sd.int),ylim = range(cross.cor.int))

	abline(a = -0.5, b = -1)
	plot(sd.int[j4], cross.cor.int[j4],pch = 16, col = 4, cex = 0.7, log = 'x', xlim = range(sd.int),ylim = range(cross.cor.int))
	abline(a = -0.5, b = -1)
	
	
	h3 = hist(sd.int[T$LRT.best.model==3], breaks = seq(0,2,by = 0.02), col = 3, border = NA, xlim = range(0,sd.int, na.rm = TRUE))
	abline(v = 0.2, col = 'gray40'); 
	h2 = hist(sd.int[T$LRT.best.model==2], breaks =  seq(0,2,by = 0.02), col = 2, border = NA, xlim = range(0,sd.int, na.rm = TRUE))
	abline(v = 0.2, col = 'gray40');
	
	
	
	plot(sd.int, sqrt(T$intron.var),pch = 16, col = T$LRT.best.model+(-sign(T$LRT.best.model-1.5)+1)*2, cex = 0.7, log = 'x')
	plot(sd.int[T$LRT.best.model==2], T$fold.change.int.m2[T$LRT.best.model==2],pch = 16, col = T$LRT.best.model+(-sign(T$LRT.best.model-1.5)+1)*2, cex = 0.7, ylim = c(0,50))

	plot(sd.int, sd.ex,pch = 16, col = T$LRT.best.model+(-sign(T$LRT.best.model-1.5)+1)*2, cex = 0.7)#, log = 'yx')
	abline(a = 0, b = 1)
	
	plot(sd.int[j2], sd.ex[j2],pch = 16, col = 2, cex = 0.7, xlim = c(0,1.5), ylim = c(0,1.5))#, log = 'yx')
	abline(a = 0, b = 1)
	
	plot(sd.int[j3], sd.ex[j3],pch = 16, col = 3, cex = 0.7, xlim = c(0,1.5), ylim = c(0,1.5))#, log = 'yx')
	abline(a = 0, b = 1)
	
	
	plot(sd.int[j4], sd.ex[j4],pch = 16, col = 4, cex = 0.7, xlim = c(0,1.5), ylim = c(0,1.5))#, log = 'yx')
	abline(a = 0, b = 1)


	par(mfcol = c(4,2));
	xlim = c(1, nrow(T));
	ylim = c(0,1)
	#xlim = c(1,2000);
	#ylim = c(0,0.1)
	FDR = 0.2
	plot(sort(T$LRT.pval.p), type = 'l', xlim = xlim, ylim = ylim); abline(a = 0, b = 1/nrow(T), col = 'gray');  abline(a = 0, b = 1/nrow(T)* FDR, col = 'orange'); abline(h = 0.05, col = 'lightblue')
	plot(sort(T$LRT.pval.d), type = 'l', xlim = xlim, ylim = ylim); abline(a = 0, b = 1/nrow(T), col = 'gray');  abline(a = 0, b = 1/nrow(T)*FDR, col = 'orange'); abline(h = 0.05, col = 'lightblue')
	plot(sort(T$LRT.pval.P), type = 'l', xlim = xlim, ylim = ylim); abline(a = 0, b = 1/nrow(T), col = 'gray');  abline(a = 0, b = 1/nrow(T)*FDR, col = 'orange'); abline(h = 0.05, col = 'lightblue')
	plot(sort(T$LRT.pval.D), type = 'l', xlim = xlim, ylim = ylim); abline(a = 0, b = 1/nrow(T), col = 'gray');  abline(a = 0, b = 1/nrow(T)*FDR, col = 'orange'); abline(h = 0.05, col = 'lightblue')
	LRT.pval.m4 = apply(cbind(T$LRT.pval.P,T$LRT.pval.D),1,max)
	plot(sort(LRT.pval.m4), type = 'l', xlim = xlim, ylim = ylim); abline(a = 0, b = 1/nrow(T), col = 'gray'); abline(a = 0, b = 1/nrow(T)*FDR, col = 'orange'); abline(h = 0.05, col = 'lightblue')
	LRT.combined = pchisq(-2*log(T$LRT.pval.P * T$LRT.pval.D), df = 4,lower.tail = FALSE)
	plot(sort(LRT.combined), type = 'l', xlim = xlim, ylim = ylim); abline(a = 0, b = 1/nrow(T), col = 'gray'); abline(a = 0, b = 1/nrow(T)*FDR, col = 'orange'); abline(h = 0.05, col = 'lightblue')

	plot(sort(T$LRT.pval.A), type = 'l', xlim = xlim, ylim = ylim); abline(a = 0, b = 1/nrow(T), col = 'gray'); abline(a = 0, b = 1/nrow(T)*FDR, col = 'orange'); abline(h = 0.05, col = 'lightblue')

	

	
	plot(T$phase.gamma.m4[j4],T$phase.int.m4[j4], col = 4+4, pch = 16 )
	abline(a = 0, b = 1); 	abline(a = 12, b = 1, lty = 2); abline(a = -12, b = 1, lty = 2);  abline(a = 24, b = 1); abline(a = -24, b = 1); 
	points(T$phase.gamma.m4[j4&(perc.var>0.4)],T$phase.int.m4[j4&(perc.var>0.4)], col = 4, pch = 16 )
	
	
	cols = rainbow(n = max(round(T$fold.change.int.m4[j4&(perc.var>0.4)]), na.rm = TRUE),s = 0.8, v = 0.8, start = 0.5); cols[1] = 'gray'
	plot(T$phase.gamma.m4[j4&(perc.var>0.4)],T$phase.int.m4[j4&(perc.var>0.4)], col = cols[round(T$fold.change.int.m4[j4&(perc.var>0.4)])], pch = 1, type  = 'p', cex = 3.5, xlim = c(0,24), ylim = c(0,24))
	text(T$phase.gamma.m4[j4&(perc.var>0.4)],T$phase.int.m4[j4&(perc.var>0.4)],round(T$fold.change.int.m4[j4&(perc.var>0.4)], digits = 1), col = cols[round(T$fold.change.int.m4[j4&(perc.var>0.4)])],)
	abline(a = 0, b = 1); 	abline(a = 12, b = 1, lty = 2); abline(a = -12, b = 1, lty = 2);  abline(a = 24, b = 1); abline(a = -24, b = 1); 
	
	rel.ampl.ex = pmin(T$rel.ampl.ex[j4&(perc.var>0.4)],0.99)
	val = round(rel.ampl.ex, digits = 1) * 10 +1
	cols = rainbow(n = max(val, na.rm = TRUE),s = 0.8, v = 0.8, start = 0.5); cols[1] = 'gray'
	plot(T$phase.gamma.m4[j4&(perc.var>0.4)],T$phase.int.m4[j4&(perc.var>0.4)], col = cols[val], pch = 1, type  = 'p', cex = 3.5, xlim = c(0,24), ylim = c(0,24))
	text(T$phase.gamma.m4[j4&(perc.var>0.4)],T$phase.int.m4[j4&(perc.var>0.4)],round((1+ rel.ampl.ex)/(1-rel.ampl.ex), digits = 1), col = cols[val])
	abline(a = 0, b = 1); 	abline(a = 12, b = 1, lty = 2); abline(a = -12, b = 1, lty = 2);  abline(a = 24, b = 1); abline(a = -24, b = 1); 
	
	
	
	
	
	cols = rainbow(n = 24, s = 0.8, v = 0.9)
	cols.val = cols[round(T$phase.int.m4[j4&(perc.var>0.4)])]
	phase.diff = (T$phase.gamma.m4[j4&(perc.var>0.4)] - T$phase.int.m4[j4&(perc.var>0.4)])%%24
	fc.int = T$fold.change.int.m4[j4&(perc.var>0.4)]
	
	x = fc.int*cos(phase.diff/24*2*pi); y = fc.int*sin(phase.diff/24*2*pi); 
	x.c = cos(seq(0,2*pi+0.1, by = 0.1)); y.c = sin(seq(0,2*pi+0.1, by = 0.1))
	lim = range(x,y,-x,-y, na.rm = TRUE)
	plot(0,0, pch = 16, cex = 0.5, col = 'gray', xlim = lim, ylim = lim, main = 'radius  = intron FC')
	points(x.c,y.c,type = 'l', col = 'gray')
	points(2*x.c,2*y.c,type = 'l', col = 'gray')
	points(3*x.c,3*y.c,type = 'l', col = 'gray')
	points(x,y, pch = 16)#, col = cols.val)
	text(max(lim)/2,0,'in phase', pos = 4, col = )
	text(-max(lim)/2,0,'anti-phase', pos = 2)
	#x.cc = cos(seq(0,2*pi, length = 25)); y.cc = sin(seq(0,2*pi, length = 25))
	#points(max(lim)+ 0.4*x.cc, max(lim)+ 0.4*y.cc,pch = 16, col = cols )
	
	
	max.ex = apply(T[j4&(perc.var>0.4), ZT.ex],1,max)
	min.ex = apply(T[j4&(perc.var>0.4), ZT.ex],1,min)
	fc.ex = max.ex/min.ex
	
	x = fc.ex*cos(phase.diff/24*2*pi); y = fc.ex*sin(phase.diff/24*2*pi); 
	x.c = cos(seq(0,2*pi+0.1, by = 0.1)); y.c = sin(seq(0,2*pi+0.1, by = 0.1))
	lim = range(x,y,-x,-y, na.rm = TRUE)
	plot(0,0, pch = 16, cex = 0.5, col = 'gray', xlim = lim, ylim = lim, main = 'radius  = exon FC')
	points(x.c,y.c,type = 'l', col = 'gray')
	points(2*x.c,2*y.c,type = 'l', col = 'gray')
	points(3*x.c,3*y.c,type = 'l', col = 'gray')
	points(x,y, pch = 16)
	text(max(lim)/2,0,'in phase', pos = 4)
	text(-max(lim)/2,0,'anti-phase', pos = 2)
	
	x = fc.ex*cos(phase.diff/24*2*pi); y = fc.ex*sin(phase.diff/24*2*pi); 
	x.c = cos(seq(0,2*pi+0.1, by = 0.1)); y.c = sin(seq(0,2*pi+0.1, by = 0.1))
	lim = c(-10,10)
	plot(0,0, pch = 16, cex = 0.5, col = 'gray', xlim = lim, ylim = lim,  main = 'radius  = exon FC')
	points(x.c,y.c,type = 'l', col = 'gray')
	points(2*x.c,2*y.c,type = 'l', col = 'gray')
	points(3*x.c,3*y.c,type = 'l', col = 'gray')
	points(x,y, pch = 16)
	text(max(lim)/2,0,'in phase', pos = 4)
	text(-max(lim)/2,0,'anti-phase', pos = 2)
	
	
	fc.ratio = fc.ex/fc.int
	
	diff.phase.ex.int = (T$phase.ex[j4&(perc.var>0.4)]- T$phase.int[j4&(perc.var>0.4)])%%24
	cols = rainbow(24, s = 0.85, v = 0.85)
	color.function = colorRampPalette(c('lightblue2','lightblue3','lightblue4', 'black', 'lightblue4','lightblue2'))
	cols = color.function(24)
	x = fc.ratio*cos(phase.diff/24*2*pi); y = fc.ratio*sin(phase.diff/24*2*pi); 
	x.c = cos(seq(0,2*pi+0.1, by = 0.1)); y.c = sin(seq(0,2*pi+0.1, by = 0.1))
	lim = range(x,y,-x,-y, na.rm = TRUE)
	plot(0,0, pch = 16, cex = 0.5, col = 'gray', xlim = lim, ylim = lim, main = 'radius  = exon FC / intron FC')
	points(x.c,y.c,type = 'l', col = 'gray')
	points(2*x.c,2*y.c,type = 'l', col = 'gray')
	points(3*x.c,3*y.c,type = 'l', col = 'gray')
	points(x,y, pch = 16, col = cols[ceiling(diff.phase.ex.int)])
	text(max(lim)/2,0,'in phase', pos = 4)
	text(-max(lim)/2,0,'anti-phase', pos = 2)
	
	
	

	# comparison with Sharova h-l
	
	shar = read.table('../complementary_data/half-lives/Sharova.txt', stringsAsFactors = FALSE, sep = '\t', header = TRUE)
	fried = read.table('../complementary_data/half-lives/Friedel.txt', stringsAsFactors = FALSE, sep = '\t', header = TRUE)
	colnames(shar) = c('gene', 'hl')
	colnames(fried) = c('gene', 'probeset.ID','hlmin','hl')
	hist(shar$hl, breaks = 100)
	
	
	T$hl = NA
	T$hl[which(T$LRT.best.model == 2)] = log(2)/T$gamma.m2[which(T$LRT.best.model == 2)]
	T$hl[which(T$LRT.best.model == 3)] = log(2)/T$gamma.m3[which(T$LRT.best.model == 3)]
	T$hl[which(T$LRT.best.model == 4)] = log(2)/T$gamma.m4[which(T$LRT.best.model == 4)]


	pdf('plots/comparison_of_half_lives_with_other_datasets.pdf', width = 6.6, height = 3)
	par(mfrow = c(1,2),  cex = 0.6, mar = c(3,3,2,0.2)+0.1, mgp = c(1.6,0.5,0),las = 1, tcl = -0.3)

	m = match(T$gene , shar$gene )
	plot(shar$hl[m], T$hl, pch = 16, col = T$LRT.best.model, cex = 0.6, log = 'xy', main = 'Comparison of Half-lives - Sharova', xlab = 'Sharova et al.', ylab = 'MicroArray') #
	abline(a = 0, b = 1)
	
	n = match(T$gene, fried$gene)
	plot(fried$hl[n], T$hl, pch = 16, col = T$LRT.best.model, cex = 0.6, log = 'xy', main = 'Comparison of Half-lives - Friedel', xlab = 'Friedel et al.', ylab = 'MicroArray') #
	abline(a = 0, b = 1)

	dev.off()
	
	
	
	T$f2.ex = apply(T[, ZT.ex], 1,f24_R2)
	T$f2.int = apply(T[,ZT.int], 1,f24_R2)
	T$f2.qval.ex = fdrtool(T$f2.ex, statistic = 'pvalue', plot = FALSE)$qval
	T$f2.qval.int = fdrtool(T$f2.int, statistic = 'pvalue', plot = FALSE)$qval
	
	hist(T$LRT.best.model[T$f2.qval.ex<=0.1], breaks = seq(0.5,4.5,by=1), col = 5:8, border= FALSE)
	hist(T$LRT.best.model[T$f2.qval.int<=0.1], breaks = seq(0.5,4.5,by=1), col = 1:4, border= FALSE, add = TRUE)

	plot(-log10(T$f2.ex), -log10(T$f2.int), pch = 16, col = T$LRT.best.model+(-sign(T$LRT.best.model-1.5)+1)*2, cex = 0.5)
	
	
	out.ex = ((T$LRT.best.model == 1) & !is.na(T$LRT.best.model ) & (T$f2.ex<= 0.000001 ) )
	sum(out.ex)
	text(-log10(T$f2.ex[out.ex]), -log10(T$f2.int[out.ex]),T$gene[out.ex], pos = 3, col = 1, cex = 0.6)
	
	out.int = ((T$LRT.best.model == 1) & !is.na(T$LRT.best.model ) & (T$f2.int<= 0.0001 ) )
	sum(out.int)
	text(-log10(T$f2.ex[out.int]), -log10(T$f2.int[out.int]),T$gene[out.int], pos = 3, col = 1, cex = 0.6)
	
	
	
	
	
	out = out.int
	i = 4
	
	
	plot(zt, unlist(T[which(out)[i],ZT.ex]), type = 'b', col = 'steelblue', ylim = range(0,2,unlist(T[which(out)[i],ZT.ex]),unlist(T[which(out)[i],ZT.int])))
	points(zt, unlist(T[which(out)[i],ZT.int]), type = 'b', col = 'green3')
	
	zt.p = seq(0,46,by = 0.1)
	#fit.ex = compute.m.sigmoid(t = zt.p, gamma = T$gamma.m3[which(out)[i]],eps.gamma = T$eps.gamma.m3[which(out)[i]], phase.gamma = T$phase.gamma.m3[which(out)[i]] , fold.change = 1)
	#points(zt.p, fit.ex, type = 'l', col = 'steelblue')

	T[which(out)[i],93:99]
	


	max.ex = apply(T[, ZT.ex],1,max)
	min.ex = apply(T[, ZT.ex],1,min)
	rel.ampl = pmin(T$rel.ampl.ex,1)
	fc.ex = apply(cbind(max.ex/min.ex, (1+ rel.ampl)/(1-rel.ampl)),1,mean)
	
	
	max.int = apply(T[, ZT.int],1,max)
	min.int = apply(T[, ZT.int],1,min)
	rel.ampl = pmin(T$rel.ampl.int,1)
	fc.int = apply(cbind(max.int/min.int, (1+ rel.ampl)/(1-rel.ampl)),1,mean)
	
	plot(T$fold.change.int.m2,fc.int)
	abline(a = 0, b = 1)
	
	par(mfrow = c(1,3))
	boxplot(fc.ex~T$LRT.best.model, col =5:8, border = 1:4, log = 'y', pch = 16, main = 'Fold change Exon')
	boxplot(fc.int ~T$LRT.best.model, col =5:8, border = 1:4, log = 'y', pch = 16, main = 'Fold change Intron')
	boxplot(T$fold.change.int.m2~T$LRT.best.model, col =5:8, border = 1:4, log = 'y', pch = 16, main = 'Fold change Intron : fold.change.int.m2')


}



###################################################################################
###################################################################################
###################################################################################



if(plot){
	cat('plot : FITS RESULTS ON REAL DATA ----- MCMC\n')
	if(load.for.plot){load(file = 'Rdata/gene_fits_results.Rdata')}
	source('functions.R')
	
	best.model = T$model.MCMC
		
	pdf('plots/fits_results_rhythmic_genes_MCMC.pdf', width = 6.6, height = 6)
	
	par(mfrow = c(3,3), cex = 0.6, mar = c(3,3,2,0.2)+0.1, mgp = c(1.6,0.5,0),las = 1, tcl = -0.3)
	
	m1 = T$rhythmic.ex.24 | T$rhythmic.int.24 | T$rhythmic.ex.12 | T$rhythmic.int.12
	m2 = rep(TRUE, nrow(T))
	M = cbind(m1,m2)
	
	for(i in 1:2){
		m = M[,i]
		
		var = apply(T[m,c(ZT.int, ZT.ex)],1,var)
		error = T$error.MCMC[m]/length(c(ZT.int, ZT.ex))
		perc.of.variance.explained = 100*(var-error)/var; perc.of.variance.explained  = apply(cbind(perc.of.variance.explained,-5),1,max)
		o = order(perc.of.variance.explained)
		plot(perc.of.variance.explained[o], ylim =c(-5,100), pch = '_', col = best.model[o], axes = FALSE, main = '% of var. explained by the fits', ylab = '', xlab = '')
		abline(h = 0); axis(2); box();if(gene.set == gene.sets[1]){axis(1, at = c(1:nrow(T)), labels = T$gene[o],las = 3)}
	

		h = hist(best.model[m],breaks = seq(0.5,4.5,by = 1), plot = FALSE)
		pie(h$counts, clockwise = TRUE, col = 1:4, border = FALSE, labels = c('CS-CD','RS-CD','CS-RD','RS-RD'))

		plot.new()


		bounds = set.bounds(model = 4)

		# half-lives 
		breaks = seq(log(2)/bounds$upper[1],log(2)/bounds$lower[1], len = 100); dgamma = breaks[2]-breaks[1]
		gamma = T$gamma.MCMC[m]
		breaks = seq(log(2)/max(gamma,na.rm = TRUE), log(2)/min(gamma,na.rm = TRUE)+1, by = dgamma)
		h = hist(log(2)/gamma[best.model[m] ==2], breaks = breaks, col = 2+3, border = 2, xlab = '[h]', ylab = '', main = 'Half-lives')
		h = hist(log(2)/gamma[best.model[m] ==3], breaks = breaks, col = 3+3, border = 3, add = TRUE)
		h = hist(log(2)/gamma[best.model[m] ==4], breaks = breaks, col = 4+3, border = 4, add = TRUE)
	
		# amplitude of degradation 
		breaks = seq(bounds$upper[2],bounds$lower[2], len = 50)
		eps.gamma = T$eps.gamma.MCMC[m]
		h = hist(eps.gamma[best.model[m] ==3], breaks = breaks, col = 3+3, border = 3, xlab = '',ylab = '', main = 'Amplitude of degradation')
		h = hist(eps.gamma[best.model[m] ==4], breaks = breaks, col = 4+3, border = 4, add = TRUE)
	
		# phase of maximal degradation
		breaks = seq(bounds$upper[3],bounds$lower[3], len = 25)
		phase.gamma = T$phase.gamma.MCMC[m]
		h = hist(phase.gamma[best.model[m] ==3], breaks = breaks, col = 3+3, border = 3, xlab = '[h]', ylab = '', main = 'Phase of maximal degradation', axes = FALSE)
		h = hist(phase.gamma[best.model[m] ==4], breaks = breaks, col = 4+3, border = 4, add = TRUE)
		axis(1, at = seq(0,24,by = 4)); axis(2)
		
		
		
		# half-lives 
		breaks = seq(log(2)/bounds$upper[1],log(2)/bounds$lower[1], len = 100); dgamma = breaks[2]-breaks[1]
		gamma = T$gamma.MCMC[m]; gamma = pmin(gamma,bounds$upper[1]); gamma = pmax(gamma,bounds$lower[1])
		h = hist(log(2)/gamma[best.model[m] ==2], breaks = breaks, col = 2+3, border = 2, xlab = '[h]', ylab = '', main = 'Half-lives')
		h = hist(log(2)/gamma[best.model[m] ==3], breaks = breaks, col = 3+3, border = 3, add = TRUE)
		h = hist(log(2)/gamma[best.model[m] ==4], breaks = breaks, col = 4+3, border = 4, add = TRUE)
	
		
		plot.new()
		plot.new()


		
		
	}

	dev.off()
}







if(plot){
	cat('plot : FITS RESULTS ON REAL DATA : comparison of phases\n')
	if(load.for.plot){
		if(gene.set == gene.sets[1]){load(file = "Rdata/gene_fit_results_example_genes_and_LRT.Rdata")}
		if(gene.set == gene.sets[2]){load(file = "Rdata/gene_fit_results_for_rhythmic_genes_and_LRT.Rdata")}
	}
	source('functions.R')
	best.model = T$LRT.best.model
	
	bounds = set.bounds(model = 4, parametrization = 'sigmoid')
	breaks = seq(bounds$upper[3],bounds$lower[3], len = 25)
		
	#if(gene.set == gene.sets[1]){pdf('plots/fits_results_example_genes_phase_of_degradation.pdf', width = 6.6, height = 5)}
	#if(gene.set == gene.sets[2]){pdf('plots/fits_results_rhythmic_genes_phase_of_degradation.pdf', width = 6.6, height = 5)}
	
	pdf('plots/fits_results_phase_of_degradation_FDRb.pdf', width = 6.6, height = 6)
	
	par(mfrow = c(3,3), cex = 0.6, mar = c(3,3,2,0.2)+0.1, mgp = c(1.6,0.5,0),las = 1, tcl = -0.3)

	Deg.phase.m3 = T$phase.gamma.m3[best.model ==3]
	Deg.phase.m4 = T$phase.gamma.m4[best.model ==4]
	
	
	h = hist(Deg.phase.m3, breaks = breaks, col = 3+4, border = 3, xlab = '[h]', ylab = '', main = 'Phase of maximal degradation', axes = FALSE)
	h = hist(Deg.phase.m4, breaks = breaks, col = 4+4, border = 4, add = TRUE)
	legend('topright', legend = c('model CS-RD', 'model RS-RD'), fill = c(7,8), border = c(3,4), bty = 'n')
	axis(1, at = seq(0,24,by = 4)); axis(2)
	
	plot.new()
	
	h = hist(c(Deg.phase.m3,Deg.phase.m4), breaks = breaks, col = 3, border = 3, xlab = '[h]', ylab = '', main = 'Phase of maximal degradation', axes = FALSE)
	legend('topright', legend = c('model CS-RD & RS-RD'), fill = c(3), border = c(3), bty = 'n')
	axis(1, at = seq(0,24,by = 4)); axis(2)

	
	##### mRNA
	

	h = hist(T$phase.ex[(best.model ==2)], breaks = breaks, col = 'transparent', border = 'steelblue', xlab = '[h]', ylab = '', main = 'Phase of mRNA', axes = FALSE)
	h = hist(T$phase.ex[(best.model ==3)], breaks = breaks, col = rgb(0.6,0,0.9,0.1), border =  rgb(0.2,0,0.9,1), add = TRUE)
	h = hist(T$phase.ex[(best.model ==4)], breaks = breaks, col = rgb(0,0,0.5,0.3), border =  rgb(0,0,0.5,1), add = TRUE)
	legend('topleft', legend = c('model RS-CD', 'model CS-RD','model RS-RD'), fill = c('transparent',rgb(0.6,0,0.9,0.1),rgb(0,0,0.5,0.3)), border = c('steelblue',rgb(0.2,0,0.9,1),rgb(0,0,0.5,1)), bty = 'n')
	axis(1, at = seq(0,24,by = 4)); axis(2)

	delay.m3 = (T$phase.ex[(best.model ==3)] - Deg.phase.m3)%%24
	delay.m4 = (T$phase.ex[(best.model ==4)] - Deg.phase.m4)%%24
	h = hist(delay.m3, breaks = breaks, col = 3+4, border = 3, xlab = '[h]', ylab = '', main = 'Delay between mRNA \nand Degradation', axes = FALSE)
	h = hist(delay.m4, breaks = breaks, col = 4+4, border = 4, add = TRUE)
	axis(1, at = seq(0,24,by = 4)); axis(2)
	
	h = hist(c(delay.m3, delay.m4), breaks = breaks, col = 3, border = 3, xlab = '[h]', ylab = '', main = 'Delay between mRNA \nand Degradation', axes = FALSE)
	axis(1, at = seq(0,24,by = 4)); axis(2)
	

	
	##### pre-mRNA
	
	h = hist(T$phase.int[(best.model ==2)|(best.model ==4)], breaks = breaks, col = 'transparent', border = 'green3', xlab = '[h]', ylab = '', main = 'Phase of pre-mRNA', axes = FALSE)
	h = hist(T$phase.int[(best.model ==4)], breaks = breaks, col = 'green3', border = 'green3', add = TRUE)
	legend('topleft', legend = c('model RS-RD & RS-CD', 'model RS-RD'), fill = c('transparent','green3'), border = 'green3', bty = 'n')
	axis(1, at = seq(0,24,by = 4)); axis(2)
	
	delay.m4 = (T$phase.int[(best.model ==4)] - Deg.phase.m4)%%24
	h = hist(delay.m4, breaks = breaks, col = 4+4, border = 4, xlab = '[h]', ylab = '', main = 'Delay between pre-mRNA \nand Degradation', axes = FALSE)
	axis(1, at = seq(0,24,by = 4)); axis(2)

	plot.new()
	
	#plot(T$phase.ex[(best.model == 3)], T$phase.gamma.m3[(best.model == 3)], pch = 16, col = 3+4);
	#points(T$phase.ex[(best.model == 4)], T$phase.gamma.m4[(best.model == 4)], pch = 16, col = 4+4);
	#abline(a = 12, b = 1, col = 'gray'); abline(a = -12, b = 1, col = 'gray')
	#abline(a = 6, b = 1, col = 'gray'); abline(a = -18, b = 1, col = 'gray')
	
	#plot((T$phase.ex[(best.model == 4)]-T$phase.gamma.m4[(best.model == 4)])%%24,(T$fold.change.int.m4[(best.model == 4)]-T$phase.gamma.m4[(best.model == 4)])%%24, pch = 16, xlab = 'diff mRNA-degradation', ylab = 'diff pre-mRNA - degradation', xlim = c(0,24), ylim = c(0,24))

	
	
	dev.off()
}





if(plot){
	cat('plot : FITS RESULTS ON REAL DATA : comparison of phases ----- MCMC\n')
	if(load.for.plot){load(file = "Rdata/gene_fits_results.Rdata")}

	source('functions.R')
	best.model = T$model.MCMC
	
	bounds = set.bounds(model = 4)
	breaks = seq(bounds$upper[3],bounds$lower[3], len = 25)
		
	pdf('plots/fits_results_phase_of_degradation.pdf', width = 6.6, height = 5)
	
	par(mfrow = c(3,3), cex = 0.6, mar = c(3,3,2,0.2)+0.1, mgp = c(1.6,0.5,0),las = 1, tcl = -0.3)
	
	m1 = T$rhythmic.ex.24 | T$rhythmic.int.24 | T$rhythmic.ex.12 | T$rhythmic.int.12
	m2 = rep(TRUE, nrow(T))
	M = cbind(m1,m2)
	
	for(i in 1:2){
		m = M[,i]
		Deg.phase.m3 = T$phase.gamma.MCMC[m & (best.model ==3)]
		Deg.phase.m4 = T$phase.gamma.MCMC[m & (best.model ==4)]
	
	
		h = hist(Deg.phase.m3, breaks = breaks, col = 3+3, border = 3, xlab = '[h]', ylab = '', main = 'Phase of maximal degradation', axes = FALSE)
		h = hist(Deg.phase.m4, breaks = breaks, col = 4+3, border = 4, add = TRUE)
		legend('topright', legend = c('model CS-RD', 'model RS-RD'), fill = c(6,7), border = c(3,4), bty = 'n')
		axis(1, at = seq(0,24,by = 4)); axis(2)
		
		plot.new()
		
		h = hist(c(Deg.phase.m3,Deg.phase.m4), breaks = breaks, col = 3, border = 3, xlab = '[h]', ylab = '', main = 'Phase of maximal degradation', axes = FALSE)
		legend('topright', legend = c('model CS-RD & RS-RD'), fill = c(3), border = c(3), bty = 'n')
		axis(1, at = seq(0,24,by = 4)); axis(2)

	
		##### mRNA
				
	
		h = hist(T$phase.ex[m&((best.model ==2)|(best.model ==3)|(best.model ==4))], breaks = breaks, col = 'transparent', border = 'steelblue', xlab = '[h]', ylab = '', main = 'Phase of mRNA', axes = FALSE)
		h = hist(T$phase.ex[m&(best.model ==3)], breaks = breaks, col = 'steelblue1', border = 'steelblue', add = TRUE)
		h = hist(T$phase.ex[m&(best.model ==4)], breaks = breaks, col = rgb(0,0.2,0.6,0.2), border = 'steelblue4', add = TRUE, lwd = 3)
		legend('topleft', legend = c('all rhythmic mRNA', 'model CS-RD', 'model RS-RD'), fill = c('transparent','steelblue1',rgb(0,0.2,0.6,0.2)), border = c('steelblue','steelblue','steelblue4'), bty = 'n')
		axis(1, at = seq(0,24,by = 4)); axis(2)

		delay.m3 = (T$phase.ex[m&(best.model ==3)] - Deg.phase.m3)%%24
		delay.m4 = (T$phase.ex[m&(best.model ==4)] - Deg.phase.m4)%%24
		h = hist(delay.m3, breaks = breaks, col = 3+3, border = 3, xlab = '[h]', ylab = '', main = 'Delay between mRNA \nand Degradation', axes = FALSE)
		h = hist(delay.m4, breaks = breaks, col = 4+3, border = 4, add = TRUE)
		axis(1, at = seq(0,24,by = 4)); axis(2)
	
		h = hist(c(delay.m3, delay.m4), breaks = breaks, col = 3, border = 3, xlab = '[h]', ylab = '', main = 'Delay between mRNA \nand Degradation', axes = FALSE)
		axis(1, at = seq(0,24,by = 4)); axis(2)
	

	
		##### pre-mRNA
	
		h = hist(T$phase.int[m&((best.model ==2)|(best.model ==4))], breaks = breaks, col = 'transparent', border = 'green3', xlab = '[h]', ylab = '', main = 'Phase of pre-mRNA', axes = FALSE)
		h = hist(T$phase.int[m&(best.model ==4)], breaks = breaks, col = 'green3', border = 'green3', add = TRUE)
		legend('topleft', legend = c('model RS-RD & RS-CD', 'model RS-RD'), fill = c('transparent','green3'), border = 'green3', bty = 'n')
		axis(1, at = seq(0,24,by = 4)); axis(2)
		
		delay.m4 = (T$phase.int[m&(best.model == 4)] - Deg.phase.m4)%%24
		h = hist(delay.m4, breaks = breaks, col = 4+3, border = 4, xlab = '[h]', ylab = '', main = 'Delay between pre-mRNA \nand Degradation', axes = FALSE)
		axis(1, at = seq(0,24,by = 4)); axis(2)

		plot.new()
	}
	
	dev.off()
}




################################################################################################################################################################
################################################################################
################################################################################ most model selection parts stop here
################################################################################
################################################################################################################################################################




############################################
############### PLOT INDIVIDUAL GENE SUMMARY WITH FITS


transparent.factor = 0.3
transparent.gray = rgb(col2rgb('gray')[1]/255,col2rgb('gray')[2]/255,col2rgb('gray')[3]/255, transparent.factor)
transparent.green3 = rgb(col2rgb('green3')[1]/255,col2rgb('green3')[2]/255,col2rgb('green3')[3]/255, transparent.factor)
transparent.tomato = rgb(col2rgb('tomato')[1]/255,col2rgb('tomato')[2]/255,col2rgb('tomato')[3]/255, transparent.factor)
transparent.black = rgb(col2rgb('black')[1]/255,col2rgb('black')[2]/255,col2rgb('black')[3]/255, transparent.factor)

palette(c('gray','green3','tomato','black', transparent.gray ,transparent.green3, transparent.tomato, transparent.black))


if(plot){
	if(load.for.plot){load('Rdata/gene_fits_results.Rdata')}
	RefSeq = read.table('../complementary_data/RefSeq_Mouse_mm9.txt', stringsAsFactors = FALSE, sep = '\t', header = TRUE)
	source('functions.R') ;
	

	#delay.deg.int = (T$phase.int.MCMC[T$model.MCMC == 4] - T$phase.gamma.MCMC[T$model.MCMC == 4])%%24
	#o = order(delay.deg.int); cbind(T$pval.int[T$model.MCMC == 4][o], T$gene[T$model.MCMC == 4][o])

	gene_list = c("Per1","Per2","Per3","Cry1","Cry2","Dbp","Tef","Hlf","Nr1d1","Nr1d2","Rora","Rorb","Rorc","Arntl","Bmal1","Clock","Npas2","Bhlhe40","Bhlhe41","Cirbp","Tfrc","Hamp","Hamp2","Nr4a2","Wee1", "Por", "Alb")
	gene_list_small = c("Tspan32","Atp5sl","Sepp1","Nsdhl","Dbi","Fuk","Cpn1_chr19","Pten_chr19","Igsf11","Jmjd8","Arl10")
	gene_list_large = c("Ahsg","Acss2","Pzp","Qdpr","B2m","Pcca","Gpr156","Mmgt2","Itih3","Atf4","Serhl")
	gene_list_other = c('Ppard',"Cdkn1a")

	gene_list_menet = c('Phb2','Wdr91','Klf15','Etnk2','Zrsr1','Eme2')
	gene_list = unique(c(gene_list, gene_list_small, gene_list_large, gene_list_other,gene_list_menet))
	
	gene_list = 'Efhd1'	
	
	gene_list = T$gene[which(T$LRT.best.model==2)]
		
	#pdf.folder = 'plots/individual_genes_full_plots/individual_genes_full_plots_with_fits_MCMC'
	#for(gene in gene_list){cat(gene, '\n');plot.this.gene.full(pdf.folder = pdf.folder,gene = gene, ylim.rel.type = 'adjusted', ylim.rel = c(0,2), log2.rel = FALSE, with.fit = TRUE, fitting.method = 'MCMC', parametrization = 'sigmoid')}
	pdf.folder = 'plots/individual_genes_full_plots/individual_genes_full_plots_with_fits_Optim_all_m2-3-4_FDRb/model2'
	for(gene in gene_list){cat(gene, '\n');plot.this.gene.full(pdf.folder = pdf.folder,gene = gene, ylim.rel.type = 'adjusted', ylim.rel = c(0,2), log2.rel = FALSE, with.fit = TRUE, fitting.method = 'Optim', parametrization = 'sigmoid')}
	
	
	gene_list = T$gene[which(T$LRT.best.model==3)]
	pdf.folder = 'plots/individual_genes_full_plots/individual_genes_full_plots_with_fits_Optim_all_m2-3-4_FDRb/model3'
	for(gene in gene_list){cat(gene, '\n');plot.this.gene.full(pdf.folder = pdf.folder,gene = gene, ylim.rel.type = 'adjusted', ylim.rel = c(0,2), log2.rel = FALSE, with.fit = TRUE, fitting.method = 'Optim', parametrization = 'sigmoid')}
	
	
	gene_list = T$gene[which(T$LRT.best.model==4)]
	pdf.folder = 'plots/individual_genes_full_plots/individual_genes_full_plots_with_fits_Optim_all_m2-3-4_FDRb/model4'
	for(gene in gene_list){cat(gene, '\n');plot.this.gene.full(pdf.folder = pdf.folder,gene = gene, ylim.rel.type = 'adjusted', ylim.rel = c(0,2), log2.rel = FALSE, with.fit = TRUE, fitting.method = 'Optim', parametrization = 'sigmoid')}
	
}


##### FOR ALL GENES!!!

if(plot){
	if(load.for.plot){load('Rdata/gene_fits_results.Rdata')}
	RefSeq = read.table('../complementary_data/RefSeq_Mouse_mm9.txt', stringsAsFactors = FALSE, sep = '\t', header = TRUE)
	source('functions.R') ;
	

	gene_list = T$gene
	
	pdf.folder = '/Users/symul/Documents/Work/_x_MicroArray_all_genes_plots/individual_genes_full_plots_with_fits_Optim_ALL'
	
	already_done = list.files(pdf.folder)
	already_done = unlist(strsplit(already_done, split = '\\.pdf'))
		
	gene_list = setdiff(gene_list, already_done)	
		
	for(gene in gene_list){cat(gene, '\n');plot.this.gene.full(pdf.folder = pdf.folder,gene = gene, ylim.rel.type = 'adjusted', ylim.rel = c(0,2), log2.rel = FALSE, with.fit = TRUE, fitting.method = 'Optim', parametrization = 'sigmoid')}
		
}

####### copy the genes I need in a folder

if(do.step){
	if(load){load('Rdata/gene_fit_results_and_LRT_name_cleaned.Rdata')}
	source('functions.R') ;
	
	
	origin_directory = '/Users/symul/Documents/Work/_x_MicroArray_all_genes_plots/individual_genes_full_plots_with_fits_Optim_ALL'
	target_directory_name = 'plots/individual_genes_full_plots/individual_genes_full_plots_with_fits_Optim_and_AIC/model3'
	dir.create(target_directory_name)
	
	gene_list = T$gene[(AIC.c$best.model == 3)]; gene_list = gene_list[!is.na(gene_list)]
	
	file_names = paste(origin_directory,'/',gene_list, '.pdf', sep = '')
		
	cmd = paste('cp ', paste(file_names, collapse = ' '),' ', target_directory_name,'/.', sep = '')
	system(cmd)
	
}




############################################
############### Comparison with Cyclix Lists


if(plot){
	if(load.for.plot){load('Rdata/gene_fits_results.Rdata')}
	Cyclix.S4 = read.table('../complementary_data/Cyclix/CycliX2012_table_S4.txt', sep = '\t', header = TRUE)
	m = match(T$gene, Cyclix.S4$name)
	Cyclix = Cyclix.S4[m,]
	Cyclix$class.Fig6[is.na(Cyclix$class.Fig6)&!is.na(Cyclix$name)] = 0
	
	
	
	cats = paste(Cyclix$class.Fig6,T$model.MCMC)
	rle = rle(sort(cats))
	rle.values = rle$values
	rle.counts = rle$length
	
	rle.counts.fake = c(max(rle.counts),0,0,0,0,0,max(rle.counts)/5,0,max(rle.counts)/60,0,0,max(rle.counts)/5,0,max(rle.counts)/30,0,0,0,max(rle.counts)/30,0,0,0,0,0,0,0)
	RLE = cbind(rle.counts.fake,rle.counts); mains = c('Theory','Results')
	
	
	pdf('plots/Comparison_Cyclix_list_with_fits_results.pdf', width = 6.6, height = 3)
	
	par(mfrow = c(1,2), cex = 0.6, mar = c(3,3,2,0.2)+0.1, mgp = c(1.6,0.5,0),las = 1, tcl = -0.3)
	
	for(i in 1:2){
		rle.counts = RLE[,i]; 
		plot(1,1,type = 'n', axes = FALSE, xlim = c(-0.5,4.2), ylim = c(0.5,4.2), xlab = 'Cyclix Class', ylab = 'Exon-Array MCMC' , main = mains[i])
		axis(1, at = 0:3, labels = c('---','1','2','3'), lwd = 0)
		axis(2, at = 1:4, labels = c('CS-CD','RS-CD','CS-RD','RS-RD'), lwd = 0, las = 3)
		legend('bottomright',legend = c('MCMC prediction','Cyclix Class') , pch = c(1,16))
		for(cc in 0:3){
			for(tt in 1:4){
				ind = grep(paste(cc,tt),rle$values)
				points(cc,tt,cex = sqrt(rle.counts[ind])/5, pch = 21, bg = c('gray80','gray40','orange','red')[cc+1], col = c('gray60','green4','tomato','black')[tt], lwd = 2)
			}
		}
		
	}
	dev.off()

}




############################################
############### Comparison with Carla Green paper : Poly(A) Rhythms.




if(plot){
	if(load.for.plot){load('Rdata/gene_fits_results.Rdata')}
	
	RESS.o = read.table('../complementary_data/Green-Poly_A/SupplementalTable1_RESSmRNAs.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
	m = match(T$gene, RESS.o$Gene)
	RESS = RESS.o[m,]
	
	PAR.o = read.table('../complementary_data/Green-Poly_A/SupplementalTable2_PARmRNAs.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
	m = match(T$gene, PAR.o$Gene)
	PAR = PAR.o[m,]
	
	
	RefSeq = read.table('../complementary_data/RefSeq_Mouse_mm9.txt', stringsAsFactors = FALSE, sep = '\t', header = TRUE)
	source('functions.R') ; intense.debug = FALSE
	pdf.folder = 'plots/comparison_with_C_Green_paper/'
	
	gene.list = c('Cyp2a4', 'Hsd17b6','Polr2h','Hspb1','Gstt2','Mif','Slc44a3','Pdzk1ip1')
		
	
	for(gene in gene.list){cat(gene, '\n');plot.this.gene.full(pdf.folder = pdf.folder,gene = gene, ylim.rel.type = 'adjusted', ylim.rel = c(0,2), log2.rel = FALSE, with.fit = TRUE, fitting.method = 'MCMC')}



	transparent.factor = 0.3
	transparent.green3 = rgb(col2rgb('green3')[1]/255,col2rgb('green3')[2]/255,col2rgb('green3')[3]/255, transparent.factor)
	transparent.tomato = rgb(col2rgb('tomato')[1]/255,col2rgb('tomato')[2]/255,col2rgb('tomato')[3]/255, transparent.factor)
	transparent.black = rgb(col2rgb('black')[1]/255,col2rgb('black')[2]/255,col2rgb('black')[3]/255, transparent.factor)

	palette(c('gray','green3','tomato','black', transparent.green3, transparent.tomato, transparent.black))




	j3 = (T$model.MCMC==3) & (!is.na(PAR$Gene)) & (!is.na(T$model.MCMC))
	j4 = (T$model.MCMC==4) & (!is.na(PAR$Gene)) & (!is.na(T$model.MCMC))
	j = j3|j4
	T.phase.gamma = T$phase.gamma.MCMC

	
	j3 = (T$LRT.best.model ==3) & (!is.na(PAR$Gene)) & (!is.na(T$LRT.best.model))
	j4 = (T$LRT.best.model ==4) & (!is.na(PAR$Gene)) & (!is.na(T$LRT.best.model))
	j = j3|j4
	T.phase.gamma = rep(NA, nrow(T))
	T.phase.gamma[j3] = T$phase.gamma.m3[j3]
	T.phase.gamma[j4] = T$phase.gamma.m4[j4]

	
	pdf('plots/comparison_with_C_Green_paper/phases_comparison.pdf', height = 9, width = 3.3)
	par(mfrow = c(3,1), cex = 0.7, mar = c(3,3,3,0.2)+0.1, mgp = c(1.6,0.5,0),las = 1, tcl = -0.3)
	
	breaks = seq(0,24,by = 2)
	hist(PAR$Peaktime[j3], breaks = breaks, col = 3+4, border = 3, xlab = 'time of max length of Poly(A) tail', main = 'Phase of Poly(A) tail', axes = FALSE)
	hist(PAR$Peaktime[j4], breaks = breaks, col = 4+4, border = 4, add = TRUE)
	legend('topleft', c('model 3 (CS-RD)', 'model 4 (RS-RD)'), fill= c(3,4)+4,  bty = 'n', border = c(3,4))
	axis(1, at = seq(0,24,by =6)); axis(2)
	
	plot(PAR$Peaktime[j3], T.phase.gamma[j3], pch = 16, col = 3, xlim = c(0,24), ylim =c(0,24), xlab = 'max length of Poly(A) tail', ylab = 'phase of max degradation')
	points(PAR$Peaktime[j4], T.phase.gamma[j4], pch = 16, col = 4)
	abline(a = 0, b = 1); abline(a = 12, b = 1, lty = 2); abline(a = -12, b = 1, lty = 2)
	legend('topright', c('model 3 (CS-RD)', 'model 4 (RS-RD)'), col= c(3,4), pch = 16, bty = 'n')
	
	phase.diff = (PAR$Peaktime-T.phase.gamma)%%24
	
	breaks = seq(0,24,by = 2)
	hist(phase.diff[j3], breaks = breaks, col = 3+4, border = 3, xlab = 'phase diff. [time of longuest Poly(A) tail - time of max degr.]', axes = FALSE)
	hist(phase.diff[j4], breaks = breaks, col = 4+4, border = 4, add = TRUE)
	legend('topright', c('model 3 (CS-RD)', 'model 4 (RS-RD)'), fill= c(3,4)+4,  bty = 'n', border = c(3,4))
	axis(1, at = seq(0,24,by =6)); axis(2)
	
	dev.off()

	
	Cyclix.S4 = read.table('../complementary_data/Cyclix/CycliX2012_table_S4.txt', sep = '\t', header = TRUE)
	m = match(T$gene, Cyclix.S4$name)
	Cyclix = Cyclix.S4[m,]
	Cyclix$class.Fig6[is.na(Cyclix$class.Fig6)&!is.na(Cyclix$name)] = 0
	
	
	
	cats = paste(Cyclix$class.Fig6,T$model.MCMC)
	rle = rle(sort(cats))
	rle.values = rle$values
	rle.counts = rle$length
	
	rle.counts.fake = c(max(rle.counts),0,0,0,0,0,max(rle.counts)/5,0,max(rle.counts)/60,0,0,max(rle.counts)/5,0,max(rle.counts)/30,0,0,0,max(rle.counts)/30,0,0,0,0,0,0,0)
	RLE = cbind(rle.counts.fake,rle.counts); mains = c('Theory','Results')
	
	
	pdf('plots/Comparison_Cyclix_list_with_fits_results.pdf', width = 6.6, height = 3)
	
	par(mfrow = c(1,2), cex = 0.6, mar = c(3,3,2,0.2)+0.1, mgp = c(1.6,0.5,0),las = 1, tcl = -0.3)
	
	for(i in 1:2){
		rle.counts = RLE[,i]; 
		plot(1,1,type = 'n', axes = FALSE, xlim = c(-0.5,4.2), ylim = c(0.5,4.2), xlab = 'Cyclix Class', ylab = 'Exon-Array MCMC' , main = mains[i])
		axis(1, at = 0:3, labels = c('---','1','2','3'), lwd = 0)
		axis(2, at = 1:4, labels = c('CS-CD','RS-CD','CS-RD','RS-RD'), lwd = 0, las = 3)
		legend('bottomright',legend = c('MCMC prediction','Cyclix Class') , pch = c(1,16))
		for(cc in 0:3){
			for(tt in 1:4){
				ind = grep(paste(cc,tt),rle$values)
				points(cc,tt,cex = sqrt(rle.counts[ind])/5, pch = 21, bg = c('gray80','gray40','orange','red')[cc+1], col = c('gray60','green4','tomato','black')[tt], lwd = 2)
			}
		}
		
	}
	dev.off()

}












plot = TRUE; load.for.plot = TRUE
do.step = TRUE; load = TRUE
stop('Normal. End of code. Everything went fine :-)')












if(do.step){
	cat("GETTING 3'UTR REGIONS FOR SELECTED GENES \n")
	if(load){load('Rdata/gene_fit_results_and_LRT.Rdata')}
	source('functions.R')
	
	gene_list = c("Per1","Per2","Per3","Cry1","Cry2","Dbp","Tef","Hlf","Nr1d1","Nr1d2","Arntl","Npas2","Tfrc","Wee1", "Por", "Alb")
	gene_list = T$gene
	
	library(biomaRt)
	ensembl = useMart("ensembl", dataset = 'mmusculus_gene_ensembl')
	#listAttributes(ensembl); listFilters(ensembl)

	
	annot = getBM(attributes = c('3utr','wikigene_name','ensembl_transcript_id'), filters = 'wikigene_name', values = gene_list, mart = ensembl) # other possible attributes : 'external_transcript_id','mgi_symbol','wikigene_name',
	colnames(annot) = c('sequence', 'gene','ensembl_ID')
	
	k = which(annot$sequence == 'Sequence unavailable')
	annot = annot[-k,]

	SEQ = gene_list
		
	for(i in 1:length(gene_list)){
		if(i%%400 == 0){cat(round(i/length(gene_list), digits = 2),'\n')}
		gene = gene_list[i]
		k = which(annot$gene == gene)	
		if(length(k)==0){cat(gene, '\t no sequence\n'); seq = ''}
		if(length(k)==1){seq = annot$sequence[k]; cat(gene,'\t only 1 sequence\n')}
		if(length(k)>1){
			seqs = annot$sequence[k]
			lengths = nchar(seqs);  cat(gene,'\t sequences length : ',sort(lengths),'\n')
			o = order(lengths); seqs = seqs[o]
			seq = ''
			for(l in 1:(length(seqs)-1)){
				seql = seqs[l];
				lim = 2560; if(nchar(seql)<lim){index = grep(seql, seqs)}else{index = l; cat('too long sequence! \n\n')}
				if(all(index<=l)){if(seq == ''){seq = seql}else{seq = paste(seq,seql, sep = '/')}}
			}
			seql = seqs[length(seqs)];
			if(seq == ''){seq = seql}else{seq = paste(seq,seql, sep = '/')}
		}
		SEQ[i] = seq
		#cat('\t',seq,'\n')	
	}
	
	SEQ = data.frame(gene = gene_list, sequence = SEQ, stringsAsFactors = FALSE)
	
	save(SEQ, file = 'Rdata/sequences_of_3UTR.Rdata')
	
}



if(do.step){
	cat("SCANNING 3'UTR REGIONS WITH RNA-BINDING-PROTEIN PWM \n")
	load('Rdata/sequences_of_3UTR.Rdata')
	load('Rdata/gene_fit_results_and_LRT_name_cleaned.Rdata')
	source('functions.R')
	library('gtools')
	N = nrow(SEQ)
	#matrix.list = read.table('../complementary_data/PWM_RNA_Binding_Proteins/matrix_list.txt', sep = '\t', header = FALSE, stringsAsFactors = FALSE) 
	matrix.list = read.table('../complementary_data/PWM_RNA_Binding_Proteins/matrix_list_PFM.txt', sep = '\t', header = FALSE, stringsAsFactors = FALSE) 
	M = nrow(matrix.list)
	
	#create a fasta file with all sequences [and an ID]
	
	FA = SEQ;
	FA$gene = paste('>',FA$gene,sep ='')
	fasta.file = '../complementary_data/scanning/fasta/fasta_all_genes_3UTR.fa'
	write.table(FA, file = fasta.file, quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE)
	
	
	# create the background weight matrix from all 3 UTR sequences
	all_bases = unlist(strsplit(SEQ[,2], split = ''))
	bases = c('A','C','G','T')
	stat = rep(0,length(bases))
	for(i in 1:length(bases)){
		b = bases[i]
		stat[i] = sum(all_bases == b)
		}
	tot = sum(stat);
	freq.bg = stat/tot
	background.weight.matrix.file = '../complementary_data/PWM_RNA_Binding_Proteins/3UTR_bg.mat'
	write.table(matrix(freq.bg,nrow = 1, ncol = length(freq.bg)), file = background.weight.matrix.file , sep = '\t', col.names = FALSE, quote = FALSE)
	
	#create fasta files with all sequences possible with X bp
#	X = c(4:21)
#	for(x in X){
#		all.possible.seq = apply(matrix(bases[permutations(4,x, repeats.allowed = TRUE)],4^x,x),1,paste, collapse = '')
#		fasta.all.possible.seq = data.frame(ID = paste('>',1:4^x,sep = ''), seq = all.possible.seq, stringsAsFactors = FALSE)
#		fasta.all.possible.seq.file = paste('../complementary_data/scanning/fasta/sequences_for_',x,'.fa',sep = '')
#		write.table(fasta.all.possible.seq, file = fasta.all.possible.seq.file, quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE)
#		}
	
	
	
	
	
	#scan the sequences with all PWM
	cutoff = -1000;
	
	score.PWM.RBP.per.gene =  matrix(NA, nrow= nrow(T), ncol = nrow(matrix.list)); colnames(score.PWM.RBP.per.gene) = matrix.list[,1]\
	best.scores = rep(0, nrow(matrix.list))
	
	for(j in 1:M){
		cat(j,'\n')
		weight.matrix.file.pwm = paste('../complementary_data/PWM_RNA_Binding_Proteins/', matrix.list[j,1],'.pwm',sep ='') 
		matrix = t(read.table(weight.matrix.file.pwm, stringsAsFactors = FALSE))
		#weight.matrix.file.freq = paste('../complementary_data/PWM_RNA_Binding_Proteins/',matrix.list[j,1],'.pfm',sep ='')
		#freq.m = t(read.table(weight.matrix.file.freq, stringsAsFactors = FALSE))
		#sum.mot = apply(freq.m, 1, sum)
		#l = nrow(freq.m)
		
		dim(matrix)
		l = nrow(matrix)
		matrix = exp(matrix)/apply(exp(matrix),1,sum)
		
		weight.matrix.file = paste('../complementary_data/PWM_RNA_Binding_Proteins/',matrix.list[j,1],'.mat',sep ='') 
		#write.table(cbind(rep(1,nrow(freq.m)),freq.m/sum.mot), file = weight.matrix.file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
		write.table(cbind(rep(1,nrow(matrix)), matrix), file = weight.matrix.file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
		
		
		
		
		#make a pdf with the logo
		cmd =  paste('perl ../complementary_scripts/mat2pdf.pl',weight.matrix.file)
		system(cmd)
		
		# find the threshold score	 = 0.3* max score possible with PWM
		best.seq.fa.file = paste('../complementary_data/PWM_RNA_Binding_Proteins/',matrix.list[j,1],'_best_seq.fa',sep = '')
		best.score = 0;
		for(ll in 1:l){p = max(matrix[ll,]); b = freq.bg[which.max(matrix[ll,])]; best.score = best.score + log2(p/b)}
		best.scores[j] = best.score
		
#		i.max = apply(matrix, 1, which.max); best.seq = bases[i.max]
#		write.table(paste('>best.seq\n',paste(best.seq, collapse = ''),sep = ''),file = best.seq.fa.file, quote = FALSE, row.names = FALSE, col.names = FALSE)
#		output.file = paste('../complementary_data/scanning/results/',matrix.list[j,1],'_best_seq.txt',sep = '')
#		cmd = paste('./../complementary_scripts/S1K',weight.matrix.file,background.weight.matrix.file,cutoff, best.seq.fa.file,'>', output.file)
#		system(cmd)
#		res.best = read.table(output.file, sep = '\t',col.names = c('gene','seq','score','pos','strand'), stringsAsFactors = FALSE)
#		best.score = max(res.best$score)

		#fasta.all.possible.seq.file = paste('../complementary_data/scanning/fasta/sequences_for_',l,'.fa',sep = '')
		
		
		#output.file = paste('../complementary_data/scanning/results/',matrix.list[j,1],'_all_seq.txt',sep = '')
		#cmd = paste('./../complementary_scripts/S1K',weight.matrix.file,background.weight.matrix.file,cutoff, fasta.all.possible.seq.file,'>', output.file)
		#system(cmd)
		#res.all = read.table(output.file, sep = '\t',col.names = c('gene','seq','score','pos','strand'), stringsAsFactors = FALSE)
		
		
		
		# scan the 3 UTR sequences
		output.file = paste('../complementary_data/scanning/results/',matrix.list[j,1],'.txt',sep = '')
		cmd = paste('./../complementary_scripts/S1K',weight.matrix.file,background.weight.matrix.file,cutoff,fasta.file,'>', output.file)
		cat(cmd,'\n')
		system(cmd)
		
		# remove the lines with "/"
		light.output.file = paste('../complementary_data/scanning/results/',matrix.list[j,1],'_light.txt',sep = '')
		cmd = paste("grep -v '/'", output.file,">", light.output.file)
		system(cmd)
		
		

		# for each gene, keep the best value of the score
		best.output.file = paste('../complementary_data/scanning/results/',matrix.list[j,1],'_best.txt',sep = '')
		cmd = paste('perl ./../complementary_scripts/process_scanning_results.pl ', light.output.file, best.output.file)
		cat(cmd,'\n')
		system(cmd)
		

		#read the results
		res = read.table(best.output.file, sep = '\t',col.names = c('gene','seq','score','pos','strand'), stringsAsFactors = FALSE)
		save(res,file = paste('../complementary_data/scanning/results/',matrix.list[j,1],'.Rdata',sep = ''))
		m = match(T$gene, res$gene)
		
		## store the results together

		score.PWM.RBP.per.gene[,j] = res$score[m]
		
		#clean the directory and rm heavy files
		cmd = paste('rm', output.file, light.output.file, best.output.file )
		system(cmd);
		
		
						
	}

	colnames(score.PWM.RBP.per.gene) = matrix.list[,3]; names(best.scores) =  matrix.list[,3]
	save(best.scores ,score.PWM.RBP.per.gene, file = 'Rdata/RBP_scanning_score.Rdata')

	
}






if(do.step){
	cat("CHECKING miRNA MAPPING \n")
	if(load){load('Rdata/gene_fit_results_and_LRT.Rdata')}
	source('functions.R')
	
	### here I used the TargetScan data base downloaded from UCSC Table Browser. http://www.targetscan.org/mmu_50/
	# There is also another tool [that I have not used so far] from the Lab of Eran Segal - Weizmann Institute : http://genie.weizmann.ac.il/pubs/mir07/mir07_prediction.html
	
	miRNA = read.table('../complementary_data/miRNA/miRNA_mouse_genomic_target.txt', header = TRUE, sep = '\t', stringsAsFactors =FALSE)
	
	
	miRNA$gene = NA; miRNA$miRNA.name = NA; miRNA$number = NA
	for(i in 1:nrow(miRNA)){
		if(i%%400 == 0){cat(round(i/nrow(miRNA), digits = 2),'\n')}
		split = unlist(strsplit(miRNA$name[i],':'))
		miRNA$gene[i] = split[1]; miRNA$miRNA.name[i] = split[2]
		if(length(split)==3){miRNA$number[i] = split[3]}
		if(length(split)>3){cat(i,'\t',split,'\n')}
	}
	
	save(miRNA, file = 'Rdata/miRNA.Rdata')
	
	
	gene_list = c("Per1","Per2","Per3","Cry1","Cry2","Dbp","Tef","Hlf","Nr1d1","Nr1d2","Arntl","Npas2","Tfrc","Wee1", "Por", "Alb")
	gene_list = c("Per1","Per2","Per3")

	m = match(miRNA$gene, gene_list)
	miRNA[!is.na(m),]
	
		
}


###### building the matrix saying which 3 UTR is bound by which RBP or miRNA


if(do.step){
	cat("building the mapping matrix\n")
	if(load){load('Rdata/gene_fit_results_and_LRT_name_cleaned.Rdata'); load('Rdata/RBP_scanning_score.Rdata'); load('Rdata/miRNA.Rdata')}
	source('functions.R')
	
	
	matching.matrix = matrix(FALSE, nrow = nrow(T), ncol = ncol(score.PWM.RBP.per.gene)  + length(unique(miRNA$miRNA.name)))
	rownames(matching.matrix) = T$gene; colnames(matching.matrix) = c(colnames(score.PWM.RBP.per.gene) , sort(unique(miRNA$miRNA.name)))
	rbps = 1:ncol(score.PWM.RBP.per.gene); mis = (ncol(score.PWM.RBP.per.gene)+1): ncol(matching.matrix)
	
	for(rbp in rbps){
		score.max = best.scores[rbp]
		cat(score.max,'\n')
		matching.matrix[,rbp] = score.PWM.RBP.per.gene[,rbp] >= (1/(score.max^2/2000+1))* score.max
		}
		
	for(mi in mis){
		cat(mi,'\t',colnames(matching.matrix)[mi] ,'\n')
		miRNA.name = colnames(matching.matrix)[mi]
		
		j = which(miRNA$miRNA.name ==  miRNA.name)
		m = match(T$gene, miRNA$gene[j])
		matching.matrix[,mi] = !is.na(m)
		}
		
	save(matching.matrix, file = 'Rdata/matching_matrix_gene_RBP_miRNA.Rdata')
	
#	n.gene = apply(matching.matrix,2,sum, na.rm = TRUE)
#	plot(n.gene)
#	plot(n.gene[rbps] , best.scores)

}




if(do.step){
	cat("GLMNET\n")
	if(load){load('Rdata/gene_fit_results_and_LRT_name_cleaned.Rdata'); load('Rdata/RBP_scanning_score.Rdata'); load('Rdata/miRNA.Rdata'); load('Rdata/matching_matrix_gene_RBP_miRNA.Rdata')}
	source('functions.R')
	
	best.model  = T$LRT.best.model
	
	
	require(glmnet)
	require(gplots)
	
	j3 = (best.model == 3)&!is.na(best.model)
	j4 = (best.model == 4)&!is.na(best.model)

	y1 = c(T$eps.gamma.m3[j3]*cos(2*pi*T$phase.gamma.m3[j3]/24), T$eps.gamma.m4[j4]*cos(2*pi*T$phase.gamma.m4[j4]/24))
	y2 = c(T$eps.gamma.m3[j3]*sin(2*pi*T$phase.gamma.m3[j3]/24), T$eps.gamma.m4[j4]*sin(2*pi*T$phase.gamma.m4[j4]/24))
	y=cbind(y1,y2)
	
	x = matching.matrix[c(which(j3),c(which(j4))),]+1-1

	fit=glmnet(x,y,alpha=0.1,family='mgaussian', type.multinomial=c("grouped"), standardize=F, standardize.response=F)



	x = matrix(c(1,0,0,0,2,0,0,0,1),3,3)
	y = c(1,0,2)
	fit=glmnet(x,y,alpha=0.1,family='gaussian', type.multinomial=c("grouped"), standardize=F, standardize.response=F)

	
	
	
	
}






#### DEADENYLASES

D = read.table('../complementary_data/deadenylases.txt', stringsAsFactors = FALSE)
colnames(D) = 'gene'
m  = match( T$gene, D$gene)
mi = which(!is.na(m))
T[mi,]



pdf('plots/deadenylases.pdf')


cols = rainbow(n = length(mi), s = 0.8, v = 0.85)

matplot(zt,t(T[mi,ZT.ex]), type = 'l', col = cols, axes = FALSE, main = 'Deadenylases - exon', lwd = 2)
axis(1, at = seq(0,48,by = 4))
axis(2)
box()
legend('topleft', legend = T$gene[mi], col = cols, lty = c(1:5))

matplot(zt,t(T[mi,ZT.int]), type = 'l', col = cols, axes = FALSE, main = 'Deadenylases - intron', lwd = 2)
axis(1, at = seq(0,48,by = 4))
axis(2)
box()



for(i in 1:length(mi)){
	
	ex = T[mi[i],ZT.ex]; int =  T[mi[i],ZT.int]
	plot(zt, int, col = 'green3', type = 'b', pch = 5, lwd = 2, lty = 2, main = T$gene[mi[i]], axes = FALSE, ylim = range(ex, int))
	points(zt, ex, col = 'steelblue', type = 'b', pch = 16 , lwd = 2, lty = 1)
	axis(1, at = seq(0,48,by = 4))
	axis(2)
	box()
	}
dev.off()




	annot2 = getBM(attributes = c('ensembl_transcript_id', '3_utr_start','3_utr_end'), filters = 'wikigene_name', values = gene_list, mart = ensembl) # other possible attributes : 'external_transcript_id','mgi_symbol','wikigene_name',
	colnames(annot2) = c('ensembl_ID', 'start_3UTR','end_3UTR')
	
	k = which(is.na(annot2$start_3UTR))
	annot2 = annot2[-k,]
	m = match(annot2$ensembl_ID, annot$ensembl_ID)
	k = which(is.na(m));
	if(length(k)>0){annot2 = annot2[-k,]}
	m = match(annot2$ensembl_ID, annot$ensembl_ID)

	
	length1 = nchar(annot$sequence);
	length2 = annot2$end_3UTR - annot2$start_3UTR
	
	plot(length1[m], length2)
	abline(a = 0, b = 1)

























