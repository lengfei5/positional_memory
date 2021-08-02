library(emdbook)
library(deSolve)
library(fdrtool)
library(circular)
library(preprocessCore)
library(gtools)
library(biomaRt)
#library(plotrix)


transparent.factor = 0.3
transparent.gray = rgb(col2rgb('gray')[1]/255,col2rgb('gray')[2]/255,col2rgb('gray')[3]/255, transparent.factor)
transparent.green3 = rgb(col2rgb('green3')[1]/255,col2rgb('green3')[2]/255,col2rgb('green3')[3]/255, transparent.factor)
transparent.tomato = rgb(col2rgb('tomato')[1]/255,col2rgb('tomato')[2]/255,col2rgb('tomato')[3]/255, transparent.factor)
transparent.black = rgb(col2rgb('black')[1]/255,col2rgb('black')[2]/255,col2rgb('black')[3]/255, transparent.factor)

palette(c('gray','green3','tomato','black', transparent.gray ,transparent.green3, transparent.tomato, transparent.black))

make.color.scale = function(xlim = c(10,20)){
	
	polygon(x = c(xlim[1],xlim[2]*2,xlim[2]*2,xlim[1]), y = c(0,0,100,100), col= 'white', border= 'white')
	
	pch = 3
	
	xl = c(xlim[1]+0.2*(xlim[2]-xlim[1]), xlim[2]); x1 = xl[1]; x2 = xl[2]
	length = xl[2]-xl[1]
	
	col = col.phases; N = length(col); y.offset = offsets[length(offsets)]; 
	text(seq(x1,x2,len=5),y = rep(y.offset,5), paste('ZT',seq(0,24,len = 5), sep = ''), pos = 3 , cex = 0.8)
	#points(seq(0,1,len=5),y = rep(y.offset +h,5), pch = pch);
	y = c(y.offset-hh,y.offset-hh,y.offset+hh,y.offset+hh);
	for(i in 1:N){x = c((i-1)/N, i/N, i/N, (i-1)/N)*length+x1;  polygon(x,y,col = col[i], border = col[i])}; polygon(c(0,1,1,0),y, col = 'transparent', border = 'black' , lwd = 0.5)

	y.min = 5; y.max = 27; text(rep(mean(xl),2), c(y.min+0.1,y.max-0.3), c('min','max'),pos = 3, cex = 1 , offset = 0); x.ex = c(0,0.3,0.3,0)*length+x1; x.int = c(0.7,1,1,0.7)*length+x1; 
	for(i in 1:Ncolors){
		y = c(y.min + 22*(i-1)/Ncolors, y.min +22*(i-1)/Ncolors, y.min + 22*(i)/Ncolors, y.min + 22*(i)/Ncolors);
		polygon(x.ex, y, col = col.rel.ex[i], border = col.rel.ex[i] )
		polygon(x.int, y, col = col.rel.int[i], border = col.rel.int[i] )
		}; text(c(0.15,0.85)*length+x1, rep(y.max,2), c('exon', 'intron'),pos = 3)
#
	col = col.pval; N = length(col); y.offset = 3;
	y = c(y.offset-hh,y.offset-hh,y.offset+hh,y.offset+hh);
	for(i in 1:N){x = c((i-1)/N, i/N, i/N, (i-1)/N)*length/2+x1+length/4; polygon(x,y,col = col[N-i+1], border = col[N-i+1])}; polygon(c(0,1,1,0),y, col = 'transparent', border = 'black', lwd = 0.5 )
	#text(seq(x1+length/2/N,x2-length/2/N,len=5),y = rep(y.offset,5), seq(1,0,len = 5), cex = 0.8, col=  'white')
	text(x1,y.offset,'high',pos = 4, offset =0); text(x2,y.offset,'low',pos = 2, offset =0); 
#	
	col = col.median.level; N = length(col); y.offset = 2;
	y = c(y.offset-hh,y.offset-hh,y.offset+hh,y.offset+hh);
	for(i in 1:N){x = c((i-1)/N, i/N, i/N, (i-1)/N)*length/2+x1+length/4; polygon(x,y,col = col[i], border = col[i])}; polygon(c(0,1,1,0),y, col = 'transparent', border = 'black' , lwd = 0.5)
	#text(seq(x1+length/2/N,x2-length/2/N,len=5),y = rep(y.offset,5), seq(3,15,len = 5), cex = 0.8, col=  'white')
	text(x1,y.offset,'low',pos = 4, offset =0); text(x2,y.offset,'high',pos = 2, offset =0); 
#	
	y.offset = 1; y = c(y.offset-hh,y.offset-hh,y.offset+hh,y.offset+hh);
	polygon(c(0,0.4,0.4,0)*length+x1, y, col = col.kept[1], border = 'NA'); polygon(c(0.6,1,1,0.6)*length+x1, y, col = col.kept[2], border = 'NA');
	text(c(0.2,0.8)*length+x1, rep(y.offset,2),c('accepted', 'rejected'), pos = 1, cex = 1) 
#		
	}
	
	
attribute.global.variable.colors = function(){
	col.ex <<- 'steelblue'; col.int <<- 'green3'; col.accepted.ex <<- 'steelblue2'; col.accepted.int <<- 'limegreen'; col.rejected <<- 'gray' ; zt <<- seq(0,46,by=2)
	col.phases <<- rainbow(n = 241,s = 0.9, v = 0.9); 
	col.deg <<- 'orangered'
	Ncolors <<- 21; col.rel.ex <<- col.ramp.ex(n = Ncolors); col.rel.int <<- col.ramp.int(n = Ncolors); 
	col.pval <<- col.ramp.pval(n = 11); col.median.level  <<- col.ramp.median(n = 15);col.kept <<- c('yellowgreen','tomato1')
	}	

plot.this.gene.full = function(pdf.folder = '',gene = 'Dbp', ylim.rel.type = c('fixed','adjusted','personalized'), ylim.rel = c(0,2), log2.rel = FALSE, with.fit = FALSE, fitting.method = c('MCMC','Optim','both'), parametrization = c('sigmoid','cosine')){
	ylim.rel.type = ylim.rel.type[1]; parametrization = parametrization[1]; fitting.method = fitting.method[1]
	
	#global variable
	ylim.default <<- c(0,2.2); ylim.default.log2 <<- c(log2(0.1),log2(4))
	attribute.global.variable.colors()
	
	# create pdf and plots
	pdfname = paste(pdf.folder,'/',gene,'.pdf',sep = '')
	nrow = 4; ncol = 6; dim = 3
	pdf(pdfname, width = ncol/3*dim *1.7 , height = nrow* dim)
	if(fitting.method == 'both'){layout(matrix(c(rep(1,6),rep(2:5,each = 3), rep(6:8,each=2)),nrow = nrow,ncol = ncol, byrow = TRUE), widths = rep(1,ncol), heights = c(1.5,1,1,1))}
	else{layout(matrix(c(rep(1,6),rep(2:7,each = 3)),nrow = nrow,ncol = ncol, byrow = TRUE), widths = rep(1,ncol), heights = c(1.5,1,1,1))}
	par(mar = c(0.5, 4, 2, 0.5) + 0.1, mgp = c(2, 0, 0)); make.UCSC.plot( gene = gene) ; cat('\tUCSC plot done \n');
	par(mar = c(3, 3, 2, 0.5) + 0.1, mgp = c(2, 0.7, 0)); show.all.probesets(gene = gene, ylim.rel.type = ylim.rel.type, ylim.rel = ylim.rel, log2.rel = log2.rel); cat('\tprobesets plots done \n')
	show.summarized.signal(gene = gene, ylim.rel.type = ylim.rel.type, ylim.rel = ylim.rel, log2.rel = log2.rel, with.fit = with.fit, fitting.method = fitting.method, parametrization = parametrization); cat('\tsummary plots done \n')
	dev.off()
}




make.UCSC.plot = function( gene = gene){
	i = t$gene == gene; if(sum(i) == 0){cat('ERROR : no gene with this name found in the table');return()}
	rf = which(RefSeq$name2 == gene); 
	if(length(rf) == 0){ gene.name = unlist(strsplit(gene,'_chr'))[1]; rf = which(RefSeq$name2 == gene); }
	if(length(rf) == 0){ gene.accession.numbers = setdiff(unique(t$mRNA.acc.num[i]),'---'); gene.accession.numbers = gsub(' ','', gene.accession.numbers);if(length(gene.accession.numbers)>0){rf = which(!is.na(match(RefSeq$name, gene.accession.numbers)))}}
	if(length(rf) == 0){ cat('WARNING : no gene with this name found in RefSeq\n')}#; plot(1,1,type = 'n', axes = FALSE, xlab = '', ylab = '');text(1,1,'no gene with this name in RefSeq');return()}
	if(length(rf)>0){xlim.gene = range(RefSeq$txStart[rf],RefSeq$txEnd[rf])}else{xlim.gene = c(min(t$start[i], na.rm = TRUE), max(t$stop[i], na.rm = TRUE))}
	gene.length = xlim.gene[2]-xlim.gene[1]; xlim = xlim.gene + c(0,0.2*(gene.length)); xlim.scale = c(xlim.gene[2]+0.05*gene.length, xlim[2])

	plot(1,1, type = 'n', xlab = '', ylab = '', axes = FALSE, xlim = xlim, ylim = c(-length(rf)-1,31), main = gene); 
	#text(1,1, 'UCSC-like plot')
	offsets <<- c(1:3,5:28,30); hh <<- 0.4
	for(h in offsets){abline(h = h, col = 'gray80')}
	for(p in 1:sum(i)){ip = which(i)[p]; cols = rep('transparent',16); 
		cols[1] = col.kept[-t$kept[ip]+2]; 
		if(t$kept[ip]){ 
			cols[2] = col.median.level[round(log2(t$mean[ip]))+1]
			cols[3] = col.pval[round(10*t$pval[ip]+1)]
			timepoint.index.value.no.cut = unlist(round(log10(t[ip,grep('.rel.ampl',colnames(t))])*(Ncolors-1)+ceiling(Ncolors/2)))
			timepoint.index.value = pmax(pmin(timepoint.index.value.no.cut, Ncolors),1)
			if(t$ex_in[ip] == 1){cols[27:4] = col.rel.ex[timepoint.index.value]}else{cols[28:5] = col.rel.int[timepoint.index.value]} 
			cols[28] = col.phases[round(10*t$phase[ip])+1]}
		plot.polygons(index.probe = ip , cols = cols)
	}
	if(length(rf)>0){plot.gene.pattern(index = rf)}
	make.color.scale(xlim = xlim.scale)
	axis(2,las = 1, at = offsets, labels = c('kept','mean exp.','F24 pval',paste('ZT',rev(zt[6:24]),sep = ''),paste('ZT0',rev(zt[1:5]),sep = ''),'phase'), tick = FALSE)
	}

col.ramp.ex <- colorRampPalette(c('steelblue4', 'steelblue1' , 'white' , 'orangered' , 'orangered3'), space = "Lab")
col.ramp.int <- colorRampPalette(c('darkturquoise', 'turquoise' , 'white' , 'gold1' , 'darkgoldenrod2'), space = "Lab")
col.ramp.pval <- colorRampPalette(c('black', 'blue' , 'royalblue' , 'steelblue1' , 'lightsteelblue1', 'white'), space = "Lab")
col.ramp.median <- colorRampPalette(c('gray90', 'firebrick1' , 'black'), space = "Lab")




	
plot.polygons = function(index.probe = 1, cols = cols){
	x0 = t$start[index.probe]; x1 = t$stop[index.probe]
	for(i in 1:length(cols)){ x = c(x0,x1,x1,x0); y = c(offsets[i]-hh,offsets[i]-hh,offsets[i]+hh,offsets[i]+hh); polygon(x,y, border = NA, col = cols[i])}
}


plot.gene.pattern = function(index = 1){
	for(i in 1:length(index)){plot.transcript.pattern(index = index[i], hline = -i)}
	}

plot.transcript.pattern = function(index = 1, hline = 0){
	if(RefSeq$strand[index] == '+'){pch = '>'}else{pch = '<'}
	points(x = c(RefSeq$txStart[index], RefSeq$txEnd[index]), y = rep(hline,2), type = 'l', lwd = 1); points(x = seq(RefSeq$txStart[index], RefSeq$txEnd[index], len = 11)[2:10], y = rep(hline,9), type = 'p', pch = pch, cex = 1.3); 
	exon.starts = unlist(strsplit(RefSeq$exonStarts[index],','));exon.ends = unlist(strsplit(RefSeq$exonEnds[index],','));
	for(i in 1:RefSeq$exonCount[index]){x = c(exon.starts[i], exon.ends[i], exon.ends[i],exon.starts[i]); y = c(hline-0.3, hline-0.3, hline +0.3, hline +0.3);polygon(x,y, border = NA, col = 'black')}
	}
	

show.all.probesets = function(gene = gene, ylim.rel.type = c('fixed','adjusted','personalized'),ylim.rel = c(0,2), log2.rel = FALSE){
	ylim.rel.type = ylim.rel.type[1]
	show.all.probesests.detail(gene = gene, probe.type = 'intron', plot.rejected = TRUE , signal.type = 'absolute', ylim.rel.type = ylim.rel.type, ylim.rel = ylim.rel, log2.rel = log2.rel)#; cat('1st probesets plots done \n')
	show.all.probesests.detail(gene = gene, probe.type = 'exon'  , plot.rejected = TRUE , signal.type = 'absolute', ylim.rel.type = ylim.rel.type, ylim.rel = ylim.rel, log2.rel = log2.rel)#; cat('2nd probesets plots done \n')
	show.all.probesests.detail(gene = gene, probe.type = 'intron', plot.rejected = FALSE, signal.type = 'relative', ylim.rel.type = ylim.rel.type, ylim.rel = ylim.rel, log2.rel = log2.rel)#; cat('3rd probesets plots done \n')
	show.all.probesests.detail(gene = gene, probe.type = 'exon'  , plot.rejected = FALSE, signal.type = 'relative', ylim.rel.type = ylim.rel.type, ylim.rel = ylim.rel, log2.rel = log2.rel)#; cat('4th probesets plots done \n')

	}	


show.all.probesests.detail = function(gene = gene, probe.type = 'intron', plot.rejected = TRUE, signal.type = c('absolute', 'relative'), ylim.rel.type = c('fixed','adjusted','personalized'), ylim.rel = c(0,2), log2.rel = FALSE){
	signal.type = signal.type[1];ylim.rel.type = ylim.rel.type[1]; #zt = seq(0,22,by = 2)
	
	i = t$gene == gene; if(sum(i) == 0){cat('ERROR : no gene with this name found in the table');return()}
	if(probe.type == 'intron'){i = i & (t$ex_in == -1); col.accepted = col.accepted.int}else{i = i & (t$ex_in == 1); col.accepted = col.accepted.ex}; if(!plot.rejected){i  = i & t$kept}
	
	
	if(signal.type == 'absolute'){
		ZT = c(1:24); ylim = c(0,15); if(sum(i)>0){MAT = log2(t[i,ZT])}else{MAT = rep(0,12)}; hline = -10; ylab = 'absolute signal [log2 scale]'}else{
		ZT = grep('.rel.ampl', colnames(t)); if(log2.rel){if(sum(i)>0){MAT = log2(t[i,ZT])}else{MAT = rep(0,12)}; hline = 0; ylab = 'relative signal [log2 scale]'}else{if(sum(i)>0){MAT = t[i,ZT]}else{MAT = rep(0,12)}; hline = 1; ylab = 'relative signal [lin scale]'}
		if(ylim.rel.type == 'fixed'){ylim = ylim.default; if(log2.rel){ylim = ylim.default.log2}}else if(ylim.rel.type == 'personalized'){ylim = ylim.rel}else{ylim = if(sum(i)>0){range(MAT)}else{ylim = ylim.default; if(log2.rel){ylim = ylim.default.log2}}}}
	
	
	
	style = (!t$kept[i]) + 1; if(ylim.rel.type == 'adjusted'){}
	palette(c(col.accepted, col.rejected))
	
	plot(1,1,type = 'n', axes = FALSE, xlab = '', ylab = ylab, main = probe.type,ylim = ylim, xlim = range(zt))
	if(sum(i)>0){
		matplot(zt,t(MAT), type = 'l',pch = 16, col = style, lwd = -style+3, lty = style,  xlab = '', ylab = ylab, main = probe.type, add = TRUE) 
		matplot(zt,t(MAT), type = 'p',pch = 16, col = style, cex = 0.6,add = TRUE)
	}	
	abline(h = hline, col = 'gray', lty = 3)
	time.axis()
	axis(2); box()
}




show.summarized.signal = function(gene = gene, ylim.rel.type = c('fixed','adjusted','personalized'), ylim.rel = c(0,2), log2.rel = FALSE, with.fit = FALSE, fitting.method = 'MCMC', parametrization = c('sigmoid','cosine'), abline = FALSE, dt = 4, title = 'default'){
	ylim.rel.type = ylim.rel.type[1]; parametrization = parametrization[1];
	show.summarized.signal.detail(gene = gene, ylim.rel.type = ylim.rel.type, ylim.rel = ylim.rel, log2.rel = log2.rel, with.cosine.fit = FALSE, with.fit = FALSE, fitting.method = 'MCMC', parametrization = parametrization, abline = abline, dt = dt, title = title[1])
	with.cosine.fit = FALSE
	#if(with.fit){with.cosine.fit = FALSE}else{with.cosine.fit = TRUE}
	if(length(title)>1){title = title[2]}
	if(fitting.method == 'both'){
		cat('\n\n\n HERE\n\n\n')
		show.summarized.signal.detail(gene = gene, ylim.rel.type = ylim.rel.type, ylim.rel = ylim.rel, log2.rel = log2.rel, with.cosine.fit = with.cosine.fit, with.fit = with.fit, fitting.method = 'Optim', parametrization = parametrization, abline = abline, dt = dt, title = title)
		show.summarized.signal.detail(gene = gene, ylim.rel.type = ylim.rel.type, ylim.rel = ylim.rel, log2.rel = log2.rel, with.cosine.fit = with.cosine.fit, with.fit = with.fit, fitting.method = 'MCMC', parametrization = parametrization, abline = abline, dt = dt, title = title)
	}else{
		show.summarized.signal.detail(gene = gene, ylim.rel.type = ylim.rel.type, ylim.rel = ylim.rel, log2.rel = log2.rel, with.cosine.fit = with.cosine.fit, with.fit = with.fit, fitting.method = fitting.method, parametrization = parametrization, abline = abline, dt = dt, title = title)
	}
}


show.summarized.signal.detail = function(gene = gene, ylim.rel.type = c('fixed','adjusted','personalized'), ylim.rel = c(0,2), log2.rel = FALSE, with.cosine.fit = FALSE, with.fit = FALSE, fitting.method = 'MCMC', parametrization = c('sigmoid','cosine'), abline = FALSE, dt = 4, title = 'default'){
	ylim.rel.type = ylim.rel.type[1]; parametrization = parametrization[1]
	i = T$gene == gene; 
	if(sum(i)>0){
		ex = unlist(T[i,grep('.rel.ampl.ex', colnames(T))]); int = unlist(T[i,grep('.rel.ampl.int', colnames(T))]);hline = 1; ylab = 'relative signal [lin scale]'
		if(log2.rel){ex = log2(ex);hline = 0; ylab = 'relative signal [log2 scale]'}; 
	
		if(ylim.rel.type == 'fixed'){ylim = ylim.default; if(log2.rel){ylim = ylim.default.log2}}else if(ylim.rel.type == 'adjusted'){ylim = range(ex,int)}else{ylim = ylim.rel}
		if(title == 'default'){title = 'summarized signal'; if(with.cosine.fit){title = paste(title,'+ cos. fits (2 comp.)')}; if(with.fit){}}
		plot(1,1, type = 'n', xlab = 'time (hours)', ylab = ylab, axes = FALSE, xlim = range(zt), ylim = ylim); box(); 
		time.axis(dt = dt); axis(2)
		if(abline){abline(h = 1, col = 'gray')}
		points(zt, ex, col = col.ex, type = 'p', lwd = 2, pch = 16); points(zt, int, col = col.int, type = 'p', lwd = 1, pch = 5)
		if(with.cosine.fit){
			zt.p = seq(0,48,by = 0.05)
			cos.ex = T$mean.ex[i] + min(1,T$rel.ampl.ex[i]) * cos(2*pi*(zt.p-T$phase.ex[i])/24) +  min(1,T$rel.ampl.12.ex[i]) * cos(2*pi*(zt.p-T$phase.12.ex[i])/12)
			cos.int = T$mean.int[i] + min(1,T$rel.ampl.int[i]) * cos(2*pi*(zt.p-T$phase.int[i])/24) +  min(1,T$rel.ampl.12.int[i]) * cos(2*pi*(zt.p-T$phase.12.int[i])/12)
			points(zt.p, cos.ex, type = 'l', lwd = 2, col = col.ex);points(zt.p, cos.int, type = 'l', lwd = 2, col = col.int)
			points(zt, ex, col = col.ex, type = 'l', lwd = 1);points(zt, int, col = col.int, type = 'l', lwd = 0.5)
		}else if(with.fit){
			zt.p = seq(0,48,by = 0.05)
			
			gamma = 0; eps.gamma = 0; phase.gamma = 0; 
			rel.ampl.int = 0; phase.int = 0; rel.ampl.12.int = 0; phase.12.int = 0
			fold.change = 1; up.time = 12; down.time = 12;

			if((fitting.method == 'MCMC')&!is.null(T$model.MCMC[i]) && !is.na(T$model.MCMC[i])){ 
				best.model = T$model.MCMC[i]
				gamma = T$gamma.MCMC[i]; eps.gamma = T$eps.gamma.MCMC[i]; phase.gamma = T$phase.gamma.MCMC[i]; 
				rel.ampl.int = T$rel.ampl.int.MCMC[i]; phase.int = T$phase.int.MCMC[i]; rel.ampl.12.int = T$rel.ampl.12.int.MCMC[i]; phase.12.int = T$phase.12.int.MCMC[i]
				if(title == 'default'){title = paste(title,'-',fitting.method)}
				if(title == 'model'){title = c('CS-CD','RS-CD','CS-RD','RS-RD')[best.model]}
			}else if(!is.na(T$LRT.best.model[i])){
				best.model = T$LRT.best.model[i]; if(!is.null(T$LRT.eff.best.model[i])){best.model = T$LRT.eff.best.model[i]}
				eval(parse(text =paste('gamma = T$gamma.m',best.model,'[i]',sep = '')));eval(parse(text =paste('eps.gamma = T$eps.gamma.m',best.model,'[i]',sep = '')));eval(parse(text =paste('phase.gamma = T$phase.gamma.m',best.model,'[i]',sep = '')));
				if(parametrization == 'cosine'){
				eval(parse(text = paste('rel.ampl.int = T$rel.ampl.int.m',best.model,'[i]; phase.int = T$phase.int.m',best.model,'[i]; rel.ampl.12.int = T$rel.ampl.12.int.m',best.model,'[i]; phase.12.int = T$phase.12.int.m',best.model,'[i]',sep = '')))}
				else{eval(parse(text = paste( 'fold.change = T$fold.change.int.m',best.model,'[i]; phase.int = T$phase.int.m',best.model,'[i]; up.time = T$up.time.int.m', best.model,'[i]; down.time = T$down.time.int.m', best.model,'[i]', sep = '')))}
				if(title == 'default'){title = paste(title,'- Optim')}
				if(title == 'model'){title = c('CS-CD','RS-CD','CS-RD','RS-RD')[best.model]}
			}
			
			cat('is.null(fold.change) :', is.null(fold.change),'\n')
			
			if(is.null(fold.change)){fold.change = 1; phase.int = 0; up.time = 12; down.time = 12}
			if(is.null(eps.gamma)){if(parametrization == 'cosine'){eps.gamma = 0; phase.gamma = 0}else{eps.gamma = 1; phase.gamma = 0}}

			
			
			cat('parametrization :', parametrization, '\n gamma = ',gamma,'\n eps.gamma = ',eps.gamma, '\n phase.gamma = ',phase.gamma,'\n fold.change = ',fold.change, '\n phase.int = ',phase.int, '\n up.time = ',up.time, '\n down.time = ', down.time, '\n' );
			real.gamma = gamma
			gamma = min(gamma, 5)
			
			profile.int = rep(1,length(zt.p)); 
			if(best.model > 1){
				if(parametrization == 'cosine'){profile.int = compute.s(t = zt.p, eps.24.S = rel.ampl.int, phase.24.S = phase.int, eps.12.S = rel.ampl.12.int, phase.12.S = phase.12.int  )}
				else{profile.int = compute.sigmoid(t = zt.p, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)}
			}
			
			profile.deg = rep(1,length(zt.p)); 
			if(best.model > 2){
				if(parametrization == 'cosine'){profile.deg = cos.deg + eps.gamma * cos(2*pi*(zt.p-phase.gamma)/24)}
				else{profile.deg = compute.sigmoid(t = zt.p, fold.change = eps.gamma, phase = phase.gamma, up.time = 12, down.time = 12)}
			} 
			
			
			if(best.model == 1){profile.ex = rep(1,length(zt.p))}
			if(best.model > 1){
				if(parametrization == 'cosine'){profile.ex = compute.m.changing(t = zt.p, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, eps.24.S = rel.ampl.int, phase.24.S = phase.int, eps.12.S = rel.ampl.12.int, phase.12.S = phase.12.int)}
				else{profile.ex = compute.m.sigmoid(t = zt.p, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)}
			}
			
			#cat('length(zt.p) = ',length(zt.p), '\t length(profile.ex) = ',length(profile.ex),'\t length(profile.int) = ',length(profile.int),'\t length(profile.deg) = ',length(profile.deg),'\n' )

			points(zt.p, profile.ex, type = 'l', lwd = 2, col = col.ex);points(zt.p, profile.int, type = 'l', lwd = 2, col = col.int); points(zt.p, profile.deg, type = 'l', lwd = 2, col = col.deg);
			
			
			if(best.model>1){legend('topright',legend = paste('half-life =',round(log(2)/real.gamma, digits = 2),'h'), text.col = col.deg, bty = 'n', lty = 0)}
			if(title == 'default'){legend('topleft',legend = paste('Model =',c('CS-CD','RS-CD','CS-RD','RS-RD')[best.model]), text.col = c('gray','green4','tomato','black')[best.model], bty = 'n', lty = 0)}
			points(zt, ex, col = col.ex, type = 'l', lwd = 1);points(zt, int, col = col.int, type = 'l', lwd = 0.5)
			
		}else{points(zt, ex, col = col.ex, type = 'l', lwd = 2);points(zt, int, col = col.int, type = 'l', lwd = 2)}
	}else{plot(1,1,type = 'n', axes = FALSE, xlab = '', ylab = ''); text(1,1, 'no summarized signal for this gene')}
	col.title = 'black'; if(any(title == c('CS-CD','RS-CD','CS-RD','RS-RD'))){col.title = best.model}
	mtext(title, side = 3, line = 1, cex = 0.8, font = 2, col = col.title)
	
}


time.axis = function(dt = 4){axis(1, at = zt, labels = rep('',length(zt))); axis(1, at = seq(0,48,by = dt))}



################################################
################################################
################################################

make.gif = function(gene = 'Dbp', gif.folder = ''){
	
	gif.filename = paste(gif.folder,'/',gene,'.gif', sep = '')
	zt <<- seq(0,46,by = 2)
	png.folder = paste(gif.folder,'/',gene,'_png', sep = '');dir.create(png.folder)
	for(i in 1:length(zt)){
	#for(i in 1:1){	
		ZT = as.character(zt[i]); if(zt[i]<10){ZT = paste('0',ZT,sep = '')}
		png.filename = paste(png.folder,'/',gene,'_ZT',ZT,'.png', sep = '')
		png(filename = png.filename, width = 1500, height = 700, bg = 'white')
		layout(matrix(c(1,2),nrow = 1, ncol = 2), widths = c(4,1))
		par(mar = c(1, 6, 1, 1) + 0.1, mgp = c(4,1,0))
		plot.time.point(ZT = ZT, gene = gene)
		par(mar = c(1, 1, 1, 1) + 0.1, mgp = c(4,1,0))
		plot.clock(time = zt[i], gene = gene)
		dev.off()
		}
	cmd = paste('convert -delay 40 ',png.folder,'/*.png ',gif.filename,sep='')
	system(cmd)
		
}


plot.time.point = function(ZT = '00', gene = 'Dbp'){
	i = (t$gene == gene)&(t$kept); if(sum(i) == 0){cat('ERROR : no gene with this name found in the table');return()}
	rf = which(RefSeq$name2 == gene); 
	if(length(rf) == 0){ gene.name = unlist(strsplit(gene,'_chr'))[1]; rf = which(RefSeq$name2 == gene); }
	if(length(rf) == 0){ gene.accession.numbers = setdiff(unique(t$mRNA.acc.num[i]),'---'); gene.accession.numbers = gsub(' ','', gene.accession.numbers);if(length(gene.accession.numbers)>0){rf = which(!is.na(match(RefSeq$name, gene.accession.numbers)))}}
	if(length(rf) == 0){ 'ERROR : no gene with this name found in RefSeq'; return()}
	xlim.gene = range(RefSeq$txStart[rf],RefSeq$txEnd[rf]); #xlim = xlim.gene + c(0,0.2*(xlim.gene[2]-xlim.gene[1]))
	plot(1,1, type = 'n', xlab = '', ylab = 'signal (log2 scale)', axes = FALSE, xlim = xlim.gene, ylim = c(-length(rf)-1,15), cex.lab = 2);
	text(mean(xlim.gene), 15, gene, cex = 2, font = 4, pos = 1)
	for(p in 1:sum(i)){
		ip = which(i)[p]	
		x = c(t$start[ip],t$stop[ip],t$stop[ip],t$start[ip]); y.max = log2(c(1,1,rep(max(t[ip,2:13]),2))); if(t$ex_in[ip] == 1){col = 'steelblue'}else{col = 'green3'}
		polygon(x,y.max, col = 'transparent', border = col, lwd = 2)
		y = log2(c(1,1,rep(t[ip,paste('ZT',ZT,sep= '')],2)));
		polygon(x,y, col = col, border = col, lwd = 2)
		}
	
	plot.gene.pattern(index = rf)
	axis(2, las = 1, at = seq(0,15,by = 3), cex.axis = 2)
}


plot.clock =  function(time = 0, gene = 'Dbp'){
	i = (t$gene == gene)&t$kept; if(sum(i) == 0){cat('ERROR : no gene with this name found in the table');return()} 
	phases = t$phase[i]; colors = c('green3', 'steelblue'); col = colors[(t$ex_in+1)/2+1]
	radius = 1; lim = c(-1,1)*radius*1.2
	plot(0,0, type = 'n', axes= FALSE, xlab = '', ylab = '', asp = 1, xlim = lim, ylim = lim)
	matplot(rbind(rep(0,sum(i)),-radius*cos(2*pi*phases/24 + pi/2)), rbind(rep(0,sum(i)),radius*sin(2*pi* phases/24 + pi/2)), type = 'l', lty = 1, col = col, lwd = 2, add = TRUE)
	col = 'orange'; if(time>22){col = 'red'};points(c(0,-radius*cos(2*pi*time/24 + pi/2)), c(0,radius*sin(2*pi*time/24 + pi/2)), type = 'l', col = col, lwd = 5)
	draw.circle(x = 0, y = 0, radius = radius, lwd = 5)
	text(-1.2*radius*cos(2*pi*zt[1:12]/24 + pi/2),1.2*radius*sin(2*pi*zt[1:12]/24 + pi/2),zt[1:12], cex = 2)
	day  = 1; if(time>22){day = 2}; text(0,1.4,paste('day',day), col = col, pos = 3, cex = 2, font = 2)
	}



################################################
################################################
################################################

do.rhythmic.analysis = function(T = T1, ZT = ZT.ex, delta.t = 2)
{
	FFT = t(apply(T[,ZT],1,fft))
	n = length(ZT); 
	if(delta.t == 2){if(n == 12){i.24 = 2; i.12 = 3}else if(n == 24){i.24 = 3; i.12 = 5}}
	if(delta.t == 4){if(n == 6){i.24 = 2; i.12 = 3}else if(n == 12){i.24 = 3; i.12 = 5}}
	mean = abs(FFT[,1])/n 
	RelAmp = 2*abs(FFT[,i.24])/n/mean 
	RelAmp.12 = 2*abs(FFT[,i.12])/n/mean 
	phase = apply(T[,ZT],1,phase.n, n = i.24-1)*24/2/pi
	phase.12 =  apply(T[,ZT],1,phase.n, n = i.12-1)*12/2/pi
	score = apply(T[,ZT],1,fs.n, n = i.24-1)
	pval = apply(T[,ZT],1,pval.n, n = i.24 -1)
	pval.12 = apply(T[,ZT],1,pval.n, n = i.12 -1)
	return(data.frame(mean = mean, rel.ampl = RelAmp, rel.ampl.12 = RelAmp.12, phase = phase, phase.12 = phase.12, fscore = score,pval = pval, pval.12 = pval.12))
}

plot.cosine.fit = function(gene = 'Dbp', replicate = 1, probe.type = 'exon', n.components = 2){
	if(replicate == 1){T= T1}else{T = T2}
	k = which(T$gene == gene); if(length(k) == 0){k = grep(gene, T$gene, ignore.case = TRUE)}; if(length(k) == 0){cat("GENE NOT FOUND! \n"); return()}
	
	if(probe.type == 'exon'){ZT = ZT.ex; mean = T$mean.ex[k]; rel.ampl = T$rel.ampl.ex[k]; phase = T$phase.ex[k]}else{ZT = ZT.int;mean = T$mean.int[k]; rel.ampl = T$rel.ampl.int[k]; phase = T$phase.int[k]}

	
	zt = seq(0,22,by = 2)

	Tr = T[k,ZT]; TR = as.vector(as.numeric(unlist(Tr)))
	FFT = fft(TR)

	sup = seq(0,22,by = 0.05)
	
	cos.fit = mean + rel.ampl*cos(2*pi*(sup-phase)/24)
	cos.fit.2 = mean + min(1,rel.ampl)*cos(2*pi*(sup-phase)/24)
	
	rel.ampl.12 = 2*abs(FFT[3])/12; phase.12 = phase12(TR)*12/2/pi;
	cos.fit.12 = mean + rel.ampl*cos(2*pi*(sup-phase)/24) + rel.ampl.12*cos(4*pi*(sup-phase.12)/24) 
	cos.fit.12.2 = mean + min(1,rel.ampl)*cos(2*pi*(sup-phase)/24) + min(1,rel.ampl.12)*cos(4*pi*(sup-phase.12)/24)
	
	title = paste(gene, '-', probe.type, ' - replicate',replicate)
	plot(1,1, type = 'n', main = title, xlab = 'Time [in hour of a day]', ylab = 'Relative signal, lin scale', xlim = c(0,22), ylim = range(Tr, cos.fit,0,2, na.rm = TRUE), axes = FALSE)
	axis(2); axis(1, at = zt); abline(h = 1, col = 'gray'); abline(h = 0, col = 'red')
	points(zt,Tr, type = 'b', col = 'black', pch = 16, lwd = 2)
	points(sup, cos.fit, type = 'l', col = 'steelblue3', lwd = 3)
	points(sup, cos.fit.2, type = 'l', col = 'darkseagreen3', lwd = 3)
	
	if(n.components>2){
		points(sup, cos.fit.12, type = 'l', col = 'coral', lwd = 3)
		points(sup, cos.fit.12.2, type = 'l', col = 'brown1', lwd = 3)
		}
	
	}
	


fs.n = function(x, n = 1){
	mx = mean(x)
	x = x-mx
	P = (Mod(fft(x)))^2 
	P[n+1]/sum(P[2:floor((length(x)+1)/2)])
	}

	
score_temporal_smoothness = function(x){
	mx = mean(x)
	x = x-mx
	P = (Mod(fft(x)))^2
	L = length(x)/2
	return(sum(P[(L-1):L])/sum(P[1:L]))
	}
	

pval.n = function(x, n = 1){
	S = fs.n(x, n = n)
	pval = (1-S)^(length(x)/2-2)
	return(pval)
	}
	

phase = function(x){
	e = fft(x)[2]
	if(is.na(e)){return(NA)}
	if(Re(e)<0) p = (atan(Im(e)/Re(e))+pi)
	else p = (atan(Im(e)/Re(e))) %% (2*pi)
	return((-p+2*pi)%%(2*pi))
	}
		
	
phase12 = function(x){
	i = 3 # index of the 12h component
	e = fft(x)[i]
	if(is.na(e)){return(NA)}
	if(Re(e)<0) p = (atan(Im(e)/Re(e))+pi)
	else p = (atan(Im(e)/Re(e))) %% (2*pi)
	return((-p+2*pi)%%(2*pi))
	}

phase.n = function(x, n = 1){
	e = fft(x)[n+1]
	if(is.na(e)){return(NA)}
	if(Re(e)<0) p = (atan(Im(e)/Re(e))+pi)
	else p = (atan(Im(e)/Re(e))) %% (2*pi)
	return((-p+2*pi)%%(2*pi))
	}





f24_R2=function(x, t=2*(0:(length(x)-1)), period=24, offset=0){
	n=length(x)
	# mu=mean(x)
	sig2=var(x)

#
	c=cos(2*pi*t/period)
	s=sin(2*pi*t/period)


	# x=x-mean(x)	
	# a=mean(x*(c-mean(c))
	# b=mean(x*(s-mean(s))

	# a=mean(x*c)
	# b=mean(x*s)
	# x.hat=mu+2*a*c+2*b*s
	# sig2.1=var(x-x.hat)
	
	fit = lm(x~c+s)

	a=fit$coef[2]
	b=fit$coef[3]
	
	R2=0
	if(sig2>0) R2 =1.0-sum(fit$residuals^2)/(n-1)/sig2


	p=3


	pv=pbeta(R2, (p-1)/2, (n-p)/2, lower.tail = FALSE, log.p = FALSE)

	lu=length(unique(x))
	amp=max(x)-min(x)
	phase=period/(2*pi)*atan2(b, a)
	if(phase<0) phase=phase+period
	if(phase>period) phase=phase-period
	phase=(phase+offset)%%period
	
	result = c(lu=lu, mean=mean(x), amp=2*sqrt(a^2+b^2), relamp=(max(x)-min(x))/(2.0*fit$coef[1]), phase=phase, pval=pv)
	return(result)
}




fdr2pval = function(P,fdr = 0.1){
	N = length(P)
	sorted_P = sort(P)
	in_fdr = which((sorted_P - c(1:N)*fdr/N)<=0)
	x = length(in_fdr)
	if(x>0){
		in_fdr = in_fdr[length(in_fdr)]; 
		p_th = sorted_P[in_fdr]
	}else{p_th = 1}
	return(p_th)
}



qvals=function(pv){
	PV=pv
	if(any(is.na(PV))) pv=pv[!is.na(PV)]

	r=rank(-pv)
	o=order(-pv)
	lpv=length(pv)
	q=double(lpv)
	
	pv.s=pv[o]

	for (n in 1:lpv){
	        q[n] = pv.s[n] * lpv / (lpv-n+1)
	        if(n==1) qmin=q[1]
	        if(n>1 & q[n]>qmin) {q[n]=qmin}
	        qmin=q[n]
	}
	Q=rep(NA,length(PV))
	if(any(is.na(PV))){Q[!is.na(PV)]=q[r]}else{Q=q[r]}
	return(Q)
}



passes.fdr = function(t = t, fdr = 0.1){
	pval.col = grep('.pval', colnames(t)); X = length(pval.col);
	pval.lim = apply(t[,pval.col],2,fdr2pval,0.1)
	pass = t( t(t[,pval.col])<= pval.lim )
	return(apply(pass,1,sum)>(X/2))
}



remove.low.probesets = function(t = t1, th = 3){
	median = apply(t[t$kept,ZT],1, median)
	ind.small = median>th
	sup = rep(FALSE,nrow(t)); ind.long = sup; ind.long[t$kept] = ind.small 
	return(ind.long)
	}

remove.high.Xhyb = function(t = t1, th = 3){
	return(t$kept & (t$cross.hyb.type<th))
	}
	
enough.probesets = function(t = t1){
	gene_ex = as.character(t$gene[(t$ex_in == 1)& (t$kept)]); gene_list_ex = unique(gene_ex[duplicated(gene_ex)])
	gene_int = as.character(t$gene[(t$ex_in == -1)& (t$kept)]); gene_list_int = unique(gene_int[duplicated(gene_int)])
	gene_list = intersect(gene_list_ex, gene_list_int);#gene_list = c(circ,gene_list[1:200])#; gene_list = c(circ)
	ind = match(t$gene, gene_list); ind = !is.na(ind); ind = ind&t$kept
	return(ind)
	}


#plot.this.gene = function(gene_table = gene_table,gene = 'Dbp'){
#	i = which(gene_table$gene == gene); if(length(i) == 0){i = grep(gene, gene_table$gene, ignore.case = TRUE)}; if(length(i) == 0){cat("GENE NOT FOUND! \n"); return()}
#	
#	EX = unlist(gene_table[i,zt_ex]); INT = unlist(gene_table[i,zt_int])
#	
#	par(mfrow = c(2,2), lwd = 2)
#	
#	
#	ex = EX; int = INT
#	plot(1,1, type = 'n', main = paste(gene,'Absolute Signal, lin scale'), xlab = 'Time [in hour of a day]', ylab = 'Signal',xlim = c(0,22), ylim = range(ex,int), axes = FALSE)
#	axis(2)
#	axis(1,at=zt)
#	points(zt,ex,type = 'b',col = exon_col, pch = 17)
#	points(zt,int,type = 'b',col = intron_col, pch = 16)
#	
#	
#	ex = EX/mean(EX); int = INT/mean(INT)
#	plot(1,1, type = 'n', main = paste(gene,'Relative Signal, lin scale'), xlab = 'Time [in hour of a day]', ylab = 'Relative Amplitude',xlim = c(0,22), ylim = range(ex,int), axes = FALSE)
#	axis(2)
#	axis(1,at=zt)
#	abline(h = 1, col = 'gray')
#	points(zt,ex,type = 'b',col = exon_col, pch = 17)
#	points(zt,int,type = 'b',col = intron_col, pch = 16)
#
#	
#
#	
#	ex = log2(EX); int = log2(INT)
#	plot(1,1, type = 'n', main = paste(gene,'Absolute Signal, log2 scale'), xlab = 'Time [in hour of a day]', ylab = 'Signal (log2 scale)',xlim = c(0,22), ylim = range(ex,int), axes = FALSE)
#	axis(2)
#	axis(1,at=zt)
#	points(zt,ex,type = 'b',col = exon_col, pch = 17)
#	points(zt,int,type = 'b',col = intron_col, pch = 16)
#	
#	
#	ex = ex/mean(ex); int = int/mean(int)
#	plot(1,1, type = 'n', main = paste(gene,'Relative Signal, log2 scale'), xlab = 'Time [in hour of a day]', ylab = 'Relative Amplitude',xlim = c(0,22), ylim = range(ex,int), axes = FALSE)
#	axis(2)
#	axis(1,at=zt)
#	abline(h = 1, col = 'gray')
#	points(zt,ex,type = 'b',col = exon_col, pch = 17)
#	points(zt,int,type = 'b',col = intron_col, pch = 16)
#
#	
#	}	
#	
	
	

plot.duplicates = function(gene = 'Dbp'){
	i = which(t1$gene == gene); if(length(i) == 0){i = grep(gene, t1$gene, ignore.case = TRUE)}; i1 = i;	i = which(t1$gene == gene); if(length(i) == 0){i = grep(gene, t1$gene, ignore.case = TRUE)}; i2 = i; if(length(i1)+length(i2) == 0){cat("GENE NOT FOUND! \n"); return()}

	probesets = unique(t1$probe.set.ID[i1], t2$probe.set.ID[i2])
	n = length(probesets); zt = seq(0,22,by=2); k = 1
	for(probeset in probesets){
		cat('\t',probeset,'\t')
		j1 = which(t1$probe.set.ID == probeset);if(length(j1)>0){s1 = t1[j1,ZT]; r1 = t1[j1,ZT+37]; ex_in = t1$ex_in[j1]}else{s1 = rep(0,length(ZT)); r1 = s1}
		j2 = which(t2$probe.set.ID == probeset);if(length(j2)>0){s2 = t2[j2,ZT]; r2 = t2[j2,ZT+37]; ex_in = t2$ex_in[j2]}else{s2 = rep(0,length(ZT)); r2 = s2}
		if(ex_in == 1){probe_type = 'exonic'}else if(ex_in == -1){probe_type = 'intronic'}else{probe_type = 'other probe type'}
		title = paste(gene,' - ', probeset,' - ',probe_type, ' - ', k, '/', n, sep = '')
				
		plot(1,1, type = 'n', main = title, xlab = 'Time [in hour of a day]', ylab = 'Absolute signal, log2 scale',xlim = c(0,22), ylim = range(0,16), axes = FALSE)
		axis(2, at = seq(0,16,by = 2)); axis(1, at = zt)
		points(zt,log2(s1), type = 'b', col = 'black', pch = 4, lwd = 2)
		points(zt,log2(s2), type = 'b', col = 'steelblue3', pch = 16, lwd = 2)
		
		
		plot(1,1, type = 'n', main = title, xlab = 'Time [in hour of a day]', ylab = 'Absolute signal, lin scale',xlim = c(0,22), ylim = range(s1,s2)*c(0.95,1.05), axes = FALSE)
		axis(2); axis(1, at = zt)
		points(zt,s1, type = 'b', col = 'black', pch = 4, lwd = 2)
		points(zt,s2, type = 'b', col = 'steelblue3', pch = 16, lwd = 2)
		
		cat(unlist(r1),"\n\t\t\t",unlist(r2), "\t")
		range(r1,r2)
		cat("b\n")		
		
		plot(1,1, type = 'n', main = title, xlab = 'Time [in hour of a day]', ylab = 'Relative signal, lin scale',xlim = c(0,22), ylim = range(r1,r2)*c(0.95,1.05), axes = FALSE)
		axis(2); axis(1, at = zt); abline(h = 1, col = 'gray')
		points(zt,r1, type = 'b', col = 'black', pch = 4, lwd = 2)
		points(zt,r2, type = 'b', col = 'steelblue3', pch = 16, lwd = 2)

		
		plot(1,1, type = 'n', main = title, xlab = 'Time [in hour of a day]', ylab = 'Relative signal, lin scale',xlim = c(0,22), ylim = range(0,2), axes = FALSE)
		axis(2); axis(1, at = zt); abline(h = 1, col = 'gray')
		points(zt,r1, type = 'b', col = 'black', pch = 4, lwd = 2)
		points(zt,r2, type = 'b', col = 'steelblue3', pch = 16, lwd = 2)
		k = k+1; cat('done\n')
		}
	}	
	
	
	
remove.useless.line = function(t = t){
	na = is.na(t$ex_in); t = t[!na,] ;  gene_is_na = is.na(t$gene) ; t = t[!gene_is_na,]; is_zero = t$ex_in == 0 ;t = t[!is_zero,] ;
	is_two = t$ex_in == 2; t = t[!is_two,]
	cat('ex_in != NA :', sum(!na),'\n gene != NA :', sum(! gene_is_na),'\n ex_in != 0 :', sum(!is_zero),'\n')
	return(t)
	}	
	
	


signal.selection = function(t = t1, gene_list = gene_list, n.days = 2){

	N = length(gene_list)	
	
	EX_SIGNAL = matrix(NA,N,12*n.days); INT_SIGNAL = EX_SIGNAL;colnames(EX_SIGNAL) = paste(colnames(t)[ZTr],'.ex',sep=''); colnames(INT_SIGNAL) = paste(colnames(t)[ZTr],'.int', sep = ''); 
		
	gene_signal = data.frame(gene = gene_list, EX_SIGNAL, INT_SIGNAL,exon.median = rep(NA,N), exon.var = rep(NA,N), intron.median = rep(NA,N), intron.var = rep(NA,N), stringsAsFactors = FALSE)
	ZT_ex = c(2:(12*n.days+1)); ZT_int = ZT_ex +12*n.days
	
	flag.ex = rep(FALSE,nrow(t)); flag.int =rep(FALSE, nrow(t))

	for(i in 1:N){
		if((i%%140)==0){cat(round(100*i/N, digits = 1),'%\n')}
		gene = gene_list[i]; 
		#cat(gene, '\n')
		i.ex = select(t = t,gene = gene,probe = 'exon')
		i.int = select(t = t,gene = gene,probe = 'intron')

		out = gene.signal.selection(t = t, i.ex = i.ex, i.int = i.int)
		
		
		gene_signal[i,ZT_ex] = out$ex
		gene_signal[i,ZT_int] = out$int
		gene_signal$exon.var[i] = out$var.ex
		gene_signal$intron.var[i] = out$var.int
		gene_signal$exon.median[i] = median(unlist(t[i.ex,ZT][out$flag.ex,]))
		gene_signal$intron.median[i] = median(unlist(t[i.int,ZT][out$flag.int,]))
		flag.ex[i.ex] = out$flag.ex; flag.int[i.int] = out$flag.int
		}
	return(list(T = gene_signal,flag.ex = flag.ex, flag.int = flag.int))
	}



gene.signal.selection = function(t = t1, i.ex = i.ex, i.int = i.int){
	out.ex = rep(NA,12); out.int = out.ex; var.int = NA; var.ex = NA
	n.ex = sum(i.ex, na.rm = TRUE); n.int = sum(i.int, na.rm = TRUE)
	if((n.ex < 2)||(n.int < 2)){return(list(ex = out.ex,int = out.int,var.ex = var.ex,var.int = var.int, flag.ex = rep(FALSE,n.ex),flag.int = rep(FALSE,n.int)))}else{		
		
		#### Selection of the intronic probes
		
		flag.int.cluster = rep(TRUE, n.int)
		flag.t.s = rep(TRUE, n.int)
		if(n.int > 2){
			#clustering
			k = kmeans(log2(t[i.int,ZT]),2);c.high = which.max(apply(k$center,1,mean)); n.high = sum(k$cluster == c.high)
			if(n.high>=round(n.int/2)){
				med.high = apply(t[i.int,ZTr][k$cluster == c.high,],2,median); med.low = apply(t[i.int,ZTr][k$cluster != c.high,],2,median);
				if(cor(med.high, med.low)<0.5){flag.int.cluster[k$cluster != c.high] = FALSE }
			}
			# temporal smoothness
			threshold = 0.3
			scores = apply(t[i.int,ZTr],1,score_temporal_smoothness); j = scores < threshold; k = j
			if(sum(!j)>0){median = apply(t[i.int,ZTr],2,median); cor = apply(t[i.int,ZTr],1,cor,median); k = cor>0.7 }
			flag.t.s = j|k;
		}
		flag.int = flag.int.cluster & flag.t.s; if(sum(flag.int) < round(n.int/2)){flag.int = flag.int.cluster}

		# compute the intronic signal
		out.int = apply(as.matrix(t[i.int,ZTr])[flag.int,],2,median); out.int = out.int/mean(out.int)
		var.int = mean(apply(as.matrix(t[i.int,ZTr])[flag.int,],2,var))
				
		###### Selection of the exonic probes
		
		
		out.ex = apply(as.matrix(t[i.ex,ZTr]),2,median); out.ex = out.ex/mean(out.ex)

		cor.int.ex = cor(out.int, out.ex)
		if(cor.int.ex<0.9){
			COR.INT = rep(0,n.ex); COR.EX = COR.INT
			COR.INT = apply(t[i.ex,ZTr],1,cor,out.int)
			COR.EX = apply(t[i.ex,ZTr],1,cor,out.ex)
			flag.ex.cor = ((COR.INT/COR.EX)<1)
		}else{flag.ex.cor = rep(TRUE, n.ex)}
		
		
		if((n.int>=5)&(sum(flag.int)>=3))
		{
			int = log2(t[i.int,ZT][flag.int,]); 
			ex  = log2(t[i.ex,ZT])
			t.test = t.test(int,ex);
			if( t.test$p.value < 0.1){		
				p = matrix(1,n.ex,12)
				mu = apply(int,2,mean); 
				sigma = apply(int,2,var)
				p = 1-t(apply(t(ex),2,pnorm, mean = mu, sd = sqrt(sigma)))
				Sp = apply(p<0.1,1,sum); 
				flag.ex.2 = (Sp>0);
			}
			else{flag.ex.2 = rep(TRUE, n.ex)}
		}
		else{flag.ex.2 = rep(TRUE, n.ex)}
		flag.ex = flag.ex.cor&flag.ex.2
		if(sum(flag.ex)<min(3,round(n.ex/2)+1)){flag.ex = flag.ex.cor|flag.ex.2 }
		out.ex = apply(as.matrix(t[i.ex,ZTr][flag.ex,]),2,median); out.ex = out.ex/mean(out.ex)
		var.ex = mean(apply(as.matrix(t[i.ex,ZTr][flag.ex,]),2,var))
		
		### PLOTS
		#plot(1,1,xlim = c(0,1),ylim = c(0,1), xlab = 'COR.INT', ylab = 'COR.EX', type = 'n')
		#text(COR.INT,COR.EX, labels = c(1:n.ex)); points(COR.INT,COR.EX,pch = 16)
		#abline(a = 0, b = 1)
		#quartz()
		#matplot(t(t[i.ex,ZTr]), type = 'l', col =rgb(0,0,1,0.3), main = cor.int.ex)
		#points(out.int, type = 'l', lwd = 2, col = 'green4')
		#points(out.ex, type = 'l', lwd = 2, col = 'blue')
		#matplot(t(t[i.ex,ZTr][c(47,20,49),]), type = 'l', col = 'black', add = TRUE)
		#matplot(t(t[i.ex,ZTr][c(36,15,3),]), type = 'l', col = 'red', add = TRUE)
		#quartz()
		#plot(COR.INT/COR.EX); text(COR.INT/COR.EX, labels = c(1:n.ex)); abline(h =1)
	}
	return(list(ex = out.ex, int = out.int, var.ex = var.ex, var.int = var.int, flag.ex = flag.ex, flag.int = flag.int))
	
}

########################## here is a modified version of signal.selection function

select = function(t = t, gene = gene, probe = 'exon', with.rejected = FALSE)
{
	js = (t$gene == gene);
	if(!with.rejected){js = js & (t$kept)}
	jna = is.na(js);
	js[jna] = FALSE
	if(probe == 'exon'){j = js&(t$ex_in == 1)}else if(probe == 'intron'){j = js&(t$ex_in == -1)}else{j = js}
	return(j)
}


signal.selection.myv = function(t = t1, gene_list = gene_list, n.days = 2) ### no modification for this function # not used in main code
{
	
	N = length(gene_list)	
	
	EX_SIGNAL = matrix(NA,N,12*n.days); ### selected exonic singal for each gene 
	INT_SIGNAL = EX_SIGNAL; ## selected intronic signal for each gene 
	colnames(EX_SIGNAL) = paste(colnames(t)[ZT],'.ex',sep=''); 
	colnames(INT_SIGNAL) = paste(colnames(t)[ZT],'.int', sep = ''); 
	
	gene_signal = data.frame(gene = gene_list, EX_SIGNAL, INT_SIGNAL,exon.median = rep(NA,N), exon.var = rep(NA,N), intron.median = rep(NA,N), intron.var = rep(NA,N), stringsAsFactors = FALSE) ## construct the final output
	ZT_ex = c(2:(12*n.days+1)); 
	ZT_int = ZT_ex +12*n.days
	
	flag.ex = rep(FALSE,nrow(t)); 
	flag.int =rep(FALSE, nrow(t));
	
	for(i in 1:N)
	{
		if((i%%140)==0){cat(round(100*i/N, digits = 1),'%\n')} ### indicating the progress
		gene = gene_list[i]; 
		#cat(gene, '\n')
		i.ex = select(t = t,gene = gene,probe = 'exon')
		i.int = select(t = t,gene = gene,probe = 'intron')
		
		out = gene.signal.selection.myv(t = t, i.ex = i.ex, i.int = i.int)
		
		
		gene_signal[i,ZT_ex] = out$ex
		gene_signal[i,ZT_int] = out$int
		gene_signal$exon.var[i] = out$var.ex
		gene_signal$intron.var[i] = out$var.int
		gene_signal$exon.median[i] = median(unlist(t[i.ex,ZT][out$flag.ex,]))
		gene_signal$intron.median[i] = median(unlist(t[i.int,ZT][out$flag.int,]))
		flag.ex[i.ex] = out$flag.ex; flag.int[i.int] = out$flag.int
	}
	return(list(T = gene_signal,flag.ex = flag.ex, flag.int = flag.int))
}



gene.signal.selection.myv = function(t = t1, i.ex = i.ex, i.int = i.int)
{
	
	#t = test.table
	#t = t[t$kept,]
	gene_list = unique(t$gene)
	
	pdf("testing_examples_exon-vs-intron.pdf",width=12,height=10)
	for(i in 1:length(gene_list))
	{
	#
	gene = gene_list[i]; 
	print(i);print(gene);
	#gene = 'Serpina3f'
	cat(gene, '\n')
	i.ex = select(t = t,gene = gene,probe = 'exon')
	i.int = select(t = t,gene = gene,probe = 'intron')
	out.ex = rep(NA,12); 
	out.int = out.ex; 
	var.int = NA; 
	var.ex = NA;
	n.ex = sum(i.ex, na.rm = TRUE); 
	n.int = sum(i.int, na.rm = TRUE);
	#matplot(c(0:23)*2,t(log2(test.table[i.int,c(1:24)])),type='l')
	#matplot(c(0:23)*2,t(log2(test.table[i.ex,c(1:24)])),type='l')
	if((n.ex < 2)||(n.int < 2)){return(list(ex = out.ex,int = out.int,var.ex = var.ex,var.int = var.int, flag.ex = rep(FALSE,n.ex),flag.int = rep(FALSE,n.int)))}else{		
		
		#### Selection of the intronic probes
		
		flag.int.cluster = rep(TRUE, n.int)
		flag.t.s = rep(TRUE, n.int)
		if(n.int > 2)
		{
			#clustering
			k = kmeans(log2(t[i.int,ZT]),2);
			c.high = which.max(apply(k$center,1,mean)); 
			n.high = sum(k$cluster == c.high)
			if(n.high>=round(n.int/2))
			{
				med.high = apply(t[i.int,ZT][k$cluster == c.high,],2,median); 
				med.low = apply(t[i.int,ZT][k$cluster != c.high,],2,median);
				if(cor(med.high, med.low)<0.5){flag.int.cluster[k$cluster != c.high] = FALSE }
			}
		# temporal smoothness
			threshold = 0.3
			scores = apply(t[i.int,ZT],1,score_temporal_smoothness); 
			j = scores < threshold; 
			k = j;
			if(sum(!j)>0){median = apply(t[i.int,ZT],2,median); cor = apply(t[i.int,ZT],1,cor,median); k = cor>0.7 }
			flag.t.s = j|k;
		}
		flag.int = flag.int.cluster & flag.t.s; 
		if(sum(flag.int) < round(n.int/2)){flag.int = flag.int.cluster}
		
		# compute the intronic signal
				
		out.int = apply(as.matrix(t[i.int,ZT])[flag.int,],2,median);
		
		#mm.int = mean(out.int)
		#out.int[1:12] = out.int[1:12]/mean(out.int[1:12])*mm.int
		#out.int[13:24] = out.int[13:24]/mean(out.int[13:24])*mm.int
		
		
		out.int2 = apply(as.matrix(t[i.int,ZT])[flag.int,],2,mean);
		#out.int3 = apply(as.matrix(t[i.int,ZT]),2,mean);
		#out.int = out.int/mean(out.int)
		var.int = mean(apply(as.matrix(t[i.int,ZT])[flag.int,],2,var))
		
		###### Selection of the exonic probes
			
		out.ex = apply(as.matrix(t[i.ex,ZT]),2,median); 
		out.ex = out.ex/mean(out.ex)
		
		cor.int.ex = cor(out.int, out.ex)
		if(cor.int.ex<0.9){
			COR.INT = rep(0,n.ex); 
			COR.EX = COR.INT
			COR.INT = apply(t[i.ex,ZT],1,cor,out.int)
			COR.EX = apply(t[i.ex,ZT],1,cor,out.ex)
			flag.ex.cor = ((COR.INT/COR.EX)<1)
		}else{flag.ex.cor = rep(TRUE, n.ex)}
		
		
		if((n.int>=5)&(sum(flag.int)>=3)){
			int = log2(t[i.int,ZT][flag.int,]); 
			ex  = log2(t[i.ex,ZT])
			t.test = t.test(int,ex);
			if( t.test$p.value < 0.1){		
				p = matrix(1,n.ex,12)
				mu = apply(int,2,mean); sigma = apply(int,2,var)
				p = 1-t(apply(t(ex),2,pnorm, mean = mu, sd = sqrt(sigma)))
				Sp = apply(p<0.1,1,sum); flag.ex.2 = (Sp>0);
			}else{flag.ex.2 = rep(TRUE, n.ex)}
		}else{flag.ex.2 = rep(TRUE, n.ex)}
		flag.ex = flag.ex.cor&flag.ex.2
		if(sum(flag.ex)<min(3,round(n.ex/2)+1)){flag.ex = flag.ex.cor|flag.ex.2 }
		
		out.ex = apply(as.matrix(t[i.ex,ZT][flag.ex,]),2,median);
		#out.ex[1:12] = apply(as.matrix(t[i.ex,ZT][flag.ex,1:12]),2,mean);
		#out.ex[13:24] = apply(as.matrix(t[i.ex,ZT][flag.ex,13:24]),2,mean);
		#mm.ex = mean(out.ex)
		#out.ex[1:12] = out.ex[1:12]/mean(out.ex[1:12])*mm.ex
		#out.ex[13:24] = out.ex[13:24]/mean(out.ex[13:24])*mm.ex
		
		#out.ex = out.ex/mean(out.ex)
		var.ex = mean(apply(as.matrix(t[i.ex,ZT][flag.ex,]),2,var))
		
		par(mfrow = c(3,2))
		ylim = range(log2(t[i.ex,c(1:24)][flag.ex,]))
		ylim2 = range(log2(t[i.int,c(1:24)][flag.int,]))
		matplot(c(0:23)*2,t(log2(t[i.int,c(1:24)])),type='b',lwd=2.0,ylim=ylim2,main='Intron before selection')
		matplot(c(0:23)*2,t(log2(t[i.ex,c(1:24)])),type='b',lwd=2.0,ylim=ylim,main=paste(gene,'Exon before selection'))
								
		
				#readline()
		 
		matplot(c(0:23)*2,t(log2(t[i.int,c(1:24)][flag.int,])),ylim=ylim2,type='b',lwd=2.0,main='Intron after selection')
		abline(v=24,col='darkred',lwd=2.0)
		abline(v=22,col='darkred',lwd=2.0)
		
		matplot(c(0:23)*2,t(log2(t[i.ex,c(1:24)][flag.ex,])), ylim=ylim,type='b',lwd=2.0, main='Exon after selection')
		abline(v=24,col='darkred',lwd=2.0)
		abline(v=22,col='darkred',lwd=2.0)
		
		out.ex2 = apply(as.matrix(t[i.ex,ZT])[flag.ex,],2,mean);
		#out.ex3 = apply(as.matrix(t[i.ex,ZT]),2,mean);
		#out.int = out.int/mean(out.int)
				
		lim = range(log2(out.int), log2(out.int2))
		plot(c(0:23)*2,log2(out.int),type='b',col='blue', ylim=lim,main=paste(gene,'pre-mRNA represented by intron'),lwd=2.0)
		points(c(0:23)*2,log2(out.int2),type='b',col='red')
		#abline(v=24,col='darkred',lwd=2.0)
		#abline(v=22,col='darkred',lwd=2.0)
		
		lim = range(log2(out.ex), log2(out.ex2), log2(out.int), log2(out.int2))
		plot(c(0:23)*2,log2(out.ex),type='b',col='blue',ylim=lim,main=paste(gene,'mRNA represented by exons'),lwd=2.0)
		points(c(0:23)*2,log2(out.ex2),type='b',col='red')
		#points(c(0:23)*2,log2(out.ex3),type='b',col='black')
		abline(v=24,col='darkred',lwd=2.0)
		abline(v=22,col='darkred',lwd=2.0)
		abline(h=mean(log2(out.int)),lwd=2.0, col='gray')
		abline(h=mean(log2(out.int2)),lwd=2.0, col='green')
		
		
		
		#readline()
		
		}
	
	}
	dev.off()
	
	
	return(list(ex = out.ex, int = out.int, var.ex = var.ex, var.int = var.int, flag.ex = flag.ex, flag.int = flag.int))
	
}

signal.selection.myv2 = function(t = t1, gene_list = gene_list, n.days = 2) ### no modification for this function
{
	
	N = length(gene_list)	
	
	EX_SIGNAL = matrix(NA,N,12*n.days); ### selected exonic singal for each gene 
	INT_SIGNAL = EX_SIGNAL; ## selected intronic signal for each gene 
	colnames(EX_SIGNAL) = paste(colnames(t)[ZT],'.ex',sep=''); 
	colnames(INT_SIGNAL) = paste(colnames(t)[ZT],'.int', sep = ''); 
	
	gene_signal = data.frame(gene = gene_list, EX_SIGNAL, INT_SIGNAL,exon.median = rep(NA,N), exon.var = rep(NA,N), intron.median = rep(NA,N), intron.var = rep(NA,N), stringsAsFactors = FALSE) ## construct the final output
	ZT_ex = c(2:(12*n.days+1)); 
	ZT_int = ZT_ex +12*n.days
	
	flag.ex = rep(FALSE,nrow(t)); 
	flag.int =rep(FALSE, nrow(t));
	
	for(i in 1:N)
	{
		if((i%%140)==0){cat(round(100*i/N, digits = 1),'%\n')} ### indicating the progress
		gene = gene_list[i]; 
		cat(i, '\n')
		i.ex = select(t = t,gene = gene,probe = 'exon')
		i.int = select(t = t,gene = gene,probe = 'intron')
		
		out = gene.signal.selection.myv2(t = t, i.ex = i.ex, i.int = i.int)
		
		
		gene_signal[i,ZT_ex] = out$ex
		gene_signal[i,ZT_int] = out$int
		gene_signal$exon.var[i] = out$var.ex
		gene_signal$intron.var[i] = out$var.int
		gene_signal$exon.median[i] = median(unlist(t[i.ex,ZT][out$flag.ex,]))
		gene_signal$intron.median[i] = median(unlist(t[i.int,ZT][out$flag.int,]))
		flag.ex[i.ex] = out$flag.ex; flag.int[i.int] = out$flag.int
	}
	return(list(T = gene_signal,flag.ex = flag.ex, flag.int = flag.int))
}


gene.signal.selection.myv2 = function(t = t1, i.ex = i.ex, i.int = i.int) ### version 2 of the modified function
{
	out.ex = rep(NA,12); 
	out.int = out.ex; 
	var.int = NA; 
	var.ex = NA;
	n.ex = sum(i.ex, na.rm = TRUE); 
	n.int = sum(i.int, na.rm = TRUE);
	if((n.ex < 2)||(n.int < 2)){return(list(ex = out.ex,int = out.int,var.ex = var.ex,var.int = var.int, flag.ex = rep(FALSE,n.ex),flag.int = rep(FALSE,n.int)))}else{		
		
		#### Selection of the intronic probes
		
		flag.int.cluster = rep(TRUE, n.int)
		flag.t.s = rep(TRUE, n.int)
		if(n.int > 2)
		{
			#clustering
			k = kmeans(log2(t[i.int,ZT]),2);
			c.high = which.max(apply(k$center,1,mean)); 
			n.high = sum(k$cluster == c.high)
			if(n.high>=round(n.int/2))
			{
				med.high = apply(t[i.int,ZT][k$cluster == c.high,],2,median); 
				med.low = apply(t[i.int,ZT][k$cluster != c.high,],2,median);
				if(cor(med.high, med.low)<0.5){flag.int.cluster[k$cluster != c.high] = FALSE }
			}
			# temporal smoothness
			threshold = 0.3
			scores = apply(t[i.int,ZT],1,score_temporal_smoothness); 
			j = scores < threshold; 
			k = j;
			if(sum(!j)>0){median = apply(t[i.int,ZT],2,median); cor = apply(t[i.int,ZT],1,cor,median); k = cor>0.7 }
			flag.t.s = j|k;
		}
		flag.int = flag.int.cluster & flag.t.s; 
		if(sum(flag.int) < round(n.int/2)){flag.int = flag.int.cluster}
		
		# compute the intronic signal
		
		#out.int = apply(as.matrix(t[i.int,ZT])[flag.int,],2,mean);
		out.int = apply(as.matrix(t[i.int,ZT])[flag.int,],2,median);
		
		var.int = mean(apply(as.matrix(t[i.int,ZT])[flag.int,],2,var))
		
		###### Selection of the exonic probes
		
		out.ex = apply(as.matrix(t[i.ex,ZT]),2,median); 
		#out.ex = out.ex
		
		cor.int.ex = cor(out.int, out.ex)
		if(cor.int.ex<0.9){
			COR.INT = rep(0,n.ex); 
			COR.EX = COR.INT
			COR.INT = apply(t[i.ex,ZT],1,cor,out.int)
			COR.EX = apply(t[i.ex,ZT],1,cor,out.ex)
			flag.ex.cor = ((COR.INT/COR.EX)<1)
		}else{flag.ex.cor = rep(TRUE, n.ex)}
		
		
		if((n.int>=5)&(sum(flag.int)>=3))
		{
			int = log2(t[i.int,ZT][flag.int,]); 
			ex  = log2(t[i.ex,ZT])
			t.test = t.test(int,ex);
			if( t.test$p.value < 0.1)
			{		
				p = matrix(1,n.ex,12)
				mu = apply(int,2,mean); sigma = apply(int,2,var)
				p = 1-t(apply(t(ex),2,pnorm, mean = mu, sd = sqrt(sigma)))
				Sp = apply(p<0.1,1,sum); flag.ex.2 = (Sp>0);
			}else{flag.ex.2 = rep(TRUE, n.ex);}
		}else{flag.ex.2 = rep(TRUE, n.ex)}
		
		flag.ex = flag.ex.cor&flag.ex.2
		if(sum(flag.ex)<min(3,round(n.ex/2)+1)){flag.ex = flag.ex.cor|flag.ex.2 }
		
		#out.ex = apply(as.matrix(t[i.ex,ZT][flag.ex,]),2,mean);
		out.ex = apply(as.matrix(t[i.ex,ZT][flag.ex,]),2,median);
				
		#out.ex = out.ex/mean(out.ex)
		var.ex = mean(apply(as.matrix(t[i.ex,ZT][flag.ex,]),2,var))
	
	}
	
	return(list(ex = out.ex, int = out.int, var.ex = var.ex, var.int = var.int, flag.ex = flag.ex, flag.int = flag.int))
	
}


################################################### finishing line of modified version of signal.selection function
	
	
	
plot.selection = function(t = t1, T = T1, gene = 'Per2',  probe.type = 'exon', n.days = 2){
	ZTa = ZT
	if(probe.type == 'exon'){ZT = ZT.ex; flag = t$kept}else{ZT = ZT.int;flag = t$kept}
	
	ind = select(t = t, gene = gene, probe = probe.type, with.rejected = TRUE)
	if(sum(ind)>0){tr = t[ind,ZTr]; ta = log2(t[ind, ZTa]); col = rgb(0,0.8,1,0.05) }else{tr = matrix(0,2,12*n.days);ta = matrix(0,2,12*n.days); col = rgb(0,0,0,0)}
	#cat('\t', tr, '\n')

	
	k = which(T$gene == gene)
	if(length(k)>0){Tr = unlist(T[k,ZT])}else{Tr = rep(0,12*n.days)}
	#cat('\t', Tr, '\t')
	#zt = seq(0,22,by = 2)
	
	
	### Absolute levels
	
	title = paste(gene, '-', probe.type); ylim = c(0,16)
	plot(1,1, type = 'n', main = title, xlab = 'Time [in hour of a day]', ylab = 'Absolute signal, log2 scale', xlim = range(zt), ylim = ylim , axes = FALSE)
	axis(2); axis(1, at = zt); 	
	palette(c('steelblue3', 'gray','orange'))
	l = -flag[ind] +2 + ((t$rejection.code == 'exon')|(t$rejection.code == 'intron'))[ind]; 
	if(sum(ind)>0){matplot(zt,t(ta), type = 'l', col = l, lty = l, lwd = 2, add = TRUE)}
	palette('default')

	
	### Relative amplitude
	
	
	title = paste(gene, '-', probe.type); ylim = range(tr,Tr, na.rm = TRUE) #range(Tr, na.rm = TRUE) 
	plot(1,1, type = 'n', main = title, xlab = 'Time [in hour of a day]', ylab = 'Relative signal, lin scale', xlim = range(zt), ylim = ylim , axes = FALSE)
	axis(2); axis(1, at = zt); abline(h = 1, col = 'gray')
	
	#matplot(zt,t(tr), type = 'l', col = col, lwd = 10, add = TRUE, lty = 1)
	
	
	palette(c('steelblue3', 'gray','orange'))
	#if(sum(ind)>0){matplot(zt,t(tr), type = 'l', col = t$cross.hyb.type[ind], lwd = 2, add = TRUE, lty = 1)}
	l = -flag[ind] +2 + ((t$rejection.code == 'exon')|(t$rejection.code == 'intron'))[ind]; 
	if(sum(ind)>0){matplot(zt,t(tr), type = 'l', col = l, lty = l, lwd = 2, add = TRUE)}
	points(zt,Tr, type = 'b', col = 'black', pch = 16, lwd = 2)
	palette('default')
	#cat(probe.type , ' done \n')
	
	
	
	}
	
	
	
remove.gene.witout.signal = function(T = T1){
	j = !is.na(apply(T[,-1],1,sum)); #j = (apply(T[,-1],1,prod)) == 0
	return(j)
	}	
	
	
	

make.beds = function(){
		
	o = order(t$chr[t$kept],t$start[t$kept])

	phase.color.vec = apply(col2rgb(rainbow(240)),2,paste,collapse = ',')
	phase.color = phase.color.vec[ceiling(10*t$phase[t$kept])]
	phase = data.frame(t[t$kept, c('chr','start','stop','probeset.ID')], score = rep(0,sum(t$kept)), t[t$kept,c('strand','start','stop')], phase = phase.color, stringsAsFactors = FALSE)
	phase = phase[o,]
	
	bed.file = paste('ucsc_BED/phase.bed', sep = '')
	track.header = paste("track name='phase' description=' ' priority=",1," visibility=dense itemRgb=On",sep = '')
	write.table(track.header, file = bed.file, quote = FALSE, row.names = FALSE, col.names = FALSE)
	write.table(phase, file = bed.file, quote = FALSE,append = TRUE, sep = '\t', row.names = FALSE, col.names = FALSE)
	
	
	
	Ncol = 100
	colorlow = 'black';  colormed1 = 'blue' ; colormed2 = 'royalblue'; colormed3  = 'steelblue1' ; colormed4 = 'lightsteelblue1'; colorhigh = 'white'
	col.ramp <- colorRampPalette(c(colorlow, colormed1 , colormed2 , colormed3 ,colormed4, colorhigh), space = "Lab")
	pval.color.vec = apply(col2rgb(col.ramp(Ncol)),2,paste,collapse = ',')
	pval.color = pval.color.vec[ceiling(Ncol*t$pval[t$kept])]
	pval = data.frame(t[t$kept, c('chr','start','stop','probeset.ID')], score = rep(0,sum(t$kept)), t[t$kept,c('strand','start','stop')], pval = pval.color, stringsAsFactors = FALSE)
	pval = pval[o,]
	
	bed.file = paste('ucsc_BED/pval.bed', sep = '')
	track.header = paste("track name='pval' description=' ' priority=",length(ZTr)+2," visibility=dense itemRgb=On",sep = '')
	write.table(track.header, file = bed.file, quote = FALSE, row.names = FALSE, col.names = FALSE)
	write.table(pval, file = bed.file, quote = FALSE,append = TRUE, sep = '\t', row.names = FALSE, col.names = FALSE)
	
	
	
	Ncol = 101
	timepoint.index.value.no.cut = as.matrix(round(log10(t[t$kept,ZTr])*(Ncol-1)+ceiling(Ncol/2)))
	timepoint.index.value = pmax(pmin(timepoint.index.value.no.cut,Ncol),1)
	#exon
	colorlow = 'steelblue4';  colormedlow = 'steelblue1' ; colormed = 'white'; colormedhigh  = 'orangered' ; colorhigh = 'orangered3'
	col.ramp <- colorRampPalette(c(colorlow, colormedlow , colormed , colormedhigh , colorhigh), space = "Lab")
	timepoint.color.vec.exon = apply(col2rgb(col.ramp(Ncol)),2,paste,collapse = ',')
	j.ex = which(t$ex_in[t$kept] == 1)
	#inton
	colorlow = 'darkturquoise';  colormedlow = 'turquoise' ; colormed = 'white'; colormedhigh  = 'gold1' ; colorhigh = 'darkgoldenrod2'
	col.ramp <- colorRampPalette(c(colorlow, colormedlow , colormed , colormedhigh , colorhigh), space = "Lab")
	timepoint.color.vec.intron = apply(col2rgb(col.ramp(Ncol)),2,paste,collapse = ',')
	j.int = which(t$ex_in[t$kept] == -1)
	
	for(i in 1:length(ZTr)){
		timepoint.color = vector(mode = 'character',length = sum(t$kept))
		timepoint.color[j.ex] = timepoint.color.vec.exon[timepoint.index.value[j.ex,i]]
		timepoint.color[j.int] = timepoint.color.vec.intron[timepoint.index.value[j.int,i]]
	
		TimePoint = data.frame(t[t$kept, c('chr','start','stop','probeset.ID')], score = rep(0,sum(t$kept)), t[t$kept,c('strand','start','stop')], timepoint = timepoint.color, stringsAsFactors = FALSE)
		TimePoint = TimePoint[o,]
		
		bed.file = paste('ucsc_BED/RelAmpl_ZT',2*(i-1),'.bed', sep = '')
		track.header = paste("track name='RelAmpl_ZT",2*(i-1),"' description=' ' priority=",i+1," visibility=dense itemRgb=On",sep = '')
		write.table(track.header, file = bed.file, quote = FALSE, row.names = FALSE, col.names = FALSE)
		write.table(TimePoint, file = bed.file, quote = FALSE,append = TRUE, sep = '\t', row.names = FALSE, col.names = FALSE)

	}
	
	
	Ncol = 100
	colorlow = 'gray90';  colormed = 'firebrick1' ;  colorhigh = 'black'
	col.ramp <- colorRampPalette(c(colorlow, colormed, colorhigh), space = "Lab")
	mean.color.vec = apply(col2rgb(col.ramp(Ncol)),2,paste,collapse = ',')
	mean.color = mean.color.vec[ceiling(log2(t$mean[t$kept]))/15*(Ncol-1)+1]
	mean = data.frame(t[t$kept, c('chr','start','stop','probeset.ID')], score = rep(0,sum(t$kept)), t[t$kept,c('strand','start','stop')], mean = mean.color, stringsAsFactors = FALSE)
	mean = mean[o,]
	
	bed.file = paste('ucsc_BED/mean.bed', sep = '')
	track.header = paste("track name='mean' description=' ' priority=",length(ZTr)+3," visibility=dense itemRgb=On",sep = '')
	write.table(track.header, file = bed.file, quote = FALSE, row.names = FALSE, col.names = FALSE)
	write.table(mean, file = bed.file, quote = FALSE,append = TRUE, sep = '\t', row.names = FALSE, col.names = FALSE)
	
	
	
	
	
	oo = order(t$chr,t$start)
	green = '0,211,55'; red = '252,86,0'
	kept.color = rep(red,nrow(t)); kept.color[t$kept[oo]] = green; name =  paste(t$probeset.ID[oo],t$rejection.code[oo]); name = gsub(' ','_', name)
	Kept = data.frame(t[oo, c('chr','start','stop')],name  = name , score = rep(0,nrow(t)), t[oo,c('strand','start','stop')], kept = kept.color, stringsAsFactors = FALSE); j = which(Kept$chr == '---'); Kept = Kept[-j,]
	
	bed.file = paste('ucsc_BED/Kept.bed', sep = '')
	track.header = paste("track name='Kept' description=' ' priority=",length(ZTr)+4," visibility=squish itemRgb=On",sep = '')
	write.table(track.header, file = bed.file, quote = FALSE, row.names = FALSE, col.names = FALSE)
	write.table(Kept, file = bed.file, quote = FALSE,append = TRUE, sep = '\t', row.names = FALSE, col.names = FALSE)


}




compute.statistics.on.probeset.numbers = function(){
	gene.list = T$gene;
	j = match(t$gene,T$gene)
	N = nrow(T); numbers = matrix(0,N,7); colnames(numbers) = c('kept.intronic','rejected.intronic', 'intronic.probeset','kept.exonic','rejected.exonic','exonic.probeset', 'total.probeset')
	for(i in 1:N){
		if((i%%300) == 0){cat(100*i/N,' % \n')}
		numbers[i,'kept.intronic'] = sum((j == i) & t$kept & (t$ex_in == -1), na.rm = TRUE)
		numbers[i,'rejected.intronic'] = sum((j == i) & !t$kept & (t$ex_in == -1), na.rm = TRUE)
		numbers[i,'intronic.probeset'] = sum((j == i)  & (t$ex_in == -1), na.rm = TRUE)
		numbers[i,'kept.exonic'] = sum((j == i) & t$kept & (t$ex_in == 1), na.rm = TRUE)
		numbers[i,'rejected.exonic'] = sum((j == i) & !t$kept & (t$ex_in == 1), na.rm = TRUE)
		numbers[i,'exonic.probeset'] = sum((j == i)  & (t$ex_in == 1), na.rm = TRUE)
		numbers[i,'total.probeset'] = sum((j == i), na.rm = TRUE)
		}
	return(numbers)
}
	
	
	
hist.fourier.compo = function(replicate = 1, type = 'exon', pval = 1,exp.lim = 5){
	if(replicate == 1){T = T1}else{T = T2}
	if(type == 'exon'){ZT = ZT.ex; j = (T$pval.ex<=pval)&(T$exon.median>exp.lim)}else{ZT = ZT.int; j = (T$pval.int<=pval)&(T$intron.median>exp.lim)}
	
	FFT = t(apply(T[j,ZT],1,fft))/length(ZT)	
	ABS = abs(FFT)
	
#	breaks = seq(0,1,by = 0.002)
#	HH  = apply(ABS, 2,hist, breaks = breaks, plot = FALSE)
#	mids = HH[[1]]$mids
#	H = matrix(0,nrow = 5, ncol = length(mids))
#	
#	for(i in 1:5){H[i,] = HH[[i+1]]$counts}
#
#	plot(1,1,type = 'n', ylim =range(H), xlim = c(0,0.5), xlab = 'Fourier component (Modulo)', ylab = 'number of genes', main = paste('Replicate ', replicate, '-',type, ' | pval.max = ',pval))
#	palette(rainbow(n = 5,alpha = 1))
#	matplot(mids, t(H), col= c(1:5), type = 'l', lwd = 2, lty = 1, add = TRUE)
#	legend('topright',legend = 24/c(1:5), lty = 1, col = c(1:5), lwd = 2)
#	
	
	ymax = 0;for(i in 1:5){d = density(ABS[,i+1]);  ymax = max(d$y,ymax)}
	plot(1,1,type= 'n', ylim = c(0,ymax),xlim = c(0,0.5), xlab = 'Fourier component (Modulo)', ylab = 'density', main = paste('Replicate ', replicate, '-',type, ' | pval.max = ',pval))
	for(i in 1:5){points(density(ABS[,i+1]), col = i, type = 'l', lwd = 3);}
	legend('topright',legend = 24/c(1:5), lty = 1, col = c(1:5), lwd = 3)
	
}
	
	
hist.fourier.compo.2 = function(replicate = 1, pval = 1,exp.lim = 5){
	if(replicate == 1){T = T1}else{T = T2}
	j.ex = (T$pval.ex<=pval)&(T$exon.median>exp.lim); j.int = (T$pval.int<=pval)&(T$intron.median>exp.lim)
	
	
	FFT.ex = t(apply(T[j.ex,ZT.ex],1,fft))/length(ZT); FFT.int = t(apply(T[j.int,ZT.int],1,fft))/length(ZT)
	ABS.ex = abs(FFT.ex); ABS.int = abs(FFT.int)
	
	col.ex = rainbow(n = 5, s = 1, v = 1)
	col.int = rainbow(n = 5, s = 1, v = 0.6)

	for(i in 1:5){
		d.ex = density(ABS.ex[,i+1])
		d.int = density(ABS.int[,i+1])

		plot(1,1,type= 'n', ylim = range(d.ex$y,d.int$y),xlim = c(0,0.5), xlab = 'Fourier component (Modulo)', ylab = 'density', main = paste('Replicate ', replicate, '- component at',24/i,'h', '- pval : ',pval));
		points(d.ex, col = col.ex[i], lty = 1, lwd = 2, type = 'l') 
		points(d.int, col = col.int[i], lty = 1, lwd = 2, type = 'l')
		legend('topright',legend = c('exon', 'intron'), lty = 1, col = c(col.ex[i], col.int[i]), lwd = 2)
	}
		
}	
	
scatter.fourier.compo = function(replicate = 1, pval = 1,exp.lim = 5){
	if(replicate == 1){T = T1}else{T = T2}
	j.ex = (T$pval.ex<=pval)&(T$exon.median>exp.lim); j.int = (T$pval.int<=pval)&(T$intron.median>exp.lim)
	
	
	FFT.ex = t(apply(T[j.ex,ZT.ex],1,fft))/length(ZT.ex); FFT.int = t(apply(T[j.int,ZT.int],1,fft))/length(ZT.int)
	ABS.ex = abs(FFT.ex); ABS.int = abs(FFT.int)
	
	col.ex = rgb(1,0.2,0,0.1)
	col.int = rgb(0.5,0.5,0,0.1)
	
	smoothScatter(ABS.ex[,2:3], xlim = range(ABS.ex[,2:3]), ylim = range(ABS.ex[,2:3]), main = paste('Replicate ',replicate,' - exon'), xlab = '24h component', ylab = '12h component');	abline(a = 0, b = 1)
	smoothScatter(ABS.int[,2:3], xlim = range(ABS.int[,2:3]), ylim = range(ABS.int[,2:3]), main = paste('Replicate ',replicate,' - intron'), xlab = '24h component', ylab = '12h component');	abline(a = 0, b = 1)
	
#	colors  <- densCols(ABS.ex[,2:3])
#	plot(ABS.ex[,2:3], col=colors, pch=16, cex = 0.7, xlim = range(ABS.ex[,2:3]), ylim = range(ABS.ex[,2:3]));	abline(a = 0, b = 1)
#	colors  <- densCols(ABS.int[,2:3])
#	plot(ABS.int[,2:3], col=colors, pch=16, cex = 0.7, xlim = range(ABS.int[,2:3]), ylim = range(ABS.int[,2:3]));	abline(a = 0, b = 1)

}	


decreasing.errors = function(replicate = 1 ,type = 'exon', pval = 1, exp.lim = 0){
	if(replicate == 1){T = T1}else{T = T2}
	if(type == 'exon'){j = (T$pval.ex<=pval)&(T$exon.median>exp.lim); ZT = ZT.ex}else{j = (T$pval.int<=pval)&(T$intron.median>exp.lim); ZT = ZT.int}
	
	T = T[j,]; N = sum(j); K = 5
	Errors = matrix(Inf,N,K)
	for(k in 1:K){
		cat(k,'\n')
		for(i in 1:N){
			if((i/N)%%0.25 == 0){cat('\t',100*round(i/N, digits = 4),'\n')}
			X = unlist(T[i,ZT])
			x = estimate.x(X,k = k-1)
			error = sum((X-x)^2)
			Errors[i,k] = error
		}
	}
	return(Errors)
}


estimate.x = function(X = rep(1,12), k = 0){
	sup = seq(0,2*pi,by = pi/6); sup = sup[1:12]
	fft = fft(X)/length(X); 
	x = mean(X) 
	if(k>0){x = x + 2*abs(fft[2])/mean(X) * cos(sup-phase.n(X,n=1))}
	if(k>1){x = x + 2*abs(fft[3])/mean(X) * cos(2*sup - phase.n(X,n = 2))}
	if(k>2){x = x + 2*abs(fft[4])/mean(X) * cos(3*sup - phase.n(X,n = 3))} 
	if(k>3){x = x + 2*abs(fft[5])/mean(X) * cos(4*sup - phase.n(X,n = 4))}
	return(x)
}	


select.rhythmic.genes = function( probe.type = 'exon', period = 24, fdr = 0.1, ampl.lim = 0.15, fact = 1)
{
	if(probe.type == 'exon')
	{
		if(period == 24){pval = T$pval.ex; rel.ampl = T$rel.ampl.ex}else{pval = T$pval.12.ex; rel.ampl = T$rel.ampl.12.ex}; var = T$exon.var
	}else{
		if(period == 24){pval = T$pval.int; rel.ampl = T$rel.ampl.int}else{pval = T$pval.12.int; rel.ampl = T$rel.ampl.12.int}; var = T$intron.var
	}
	pval.lim = fdr2pval(P = pval, fdr = fdr); if(pval.lim == 1){pval.lim = 0}
	j = (pval <= pval.lim)&(rel.ampl > ampl.lim)&(rel.ampl > fact*sqrt(var))
	return(j)
}






check.for.other.pattern = function(probe.type = 'exon'){
	N = nrow(T)
	if(probe.type == 'exon'){ZT = grep('.rel.ampl.ex', colnames(T)); var = T$exon.var; r24 = T$rhythmic.ex.24; median = log2(T$exon.median)}else{ZT = grep('.rel.ampl.int', colnames(T)); var = T$intron.var; r24 = T$rhythmic.int.24;  median = log2(T$intron.median)}
	
	maximum = apply(T[,ZT],1,max);minimum = apply(T[,ZT],1,min); FC = maximum - minimum
	
	with.pattern = (FC > 2*ampl.lim)&((FC/2) > fact*sqrt(var)) #
	n.no.pattern = nrow(T) - sum(with.pattern)
	other.pattern = (with.pattern & !r24)
	n.other.pattern = sum(other.pattern)
	n.24 = sum(r24)
	
	if((n.no.pattern + n.other.pattern + n.24) != nrow(T)){cat('OUPS\n')}
	
	stat = data.frame(n.no.pattern = n.no.pattern,n.24 = n.24, n.other.pattern= n.other.pattern)
	indexes = data.frame(no.pattern = rep(TRUE,nrow(T))&!with.pattern, r24 = r24, other.pattern= other.pattern)
	col = c('gray','tomato','slategray1'); alpha = 0.5; transp.col = c(rgb(col2rgb('gray')[1]/255,col2rgb('gray')[2]/255,col2rgb('gray')[3]/255,alpha), rgb(col2rgb('tomato')[1]/255,col2rgb('tomato')[2]/255,col2rgb('tomato')[3]/255, alpha),rgb(col2rgb('slategray1')[1]/255,col2rgb('slategray1')[2]/255,col2rgb('slategray1')[3]/255, alpha))
	col.border = c('gray30','tomato4','slategray4') 
	labels = c('No Pattern','Circadian [24h]','Other patterns')
	
	n.cat = length(stat); n.comp = length(ZT)/2-1
	H = matrix(NA,nrow = n.cat, ncol = n.comp); rownames(H) = labels; colnames(H) =  paste(round(rep(48,n.comp)/c(2:(n.comp+1)),digits = 1),'h', sep = '')# histogram of components
	COMP = matrix(NA, nrow = N, ncol = n.cat* n.comp)
	F.COMPO = t(apply(T[,ZT],1,fft)); F.COMPO = F.COMPO[,3:(n.comp+2)]; F.COMPO = Mod(F.COMPO); 
	
	for(i in 1:n.cat){
		f.compo = t(apply(T[indexes[,i],ZT],1,fft)); f.compo = f.compo[,3:(n.comp+2)];f.compo = Mod(f.compo); 
		stronger.compo = apply(f.compo, 1, which.max)
		h = hist(stronger.compo, breaks = seq(0.5,n.comp+0.5,by = 1), plot = FALSE)
		H[i,] = h$counts
		COMP[indexes[,i],c(1:n.comp)+n.comp*(i-1)] = f.compo #(c(0:5)*4 + i)
	}
	
	#par(mfrow = c(1,4))
	layout(matrix(c(1,2,6,7,3,3,6,7,4,4,6,7,5,5,6,7),4,4), heights = c(0.5,0.5,1,1))
	
	par(mar = c(0,0,2,0))
	pie(unlist(stat), init.angle = 90, clockwise = TRUE, col = col , border= NA, labels = labels, col.lab = 'deeppink', main = paste(probe.type)) #labels = labels,

	par(mar = c(5, 4, 4, 2) + 0.1)
	boxplot( median[indexes[,1]],median[indexes[,3]],median[indexes[,2]], pch = 16, cex = 0.5,col = col[c(1,3,2)],border = col.border[c(1,3,2)], names = labels[c(1,3,2)], ylab = 'median level of expression (log2)')
	
	barplot(H, col = col, border = NA, names = colnames(H), ylab = 'number of genes' ,main = paste(probe.type, '\n histogram of the strongest fourier component'));
	legend('topright', legend = labels, fill = col, border = NA, lty = 0, bty = 'n')
	
	ind = c(1,3)
	barplot(H[ind,]/c(n.no.pattern, n.other.pattern), beside = TRUE, col = col[ind], border = NA, names = colnames(H), ylab = 'fraction of genes', main = paste(probe.type, '\n histogram of the strongest fourier component'))
	legend('topright', legend = labels[ind], fill = col[ind], border = NA, lty = 0, bty = 'n')


	breaks = 2^seq(-3.2,3.2,by = 0.1)#seq(0,6.5,by = 0.1)
	h1 = hist(FC[indexes[,1]], breaks = breaks, plot = FALSE)
	h2 = hist(FC[indexes[,2]], breaks = breaks, plot = FALSE)
	h3 = hist(FC[indexes[,3]], breaks = breaks, plot = FALSE)


	plot(1,1,type = 'n', xlim = range(h1$mids), axes = FALSE, log = 'x', ylim = range(h1$count/sum(h1$count),h2$count/sum(h2$count),h3$count/sum(h3$count)), xlab = 'Fold Change', ylab = '% of genes')
	polygon(h1$mids, h1$count/sum(h1$count), col = transp.col[1], border = col.border[1], lwd = 1)
	polygon(h2$mids, h2$count/sum(h2$count), col = transp.col[2], border = col.border[2], lwd = 1)
	polygon(h3$mids, h3$count/sum(h3$count), col = transp.col[3], border = col.border[3], lwd = 1)
	abline(v = 2*ampl.lim-0.01, col = 'green3', lty = 2, lwd = 2)
	axis(1); axis(2, las =1, at = c(0:100)/100, labels = c(0:100))
	legend('topright', legend = labels, fill = col, border = NA, lty = 0, bty = 'n')



	#par(mfrow = c(1,1))

	boxplot(COMP, col = rep(col, each = n.comp), cex = 0.3, pch = 16, axes = FALSE, at = c(1:(n.cat*(n.comp+1)))[-(n.comp+1)*c(1:n.cat)], xlim = c(0.5,n.cat*(n.comp+1)-0.5), border = rep(col.border, each = n.comp), main = paste(probe.type, '\n Fourier Components by Pattern Categories'))
	axis(2, las = 1); box()
	axis(1, at = c(1:(n.cat*(n.comp+1)))[-(n.comp+1)*c(1:n.cat)], labels = rep(colnames(H), n.cat))
	legend('topright', legend = labels, fill = col, border = NA, lty = 0, bty = 'n')


	boxplot(COMP[,rep(c(0:(n.cat-1))*n.comp+1,n.comp) + rep(c(0:(n.comp-1)), each = n.cat)], col = col, cex = 0.3, pch = 16, axes = FALSE, at =  c(1:(n.comp*(n.cat+1)))[-(n.cat +1)*c(1:n.comp)] , xlim = c(0.5,n.comp*(n.cat+1)-0.5), border = col.border, main = paste(probe.type, '\n Fourier Components by Pattern Categories'))
	axis(2, las = 1); box()
	axis(1, at = c(1:(n.comp*(n.cat+1)))[-(n.cat +1)*c(1:n.comp)], labels = rep(colnames(H), each = n.cat))
	legend('topright', legend = labels, fill = col, border = NA, lty = 0, bty = 'n')
}		






select.fit.candidates = function(T = T1, pval.lim = pval.lim, fdr = 0.1, ampl.lim = ampl.lim, fact = 1){
	
	if(fdr>0){pval.lim = fdr2pval(T$pval.ex, fdr = fdr);cat('exon\n')}
	j1 = (T$pval.ex <= pval.lim)&(T$rel.ampl.ex > ampl.lim)&(T$rel.ampl.ex > fact*sqrt(T$exon.var))
	cat('\t',pval.lim,'\t', sum(j1),'\n')
	
	if(fdr>0){pval.lim = fdr2pval(T$pval.int, fdr = fdr);cat('intron\n')}
	j2 = (T$pval.int <= pval.lim)&(T$rel.ampl.int > ampl.lim)&(T$rel.ampl.int > fact*sqrt(T$intron.var))
	cat('\t',pval.lim,'\t', sum(j2),'\n')
	
	j = j1|j2
	return(j)
	}
	


######## FUNCTIONS FOR SIMULATED DATA

set.global.variable = function()
{
	splicing.rate <<- log(2)/8*60
}

set.global.sigma = function() ## this function is to used to generate fake data with different standard deviation
{
	sigma.s <<- 0.2733805
	sigma.m <<- 0.1962412
}



generate.fake.data = function(T = T[1:100,],X = 100, sd = 0.05, model = 1, parametrization = c('cosine','sigmoid'), sum.species = TRUE, random.splicing=FALSE)
{
	parametrization = parametrization[1]
	set.seed(42)
	
	set.global.variable(); ### to set parameter of splicing rate 'k' which is 8 minutes 
	k = splicing.rate;
	if(random.splicing) {kk=sample(seq(4.1,8.3,by=0.042), X, replace=TRUE)}
	
	
	T$gene = paste('fake_m',model,'_',1:X,sep = '')
	T$intron.median = 2^rnorm(X,mean = 5.7, sd = 1.4) ### make the median of intronic signals
	T$exon.var = 0; 
	T$intron.var = 0;
	bounds = set.bounds(model = 4, parametrization = parametrization); 
	lower = bounds$lower; 
	upper = bounds$upper; 
	gamma.min = lower[1]; 
	gamma.max = upper[1];
	eps.gamma.min = 0.1; 
	eps.gamma.max = upper[2]; 
	if(parametrization == 'sigmoid'){eps.gamma.min = 1.2};
	phase.gamma.min = lower[3]; 
	phase.gamma.max = upper[3];
	
	if(parametrization == 'cosine')
	{
		rel.ampl.int.min = 0.1; rel.ampl.int.max = upper[4];
		phase.int.min = lower[5]; phase.int.max = upper[5];
		rel.ampl.12.int.min = lower[6]; rel.ampl.12.int.max = upper[6];
		phase.12.int.min = lower[7]; phase.12.int.max = upper[7];
	}else{
		fold.change.min = 1.25; 
		fold.change.max = upper[4];
		phase.int.min = lower[5]; 
		phase.int.max = upper[5];
		up.time.min = lower[6]; 
		up.time.max = upper[6];
		down.time.min = lower[7]; 
		up.time.max = upper[7];
	}
	
	if((model == 2)|(model == 4))
	{
		if(parametrization == 'cosine')
		{
			T$rel.ampl.int = sample(seq(rel.ampl.int.min, rel.ampl.int.max,by = 0.01),X, replace = TRUE)
			T$rel.ampl.12.int = apply(cbind(1-T$rel.ampl.int, sample(seq(rel.ampl.12.int.min,rel.ampl.12.int.max,by = 0.01),X, replace = TRUE)),1,min) 
			T$phase.int = sample(seq(phase.int.min ,phase.int.max,by = 0.1),X, replace = TRUE)
			T$phase.12.int = sample(seq(phase.12.int.min,phase.12.int.max,by = 0.1),X, replace = TRUE)
		}else{
					
			T$rel.ampl.int = sample(seq((fold.change.min-1)/(fold.change.min+1), (fold.change.max-1)/(fold.change.max+1),by = 0.01),X, replace = TRUE)
			T$fold.change.int  = (1+T$rel.ampl.int)/(1-T$rel.ampl.int)
			T$phase.int =  sample(seq(phase.int.min ,phase.int.max,by = 0.1),X, replace = TRUE)
			T$up.time.int = pmax(pmin(rnorm(n = X, mean = 12, sd = 4.5), up.time.max), up.time.min)
			T$down.time.int = pmax(apply(cbind(rnorm(n = X, mean = 12, sd = 4.5),24-T$up.time),1,min),down.time.min)
			T$rel.ampl.12.int = rep(0,X)
			T$phase.12.int = rep(0,X)
		}
	}else{
		T$rel.ampl.int = rep(0,X);
		T$rel.ampl.12.int = rep(0,X);
		T$phase.int = rep(0,X); 
		T$phase.12.int = rep(0,X)
		if(parametrization == 'sigmoid'){T$fold.change.int = rep(1,X); T$phase.int = 0; T$up.time.int = rep(12,X); T$down.time.int = rep(12,X)}
	}
		
	if(model>2)
	{
		T$gamma = pmin(pmax(log(2)/exp(rnorm(X, mean = log(3), sd = 1)), gamma.min), gamma.max)
		T$eps.gamma = sample(seq(eps.gamma.min, eps.gamma.max,by = 0.01),X, replace = TRUE)
		if(parametrization == 'sigmoid'){T$eps.gamma = sample(seq((eps.gamma.min-1)/(eps.gamma.min+1), (eps.gamma.max-1)/(eps.gamma.max+1),by = 0.01),X, replace = TRUE); }
		T$phase.gamma = sample(seq(phase.gamma.min, phase.gamma.max,by = 0.1),X, replace = TRUE)
	}else if(model ==2){
		T$gamma = pmin(pmax(log(2)/exp(rnorm(X, mean = log(3), sd = 1)), gamma.min), gamma.max)
		T$eps.gamma = rep(0,X)
		T$phase.gamma = rep(0,X)
	}else{
		T$gamma = rep(1,X)
		T$eps.gamma = rep(0,X)
		T$phase.gamma = rep(0,X)
	}	

	if(parametrization == 'sigmoid'){T$eps.gamma = (1+T$eps.gamma)/(1-T$eps.gamma)}

	if(parametrization == 'cosine')
	{
		for(i in 1:X){
			T[i, ZT.int] = compute.s(t = zt, eps.24.S = T$rel.ampl.int[i] , phase.24.S = T$phase.int[i], eps.12.S = T$rel.ampl.12.int[i], phase.12.S = T$phase.12.int[i])
			if(model == 1){ T[i, ZT.ex] = rep(1,length(zt))}
			if(model == 2){ T[i, ZT.ex] = compute.m(t = zt, gamma = T$gamma[i] ,eps.24.S = T$rel.ampl.int[i] , phase.24.S = T$phase.int[i], eps.12.S = T$rel.ampl.12.int[i], phase.12.S = T$phase.12.int[i])}
			if(model>2){T[i, ZT.ex] = compute.m.changing(t = zt, gamma = T$gamma[i] , eps.gamma = T$eps.gamma[i], phase.gamma = T$phase.gamma[i],eps.24.S = T$rel.ampl.int[i] , phase.24.S = T$phase.int[i], eps.12.S = T$rel.ampl.12.int[i], phase.12.S = T$phase.12.int[i])}
		}
		
		fact.ex = k * T$intron.median/ T$gamma
		if(model>2){ RA.ex = do.rhythmic.analysis(T = T, ZT = ZT.ex); fact.ex  = fact.ex /(1+0.5*T$eps.gamma*RA.ex$rel.ampl * cos(RA.ex$phase-T$phase.gamma))}
	
		abs.int = T[,ZT.int] * T$intron.median
		abs.ex = T[,ZT.ex]* fact.ex
	
	}else{
		abs.ex = T[,ZT.ex]*0
		for(i in 1:X)
		{
			#cat(i,'\n')
			## generate relative signals for introns with average of 1.0; the median of intronic signals is generated at the beginning.
			T[i,ZT.int] = compute.sigmoid(t = zt, fold.change = T$fold.change.int[i], phase = T$phase.int[i], up.time = T$up.time.int[i], down.time = T$down.time.int[i])
			
			## generate absolute singals of exons
			### very weired why multipled again by log(2)
			#if(model == 1){abs.ex[i,] = k * T$intron.median[i] * rep(1,length(zt)) / log(2)}
			if(model == 1)
			{
				abs.ex[i,] = k * T$intron.median[i] * rep(1,length(zt));
				if(random.splicing) abs.ex[i,] = kk[i] * T$intron.median[i] * rep(1,length(zt));
			} 
			if(model>1)
			{
				#cat('gamma =', T$gamma[i],'eps.gamma = ',(1+T$eps.gamma[i])/(1-T$eps.gamma[i]),'phase.gamma = ',T$phase.gamma[i],' fold.change =', T$fold.change[i], 'phase = ',T$phase.int[i], 'up.time = ',T$up.time[i],' down.time = ',T$down.time[i],'\n')
				if(random.splicing){
					abs.ex[i,] = compute.m.sigmoid(t = zt, gamma = T$gamma[i],eps.gamma = T$eps.gamma[i],phase.gamma = T$phase.gamma[i], 
												   fold.change = T$fold.change.int[i], phase = T$phase.int[i], up.time = T$up.time.int[i],
												   down.time = T$down.time.int[i], rescale = FALSE, synthesis.factor = kk[i] * T$intron.median[i])
				}else{
					abs.ex[i,] = compute.m.sigmoid(t = zt, gamma = T$gamma[i],eps.gamma = T$eps.gamma[i],phase.gamma = T$phase.gamma[i], 
												   fold.change = T$fold.change.int[i], phase = T$phase.int[i], up.time = T$up.time.int[i],
												   down.time = T$down.time.int[i], rescale = FALSE, synthesis.factor = k * T$intron.median[i])
				}
								
				#cat('finished!\n')
			}
		}	
		abs.int = T[,ZT.int] * T$intron.median ## absolute singals for introns
		
	}
	
	if(sum.species){abs.ex = abs.ex+abs.int}   ## here to sum up the exonic signals from premrnas and mrnas
	
	log.abs.int = log(abs.int)	
	log.abs.ex = log(abs.ex)
	
	set.global.sigma();
	sd.noise = c(0.1,0.25,0.5);
	if(sd==0.1) {sigma.s = sigma.s/2; sigma.m = sigma.m/2;}
	if(sd==0.5) {sigma.s = sigma.s*2; sigma.m = sigma.m*2;}
	
	noise.int = matrix(rnorm(X*length(ZT.int), mean = 0, sd = sigma.s),nrow = X, ncol = length(ZT.int))
	noise.ex = matrix(rnorm(X*length(ZT.ex), mean = 0, sd = sigma.m),nrow = X, ncol = length(ZT.ex))

	log.noise.abs.int = log.abs.int + noise.int
	log.noise.abs.ex = log.abs.ex + noise.ex
	
	noise.abs.int = exp(log.noise.abs.int)
	noise.abs.ex = exp(log.noise.abs.ex)

	T[,ZT.int] = noise.abs.int/apply(noise.abs.int, 1, mean)
	T[,ZT.ex] =  noise.abs.ex/apply(noise.abs.ex, 1, mean)
	T$exon.median = apply(noise.abs.ex,1,median)		
	return(T)

}


compare.parameters = function(T = T, fitting.method = c('Optim','MCMC'), quality = FALSE, parametrization = c('sigmoid','cosine'))
{
	#fitting.method = c('Optim','MCMC'); quality = FALSE; parametrization = c('sigmoid','cosine');
	parametrization = parametrization[1]
	fitting.method = fitting.method[1]
	if(fitting.method == 'Optim'){best.model = T$LRT.best.model;true.model=T$true.model;}else{best.model = T$Model.MCMC;true.model=T$true.model;}
	
	if(parametrization == 'cosine') # not used
	{
		mains = c('Degr. Rate','Degr. Amplitude','Degr. Phase','Rel. Ampl. - Intron','Phase 24h - Intron','Rel. Ampl. 12h - Intron','Phase 12h - Intron')
		parameter.list = c('gamma','eps.gamma','phase.gamma','rel.ampl.int','phase.int','rel.ampl.12.int','phase.12.int')
	}else{
		mains = c('Degr. Rate','Degr. FC','Degr. Phase','FC - Intron','Phase - Intron','Up-time - Intron','Down-time - Intron')
		parameter.list = c('gamma','eps.gamma','phase.gamma','fold.change.int','phase.int','up.time.int','down.time.int')
	}
	models.p = cbind(c(0,2, 3, 4),c(0,0, 3, 4),c(0,0, 3, 4),c(0, 2,0, 4),c(0, 2,0, 4),c(0, 2,0, 4),c(0, 2,0, 4))
	bounds = set.bounds(model = 4, parametrization = parametrization); 
	lower = bounds$lower; 
	upper = bounds$upper;
	# 7 parpameters and 7 plots
	for(p in 1:7)
	{
		log = ''; 
		if((p==1)|(p==2)|(p==4)){log ='xy'}
		plot(c(lower[p], upper[p]),c(lower[p],upper[p]),type = 'n', xlab = 'True values', ylab = 'Estimated Values', main = mains[p], log = log)
		
		#for(m in c(1:4)[models.p[,p]]){
		for(m in c(2:4)) # for each model 2,3,4
		{	
			eval(parse(text = paste('ttrue = T$', parameter.list[p],'[best.model == m & true.model==m]', sep ='')))
			eval(parse(text = paste('ftrue = T$', parameter.list[p],'[best.model == m & true.model!=m]', sep ='')))
			
			if(fitting.method == 'Optim')
			{
				if(models.p[m,p] == 0){
					estimated = rep(0, sum(best.model == m));
					#print('here')
				}else{
					eval(parse(text = paste('estimated.t = T$', parameter.list[p],'.m',m,'[best.model == m & true.model==m ]', sep =''))); 
					eval(parse(text = paste('estimated.f = T$', parameter.list[p],'.m',m,'[best.model == m & true.model!=m ]', sep ='')));
					
				}
			}else{
				#eval(parse(text = paste('estimated.t = T$', parameter.list[p],'.MCMC',m,'[best.model == m & true.model==m ]', sep =''))); 
				#eval(parse(text = paste('estimated.f = T$', parameter.list[p],'.MCMC',m,'[best.model == m & true.model!=m ]', sep ='')));
				eval(parse(text = paste('estimated.t = T$', parameter.list[p],'.MCMC[best.model == m & true.model==m]', sep ='')))
				eval(parse(text = paste('estimated.f = T$', parameter.list[p],'.MCMC[best.model == m & true.model!=m]', sep ='')))
				
			}
			
			if(quality){if(fitting.method == 'Optim'){cex = T$LRT.quality[best.model == m]/10}else{cex = T$MCMC.quality[best.model == m]}}else{cex = 1}
			if(models.p[m,p]>0)
			{
				points(ttrue, estimated.t, pch = 21, col = m, bg = m, cex = cex);
				points(ftrue, estimated.f, pch = 2, col = m, bg = m, cex = 1.5);
				
			}
		}	
		
		abline(a =0, b = 1, col = 'gray')
	}
	
#	plot(c(0,1),c(0,1),type = 'n', xlab = 'True values', ylab = 'Estimated Values', main = 'Relative Amplitude - Intron')
#	points(T$rel.ampl.int[T$LRT.best.model == 2], T$rel.ampl.int.m2[T$LRT.best.model == 2], pch = 21, col = rep(1:4, each = X)[T$LRT.best.model == 2], bg = 2)
#	points(T$rel.ampl.int[T$LRT.best.model == 4], T$rel.ampl.int.m2[T$LRT.best.model == 4], pch = 21, col = rep(1:4, each = X)[T$LRT.best.model == 4], bg = 4)
#	abline(a =0, b = 1, col = 'gray')
	
}




######################
# Functions for the fits
######################

#### not used now
make.fit.with.all.models.for.all.genes = function(T = T, debug = FALSE, parametrization = c('cosine','sigmoid'), method = c('integration','simulation'), sum.species = FALSE, zt = seq(0,46,by = 2), i.ex = ZT.ex, i.int = ZT.int, absolute.signal = FALSE){
	parametrization = parametrization[1]; method = method[1]
	N = nrow(T); steps = max(10,round(N/100))
	for(i in 1:N){
		if(i%%steps == 0){cat(round(100*i/N),'%   - ',as.character(Sys.time()),'\n')}
		if(debug){cat('starting gene ',T$gene[i],'\n')}
		Param.fit.for.gene = make.fits.with.all.models.for.one.gene(T = T, gene.index = i, debug = debug, parametrization = parametrization, method = method, sum.species = sum.species, zt = zt, i.ex = i.ex, i.int = i.int, absolute.signal = absolute.signal)
		if(i == 1){PARAM.FIT = Param.fit.for.gene}else{PARAM.FIT = rbind(PARAM.FIT, Param.fit.for.gene)}
	}
	return(PARAM.FIT)
}

### Main Function for model selection and fitting

make.fits.with.all.models.for.one.gene = function(T = T, gene.index = 1, debug = FALSE, parametrization = c('cosine','sigmoid'), method =  c('integration','simulation'), sum.species = FALSE, zt = seq(0,46,by = 2), i.ex = ZT.ex, i.int = ZT.int, absolute.signal = FALSE)
{
	#T = T; gene.index = j; debug = TRUE; parametrization = 'sigmoid'; method = 'integration'; sum.species = TRUE; zt = zt; i.ex = ZT.ex; i.int = ZT.int; absolute.signal = TRUE;
	parametrization = parametrization[1]; 
	method = method[1];
	for(model in 1:4)
	{
		if(debug){cat('\t starting model ',model,'\n');}
		param.fit = make.fit.spec.model(T = T, gene.index = gene.index, model = model, debug = debug, parametrization = parametrization, method = method, sum.species = sum.species, zt = zt, i.ex = i.ex, i.int = i.int, absolute.signal = absolute.signal);
		
		if(model == 1)
		{
			Param.fit.for.gene = param.fit
		}else{
			Param.fit.for.gene = c(Param.fit.for.gene, param.fit)
		}
		if(debug){cat('\t model ',model,' finished \n')};
	}
	return(Param.fit.for.gene)
}


make.fit.spec.model = function(T = T, gene.index = 1, model = 1, debug = FALSE, parametrization = c('cosine','sigmoid'),method =  c('integration','simulation'), sum.species = FALSE, zt = seq(0,46,by = 2), i.ex = ZT.ex, i.int = ZT.int, absolute.signal = FALSE)
{
	parametrization = parametrization[1];
	method = method[1];
	param.fit = NA
	S = unlist(T[gene.index, i.int])
	M = unlist(T[gene.index, i.ex])
	if(model == 1)
	{
		s = rep(1,length(S)); 
		m =  rep(1,length(M)); 
		if(absolute.signal)
		{
			set.global.variable(); 
			gamma.m1 = splicing.rate * mean(S)/mean(M);
			names(gamma.m1) = 'gamma.m1';
			m = m*mean(M);
		}
		err = error(S,s,M,m); 
		param.fit = err; 
		names(param.fit) = 'error.m1';
		if(absolute.signal){param.fit = c(param.fit, gamma.m1)}
	}
	if(model > 1)
	{
		param.fit = make.optimization(T = T, i = gene.index, model = model, Nfit = NA, debug = debug,  parametrization = parametrization,method = method, sum.species = sum.species, zt = zt, i.ex = i.ex, i.int = i.int, absolute.signal = absolute.signal)
	}
	return(param.fit)
}


error = function(S,s,M,m, log = TRUE)
{

	S = as.numeric(unlist(S)); 
	s = as.numeric(unlist(s)); 
	M =  as.numeric(unlist(M)); 
	m =  as.numeric(unlist(m));
	## here we have to use the estimated values of sigma.ss and sigma.mm
	sigma.ss = 0.2149748;
	sigma.mm = 0.1760124;
	
	if(log)
	{
		if(any(S<=0)|any(s<=0)){
			error.S = 10^8;
		}else{
			error.S = sum((log(S)-log(s))^2);
			error.S = error.S/(sigma.ss^2);
			
		}
	}else{
		error.S = sum((S-s)^2)
	}
	if(log)
	{
		if(any(M<=0)|any(m<=0))
		{
			error.M = 10^8
		}else{
			error.M = sum((log(M)-log(m))^2);
			error.M = error.M/(sigma.mm^2);
		}
	}else{
		error.M = sum((M-m)^2)
	}
	error = error.S + error.M 
	
	if(intense.debug){cat('______________________________-------- error.S = ',error.S,'\n')}
	if(intense.debug){cat('______________________________-------- error.M = ',error.M,'\n')}

	return(error)  
}

errors = function(S = rep(1,24),s = rep(1,24),M = rep(1,24),m = rep(1,24),param = NA,model = 1, parametrization =c('cosine','sigmoid'),  log = TRUE)
{
	parametrization = parametrization[1]; 
	
	alpha = 50;
	if(parametrization == 'cosine')
	{
		constrain.s.positive = min(sum(exp(-alpha*s)),10^8)
	}else{
		constrain.s.positive = 0
	}
	if(intense.debug){cat('______________________________-------- constrain.s.positive = ', constrain.s.positive,'\n')}
	
	## error function
	if((parametrization == 'sigmoid')&((model == 2)|(model == 4)))
		{
			up.time = param[model+2]; 
			down.time = param[model+3];
			## here is to constrain the sum of uptime and downtime <24.1h
			constrain.param = min(sum(exp(alpha*(up.time+down.time-24.1))),10^8)
		}else{
			constrain.param = 0
		}
	if(intense.debug){cat('______________________________-------- constrain.param = ', constrain.param,'\n')}

	error = error(S,s,M,m,log = log)
	errors = error + constrain.s.positive + constrain.param
	return(errors)	
}



make.optimization = function(T = T, i = 1, model = 2, Nfit = NA, debug = FALSE, parametrization =c('cosine','sigmoid'), method =  c('integration','simulation'), sum.species = FALSE, zt =  seq(0,46,by = 2), i.ex = ZT.ex, i.int = ZT.int, absolute.signal = FALSE)
{
	parametrization = parametrization[1];
	method = method[1]
	
	param.fit = NA
	S = unlist(T[i, i.int])
	M = unlist(T[i, i.ex])
	
	bounds = set.bounds(model = model, parametrization = parametrization, absolute.signal = absolute.signal); 
	upper = bounds$upper; 
	lower = bounds$lower;
	
	# Nfit the number of optimization which is twice of nb of parameters and each optimization is from different initial conditions
	if(is.na(Nfit)){Nfit = 2*length(upper)}
	if(model==3) {set.global.variable(); k.splicing=splicing.rate; index.exon = grep('exon.median', colnames(T)); index.intron = grep('intron.median', colnames(T)); prior.factor.gamma=2*k.splicing*T[i,index.intron]/T[i, index.exon];} 
	
	#Nfit times of initialization
	PAR.INIT = do.parameter.initialization(T = T, S = S, M = M, model = model, index = i, Nfit = Nfit, parametrization = parametrization, absolute.signal = absolute.signal, zt = zt)
	
	#if(model==3) Nfit = Nfit*2;
	
	errors.fit = rep(NA, Nfit)
	debug = TRUE
	for(fit.number in 1:Nfit)
	{
		par.init = PAR.INIT[fit.number,]
		if(debug){cat('\t\t starting optimization # ',fit.number,' for model ', model,':',par.init,'\n')}
		
		## the optimization function used in the fitting
		## f2min is the likelihood function to minimize
		opt = optim(par.init, f2min, M = M, S = S, prior.factor.gamma= prior.factor.gamma, model = model, parametrization = parametrization, method.intern = method, debug = debug, sum.species = sum.species, zt = zt, absolute.signal = absolute.signal,  method = 'L-BFGS-B', lower = lower, upper = upper)
		
		## extract the parameter of optimization results and errors
		res.fit = opt$par
		if(model==3) {res.fit[1]=prior.factor.gamma/(1+res.fit[2])}
		errors.fit[fit.number] = opt$value
		eval(parse(text = paste('res.fit.', fit.number, ' = res.fit', sep = '')))
		if(debug){cat('\t\t optimization # ',fit.number,' finished : \t\t\t',res.fit,'\t',opt$value,'\n')}
	}
	
	imin = which.min(errors.fit); 
	eval(parse(text = paste('res.fit = res.fit.', imin, sep = '')))
	## save the minimal error and corresponding parameters
	param.fit = c(errors.fit[imin], res.fit);
	names(param.fit) = paste(c('error',colnames(PAR.INIT)),'.m',model,sep = '')
	return(param.fit);
	if(debug){cat('\t\t\t optimization finished\n')}
}




set.bounds = function(model = 2, parametrization =c('cosine','sigmoid'), absolute.signal = FALSE)
{
	parametrization = parametrization[1]

	gamma.max = log(2)/10*60 ; 
	gamma.min = log(2)/24
	phase.gamma.max = 24; 
	phase.gamma.min = 0
	
	if(parametrization == 'cosine')
	{
		eps.gamma.max = 1; eps.gamma.min = 0; 
		eps.24.S.max = 1; eps.24.S.min = 0
		phase.24.S.max = 24; phase.24.S.min = 0;
		eps.12.S.max = 1; eps.12.S.min = 0
		phase.12.S.max = 12; phase.12.S.min = 0; 
		param.synthesis.upper = c(eps.24.S.max, phase.24.S.max, eps.12.S.max, phase.12.S.max)
		param.synthesis.lower = c(eps.24.S.min, phase.24.S.min, eps.12.S.min, phase.12.S.min)
	}else
	{
		eps.gamma.max = 100; 
		eps.gamma.min = 1; 
		fold.change.max = 1000; 
		fold.change.min = 1;
		phase.max = 24; 
		phase.min = 0;
		up.time.max = 22; 
		up.time.min = 2; 
		if(absolute.signal){up.time.max = 20; up.time.min = 4}
		down.time.max = 22;
		down.time.min = 2; 
		if(absolute.signal){down.time.max = 20; down.time.min = 4}
		param.synthesis.upper = c(fold.change.max, phase.max, up.time.max, down.time.max)
		param.synthesis.lower = c(fold.change.min, phase.min, up.time.min, down.time.min)
	}
	
	if(model == 1){upper = c(0); lower = c(0)}
	if(model == 2)
	{
		upper = c(gamma.max, param.synthesis.upper)
		lower = c(gamma.min, param.synthesis.lower)
	}
	if(model == 3)
	{
		upper = c(gamma.max, eps.gamma.max, phase.gamma.max)
		lower = c(gamma.min, eps.gamma.min, phase.gamma.min)
	}
	if(model == 4)
	{
		upper = c(gamma.max, eps.gamma.max, phase.gamma.max, param.synthesis.upper)
		lower = c(gamma.min, eps.gamma.min, phase.gamma.min, param.synthesis.lower)
	}
	return(list(lower = lower, upper = upper))
}

do.parameter.initialization = function(T = T, S = rep(1,24),M = rep(1,24),model = 2, index = 1, Nfit = 5,parametrization =c('cosine','sigmoid'), absolute.signal = FALSE, zt = seq(0,46,by = 2))
{
	parametrization = parametrization[1]
	
	bounds = set.bounds(model = model, parametrization = parametrization, absolute.signal = absolute.signal); 
	upper = bounds$upper; 
	lower = bounds$lower;
	
	#gamma.init = rep(log(2)/120*60, Nfit)# initial half-life = 120 min
	gamma.init = c(rep(log(2)/lseq(log(2)/upper[1],log(2)/lower[1], length = 3), each = Nfit%/%3),rep(log(2)/5,Nfit%%3))
	eps.gamma.init = rep(0.2, Nfit); if(parametrization == 'sigmoid'){eps.gamma.init = (1+ eps.gamma.init)/(1-eps.gamma.init); eps.gamma.init = eps.gamma.init*sample(seq(1/max(eps.gamma.init),2-(1/max(eps.gamma.init)),len = 1000),Nfit)}
	phases.gamma.init = c(rep(c(6,18), round(Nfit/2)), rep(12, Nfit%%2))
	
	if(parametrization == 'cosine')
	{
		## cosine function
		eps.24.S.init = rep(T$rel.ampl.int[index],Nfit)
		phase.24.S.init = (rep(T$phase.int[index],Nfit)+rnorm(Nfit,sd = 2))%%24
		eps.12.S.init = rep(T$rel.ampl.12.int[index],Nfit) 
		phase.12.S.init = rep(T$phase.12.int[index],Nfit)
		param.init.synthesis = cbind(eps.24.S.init,phase.24.S.init,eps.12.S.init,phase.12.S.init)
		param.init.synthesis.name = c('rel.ampl.int','phase.int','rel.ampl.12.int','phase.12.int')
	}else{
		# sigmoid function
		# here the column of rel.ampl.int of T is used
		if(!is.null(T$rel.ampl.int[index])){rel.ampl.int = T$rel.ampl.int[index]; phase.int = T$phase.int[index]}else{rel.ampl.int = (max(S)-min(S))/2; phase.int = zt[which.max(S)]}
		fold.change.init = rep(min((1+ rel.ampl.int)/(1-rel.ampl.int),1000),Nfit); fold.change.init = fold.change.init*sample(seq(1/max(fold.change.init),2-(1/max(fold.change.init)),len = 1000),Nfit)
		phase.init = (rep(phase.int,Nfit)+rnorm(Nfit,sd = 2)+c(rep(0,(Nfit%/%3)*3), rep(12,Nfit%%3)))%%24
		if(model!=3){times.min = lower[model+2];times.max = upper[model+2]}else{times.min = 2; times.max = 22}
		up.time.init = c(rep(seq(times.min, times.max, length = Nfit%/%3),3),rep(mean(c(times.min,times.max)),Nfit%%3));  down.time.init =  c(rep(seq(times.max , times.min, length = Nfit%/%3),3),rep(mean(c(times.min,times.max)),Nfit%%3))
		param.init.synthesis = cbind(fold.change.init, phase.init, up.time.init, down.time.init)
		param.init.synthesis.name = c('fold.change.int','phase.int','up.time.int','down.time.int')
	}
	
	cat('parameter initialized\n')

	if(model == 2){PAR.INIT = cbind(gamma.init, param.init.synthesis); colnames(PAR.INIT) = c('gamma', param.init.synthesis.name)}
	if(model == 3){PAR.INIT = cbind(gamma.init, eps.gamma.init, phases.gamma.init); colnames(PAR.INIT) = c('gamma','eps.gamma','phase.gamma');}
	if(model == 4){PAR.INIT = cbind(gamma.init, eps.gamma.init, phases.gamma.init, param.init.synthesis); colnames(PAR.INIT) = c('gamma','eps.gamma','phase.gamma', param.init.synthesis.name)}
	return(PAR.INIT)	
}

### This is the likelihood function (or error function) to minimize in the optimization
### also in this function the sum or not sum is considered.

f2min = function(par, M, S, prior.factor.gamma= parior.factor.gamma, model, parametrization =c('cosine','sigmoid'), method.intern =  c('integration','simulation'), debug = FALSE, sum.species = FALSE,  zt = seq(0,46,by = 2), absolute.signal = FALSE)
{
	
	parametrization = parametrization[1]; 
	method = method.intern[1]
	if(intense.debug){cat('______________________________in f2min - model',model,' \n')}
	if(intense.debug){cat(par,'\n')}
	
	# initialization of the parameters
	
	gamma = par[1]; 
	eps.gamma =0 ; 
	phase.gamma = 0
	eps.24.S = 0; 
	phase.24.S = 0
	eps.12.S = 0; 
	phase.12.S = 0
	param.synthesis.1 = 0; 
	param.synthesis.2 = 0; 
	param.synthesis.3 = 0; 
	param.synthesis.4 = 0;
	
	if(parametrization == 'sigmoid'){eps.gamma =1; param.synthesis.1 = 1;  param.synthesis.3 = 1; param.synthesis.4 = 1}
	
	if(model > 2){eps.gamma = par[2]; phase.gamma = par[3];}
	if((model == 2)|(model == 4)){j = model; param.synthesis.1 = par[j]; param.synthesis.2 = par[j+1]; param.synthesis.3 = par[j+2]; param.synthesis.4 = par[j+3]}
	if(model==3) {gamma = prior.factor.gamma/(1+eps.gamma)}
	#pre-mRNA profile
	#caculated directly with parameters of sigmoid function for different time points (zt)
	if(parametrization == 'cosine')
	{
		s = compute.s(t = zt, eps.24.S = param.synthesis.1, phase.24.S = param.synthesis.2, eps.12.S = param.synthesis.3, phase.12.S =  param.synthesis.4)
	}else{
		s = compute.sigmoid(t = zt, fold.change = param.synthesis.1, phase = param.synthesis.2, up.time = param.synthesis.3, down.time =  param.synthesis.4)
	}

	#mRNA profile
	#caculated by simulation or integration
	if(method == 'simulation')
	{
		## 'simulation' method to find mrna profile
		m = simulate.m(t = zt, par = par, model = model, parametrization = parametrization)
		
		if(sum.species|absolute.signal)
		{
			set.global.variable();
			syn.f = splicing.rate; 
			rescale = FALSE
		}else{
			syn.f  = 1; 
			rescale = TRUE
		}
		if(rescale){mean = mean(m); m = m/mean;}else{m = syn.f * m}
		
		
	}else{
		if(parametrization == 'cosine')
		{
			if(model == 2){
				if(intense.debug){cat('______________________________-------- should be in model 2\n')};
				m = compute.m(t = zt, gamma,  eps.24.S = param.synthesis.1, phase.24.S = param.synthesis.2, eps.12.S = param.synthesis.3, phase.12.S =  param.synthesis.4);
			}else{
				if(intense.debug){cat('______________________________-------- should be in model 3 or 4\n')};
				m = compute.m.changing(t = zt, gamma = gamma,eps.gamma = eps.gamma, phase.gamma = phase.gamma, eps.24.S = param.synthesis.1, phase.24.S = param.synthesis.2, eps.12.S = param.synthesis.3, phase.12.S =  param.synthesis.4)}
		}else{
			## integration method to find mrna profile
			zt.for.sigmoid = zt; 
			if(((max(zt)-min(zt))>24) & (max(zt)%%24 == 24-zt[2]+zt[1])){zt.for.sigmoid = zt[1:(length(zt)/2)]}
			
			if(sum.species|absolute.signal)
			{
				set.global.variable();
				syn.f = splicing.rate; 
				rescale = FALSE
			}else{
				syn.f  = 1; 
				rescale = TRUE
			}
			## HERE DEAL WITH SUM OR NOT SUM
			m = compute.m.sigmoid(t = zt.for.sigmoid, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, 
								  fold.change = param.synthesis.1, phase = param.synthesis.2, up.time = param.synthesis.3, down.time =  param.synthesis.4, 
								  rescale = rescale, synthesis.factor = syn.f)
			
			if(length(zt.for.sigmoid)!=length(zt)){m = rep(m,max(zt)%/%24+1)}
			if(sum.species&!absolute.signal){m = (m+s)/mean(m+s)}else if(absolute.signal){m = m+s}
			
			################################################################################################################ finishing line of f2min function
			## test the indetermination of model 3
			indetermination.test.model3 = FALSE
			if(indetermination.test.model3)
			{
			test = 201
			if(parametrization == 'sigmoid'){eps.gamma =1; param.synthesis.1 = 1;  param.synthesis.3 = 1; param.synthesis.4 = 1}
			# true parameters
			gamma = T$gamma[test]; eps.gamma = T$eps.gamma[test];phase.gamma=T$phase.gamma[test];
			print(c(gamma, eps.gamma, phase.gamma))
			
			zt.for.sigmoid = zt; 
			if(((max(zt)-min(zt))>24) & (max(zt)%%24 == 24-zt[2]+zt[1])){zt.for.sigmoid = zt[1:(length(zt)/2)]}
			
			if(sum.species|absolute.signal)
			{
				set.global.variable();
				syn.f = splicing.rate; 
				rescale = FALSE
			}else{
				syn.f  = 1; 
				rescale = TRUE
			}
			
			m = compute.m.sigmoid(t = zt.for.sigmoid, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, 
								  fold.change = param.synthesis.1, phase = param.synthesis.2, up.time = param.synthesis.3, down.time =  param.synthesis.4, 
								  rescale = rescale, synthesis.factor = syn.f)
			
			if(length(zt.for.sigmoid)!=length(zt)){m = rep(m,max(zt)%/%24+1)}
			if(sum.species&!absolute.signal){m = (m+s)/mean(m+s)}else if(absolute.signal){m = m+s}
			
			xx = m
			
			## inferred paramters
			gamma = T$gamma.m3[test]; eps.gamma = T$eps.gamma.m3[test];phase.gamma=T$phase.gamma.m3[test];
			print(c(gamma, eps.gamma, phase.gamma))
			
			zt.for.sigmoid = zt; 
			if(((max(zt)-min(zt))>24) & (max(zt)%%24 == 24-zt[2]+zt[1])){zt.for.sigmoid = zt[1:(length(zt)/2)]}
			
			if(sum.species|absolute.signal)
			{
				set.global.variable();
				syn.f = splicing.rate; 
				rescale = FALSE
			}else{
				syn.f  = 1; 
				rescale = TRUE
			}

			m = compute.m.sigmoid(t = zt.for.sigmoid, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, 
								  fold.change = param.synthesis.1, phase = param.synthesis.2, up.time = param.synthesis.3, down.time =  param.synthesis.4, 
								  rescale = rescale, synthesis.factor = syn.f)
			
			if(length(zt.for.sigmoid)!=length(zt)){m = rep(m,max(zt)%/%24+1)}
			if(sum.species&!absolute.signal){m = (m+s)/mean(m+s)}else if(absolute.signal){m = m+s}
			
			yy = m
			M = T[test, c(2:25)]
			
			plot(c(0:23)*2, M, type='b',col='black')
			points(c(0:23)*2, xx, type='b',col='blue')
			points(c(0:23)*2, yy, type='b',col='red')
			abline(h=1, col='gray',lwd=2.0)
			
			kk = c(201:300)
			gamma.true = T$gamma[kk]
			
			k.splicing = log(2)/8*60
			prior.gamma=k.splicing*T[kk,index.intron]/T[kk,index.exon]*2/(1+T[kk,]$eps.gamma.m3)
			plot(gamma.true, prior.gamma, log='xy')
			abline(0,1,col='red')
			}
			
			
			
		}
	}
		
	err.fit = errors(S = S,s = s,M = M,m = m,param = par, model = model,parametrization = parametrization, log = TRUE); 
	if(intense.debug){cat('______________________________-------- error = ',error,'\n')}
	return(err.fit)
}



compute.s = function(t = zt, eps.24.S = 1, phase.24.S = 0, eps.12.S = 0, phase.12.S = 0){
	w = 2*pi/24
	s = 1 + eps.24.S * cos((t-phase.24.S)*w) + eps.12.S * cos((t-phase.12.S)*2*w)
	return(s)
	}
	
compute.m = function(t = zt, gamma = log(2), eps.24.S = 1, phase.24.S = 0, eps.12.S = 0, phase.12.S = 0){
	w = 2*pi/24
	m = 1 + eps.24.S * gamma /sqrt(w^2 + gamma^2) * cos((t-phase.24.S)*w - atan(w/gamma)) + eps.12.S * gamma /sqrt((2*w)^2 + gamma^2) * cos((t-phase.12.S)*2*w - atan((2*w)/gamma))
	return(m)
	}


compute.m.changing = function(t = zt, gamma = log(2), eps.gamma = 0, phase.gamma = 0, eps.24.S = 1, phase.24.S = 0, eps.12.S = 0, phase.12.S = 0){
	w <<- 2*pi/24
	m = integrate.m(t = t, w = w, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, eps.24.S = eps.24.S, phase.24.S = phase.24.S, eps.12.S = eps.12.S, phase.12.S = phase.12.S)
	return(m)
	}
	
integrate.m = function(t = zt, w = w, gamma = log(2), eps.gamma = 0, phase.gamma = 0, eps.24.S = 1, phase.24.S = 0, eps.12.S = 0, phase.12.S = 0){
	m = rep(0, length(t))
	for(i in 1:length(t)){
		time = t[i]; Tstable = Tstable =  24*(ceiling(log(2000)/gamma/24) +1)
		#if(intense.debug){cat('time :',time,' \t\t')}
		m[i] = exp(-Gamma(t = Tstable+time,w = w, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma)) * integrate(f2integrate,lower = 0, upper = Tstable+time, par = c(gamma, eps.gamma, phase.gamma, eps.24.S, phase.24.S, eps.12.S, phase.12.S))$value	}
	mean = mean(m); m = m/mean
}

f2integrate = function(x, par)
{
	if(intense.debug){cat('in f2integrate\n')}
	gamma = par[1]; eps.gamma= par[2]; phase.gamma= par[3]; eps.24.S= par[4]; phase.24.S= par[5]; eps.12.S= par[6]; phase.12.S = par[7] 
	return(exp(Gamma(t = x,w = w, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma)) * compute.s(t = x, eps.24.S = eps.24.S, phase.24.S = phase.24.S, eps.12.S = eps.12.S, phase.12.S = phase.12.S))
}

Gamma = function(t = 0, w = 2*pi/24, gamma = log(2), eps.gamma = 0, phase.gamma = 0)
{
	Gamma = gamma*(t+ eps.gamma/w * sin(w*(t-phase.gamma)))
	return(Gamma)
}	
	
	
compute.sigmoid = function(t = seq(0,46,by = 2), fold.change = 2, phase = 12, up.time = 8, down.time = 8, rescale = TRUE)
{
	profile = compute.sigmoid.intern(t = t%%48, fold.change = fold.change, phase = phase, up.time = up.time, down.time = down.time)
	if(rescale)
	{
		#profile.for.mean = compute.sigmoid.intern(t = seq(0,23,by =1), fold.change = fold.change, phase = phase, up.time = up.time, down.time = down.time)
		#mean = mean(profile.for.mean)
		mean = (24 + (fold.change-1)/2*(up.time+down.time))/24; ## rescale the profile so that the average of profile is 1 (ex. for intronic singals)
		x = profile/mean
	}else{x = profile}
	return(x)
}

compute.sigmoid.intern = function(t = seq(0,46,by = 2), fold.change = 2, phase = 12, up.time = 8, down.time = 8)
{
	
	up.time = up.time/8;
	down.time = down.time/8
	
	#P = -4
	#for(i in -1:2){P = P + 1/(1+exp(-(t-phase - i*24 +4* up.time)/up.time)) + 1/(1+exp((t-phase - i*24 -4* down.time)/down.time))}
	
	#P = -4 + 1/(1+exp(-(t-phase + 24 +4* up.time)/up.time)) + 1/(1+exp((t-phase + 24 -4* down.time)/down.time)) + 1/(1+exp(-(t-phase +4* up.time)/up.time)) + 1/(1+exp((t-phase -4* down.time)/down.time)) + 1/(1+exp(-(t-phase - 24 +4* up.time)/up.time)) + 1/(1+exp((t-phase - 24 -4* down.time)/down.time)) + 1/(1+exp(-(t-phase - 48 +4* up.time)/up.time)) + 1/(1+exp((t-phase - 48 -4* down.time)/down.time))
	
	x = 1+(fold.change-1)*(-4 + 1/(1+exp(-(t-phase + 24 +4* up.time)/up.time)) + 1/(1+exp((t-phase + 24 -4* down.time)/down.time)) + 1/(1+exp(-(t-phase +4* up.time)/up.time)) + 1/(1+exp((t-phase -4* down.time)/down.time)) + 1/(1+exp(-(t-phase - 24 +4* up.time)/up.time)) + 1/(1+exp((t-phase - 24 -4* down.time)/down.time)) + 1/(1+exp(-(t-phase - 48 +4* up.time)/up.time)) + 1/(1+exp((t-phase - 48 -4* down.time)/down.time)))
	return(x)
	
}
	
	
compute.m.sigmoid = function(t = zt, gamma = 2, eps.gamma = 3, phase.gamma = 18, fold.change = 2, phase = 12, up.time = 8, down.time = 8, rescale = TRUE, synthesis.factor = 1, integrate.function = TRUE, simulate = TRUE, simulate.only = FALSE)
{
	
	m = rep(0, length(t)); 	
	w <<- 2*pi/24;
	Tstable =  24*(ceiling(log(2000)/gamma/24) +1)
	
	t = t%%24

	if(!simulate.only){
	if(!integrate.function){
		#cat('alternative integration \n')
		G1 = Gamma.sigmoid(t = Tstable+t, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma)
		eG1 = exp(-G1)
		dt = 0.004; fact = 500; if(length(t)>1){dt.original = t[2]-t[1]; dt = dt.original/fact};
		time.sup = seq(t[1],Tstable+t[length(t)], by = dt)
		prod = exp(Gamma.sigmoid(t = time.sup,gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma)) * compute.sigmoid(t = time.sup,fold.change = fold.change, phase = phase, up.time = up.time, down.time = down.time, rescale = !rescale)
		integral = (cumsum(prod)-prod[1])*dt; integral.t = integral[seq(length(integral)-fact*(length(t)-1),length(integral), by = fact)]
		m = eG1*integral.t
		cat(m, '\n')
		
		prod = exp(Gamma.sigmoid(t = time.sup,gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma)-Gamma.max) * compute.sigmoid(t = time.sup,fold.change = fold.change, phase = phase, up.time = up.time, down.time = down.time, rescale = !rescale)
		integral = (cumsum(prod)-prod[1])*dt; integral.t = integral[seq(length(integral)-fact*(length(t)-1),length(integral), by = fact)]
		m = exp(-G1 + Gamma.max + log(integral.t))
		cat(m, '\n')				
		}
	else{
		#cat('integrate function \n')	
		Gamma.T = m
		for(i in 1:length(t))
		{
			time = t[i]; 	
			Gamma.t = Gamma.sigmoid(t = Tstable+time,gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma) #+max(t)
			Gamma.T[i] = Gamma.t
			#if(intense.debug){cat('time :',time,' == ',Tstable+time, '\t\t')}
			m[i] = tryCatch({
				mi =integrate(f2integrate.sigmoid.2,lower = 0, upper = Tstable+time, par = c(gamma, eps.gamma, phase.gamma, fold.change, phase, up.time, down.time), rescale = !rescale, Gamma.t = Gamma.t, abs.tol = 10^-16);
				mi = mi$value
				},warning = function(war){},error = function(err){cat('error in integration\n');mi = 'NaN'},finally = {})
				#print(m[i])
			if(any((m == 'Inf')|(m == 'NaN')|is.na(m))){break}	
		}
	}}
	
	#mm = m
	#simulate = TRUE
	#if(simulate)
	if((any((m == 'Inf')|(m == 'NaN')|is.na(m)|(m<10^-100))) & simulate)
	{ 
		#(any((m == 'Inf')|(m == 'NaN')|is.na(m)|(m==0)))&
		#cat('simulate\n')
		m = simulate.m(t = t, par = c(gamma, eps.gamma, phase.gamma, fold.change, phase, up.time, down.time), parametrization = 'sigmoid');
		#print('here')
	}
	
	if(rescale){mean = mean(m); m = m/mean;}
	else{m = synthesis.factor * m}
	return(m)
}


f2integrate.sigmoid = function(x, par, rescale = FALSE)
{
	# if(intense.debug){cat('in f2integrate.sigmoid --- x :')}
	gamma = par[1]; 
	eps.gamma= par[2]; 
	phase.gamma= par[3]; 
	fold.change=  par[4]; 
	phase = par[5]; 
	up.time = par[6]; 
	down.time = par[7] 
	value = exp(Gamma.sigmoid(t = x,gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma)) * compute.sigmoid(t = x,fold.change = fold.change, phase = phase, up.time = up.time, down.time = down.time, rescale = rescale)
	# if(intense.debug){cat(x,'\n')}

	return(value)
}
	
f2integrate.sigmoid.2 = function(x, par, rescale = FALSE, Gamma.t = 0)
{
	# if(intense.debug){cat('in f2integrate.sigmoid --- x :')}
	gamma = par[1]; 
	eps.gamma= par[2]; 
	phase.gamma= par[3]; 
	fold.change=  par[4]; 
	phase = par[5];
	up.time = par[6]; 
	down.time = par[7] 
	value = exp(Gamma.sigmoid(t = x,gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma)-Gamma.t) * compute.sigmoid(t = x,fold.change = fold.change, phase = phase, up.time = up.time, down.time = down.time, rescale = rescale)
	# if(intense.debug){cat(x,'\n')}

	return(value)
}	

Gamma.sigmoid = function(t = 0, gamma = log(2), eps.gamma = 0, phase.gamma = 0)
{
	T = t%/%24; 
	tt = t%%24
	Gamma = gamma * ( t + 24* T * (eps.gamma-1)/2 +  (eps.gamma-1) * 3/2 * (Log(t = tt, i = 0, p = phase.gamma)+ Log(t = tt, i = 1, p = phase.gamma) + Log(t = tt, i = -1, p = phase.gamma)) )
	return(Gamma)
}

Log = function(t = 0, i = 0, p = 0)
{
	Log = log((exp(-2/3*(p+i*24+6)) + 1)/(exp(-2/3*(p+i*24-6)) + 1) * (exp(-2/3*(p+i*24-6 - t)) + 1) / (exp(-2/3*(p+i*24+6- t)) + 1))
	return(Log)	
}
	
	
#Log.alt = function(t = 0, i = 0, p = 0){
#	Log = log(exp(-2/3*(p+i*24+6)) + 1) - log(exp(-2/3*(p+i*24-6)) + 1) + log(exp(-2/3*(p+i*24-6 - t)) + 1) -log(exp(-2/3*(p+i*24+6- t)) + 1) 
#	return(Log)	
#	}	
#	
#Gamma.sigmoid.alt =function(t = 0, gamma = log(2), eps.gamma = 0, phase.gamma = 0){
#	time.res = 0.1; factor = 10
#	if(length(t) != 1){time.res.t = t[2]-t[1]; time.res = time.res.t/factor}
#	gamma.sigm = compute.sigmoid(t = seq(0,max(t),by = time.res), fold.change = eps.gamma, phase = phase.gamma, up.time = 12, down.time = 12, rescale = FALSE)
#	if(length(t) == 1){ind = length(gamma.sigm)}else{ind = seq(1,length(gamma.sigm),by = factor)}
#	
#	Gamma = gamma * (cumsum(gamma.sigm) - gamma.sigm[1])[ind]*time.res
#	return(Gamma)
#}	


dmdt = function(t, y, par, parametrization = c('cosine','sigmoid'))
{
	m = y
	parametrization = parametrization[1]
	gamma = par[1];
	eps.gamma = par[2]; 
	phase.gamma = par[3];
	param.synthesis.1 = par[4]; 
	param.synthesis.2 = par[5]; 
	param.synthesis.3 = par[6]; 
	param.synthesis.4 = par[7]
	
	if(parametrization == 'cosine'){s.t = compute.s(t = t, eps.24.S = param.synthesis.1, phase.24.S = param.synthesis.2, eps.12.S = param.synthesis.3, phase.12.S =  param.synthesis.4)}
	else{s.t = compute.sigmoid(t = t, fold.change = param.synthesis.1, phase = param.synthesis.2, up.time = param.synthesis.3, down.time =  param.synthesis.4 )}#cat('here',"param.synthesis.1", param.synthesis.1,"param.synthesis.2", param.synthesis.2,"param.synthesis.3", param.synthesis.3,"param.synthesis.4", param.synthesis.4,'\n'); 
	
	#points(t,s.t)
	
	#cat('time :',t, '\t s(t) :',s.t,'\n')
	if(parametrization == 'cosine'){gamma.t = gamma * (1 + eps.gamma*cos((t-phase.gamma)/24*2*pi))}
	else{gamma.t = gamma*compute.sigmoid(t = t, fold.change= eps.gamma, phase = phase.gamma, up.time = 12, down.time = 12, rescale = FALSE)}
	dmdt = s.t - gamma.t * m
	#cat("t:",t," dmdt",dmdt," s.t",s.t," m",m," gamma.t",gamma.t,'\n')
	list(dmdt,NULL)
}


simulate.m = function(t, par, model=4, parametrization)
{
	gamma = par[1];
	if(parametrization == 'sigmoid'){eps.gamma =1; param.synthesis.1 = 1;  param.synthesis.3 = 1; param.synthesis.4 = 1}
	
	if(model > 2){eps.gamma = par[2]; phase.gamma = par[3];}
	if((model == 2)|(model == 4)){j = model; param.synthesis.1 = par[j]; param.synthesis.2 = par[j+1]; param.synthesis.3 = par[j+2]; param.synthesis.4 = par[j+3]}
	par = c(gamma, eps.gamma, phase.gamma, param.synthesis.1, param.synthesis.2, param.synthesis.3,param.synthesis.4)
	
	Tstable =  24*(ceiling(log(2000)/gamma/24) + 1)
	t.res = 2; 
	#if(length(t)!=1){t.res = t[2]-t[1]}
	t.sup = seq(0, Tstable+max(t) ,by= t.res)
	
	#plot(t.sup,t.sup*0+1, type = 'n')
	soln = lsoda(dmdt,
			times= t.sup, ## times
			y = 0, #init.conditions
			par=par,parametrization = parametrization) ## parameter values
	
	#soln[match(t+48, soln[,1]),2]
	
	i.last = nrow(soln); 
	i.keep = seq(i.last - length(t)+1,i.last,by = 1)
	#plot(soln[,1],soln[,2], type = 'l')
	m = soln[i.keep,2]; 
	#mean = mean(m);m = m/mean
	#cat(soln[,2],'\n')
	#cat(parametrization,'\n')
	return(m)
}

dpdmdt = function(t, y, par, parametrization = c('cosine','sigmoid'))
{
	p = y[1];m = y[2];
	parametrization = parametrization[1]
	s = par[1]; k = par[2]; k.fc = par[3]; k.phase = par[4]
	gamma = par[1]; eps.gamma = par[2]; phase.gamma = par[3];param.synthesis.1 = par[4]; param.synthesis.2 = par[5]; param.synthesis.3 = par[6]; param.synthesis.4 = par[7]
	
	if(parametrization == 'cosine'){s.t = compute.s(t = t, eps.24.S = param.synthesis.1, phase.24.S = param.synthesis.2, eps.12.S = param.synthesis.3, phase.12.S =  param.synthesis.4)}
	else{s.t = compute.sigmoid(t = t, fold.change = param.synthesis.1, phase = param.synthesis.2, up.time = param.synthesis.3, down.time =  param.synthesis.4 )}
	
	#points(t,s.t)
	
	#cat('time :',t, '\t s(t) :',s.t,'\n')
	if(parametrization == 'cosine'){gamma.t = gamma * (1 + eps.gamma*cos((t-phase.gamma)/24*2*pi))}
	else{gamma.t = gamma*compute.sigmoid(t = t, fold.change= eps.gamma, phase = phase.gamma, up.time = 12, down.time = 12, rescale = FALSE)}
	dmdt = s.t - gamma.t * m
	list(dmdt,NULL)
}
	
	

simulate.m.splicing = function(t,par, parametrization)
{
	gamma = par[1];
	
	Tstable =  24*(ceiling(log(2000)/gamma/24) + 1)
	t.res = 2; if(length(t)!=1){t.res = t[2]-t[1]}
	t.sup = seq(0, Tstable+max(t) ,by= t.res)
	
	#plot(t.sup,t.sup*0+1, type = 'n')
	soln = lsoda(dmdt,
			times= t.sup, ## times
			y = 0, #init.conditions
			par=par,parametrization = parametrization) ## parameter values
	
	i.last = nrow(soln); i.keep = seq(i.last - length(t)+1,i.last,by = 1)
	#plot(soln[,1],soln[,2], type = 'l')
	m = soln[i.keep,2]; #mean = mean(m); m = m/mean
	return(m)
}








######################
# Functions for the Likelihood Ratio Test
######################
vuong.test.nonest = function(xx, parametrization =c('cosine','sigmoid'), sigma=0.1,method.intern =  c('integration','simulation'), debug = FALSE, sum.species = FALSE,  zt = seq(0,46,by = 2), absolute.signal = FALSE) 
{
	#T=T; parametrization =c('cosine','sigmoid');sigma=0.1;method.intern =  c('integration','simulation'); debug = FALSE; sum.species = FALSE;  zt = seq(0,46,by = 2); absolute.signal = FALSE;
	parametrization = parametrization[1];
	method = method.intern[1]
	#xx = TT[1,]
	#if(intense.debug){cat('______________________________in f2min - model',model,' \n')}
	#if(intense.debug){cat(par,'\n')}
	
	# loop for each gene
	#first for model 2
	#model = 2
		n.par = 48
		gamma = xx[n.par+1]; 
		eps.gamma =1;
		phase.gamma = 0;
		param.synthesis.1 = xx[n.par+2]; param.synthesis.2 = xx[n.par+3]; param.synthesis.3 = xx[n.par+4]; param.synthesis.4 = xx[n.par+5];
			
		s2 = compute.sigmoid(t = zt, fold.change = param.synthesis.1, phase = param.synthesis.2, up.time = param.synthesis.3, down.time =  param.synthesis.4)
		
		zt.for.sigmoid = zt; 
		if(((max(zt)-min(zt))>24) & (max(zt)%%24 == 24-zt[2]+zt[1])){zt.for.sigmoid = zt[1:(length(zt)/2)]}
		if(sum.species|absolute.signal)
		{
			set.global.variable();
			syn.f = splicing.rate; 
			rescale = FALSE
		}else{
			syn.f  = 1; 
			rescale = TRUE
		}
		## HERE DEAL WITH SUM OR NOT SUM again
		m = compute.m.sigmoid(t = zt.for.sigmoid, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, 
							  fold.change = param.synthesis.1, phase = param.synthesis.2, up.time = param.synthesis.3, down.time =  param.synthesis.4, 
							  rescale = rescale, synthesis.factor = syn.f)
		
		if(length(zt.for.sigmoid)!=length(zt)){m = rep(m,max(zt)%/%24+1)}
		if(sum.species&!absolute.signal){m = (m+s)/mean(m+s)}else if(absolute.signal){m = m+s}
		m2 = m
		
		#model = 3
		gamma = xx[n.par+6]; 
		eps.gamma =xx[n.par+7]; 
		phase.gamma = xx[n.par+8]; 
		param.synthesis.1 = 1;  param.synthesis.2 = 0;param.synthesis.3 = 1; param.synthesis.4 = 1
		#param.synthesis.1 = T$fold.change.int.m2[n]; param.synthesis.2 = T$phase.int.m2[n]; param.synthesis.3 = T$up.time.int.m2[n]; param.synthesis.4 = T$down.time.int.m2[n]
		
		s3 = compute.sigmoid(t = zt, fold.change = param.synthesis.1, phase = param.synthesis.2, up.time = param.synthesis.3, down.time =  param.synthesis.4)
		
		zt.for.sigmoid = zt; 
		if(((max(zt)-min(zt))>24) & (max(zt)%%24 == 24-zt[2]+zt[1])){zt.for.sigmoid = zt[1:(length(zt)/2)]}
		if(sum.species|absolute.signal)
		{
			set.global.variable();
			syn.f = splicing.rate; 
			rescale = FALSE
		}else{
			syn.f  = 1; 
			rescale = TRUE
		}
		## HERE DEAL WITH SUM OR NOT SUM again
		m = compute.m.sigmoid(t = zt.for.sigmoid, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, 
							  fold.change = param.synthesis.1, phase = param.synthesis.2, up.time = param.synthesis.3, down.time =  param.synthesis.4, 
							  rescale = rescale, synthesis.factor = syn.f)
		
		if(length(zt.for.sigmoid)!=length(zt)){m = rep(m,max(zt)%/%24+1)}
		if(sum.species&!absolute.signal){m = (m+s)/mean(m+s)}else if(absolute.signal){m = m+s}
		m3 = m
		
		S = xx[25:48]
		M = xx[1:24]
		err.s2 = (log(s2)-log(S))^2
		err.m2 = (log(m2)-log(M))^2
		err2 = err.s2+err.m2
		
		err.s3 = (log(s3)-log(S))^2
		err.m3 = (log(m3)-log(M))^2
		err3 = err.s3+err.m3
		
		## here model 3 is the null model
		LR.t = 1/sigma^2*(err2-err3);
		LR.v = sd(LR.t)
		k3 = 3
		k2 = 5
		v = (sum(LR.t)-(k2/2*log(48)-k3/2*log(24)))/(sqrt(24)*LR.v)
		pvals.nonest = c(1-pnorm(v, mean = 0, sd = 1, lower.tail = FALSE), pnorm(v, mean = 0, sd = 1, lower.tail = FALSE))
		#print(pvals.nonest)
	return(pvals.nonest)
	#pvals.nonest = data.frame(pvals.nonest, stringsAsFactors=FALSE)
}

set.nb.param = function()
{
	n.param <<- c(1,5,3,7);
}


my.Likelihood.Ratio.Test = function(T = T, pval = 0.05, pval.stringent = 10^(-15), V1 = FALSE, combined = FALSE, FDR = FALSE, fdr.value = 0.2, direct = FALSE, zt = seq(0,46,by = 2), diff.sigma=TRUE)
{
	#pval.stringent = 10^(-15); V1 = FALSE; combined = FALSE; FDR = FALSE; fdr.value = 0.2; direct = FALSE; zt = seq(0,46,by = 2);Mine=TRUE
	model.orders = list(prod = c(1,2,4), deg = c(1,3,4), all = c(1,4));
	set.nb.param();
	LRP = matrix(NA, nrow = nrow(T), ncol = length(unlist(model.orders)) - length(model.orders))
	i = 1
	
	for(P in 1:length(model.orders))
	{
		order = unlist(model.orders[P])
		for(m in 1:(length(order)-1))
		{
			m1 = order[m]; 
			m2 = order[m+1]
			eval(parse(text = paste('err1 = T$error.m',m1,sep = '')))
			eval(parse(text = paste('err2 = T$error.m',m2,sep = '')))
			
			#LR = 1/sigma^2*(err1-err2);
			LR = err1 - err2;
			## LAURA old manner to calculate likelihood ratio
			#if(Old) LR = length(zt)*(log(err1)-log(err2));
			
			df = n.param[m2]-n.param[m1]
			pvals.LR = pchisq(LR,df = df ,lower.tail = FALSE);
			LRP[,i] = pvals.LR; 
			i = i+1 	
			
			#sup = seq(0,max(LR),by =0.1)
			#plot(sup, pchisq(sup, df = df), type = 'l')
			#abline(v = LR, col = rep(1:4,each = 10))
			#			
			#plot(pvals.LR, pch = 16, col = rep(1:4,each = 10))
			#plot(sort(pvals.LR), pch = 16, col = rep(1:4,each = 10), ylim = c(0,1))
		}
	
	}
	LRP = data.frame(LRP, stringsAsFactors = FALSE); 
	
	Mine = FALSE;
	if(Mine)
	{
		ZT.ex = grep('.rel.ampl.ex', colnames(T)); 
		ZT.int = grep('.rel.ampl.int', colnames(T));
		par.index.m2 = grep('.m2', colnames(T))
		par.index.m3 = grep('.m3', colnames(T))
		par.index.m2 = par.index.m2[-1]
		par.index.m3 = par.index.m3[-1]
		TT = as.matrix(T[,c(ZT.ex, ZT.int, par.index.m2, par.index.m3)])
		pvals.vuong = t(apply(TT, 1, vuong.test.nonest, parametrization ='sigmoid',sigma=sigma, method.intern = 'integration', debug = FALSE, sum.species = FALSE,  zt = seq(0,46,by = 2), absolute.signal = FALSE))
		#pvals.vuong = vuong.test.nonest(T,parametrization ='sigmoid',sigma=sigma, method.intern = 'simulation', debug = FALSE, sum.species = FALSE,  zt = seq(0,46,by = 2), absolute.signal = FALSE)
		pvals.vuong = data.frame(pvals.vuong, stringsAsFactors=FALSE)
		
		LRP$b = pvals.vuong[,1]
		LRP$B = pvals.vuong[,2]
	
		colnames(LRP) = paste('LRT.pval',c('p','P','d','D','A','nonest.m2','nonest.m3'),sep = '.')
		LRP$LRT.pval.C = pchisq(-2*log(LRP$LRT.pval.P * LRP$LRT.pval.D), df = 4, lower.tail = FALSE)
		LOOSE = LRP<pval; 
		colnames(LOOSE) =  c('p','P','d','D','A','nonest.m2','nonest.m3','C'); 
		LOOSE = data.frame(LOOSE, stringsAsFactors = FALSE)
		STRINGENT = LRP<pval.stringent;
		colnames(STRINGENT) = c('p','P','d','D','A','nonest.m2','nonest.m3','C'); 
		STRINGENT = data.frame(STRINGENT, stringsAsFactors = FALSE)
		# my method to select models
		LRP$LRT.best.model = NA
		#LRP$LRT.best.model = 1 # if it is not possible to select the model, we chose the simplest one 
		## selection with stringent criterion
		LRP$LRT.best.model[LOOSE$P & LOOSE$D & LOOSE$A] = 4
		#LRP$LRT.best.model[LOOSE$P & LOOSE$D ] = 4
		LRP$LRT.best.model[!LOOSE$p & !LOOSE$d & !LOOSE$A] = 1
		#LRP$LRT.newmodel[LOOSE$p & !LOOSE$P &(LRP$LRT.pval.nonest.m2<LRP$LRT.pval.nonest.m3)] = 2
		#LRP$LRT.newmodel[LOOSE$d & !LOOSE$D &(LRP$LRT.pval.nonest.m2>LRP$LRT.pval.nonest.m3)] = 3
		LRP$LRT.best.model[LOOSE$p & !LOOSE$P &(LRP$LRT.pval.nonest.m2<pval)] = 2
		LRP$LRT.best.model[LOOSE$d & !LOOSE$D &(pval>LRP$LRT.pval.nonest.m3)] = 3
		
		#kk = which((LRP$LRT.pval.p>=pval & LRP$LRT.pval.P<pval & LRP$LRT.pval.nonest.m2>=pval)|LRP$LRT.pval.nonest.m3>=pval)
		#LRP$LRT.best.model[is.na(LRP$LRT.best.model) & (LRP$LRT.pval.nonest.m2>=pval|LRP$LRT.pval.nonest.m3>=pval) & LOOSE$p & !LOOSE$P & LRP$LRT.pval.P>LRP$LRT.pval.D] = 2 
		#LRP$LRT.best.model[is.na(LRP$LRT.best.model) & (LRP$LRT.pval.nonest.m2>=pval|LRP$LRT.pval.nonest.m3>=pval) & LOOSE$d & !LOOSE$D & LRP$LRT.pval.D>LRP$LRT.pval.P] = 3
		#index = c(1:400)
		#LRP$true.model = ceiling(index/100)
		#ii = which(LRP$LRT.newmodel!=LRP$true.model)
		#jj = which(is.na(LRP$LRT.newmodel)==TRUE)
		
		
	}else{
		
		colnames(LRP) = paste('LRT.pval',c('p','P','d','D','A'),sep = '.')
		LRP$LRT.pval.C = pchisq(-2*log(LRP$LRT.pval.P * LRP$LRT.pval.D), df = 4, lower.tail = FALSE)
		LOOSE = LRP<pval; 
		colnames(LOOSE) =  c('p','P','d','D','A','C'); 
		LOOSE = data.frame(LOOSE, stringsAsFactors = FALSE)
		STRINGENT = LRP<pval.stringent;  
		colnames(STRINGENT) = c('p','P','d','D','A','C'); 
		STRINGENT = data.frame(STRINGENT, stringsAsFactors = FALSE)
		
	
		if(V1)
		{
			LRP$LRT.best.model = NA
			LRP$LRT.best.model[!LOOSE$P &!LOOSE$D] = 1
			LRP$LRT.best.model[!LOOSE$P & LOOSE$D] = 2
			LRP$LRT.best.model[LOOSE$P & !LOOSE$D] = 3
			LRP$LRT.best.model[LOOSE$P & LOOSE$D] = 4
		}else{
			#print('here')
			# chose models according to the pval of likelihood ratio
			LRP$LRT.best.model = 1
			LRP$LRT.best.model[(LOOSE$p|LOOSE$d)&!is.na(LOOSE$p)] =apply(cbind(LRP$LRT.pval.p[(LOOSE$p|LOOSE$d)&!is.na(LOOSE$p)], LRP$LRT.pval.d[(LOOSE$p|LOOSE$d)&!is.na(LOOSE$p)]),1, which.min)+1
			#LRP$LRT.best.model[LOOSE$P & LOOSE$D ] = 4
			LRP$LRT.best.model[LOOSE$P & LOOSE$D & LOOSE$A] = 4
			LRP$LRT.best.model[is.na(LRP$LRT.pval.p)] = NA
		
		}
		if(direct)
		{
			LRP$LRT.best.model = 1
			LRP$LRT.best.model[(LOOSE$p|LOOSE$d|LOOSE$A)&!is.na(LOOSE$p)] = apply(cbind(LRP$LRT.pval.p[(LOOSE$p|LOOSE$d|LOOSE$A)&!is.na(LOOSE$p)], LRP$LRT.pval.d[(LOOSE$p|LOOSE$d|LOOSE$A)&!is.na(LOOSE$p)],LRP$LRT.pval.A[(LOOSE$p|LOOSE$d|LOOSE$A)&!is.na(LOOSE$p)]),1, which.min)+1
			LRP$LRT.best.model[is.na(LRP$LRT.pval.p)] = NA
		}
		if(combined)
		{
			LRP$LRT.best.model = 1
			LRP$LRT.best.model[(LOOSE$p|LOOSE$d|LOOSE$C)&!is.na(LOOSE$p)] = apply(cbind(LRP$LRT.pval.p[(LOOSE$p|LOOSE$d|LOOSE$C)&!is.na(LOOSE$p)], LRP$LRT.pval.d[(LOOSE$p|LOOSE$d|LOOSE$C)&!is.na(LOOSE$p)],LRP$LRT.pval.C[(LOOSE$p|LOOSE$d|LOOSE$C)&!is.na(LOOSE$p)]),1, which.min)+1
			LRP$LRT.best.model[is.na(LRP$LRT.pval.p)] = NA
		}
		if(FDR)
		{
			pval.m2 = fdr2pval(LRP$LRT.pval.p[!is.na(LRP$LRT.pval.p)], fdr = fdr.value); MODEL2 = LRP$LRT.pval.p<=pval.m2
			pval.m3 = fdr2pval(LRP$LRT.pval.d[!is.na(LRP$LRT.pval.d)], fdr = fdr.value); MODEL3 = LRP$LRT.pval.d<=pval.m3
			pval.m4 = fdr2pval(LRP$LRT.pval.A[!is.na(LRP$LRT.pval.A)], fdr = fdr.value); MODEL4 = LRP$LRT.pval.A<=pval.m4
		
			fdr.m2 = fdrtool(LRP$LRT.pval.p[!is.na(LRP$LRT.pval.p)], statistic = 'pvalue')
			fdr.m3 = fdrtool(LRP$LRT.pval.d[!is.na(LRP$LRT.pval.d)], statistic = 'pvalue')
			fdr.m4 = fdrtool(LRP$LRT.pval.A[!is.na(LRP$LRT.pval.A)], statistic = 'pvalue')
			FDR.M2 = LRP$LRT.pval.p; FDR.M2[!is.na(LRP$LRT.pval.p)] = fdr.m2$qval
			FDR.M3 = LRP$LRT.pval.p; FDR.M3[!is.na(LRP$LRT.pval.p)] = fdr.m3$qval
			FDR.M4 = LRP$LRT.pval.p; FDR.M4[!is.na(LRP$LRT.pval.p)] = fdr.m4$qval
		
			MODEL2 = FDR.M2<= fdr.value
			MODEL3 = FDR.M3<= fdr.value
			MODEL4 = FDR.M4<= fdr.value

			LRP$LRT.best.model = 1
			LRP$LRT.best.model[(MODEL2| MODEL3 | MODEL4)&!is.na(MODEL2)] = apply(cbind(LRP$LRT.pval.p[(MODEL2| MODEL3 | MODEL4)&!is.na(MODEL2)], LRP$LRT.pval.d[(MODEL2| MODEL3 | MODEL4)&!is.na(MODEL2)],LRP$LRT.pval.A[(MODEL2| MODEL3 | MODEL4)&!is.na(MODEL2)]),1, which.min)+1
			LRP$LRT.best.model[is.na(LRP$LRT.pval.p)] = NA
		}	
	}
	ok = !is.na(LRP$LRT.best.model)
	
	LRP$LRT.quality = NA;
	LRP$LRT.quality[ok & (LRP$LRT.best.model == 1)] = -log10(1-apply(cbind(LRP$LRT.pval.p[ok & (LRP$LRT.best.model == 1)], LRP$LRT.pval.d[ok & (LRP$LRT.best.model == 1)]),1,min))
	LRP$LRT.quality[ok & (LRP$LRT.best.model == 2)] = -log10(LRP$LRT.pval.p[ok & (LRP$LRT.best.model == 2)])
	LRP$LRT.quality[ok & (LRP$LRT.best.model == 3)] = -log10(LRP$LRT.pval.d[ok & (LRP$LRT.best.model == 3)])
	LRP$LRT.quality[ok & (LRP$LRT.best.model == 4)] = -log10(LRP$LRT.pval.A[ok & (LRP$LRT.best.model == 4)])
	
	#	LRP$LRT.quality = NA;
	#	LRP$LRT.quality[LRP$LRT.best.model == 1] = apply(cbind(LRP$LRT.pval.p[LRP$LRT.best.model == 1], LRP$LRT.pval.d[LRP$LRT.best.model == 1]),1,min)
	#	LRP$LRT.quality[LRP$LRT.best.model == 2] = 1-LRP$LRT.pval.p[LRP$LRT.best.model == 2]
	#	LRP$LRT.quality[LRP$LRT.best.model == 3] = 1-LRP$LRT.pval.d[LRP$LRT.best.model == 3]
	#	LRP$LRT.quality[LRP$LRT.best.model == 4] = 1-apply(cbind(LRP$LRT.pval.P[LRP$LRT.best.model == 4], LRP$LRT.pval.D[LRP$LRT.best.model == 4]),1,max)

	return(LRP)
}

#### Akaike Information Criteria (AIC)
my.AIC = function(T = T, correction = FALSE, number.of.data.point = 48, diff.sigma=TRUE)
{
	set.nb.param();
	N = number.of.data.point;
	
	if(correction){corr = 1}else{corr = 0}
	
	## formula of AIC n*ln(error/n)+2*k+Constant == n*ln(error)+2*k+Constant and the sigma of noise is not used
	AIC.m1 = 2*n.param[1] + T$error.m1 + corr*2*n.param[1]*(n.param[1]+1)/(N-n.param[1]-1) # + N*log(2*pi*exp(1)*T$error.m1/N)
	AIC.m2 = 2*n.param[2] + T$error.m2 + corr*2*n.param[2]*(n.param[2]+1)/(N-n.param[2]-1) # + N*log(2*pi*exp(1)*T$error.m2/N)
	AIC.m3 = 2*n.param[3] + T$error.m3 + corr*2*n.param[3]*(n.param[3]+1)/(N-n.param[3]-1) # + N*log(2*pi*exp(1)*T$error.m3/N)
	AIC.m4 = 2*n.param[4] + T$error.m4 + corr*2*n.param[4]*(n.param[4]+1)/(N-n.param[4]-1) # + N*log(2*pi*exp(1)*T$error.m4/N) 
		
		
	AIC = data.frame(AIC.m1, AIC.m2, AIC.m3, AIC.m4, stringsAsFactors = FALSE)
	AIC$best.model = as.numeric(apply(AIC,1,which.min))
	
	if(correction){
		colnames(AIC) = c('AICc.m1', 'AICc.m2','AICc.m3','AICc.m4', 'AICc.best.model')
	}else{
		
		colnames(AIC)[5] = 'AIC.best.model';
	}
	
	return(AIC)								
}

my.AIC.prior = function(T = T, correction = FALSE, number.of.data.point = 48, fact = 2)
{
	set.nb.param();
	N = number.of.data.point;
	if(correction){corr = 1}else{corr = 0}
	AIC.m1 = 2*n.param[1] + N*log(T$error.m1) + corr*2*n.param[1]*(n.param[1]+1)/(N-n.param[1]-1) # + N*log(2*pi*exp(1)*T$error.m1/N)
	AIC.m2 = 2*n.param[2] + N*log(T$error.m2) + corr*2*n.param[2]*(n.param[2]+1)/(N-n.param[2]-1) # + N*log(2*pi*exp(1)*T$error.m2/N)
	AIC.m3 = 2*n.param[3] + N*log(T$error.m3) + corr*2*n.param[3]*(n.param[3]+1)/(N-n.param[3]-1) # + N*log(2*pi*exp(1)*T$error.m3/N)
	AIC.m4 = 2*n.param[4] + N*log(T$error.m4) + corr*2*n.param[4]*(n.param[4]+1)/(N-n.param[4]-1) # + N*log(2*pi*exp(1)*T$error.m4/N) 
	
	AIC = data.frame(AIC.m1, AIC.m2, AIC.m3, AIC.m4, stringsAsFactors = FALSE)
	AIC$best.model = as.numeric(apply(AIC,1,which.min))
	return(AIC)								
}

### BIC function
my.BIC = function(T = T, diff.sigma=TRUE)
{
	set.nb.param();
	
	## the formula of BIC used here is chi-square+k*ln(n)==error/sigma^2+k(ln(n)) in which sigma of noise is supported to be known.
	BIC.m1 = log(48)*n.param[1] + T$error.m1
	BIC.m2 = log(48)*n.param[2] + T$error.m2
	BIC.m3 = log(48)*n.param[3] + T$error.m3
	BIC.m4 = log(48)*n.param[4] + T$error.m4
	BIC = data.frame(BIC.m1, BIC.m2, BIC.m3, BIC.m4, stringsAsFactors = FALSE)
	
	BIC$BIC.best.model = apply(BIC,1,which.min)
	
	#colnames(BIC)[5] = 'BIC.best.model'
	
	return(BIC)								
}	

######### My model selction methods and specify only one sigma for exon and intron
my.Likelihood.Ratio.Test.Old = function(T = T, pval = 0.05,fact =0.25, pval.stringent = 10^(-15), V1 = FALSE, combined = FALSE, FDR = FALSE, fdr.value = 0.2, direct = FALSE, zt = seq(0,46,by = 2))
{
	#pval.stringent = 10^(-15); V1 = FALSE; combined = FALSE; FDR = FALSE; fdr.value = 0.2; direct = FALSE; zt = seq(0,46,by = 2);Mine=TRUE
	model.orders = list(prod = c(1,2,4), deg = c(1,3,4), all = c(1,4));
	set.nb.param();
	LRP = matrix(NA, nrow = nrow(T), ncol = length(unlist(model.orders)) - length(model.orders))
	i = 1
	sigma = fact
	for(P in 1:length(model.orders))
	{
		order = unlist(model.orders[P])
		for(m in 1:(length(order)-1))
		{
			m1 = order[m]; 
			m2 = order[m+1]
			eval(parse(text = paste('err1 = T$error.m',m1,sep = '')))
			eval(parse(text = paste('err2 = T$error.m',m2,sep = '')))
			
			## only one sigma considered for exon and intron
			LR = 1/sigma^2*(err1-err2);
			
			df = n.param[m2]-n.param[m1];
			pvals.LR = pchisq(LR,df = df ,lower.tail = FALSE);
			LRP[,i] = pvals.LR; 
			i = i+1 	
		}
		
	}
	LRP = data.frame(LRP, stringsAsFactors = FALSE); 
	
	Mine = FALSE
	if(Mine)
	{
		ZT.ex = grep('.rel.ampl.ex', colnames(T)); 
		ZT.int = grep('.rel.ampl.int', colnames(T));
		par.index.m2 = grep('.m2', colnames(T))
		par.index.m3 = grep('.m3', colnames(T))
		par.index.m2 = par.index.m2[-1]
		par.index.m3 = par.index.m3[-1]
		TT = as.matrix(T[,c(ZT.ex, ZT.int, par.index.m2, par.index.m3)])
		pvals.vuong = t(apply(TT, 1, vuong.test.nonest, parametrization ='sigmoid',sigma=sigma, method.intern = 'integration', debug = FALSE, sum.species = FALSE,  zt = seq(0,46,by = 2), absolute.signal = FALSE))
		#pvals.vuong = vuong.test.nonest(T,parametrization ='sigmoid',sigma=sigma, method.intern = 'simulation', debug = FALSE, sum.species = FALSE,  zt = seq(0,46,by = 2), absolute.signal = FALSE)
		pvals.vuong = data.frame(pvals.vuong, stringsAsFactors=FALSE)
		
		LRP$b = pvals.vuong[,1]
		LRP$B = pvals.vuong[,2]
		
		colnames(LRP) = paste('LRT.pval',c('p','P','d','D','A','nonest.m2','nonest.m3'),sep = '.')
		LRP$LRT.pval.C = pchisq(-2*log(LRP$LRT.pval.P * LRP$LRT.pval.D), df = 4, lower.tail = FALSE)
		LOOSE = LRP<pval; 
		colnames(LOOSE) =  c('p','P','d','D','A','nonest.m2','nonest.m3','C'); 
		LOOSE = data.frame(LOOSE, stringsAsFactors = FALSE)
		STRINGENT = LRP<pval.stringent;
		colnames(STRINGENT) = c('p','P','d','D','A','nonest.m2','nonest.m3','C'); 
		STRINGENT = data.frame(STRINGENT, stringsAsFactors = FALSE)
		# my method to select models
		LRP$LRT.best.model = NA
		#LRP$LRT.best.model = 1 # if it is not possible to select the model, we chose the simplest one 
		## selection with stringent criterion
		LRP$LRT.best.model[LOOSE$P & LOOSE$D & LOOSE$A] = 4
		#LRP$LRT.best.model[LOOSE$P & LOOSE$D ] = 4
		LRP$LRT.best.model[!LOOSE$p & !LOOSE$d & !LOOSE$A] = 1
		#LRP$LRT.newmodel[LOOSE$p & !LOOSE$P &(LRP$LRT.pval.nonest.m2<LRP$LRT.pval.nonest.m3)] = 2
		#LRP$LRT.newmodel[LOOSE$d & !LOOSE$D &(LRP$LRT.pval.nonest.m2>LRP$LRT.pval.nonest.m3)] = 3
		LRP$LRT.best.model[LOOSE$p & !LOOSE$P &(LRP$LRT.pval.nonest.m2<pval)] = 2
		LRP$LRT.best.model[LOOSE$d & !LOOSE$D &(pval>LRP$LRT.pval.nonest.m3)] = 3
		
	}else{
		
		colnames(LRP) = paste('LRT.pval',c('p','P','d','D','A'),sep = '.')
		LRP$LRT.pval.C = pchisq(-2*log(LRP$LRT.pval.P * LRP$LRT.pval.D), df = 4, lower.tail = FALSE)
		LOOSE = LRP<pval; 
		colnames(LOOSE) =  c('p','P','d','D','A','C'); 
		LOOSE = data.frame(LOOSE, stringsAsFactors = FALSE)
		STRINGENT = LRP<pval.stringent;  
		colnames(STRINGENT) = c('p','P','d','D','A','C'); 
		STRINGENT = data.frame(STRINGENT, stringsAsFactors = FALSE)
		
		
		if(V1)
		{
			LRP$LRT.best.model = NA
			LRP$LRT.best.model[!LOOSE$P &!LOOSE$D] = 1
			LRP$LRT.best.model[!LOOSE$P & LOOSE$D] = 2
			LRP$LRT.best.model[LOOSE$P & !LOOSE$D] = 3
			LRP$LRT.best.model[LOOSE$P & LOOSE$D] = 4
		}else{
			#print('here')
			# chose models according to the pval of likelihood ratio
			LRP$LRT.best.model = 1
			LRP$LRT.best.model[(LOOSE$p|LOOSE$d)&!is.na(LOOSE$p)] =apply(cbind(LRP$LRT.pval.p[(LOOSE$p|LOOSE$d)&!is.na(LOOSE$p)], LRP$LRT.pval.d[(LOOSE$p|LOOSE$d)&!is.na(LOOSE$p)]),1, which.min)+1
			#LRP$LRT.best.model[LOOSE$P & LOOSE$D ] = 4
			LRP$LRT.best.model[LOOSE$P & LOOSE$D & LOOSE$A] = 4
			LRP$LRT.best.model[is.na(LRP$LRT.pval.p)] = NA
			
		}
		if(direct)
		{
			LRP$LRT.best.model = 1
			LRP$LRT.best.model[(LOOSE$p|LOOSE$d|LOOSE$A)&!is.na(LOOSE$p)] = apply(cbind(LRP$LRT.pval.p[(LOOSE$p|LOOSE$d|LOOSE$A)&!is.na(LOOSE$p)], LRP$LRT.pval.d[(LOOSE$p|LOOSE$d|LOOSE$A)&!is.na(LOOSE$p)],LRP$LRT.pval.A[(LOOSE$p|LOOSE$d|LOOSE$A)&!is.na(LOOSE$p)]),1, which.min)+1
			LRP$LRT.best.model[is.na(LRP$LRT.pval.p)] = NA
		}
		if(combined)
		{
			LRP$LRT.best.model = 1
			LRP$LRT.best.model[(LOOSE$p|LOOSE$d|LOOSE$C)&!is.na(LOOSE$p)] = apply(cbind(LRP$LRT.pval.p[(LOOSE$p|LOOSE$d|LOOSE$C)&!is.na(LOOSE$p)], LRP$LRT.pval.d[(LOOSE$p|LOOSE$d|LOOSE$C)&!is.na(LOOSE$p)],LRP$LRT.pval.C[(LOOSE$p|LOOSE$d|LOOSE$C)&!is.na(LOOSE$p)]),1, which.min)+1
			LRP$LRT.best.model[is.na(LRP$LRT.pval.p)] = NA
		}
		if(FDR)
		{
			pval.m2 = fdr2pval(LRP$LRT.pval.p[!is.na(LRP$LRT.pval.p)], fdr = fdr.value); MODEL2 = LRP$LRT.pval.p<=pval.m2
			pval.m3 = fdr2pval(LRP$LRT.pval.d[!is.na(LRP$LRT.pval.d)], fdr = fdr.value); MODEL3 = LRP$LRT.pval.d<=pval.m3
			pval.m4 = fdr2pval(LRP$LRT.pval.A[!is.na(LRP$LRT.pval.A)], fdr = fdr.value); MODEL4 = LRP$LRT.pval.A<=pval.m4
			
			fdr.m2 = fdrtool(LRP$LRT.pval.p[!is.na(LRP$LRT.pval.p)], statistic = 'pvalue')
			fdr.m3 = fdrtool(LRP$LRT.pval.d[!is.na(LRP$LRT.pval.d)], statistic = 'pvalue')
			fdr.m4 = fdrtool(LRP$LRT.pval.A[!is.na(LRP$LRT.pval.A)], statistic = 'pvalue')
			FDR.M2 = LRP$LRT.pval.p; FDR.M2[!is.na(LRP$LRT.pval.p)] = fdr.m2$qval
			FDR.M3 = LRP$LRT.pval.p; FDR.M3[!is.na(LRP$LRT.pval.p)] = fdr.m3$qval
			FDR.M4 = LRP$LRT.pval.p; FDR.M4[!is.na(LRP$LRT.pval.p)] = fdr.m4$qval
			
			MODEL2 = FDR.M2<= fdr.value
			MODEL3 = FDR.M3<= fdr.value
			MODEL4 = FDR.M4<= fdr.value
			
			LRP$LRT.best.model = 1
			LRP$LRT.best.model[(MODEL2| MODEL3 | MODEL4)&!is.na(MODEL2)] = apply(cbind(LRP$LRT.pval.p[(MODEL2| MODEL3 | MODEL4)&!is.na(MODEL2)], LRP$LRT.pval.d[(MODEL2| MODEL3 | MODEL4)&!is.na(MODEL2)],LRP$LRT.pval.A[(MODEL2| MODEL3 | MODEL4)&!is.na(MODEL2)]),1, which.min)+1
			LRP$LRT.best.model[is.na(LRP$LRT.pval.p)] = NA
		}	
	}
	ok = !is.na(LRP$LRT.best.model)
	
	LRP$LRT.quality = NA;
	LRP$LRT.quality[ok & (LRP$LRT.best.model == 1)] = -log10(1-apply(cbind(LRP$LRT.pval.p[ok & (LRP$LRT.best.model == 1)], LRP$LRT.pval.d[ok & (LRP$LRT.best.model == 1)]),1,min))
	LRP$LRT.quality[ok & (LRP$LRT.best.model == 2)] = -log10(LRP$LRT.pval.p[ok & (LRP$LRT.best.model == 2)])
	LRP$LRT.quality[ok & (LRP$LRT.best.model == 3)] = -log10(LRP$LRT.pval.d[ok & (LRP$LRT.best.model == 3)])
	LRP$LRT.quality[ok & (LRP$LRT.best.model == 4)] = -log10(LRP$LRT.pval.A[ok & (LRP$LRT.best.model == 4)])
	
#	LRP$LRT.quality = NA;
#	LRP$LRT.quality[LRP$LRT.best.model == 1] = apply(cbind(LRP$LRT.pval.p[LRP$LRT.best.model == 1], LRP$LRT.pval.d[LRP$LRT.best.model == 1]),1,min)
#	LRP$LRT.quality[LRP$LRT.best.model == 2] = 1-LRP$LRT.pval.p[LRP$LRT.best.model == 2]
#	LRP$LRT.quality[LRP$LRT.best.model == 3] = 1-LRP$LRT.pval.d[LRP$LRT.best.model == 3]
#	LRP$LRT.quality[LRP$LRT.best.model == 4] = 1-apply(cbind(LRP$LRT.pval.P[LRP$LRT.best.model == 4], LRP$LRT.pval.D[LRP$LRT.best.model == 4]),1,max)
	
	colnames(LRP) = paste(colnames(LRP), '.Old', sep='')
	return(LRP)
}

my.AIC.Old = function(T = T, correction = FALSE, number.of.data.point = 48, fact = 2)
{
	set.nb.param();
	N = number.of.data.point;
	sigma = fact
	if(correction){corr = 1}else{corr = 0}
	
	## formula of AIC n*ln(error/n)+2*k+Constant == n*ln(error)+2*k+Constant and the sigma of noise is not used
	AIC.m1 = 2*n.param[1] + 1/sigma^2*T$error.m1 + corr*2*n.param[1]*(n.param[1]+1)/(N-n.param[1]-1) # + N*log(2*pi*exp(1)*T$error.m1/N)
	AIC.m2 = 2*n.param[2] + 1/sigma^2*T$error.m2 + corr*2*n.param[2]*(n.param[2]+1)/(N-n.param[2]-1) # + N*log(2*pi*exp(1)*T$error.m2/N)
	AIC.m3 = 2*n.param[3] + 1/sigma^2*T$error.m3 + corr*2*n.param[3]*(n.param[3]+1)/(N-n.param[3]-1) # + N*log(2*pi*exp(1)*T$error.m3/N)
	AIC.m4 = 2*n.param[4] + 1/sigma^2*T$error.m4 + corr*2*n.param[4]*(n.param[4]+1)/(N-n.param[4]-1) # + N*log(2*pi*exp(1)*T$error.m4/N) 
		
	AIC = data.frame(AIC.m1, AIC.m2, AIC.m3, AIC.m4, stringsAsFactors = FALSE)
	AIC$best.model = as.numeric(apply(AIC,1,which.min))
	
	if(correction){
		colnames(AIC) = c('AICc.m1', 'AICc.m2','AICc.m3','AICc.m4', 'AICc.best.model');
		colnames(AIC) = paste(colnames(AIC), '.Old', sep='');
		
	}else{
		colnames(AIC)[5] = 'AIC.best.model';
		colnames(AIC) = paste(colnames(AIC), '.Old', sep='');
	}
	
	return(AIC)								
}

my.BIC.Old = function(T = T, fact = 2)
{
	sigma =fact
	set.nb.param();
	## the formula of BIC used here is chi-square+k*ln(n)==error/sigma^2+k(ln(n)) in which sigma of noise is supported to be known.
	BIC.m1 = log(48)*n.param[1] + 1/sigma^2*T$error.m1
	BIC.m2 = log(48)*n.param[2] + 1/sigma^2*T$error.m2
	BIC.m3 = log(48)*n.param[3] + 1/sigma^2*T$error.m3
	BIC.m4 = log(48)*n.param[4] + 1/sigma^2*T$error.m4
	
	BIC = data.frame(BIC.m1, BIC.m2, BIC.m3, BIC.m4, stringsAsFactors = FALSE)
	BIC$best.model = apply(BIC,1,which.min)
	
	colnames(BIC)[5] = 'BIC.best.model'
	colnames(BIC) = paste(colnames(BIC), '.Old', sep='')
	
	return(BIC)								
}	


########## LAURA OLD codes for the model selection parts
Likelihood.Ratio.Test.LAURA = function(T = T, pval = 0.05, pval.stringent = 10^(-15), V1 = FALSE, combined = FALSE, FDR = FALSE, fdr.value = 0.2, direct = FALSE, zt = seq(0,46,by = 2))
{
	model.orders = list(prod = c(1,2,4), deg = c(1,3,4), all = c(1,4));
	set.nb.param();
	LRP = matrix(NA, nrow = nrow(T), ncol = length(unlist(model.orders)) - length(model.orders))
	i = 1
	
	for(P in 1:length(model.orders)){
		order = unlist(model.orders[P])
		for(m in 1:(length(order)-1)){
			m1 = order[m]; m2 = order[m+1]
			eval(parse(text = paste('err1 = T$error.m',m1,sep = '')))
			eval(parse(text = paste('err2 = T$error.m',m2,sep = '')))
			LR = length(zt)*(log(err1)-log(err2));#
			df = n.param[m2]-n.param[m1]
			pvals.LR = pchisq(LR,df = df ,lower.tail = FALSE);
			LRP[,i] = pvals.LR; i = i+1 	
			
		}
		
	}
	LRP = data.frame(LRP, stringsAsFactors = FALSE); colnames(LRP) = paste('LRT.pval',c('p','P','d','D','A'),sep = '.')
	LRP$LRT.pval.C = pchisq(-2*log(LRP$LRT.pval.P * LRP$LRT.pval.D), df = 4, lower.tail = FALSE)
	LOOSE = LRP<pval; colnames(LOOSE) =  c('p','P','d','D','A','C'); LOOSE = data.frame(LOOSE, stringsAsFactors = FALSE)
	STRINGENT = LRP<pval.stringent;  colnames(STRINGENT) = c('p','P','d','D','A','C'); STRINGENT = data.frame(STRINGENT, stringsAsFactors = FALSE)
	
	
	if(V1)
	{
		LRP$LRT.best.model = NA
		LRP$LRT.best.model[!LOOSE$P &!LOOSE$D] = 1
		LRP$LRT.best.model[!LOOSE$P & LOOSE$D] = 2
		LRP$LRT.best.model[LOOSE$P & !LOOSE$D] = 3
		LRP$LRT.best.model[LOOSE$P & LOOSE$D] = 4
	}else{
		LRP$LRT.best.model = 1
		LRP$LRT.best.model[(LOOSE$p|LOOSE$d)&!is.na(LOOSE$p)] =apply(cbind(LRP$LRT.pval.p[(LOOSE$p|LOOSE$d)&!is.na(LOOSE$p)], LRP$LRT.pval.d[(LOOSE$p|LOOSE$d)&!is.na(LOOSE$p)]),1, which.min)+1
		LRP$LRT.best.model[LOOSE$P & LOOSE$D] = 4
		LRP$LRT.best.model[is.na(LRP$LRT.pval.p)] = NA
	}
	if(direct){
		LRP$LRT.best.model = 1
		LRP$LRT.best.model[(LOOSE$p|LOOSE$d|LOOSE$A)&!is.na(LOOSE$p)] = apply(cbind(LRP$LRT.pval.p[(LOOSE$p|LOOSE$d|LOOSE$A)&!is.na(LOOSE$p)], LRP$LRT.pval.d[(LOOSE$p|LOOSE$d|LOOSE$A)&!is.na(LOOSE$p)],LRP$LRT.pval.A[(LOOSE$p|LOOSE$d|LOOSE$A)&!is.na(LOOSE$p)]),1, which.min)+1
		LRP$LRT.best.model[is.na(LRP$LRT.pval.p)] = NA
	}
	if(combined){
		LRP$LRT.best.model = 1
		LRP$LRT.best.model[(LOOSE$p|LOOSE$d|LOOSE$C)&!is.na(LOOSE$p)] = apply(cbind(LRP$LRT.pval.p[(LOOSE$p|LOOSE$d|LOOSE$C)&!is.na(LOOSE$p)], LRP$LRT.pval.d[(LOOSE$p|LOOSE$d|LOOSE$C)&!is.na(LOOSE$p)],LRP$LRT.pval.C[(LOOSE$p|LOOSE$d|LOOSE$C)&!is.na(LOOSE$p)]),1, which.min)+1
		LRP$LRT.best.model[is.na(LRP$LRT.pval.p)] = NA
	}
	if(FDR){
		pval.m2 = fdr2pval(LRP$LRT.pval.p[!is.na(LRP$LRT.pval.p)], fdr = fdr.value); MODEL2 = LRP$LRT.pval.p<=pval.m2
		pval.m3 = fdr2pval(LRP$LRT.pval.d[!is.na(LRP$LRT.pval.d)], fdr = fdr.value); MODEL3 = LRP$LRT.pval.d<=pval.m3
		pval.m4 = fdr2pval(LRP$LRT.pval.A[!is.na(LRP$LRT.pval.A)], fdr = fdr.value); MODEL4 = LRP$LRT.pval.A<=pval.m4
		
		fdr.m2 = fdrtool(LRP$LRT.pval.p[!is.na(LRP$LRT.pval.p)], statistic = 'pvalue')
		fdr.m3 = fdrtool(LRP$LRT.pval.d[!is.na(LRP$LRT.pval.d)], statistic = 'pvalue')
		fdr.m4 = fdrtool(LRP$LRT.pval.A[!is.na(LRP$LRT.pval.A)], statistic = 'pvalue')
		FDR.M2 = LRP$LRT.pval.p; FDR.M2[!is.na(LRP$LRT.pval.p)] = fdr.m2$qval
		FDR.M3 = LRP$LRT.pval.p; FDR.M3[!is.na(LRP$LRT.pval.p)] = fdr.m3$qval
		FDR.M4 = LRP$LRT.pval.p; FDR.M4[!is.na(LRP$LRT.pval.p)] = fdr.m4$qval
		
		MODEL2 = FDR.M2<= fdr.value
		MODEL3 = FDR.M3<= fdr.value
		MODEL4 = FDR.M4<= fdr.value
		
		LRP$LRT.best.model = 1
		LRP$LRT.best.model[(MODEL2| MODEL3 | MODEL4)&!is.na(MODEL2)] = apply(cbind(LRP$LRT.pval.p[(MODEL2| MODEL3 | MODEL4)&!is.na(MODEL2)], LRP$LRT.pval.d[(MODEL2| MODEL3 | MODEL4)&!is.na(MODEL2)],LRP$LRT.pval.A[(MODEL2| MODEL3 | MODEL4)&!is.na(MODEL2)]),1, which.min)+1
		LRP$LRT.best.model[is.na(LRP$LRT.pval.p)] = NA
	}	
	
	ok = !is.na(LRP$LRT.best.model)
	
	LRP$LRT.quality = NA;
	LRP$LRT.quality[ok & (LRP$LRT.best.model == 1)] = -log10(1-apply(cbind(LRP$LRT.pval.p[ok & (LRP$LRT.best.model == 1)], LRP$LRT.pval.d[ok & (LRP$LRT.best.model == 1)]),1,min))
	LRP$LRT.quality[ok & (LRP$LRT.best.model == 2)] = -log10(LRP$LRT.pval.p[ok & (LRP$LRT.best.model == 2)])
	LRP$LRT.quality[ok & (LRP$LRT.best.model == 3)] = -log10(LRP$LRT.pval.d[ok & (LRP$LRT.best.model == 3)])
	LRP$LRT.quality[ok & (LRP$LRT.best.model == 4)] = -log10(LRP$LRT.pval.A[ok & (LRP$LRT.best.model == 4)])
	
#	LRP$LRT.quality = NA;
#	LRP$LRT.quality[LRP$LRT.best.model == 1] = apply(cbind(LRP$LRT.pval.p[LRP$LRT.best.model == 1], LRP$LRT.pval.d[LRP$LRT.best.model == 1]),1,min)
#	LRP$LRT.quality[LRP$LRT.best.model == 2] = 1-LRP$LRT.pval.p[LRP$LRT.best.model == 2]
#	LRP$LRT.quality[LRP$LRT.best.model == 3] = 1-LRP$LRT.pval.d[LRP$LRT.best.model == 3]
#	LRP$LRT.quality[LRP$LRT.best.model == 4] = 1-apply(cbind(LRP$LRT.pval.P[LRP$LRT.best.model == 4], LRP$LRT.pval.D[LRP$LRT.best.model == 4]),1,max)
	
	
	colnames(LRP) = paste(colnames(LRP), '.Laura', sep='')
	return(LRP)
}

AIC.LAURA = function(T = T, correction = FALSE, number.of.data.point = 48, fact = 2)
{
	set.nb.param();
	N = number.of.data.point;
	if(correction){corr = 1}else{corr = 0}
	AIC.m1 = 2*n.param[1] + N*log(T$error.m1) + corr*2*n.param[1]*(n.param[1]+1)/(N-n.param[1]-1) # + N*log(2*pi*exp(1)*T$error.m1/N)
	AIC.m2 = 2*n.param[2] + N*log(T$error.m2) + corr*2*n.param[2]*(n.param[2]+1)/(N-n.param[2]-1) # + N*log(2*pi*exp(1)*T$error.m2/N)
	AIC.m3 = 2*n.param[3] + N*log(T$error.m3) + corr*2*n.param[3]*(n.param[3]+1)/(N-n.param[3]-1) # + N*log(2*pi*exp(1)*T$error.m3/N)
	AIC.m4 = 2*n.param[4] + N*log(T$error.m4) + corr*2*n.param[4]*(n.param[4]+1)/(N-n.param[4]-1) # + N*log(2*pi*exp(1)*T$error.m4/N) 
	
	AIC = data.frame(AIC.m1, AIC.m2, AIC.m3, AIC.m4, stringsAsFactors = FALSE)
	AIC$best.model = as.numeric(apply(AIC,1,which.min))
	
	if(correction){
		colnames(AIC) = c('AICc.m1', 'AICc.m2','AICc.m3','AICc.m4', 'AICc.best.model');
		colnames(AIC) = paste(colnames(AIC), '.Laura', sep='');
		
	}else{
		colnames(AIC)[5] = 'AIC.best.model';
		colnames(AIC) = paste(colnames(AIC), '.Laura', sep='');
	}	
	
	return(AIC)								
}

BIC.LAURA = function(T = T, fact = 2)
{
	set.nb.param();
	BIC.m1 = log(48)*n.param[1] + 48*log(T$error.m1)
	BIC.m2 = log(48)*n.param[2] + 48*log(T$error.m2)
	BIC.m3 = log(48)*n.param[3] + 48*log(T$error.m3)
	BIC.m4 = log(48)*n.param[4] + 48*log(T$error.m4)
	
	BIC = data.frame(BIC.m1, BIC.m2, BIC.m3, BIC.m4, stringsAsFactors = FALSE)
	BIC$best.model = apply(BIC,1,which.min)
	
	colnames(BIC)[5] = 'BIC.best.model'
	colnames(BIC) = paste(colnames(BIC), '.Laura', sep='')
	
	return(BIC)								
}	

################################################################################### finishing line for the model selection methods

## clean the results of model selection:

cleaning.model.selection.results = function(T=Tt, parametrization =c('cosine','sigmoid'),absolute.signal = FALSE, model=4)
{
	## T=Tt; parametrization ='sigmoid';absolute.signal = FALSE; model=4
	
	parametrization = parametrization[1]
	
	bounds = set.bounds(model = model, parametrization = parametrization, absolute.signal = absolute.signal); 
	upper = bounds$upper; 
	lower = bounds$lower;
	model = c(2,3,4)
		
	LRT.best.model = T$LRT.best.model
	AIC.best.model = T$AIC.best.model
	AICc.best.model =T$AICc.best.model
	
	cutoff = 0.00001
	
	for(n in 1:nrow(T))
	{
		cat(n, '.....\n')
		bool.2 =	abs(T$gamma.m2[n]-upper[1])<cutoff|abs(T$gamma.m2[n]-lower[1])<cutoff|
				abs(T$fold.change.int.m2[n]-upper[4])<cutoff|abs(T$fold.change.int.m2[n]-lower[4])<cutoff|
				abs(T$phase.int.m2[n]-upper[5])<cutoff|	abs(T$phase.int.m2[n]-lower[5])<cutoff|
				abs(T$up.time.int.m2[n]-upper[6])<cutoff|abs(T$up.time.int.m2[n]-lower[6])<cutoff|
				abs(T$down.time.int.m2[n]-lower[7])<cutoff|abs(T$down.time.int.m2[n]-upper[7])<cutoff
		
		bool.3 =	abs(T$gamma.m3[n]-upper[1])<cutoff|abs(T$gamma.m3[n]-lower[1])<cutoff|
				abs(T$eps.gamma.m3[n]-upper[2])<cutoff|abs(T$eps.gamma.m3[n]-lower[2])<cutoff|
				abs(T$phase.gamma.m3[n]-upper[3])<cutoff|abs(T$phase.gamma.m3[n]-lower[3])<cutoff	
		
		bool.4 =abs(T$gamma.m4[n]-upper[1])<cutoff|abs(T$gamma.m4[n]-lower[1])<cutoff|
				abs(T$eps.gamma.m4[n]-upper[2])<cutoff|abs(T$eps.gamma.m4[n]-lower[2])<cutoff|
				abs(T$phase.gamma.m4[n]-upper[3])<cutoff|abs(T$phase.gamma.m4[n]-lower[3])<cutoff|
				abs(T$fold.change.int.m4[n]-upper[4])<cutoff|abs(T$fold.change.int.m4[n]-lower[4])<cutoff|
				abs(T$phase.int.m4[n]-upper[5])<cutoff|	abs(T$phase.int.m4[n]-lower[5])<cutoff|
				abs(T$up.time.int.m4[n]-upper[6])<cutoff|abs(T$up.time.int.m4[n]-lower[6])<cutoff|abs(T$down.time.int.m4[n]-lower[7])<cutoff|abs(T$down.time.int.m4[n]-upper[7])<cutoff
		
		#parameters.model = c(c('gamma.m2', 'fold.change.int.m2', 'phase.int.m2', 'up.time.int.m2', 'down.time.int.m2'), c('gamma.m3', 'eps.gamma.m3', 'phase.gamma.m3'), 
		#					 c('gamma.m4', 'eps.gamma.m4', 'phase.gamma.m4', 'fold.change.int.m4', 'phase.int.m4', 'up.time.int.m4', 'down.time.int.m4'))
		#index = match(parameters.model, colnames(T))
		
		#bool.2;bool.3;bool.4
		#T[n,index];
		#upper;lower
		#LRT.best.model[n];AIC.best.model[n];AICc.best.model[n]
		
		if(bool.2)
		{
			if(!is.na(LRT.best.model[n])) 
			{
				if(LRT.best.model[n]==2) LRT.best.model[n]=NA
			}
			if(!is.na(AIC.best.model[n])) 
			{
				if(AIC.best.model[n]==2) AIC.best.model[n]=NA
			}
			if(!is.na(AICc.best.model[n]))
			{
				if(AICc.best.model[n]==2) AICc.best.model[n]=NA
			}
		}
		
		
		
		if(bool.3)
		{
			if(!is.na(LRT.best.model[n])) 
			{
				if(LRT.best.model[n]==3) LRT.best.model[n]=NA
			}
			if(!is.na(AIC.best.model[n])) 
			{
				if(AIC.best.model[n]==3) AIC.best.model[n]=NA
			}
			if(!is.na(AICc.best.model[n]))
			{
				if(AICc.best.model[n]==3) AICc.best.model[n]=NA
			}
			
		}
		
				
		if(bool.4)
		{
			if(!is.na(LRT.best.model[n])) 
			{
				if(LRT.best.model[n]==4) LRT.best.model[n]=NA
			}
			if(!is.na(AIC.best.model[n])) 
			{
				if(AIC.best.model[n]==4) AIC.best.model[n]=NA
			}
			if(!is.na(AICc.best.model[n]))
			{
				if(AICc.best.model[n]==4) AICc.best.model[n]=NA
			}
		}
		
	}
	T$LRT.best.model = LRT.best.model
	T$AIC.best.model = AIC.best.model
	T$AICc.best.model = AICc.best.model
	return(T)
	
	
}


############## NOISE ON THE DATA
noise.on.data = function(probetype = c('exon', 'intron'), method = c('time', 'replicates'), table = c('t','T'), log = 'FALSE')
{
	probetype = probetype[1]; 
	method = method[1];
	table = table[1]
	eval(parse(text = paste('X = ',table)))
	if(table == 't')
	{
		if(method == 'time'){i = 1:24}else{i = c(1,12)}
		if(probetype == 'exon'){ j = t$ex_in == 1} else{j = t$ex_in == -1}
		if(!log)
		{
			means = apply(X[j,i],1,mean); 
			var = apply(X[j,i],1,var)
		}else{
			means = apply(log(X[j,i]),1,mean); 
			var = apply(log(X[j,i]),1,var)
		}
		
	}else{
		if(method == 'time'){ if(probetype == 'exon'){ i = grep('.rel.ampl.ex',colnames(X))}else{i = grep('.rel.ampl.int',colnames(X))}}
		else{if(probetype == 'exon'){ i = match(c('ZT00.rel.ampl.ex','ZT26.rel.ampl.ex'), colnames(X))}else{i = match(c('ZT00.rel.ampl.int','ZT26.rel.ampl.int'), colnames(X))}}
		if(probetype == 'exon'){if(!log){means = X$exon.median}else{means = log(X$exon.median)}}else{if(!log){means = X$intron.median}else{means = log(X$intron.median)}}
		if(!log){var = means^2*apply(X[,i],1,var)}else{var = means^2*apply(log(X[,i]),1,var)}
	}
	
	
	par(mfrow = c(5,2))
	plot(means, var, pch = 16, cex = 0.3, col = rgb(0,0,0,0.5), main = paste(table,probetype, method, sep = '  -  '))
	plot(means, sqrt(var), pch = 16, cex = 0.3, col = rgb(0,0,0,0.5), main = paste(table,probetype, method, sep = '  -  '))
	plot(means, var, pch = 16, cex = 0.3, col = rgb(0,0,0,0.5), log = 'y', main = paste(table,probetype, method, sep = '  -  '))
	plot(means, sqrt(var), pch = 16, cex = 0.3, col = rgb(0,0,0,0.5), log = 'y', main = paste(table,probetype, method, sep = '  -  '))
	plot(means, var, pch = 16, cex = 0.3, col = rgb(0,0,0,0.5), log = 'x', main = paste(table,probetype, method, sep = '  -  '))
	plot(means, sqrt(var), pch = 16, cex = 0.3, col = rgb(0,0,0,0.5), log = 'x', main = paste(table,probetype, method, sep = '  -  '))
	plot(means, var, pch = 16, cex = 0.3, col = rgb(0,0,0,0.5), log = 'xy', main = paste(table,probetype, method, sep = '  -  '))
	plot(means, sqrt(var), pch = 16, cex = 0.3, col = rgb(0,0,0,0.5), log = 'xy', main = paste(table,probetype, method, sep = '  -  '))

	diff = log2(sqrt(var)) - log2(means);
	mm = mean(diff[diff!=-Inf] )
	plot(log2(means), log2(sqrt(var)), pch = 16, cex = 0.3, col = rgb(0,0,0,0.5), main = paste(table,probetype, method, sep = '  -  '))
	abline(a = 0, b = 1, col = 'gray')
	abline(a = mm, b = 1, col = 'green')
	
	plot(log2(means), diff , pch = 16, cex = 0.3, col = rgb(0,0,0,0.5), main = paste('a = ', round(mm, digits = 3),'||| sigma^2 = ', round(log(2^(2*mm)+1), digits = 3)))
	abline(h = mm, col = 'green')
	
}


plot.figure.fit.results = function(T = T,best.model = T$model.MCMC, gamma = T$gamma.MCMC , gamma.fc = T$eps.gamma.MCMC, phase.gamma = T$phase.gamma.MCMC, phase.int = T$phase.int.MCMC, P = P){

	layout(mat = matrix(c(1,1:9),nrow = 2, ncol = 5, byrow = TRUE), widths =c(0.65,0.55,1,1,1.05))
	#layout(mat = matrix(c(1:8),nrow = 2, ncol = 4, byrow = TRUE), widths =c(1.2,1,1,1.05))
	
	par(cex = 0.6, mgp = c(0,0.5,0),las = 1, tcl = -0.3, cex.main = 1, cex.axis = 0.9)
	
	par(mar = c(2,0,2,2))

	ylim = c(-0.1,1.15)

	h = hist(best.model,breaks = seq(0.5,4.5,by = 1), plot = FALSE)
	H = matrix(c(h$density, 0,h$counts[2:4]/sum(h$counts[2:4])), 4,2)
	space = 1
	b = barplot(H, beside = FALSE, border = FALSE, col = 1:4, space = space, axes = FALSE, xlim = c(-6.2,4), ylim = ylim)
	lines(x = b + c(space/2,-space/2), y = H[1,], col = 'gray', lty = '22')
	lines(x = b + c(space/2,-space/2), y = rep(1,2), col = 'gray', lty = '22')
	axis(4, las = 1)
	text(x = rep(0.7,4), y = c(0.25,0.7,0.85,0.96), c('Constitutive transcription\n& constant degradation', 'Rhythmic transcription', 'Rhythmic degradation','Rhythmic transcr. & degr.'), col = c('gray50',2:4), pos = 2, cex = 0.9, font = 2)
	lines(x = c(0.7 , b[1]-space/2), y =c(0.7, 0.9), col = 2)
	lines(x = c(0.7 , b[1]-space/2), y =c(0.85,0.98), col = 3)
	lines(x = c(0.7 , b[1]-space/2), y =c(0.97,1), col = 4)

	mtext('A',side =3, at = -6, font = 2, cex = 0.8, line = 0.8)



	fractions.ex = matrix(NA, nrow = 14, ncol = 4); 
	numbers.ex = rep(0,14)
	fractions.int = matrix(NA, nrow = 14, ncol = 4); 
	numbers.int = rep(0,14)
	for(r in 1:14){
		ok = (log2(T$exon.median)>(r-1))&(log2(T$exon.median)<=r)
		fractions.ex[r,] =  hist( best.model[ok], breaks = seq(0.5,4.5,by=1),plot = FALSE)$density
		numbers.ex[r] = sum(ok, na.rm = TRUE)
		ok = (log2(T$intron.median)>(r-1))&(log2(T$intron.median)<=r)
		fractions.int[r,] =  hist( best.model[ok], breaks = seq(0.5,4.5,by=1),plot = FALSE)$density
		numbers.int[r] = sum(ok, na.rm = TRUE)
		}
		
	par(mgp = c(0,0.5,0), mar = c(2,0,2,0))
	
	b = barplot(t(fractions.ex), col = c(1:4), border = c(1:4), ylim = ylim, space = 0.0, xlab = 'exonic probes expression level', axes = FALSE)
	#text(b,rep(1.005,14), numbers.ex, pos = 4, cex = 0.7, srt = 90, offset = 0.1)
	pos = c(b[1]-diff(b)[1],b)+diff(b)[1]/2
	text(pos[seq(1,15,by = 2)],rep(0,15)[seq(1,15,by = 2)], c(0:14)[seq(1,15,by = 2)], pos = 1)
	
	numbers.ex
	frac.ex = numbers.ex/max(numbers.ex)*0.13
	frac.ex = frac.ex+1.02
	polygon(rep(seq(b[1]-(b[2]-b[1])/2,b[length(b)]+(b[2]-b[1])/2,by = b[2]-b[1] ),each = 2), c(1.02,rep(frac.ex,each=2),1.02), border = NA, col = 'gray40')
	#text(0,1.08,'# of mRNA per\nexpr. level', col = 'gray40',pos = 4, cex = 0.8)	
	mtext('# of mRNA per\nexpr. level', line = -0.4 ,side = 3, at = 3,col = 'gray40', cex = 0.5 )
	mtext('B',side =3, at = -0.5, font = 2, cex = 0.8, line = 0.8)

	
	
	
	j2 = best.model == 2
	j3 = best.model == 3
	j4 = best.model == 4

	
	bounds = set.bounds(model = 4, parametrization = 'sigmoid')

	# half-lives 
	
	breaks = lseq(log(2)/max(gamma,na.rm = TRUE),log(2)/min(gamma,na.rm = TRUE), len = 31)
	
	
	par(mgp = c(1.6,0.5,0), mar = c(3,2,2,0))

	h2 = hist(log(2)/gamma[j2], breaks = breaks, plot = FALSE)
	h3 = hist(log(2)/gamma[j3], breaks = breaks, plot = FALSE)
	h4 = hist(log(2)/gamma[j4], breaks = breaks, plot = FALSE)
	
	plot(c(h2$breaks[1], h2$breaks, h2$breaks[length(h2$breaks)]),c(0, h2$counts, h2$counts[length(h2$counts)],0), log = 'x', type = 'n', main = 'Maximal half-lives',  xlab = 'half-lives [h]', ylab = '', ylim = range(h2$counts, h3$count, h4$counts))
	polygon(c(rep(h2$breaks, each = 2)), c(0, rep(h2$counts, each=2),0), col = 2+4, border = 2, lwd = 1)
	polygon(c(rep(h3$breaks, each = 2)), c(0, rep(h3$counts, each=2),0), col = 3+4, border = 3, lwd = 1)
	polygon(c(rep(h4$breaks, each = 2)), c(0, rep(h4$counts, each=2),0), col = 4+4, border = 4, lwd = 1)

	mtext('C',side =3, at = 0.12, font = 2, cex = 0.8, line = 0.8)


	
	# MINIMAL half-lives
	
	
	minimal.h.l = pmax(log(2)/gamma/gamma.fc,5/60)
	 
	breaks = lseq(min(minimal.h.l[j3|j4],na.rm = TRUE),max(minimal.h.l[j3|j4],na.rm = TRUE), len = 31)
	
	par(mgp = c(1.6,0.5,0), mar = c(3,2,2,0.5))

	h3 = hist(minimal.h.l[j3], breaks = breaks, plot = FALSE)
	h4 = hist(minimal.h.l[j4], breaks = breaks, plot = FALSE)

	plot(c(h3$breaks[1],h3$breaks,h3$breaks[length(h3$breaks)]),c(0,h3$counts,h3$counts[length(h3$counts)],0), log = 'x', type = 'n', main = 'Minimal half-lives',  xlab = 'half-lives [h]', ylab = '', ylim = range(h3$count, h4$counts))
	polygon(c(rep(h3$breaks, each = 2)), c(0, rep(h3$counts, each=2),0), col = 3+4, border = 3, lwd = 1)
	polygon(c(rep(h4$breaks, each = 2)), c(0, rep(h4$counts, each=2),0), col = 4+4, border = 4, lwd = 1)
	mtext('D',side =3, at = 0.05, font = 2, cex = 0.8, line = 0.8)


	
	# degradation fold changes
	
	par(mar = c(3,2,3,0.5))

	
	breaks = lseq(1,max(gamma.fc,na.rm = TRUE), len = 31)

	h3 = hist(gamma.fc[j3], breaks = breaks, plot = FALSE)
	h4 = hist(gamma.fc[j4], breaks = breaks, plot = FALSE)

	plot(c(h3$breaks[1],h3$breaks,h3$breaks[length(h3$breaks)]),c(0,h3$counts,h3$counts[length(h3$counts)],0), log = 'x', type = 'n', main = 'Degr. FC',  xlab = 'fold changes', ylab = '', ylim = range(0,h3$counts, h4$counts))
	polygon(c(rep(h3$breaks, each = 2)), c(0, rep(h3$counts, each=2),0), col = 3+4, border = 3, lwd = 1)
	polygon(c(rep(h4$breaks, each = 2)), c(0, rep(h4$counts, each=2),0), col = 4+4, border = 4, lwd = 1)
	
	mtext('E',side =3, at = 0.6, font = 2, cex = 0.8, line = 1.7)

	
	# fold.changes in exons  ## it does not work at all with Nacho results....!
	
	###
	ind = c(1:(ncol(P)-1))
	
	FC = apply(P[,ind],1,max)/apply(P[,ind],1,min); FC = pmax(FC,1)
	boxplot(FC~best.model,log = 'y', pch = 18, col = 5:8, border = 1:4, xlim = c(1.5,4.5), main = 'mRNA FC', axes = FALSE, names = c('CS-CD','RS-CD','CS-RD','RS-RD')) # , ylim = c(1,10)
	box(); axis(2)
	axis(1, at = 2, labels = 'RS-CD', col.axis =2, las = 2);
	axis(1, at = 3, labels = 'CS-RD', col.axis =3, las = 2);
	axis(1, at = 4, labels = 'RS-RD', col.axis =4, las = 2);
	mtext('F',side =3, at = 1, font = 2, cex = 0.8, line = 1.7)

	
	# degradation phases
	

	##par(mgp = c(1.6,0.5,0), mar = c(1,1,3,0.5))

	breaks = seq(bounds$upper[3],bounds$lower[3], len = 25)
	h3 = hist(phase.gamma[j3], breaks = breaks, plot = FALSE);
	h4 = hist(phase.gamma[j4], breaks = breaks, plot = FALSE)
	h3 = hist(phase.gamma[j3], breaks = breaks, col = 3+4, border = 3, xlab = '[h]', ylab = '', main = 'Phase of max. degr.', axes = FALSE, plot = TRUE, ylim = range(0,h3$counts, h4$counts))
	h4 = hist(phase.gamma[j4], breaks = breaks, col = 4+4, border = 4, add = TRUE, plot = TRUE)
	axis(1, at = seq(0,24,by = 6)); axis(2); box()
	
	# factor.rose.diag = 0.75
	# h3 = hist(phase.gamma[j3], breaks = breaks, plot = FALSE)
	# h4 = hist(phase.gamma[j4], breaks = breaks, plot = FALSE)
	
	# pg3 = circular(phase.gamma[j3], type = 'angles',units = 'hours',template = 'clock24')
	# pg4 = circular(phase.gamma[j4], type = 'angles',units = 'hours',template = 'clock24')
	# rose.diag(pg3, bins = 24, col = 3+4, border = 3, prop = 1/sqrt(max(h3$density))* factor.rose.diag, main =  'Phase of max. degr.', asp = 1)
	# points.circular(pg3, bins = 24, col = 3+4)
	# rose.diag(pg4, bins = 24, col = 4+4, border = 4, prop = 1/sqrt(max(h3$density))* factor.rose.diag, add = TRUE)
	# points.circular(pg4, bins = 24, col = 4+4)
	
	
	mtext('G',side =3, at = -1, font = 2, cex = 0.8, line = 1.7)

	

	par(mar = c(3,2,3,0.5))
	
	phase.ex = apply(P[,ind],1,which.max)/10
	diff.phase = (phase.ex-phase.gamma)%%24
	#diff.phase[which(diff.phase>12)] = diff.phase[which(diff.phase>12)]-24
	h = hist(diff.phase[j3], breaks = breaks, col = 3+4, border = 3, xlab = '[h]', ylab = '', main = 'Delay between\nmax. mRNA & degr.', axes = FALSE)
	h = hist(diff.phase[j4], breaks = breaks, col = 4+4, border = 4, add = TRUE)
	axis(1, at = seq(0,24,by = 6)); axis(2); box()
	
	# h3 = hist(diff.phase[j3], breaks = breaks,  plot = FALSE)
	# h4 = hist(diff.phase[j4], breaks = breaks, plot = FALSE)
	# diff3 = circular(diff.phase[j3], type = 'angles',units = 'hours',template = 'clock24')
	# diff4 = circular(diff.phase[j4], type = 'angles',units = 'hours',template = 'clock24')
	# rose.diag(diff3, bins = 24, col = 3+4, border = 3, prop = 1/sqrt(max(h3$density))* factor.rose.diag, main = 'Delay between\nmax. mRNA & degr.', asp = 1)
	# points.circular(diff3, bins = 24, col = 3+4)
	# rose.diag(diff4, bins = 24, col = 4+4, border = 4, prop = 1/sqrt(max(h3$density))* factor.rose.diag, add = TRUE)
	# points.circular(diff4, bins = 24, col = 4+4)
	
	mtext('H',side =3, at = -1, font = 2, cex = 0.8, line = 1.7)

	diff.phase = (phase.int-phase.gamma)%%24
	#diff.phase[which(diff.phase>12)] = diff.phase[which(diff.phase>12)]-24
	h = hist(diff.phase[j4], breaks = breaks, col = 4+4, border = 4, xlab = '[h]', ylab = '', main = 'Delay between\nmax. pre-mRNA & degr.', axes = FALSE)
	axis(1, at = seq(0,24,by = 6)); axis(2); box()
	mtext('I',side =3, at = -1, font = 2, cex = 0.8, line = 1.7)

}

























































#do.lower.and.upper.bound = function(n.comp = 2, changing = FALSE){
#	gamma.max = log(2)/10*60 ; gamma.min = log(2)/24
#	if(changing){eps.gamma.max = 1; eps.gamma.min = 0; phase.gamma.max = 24; phase.gamma.min = 0}
#	eps.24.S.max = 1; eps.24.S.min = 0
#	phase.24.S.max = 24; phase.24.S.min = 0; 
#	upper = c(gamma.max , eps.24.S.max, phase.24.S.max);
#	lower = c(gamma.min , eps.24.S.min, phase.24.S.min);
#	if(n.comp>1){
#		eps.12.S.max = 1; eps.12.S.min = 0
#		phase.12.S.max = 12; phase.12.S.min = 0; 
#		upper = c(upper,eps.12.S.max,phase.12.S.max)
#		lower = c(lower,eps.12.S.min,phase.12.S.min)
#	}
#	if(!changing){lower  <<- lower; upper <<- upper}else{
#		lower <<- c(gamma.min , eps.gamma.min, phase.gamma.min,  eps.24.S.min, phase.24.S.min,eps.12.S.min,phase.12.S.min);
#		upper <<- c(gamma.max , eps.gamma.max, phase.gamma.max,  eps.24.S.max, phase.24.S.max,eps.12.S.max,phase.12.S.max);}
#}
#
#
#
#
#
#
#
#do.fit.constant.degradation = function(index = which(T$gene == 'Dbp'), n.comp = 2){
#	gene = T$gene[index]; cat(gene, '\n')
#	M = unlist(T[index, ZT.ex]); S = unlist(T[index,ZT.int])
#
#	do.parameter.initialization(n.comp = n.comp, index = index)
#	do.lower.and.upper.bound(n.comp = n.comp)
#
#	opt = optim(par.init, f2min.constant, var = rbind(M,S), method = 'L-BFGS-B', lower = lower, upper = upper)
#	res.fit = list(gamma = opt$par[1], error = opt$value, rel.ampl.int = opt$par[2], phase.int = opt$par[3])
#	if(n.comp>1){res.fit$rel.ampl.12.int = opt$par[4]; res.fit$phase.12.int = opt$par[5]}
#	return(res.fit)
#	}
#
#plot.fit.constant.degradation = function(index = which(T$gene == 'Dbp'), n.comp = 2){
#	gene = T$gene[index]; cat(gene, '\n')
#	M = unlist(T[index, ZT.ex]); S = unlist(T[index,ZT.int])
#
#	plot(1,1,type = 'n', xlim = range(zt), ylim = range(M,S,0,2), xlab = 'time', ylab = '', main = paste(gene, '- n.comp :', n.comp))
#	points(zt,M, type = 'b', col = 'steelblue', lwd = 2, pch = 16)
#	points(zt,S, type = 'b', col = 'green3', lwd = 2, pch = 5)
#	
#	if(n.comp == 1){err = T$err.fit.1.comp[index]; gamma = T$DEG.const.deg.fit.1.comp[index]; eps.24.S = T$rel.ampl.int.fit.1.comp[index]; phase.24.S = T$phase.int.fit.1.comp[index]; eps.12.S =0 ; phase.12.S = 0
#		}else{err = T$err.fit.2.comp[index]; gamma = T$DEG.const.deg.fit.2.comp[index]; eps.24.S = T$rel.ampl.int.fit.2.comp[index]; phase.24.S = T$phase.int.fit.2.comp[index]; eps.12.S = T$rel.ampl.12.int.fit.2.comp[index] ; phase.12.S = T$phase.12.int.fit.2.comp[index]}
#	VAR = var(M)+var(S); perc = (VAR-err/24)/VAR
#	t = seq(min(zt), max(zt), by = 0.1)
#	points(t, compute.m(t = t,gamma = gamma, eps.24.S = eps.24.S, phase.24.S = phase.24.S, eps.12.S = eps.12.S, phase.12.S = phase.12.S), type = 'l', col = 'steelblue')
#	points(t, compute.s(t = t, eps.24.S = eps.24.S, phase.24.S = phase.24.S, eps.12.S = eps.12.S, phase.12.S = phase.12.S), type = 'l', col = 'green3')
#	legend('topright', legend = c(paste('half-life =',round(log(2)/gamma*60),'min =~', round(log(2)/gamma, digits = 2),'h'), paste('err =',round(err,digits = 3),'  %var explained =',round(perc, digits=2))) , bty = 'n', lty = 0)
#	}	
#
#
#
#
#f2min.constant = function(par, var){
#	M = var[1,]; S = var[2,]
#	n.comp = floor(length(par)/2)
#	gamma = par[1]
#	eps.24.S = par[2]; phase.24.S = par[3]; eps.12.S = 0; phase.12.S = 0
#	if(n.comp>1){eps.12.S = par[4]; phase.12.S = par[5]}
#
#	s = compute.s(t = zt, eps.24.S = eps.24.S, phase.24.S = phase.24.S, eps.12.S = eps.12.S, phase.12.S = phase.12.S)
#	m = compute.m(t = zt, gamma, eps.24.S = eps.24.S, phase.24.S = phase.24.S, eps.12.S = eps.12.S, phase.12.S = phase.12.S)
#		
#	error = sum((S-s)^2)+sum((M-m)^2)
#	return(error)
#	}
#
####### changing half-lives
#
#do.fit.changing.degradation = function(index = which(T$gene == 'Dbp')){
#	gene = T$gene[index]; cat(gene, '\n')
#	M = unlist(T[index, ZT.ex]); S = unlist(T[index,ZT.int])
#
#	do.lower.and.upper.bound(changing = TRUE)
#	Nfit = 5; errors = rep(NA, Nfit)
#	
#	for(fit.number in 1:Nfit){
#		do.parameter.initialization.changing(index = index, k = fit.number, Nfit = Nfit)
#		opt = optim(par.init, f2min.changing, var = rbind(M,S), method = 'L-BFGS-B', lower = lower, upper = upper)
#		res.fit = list(gamma = opt$par[1], eps.gamma = opt$par[2], phase.gamma = opt$par[3],  error = opt$value, rel.ampl.int = opt$par[4], phase.int = opt$par[5], rel.ampl.12.int = opt$par[6], phase.12.int = opt$par[7] )
#		errors[fit.number] = opt$value
#		eval(parse(text = paste('res.fit.', fit.number, ' = res.fit', sep = '')))
#	}
#	imin = which.min(errors); eval(parse(text = paste('res.fit = res.fit.', imin, sep = '')))
#	return(res.fit)
#	}
#
#plot.fit.changing.degradation = function(index = which(T$gene == 'Dbp')){
#	gene = T$gene[index]; cat(gene, '\n')
#	M = unlist(T[index, ZT.ex]); S = unlist(T[index,ZT.int])
#
#	plot(1,1,type = 'n', xlim = range(zt), ylim = range(M,S,0,2), xlab = 'time', ylab = '', main = paste(gene, '- changing half-life'))
#	points(zt,M, type = 'b', col = 'steelblue', lwd = 2, pch = 16)
#	points(zt,S, type = 'b', col = 'green3', lwd = 2, pch = 5)
#	
#	err = T$err.fit[index]; gamma = T$DEG.chang.deg.fit[index]; eps.gamma = T$DEG.REL.AMPL.chang.deg.fit[index]; phase.gamma = T$DEG.PHASE.chang.deg.fit[index]; eps.24.S = T$rel.ampl.int.fit[index]; phase.24.S = T$phase.int.fit[index]; eps.12.S = T$rel.ampl.12.int.fit[index] ; phase.12.S = T$phase.12.int.fit[index]
#	VAR = var(M)+var(S); perc = (VAR-err/24)/VAR
#	t = seq(min(zt), max(zt), by = 0.1)
#	points(t, compute.m.changing(t = t,gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma,  eps.24.S = eps.24.S, phase.24.S = phase.24.S, eps.12.S = eps.12.S, phase.12.S = phase.12.S), type = 'l', col = 'steelblue')
#	points(t, compute.s(t = t, eps.24.S = eps.24.S, phase.24.S = phase.24.S, eps.12.S = eps.12.S, phase.12.S = phase.12.S), type = 'l', col = 'green3')
#	points(t, compute.s(t = t, eps.24.S = eps.gamma, phase.24.S = phase.gamma), type = 'l', col = 'orangered')
#	legend('topright', legend = c(paste('half-life =',round(log(2)/gamma*60),'min =~', round(log(2)/gamma, digits = 2),'h'), paste('degradation rel.ampl =', round(eps.gamma, digits =2)), paste('degradation phase =', round(phase.gamma,digits= 1)), paste('err =',round(err,digits = 3),'  %var explained =',round(perc, digits=2))) , bty = 'n', lty = 0)
#	}	
#
#
#
#	
#do.parameter.initialization.changing = function(index = index, k = 1, Nfit = 5){
#	phases.gamma.init = seq(0,24,len = Nfit)
#	
#	gamma.init = log(2)/120*60 # initial half-life = 120 min
#	eps.gamma.init = 0.2 ; phase.gamma.init = phases.gamma.init[k];
#	eps.24.S.init = T$rel.ampl.int[index] ; phase.24.S.init = T$phase.int[index]
#	eps.12.S.init = T$rel.ampl.12.int[index] ; phase.12.S.init = T$phase.12.int[index]
#	par.init <<- c(gamma.init, eps.gamma.init, phase.gamma.init, eps.24.S.init, phase.24.S.init, eps.12.S.init,phase.12.S.init)
#}
#
#
#





#set.fake.values = function(i = 1, sd.noise = 0.1){
#	T$gene[i] <<- paste('fake',i, sep = '')
#	T$rel.ampl.int[i] <<- sample(seq(0,1,by = 0.01),1)
#	T$rel.ampl.12.int[i] <<- sample(seq(0,1-T$rel.ampl.int[i],by = 0.01),1)
#	T$phase.int[i] <<- sample(seq(0,24,by = 0.1),1)
#	T$phase.12.int[i] <<- sample(seq(0,12,by = 0.1),1)
#	T$gamma[i] <<- sample(seq(log(2)/24,log(2)/10*60, by = 0.1),1)
#	T$eps.gamma[i] <<- sample(seq(0,1,by = 0.01),1)
#	T$phase.gamma[i] <<- sample(seq(0,24,by = 0.1),1)
#		
#	T[i, ZT.int] <<- compute.s(t = zt, eps.24.S = T$rel.ampl.int[i] , phase.24.S = T$phase.int[i], eps.12.S = T$rel.ampl.12.int[i], phase.12.S = T$phase.12.int[i])
#	T[i, ZT.ex] <<- compute.m.changing(t = zt, gamma = T$gamma[i] , eps.gamma = T$eps.gamma[i], phase.gamma = T$phase.gamma[i],eps.24.S = T$rel.ampl.int[i] , phase.24.S = T$phase.int[i], eps.12.S = T$rel.ampl.12.int[i], phase.12.S = T$phase.12.int[i])
#		
#	noise.int = rnorm(length(zt),mean = 0, sd = sd.noise)
#	noise.ex = rnorm(length(zt),mean = 0, sd = sd.noise)
#	T[i, ZT.int] <<- T[i, ZT.int]+ noise.int
#	T[i, ZT.ex] <<- T[i, ZT.ex] + noise.ex
#	}



#times = seq(0,2*pi,by = 0.1)
#integrale = 0*times
#for(i in 1:length(times)){
#	integrale[i] = integrate(sin,lower = times[1], upper = times[i])$value
#	}
#
#plot(integrale, type = 'l')
#points(-cos(times), type = 'l', col = 'red')
#
#
#
#cosine.fit.24h = function(par = c(1,0), sup = sup){
#	eps = par[1]; psi.24 = par[2];
#	cosine = 1 + eps*cos(2*pi/24*sup - psi.24) 	
#	return(cosine)
#	}
#	
#cosine.fit.12h = function(par = c(1,0,1,0), sup = sup){
#	eps = par[1]; psi.24 = par[2]; nu = par[3]; psi.12 = par[4];
#	cosine = cosine.fit.24h(par = par[1:2], sup = sup) + nu * cos(2*2*pi/24*sup - psi.12)
#	return(cosine)
#	}
#	
#	
#cosine.fit.8h = function(par = c(1,0,1,0,1,0), sup = sup){
#	eps = par[1]; psi.24 = par[2]; nu = par[3]; psi.12 = par[4]; mu = par[5]; psi.8 = par[6]
#	cosine = cosine.fit.12h(par = par[1:4],sup = sup) + mu * cos(3*2*pi/24*sup - psi.8)
#	return(cosine)
#	}	
#
#
#
#
#f2min.for.cosine.12h = function(par = c(1,0,1,0),var = rep(1,12)){
#
#	M = unlist(var);# cat(M,'\n')
#	m.est = cosine.fit.12h(par = par, sup = sup)
#	err.m = sum((m.est-M)^2)
#	
#	alpha = 50
#	constrain.positive = min(sum(exp(-alpha*m.est)),10000000)
#	
#	err = err.m + constrain.positive
#	return(err)
#	}	
#	
#f2min.for.cosine.8h = function(par = c(1,0,1,0,1,0),var = rep(1,12)){
#
#	S = unlist(var);# cat(M,'\n')
#	s.est = cosine.fit.8h(par = par, sup = sup)
#	err.s = sum((s.est-S)^2)
#	
#	alpha = 50
#	constrain.positive = min(sum(exp(-alpha*s.est)),10000000)
#	
#	err = err.s + constrain.positive
#	return(err)
#	}	
#
#
#
#f2min = function(par = c(1,1,0),var = rep(0,25), add.penalty = FALSE){	
#	fft.m.est = unlist(var[1:12]);
#	S = unlist(var[13:24])
#	n = var[25]
#	
#	gamma = par[1]; eps_g = par[2]; psi_g. = par[3]
#	hat.s = estimate.s(gamma = gamma , eps_g = eps_g , psi_g. = psi_g. , fft.m.est = fft.m.est)
#	err.s = sum((hat.s-S)^2)
#	
#	alpha = 50
#	constrain.s.positive = min(sum(exp(-alpha*hat.s)),10000000)
#	#cat(constrain.s.positive,'\n')
#	
#	
#	err = err.s + constrain.s.positive
#	if(add.penalty){err = err+ (eps_g)^2}
#	return(Re(err))
#	}	
#
#	
#	
#	
#estimate.s = function(gamma = 1 , eps_g  = 1, psi_g. = 0 , fft.m.est  = rep(0,12), normalize.s = TRUE){
#	
#	N = length(fft.m.est); mod.fft.m.est = fft.m.est[c((N/2+2):N,1:(N/2))]
#	Nmod = length(mod.fft.m.est)
#	X = compute.X(gamma = gamma, eps_g = eps_g, psi_g. = psi_g., N = Nmod)
#	
#	#cat('dim X : ', dim(X), '\n')
#	#cat('length(mod.fft.m.est) : ',length(mod.fft.m.est),'\n' )
#	
#	fft.s = X %*% mod.fft.m.est
#	
#	fft.s = c(fft.s[ceiling(Nmod/2):Nmod],0,fft.s[1:floor(Nmod/2)])
#
#	s = Re( fft(fft.s, inverse = TRUE))/length(fft.s); 
#	if(normalize.s){s = s/mean(s)}
#	return(s)
#}	
#	
#	
#compute.X = function(gamma = gamma, eps_g = eps_g, psi_g. = psi_g., N = 11){
#	n = floor(N/2)
#	diag.vec = gamma + c(-n:n)*1i*2*pi/24; nX = length(diag.vec)
#	A = 0.5*gamma*eps_g*exp(1i*psi_g.)
#	A. = 0.5*gamma*eps_g*exp(-1i*psi_g.)
#	X = diag(diag.vec) + A * diag(1,nX+1,nX+1)[2:(nX+1),1:nX] + A.* diag(1,nX+1,nX+1)[1:nX,2:(nX+1)]
#	return(X)
#	}	
#	
#	
#make.fit = function(rep = 1, n = 2, i = 1, plot = TRUE, add.penalty = FALSE, constant.degradation = FALSE){
#	if(rep == 1){T = T1}else{T = T2}
#	
#	
#	cat(T$gene[i], '\n'); 
#	########### 
#	#Estimation of m(t)
#	########### 	
#
#	M = unlist(T[i,ZT.ex])
#	S = unlist(T[i, ZT.int])
#	
#	if(plot){
#		title = paste('Replicate', rep ,' - ',T$gene[i],'\n M - pval : ',T$pval.ex[i], ' ||  S - pval : ', T$pval.int[i])
#		ylim = range(0,S,M)
#		plot(1,1, type = 'n', xlim = c(0,24), ylim = ylim , axes = FALSE, xlab = 'time [h of a day]', ylab = '', main = title)
#		abline(h = 1, col = 'gray')
#		points(sup, M, pch = 16, col = 'steelblue3', type = 'b')
#		points(sup, S, pch = 5, col = 'green3', type = 'b')
#		axis(1,at = sup); axis(2)
#	} 
#	
#	if(n == 3){	
#		par.i = c(min(1,T$rel.ampl.ex[i]),T$phase.ex[i]*2*pi/24, 0,0)
#		out.optim = optim(par.i, f2min.for.cosine.12h, var = M, method = 'L-BFGS-B', lower = c(0,-2*pi,0,-2*pi), upper = c(1,2*pi,1,2*pi))
#		par.f = out.optim$par; eps.m = par.f[1];psi.m = par.f[2]; nu.m = par.f[3];  psi.12.m = par.f[4]
#		m.t.i = cosine.fit.12h(par = par.i, sup = sup); m.t.i.p = cosine.fit.12h(par = par.i, sup = sup.p)
#		m.t.f = cosine.fit.12h(par = par.f, sup = sup); m.t = m.t.f; m.t.p = cosine.fit.12h(par = par.f,sup = sup.p)
#		fft.m.est = fft(m.t)
#	}else{
#		eps.m = min(1,T$rel.ampl.ex[i]); psi.m = T$phase.ex[i]*2*pi/24; nu.m = 0; psi.12.m =0;
#		m.t = cosine.fit.24h(par = c(eps.m, psi.m), sup = sup); m.t.p = cosine.fit.24h(par = c(eps.m, psi.m), sup = sup.p);
#		fft.m.est = fft(m.t); m.t.i = m.t; m.t.i.p = m.t.p
#	}	
#	
#	if(plot){
#		title = paste('Replicate', rep ,' - ',T$gene[i], '\n eps.m : ' , round(eps.m, digits = 2),' ||  psi.m : ', round(((psi.m)%%(2*pi))*24/2/pi, digits = 1)  ,'h ||  nu.m : ', round(nu.m, digits = 2),' || psi.12.m : ',round(((psi.12.m)%%pi)*24/2/pi, digits = 1))
#		plot(sup,M, type = 'b',pch = 16, xlab = 'time', ylab = 'm(t)', main = title, ylim = range(0,M,m.t.i.p,m.t.p), axes = FALSE, col = 'steelblue'); 
#		abline(h = 1, col = 'gray')
#		points(sup.p, m.t.i.p, type = 'l', col = 'gray', lwd =  1) ; #points(sup, m.t.i, type = 'l', col = 'gray');
#		points(sup.p, m.t.p, type = 'l', col = 'steelblue', lwd =  3); #points(sup,m.t, type = 'l',col = 'steelblue'); 
#		axis(1,at = sup); axis(2)
#	}
#	
#	error.m = sum((M-m.t)^2)
#
#	gamma.min = log(2)/100; gamma.max = log(2)/(20/60); 
#	eps_g.min = 0 ; 		eps_g.max = 1 
#	psi_g.min = -pi; 		psi_g.max = pi 
#	if(constant.degradation){eps_g.max = .Machine$double.eps; psi_g.max = -pi + 10*.Machine$double.eps }
#
#	var = c(fft.m.est,S,n)
#
#	# initial parameters for g(t)
#	gammas = c(log(2)/1, log(2)/10)
#	eps_gs = c(0,0.5,1)
#	psi_g.s = c( -pi,-pi/2, 0 ,pi/2, pi)
#	
#	ERROR = 10^10; errors = 10^10; N.PAR = c(0,0,0); n.par = rep(0,3)
#	for(gamma in gammas){
#		#cat('\t gamma : ',gamma, '\n')
#		for(eps_g in eps_gs){
#			#cat('\t\t eps_g : ', eps_g, '\n')
#			for(psi_g. in psi_g.s){
#				#cat('\t\t\t psi_g. : ', if(psi_g. == 0){'0.00'}else{round(psi_g.,digits = 2)}, '\t\t')
#				
#				par =  c(gamma, eps_g, psi_g.)
#				res = optim(par = par, f2min , gr = NULL, var = var , add.penalty = add.penalty, method = "L-BFGS-B", lower = c(gamma.min, eps_g.min, 	psi_g.min), upper = c(gamma.max, eps_g.max, psi_g.max) )
#				error = res$value; tmp.par = res$par; 
#				#cat('error = ',error,'\n')
#
#				if(error<ERROR){ERROR = error; n.par = tmp.par}
#				errors = c(errors, error); N.PAR = rbind(N.PAR, tmp.par)
#				
#			}	
#		}
#	}
#	
#
#
#	hl = log(2)/n.par[1]; hours = floor(hl); minutes = round((hl-hours)*60); if(hours > 0){halflife = paste(hours, 'h', minutes, 'min')}else{halflife = paste(minutes, 'min')}
#	phase.in.hours =  round(((n.par[3])%%(2*pi))*24/2/pi, digits = 1)
#	
#	text = paste('\t Best fit: \n\t\t half-life : ' , halflife,'\n\t\t eps_g : ', round(n.par[2], digits = 2) ,'\n\t\t phase psi_g. : ', phase.in.hours,'h \n\n\n' );cat(text) 
#
#	title = paste('Replicate', rep ,' - ',T$gene[i], '\n half-life : ' , halflife,' ||  eps_g : ', round(n.par[2], digits = 2),' ||  phase psi_g. : ', phase.in.hours,'h')
#
#	s.t = estimate.s(gamma = n.par[1],eps_g = n.par[2], psi_g.  = n.par[3], fft.m.est = fft.m.est);
#	fft.s.est = fft(s.t)/length(s.t); eps.s = 2*abs(fft.s.est[2]); phase.s.24 = min(4*pi,phase.n(s.t, n = 1),na.rm = TRUE); nu.s = 2*abs(fft.s.est[3]); phase.s.12 = min(4*pi,phase.n(s.t, n = 2),na.rm = TRUE); mu.s = 2*abs(fft.s.est[4]); phase.s.8 = min(4*pi,phase.n(s.t, n = 3),na.rm = TRUE);
#
#	if(n == 3){ s.t.p = cosine.fit.8h(par = c(eps.s,phase.s.24, nu.s, phase.s.12, mu.s, phase.s.8),sup = sup.p)
#	}else{s.t.p = cosine.fit.12h(par = c(eps.s,phase.s.24,nu.s,phase.s.12), sup = sup.p)}
#
#
#	if(n == 3){
#		par.i = c(T$rel.ampl.int[i],T$phase.int[i]*2*pi/24,T$rel.ampl.12.int[i],T$phase.12.int[i]*2*pi/12,0,0)
#		out.optim = optim(par.i, f2min.for.cosine.8h, var = S, method = 'L-BFGS-B', lower = c(0,-2*pi,0,-2*pi,0,-2*pi), upper = c(1,2*pi,1,2*pi,1,2*pi))
#		par.s.f = out.optim$par; 
#		best.fit.on.s = cosine.fit.8h(par = par.s.f, sup = sup.p); err.best = out.optim$value
#	}else{
#		par.i = c(T$rel.ampl.int[i],T$phase.int[i]*2*pi/24,T$rel.ampl.12.int[i],T$phase.12.int[i]*2*pi/12)
#		out.optim = optim(par.i, f2min.for.cosine.12h, var = S, method = 'L-BFGS-B', lower = c(0,-2*pi,0,-2*pi), upper = c(1,2*pi,1,2*pi))
#		par.s.f = out.optim$par; 
#		best.fit.on.s = cosine.fit.12h(par = par.s.f, sup = sup.p);err.best = out.optim$value
#	}
#	
#	degradation = cosine.fit.24h(par = c(n.par[2],n.par[3]), sup = sup.p)
#
#	if(plot){
#		ylim = range(0,s.t,m.t,S,M)
#		plot(1,1, type = 'n', xlim = c(0,24), ylim = ylim , axes = FALSE, xlab = 'time [h of a day]', ylab = '', main = title)
#		abline(h = 1, col = 'gray')
#		points(sup.p,degradation, col = 'red', type = 'l'); if(!constant.degradation){abline(v = phase.in.hours, col = 'red', lwd = 0.2)}
#		points(sup, M, pch = 16, col = 'steelblue3', type = 'b')
#		points(sup, S, pch = 5, col = 'green3', type = 'b')
#		#points(sup.p, best.fit.on.s, type = 'l', col = 'gray', lwd = 1)
#		points(sup.p, m.t.p, type = 'l', col = 'steelblue3', lwd = 3)
#		points(sup.p, s.t.p, type = 'l', col = 'green3', lwd = 3)
#		#points(sup.p, s.t.p.a, type = 'l', col = 'pink', lwd = 2)
#		#points(sup, s.t, type = 'l', col = 'red', lwd = 1)
#
#		abline(h = 0, col = 'red')
#		axis(1, at = sup);axis(2)
#		
#		
#		
#		title = paste('Replicate', rep ,' - ',T$gene[i], '\n error on s using the model : ' , round(f2min.for.cosine.8h(par = c(eps.s,phase.s.24, nu.s, phase.s.12, mu.s, phase.s.8), var= S),digits = 2),' ||  error for best fit for s : ', round(err.best,digits = 2))
#		
#		ylim = range(0,s.t,m.t,S,M)
#		plot(1,1, type = 'n', xlim = c(0,24), ylim = ylim , axes = FALSE, xlab = 'time [h of a day]', ylab = 's(t)', main = title)
#		abline(h = 1, col = 'gray')
#		points(sup, S, pch = 5, col = 'green3', type = 'b')
#		points(sup.p, best.fit.on.s, type = 'l', col = 'gray', lwd = 1)
#		points(sup.p, s.t.p, type = 'l', col = 'green3', lwd = 3)
#		abline(h = 0, col = 'red')
#		axis(1, at = sup);axis(2)
#	}
#	return(list(opt.par = n.par, error.s = ERROR, error.m = error.m, errors = errors[-1], N.PAR = N.PAR[-1,]))
#}
#
#
#grid.it = function(rep = 1,n = 3, add.penalty = TRUE){
#	if(rep == 1){T = T1}else{T = T2}
#
#	S = unlist(T[i, ZT.int]) 
#	M = unlist(T[i,ZT.ex])
#
#	if(n == 3){	
#		par.i = c(min(1,T$rel.ampl.ex[i]),T$phase.ex[i]*2*pi/24, 0,0)
#		out.optim = optim(par.i, f2min.for.cosine.12h, var = M, method = 'L-BFGS-B', lower = c(0,-2*pi,0,-2*pi), upper = c(1,2*pi,1,2*pi))
#		par.f = out.optim$par; eps.m = par.f[1];psi.m = par.f[2]; nu.m = par.f[3];  psi.12.m = par.f[4]
#		m.t.i = cosine.fit.12h(par = par.i, sup = sup); m.t.i.p = cosine.fit.12h(par = par.i, sup = sup.p)
#		m.t.f = cosine.fit.12h(par = par.f, sup = sup); m.t = m.t.f; m.t.p = cosine.fit.12h(par = par.f,sup = sup.p)
#		fft.m.est = fft(m.t)
#	}else{
#		eps.m = min(1,T$rel.ampl.ex[i]); psi.m = T$phase.ex[i]*2*pi/24; nu.m = 0; psi.12.m =0;
#		m.t = cosine.fit.24h(par = c(eps.m, psi.m), sup = sup); m.t.p = cosine.fit.24h(par = c(eps.m, psi.m), sup = sup.p);
#		fft.m.est = fft(m.t); m.t.i = m.t; m.t.i.p = m.t.p
#	}	
#	
#	var = c(fft.m.est,S,n)
#	
#	n.gamma = 21; n.eps = 21; n.psi = 21
#	gammas = log(2)/exp(seq(log(10/60),log(15),len = n.gamma))
#	eps_gs = seq(0,1,len = n.eps)
#	psi_g.s = seq(0,2*pi, len = n.psi)
#	
#	ERROR = array(NA,dim = c(n.gamma, n.eps, n.psi)); dimnames(ERROR)[[1]] = gammas; dimnames(ERROR)[[2]] = eps_gs; dimnames(ERROR)[[3]] = psi_g.s; 
#	for(gamma in gammas){
#		cat('\t gamma : ',gamma, '\n')
#		for(eps_g in eps_gs){
#			cat('\t\t eps_g : ', eps_g, '\n')
#			for(psi_g. in psi_g.s){
#				cat('\t\t\t psi_g. : ', if(psi_g. == 0){'0.00'}else{round(psi_g.,digits = 2)}, '\t\n')
#				
#				par =  c(gamma, eps_g, psi_g.)
#				error =  f2min(par = par, var = var, add.penalty = add.penalty)
#				ERROR[as.character(gamma),as.character(eps_g),as.character(psi_g.)] = error
#				
#			}	
#		}
#	}
#	return(ERROR)
#}	
#
#
#plot.landscape = function(ERROR, opt.par){
#	
#	colorlow = 'black'; colormed0 = 'red'; colormed1 = 'deeppink' ; colormed2 = 'pink';colormed3 = 'white';colorhigh  = 'white' #colorlow = rgb(1,1,0.7,0.5)
#	col.ramp <- colorRampPalette(c(colorlow, rep(c(colormed0,colormed1,colormed2,colormed3), each = 3), colorhigh), space = "Lab")
#	ERR = log10(ERROR)
#	zlim = range(log10(ERROR)); 
#	dims = dim(ERROR); n.gamma = dims[1]; n.eps = dims[2]; n.psi = dims[3]
#	gammas = as.numeric(dimnames(ERROR)[[1]]); eps_gs = as.numeric(dimnames(ERROR)[[2]]); psi_g.s = as.numeric(dimnames(ERROR)[[3]]); 
#
#	for(j in 1:max(n.gamma,n.eps, n.psi)){
#		if(j<=n.gamma){
#		image(ERR[j,,], zlim = zlim, col = col.ramp(100), ylab = 'psi_g', xlab = 'eps_g', main = paste('gamma : ',round(gammas[j],digits = 2)), axes = FALSE)
#		axis(1, at = seq(0,1,len = n.eps)); axis(2, at = seq(0,1,len = n.psi), labels = round(psi_g.s,digits = 2)) ;
#		dist = 1-abs(log(opt.par[1])-log(gammas[j]))/max(abs(log(opt.par[1])-log(gammas)))
#		points(opt.par[2],(opt.par[3]%%(2*pi))/(2*pi), pch = 16, cex = 2*dist, col = rgb(0,1,1,dist))
#		}else{plot.new()}
#
#		if(j<=n.eps){
#		image(ERR[,j,], zlim = zlim, col = col.ramp(100), ylab = 'psi_g', xlab = 'gamma', main = paste('eps_g : ',round(eps_gs[j],digits = 2)), axes = FALSE)
#		axis(1, at = seq(0,1,len = n.gamma), labels = round(gammas,digits = 2)); axis(2, at = seq(0,1,len = n.psi), labels = round(psi_g.s,digits = 2)) 
#		dist = 1-abs(opt.par[2]-eps_gs[j])/max(abs(opt.par[2]-eps_gs))
#		points(approx(gammas,y = seq(0,1,len = n.gamma), xout = opt.par[1])$y,(opt.par[3]%%(2*pi))/(2*pi), pch = 16, cex = 2*dist, col = rgb(0,1,1,dist))
#		}else{plot.new()}
#		
#		if(j<=n.psi){
#		image(ERR[,,j], zlim = zlim, col = col.ramp(100), ylab = 'eps_g', xlab = 'gamma', main = paste('psi_g. : ',round(psi_g.s[j],digits = 2)), axes = FALSE)
#		axis(1, at = seq(0,1,len = n.gamma), labels = round(gammas,digits = 2)); axis(2, at = seq(0,1,len = n.eps)); 
#		opt.phase = opt.par[3]; diff = abs(opt.phase-psi_g.s)%%(2*pi); diff[diff>pi] = 2*pi - diff[diff>pi]
#		dist = 1-diff[j]/max(diff)
#		points(approx(gammas,y = seq(0,1,len = n.gamma), xout = opt.par[1])$y,opt.par[2], pch = 16, cex = 2*dist, col = rgb(0,1,1,dist))
#		}else{plot.new()}
#		
#		}
#	
#}
#


##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
# OLDIES
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################




signal.selection.v2 = function(t = t1, gene_list = gene_list){

	N = length(gene_list)	
	
	EX_SIGNAL = matrix(NA,N,12); INT_SIGNAL = EX_SIGNAL;colnames(EX_SIGNAL) = paste(colnames(t)[ZTr],'.ex',sep=''); colnames(INT_SIGNAL) = paste(colnames(t)[ZTr],'.int', sep = ''); 
		
	gene_signal = data.frame(gene = gene_list, EX_SIGNAL, INT_SIGNAL,exon.median = rep(NA,N), intron.median = rep(NA,N), stringsAsFactors = FALSE)
	ZT_ex = ZT; ZT_int = ZT+12

	for(i in 1:N){
		if((i%%300)==0){cat(round(100*i/N, digits = 1),'%\n')}
		gene = gene_list[i]; #cat(gene, '\n')
		ind = select(t = t,gene = gene,probe = 'exon')
		ex.out = gene.signal.selection.v2(t = t,ind = ind)
		gene_signal$exon.median[i] = median(unlist(t[ind,ZT]))
		ind = select(t = t,gene = gene,probe = 'intron')
		int.out = gene.signal.selection.v2(t = t,ind = ind)
		gene_signal$intron.median[i] = median(unlist(t[ind,ZT]))
		gene_signal[i,ZT_ex] = ex.out
		gene_signal[i,ZT_int] = int.out
		}
	return(gene_signal)
	}

gene.signal.selection.v2 = function(t = t,ind = ind){
	out = rep(NA,12)
	n = sum(ind,na.rm = TRUE)
	if(n == 0){return(out)}else if(n == 1){return(as.vector(t[ind,ZTr]))}else if(n >= 2){out = apply(as.matrix(t[ind,ZTr]),2,median); out = out/mean(out); return(out)}
	}		
			


gene.signal.selection.v1 = function(t = t,ind = ind){
	out = rep(NA,12)
	n = sum(ind,na.rm = TRUE)
	if(n == 0){return(out)}else if(n == 1){return(as.vector(t[ind,ZT]))}else if(n == 2){return(apply(as.matrix(t[ind,ZT]),2,median))}else{
		T = as.matrix(t[ind,ZT])
		median = apply(T,2,median, na.rm=TRUE)
		COR = vector(mode='numeric', length = n);for(k in 1:n){COR[k] = cor(T[k,],median)}; ok = (COR >= 0.9); nok = sum(ok, na.rm = TRUE)
		if(((nok/n)>=0.25)&(nok>1)){out = apply(T[ok,],2,median,na.rm =TRUE)} #at least 25 % of the probes and more than 2 probes are like the median
		else{out = median}
		}
	return(out)	
	}		
				
	

	
	
	
#		e_m = min(1,T$rel.ampl.ex[i]); psi_m =T$phase.ex[i]*2*pi/24 ; w = 2*pi/24; g = n.par[1]; e_g = n.par[2]; 
#		psi_g.f = n.par[3]; if(psi_g.f<psi_m){psi_g.f = psi_g.f + 2*pi} 
#		psi_g = (psi_g.f - psi_m)
#		e_s = e_m * (sqrt(w^2 + g^2))/(g*(1+0.5*e_m*e_g*cos(psi_g)));phi = atan(w/g)
#		F = 1/(1+ 0.5*e_m*e_g*cos(psi_g) )
#		s.t.p.a = 1 +  F * ( e_m * (sqrt(w^2 + g^2))/g * cos(w*sup.p  - psi_m + phi) + e_g * cos(w*sup.p - psi_m - psi_g) + 0.5 * e_g * e_m * cos(2*w*sup.p - 2*psi_m - psi_g) )
	
	
	
	
	
plot.distrib.t = function(t = t, zt = 1:24, zt.rel = zt.rel, main = "just out from RMA"){
	
	i.all = 1:nrow(t)
	i.int = which((t$ex_in == -1)&(t$kept))
	i.ex = which((t$ex_in == 1)&(t$kept))
	i.other = which((!((t$ex_in == -1)|(t$ex_in == 1)))|(!t$kept))
	
	
	i.spec = 1:nrow(t);
	plot.densities(t = t, zt = zt, i.spec = 1:nrow(t), main = paste(main,'- all probesets'), type = 'abs')
	plot.densities(t = t, zt = zt, i.spec = i.int, main = paste(main,'- intronic probesets'), type = 'abs')
	plot.densities(t = t, zt = zt, i.spec = i.ex, main = paste(main,'- exonic probesets'), type = 'abs')
	plot.densities(t = t, zt = zt, i.spec = i.other, main = paste(main,'- other probesets'), type = 'abs')
	
	plot.densities(t = t, zt = zt.rel, i.spec = 1:nrow(t), main = paste(main,'- all probesets'), type = 'rel')
	plot.densities(t = t, zt = zt.rel, i.spec = i.int, main = paste(main,'- intronic probesets'), type = 'rel')
	plot.densities(t = t, zt = zt.rel, i.spec = i.ex, main = paste(main,'- exonic probesets'), type = 'rel')
	plot.densities(t = t, zt = zt.rel, i.spec = i.other, main = paste(main,'- other probesets'), type = 'rel')
	
	
	plot(1,1,type = 'n', xlab = '', ylab = '', axes = FALSE)

	plot.abs.rel(t = t, zt = zt, zt.rel = zt.rel, i.spec = i.int, main = paste(main,'- intronic probesets'))
	plot.abs.rel(t = t, zt = zt, zt.rel = zt.rel, i.spec = i.ex, main = paste(main,'- exonic probesets'))	
	plot.abs.rel(t = t, zt = zt, zt.rel = zt.rel, i.spec = i.other, main = paste(main,'- other probesets'))	
		
}
	

plot.densities = function(t = t, zt = zt, i.spec =  1:nrow(t), main = 'title', type = c('abs','rel')){
	
	
	type = type[1]
	cat('plot densities -',type,'-', main,'\n')

	rainbow = rainbow(12,s = 0.85, v = 0.85)

	
	d = density(t[i.spec,zt[1]])

	if(type == 'abs'){xlab =  'abs. signal (log2)'; xlim = c(0,15)}else{xlab =  'rel. signal'; xlim = c(0,2.2)}
	ylim = c(0,1.05*max(d$y))
	plot(1,1,type = 'n', xlab = xlab, ylab = 'density', xlim = xlim, ylim = ylim, main = main)
	abline(v = 1)
	for(i in zt){
		i.rel = i-zt[1]+1; cat(i.rel,"\n")
		d = density(t[i.spec,i])
		points(d, type = 'l', col = rainbow[(i.rel-1)%%12+1], lty = (i.rel-1)%/%12+1)
		legend(x = 0.8*xlim[2],y = ylim[2]*((i.rel+1)/length(zt)), legend = paste('ZT',2*(i.rel-1)), col =  rainbow[(i.rel-1)%%12+1], lty = (i.rel-1)%/%12+1, bty = 'n' )
	}
	
}	
	
	
plot.abs.rel = function(t = t, zt = zt, zt.rel = zt.rel, i.spec = i.int, main = 'title'){
	
	cat('plot abs. vs rel. -', main,'\n')

	
	rainbow = rainbow(12,s = 0.85, v = 0.85)

	plot(1,1,type = 'n', xlab = 'abs. signal (log2)', ylab = 'rel. signal', xlim = c(1,12), ylim = c(0.3,1.8), main = main)
	
	for(i in 1:length(zt)){
		j = zt[i]; k = zt.rel[i];
		#points(t[i.spec,j],t[i.spec,k], cex = 0.3, pch = i%/%12+16, col = rainbow[(i-1)%%12+1])
		d = kde2d(x = t[i.spec,j], y = t[i.spec,k], h = c(3,0.5))
		contour(d, col = rainbow[(i-1)%%12+1],lty = (i-1)%/%12+1, nlevels = 5, label = paste('ZT',2*(i-1), sep = ''),drawlabels = TRUE, add = TRUE)
	}
}
	
