library(emdbook)
library(deSolve)
library(fdrtool)
library(circular)
library(preprocessCore)
library(gtools)
library(biomaRt)
#library(plotrix)

######## FUNCTIONS FOR SIMULATED DATA and for the optimization and model selection

set.global.variable = function()
{
	splicing.rate <<- log(2)/8*60
}

set.global.sigma = function() ## this function is to used to generate fake data with different standard deviation
{
	sigma.s <<- 0.2149748
	sigma.m <<- 0.1760124
}


#####
##### fake data generating
#####
generate.fake.data = function(T = T[1:100,],X = 100, sd = 0.05, model = 1, parametrization = c('cosine','sigmoid'), sum.species = TRUE, random.splicing=FALSE)
{
	#parametrization = c('sigmoid');random.splicing=FALSE;sum.species = FALSE;X=100;
	parametrization = parametrization[1]
	set.seed(42)
	
	## Set parameter of splicing rate 'k' which is 8 minutes 
	set.global.variable(); 
	k = splicing.rate;
	if(random.splicing) {kk=sample(seq(5,8.3,by=0.042), X, replace=TRUE)}
	
	
	T$gene = paste('fake_m',model,'_',1:X,sep = '')
	T$intron.median = 2^rnorm(X,mean = 5.7, sd = 1.4) ### first generate the median of intronic signals
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
		## sigmoid function
		fold.change.min = 1.25; 
		fold.change.max = upper[4];
		phase.int.min = lower[5]; 
		phase.int.max = upper[5];
		up.time.min = lower[6]; 
		up.time.max = upper[6];
		down.time.min = lower[6]; 
		down.time.max = upper[6];
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
			#T$down.time.int = pmax(apply(cbind(rnorm(n = X, mean = 12, sd = 4.5),24-T$up.time),1,min),down.time.min)
			T$down.time.int = T$up.time.int
			T$rel.ampl.12.int = rep(0,X)
			T$phase.12.int = rep(0,X)
		}
	}else{
		T$rel.ampl.int = rep(0,X);
		T$rel.ampl.12.int = rep(0,X);
		T$phase.int = rep(0,X); 
		T$phase.12.int = rep(0,X)
		if(parametrization == 'sigmoid'){
			T$fold.change.int = rep(1,X); 
			T$phase.int = 0; 
			T$up.time.int = rep(12,X); 
			T$down.time.int = T$up.time.int
		}
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
			T[i,ZT.int] = compute.sigmoid(t = zt, fold.change = T$fold.change.int[i], phase = T$phase.int[i], up.time = T$up.time.int[i], down.time = T$up.time.int[i])
			
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
												   down.time = T$up.time.int[i], rescale = FALSE, synthesis.factor = kk[i] * T$intron.median[i])
				}else{
					abs.ex[i,] = compute.m.sigmoid(t = zt, gamma = T$gamma[i],eps.gamma = T$eps.gamma[i],phase.gamma = T$phase.gamma[i], 
												   fold.change = T$fold.change.int[i], phase = T$phase.int[i], up.time = T$up.time.int[i],
												   down.time = T$up.time.int[i], rescale = FALSE, synthesis.factor = k * T$intron.median[i])
				}
				
				#cat('finished!\n')
			}
		}	
		abs.int = T[,ZT.int] * T$intron.median ## absolute singals for introns
		
	}
	
	if(sum.species){abs.ex = abs.ex+abs.int}   ## here to sum up the exonic signals from premrnas and mrnas
	
	log.abs.int = log(abs.int);
	log.abs.ex = log(abs.ex);
	
	## add noise
	set.global.sigma();
	sd.noise = c(0.1,0.25,0.5);
	if(sd==0.1) {sigma.s = sigma.s/5; sigma.m = sigma.m/5;}
	if(sd==0.25) {sigma.s = sigma.s/2; sigma.m = sigma.m/2;}
	if(sd==0.5) {sigma.s = sigma.s; sigma.m = sigma.m;}
	
	noise.int = matrix(rnorm(X*length(ZT.int), mean = 0, sd = sigma.s),nrow = X, ncol = length(ZT.int))
	noise.ex = matrix(rnorm(X*length(ZT.ex), mean = 0, sd = sigma.m),nrow = X, ncol = length(ZT.ex))
	
	## add different background for introns and exons
	cc.int =  matrix(rep(rnorm(X, mean = -2, sd = 1.25), length(ZT.int)),nrow = X, byrow = FALSE)
	
	log.noise.abs.int = log.abs.int + noise.int
	#log.noise.abs.int = log.abs.int 
	log.noise.abs.ex = log.abs.ex + cc.int + noise.ex
	#log.noise.abs.ex = log.abs.ex + noise.ex 
	#log.noise.abs.ex = log.abs.ex
	
	noise.abs.int = exp(log.noise.abs.int)
	noise.abs.ex = exp(log.noise.abs.ex)
	
	T[,ZT.int] = noise.abs.int/apply(noise.abs.int, 1, mean)
	T[,ZT.ex] =  noise.abs.ex/apply(noise.abs.ex, 1, mean)
	T$exon.median = apply(noise.abs.ex,1,median)		
	return(T)
	
}


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

fs.n = function(x, n = 1){
	mx = mean(x)
	x = x-mx
	P = (Mod(fft(x)))^2 
	P[n+1]/sum(P[2:floor((length(x)+1)/2)])
}

pval.n = function(x, n = 1){
	S = fs.n(x, n = n)
	pval = (1-S)^(length(x)/2-2)
	return(pval)
}

compare.parameters = function(T = T, fitting.method = c('Optim','MCMC'), quality = FALSE, parametrization = c('sigmoid','cosine'))
{
	#fitting.method = c('Optim','MCMC'); quality = FALSE; parametrization = c('sigmoid','cosine');
	parametrization = parametrization[1]
	fitting.method = fitting.method[1]
	if(fitting.method == 'Optim'){
		best.model = T$LRT.best.model;
		true.model=T$true.model;
	}else{
		best.model = T$Model.MCMC;
		true.model=T$true.model;
	}
	
	if(parametrization == 'cosine') # not used
	{
		mains = c('Degr. Rate','Degr. Amplitude','Degr. Phase','Rel. Ampl. - Intron','Phase 24h - Intron','Rel. Ampl. 12h - Intron','Phase 12h - Intron')
		parameter.list = c('gamma','eps.gamma','phase.gamma','rel.ampl.int','phase.int','rel.ampl.12.int','phase.12.int')
	}else{
		mains = c('Degr. Rate','Degr. FC','Degr. Phase','FC - Intron','Phase - Intron','Up-time - Intron','Down-time - Intron')
		parameter.list = c('gamma','eps.gamma','phase.gamma','fold.change.int','phase.int','up.time.int','up.time.int')
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
		if(p==7)
		{
			p = 6;
			plot(c(lower[p], upper[p]),c(lower[p],upper[p]),type = 'n', xlab = 'True values', ylab = 'Estimated Values', main = mains[p+1], log = log)
		}else{
			plot(c(lower[p], upper[p]),c(lower[p],upper[p]),type = 'n', xlab = 'True values', ylab = 'Estimated Values', main = mains[p], log = log)
		}
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


##################################################################################################################
##################################################################################################################
################################################################################################################## Functions for the optimization
##################################################################################################################
##################################################################################################################

### Main Function for model selection and fitting
make.fits.with.all.models.for.one.gene = function(T = T, gene.index = 1, debug = FALSE, parametrization = c('cosine','sigmoid'), method =  c('integration','simulation'), underdetermination.model3 = FALSE, sum.species = FALSE, zt = seq(0,46,by = 2), i.ex = ZT.ex, i.int = ZT.int, absolute.signal = FALSE)
{
	#T = T; gene.index = j; debug = TRUE; parametrization = 'sigmoid'; method = 'integration'; sum.species = TRUE; zt = zt; i.ex = ZT.ex; i.int = ZT.int; absolute.signal = FALSE;
	parametrization = parametrization[1]; 
	method = method[1];
	for(model in 1:4)
	{
		if(debug){cat('\t starting model ',model,'\n');}
		param.fit = make.fit.spec.model(T = T, gene.index = gene.index, model = model, debug = debug, parametrization = parametrization, method = method, underdetermination.model3=underdetermination.model3, sum.species = sum.species, zt = zt, i.ex = i.ex, i.int = i.int, absolute.signal = absolute.signal);
		
		if(model == 1)
		{
			Param.fit.for.gene = param.fit
		}else{
			Param.fit.for.gene = c(Param.fit.for.gene, param.fit)
		}
		if(debug){cat('\t model ',model,' finished \n')};
	}
	## model selction for one gene
	#BIC = my.BIC.this.gene(Param.fit.for.gene)
	#Param.fit.for.gene = c(Param.fit.for.gene, BIC)
	return(Param.fit.for.gene)
}

## BIC model-selection for one gene
my.BIC.this.gene = function(param.fits.results, diff.sigma=TRUE)
{
	set.nb.param();
	index = match(c('error.m1', 'error.m2', 'error.m3', 'error.m4'), names(param.fits.results))
	error.m1 = param.fits.results[index[1]]
	error.m2 = param.fits.results[index[2]]
	error.m3 = param.fits.results[index[3]]
	error.m4 = param.fits.results[index[4]]
	## the formula of BIC used here is chi-square+k*ln(n)==error/sigma^2+k(ln(n)) in which sigma of noise is supported to be known.
	BIC.m1 = log(48)*n.param[1] + 48*log(error.m1)
	BIC.m2 = log(48)*n.param[2] + 48*log(error.m2)
	BIC.m3 = log(48)*n.param[3] + 48*log(error.m3)
	BIC.m4 = log(48)*n.param[4] + 48*log(error.m4)
	
	BIC = c(BIC.m1, BIC.m2, BIC.m3, BIC.m4)
	BIC = c(BIC, which.min(BIC))
	names(BIC) = c('BIC.m1', 'BIC.m2', 'BIC.m3', 'BIC.m4', 'BIC.best.model')
	
	return(BIC)								
}	


make.fit.spec.model = function(T = T, gene.index = 1, model = 1, debug = FALSE, parametrization = c('cosine','sigmoid'),method =  c('integration','simulation'), underdetermination.model3 = FALSE, sum.species = FALSE, zt = seq(0,46,by = 2), i.ex = ZT.ex, i.int = ZT.int, absolute.signal = FALSE)
{
	parametrization = parametrization[1];
	method = method[1];
	param.fit = NA
	S = unlist(T[gene.index, i.int])
	M = unlist(T[gene.index, i.ex])
	## Model 1 premrna and mran are both flat
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
	## parameter estimations for model 2,3,4
	if(model > 1)
	{
		param.fit = make.optimization(T = T, i = gene.index, model = model, Nfit = NA, debug = debug,  parametrization = parametrization,method = method, underdetermination.model3=underdetermination.model3, sum.species = sum.species, zt = zt, i.ex = i.ex, i.int = i.int, absolute.signal = absolute.signal)
	}
	return(param.fit)
}

## use two replicates to estimate gene-specific noise, but does not work well
noise.estimate = function(data)
{
	nb = length(data)
	half = nb/2
	xx = data[1:half]-data[(half+1):nb]
	return(sd(xx)/sqrt(2))		
}

error = function(S,s,M,m, log = TRUE)
{
	S = as.numeric(unlist(S)); 
	s = as.numeric(unlist(s)); 
	M =  as.numeric(unlist(M)); 
	m =  as.numeric(unlist(m));
	
	## here we have to use the estimated values of sigma.ss and sigma.mm
	set.global.sigma();
	sigma.ss = sigma.s;
	sigma.mm = sigma.m;
	
	if(log)
	{
		if(any(S<=0)|any(s<=0)){
			error.S = 10^8;
		}else{
			error.S = sum((log(S)-log(s))^2);
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
		}
	}else{
		error.M = sum((M-m)^2)
	}
	
	error = error.M + error.S*sigma.mm^2/sigma.ss^2
	#error = error.M + error.S
	
	#error = error.M*sigma.ss^2/sigma.mm^2 + error.S
	#cat('error of S is ', error.S,'\n');
	#cat('error of M is ', error.M,'\n');
	#cat('error of total is ', error,'\n');
	
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
		#down.time = param[model+3];
		down.time = up.time;
		
		## here is to constrain the sum of uptime and downtime <24.1h
		#constrain.param = min(sum(exp(alpha*(up.time+down.time-24.1))),10^8)
		constrain.param = 0.0;
		
	}else{
		constrain.param = 0.0;
	}
	if(intense.debug){cat('______________________________-------- constrain.param = ', constrain.param,'\n')}
	
	error = error(S,s,M,m,log = log)
	errors = error + constrain.s.positive + constrain.param
	return(errors)	
}

## function to estimat parameters for model 2, 3, 4
make.optimization = function(T = T, i = 1, model = 2, Nfit = NA, debug = FALSE, parametrization =c('cosine','sigmoid'), method =  c('integration','simulation'), underdetermination.model3=FALSE, sum.species = FALSE, zt =  seq(0,46,by = 2), i.ex = ZT.ex, i.int = ZT.int, absolute.signal = FALSE)
{
	#i = j; zt =  seq(0,46,by = 2); i.ex = ZT.ex; i.int = ZT.int;absolute.signal = FALSE; Nfit=NA; parametrization =c('sigmoid'); method =  c('integration','simulation');sum.species = FALSE;
	
	parametrization = parametrization[1];
	method = method[1]
	
	param.fit = NA
	S = unlist(T[i, i.int])
	M = unlist(T[i, i.ex])
	
	bounds = set.bounds(model = model, parametrization = parametrization, absolute.signal = absolute.signal); 
	upper = bounds$upper; 
	lower = bounds$lower;
	
	fit.separate = TRUE
	
	if((model==2|model==4) & fit.separate)
	{
		## choose different initial values of parameters for S fitting
		#Nfit.S = 10
		Nfit.S = 6
		index = i
		#gamma.init = c(rep(log(2)/lseq(log(2)/lower[1],log(2)/upper[1], length = 5), Nfit%/%5), rep(log(2)/5,Nfit%%5))
		
		if(!is.null(T$rel.ampl.int[index]))
		{
			rel.ampl.int = T$rel.ampl.int[index]; 
			phase.int = T$phase.int[index]
		}else{
			rel.ampl.int = (max(S)-min(S))/2;
			phase.int = zt[which.max(S)]
		}
		
		fold.change.init = rep(min((1+ rel.ampl.int)/(1-rel.ampl.int),1000),Nfit.S);
		fold.change.init = fold.change.init*sample(seq(1/max(fold.change.init),2-(1/max(fold.change.init)),len = 1000),Nfit.S)
		phase.init = (rep(phase.int,Nfit.S)+rnorm(Nfit.S,sd = 2))%%24
		
		times.min = 2
		times.max = 12
		up.time.init = c(lseq(times.min,times.max, length = Nfit.S))
		PAR.INIT.S = cbind(fold.change.init, phase.init, up.time.init)
		
		## fitting S 
		errors.fit.s = rep(NA, Nfit.S)
		for(fit.nb.s in 1:Nfit.S)
		{
			par.init.s = PAR.INIT.S[fit.nb.s,]
			opt.s = optim(par.init.s, f2min.int, S=S, absolute.signal = absolute.signal, zt = zt, method = 'L-BFGS-B', lower = c(1,0,2), upper = c(1000, 24, 12))
			res.fit.s = opt.s$par
			errors.fit.s[fit.nb.s] = opt.s$value
			eval(parse(text = paste('res.fit.s.', fit.nb.s, ' = res.fit.s', sep = '')))
			
		}
		
		## choose the best-fitting parameters for S and select the best model with BIC
		imin.s = which.min(errors.fit.s); 
		eval(parse(text = paste('res.fit.s = res.fit.s.', imin.s, sep = '')))
		
		s = rep(1,length(S)); 
		m =  M;
		err.s = error(S,s,M,m); 
		
		if(model==2)
		{
			param.fit = c(err.s, errors.fit.s[imin.s], res.fit.s)
			BIC.s13 = log(24)*0 + 24*log(param.fit[1])
			BIC.s24= log(24)*3 + 24*log(param.fit[2])
		
			param.fit = c(param.fit, BIC.s13, BIC.s24)
			names(param.fit) = c('error.s13', 'error.s24','fold.change.int.s24', 'phase.int.s24', 'up.time.int.s24', 'BIC.s13', 'BIC.s24')
		}
		## choose different initial conditions to fit S and M according to above results
		if(is.na(Nfit))
		{
			#Nfit = 10
			Nfit = 6
		}
		
		gamma.init = c(rep(log(2)/lseq(log(2)/lower[1],log(2)/upper[1], length = Nfit), Nfit%/%Nfit), rep(log(2)/5,Nfit%%Nfit))
		
		eps.gamma.init = rep(0.2, Nfit); 
		if(parametrization == 'sigmoid')
		{
			eps.gamma.init = (1+ eps.gamma.init)/(1-eps.gamma.init); 
			eps.gamma.init = eps.gamma.init*sample(seq(1/max(eps.gamma.init),2-(1/max(eps.gamma.init)),len = 1000),Nfit)
		}
		
		phases.gamma.init = c(rep(c(6,12,18), round(Nfit/3)), rep(12, Nfit%%3))
		
		fold.change.init = lseq(max(1, res.fit.s[1]/2), min(res.fit.s[1]*2, 1000), length=Nfit)
		phase.init = (rep(res.fit.s[2],Nfit)+rnorm(Nfit,sd = 2))%%24
		up.time.init = lseq(max(2, (res.fit.s[3]-2)), min((res.fit.s[3]+2), 12), length=Nfit)
		
		if(model==2) {PAR.INIT = cbind(gamma.init, fold.change.init, phase.init, up.time.init); colnames(PAR.INIT) = c('gamma','fold.change.int','phase.int','up.time.int');}
		if(model==4) {
			PAR.INIT = cbind(gamma.init, eps.gamma.init, phases.gamma.init, fold.change.init, phase.init, up.time.init); 
			colnames(PAR.INIT) = c('gamma','eps.gamma','phase.gamma','fold.change.int','phase.int','up.time.int')
		}
		
		## fit the S and M 
		prior.factor.gamma = 1;
		errors.fit = rep(NA, Nfit)
		
		for(fit.number in 1:Nfit)
		{
			##cat(fit.number, '\n');
			par.init = PAR.INIT[fit.number,]
			if(debug){cat('\t\t starting optimization # ',fit.number,' for model ', model,':',par.init,'\n')}
			
			## the optimization function used in the fitting
			## f2min is the likelihood function to minimize
			opt = optim(par.init, f2min, M = M, S = S, prior.factor.gamma= prior.factor.gamma, model = model, parametrization = parametrization, method.intern = method, debug = debug, underdetermination.model3 = underdetermination.model3, sum.species = sum.species, zt = zt, absolute.signal = absolute.signal,  
						method = 'L-BFGS-B', lower = lower, upper = upper, hessian = TRUE)
			
			ttry = try(sqrt(diag(solve(opt$hessian))), silent = TRUE)
			if(!inherits(ttry, "try-error")){
				opt$stderr = sqrt(diag(solve(opt$hessian)))
			}else{
				opt$stderr = rep(NA, length(opt$par))
			}
			
			## extract the parameter of optimization results and errors
			res.fit = opt$par
			res.fit.stderr = opt$stderr
			#if(model==3 & !sum.species & underdetermination.model3) {print('HERE');res.fit[1]=prior.factor.gamma/(1+res.fit[2]); res.fit.stderr[1]=res.fit.stderr[2]*prior.factor.gamma/(1+res.fit[2])^2;}
			errors.fit[fit.number] = opt$value
			eval(parse(text = paste('res.fit.', fit.number, ' = res.fit', sep = '')))
			eval(parse(text = paste('res.fit.stderr.', fit.number, ' = res.fit.stderr', sep = '')))
			if(debug){cat('\t\t optimization # ',fit.number,' finished : \t\t\t',res.fit,'\t',opt$value,'\n')}
		}
		
		## choose the best estimated parameters
		imin = which.min(errors.fit); 
		eval(parse(text = paste('res.fit = res.fit.', imin, sep = '')))
		eval(parse(text = paste('res.fit.stderr = res.fit.stderr.', imin, sep = '')))
		
		param.fit1 = c(errors.fit[imin], res.fit, res.fit.stderr);
		names(param.fit1) = paste(c('error',colnames(PAR.INIT), paste(colnames(PAR.INIT), '.stderr', sep='')),'.m',model,sep = '')
		
		if(model==2) param.fit = c(param.fit, param.fit1)
		if(model==4) param.fit = param.fit1
		
		
	}else{
	# Nfit the number of optimization which is twice of nb of parameters and each optimization is from different initial conditions
	#if(is.na(Nfit)){Nfit = 2*length(upper)}
	if(is.na(Nfit))
	{
		if(model==1) Nfit = 2*length(upper)
		if(model==2) Nfit = 5*5
		if(model==3) Nfit = 5*2
		if(model==4) Nfit = 5*5
		
	}
	if(model==3 & !sum.species & underdetermination.model3)
	{
		#print('HERE 1')
		set.global.variable(); 
		k.splicing=splicing.rate; 
		index.exon = grep('exon.median', colnames(T));
		index.intron = grep('intron.median', colnames(T));
		prior.factor.gamma=2*k.splicing*T[i,index.intron]/T[i, index.exon];
	}else{
		prior.factor.gamma = 1;
	}
	
	#Nfit times of initialization
	PAR.INIT = do.parameter.initialization(T = T, S = S, M = M, model = model, index = i, Nfit = Nfit, parametrization = parametrization, absolute.signal = absolute.signal, zt = zt)
	
	
	errors.fit = rep(NA, Nfit)
	debug = TRUE
	for(fit.number in 1:Nfit)
	{
		##cat(fit.number, '\n');
		par.init = PAR.INIT[fit.number,]
		if(debug){cat('\t\t starting optimization # ',fit.number,' for model ', model,':',par.init,'\n')}
		
		## the optimization function used in the fitting
		## f2min is the likelihood function to minimize
		opt = optim(par.init, f2min, M = M, S = S, prior.factor.gamma= prior.factor.gamma, model = model, parametrization = parametrization, method.intern = method, debug = debug, 
					underdetermination.model3 = underdetermination.model3, sum.species = sum.species, zt = zt, absolute.signal = absolute.signal,  
					method = 'L-BFGS-B', lower = lower, upper = upper, hessian = TRUE)
		
		ttry = try(sqrt(diag(solve(opt$hessian))), silent = TRUE)
		if(!inherits(ttry, "try-error")){
			opt$stderr = sqrt(diag(solve(opt$hessian)))
		}else{
			opt$stderr = rep(NA, length(opt$par))
		}
		
		## extract the parameter of optimization results and errors
		res.fit = opt$par
		res.fit.stderr = opt$stderr
		if(model==3 & !sum.species & underdetermination.model3) {res.fit[1]=prior.factor.gamma/(1+res.fit[2]); res.fit.stderr[1]=res.fit.stderr[2]*prior.factor.gamma/(1+res.fit[2])^2;}
		errors.fit[fit.number] = opt$value
		eval(parse(text = paste('res.fit.', fit.number, ' = res.fit', sep = '')))
		eval(parse(text = paste('res.fit.stderr.', fit.number, ' = res.fit.stderr', sep = '')))
		if(debug){cat('\t\t optimization # ',fit.number,' finished : \t\t\t',res.fit,'\t',opt$value,'\n')}
	}
	
	imin = which.min(errors.fit); 
	eval(parse(text = paste('res.fit = res.fit.', imin, sep = '')))
	eval(parse(text = paste('res.fit.stderr = res.fit.stderr.', imin, sep = '')))
	
	## save the minimal error and corresponding parameters
	param.fit = c(errors.fit[imin], res.fit, res.fit.stderr);
	names(param.fit) = paste(c('error',colnames(PAR.INIT), paste(colnames(PAR.INIT), '.stderr', sep='')),'.m',model,sep = '')
	
	}
	
	return(param.fit);
	if(debug){cat('\t\t\t optimization finished\n')}
}


## set the boundary of parameters
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
		## if we consider the same uptime and down time, the boundary is between 6 and 12 h
		#up.time.max = 22; 
		up.time.max = 12; 
		up.time.min = 2; 
		if(absolute.signal){up.time.max = 12; up.time.min = 4}
		down.time.max = 12;
		down.time.min = 2; 
		if(absolute.signal){down.time.max = 12; down.time.min = 6}
		#param.synthesis.upper = c(fold.change.max, phase.max, up.time.max, down.time.max)
		#param.synthesis.lower = c(fold.change.min, phase.min, up.time.min, down.time.min)
		## here we use the same uptime and downtime
		param.synthesis.upper = c(fold.change.max, phase.max, up.time.max)
		param.synthesis.lower = c(fold.change.min, phase.min, up.time.min)
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

## choose the different initial values for parameters
do.parameter.initialization = function(T = T, S = rep(1,24),M = rep(1,24),model = 2, index = 1, Nfit = 5,parametrization =c('cosine','sigmoid'), absolute.signal = FALSE, zt = seq(0,46,by = 2))
{
	parametrization = parametrization[1]
	bounds = set.bounds(model = model, parametrization = parametrization, absolute.signal = absolute.signal); 
	upper = bounds$upper; 
	lower = bounds$lower;
	
	#gamma.init = rep(log(2)/120*60, Nfit)# initial half-life = 120 min
	gamma.init = c(rep(log(2)/lseq(log(2)/lower[1],log(2)/upper[1], length = 5), Nfit%/%5), rep(log(2)/5,Nfit%%5))
	#gamma.init = c(rep(log(2)/lseq(log(2)/upper[1],log(2)/lower[1], length = 3), each = Nfit%/%3),rep(log(2)/5,Nfit%%3))
	
	eps.gamma.init = rep(0.2, Nfit); 
	if(parametrization == 'sigmoid'){
		eps.gamma.init = (1+ eps.gamma.init)/(1-eps.gamma.init); 
		eps.gamma.init = eps.gamma.init*sample(seq(1/max(eps.gamma.init),2-(1/max(eps.gamma.init)),len = 1000),Nfit)
	}
	
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
		
		if(!is.null(T$rel.ampl.int[index])){
			rel.ampl.int = T$rel.ampl.int[index]; 
			phase.int = T$phase.int[index]}
		else{
			rel.ampl.int = (max(S)-min(S))/2;
			phase.int = zt[which.max(S)]
		}
		
		fold.change.init = rep(min((1+ rel.ampl.int)/(1-rel.ampl.int),1000),Nfit); 
		fold.change.init = fold.change.init*sample(seq(1/max(fold.change.init),2-(1/max(fold.change.init)),len = 1000),Nfit)
		phase.init = (rep(phase.int,Nfit)+rnorm(Nfit,sd = 2)+c(rep(0,(Nfit%/%3)*3), rep(12,Nfit%%3)))%%24
		if(model!=3){
			times.min = lower[model+2];
			times.max = upper[model+2];
		}else{
			times.min = 7; 
			times.max = 11
		}
		up.time.init = c(rep(3,Nfit%/%5), rep(5,Nfit%/%5), rep(7,Nfit%/%5),rep(9,Nfit%/%5),rep(11,Nfit%/%5), rep(8, Nfit%%5))
		#up.time.init = c(rep(seq(times.min, times.max, length = Nfit%/%3),3),rep(mean(c(times.min,times.max)),Nfit%%3));  
		
		## the same uptime and downtime for the synthesis
		#down.time.init =  c(rep(seq(times.max , times.min, length = Nfit%/%3),3),rep(mean(c(times.min,times.max)),Nfit%%3))
		#param.init.synthesis = cbind(fold.change.init, phase.init, up.time.init, down.time.init)
		param.init.synthesis = cbind(fold.change.init, phase.init, up.time.init)
		#param.init.synthesis.name = c('fold.change.int','phase.int','up.time.int','down.time.int')
		param.init.synthesis.name = c('fold.change.int','phase.int','up.time.int')
	}
	
	cat('parameter initialized\n')
	
	if(model == 2){PAR.INIT = cbind(gamma.init, param.init.synthesis); colnames(PAR.INIT) = c('gamma', param.init.synthesis.name)}
	if(model == 3){PAR.INIT = cbind(gamma.init, eps.gamma.init, phases.gamma.init); colnames(PAR.INIT) = c('gamma','eps.gamma','phase.gamma');}
	if(model == 4){PAR.INIT = cbind(gamma.init, eps.gamma.init, phases.gamma.init, param.init.synthesis); colnames(PAR.INIT) = c('gamma','eps.gamma','phase.gamma', param.init.synthesis.name)}
	return(PAR.INIT)	
}


### also in this function the sum or not sum is considered.
f2min.int = function(par.int, S, parametrization =c('sigmoid'), debug = FALSE, zt = seq(0,46,by = 2), absolute.signal = FALSE)
{
	S = as.numeric(unlist(S)); 
	#s = as.numeric(unlist(s)); 
	#M =  as.numeric(unlist(M)); 
	#m =  as.numeric(unlist(m));
	
	parametrization = parametrization[1]; 
	param.synthesis.1 = par.int[1]; 
	param.synthesis.2 = par.int[2]; 
	param.synthesis.3 = par.int[3]; 
	param.synthesis.4 = par.int[3]; 
	s = compute.sigmoid(t = zt, fold.change = param.synthesis.1, phase = param.synthesis.2, up.time = param.synthesis.3, down.time =  param.synthesis.4)
	err = sum((log(S)-log(s))^2);
	#err = sum((S-s)^2)
	#print(s[1:12]);
	#print(s[13:24]);
	#print(err);
	#print('>>>>>>>>>>')
	return(err)

}

f2min = function(par, M, S, prior.factor.gamma= parior.factor.gamma, model, parametrization =c('cosine','sigmoid'), method.intern =  c('integration','simulation'), debug = FALSE, underdetermination.model3=FALSE, sum.species = FALSE,  zt = seq(0,46,by = 2), absolute.signal = FALSE)
{
	parametrization = parametrization[1]; 
	method = method.intern[1]
	#print('here')
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
	
	## the same uptime and downtime
	if((model == 2)|(model == 4)){j = model; param.synthesis.1 = par[j]; param.synthesis.2 = par[j+1]; param.synthesis.3 = par[j+2]; param.synthesis.4 = par[j+2]}
	#if((model == 2)|(model == 4)){j = model; param.synthesis.1 = par[j]; param.synthesis.2 = par[j+1]; param.synthesis.3 = par[j+2]; param.synthesis.4 = par[j+3]}
	if(model==3 & !sum.species & underdetermination.model3) {gamma = prior.factor.gamma/(1+eps.gamma);}
	
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
								  fold.change = param.synthesis.1, phase = param.synthesis.2, up.time = param.synthesis.3, down.time = param.synthesis.4, 
								  rescale = rescale, synthesis.factor = syn.f)
			
			if(length(zt.for.sigmoid)!=length(zt)){m = rep(m,max(zt)%/%24+1)}
			
			if(sum.species & !absolute.signal){m = (m+s)/mean(m+s)}else if(absolute.signal){m = m+s}
			
		}
	}
	
	#err.fit = errors(S = S,s = s,M = M,m = m, param = par, model = model, parametrization = parametrization, log = TRUE); 
	err.fit = error(S=S, s=s, M=M, m=m, log=TRUE)
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
		
		#mean = mean(profile) ## this is wrong !!!
		#print(c(mean, mean(profile)));
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
	
	x = 1+(fold.change-1)*(-1+1/(1+exp(-(t-phase + 24 +4* up.time)/up.time)) + 1/(1+exp((t-phase + 24 -4* down.time)/down.time)) 
						   -1+1/(1+exp(-(t-phase +4* up.time)/up.time)) + 1/(1+exp((t-phase -4* down.time)/down.time)) 
						   -1+1/(1+exp(-(t-phase - 24 +4* up.time)/up.time)) + 1/(1+exp((t-phase - 24 -4* down.time)/down.time)) 
						   -1+1/(1+exp(-(t-phase - 48 +4* up.time)/up.time)) + 1/(1+exp((t-phase - 48 -4* down.time)/down.time))
						   -1+1/(1+exp(-(t-phase + 48 +4* up.time)/up.time)) + 1/(1+exp((t-phase + 48 -4* down.time)/down.time))
						   -1+1/(1+exp(-(t-phase - 72 +4* up.time)/up.time)) + 1/(1+exp((t-phase - 72 -4* down.time)/down.time))
						   -1+1/(1+exp(-(t-phase + 72 +4* up.time)/up.time)) + 1/(1+exp((t-phase + 72 -4* down.time)/down.time))
						   -1+1/(1+exp(-(t-phase - 96 +4* up.time)/up.time)) + 1/(1+exp((t-phase - 96 -4* down.time)/down.time))
						   -1+1/(1+exp(-(t-phase + 96 +4* up.time)/up.time)) + 1/(1+exp((t-phase + 96 -4* down.time)/down.time))
						   )
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
				
				if(any((m == 'Inf')|(m == 'NaN')|is.na(m))){break}	
			}
		}}
	
	if((any((m == 'Inf')|(m == 'NaN')|is.na(m)|(m<10^-100))) & simulate)
	{ 
		m = simulate.m(t = t, par = c(gamma, eps.gamma, phase.gamma, fold.change, phase, up.time, down.time), parametrization = 'sigmoid');
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


##################################################################################################################
################################################################################################################## finishing line of optimization
##################################################################################################################

######################
# Functions for model selecion
######################
set.nb.param = function()
{
	n.param <<- c(0,4,3,6);
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
	
	
	#colnames(LRP) = paste(colnames(LRP), '.Laura', sep='')
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
		#colnames(AIC) = paste(colnames(AIC), '.Laura', sep='');
		
	}else{
		colnames(AIC)[5] = 'AIC.best.model';
		#colnames(AIC) = paste(colnames(AIC), '.Laura', sep='');
	}	
	
	return(AIC)								
}

BIC.LAURA = function(T = T, fact = 2)
{
	set.nb.param();
	nb.data = 48
	BIC.m1 = log(nb.data)*n.param[1] + nb.data*log(T$error.m1)
	BIC.m2 = log(nb.data)*n.param[2] + nb.data*log(T$error.m2)
	BIC.m3 = log(nb.data)*n.param[3] + nb.data*log(T$error.m3)
	BIC.m4 = log(nb.data)*n.param[4] + nb.data*log(T$error.m4)
	
	BIC = data.frame(BIC.m1, BIC.m2, BIC.m3, BIC.m4, stringsAsFactors = FALSE)
	BIC$best.model = as.numeric(apply(BIC,1,which.min))
	
	colnames(BIC)[5] = 'BIC.best.model'
	#colnames(BIC) = paste(colnames(BIC), '.Laura', sep='')
	
	return(BIC)								
}	

################################################################################### finishing line for the model selection methods

## clean the results of model selection:
cleaning.model.selection.results = function(T=Tt, parametrization =c('sigmoid','cosine'),absolute.signal = FALSE, model=4)
{
## T=Tt; absolute.signal = FALSE; model=4; i = 1;Nfit = NA; debug = FALSE; parametrization =c('sigmoid');method =  c('integration','simulation'); ZT.int = grep('.rel.ampl.int', colnames(T));ZT.ex = grep('.rel.ampl.ex', colnames(T));	zt = seq(0,46,by = 2)
	
	parametrization = parametrization[1];
	
	bounds = set.bounds(model = model, parametrization = parametrization, absolute.signal = absolute.signal); 
	upper = bounds$upper; 
	lower = bounds$lower;
	
	model = c(2,3,4)
	
	LRT.best.model = T$LRT.best.model
	AIC.best.model = T$AIC.best.model
	AICc.best.model =T$AICc.best.model
	BIC.best.model = T$BIC.best.model
	
	cutoff = 0.000001
	
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
			if(!is.na(BIC.best.model[n])) 
			{
				if(BIC.best.model[n]==2) BIC.best.model[n]=NA
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
			if(!is.na(BIC.best.model[n])) 
			{
				if(BIC.best.model[n]==3) BIC.best.model[n]=NA
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
			if(!is.na(BIC.best.model[n])) 
			{
				if(BIC.best.model[n]==4) BIC.best.model[n]=NA
			}
		}
		
	}
	T$LRT.best.model.eff = LRT.best.model
	T$AIC.best.model.eff = AIC.best.model
	T$AICc.best.model.eff = AICc.best.model
	T$BIC.best.model.eff = BIC.best.model
	return(T)
	
	
}

### using Z score to make some selection when extracting estimated parameters

plot.figure.fit.results.zscore = function(T = T, best.model = best.model, gamma = gamma, gamma.zscore = gamma.zscore, gamma.fc = gamma.fc, gamma.fc.zscore = gamma.fc.zscore, phase.gamma = phase.gamma,phase.gamma.score = phase.gamma.score, phase.int = phase.int, phase.int.score = phase.int.score, 
up.time.int = up.time.int, up.time.int.score = up.time.int.score, P = P)
{
	#T =T; best.model = best.model; gamma = gamma; gamma.fc = gamma.fc; phase.gamma = phase.gamma; phase.int = phase.int; P = P;
	
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
	for(r in 1:14)
	{
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
	
		
	bounds = set.bounds(model = 4, parametrization = 'sigmoid')
	
	# half-lives 
	j2 = best.model == 2
	j3 = best.model == 3
	j4 = best.model == 4

	#j2 = which(best.model == 2 & log(2)/T$gamma.m2^2*T$gamma.stderr.m2<2)
	#j3 = which(best.model == 3 & log(2)/T$gamma.m3^2*T$gamma.stderr.m3<2)
	#j4 = which(best.model == 4 & log(2)/T$gamma.m4^2*T$gamma.stderr.m4<2)
	
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


plot.figure.fit.results = function(T = T, best.model = T$model.MCMC, gamma = T$gamma.MCMC, gamma.fc = T$eps.gamma.MCMC, phase.gamma = T$phase.gamma.MCMC, phase.int = T$phase.int.MCMC, up.time.int = up.time.int, P = P)
{
	#T =T; best.model = best.model; gamma = gamma; gamma.fc = gamma.fc; phase.gamma = phase.gamma; phase.int = phase.int; P = P;
	
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
	for(r in 1:14)
	{
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
	
	
	bounds = set.bounds(model = 4, parametrization = 'sigmoid')
	
	cutoff = 1.036
	
	# half-lives 
	j2 = best.model == 2
	j3 = best.model == 3
	j4 = best.model == 4
	
	#j2 = which(best.model == 2 & log(2)/T$gamma.m2^2*T$gamma.stderr.m2<2)
	#j3 = which(best.model == 3 & log(2)/T$gamma.m3^2*T$gamma.stderr.m3<2)
	#j4 = which(best.model == 4 & log(2)/T$gamma.m4^2*T$gamma.stderr.m4<2)
	
	
	breaks = lseq(log(2)/max(gamma,na.rm = TRUE),log(2)/min(gamma,na.rm = TRUE), len = 31)
	
	par(mgp = c(1.6,0.5,0), mar = c(3,2,2,0))
	
	h2 = hist(log(2)/gamma[j2[which(gamma.zscore[j2]>cutoff)]], breaks = breaks, plot = FALSE)
	h3 = hist(log(2)/gamma[j3[which(gamma.zscore[j3]>cutoff)]], breaks = breaks, plot = FALSE)
	h4 = hist(log(2)/gamma[j4[which(gamma.zscore[j4]>cutoff)]], breaks = breaks, plot = FALSE)
	
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





