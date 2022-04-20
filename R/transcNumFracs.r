#  formulas for contfrac of "special" numbers
#brouncker 4/pi = 1 +  numseq{1,inf}((2k-1)^2);denomseq(2)
#stern pi/2 = 1 - numseq(1, 2*3,1*2,4*5, 3*4, 6*7, ....); denomseq3,1,3,1,3,1...); make sure all denoms are N - (term)
# 2008, Coleman and Pickett
#pi/2 =  {1; series(2,N,1/j)}  i.e  1 +  denoms 1, 1/2, 1/3, 1/4 .....

pi2cfrac <- function(nterms, method = c('brouncker','stern','coleman'),...) {
method <- tolower(method[1])
goody <- c('brouncker','stern', 'coleman')
method <- goody[startsWith(goody,'col')]
if(!length(method)) stop ('unknown method %s ', method)
switch(method,
	'brouncker' =
	{
# change these so it's already reciprocal and multiply first num  by 4 (4/pi = x --> pi = 4/x) and stick that zero in denom. 
# 4/pi is this:
	num <- (2*(1:nterms)-1)^2
	denom <- c(1, rep(2,nterms))
# so pi is this:
	num <- c(4,num)
	denom <- c(0,denom)
	},
	'stern' = 
	{
	jj = 1:nterms
# just make first num == 2 and first denom == 2 
	num <- -c(2,rbind(((2*jj)*(2*jj+1)), (2*jj)*(2*jj-1))) 
	denom <- c(2,rep(c(3,1),length=nterms-1))
	},
	'coleman' = 
	{
	nterms = 1:nterms
#just make first num == 2 and first denom == 2
	denom = c(as.bigz(2), as.bigq(rep(1,times = length(nterms)), nterms) )
	num = c(2,rep(1,times = length(denom)-2 ))
	},
	stop('unknown method %s ',method) #redundant but keep for later 
)
return(invisible(list( nterms=nterms, method=method, num=num,denom=denom )) )
}


# Euler number 'e'   
# simplest one is Euler  2; 1,2,1,1,4,1,1,6,1,1,8,...
# Wall, 1948:   2 + numseq(1,1,2,3,4...) ; denomseq(1,2,3,4,...) 
# Euler, or possibly Gauss via hyperbolic functions: 
# e^(x/y) = 1 +numseq(2x,x^2,x^2,x^2,...) ; denomseq(2y-x,6y,10y,14y, 18y, 22y,..)
# which, annoyingly, for y == x == 1, gives denoms 2,6,10,14,.. because the nums aren't all 1

e2cfrac <- function(nterms, pownum = 1, powdenom = 1,  method = c('euler','wall','gauss'), ...){
method <- tolower(method[1])
goody <- c('euler','wall','gauss')
method <- goody[startsWith(goody,method)]
if(!length(method)) stop ('unknown method ', method)
pownum <- floor(pownum[1])
powdenom <- floor(powdenom[1])
# if num/denom != 1, ignore "method" and use the euler/gauss formula
if (pownum/powdenom != 1 && !(method == 'gauss')) {
	warning('exponent is not 1; using Euler-Gauss formula for e^(pownum/powdenom)')
	method = 'gauss'
	}
#turn into switch TODO
if (method == 'gauss') {
	num = c(2*pownum, rep(pownum^2, times = nterms-2))
	denom = c(1, 2*powdenom-pownum, powdenom*(2 + 4*(1:(nterms-2)))  )
	method = 'EulerGauss'
	} else {
if (method == 'euler'){
	num = 1
	jd = ceiling((nterms-2)/3)
	tmpd <- rbind(2*(1:jd),rep(1,times=jd),rep(1,times=jd))
	denom = c(2,1,tmpd)[1:nterms]
	} else {
if (method =='wall'){
	num <- c(1,1:(nterms-2))
	denom <- c(2,1:(nterms-1))
	} else stop('unknown method %s ', method)
		}
	} #end of else 
return(invisible(list( nterms = nterms, method = method, num=num, denom = denom, pownum = pownum, powdenom = powdenom)))
}

#
# phi, (1+sqrt(5))/2,  is ,cutely enough  1; 1,1,1,1,1,1
# phi ^3 is  4; 4,4,4,4,4  and in fact ...
# The Lucas numbers (Ln) are a sequence defined by the recurrence Ln+2 =
#Ln+1 + Ln, where L0 = 2 and L1 = 1.
# and , proved in various places, for n odd,  phi^n = Ln; Ln,Ln...  and
# for n even,  phi^n = Ln-1 ; 1,Ln-2,1,Ln-2,...

phi2cfrac <- function( nterms = 10, exponent = 1,  ... ) {
exponent <- floor(exponent[1])
nterms <- floor(nterms[1])
lucas <- vector()
lucas[1:2] <- c(2,1)
for (jn in 3:(max(3,exponent+1))){
	lucas[jn] <- lucas[jn-1] + lucas[jn-2]
	}
# be warned - standard notation starts with lucas[0], so shift things here
lucval <- lucas[exponent+1]
if (exponent%%2) {
	phidenom <- rep(lucval,times=nterms)
	}else{  
	# it's an even power
		phidenom <- c(lucval-1, rep(c(1,lucval-2),times = ceiling(nterms/2)))
	}
return(invisible(list(denom = phidenom, nterms = nterms, exponent = exponent)))
}