#  algebraic calculation of square root simple continued fracs
# This function finds the periodic expansion. Use root2cfrac for nonperiodic but
# faster convergence. 
# input must be int or bigz, or a charstring which will be floor-ed to bigz
# Default limit on nterms just to avoid monstrous output
#  sqrt(n/m) =  sqrt(n*m)/m ,so can do this for any rational number. 
# irrationals we can only approximate anyway. 
# reference:
#  rt(67) is 8;   5  2  1  1  7  1  1  2  5 16
# 
#output definition:
#	repeated: TRUE if repeat sequence found. FALSE if exceed nterms w/o repeat
#	input: echo the inputs

sqrt2periodicCfrac <- function (num, denom = 1,  nterms = 50, ...){
if ( length(num) >1 || length(denom) >1 ){
	warning('Only first element of num or denom will be used')
	num <- num[1]
	denom <- denom[1]
	}
input <- c(num,denom) #save to dump into output
#convert to a single bigz and save the denom to re-apply at end of calculation
#first, reduce to minimum fraction.
frac <- abs(as.bigz(num) /as.bigz(denom))
num <- numerator(frac) * denominator(frac)
# This is a final multiplier to get the desired sqrt value
thedivisor <- denominator(frac)
# this will store the cont frac denoms 
#TODO: replact with initialization thedenom <- as.bigz(NULL)
# and see about simplifying the converstion to numerics
thedenom <- as.bigz(NULL)  #list()
numericdenom <- vector()
#follow math.stackexchange question 2215918
# get the greatest square less than x; make sure mpfr has reasonable precision
kbits <- max(5 * ceiling(log10(num)), 500)  
mpnum <- mpfr(num,kbits)
#start the engine
m <- floor(sqrt(mpnum))  # the first term in standard contfrac notation
thedenom <- c(thedenom, .mpfr2bigz(m) )#this is the "integer" part of contfrac form  
pq <- c(m,1)  # this is mpfrs
newq <- (num - pq[1]^2) / pq[2] 
#if newq is zero, terminate, perfect square was input
if (newq == 0 ) {
	return(invisible(list(thefrac = thedenom, numericfrac = numericdenom,  thedivisor=thedivisor, repeated='perfect_square', input=input) ))
}
mpdenom <- floor( (pq[1] + m)/newq)
thedenom <- c(thedenom, .mpfr2bigz( (pq[1] + m)/newq) )
# numericdenom[2] <- as.integer(thedenom[[2]])
newp <-  mpdenom * newq - pq[[1] ]
pq <- c(newp,newq)
#loop on that for either nterms or until  mpdenom == 2*m 
# Since we can only feed rationals to the function, and we convert
# to an integer value, there can't be a 'preamble' 
idx <- 2
#initialize
repeated = FALSE # indicating no repeat sequence found
while ( (idx <= nterms) && !repeated  ) {
	newq <- (num - pq[[1]]^2) / pq[[2]] 
	mpdenom <- floor( (pq[1] + m)/newq)
	newp <-  mpdenom *  newq - pq[[1]]  
	pq <- c(newp,newq)
	thedenom <- c(thedenom, .mpfr2bigz(mpdenom) )	# watch for giant integers... 
	if (mpdenom == 2*m) {repeated = TRUE}
	idx <- idx + 1
}
# now stick that 1/m term back in
thedenom[1] <- thedenom[1]/thedivisor
numericdenom <- as.numeric(thedenom)
numericnums <- vector()
numericnums[1] <- 1/as.integer(thedivisor)
numericnums[2:(length(numericdenom) -1)] <- 1
thenum <- as.bigq(numericnums)
thenum[1] <- 1/thedivisor  #removes possible precision foulup 
return(invisible(list(denom = thedenom, num = thenum, numericdenom = numericdenom, numericnum =numericnums,  repeated=repeated, input=input) ))

}

