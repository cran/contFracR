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

# new Nov 2022:   input validation, 
#   add a warning if num or denom is not integer, as I do force them to be bigz
#  calculate convergents as in latest num2cfrac

sqrt2periodicCfrac <- function (num, denom = 1,  nterms = 50, ...){
if ( length(num) >1 || length(denom) >1 ){
	warning('Only first element of num or denom will be used')
	num <- num[1]
	denom <- denom[1]
	}
# input manipulation as in num2cfrac
# new sortorific bit of code to catch bigq
if (is(num,'bigq') ) {
	denom <- denominator(num)
	num <- numerator(num)
}
numclass <- class(num)
denomclass <- class(denom) 
classlist <- c('numeric','character', 'bigz')
if (!(numclass %in% classlist) || !(denomclass %in% classlist)){
	stop('Inputs must be bigz,bigq, numeric, or character')
	}
# do go2bigq on both 
if(!is.bigz(num)) num <- numerator(go2bigq(num)[[1]])
if(!is.bigz(denom)) denom <- numerator(go2bigq(denom)[[1]])
	
#	browser()
if (floor(num) != num || floor(denom) != denom) {
		warning('num, denom must be integers. Will truncate')
		num <- floor(num)
		denom <- floor(denom)
	}
if (num * denom <  0 ) {
	num <- abs(num)
	denom <- abs(denom)
	warning(' negative values not allowed; taking abs val')
}
if ( denom == 0) {stop(' denom of zero disallowed')}
input <- c(num,denom) #save to dump into output
if (nterms < 1) {stop('nterms must be a positive integer')}
nterms <- floor(nterms)
#convert to a single bigz and save the denom to re-apply at end of calculation
#first, reduce to minimum fraction.
frac <- abs(as.bigz(num) /as.bigz(denom))
#  this means we need to report the actual value is the returned value divided by 
# the denominator provided as input 
num <- numerator(frac) * denominator(frac)
# This is a final multiplier to get the desired sqrt value
thedivisor <- denominator(frac)
# this will store the cont frac denoms 
thedenom <- as.bigz(NULL)  
numericdenom <- vector()
#follow math.stackexchange question 2215918
# get the greatest square less than x; make sure mpfr has reasonable precision
kbits <- max(5 * ceiling(log10(num)), 500)  
mpnum <- mpfr(num,kbits)
#start the engine
m <- floor(sqrt(mpnum))  # the first term in standard contfrac notation, aka b[0]
mz <- .mpfr2bigz(m)
# initialize convergents. Notice that first elements are the "minus 1" indices
# since all nums are 1 at this point, can substitute that in for a[j]
A <- c(as.bigz(1), mz )  #
B <- as.bigz(c(0 ,1)) 

thedenom <- c(thedenom, mz )#this is the "integer" part of contfrac form  
pq <- c(m,1)  # this is mpfrs pair (p,q)
#  do I really want the (numerator * denominator) yes
  newq <- (num - pq[1]^2) / pq[2] 
#newq <- (m - pq[1]^2) / pq[2] 
#if newq is zero, terminate, perfect square was input
if (newq == 0 ) {
	return(invisible(list(thefrac = thedenom, numericfrac = numericdenom,  thedivisor=thedivisor, repeated='perfect_square', input=input) ))
}
mpdenom <- floor( (pq[1] + m)/newq)
thedenom <- c(thedenom, .mpfr2bigz( (pq[1] + m)/newq) )
newp <-  mpdenom * newq - pq[[1]]
pq <- c(newp,newq)
# repeating first num calc from below; clean up later
tmpden <- .mpfr2bigz( (pq[1] + m)/newq)
# this is not right!  Keep the divisor out of my calc'n
# tmpnum <- 1/thedivisor  # i.e.  a[1]  no, it isn't - the formula will have all a[j] = 1
idxab <- 3  # while num[1] may be !=1, it doesn't get called in these formulas

# want a[1] and b[1] here; remember A,B have index off by 2 for [-1,0] indices
# AND that b starts with b0 ; want b[1]    
A[idxab] <- thedenom[idxab-1] *A[idxab-1] + 1 * A[idxab-2]
B[idxab] <- thedenom[idxab-1] * B[idxab-1] + 1 * B[idxab-2]


#loop on that for either nterms or until  mpdenom == 2*m 
# Since we can only feed rationals to the function, and we convert
# to an integer value, there can't be a 'preamble' 
idx <- 2
#initialize
repeated = FALSE # indicating no repeat sequence found

while ( (idx <= nterms) && !repeated  ) {
	idxab <- idxab + 1 
	newq <- (num - pq[[1]]^2) / pq[[2]] 
	mpdenom <- floor( (pq[1] + m)/newq)
	mpdenomz <- .mpfr2bigz(mpdenom)
	newp <-  mpdenom *  newq - pq[[1]]  
	pq <- c(newp,newq)
	thedenom <- c(thedenom, mpdenomz )	# watch for giant integers... 
	if (mpdenom == 2*m) {repeated = TRUE}
	idx <- idx + 1
	# update A, B 
	A[idxab] <- mpdenomz * A[idxab-1] + A[idxab-2]
	B[idxab] <- mpdenomz * B[idxab-1] + B[idxab-2]


}
# 
# thedenom[1] <- thedenom[1]/thedivisor
numericdenom <- as.numeric(thedenom)
numericnums <- vector()
numericnums[1] <- 1  # no, don't /as.integer(thedivisor)
numericnums[2:(length(numericdenom) -1)] <- 1
thenum <- as.bigq(numericnums)
#  thenum[1] <- 1  #  NO /thedivisor   
return(invisible(list(denom = thedenom, num = thenum, numericdenom = numericdenom, numericnum =numericnums,  repeated=repeated, input=input,convA = A, convB = B) ))

}

