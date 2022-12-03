# requires go2bigq; I do want that to be a separate package, though. 
# Function  num2cfrac

#Inputs can be numeric, bigz, or reasonable character strings representing numbers. Or feed a bigq to "numerator"  
# char strings can include '-' or '+' leading sign and 'eX' exponential notation
#  you can calculate irrationals to arbitrary precision by setting the inputs to  num = (irat_num)*1eN , denom = 1eN for N some nice integer
# mpfr values not allowed because teensy rounding errors, e.g. "3.2400000000001" make  a hash of the fraction

#TODO revision Nov 2022: think about any other possible input validation to add
#	 can I calc convergents along the way?  
# to get convergents WHILE calculating the cfrac, use recursion formulas. bj are denoms with b0 the integer part,  aj are nums.

# "starter values A[-1] = 1, B[-1] = 0

#  A0 = b0 
#  B0 = 1
# Aj = bj*A(j-1) + aj*A[j-2]
# Bj =  bj*B[n-1] + aj*B[n-2] 
#
#
num2cfrac <-function(num, denom = 1, ...) {
if(length(num) >1 ) {
	warning('only using first element of num')
	num <- num[1]
	}
if(length(denom) >1) {
	warning('only using first element of denom')
	denom <- denom[1]
	}
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
# do go2bigq on both and 'flip' the division in gmp space
if(!is.bigz(num)) num <- go2bigq(num)
if(!is.bigz(denom)) denom <- go2bigq(denom)
# ok, now can reduce from bigz or bigq form. In addition, this step
# produces a bigq in reduced form. 
thebigq <- num[[1]]/denom[[1]]
numdenom  <- c( abs(numerator(thebigq)), denominator(thebigq))
thesign <- sign(thebigq)
ND <- numdenom 
if (numdenom[2] == 0 || numdenom[1] == 0 ) stop('zeros not allowed')
if (numdenom[1] == 1) return (invisible(list(intvec = c(0,denom), numdenom = ND) ) )
# intvec = as.bigz(NULL)  #NOT a R-vector, but NULL will be "replaced" with first value assigned.  
jlev = 2   # start vector index, with offset for A[-1] stuff

# initialize convergents. Notice that first elements are the "minus 1" indices
# since all nums are 1 at this point, can substitute that in for a[j]
#  browser()
 intvec <- bzero <- floor(numdenom[1]/numdenom[2] )
 
	num = numdenom[1] - (intvec * numdenom[2] )
	numdenom = c(numdenom[2],num) # reciprocalling
A <- c( as.bigz(1) , bzero )  #otherwise coercion screws up bigz value 
B <- as.bigz(c(0 ,1)) 
# numerator == 0 means there was a common factor; this is easiest way to catch it.
# rats: is.element, %in% don't like bigq
#while (!is.element(numdenom[1], c(0,1)) && numdenom[2] != 0 ) {
while(numdenom[1] !=1 && numdenom[1] != 0 && numdenom[2]!= 0){
	intvec = c(intvec,  floor(numdenom[1]/numdenom[2] ) )
	num = numdenom[1] - (intvec[jlev] * numdenom[2] )
	numdenom = c(numdenom[2],num) # reciprocalling
# numdenom is "internal" calc. The next denom is latest intvec
	A[jlev + 1] <- intvec[jlev] * A[jlev ] +  A[jlev -1]
	B[jlev + 1] <- intvec[jlev] * B[jlev ] +  B[jlev -1]
	 jlev = jlev + 1
	}
intvec <- intvec * thesign
return(invisible(list(denom = intvec , numdenom = ND ,convA = A, convB = B )) )
}


# FUNCTION cf2latex -- build the latex formula and the inline equation
# the input 'denomvals' is the 'intvec' output of num2cfrac, or any other legal Continued Fraction form vector of integers. Similar for numvals if setting up a non-simple continued fraction. 
# Beware: R interprets backslashes in a char string as escapers in most cases. Example:
 # foo <- "this string has \f backslashes \\ and stuff"
# > unlist(strsplit(foo,''))
 # [1] "t"  "h"  "i"  "s"  " "  "s"  "t"  "r"  "i"  "n"  "g"  " "  "h"  "a"  "s" 
# [16] " "  "\f" " "  "b"  "a"  "c"  "k"  "s"  "l"  "a"  "s"  "h"  "e"  "s"  " " 
# [31] "\\" " "  "a"  "n"  "d"  " "  "s"  "t"  "u"  "f"  "f"   

# So strings must be submitted in  double-escaped format
 
cfrac2latex <-function(denomvals, numvals = 1, denrepeats = 1, doellipsis = FALSE,  ...) {
vrep <- denomvals[-1] # don't repeat the lead integer
denomvals <- c(denomvals[1], rep(vrep, times= denrepeats) )
# Do some legerdemain to get length(numvals) to match length(denomvals)
# and some new variable PM , pm
numvals <- rep(numvals, length = length(denomvals)-1)
PM = c(' - ',' + ',' + ')
pm <- PM[sign(numvals) + 2] 
# now remove signs from numeric numvals
numvals <- abs(numvals)
eqn <- paste0(denomvals[1], pm[1])
#  stop the loop one short to avoid extra "+"
for (jj in 2:(length(denomvals) -1) ) {
#first denomval is integer, hence teh index shift for nums
	eqn <- paste0(eqn, numvals[jj-1],'/(',denomvals[jj], pm[jj], collapse='')
	}
#now the last integer, no plus sign
if (doellipsis){
	eqn <- paste0(eqn, numvals[jj],'/', denomvals[jj+1],' ...',collapse='')
	}else{
		eqn <- paste0(eqn, numvals[jj],'/', denomvals[jj+1],collapse='')
	}
#finish the inline version 
clospar <-paste0(rep(')',jj-1),collapse='')
eqn <- paste0(eqn, clospar,collapse='' )
# Now a Q&D conversion to LaTeX, skipping markdown "$"
# tex2copy is useful only for copy/paste, as the \ will be parsed by console 
# replace '1' with num value.  ([0-9]{1,}) 
tex2copy <- gsub('([0-9]{1,})/[(]', '\\\frac{\\1}{' , eqn)
tex2copy <- gsub('[)]', '}', tex2copy)
#cleanup
tex2copy <- gsub('([0-9]{1,})/', '\\\frac{\\1}{', tex2copy)
tex2copy <- paste0(tex2copy,'}',collapse='')
# convert the single actual character "\f" to two characters "\\" and "f" 
texeqn <- gsub('\\f', '\\\\\\f', tex2copy)
return(invisible(list(eqn=eqn, tex2copy = tex2copy, texeqn = texeqn)) )
}


# calculating numeric value of a continued fraction representation ,
# where K is either the Continued Fraction form sequence of integers, or a single value
# if the particular source is that value repeated 'numterms' times. 
# Input 'nreps' is how many times to repeat a single-value denom input
#
# outputs
#	cgmp is the gmp num&denom
#	cmpfr is the mpfr conversion of cgmp
#	nums, denoms returns the input vectors used

# TODO: check for list inputs and either ban them or unlist(carefully)
# Revision TODO Nov 2022.  Per Hans Borcher,  
# cfrac2num(c(0, 1,-1,1)); x
# Error in `/.bigq`(num[[kk]], spart) : division by zero

#  One can certainly argue whether negative entries shall be allowed for continued
#  fractions, but I think it will be better if the function will check the input instead of 
#  the system.

#  store every convergent , not just the final one (cgmp)  for people to play with  

#BUGFIX: nreps fouls up when length of denom is bigger, as it will be for a periodic 
#  contfrac resulting from a square root, for example. 
cfrac2num  <- function(denom, num = 1, init = 0, nreps= 1, msgrate = NULL ) {
if (length(denom) == 1){
#  "exit strategy" when denoms == 0.  
# Zero is only legal in the first position of  the contfrac sequence anyway. 
        if( denom == 0 ) {
                return(invisible(list(cgmp = as.bigq(0), mdec=mpfr(0,10),denom = denom, num = num)))
                }
# bug-ish in tail(gmp), so check nterms value
        denom <- rep(denom, max(nreps,1) )
# One More Thing(TM) to check: if all we have is the 'a' term:
        if (length(denom) == 1) {
                return(invisible(list(cgmp = as.bigq(denom), mdec=mpfr(denom,10),denom = denom, num = num)))
                }
        } else {
      # Now apply nreps to a putatively periodic denominator
   			denom <- c(denom[1], rep(denom[-1],times=nreps))
        }

#  expand num to be  length(denom) -1  if it isn't already (first denom is the integer)
#  and if length(num) == 1 , repeat it by nterms-1
if( !is(num[[1]],'bigq')){
	num <- as.bigq(rep(num,length=length(denom)-1 ) )
}else{
	num <- rep(num,length = length(denom)-1)
}
if (!is(denom[[1]],'bigq')) {
	Kq <- as.bigq(denom)
}else{
	Kq <- denom
}
# initalize
spart <- ( tail(Kq, 1) )[[1]] + init 
# new: save all convergents
#  NOPE:  convergents have to be calc'd from start, not end 
#  conv <- spart
# add a status msg. Also notice reverse order as befits unwinding a cont frac
loopvals <- (length(denom) - 1):1
if (!is.null(msgrate)){
	msgrate <- floor(msgrate)
	msgnums <- loopvals[!(loopvals %% msgrate)]
	for (kk in loopvals) {
		if ( kk %in% msgnums) message('kk is now ', kk)
        spart <- Kq[kk] + num[kk] /spart  #1 / spart
 #       conv <- c(conv, spart)
        }
	
	} else {
		for (kk in loopvals) {
			spart <- Kq[[kk]] + num[[kk]] /spart  #1 / spart
#		conv <- c(conv, spart)
        }
		}
	
mdec <- .bigq2mpfr(spart) # no reason to apply nonoptimal nbr bits
fracpart <- spart - Kq[[kk]]
return(invisible(list(cgmp = spart, cmpfr = mdec, fracpart = fracpart, denom = denom, num = num)))
}
