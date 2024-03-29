# wikipedia  suggests that making x and y big (  nth root of z^m  with arbitrary z  = x^n + y . Maybe choose largest possible x? In simplest
# setup,let m == 1  . Then find the largest n-th cube less than  z.   
# lead term is a^m 
# nums are my, (n-m)b, (n+m)b, (2n-m)b, (2n+m)b, 3,3,4,4,.....
# denoms are  1na^(n-m), 2a^m, 3na^(n-m), 2a^m, 5na^(n-m), alt 2term and 7,9,11,...term

# Calculating  x^(m/n)  

# Nov 2022: improve input validation work.  

root2cfrac <- function(x, n, m = 1, nterms = 10, ...){	
#sanitize
n <- floor(n[1])
m <- floor(m[1])
#TODO: check to see if can replace x, m with  some root K(x) and m*K . Probably 
# doesn't guarantee faster convergence. 
if(length(x) > 1) warning('only first element of x will be used') 
x <- as.bigq(x[1]) #(just in case x is not an integer)
if (x < 0 ) {
	x <- abs(x)
	warning('negative x not allowed. Taking abs value')
}
if ( m * n < 0)  {
	stop('negative exponents not allowed. Consider entering 1/x in bigq form')
}
if (n == 0) stop('infinite exponent not allowed (n > 0 required)')
# want  a^n + b = x 
# gmp doesn't allow fraction powers (but .bigq2mpfr "guarantees" sufficient digits are used to be exact)
xm <- .bigq2mpfr(x)
a <- go2bigq(floor( xm^(1/n)))[[1]]
b = x - pow.bigq(a, n )
denom <- as.bigz(NULL)  
num <- as.bigz(NULL)  
numericdenom <- vector()
numericnum <- vector()
denom <- c(denom,pow.bigq(a, m))
num <- c(num,m * b)
jnum = 1  #indexer
for (jt in seq(1,nterms, by = 2)) {
#that skipseq works  but take care
	num <- c(num, (jnum *n-m)*b,(jnum *n + m)*b )
	denom <- c(denom,jt*n* pow.bigq(a, (n-m)), 2*pow.bigq(a,m) )
# bump num index
	jnum = jnum + 1
}
#finish off denom  but watch for oddeven jt as it always ends with odd regardless of nterms
denom <- c(denom, (jt+3) *n* pow.bigq(a, (n-m)))
numericdenom <- as.numeric(denom)
numericnum <- as.numeric(num)
return(invisible(list(num = num, denom=denom, numericnum = numericnum, numericdenom = numericdenom, x=x, n= n, m= m)) ) 
}