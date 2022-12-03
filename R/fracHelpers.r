# simple helpers
#  see https://math.stackexchange.com/questions/86043/
#   move original b1 inside the a1....an sequence and setting the new b1 == 0
# and for numerator, prepend a "1" 
# generate terms for 1/X 

#TODO  even tho' these are primarily helper funcs, add input validation 

recipCfrac <- function( denom, num = 1, ...){
# ensure proper lengths
num <- rep(num,length=length(denom) -1)
denom <- c(0,denom)
num <- c(1,num)
return(invisible( list(denom=denom, num = num)) )
}

# generate equivalent simple denom for any num,denom vectors
cfrac2simple <- function(denom, num=1, ...){
foo <- cfrac2num(denom = denom, num = num)$cgmp
thenum <- numerator(foo)
thedenom <- denominator(foo)
denom <- num2cfrac(thenum, thedenom)$denom
return(invisible(list(intfrac = foo, denom = denom)) )
}
