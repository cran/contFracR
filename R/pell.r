# pell  x^2 - d * y^2 = 1 , d a non-square integer
#wiki info:
# Let   h_i / k_i denote the sequence of convergents to the regular continued fraction for sqrt(n). This sequence is unique. Then the pair  solving Pell's equation and minimizing x satisfies x1 = hi and y1 = ki for some i. This pair is called the fundamental solution. Thus, the fundamental solution may be found by performing the continued fraction expansion and testing each successive convergent until a solution to Pell's equation is found.
# h_i are nums, k_i denoms of each successive layer in contfrac  

# test example  for d ==7
# The sequence of convergents for the square root of seven are
# h/k (convergent)	h^2   7k^2 (Pell-type approximation)
# 2/1			 3
# 3/1			+2
# 5/2			 3
# 8/3			+1
# so 8 / 3 is the one we want


# rev after 1.2: more input validation, done better
#  generate more convergents if no "win" found.  
#  allow A,B inputs to continue  a previous run
#  make sure $repeated == TRUE before running with the results of
#  sqrt2periodicCfrac

pell <- function (d , maxrep = 10, A = NULL, B = NULL){
	# 	wins <- NULL  # NO!  gmp really hates NULL things
	wins <- as.bigz( matrix(0, nrow= 1 ,ncol = 2)) 
	d <- as.bigz(d[1])
	# get the greatest square less than x; make sure mpfr has reasonable precision
	dm  <- .bigz2mpfr(d, max(5 * ceiling(log10(d)), 500))
	#validate - check that sqrt(d) is not an integer
	if (floor(sqrt(dm)) == sqrt(dm)) {
		stop('Pell is pointless when d is a square integer')
	}
	foo <- sqrt2periodicCfrac(d)
	# check for success
	nterm = 100  # default is 50
	while (!foo$repeated) {
		foo <- sqrt2periodicCfrac(d, nterms = nterm)
		nterm = nterm *2 
	}
	repden <- foo$denom[-1]
	replen <- length(repden)
	jrep = 0
	if (is.null(A)) {
	# build initial stuff
		# conv{A,B}[1] is the "starter" value A[-1], so ignore
		for (jc in 2:length(foo$convA)) {
			if (foo$convA[jc]^2 -d * foo$convB[jc]^2 -1 == 0) {
				wins <-	rbind(wins,c(foo$convA[jc],foo$convB[jc]) )
				}
		}
		jconv = jc+1
		A = foo$convA
		B = foo$convB
		# end of if isnull convA
	}  else {
		A <- A
		B <- B
		jconv = length(A) + 1
	}
#   (wins < 4 because of dummy first row)

	while(length(wins ) < 4 && jconv <= maxrep * replen) {
	# need more convergents. use the repeated section of foo$denom	
		A[jconv] <- repden[jrep + 1] *A[jconv-1] + 1 * A[jconv-2]
		B[jconv] <- repden[jrep + 1] * B[jconv-1] + 1 * B[jconv-2]
		if (A[jconv]^2 - d * B[jconv]^2 -1 == 0) wins <- rbind(wins,c(A[jconv],B[jconv]) )
		# increment jrep mod(length(repden))
		jrep <- (jrep+1) %% replen
		jconv = jconv + 1 
	}
	return(invisible(list(d = d, wins=wins[-1,], cfracdata = foo, A= A, B= B)))
}