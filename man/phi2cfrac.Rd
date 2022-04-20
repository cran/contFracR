\name{phi2cfrac}
\alias{phi2cfrac}

\title{
Function to Calculate Phi And Powers of Phi In Continued Fraction Form.  
}
\description{
This function generates the continued fraction form of the "golden ratio",  \code{phi^N} for integer powers N.   
}
\usage{
phi2cfrac( nterms = 10, exponent = 1,  ...)
}

\arguments{
  \item{ nterms}{
How many denominators to calculate.   
}
  \item{exponent}{	
An positive integer indicating the power of \code{phi} desired. The default is 1.
}
  \item{...}{	
Reserved for future use.
}
 
}
\details{
The 'golden ratio' ,  equal to \code{(1 + sqrt(5))/2}, is the ratio of two sides x < y of a rectangle such that, by removing a square of side x, the remaining rectangle has the same ratio.  

It turns out, in one of those mathematical curiosities, the denominators of the continued fraction form of \code{phi} are all equal to one.  Some people use this to state, humorously, that this makes \code{phi} "the most irrational irrational number."  It also happens that the continued fraction form for powers of \code{phi} consist of Lucas Numbers (see References). 
 
}
\value{
The continued fraction denominators are provided in \code{denom}.  The inputs \code{nterms} and \code{exponent} are echoed back for reference.  
}

\author{
Carl Witthoft,  \email{carl@witthoft.com} 
}
\references{
\url{https://en.wikipedia.org/wiki/Lucas_number}
\url{https://en.wikipedia.org/wiki/Golden_ratio}

}
\seealso{
\code{\link[contFracR]{num2cfrac}}
}

\examples{
phi2cfrac(nterms = 10)
phi2cfrac(exponent = 3)
foop <- phi2cfrac(nterms = 20)
cfrac2num(denom  = foop$denom)
# compare with:
library(Rmpfr)
(1 + sqrt(mpfr(5,1000)))/2 

}
