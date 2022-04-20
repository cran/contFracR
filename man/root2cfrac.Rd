\name{root2cfrac}
\alias{root2cfrac}

\title{
Function To  Generate Continued Fraction For Arbitrary Roots
}
\description{
This function generates the generalized continued fraction  for any input value \code{x^(m/n)}.   
}
\usage{
root2cfrac(x, n, m = 1, nterms = 10, ...)
}

\arguments{
  \item{ x}{
The number itself.  Integers, doubles, bigz, and bigq classes are allowed.   
}
  \item{n}{	
The integer denominator of the power to which x is raised. That is, when \code{m} is 1, the n-th root of x is generated. 
}
  \item{m}{	
The integer numerator of the power to which x is raised. The default is 1. 
}
 
 \item{nterms}{
 How many terms (denominators) to calculate. 
 }
 \item{...}{
 Reserved for future use
 }
}
\details{
The generalized continued fraction for arbitrary roots will not be periodic, and may not even show a pattern in the denominator values. By comparison, \code{\link[contFracR]{sqrt2periodicCfrac}} generates a simple continued fraction with a periodic sequence for square roots only.  That periodic sequence tends to converge more slowly than the aperiodic sequence produced here. 

 
}
\value{
A list, containing:
The continued fraction numerators and denominators in \code{bigz} form \code{num} , \code{denom} .
The continued fraction numerators and denominators in numeric form \code{numericnum}, \code{numericdenom} .  In the extreme case that a value exceeds the machine size of a numeric, \code{NA} is returned. 
The inputs \code{x, n, m} are echoed back. 

}

\references{

\url{https://en.wikipedia.org/wiki/Generalized_continued_fraction }
  

}
\author{
Carl Witthoft,  \email{carl@witthoft.com} 
}

\seealso{
\code{\link[contFracR]{sqrt2periodicCfrac}},  \code{\link[contFracR]{cfrac2num}}
}

\examples{
root2cfrac(x = 2, n = 3)
root2cfrac(x=17, n= 5, m= 2)
root2cfrac(x = 2, n = 2, nterms = 20)
#compare with 
sqrt2periodicCfrac(num = 2, nterms = 20)


}
