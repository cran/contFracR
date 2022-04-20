\name{recipCfrac}
\alias{recipCfrac}

\title{
Function to Generate the Reciprocal of A Continued Fraction 
}
\description{
Given the numerators and denominators of a continued fraction, generate the numerators and denominators of the reciprocal of that fraction.   
}
\usage{
recipCfrac( denom, num = 1, ...)
}

\arguments{
  \item{ denom}{
  A vector of the continued fraction denominators.
}
  \item{num}{	
A vector of the continued fraction numerators.  Default is 1, indicating that all numerators have the value 1. 
}
  \item{...}{	
Reserved for future use. 
}
 
}
\details{
This is a Q&D tool to perform a simple operation.  It can easily be shown that the reciprocal of a continued fraction is obtained by prepending \code{0} to the current denominator sequence and prepending \code{1} to the current numerator sequence.  To be explicit, given the standard denominator notation \code{[b0; a1,a2,a3...]} , the new denominator vector is \code{[ 0; b0, a1,a2,a3...]}  
 
}
\value{
A list containing \code{denom} and \code{num} (numerators) for the reciprocal value.  
}

\author{
Carl Witthoft,  \email{carl@witthoft.com} 
}




\examples{
foon <- c(1,3,5,7,9)
food <- c(1,2,4,6,8,10)
foor <- recipCfrac(num = foon, den=food)
# compare:
cfrac2num(num=foor$num, denom = foor$denom)$cgmp
original <- cfrac2num(num=foon, denom=food)
original$cgmp *cfrac2num(num=foor$num, denom = foor$denom)$cgmp


}
