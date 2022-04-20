\name{cfrac2simple}
\alias{cfrac2simple}

\title{
Function to Convert Continued Fraction to Simple Fraction. 
}
\description{
Given a vector of denominators and optionally numerators for a continued fraction, return the "simple" continued fraction where all numerators are 1. In addition, return the exact "simple" fraction x/y where both are integers.   
}
\usage{
cfrac2simple(denom, num =1, ...)
}

\arguments{
\item{denom}{
  A vector of values representing the standard \code{[b0; a1,a2,a3...]}  continued fraction denominators, starting with the integer part.
}
\item{num}{
A vector of values representing the numerators in a continued fraction. The default is \code{ nums = 1} as is the case for \code{\link[contFracR]{cfrac2num}} . 
}

\item{...}{
Reserved for future use
}

}

\value{
A list containing:
\code{intfrac} ,the integer fraction  as a \code{bigq} value. 
\code{denom}, the denominator sequence for the continued fraction with all numerators equal to one. 
}

\author{
Carl Witthoft,  \email{carl@witthoft.com} 
}

\seealso{
\code{\link[contFracR]{cfrac2num}}, \code{\link[contFracR]{num2cfrac}}
}

\examples{
foon <- c(1,3,5,7,9)
food <- c(1,2,4,6,8,10)
foosimple <- cfrac2simple(denom =food, num=foon)
# compare with these two:
cfrac2num(den=food, num= foon)
cfrac2num(den=foosimple$denom)
}
