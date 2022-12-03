\name{cfrac2num}
\alias{cfrac2num}

\title{
Function to convert a defined set of numerators and denominators into the equivalent number.  
}
\description{
Given a vector of denominators and optionally a vector of numerators, calculate the \code{\link[gmp]{bigq}} fraction and the \code{mpfr} extended precision decimal value.  
}
\usage{
cfrac2num(denom, num = 1, init = 0, nreps= 1, msgrate = NULL )
}

\arguments{
  \item{ denom}{
 A vector in standard continued fraction form \code{[b0; a1,a2,a3,...]} where \code{b0} represents the integer part of the number and the \code{a_j} are the denominator terms. If all terms, including \code{b0} , are the same, optionally enter a single value and specify the length with the \code{numterms} value
}
  \item{num}{	
A vector of numerator values. If all numerators are the same, a single value will suffice.  The default is \code{1} 
}
\item{init}{
In those cases where the denominator values repeat (e.g., for continued fraction form of square roots), it is possible to take the output of a previous run and use it to initialize a new run. This input must be a \code{bigq} fraction. Use with care. 

}
  \item{nreps}{	
If \code{denoms} is a single value, repeat that value \code{nreps} times to set up the continued fraction. Otherwise, assume a periodic denominator sequence was submitted, and repeat the \code{denom[-1] } section \code{nreps} times (the first denominator is the integer term) . 
}

\item{msgrate}{
If desired enter an integer indicating how often a status message should be displayed. Leave as \code{NULL} to suppress messages.  For example, setting thi to 10 will generate a message every 10 sub-fraction reductions. 

}

}
\details{
All calculations are done with \code{\link[gmp]{bigq}} fractions to preserve full precision.  Naturally, a finite input will not yield the exact value for any irrational number.  
 
}
\value{
The exact numeric value is provided in \code{cgmp}, a \code{bigq} fraction.  The \code{mpfr} value is calculated using \code{\link[Rmpfr]{.bigq2mpfr}}, which generates a decimal number of sufficient precision to match the implied precision of the fraction. 
 \code{fracpart} equals \code{cgmp} minus the integer part, if any. 
 \code{denom} and \code{num} provide the vectors of denominator and numerator values used in the calculation. 
 
}

\author{
Carl Witthoft,  \email{carl@witthoft.com} 
}

\seealso{
\code{\link[contFracR]{num2cfrac}}
}

\examples{

cfrac2num(rep(1,10)) # approximate phi

frac2 <- sqrt2periodicCfrac(2) 
cfrac2num(frac2$numericdenom)
#simple cases
cfrac2num(denom=1)
cfrac2num(denom = c(0,2),num=1)
}
