\name{e2cfrac}
\alias{e2cfrac}

\title{
Function to Calculate e, Euler's Number With Continued Fractions
}
\description{
Select a particular formula to get an approximate value for Euler's number or any power \code{e^(x/y)} . 
}
\usage{
e2cfrac(nterms, pownum = 1, powdenom = 1,  method = c('euler','wall','gauss'), ...)
}

\arguments{
  \item{ nterms}{
  How many terms of the continued fraction formula to use.   
}
  \item{pownum}{	
An integer identifying the numerator of the (possibly fractional) power to raise \code{e} to.  See Details.
}
  \item{powdenom}{	
An integer identifying the denominator of the (possibly fractional) power to raise \code{e} to.  See Details.
}
  \item{method}{	
Select the method, i.e. the continued fraction formula, for the calculation. If \code{pownum/powdenom != 1}, the method 'gauss' is automatically chosen.  See Details. 
}
\item{...}{
Reserved for future use
}

}
\details{
The two methods which  calculate \code{e} but not to any power, use the following formulas:
"euler" denominators = 2; 1,2,1,1,4,1,1,6,1,1,8,...
"wall"  denominators = 2; 1,2,3,4,...  and numerators = 1,1,2,3,4,... 
The third method, which can calculate powers of \code{e} is variously listed as derived by Euler or proved by Gauss based on hypergeometric functions, is automatically invoked if \code{pownum/powdenom != 1 } . 
"gauss" denominators = 1; 2y-x,6y,10y,14y, 18y, 22y,.. and numerators = 2x,x^2,x^2,x^2,... where the exponent is x/y . 

Due to the cyclic formula for the "euler" case, the exact number of terms used in the calculation may differ from the input \code{nterms} by 1 or 2 counts. 

}
\value{
\code{edec} contains the \code{mpfr} extended decimal representation.  \code{egmp} contains the fraction representation.
\code{nterms, pownum, powdenom} echo back the input arguments used.
\code{method, nums, denoms} report the actual values used in the calculation. 

}

\author{
Carl Witthoft,  \email{carl@witthoft.com} 
}



\references{

\url{https://mathworld.wolfram.com/eContinuedFraction.html}
\url{https://en.wikipedia.org/wiki/Euler's_continued_fraction_formula}
\url{https://en.wikipedia.org/wiki/Gauss's_continued_fraction}

}

\examples{
e2cfrac(nterms = 10)
e2cfrac(nterms = 10, method = 'wall')
e2cfrac(nterms = 10, pownum = 2)
e2cfrac(nterms  = 10, powdenom = 2)

}
