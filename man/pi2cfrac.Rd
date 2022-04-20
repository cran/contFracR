\name{pi2cfrac}
\alias{pi2cfrac}

\title{
Function to Calculate Pi In Continued Fraction Form
}
\description{
This function generates the continued fraction coefficients for pi, with three different algorithms to choose from.    
}
\usage{
pi2cfrac(nterms, method = c('brouncker','stern','coleman'),...)
}

\arguments{
  \item{ nterms}{
  The number of terms (basically depth of denominators) to generate.   
}
  \item{method}{	
The method to use, entered as a string. The choice can be abbreviated so long as it's unambiguous. See the Details section for discussion of the methods and their discoverers. 
}
  \item{...}{	
Reserved for future use. 
}
 
}
\details{
The three methods are named for their discoverers: Brouncker in 1655, Stern in  1833, and Coleman & Pickett in 2008.  There are others which can be found in various papers; these may be added in a future release. 
 
}
\value{
A list containing:
The numerators and denominators of the continued fraction in \code{nums} and \code{denoms} , as \code{bigz} integers,
\code{nterms} and \code{method} echoed back for reference
}

\author{
Carl Witthoft,  \email{carl@witthoft.com} 
}

\seealso{
\code{\link[contFracR]{cfrac2num}}
}

\references{
\url{http://scihi.org/william-brouncker-approximation-pi/}
\url{https://divisbyzero.com/2008/12/09/a-new-continued-fraction-for-pi/}
\url{https://mathworld.wolfram.com/PiContinuedFraction.html}

}

\examples{
 pi2cfrac(nterms = 100, method='stern')
 pi2cfrac(nterms = 100, method = 'coleman')


}
