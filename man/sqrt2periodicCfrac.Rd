\name{sqrt2periodicCfrac}
\alias{sqrt2periodicCfrac}

\title{
Function To  Generate Periodic Continued Fraction For Square Roots
}
\description{
This function produces the denominators of the continued fraction form of any square root.   
}
\usage{
sqrt2periodicCfrac(num, denom = 1,  nterms = 50, ...)
}

\arguments{
  \item{ num}{
 An integer, \code{bigz} integer, or a character string form of an integer,  representing the number or the numerator of the number for which the continued fraction form of the square root is to be generated. 
}
  \item{denom}{	
 An integer ,\code{bigz} integer, or a character string form of an integer,representing  the denominator of the number for which the continued fraction form of the square root is to be generated. The default value is 1
 }
  \item{nterms}{	
The maximum number of terms (denominators) to calculate. This is a "safety" limit in case the denominator repeat pattern turns out to be extremely long.  See the Details section. 
}
\item{...}{
Reserved for future use.
}
 
}
\details{
As discussed in the references, this algorithm will produce a periodic sequence of denominator values for any rational input value. 
If the returned value of \code{repeated} is "FALSE" then increase the input argument \code{nterms}.  The default value (50) exists so that the function can terminate rather than spend (possibly undesired) time calculating extremely long denominator repeat sequences. 
 
 
}
\value{
A list, with:
The continued fraction numerators and denominators in \code{bigz} form \code{num} , \code{denom} .
The continued fraction numerators and denominators in numeric form \code{numericnum}, \code{numericdenom} .  In the extreme case that a value exceeds the machine size of a numeric, \code{NA} is returned.
\code{repeated} returns TRUE if the denominator repeat argument is found, and FALSE if not. 
\code{input} echoes back the input num and denom arguments. 

}

\references{
\url{https://r-knott.surrey.ac.uk/Fibonacci/cfINTRO.html}  section6.2
\url{https://math.stackexchange.com/questions/2215918} 

Proof of periodicity:
\url{https://web.math.princeton.edu/mathlab/jr02fall/Periodicity/mariusjp.pdf} 

}
\author{
Carl Witthoft,  \email{carl@witthoft.com} 
}

\seealso{
\code{\link[base]{strtoi}}
}

\examples{
sqrt2periodicCfrac('101011',2,10)
sqrt2periodicCfrac(101011,2,10)
sqrt2periodicCfrac('4r8s',30,10)
sqrt2periodicCfrac('101011',12,5)


}
