\name{pell}
\alias{pell}

\title{
Function to Solve Pell's Equation With Continued Fractions
}
\description{
Solve Pell's equation, \eqn{x^2 - d*y^2 = 1} ,  using continued fractions. 
}
\usage{
pell(d, maxrep=10 , A = NULL, B = NULL)
}

\arguments{
  \item{ d}{
  The integer which multiplies y in the equation. For obvious reasons, there is no value to selecting a perfect square.  
}
\item{maxrep}{
The maximum number of times to repeat the periodic section of the continued fraction while generating convergents.  This exists only to avoid the possibility of an infinite "while" loop should something have gone wrong. 
}
\item{A,B}{
Should a  run of \code{pell} fail to find a solution (due, most likely, to a value of the argument \code{maxrep}), enter the convergents \code{A} and \code{B} from the previous run. This simply reduces runtime by avoiding the need to recalculate these convergent pairs.
}
 
}

\details{
As is nicely presented in the Wikipedia entry for Pell's Equation, 
Let h[i] / k[i] denote the sequence ofconvergents to the regular continued fraction for sqrt(n). This sequence is unique. Then the pairsolving Pell's equation and minimizing x satisfies x1 = h[i] and y1 = k[i] for some i. This pair is called thefundamental solution. Thus, the fundamental solution may be found by performing the continued fraction expansion and testing each successive convergent until a solution to Pell's equation is found.
The sequence of convergents for the square root of seven are
 h/k(convergent)	h^2-7*k^2(Pell-type approximation)
 2/1			3
 3/1			+2
 5/2			3
 8/3			+1
 so 8 / 3 is the one we want 

}
\value{
A list.
\code{wins} contains the solutions found
\code{d} echo back the input argument used.
\code{cfracdata} contains the output generated with \code{\link[contFracR]{sqrt2periodicCfrac} }
\code{A,B} the numerators and denominators of the convergents calculated. 

}

\author{
Carl Witthoft,  \email{carl@witthoft.com} 
}



\references{

\url{https://mathworld.wolfram.com/eContinuedFraction.html}
\url{https://en.wikipedia.org/wiki/Pell's_equation}

}
