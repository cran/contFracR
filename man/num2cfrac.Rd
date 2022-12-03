\name{num2cfrac}
\alias{num2cfrac}

\title{
Function to Generate The Continued Fraction For a Number
}
\description{
This function takes a numeric value or a character string version of a number and returns the continued fraction equivalent representation.  
}
\usage{
num2cfrac(num, denom = 1, ...)
}

\arguments{
  \item{ num}{
  A value of class \code{numeric, bigz, bigq}, or \code{character} representing the numerator of the value to be converted. If the class is \code{bigq} then the input argument \code{denom} is ignored. 
  Warning: if a numeric value with a fractional part is entered, errors will occur due to decimal - binary precision errors. Use either the character string form or enter integers for both numerator and denominator. 
   A character string must be of the form [+,-][X][.][Y][[E,e][+,-]Z]  where X,Y, and Z are strings of numerals. See the Details section for more information. 
}
  \item{denom}{	
Optional argument representing the denominator of the value to be converted. The default is \code{denom = 1} .  A character string must be of the form [+,][X][.][Y][[E,e][+,-]Z]  where X,Y, and Z are strings of numerals. See the Details section for more information.
}
  \item{...}{	
Reserved for future use. 
}
 
}
\details{
The use of \code{mpfr} - class numbers is deliberately disallowed, because the small but generally present binary vs. decimal rounding error will lead to horrific continued fraction representations.
As noted above, character strings must be in the form [+,-][X][.][Y][[E,e][+,-]Z]  where X,Y, and Z are strings of numerals. To be explicit: if an exponent is desired, the "+, -" sign is optional but the numeric value Z is required. 

While irrationals can't be expressed as a finite-digit number, by entering as many digits of a known approximation as desired, the equivalent close approximation as a finite continued fraction can be generated.  See the functions \code{\link[contFracR]{pi2cfrac}}, \code{\link[contFracR]{phi2cfrac}}, \code{\link[contFracR]{e2cfrac}} to generate as many elements of these known generators as desired. 
 
}
\value{
The denominators are returned as a \code{bigz} integer vector in \code{denom}.  The actual numerator and denominator values used in the continued fraction creation are echoed back as \code{bigz} integers in \code{numdenom} .  For example, if the input numerator and denominator were 6 and 15, the values 2 and 5 (reduced fraction) would be returned.  Or, if the numerator input is the character string "1.25" , then the values 5 and 4 will be returned. 
The numerators (A) and denominators (B) of the convergents are provided in \code{convA} and \code{convB}
}

\author{
Carl Witthoft,  \email{carl@witthoft.com} 
}

\seealso{
\code{\link[go2bigq]{go2bigq}}, \code{\link[contFracR]{cfrac2num}}
}

\examples{
num2cfrac(num = 35, denom = 11)
num2cfrac(num = "-5.4967E+2")
num2cfrac(num = "1.3579")
# but this will be "wrong" :
num2cfrac(num = 1.3579)


}
