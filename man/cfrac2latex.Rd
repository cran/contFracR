\name{cfrac2latex}
\alias{cfrac2latex}

\title{
Function which builds the LaTeX formula to display the continued fraction 
}
\description{
Take the vectors of numerators and denominators and create the LaTeX formula which will draw the "cascading" continued fraction.   
}
\usage{
cfrac2latex(denomvals, numvals = 1, denrepeats = 1, doellipsis = FALSE,  ...)
}

\arguments{
  \item{denomvals}{
 A vector containing the denominator values, starting with the not-actually-denominator integer part of the continued fraction. If the denominator has a repeat sequence, specify the number of repeats using the argument \code{denrepeats} . See the Details section for more information.  
    
}
  \item{numvals}{	
A vector containing the numerator values.  If the length of this vector is less than the length of the (perhaps repeated) \code{denomvals}, it will be recycled as necessary. 

}
  \item{denrepeats}{	
If the denominator sequence has a repeat pattern, repeat that pattern this many times.  See the Details section for more information. 


}
\item{doellipsis}{
If TRUE, an ellipsis is printed in the deepest denominator, indicating an infinite sequence continues.  If FALSE, not added. 

}
\item{...}{
Reserved for future upgrades


}
}
\details{
The standard notation for a continued fraction defines a sequence [a0; b1, b2, b3,...bn] indicating the formula is  x = a0 + 1/(b1 + 1/(b2 + 1/(b3...))) (or replace the numerators with a specified sequence of values of length n).  To save input effort, if there's a repeat pattern, e.g. [a0,b1,b2,b3,b1,b2,b3,...] then the user can enter the vector \code{denomvals = c(a0,b1,b2,b3)} and enter the desired number of repeats with the argument \code{denrepeats}. 

}
\value{
Three versions are returned.  
 \code{texeqn} contains a text string compatible with the R-console's parser. Use this to write the correct LaTeX string to a file, or to pass it to tools such as \code{\link[graphics]{title}}. 
 \code{tex2copy} can be used if you wish to use Select/Copy/Paste operations on the text displayed in the console. Do not enter this into any R command, as the parser will interpret the backslashes as escape characters, e.g. "backslash-f" turns into a newline.  
  Finally,  \code{eqn} returns an ASCII-only string along the lines of "1/(1+1/(2+1/(3+...)))"  .  

}

\author{
Carl Witthoft,  \email{carl@witthoft.com} 
}



\examples{
cfrac2latex( 1:5,1, doellipsis= FALSE)
#$eqn
#[1] "1 + 1/(2 + 1/(3 + 1/(4 + 1/5 ...)))"
#$tex2copy
#[1] "1 + \frac{1}{2 + \frac{1}{3 + \frac{1}{4 + \frac{1}{5 ...}}}}"

# Notice the additional backslashes to make the console parser happy. 
#$texeqn
#[1] "1 + \\frac{1}{2 + \\frac{1}{3 + \\frac{1}{4 + \\frac{1}{5 ...}}}}"

}
