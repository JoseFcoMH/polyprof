\name{peak_finder}
\alias{peak_finder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find n peaks in a given dataframe.
}
\description{
Takes a dataframe as input and uses the pracma::findpeaks function to find n peaks in the profile within the specified range. Outputs a dataframe consisting of the rows corresponding to peaks in the input dataframe.
}
\usage{
peak_finder(dtf, npeaks=4, show_peaks=FALSE, minPeakPos=15, maxPeakPos=70)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dtf}{
The profile you want to scan for peaks
}
  \item{npeaks}{
The number of peaks you want to keep, being the first n peaks found in the specified range.
}
  \item{show_peaks}{
If TRUE, prints a plot which shows an overlay of the peaks found in the given profile. Also specifies which peaks are present in the output (good == TRUE).
}
  \item{minPeakPos}{
Start of the range you want to scan for peaks in (in Position(mm)).
}
  \item{maxPeakPos}{
End of the range you want to scan for peaks in (in Position(mm)).
}
}
\details{
Works best when peaks are clear (i.e 40S, 60S, 80S, polysomes), but can handle a fair amount of noise.
}
\value{
Returns a dataframe consisting of the rows corresponding to peaks in the input dataframe.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
peak_finder(dtf, npeaks = 4, show_peaks = TRUE, minPeakPos = 15, maxPeakPos = 70)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
