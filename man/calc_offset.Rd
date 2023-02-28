\name{calc_offset}
\alias{calc_offset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates the offset between two given profiles.
}
\description{
Calculates the total  pair-wise vertical distance between a reference and another profile over a range of different x-axis offsets and returns the x-axis offset which results in the minimum vertical distance.
}
\usage{
calc_offset(ref, other, span = seq(-250, 250, 1), return_plot=FALSE, ref_start = 270, ref_end = 940)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ref}{
Reference profile.
}
  \item{other}{
Profile to align to ref.
}
  \item{span}{
List of values to try and offset the profile. 1 measurement = 1 position value.
}
  \item{return_plot}{
Whether to return a plot of the calculated total distance vs. given offset.
}
  \item{ref_start}{
Position in the profile (by measurement #) of the first absorbance value to consider in the aligning range.
}
  \item{ref_end}{
Position in the profile (by measurement #) of the last absorbance value to consider in the aligning range.
}
}
\details{

}
\value{
Returns the offset that minimizes the vertical pairwise distance of the ref and other profiles, in the given range.
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
lapply(c('tidyverse', 'ggthemes', 'pracma', 'ggprism', 'data.table', 'polyprof', 'readxl'), require, character.only = TRUE)

tst1 <- ProfLoader(excel_path = system.file("extdata", "Polysome_profs.xlsx", package = "polyprof"), excel_sheet = '20230207', excel_range = 'A1:L7', prof_path = system.file("extdata", package = "polyprof"), prof_pattern = '^230207', output_as = 'list')

calc_offset(tst1[[2]], tst1[[4]], return_plot = TRUE)
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
