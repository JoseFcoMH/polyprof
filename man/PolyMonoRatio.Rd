\name{PolyMonoRatio}
\alias{PolyMonoRatio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the polysome to monosome ratio.
}
\description{
Calculates the AUC in two ranges for a profile and returns a ratio of the second range to the first.
}
\usage{
PolyMonoRatio(dtf, mono_start=29, mono_end=37, poly_start=mono_end, poly_end=60, show_cutoffs=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dtf}{
The profile to operate on.
}
  \item{mono_start}{
Start of the monosome range.
}
  \item{mono_end}{
End of the monosome range.
}
  \item{poly_start}{
Start of the polysome range.
}
  \item{poly_end}{
End of the polysome range.
}
  \item{show_cutoffs}{
Whether to return a plot overlaying the start and end of the mono- and polysome areas over the profile.
}
}
\details{

}
\value{
Returns the ratio between the areas of the specified ranges.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
Better to use in normalized profiles to avoid negative values. May also use with lapply after aligning. Of course, could also do other ratios like 60S/40S by juggling the ranges around.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
lapply(c('tidyverse', 'ggthemes', 'pracma', 'ggprism', 'data.table', 'polyprof', 'readxl'), require, character.only = TRUE)

tst1 <- ProfLoader(excel_path = system.file("extdata", "Polysome_profs.xlsx", package = "polyprof"), excel_sheet = '20230207', excel_range = 'A1:L7', prof_path = system.file("extdata", package = "polyprof"), prof_pattern = '^230207', output_as = 'list')

PolyMonoRatio(normalize(tst1[[1]]), mono_start = 33, mono_end = 37, poly_start = 42, show_cutoffs = TRUE)
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
