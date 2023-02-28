\name{align}
\alias{align}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Aligns two given profiles.
}
\description{
Attemps to align two profiles either by offseting one of them by the average pair-wise distance of chosen peaks or minimizing total pair-wise vertical distance. Usually succeeds.
}
\usage{
align(ref, other, span = seq(-250, 250, 1), return_plot=FALSE, ref_start = 270, ref_end = 940, by_peaks = TRUE, npeaks = 4, show_peaks = FALSE, minPeakPos = 15, maxPeakPos = 70)
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
List of values to try and offset the profile, if not aligned by peaks. 1 measurement = 1 position.
}
  \item{return_plot}{
Whether to return a plot of the total calculated distance vs. given offset, if not aligning by peaks.
}
  \item{ref_start}{
Position in the profile (by measurement #) of the first absorbance value to consider in the aligning range, if not aligning by peaks.
}
  \item{ref_end}{
Position in the profile (by measurement #) of the last absorbance value to consider in the aligning range, if not aligning by peaks.
}
  \item{by_peaks}{
Whether to align by considering peaks or total vertical distance in the given range.
}
  \item{npeaks}{
Number of peaks to keep if aligning by peaks. Will keep the first n peaks in the given range.
}
  \item{show_peaks}{
Whether to return a plot of the peaks found and used to align the profiles.
}
  \item{minPeakPos}{
Minimun position value (in mm) at with peaks are searched.
}
  \item{maxPeakPos}{
Maximun position value (in mm) at with peaks are searched.
}
}
\details{
Usually better to try and align by peaks, even if you end up only using 1. But if there are no clear peaks common to all profiles (i. e. clear 80S/40S/etc. peaks) and the profiles have similar overall shape, I would not align by peaks.
}
\value{
Returns the other dataframe, offset to be aligned to the ref dataframe.
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

tst1_aligned <- lapply(tst1, align, ref = tst1[[1]], by_peaks = TRUE, show_peaks  = TRUE)
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
