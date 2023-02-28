\name{normalize}
\alias{normalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalizes the given profile to various possible outcomes.
}
\description{
Takes a profile and applies the chosen operation to it. This includes normalizing to total AUC or estimated ribosome range AUC. Can also make the profile 'smoother', removing some noise, by omiting certain values or overfitting a curve to the profile.
}
\usage{
normalize(dtf, max_abs = Inf, pos_start = 7, pos_end = 60, pos_offset = 0, to = '', max_jump = Inf, smoothen = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dtf}{
Profile to normalize.
}
  \item{max_abs}{
Max absorbance for the profile. Higher values will be substituted by the max_abs value.
}
  \item{pos_start}{
Start of the range of the profile to keep (in Position(mm)).
}
  \item{pos_end}{
End of the range of the profile to keep (in Position(mm)).
}
  \item{pos_offset}{
Offset to apply to the profile (in Position(mm)), if for any reason it has to be done manually.
}
  \item{to}{
Whether to normalize to total AUC in the kept range (to = 'AUC'), to total AUC in the estimated ribosome range (to = 'AUC_Rib') or 80S peak (to = '80S', not recommended).
}
  \item{max_jump}{
Max consecutive increase in value to allow. May help with certain types of noise.
}
  \item{smoothen}{
Whether to overfit a curve to the profile and substitute its Absorbance values for the values predicted by the overfitted curve. Usually redudes noise drastically while preserving the shape of the curve.
}
}
\details{
Use at your own risk.
}
\value{
Returns a profile with the specified modifications.
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

tst1_normalized <- lapply(tst1, normalize, max_abs = 0.5, to = 'AUC')
tst1_normalized
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
