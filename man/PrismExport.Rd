\name{PrismExport}
\alias{PrismExport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exports the given profile list in a format that GraphPad Prism likes
}
\description{
Takes a dataframe and formats it to something GraphPad Prism can easily take as input. That is, one column for x axis values (Position(mm)) and another one for each group with the respective Absorbance values.
}
\usage{
PrismExport(dtf, wider_names = c('Sample_ID'), wider_vals = c('Absorbance'), mode = 'percentage')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dtf}{
The dataframe to format.
}
  \item{wider_names}{
The columns from where to extract the names of the new columns.
}
  \item{wider_vals}{
The column from where to extract the values of the new columns.
}
  \item{mode}{
If 'percentage', wil scale the values of each dtf from 0 (min) to 100 (max).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a dataframe more palatable for GraphPad Prism
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

PrismExport(tst1, mode = 'percentage')
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
