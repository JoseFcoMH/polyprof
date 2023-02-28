\name{ProfLoader}
\alias{ProfLoader}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load polysome profiles
}
\description{
Lets you load and label all the profiles in a run at once - no more painful one by one labelling. Relies on a companion excel sheet to load and label everything (more on that below).
}
\usage{
ProfLoader(excel_path, excel_sheet, excel_range, prof_path, prof_pattern, output_as = 'list')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{excel_path}{
The path to the reference excel file (a string).
}
  \item{excel_sheet}{
The name of the reference sheet within the reference excel file(a string).
}
  \item{excel_range}{
The range to import in the reference sheet, in Excel syntax (a string).
}
  \item{prof_path}{
Path to the folder with the profiles to import (a string).
}
  \item{prof_pattern}{
Pattern shared among the name file of the profiles to import.
}
  \item{output_as}{
If TRUE, returns a list of the selected profiles. If FALSE, returns a single dataframe containing all the profiles.
}
}
\details{
This function assumes that the excel sheet contains:
  - A column named Sample_ID with unique values for each profile.
  - A column named Sample_name containing the exact name of the .csv file for each profile.
  - A row for each profile, each with unique Sample_ID and Sample_name values.
}
\value{
Outputs the selected profiles either as a list of dataframes or merged in one dataframe.
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

tst1

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
