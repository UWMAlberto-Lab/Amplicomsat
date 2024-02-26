\name{getpoints}
\alias{getpoints}
\title{
Interacts with an open plot to extract sample names
}
\description{
A tool to extract sample name information
}
\usage{
getpoints(GenotypeA.DF,locus)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{GenotypeA.DF}{
A data.frame as produced by \code{\link{alleleFinder}}
}
\item{locus}{
A string with the name of the locus to plot. The name needs to be written as in the file provided to locusInfo argument to
other functions, e.g., \code{\link{read.size.freqs}}
}
}
\details{
Interacts with the plots produced by \code{\link{AlleleCum}} to extract sample name information for specific allele sizes. Upon running the 
command the R locator function is called for plot interaction. The user needs to perform two clicks on the plot, first by moving the cursor 
below the targeted allele points, the second click above them. The list of samples in GenotypeA.DF that contain alleles sizes within the y range 
(size in bp) of the two clicks is printed to the screen.
}
\value{
The sample names containing the selected allele sizes are printed to screen.
}
\references{
Alberto F (in prep)...
}
\author{
F. Alberto
}
\seealso{
\code{\link{AlleleCum}}
}
\keyword{ hplot }
