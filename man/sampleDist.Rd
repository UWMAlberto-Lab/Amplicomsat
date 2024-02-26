\name{sampleDist}
\alias{sampleDist}
\title{
Plots  barplots with filtered read frequency (coverage) for a selected sample and given locus.
}
\description{
Plots coverage barplots of a selected sample for a given locus
}
\usage{
sampleDist(sample,xlim=c(75,400))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample}{
The name of the sample to plot including the name of the loci. The format is exactly as shown in the main tittle of the plots saved in the Plots 
directory, thus copy them from here.
}
\item{xlim}{
A vector with two values for lower and upper limits of the x-axis
}
}
\details{
A plotting function to zoom in to specific sample locus filtered read plots. Produces a plot similar to those printed to pdfs in the Plots directory
 for a single sample.
}
\value{
A barplot  with the read size coverage of a selected sample for a given locus.
}
\references{
Alberto F (in prep)...
}
\author{
F. Alberto
}
\seealso{
\code{\link{reads.plots}}
}
\keyword{ hplot }
