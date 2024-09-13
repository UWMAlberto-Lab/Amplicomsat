\name{read.quality.filter}
\alias{read.quality.filter}
\title{
Filters the fastq reads for sequencing quality
}
\description{
Filters the reads so that if a single nucleotide is below a phred score of 30 the read is removed. Future versions will allow  
user defined phred scores.
}
\usage{
read.quality.filter(sampleDir="./samples",filterDir="samplesF")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleDir}{
A string for the path of the directory containing the samples. This should contain only the fastq files and nothing else.
The default is a folder named samples
}
\item{filterDir}{
A string with the directory path where the filtered fastq files will be stored. The default is a folder named samplesF.
This folder is created in the session working directory when running the function. When running multiple times, users might want to rename 
this folder.
}
}
\details{
A read is filtered out if a single nucleotide is below a phred score 30. The input fastq files should be stored in a directory
named samples, and the output filtered files will be written in a folder named samplesF created when the function is run. Users can change these names.
}
\value{
Fstq read files filtered by sequence quality. A Filterlog file is written to the working directory with the number and percentage of reads that were 
removed for each sample.
}
\references{
Alberto F (in prep)...
}
\author{
F. Alberto
}
\seealso{
\code{\link{read.size.freqs}}
}
\keyword{ IO }
