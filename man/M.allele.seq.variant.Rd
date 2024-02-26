\name{M.allele.seq.variant}
\alias{M.allele.seq.variant}
\title{
Identifies multiple sequence based alleles in eDNA or forensic samples.
}
\description{
The function searches samples for multiple sequenced based alleles for alleles of the same size in each loci. Examples of such samples 
come from forensic analysis, where DNA from multiple individuals might be present, or eDNA samples containing alleles from several
individuals. Sequenced based AllelesSample are searched from a list containing the multiple size-based alleles found
by \code{\link{M.allele.Finder}}.
}
\usage{
M.allele.seq.variant(sampleDir="./samplesF.MA",AlleleList,locusInfo,alleleDB, nmotif=2,
                    HeteroThreshold=0.4,coverageF=25)
}
\arguments{
\item{sampleDir}{
The directory where the filtered reads are stored. This directory should only contain the read files and 
nothing else. By default this is a directory named samplesF written in the session working directory, 
as produced by \code{\link{read.quality.filter}}
}
\item{AlleleList}{
A list with the alleles found for each loci and sample as created by \code{\link{M.allele.Finder}}.
}
\item{locusInfo}{
A tab delimited text file containing the information for the microsatellite loci used. Should contain four columns
with the following header names:\cr
\cr
Locus: The loci names\cr
\cr
Fprimer: The sequence of the FWD primer in upper case\cr
\cr
Rprimer: The sequence of the RVS primer in upper case\cr
\cr
Motif: The sequence for the microsatellite motif (e.g., CTTG for tetranucleotide). Note that for some loci you might need to use the reverse
complementar sequence of what your reference value is. Thus, inspect your reads first or if you are not getting data for some loci try this out.
For example, for a motif CAGG, you might need to use CCTG.
}
\item{alleleDB}{
A data.frame containing the database of allele size, sequence variation, new allele code, and nucleotide sequence, as created by 
\code{\link{score.seq.variants}}. If necessary, the function will update the database with new allele sizes or sequence variants or both, 
in which case the updated database is written to the working directory. The format is a tab-delimited table with different sequence variants 
(rows) and the following required columns:\cr
\cr
Locus: The locus name for the sequence-based allele. Needs to match the loci names supplied in argument locusInfo.\cr
\cr	
SIZE: The size in bp for the sequence-based allele.\cr
\cr
VARIANT: A low case letter code to distinguish the different sequence-based alleles that share the same sequence size.\cr
\cr
SEQ: The nucleotide sequence for each sequence-based allele\cr
\cr
AlleleCode: A numeric code used to write alleles in format required by most population genetics software.\cr
}
\item{nmotif}{
The number of microsatellite repeat motif repetitions to search for. It takes either a single value to be applied to all loci or a 
vector of values. We highly suggest you use the same values for this function as used for \code{\link{read.size.freqs}}.
}
\item{HeteroThreshold}{
When sequence variation exists within reads of the same size for the same loci, two or more sequences will dominate the read frequencies, with many
more in low frequency due to sequencing errors. This argument refers to the proportion of reads in the second and following most frequent variants 
compared to the first. The value sets the minimum allowed proportion required to score these reads as sequence-based variant alleles, i.e., not
sequencing errors. It defaults to 0.4.
}
\item{coverageF}{
Reads scored as sequence-based alleles need to have coverage larger than this value.  Defaults to coverage of 25. It is possible to send a vector with 
different values for each locus. You can input a single value if you want to use the same coverage threshold filter for all loci. We highly recommend using values similar
to what was used in \code{\link{alleleFinder}} or \code{\link{score.seq.variants}}. 
}
}
\details{
The function searches for the presence of sequenced-based alleles in  list containing multiple size-based alleles supplied by the user (as produced 
by \code{\link{M.allele.Finder}}). The function requires a data base of allele variants, as produced by \code{\link{score.seq.variants}}. If new variants are found, 
the allele database will be updated and written to working environment (Allele DataBase UPDATED.txt). A log file (DataBase new entries.txt) is also written to the working environment with the new alleles added to the database. The fifth column in this file informs if the new allele represents a new size or a new sequence variant.
}
\value{
	A list with the sequence-based alleles found for each locus and individual.
}
\references{
	Alberto F (in prep)...
}
\author{
	F. Alberto
}
\seealso{
	\code{\link{score.seq.variants}}, \code{\link{alleleFinder}}
}
\keyword{ IO }