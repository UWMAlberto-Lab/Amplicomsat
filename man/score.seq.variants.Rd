\name{score.seq.variants}
\alias{score.seq.variants}
\title{
Scores alleles based on sequence polymorphism
}
\description{
The function searches for sequence based polymorphism in a data set of size-only scored alleles. 
}
\usage{
score.seq.variants(locusInfo="primers",RefGenotypes,sampleDir="./samplesF",
              nmotif=2,HeteroThreshold=0.4,alleleDB=NULL,coverageF=25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
\item{RefGenotypes}{
A tab delimited file with two columns per locus containing the size only based alleles to be inspected for sequence based variation. The 
format is of file "/GenotypesOut/Genotypes allele table.txt" created by function \code{\link{alleleFinder}}.
}
\item{sampleDir}{
The directory where the filtered reads are stored. This directory should only contain the read files and 
nothing else. By default this is a directory named samplesF written in the session working directory, 
as produced by \code{\link{read.quality.filter}}
}
\item{nmotif}{
Reads are retained if they have at least these many repetitions of the repeat motif. It takes either a single value to be 
applied to all loci or a vector of values
}
\item{HeteroThreshold}{
When sequence variation exists within reads of the same size for the same loci, two variants will dominate the read frequencies, with many more in low frequency due to sequencing errors. This argument refers to the proportion of reads in the second most frequent variant compared to the first. The value sets the minimum allowed proportion required to score the second allele as a sequence-based variant allele, i.e., not a sequencing error. It defaults to 0.4.
}
\item{alleleDB}{
A database with allele size, sequence variation, new allele code, and nucleotide sequence. The default is NULL, in which case the function
will produce the database from the data supplied. When a database is supplied, the function will update it with any new variants 
that might be identified. The format is a tab-delimited table with different sequence variants (rows) and the following required columns:\cr
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
\item{coverageF}{
Reads scored as sequence-based alleles need to have coverage above this frequency otherwise are written as missing data NA. Defaults to coverage of 25.
}
}
\details{
The function searches for sequence-based polymorphism within a size-only genotype table supplied by the user (refGenotypes argument to function). The filtered reads are reanalyzed and searched for sequence-based polymorphism. The analysis of sequence-based variation within loci reads of the same size reveals many variants in low frequency that represent sequencing errors. Two high-frequency reads, which would otherwise be cryptic (by homoplasy) alleles, are often revealed. The function is intended for diploid data, and therefore, if three or more high-frequency sequence-based alleles are identified, the genotype is scored as NA, and a message is logged to a WarningLog file written to the working directory.
}
\value{
Outputs three files written in a folder created in the working directory named AlleleBySeq. The three files are Allele DataBase.txt, containing the reference database of alllele sequence variants found; Genotypes.by.Sequence.txt, a population genetics file with two columns per locus and alleles represented by the concatenation of allele size and variant letter code; and Genotype Codes.txt contains new numeric allele codes, one column per locus, that are required for most population genetics software. The allele database file contains the correspondence between the different allele code formats. A warningLog file is also written directly to the working directory reporting the cases where more than two alleles were found, in which case NAs were written to the output file.
}
\references{
Alberto F (in prep)...
}
\author{
F. Alberto
}
\seealso{
\code{\link{GenotypeModel}}, \code{\link{alleleFinder}}
}
\keyword{ IO }
