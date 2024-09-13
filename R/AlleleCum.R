AlleleCum<-function(GenotypeA.DF="./GenotypesOut/Genotypes allele table.txt",locus,locusInfo="primers",ymin=NULL,ymax=NULL){
  
LocusData<-GenotypeA.DF[,grep(names(GenotypeA.DF),pattern=locus)]
primers<-read.delim(locusInfo)
allelesPresent<-c(LocusData[,1],LocusData[,2])
allelesPresent<-allelesPresent[!is.na(allelesPresent)]
totalAlleles<-length(allelesPresent)
range.alleles<-range(allelesPresent)
if(is.null(ymin) & !is.null(ymax)) {ymin<-range.alleles[1]-2}
if(!is.null(ymin) & is.null(ymax)) {ymax<-range.alleles[2]+2}

allelesFreq<-table(allelesPresent)
colVect<-NULL
colors<-c("black","gray50")
for(a in 1:length(allelesFreq)){
  if(a%%2==0){
    colVect<-c(colVect,rep(colors[1],allelesFreq[a]))}else{
    colVect<-c(colVect,rep(colors[2],allelesFreq[a]))
  }
  
}

refAllele<-as.numeric(names(allelesFreq)[which.max(allelesFreq)])
#x11()
refgrid<-c(seq(refAllele,range.alleles[1],-nchar(primers$Motif[primers$Locus==locus])),
  seq(refAllele,range.alleles[2],nchar(primers$Motif[primers$Locus==locus]))[-1])

plot(1:totalAlleles,sort(allelesPresent),xlab="Allele index",ylab="Allele size (bp)",
     ylim=c(ymin,ymax),type="n",main=locus,axes=F)
axis(1)
axis(2, at=refgrid,las=2)
     
abline(h=refgrid,col="grey80"  )

points(1:totalAlleles,sort(allelesPresent),col=colVect)

}     

