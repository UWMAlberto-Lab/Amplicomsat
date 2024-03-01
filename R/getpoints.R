getpoints<-function(GenotypeA.DF,locus){
    
    y.limits<-round(locator(2)$y,0)
    LocusColumns<-grep(names(GenotypeA.DF),pattern=locus)
    LocusData<-GenotypeA.DF[,LocusColumns]
    c1<-LocusData[,1]>=y.limits[1] & LocusData[,1]<=y.limits[2] 
    c2<-LocusData[,2]>=y.limits[1] & LocusData[,2]<=y.limits[2] 
    c3<-c1 | c2
    GenotypeA.DF<-GenotypeA.DF[c3,c(1,LocusColumns)]
    GenotypeA.DF<-GenotypeA.DF<-GenotypeA.DF[!is.na(GenotypeA.DF[,2]),]
    GenotypeA.DF
  }

