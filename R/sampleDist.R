sampleDist<-function(sample,xlim=c(75,400)){
  
  DF.sample.read.msat.sizes<-readRDS(file=paste0("./ReadSizesRepeatNumbers/",sample,".Rds"))
  hist(  DF.sample.read.msat.sizes$ReadSize,xlim=xlim,
         breaks=seq(min( DF.sample.read.msat.sizes$ReadSize),
                    max( DF.sample.read.msat.sizes$ReadSize),1),
         axes=F,xlab="Read size (bp)",main=sample)
  axis(1,at=seq(xlim[1],xlim[2],1)-0.5,
       labels=seq(xlim[1],xlim[2],1), las=2)
  axis(2)
}

