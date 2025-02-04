reads.plots<-function(sampleDir="./samplesF",locusInfo="primers",x.range=c(90,350),x.label.cex=0.8,RepNumberPlot=FALSE,G.DB=NULL ){
 
  system2(command="mkdir",args="Plots")
  files<-list.files(sampleDir)
  nfiles<-(length(files))
  primers<-read.delim(locusInfo) 
  
  if(!is.null(G.DB)){
  	system2(command="mkdir",args="Plots.Scored")
  	
  	g.DB<-read.delim(G.DB)
  	i.a1<-seq(2,ncol(g.DB)-1,2)
  	i.a2<-seq(3,ncol(g.DB),2)
  	
  }
  
  for(p in 1:nrow(primers)){
  print(paste0("Processing sample plots for locus number ",p," _ ",primers$Locus[p]))
  	if(!is.null(G.DB)){
  		pdf(file=paste0("./Plots.Scored/Plot.",primers[p,1],".pdf"),width=11,height=5)
  	}else{
    pdf(file=paste0("./Plots/Plot.",primers[p,1],".pdf"),width=11,height=5)
  	}
   for(f in 1:nfiles){   
     DF.sample.read.msat.sizes<-readRDS(paste0("./ReadSizesRepeatNumbers/",files[f],"_",primers[p,1],".Rds"))
     if(nrow(DF.sample.read.msat.sizes)==0){next}
     if(nrow(DF.sample.read.msat.sizes)==1 & is.na(DF.sample.read.msat.sizes$ReadSize[1])){next}
     xlimL<-NULL
     xlimU<-NULL
     if( max(DF.sample.read.msat.sizes$ReadSize)>x.range[2]){xlimU<-max(DF.sample.read.msat.sizes$ReadSize)}else{
       xlimU<-x.range[2]}
     if( min(DF.sample.read.msat.sizes$ReadSize)<x.range[1]){xlimL[1]<-min(DF.sample.read.msat.sizes$ReadSize)}else{
       xlimL<-x.range[1]}
        
     #processing genotypes
     if(!is.null(G.DB)){
      A1<-g.DB[g.DB$Sample==files[f],i.a1[p]]
      A2<-g.DB[g.DB$Sample==files[f],i.a2[p]]
      A.u<-unique(c(A1,A2))
      hist.breaks<-seq(0,xlimU,1)
      colV<-rep("grey60",length( hist.breaks))
      colV[A.u]<-"palegreen3"
      hist(DF.sample.read.msat.sizes$ReadSize,breaks=hist.breaks,xlim=c(xlimL,xlimU),axes=F,
      		 main=paste0(files[f],"_",primers[p,1]),xlab="Read size (bp)",col=colV)
      par(cex=x.label.cex)
      axis(1,at=seq(xlimL,xlimU,1)-0.5, labels=seq(xlimL,xlimU,1), las=2)
      par(cex=1)
      axis(2)
     }else{
          
     hist(DF.sample.read.msat.sizes$ReadSize,breaks=seq(0,xlimU,1),xlim=c(xlimL,xlimU),axes=F,
    	main=paste0(files[f],"_",primers[p,1]),xlab="Read size (bp)")
      par(cex=x.label.cex)
      axis(1,at=seq(xlimL,xlimU,1)-0.5, labels=seq(xlimL,xlimU,1), las=2)
      par(cex=1)
      axis(2)
     }
     
    if(RepNumberPlot){
	plot(DF.sample.read.msat.sizes$ReadSize,DF.sample.read.msat.sizes$RepeatNumber,pch=19,
	col=rgb(214,33,0,alpha=5,maxColorValue =255),cex=2,
	xlim=c(xlimL,xlimU),xlab="Read size (bp)",ylab="Repeat number")
	}
	}
  dev.off()
  
  }
}
