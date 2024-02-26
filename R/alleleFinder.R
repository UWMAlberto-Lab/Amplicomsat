alleleFinder<-function(sampleDir="./samplesF",locusInfo="primers",coverageF=25){
  
  system2(command="mkdir",args="GenotypesOut")
  files<-list.files(sampleDir)
  nfiles<-(length(files))
  primers<-read.delim(locusInfo) 
  nprimers<-nrow(primers)
  if(length(coverageF)==1){coverageF<-rep(coverageF,nprimers)}
  
  GenotypeA.DF<-data.frame(matrix(nrow=nfiles,ncol=(nprimers*2)+1))
  names(GenotypeA.DF)<-c("Sample",paste0(rep(primers[,1],each=2),c(".a1",".a2")))
  GenotypeA.DF[,1]<-files
  
  GenotypeL.DF<-data.frame(matrix(nrow=nfiles,ncol=(nprimers+1)))
  names(GenotypeL.DF)<-c("Sample",primers[,1])
  GenotypeL.DF[,1]<-files
  
#Calling alleles
for(p in 1:nprimers){
 print(paste0("Starting allele calling for locus: ", primers[p,1]))
 Model.Table<-readRDS(paste0("./ModelTables/GenotypeModel.",primers[p,1],".Rds"))
 for(f in 1:nfiles){
	DF.sample.read.msat.sizes<-readRDS(paste0("./ReadSizesRepeatNumbers/",files[f],":",primers[p,1],".Rds")) 
	
	tableFreqs<-table(DF.sample.read.msat.sizes$ReadSize)
	
	t3<-rep(0,400)
	indexTemp<-as.numeric(names(tableFreqs))
	c1<-indexTemp<=400
	indexTemp<-indexTemp[c1]
	t3[indexTemp]<-as.numeric(tableFreqs/sum(tableFreqs))[c1]
	
	SS.model<-apply(Model.Table[,3:ncol(Model.Table)],1,function(x){sum(sqrt((x-t3)^2))})
	genotype<-as.numeric(Model.Table[which.min(SS.model),1:2])
	
	coverage1<-as.numeric(tableFreqs[names(tableFreqs)%in%genotype[1]])
	coverage2<-as.numeric(tableFreqs[names(tableFreqs)%in%genotype[2]])
	genotype<-genotype[c(coverage1>coverageF[p],coverage2>coverageF[p])]
	
	if(length(genotype)==0){genotype<-c("NA","NA")}else{
	  if(length(genotype)==1){genotype<-rep(genotype,2)}
	  }
	
	  GenotypeA.DF[f,c(2*p,((2*p)+1))]<-genotype
	  if(genotype[1]=="NA"){  
	    GenotypeL.DF[f,p+1]<-"NA" }else{
	    GenotypeL.DF[f,p+1]<-paste0(genotype,collapse="")
	  }
 
	if(f%in%seq(10,length(files),10))print(paste0("Samples completed ",f))
  }
 print(paste0("Locus ",primers[p,1]," is processed"))
 }
  
#write this to a allefinder specific folder
write.table(GenotypeA.DF,"./GenotypesOut/Genotypes allele table.txt",sep="\t",row.names=F,quote=F)
write.table(GenotypeL.DF,"./GenotypesOut/Genotypes locus table.txt",sep="\t",row.names=F,quote=F)
GenotypeA.DF
}

