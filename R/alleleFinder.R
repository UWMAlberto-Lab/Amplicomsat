alleleFinder<-function(sampleDir="./samplesF",locusInfo="primers",calibrationFiles="CalibrationSamples.txt"){
  
  system2(command="mkdir",args="GenotypesOut")
  files<-list.files(sampleDir)
  nfiles<-(length(files))
  primers<-read.delim(locusInfo) 
  nprimers<-nrow(primers)
  #if(length(coverageF)==1){coverageF<-rep(coverageF,nprimers)}
  CalibrationFile<-read.delim(calibrationFiles)
  coverageF<-as.numeric(CalibrationFile[nrow(CalibrationFile),])
  
  
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
 
 rangeA<-as.numeric(unlist(strsplit(CalibrationFile[1,p],split="-")))
 minA<-rangeA[1]
 maxA<-rangeA[2]
 
 for(f in 1:nfiles){
	DF.sample.read.msat.sizes<-readRDS(paste0("./ReadSizesRepeatNumbers/",files[f],"_",primers[p,1],".Rds")) 
	U.ReadSizes<-unique(DF.sample.read.msat.sizes$ReadSize)
	if(!any(U.ReadSizes%in%(minA:maxA)) ){
		                             GenotypeL.DF[f,p+1]<-"NA"
		                             GenotypeA.DF[f,c(2*p,((2*p)+1))]<-"NA"
	                               next}
	
	tableFreqs<-table(DF.sample.read.msat.sizes$ReadSize)
	tableFreqs<-tableFreqs[	names(tableFreqs)%in%(minA:maxA) ]
	
	t3<-rep(0,400)
	indexTemp<-as.numeric(names(tableFreqs))
	c1<-indexTemp<=400
	indexTemp<-indexTemp[c1]
	t3[indexTemp]<-as.numeric(tableFreqs/sum(tableFreqs))[c1]
	
	NoReadIndex<-which(t3==0)
	Model.Table.S<-Model.Table[(!Model.Table$a1%in%NoReadIndex & !Model.Table$a2%in%NoReadIndex),] #removing rows for possible alleles in the range with zero coverage
	Model.Table.S<-Model.Table.S[,-(NoReadIndex+2)] #removing rows for possible alleles in the range with zero coverage
	t4<-t3[-NoReadIndex]
	
	if(nrow(Model.Table.S)==1){genotype<-as.numeric(Model.Table.S[1:2])}else{
	SS.model<-apply(Model.Table.S[,3:ncol(Model.Table.S)],1,function(x){
		        	sum( ((x-t4)^2)/(x+0.00001))/(1/sum(x+0.00001))
		}
		)
	
	genotype<-as.numeric(Model.Table.S[which.min(SS.model),1:2])
	}
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

