GenotypeModel<-function(sampleDir="./samplesF",locusInfo="primers",calibrationFiles="CalibrationSamples.txt"){
  
  system2(command="mkdir",args="ModelTables")
  files<-list.files(sampleDir)
  nfiles<-(length(files))
  primers<-read.delim(locusInfo) 
  nprimers<-nrow(primers)
  CalibrationFile<-read.delim(calibrationFiles)
  nrowsCalibration<-nrow(CalibrationFile)
  rangesV<-CalibrationFile[1,]
  Homo.samples<-CalibrationFile[-1,]
  
    
  for(p in 1:nprimers){
  print(paste0("starting locus: ", primers[p,1]))
  Homo.samples<-CalibrationFile[2:nrowsCalibration,p]
  #Creating a reference relative allele frequency for an homozygous (using data from 6 homozygous samples, could be more)
  DF.sample.read.msat.sizes<-readRDS(paste0("./ReadSizesRepeatNumbers/",Homo.samples[1],":",primers[p,1],".Rds")) 
  tableRef<-table(DF.sample.read.msat.sizes$ReadSize)
  maxRef<-as.numeric(names(tableRef)[which.max(tableRef)])
  CombFreqs<-DF.sample.read.msat.sizes$ReadSize

 if(length(Homo.samples)>1){
 for(h in 2:length(Homo.samples)){
   DF.sample.read.msat.sizes<-readRDS(paste0("./ReadSizesRepeatNumbers/",Homo.samples[h],":",primers[p,1],".Rds"))
   tableTemp<-table(DF.sample.read.msat.sizes$ReadSize)
   maxTemp<-as.numeric(names(tableTemp)[which.max(tableTemp)])
   CombFreqs<-c(CombFreqs,(DF.sample.read.msat.sizes$ReadSize+(maxRef-maxTemp)))
}
}

#restrict by range
range<-as.numeric(unlist(strsplit(as.character(rangesV[p]),split="-") ))
  
CombFreqs<-CombFreqs[CombFreqs>range[1] & CombFreqs<range[2]]
RefFreq.HomoPattern<-table(CombFreqs)/length(CombFreqs)

AlleleRange<-range[1]:range[2]
PossibleGenotypes<-cbind(data.frame(combn(AlleleRange, 2)),t(data.frame(a1=AlleleRange, a2=AlleleRange)))

Model.Table<-data.frame(matrix(nrow=ncol(PossibleGenotypes),ncol=length(1:400)+2 ))

Model.Table[,1:2]<-t(PossibleGenotypes)
names(Model.Table)<-c("a1","a2",1:400)


#Creating a Model table with the expected relative  frequencies for each peak in all different combinations of 2 out of N ,
# N defined in steps of one and limited by the observe allele range for each loci

Model.Table[,3:ncol(Model.Table)]<-0

for(r in 1:nrow(Model.Table)){
       a1<-Model.Table[r,1]
       a2<-Model.Table[r,2]
       t1<-rep(0,400)
       t2<-rep(0,400)
       Peaks1<-as.numeric(names(RefFreq.HomoPattern))+(a1-maxRef)
       Peaks2<-as.numeric(names(RefFreq.HomoPattern))+(a2-maxRef)

       t1[Peaks1]<-as.numeric( RefFreq.HomoPattern)
       t2[Peaks2]<-as.numeric( RefFreq.HomoPattern)
      Model.Table[r,3:ncol(Model.Table)]<-(t1+t2)/2
 if(r%in%seq(1,nrow(Model.Table),250)) print(paste0("Processed ", r, " possible genotypes"))
}
print(paste0("Writing model table for Locus ",primers[p,1]))
saveRDS(Model.Table,file=paste0("./ModelTables/GenotypeModel.",primers[p,1],".Rds"))
  }
  
}

