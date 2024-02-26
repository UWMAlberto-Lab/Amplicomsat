read.size.freqs<-function(locusInfo="primers",sampleDir="./samplesF",nmotif=2){

primers<-read.delim(locusInfo) #This will be a required file
nprimers<-nrow(primers)
if(length(nmotif)==1){nmotif<-rep(nmotif,nprimers)}

files<-list.files(sampleDir)
nfiles<-(length(files))
system2(command = "mkdir",args="ReadSizesRepeatNumbers")
  
#Finding reads with fwd, rvs primers and motif, notice that some motifs were the 
#reverse complementar of what was on QDD, that is to be expected. These were Nl30 Nl40, Nl50, and Nl59

  for(p in 1:nprimers){
   print(paste0("Processing primer ", primers[p,1])) #progress report primer
   p1<-paste0(rep(primers$Motif[p],nmotif[p]),collapse="")

   for(f in 1:nfiles){
    sampleTest<-read_lines(paste0(sampleDir,"/",files[f]),progress=FALSE)
   
   if(f%in%seq(1,nfiles,50)){ print(paste0("Processing sample ","number ",f," : ", files[f]))}
    #using CLI grep to find sequences that match fwd primer 
    Reads.FPmatch<-grepl(pattern=primers[p,2],x=sampleTest,fixed=TRUE)
    #Rvs primer
    Reads.RPmatch<-grepl(pattern=as.character(reverseComplement(Biostrings::DNAString(primers[p,3]))),x=sampleTest,fixed=TRUE)
    #find the intersect of fwd and rvs lines
    Reads.bothPrimers.i<-(Reads.FPmatch & Reads.RPmatch)
    Reads.bothPrimers<-sampleTest[Reads.bothPrimers.i]
 
    #finding the primer
    motif.i<-grepl(p1,x=Reads.bothPrimers,fixed=TRUE)
    Reads.bothPrimers.motif<-Reads.bothPrimers[motif.i]
    ReadSize<-nchar(Reads.bothPrimers.motif) #Main result sequence sizes
    motif<-primers$Motif[p]
    RepeatNumber<-sapply(str_extract_all(Reads.bothPrimers.motif,pattern=motif),length) #Main result number of repeats 
    saveRDS(data.frame(ReadSize,RepeatNumber),
              file=paste0("./ReadSizesRepeatNumbers/",files[f],":",primers[p,1],".Rds")
              )
	 }
  }
}

