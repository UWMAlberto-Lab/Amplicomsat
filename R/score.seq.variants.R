score.seq.variants<-function(locusInfo="primers",RefGenotypes,sampleDir="./samplesF",nmotif=2,
														 HeteroThreshold=0.4,alleleDB=NULL,coverageF=25){

	system2("mkdir","AlleleBySeq")
	Genotypes<-read.delim(RefGenotypes)
	Genotypes.by.SEQ<-Genotypes
	Genotypes.by.SEQ[,-1]<-NA

primers<-read.delim(locusInfo) #This will be a required file
nprimers<-nrow(primers)
if(length(nmotif)==1){nmotif<-rep(nmotif,nprimers)}

files<-Genotypes[,1]
nfiles<-(length(files))

#Finding reads with fwd, rvs primers and motif, notice that some motifs were the 
#reverse complementar of what was on QDD, that is to be expected. These were Nl30 Nl40, Nl50, and Nl59
ALLELE.DATABASE<-alleleDB  	

  for(p in 1:nprimers){
   print(paste0("Processing primer ", primers[p,1])) #progress report primer
   p1<-paste0(rep(primers$Motif[p],nmotif[p]),collapse="")
   if(is.null(alleleDB)){REF.ALLELE<-NULL}else
   	                         { REF.ALLELE<-ALLELE.DATABASE[ALLELE.DATABASE$Locus==primers[p,1],]  }
    
   	for(f in 1:nfiles){
   	if(f%in%seq(1,nfiles,50)){ print(paste0("Processing sample ","number ",f," : ", files[f]))}
   	
   	#finding the alleles for this loci and sample combination
   	sampleTest<-read_lines(paste0(sampleDir,"/",files[f]),progress=FALSE)
   	#Find read sizes to filter only the ones for the alleles 

  
   	Genes1<-Genotypes[f,(p*2)]
   	if(is.na(Genes1))next
   	##################################   
   	Genes2<-Genotypes[f,(p*2)+1]
   	Genes<-c(Genes1,Genes2)
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
   
    for(g in 1:2){
      ReadsGene1<-Reads.bothPrimers.motif[ReadSize==Genes[g]]
      ReadsCount1<-table( ReadsGene1)
      ReadsCount1.MAX<- max(ReadsCount1)
      if(ReadsCount1.MAX<coverageF)next
      index.HighFreq<-as.numeric(which((ReadsCount1/ReadsCount1.MAX)>HeteroThreshold))# the length of this object is the number of different sequences that have sufficient counts to be cosidered true homoplasy alleles
     Nvariants<-length(index.HighFreq)
   
        #vWarning message if there are more than 2 common alleles in any case
    
     if(Nvariants>2){print(paste0("Warning: Homozygous genotype with ",Nvariants," sequence variants detected for locus ",primers[p,1],", allele size ",Genes[g], ", for sample ",files[f]))
    	              cat(paste0("Warning: Homozygous genotype with ",Nvariants," sequence variants detected for locus ",primers[p,1],", allele size ",Genes[g], ", for sample ",files[f]),
    	                  file="Warning Log",append=T,fill=T)
    	              Genotypes.by.SEQ[f,(p*2)]<-NA
    	              Genotypes.by.SEQ[f,(p*2)+1]<-NA
    	              next  #no need to look at second allele if this is true
  								}
    #If there are two Variants for the first size allele and the sample is heterozygous for size, than we have a "three or more peaks" should write NA and print warning to log
    if(Nvariants==2 & Genes1!=Genes2){print(paste0("Warning: Heterozygous genotype (by size) with ",Nvariants," sequence variants detected for locus ",primers[p,1],", allele size ",Genes[g], ", for sample ",files[f]))
    	                              cat(paste0("Warning: Heterozygous genotype with ",Nvariants," sequence variants detected for locus ",primers[p,1],", allele size ",Genes[g], ", for sample ",files[f]),
    			                          file="Warning Log",append=T,fill=T)
    	                              Genotypes.by.SEQ[f,(p*2)]<-NA
                                  	Genotypes.by.SEQ[f,(p*2)+1]<-NA
                                  	next #no need to look at second (size) allele if this is true
                                    }
    
    Sequence=names(ReadsCount1)[which.max(ReadsCount1)] #the sequence of the dominant read 
    SequenceALL=names(ReadsCount1)[index.HighFreq]     #the sequences of the dominant reads, given the IFs above can only be of size 2 at most 
   
   ########## Start building the allele database and write a sequence bases allele file 
    
    if(is.null(REF.ALLELE)){   #if the database is still empty (ie NULL)
    	REF.ALLELE<-data.frame(Locus=rep(primers[p,1],Nvariants),SIZE=rep(Genes[g],Nvariants),
    												 VARIANT=letters[1:Nvariants],AlleleCode=rep(NA,Nvariants),SEQ=SequenceALL)     #write one or two new alleles into a new database
    	
    	Genotypes.by.SEQ[f,(p*2)]<-paste0(REF.ALLELE$SIZE[1],REF.ALLELE$VARIANT[1]) #writing the new format of allele, basically allele size and a letter for variant
    
    		if(	nrow(REF.ALLELE)==1){Genotypes.by.SEQ[f,(p*2)+1]<-Genotypes.by.SEQ[f,(p*2)]}else{
    		   Genotypes.by.SEQ[f,(p*2)+1]<-paste0(REF.ALLELE$SIZE[2],REF.ALLELE$VARIANT[2])
         	}
    	  break
    }else{
    	
    	  if(any(REF.ALLELE$SIZE==Genes[g])){  # if the allele size already exists
    	  	
    	  	for(v in 1:Nvariants){
    	  	 AlleleSizeSet<-c(REF.ALLELE$SEQ[REF.ALLELE$SIZE==Genes[g]],SequenceALL[v]) #merge sequence of allele to  check with existing sequences
    	  	 if(!duplicated(AlleleSizeSet)[length(AlleleSizeSet)]){ #if the sequence variant is new to the set of the same size sequences
    	  		
    	  		 REF.ALLELE<-rbind(REF.ALLELE,data.frame(
    	  		 	Locus=primers[p,1], 
    	  			SIZE=Genes[g],
    	  			VARIANT=letters[sum(REF.ALLELE$SIZE==Genes[g])+1], #adding a new variant
    	  			AlleleCode=NA,
    	  			SEQ=SequenceALL[v])# HEREEEEE
    	  		  )
    	  	    Genotypes.by.SEQ[f,(p*2)+(g-1)]<-paste0(Genes[g],REF.ALLELE$VARIANT[nrow(REF.ALLELE)])	
    	  	 }else{
    	  	    variantTemp<-REF.ALLELE$VARIANT[ REF.ALLELE$SIZE==Genes[g] & REF.ALLELE$SEQ==SequenceALL[v]] #find the variant letter in the database that matches the 2nd variant on the read
    	  		 	Genotypes.by.SEQ[f,(p*2)+(g-1)]<-paste0(Genes[g],variantTemp)
    	  	 }
    	  	}
    	  }else{ #new allele size into database
    	     	REF.ALLELE<-rbind(REF.ALLELE,data.frame(
    	  	      	Locus=primers[p,1],
    	  	      	SIZE=Genes[g],
    	  	      	VARIANT="a", #adding a new variant
    	  	      	AlleleCode=NA,
    	  	      	SEQ=SequenceALL[1]
    	  	      	)
    	     	)
    	  	  Genotypes.by.SEQ[f,(p*2)+(g-1)]<-paste0(Genes[g],"a")
    	  	
    	  }  
    	 }
    }
 	}	
  REF.ALLELE<-REF.ALLELE[order(REF.ALLELE$SIZE,REF.ALLELE$VARIANT),]
  REF.ALLELE$AlleleCode<-100:(99+nrow(REF.ALLELE))
  base::assign(value=REF.ALLELE,x=paste0("REF.ALLELE-",primers[p,1]))
  ALLELE.DATABASE<-rbind(ALLELE.DATABASE,base::get(x=paste0("REF.ALLELE-",primers[p,1])))
  }  

#outputing a genotype table with numeric codes only to be used by most pop gen software
GenotypeCodes<-data.frame(matrix(nrow=nfiles,ncol=nprimers+1))
GenotypeCodes[,1]<-files
GenotypeCodes[,2:ncol(GenotypeCodes)]<-NA
for(p in 1:nprimers){
	for(f in 1:nfiles ){
		
		a1<-Genotypes.by.SEQ[f,p*2]
		if(is.na(a1))next
		a2<-Genotypes.by.SEQ[f,(p*2)+1]
		if(is.na(a2))next
		
		a1.var<-str_extract(a1,regex("[a-z]"))
		a1.size<-str_extract(a1,regex("[0-9]+"))
		a2.var<-str_extract(a2,regex("[a-z]"))
		a2.size<-str_extract(a2,regex("[0-9]+"))
		
		a1code<-ALLELE.DATABASE$AlleleCode[ALLELE.DATABASE$Locus==primers[p,1] & ALLELE.DATABASE$SIZE==a1.size & ALLELE.DATABASE$VARIANT==a1.var ]
		a2code<-ALLELE.DATABASE$AlleleCode[ALLELE.DATABASE$Locus==primers[p,1] & ALLELE.DATABASE$SIZE==a2.size & ALLELE.DATABASE$VARIANT==a2.var ]
		
		GenotypeCodes[f,p+1]<-as.numeric(paste0(a1code,a2code))
	}
}

 write.table(GenotypeCodes,"./AlleleBySeq/Genotype Codes.txt",sep="\t",row.names=F,quote=F)

 #Output data.base
 write.table(ALLELE.DATABASE,"./AlleleBySeq/Allele DataBase.txt",sep="\t",row.names=F,quote=F)
 
 #Output the new genotype table
 write.table(Genotypes.by.SEQ,file="./AlleleBySeq/Genotypes.by.Sequence.txt",row.names=F,quote=F,sep="\t")
 print("Don't forget to check the Warning log file for problematic data")
 print("New data base written to file Allele DataBase.txt")
}