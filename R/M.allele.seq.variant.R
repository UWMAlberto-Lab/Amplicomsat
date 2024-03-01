M.allele.seq.variant<-function(sampleDir="./samplesF.MA",AlleleList,locusInfo,alleleDB, nmotif=2, HeteroThreshold=0.4,coverageF=25){
	
	primers<-read.delim(locusInfo)
	nprimers<-nrow(primers)
	if(length(coverageF)==1){coverageF<-rep(coverageF,nprimers)}
	
	allele.DB<-read.delim(alleleDB)
	allele.DB$Locus<-factor(allele.DB$Locus,levels=primers[,1])
	Original.DB.size<-nrow(allele.DB)
	
	if(length(nmotif)==1){nmotif<-rep(nmotif,nprimers)}
	
	samples<-names(AlleleList)
	nsamples<-length(samples)
	
	SEQ.ALLELES.Samples<-NULL
	
	for(s in 1:nsamples){
		
		SEQ.ALLELES.Loci<-NULL
		
		for(p in 1:nprimers){
			SEQ.ALLELES<-NULL
			p1<-paste0(rep(primers$Motif[p],nmotif[p]),collapse="")
			Alleles.by.Size<-AlleleList[samples[s]][[1]][[p]] # retrieve the alleles found for locus p i sample s
			nalleles<-length(Alleles.by.Size)
			if(nalleles==0){
				SEQ.ALLELES.Loci<-append(SEQ.ALLELES.Loci,list(NA))
				next
			}
			
			#read the fragments for sample s locus p
			#finding the alleles for this loci and sample combination
			sampleTest<-read_lines(paste0(sampleDir,"/",samples[s]),progress=FALSE)
			#Find read sizes to filter only the ones for the alleles 
			
			Reads.FPmatch<-grepl(pattern=primers[p,2],x=sampleTest,fixed=TRUE)
			#Rvs primer
			Reads.RPmatch<-grepl(pattern=as.character(reverseComplement(Biostrings::DNAString(primers[p,3]))),x=sampleTest,fixed=TRUE)
			#find the intersect of fwd and rvs lines
			Reads.bothPrimers.i<-(Reads.FPmatch & Reads.RPmatch)
			Reads.bothPrimers<-sampleTest[Reads.bothPrimers.i]
			#finding the primer
			motif.i<-grepl(p1,x=Reads.bothPrimers,fixed=TRUE)
			Reads.bothPrimers.motif<-Reads.bothPrimers[motif.i]
			ReadSize<-nchar(Reads.bothPrimers.motif) #Ma
			
			for(a in 1:nalleles){
				reads.A<-Reads.bothPrimers.motif[ReadSize==Alleles.by.Size[a]]
				Table.reads.A<-table(reads.A)
				Table.reads.A.max<-max(Table.reads.A)
				if(sum(Table.reads.A)<coverageF[p])next
				#if(Table.reads.A.max<coverageF)next
				index.HighFreq<-as.numeric(which((Table.reads.A/Table.reads.A.max)>HeteroThreshold))# the length of this object is the number of different sequences that have sufficient counts to be cosidered true homoplasy alleles
				Nvariants<-length(index.HighFreq) #how many variants
				MainSeqs<-reads.A[index.HighFreq] #this could be one or more
				alleleDB.Locus.Size<-allele.DB[allele.DB$Locus==primers[p,1] &  allele.DB$SIZE==Alleles.by.Size[a], ] #retrieving the sequences in the database for this locus and allele size
				AlleleCodeMax<-max(allele.DB$AlleleCode[allele.DB$Locus==primers[p,1]])
				
				# if no sequences are retrieved for that size in that loci write new variant and write to log file
				if(nrow(alleleDB.Locus.Size)==0){
					New.DB.entries<-data.frame(Locus=rep(primers[p,1],Nvariants),
																		 SIZE=rep(Alleles.by.Size[a],Nvariants),
																		 VARIANT=letters[(nrow(alleleDB.Locus.Size)+1):(nrow(alleleDB.Locus.Size)+Nvariants)],
																		 AlleleCode=(AlleleCodeMax+1):(AlleleCodeMax+Nvariants),
																		 Type=rep("New allele size",Nvariants),
																		 Sample=rep(samples[s],Nvariants),
																		 SEQ=MainSeqs)
					allele.DB<-rbind(allele.DB,New.DB.entries[,c(1:4,7)])
					
					o1<-order(allele.DB$Locus,allele.DB$SIZE,allele.DB$VARIANT) #reordering by locus and size, Allele codes won't be ordered in a logical way so that they are kept among different databases
					allele.DB<-allele.DB[o1,]
					write.table(New.DB.entries,file="DataBase new entries.txt",append=TRUE,sep="\t",row.names=FALSE,quote=F,col.names=F)
					next
				}
				letterINDEX<-nrow(alleleDB.Locus.Size)+1
				for(v in 1:Nvariants){
					DB.index<-which(alleleDB.Locus.Size$SEQ%in%MainSeqs[v])
					
					if(length(DB.index)==0){        #if allele size exists but not the variant
						AlleleCodeMax<-AlleleCodeMax+1
						New.DB.entries<-data.frame(Locus=primers[p,1],
																			 SIZE=Alleles.by.Size[a],
																			 VARIANT=letters[letterINDEX],
																			 AlleleCode=AlleleCodeMax,
																			 Type=rep("New Sequence variant"),
																			 Sample=samples[s],
																			 SEQ=MainSeqs[v])
						allele.DB<-rbind(allele.DB,New.DB.entries[,c(1:4,7)])
						letterINDEX<-(letterINDEX+1)
						
						o1<-order(allele.DB$Locus,allele.DB$SIZE,allele.DB$VARIANT) #reordering by locus and size, Allele codes won't be ordered in a logical way so that they are kept among different databases
						allele.DB<-allele.DB[o1,]
						write.table(New.DB.entries,file="DataBase new entries.txt",append=TRUE,sep="\t",row.names=FALSE,quote=F,col.names=F)
						SEQ.ALLELES<-c(SEQ.ALLELES, paste0(Alleles.by.Size[a],	letters[letterINDEX]))
						}else{
					    SEQ.ALLELES<-c(SEQ.ALLELES, paste0(alleleDB.Locus.Size$SIZE[DB.index],	alleleDB.Locus.Size$VARIANT[DB.index]))
					}
				}
				SEQ.ALLELES<-unique(SEQ.ALLELES)
			}
			SEQ.ALLELES.Loci<-append(SEQ.ALLELES.Loci,list(SEQ.ALLELES))
		}
		names(SEQ.ALLELES.Loci)<-primers[,1]
		SEQ.ALLELES.Samples<-append(SEQ.ALLELES.Samples,list(SEQ.ALLELES.Loci))
		print(paste0("Completed sample: ",samples[s]))
	}
	names(SEQ.ALLELES.Samples)<-samples
	
	Final.DB.size<-nrow(allele.DB)
	if(Final.DB.size>Original.DB.size){
		write.table(allele.DB,file="Allele DataBase UPDATED.txt",row.names=F,quote=F,sep="\t")
		print("WARNING: Allele database updated and written to file: Allele DataBase UPDATED.txt")
		print("List of new allele entries written to file DataBase new entries.txt")
	}
	
	SEQ.ALLELES.Samples
}