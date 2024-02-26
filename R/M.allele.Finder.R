M.allele.Finder<-function(sampleDir="./samplesF.MA", locusInfo="primers",calibrationFiles="CalibrationSamples.txt",
													coverageF=20,RefGenotypes="./GenotypesOut/Genotypes allele table.txt"){
	
	files<-list.files(sampleDir)
	nfiles<-length(files)
	
	primers<-read.delim(locusInfo) 
	nprimers<-nrow(primers)
	
	if(length(coverageF)==1){coverageF<-rep(coverageF,nprimers)}
	if(nprimers!=length(coverageF))stop("coverageF and the number of loci have different sizes")

	CalibrationFile<-read.delim(calibrationFiles)
	nrowsCalibration<-nrow(CalibrationFile)
	rangesV<-CalibrationFile[1,]
	Homo.samples<-CalibrationFile[-1,]
	
	SlideWindowSize<-8 #function argument, could be a vector for each locus I will keep this non-editable now. If this is too wide the code will too slow due to possible high number of combinations
	
	#finding the alleles in the ref population
	Genotypes<-read.delim(RefGenotypes)
	GenotypesA<-apply(
		rbind(as.matrix(Genotypes[,seq(2,ncol(Genotypes),2)]),
					as.matrix(Genotypes[,seq(3,ncol(Genotypes),2)])),
		2, function(x)sort(unique(x)))

	All.Samples.List<-list()
	
	for(f in 1: nfiles){
		All.loci.alleles.List<-list()
		
		for(p in 1:nprimers){
			print(paste0("Starting locus: ",primers[p,1]))
			ref.alleles<-GenotypesA[[p]]
			ref.alleles.LowR<-min(ref.alleles)-4
			ref.alleles.UppR<-max(ref.alleles)+4
			
			WindowBreaks<-seq(ref.alleles.LowR-4,ref.alleles.UppR+4,SlideWindowSize)
			
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
			Freq.HomoPattern<-table(CombFreqs)
			#RefFreq.HomoPattern<-table(CombFreqs)/length(CombFreqs) #relative allele freqs
			
			RefFreq.DF<-data.frame(Size=75:350,Freq=rep(0,length(75:350)))
			RefFreq.DF$Freq[ RefFreq.DF$Size%in%as.numeric(names(Freq.HomoPattern))]<-as.numeric(Freq.HomoPattern)
			RefMax.Index<-which.max(RefFreq.DF$Freq)
			RefFreq.DF$Rel.Freq<-RefFreq.DF$Freq/sum(RefFreq.DF$Freq)
			
			DF.sample.read.msat.sizes<-readRDS(paste0("./ReadSizesRepeatNumbers/",files[f],":",primers[p,1],".Rds")) 
			
			tableFreqs<-table(DF.sample.read.msat.sizes$ReadSize)
			all.Obs.TableFreq<-as.numeric(names(tableFreqs))
			ObsRF.DF<-data.frame(Size=75:350,Freq=rep(0,length(75:350)))  #the 75:350 could be an argument
			ObsRF.DF$Freq[ObsRF.DF$Size%in%as.numeric(names(tableFreqs))]<-as.numeric(tableFreqs)
			ObsRF.DF$Rel.Freq<-ObsRF.DF$Freq/sum(ObsRF.DF$Freq)
			
			#totalCombs<-0
			alleles.Found.All<-NULL	
			for(w in 1:(length(WindowBreaks)-1)){
				print(paste0("starting window from ",WindowBreaks[w]-4, " to ", WindowBreaks[w+1]+4))
				all.in.window<-ref.alleles[ref.alleles>=(WindowBreaks[w]-4) & ref.alleles<=WindowBreaks[w+1]+4] #-4 so that there is overlap between consecutive windows
				print(paste0(c("Alleles in window: ",all.in.window),collapse=" "))
				nall<-length(all.in.window)
				Wind.All<-seq(WindowBreaks[w]-4,WindowBreaks[w+1]+4,1)
				if(nall==0)next
				
				All.Combs<-NULL
				if(nall==1){All.Combs<-list(all.in.window)}else{
					for(x in 1:nall ){		All.Combs<-c(All.Combs,combn(all.in.window,x,simplify=FALSE))} #all possible combinations
				}
				testArray<-data.frame(matrix(ncol=length(75:350)+1,nrow=length(All.Combs))) #the 75:350 could be an argument
				names(testArray)<-c("Comb",75:350)
				testArray[,]<-0
				
				for(c in 1:length(All.Combs)){
					testArray[c,1]<-paste0(All.Combs[[c]],collapse=":")
					all.in.comb<-length(All.Combs[[c]])
					
					if(all.in.comb==1){
						Index.Comb.Size<-which(75:350==All.Combs[[c]])
						Slide.DIF<-RefMax.Index-Index.Comb.Size
						SlideDF<-data.frame(Size=RefFreq.DF$Size-Slide.DIF,Rel.Freq=RefFreq.DF$Rel.Freq)
						
						SlideDF<-	SlideDF[	SlideDF$Size>=75 & SlideDF$Size<=350,]
						testArray[c, which(names(testArray)%in%SlideDF$Size) ]<-SlideDF$Rel.Freq                  
					}else{	
						temp.RF.V<-data.frame(matrix(nrow=length(All.Combs[[c]]),ncol=length(75:350)))
						temp.RF.V[,]<-0
						names(temp.RF.V)<-75:350
						for(a in 1:all.in.comb){  #loop for multiple alleles in same combunation need to slide the ref relative freq for each sum and average in the end
							Index.Comb.Size<-which(75:350==All.Combs[[c]][a])
							Slide.DIF<-RefMax.Index-Index.Comb.Size
							SlideDF<-data.frame(Size=RefFreq.DF$Size-Slide.DIF,Rel.Freq=RefFreq.DF$Rel.Freq)
							SlideDF<-	SlideDF[	SlideDF$Size>=75 & SlideDF$Size<=350,]
							temp.RF.V[a,which(names(temp.RF.V)%in%SlideDF$Size) ]<-SlideDF$Rel.Freq 
						}
						testArray[c,2:ncol(testArray)]<-apply(temp.RF.V,2,sum)/length(All.Combs[[c]])
					}
				}

				SS.Model.Obs<- apply(testArray[,-1],1, function(x){
					All.numbers<-as.numeric(names(x))
					expected<-x[All.numbers>(WindowBreaks[w]-4) & All.numbers<=(WindowBreaks[w+1]+4)]
					Obs<-ObsRF.DF$Rel.Freq[All.numbers>(WindowBreaks[w]-4) & All.numbers<=(WindowBreaks[w+1]+4)]
					sum( (((expected-Obs)^2)/(expected+0.0000001)  ))
				}
				)
				GenotypeSS.Index<- which.min(SS.Model.Obs )
				GenotypeSS<-testArray[GenotypeSS.Index,1]
				L1<-grep(":",GenotypeSS)
				if(length(L1)!=0){ alleles.Found<-as.numeric(unlist(strsplit(GenotypeSS,split=":")))}else{
					alleles.Found<-as.numeric(GenotypeSS)
				}
				alleles.Found.All<-unique(c(alleles.Found.All,alleles.Found))
				print(paste0(c("Alleles found:",alleles.Found),collapse=" "))
			}  # end of finding alleles loop w counter (for sliding window algorithm)
			
			alllesFilteredCov<-NULL
			if(length(alleles.Found.All)>0){
				alllesFilteredCov<-NULL
				for(a in 1:length(alleles.Found.All)){
					if(ObsRF.DF$Freq[ObsRF.DF$Size==alleles.Found.All[a]]>=coverageF[p]){
						alllesFilteredCov<-c(alllesFilteredCov,alleles.Found.All[a])} #replace 20 by coverageF
				}
			}else{alllesFilteredCov<-NA}
			
			All.loci.alleles.List<-append(All.loci.alleles.List,list(sort(alllesFilteredCov)))
		}#end of loop for primers----
		
		names(All.loci.alleles.List)<-primers[,1]
		
		
		#plots
		pdf(file=paste0("Multiallelic plots for sample ",files[f],".pdf"),height=6,width=12)
		for(p in 1:nprimers){
			DF.sample.read.msat.sizes<-readRDS(paste0("./ReadSizesRepeatNumbers/",files[f],":",primers[p,1],".Rds")) 
			xlimits<-range(DF.sample.read.msat.sizes$ReadSize)
			
			hist(DF.sample.read.msat.sizes$ReadSize,breaks=seq(xlimits[1],xlimits[2],1),xlim=c(xlimits[1],xlimits[2]),axes=F,
					 main=paste0(primers[p,1]),xlab="Read size (bp)")
			par(cex=1.3)  #x.label.cex possible parameter name
			axis(1,at=seq(xlimits[1],xlimits[2],1)-0.5, labels=seq(xlimits[1],xlimits[2],1), las=2)
			par(cex=1)
			axis(2)
			scored<-All.loci.alleles.List[[p]]
			nscored<-length(scored)
			#if(nscored==0)next
			YMAX=max(table(DF.sample.read.msat.sizes$ReadSize))
			if(length(nscored)>0){
				for(l in 1:nscored){abline(v=scored[l]-0.5,col="red") }
			}
		}
		dev.off()
		print(paste0("Completed sample ",f," : ",files[f]))
		All.Samples.List[[paste0(files[f])]]<-All.loci.alleles.List
		
	} #end of loop for files----
	All.Samples.List
}