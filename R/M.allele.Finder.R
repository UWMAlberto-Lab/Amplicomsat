M.allele.Finder<-function(sampleDir="./samplesF.MA", locusInfo="primers",calibrationFiles="CalibrationSamples.txt",
		          ReadSizesRepeatNumbers="./ReadSizesRepeatNumbers/",RefGenotypes="./GenotypesOut/Genotypes allele table.txt"){
	
	files<-list.files(sampleDir)
	nfiles<-length(files)
	system2(command="mkdir",args="M.plots")
	primers<-read.delim(locusInfo) 
	nprimers<-nrow(primers)
	
	
	CalibrationFile<-read.delim(calibrationFiles)
	nrowsCalibration<-nrow(CalibrationFile)
	rangesV<-CalibrationFile[1,]
	Homo.samples<-CalibrationFile[2:(nrowsCalibration-1),]
	coverageF<-as.numeric(CalibrationFile[nrowsCalibration,])
	
	
	SlideWindowSize<-8 #function argument, could be a vector for each locus I will keep this non-editable now. If this is too wide the code will be too slow due to possible high number of combinations
	
	#finding the alleles in the ref population
	#I believe that only the alleles are allowed to be found
	Genotypes<-read.delim(RefGenotypes)
	GenotypesA<-apply(
		rbind(as.matrix(Genotypes[,seq(2,ncol(Genotypes),2)]),
					as.matrix(Genotypes[,seq(3,ncol(Genotypes),2)])),
		2, function(x)sort(unique(x)))
	names(GenotypesA)<-primers$Locus
	
	All.Samples.List<-list()
	
	for(f in 1: nfiles){
		All.loci.alleles.List<-list()
		print(paste0("Starting sample: ",files[f] ))
		for(p in 1:nprimers){
			print(paste0("Starting locus: ",primers[p,1]))
			ref.alleles<-GenotypesA[[p]]                    #I think this means that alleles not in the ref won't be called. I think that's OK because we don't have their freq in the ref population but it won't allow for new alleles to be found
			ref.alleles.LowR<-min(ref.alleles)-4
			ref.alleles.UppR<-max(ref.alleles)+4
			
			WindowBreaks<-seq(ref.alleles.LowR-4,ref.alleles.UppR+4,SlideWindowSize)
			
			Homo.samples<-CalibrationFile[2:(nrowsCalibration-1),p]   #TODO This should not be here, it's being repeated n.samples times instead of n.loci times!!!
			#Creating a reference relative allele frequency for an homozygous (using data from 6 homozygous samples, could be more)
			DF.sample.read.msat.sizes<-readRDS(paste0(ReadSizesRepeatNumbers,Homo.samples[1],"_",primers[p,1],".Rds")) 
			tableRef<-table(DF.sample.read.msat.sizes$ReadSize)
			maxRef<-as.numeric(names(tableRef)[which.max(tableRef)])
			CombFreqs<-DF.sample.read.msat.sizes$ReadSize
			
			if(length(Homo.samples)>1){    #averaging across more than one homozygous ref sample
				for(h in 2:length(Homo.samples)){
					DF.sample.read.msat.sizes<-readRDS(paste0(ReadSizesRepeatNumbers,Homo.samples[h],"_",primers[p,1],".Rds"))
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
			
			#reading observed data
			DF.sample.read.msat.sizes<-readRDS(paste0(ReadSizesRepeatNumbers,files[f],"_",primers[p,1],".Rds")) 
			#if(nrow(DF.sample.read.msat.sizes)==0)next #Likely some operations missing before next----
			tableFreqs<-table(DF.sample.read.msat.sizes$ReadSize)
			all.Obs.TableFreq<-as.numeric(names(tableFreqs))
			F1<-all.Obs.TableFreq%in%75:350
			all.Obs.TableFreq<-all.Obs.TableFreq[F1]
			tableFreqs<-tableFreqs[F1]
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
				if(nall==0)next                                         #clear here, and 4 rows above here, that only alleles present in the ref pop are allowed
				
				All.Combs<-NULL
				if(nall==1){All.Combs<-list(all.in.window)}else{
					for(x in 1:nall ){		All.Combs<-c(All.Combs,combn(all.in.window,x,simplify=FALSE))} #all possible combinations
				}
				testArray<-data.frame(matrix(ncol=length(75:350)+1,nrow=length(All.Combs))) #the 75:350 could be an argument
				names(testArray)<-c("Comb",75:350)
				testArray[,]<-0
				
				#creating a matrix (testArray) with the  allele relative frequencies expected for all combinations of alleles in the slide window
				# the relative frequencies are estimated from their proportions in the homozygous reference population, when several alleles are present
				# in a combination this produces sums that include any possible overlaps of secondary peaks
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
						for(a in 1:all.in.comb){  #loop for multiple alleles in the combination need to slide the ref relative freq for each sum and average in the end
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
					sum(((expected-Obs)^2) /(expected+0.00001) )/(1/sum(expected+0.00001))
				}
				)
				GenotypeSS.Index<- which.min(SS.Model.Obs )
				#print(paste0("N all in comb: ",all.in.comb," ; Min SS: ", min(SS.Model.Obs))) #Just to see which Sum of squares are produced
				GenotypeSS<-testArray[GenotypeSS.Index,1]
				L1<-grep(":",GenotypeSS)
				if(length(L1)!=0){ alleles.Found<-as.numeric(unlist(strsplit(GenotypeSS,split=":")))}else{
					alleles.Found<-as.numeric(GenotypeSS)
				}
				alleles.Found.All<-unique(c(alleles.Found.All,alleles.Found))
				print(paste0(c("Alleles found:",alleles.Found),collapse=" "))
			}  # end of finding alleles loop w counter (for sliding window algorithm)
			#The code above always finds some allele, I think below we will filter based on coverage 
			
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
		pdf(file=paste0("./M.plots/Multiallelic plots for sample ",files[f],".pdf"),height=6,width=12)
		for(p in 1:nprimers){
			DF.sample.read.msat.sizes<-readRDS(paste0(ReadSizesRepeatNumbers,files[f],"_",primers[p,1],".Rds")) 
			if(nrow(DF.sample.read.msat.sizes)==0)next
			xlimits<-range(DF.sample.read.msat.sizes$ReadSize)
			scored<-All.loci.alleles.List[[p]]
			breaks.hist<-seq(xlimits[1],xlimits[2],1)
			cols.hist<-c("grey60","royalblue3")[(as.numeric(breaks.hist[-1]%in%scored)+1)]
			hist(DF.sample.read.msat.sizes$ReadSize,breaks=breaks.hist,xlim=c(xlimits[1],xlimits[2]),axes=F,
					 main=paste0(primers[p,1]),xlab="Read size (bp)",col=cols.hist)
			par(cex=1.3)  #x.label.cex possible parameter name
			axis(1,at=seq(xlimits[1],xlimits[2],1)-0.5, labels=seq(xlimits[1],xlimits[2],1), las=2)
			par(cex=1)
			axis(2)
			
		}
		dev.off()
		print(paste0("Completed sample ",f," : ",files[f]))
		All.Samples.List[[paste0(files[f])]]<-All.loci.alleles.List
		
	} #end of loop for files----
	All.Samples.List
}
