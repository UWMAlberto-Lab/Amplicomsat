Edit.Primer.File<-function(sampleDir="./samplesF",locusInfo="primers",calibN=2,mismatch=1){
primers<-read.delim(locusInfo) #This will be a required file
nprimers<-nrow(primers)

files<-list.files(sampleDir)
nfiles<-length(files)

f.r<-sample(1:nfiles,calibN)

Primers.Edited<-data.frame(matrix(ncol=3,nrow=nprimers))
names(Primers.Edited)<-names(primers)[1:3]
Primers.Edited$Locus<-primers$Locus

LogTable<-data.frame(matrix(ncol=9,nrow=nprimers))
names(LogTable)<-c("Locus","FprimerKept","RprimerKept","Old.F.primer","New.F.primer",
									 "Old.R.primer","New.R.primer","add.F.primer","add.R.primer")
LogTable$Locus<-primers$Locus

cut.these<-NULL

for(p in 1:nprimers){
	
	matchesR.all<-NULL
	matchesF.all<-NULL
	
	fcount<-ifelse(is.numeric(calibN),calibN,length(calibN))
	
	for(f in 1:fcount){
		
		if(is.numeric(calibN)){
			        sampleTest<-read_lines(paste0(sampleDir,"/",files[f.r]),progress=FALSE)
	         	}else{
	         		sampleTest<-read_lines(paste0(sampleDir,"/",calibN[f]),progress=FALSE)	
	          	}
		SEQS<-sampleTest[seq(2,length(sampleTest),4)]
		
		#Searching for the primer allowing for a mismatch of one
		Search.F.PrimerRES<-sapply(SEQS,function(x){matchPattern(pattern=primers[p,2],subject=x,max.mismatch=mismatch)}) 
		#indexes for reads where a match was found
		indexFWD.P<-as.numeric(which(unlist(sapply(Search.F.PrimerRES,function(x){length(x@ranges@start)!=0})))) 
		if(length(indexFWD.P)==0)next
		# a loop to extract all the matches sequences, some might have up to 1 bp mismatched (I checked that this allows for one fewer bp too)
		matchesF<-1:length(indexFWD.P)
		for(s in 1:length(indexFWD.P)){
			matchesF[s]<-as.character(DNAStringSet(Search.F.PrimerRES[[indexFWD.P[1]]]))
		}
		matchesF.all<-c(matchesF.all,matchesF)
		
		
		#For the reverse primer 
		RevCompRVS.P<-as.character(reverseComplement(DNAString(primers[p,3])))
		Search.R.PrimerRES<-sapply(SEQS,function(x){matchPattern(pattern=RevCompRVS.P,subject=x,max.mismatch=mismatch)}) 
		#indexes for reads where a match was found
		indexRVS.P<-as.numeric(which(unlist(sapply(Search.R.PrimerRES,function(x){length(x@ranges@start)!=0})))) 
		if(length(indexRVS.P)==0)next
		
		# a loop to extract all the matches sequences, some might have up to 1 bp mismatched (I checked that this allows for one fewer bp too)
		matchesR<-1:length(indexRVS.P)
		for(s in 1:length(indexRVS.P)){
			matchesR[s]<-as.character(DNAStringSet(Search.R.PrimerRES[[indexRVS.P[1]]]))
		}
		matchesR.all<-c(matchesR.all,matchesR)
		print(paste0("Completed sample ", f, " for locus ", primers[p,1]))
	}
	#### CONDTION TO FILL THE TABLES IF NO PRIMERS WERE FOUND
	if(is.null(matchesF.all) | is.null(matchesR.all)){
		cut.these<-c(cut.these,p)
		next}
	#Counting how many times each match showed up for FWD Primer
	T.FWD<-table(matchesF.all)
	Primers.Edited[p,2]<-names(T.FWD)[which.max(T.FWD)]  #extraction the most common match. What about if there is more?
	
	if(length(T.FWD)>1){
		sortF<-sort(T.FWD,decreasing=T)
		if((sortF[2]/sort[1])>0.1){add.Fprimer<-names(sortF)[2]}
	}else{
		add.Fprimer<-"NA"
	}
	
	if(Primers.Edited[p,2]==primers[p,2]){
		LogTable$FprimerKept[p]<-TRUE
		LogTable$Old.F.primer[p]<-primers$Fprimer[p]
		LogTable$New.F.primer[p]<-"NA"
		LogTable$add.F.primer[p]<-add.Fprimer
	}else{
		LogTable$FprimerKept[p]<-FALSE
		LogTable$Old.F.primer[p]<-primers$Fprimer[p]
		LogTable$New.F.primer[p]<-Primers.Edited[p,2]
		LogTable$add.F.primer[p]<-add.Fprimer
	}
	
	#Counting how many times each match showed up for RVS primer
	T.RVS<-table(matchesR.all)
	Primers.Edited[p,3]<-as.character(reverseComplement(DNAString(names(T.RVS)[which.max(T.RVS)])))  #extraction the most common match. What about if there is more?
	
	if(length(T.RVS)>1){
		sortR<-sort(T.RVS,decreasing=T)
		if((sortR[2]/sort[1])>0.1){add.Rprimer<-as.character(reverseComplement(DNAString(names(sortR)[2])))}
	}else{
		add.Rprimer<-"NA"
	}
	
	if(Primers.Edited[p,3]==primers[p,3]){
		LogTable$RprimerKept[p]<-TRUE
		LogTable$Old.R.primer[p]<-primers$Rprimer[p]
		LogTable$New.R.primer[p]<-"NA"
		LogTable$add.R.primer[p]<-add.Rprimer
	}else{
		LogTable$RprimerKept[p]<-FALSE
		LogTable$Old.R.primer[p]<-primers$Rprimer[p]
		LogTable$New.R.primer[p]<-Primers.Edited[p,3]
		LogTable$add.R.primer[p]<-add.Rprimer
	}
	print(paste0(p ," markers done, out of ", nprimers))
}


Primers.Edited<-cbind(Primers.Edited,primers[,4:5])
#remove rows where one of the primers was not f
if(!is.null(cut.these)){Primers.Edited<-Primers.Edited[-cut.these,]}
write.table(Primers.Edited,"primers.edited",sep="\t",row.names=F,quote=F)
write.table(LogTable,"Log of changes to primers file",sep="\t",row.names=F,quote=F)
Primers.Edited
}