read.quality.filter<-function(sampleDir="./samples",Phred.threshold=30,All.below.Q=5,filterDir="samplesF"){

PhredString<-"+,-./0123456789:;<=>?@ABCDEFGHI" #Minimum quality allowed is 10

PhredString.T<-paste0("[^",substr(PhredString,start=Phred.threshold-9, stop=nchar(PhredString)),"]+")


system2(command="mkdir",args=filterDir)
files<-list.files(sampleDir)
pdf(file="Distribution of bases below Q threshold quality.pdf")
for(f in 1:length(files)){
  sampleTest<-readLines(paste0(sampleDir,"/",files[f]))
  if(length(sampleTest)==0)next

  QualitySearch<-lapply(strsplit(sampleTest[seq(4,length(sampleTest),4)],split=""),
  		function(x){grepl(x,pattern=PhredString.T)})
      
  QualitySearch.SUM<-sapply(QualitySearch,sum)
  QualitySearch.R<-which(QualitySearch.SUM>All.below.Q)
  
  if(identical(QualitySearch.R, integer(0))){
  	readr::write_lines(sampleTest,file=paste0("./",filterDir,"/",sub(".fastq",".F.fastq",files[f])))
  	print(paste0("Filtered file ", files[f]," - ",length(QualitySearch.R)," reads removed"))
  	readr::write_lines(paste0(files[f]," - ",length(QualitySearch.R),"(",round((length(QualitySearch.R)/(length(sampleTest)/4)*100),2),"%) reads removed")
  										 ,append = T,file="FilterLog")
  	hist(QualitySearch.SUM,	xlab="Number of bases below Q threshold quality",
  			 main=paste0(files[f]," Q threshold: ",Phred.threshold))
  }else{
   
  QR2remove<-QualitySearch.R*4
  QR2removeALL<-c(QR2remove-3,QR2remove-2,QR2remove-1,QR2remove)
  sampleTestF<-sampleTest[-QR2removeALL]
  readr::write_lines(sampleTestF,file=paste0("./",filterDir,"/",sub(".fastq",".F.fastq",files[f])))
  print(paste0("Filtered file ", files[f]," - ",length(QualitySearch.R)," reads removed"))
  readr::write_lines(paste0(files[f]," - ",length(QualitySearch.R),"(",round((length(QualitySearch.R)/(length(sampleTest)/4)*100),2),"%) reads removed")
              ,append = T,file="FilterLog")
  hist(QualitySearch.SUM,	xlab="Number of bases below Q threshold quality",
  	main=paste0(files[f]," Q threshold: ",Phred.threshold))
  }
  }
  dev.off()
}

