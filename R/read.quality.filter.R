read.quality.filter<-function(sampleDir="./samples",filterDir="samplesF"){
  
system2(command="mkdir",args=filterDir)
files<-list.files(sampleDir)
for(f in 1:length(files)){
  sampleTest<-readLines(paste0(sampleDir,"/",files[f]))
  if(length(sampleTest)==0)next
  sampleTest[seq(4,length(sampleTest),4)]
  QualitySearch<-grep(sampleTest[seq(4,length(sampleTest),4)],pattern="[^?@ABCDEFGHI]")
  QR2remove<-QualitySearch*4
  QR2removeALL<-c(QR2remove-3,QR2remove-2,QR2remove-1,QR2remove)
  sampleTestF<-sampleTest[-QR2removeALL]
  readr::write_lines(sampleTestF,file=paste0("./",filterDir,"/",sub(".fastq",".F.fastq",files[f])))
  print(paste0("Filtered file ", files[f]," - ",length(QualitySearch)," reads removed"))
  readr::write_lines(paste0(files[f]," - ",length(QualitySearch),"(",round((length(QualitySearch)/(length(sampleTest)/4)*100),2),"%) reads removed")
              ,append = T,file="FilterLog")
  }
}

