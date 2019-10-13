setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")

manifest2barcode(gdc_manifest.2019-10-12.txt)

files=list.files(pattern="*.FPKM-UQ.txt$",recursive = T)
rnaseqdata<-c()
for(i in 1:length(files)){
  temp<-read.table(files[i],head=F,sep="\t",row.names = 1)
  rnaseqdata<-cbind(rnaseqdata,temp[,1])
  print(i)
}
colnames(rnaseqdata)<-files
save(rnaseqdata,file="rnaseqdata.pancancer.RData")
save.image("rnaseqdata.pancancer.env.RData")



library("randomForest")
library("arm")
library("plyr") 
library("PredictABEL")
library("neuralnet")
library("caret")

manifest2barcode<-function(manifest){
x=read.table(manifest,header = T)
manifest_length= nrow(x)
id= toString(sprintf('"%s"', x$id))
Part1= '{"filters":{"op":"in","content":{"field":"files.file_id","value":[ '
Part2= '] }},"format":"TSV","fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id","size":'
Part3= paste0("\"",manifest_length, "\"", "}")
Sentence= paste(Part1,id,Part2,Part3, collapse=" ")
write.table(Sentence,"Payload.txt",quote=F,col.names=F,row.names=F)
system("curl --request POST --header \"Content-Type: application/json\" --data @Payload.txt \"https://api.gdc.cancer.gov/files\" > barcode.txt")
}

source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/bin/id2phen4.R")

file=list.files(pattern="*FPKM-UQ.txt$",recursive = TRUE)
manifest2barcode("gdc_manifest.2019-10-12.txt")
barcode<-read.table("barcode.txt",sep="\t",head=T)
data<-c()
for(i in 1:length(file)){
  tmp<-read.table(file[i],head=F,sep="\t",as.is=F)  
  data<-cbind(data,tmp[,2])
  print(paste(i,"in",length(file),file[i],sep=" "))
}
rownames(data)<-tmp[,1]
barcode$file_name<-gsub(".gz","",barcode$file_name)
colnames(data)<-id2phen4(barcode[match(unlist(lapply(file,function(x) unlist(strsplit(x,"[/]"))[2])),barcode$file_name),]$cases.0.samples.0.submitter_id)
data<-data[,grep("TCGA",colnames(data))]
data<-data[,match(unique(colnames(data)),colnames(data))]
save(data,file="../TCGA-Pancancer.mRNAseq.RData")
                                          
