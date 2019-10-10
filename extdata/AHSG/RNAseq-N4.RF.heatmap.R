################
setwd("/home/guosa/hpc/project/TCGA")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
phen<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/drug_response/pancancer.chemotherapy.response.txt",head=T,sep="\t")
barcode<-read.table("~/hpc/project/TCGA/pancancer/FPKM/barcode.txt",head=T,sep="\t")
load("~/hpc/methylation/Pancancer/RNA-seq/rnaseqdata.pancancer.env.RData")
source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/bin/id2phen4.R")
ncn<-barcode[match(unlist(lapply(colnames(rnaseqdata),function(x) unlist(strsplit(x,"[/]"))[2])),barcode$file_name<-gsub(".gz","",barcode$file_name)),]
ncol<-match("cases.0.samples.0.submitter_id",colnames(ncn))
colnames(rnaseqdata)<-ncn[,ncol]
phen$ID<-paste(phen$bcr_patient_barcode,"-01",sep="")
rnaseq<-rnaseqdata[,na.omit(match(unique(phen$ID),id2phen4(colnames(rnaseqdata))))]
rnaseq<-rnaseq[which(unlist(apply(rnaseq,1,function(x) sd(x)>0))),]
colnames(rnaseq)<-id2phen4(colnames(rnaseq))
newphen<-phen[unlist(lapply(colnames(rnaseq),function(x) match(x,phen$ID)[1])),]
sort(table(newphen$bcr_patient_barcode))
table(newphen$measure_of_response)
input<-data.frame(phen=newphen$measure_of_response,log(t(rnaseq)+1,2))
P=apply(input[,2:ncol(input)],2,function(x) summary(glm(as.factor(input[,1])~x,family=binomial))$coefficients[2,4])
input<-input[,c(1,which(P<0.05/length(P))+1)]
dim(input)
RF <- randomForest(as.factor(phen) ~ ., data=input, importance=TRUE,proximity=T)
imp<-RF$importance
head(imp,n=5)
imp<-imp[order(imp[,4],decreasing = T),]
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/HeatMap.R")
newinput<-t(input[,match(rownames(imp)[1:50],colnames(input))])
colnames(newinput)<-input[,1]
pdf("heatmap.randomForest.pdf")
HeatMap(newinput)
dev.off()
save.image("RNAseq-N2.RF.heatmap.RData")
