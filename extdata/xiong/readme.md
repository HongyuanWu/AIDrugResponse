```R
sgene<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AIDrugResponse/master/extdata/xiong/sgene.txt",head=T)
ENSG2Symbol<-function(ENSG){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
  ENSG<-unlist(lapply(strsplit(ENSG,split="[.]"),function(x) x[1]))
  Symbol<-db[match(as.character(ENSG),db$V8),4]
  return(Symbol)
}

out<-data.frame(sgene,symbol=ENSG2Symbol(as.character(sgene[,2])))
write.table(out,file="TCGA_DATA_TTEST_10011409_Linear.symbol.txt",sep="\t",quote=F,col.names = NA,row.names = T)


data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AIDrugResponse/master/extdata/survivaltime.txt",head=T,sep="\t")
d1109<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AIDrugResponse/master/extdata/overlapid/id.1109.txt",head=T,sep="\t")
out<-data.frame(d1109,data[match(d1109$bcr_patient_barcode,data$submitter_id),])
write.table(out,file="S1109.SurvivalTime.txt",sep="\t",quote=F,col.names = NA,row.names = T)
```
