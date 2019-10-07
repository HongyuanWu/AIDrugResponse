
OS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AIDrugResponse/master/extdata/survivaltime.txt",head=T)
radio<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AIDrugResponse/master/extdata/xiong/radiotherapy.id.txt")
head(OS)
head(radio)
out<-na.omit(data.frame(radio,OS[match(radio$V1,OS$submitter_id),]))
write.table(out[,2:5],file="radiotherapy.survivaltime.txt",sep="\t",col.names = T,row.names = F,quote=F)
