
source("GscTools.R")
id1187<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AIDrugResponse/master/extdata/overlapid/id.1187.txt",head=T,sep="\t")
image<-read.table("gdc_manifest.2019-09-12_image.txt",head=T)
id1109<-id1187[id1187[,1] %in% id2phen3(image$filename),]
write.table(id1109,file="id1109.txt",sep="\t",quote=F,col.names = T,row.names = F)
