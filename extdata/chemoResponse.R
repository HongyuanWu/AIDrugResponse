library("TCGAbiolinks")
# 1) receive pid from TCGAbiolinks
pid<-TCGAbiolinks:::getGDCprojects()$project_id
pid<-pid[grep("TCGA",pid)]
# 1) receive pid from github
pid<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/drug_response/Pid.drugResponse.txt",head=F,sep="\t")
pid<-as.character(pid[,1])

drug2csv<-function(clinical.drug){
  bcr_patient_barcode<-clinical.drug$bcr_patient_barcode
  therapy_types<-clinical.drug$therapy_types
  drug_name<-clinical.drug$drug_name
  measure_of_response<-clinical.drug$measure_of_response
  days_to_drug_therapy_start<-clinical.drug$days_to_drug_therapy_start
  days_to_drug_therapy_end<-clinical.drug$days_to_drug_therapy_end
  therapy_ongoing<-clinical.drug$therapy_ongoing
  new.clinical.drug<-data.frame(bcr_patient_barcode,therapy_types,drug_name,measure_of_response,days_to_drug_therapy_start,days_to_drug_therapy_end,therapy_ongoing)
  return(new.clinical.drug)
}

rlt<-c()
for(i in pid){
  query <- GDCquery(project=i,data.category = "Clinical",file.type = "xml")
  GDCdownload(query)
  clinical.drug <- GDCprepare_clinic(query,"drug")
  drugResponse<-drug2csv(clinical.drug)
  rlt<-rbind(rlt,drugResponse)
}
write.table(rlt,file="Pancancer.drugresponse.txt",sep="\t",col.names=NA,row.names=T,quote=F)
