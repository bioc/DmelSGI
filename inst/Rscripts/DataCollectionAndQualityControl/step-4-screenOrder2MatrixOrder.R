
library(abind)
library(rhdf5)

host = system("echo $HOSTNAME",intern=TRUE)
if (host != "schroedinger.embl.de") {
  stop("This script has to be executed on schroedinger")
}

fileScreenData = "../h5/screenData.h5"
fileMatrixDataRaw = "../h5/matrixDataRaw.h5"

if(file.exists(fileMatrixDataRaw)) {
  stop("file ",fileMatrixDataRaw," already exists. Remove file before running this script!")
}

Anno = h5read("Anno",file = fileScreenData)
Phenotype = h5read("Phenotype",file = fileScreenData)
QueryLayout = h5read("QueryLayout",file = fileScreenData)
TemplateLayout = h5read("TemplateLayout",file = fileScreenData)

h5createFile(file=fileMatrixDataRaw)

# screenData = h5read("screenData", file="data/screenData.h5")
# attach(screenData)

##########################
## Annotation
##########################

AnnoTemplate = TemplateLayout[ (TemplateLayout$group == "sample") & (TemplateLayout$TemplateDesign == 1),
                               c("TID","TemplatePlate","group","Symbol","AnnotationSymbol","Name")]
AnnoTemplate = unique(AnnoTemplate)

AnnoTemplateDesign = data.frame(design=1:2)

AnnoQuery = QueryLayout[QueryLayout$group=="sample",c("TID","Batch","Symbol","AnnotationSymbol","Name")]
AnnoQuery = unique(AnnoQuery)

AnnoQueryDesign = data.frame(design=1:2)

AnnoBatch = data.frame(batch=1:12)

AnnoPhenotype = data.frame(phenotype = Phenotype$Phenotype,stringsAsFactors=FALSE)

CTRLnr = rep(NA,nrow(TemplateLayout))
I = which((TemplateLayout$group %in% c("pos","neg"))
          & (TemplateLayout$TemplateDesign == 1))
CTRLnr[I] = rank(paste(TemplateLayout$TemplatePlate[I],TemplateLayout$TID[I],sep="-"),ties.method="first")
I = which((TemplateLayout$group %in% c("pos","neg"))
          & (TemplateLayout$TemplateDesign == 2))
CTRLnr[I] = rank(paste(TemplateLayout$TemplatePlate[I],TemplateLayout$TID[I],sep="-"),ties.method="first")
AnnoCTRL = data.frame(TID = rep("",length(I)), TemplatePlate = NA, group = "",stringsAsFactors=FALSE)
AnnoCTRL$TID[CTRLnr[I]] = TemplateLayout$TID[I]
AnnoCTRL$TemplatePlate[CTRLnr[I]] = TemplateLayout$TemplatePlate[I]
AnnoCTRL$group[CTRLnr[I]] = TemplateLayout$group[I]
AnnoCTRL$Symbol[CTRLnr[I]] = TemplateLayout$Symbol[I]
AnnoCTRL$AnnotationSymbol[CTRLnr[I]] = TemplateLayout$AnnotationSymbol[I]
AnnoCTRL$Name[CTRLnr[I]] = TemplateLayout$Name[I]

I = which((TemplateLayout$group %in% c("posSKD","negSKD"))
          & (TemplateLayout$TemplateDesign == 1))
CTRLnr[I] = rank(paste(TemplateLayout$TemplatePlate[I],TemplateLayout$TID[I],sep="-"),ties.method="first")
I = which((TemplateLayout$group %in% c("posSKD","negSKD"))
          & (TemplateLayout$TemplateDesign == 2))
CTRLnr[I] = rank(paste(TemplateLayout$TemplatePlate[I],TemplateLayout$TID[I],sep="-"),ties.method="first")
AnnoCTRLSKD = data.frame(TID = rep("",length(I)),TemplatePlate = NA, group = "",stringsAsFactors=FALSE)
AnnoCTRLSKD$TID[CTRLnr[I]] = TemplateLayout$TID[I]
AnnoCTRLSKD$TemplatePlate[CTRLnr[I]] = TemplateLayout$TemplatePlate[I]
AnnoCTRLSKD$group[CTRLnr[I]] = TemplateLayout$group[I]
AnnoCTRLSKD$Symbol[CTRLnr[I]] = TemplateLayout$Symbol[I]
AnnoCTRLSKD$AnnotationSymbol[CTRLnr[I]] = TemplateLayout$AnnotationSymbol[I]
AnnoCTRLSKD$Name[CTRLnr[I]] = TemplateLayout$Name[I]

##########################
## Experiment layout
##########################

ExperimentLayout = merge(cbind(TemplateLayout[,c("TemplatePlate", "TemplateDesign",
                                                 "Well", "QueryNr", "TID",
                                                 "group")],CTRLnr),
                         QueryLayout[,c("Plate","TemplatePlate",
                                        "TemplateDesign","QueryDesign","Batch",
                                        "TID", "group")],
                         by=c("TemplatePlate","TemplateDesign"),
                         suffixes = c(".template",".query"))

m = match(sprintf("%s__%s",Anno$Plate, Anno$Well),
      sprintf("%s__%s",ExperimentLayout$Plate, ExperimentLayout$Well))
ExperimentLayout = ExperimentLayout[m,]

##########################
## DKD
##########################

I = which((ExperimentLayout$group.template == "sample") & (ExperimentLayout$group.query == "sample"))
P = cbind(match(ExperimentLayout$TID.template[I], AnnoTemplate$TID),
  ExperimentLayout$TemplateDesign[I],
  match(ExperimentLayout$TID.query[I], AnnoQuery$TID),
  ExperimentLayout$QueryDesign[I],
  rep(1,length(I)))
d = c(nrow(AnnoTemplate),2,nrow(AnnoQuery),2,nrow(AnnoPhenotype))

Ictrl = which((ExperimentLayout$group.template %in% c("pos","neg")) 
          & (ExperimentLayout$group.query == "sample"))
Pctrl = cbind(ExperimentLayout$CTRLnr[Ictrl],
  ExperimentLayout$TemplateDesign[Ictrl],
  match(ExperimentLayout$TID.query[Ictrl], AnnoQuery$TID),
  ExperimentLayout$QueryDesign[Ictrl],
  rep(1,length(Ictrl)))
dctrl = c(nrow(AnnoCTRL),2,nrow(AnnoQuery),2,nrow(AnnoPhenotype))

I = c(I,Ictrl)
Pctrl[,1] = Pctrl[,1] + d[1]
P = rbind(P,Pctrl)
d[1] = d[1] + dctrl[1]

A = list(template = rbind(AnnoTemplate, AnnoCTRL),
         templateDesign = AnnoTemplateDesign,
         query = AnnoQuery,
         queryDesign = AnnoQueryDesign,
         phenotype = AnnoPhenotype)

h5createGroup(file=fileMatrixDataRaw, group="DKD")
h5write(A, file=fileMatrixDataRaw, name="DKD/A")
h5createDataset(file=fileMatrixDataRaw,dataset="DKD/D",dims=d,storage.mode="double",level=0)

d[5] = 1
for (i in 1:nrow(Phenotype)) {
  cat("write DKD phenotype ",i," of ",nrow(Phenotype),"\r")
  Data = h5read(file=fileScreenData,name="Data",index=list(NULL,i))
  D = array(NA, dim=d)
  D[P] = Data[I,]
  h5write(D, file=fileMatrixDataRaw,name="DKD/D",index=list(NULL,NULL,NULL,NULL,i))
}
cat("\n")

##########################
## SKD
##########################

I = which((ExperimentLayout$group.template == "sample") & (ExperimentLayout$group.query == "neg"))
P = cbind(match(ExperimentLayout$TID.template[I], AnnoTemplate$TID),
  ExperimentLayout$TemplateDesign[I],
  ExperimentLayout$Batch[I],
  rep(1,length(I)))
d = c(nrow(AnnoTemplate),2,max(QueryLayout$Batch),nrow(AnnoPhenotype))

Ictrl = which((ExperimentLayout$group.template %in% c("pos","neg")) 
          & (ExperimentLayout$group.query == "neg"))
Pctrl = cbind(ExperimentLayout$CTRLnr[Ictrl],
  ExperimentLayout$TemplateDesign[Ictrl],
  ExperimentLayout$Batch[Ictrl],
  rep(1,length(Ictrl)))
dctrl = c(nrow(AnnoCTRL),2,max(QueryLayout$Batch),nrow(AnnoPhenotype))

I = c(I,Ictrl)
Pctrl[,1] = Pctrl[,1] + d[1]
P = rbind(P,Pctrl)
d[1] = d[1] + dctrl[1]

A = list(template = rbind(AnnoTemplate, AnnoCTRL),
         templateDesign = AnnoTemplateDesign,
         batch = AnnoBatch,
         phenotype = AnnoPhenotype)

h5createGroup(file=fileMatrixDataRaw, group="SKD")
h5write(A, file=fileMatrixDataRaw, name="SKD/A")
h5createDataset(file=fileMatrixDataRaw,dataset="SKD/D",dims=d,storage.mode="double",level=0)

d[4] = 1
D = array(NA, dim=d)
for (i in 1:nrow(Phenotype)) {
  cat("write SKD phenotype ",i," of ",nrow(Phenotype),"\r")
  Data = h5read(file=fileScreenData,name="Data",index=list(NULL,i))
  D[P] = Data[I,]
  h5write(D, file=fileMatrixDataRaw,name="SKD/D",index=list(NULL,NULL,NULL,i))
}
cat("\n")

##########################
## plateControlDKD
##########################

I = which((ExperimentLayout$group.template %in% c("posSKD","negSKD"))
          & (ExperimentLayout$group.query == "sample"))
P = cbind(ExperimentLayout$CTRLnr[I],
  ExperimentLayout$TemplateDesign[I],
  match(ExperimentLayout$TID.query[I], AnnoQuery$TID),
  ExperimentLayout$QueryDesign[I],
  rep(1,each=length(I)))
d = c(nrow(AnnoCTRLSKD),2,nrow(AnnoQuery),2,nrow(AnnoPhenotype))

A = list(template = AnnoCTRLSKD,
             templateDesign = AnnoTemplateDesign,
             query = AnnoQuery,
             queryDesign = AnnoQueryDesign, 
             phenotype= AnnoPhenotype)

h5createGroup(file=fileMatrixDataRaw, group="plateControlDKD")
h5write(A, file=fileMatrixDataRaw, name="plateControlDKD/A")
h5createDataset(file=fileMatrixDataRaw,dataset="plateControlDKD/D",dims=d,
                storage.mode="double",level=0)

d[5] = 1
D = array(NA, dim=d)
for (i in 1:nrow(Phenotype)) {
  cat("write plateControlDKD phenotype ",i," of ",nrow(Phenotype),"\r")
  Data = h5read(file=fileScreenData,name="Data",index=list(NULL,i))
  D[P] = Data[I,]
  h5write(D, file=fileMatrixDataRaw,name="plateControlDKD/D",index=list(NULL,NULL,NULL,NULL,i))
}
cat("\n")

##########################
## plateControlSKD
##########################

I = which((ExperimentLayout$group.template %in% c("posSKD","negSKD"))
          & (ExperimentLayout$group.query == "neg"))
P = cbind(ExperimentLayout$CTRLnr[I],
  ExperimentLayout$TemplateDesign[I],
  ExperimentLayout$Batch[I],
  rep(1,each=length(I)))
d = c(nrow(AnnoCTRLSKD),2,max(QueryLayout$Batch),nrow(AnnoPhenotype))
A = list(template = AnnoCTRLSKD,
             templateDesign = AnnoTemplateDesign,
             batch = AnnoBatch,
             phenotype = AnnoPhenotype)

h5createGroup(file=fileMatrixDataRaw, group="plateControlSKD")
h5write(A, file=fileMatrixDataRaw, name="plateControlSKD/A")
h5createDataset(file=fileMatrixDataRaw,dataset="plateControlSKD/D",dims=d,
                storage.mode="double",level=0)

d[4] = 1
D = array(NA, dim=d)
for (i in 1:nrow(Phenotype)) {
  cat("write plateControlSKD phenotype ",i," of ",nrow(Phenotype),"\r")
  Data = h5read(file=fileScreenData,name="Data",index=list(NULL,i))
  D[P] = Data[I,]
  h5write(D, file=fileMatrixDataRaw,name="plateControlSKD/D",index=list(NULL,NULL,NULL,i))
}
cat("\n")



# inputData = list(DKD = DKD, SKD = SKD, plateControlDKD = plateControlDKD, 
#                  plateControlSKD = plateControlSKD)
# save(inputData, file=file.path("data","inputData.rda"))
# 






# 
# screenData = list(Data = data, Anno = Anno, TemplateLayout = TemplateLayout, QueryLayout = QueryLayout)
# h5createFile(file="data/screenData.h5")
# h5write(screenData, name = "screenData", file="data/screenData.h5")
# 
# 
# 
# 
# 
# 
# h5write(TemplateLayout, name = "TemplateLayout", file="data/screenData.h5")
# h5write(screenData, name = "screenData", file="data/screenData.h5")
# 
# 
# library(rhdf5)
# file.remove("data/screenData.h5")
# load("data/Anno.rda")
# h5createFile(file="data/screenData.h5")
# h5write(Anno, name = "Anno", file="data/screenData.h5")
# 
# 
# f = H5Fopen("data/screenData.h5")
# d = H5Dopen(h5loc=f,name="screenData/Data")
# 
# 
# 
