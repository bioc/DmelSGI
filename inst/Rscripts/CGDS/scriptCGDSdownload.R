
library(ChromatinSet)
data(Interactions, package="ChromatinSet")
Fly2Human = read.table("Fly2Human.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
ENSG2HGNC = read.table("ENSG2HGNC.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
Fly2HUGO = merge(Fly2Human, ENSG2HGNC, 
                 by.x=2, by.y=1,all.x=TRUE)[,
                   c("Ensembl.Gene.ID","HGNC.symbol","Human...Identity")]

#Fly2HUGO = read.table("FB2Human.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)[,c(3,7,9)]
#s=t(sapply(split(Fly2HUGO,Fly2HUGO[[2]]),function(x) { x[which.max(x[[3]]),,drop=FALSE] }  ))

Fly2HUGO = unique(Fly2HUGO)
Fly2HUGO = Fly2HUGO[nchar(Fly2HUGO[,2]) > 0,]

SP = split(Fly2HUGO[,c(1,3)], Fly2HUGO[[2]])
SP = sapply(SP, function(x) { x[which.max(x[,2]),1] } )
Fly2HUGO = data.frame(SP,names(SP), stringsAsFactors=FALSE)

I = which(Fly2HUGO[,1] %in% Interactions$Anno$template$TID)
Fly2HUGO = Fly2HUGO[I,]

TID2HUGO = tapply(Fly2HUGO[,2], 
                  factor(Fly2HUGO[,1],levels=Interactions$Anno$template$TID),
                  function(x) { unique(x) } )

# number of orthologous genes in human for each fly gene
n = sapply(TID2HUGO,length)
table(n)

# number of orthologous genes in fly for each human gene
table(table(unlist(TID2HUGO)))

# list of human genes with more than one fly orthologue
I = which(table(unlist(TID2HUGO)) > 1)
table(unlist(TID2HUGO))[I]
save(TID2HUGO, file="TID2HUGO.rda")

HUGOgenes = unique(unlist(TID2HUGO))


library(cgdsr)
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

# Test the CGDS endpoint URL using a few simple API tests
test(mycgds) 

# Get list of cancer studies at server
Studies = getCancerStudies(mycgds)

K = c()
# Get available case lists (collection of samples) for a given cancer study  
for (i in 1:nrow(Studies)) {
  mycancerstudy = Studies[i,1]
  mycaselist = sprintf("%s_all",mycancerstudy)

  print("########################")
  P = getGeneticProfiles(mycgds,mycancerstudy)
  I = which(P$genetic_profile_name == "Mutations")
  if (length(I) > 0) {
    K = c(K,i)
  }
}
StudiesWithMutations = Studies[K,]

Mutations = list()
Bad = list()
for (i in 1:nrow(StudiesWithMutations)) {
  mycancerstudy = StudiesWithMutations[i,1]
  mycaselist = sprintf("%s_all",mycancerstudy)
  print("########################")
  P = getGeneticProfiles(mycgds,mycancerstudy)
  I = which(P$genetic_profile_name == "Mutations")
  if (length(I) > 0) {
    mygeneticprofile = P[I,1]
    x = getProfileData(mycgds,HUGOgenes[1],mygeneticprofile,mycaselist)
    Bad[[mycancerstudy]] = c()
    Mutations[[mycancerstudy]] = matrix(NA_character_,nr=nrow(x),nc=length(HUGOgenes))
    row.names(Mutations[[mycancerstudy]]) = row.names(x)
    colnames(Mutations[[mycancerstudy]]) = HUGOgenes
    for (j in seq_len(length(HUGOgenes))) {
      cat("study = ", mycancerstudy, " j=",j," nbad=",length(Bad[[mycancerstudy]]),"           \r")
      e = try( {
        y = getProfileData(mycgds,HUGOgenes[j],mygeneticprofile,mycaselist)
      } )
      k=1
      while(class(e) == "try-error") {
        k=k+1
        e = try( {
          cat("study = ", mycancerstudy, " j=",j," nbad=",length(Bad[[mycancerstudy]])," attempt=",k,"           \r")
          Sys.sleep(time=3)
          y = getProfileData(mycgds,HUGOgenes[j],mygeneticprofile,mycaselist)
        } )
      }
      if (nrow(y) != nrow(Mutations[[mycancerstudy]])) {
        Bad[[mycancerstudy]] = c(Bad[[mycancerstudy]],j)
      } else{
        x = as.character(y[,1])
        names(x) = row.names(y)
        if (length(x) != nrow(Mutations[[mycancerstudy]])) {
          stop(mycancerstudy," j=",j," length of result and nrow of matrix differ")        
        }
        K = match(row.names(Mutations[[mycancerstudy]]),names(x))
        if (any(is.na(K))) {
          stop(mycancerstudy," j=",j," not all cases in result.")
        }
        x = x[K]
        Mutations[[mycancerstudy]][,HUGOgenes[j]] = x
      }
    }
    cat("study = ", mycancerstudy, " completed.\n")
  }
}
for (i in 1:nrow(StudiesWithMutations)) {
  Mutations[[i]][Mutations[[i]] == "NaN"] = NA_character_
}

##########################################

# K = c()
# for (i in 1:nrow(Studies)) {
#   mycancerstudy = Studies[i,1]
#   mycaselist = sprintf("%s_all",mycancerstudy)
#   
#   print("########################")
#   P = getGeneticProfiles(mycgds,mycancerstudy)
#   I = which(P$genetic_profile_name == "mRNA Expression z-Scores (microarray)")
#   if (length(I) > 0) {
#     K = c(K,i)
#   }
# }
# StudiesWithMicroArray = Studies[K,]
# 
# mRNAexpressionZ = list()
# Bad = list()
# for (i in 6:nrow(StudiesWithMicroArray)) {
#   mycancerstudy = StudiesWithMicroArray[i,1]
#   mycaselist = sprintf("%s_all",mycancerstudy)
#   print("########################")
#   P = getGeneticProfiles(mycgds,mycancerstudy)
#   I = which(P$genetic_profile_name == "mRNA Expression z-Scores (microarray)")
#   if (length(I) > 0) {
#     mygeneticprofile = P[I,1]
#     x = getProfileData(mycgds,HUGOgenes[1],mygeneticprofile,mycaselist)
#     Bad[[mycancerstudy]] = c()
#     mRNAexpressionZ[[mycancerstudy]] = matrix(NA_character_,nr=nrow(x),nc=length(HUGOgenes))
#     row.names(mRNAexpressionZ[[mycancerstudy]]) = row.names(x)
#     colnames(mRNAexpressionZ[[mycancerstudy]]) = HUGOgenes
#     for (j in seq_len(length(HUGOgenes))) {
#       cat("study = ", mycancerstudy, " j=",j," nbad=",length(Bad[[mycancerstudy]]),"           \r")
#       e = try( {
#         y = getProfileData(mycgds,HUGOgenes[j],mygeneticprofile,mycaselist)
#       } )
#       k=1
#       while(class(e) == "try-error") {
#         k=k+1
#         e = try( {
#           cat("study = ", mycancerstudy, " j=",j," nbad=",length(Bad[[mycancerstudy]])," attempt=",k,"           \r")
#           Sys.sleep(time=3)
#           y = getProfileData(mycgds,HUGOgenes[j],mygeneticprofile,mycaselist)
#         } )
#       }
#       if (nrow(y) != nrow(mRNAexpressionZ[[mycancerstudy]])) {
#         Bad[[mycancerstudy]] = c(Bad[[mycancerstudy]],j)
#       } else{
#         x = as.character(y[,1])
#         names(x) = row.names(y)
#         if (length(x) != nrow(mRNAexpressionZ[[mycancerstudy]])) {
#           stop(mycancerstudy," j=",j," length of result and nrow of matrix differ")        
#         }
#         K = match(row.names(mRNAexpressionZ[[mycancerstudy]]),names(x))
#         if (any(is.na(K))) {
#           stop(mycancerstudy," j=",j," not all cases in result.")
#         }
#         x = x[K]
#         mRNAexpressionZ[[mycancerstudy]][,HUGOgenes[j]] = x
#       }
#     }
#     cat("study = ", mycancerstudy, " completed.\n")
#   }
# }

##########################################

ClinicalData = list()
for (i in 1:nrow(StudiesWithMutations)) {
  mycancerstudy = StudiesWithMutations[i,1]
  mycaselist = sprintf("%s_all",mycancerstudy)
  ClinicalData[[mycancerstudy]] = getClinicalData(mycgds, mycaselist)
}

CaseLists = list()
for (i in 1:nrow(StudiesWithMutations)) {
  mycancerstudy = StudiesWithMutations[i,1]
  CaseLists[[mycancerstudy]] = getCaseLists(mycgds,mycancerstudy)
}

CGDS = list(Mutations = Mutations, 
            CaseLists = CaseLists,
            ClinicalData = ClinicalData)
save(CGDS, file="CGDS.rda")

