# R script to download selected samples
# Copy code and run on a local machine to initiate download

#https://docs.google.com/document/d/1jm1sJtHgJkSV7hEjMC_VK7l7gzshMQzr4T7NdgY9rGE/edit#heading=h.skn982pchl79
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4930834/

# Check for dependencies and install if missing
packages <- c("rhdf5","pkgconfig")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  BiocManager::install("rhdf5")
}
#make ~/.R/Makevars  and add the following (remove # in the beginning of each line)
# LDFLAGS= -L/usr/local/clang4/lib
# FLIBS=-L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin16/6.3.0 -L/usr/local/gfortran/lib -lgfortran -lquadmath -lm
# CC=/usr/local/clang4/bin/clang
# CXX=/usr/local/gfortran/bin/g++
#   CXX1X=/usr/local/gfortran/bin/g++
#   CXX98=/usr/local/gfortran/bin/g++
#   CXX11=/usr/local/gfortran/bin/g++
#   CXX14=/usr/local/gfortran/bin/g++
#   CXX17=/usr/local/gfortran/bin/g++
#   
#   LLVM_LOC = /usr/local/opt/llvm
# CC=/usr/local/gfortran/bin/gcc -fopenmp
# CXX=/usr/local/gfortran/bin/g++ -fopenmp
# # -O3 should be faster than -O2 (default) level optimisation ..
# CFLAGS=-g -O3 -Wall -pedantic -std=gnu99 -mtune=native -pipe
# CXXFLAGS=-g -O3 -Wall -pedantic -std=c++11 -mtune=native -pipe
# LDFLAGS=-L/usr/local/opt/gettext/lib -L$(LLVM_LOC)/lib -Wl,-rpath,$(LLVM_LOC)/lib
# CPPFLAGS=-I/usr/local/opt/gettext/include -I$(LLVM_LOC)/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include

#add this path to your mac:  /usr/local/gfortran/bin
#If you have data.table error: ‘datatable.so’ not found
#https://github.com/Rdatatable/data.table/wiki/Installation#openmp-enabled-compiler-for-mac

# sample size
# platform: not array
# year
# organism


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("GEOquery", quietly = TRUE))
  install.packages("GEOquery")

library("rhdf5")
library("tools")
library("GEOquery")

setwd("~/Dropbox (Gladstone)/Pathway-Project/human_matrix_download.h5")
destination_file = "human_matrix_download.h5"

data <- getGEO("GSE67492")
data$GSE67492_series_matrix.txt.gz@experimentData@title
data$GSE67492_series_matrix.txt.gz@experimentData@abstract

gse = h5read(destination_file, "meta/Sample_series_id")
head(gse)
length(gse)
238522
length(unique(gse))
7909

unique_gse=unique(gse)



###gse issue example: "GSE27452Xx-xXGSE32303"
##removing
library(stringr)
str_count(unique_gse[3],"GSE")
GSE_counts=apply(unique_gse,1,function(x) str_count(x,"GSE"))
max(GSE_counts)
3
length(which(GSE_counts>1))
1607
unique_gse=unique_gse[-which(GSE_counts>1)]

dieases=c("Renal cell carcinoma","Alzheimer’s disease","Alzheimer","Thyroid cancer","Dilated cardiomyopathy",
          "Colorectal cancer","Prostate cancer","Huntington’s disease","Huntington","Acute Myeloid Leukemia",
          "Pancreatic cancer","Non-small cell lung cancer","Glioma","Parkinson’s disease","Parkinson",
          "Type II diabetes mellitus","Type 2 diabetes mellitus","Chronic Myeloid Leukemia","Endometrial cancer")

match_gse_disease<-function(gse_id){
  print(gse_id)
  
 
  
   tryCatch({
    
    data <- getGEO(gse_id)
    
  }, error = function(err) {
    
    
   
    return(-1)
    
  }) # END tryCatch
  
  
  title=data[[1]]@experimentData@title
  abstract=data[[1]]@experimentData@abstract
  
  matched_disease_title=apply(as.matrix(dieases),1,function(x) grep(x,title,ignore.case=TRUE) )
  matched_disease_abstract=apply(as.matrix(dieases),1,function(x) grep(x,abstract,ignore.case=TRUE) )
  
  if(length(matched_disease_title)>0 | length(matched_disease_abstract)>0){
    if(length(matched_disease_title)>0){
      which_disease=which(matched_disease_title>0)
    }else{
      which_disease=which(matched_disease_abstract>0)
    }
    return(which_disease) 
  }else{
    return(0)
  }
  if(which_disease>0){print(disease[which_disease])}
  
}

matched_diseases=apply(unique_gse,1,match_gse_disease)

matched_diseases_unlisted=as.matrix(unlist(lapply(matched_diseases, `[[`, 1)))

length(which(matched_diseases_unlisted>0))
507

tb=table(matched_diseases_unlisted)
names(tb)=c("NA",dieases)
tb

# NA       Renal cell carcinoma        Alzheimer’s disease                  Alzheimer             Thyroid cancer 
# 5795                         22                         10                         13                         13 
# Dilated cardiomyopathy          Colorectal cancer            Prostate cancer       Huntington’s disease                 Huntington 
# 4                         68                        133                          5                          9 
# Acute Myeloid Leukemia          Pancreatic cancer Non-small cell lung cancer                     Glioma        Parkinson’s disease 
# 70                         37                         27                         61                         12 
# Parkinson  Type II diabetes mellitus   Type 2 diabetes mellitus   Chronic Myeloid Leukemia         Endometrial cancer 
# 8                          1                          6                          4                          4 

#distinguish sample names

save.image(file="15disease_parsed.RData")

data <- getGEO("GSE67492")
data$GSE67492_series_matrix.txt.gz@experimentData@title
data$GSE67492_series_matrix.txt.gz@experimentData@abstract
data[[1]]@phenoData@data$title
data[[1]]@phenoData@data$description
# 
# > data[[1]]@phenoData@data$title
# V2                                  V3                                  V4                                  V5 
# Idiopathic Dilated Cardiomyopathy 1 Idiopathic Dilated Cardiomyopathy 2           PAH with BMPR2 mutation 1           PAH with BMPR2 mutation 2 
# V6                                  V7 
# Control RV 1                        Control RV 2 
# 6 Levels: Control RV 1 Control RV 2 Idiopathic Dilated Cardiomyopathy 1 Idiopathic Dilated Cardiomyopathy 2 ... PAH with BMPR2 mutation 2
# > data[[1]]@phenoData@data$description
# V2        V3        V4        V5        V6        V7 
# IDCC 1    IDCC 2       PPH       SPH Control 1 Control 2 
# Levels: Control 1 Control 2 IDCC 1 IDCC 2 PPH SPH

data <- getGEO(unique_gse[which(matched_diseases_unlisted>0)[2]])
data[[1]]@phenoData@data$title
data[[1]]@phenoData@data$description
# 
# > data[[1]]@phenoData@data$title
# V2                                                       V3 
# PSC from PDA patient [CA2]                               PSC from PDA patient [CA3] 
# V4                                                       V5 
# PSC from PDA patient [CA6]                               PSC from PDA patient [CA7] 
# V6                                                       V7 
# PSC from non-cancer control [DI1]                        PSC from non-cancer control [DI2] 
# V8                                                       V9 
# PSC from non-cancer control [DI3]                        PSC from non-cancer control [DI4] 
# V10                                                      V11 
# PSC from non-cancer control [DI5]            MiaPaCa-2 cells treated with DMEM  [Mia-DMEM] 
# V12                                                      V13 
# MiaPaCa-2 cells treated with calcipotriol  [Mia-CAL] MiaPaCa-2 cells treated with conditioned media  [Mia-CM] 
# V14 
# MiaPaCa-2 cells treated with conditioned media [Mia-CMS] 
# 13 Levels: MiaPaCa-2 cells treated with calcipotriol  [Mia-CAL] ... PSC from PDA patient [CA7]
# > data[[1]]@phenoData@data$description
# V2                                                    V3 
# processed data file available on Series record: CA.bw processed data file available on Series record: CA.bw 
# V4                                                    V5 
# processed data file available on Series record: CA.bw processed data file available on Series record: CA.bw 
# V6                                                    V7 
# processed data file available on Series record: DI.bw processed data file available on Series record: DI.bw 
# V8                                                    V9 
# processed data file available on Series record: DI.bw processed data file available on Series record: DI.bw 
# V10                                                   V11 
# processed data file available on Series record: DI.bw                     MiaPaCa-2 cells treated with DMEM 
# V12                                                   V13 
# MiaPaCa-2 cells treated with calcipotriol        MiaPaCa-2 cells treated with conditioned media 
# V14 
# MiaPaCa-2 cells treated with conditioned media 
# 5 Levels: MiaPaCa-2 cells treated with calcipotriol MiaPaCa-2 cells treated with conditioned media ... processed data file available on Series record: DI.bw




#code that retrieves samples
#https://github.com/ctlab/phantasus/blob/master/R/loadGEO.R#L708


#################sort titles



group_samples<-function(titles){

library(stringdist)
library(igraph)
library(raster)

controls=NULL
controls=c(controls,grep("control",titles,ignore.case = TRUE))
controls=c(controls,grep("wild",titles,ignore.case = TRUE))

controls=c(controls,which(regexpr("[^a-z]wt[^a-z]",titles,ignore.case = TRUE) >0))

titles[controls]

if(length(controls)==0){
  print("No control samples")
  return(list())
}

###case

distm=as.matrix(stringdist::stringdistmatrix(titles[-controls],method="lv"))
colnames(distm)=titles[-controls]
rownames(distm)=titles[-controls]
distm


# s.a <- strsplit(titles[1], "")[[1]]
# s.b <- strsplit(titles[2], "")[[1]]
# paste(s.a[s.a != s.b], collapse = "")


distm[which(distm>1)]=0

g  <- graph.adjacency(distm) 
#plot(g)
clu <- components(g)
gr=groups(clu)

gr=gr[which(lengths(gr)>1)]

if(length(gr)==0){
  print("No case samples")
  return(list())
}


sample_groups=list()

for(i in 1:length(gr)){
  #print(gr[[i]])
  sample_groups<-c(sample_groups,list(gr[[i]]))
  names(sample_groups)[i]=paste("case",i,sep="")
  last_case=i
}



###control

distm=as.matrix(stringdist::stringdistmatrix(titles[controls],method="lv"))
colnames(distm)=titles[controls]
rownames(distm)=titles[controls]
distm


distm[which(distm>1)]=0

g  <- graph.adjacency(distm) 
#plot(g)
clu <- components(g)
gr=groups(clu)

gr=gr[which(lengths(gr)>1)]

if(length(gr)==0){
  print("No control samples")
  return(list())
}

for(i in 1:length(gr)){
  #print(gr[[i]])
  sample_groups<-c(sample_groups,list(gr[[i]]))
  names(sample_groups)[last_case+i]=paste("control",i,sep="")
}

return(sample_groups)
}


#######################run the function
#result_1: successed datasets
#result_2: failed results
#results: successful results

result_1=list()
result_0=list()
results=list()
for(i in 1:length(unique_gse)){
print(i)
  data <- getGEO(unique_gse[which(matched_diseases_unlisted>0)[i]])
  titles=as.character(data[[1]]@phenoData@data$title)
  titles
  group_samples_result=group_samples(titles)
  

  if(length(group_samples_result)>0){
    results<-c(results,list(group_samples_result))
    result_1<-c(result_1,list(list(gse_id=unique_gse[which(matched_diseases_unlisted>0)[3]],titles=titles)))
    
  }else{
    
    result_0<-c(result_0,list(list(gse_id=unique_gse[which(matched_diseases_unlisted>0)[3]],titles=titles)))
  }
  
}


setwd("/Users/mingyoungshin/git/PFOCRInPathwayAnalyses")
save.image("title_parsing.RData")




#########Compare archs4 and hormonizome
harmonizome<-read.csv("~/git/PFOCRInPathwayAnalyses/target-datasets.csv")
length(which(as.character(harmonizome$dataset_id)%in%unique_gse))
0
#########

extracted_expression_file = "Hair Follicle_expression_matrix.tsv"
#url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
# Check if gene expression file was already downloaded and check integrity, if not in current directory download file form repository
# if(!file.exists(destination_file)){
#   print("Downloading compressed gene expression matrix.")
#   download.file(url, destination_file, quiet = FALSE)
# }

# Selected samples to be extracted
samp = c("GSM2072459","GSM2072458","GSM2703811","GSM2703812","GSM2703813","GSM2703814","GSM2703815","GSM2703816","GSM3717034","GSM3717035","GSM3717036","GSM3717037","GSM3717038","GSM3731086","GSM3731087","GSM3731088","GSM3731089","GSM3731090","GSM3731091","GSM3731092","GSM3731093","GSM3731094","")

# Retrieve information from compressed data
samples = h5read(destination_file, "meta/Sample_geo_accession")
tissue = h5read(destination_file, "meta/Sample_source_name_ch1")
genes = h5read(destination_file, "meta/genes")
datah5 = h5read(destination_file, "meta")
gse = h5read(destination_file, "meta/Sample_series_id")
# Identify columns to be extracted
sample_locations = which(samples %in% samp)

# extract gene expression from compressed data
expression = h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))
H5close()
rownames(expression) = genes
colnames(expression) = samples[sample_locations]

# Print file
write.table(expression, file=extracted_expression_file, sep="\t", quote=FALSE)
print(paste0("Expression file was created at ", getwd(), "/", extracted_expression_file))
