load("groups_filtered_merged.RData")

packages <- c("rhdf5")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  print("Install required packages")
  source("https://bioconductor.org/biocLite.R")
  biocLite("rhdf5")
}
library("rhdf5")
library("tools")
library('edgeR')

destination_file = "~/Dropbox (Gladstone)/Pathway-Project/human_matrix_download.h5"
extracted_expression_file = "Hair Follicle_expression_matrix.tsv"
url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"

# Check if gene expression file was already downloaded and check integrity, if not in current directory download file form repository
if(!file.exists(destination_file)){
  print("Downloading compressed gene expression matrix.")
  download.file(url, destination_file, quiet = FALSE)
}


# Retrieve information from compressed data
samples = h5read(destination_file, "meta/Sample_geo_accession")
tissue = h5read(destination_file, "meta/Sample_source_name_ch1")
genes = h5read(destination_file, "meta/genes")
series = h5read(destination_file, "meta/Sample_series_id")

start_t=Sys.time()
passed=NULL

sigGene_results=array(0,dim=c(1,length(merged_filtered_gses)))
pvalue_results=array(10,dim=c(length(genes),length(merged_filtered_gses)))
dim(pvalue_results)
rownames(pvalue_results)=genes
colnames(pvalue_results)=unlist(lapply(merged_filtered_gses , "[[" , "GSE" ))
colnames(sigGene_results)=unlist(lapply(merged_filtered_gses , "[[" , "GSE" ))

#35238   536
###########loop
for(i in 1:length(merged_filtered_gses)){

  print(i)
geo_data<-getGEO(merged_filtered_gses[[i]]$GSE,getGPL = FALSE)
merged_filtered_gses[[i]]$GSE
show(geo_data)

pheno=pData(geo_data[[1]])
dim(pheno)
colnames(pheno)
gsm_title=cbind(rownames(pheno),as.character(pheno$title))
head(gsm_title)
merged_filtered_gses[[i]]$group1

# Selected samples to be extracted
samp = gsm_title[,1]

# Identify columns to be extracted
sample_locations = which(samples %in% samp)
sample_locations

# extract gene expression from compressed data
expression = h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))
H5close()
rownames(expression) = genes
colnames(expression) = samples[sample_locations]
# Print file
#write.table(expression, file=paste0(merged_filtered_gses[[i]]$GSE,"_expression.txt"), sep="\t", quote=FALSE)


###DEG

filter <- apply(expression, 1, function(x) length(x[x>5])>=2) #minimum count (4) is satified in at least two conditions
FilterGeneCountData <- expression[filter,]


###group GSMs
Group=gsm_title[,2]
group_index=grep("group",names(merged_filtered_gses[[i]]))
selected_all=NULL
for(gi in group_index){
  selected=which( (Group%in%unlist(merged_filtered_gses[[i]][gi])) ==TRUE)
  Group[selected]=paste0("group",gi)
  selected_all=c(selected_all,selected)
}
gsm_title2=cbind(gsm_title,Group)
gsm_title2=gsm_title2[selected_all,]
colnames(gsm_title2)=c("GSM","Pheno","Group")

gsm_title2=as.data.frame(gsm_title2)
gsm_title2$Group <- as.factor(gsm_title2$Group)

TempIndices <- match(gsm_title2[,1], colnames(FilterGeneCountData))

colN=colnames(FilterGeneCountData)
if(length(gsm_title2[,1])!=length(TempIndices[!is.na(TempIndices)])){
  print("sample size mismatch in pheno data and expression data")
  passed=c(passed,i)
  next
}

FilterGeneCountData=FilterGeneCountData[,TempIndices]
cbind(gsm_title2,colnames(FilterGeneCountData))

Group_p <- gsm_title2$Group
print(levels(Group_p))

design <- model.matrix(~Group_p)
design

y <- DGEList(counts=FilterGeneCountData, group=NULL) 
y <- calcNormFactors(y)

y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
Coef <- fit$coefficients
head(Coef)

lrt <- glmQLFTest(fit, coef=c(2))

pvalue_results[match(rownames(lrt),rownames(pvalue_results)),merged_filtered_gses[[i]]$GSE]=lrt$table$PValue
sigGene_results[1,merged_filtered_gses[[i]]$GSE]=sum(p.adjust(lrt$table$PValue,method="BH") < 0.05)

}

Sys.time()

#i5 gsm don't exist in h5
#i39 GSM1895372 don't exist in h5
#i45 sample size mismatch
#i46 sample size mismatch
#########end of loop
sigGene_results_bk=sigGene_results
pvalue_results_bk=pvalue_results
length(passed)
54

sigGene_results[passed]=NA
pvalue_results[which(pvalue_results==10)]=NA

jpeg("hist_number_sig_genes_human.jpeg")
hist(sigGene_results,breaks=50)
dev.off()

table(sigGene_results)

colnames(sigGene_results)[which(sigGene_results==424)]
save.image("human_DEG_result.RData")



geo_data<-getGEO(merged_filtered_gses[[i]]$GSE,getGPL = FALSE)
merged_filtered_gses[[i]]$GSE
show(geo_data)

pheno=pData(geo_data[[1]])
dim(pheno)
colnames(pheno)
gsm_title=cbind(rownames(pheno),as.character(pheno$title))
head(gsm_title)
merged_filtered_gses[[i]]$group1

# Selected samples to be extracted
samp = gsm_title[,1]

# Identify columns to be extracted
sample_locations = which(samples %in% samp)
sample_locations

# extract gene expression from compressed data
expression = h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))
H5close()
rownames(expression) = genes
colnames(expression) = samples[sample_locations]
# Print file
#write.table(expression, file=paste0(merged_filtered_gses[[i]]$GSE,"_expression.txt"), sep="\t", quote=FALSE)


###DEG

filter <- apply(expression, 1, function(x) length(x[x>5])>=2) #minimum count (4) is satified in at least two conditions
FilterGeneCountData <- expression[filter,]


###group GSMs
Group=gsm_title[,2]
group_index=grep("group",names(merged_filtered_gses[[i]]))
selected_all=NULL
for(gi in group_index){
  selected=which( (Group%in%unlist(merged_filtered_gses[[i]][gi])) ==TRUE)
  Group[selected]=paste0("group",gi)
  selected_all=c(selected_all,selected)
}
gsm_title2=cbind(gsm_title,Group)
gsm_title2=gsm_title2[selected_all,]
colnames(gsm_title2)=c("GSM","Pheno","Group")

gsm_title2=as.data.frame(gsm_title2)
gsm_title2$Group <- as.factor(gsm_title2$Group)

TempIndices <- match(gsm_title2[,1], colnames(FilterGeneCountData))

colN=colnames(FilterGeneCountData)
if(length(gsm_title2[,1])!=length(TempIndices[!is.na(TempIndices)])){
  print("sample size mismatch in pheno data and expression data")
  passed=c(passed,i)
  next
}

FilterGeneCountData=FilterGeneCountData[,TempIndices]
cbind(gsm_title2,colnames(FilterGeneCountData))

Group_p <- gsm_title2$Group
print(levels(Group_p))

design <- model.matrix(~Group_p)
design

y <- DGEList(counts=FilterGeneCountData, group=NULL) 
y <- calcNormFactors(y)

y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
Coef <- fit$coefficients
head(Coef)

lrt <- glmQLFTest(fit, coef=c(2))

pvalue_results[match(rownames(lrt),rownames(pvalue_results)),merged_filtered_gses[[i]]$GSE]=lrt$table$PValue
sigGene_results[1,merged_filtered_gses[[i]]$GSE]=sum(p.adjust(lrt$table$PValue,method="BH") < 0.05)

#######################1000 permutation


above_threshold=colnames(sigGene_results)[which(sigGene_results<7000)]
383
group1_counts=apply(as.matrix(1:length(merged_filtered_gses)),1,function(i) length(merged_filtered_gses[[i]]$group1))
hist(group1_counts,breaks = 50)  
head(sort(group1_counts,decreasing = TRUE))
which(group1_counts==39)
128
497
merged_filtered_gses[[128]]
merged_filtered_gses[[497]]

i=128
perm_test_set1=pvalue_results[,i]
colnames(pvalue_results)[i]
#"GSE102741"
sigGene_results[i]
4


pvalue_results_perm1=array(10,dim=c(length(genes),1000))
dim(pvalue_results_perm1)
rownames(pvalue_results_perm1)=genes
colnames(pvalue_results_perm1)=paste0("set",1:1000)
sigGene_results_perm1=array(0,dim=c(length(merged_filtered_gses)))



geo_data<-getGEO("GSE102741",getGPL = FALSE)
show(geo_data)

pheno=pData(geo_data[[1]])
dim(pheno)
colnames(pheno)
gsm_title=cbind(rownames(pheno),as.character(pheno$title))
head(gsm_title)
merged_filtered_gses[[i]]

# Selected samples to be extracted
samp = gsm_title[,1]

# Identify columns to be extracted
sample_locations = which(samples %in% samp)
sample_locations

# extract gene expression from compressed data
expression = h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))
H5close()
rownames(expression) = genes
colnames(expression) = samples[sample_locations]

start_t=Sys.time()

set.seed(511)
for(k in 1:1000){
  
  ######permute samples
  
  expression_perm1=expression[,sample(1:ncol(expression),ncol(expression),replace=FALSE)]
  colnames(expression_perm1)=colnames(expression)
  ###DEG
  
  filter <- apply(expression_perm1, 1, function(x) length(x[x>5])>=2) #minimum count (4) is satified in at least two conditions
  FilterGeneCountData <- expression_perm1[filter,]
  
  
  ###group GSMs
  Group=gsm_title[,2]
  group_index=grep("group",names(merged_filtered_gses[[i]]))
  selected_all=NULL
  for(gi in group_index){
    selected=which( (Group%in%unlist(merged_filtered_gses[[i]][gi])) ==TRUE)
    Group[selected]=paste0("group",gi)
    selected_all=c(selected_all,selected)
  }
  gsm_title2=cbind(gsm_title,Group)
  gsm_title2=gsm_title2[selected_all,]
  colnames(gsm_title2)=c("GSM","Pheno","Group")
  
  gsm_title2=as.data.frame(gsm_title2)
  gsm_title2$Group <- as.factor(gsm_title2$Group)
  
  TempIndices <- match(gsm_title2[,1], colnames(FilterGeneCountData))
  
  colN=colnames(FilterGeneCountData)
  
  FilterGeneCountData=FilterGeneCountData[,TempIndices]
  cbind(gsm_title2,colnames(FilterGeneCountData))
  
  Group_p <- gsm_title2$Group
  print(levels(Group_p))
  
  design <- model.matrix(~Group_p)
  design
  
  y <- DGEList(counts=FilterGeneCountData, group=NULL) 
  y <- calcNormFactors(y)
  
  y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmQLFit(y, design, robust=TRUE)
  Coef <- fit$coefficients
  head(Coef)
  
  lrt <- glmQLFTest(fit, coef=c(2))
  
  pvalue_results_perm1[match(rownames(lrt),rownames(pvalue_results_perm1)),k]=lrt$table$PValue
  sigGene_results_perm1[k]=sum(p.adjust(lrt$table$PValue,method="BH") < 0.05)
  
  print(k)
  
}

end_t=Sys.time()#started at 5:55 pm






i=497
perm_test_set2=pvalue_results[,i]
colnames(pvalue_results)[i]
#"GSE139061"
sigGene_results[i]
4


pvalue_results_perm2=array(10,dim=c(length(genes),1000))
dim(pvalue_results_perm2)
rownames(pvalue_results_perm2)=genes
colnames(pvalue_results_perm2)=paste0("set",1:1000)
sigGene_results_perm2=array(0,dim=c(length(merged_filtered_gses)))



geo_data<-getGEO("GSE139061",getGPL = FALSE)
show(geo_data)

pheno=pData(geo_data[[1]])
dim(pheno)
colnames(pheno)
gsm_title=cbind(rownames(pheno),as.character(pheno$title))
head(gsm_title)
merged_filtered_gses[[i]]

# Selected samples to be extracted
samp = gsm_title[,1]

# Identify columns to be extracted
sample_locations = which(samples %in% samp)
sample_locations

# extract gene expression from compressed data
expression = h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))
H5close()
rownames(expression) = genes
colnames(expression) = samples[sample_locations]

start_t=Sys.time()

set.seed(511)
for(k in 1:1000){
  
  ######permute samples
  
  expression_perm1=expression[,sample(1:ncol(expression),ncol(expression),replace=FALSE)]
  colnames(expression_perm1)=colnames(expression)
  ###DEG
  
  filter <- apply(expression_perm1, 1, function(x) length(x[x>5])>=2) #minimum count (4) is satified in at least two conditions
  FilterGeneCountData <- expression_perm1[filter,]
  
  
  ###group GSMs
  Group=gsm_title[,2]
  group_index=grep("group",names(merged_filtered_gses[[i]]))
  selected_all=NULL
  for(gi in group_index){
    selected=which( (Group%in%unlist(merged_filtered_gses[[i]][gi])) ==TRUE)
    Group[selected]=paste0("group",gi)
    selected_all=c(selected_all,selected)
  }
  gsm_title2=cbind(gsm_title,Group)
  gsm_title2=gsm_title2[selected_all,]
  colnames(gsm_title2)=c("GSM","Pheno","Group")
  
  gsm_title2=as.data.frame(gsm_title2)
  gsm_title2$Group <- as.factor(gsm_title2$Group)
  
  TempIndices <- match(gsm_title2[,1], colnames(FilterGeneCountData))
  
  colN=colnames(FilterGeneCountData)
  
  FilterGeneCountData=FilterGeneCountData[,TempIndices]
  cbind(gsm_title2,colnames(FilterGeneCountData))
  
  Group_p <- gsm_title2$Group
  print(levels(Group_p))
  
  design <- model.matrix(~Group_p)
  design
  
  y <- DGEList(counts=FilterGeneCountData, group=NULL) 
  y <- calcNormFactors(y)
  
  y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmQLFit(y, design, robust=TRUE)
  Coef <- fit$coefficients
  head(Coef)
  
  lrt <- glmQLFTest(fit, coef=c(2))
  
  pvalue_results_perm2[match(rownames(lrt),rownames(pvalue_results_perm2)),k]=lrt$table$PValue
  sigGene_results_perm2[k]=sum(p.adjust(lrt$table$PValue,method="BH") < 0.05)
  
  print(k)
  
}

end_t=Sys.time()#started at 5:55 pm
save.image("two_permutation_tests.RData")
###################################################################################################################rSEA
library(rSEA)
library("clusterProfiler")
library("org.Hs.eg.db")
library(dplyr)

wp <- read.gmt("wikipathways-20200210-gmt-Homo_sapiens.gmt")
wp <- wp %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wp=wp[,c("name","wpid","gene")]
wp=unique(wp)

wp_list=split(wp$gene,wp$wpid)


   gene_entrez <- bitr(names(perm_test_set1), fromType = "SYMBOL",
                      toType = c("ENTREZID"),
                      OrgDb = org.Hs.eg.db, drop=FALSE)
  head(gene_entrez)
  
  colnames(gene_entrez)[1]="Gene"
  gene_entrez=unique(gene_entrez)
  
  
  gene_entrez$Gene[setdiff(1:length(names(perm_test_set1)),which(gene_entrez$Gene%in%names(perm_test_set1)))]
  tb=table(gene_entrez$Gene)
  tb[which(tb>1)]
  
  # HBD MEMO1  MMD2   TEC 
  # 2     2     2     2 
  
  which(gene_entrez$Gene=="HBD")
  which(gene_entrez$Gene=="MEMO1")
  which(gene_entrez$Gene=="MMD2")
  which(gene_entrez$Gene=="TEC")
  
  gene_entrez=gene_entrez[-c(11401,15467,15707, 31349),]
  
  perm_test_set1_m=as.data.frame(cbind(names(perm_test_set1),as.matrix(perm_test_set1)))
  colnames(perm_test_set1_m)=c("Gene","pvalue")
  merged=merge(perm_test_set1_m,gene_entrez,by="Gene",all = "L",all.x=TRUE,all.y=FALSE)
  
  head(merged)
  
  sea_result<-SEA(head(as.numeric(as.character(merged$pvalue))), merged$ENTREZID, pathlist = wp_list)
  colnames(sea_result)[2]="wpid"
  head(sea_result)
  wiki_result=merge(sea_result,unique(wp[,-3]),by="wpid")
  head(wiki_result)
  
  dim(sea_result)
  dim(wp)
  sum( sea_result$wpid%in% wp$wpid)
  
  wiki_result_sorted<-topSEA(wiki_result, by=Comp.adjP, descending=FALSE, n=length(which(wiki_result$Comp.adjP<0.2))) #sorted by smallest competitive adj.p-values 
  head(wiki_result_sorted,50)
  
  write.xlsx(wiki_result_sorted,file=paste(substring(fileNs[filei],1,nchar(fileNs[filei])-3),"_pathwayEnrichment_newWiki.xlsx",sep=""))
  

# 
# 
# 
# 
# library(GEOquery)
# library(R.utils)
# i=6
# 
# 
# geo_data<-getGEO(merged_filtered_gses[[i]]$GSE,getGPL = FALSE)
# show(geo_data)
# 
# if(!dir.exists(merged_filtered_gses[[i]]$GSE)){
#   filePaths=getGEOSuppFiles(merged_filtered_gses[[i]]$GSE)
# }
# 
# 
# pheno=pData(geo_data[[1]])
# dim(pheno)
# colnames(pheno)
# gsm_title=cbind(rownames(pheno),as.character(pheno$title))
# head(gsm_title)
# merged_filtered_gses[[i]]$group1
# 
# fileN=list.files(path=paste0("./",merged_filtered_gses[[i]]$GSE) ,pattern="*.txt$")
# 
# if(length(fileN)>0){
#   data<-read.table(paste0("./",merged_filtered_gses[[i]]$GSE,"/",fileN),header=TRUE)
# }else{
#   
# fileN=list.files(path=paste0("./",merged_filtered_gses[[i]]$GSE) ,pattern="*.xlsx$")
# 
# if(length(fileN)>0){
#   data<-read_excel(paste0("./",merged_filtered_gses[[i]]$GSE,"/",fileN))
# }else{
#   
#   fileN=list.files(path=paste0("./",merged_filtered_gses[[i]]$GSE) ,pattern="*.csv$")
#   
# if(length(fileN)>0){
#   data<-read.csv(paste0("./",merged_filtered_gses[[i]]$GSE,"/",fileN))
# }else{
#   fileN=list.files(path=paste0("./",merged_filtered_gses[[i]]$GSE) ,pattern="*.gz$")
#   
# 
#   if(length(fileN)>0){
#     if(length(grep("txt",fileN)>0)){
#       data<-read.table(gzfile(paste0("./",merged_filtered_gses[[i]]$GSE,"/",fileN),'rt'),header=TRUE)
#     }else if(length(grep("xlsx",fileN)>0)){
#       R.utils:gunzip(paste0("./",merged_filtered_gses[[i]]$GSE,"/",fileN))
# 
#       fileN=list.files(path=paste0("./",merged_filtered_gses[[i]]$GSE) ,pattern="*.xlsx$")
#       data<-read_excel(paste0("./",merged_filtered_gses[[i]]$GSE,"/",fileN))
#     }else if(length(grep("csv",fileN)>0)){
#       data<-read.csv(gzfile(paste0("./",merged_filtered_gses[[i]]$GSE,"/",fileN),'rt'))
#     }
#     
#   }else{
#     print("Cannot read")
#     print(merged_filtered_gses[[i]]$GSE)
#     print(fileN)
#     print("is not readable")
#   }
#   
#   
#   
# }
# }
# }
# 
# fileN
# dim(data)
# colnames(data)
# head(gsm_title)
# 
# data<-read.csv(gzfile(paste0("./",merged_filtered_gses[[i]]$GSE,"/",fileN),'rt'),header=TRUE,sep="\t")
# 
# 
# 
# 
# library(readxl)
# untar(paste0("./",merged_filtered_gses[[i]]$GSE,"/",fileN))
# data<-read_excel(gzfile(paste0("./",merged_filtered_gses[[i]]$GSE,"/",fileN),'rt'))
#                  
# library(data.table)
# dt = fread(paste0("./",merged_filtered_gses[[i]]$GSE,"/",fileN))
# 
# 
# 
# 
# colnames(data)
