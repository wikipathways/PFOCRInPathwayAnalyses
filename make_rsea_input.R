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
passed=NULL ## GSEs with errors

sigGene_results=array(0,dim=c(1,length(merged_filtered_gses)))
pvalue_results=array(10,dim=c(length(genes),length(merged_filtered_gses)))
logFC_results=array(10,dim=c(length(genes),length(merged_filtered_gses)))

dim(pvalue_results)
rownames(pvalue_results)=genes
colnames(pvalue_results)=unlist(lapply(merged_filtered_gses , "[[" , "GSE" ))
colnames(sigGene_results)=colnames(pvalue_results)
colnames(logFC_results)=colnames(pvalue_results)

#gene count / GSE count
#35238   536


##################################################################################################run DEG analysis for all GSE
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
logFC_results[match(rownames(lrt),rownames(pvalue_results)),merged_filtered_gses[[i]]$GSE]=lrt$table$logFC
}

Sys.time()

#i5 gsm don't exist in h5
#i39 GSM1895372 don't exist in h5
#i45 sample size mismatch
#i46 sample size mismatch
##################################################################################################'run DEG analysis for all GSE' ends here' 

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

save.image("human_DEG_result.RData")

##################################################################################################'1000 permutation

above_threshold=colnames(sigGene_results)[which(sigGene_results<7000)]
383
group1_counts=apply(as.matrix(1:length(merged_filtered_gses)),1,function(i) length(merged_filtered_gses[[i]]$group1))
hist(group1_counts,breaks = 50)  
head(sort(group1_counts,decreasing = TRUE))
which(group1_counts==39)
sigGene_results[which(group1_counts==39)]
128#just 4 sigificant genes
497
merged_filtered_gses[[497]]


which(group1_counts==29)
353
sigGene_results[which(group1_counts==29)]
601
merged_filtered_gses[[353]]


##################################################################################################'data 1


i=353
perm_test_set1=as.matrix(pvalue_results[,i])
colnames(pvalue_results)[i]
#"GSE107868"

perm_test_set1[which(perm_test_set1==10)]=NA

pvalue_results_perm1=array(10,dim=c(length(genes),1000))
dim(pvalue_results_perm1)
rownames(pvalue_results_perm1)=genes
colnames(pvalue_results_perm1)=paste0("set",1:1000)
sigGene_results_perm1=array(0,dim=c(length(merged_filtered_gses)))
logFC_results_perm1=array(0,dim=c(length(merged_filtered_gses)))

geo_data<-getGEO("GSE107868",getGPL = FALSE)
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
  gsm_title1=cbind(gsm_title,Group)
  gsm_title1=gsm_title1[selected_all,]
  colnames(gsm_title1)=c("GSM","Pheno","Group")
  
  gsm_title1=as.data.frame(gsm_title1)
  gsm_title1$Group <- as.factor(gsm_title1$Group)
  
  TempIndices <- match(gsm_title1[,1], colnames(FilterGeneCountData))
  
  colN=colnames(FilterGeneCountData)
  
  FilterGeneCountData=FilterGeneCountData[,TempIndices]
  cbind(gsm_title1,colnames(FilterGeneCountData))
  
  Group_p <- gsm_title1$Group
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
  logFC_results_perm1[match(rownames(lrt),rownames(pvalue_results_perm1)),k]=lrt$table$logFC
  
  print(k)
  
}

pvalue_results_perm1[which(pvalue_results_perm1==10)]=NA

end_t=Sys.time()#started at 5:55 pm
save.image("two_permutation_tests.RData")



##################################################################################################'data 2

i=497
perm_test_set2=as.matrix(pvalue_results[,i])
colnames(pvalue_results)[i]
#"GSE139061"

perm_test_set2[which(perm_test_set2==10)]=NA

pvalue_results_perm2=array(10,dim=c(length(genes),1000))
dim(pvalue_results_perm2)
rownames(pvalue_results_perm2)=genes
colnames(pvalue_results_perm2)=paste0("set",1:1000)
sigGene_results_perm2=array(0,dim=c(length(merged_filtered_gses)))
logFC_results_perm2=array(0,dim=c(length(merged_filtered_gses)))

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
  logFC_results_perm2[match(rownames(lrt),rownames(pvalue_results_perm2)),k]=lrt$table$logFC
  
  print(k)
  
}

pvalue_results_perm2[which(pvalue_results_perm2==10)]=NA

end_t=Sys.time()#started at 5:55 pm
save.image("two_permutation_tests.RData")

##################data 2
# 
# 
# i=497
# perm_test_set2=as.matrix(pvalue_results[,i])
# colnames(pvalue_results)[i]
# #"GSE139061"
# sigGene_results[i]
# 4
# perm_test_set2[which(perm_test_set2==10)]=NA
# 
# 
# pvalue_results_perm2=array(10,dim=c(length(genes),1000))
# dim(pvalue_results_perm2)
# rownames(pvalue_results_perm2)=genes
# colnames(pvalue_results_perm2)=paste0("set",1:1000)
# sigGene_results_perm2=array(0,dim=c(length(merged_filtered_gses)))
# 
# 
# 
# geo_data<-getGEO("GSE139061",getGPL = FALSE)
# show(geo_data)
# 
# pheno=pData(geo_data[[1]])
# dim(pheno)
# colnames(pheno)
# gsm_title=cbind(rownames(pheno),as.character(pheno$title))
# head(gsm_title)
# merged_filtered_gses[[i]]
# 
# # Selected samples to be extracted
# samp = gsm_title[,1]
# 
# # Identify columns to be extracted
# sample_locations = which(samples %in% samp)
# sample_locations
# 
# # extract gene expression from compressed data
# expression = h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))
# H5close()
# rownames(expression) = genes
# colnames(expression) = samples[sample_locations]
# 
# start_t=Sys.time()
# 
# set.seed(511)
# for(k in 1:1000){
#   
#   ######permute samples
#   
#   expression_perm1=expression[,sample(1:ncol(expression),ncol(expression),replace=FALSE)]
#   colnames(expression_perm1)=colnames(expression)
#   ###DEG
#   
#   filter <- apply(expression_perm1, 1, function(x) length(x[x>5])>=2) #minimum count (4) is satified in at least two conditions
#   FilterGeneCountData <- expression_perm1[filter,]
#   
#   
#   ###group GSMs
#   Group=gsm_title[,2]
#   group_index=grep("group",names(merged_filtered_gses[[i]]))
#   selected_all=NULL
#   for(gi in group_index){
#     selected=which( (Group%in%unlist(merged_filtered_gses[[i]][gi])) ==TRUE)
#     Group[selected]=paste0("group",gi)
#     selected_all=c(selected_all,selected)
#   }
#   gsm_title2=cbind(gsm_title,Group)
#   gsm_title2=gsm_title2[selected_all,]
#   colnames(gsm_title2)=c("GSM","Pheno","Group")
#   
#   gsm_title2=as.data.frame(gsm_title2)
#   gsm_title2$Group <- as.factor(gsm_title2$Group)
#   
#   TempIndices <- match(gsm_title2[,1], colnames(FilterGeneCountData))
#   
#   colN=colnames(FilterGeneCountData)
#   
#   FilterGeneCountData=FilterGeneCountData[,TempIndices]
#   cbind(gsm_title2,colnames(FilterGeneCountData))
#   
#   Group_p <- gsm_title2$Group
#   print(levels(Group_p))
#   
#   design <- model.matrix(~Group_p)
#   design
#   
#   y <- DGEList(counts=FilterGeneCountData, group=NULL) 
#   y <- calcNormFactors(y)
#   
#   y <- estimateDisp(y, design, robust=TRUE)
#   fit <- glmQLFit(y, design, robust=TRUE)
#   Coef <- fit$coefficients
#   head(Coef)
#   
#   lrt <- glmQLFTest(fit, coef=c(2))
#   
#   pvalue_results_perm2[match(rownames(lrt),rownames(pvalue_results_perm2)),k]=lrt$table$PValue
#   sigGene_results_perm2[k]=sum(p.adjust(lrt$table$PValue,method="BH") < 0.05)
#   
#   print(k)
#   
# }
# 
# pvalue_results_perm2[which(pvalue_results_perm2==10)]=NA
# 
# end_t=Sys.time()#started at 5:55 pm
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
  

  
  
  
  ###################################################################################################################run rSEA function  
run_rSEA<-function(data,list_result){
  
  na_row=which(is.na(data)==TRUE)
  data=as.matrix(data[-na_row,1])
  
  data_m=as.data.frame(cbind(rownames(data),data))
  colnames(data_m)=c("Gene","pvalue")
  merged=merge(data_m,gene_entrez,by="Gene",all = "L",all.x=TRUE,all.y=FALSE)
  
 # head(merged)
  
  sea_result<-SEA(as.numeric(as.character(merged$pvalue)), merged$ENTREZID, pathlist = wp_list)
  colnames(sea_result)[2]="wpid"
  #head(sea_result)
  wiki_result=merge(sea_result,unique(wp[,-3]),by="wpid")
 # head(wiki_result)
  

  
  list_result=modifyList(list_result, list(wpid = cbind(list_result$wpid,wiki_result$wpid), Coverage = cbind(list_result$Coverage,wiki_result$Coverage) ,
                         TDP.bound = cbind(list_result$TDP.bound,wiki_result$TDP.bound), TDP.estimate = cbind(list_result$TDP.estimate,wiki_result$TDP.estimate),  
                         SC.adjP = cbind(list_result$SC.adjP,wiki_result$SC.adjP), Comp.adjP = cbind(list_result$Comp.adjP,wiki_result$Comp.adjP),
                         ID = cbind(list_result$ID), Size = cbind(list_result$Size), name=list_result$name) )
  
  return(list_result)
                         
  
}

perm_test_set1_list <- vector("list", 9)
names(perm_test_set1_list) <- c("wpid", "ID", "Size", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP", "name")


perm_test_set2_list <- vector("list", 9)
names(perm_test_set2_list) <- c("wpid", "ID", "Size", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP", "name")

start_t=Sys.time()
perm_test_set2_list=apply(as.matrix(1:ncol(pvalue_results_perm2)),1,function(x) run_rSEA(as.matrix(pvalue_results_perm2[,x]),perm_test_set2_list))
end_t=Sys.time()

save.image("../PFOCRInPathwayAnalyses_RData/permuted_rsea.RData")

start_t=Sys.time()
perm_test_set1_list=apply(as.matrix(1:ncol(pvalue_results_perm1)),1,function(x) run_rSEA(as.matrix(pvalue_results_perm1[,x]),perm_test_set1_list))
end_t=Sys.time()

save.image("../PFOCRInPathwayAnalyses_RData/permuted_rsea.RData")

head(perm2_rsea_result$wpid)
head(perm2_rsea_result$Comp.adjP)

###################################################################################################################run ORA function#20 per GSE
run_ORA<-function(data, list_result){
  
  
  GeneRatio=matrix(NA,nrow=length(unique(wp$wpid)),ncol=1)
  BgRatio=matrix(NA,nrow=length(unique(wp$wpid)),ncol=1)
  pvalue=matrix(NA,nrow=length(unique(wp$wpid)),ncol=1)
  p.adjust=matrix(NA,nrow=length(unique(wp$wpid)),ncol=1)
  qvalue=matrix(NA,nrow=length(unique(wp$wpid)),ncol=1)
  Count=matrix(NA,nrow=length(unique(wp$wpid)),ncol=1)
  rownames(GeneRatio)=unique(wp$wpid)
  rownames(BgRatio)=unique(wp$wpid)
  rownames(pvalue)=unique(wp$wpid)
  rownames(p.adjust)=unique(wp$wpid)
  rownames(qvalue)=unique(wp$wpid)
  rownames(Count)=unique(wp$wpid)
  
  na_row=which(is.na(data)==TRUE)
  data=as.matrix(data[-na_row,1])
  
  data_m=as.data.frame(cbind(rownames(data),data))
  colnames(data_m)=c("Gene","pvalue")
  merged=merge(data_m,gene_entrez,by="Gene",all = "L",all.x=TRUE,all.y=FALSE)
  

  
  enrichment_result <- clusterProfiler::enricher(
    merged$ENTREZID,
    #universe = all genes,
    pAdjustMethod = "fdr",
    pvalueCutoff = 1, #p.adjust cutoff
    # minGSSize = 1,
    # maxGSSize = 200,
    TERM2GENE = wp[,2:3],
    TERM2NAME = wp[,c(2,1)])
  
  #enrichment_result <- DOSE::setReadable(enrichment_result, org.Hs.eg.db, keyType = "ENTREZID")
  
  enrichment_result=as.data.frame(enrichment_result)
  enrichment_result=enrichment_result[,-match("geneID",colnames(enrichment_result))]
  
  GeneRatio[match(enrichment_result$ID,rownames(GeneRatio))]=enrichment_result$GeneRatio
  BgRatio[match(enrichment_result$ID,rownames(GeneRatio))]=enrichment_result$BgRatio
  pvalue[match(enrichment_result$ID,rownames(GeneRatio))]=enrichment_result$pvalue
  p.adjust[match(enrichment_result$ID,rownames(GeneRatio))]=enrichment_result$p.adjust
  qvalue[match(enrichment_result$ID,rownames(GeneRatio))]=enrichment_result$qvalue
  Count[match(enrichment_result$ID,rownames(GeneRatio))]=enrichment_result$Count
  
  list_result=modifyList(list_result, list(GeneRatio = cbind(list_result$GeneRatio,GeneRatio), BgRatio = cbind(list_result$BgRatio,BgRatio) ,
                                           pvalue = cbind(list_result$pvalue,pvalue), p.adjust = cbind(list_result$p.adjust,p.adjust),  
                                           qvalue = cbind(list_result$qvalue,qvalue), Count = cbind(list_result$Count,Count)) )
  
  return(list_result)
  

  
 
}




perm_test_set2_list_ora <- vector("list", 6)
names(perm_test_set2_list_ora) <- c("GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue","Count")

start_t=Sys.time()
perm_test_set2_list_ora=apply(as.matrix(1:ncol(pvalue_results_perm2)),1,function(x) run_ORA(as.matrix(pvalue_results_perm2[,x]),perm_test_set2_list_ora))
end_t=Sys.time()

###################################################################################################################run GSEA function
run_gsea<-function(data){
  
  ### geneList prep
  h1d.data.hs$rank<-with(h1d.data.hs, sign(h1d.data.hs$log2FC.GL05520) * - log10(h1d.data.hs$P.Value))
  h1d.genelist.hs<-h1d.data.hs$rank
  names(h1d.genelist.hs)<-h1d.data.hs$entrez.hs
  h1d.genelist.hs = sort(h1d.genelist.hs[unique(names(h1d.genelist.hs))], decreasing = TRUE)

    gsewp.p <- GSEA(
    h1d.genelist.hs,
    TERM2GENE = wp[,2:3],
    TERM2NAME = wp[,c(2,1)],
    nPerm        = 1000,
    minGSSize    = 10,
    maxGSSize    = 500,
    pvalueCutoff = 0.5,
    verbose=FALSE)
    
  gsewp.p <- DOSE::setReadable(gsewp.p, org.Hs.eg.db, keyType = "ENTREZID")
  head(gsewp.p, 20)


  
}