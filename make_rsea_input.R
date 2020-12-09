load("groups_filtered_merged.RData")

packages <- c("rhdf5")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  print("Install required packages")
  source("https://bioconductor.org/biocLite.R")
  biocLite("rhdf5")
}

require("rhdf5")
require("tools")
require('edgeR')

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

##################################################################################################GSE key word match matrix 
require(GEOquery)
gse_keyword_A=matrix(NA,nrow=length(merged_filtered_gses),ncol=5)

for(i in 1:length(merged_filtered_gses)){
  
  print(i)
  gse_keyword_A[i,1]=merged_filtered_gses[[i]]$GSE

  geo_data<-getGEO(merged_filtered_gses[[i]]$GSE,getGPL = FALSE)
  description=geo_data[[1]]@experimentData@title
  abstract=geo_data[[1]]@experimentData@abstract
  overall_design=geo_data[[1]]@experimentData@other$overall_design
  
  gse_keyword_A[i,2]=description
  gse_keyword_A[i,3]=abstract
  gse_keyword_A[i,5]=overall_design
  if(length(merged_filtered_gses[[i]]$matched_term)>0) gse_keyword_A[i,4]=paste(merged_filtered_gses[[i]]$matched_term,collapse = '/')
 
  
}
gse_keyword_A_bk=gse_keyword_A
gse_keyword_A[,3]=gsub("\t"," ",gse_keyword_A[,3])
gse_keyword_A[,3]=gsub("\n"," ",gse_keyword_A[,3])

gse_keyword_A[,5]=gsub("\t"," ",gse_keyword_A[,3])
gse_keyword_A[,5]=gsub("\n"," ",gse_keyword_A[,3])
write.table(gse_keyword_A,file="GSE_keyword_description_human.txt",col.names=c("GSE","description","abstract","keyword","overall_design"),row.names=FALSE,sep="\t",quote=FALSE)

##########2 group only gse & gse sample names
gse_keyword_A2=matrix(NA,nrow=length(merged_filtered_gses),ncol=4)

for(i in 1:length(merged_filtered_gses)){
  
  gse_keyword_A2[i,1]=merged_filtered_gses[[i]]$GSE
  gse_keyword_A2[i,2]=merged_filtered_gses[[i]]$failed_sample_count
  gse_keyword_A2[i,3]=paste(merged_filtered_gses[[i]]$group1,collapse = '/')
  gse_keyword_A2[i,4]=paste(merged_filtered_gses[[i]]$group2,collapse = '/')
  
}


write.table(gse_keyword_A2,file="GSE_keyword_description_human_2groupOnly.txt",col.names=c("GSE","failed_sample_count","group1","group2"),row.names=FALSE,sep="\t",quote=FALSE)
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
Coef <- fit$coefficientsf
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
logFC_results[which(pvalue_results==10)]=NA
pvalue_results[which(pvalue_results==10)]=NA


jpeg("hist_number_sig_genes_human.jpeg")
hist(sigGene_results,breaks=50)
dev.off()

table(sigGene_results)

save.image("human_DEG_result.RData")
##################################################################################################run voom DEG analysis for all GSE

start_t=Sys.time()
pvalue_results_human_voom=array(NA,dim=c(length(genes),length(merged_filtered_gses)))
dim(pvalue_results_human_voom)
rownames(pvalue_results_human_voom)=genes
colnames(pvalue_results_human_voom)=unlist(lapply(merged_filtered_gses , "[[" , "GSE" ))
sigGene_results_human_voom=array(NA,dim=c(length(merged_filtered_gses)))
logFC_results_human_voom=array(NA,dim=c(length(genes),length(merged_filtered_gses)))
rownames(logFC_results_human_voom)=genes
colnames(logFC_results_human_voom)=unlist(lapply(merged_filtered_gses , "[[" , "GSE" ))
colnames(sigGene_results_human_voom)=colnames(pvalue_results_human_voom)


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
  
  design <- model.matrix(~0+Group_p)
  design
  
  y <- DGEList(counts=FilterGeneCountData, group=as.character(gsm_title2$Group))
  y <- calcNormFactors(y)
  y <- estimateDisp(y)
  
  y <- voom(y, design, plot = F)
  
  fit <- lmFit(y, design)
  head(coef(fit))
  
  contr <- makeContrasts(Group_pgroup2-Group_pgroup1, levels = colnames(coef(fit)))
  contr
  
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  voom_result=topTable(tmp, sort.by = "P", n = Inf)
  
  pvalue_results_human_voom[match(rownames(voom_result),rownames(pvalue_results_human_voom)),merged_filtered_gses[[i]]$GSE]=voom_result$P.Value
  sigGene_results_human_voom[merged_filtered_gses[[i]]$GSE]=sum(voom_result$adj.P.Val < 0.05)
  logFC_results_human_voom[match(rownames(voom_result),rownames(logFC_results_human_voom)),merged_filtered_gses[[i]]$GSE]=voom_result$logFC
  
}

end_t=Sys.time()
save.image("human_DEG_result.RData")

sigGene_results_human_voom=array(NA,dim=c(length(merged_filtered_gses)))
sigGene_results_human_voom=apply(as.matrix(1:ncol(pvalue_results_human_voom)),1,function(x) sum(p.adjust(pvalue_results_human_voom[,x],method="BH") < 0.05,na.rm=TRUE))
na_gses=apply(as.matrix(1:ncol(pvalue_results_human_voom)),1,function(x) length(is.na(pvalue_results_human_voom[,x])==FALSE)==0)
##################################################################################################'run voom DEG analysis for all GSE' ends here' 
##################################################################################################'1000 permutation

above_threshold=colnames(sigGene_results)[which(sigGene_results<7000)]
383
group1_counts=apply(as.matrix(1:length(merged_filtered_gses)),1,function(i) length(merged_filtered_gses[[i]]$group1))
hist(group1_counts,breaks = 50)  
head(sort(group1_counts,decreasing = TRUE))
which(group1_counts==39)
sigGene_results[which(group1_counts==39)]
128#just 4 sigificant genes, ignore
497
merged_filtered_gses[[497]]


which(group1_counts==29)
353
sigGene_results[which(group1_counts==29)]
601
merged_filtered_gses[[353]]


#####################################################################################PCA plot
#################data 1
i=353

geo_data<-getGEO("GSE107868",getGPL = FALSE)
show(geo_data)

pheno=pData(geo_data[[1]])
dim(pheno)
colnames(pheno)
gsm_title=cbind(rownames(pheno),as.character(pheno$title))
head(gsm_title)

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
gsm_title1=cbind(gsm_title,Group)
gsm_title1=gsm_title1[selected_all,]
colnames(gsm_title1)=c("GSM","Pheno","Group")

gsm_title1=as.data.frame(gsm_title1)
gsm_title1$Group <- as.factor(gsm_title1$Group)

TempIndices <- match(gsm_title1[,1], colnames(FilterGeneCountData))

colN=colnames(FilterGeneCountData)

FilterGeneCountData=FilterGeneCountData[,TempIndices]

y <- DGEList(counts=FilterGeneCountData, group=NULL)
y <- calcNormFactors(y)

cpm <- cpm(y, log = TRUE, prior.count = 0.01)
rv <- apply(cpm,1,var) 

#Select genes with highest variance.
keep <- order(rv, decreasing = TRUE)[1:500]
selected <- cpm[keep, ] %>% t()

#Transpose is needed to ensure that each row is a vector.
pca <- prcomp(selected, scale=T, center = T)

stddev <- pca$sdev
pc1_var <- round(100*stddev[1]^2/sum(stddev^2))
pc2_var <- round(100*stddev[2]^2/sum(stddev^2))
pc3_var <- round(100*stddev[3]^2/sum(stddev^2))
PlotData <- data.frame(cbind(PC1 = pca$x[,1], PC2 = pca$x[,2],PC3 = pca$x[,3]))
PlotData <- gsm_title1[, c("Group")] %>%
  cbind(PlotData, .)
colnames(PlotData)[4]="Group"
gp1=ggplot(PlotData, aes(x=PC1, y=PC2, color=Group, shape=Group)) + 
  geom_point(size=4.5) +  
  xlab(paste("PC1:", pc1_var, "% variance")) + 
  ylab(paste("PC2:", pc2_var, "% variance"))

gp1
jpeg("PCA_GSE107868_permSet1.jpeg")
print(gp1)
dev.off()

#################data 2
i=497

geo_data<-getGEO("GSE139061",getGPL = FALSE)
show(geo_data)

pheno=pData(geo_data[[1]])
dim(pheno)
colnames(pheno)
gsm_title=cbind(rownames(pheno),as.character(pheno$title))
head(gsm_title)

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
  gsm_title1=cbind(gsm_title,Group)
  gsm_title1=gsm_title1[selected_all,]
  colnames(gsm_title1)=c("GSM","Pheno","Group")
  
  gsm_title1=as.data.frame(gsm_title1)
  gsm_title1$Group <- as.factor(gsm_title1$Group)
  
  TempIndices <- match(gsm_title1[,1], colnames(FilterGeneCountData))
  
  colN=colnames(FilterGeneCountData)
  
  FilterGeneCountData=FilterGeneCountData[,TempIndices]
  
y <- DGEList(counts=FilterGeneCountData, group=NULL)
y <- calcNormFactors(y)

cpm <- cpm(y, log = TRUE, prior.count = 0.01)
rv <- apply(cpm,1,var) 

#Select genes with highest variance.
keep <- order(rv, decreasing = TRUE)[1:500]
selected <- cpm[keep, ] %>% t()

#Transpose is needed to ensure that each row is a vector.
pca <- prcomp(selected, scale=T, center = T)

stddev <- pca$sdev
pc1_var <- round(100*stddev[1]^2/sum(stddev^2))
pc2_var <- round(100*stddev[2]^2/sum(stddev^2))
pc3_var <- round(100*stddev[3]^2/sum(stddev^2))
PlotData <- data.frame(cbind(PC1 = pca$x[,1], PC2 = pca$x[,2],PC3 = pca$x[,3]))
PlotData <- gsm_title1[, c("Group")] %>%
  cbind(PlotData, .)
colnames(PlotData)[4]="Group"
rm(Group)
gp1=ggplot(PlotData, aes(x=PC1, y=PC2, color=Group, shape=Group)) + 
  geom_point(size=4.5) +  
  xlab(paste("PC1:", pc1_var, "% variance")) + 
  ylab(paste("PC2:", pc2_var, "% variance"))

gp1

jpeg("PCA_GSE139061_permSet2.jpeg")
print(gp1)
dev.off()
##################################################################################################'data 1 only edgeR


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
logFC_results_perm1=array(10,dim=c(length(genes),1000))
rownames(logFC_results_perm1)=genes
colnames(logFC_results_perm1)=paste0("set",1:1000)

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

logFC_results_perm1[which(pvalue_results_perm1==10)]=NA
pvalue_results_perm1[which(pvalue_results_perm1==10)]=NA

sum(pvalue_results_perm1==10,na.rm=TRUE)
0
sum(logFC_results_perm1==10,na.rm=TRUE)
0

end_t=Sys.time()#6 hours
save.image("../PFOCRInPathwayAnalyses_RData/two_permutation_tests_2ndset.RData")

#https://www.r-bloggers.com/2020/09/exact-tests-and-plots-with-aedger-basic-differential-expression-analysis/
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
##################################################################################################'data 1 only voom and exact


i=353
colnames(pvalue_results)[i]
#"GSE107868"

pvalue_results_perm1_voom=array(NA,dim=c(length(genes),1000))
dim(pvalue_results_perm1_voom)
rownames(pvalue_results_perm1_voom)=genes
colnames(pvalue_results_perm1_voom)=paste0("set",1:1000)
sigGene_results_perm1_voom=array(NA,dim=c(length(merged_filtered_gses)))
logFC_results_perm1_voom=array(NA,dim=c(length(genes),1000))
rownames(logFC_results_perm1_voom)=genes
colnames(logFC_results_perm1_voom)=paste0("set",1:1000)

pvalue_results_perm1_exact=array(NA,dim=c(length(genes),1000))
dim(pvalue_results_perm1_exact)
rownames(pvalue_results_perm1_exact)=genes
colnames(pvalue_results_perm1_exact)=paste0("set",1:1000)
sigGene_results_perm1_exact=array(NA,dim=c(length(merged_filtered_gses)))
logFC_results_perm1_exact=array(NA,dim=c(length(genes),1000))
rownames(logFC_results_perm1_exact)=genes
colnames(logFC_results_perm1_exact)=paste0("set",1:1000)

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
  
  FilterGeneCountData=FilterGeneCountData[,TempIndices]
  cbind(gsm_title1[,1],colnames(FilterGeneCountData))
  
  Group_p <- gsm_title1$Group
  print(levels(Group_p))
  
  ###Create DGEList object
  d0 <- DGEList(FilterGeneCountData,group=as.character(gsm_title1$Group))
  
  ###Calculate normalization factors
  d0 <- calcNormFactors(d0)
  d0 <- estimateDisp(d0)
  
  #Perform an exact test for treat vs ctrl
  exact_result <- exactTest(d0)
  dim(exact_result)
  dim(FilterGeneCountData)
  exact_result=exact_result$table
  
  pvalue_results_perm1_exact[match(rownames(exact_result),rownames(pvalue_results_perm1)),k]=exact_result$PValue
  sigGene_results_perm1_exact[k]=sum(p.adjust(exact_result$PValue,method="BH") < 0.05)
  logFC_results_perm1_exact[match(rownames(exact_result),rownames(pvalue_results_perm1)),k]=exact_result$logFC
  
  
  design <- model.matrix(~0+Group_p)
  design
  
  y <- voom(d0, design, plot = F)
  
  fit <- lmFit(y, design)
  head(coef(fit))
  
  contr <- makeContrasts(Group_pgroup2-Group_pgroup1, levels = colnames(coef(fit)))
  contr
  
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  voom_result=topTable(tmp, sort.by = "P", n = Inf)
  
  pvalue_results_perm1_voom[match(rownames(voom_result),rownames(pvalue_results_perm1)),k]=voom_result$P.Value
  sigGene_results_perm1_voom[k]=sum(voom_result$adj.P.Val < 0.05)
  logFC_results_perm1_voom[match(rownames(voom_result),rownames(pvalue_results_perm1)),k]=voom_result$logFC
  
  print(k)
  
}

#logFC_results_perm1_exact[which(logFC_results_perm1_exact==10)]=NA
#pvalue_results_perm1_exact[which(pvalue_results_perm1_exact==10)]=NA

#logFC_results_perm1_voom[which(logFC_results_perm1_voom==10)]=NA
#pvalue_results_perm1_voom[which(pvalue_results_perm1_voom==10)]=NA

sum(sigGene_results_perm1_voom[1:20]>0)
0
sum(sigGene_results_perm1_exact[1:20]>0)
17

end_t=Sys.time()#6 hours
save.image("../PFOCRInPathwayAnalyses_RData/two_permutation_tests_voom_exact.RData")

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
logFC_results_perm2=array(10,dim=c(length(genes),1000))
rownames(logFC_results_perm2)=genes
colnames(logFC_results_perm2)=paste0("set",1:1000)

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
  
  expression_perm2=expression[,sample(1:ncol(expression),ncol(expression),replace=FALSE)]
  colnames(expression_perm2)=colnames(expression)
  ###DEG
  
  filter <- apply(expression_perm2, 1, function(x) length(x[x>5])>=2) #minimum count (4) is satified in at least two conditions
  FilterGeneCountData <- expression_perm2[filter,]
  
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
logFC_results_perm2[which(is.na(pvalue_results_perm2)==TRUE)]=NA

end_t=Sys.time()#started at 5:55 pm
save.image("../PFOCRInPathwayAnalyses_RData/two_permutation_tests_2ndset.RData")

##################################################################################################'data 2 only voom and exact


i=497
perm_test_set2=as.matrix(pvalue_results[,i])
colnames(pvalue_results)[i]
#"GSE139061"

pvalue_results_perm2_voom=array(NA,dim=c(length(genes),1000))
dim(pvalue_results_perm2_voom)
rownames(pvalue_results_perm2_voom)=genes
colnames(pvalue_results_perm2_voom)=paste0("set",1:1000)
sigGene_results_perm2_voom=array(NA,dim=c(length(merged_filtered_gses)))
logFC_results_perm2_voom=array(NA,dim=c(length(genes),1000))
rownames(logFC_results_perm2_voom)=genes
colnames(logFC_results_perm2_voom)=paste0("set",1:1000)

pvalue_results_perm2_exact=array(NA,dim=c(length(genes),1000))
dim(pvalue_results_perm2_exact)
rownames(pvalue_results_perm2_exact)=genes
colnames(pvalue_results_perm2_exact)=paste0("set",1:1000)
sigGene_results_perm2_exact=array(NA,dim=c(length(merged_filtered_gses)))
logFC_results_perm2_exact=array(NA,dim=c(length(genes),1000))
rownames(logFC_results_perm2_exact)=genes
colnames(logFC_results_perm2_exact)=paste0("set",1:1000)

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
  
  FilterGeneCountData=FilterGeneCountData[,TempIndices]
  cbind(gsm_title1[,1],colnames(FilterGeneCountData))
  
  Group_p <- gsm_title1$Group
  print(levels(Group_p))
  
  ###Create DGEList object
  d0 <- DGEList(FilterGeneCountData,group=as.character(gsm_title1$Group))
  
  ###Calculate normalization factors
  d0 <- calcNormFactors(d0)
  d0 <- estimateDisp(d0)
  
  #Perform an exact test for treat vs ctrl
  exact_result <- exactTest(d0)
  dim(exact_result)
  dim(FilterGeneCountData)
  exact_result=exact_result$table
  
  pvalue_results_perm2_exact[match(rownames(exact_result),rownames(pvalue_results_perm2)),k]=exact_result$PValue
  sigGene_results_perm2_exact[k]=sum(p.adjust(exact_result$PValue,method="BH") < 0.05)
  logFC_results_perm2_exact[match(rownames(exact_result),rownames(pvalue_results_perm2)),k]=exact_result$logFC
  
  
  design <- model.matrix(~0+Group_p)
  design
  
  y <- voom(d0, design, plot = F)
  
  fit <- lmFit(y, design)
  head(coef(fit))
  
  contr <- makeContrasts(Group_pgroup2-Group_pgroup1, levels = colnames(coef(fit)))
  contr
  
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  voom_result=topTable(tmp, sort.by = "P", n = Inf)
  
  pvalue_results_perm2_voom[match(rownames(voom_result),rownames(pvalue_results_perm2)),k]=voom_result$P.Value
  sigGene_results_perm2_voom[k]=sum(voom_result$adj.P.Val < 0.05)
  logFC_results_perm2_voom[match(rownames(voom_result),rownames(pvalue_results_perm2)),k]=voom_result$logFC
  
  print(k)
  
}

#logFC_results_perm2_exact[which(logFC_results_perm2_exact==10)]=NA
#pvalue_results_perm2_exact[which(pvalue_results_perm2_exact==10)]=NA

#logFC_results_perm2_voom[which(logFC_results_perm2_voom==10)]=NA
#pvalue_results_perm2_voom[which(pvalue_results_perm2_voom==10)]=NA

sum(sigGene_results_perm2_voom[1:20]>0)
1
sum(sigGene_results_perm2_exact[1:20]>0)
20

end_t=Sys.time()#6 hours
save.image("../PFOCRInPathwayAnalyses_RData/two_permutation_tests_voom_exact.RData")

################################################################################################################### number of significant genes

how_many_sig1_voom=apply(as.matrix(1:1000),1,function(i) sum(p.adjust(pvalue_results_perm1_voom[,i],method="BH") < 0.05, na.rm=TRUE) )
how_many_sig1_exact=apply(as.matrix(1:1000),1,function(i) sum(p.adjust(pvalue_results_perm1_exact[,i],method="BH") < 0.05, na.rm=TRUE) )
how_many_sig1=apply(as.matrix(1:1000),1,function(i) sum(p.adjust(pvalue_results_perm1[,i],method="BH") < 0.05, na.rm=TRUE) )

how_many_sig2_voom=apply(as.matrix(1:1000),1,function(i) sum(p.adjust(pvalue_results_perm2_voom[,i],method="BH") < 0.05, na.rm=TRUE) )
how_many_sig2_exact=apply(as.matrix(1:1000),1,function(i) sum(p.adjust(pvalue_results_perm2_exact[,i],method="BH") < 0.05, na.rm=TRUE) )
how_many_sig2=apply(as.matrix(1:1000),1,function(i) sum(p.adjust(pvalue_results_perm2[,i],method="BH") < 0.05, na.rm=TRUE) )

jpeg("voom_permutation_significant_gene_hist_set1.jpeg")
hist(how_many_sig1_voom,breaks=50)
dev.off()

jpeg("exact_permutation_significant_gene_hist_set1.jpeg")
hist(how_many_sig1_exact,breaks=50)
dev.off()

jpeg("permutation_significant_gene_hist_set1.jpeg")
hist(how_many_sig1,breaks=50)
dev.off()

jpeg("voom_permutation_significant_gene_hist_set2.jpeg")
hist(how_many_sig2_voom,breaks=50)
dev.off()

jpeg("exact_permutation_significant_gene_hist_set2.jpeg")
hist(how_many_sig2_exact,breaks=50)
dev.off()

jpeg("permutation_significant_gene_hist_set2.jpeg")
hist(how_many_sig2,breaks=50)
dev.off()


################################################################################################################### find outliers and remove them
#########################################################################################set1

i=353
colnames(pvalue_results)[i]
#"GSE107868"

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
  
  outliers1=gof(fit)
  head(outliers1$outlier)
  sum(outliers1$outlier==TRUE) #641
  oustliers_list1=names(outliers1$outlier[which(outliers1$outlier==TRUE)])
  
  pvalue_results_perm1_noOutlier=pvalue_results_perm1
  sum(rownames(pvalue_results_perm1_noOutlier)!=names(outliers1$outlier))
  pvalue_results_perm1_noOutlier[oustliers_list1,]=NA
  
  logFC_results_perm1_noOutlier=logFC_results_perm1
  logFC_results_perm1_noOutlier[oustliers_list1,]=NA

  
  #########################################################################################set2
  
  i=497
  perm_test_set2=as.matrix(pvalue_results[,i])
  #"GSE139061"
  
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
  
  outliers2=gof(fit)
  head(outliers2$outlier)
  sum(outliers2$outlier==TRUE)#1213
  oustliers_list2=names(outliers2$outlier[which(outliers2$outlier==TRUE)])
  
  pvalue_results_perm2_noOutlier=pvalue_results_perm2
  pvalue_results_perm2_noOutlier[oustliers_list2,]=NA
  
  logFC_results_perm2_noOutlier=logFC_results_perm2
  logFC_results_perm2_noOutlier[oustliers_list2,]=NA
  
###################################################################################################################call database , Entrez ID match
require(rSEA)
require("clusterProfiler")
require("org.Hs.eg.db")
require(dplyr)

wp <- read.gmt("wikipathways-20200210-gmt-Homo_sapiens.gmt")
wp <- wp %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wp=wp[,c("name","wpid","gene")]
wp=unique(wp)

wp_list=split(wp$gene,wp$wpid)

wp_annotation=wp
colnames(wp_annotation)[2]="set_id"



pfocr <- read.gmt("/Users/mingyoungshin/git/PFOCRInPathwayAnalyses/pfocr20200224.gmt")
pfocr <- pfocr %>% tidyr::separate(ont, c("name","file","pfocrid","org"), "%")
pfocr=pfocr[,c("name","file","gene")]
pfocr=unique(pfocr)

pfocr_list=split(pfocr$gene,pfocr$file)

pfocr_annotation=pfocr
colnames(pfocr_annotation)[2]="set_id"



require(GO.db)
require(org.Hs.eg.db)


go_to_entrez <- as.data.frame(org.Hs.egGO2EG)
goegbp=go_to_entrez[go_to_entrez$Ontology=="BP",]
goegbp=goegbp[,1:2]

go_list=split(goegbp$gene_id,goegbp$go_id)
go_dataframe <- data.frame(set_id = rep(names(go_list), sapply(go_list, length)),
                 gene = unlist(go_list))


go_annotation=  AnnotationDbi::select(GO.db, keys=as.character(unique(go_dataframe$set_id)), columns=c("TERM","ONTOLOGY"),
                           keytype="GOID")

colnames(go_dataframe)[1]="GOID"
go_annotation=merge(go_dataframe,go_annotation[,c(1,2)],by="GOID")
colnames(go_annotation)=c("set_id","gene","name")
go_annotation=go_annotation[,c(3,1,2)]

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
  
  #################################################################################################################### run rSEA function_ database as a parameter 
  run_rSEA2<-function(data,list_result,database,database_annotation){#limma voom + exact test + and PCA for the 2nd #wilcoxon but not for small sample size # remove patients # exact test # simulate RNAseq data
    
    na_row=which(is.na(data)==TRUE)
    if(length(na_row)>0){
      data=as.matrix(data[-na_row,1]) 
    }
    
    if(nrow(data)==0){
      list_result=modifyList(list_result, list(set_id = cbind(list_result$set_id,NA), Coverage = cbind(list_result$Coverage,NA) ,
                                               TDP.bound = cbind(list_result$TDP.bound,NA), TDP.estimate = cbind(list_result$TDP.estimate,NA),  
                                               SC.adjP = cbind(list_result$SC.adjP,NA), Comp.adjP = cbind(list_result$Comp.adjP,NA),
                                               ID = cbind(list_result$ID), Size = cbind(list_result$Size), name=list_result$name) )
      return(list_result)
    }
    
    data_m=as.data.frame(cbind(rownames(data),data))
    colnames(data_m)=c("Gene","pvalue")
    merged=merge(data_m,gene_entrez,by="Gene",all = "L",all.x=TRUE,all.y=FALSE)
    
    # head(merged)
    
    sea_result<-SEA(as.numeric(as.character(merged$pvalue)), merged$ENTREZID, pathlist = database)
    colnames(sea_result)[2]="set_id"
    #head(sea_result)
    set_result=merge(sea_result,unique(database_annotation[,-3]),by="set_id")
    # head(set_result)
    
    
    
    list_result=modifyList(list_result, list(set_id = cbind(list_result$set_id,set_result$set_id), Coverage = cbind(list_result$Coverage,set_result$Coverage) ,
                                             TDP.bound = cbind(list_result$TDP.bound,set_result$TDP.bound), TDP.estimate = cbind(list_result$TDP.estimate,set_result$TDP.estimate),  
                                             SC.adjP = cbind(list_result$SC.adjP,set_result$SC.adjP), Comp.adjP = cbind(list_result$Comp.adjP,set_result$Comp.adjP),
                                             ID = cbind(list_result$ID), Size = cbind(list_result$Size), name=list_result$name) )
    
    return(list_result)
    
    
  }
  
  
####################################################################################################################run rSEA function  
run_rSEA<-function(data,list_result){#limma voom + exact test + and PCA for the 2nd #wilcoxon but not for small sample size # remove patients # exact test # simulate RNAseq data
  
  na_row=which(is.na(data)==TRUE)
  if(length(na_row)>0){
    data=as.matrix(data[-na_row,1]) 
  }
  
  if(nrow(data)==0){
    list_result=modifyList(list_result, list(wpid = cbind(list_result$wpid,NA), Coverage = cbind(list_result$Coverage,NA) ,
                                             TDP.bound = cbind(list_result$TDP.bound,NA), TDP.estimate = cbind(list_result$TDP.estimate,NA),  
                                             SC.adjP = cbind(list_result$SC.adjP,NA), Comp.adjP = cbind(list_result$Comp.adjP,NA),
                                             ID = cbind(list_result$ID), Size = cbind(list_result$Size), name=list_result$name) )
    return(list_result)
  }
  
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

#######NA issue

which(is.na(perm_test_set1_list[[1]]$SC.adjP)==TRUE)
541
perm_test_set1_list[[1]]$wpid[541]
# "WP661"
#wp[which(wp$wpid=="WP661"),]
#name  wpid gene
#14032 Glucose Homeostasis WP661 3630

#gene_entrez[13152,]
#Gene ENTREZID
#13153  INS     3630
pvalue_results_perm2[13152,]

#check whether wpids are reported in the same order
for(i in 1:length(perm_test_set2_list[[i]]$wpid)){
  wpids=NULL
  for(j in 1:1000){
    wpids=c(wpids,perm_test_set2_list[[j]]$wpid[i])
  }
  if(length(table(wpids))!=1){
    print(paste0(i," ",j))
  }
}



perm_test_set1_list <- vector("list", 9)
names(perm_test_set1_list) <- c("wpid", "ID", "Size", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP", "name")

perm_test_set2_list <- vector("list", 9)
names(perm_test_set2_list) <- c("wpid", "ID", "Size", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP", "name")

perm_test_set1_list_voom <- vector("list", 9)
names(perm_test_set1_list_voom) <- c("wpid", "ID", "Size", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP", "name")

perm_test_set2_list_voom <- vector("list", 9)
names(perm_test_set2_list_voom) <- c("wpid", "ID", "Size", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP", "name")

perm_test_set1_list_exact <- vector("list", 9)
names(perm_test_set1_list_exact) <- c("wpid", "ID", "Size", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP", "name")

perm_test_set2_list_exact <- vector("list", 9)
names(perm_test_set2_list_exact) <- c("wpid", "ID", "Size", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP", "name")

perm_test_set1_list_noOutlier <- vector("list", 9)
names(perm_test_set1_list_noOutlier) <- c("wpid", "ID", "Size", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP", "name")

perm_test_set2_list_noOutlier <- vector("list", 9)
names(perm_test_set2_list_noOutlier) <- c("wpid", "ID", "Size", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP", "name")

rsea_results_human_voom <- vector("list", 9)
names(rsea_results_human_voom) <- c("wpid", "ID", "Size", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP", "name")


rsea_results_human_voom_pfocr <- vector("list", 9)
names(rsea_results_human_voom_pfocr) <- c("set_id", "ID", "Size", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP", "name")

rsea_results_human_voom_go <- vector("list", 9)
names(rsea_results_human_voom_go) <- c("set_id", "ID", "Size", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP", "name")

rsea_results_human_voom_wp <- vector("list", 9)
names(rsea_results_human_voom_wp) <- c("set_id", "ID", "Size", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP", "name")


start_t=Sys.time()
rsea_results_human_voom_wp=apply(as.matrix(1:ncol(pvalue_results_human_voom)),1,function(x) run_rSEA2(as.matrix(pvalue_results_human_voom[,x]),rsea_results_human_voom_wp,wp_list,wp_annotation))
end_t=Sys.time() # 20 min

start_t=Sys.time()
rsea_results_human_voom_go=apply(as.matrix(1:ncol(pvalue_results_human_voom)),1,function(x) run_rSEA2(as.matrix(pvalue_results_human_voom[,x]),rsea_results_human_voom_go,go_list,go_annotation))
end_t=Sys.time() # 20 min

start_t=Sys.time()#started at 11:20 am 
rsea_results_human_voom_pfocr=apply(as.matrix(1:ncol(pvalue_results_human_voom)),1,function(x){ run_rSEA2(as.matrix(pvalue_results_human_voom[,x]),rsea_results_human_voom_pfocr,pfocr_list,pfocr_annotation);print(x)})
end_t=Sys.time() # 20 min


save.image("~/git/PFOCRInPathwayAnalyses_RData/human_DEG_result.RData")

########printout results
apply(as.matrix(1:length(merged_filtered_gses)),1,function(x){ dir.create(file.path(paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/", merged_filtered_gses[[x]]$GSE)), showWarnings = FALSE);
      dir.create(file.path(paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/", merged_filtered_gses[[x]]$GSE),"/rSEA"), showWarnings = FALSE);
  dir.create(file.path(paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/", merged_filtered_gses[[x]]$GSE),"/rSEA/WP"), showWarnings = FALSE);
  dir.create(file.path(paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/", merged_filtered_gses[[x]]$GSE),"/rSEA/PFOCR"), showWarnings = FALSE);
  dir.create(file.path(paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/", merged_filtered_gses[[x]]$GSE),"/rSEA/GO"), showWarnings = FALSE);
      dir.create(file.path(paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/", merged_filtered_gses[[x]]$GSE),"/ORA"), showWarnings = FALSE);
  dir.create(file.path(paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/", merged_filtered_gses[[x]]$GSE),"/ORA/WP"), showWarnings = FALSE);
  dir.create(file.path(paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/", merged_filtered_gses[[x]]$GSE),"/ORA/PFOCR"), showWarnings = FALSE);
  dir.create(file.path(paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/", merged_filtered_gses[[x]]$GSE),"/ORA/GO"), showWarnings = FALSE);
      dir.create(file.path(paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/", merged_filtered_gses[[x]]$GSE),"/GSEA"), showWarnings = FALSE);
  dir.create(file.path(paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/", merged_filtered_gses[[x]]$GSE),"/GSEA/WP"), showWarnings = FALSE);
  dir.create(file.path(paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/", merged_filtered_gses[[x]]$GSE),"/GSEA/GO"), showWarnings = FALSE);
  dir.create(file.path(paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/", merged_filtered_gses[[x]]$GSE),"/GSEA/GSEA"), showWarnings = FALSE)})

apply(as.matrix(1:length(merged_filtered_gses)),1,function(x){ 
  dir.create(file.path(paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/", merged_filtered_gses[[x]]$GSE),"/GSEA/PFOCR"), showWarnings = FALSE)})

apply(as.matrix(1:length(merged_filtered_gses)),1, function(x) write.table(as.data.frame(rsea_results_human_voom_wp[[x]]),
                                                      file=paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/",  merged_filtered_gses[[x]]$GSE, "/rSEA/WP/result.txt",
                                                      colnames(merged_filtered_gses)[x]),col.names=names(rsea_results_human_voom_wp[[x]]),row.names=FALSE,quote=FALSE))
apply(as.matrix(1:length(merged_filtered_gses)),1, function(x) write.table(as.data.frame(rsea_results_human_voom_go[[x]]),
                                                     file=paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/",  merged_filtered_gses[[x]]$GSE, "/rSEA/GO/result.txt",
                                                         colnames(merged_filtered_gses)[x]),col.names=names(rsea_results_human_voom_go[[x]]),row.names=FALSE,quote=FALSE))
apply(as.matrix(1:length(merged_filtered_gses)),1, function(x) write.table(as.data.frame(rsea_results_human_voom_pfocr[[x]]),
                                                     file=paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/",  merged_filtered_gses[[x]]$GSE, "/rSEA/PFOCR/result.txt",
                                                  colnames(merged_filtered_gses)[x]),col.names=names(rsea_results_human_voom_pfocr[[x]]),row.names=FALSE,quote=FALSE))

start_t=Sys.time()
rsea_results_human_voom=apply(as.matrix(1:ncol(pvalue_results_human_voom)),1,function(x) run_rSEA(as.matrix(pvalue_results_human_voom[,x]),rsea_results_human_voom))
end_t=Sys.time() # 20 mins

start_t=Sys.time()
perm_test_set2_list=apply(as.matrix(1:ncol(pvalue_results_perm2)),1,function(x) run_rSEA(as.matrix(pvalue_results_perm2[,x]),perm_test_set2_list))
end_t=Sys.time()

start_t=Sys.time()
perm_test_set1_list=apply(as.matrix(1:ncol(pvalue_results_perm1)),1,function(x) run_rSEA(as.matrix(pvalue_results_perm1[,x]),perm_test_set1_list))
end_t=Sys.time()

save.image("../PFOCRInPathwayAnalyses_RData/two_permutation_tests_2ndset.RData")

start_t=Sys.time()
perm_test_set2_list_voom=apply(as.matrix(1:ncol(pvalue_results_perm2_voom)),1,function(x) run_rSEA(as.matrix(pvalue_results_perm2_voom[,x]),perm_test_set2_list_voom))
end_t=Sys.time()

start_t=Sys.time()
perm_test_set1_list_voom=apply(as.matrix(1:ncol(pvalue_results_perm1_voom)),1,function(x) run_rSEA(as.matrix(pvalue_results_perm1_voom[,x]),perm_test_set1_list_voom))
end_t=Sys.time()

start_t=Sys.time()
perm_test_set2_list_exact=apply(as.matrix(1:ncol(pvalue_results_perm2_exact)),1,function(x) run_rSEA(as.matrix(pvalue_results_perm2_exact[,x]),perm_test_set2_list_exact))
end_t=Sys.time()

start_t=Sys.time()
perm_test_set1_list_exact=apply(as.matrix(1:ncol(pvalue_results_perm1_exact)),1,function(x) run_rSEA(as.matrix(pvalue_results_perm1_exact[,x]),perm_test_set1_list_exact))
end_t=Sys.time()#3:20 hours

save.image("../PFOCRInPathwayAnalyses_RData/two_permutation_tests_voom_exact.RData")

start_t=Sys.time()
perm_test_set2_list_noOutlier=apply(as.matrix(1:ncol(pvalue_results_perm2_noOutlier)),1,function(x) run_rSEA(as.matrix(pvalue_results_perm2_noOutlier[,x]),perm_test_set2_list_noOutlier))
end_t=Sys.time()

start_t=Sys.time()
perm_test_set1_list_noOutlier=apply(as.matrix(1:ncol(pvalue_results_perm1_noOutlier)),1,function(x) run_rSEA(as.matrix(pvalue_results_perm1_noOutlier[,x]),perm_test_set1_list_noOutlier))
end_t=Sys.time()

save.image("../PFOCRInPathwayAnalyses_RData/two_permutation_tests_voom_exact_noOutlier.RData")
# 
# SC_sig=array(0,dim=c(nrow(perm_test_set1_list[[1]]$wpid),1))
# for(i in 1:1000){
#   
#   sigs=which(perm_test_set1_list[[i]]$SC.adjP< 0.05)
#   if(length(sigs)>0){
#     SC_sig[sigs]=SC_sig[sigs]+1
#   }
# }  
# 
# 
# Comp_sig=array(0,dim=c(nrow(perm_test_set1_list[[1]]$wpid),1))
# for(i in 1:1000){
#   
#   sigs=which(perm_test_set1_list[[i]]$Comp.adjP< 0.05)
#   if(length(sigs)>0){
#     Comp_sig[sigs]=Comp_sig[sigs]+1
#   }
# }  
# 
# length(which(Comp_sig >= 50))
# #View(as.matrix(Comp_sig))
# jpeg("Comp.adjP_permutation_GSE107868.jpeg")
# print(hist(SC_sig,breaks=100))
# dev.off()


Comp_sig_human=array(0,dim=length(rsea_results_human_voom))
for(i in 1:length(rsea_results_human_voom)){
  
  Comp_sig_human[i]=sum(rsea_results_human_voom[[i]]$Comp.adjP< 0.05,na.rm=TRUE)
  
}  
#plot.ecdf(Comp_sig_human)
sum(Comp_sig_human>0)
206

SC_sig_human=array(0,dim=length(rsea_results_human_voom))
for(i in 1:length(rsea_results_human_voom)){
  
  SC_sig_human[i]=sum(rsea_results_human_voom[[i]]$SC.adjP< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(SC_sig_human)
hist(SC_sig_human,breaks=100)
sum(SC_sig_human>0)
206



Comp_sig_perm=array(0,dim=c(1000))
for(i in 1:1000){
  
  Comp_sig_perm[i]=sum(perm_test_set1_list[[i]]$Comp.adjP< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(Comp_sig_perm)
sum(Comp_sig_perm>0)
138

SC_sig_perm=array(0,dim=c(1000))
for(i in 1:1000){
  
  SC_sig_perm[i]=sum(perm_test_set1_list[[i]]$SC.adjP< 0.05,na.rm=TRUE)

}  
plot.ecdf(SC_sig_perm)
sum(SC_sig_perm>0)
138
#View(as.matrix(SC_sig))

Comp_sig_perm_voom=array(0,dim=c(1000))
for(i in 1:1000){
  
  Comp_sig_perm_voom[i]=sum(perm_test_set1_list_voom[[i]]$Comp.adjP< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(Comp_sig_perm_voom)
sum(Comp_sig_perm_voom>0)
11

SC_sig_perm_voom=array(0,dim=c(1000))
for(i in 1:1000){
  
  SC_sig_perm_voom[i]=sum(perm_test_set1_list_voom[[i]]$SC.adjP< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(SC_sig_perm_voom)
sum(SC_sig_perm_voom>0)
12

Comp_sig_perm_exact=array(0,dim=c(1000))
for(i in 1:1000){
  
  Comp_sig_perm_exact[i]=sum(perm_test_set1_list_exact[[i]]$Comp.adjP< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(Comp_sig_perm_exact)
sum(Comp_sig_perm_exact>0)
171

SC_sig_perm_exact=array(0,dim=c(1000))
for(i in 1:1000){
  
  SC_sig_perm_exact[i]=sum(perm_test_set1_list_exact[[i]]$SC.adjP< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(SC_sig_perm_exact)
sum(SC_sig_perm_exact>0)
171

Comp_sig_perm_noOutlier=array(0,dim=c(1000))
for(i in 1:1000){
  
  Comp_sig_perm_noOutlier[i]=sum(perm_test_set1_list_noOutlier[[i]]$Comp.adjP< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(Comp_sig_perm_noOutlier)
sum(Comp_sig_perm_noOutlier>0)
121

SC_sig_perm_noOutlier=array(0,dim=c(1000))
for(i in 1:1000){
  
  SC_sig_perm_noOutlier[i]=sum(perm_test_set1_list_noOutlier[[i]]$SC.adjP< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(SC_sig_perm_noOutlier)
sum(SC_sig_perm_noOutlier>0)
121




Comp_sig_perm2=array(0,dim=c(1000))
for(i in 1:1000){
  
  Comp_sig_perm2[i]=sum(perm_test_set2_list[[i]]$Comp.adjP< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(Comp_sig_perm2)
sum(Comp_sig_perm2>0)
487

SC_sig_perm2=array(0,dim=c(1000))
for(i in 1:1000){
  
  SC_sig_perm2[i]=sum(perm_test_set2_list[[i]]$SC.adjP< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(SC_sig_perm2)
sum(SC_sig_perm2>0)
489
# 
# SC_sig2=array(0,dim=c(nrow(perm_test_set2_list[[1]]$wpid),1))
# for(i in 1:1000){
#   
#   sigs=which(perm_test_set2_list[[i]]$SC.adjP< 0.05)
#   if(length(sigs)>0){
#     SC_sig2[sigs]=SC_sig2[sigs]+1
#   }
# }  
# 
# length(which(SC_sig2 >= 50))
# 77
#View(as.matrix(SC_sig2))
# jpeg("SC.adjP_permutation_GSE139061.jpeg")
# print(hist(SC_sig2,breaks=100))
# dev.off()
# 
# Comp_sig2=array(0,dim=c(nrow(perm_test_set2_list[[1]]$wpid),1))
# for(i in 1:1000){
#   
#   sigs=which(perm_test_set2_list[[i]]$Comp.adjP< 0.05)
#   if(length(sigs)>0){
#     Comp_sig2[sigs]=Comp_sig2[sigs]+1
#   }
# }  
# 
# length(which(Comp_sig2 >= 50))
# 76
# #View(as.matrix(Comp_sig2))
# jpeg("Comp.adjP_permutation_GSE139061.jpeg")
# print(hist(Comp_sig2,breaks=100))
# dev.off()



Comp_sig_perm_voom2=array(0,dim=c(1000))
for(i in 1:1000){
  
  Comp_sig_perm_voom2[i]=sum(perm_test_set2_list_voom[[i]]$Comp.adjP< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(Comp_sig_perm_voom2)
sum(Comp_sig_perm_voom2>0)
44

SC_sig_perm_voom2=array(0,dim=c(1000))
for(i in 1:1000){
  
  SC_sig_perm_voom2[i]=sum(perm_test_set2_list_voom[[i]]$SC.adjP< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(SC_sig_perm_voom2)
sum(SC_sig_perm_voom2>0)
46

Comp_sig_perm_exact2=array(0,dim=c(1000))
for(i in 1:1000){
  
  Comp_sig_perm_exact2[i]=sum(perm_test_set2_list_exact[[i]]$Comp.adjP< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(Comp_sig_perm_exact2)
sum(Comp_sig_perm_exact2>0)
734

SC_sig_perm_exact2=array(0,dim=c(1000))
for(i in 1:1000){
  
  SC_sig_perm_exact2[i]=sum(perm_test_set2_list_exact[[i]]$SC.adjP< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(SC_sig_perm_exact2)
sum(SC_sig_perm_exact2>0)
734

Comp_sig_perm_noOutlier2=array(0,dim=c(1000))
for(i in 1:1000){
  
  Comp_sig_perm_noOutlier2[i]=sum(perm_test_set2_list_noOutlier[[i]]$Comp.adjP< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(Comp_sig_perm_noOutlier2)
sum(Comp_sig_perm_noOutlier2>0)
290

SC_sig_perm_noOutlier2=array(0,dim=c(1000))
for(i in 1:1000){
  
  SC_sig_perm_noOutlier2[i]=sum(perm_test_set2_list_noOutlier[[i]]$SC.adjP< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(SC_sig_perm_noOutlier2)
sum(SC_sig_perm_noOutlier2>0)
291



######################################################pathway intersection percentage
require(bayesbio)

set1_sig_pathway_jaccard=array(NA,dim=c(1000))

for(i in 1:1000){
  if( sum(perm_test_set1_list_voom[[i]]$SC.adjP< 0.05,na.rm=TRUE)==0 |  sum(perm_test_set1_list[[i]]$SC.adjP< 0.0,na.rm=TRUE)==0  ){
    if( sum(perm_test_set1_list_voom[[i]]$SC.adjP< 0.05,na.rm=TRUE) +  sum(perm_test_set1_list[[i]]$SC.adjP< 0.05,na.rm=TRUE) >0 ){
      set1_sig_pathway_jaccard[i]=0
    }
  }
    set1_sig_pathway_jaccard[i]=jaccardSets(perm_test_set1_list_voom[[i]]$wpid[perm_test_set1_list_voom[[i]]$SC.adjP< 0.05] , perm_test_set1_list[[i]]$wpid[perm_test_set1_list[[i]]$SC.adjP< 0.05])
  
}
jpeg("permuted_set1_signiciant_rSEA_pathway_jaccard.jpeg")
hist(set1_sig_pathway_jaccard,breaks=100)
dev.off()

set2_sig_pathway_jaccard=array(NA,dim=c(1000))

for(i in 1:1000){
 if( sum(perm_test_set2_list_voom[[i]]$SC.adjP< 0.05,na.rm=TRUE)==0 |  sum(perm_test_set2_list[[i]]$SC.adjP< 0.05,na.rm=TRUE)==0  ){
    if( sum(perm_test_set2_list_voom[[i]]$SC.adjP< 0.05,na.rm=TRUE) +  sum(perm_test_set2_list[[i]]$SC.adjP< 0.05,na.rm=TRUE) >0 ){
      set2_sig_pathway_jaccard[i]=0
    }
  }
    set2_sig_pathway_jaccard[i]=jaccardSets(perm_test_set2_list_voom[[i]]$wpid[perm_test_set2_list_voom[[i]]$SC.adjP< 0.05] , perm_test_set2_list[[i]]$wpid[perm_test_set2_list[[i]]$SC.adjP< 0.05])
  
}


jpeg("permuted_set2_signiciant_rSEA_pathway_jaccard.jpeg")
hist(set2_sig_pathway_jaccard,breaks=100)
dev.off()
for(i in 1:100){
  print(i)
  print(c(perm_test_set2_list_voom[[i]]$wpid[perm_test_set2_list_voom[[i]]$SC.adjP< 0.05] ))
  print((perm_test_set2_list[[i]]$wpid[perm_test_set2_list[[i]]$SC.adjP< 0.05]))
   
}
save.image("../PFOCRInPathwayAnalyses_RData/Nov17.RData")

###################################################################################################################run ORA function
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
  if(length(na_row)>0){
    data=as.matrix(data[-na_row,1]) 
  }
  
  data_m=as.data.frame(cbind(rownames(data),data))
  colnames(data_m)=c("Gene","pvalue")
  merged=merge(data_m,gene_entrez,by="Gene",all = "L",all.x=TRUE,all.y=FALSE)
  
  enrichment_result <- clusterProfiler::enricher(
    merged$ENTREZID[which(p.adjust(as.numeric(as.character(merged$pvalue)),method="BH") < 0.05)],
    universe = merged$ENTREZID,
    pAdjustMethod = "holm",
    pvalueCutoff = 1, #p.adjust cutoff
    minGSSize = 1,
    maxGSSize = 100000,
    TERM2GENE = wp[,2:3],
    TERM2NAME = wp[,c(2,1)])
  
  #enrichment_result <- DOSE::setReadable(enrichment_result, org.Hs.eg.db, keyType = "ENTREZID")
  if(length(enrichment_result)>0){
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
  
  }
  return(list_result)
  

  
 
}


###################################################################################################################run ORA function, databse parameter
run_ORA2<-function(data, list_result, database, database_annotation){
  #data=as.matrix(pvalue_results_human_voom[,1]);list_result=ora_results_human_voom_wp;database=wp_list;database_annotation=wp_annotation
  
  GeneRatio=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  BgRatio=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  pvalue=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  p.adjust=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  qvalue=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  Count=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  rownames(GeneRatio)=unique(database_annotation$set_id)
  rownames(BgRatio)=unique(database_annotation$set_id)
  rownames(pvalue)=unique(database_annotation$set_id)
  rownames(p.adjust)=unique(database_annotation$set_id)
  rownames(qvalue)=unique(database_annotation$set_id)
  rownames(Count)=unique(database_annotation$set_id)
  
  na_row=which(is.na(data)==TRUE)
  if(length(na_row)>0){
    data=as.matrix(data[-na_row,1]) 
  }
  
  data_m=as.data.frame(cbind(rownames(data),data))
  colnames(data_m)=c("Gene","pvalue")
  merged=merge(data_m,gene_entrez,by="Gene",all = "L",all.x=TRUE,all.y=FALSE)
  
  enrichment_result <- clusterProfiler::enricher(
    merged$ENTREZID[which(p.adjust(as.numeric(as.character(merged$pvalue)),method="BH") < 0.05)],
    universe = merged$ENTREZID,
    pAdjustMethod = "holm",
    pvalueCutoff = 1, #p.adjust cutoff
    minGSSize = 1,
    maxGSSize = 100000,
    TERM2GENE = database_annotation[,2:3],
    TERM2NAME = database_annotation[,c(2,1)])
  
  #enrichment_result <- DOSE::setReadable(enrichment_result, org.Hs.eg.db, keyType = "ENTREZID")
  if(length(enrichment_result)>0){
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
    
  }
  return(list_result)
  
  
  
  
}


########################################################################number of significant gens
sig_c_perm1=apply(as.matrix(1:1000),1,function(x) { sum(p.adjust(pvalue_results_perm1[,x],method="BH") < 0.05,na.rm=TRUE) } )
sum(sig_c_perm1==0)
157

sig_c_perm2=apply(as.matrix(1:1000),1,function(x) { sum(p.adjust(pvalue_results_perm2[,x],method="BH") < 0.05,na.rm=TRUE) } )
sum(sig_c_perm2==0)
0

perm_test_set2_list_ora <- vector("list", 6)
names(perm_test_set2_list_ora) <- c("GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue","Count")

start_t=Sys.time()
perm_test_set2_list_ora=apply(as.matrix(1:ncol(pvalue_results_perm2)),1,function(x) run_ORA(as.matrix(pvalue_results_perm2[,x]),perm_test_set2_list_ora))
end_t=Sys.time()

#start_t=Sys.time()
#for(x in 1:1000){
#  print(x)
#  run_ORA(as.matrix(pvalue_results_perm2[,x]),perm_test_set2_list_ora)
#}
#end_t=Sys.time()

perm_test_set1_list_ora <- vector("list", 6)
names(perm_test_set1_list_ora) <- c("GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue","Count")

start_t=Sys.time()
perm_test_set1_list_ora=apply(as.matrix(1:ncol(pvalue_results_perm1)),1,function(x) run_ORA(as.matrix(pvalue_results_perm1[,x]),perm_test_set1_list_ora))
end_t=Sys.time()



perm_test_set1_list_ora_voom <- vector("list", 6)
names(perm_test_set1_list_ora_voom) <- c("GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue","Count")

start_t=Sys.time()
perm_test_set1_list_ora_voom=apply(as.matrix(1:ncol(pvalue_results_perm1_voom)),1,function(x) run_ORA(as.matrix(pvalue_results_perm1_voom[,x]),perm_test_set1_list_ora_voom))
end_t=Sys.time()

perm_test_set2_list_ora_voom <- vector("list", 6)
names(perm_test_set2_list_ora_voom) <- c("GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue","Count")

ora_results_human_voom_wp <- vector("list", 6)
names(ora_results_human_voom_wp) <- c("GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue","Count")


ora_results_human_voom_pfocr <- vector("list", 6)
names(ora_results_human_voom_pfocr) <- c("GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue","Count")


ora_results_human_voom_go <- vector("list", 6)
names(ora_results_human_voom_go) <- c("GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue","Count")

start_t=Sys.time()
ora_results_human_voom_go=apply(as.matrix(1:ncol(pvalue_results_human_voom)),1,function(x) run_ORA2(as.matrix(pvalue_results_human_voom[,x]),ora_results_human_voom_go,go_list,go_annotation))
end_t=Sys.time()

saveRDS(ora_results_human_voom_go, file = "../PFOCRInPathwayAnalyses_RData/ora_results_human_voom_go")

start_t=Sys.time()
ora_results_human_voom_pfocr=apply(as.matrix(1:ncol(pvalue_results_human_voom)),1,function(x) run_ORA2(as.matrix(pvalue_results_human_voom[,x]),ora_results_human_voom_pfocr,pfocr_list,pfocr_annotation))
end_t=Sys.time()

saveRDS(ora_results_human_voom_pfocr, file = "../PFOCRInPathwayAnalyses_RData/ora_results_human_voom_pfocr")

start_t=Sys.time()
ora_results_human_voom_wp=apply(as.matrix(1:ncol(pvalue_results_human_voom)),1,function(x) run_ORA2(as.matrix(pvalue_results_human_voom[,x]),ora_results_human_voom_wp,wp_list,wp_annotation))
end_t=Sys.time()

saveRDS(ora_results_human_voom_wp, file = "../PFOCRInPathwayAnalyses_RData/ora_results_human_voom_wp.rds")

#########################print out results

apply(as.matrix(1:length(merged_filtered_gses)),1, function(x) write.table(as.data.frame(rsea_results_human_voom_wp[[x]]),
                                                                           file=paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/",  merged_filtered_gses[[x]]$GSE, "/ORA/WP/result.txt",
                                                                                       colnames(merged_filtered_gses)[x]),col.names=names(rsea_results_human_voom_wp[[x]]),row.names=FALSE,quote=FALSE))
apply(as.matrix(1:length(merged_filtered_gses)),1, function(x) write.table(as.data.frame(rsea_results_human_voom_go[[x]]),
                                                                           file=paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/",  merged_filtered_gses[[x]]$GSE, "/ORA/GO/result.txt",
                                                                                       colnames(merged_filtered_gses)[x]),col.names=names(rsea_results_human_voom_go[[x]]),row.names=FALSE,quote=FALSE))
apply(as.matrix(1:length(merged_filtered_gses)),1, function(x) write.table(as.data.frame(rsea_results_human_voom_pfocr[[x]]),
                                                                           file=paste0("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/GSE/",  merged_filtered_gses[[x]]$GSE, "/ORA/PFOCR/result.txt",
                                                                                       colnames(merged_filtered_gses)[x]),col.names=names(rsea_results_human_voom_pfocr[[x]]),row.names=FALSE,quote=FALSE))


start_t=Sys.time()
perm_test_set2_list_ora_voom=apply(as.matrix(1:ncol(pvalue_results_perm2_voom)),1,function(x) run_ORA(as.matrix(pvalue_results_perm2_voom[,x]),perm_test_set2_list_ora_voom))
end_t=Sys.time()

save.image("../PFOCRInPathwayAnalyses_RData/two_permutation_tests_voom_exact.RData")

save.image("../PFOCRInPathwayAnalyses_RData/two_permutation_tests_holm.RData")

ora_sig_perm=array(0,dim=c(1000))
for(i in 1:1000){
  
  ora_sig_perm[i]=sum(perm_test_set1_list_ora[[i]]$p.adjust< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(ora_sig_perm)
sum(ora_sig_perm>0)
364
#holm
314

ora_sig_perm2=array(0,dim=c(1000))
for(i in 1:1000){
  
  ora_sig_perm2[i]=sum(perm_test_set2_list_ora[[i]]$p.adjust< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(ora_sig_perm2)
sum(ora_sig_perm2>0)
654
#holm
573

ora_sig_perm_voom=array(0,dim=c(1000))
for(i in 1:1000){
  
  ora_sig_perm_voom[i]=sum(perm_test_set1_list_ora_voom[[i]]$p.adjust< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(ora_sig_perm_voom)
sum(ora_sig_perm_voom>0)
16
#holm
14

ora_sig_perm_voom2=array(0,dim=c(1000))
for(i in 1:1000){
  
  SC_sig_perm_voom2[i]=sum(perm_test_set2_list_ora_voom[[i]]$p.adjust< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(SC_sig_perm_voom2)
sum(SC_sig_perm_voom2>0)
64
#holm
60
# 
# 
# SC_sig=array(0,dim=c(nrow(perm_test_set1_list_ora[[1]]$p.adjust),1))
# for(i in 1:1000){
#   
#   sigs=which(perm_test_set1_list_ora[[i]]$p.adjust< 0.05)
#   if(length(sigs)>0){
#     SC_sig[sigs]=SC_sig[sigs]+1
#   }
# }  
# 
# length(which(SC_sig >= 50))
# 93
# #View(as.matrix(SC_sig))
# jpeg("ora_permutation_GSE107868.jpeg")
# print(hist(SC_sig,breaks=100))
# dev.off()
# 
# 
# SC_sig=array(0,dim=c(nrow(perm_test_set2_list_ora[[1]]$wpid),1))
# for(i in 1:1000){
#   
#   sigs=which(perm_test_set2_list_ora[[i]]$SC.adjP< 0.05)
#   if(length(sigs)>0){
#     SC_sig[sigs]=SC_sig[sigs]+1
#   }
# }  
# 
# length(which(SC_sig >= 50))
# 
# #View(as.matrix(SC_sig))
# jpeg("ora_permutation_GSE139061.jpeg")
# print(hist(SC_sig,breaks=100))
# dev.off()

###################################################################################################################run GSEA function
run_gsea<-function(data,logFC_data, list_result){
  
   
  na_row=which(is.na(data)==TRUE)
  if(length(na_row)>0){
    data=as.matrix(data[-na_row,1]) 
    logFC_data=as.matrix(logFC_data[-na_row,1])
  }

  
  data_m=as.data.frame(cbind(rownames(data),data))
  colnames(data_m)=c("Gene","pvalue")
  merged=merge(data_m,gene_entrez,by="Gene",all = "L",all.x=TRUE,all.y=FALSE)
  

  ### geneList prep
  gene_list<-sign(logFC_data) * - log10(as.numeric(as.character(merged[,2])))
  names(gene_list)=merged$ENTREZID
  gene_list = sort(gene_list[unique(names(gene_list))], decreasing = TRUE)
    
    gsewp.p <- GSEA(
      gene_list,
      pAdjustMethod="holm",
    TERM2GENE = wp[,2:3],
    TERM2NAME = wp[,c(2,1)],
    nPerm        = 1000,
    minGSSize = 1,
    maxGSSize = 100000,
    pvalueCutoff = 1,
    verbose=FALSE)
    
#  gsewp.p <- DOSE::setReadable(gsewp.p, org.Hs.eg.db, keyType = "ENTREZID")
 # head(gsewp.p, 20)
  enrichment_result=as.data.frame(gsewp.p)
  enrichment_result=enrichment_result[,-match("core_enrichment",colnames(enrichment_result))]
  
  enrichmentScore=matrix(NA,nrow=length(unique(wp$wpid)),ncol=1)
  NES=matrix(NA,nrow=length(unique(wp$wpid)),ncol=1)
  pvalue=matrix(NA,nrow=length(unique(wp$wpid)),ncol=1)
  p.adjust=matrix(NA,nrow=length(unique(wp$wpid)),ncol=1)
  qvalues=matrix(NA,nrow=length(unique(wp$wpid)),ncol=1)
  rank=matrix(NA,nrow=length(unique(wp$wpid)),ncol=1)
  rownames(enrichmentScore)=unique(wp$wpid)
  rownames(NES)=unique(wp$wpid)
  rownames(pvalue)=unique(wp$wpid)
  rownames(p.adjust)=unique(wp$wpid)
  rownames(qvalues)=unique(wp$wpid)
  rownames(rank)=unique(wp$wpid)
  
  enrichmentScore[match(enrichment_result$ID,rownames(enrichmentScore))]=enrichment_result$enrichmentScore
  NES[match(enrichment_result$ID,rownames(NES))]=enrichment_result$NES
  pvalue[match(enrichment_result$ID,rownames(pvalue))]=enrichment_result$pvalue
  p.adjust[match(enrichment_result$ID,rownames(p.adjust))]=enrichment_result$p.adjust
  qvalues[match(enrichment_result$ID,rownames(qvalues))]=enrichment_result$qvalues
  rank[match(enrichment_result$ID,rownames(rank))]=enrichment_result$rank
  
  list_result=modifyList(list_result, list(enrichmentScore = cbind(list_result$enrichmentScore,enrichmentScore), NES = cbind(list_result$NES,NES) ,
                                           pvalue = cbind(list_result$pvalue,pvalue), p.adjust = cbind(list_result$p.adjust,p.adjust),  
                                           qvalues = cbind(list_result$qvalues,qvalues), rank = cbind(list_result$rank,rank)) )
  
  return(list_result)
  
  

  
}


###################################################################################################################run GSEA function, database paramter
run_gsea1<-function(data,logFC_data, list_result, database, database_annotation){
  
  
  na_row=which(is.na(data)==TRUE)
  if(length(na_row)>0){
    data=as.matrix(data[-na_row,1]) 
    logFC_data=as.matrix(logFC_data[-na_row,1])
  }
  
  
  data_m=as.data.frame(cbind(rownames(data),data))
  colnames(data_m)=c("Gene","pvalue")
  merged=merge(data_m,gene_entrez,by="Gene",all = "L",all.x=TRUE,all.y=FALSE)
  
  
  ### geneList prep
  gene_list<-sign(logFC_data) * - log10(as.numeric(as.character(merged[,2])))
  names(gene_list)=merged$ENTREZID
  gene_list = sort(gene_list[unique(names(gene_list))], decreasing = TRUE)
  
  gsedatabase_annotation.p <- GSEA(
    gene_list,
    pAdjustMethod="holm",
    TERM2GENE = database_annotation[,2:3],
    TERM2NAME = database_annotation[,c(2,1)],
    nPerm        = 1000,
    minGSSize = 1,
    maxGSSize = 100000,
    pvalueCutoff = 1,
    verbose=FALSE)
  
  #  gsedatabase_annotation.p <- DOSE::setReadable(gsedatabase_annotation.p, org.Hs.eg.db, keyType = "ENTREZID")
  # head(gsedatabase_annotation.p, 20)
  enrichment_result=as.data.frame(gsedatabase_annotation.p)
  enrichment_result=enrichment_result[,-match("core_enrichment",colnames(enrichment_result))]
  
  enrichmentScore=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  NES=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  pvalue=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  p.adjust=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  qvalues=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  rank=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  rownames(enrichmentScore)=unique(database_annotation$set_id)
  rownames(NES)=unique(database_annotation$set_id)
  rownames(pvalue)=unique(database_annotation$set_id)
  rownames(p.adjust)=unique(database_annotation$set_id)
  rownames(qvalues)=unique(database_annotation$set_id)
  rownames(rank)=unique(database_annotation$set_id)
  
  enrichmentScore[match(enrichment_result$ID,rownames(enrichmentScore))]=enrichment_result$enrichmentScore
  NES[match(enrichment_result$ID,rownames(NES))]=enrichment_result$NES
  pvalue[match(enrichment_result$ID,rownames(pvalue))]=enrichment_result$pvalue
  p.adjust[match(enrichment_result$ID,rownames(p.adjust))]=enrichment_result$p.adjust
  qvalues[match(enrichment_result$ID,rownames(qvalues))]=enrichment_result$qvalues
  rank[match(enrichment_result$ID,rownames(rank))]=enrichment_result$rank
  
  list_result=modifyList(list_result, list(enrichmentScore = cbind(list_result$enrichmentScore,enrichmentScore), NES = cbind(list_result$NES,NES) ,
                                           pvalue = cbind(list_result$pvalue,pvalue), p.adjust = cbind(list_result$p.adjust,p.adjust),  
                                           qvalues = cbind(list_result$qvalues,qvalues), rank = cbind(list_result$rank,rank)) )
  
  return(list_result)
  
  
  
  
}



perm_test_set1_list_gsea <- vector("list", 6)
names(perm_test_set1_list_gsea) <- c("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalues","rank")

start_t=Sys.time()
perm_test_set1_list_gsea=apply(as.matrix(1:ncol(pvalue_results_perm1)),1,function(x) run_gsea(as.matrix(pvalue_results_perm1[,x]),as.matrix(logFC_results_perm1[,x]),perm_test_set1_list_gsea))
end_t=Sys.time()


perm_test_set2_list_gsea <- vector("list", 6)
names(perm_test_set2_list_gsea) <- c("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalues","rank")

start_t=Sys.time()
perm_test_set2_list_gsea=apply(as.matrix(1:ncol(pvalue_results_perm2)),1,function(x) run_gsea(as.matrix(pvalue_results_perm2[,x]),as.matrix(logFC_results_perm2[,x]),perm_test_set2_list_gsea))
end_t=Sys.time()

save.image("../PFOCRInPathwayAnalyses_RData/two_permutation_tests_2ndset_uptoGSEA.RData")



perm_test_set1_list_gsea_voom <- vector("list", 6)
names(perm_test_set1_list_gsea_voom) <- c("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalues","rank")

start_t=Sys.time()
perm_test_set1_list_gsea_voom=apply(as.matrix(1:ncol(pvalue_results_perm1_voom)),1,function(x) run_gsea(as.matrix(pvalue_results_perm1_voom[,x]),as.matrix(logFC_results_perm1_voom[,x]),perm_test_set1_list_gsea_voom))
end_t=Sys.time()


perm_test_set2_list_gsea_voom <- vector("list", 6)
names(perm_test_set2_list_gsea_voom) <- c("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalues","rank")

start_t=Sys.time()
perm_test_set2_list_gsea_voom=apply(as.matrix(1:ncol(pvalue_results_perm2_voom)),1,function(x) run_gsea(as.matrix(pvalue_results_perm2_voom[,x]),as.matrix(logFC_results_perm2_voom[,x]),perm_test_set2_list_gsea_voom))
end_t=Sys.time()

save.image("../PFOCRInPathwayAnalyses_RData/two_permutation_tests_2ndset_uptoGSEA_voom.RData")

save.image("../PFOCRInPathwayAnalyses_RData/two_permutation_tests_holm.RData")

gsea_sig_perm_voom1=array(0,dim=c(1000))
for(i in 1:1000){
  
  gsea_sig_perm_voom1[i]=sum(perm_test_set1_list_gsea_voom[[i]]$p.adjust< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(gsea_sig_perm_voom1)
sum(gsea_sig_perm_voom1>0)
0
#holm
0

gsea_sig_perm_voom2=array(0,dim=c(1000))
for(i in 1:1000){
  
  gsea_sig_perm_voom2[i]=sum(perm_test_set2_list_gsea_voom[[i]]$p.adjust< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(gsea_sig_perm_voom2)
sum(gsea_sig_perm_voom2>0)
57
#holm
57





gsea_sig_perm=array(0,dim=c(1000))
for(i in 1:1000){
  
  gsea_sig_perm[i]=sum(perm_test_set1_list_gsea[[i]]$p.adjust< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(gsea_sig_perm)
sum(gsea_sig_perm>0)
0
#holm
0

gsea_sig_perm2=array(0,dim=c(1000))
for(i in 1:1000){
  
  gsea_sig_perm2[i]=sum(perm_test_set2_list_gsea[[i]]$p.adjust< 0.05,na.rm=TRUE)
  
}  
plot.ecdf(gsea_sig_perm2)
sum(gsea_sig_perm2>0)
6
#holm
0


##############################################################################knockout keyword match

require("clusterProfiler")
require("org.Hs.eg.db")
require(dplyr)
require(readxl)




data<-read_excel("/Users/mingyoungshin/git/PFOCRInPathwayAnalyses/human_knockout_genelist.xlsx", col_names = F)
data_gene=as.character(unlist(data[,2]))


wp <- read.gmt("/Users/mingyoungshin/git/PFOCRInPathwayAnalyses/wikipathways-20200210-gmt-Homo_sapiens.gmt")
wp <- wp %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wp=wp[,c("name","wpid","gene")]
wp=unique(wp)

wp_list=split(wp$gene,wp$wpid)


gene_entrez <- bitr(data_gene, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db, drop=FALSE)
#17.65% of input gene IDs are fail to map
head(gene_entrez)
# 
#1	A20	NA
# 8	CBAP	NA
# 9	CD73	NA
# 22	E6AP	NA
# 29	Exoc5	NA
# 32	GRP78	NA
# 34	HBO1	NA
# 37	HSP90	NA
# 42	KT80	NA
# 45	miR-424	NA
# 47	MTR4	NA
# 53	PARP	NA
# 67	Rorc	NA
# 81	UBC9	NA
# 85	ZRF1	NA

colnames(gene_entrez)[1]="Gene"
gene_entrez=unique(gene_entrez)


tb=table(gene_entrez$Gene)
tb[which(tb>1)]


result=matrix(NA,nrow=nrow(gene_entrez),ncol=3)
for(i in 1:nrow(result)){
  result[i,1]=gene_entrez$Gene[i]
  result[i,2]=gene_entrez$ENTREZID[i]
  temp=NULL
  for(j in 1:length(wp_list)){
    
    if( result[i,2]%in%wp_list[[j]]){
      wp_name=wp$name[which(wp$wpid==names(wp_list)[j])[1]]
      #temp=paste(temp,wp_name,collapse =';') 
      temp=names(wp_list)[j]
    }
    if(length(temp)>0){result[i,3]=paste(result[i,3],temp,sep=";")}
    
  }
}


nrow(result)-sum(is.na(result[,3]))
40
#
41

result=gsub("NA;","",result)

colnames(result)=c("Symbol","Entrez","WP")
colnames(data)=c("GSE","Symbol")
merged_wp=merge(data,result,by="Symbol")

write.table(merged_wp,"human_knockout_wp.txt",col.names = colnames(merged_wp),quote=FALSE)




result_pfocr=matrix(NA,nrow=nrow(gene_entrez),ncol=3)
for(i in 1:nrow(result_pfocr)){
  result_pfocr[i,1]=gene_entrez$Gene[i]
  result_pfocr[i,2]=gene_entrez$ENTREZID[i]
  
  for(j in 1:length(pfocr_list)){
    temp=NULL
    if( result_pfocr[i,2]%in%pfocr_list[[j]]){
      pfocr_name=pfocr$name[which(pfocr$pfocrid==names(pfocr_list)[j])[1]]
      #temp=paste(temp,pfocr_name,collapse =';') 
      temp=names(pfocr_list)[j]
    }
    if(length(temp)>0){result_pfocr[i,3]=paste(result_pfocr[i,3],temp,sep=";")}
    
  }
}

result_pfocr=gsub("NA;","",result_pfocr)

nrow(result_pfocr)-sum(is.na(result_pfocr[,3]))
53
#
55


colnames(result_pfocr)=c("Symbol","Entrez","PFOCR")
merged_pfocr=merge(data,result_pfocr,by="Symbol")
write.table(merged_pfocr,"human_knockout_pfocr.txt",col.names = colnames(merged_pfocr),quote=FALSE)

dim(data)
dim(merged_pfocr)

save.image("../PFOCRInPathwayAnalyses_RData/updating_pfocr_match.RData")




result_go=matrix(NA,nrow=nrow(gene_entrez),ncol=3)
for(i in 1:nrow(result_go)){
  result_go[i,1]=gene_entrez$Gene[i]
  result_go[i,2]=gene_entrez$ENTREZID[i]
  
  for(j in 1:length(go_list)){
    temp=NULL
    if( result_go[i,2]%in%go_list[[j]]){
      
      #temp=paste(temp,go_name,collapse =';') 
      temp=names(go_list)[j]
    }
    if(length(temp)>0){result_go[i,3]=paste(result_go[i,3],temp,sep=";")}
    
  }
}

result_go=gsub("NA;","",result_go)

nrow(result_go)-sum(is.na(result_go[,3]))
65

colnames(result_go)=c("Symbol","Entrez","GO")
merged_go=merge(data,result_go,by="Symbol")
write.table(merged_go,"human_knockout_go.txt",col.names = colnames(merged_go),quote=FALSE)

length(intersect(intersect(which(is.na(result[,3])==FALSE),which(is.na(result_pfocr[,3])==FALSE)),which(is.na(result_go[,3])==FALSE)))
38

dim(data)
dim(merged_go)
##########################################################################################knockout benchmark
require(pROC)

targets=which(is.na(merged_wp$WP)==FALSE)
target_gses=merged_wp[targets,2]
target_index=match(target_gses,colnames(pvalue_results_human_voom))

for(i in target_index){
  
  true_sig=merged_wp$set_id[which(merged_wp$GSE==target_gses[i])]
  
  response=(rsea_results_human_voom[[i]]$set_id%in%true_sig)
  roc_r=roc(rsea_results_human_voom[[i]]$Comp.adjP)
  auc(roc_r)
 
}







