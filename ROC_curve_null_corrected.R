##########################################################################################knockout benchmark
require(pROC)

setwd("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/")
pvalue_results_human_voom<-readRDS("pvalue_results_human_voom.rds")

library(dplyr)
library(magrittr)

non0gse=which(colSums(pvalue_results_human_voom,na.rm=TRUE)!=0)
length(non0gse)

rsea_results_human_voom_go=list()

## process result.txt files for counts
for(i in non0gse){
  rsea_results_human_voom_go[[i]]<-read.table(paste0("GSE/",colnames(pvalue_results_human_voom)[i],"/rSEA/GO/result.txt"),header=TRUE)
}

rsea_results_human_voom_wp=list()

## process result.txt files for counts
for(i in non0gse){
  rsea_results_human_voom_wp[[i]]<-read.table(paste0("GSE/",colnames(pvalue_results_human_voom)[i],"/rSEA/WP/result.txt"),header=TRUE)
}

rsea_results_human_voom_pfocr=list()

## process result.txt files for counts
for(i in non0gse){
  rsea_results_human_voom_pfocr[[i]]<-read.table(paste0("GSE/",colnames(pvalue_results_human_voom)[i],"/rSEA/PFOCR/result.txt"),header=TRUE)
}

rsea_results_human_voom_pfocr_3sets=list()

## process result.txt files for counts
for(i in non0gse){
  rsea_results_human_voom_pfocr_3sets[[i]]<-read.table(paste0("GSE/",colnames(pvalue_results_human_voom)[i],"/rSEA/PFOCR_3sets/result.txt"),header=TRUE)
}

rsea_results_human_voom_pfocr_miny=list()

## process result.txt files for counts
for(i in non0gse){
  rsea_results_human_voom_pfocr_miny[[i]]<-read.table(paste0("GSE/",colnames(pvalue_results_human_voom)[i],"/rSEA/PFOCR_miny/result.txt"),header=TRUE)
}

knockout_target_gse<-read.table("knockout_target_gse",header=FALSE)
View(which(colnames(pvalue_results_human_voom)%in%as.character(unlist(knockout_target_gse))))

database_lists2<-load("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/pfocr_miny_April2021.RData")
database_lists1<-load("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/databases.RData")#has wp, pfocr, go
#database_lists3<-load("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/databases_pfocr_3sets.RData")
database_lists3<-load("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/databases_pfocr_3intersect_v2.RData")
min_set_size <- 3
max_set_size <- 500 

database_lists2<- unname(unlist(sapply(database_lists2, grep, pattern="_list$", value = T, perl = T)))
for (db in database_lists2) {
  eval(call("<-", as.name(db),  Filter(Negate(is.null), lapply(get(db), function(x){
    if(length(x) < min_set_size | length(x) > max_set_size)
      NULL
    else
      x
  }))
  ))
}

database_lists1 <- unname(unlist(sapply(database_lists1, grep, pattern="_list$", value = T, perl = T)))
for (db in database_lists1) {
  eval(call("<-", as.name(db),  Filter(Negate(is.null), lapply(get(db), function(x){
    if(length(x) < min_set_size | length(x) > max_set_size)
      NULL
    else
      x
  }))
  ))
}

database_lists3<- unname(unlist(sapply(database_lists3, grep, pattern="_list$", value = T, perl = T)))
for (db in database_lists3) {
  eval(call("<-", as.name(db),  Filter(Negate(is.null), lapply(get(db), function(x){
    if(length(x) < min_set_size | length(x) > max_set_size)
      NULL
    else
      x
  }))
  ))
}

########GO
#go match

data<-read_excel("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/human_knockout_genelist_alias_updated_April2021.xlsx", col_names = F)
data=as.data.frame(data)
colnames(data)=c("GSE","Symbol")

data_gene=as.character(unlist(data[,2]))


gene_entrez <- bitr(data_gene, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db, drop=FALSE)


colnames(gene_entrez)[1]="Gene"
gene_entrez=unique(gene_entrez)


tb=table(gene_entrez$Gene)
tb[which(tb>1)]


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
44

colnames(result_go)=c("Symbol","Entrez","GO")
merged_go=merge(data,result_go,by="Symbol")

nrow(merged_go)-sum(is.na(merged_go[,3]))
53
write.table(merged_go,"human_knockout_go_April2021version.txt",col.names = colnames(merged_go),quote=FALSE)

#go match #end
merged_go<-read.table("/Users/mingyoungshin/git/PFOCRInPathwayAnalyses/human_knockout_go.txt",header=TRUE)
merged_go<-read.table("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/human_knockout_go_April2021version.txt",header=TRUE)
targets=which(is.na(merged_go$GO)==FALSE)#get GSEs matched with knock genes
target_gses=merged_go[targets,2]
length(target_gses)
target_gses=target_gses[-which(target_gses%in%c("GSE75291","GSE78506","GSE141126","GSE76421"))]
length(target_gses)
target_index=match(target_gses,colnames(pvalue_results_human_voom))#retrieve their index
target_index_go=target_index

rsea_go_sig_gene_count=NULL
rsea_go_aucs=NULL
rsea_go_aucs_ci=NULL
rsea_go_50lower=NULL
rsea_go_TRUEs=NULL
for(i in unique(target_index)){
  print(i)
  rsea_go_sig_gene_count=c(rsea_go_sig_gene_count,sum(rsea_results_human_voom_go[[i]]$Com.adjP <0.05,na.rm=TRUE))
  
  if( sum(is.na(rsea_results_human_voom_go[[i]])) ==sum(lengths(rsea_results_human_voom_go[[i]])) ){next}#pass null cases
  gse_match=which(as.character(merged_go$GSE)==colnames(pvalue_results_human_voom)[i])
  for(j in gse_match){
    print(merged_go$GO[j])
    if(is.na(merged_go$GO[j])){next}
    
    true_sig=merged_go$GO[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    true_sig=true_sig[true_sig%in%(names(go_list))]
    print(length(true_sig))
    response=apply(as.matrix(1:length(rsea_results_human_voom_go[[i]]$set_id)),1,function(x) rsea_results_human_voom_go[[i]]$set_id[x]%in%true_sig)
    print(sum(response))
    print(table(response))
    rsea_go_TRUEs=rbind(rsea_go_TRUEs,c(i,merged_go$Symbol[j],as.numeric(table(response)['TRUE'])))
    
    if(length(table(response))==1) next
    
   
    
    roc_r=roc(response,as.numeric(rsea_results_human_voom_go[[i]]$Comp.adjP))
    auc_r=auc(roc_r)
    
  
   
    
    rsea_go_aucs=rbind(rsea_go_aucs,c(i,auc_r[1]) )
  }
  
  
}


rsea_go_aucs=rsea_go_aucs[order(rsea_go_aucs[,1]),]




########wp
#####wp match
library(readxl)
library(org.Hs.eg.db)
library(clusterProfiler)

data<-read_excel("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/human_knockout_genelist_alias_updated_April2021.xlsx", col_names = F)


data_gene=as.character(unlist(data[,2]))


gene_entrez <- bitr(data_gene, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db, drop=FALSE)
# 13.56%% of input gene IDs are fail to map
head(gene_entrez)


colnames(gene_entrez)[1]="Gene"
gene_entrez=unique(gene_entrez)


tb=table(gene_entrez$Gene)
tb[which(tb>1)]


result_wp=matrix(NA,nrow=nrow(gene_entrez),ncol=3)
for(i in 1:nrow(result_wp)){
  result_wp[i,1]=gene_entrez$Gene[i]
  result_wp[i,2]=gene_entrez$ENTREZID[i]
  
  for(j in 1:length(wp_list)){
    temp=NULL
    if( result_wp[i,2]%in%wp_list[[j]]){
      
      #temp=paste(temp,wp_name,collapse =';') 
      temp=names(wp_list)[j]
    }
    if(length(temp)>0){result_wp[i,3]=paste(result_wp[i,3],temp,sep=";")}
    
  }
}


nrow(result_wp)-sum(is.na(result_wp[,3]))
29

result_wp=gsub("NA;","",result_wp)

colnames(result_wp)=c("Symbol","Entrez","WP")
colnames(data)=c("GSE","Symbol")
merged_wp=merge(data,result_wp,by="Symbol")
nrow(merged_wp)-sum(is.na(merged_wp[,3]))
53

write.table(merged_wp,"human_knockout_wp_April2021version.txt",col.names = colnames(merged_wp),quote=FALSE)



#####wp match #end


merged_wp<-read.table("/Users/mingyoungshin/git/PFOCRInPathwayAnalyses/human_knockout_wp.txt",header=TRUE)
merged_wp<-read.table("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/human_knockout_wp_April2021version.txt",header=TRUE)

targets=which(is.na(merged_wp$WP)==FALSE)#get GSEs matched with knock genes
target_gses=merged_wp[targets,2]
length(target_gses)
target_gses=target_gses[-which(target_gses%in%c("GSE75291","GSE78506","GSE141126","GSE76421"))]
length(target_gses)

target_index=match(target_gses,colnames(pvalue_results_human_voom))#retrieve their index
target_index_wp=target_index

rsea_wp_sig_gene_count=NULL
rsea_wp_aucs=NULL
rsea_wp_aucs_ci=NULL
rsea_wp_50lower=NULL
rsea_wp_TRUEs=NULL
for(i in unique(target_index)){
  print(i)
  rsea_wp_sig_gene_count=c(rsea_wp_sig_gene_count,sum(rsea_results_human_voom_wp[[i]]$Com.adjP <0.05,na.rm=TRUE))
  
  if( sum(is.na(rsea_results_human_voom_wp[[i]])) ==sum(lengths(rsea_results_human_voom_wp[[i]])) ){next}#pass null cases
  gse_match=which(as.character(merged_wp$GSE)==colnames(pvalue_results_human_voom)[i])
  for(j in gse_match){
    print(merged_wp$WP[j])
    if(is.na(merged_wp$WP[j])){next}
    
    true_sig=merged_wp$WP[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    true_sig=true_sig[true_sig%in%(names(wp_list))]
    print(length(true_sig))
    response=apply(as.matrix(1:length(rsea_results_human_voom_wp[[i]]$set_id)),1,function(x) rsea_results_human_voom_wp[[i]]$set_id[x]%in%true_sig)
    print(sum(response))
    print(table(response))
    
    rsea_wp_TRUEs=rbind(rsea_wp_TRUEs,c(i,merged_wp$Symbol[j],as.numeric(table(response)['TRUE'])))
    
    if(length(table(response))==1) next
    
    roc_r=roc(response,as.numeric(rsea_results_human_voom_wp[[i]]$Comp.adjP))
    auc_r=auc(roc_r)
    if(auc_r<=0.5){rsea_wp_50lower=rbind(rsea_wp_50lower,c(i,j))}
    
    auc_r_ci=ci(auc_r)
    if(!is.na(auc_r_ci[1])){    if(auc_r_ci[1] < roc_r$auc[1] & auc_r_ci[2] > roc_r$auc[1]){
      rsea_wp_aucs_ci=c(rsea_wp_aucs_ci,auc_r[1])
    }}

 
    rsea_wp_aucs=rbind(rsea_wp_aucs,c(i,auc_r[1]) )
  }
  
  
}


boxplot(rsea_wp_aucs)



########pfocr
#pfocr match


result_pfocr=matrix(NA,nrow=nrow(gene_entrez),ncol=3)
for(i in 1:nrow(result_pfocr)){
  result_pfocr[i,1]=gene_entrez$Gene[i]
  result_pfocr[i,2]=gene_entrez$ENTREZID[i]
  
  for(j in 1:length(pfocr_list)){
    temp=NULL
    if( result_pfocr[i,2]%in%pfocr_list[[j]]){
      
      #temp=paste(temp,pfocr_name,collapse =';') 
      temp=names(pfocr_list)[j]
    }
    if(length(temp)>0){result_pfocr[i,3]=paste(result_pfocr[i,3],temp,sep=";")}
    
  }
}

result_pfocr=gsub("NA;","",result_pfocr)

nrow(result_pfocr)-sum(is.na(result_pfocr[,3]))
40

colnames(result_pfocr)=c("Symbol","Entrez","PFOCR")
merged_pfocr=merge(data,result_pfocr,by="Symbol")

nrow(merged_pfocr)-sum(is.na(merged_pfocr[,3]))
53
write.table(merged_pfocr,"human_knockout_pfocr_April2021version.txt",col.names = colnames(merged_pfocr),quote=FALSE)

#pfocr match end



merged_pfocr<-read.table("/Users/mingyoungshin/git/PFOCRInPathwayAnalyses/human_knockout_pfocr.txt",header=TRUE)
merged_pfocr<-read.table("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/human_knockout_pfocr_April2021version.txt",header=TRUE)

targets=which(is.na(merged_pfocr$PFOCR)==FALSE)#get GSEs matched with knock genes
target_gses=merged_pfocr[targets,2]
length(target_gses)
target_gses=target_gses[-which(target_gses%in%c("GSE75291","GSE78506","GSE141126","GSE76421"))]
length(target_gses)

target_index=match(target_gses,colnames(pvalue_results_human_voom))#retrieve their index
target_index_pfocr=target_index

rsea_pfocr_sig_gene_count=NULL
rsea_pfocr_aucs=NULL
rsea_pfocr_aucs_ci=NULL
rsea_pfocr_50lower=NULL
rsea_pfocr_TRUEs=NULL
for(i in unique(target_index)){
  print(i)
  rsea_pfocr_sig_gene_count=c(rsea_pfocr_sig_gene_count,sum(rsea_results_human_voom_pfocr[[i]]$Com.adjP <0.05,na.rm=TRUE))
  
  if( sum(is.na(rsea_results_human_voom_pfocr[[i]])) ==sum(lengths(rsea_results_human_voom_pfocr[[i]])) ){next}#pass null cases
  gse_match=which(as.character(merged_pfocr$GSE)==colnames(pvalue_results_human_voom)[i])
  for(j in gse_match){
    print(merged_pfocr$PFOCR[j])
    if(is.na(merged_pfocr$PFOCR[j])){next}
    
    true_sig=merged_pfocr$PFOCR[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    true_sig=true_sig[true_sig%in%(names(pfocr_list))]
    print(length(true_sig))
    response=apply(as.matrix(1:length(rsea_results_human_voom_pfocr[[i]]$set_id)),1,function(x) rsea_results_human_voom_pfocr[[i]]$set_id[x]%in%true_sig)
    print(sum(response))
    print(table(response))
    rsea_pfocr_TRUEs=rbind(rsea_pfocr_TRUEs,c(i,merged_pfocr$Symbol[j],as.numeric(table(response)['TRUE'])))
    
    if(length(table(response))==1) next
    
    roc_r=roc(response,as.numeric(rsea_results_human_voom_pfocr[[i]]$Comp.adjP))
    auc_r=auc(roc_r)
    if(auc_r<=0.5){rsea_pfocr_50lower=rbind(rsea_pfocr_50lower,c(i,j))}
    
    auc_r_ci=ci(auc_r)
    if(!is.na(auc_r_ci[1])){    if(auc_r_ci[1] < roc_r$auc[1] & auc_r_ci[2] > roc_r$auc[1]){
      rsea_pfocr_aucs_ci=c(rsea_pfocr_aucs_ci,auc_r[1])
    }}

   
    rsea_pfocr_aucs=rbind(rsea_pfocr_aucs,c(i,auc_r[1]))
  }
  
  
}


boxplot(rsea_pfocr_aucs)


#######pfocr 3set match
library(readxl)
library(org.Hs.eg.db)

data<-read_excel("/Users/mingyoungshin/git/PFOCRInPathwayAnalyses/human_knockout_genelist_alias_updated.xlsx", col_names = F)
data_gene=as.character(unlist(data[,2]))

gene_entrez <- bitr(data_gene, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db, drop=FALSE)


colnames(gene_entrez)[1]="Gene"
gene_entrez=unique(gene_entrez)


tb=table(gene_entrez$Gene)
tb[which(tb>1)]


result_pfocr=matrix(NA,nrow=nrow(gene_entrez),ncol=3)
for(i in 1:nrow(result_pfocr)){
  result_pfocr[i,1]=gene_entrez$Gene[i]
  result_pfocr[i,2]=gene_entrez$ENTREZID[i]
  
  for(j in 1:length(pfocr_3sets_list)){
    temp=NULL
    if( result_pfocr[i,2]%in%pfocr_3sets_list[[j]]){
      temp=names(pfocr_3sets_list)[j]
    }
    if(length(temp)>0){result_pfocr[i,3]=paste(result_pfocr[i,3],temp,sep=";")}
    
  }
}

result_pfocr=gsub("NA;","",result_pfocr)

nrow(result_pfocr)-sum(is.na(result_pfocr[,3]))
65
colnames(result_pfocr)=c("Symbol","Entrez","PFOCR")
colnames(data)=c("GSE","Symbol")
merged_pfocr=merge(data,result_pfocr,by="Symbol")

nrow(merged_pfocr)-sum(is.na(merged_pfocr[,3]))
85
write.table(merged_pfocr,"human_knockout_pfocr_3sets.txt",col.names = colnames(merged_pfocr),quote=FALSE)




library(readxl)
library(org.Hs.eg.db)

data<-read_excel("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/human_knockout_genelist_alias_updated_April2021.xlsx", col_names = F)
data_gene=as.character(unlist(data[,2]))

gene_entrez <- bitr(data_gene, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db, drop=FALSE)


colnames(gene_entrez)[1]="Gene"
gene_entrez=unique(gene_entrez)


tb=table(gene_entrez$Gene)
tb[which(tb>1)]


result_pfocr=matrix(NA,nrow=nrow(gene_entrez),ncol=3)
for(i in 1:nrow(result_pfocr)){
  result_pfocr[i,1]=gene_entrez$Gene[i]
  result_pfocr[i,2]=gene_entrez$ENTREZID[i]
  
  for(j in 1:length(pfocr_3sets_list)){
    temp=NULL
    if( result_pfocr[i,2]%in%pfocr_3sets_list[[j]]){
      temp=names(pfocr_3sets_list)[j]
    }
    if(length(temp)>0){result_pfocr[i,3]=paste(result_pfocr[i,3],temp,sep=";")}
    
  }
}

result_pfocr=gsub("NA;","",result_pfocr)

nrow(result_pfocr)-sum(is.na(result_pfocr[,3]))
41
colnames(result_pfocr)=c("Symbol","Entrez","PFOCR")
colnames(data)=c("GSE","Symbol")
merged_pfocr=merge(data,result_pfocr,by="Symbol")

nrow(merged_pfocr)-sum(is.na(merged_pfocr[,3]))
53
write.table(merged_pfocr,"human_knockout_pfocr_3sets_April2021version.txt",col.names = colnames(merged_pfocr),quote=FALSE)

#######pfocr 3set match#end




merged_pfocr_3sets<-read.table("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/human_knockout_pfocr_3sets.txt",header=TRUE)
merged_pfocr_3sets<-read.table("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/human_knockout_pfocr_3sets_April2021version.txt",header=TRUE)

targets=which(is.na(merged_pfocr_3sets$PFOCR)==FALSE)#get GSEs matched with knock genes
target_gses=merged_pfocr_3sets[targets,2]
length(target_gses)
target_gses=target_gses[-which(target_gses%in%c("GSE75291","GSE78506","GSE141126","GSE76421"))]
length(target_gses)

target_index=match(target_gses,colnames(pvalue_results_human_voom))#retrieve their index
target_index_pfocr_3sets=target_index

rsea_pfocr_3sets_sig_gene_count=NULL
rsea_pfocr_3sets_aucs=NULL
rsea_pfocr_3sets_aucs_ci=NULL
rsea_pfocr_3set_50lower=NULL
rsea_pfocr_3rsets_TRUEs=NULL
for(i in unique(target_index)){
  print(i)
  rsea_pfocr_3sets_sig_gene_count=c(rsea_pfocr_3sets_sig_gene_count,sum(rsea_results_human_voom_pfocr_3sets[[i]]$Com.adjP <0.05,na.rm=TRUE))
  
  if( sum(is.na(rsea_results_human_voom_pfocr_3sets[[i]])) ==sum(lengths(rsea_results_human_voom_pfocr_3sets[[i]])) ){next}#pass null cases
  gse_match=which(as.character(merged_pfocr_3sets$GSE)==colnames(pvalue_results_human_voom)[i])
  for(j in gse_match){
    print(merged_pfocr_3sets$PFOCR[j])
    if(is.na(merged_pfocr_3sets$PFOCR[j])){next}
    
    true_sig=merged_pfocr_3sets$PFOCR[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    true_sig=true_sig[true_sig%in%(names(pfocr_3sets_list))]
    print(length(true_sig))
    response=apply(as.matrix(1:length(rsea_results_human_voom_pfocr_3sets[[i]]$set_id)),1,function(x) rsea_results_human_voom_pfocr_3sets[[i]]$set_id[x]%in%true_sig)
    print(sum(response))
    print(table(response))
    rsea_pfocr_3rsets_TRUEs=rbind(rsea_pfocr_3rsets_TRUEs,c(i,merged_pfocr_3sets$Symbol[j],as.numeric(table(response)['TRUE'])))
    
    if(length(table(response))==1) next
    
    roc_r=roc(response,as.numeric(rsea_results_human_voom_pfocr_3sets[[i]]$Comp.adjP))
    auc_r=auc(roc_r)
    if(auc_r<=0.5){rsea_pfocr_3set_50lower=rbind(rsea_pfocr_3set_50lower,c(i,j))}
    
    auc_r_ci=ci(auc_r)
    if(!is.na(auc_r_ci[1])){
      if(auc_r_ci[1] < roc_r$auc[1] & auc_r_ci[2] > roc_r$auc[1]){
        rsea_pfocr_3sets_aucs_ci=c(rsea_pfocr_3sets_aucs_ci,auc_r[1])
      }
    }

   
    rsea_pfocr_3sets_aucs=rbind(rsea_pfocr_3sets_aucs,c(i,auc_r[1]) )
  }
  
  
}




#######pfocr miny match
library(readxl)
library(org.Hs.eg.db)

data<-read_excel("/Users/mingyoungshin/git/PFOCRInPathwayAnalyses/human_knockout_genelist_alias_updated.xlsx", col_names = F)
data_gene=as.character(unlist(data[,2]))

gene_entrez <- bitr(data_gene, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db, drop=FALSE)


colnames(gene_entrez)[1]="Gene"
gene_entrez=unique(gene_entrez)


tb=table(gene_entrez$Gene)
tb[which(tb>1)]


result_pfocr=matrix(NA,nrow=nrow(gene_entrez),ncol=3)
for(i in 1:nrow(result_pfocr)){
  result_pfocr[i,1]=gene_entrez$Gene[i]
  result_pfocr[i,2]=gene_entrez$ENTREZID[i]
  
  for(j in 1:length(pfocr_miny_list)){
    temp=NULL
    if( result_pfocr[i,2]%in%pfocr_miny_list[[j]]){
      temp=names(pfocr_miny_list)[j]
    }
    if(length(temp)>0){result_pfocr[i,3]=paste(result_pfocr[i,3],temp,sep=";")}
    
  }
}

result_pfocr=gsub("NA;","",result_pfocr)

nrow(result_pfocr)-sum(is.na(result_pfocr[,3]))
41
colnames(result_pfocr)=c("Symbol","Entrez","PFOCR")
colnames(data)=c("GSE","Symbol")
merged_pfocr=merge(data,result_pfocr,by="Symbol")

nrow(merged_pfocr)-sum(is.na(merged_pfocr[,3]))
53

write.table(merged_pfocr,"human_knockout_pfocr_miny_April2021version.txt",col.names = colnames(merged_pfocr),quote=FALSE)

#######pfocr miny match#end

merged_pfocr_miny<-read.table("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/human_knockout_pfocr_miny.txt",header=TRUE)
merged_pfocr_miny<-read.table("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/human_knockout_pfocr_miny_April2021version.txt",header=TRUE)

targets=which(is.na(merged_pfocr_miny$PFOCR)==FALSE)#get GSEs matched with knock genes
target_gses=merged_pfocr_miny[targets,2]
length(target_gses)
target_gses=target_gses[-which(target_gses%in%c("GSE75291","GSE78506","GSE141126","GSE76421"))]
length(target_gses)

target_index=match(target_gses,colnames(pvalue_results_human_voom))#retrieve their index
target_index_pfocr_miny=target_index

rsea_pfocr_miny_sig_gene_count=NULL
rsea_pfocr_miny_aucs=NULL
rsea_pfocr_miny_aucs_ci=NULL
rsea_pfocr_3set_50lower=NULL
rsea_pfocr_3rsets_TRUEs=NULL
for(i in unique(target_index)){
  print(i)
  rsea_pfocr_miny_sig_gene_count=c(rsea_pfocr_miny_sig_gene_count,sum(rsea_results_human_voom_pfocr_miny[[i]]$Com.adjP <0.05,na.rm=TRUE))
  
  if( sum(is.na(rsea_results_human_voom_pfocr_miny[[i]])) ==sum(lengths(rsea_results_human_voom_pfocr_miny[[i]])) ){next}#pass null cases
  gse_match=which(as.character(merged_pfocr_miny$GSE)==colnames(pvalue_results_human_voom)[i])
  for(j in gse_match){
    print(merged_pfocr_miny$PFOCR[j])
    if(is.na(merged_pfocr_miny$PFOCR[j])){next}
    
    true_sig=merged_pfocr_miny$PFOCR[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    true_sig=true_sig[true_sig%in%(names(pfocr_miny_list))]
    print(length(true_sig))
    response=apply(as.matrix(1:length(rsea_results_human_voom_pfocr_miny[[i]]$set_id)),1,function(x) rsea_results_human_voom_pfocr_miny[[i]]$set_id[x]%in%true_sig)
    print(sum(response))
    print(table(response))
    rsea_pfocr_3rsets_TRUEs=rbind(rsea_pfocr_3rsets_TRUEs,c(i,merged_pfocr_miny$Symbol[j],as.numeric(table(response)['TRUE'])))
    
    if(length(table(response))==1) next
    
    roc_r=roc(response,as.numeric(rsea_results_human_voom_pfocr_miny[[i]]$Comp.adjP))
    auc_r=auc(roc_r)
    if(auc_r<=0.5){rsea_pfocr_3set_50lower=rbind(rsea_pfocr_3set_50lower,c(i,j))}
    
    auc_r_ci=ci(auc_r)
    if(!is.na(auc_r_ci[1])){
      if(auc_r_ci[1] < roc_r$auc[1] & auc_r_ci[2] > roc_r$auc[1]){
        rsea_pfocr_miny_aucs_ci=c(rsea_pfocr_miny_aucs_ci,auc_r[1])
      }
    }
    
   
    rsea_pfocr_miny_aucs=rbind(rsea_pfocr_miny_aucs,c(i,auc_r[1]) )
  }
  
  
}


save.image("ROC.RData")


#################compare with null roc
library(dplyr)


permuted_auc_go.result<-read.table("permuted_auc_go.result",header=FALSE)
colnames(permuted_auc_go.result)=c("index","mean","dev")

rsea_go_aucs
colnames(rsea_go_aucs)=c("index","auc")

rsea_go_merged=merge(rsea_go_aucs,permuted_auc_go.result,by="index")
rsea_go_merged %<>% mutate(corrected_auc= (auc-mean)/dev)

permuted_auc_wp.result<-read.table("permuted_auc_wp.result",header=FALSE)
colnames(permuted_auc_wp.result)=c("index","mean","dev")

rsea_wp_aucs
colnames(rsea_wp_aucs)=c("index","auc")

rsea_wp_merged=merge(rsea_wp_aucs,permuted_auc_wp.result,by="index")
rsea_wp_merged %<>% mutate(corrected_auc= (auc-mean)/dev)

permuted_auc_pfocr.result<-read.table("permuted_auc_pfocr.result",header=FALSE)
colnames(permuted_auc_pfocr.result)=c("index","mean","dev")

rsea_pfocr_aucs
colnames(rsea_pfocr_aucs)=c("index","auc")

rsea_pfocr_merged=merge(rsea_pfocr_aucs,permuted_auc_pfocr.result,by="index")
rsea_pfocr_merged %<>% mutate(corrected_auc= (auc-mean)/dev)

permuted_auc_pfocr_3sets.result<-read.table("permuted_auc_pfocr_3sets.result",header=FALSE)
colnames(permuted_auc_pfocr_3sets.result)=c("index","mean","dev")

rsea_pfocr_3sets_aucs
colnames(rsea_pfocr_3sets_aucs)=c("index","auc")

rsea_pfocr_3sets_merged=merge(rsea_pfocr_3sets_aucs,permuted_auc_pfocr_3sets.result,by="index")
rsea_pfocr_3sets_merged %<>% mutate(corrected_auc= (auc-mean)/dev)


permuted_auc_pfocr_miny.result<-read.table("permuted_auc_pfocr_miny.result",header=FALSE)
colnames(permuted_auc_pfocr_miny.result)=c("index","mean","dev")

rsea_pfocr_miny_aucs
colnames(rsea_pfocr_miny_aucs)=c("index","auc")

rsea_pfocr_miny_merged=merge(rsea_pfocr_miny_aucs,permuted_auc_pfocr_miny.result,by="index")
rsea_pfocr_miny_merged %<>% mutate(corrected_auc= (auc-mean)/dev)
#################

library(ggplot2)

data=rbind( cbind(rsea_pfocr_miny_merged$corrected_auc,rep("PFOCR_miny",length(rsea_pfocr_miny_merged$corrected_auc))),
  cbind(rsea_pfocr_3sets_merged$corrected_auc,rep("PFOCR_3sets",length(rsea_pfocr_3sets_merged$corrected_auc))),
  cbind(rsea_pfocr_merged$corrected_auc,rep("PFOCR",length(rsea_pfocr_merged$corrected_auc))),
            cbind(rsea_wp_merged$corrected_auc,rep("WP",length(rsea_wp_merged$corrected_auc))),
            cbind(rsea_go_merged$corrected_auc,rep("GO",length(rsea_go_merged$corrected_auc))) )
colnames(data)=c("AUC","Data")

data=data.frame(data)
data$AUC=as.numeric(as.character(data$AUC))


# Change violin plot colors by groups
p<-ggplot(data, aes(x=Data, y=AUC, fill=Data)) +
  geom_boxplot() + scale_fill_manual(values=c("antiquewhite3", "darkseagreen", "darksalmon","lightblue","pink")) + theme_classic()
p

tiff("corrected_auc.tiff")
print(p)
dev.off()



p<-ggplot(subset(data,abs(AUC)<10), aes(x=Data, y=AUC, fill=Data)) +
  geom_boxplot() + scale_fill_manual(values=c("antiquewhite3", "darkseagreen", "darksalmon","lightblue","pink")) + theme_classic()
p

tiff("corrected_auc_abs_10.tiff")
print(p)
dev.off()

