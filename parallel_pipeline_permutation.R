#qsub -v start=1,end=5,method="ora",database="pfocr",set="set1" enrichment_parallel_permutation.sh 
require("rhdf5")
require("tools")
require('edgeR')
require(GEOquery)
require(rSEA)
require("clusterProfiler")
require("org.Hs.eg.db")
require(dplyr)
require(GO.db)
require(bayesbio)
require(readxl)
#require(pROC)

#pvalue_results_human_voom=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/pvalue_results_human_voom.rds")
#logFC_results_human_voom=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/logFC_results_human_voom.rds")
#merged_filtered_gses=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/merged_filtered_gses.rds")
pvalue_results_perm1_voom=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/pvalue_results_perm1_voom.rds")
logFC_results_perm1_voom=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/logFC_results_perm1_voom.rds")

pvalue_results_perm2_voom=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/pvalue_results_perm2_voom.rds")
logFC_results_perm2_voom=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/logFC_results_perm2_voom.rds")

print("loaded data")


database_lists1<-load("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/databases.RData")#has wp, pfocr, go
#database_lists2<-load("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/databases_new_pfocrs.RData")#has wp, pfocr, go

min_set_size <- 10
max_set_size <- 500 

 
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
# 
# database_lists2 <- unname(unlist(sapply(database_lists2, grep, pattern="_list$", value = T, perl = T)))
# for (db in database_lists2) {
#   eval(call("<-", as.name(db),  Filter(Negate(is.null), lapply(get(db), function(x){
#     if(length(x) < min_set_size | length(x) > max_set_size)
#       NULL
#     else
#       x
#   }))
#   ))
# }

print("loaded databses")

gene_entrez=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/gene_entrez.rds")

print("loaded entrez")

source("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/run_rSEA4.r")
source("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/run_gsea4.r")
source("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/run_ORA5.r")

print("loaded functions")
#args=c(1,536)
args <- commandArgs(trailingOnly=TRUE)

print(args)



rsea_results <- vector("list", 7)
names(rsea_results) <- c("set_id", "ID", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP")


ora_results <- vector("list", 6)
names(ora_results) <- c("GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue","Count")


gsea_results <- vector("list", 6)
names(gsea_results) <- c("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalues","rank")


print("starting process")
print(Sys.time())
range=as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2])))

if(args[5]=="set1"){
  
  if(args[3]=="gsea"){
    if(args[4]=="pfocr"){
      gsea_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_gsea4(as.matrix(pvalue_results_perm1_voom[,x]),as.matrix(logFC_results_perm1_voomm[,x]),gsea_results,pfocr_list,pfocr_annotation)})
      saveRDS(gsea_results,paste0("rsea_results_set1_pfocr","_",args[1],"_",args[2],".rds"))#11hours in total
      apply( as.matrix(1:length(gsea_results)),1, function(x) write.table(as.data.frame(gsea_results[[x]]),
                                                                          file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set1/perm",  as.numeric(as.character(args[1]))+(x-1), "/GSEA/PFOCR/result.txt"),
                                                                          col.names=names(gsea_results[[x]]),row.names=rownames(as.data.frame(gsea_results[[x]])),quote=FALSE))
    }else if(args[4]=="wp"){
      gsea_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_gsea4(as.matrix(pvalue_results_perm1_voom[,x]),as.matrix(logFC_results_perm1_voom[,x]),gsea_results,wp_list,wp_annotation)})
      saveRDS(gsea_results,paste0("gsea_results_set1_wp","_",args[1],"_",args[2],".rds"))
      apply( as.matrix(1:length(gsea_results)),1, function(x) write.table(as.data.frame(gsea_results[[x]]),
                                                                          file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set1/perm", as.numeric(as.character(args[1]))+(x-1), "/GSEA/WP/result.txt"),
                                                                          col.names=names(gsea_results[[x]]),row.names=rownames(as.data.frame(gsea_results[[x]])),quote=FALSE))
      
    }else if(args[4]=="go"){
      gsea_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_gsea4(as.matrix(pvalue_results_perm1_voom[,x]),as.matrix(logFC_results_perm1_voom[,x]),gsea_results,go_list,go_annotation)})
      saveRDS(gsea_results,paste0("gsea_results_set1_go","_",args[1],"_",args[2],".rds"))
      apply( as.matrix(1:length(gsea_results)),1, function(x) write.table(as.data.frame(gsea_results[[x]]),
                                                                          file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set1/perm", as.numeric(as.character(args[1]))+(x-1), "/GSEA/WP/result.txt"),
                                                                          col.names=names(gsea_results[[x]]),row.names=rownames(as.data.frame(gsea_results[[x]])),quote=FALSE))
    }
    
  }else if(args[3]=="rsea"){
    if(args[4]=="pfocr"){
      rsea_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_rSEA4(as.matrix(pvalue_results_perm1_voom[,x]),rsea_results,pfocr_list,pfocr_annotation)})
      saveRDS(rsea_results,paste0("rsea_results_set1_pfocr","_",args[1],"_",args[2],".rds"))#28hours
      apply( as.matrix(1:length(rsea_results)),1, function(x) write.table(as.data.frame(rsea_results[[x]]),
                                                                          file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set1/perm",   as.numeric(as.character(args[1]))+(x-1), "/rSEA/PFOCR/result.txt"),
                                                                          col.names=names(rsea_results[[x]]),row.names=rownames(as.data.frame(rsea_results[[x]])),quote=FALSE))
    }else if(args[4]=="wp"){
      rsea_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_rSEA4(as.matrix(pvalue_results_perm1_voom[,x]),rsea_results_human_voom_wp,wp_list,wp_annotation)})
      saveRDS(rsea_results,paste0("rsea_results_set1_wp","_",args[1],"_",args[2],".rds"))
      apply( as.matrix(1:length(rsea_results)),1, function(x) write.table(as.data.frame(rsea_results[[x]]),
                                                                          file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set1/perm",  as.numeric(as.character(args[1]))+(x-1), "/rSEA/WP/result.txt"),
                                                                          col.names=names(rsea_results[[x]]),row.names=rownames(as.data.frame(rsea_results[[x]])),quote=FALSE))
      
    }else if(args[4]=="go"){
      rsea_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_rSEA4(as.matrix(pvalue_results_perm1_voom[,x]),rsea_results,go_list,go_annotation)})
      saveRDS(rsea_results,paste0("rsea_results_set1_go","_",args[1],"_",args[2],".rds"))#20-21 hours
      apply( as.matrix(1:length(rsea_results)),1, function(x) write.table(as.data.frame(rsea_results[[x]]),
                                                                          file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set1/perm",  as.numeric(as.character(args[1]))+(x-1), "/rSEA/GO/result.txt"),
                                                                          col.names=names(rsea_results[[x]]),row.names=rownames(as.data.frame(rsea_results[[x]])),quote=FALSE))
      
    }
  }else if(args[3]=="ora"){
    if(args[4]=="pfocr"){
      ora_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_ORA5(as.matrix(pvalue_results_perm1_voom[,x]),as.matrix(logFC_results_perm1_voom[,x]),ora_results,pfocr_list,pfocr_annotation)})
      saveRDS(ora_results,paste0("ora_results_set1_pfocr","_",args[1],"_",args[2],".rds"))
      apply( as.matrix(1:length(ora_results)),1, function(x) write.table(as.data.frame(ora_results[[x]]),
                                                                         file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set1/perm",   as.numeric(as.character(args[1]))+(x-1), "/ORA/PFOCR/result.txt"),
                                                                         col.names=names(ora_results[[x]]),row.names=rownames(as.data.frame(ora_results[[x]])),quote=FALSE))
      
    }else if(args[4]=="wp"){
      ora_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_ORA5(as.matrix(pvalue_results_perm1_voom[,x]),as.matrix(logFC_results_perm1_voom[,x]),ora_results,wp_list,wp_annotation)})
      saveRDS(ora_results,paste0("ora_results_set1_wp","_",args[1],"_",args[2],".rds"))
      apply( as.matrix(1:length(ora_results)),1, function(x) write.table(as.data.frame(ora_results[[x]]),
                                                                         file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set1/perm",   as.numeric(as.character(args[1]))+(x-1), "/ORA/WP/result.txt"),
                                                                         col.names=names(ora_results[[x]]),row.names=rownames(as.data.frame(ora_results[[x]])),quote=FALSE))
      
    }else if(args[4]=="go"){
      ora_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_ORA5(as.matrix(pvalue_results_perm1_voom[,x]),as.matrix(logFC_results_perm1_voom[,x]),ora_results,go_list,go_annotation)})
      saveRDS(ora_results,paste0("ora_results_set1_go","_",args[1],"_",args[2],".rds"))
      apply( as.matrix(1:length(ora_results)),1, function(x) write.table(as.data.frame(ora_results[[x]]),
                                                                         file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set1/perm",   as.numeric(as.character(args[1]))+(x-1), "/ORA/GO/result.txt"),
                                                                         col.names=names(ora_results[[x]]),row.names=rownames(as.data.frame(ora_results[[x]])),quote=FALSE))
      
    }
  }
  
}else if(args[5]=="set2"){
  
  if(args[3]=="gsea"){
    if(args[4]=="pfocr"){
      gsea_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_gsea3(as.matrix(pvalue_results_perm2_voom[,x]),as.matrix(logFC_results_perm2_voomm[,x]),gsea_results,pfocr_list,pfocr_annotation)})
      saveRDS(gsea_results,paste0("rsea_results_set2_pfocr","_",args[1],"_",args[2],".rds"))#11hours in total
      apply( as.matrix(1:length(gsea_results)),1, function(x) write.table(as.data.frame(gsea_results[[x]]),
                                                                          file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set2/perm",  as.numeric(as.character(args[1]))+(x-1), "/GSEA/PFOCR/result.txt"),
                                                                          col.names=names(gsea_results[[x]]),row.names=rownames(as.data.frame(gsea_results[[x]])),quote=FALSE))
    }else if(args[4]=="wp"){
      gsea_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_gsea3(as.matrix(pvalue_results_perm2_voom[,x]),as.matrix(logFC_results_perm2_voom[,x]),gsea_results,wp_list,wp_annotation)})
      saveRDS(gsea_results,paste0("gsea_results_set2_wp","_",args[1],"_",args[2],".rds"))
      apply( as.matrix(1:length(gsea_results)),1, function(x) write.table(as.data.frame(gsea_results[[x]]),
                                                                          file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set2/perm", as.numeric(as.character(args[1]))+(x-1), "/GSEA/WP/result.txt"),
                                                                          col.names=names(gsea_results[[x]]),row.names=rownames(as.data.frame(gsea_results[[x]])),quote=FALSE))
      
    }else if(args[4]=="go"){
      gsea_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_gsea3(as.matrix(pvalue_results_perm2_voom[,x]),as.matrix(logFC_results_perm2_voom[,x]),gsea_results,go_list,go_annotation)})
      saveRDS(gsea_results,paste0("gsea_results_set2_go","_",args[1],"_",args[2],".rds"))
      apply( as.matrix(1:length(gsea_results)),1, function(x) write.table(as.data.frame(gsea_results[[x]]),
                                                                          file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set2/perm", as.numeric(as.character(args[1]))+(x-1), "/GSEA/WP/result.txt"),
                                                                          col.names=names(gsea_results[[x]]),row.names=rownames(as.data.frame(gsea_results[[x]])),quote=FALSE))
    }
    
  }else if(args[3]=="rsea"){
    if(args[4]=="pfocr"){
      rsea_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_rSEA4(as.matrix(pvalue_results_perm2_voom[,x]),rsea_results,pfocr_list,pfocr_annotation)})
      saveRDS(rsea_results,paste0("rsea_results_set2_pfocr","_",args[1],"_",args[2],".rds"))#28hours
      apply( as.matrix(1:length(rsea_results)),1, function(x) write.table(as.data.frame(rsea_results[[x]]),
                                                                          file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set2/perm",   as.numeric(as.character(args[1]))+(x-1), "/rSEA/PFOCR/result.txt"),
                                                                          col.names=names(rsea_results[[x]]),row.names=rownames(as.data.frame(rsea_results[[x]])),quote=FALSE))
    }else if(args[4]=="wp"){
      rsea_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_rSEA4(as.matrix(pvalue_results_perm2_voom[,x]),rsea_results_human_voom_wp,wp_list,wp_annotation)})
      saveRDS(rsea_results,paste0("rsea_results_set2_wp","_",args[1],"_",args[2],".rds"))
      apply( as.matrix(1:length(rsea_results)),1, function(x) write.table(as.data.frame(rsea_results[[x]]),
                                                                          file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set2/perm",  as.numeric(as.character(args[1]))+(x-1), "/rSEA/WP/result.txt"),
                                                                          col.names=names(rsea_results[[x]]),row.names=rownames(as.data.frame(rsea_results[[x]])),quote=FALSE))
      
    }else if(args[4]=="go"){
      rsea_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_rSEA4(as.matrix(pvalue_results_perm2_voom[,x]),rsea_results,go_list,go_annotation)})
      saveRDS(rsea_results,paste0("rsea_results_set2_go","_",args[1],"_",args[2],".rds"))#20-21 hours
      apply( as.matrix(1:length(rsea_results)),1, function(x) write.table(as.data.frame(rsea_results[[x]]),
                                                                          file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set2/perm",  as.numeric(as.character(args[1]))+(x-1), "/rSEA/GO/result.txt"),
                                                                          col.names=names(rsea_results[[x]]),row.names=rownames(as.data.frame(rsea_results[[x]])),quote=FALSE))
      
    }
  }else if(args[3]=="ora"){
    if(args[4]=="pfocr"){
      ora_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_ORA4(as.matrix(pvalue_results_perm2_voom[,x]),as.matrix(logFC_results_perm2_voom[,x]),ora_results,pfocr_list,pfocr_annotation)})
      saveRDS(ora_results,paste0("ora_results_set2_pfocr","_",args[1],"_",args[2],".rds"))
      apply( as.matrix(1:length(ora_results)),1, function(x) write.table(as.data.frame(ora_results[[x]]),
                                                                         file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set2/perm",   as.numeric(as.character(args[1]))+(x-1), "/ORA/PFOCR/result.txt"),
                                                                         col.names=names(ora_results[[x]]),row.names=rownames(as.data.frame(ora_results[[x]])),quote=FALSE))
      
    }else if(args[4]=="wp"){
      ora_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_ORA4(as.matrix(pvalue_results_perm2_voom[,x]),as.matrix(logFC_results_perm2_voom[,x]),ora_results,wp_list,wp_annotation)})
      saveRDS(ora_results,paste0("ora_results_set2_wp","_",args[1],"_",args[2],".rds"))
      apply( as.matrix(1:length(ora_results)),1, function(x) write.table(as.data.frame(ora_results[[x]]),
                                                                         file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set2/perm",   as.numeric(as.character(args[1]))+(x-1), "/ORA/WP/result.txt"),
                                                                         col.names=names(ora_results[[x]]),row.names=rownames(as.data.frame(ora_results[[x]])),quote=FALSE))
      
    }else if(args[4]=="go"){
      ora_results=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){print(x);run_ORA4(as.matrix(pvalue_results_perm2_voom[,x]),as.matrix(logFC_results_perm2_voom[,x]),ora_results,go_list,go_annotation)})
      saveRDS(ora_results,paste0("ora_results_set2_go","_",args[1],"_",args[2],".rds"))
      apply( as.matrix(1:length(ora_results)),1, function(x) write.table(as.data.frame(ora_results[[x]]),
                                                                         file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/permuted_GSE/set2/perm",   as.numeric(as.character(args[1]))+(x-1), "/ORA/GO/result.txt"),
                                                                         col.names=names(ora_results[[x]]),row.names=rownames(as.data.frame(ora_results[[x]])),quote=FALSE))
      
    }
  }
  
}


print("ending process")
print(Sys.time())
q(save='no')

##############databses

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






apply(as.matrix(1:length(merged_filtered_gses)),1,function(x){ dir.create(file.path(paste0("./GSE/", merged_filtered_gses[[x]]$GSE)), showWarnings = FALSE);
  dir.create(file.path(paste0("./GSE/", merged_filtered_gses[[x]]$GSE),"/rSEA"), showWarnings = FALSE);
  dir.create(file.path(paste0("./GSE/", merged_filtered_gses[[x]]$GSE),"/rSEA/WP"), showWarnings = FALSE);
  dir.create(file.path(paste0("./GSE/", merged_filtered_gses[[x]]$GSE),"/rSEA/PFOCR"), showWarnings = FALSE);
  dir.create(file.path(paste0("./GSE/", merged_filtered_gses[[x]]$GSE),"/rSEA/GO"), showWarnings = FALSE);
  dir.create(file.path(paste0("./GSE/", merged_filtered_gses[[x]]$GSE),"/ORA"), showWarnings = FALSE);
  dir.create(file.path(paste0("./GSE/", merged_filtered_gses[[x]]$GSE),"/ORA/WP"), showWarnings = FALSE);
  dir.create(file.path(paste0("./GSE/", merged_filtered_gses[[x]]$GSE),"/ORA/PFOCR"), showWarnings = FALSE);
  dir.create(file.path(paste0("./GSE/", merged_filtered_gses[[x]]$GSE),"/ORA/GO"), showWarnings = FALSE);
  dir.create(file.path(paste0("./GSE/", merged_filtered_gses[[x]]$GSE),"/GSEA"), showWarnings = FALSE);
  dir.create(file.path(paste0("./GSE/", merged_filtered_gses[[x]]$GSE),"/GSEA/WP"), showWarnings = FALSE);
  dir.create(file.path(paste0("./GSE/", merged_filtered_gses[[x]]$GSE),"/GSEA/GO"), showWarnings = FALSE);
  dir.create(file.path(paste0("./GSE/", merged_filtered_gses[[x]]$GSE),"/GSEA/PFOCR"), showWarnings = FALSE)})


apply(as.matrix(1:length(merged_filtered_gses)),1, function(x) write.table(as.data.frame(ora_results_human_voom_wp[[x]])[c("GeneRatio","BgRatio","pvalue","p.adjust","qvalue","Count")],
                                                                           file=paste0("./GSE/",  merged_filtered_gses[[x]]$GSE, "/ORA/WP/result.txt"),
                                                                                       colnames(merged_filtered_gses)[x]),col.names=c("GeneRatio","BgRatio","pvalue","p.adjust","qvalue","Count"),row.names=rownames(as.data.frame(ora_results_human_voom_wp[[x]])),quote=FALSE))
apply(as.matrix(1:length(merged_filtered_gses)),1, function(x) write.table(as.data.frame(ora_results_human_voom_go[[x]])[c("GeneRatio","BgRatio","pvalue","p.adjust","qvalue","Count")],
                                                                           file=paste0("./GSE/",  merged_filtered_gses[[x]]$GSE, "/ORA/GO/result.txt"),
                                                                                       colnames(merged_filtered_gses)[x]),col.names=c("GeneRatio","BgRatio","pvalue","p.adjust","qvalue","Count"),row.names=rownames(as.data.frame(ora_results_human_voom_go[[x]])),quote=FALSE))
apply(as.matrix(1:length(merged_filtered_gses)),1, function(x){ write.table(as.data.frame(ora_results_human_voom_pfocr[[x]]),
                                                                            file=paste0("./GSE/",  merged_filtered_gses[[x]]$GSE, "/ORA/PFOCR/result.txt"),
                                                                            col.names=names(ora_results_human_voom_pfocr[[x]]),row.names=rownames(as.data.frame(ora_results_human_voom_pfocr[[x]])),quote=FALSE));print(x)})


apply(as.matrix(1:length(merged_filtered_gses)),1, function(x){ print(x);write.table(  as.data.frame(rsea_results_human_voom_wp[[x]][c("set_id","ID","Coverage","TDP.bound","TDP.estimate","SC.adjP","Comp.adjP")]),
                                                                                       file=paste0("./GSE/",  merged_filtered_gses[[x]]$GSE, "/rSEA/WP/result.txt"),
                                                                                       col.names=c("set_id","ID","Coverage","TDP.bound","TDP.estimate","SC.adjP","Comp.adjP"),row.names=rownames(as.data.frame(rsea_results_human_voom_wp[[x]][c("set_id","ID","Coverage","TDP.bound","TDP.estimate","SC.adjP","Comp.adjP")])),quote=FALSE)})
apply(as.matrix(1:length(merged_filtered_gses)),1, function(x){ print(x);write.table(  as.data.frame(rsea_results_human_voom_go[[x]][c("set_id","ID","Coverage","TDP.bound","TDP.estimate","SC.adjP","Comp.adjP")]),
                                                                           file=paste0("./GSE/",  merged_filtered_gses[[x]]$GSE, "/rSEA/GO/result.txt"),
                                                                                       col.names=c("set_id","ID","Coverage","TDP.bound","TDP.estimate","SC.adjP","Comp.adjP"),row.names=rownames(as.data.frame(rsea_results_human_voom_go[[x]][c("set_id","ID","Coverage","TDP.bound","TDP.estimate","SC.adjP","Comp.adjP")])),quote=FALSE)})
apply(as.matrix(1:length(merged_filtered_gses)),1, function(x){ print(x);write.table(  as.data.frame(rsea_results_human_voom_pfocr[[x]][c("set_id","ID","Coverage","TDP.bound","TDP.estimate","SC.adjP","Comp.adjP")]),
                                                                                       file=paste0("./GSE/",  merged_filtered_gses[[x]]$GSE, "/rSEA/PFOCR/result.txt"),
                                                                                       col.names=c("set_id","ID","Coverage","TDP.bound","TDP.estimate","SC.adjP","Comp.adjP"),row.names=rownames(as.data.frame(rsea_results_human_voom_pfocr[[x]][c("set_id","ID","Coverage","TDP.bound","TDP.estimate","SC.adjP","Comp.adjP")])),quote=FALSE)})


apply(as.matrix(1:length(merged_filtered_gses)),1, function(x){ print(x);write.table(  as.data.frame(gsea_results_human_voom_wp[[x]]),
                                                                                       file=paste0("./GSE/",  merged_filtered_gses[[x]]$GSE, "/GSEA/WP/result.txt"),
                                                                                       col.names=colnames(as.data.frame(gsea_results_human_voom_wp[[x]])),row.names=rownames(as.data.frame(gsea_results_human_voom_wp[[x]])),quote=FALSE))})
apply(as.matrix(1:length(merged_filtered_gses)),1, function(x){ print(x);write.table(  as.data.frame(gsea_results_human_voom_go[[x]]),
                                                                                       file=paste0("./GSE/",  merged_filtered_gses[[x]]$GSE, "/GSEA/GO/result.txt"),
                                                                                       col.names=colnames(as.data.frame(gsea_results_human_voom_go[[x]])),row.names=rownames(as.data.frame(gsea_results_human_voom_go[[x]])),quote=FALSE))})
apply(as.matrix(1:length(merged_filtered_gses)),1, function(x){ print(x);write.table(  as.data.frame(gsea_results_human_voom_pfocr[[x]]),
                                                                                       file=paste0("./GSE/",  merged_filtered_gses[[x]]$GSE, "/GSEA/PFOCR/result.txt"),
                                                                                       col.names=colnames(as.data.frame(gsea_results_human_voom_pfocr[[x]])),row.names=rownames(as.data.frame(gsea_results_human_voom_pfocr[[x]])),quote=FALSE))})


apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum(is.na(gsea_results_human_voom_go[[x]]$p.adjust))) } )#max <2000 among 12272
apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum(is.na(gsea_results_human_voom_wp[[x]]$p.adjust))) } )#max <40 among 568
apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum(is.na(gsea_results_human_voom_pfocr[[x]]$p.adjust))) } )#max <400 among 32277

apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum(is.na(rsea_results_human_voom_go[[x]]$Comp.adjP))) } )#max 4090 among 12272
apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum(is.na(rsea_results_human_voom_wp[[x]]$Comp.adjP))) } )#max among 568
apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum(is.na(rsea_results_human_voom_pfocr[[x]]$Comp.adjP))) } )#max  among 32277

apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum(is.na(ora_results_human_voom_go[[x]]$p.adjust))) } )#max  among 12272
apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum(is.na(ora_results_human_voom_wp[[x]]$p.adjust))) } )#max 568 among 568
apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum(is.na(ora_results_human_voom_pfocr[[x]]$p.adjust))) } )#max 32277 among 32277


knockout_target_gse<-read.table("knockout_target_gse",header=FALSE)#intersect of three databases
knockout_target_gse=as.matrix(knockout_target_gse)

all_gse_titles=apply(as.matrix(1:length( merged_filtered_gses)), 1, function(x) merged_filtered_gses[[x]]$GSE)

target_index=match(knockout_target_gse,(all_gse_titles))#retrieve their index


#########ora
temp=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum( (ora_results_human_voom_go[[x]]$p.adjust) <0.05 ,na.rm=TRUE )) } )#0 significance
png("ora_results_human_voom_go_sig.png")
hist(temp,breaks=100)
dev.off()

temp=apply(as.matrix(target_index),1,function(x){ print(sum( (ora_results_human_voom_go[[x]]$p.adjust) <0.05 ,na.rm=TRUE )) } )#0 significance

png("ora_results_human_voom_go_sig_knockout.png")
hist(temp,breaks=100)
dev.off()



temp=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum( (ora_results_human_voom_wp[[x]]$p.adjust) <0.05 ,na.rm=TRUE )) } )#0 significance
png("ora_results_human_voom_wp_sig.png")
hist(temp,breaks=100)
dev.off()

temp=apply(as.matrix(target_index),1,function(x){ print(sum( (ora_results_human_voom_wp[[x]]$p.adjust) <0.05 ,na.rm=TRUE )) } )#0 significance

png("ora_results_human_voom_wp_sig_knockout.png")
hist(temp,breaks=100)
dev.off()


temp=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum( (ora_results_human_voom_pfocr[[x]]$p.adjust) <0.05 ,na.rm=TRUE )) } )#0 significance
png("ora_results_human_voom_pfocr_sig.png")
hist(temp,breaks=100)
dev.off()

temp=apply(as.matrix(target_index),1,function(x){ print(sum( (ora_results_human_voom_pfocr[[x]]$p.adjust) <0.05 ,na.rm=TRUE )) } )#0 significance

png("ora_results_human_voom_pfocr_sig_knockout.png")
hist(temp,breaks=100)
dev.off()



#########gsea
temp=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum( (gsea_results_human_voom_go[[x]]$p.adjust) <0.05 ,na.rm=TRUE )) } )#0 significance
png("gsea_results_human_voom_go_sig.png")
hist(temp,breaks=100)
dev.off()

temp=apply(as.matrix(target_index),1,function(x){ print(sum( (gsea_results_human_voom_go[[x]]$p.adjust) <0.05 ,na.rm=TRUE )) } )#0 significance

png("gsea_results_human_voom_go_sig_knockout.png")
hist(temp,breaks=100)
dev.off()



temp=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum( (gsea_results_human_voom_wp[[x]]$p.adjust) <0.05 ,na.rm=TRUE )) } )#0 significance
png("gsea_results_human_voom_wp_sig.png")
hist(temp,breaks=100)
dev.off()

temp=apply(as.matrix(target_index),1,function(x){ print(sum( (gsea_results_human_voom_wp[[x]]$p.adjust) <0.05 ,na.rm=TRUE )) } )#0 significance

png("gsea_results_human_voom_wp_sig_knockout.png")
hist(temp,breaks=100)
dev.off()


temp=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum( (gsea_results_human_voom_pfocr[[x]]$p.adjust) <0.05 ,na.rm=TRUE )) } )#0 significance
png("gsea_results_human_voom_pfocr_sig.png")
hist(temp,breaks=100)
dev.off()

temp=apply(as.matrix(target_index),1,function(x){ print(sum( (gsea_results_human_voom_pfocr[[x]]$p.adjust) <0.05 ,na.rm=TRUE )) } )#0 significance

png("gsea_results_human_voom_pfocr_sig_knockout.png")
hist(temp,breaks=100)
dev.off()



#########rsea
temp=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum( (rsea_results_human_voom_go[[x]]$Comp.adjP) <0.05 ,na.rm=TRUE )) } )#0 significance
png("rsea_results_human_voom_go_sig.png")
hist(temp,breaks=100)
dev.off()

temp=apply(as.matrix(target_index),1,function(x){ print(sum( (rsea_results_human_voom_go[[x]]$Comp.adjP) <0.05 ,na.rm=TRUE )) } )#0 significance

png("rsea_results_human_voom_go_sig_knockout.png")
hist(temp,breaks=100)
dev.off()



temp=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum( (rsea_results_human_voom_wp[[x]]$Comp.adjP) <0.05 ,na.rm=TRUE )) } )#0 significance
png("rsea_results_human_voom_wp_sig.png")
hist(temp,breaks=100)
dev.off()

temp=apply(as.matrix(target_index),1,function(x){ print(sum( (rsea_results_human_voom_wp[[x]]$Comp.adjP) <0.05 ,na.rm=TRUE )) } )#0 significance

png("rsea_results_human_voom_wp_sig_knockout.png")
hist(temp,breaks=100)
dev.off()


temp=apply(as.matrix(as.numeric(as.character(args[1])):as.numeric(as.character(args[2]))),1,function(x){ print(sum( (rsea_results_human_voom_pfocr[[x]]$Comp.adjP) <0.05 ,na.rm=TRUE )) } )#0 significance
png("rsea_results_human_voom_pfocr_sig.png")
hist(temp,breaks=100)
dev.off()

temp=apply(as.matrix(target_index),1,function(x){ print(sum( (rsea_results_human_voom_pfocr[[x]]$Comp.adjP) <0.05 ,na.rm=TRUE )) } )#0 significance

png("rsea_results_human_voom_pfocr_sig_knockout.png")
hist(temp,breaks=100)
dev.off()




knockout_target_gse<-read.table("knockout_target_gse",header=FALSE)#intersect of three databases
knockout_target_gse=as.matrix(knockout_target_gse)

############################################################################################################################################gsea wp
merged_wp<-read.table("human_knockout_wp.txt",header=TRUE)

targets=which(is.na(merged_wp$WP)==FALSE)#get GSEs matched with knock genes
target_gses=merged_wp[targets,2]

target_index1=match(target_gses,(knockout_target_gse))#retrieve their index

all_gse_titles=apply(as.matrix(1:length( merged_filtered_gses)), 1, function(x) merged_filtered_gses[[x]]$GSE)

target_index2=match(target_gses,(all_gse_titles))#retrieve their index


gsea_wp_sig_gene_count=NULL
gsea_wp_aucs=NULL
for(i in target_index2){
  print(i)
  gsea_wp_sig_gene_count=c(gsea_wp_sig_gene_count,sum(gsea_results_human_voom_wp[[i]]$p.adjust <0.05,na.rm=TRUE))
  
  if( sum(is.na(gsea_results_human_voom_wp[[i]]$p.adjust)) ==(length(gsea_results_human_voom_wp[[i]]$p.adjust)) || sum(!is.na(gsea_results_human_voom_wp[[i]]$p.adjust)) <3 ){next}#pass null cases
  gse_match=which(as.character(merged_wp$GSE)==all_gse_titles[i])#match with the current gse
 
   for(j in gse_match){# there are cases where more than one same GSE exist
    print(merged_wp$WP[j])
    if(is.na(merged_wp$WP[j])){next}
    
    true_sig=merged_wp$WP[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    print(length(true_sig))
    response=apply(as.matrix(1:length(rownames(as.data.frame(gsea_results_human_voom_wp[[i]])))),1,function(x) rownames(as.data.frame(gsea_results_human_voom_wp[[i]]))[x]%in%true_sig)
    print(sum(response))
    print(table(response))
    roc_r=roc(response,as.numeric(gsea_results_human_voom_wp[[i]]$p.adjust))
    auc_r=auc(roc_r)
    gsea_wp_aucs=c(gsea_wp_aucs,auc_r[1])
  }
  
  
}
#i 452 j32  'response' must have two levels

png("gsea_wp_aucs.png")
boxplot(gsea_wp_aucs)
dev.off()

save(gsea_wp_sig_gene_count,gsea_wp_aucs,file="gsea_wp_auc.RData")

############################################################################################################################################gsea pfocr
merged_pfocr<-read.table("human_knockout_pfocr.txt",header=TRUE)

targets=which(is.na(merged_pfocr$PFOCR)==FALSE)#get GSEs matched with knock genes
target_gses=merged_pfocr[targets,2]

target_index1=match(target_gses,(knockout_target_gse))#retrieve their index

all_gse_titles=apply(as.matrix(1:length( merged_filtered_gses)), 1, function(x) merged_filtered_gses[[x]]$GSE)

target_index2=match(target_gses,(all_gse_titles))#retrieve their index

gsea_pfocr_sig_gene_count=NULL
gsea_pfocr_aucs=NULL
for(i in target_index2){
  print(i)
  gsea_pfocr_sig_gene_count=c(gsea_pfocr_sig_gene_count,sum(gsea_results_human_voom_pfocr[[i]]$p.adjust <0.05,na.rm=TRUE))
  
  if( sum(is.na(gsea_results_human_voom_pfocr[[i]]$p.adjust)) ==(length(gsea_results_human_voom_pfocr[[i]]$p.adjust)) || sum(!is.na(gsea_results_human_voom_pfocr[[i]]$p.adjust)) <3 ){next}#pass null cases
  gse_match=which(as.character(merged_pfocr$GSE)==all_gse_titles[i])#match with the current gse
  
  for(j in gse_match){# there are cases where more than one same GSE exist
    print(merged_pfocr$PFOCR[j])
    if(is.na(merged_pfocr$PFOCR[j])){next}
    
    true_sig=merged_pfocr$PFOCR[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    print(length(true_sig))
    response=apply(as.matrix(1:length(rownames(as.data.frame(gsea_results_human_voom_pfocr[[i]])))),1,function(x) rownames(as.data.frame(gsea_results_human_voom_pfocr[[i]]))[x]%in%true_sig)
    print(sum(response))
    print(table(response))
    roc_r=roc(response,as.numeric(gsea_results_human_voom_pfocr[[i]]$p.adjust))
    auc_r=auc(roc_r)
    gsea_pfocr_aucs=c(gsea_pfocr_aucs,auc_r[1])
  }
  
  
}


png("gsea_pfocr_aucs.png")
boxplot(gsea_pfocr_aucs)
dev.off()

save(gsea_pfocr_sig_gene_count,gsea_pfocr_aucs,file="gsea_pfocr_auc.RData")

#response
#FALSE  TRUE 
#32258    19 




############################################################################################################################################gsea go
merged_go<-read.table("human_knockout_go.txt",header=TRUE)

targets=which(is.na(merged_go$GO)==FALSE)#get GSEs matched with knock genes
target_gses=merged_go[targets,2]

target_index1=match(target_gses,(knockout_target_gse))#retrieve their index

all_gse_titles=apply(as.matrix(1:length( merged_filtered_gses)), 1, function(x) merged_filtered_gses[[x]]$GSE)

target_index2=match(target_gses,(all_gse_titles))#retrieve their index

gsea_go_sig_gene_count=NULL
gsea_go_aucs=NULL
for(i in target_index2){
  print(i)
  gsea_go_sig_gene_count=c(gsea_go_sig_gene_count,sum(gsea_results_human_voom_go[[i]]$p.adjust <0.05,na.rm=TRUE))
  
  
  if( sum(is.na(gsea_results_human_voom_go[[i]]$p.adjust)) ==(length(gsea_results_human_voom_go[[i]]$p.adjust)) || sum(!is.na(gsea_results_human_voom_go[[i]]$p.adjust)) <3 ){next}#pass null cases
  gse_match=which(as.character(merged_go$GSE)==all_gse_titles[i])#match with the current gse
  
  for(j in gse_match){# there are cases where more than one same GSE exist
    print(merged_go$GO[j])
    if(is.na(merged_go$GO[j])){next}
    
    true_sig=merged_go$GO[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    print(length(true_sig))
    response=apply(as.matrix(1:length(rownames(as.data.frame(gsea_results_human_voom_go[[i]])))),1,function(x) rownames(as.data.frame(gsea_results_human_voom_go[[i]]))[x]%in%true_sig)
    print(sum(response))
    print(table(response))
    roc_r=roc(response,as.numeric(gsea_results_human_voom_go[[i]]$p.adjust))
    auc_r=auc(roc_r)
    gsea_go_aucs=c(gsea_go_aucs,auc_r[1])
  }
  
  
}


png("gsea_go_aucs.png")
boxplot(gsea_go_aucs)
dev.off()

save(gsea_go_sig_gene_count,gsea_go_aucs,file="gsea_go_auc.RData")




############################################################################################################################################rsea go
merged_go<-read.table("human_knockout_go.txt",header=TRUE)

targets=which(is.na(merged_go$GO)==FALSE)#get GSEs matched with knock genes
target_gses=merged_go[targets,2]

target_index1=match(target_gses,(knockout_target_gse))#retrieve their index

all_gse_titles=apply(as.matrix(1:length( merged_filtered_gses)), 1, function(x) merged_filtered_gses[[x]]$GSE)

target_index2=match(target_gses,(all_gse_titles))#retrieve their index

rsea_go_sig_gene_count=NULL
rsea_go_aucs=NULL
for(i in target_index2){
  print(i)
  rsea_go_sig_gene_count=c(rsea_go_sig_gene_count,sum(rsea_results_human_voom_go[[i]]$Comp.adjP <0.05,na.rm=TRUE))
  
  if(sum(is.na(rsea_results_human_voom_go[[i]]$Comp.adjP))>0){
    if( sum(is.na(rsea_results_human_voom_go[[i]]$Comp.adjP)) ==(length(rsea_results_human_voom_go[[i]]$Comp.adjP)) || sum(!is.na(rsea_results_human_voom_go[[i]]$Comp.adjP)) <3 ){next}#pass null cases
  }
    gse_match=which(as.character(merged_go$GSE)==all_gse_titles[i])#match with the current gse
    

  
  for(j in gse_match){# there are cases where more than one same GSE exist
    print(merged_go$GO[j])
    if(is.na(merged_go$GO[j])){next}
    
    true_sig=merged_go$GO[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    print(length(true_sig))
    response=apply(as.matrix(1:length(rownames( as.data.frame(rsea_results_human_voom_go[[i]][c("set_id","ID","Coverage","TDP.bound","TDP.estimate","SC.adjP","Comp.adjP")])  ))),1,function(x)  unlist(as.data.frame(rsea_results_human_voom_go[[i]][c("set_id")]))[x] %in%true_sig)
    print(sum(response))
    print(table(response))
    roc_r=roc(response,as.numeric(rsea_results_human_voom_go[[i]]$Comp.adjP))
    auc_r=auc(roc_r)
    rsea_go_aucs=c(rsea_go_aucs,auc_r[1])
  }
  
  
}


png("rsea_go_aucs.png")
boxplot(rsea_go_aucs)
dev.off()

save(rsea_go_sig_gene_count,rsea_go_aucs,file="rsea_go_auc.RData")




############################################################################################################################################rsea pfocr
knockout_target_gse<-read.table("knockout_target_gse",header=FALSE)#intersect of three databases
knockout_target_gse=as.matrix(knockout_target_gse)


merged_pfocr<-read.table("human_knockout_pfocr.txt",header=TRUE)

targets=which(is.na(merged_pfocr$PFOCR)==FALSE)#get GSEs matched with knock genes
target_gses=merged_pfocr[targets,2]

target_index1=match(target_gses,(knockout_target_gse))#retrieve their index

all_gse_titles=apply(as.matrix(1:length( merged_filtered_gses)), 1, function(x) merged_filtered_gses[[x]]$GSE)

target_index2=match(target_gses,(all_gse_titles))#retrieve their index

rsea_pfocr_sig_gene_count=NULL
rsea_pfocr_aucs=NULL
for(i in target_index2){
  print(i)
  rsea_pfocr_sig_gene_count=c(rsea_pfocr_sig_gene_count,sum(rsea_results_human_voom_pfocr[[i]]$Comp.adjP <0.05,na.rm=TRUE))
  
  if(sum(is.na(rsea_results_human_voom_pfocr[[i]]$Comp.adjP))>0){
    if( sum(is.na(rsea_results_human_voom_pfocr[[i]]$Comp.adjP)) ==(length(rsea_results_human_voom_pfocr[[i]]$Comp.adjP)) || sum(!is.na(rsea_results_human_voom_pfocr[[i]]$Comp.adjP)) <3 ){next}#pass null cases
  }
  gse_match=which(as.character(merged_pfocr$GSE)==all_gse_titles[i])#match with the current gse
  
  
  
  for(j in gse_match){# there are cases where more than one same GSE exist
    print(merged_pfocr$PFOCR[j])
    if(is.na(merged_pfocr$PFOCR[j])){next}
    
    true_sig=merged_pfocr$PFOCR[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    print(length(true_sig))
    response=apply(as.matrix(1:length(rownames( as.data.frame(rsea_results_human_voom_pfocr[[i]][c("set_id","ID","Coverage","TDP.bound","TDP.estimate","SC.adjP","Comp.adjP")])  ))),1,function(x)  unlist(as.data.frame(rsea_results_human_voom_pfocr[[i]][c("set_id")]))[x] %in%true_sig)
    print(sum(response))
    print(table(response))
    roc_r=roc(response,as.numeric(rsea_results_human_voom_pfocr[[i]]$Comp.adjP))
    auc_r=auc(roc_r)
    rsea_pfocr_aucs=c(rsea_pfocr_aucs,auc_r[1])
  }
  
  
}


png("rsea_pfocr_aucs.png")
boxplot(rsea_pfocr_aucs)
dev.off()

save(rsea_pfocr_sig_gene_count,rsea_pfocr_aucs,file="rsea_pfocr_auc.RData")


############################################################################################################################################rsea wp

knockout_target_gse<-read.table("knockout_target_gse",header=FALSE)#intersect of three databases
knockout_target_gse=as.matrix(knockout_target_gse)

merged_wp<-read.table("human_knockout_wp.txt",header=TRUE)

targets=which(is.na(merged_wp$WP)==FALSE)#get GSEs matched with knock genes
target_gses=merged_wp[targets,2]

target_index1=match(target_gses,(knockout_target_gse))#retrieve their index

all_gse_titles=apply(as.matrix(1:length( merged_filtered_gses)), 1, function(x) merged_filtered_gses[[x]]$GSE)

target_index2=match(target_gses,(all_gse_titles))#retrieve their index

rsea_wp_sig_gene_count=NULL
rsea_wp_aucs=NULL
for(i in target_index2){
  print(i)
  rsea_wp_sig_gene_count=c(rsea_wp_sig_gene_count,sum(rsea_results_human_voom_wp[[i]]$Comp.adjP <0.05,na.rm=TRUE))
  
  if(sum(is.na(rsea_results_human_voom_wp[[i]]$Comp.adjP))>0){
    if( sum(is.na(rsea_results_human_voom_wp[[i]]$Comp.adjP)) ==(length(rsea_results_human_voom_wp[[i]]$Comp.adjP)) || sum(!is.na(rsea_results_human_voom_wp[[i]]$Comp.adjP)) <3 ){next}#pass null cases
  }
  gse_match=which(as.character(merged_wp$GSE)==all_gse_titles[i])#match with the current gse
  
  
  
  for(j in gse_match){# there are cases where more than one same GSE exist
    print(merged_wp$WP[j])
    if(is.na(merged_wp$WP[j])){next}
    
    true_sig=merged_wp$WP[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    print(length(true_sig))
    response=apply(as.matrix(1:length(rownames( as.data.frame(rsea_results_human_voom_wp[[i]][c("set_id","ID","Coverage","TDP.bound","TDP.estimate","SC.adjP","Comp.adjP")])  ))),1,function(x)  unlist(as.data.frame(rsea_results_human_voom_wp[[i]][c("set_id")]))[x] %in%true_sig)
    print(sum(response))
    print(table(response))
    roc_r=roc(response,as.numeric(rsea_results_human_voom_wp[[i]]$Comp.adjP))
    auc_r=auc(roc_r)
    rsea_wp_aucs=c(rsea_wp_aucs,auc_r[1])
  }
  
  
}


png("rsea_wp_aucs.png")
boxplot(rsea_wp_aucs)
dev.off()

save(rsea_wp_sig_gene_count,rsea_wp_aucs,file="rsea_wp_auc.RData")






############################################################################################################################################ora wp
sink("ora_results_human_voom_go_auc.log")
knockout_target_gse<-read.table("knockout_target_gse",header=FALSE)#intersect of three databases
knockout_target_gse=as.matrix(knockout_target_gse)


merged_wp<-read.table("human_knockout_wp.txt",header=TRUE)

targets=which(is.na(merged_wp$WP)==FALSE)#get GSEs matched with knock genes
target_gses=merged_wp[targets,2]

target_index1=match(target_gses,(knockout_target_gse))#retrieve their index

all_gse_titles=apply(as.matrix(1:length( merged_filtered_gses)), 1, function(x) merged_filtered_gses[[x]]$GSE)

target_index2=match(target_gses,(all_gse_titles))#retrieve their index


ora_wp_sig_gene_count=NULL
ora_wp_aucs=NULL
for(i in target_index2){
  print(i)
  ora_wp_sig_gene_count=c(ora_wp_sig_gene_count,sum(ora_results_human_voom_wp[[i]]$p.adjust <0.05,na.rm=TRUE))
  
  if( sum(is.na(ora_results_human_voom_wp[[i]]$p.adjust)) ==(length(ora_results_human_voom_wp[[i]]$p.adjust)) || sum(!is.na(ora_results_human_voom_wp[[i]]$p.adjust)) <3 ){next}#pass null cases
  gse_match=which(as.character(merged_wp$GSE)==all_gse_titles[i])#match with the current gse
  
  for(j in gse_match){# there are cases where more than one same GSE exist
    print(merged_wp$WP[j])
    if(is.na(merged_wp$WP[j])){next}
    
    true_sig=merged_wp$WP[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    print(length(true_sig))
    response=apply(as.matrix(1:length(rownames(as.data.frame(ora_results_human_voom_wp[[i]])))),1,function(x) rownames(as.data.frame(ora_results_human_voom_wp[[i]]))[x]%in%true_sig)
    print(sum(response))
    print(table(response))
    
    
    print("table(temp[which(response==TRUE),2])")
    print(table(temp[which(response==TRUE),2]))
    print("table(temp[which(response==FALSE),2])")
    print(table(temp[which(response==FALSE),2]))
    
    
    temp=cbind(response,as.numeric(ora_results_human_voom_wp[[i]]$p.adjust))                                                                                                                                   
    if( sum(temp[which(response==TRUE),2], na.rm=TRUE) ==0 || sum(temp[which(response==FALSE),2], na.rm=TRUE) ==0){next}
    
    roc_r=roc(response,as.numeric(ora_results_human_voom_wp[[i]]$p.adjust))
    auc_r=auc(roc_r)
    ora_wp_aucs=c(ora_wp_aucs,auc_r[1])
  }
  
  
}
sink()

png("ora_wp_aucs.png")
boxplot(ora_wp_aucs)
dev.off()

save(ora_wp_sig_gene_count,ora_wp_aucs,file="ora_wp_auc.RData")

############################################################################################################################################ora pfocr
sink("ora_results_human_voom_pfocr_auc.log")
knockout_target_gse<-read.table("knockout_target_gse",header=FALSE)#intersect of three databases
knockout_target_gse=as.matrix(knockout_target_gse)

merged_pfocr<-read.table("human_knockout_pfocr.txt",header=TRUE)

targets=which(is.na(merged_pfocr$PFOCR)==FALSE)#get GSEs matched with knock genes
target_gses=merged_pfocr[targets,2]

target_index1=match(target_gses,(knockout_target_gse))#retrieve their index

all_gse_titles=apply(as.matrix(1:length( merged_filtered_gses)), 1, function(x) merged_filtered_gses[[x]]$GSE)

target_index2=match(target_gses,(all_gse_titles))#retrieve their index

ora_pfocr_sig_gene_count=NULL
ora_pfocr_aucs=NULL
for(i in target_index2){
  print(i)
  ora_pfocr_sig_gene_count=c(ora_pfocr_sig_gene_count,sum(ora_results_human_voom_pfocr[[i]]$p.adjust <0.05,na.rm=TRUE))
  
  if( sum(is.na(ora_results_human_voom_pfocr[[i]]$p.adjust)) ==(length(ora_results_human_voom_pfocr[[i]]$p.adjust)) || sum(!is.na(ora_results_human_voom_pfocr[[i]]$p.adjust)) <3 ){next}#pass null cases
  gse_match=which(as.character(merged_pfocr$GSE)==all_gse_titles[i])#match with the current gse
  
  for(j in gse_match){# there are cases where more than one same GSE exist
    print(merged_pfocr$PFOCR[j])
    if(is.na(merged_pfocr$PFOCR[j])){next}
    
    true_sig=merged_pfocr$PFOCR[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    print(length(true_sig))
    response=apply(as.matrix(1:length(rownames(as.data.frame(ora_results_human_voom_pfocr[[i]])))),1,function(x) rownames(as.data.frame(ora_results_human_voom_pfocr[[i]]))[x]%in%true_sig)
    print(sum(response))
    print(table(response))
    
    print("table(temp[which(response==TRUE),2])")
    print(table(temp[which(response==TRUE),2]))
    print("table(temp[which(response==FALSE),2])")
    print(table(temp[which(response==FALSE),2]))
    
    temp=cbind(response,as.numeric(ora_results_human_voom_pfocr[[i]]$p.adjust))                                                                                                                                   
    if( sum(temp[which(response==TRUE),2], na.rm=TRUE) ==0 ||  sum(temp[which(response==FALSE),2], na.rm=TRUE)==0){next}
    
    
    roc_r=roc(response,as.numeric(ora_results_human_voom_pfocr[[i]]$p.adjust))
    auc_r=auc(roc_r)
    ora_pfocr_aucs=c(ora_pfocr_aucs,auc_r[1])
  }
  
  
}
sink()

png("ora_pfocr_aucs.png")
boxplot(ora_pfocr_aucs)
dev.off()

save(ora_pfocr_sig_gene_count,ora_pfocr_aucs,file="ora_pfocr_auc.RData")

#response
#FALSE  TRUE 
#32258    19 




############################################################################################################################################ora go
sink("ora_results_human_voom_go_auc.log")
merged_go<-read.table("human_knockout_go.txt",header=TRUE)

targets=which(is.na(merged_go$GO)==FALSE)#get GSEs matched with knock genes
target_gses=merged_go[targets,2]

target_index1=match(target_gses,(knockout_target_gse))#retrieve their index

all_gse_titles=apply(as.matrix(1:length( merged_filtered_gses)), 1, function(x) merged_filtered_gses[[x]]$GSE)

target_index2=match(target_gses,(all_gse_titles))#retrieve their index

ora_go_sig_gene_count=NULL
ora_go_aucs=NULL
for(i in target_index2){
  print(i)
  ora_go_sig_gene_count=c(ora_go_sig_gene_count,sum(ora_results_human_voom_go[[i]]$p.adjust <0.05,na.rm=TRUE))
  
  
  if( sum(is.na(ora_results_human_voom_go[[i]]$p.adjust)) ==(length(ora_results_human_voom_go[[i]]$p.adjust)) || sum(!is.na(ora_results_human_voom_go[[i]]$p.adjust)) <3 ){next}#pass null cases
  gse_match=which(as.character(merged_go$GSE)==all_gse_titles[i])#match with the current gse
  
  for(j in gse_match){# there are cases where more than one same GSE exist
    print(merged_go$GO[j])
    if(is.na(merged_go$GO[j])){next}
    
    true_sig=merged_go$GO[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    print(length(true_sig))
    response=apply(as.matrix(1:length(rownames(as.data.frame(ora_results_human_voom_go[[i]])))),1,function(x) rownames(as.data.frame(ora_results_human_voom_go[[i]]))[x]%in%true_sig)
    print(sum(response))
    print(table(response))
    
    temp=cbind(response,as.numeric(ora_results_human_voom_go[[i]]$p.adjust))
    
    print("table(temp[which(response==TRUE),2])")
    print(table(temp[which(response==TRUE),2]))
    print("table(temp[which(response==FALSE),2])")
    print(table(temp[which(response==FALSE),2]))
    
    if( sum(temp[which(response==TRUE),2], na.rm=TRUE) ==0 ){next}
    if(sum(temp[which(response==FALSE),2], na.rm=TRUE) ==0){next} 
    
    
    roc_r=roc(response,as.numeric(ora_results_human_voom_go[[i]]$p.adjust))
    auc_r=auc(roc_r)
    ora_go_aucs=c(ora_go_aucs,auc_r[1])
  }
  
  
}
sink()

png("ora_go_aucs.png")
boxplot(ora_go_aucs)
dev.off()

save(ora_go_sig_gene_count,ora_go_aucs,file="ora_go_auc.RData")



############################################################################################################################################check correlation between gsea and rsea results
x=1


pvalR=apply(as.matrix(1:length(merged_filtered_gses)),1,function(x){
  
  rsea_temp=rsea_results_human_voom_go[[x]]$Comp.adjP
  rownames(rsea_temp)=(rsea_results_human_voom_go[[x]])$set_id
  
  gsea_temp=gsea_results_human_voom_go[[x]]$pvalue
  rownames(gsea_temp)=rownames(as.data.frame(gsea_results_human_voom_go[[x]]))
  
  index=match(rownames(rsea_temp),rownames(gsea_temp))
  gsea_temp=gsea_temp[index]
  print(head(gsea_temp))
  print(head(rsea_temp))
  
  if(dim(as.data.frame(gsea_results_human_voom_go[[x]]))[1]==1){return(NULL)}
  
  corR=cor.test(gsea_temp,rsea_temp,use="pairwise.complete.obs")
  corR$p.value
  
})

jpeg("gsea_rsea_pval_cor_pval.jpeg")
hist(-log10(unlist(pvalR)),breaks=100)
dev.off()