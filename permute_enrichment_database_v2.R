# scriptDir=/wynton/group/gladstone/biocore/projects/pfocr_pathway_enrichment_evaluation/PFOCRInPathwayAnalyses/
# containerDir=/wynton/group/gladstone/biocore/containers
# dataDir=/wynton/group/gladstone/biocore/projects/PFOCR/
# export SINGULARITY_BINDPATH="$containerDir,$scriptDir,$dataDir"
# singularity exec $containerDir/pathway_enrichment_pathway_databases_evaluation_latest.sif R

##################################required libraries
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
require(pROC)
args <- commandArgs(trailingOnly=TRUE)

GSE_index <- as.integer(args[1])
print(GSE_index)

min_set_size <- 3
max_set_size <- 500 
dataDir <- "/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/"
database_lists1 <- load("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/databases_pfocr_3intersect_v2.RData")#pfocr_3sets
database_lists2 <- load("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/databases.RData")#has wp, pfocr, go
database_lists3<-load("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/pfocr_miny_April2021.RData")

database_lists <- c(database_lists1, database_lists2,database_lists3)
database_lists <- unname(unlist(sapply(database_lists, grep, pattern="_list$", value = T, perl = T)))
for (db in database_lists) {
  eval(call("<-", as.name(db),  Filter(Negate(is.null), lapply(get(db), function(x){
    if(length(x) < min_set_size | length(x) > max_set_size)
      NULL
    else
      x
  }))
  ))
}



##################################required data sets
pvalue_results_human_voom=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/pvalue_results_human_voom.rds")
logFC_results_human_voom=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/logFC_results_human_voom.rds")
merged_filtered_gses=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/merged_filtered_gses.rds")

gene_entrez=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/gene_entrez.rds")

##################################required functions
source("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/run_rSEA5.r")
source("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/run_gsea4.r")
source("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/run_ORA5.r")


set.seed(1234)

Nperm <- 100

PermuteDatabase <- function(p, GSE_index, path_list, path_annotation,pvalue_results_human_voom, gene_entrez) {
  

  #print(p)
  #path_list=wp_list;path_annotation=wp_annotation
  temp_list <- path_list
  temp_annotation <- path_annotation
  temp_annotation$gene <- path_annotation$gene[sample(1:nrow(path_annotation), nrow(path_annotation), replace = FALSE)]
  for(index in 1:length(temp_list)) {
    temp_list[[index]] <- temp_annotation %>%
      subset(set_id == names(temp_list)[index]) %>%
      pull(gene) %>%
      as.character()
  }
  
  data_m=data.frame(rownames(pvalue_results_human_voom),pvalue_results_human_voom[,GSE_index])
  colnames(data_m)=c("Gene","pvalue")
  merged=merge(data_m,gene_entrez,by="Gene",all = "L",all.x=TRUE,all.y=FALSE)
  merged <- merged[!is.na(merged$pvalue) & !is.na(merged$ENTREZID),]
  
  sea_result <- SEA(merged$pvalue, merged$ENTREZID, pathlist = temp_list)

  return(list(temp_list,sea_result))
  
}

run_permueDatabase<-function(p,GSE_index){

rNsig_wp<-PermuteDatabase(1, GSE_index, wp_list, wp_annotation,pvalue_results_human_voom, gene_entrez)
rNsig_go <- PermuteDatabase(1, GSE_index, go_list, go_annotation,pvalue_results_human_voom, gene_entrez)
rNsig_pfocr <- PermuteDatabase(1, GSE_index, pfocr_list, pfocr_annotation,pvalue_results_human_voom, gene_entrez)
rNsig_pfocr_3sets <- PermuteDatabase(1, GSE_index, pfocr_3sets_list, pfocr_annotation_3sets,pvalue_results_human_voom, gene_entrez)
rNsig_pfocr_miny <- PermuteDatabase(1, GSE_index, pfocr_miny_list, pfocr_miny_annotation,pvalue_results_human_voom, gene_entrez)

if(!dir.exists(paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE))){
  dir.create(file.path(paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE)), showWarnings = FALSE);
}


apply(as.matrix(GSE_index:GSE_index),1,function(x){ 
  dir.create(file.path(paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE,"/rSEA")), showWarnings = FALSE);
  dir.create(file.path(paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE,"/rSEA/WP")), showWarnings = FALSE);
  dir.create(file.path(paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE,"/rSEA/PFOCR")), showWarnings = FALSE);
  dir.create(file.path(paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE,"/rSEA/PFOCR_3sets")), showWarnings = FALSE);
  dir.create(file.path(paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE,"/rSEA/PFOCR_miny")), showWarnings = FALSE);
  dir.create(file.path(paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE,"/rSEA/GO")), showWarnings = FALSE);

}
)

write.table(rNsig_wp[[2]],paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE,"/rSEA/WP/result",p,".txt"),
            col.names=colnames(rNsig_wp[[2]]),row.names=rownames(rNsig_wp[[2]]),quote=FALSE)
write.table(rNsig_go[[2]],paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE,"/rSEA/GO/result",p,".txt"),
            col.names=colnames(rNsig_go[[2]]),row.names=rownames(rNsig_go[[2]]),quote=FALSE)
write.table(rNsig_pfocr[[2]],paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE,"/rSEA/PFOCR/result",p,".txt"),
            col.names=colnames(rNsig_pfocr[[2]]),row.names=rownames(rNsig_pfocr[[2]]),quote=FALSE)
write.table(rNsig_pfocr_3sets[[2]],paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE,"/rSEA/PFOCR_3sets/result",p,".txt"),
            col.names=colnames(rNsig_pfocr_3sets[[2]]),row.names=rownames(rNsig_pfocr_3sets[[2]]),quote=FALSE)
write.table(rNsig_pfocr_miny[[2]],paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE,"/rSEA/PFOCR_miny/result",p,".txt"),
            col.names=colnames(rNsig_pfocr_miny[[2]]),row.names=rownames(rNsig_pfocr_miny[[2]]),quote=FALSE)

saveRDS(rNsig_wp[[1]],file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE,"/rSEA/WP/result",p,".rds"))
saveRDS(rNsig_go[[1]],file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE,"/rSEA/GO/result",p,".rds"))
saveRDS(rNsig_pfocr[[1]],file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE,"/rSEA/PFOCR/result",p,".rds"))
saveRDS(rNsig_pfocr_3sets[[1]],file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE,"/rSEA/PFOCR_3sets/result",p,".rds"))
saveRDS(rNsig_pfocr_miny[[1]],file=paste0("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/GSE_db_permuted/", merged_filtered_gses[[GSE_index]]$GSE,"/rSEA/PFOCR_miny/result",p,".rds"))

}

sapply(args[2]:args[3], run_permueDatabase, GSE_index)

#rNsig_wp <- t(sapply(1:Nperm, run_permueDatabase, GSE_index, wp_list, wp_annotation,pvalue_results_human_voom, gene_entrez))
#rNsig_go <- t(sapply(1:Nperm, run_permueDatabase, GSE_index, go_list, go_annotation,pvalue_results_human_voom, gene_entrez))
#rNsig_pfocr <- t(sapply(1:Nperm, run_permueDatabase, GSE_index, pfocr_list, pfocr_annotation,pvalue_results_human_voom, gene_entrez))
#rNsig_pfocr_3sets <- t(sapply(1:Nperm, run_permueDatabase, GSE_index, pfocr_3sets_list, pfocr_annotation_3sets,pvalue_results_human_voom, gene_entrez))
#rNsig_pfocr_miny <- t(sapply(1:Nperm, run_permueDatabase, GSE_index, pfocr_miny_list, pfocr_miny_annotation,pvalue_results_human_voom, gene_entrez))

