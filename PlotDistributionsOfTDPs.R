rm(list = ls())
setwd("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses")

load(paste0(getwd(), "/PFOCRInPathwayAnalyses/databases.RData"))
ocr_length <- sapply(pfocr_list, length)
go_length <- sapply(go_list, length)
wp_length <- sapply(wp_list, length)
GSEid <- list.files(path = "./PFOCRInPathwayAnalyses/GSE/")
minSig <- 5
minSize <- 10
theta <- 0.9
thetaL <- as.character(100*theta)
ChooseGSE <- 1
getTDPQuantiles <- function(ChooseGSE, GSEid, ocr_length, go_length, wp_length, minSig, minSize, theta) {
  print(ChooseGSE)
  q_ocr <- NA
  q_wp <- NA
  q_go <- NA
  ocr_result <- read.table(paste0(getwd(), "/PFOCRInPathwayAnalyses/GSE/",GSEid[ChooseGSE],"/rSEA/PFOCR_3sets/result.txt"), header = T)
  tempIndices <- match(ocr_result$set_id, names(ocr_length))
  temp_ocr_length <- ocr_length[tempIndices]
  Nocr <- sum(!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 & temp_ocr_length >= minSize, na.rm = T)
  if(Nocr > minSig) {
    sig_ocr_result <- data.frame(ocr_result[!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 & temp_ocr_length >= minSize, ], database=rep("ocr", Nocr))
    q_ocr <- quantile(sig_ocr_result$TDP.estimate, theta, na.rm = T)
  }
  
  wp_result <- read.table(paste0(getwd(), "/PFOCRInPathwayAnalyses/GSE/",GSEid[ChooseGSE],"/rSEA/WP/result.txt"), header = T)
  tempIndices <- match(wp_result$set_id, names(wp_length))
  temp_wp_length <- wp_length[tempIndices]
  
  Nwp <- sum(!is.na(wp_result$Comp.adjP) & wp_result$Comp.adjP < 0.05 & temp_wp_length >= minSize, na.rm = T)
  if(Nwp > minSig) {
    sig_wp_result <- data.frame(wp_result[!is.na(wp_result$Comp.adjP) & wp_result$Comp.adjP < 0.05 & temp_wp_length >= minSize, ], database=rep("wp", Nwp))
    q_wp <- quantile(sig_wp_result$TDP.estimate, theta, na.rm = T)
  }
  
  go_result <- read.table(paste0(getwd(), "/PFOCRInPathwayAnalyses/GSE/",GSEid[ChooseGSE],"/rSEA/GO/result.txt"), header = T)
  tempIndices <- match(go_result$set_id, names(go_length))
  temp_go_length <- go_length[tempIndices]
  
  Ngo <- sum(!is.na(go_result$Comp.adjP) & go_result$Comp.adjP < 0.05 & temp_go_length >= minSize & temp_go_length < 500, na.rm = T)
  if(Ngo > minSig) {
    sig_go_result <- data.frame(go_result[!is.na(go_result$Comp.adjP) & go_result$Comp.adjP < 0.05 & temp_go_length >= minSize & temp_go_length < 500, ], database=rep("go", Ngo))
    q_go <- quantile(sig_go_result$TDP.estimate, theta, na.rm = T)
  }
  
  # PlotData <- rbind(sig_ocr_result, sig_wp_result)
  # PlotData <- data.frame(rbind(PlotData, sig_go_result))
  # ggplot(PlotData, aes(x=TDP.estimate, color=database)) + stat_ecdf() + geom_vline(xintercept = q_go,lty=2)
  # ggplot(PlotData, aes(y=TDP.estimate, x=database)) + geom_boxplot()
  
  return(c(q_ocr, q_wp, q_go))
}

q_results <- t(sapply(1:length(GSEid), getTDPQuantiles, GSEid, ocr_length, go_length, wp_length, minSig, minSize, theta))
q_results <- data.frame(q_results)
colnames(q_results) <- c("median_tdp_ocr_3sets", "median_tdp_wp", "median_tdp_go")

require(ggplot2)

pdf(paste0("pathway_databases_w_pfocr_3sets_comparison_theta_level_", thetaL,"_percent.pdf"))
for(i in 1:(ncol(q_results) - 1)) {
  for(j in (i+1):ncol(q_results)) {
    databases <- names(q_results)
    x <- databases[i]
    y <- databases[j]
    tempRes <- t.test(q_results[,i], q_results[,j], paired = T)
    print(ggplot(q_results, aes(x=q_results[,i], y=q_results[,j])) 
          + geom_point() 
          + geom_abline(slope = 1, intercept = 0, lty=2, col="red") 
          + xlab(x) + ylab(y)
          + ggtitle(paste0("Mean difference = ", round(1e4*tempRes$estimate)/1e4, ";", "pvalue = ", round(1e4*tempRes$p.value)/1e4)))
  }
}
dev.off()

get_pfocr_TDPQuantiles <- function(ChooseGSE, GSEid, ocr_length, go_length, wp_length, minSig, minSize, theta) {
  print(ChooseGSE)
  q_ocr_orig <- NA
  q_ocr_noalias <- NA
  q_ocr_nobe0 <- NA
  q_ocr_nobe2 <- NA
  q_ocr_nobe3 <- NA
  q_ocr_noprev <- NA
 
  ocr_result <- read.table(paste0(getwd(), "/GSE/",GSEid[ChooseGSE],"/rSEA/PFOCR/result.txt"), header = T)
  tempIndices <- match(ocr_result$set_id, names(ocr_length))
  temp_ocr_length <- ocr_length[tempIndices]
  Nocr <- sum(!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , na.rm = T)
  if(Nocr > minSig) {
    sig_ocr_result <- data.frame(ocr_result[!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , ], database=rep("ocr_noalias", Nocr))
    q_ocr_orig <- quantile(sig_ocr_result$TDP.estimate, theta, na.rm = T)
  }
  
  ocr_result <- read.table(paste0(getwd(), "/GSE/",GSEid[ChooseGSE],"/rSEA/PFOCR_noalias/result.txt"), header = T)
  tempIndices <- match(ocr_result$set_id, names(ocr_length))
  temp_ocr_length <- ocr_length[tempIndices]
  Nocr <- sum(!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , na.rm = T)
  if(Nocr > minSig) {
    sig_ocr_result <- data.frame(ocr_result[!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , ], database=rep("ocr_noalias", Nocr))
    q_ocr_noalias <- quantile(sig_ocr_result$TDP.estimate, theta, na.rm = T)
  }
  
  ocr_result <- read.table(paste0(getwd(), "/GSE/",GSEid[ChooseGSE],"/rSEA/PFOCR_nobe0/result.txt"), header = T)
  tempIndices <- match(ocr_result$set_id, names(ocr_length))
  temp_ocr_length <- ocr_length[tempIndices]
  Nocr <- sum(!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , na.rm = T)
  if(Nocr > minSig) {
    sig_ocr_result <- data.frame(ocr_result[!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , ], database=rep("ocr_nobe0", Nocr))
    q_ocr_nobe0 <- quantile(sig_ocr_result$TDP.estimate, theta, na.rm = T)
  }
  
  
  ocr_result <- read.table(paste0(getwd(), "/GSE/",GSEid[ChooseGSE],"/rSEA/PFOCR_nobe2/result.txt"), header = T)
  tempIndices <- match(ocr_result$set_id, names(ocr_length))
  temp_ocr_length <- ocr_length[tempIndices]
  Nocr <- sum(!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , na.rm = T)
  if(Nocr > minSig) {
    sig_ocr_result <- data.frame(ocr_result[!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , ], database=rep("ocr_nobe2", Nocr))
    q_ocr_nobe2 <- quantile(sig_ocr_result$TDP.estimate, theta, na.rm = T)
  }
  
  ocr_result <- read.table(paste0(getwd(), "/GSE/",GSEid[ChooseGSE],"/rSEA/PFOCR_nobe3/result.txt"), header = T)
  tempIndices <- match(ocr_result$set_id, names(ocr_length))
  temp_ocr_length <- ocr_length[tempIndices]
  Nocr <- sum(!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , na.rm = T)
  if(Nocr > minSig) {
    sig_ocr_result <- data.frame(ocr_result[!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , ], database=rep("ocr_nobe3", Nocr))
    q_ocr_nobe3 <- quantile(sig_ocr_result$TDP.estimate, theta, na.rm = T)
  }
  
  ocr_result <- read.table(paste0(getwd(), "/GSE/",GSEid[ChooseGSE],"/rSEA/PFOCR_noprev/result.txt"), header = T)
  tempIndices <- match(ocr_result$set_id, names(ocr_length))
  temp_ocr_length <- ocr_length[tempIndices]
  Nocr <- sum(!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , na.rm = T)
  if(Nocr > minSig) {
    sig_ocr_result <- data.frame(ocr_result[!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , ], database=rep("ocr_noprev", Nocr))
    q_ocr_noprev <- quantile(sig_ocr_result$TDP.estimate, theta, na.rm = T)
  }
  
  # PlotData <- rbind(sig_ocr_result, sig_wp_result)
  # PlotData <- data.frame(rbind(PlotData, sig_go_result))
  # ggplot(PlotData, aes(x=TDP.estimate, color=database)) + stat_ecdf() + geom_vline(xintercept = q_go,lty=2)
  # ggplot(PlotData, aes(y=TDP.estimate, x=database)) + geom_boxplot()
  
  return(c(q_ocr_orig, q_ocr_noalias, q_ocr_nobe0, q_ocr_nobe2, q_ocr_nobe3, q_ocr_noprev))
}

ocr_q_results <- t(sapply(1:length(GSEid), get_pfocr_TDPQuantiles, GSEid, ocr_length, go_length, wp_length, minSig, minSize, theta))
ocr_q_results <- data.frame(ocr_q_results)
colnames(ocr_q_results) <- c("median_tdp_ocr_orig","median_tdp_ocr_noalias", "median_tdp_ocr_nobe0", "median_tdp_ocr_nobe2", "median_tdp_ocr_nobe3", "median_tdp_ocr_noprev")

pdf(paste0("ocr_databases_comparison_theta_level_", thetaL,"_percent.pdf"))
for(i in 1:(ncol(ocr_q_results) - 1)) {
  for(j in (i+1):ncol(ocr_q_results)) {
    databases <- names(ocr_q_results)
    x <- databases[i]
    y <- databases[j]
    tempRes <- t.test(ocr_q_results[,i], ocr_q_results[,j], paired = T)
    print(ggplot(ocr_q_results, aes(x=ocr_q_results[,i], y=ocr_q_results[,j])) 
          + geom_point() 
          + geom_abline(slope = 1, intercept = 0, lty=2, col="red") 
          + xlab(x) + ylab(y)
          + ggtitle(paste0("Mean difference = ", round(1e4*tempRes$estimate)/1e4, ";", "pvalue = ", round(1e4*tempRes$p.value)/1e4)))
  }
}
dev.off()
