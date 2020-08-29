# Identify disease ontology terms tagging WikiPathways and PFOCR content,
# then intersect with Harmonizome datasets

library(tidyverse)

###############################################################################

# ## Identify disease terms used for WikiPathways
# library(rWikiPathways)
# 
# # get all human pathways tagged with any disease ontology term
# disease.pathways <- rWikiPathways::getPathwayIdsByParentOntologyTerm("DOID:4")
# hs.pathways <- rWikiPathways::listPathwayIds("Homo sapiens")
# hs.disease.pathways <- intersect(disease.pathways, hs.pathways)
# 
# # examine disease terms among these pathways
# all.wp.terms <- lapply(hs.disease.pathways, function(p){
#   ont <- getOntologyTerms(p)
#   # transform into desired df
#   df <- as.data.frame(unlist(ont), stringsAsFactors = F)
#   names(df) <- 'terms'
#   df2 <- df %>%
#     mutate(key = rep(names(ont[[1]]), n() / 3), 
#          id = cumsum(key == names(ont[[1]])[[1]])) %>% 
#     spread(key, terms)
#   
#   df2 %>%
#     filter(ontology == 'Disease') %>%
#     mutate(wpid = p, .before = 1)
# 
# 
# })
# # combine into single df of WikiPathways-associated disease terms
# all.wp.terms <- bind_rows(all.wp.terms)
# 
# # 56 unique disease terms annotating WikiPathways
# unique(all.wp.terms$name)

####################

## Alternative source for disease terms for WikiPathways,
##  using same enrichment method used for PFOCR

gmt.wp.overlaps <- readRDS("gmt-wp-overlaps.RDS")

gmt.wp.overlaps.5 <- filter(gmt.wp.overlaps, wp.overlap.cnt >= 5)
max.wp.terms.5 <- gmt.wp.overlaps.5 %>%
  group_by(term) %>%
  top_n(1,wp.overlap.cnt) %>%
  select(c(term, wpid, wp.overlap.cnt)) %>%
  mutate(ontology = "Disease") %>%
  droplevels()
unique(max.wp.terms.5$term)

## check against pathway names
# wp.check <- lapply(1:nrow(max.wp.terms.5), function(i){
#   p <- max.wp.terms.5$wpid[[i]]
#   n <- getPathwayInfo(p)$name
#   d <- max.wp.terms.5$term[[i]]
#   data_frame(disease=d, name=n, wpid=p, cnt=max.wp.terms.5$wp.overlap.cnt[[i]])
# })
# wp.check <- as.data.frame(bind_rows(wp.check))

cur.wp.rows <- c(1:3,6:11,13:17,19,20,22,25,29,31,32,36,37,42,43,44,45,46,47,49,50,
                 51,53,54,55,58,61,64,65,66,67,69,70,71,72,73,75,76,77,80,81,
                 82,85,88,89,90,91,93)

cur.wp.terms <- max.wp.terms.5[cur.wp.rows,]

saveRDS(cur.wp.terms, "cur-wp-terms.RDS")
###############################################################################

## Identify disease terms used for PFOCR

gmt.pfocr.overlaps <- readRDS("gmt-pfocr-overlaps.RDS")
gmt.pfocr.overlaps.5 <- filter(gmt.pfocr.overlaps, pf.overlap.cnt >= 5)
max.pf.terms.5 <- gmt.pfocr.overlaps.5 %>%
  group_by(term) %>%
  top_n(1,pf.overlap.cnt) %>%
  select(c(term, figid, pf.overlap.cnt)) %>%
  mutate(ontology = "Disease") %>%
  droplevels()
unique(max.pf.terms.5$term)
max.pf.terms.5 <- as.data.frame(max.pf.terms.5)

## Check and curate annotations
# pfocr.figs <- readRDS("~/Dropbox (Gladstone)/PFOCR_25Years/pfocr_figures.rds")
# pf.check <- lapply(1:nrow(max.pf.terms.5), function(i){
#   f <- max.pf.terms.5$figid[[i]]
#   ft <- pfocr.figs %>%
#     filter(figid == f) %>%
#     select(figtitle)
#   pt <- pfocr.figs %>%
#     filter(figid == f) %>%
#     select(papertitle)
#   d <- max.pf.terms.5$term[[i]]
#   data_frame(disease=d, paper=pt, fig=ft, figid=f, cnt=max.pf.terms.5$pf.overlap.cnt[[i]])
# })
# pf.check <- as.data.frame(bind_rows(pf.check))

cur.pf.rows <- c(1,2,3,6,8,10,11,12,14,15,16,19,20,22,24,26,28,29,30,121,122,123,124,125,126,
                 131,132,135,136,137,138,139,140,141,143,144,145,151,152,153,154,155,157,158,
                 159,160,162,167,179,183,184,186,187,200,201,202,203,205,206,207,210,
                 211,212,213,214,215,216,217,220,221,222,223,244,246,247,255,256,257)

cur.pf.terms <- max.pf.terms.5[cur.pf.rows,]

saveRDS(cur.pf.terms, "cur-pf-terms.RDS")

###############################################################################

# Identify disease annotations for GEO datasets

## NOTE: First clone the gladstone-institutes/bioinformatics-datasets repo in 
##       the same parent dir as this repo, such that the following line works.

setwd("../bioinformatics-datasets/")
source("query-gene-sets.R")

ds <- listHarmonizomeSets()

#163 unique disease terms annotating harmonizome datasets
unique(ds$dataset_annot)

###############################################################################

# Intersection of disease terms

## WP and PFOCR
wp.pf.terms.1 <- lapply(unique(cur.wp.terms$term), function(p){
  hits <- agrep(p, unique(cur.pf.terms$term), max = 1, ignore.case = T)
  if(length(hits) > 0){
  data.frame(wp.term = p, pf.term.index = hits) %>%
    mutate(pf.term.name =  unique(cur.pf.terms$term)[pf.term.index]) %>%
    select(c(1,3))
  }
})
wp.pf.terms.1 <- bind_rows(wp.pf.terms.1)

wp.pf.terms.selected <- wp.pf.terms.1 %>%
  filter(!wp.term %in% c("Cancer"))

## WP and Data
wp.data.terms.1 <- lapply(unique(cur.wp.terms$term), function(p){
  hits <- agrep(p, unique(ds$dataset_annot), max = 1, ignore.case = T)
  if(length(hits) > 0){
    data.frame(wp.term = p, data.term.index = hits) %>%
      mutate(data.term.name =  unique(ds$dataset_annot)[data.term.index]) %>%
      select(c(1,3))
  }
})
wp.data.terms.1 <- bind_rows(wp.data.terms.1)

# 10 unique wp.terms overlap with 13 unique data.terms
wp.data.terms.selected <- wp.data.terms.1 %>%
  filter(!wp.term %in% c("Cancer"))

## PFOCR and Data
pf.data.terms.1 <- lapply(unique(cur.pf.terms$term), function(p){
  hits <- agrep(p, unique(ds$dataset_annot), max = 1, ignore.case = T)
  if(length(hits) > 0){
    data.frame(pf.term = p, data.term.index = hits) %>%
      mutate(data.term.name =  unique(ds$dataset_annot)[data.term.index]) %>%
      select(c(1,3))
  }
})
pf.data.terms.1 <- bind_rows(pf.data.terms.1)

# 16 unique pf.terms overlap with 19 unique data.terms
pf.data.terms.selected <- pf.data.terms.1


## Intersection of data.terms: 11 unique disease terms
intersect(pf.data.terms.selected$data.term.name, wp.data.terms.selected$data.term.name)

# [1] "Alzheimer's Disease"                             "Arrhythmogenic Right Ventricular Cardiomyopathy" "Breast Cancer"                                  
# [4] "Cardiomyopathy"                                  "Cardiomyopathy, Dilated"                         "Viral cardiomyopathy"                           
# [7] "JRA - Juvenile rheumatoid arthritis"             "Leber congenital amaurosis"                      "Diabetic Neuropathy"                            
# [10] "Peripheral motor neuropathy"                     "Parkinson's Disease"               

## Union of data.terms: 21 unique disease terms
all.overlapping.terms <- union(pf.data.terms.selected$data.term.name, wp.data.terms.selected$data.term.name)

# [1] "Alzheimer's Disease"                             "Arrhythmogenic Right Ventricular Cardiomyopathy" "Breast Cancer"                                  
# [4] "Cardiomyopathy"                                  "Cardiomyopathy, Dilated"                         "Viral cardiomyopathy"                           
# [7] "Crohn's disease"                                 "Type 1 diabetes mellitus"                        "Type 2 diabetes mellitus"                       
# [10] "JRA - Juvenile rheumatoid arthritis"             "Leber congenital amaurosis"                      "In situ melanoma of skin"                       
# [13] "Diabetic Neuropathy"                             "Peripheral motor neuropathy"                     "Obesity"                                        
# [16] "Parkinson's Disease"                             "Polycystic Ovary Syndrome"                       "Retinitis Pigmentosa"                           
# [19] "SCID - Severe combined immunodeficiency"         "Epilepsy"                                        "Schizophrenia"                  

###############################################################################

# Select GEO datasets based on disease terms overlapping with WP and PFOCR
target.ds <- as.data.frame(ds %>%
  filter(dataset_annot %in% all.overlapping.terms) 
  )

saveRDS(target.ds, "target-datasets.RDS")
write.csv(target.ds, "target-datasets.csv", sep=",",row.names = F)
