# Identify disease ontology terms tagging WikiPathways and PFOCR content,
# then intersect with Harmonizome datasets

library(tidyverse)

###############################################################################

## Identify disease terms used for WikiPathways
library(rWikiPathways)

# get all human pathways tagged with any disease ontology term
disease.pathways <- rWikiPathways::getPathwayIdsByParentOntologyTerm("DOID:4")
hs.pathways <- rWikiPathways::listPathwayIds("Homo sapiens")
hs.disease.pathways <- intersect(disease.pathways, hs.pathways)

# examine disease terms among these pathways
all.wp.terms <- lapply(hs.disease.pathways, function(p){
  ont <- getOntologyTerms(p)
  # transform into desired df
  df <- as.data.frame(unlist(ont), stringsAsFactors = F)
  names(df) <- 'terms'
  df2 <- df %>%
    mutate(key = rep(names(ont[[1]]), n() / 3), 
         id = cumsum(key == names(ont[[1]])[[1]])) %>% 
    spread(key, terms)
  
  df2 %>%
    filter(ontology == 'Disease') %>%
    mutate(wpid = p, .before = 1)


})
# combine into single df of WikiPathways-associated disease terms
all.wp.terms <- bind_rows(all.wp.terms)

# 56 unique disease terms annotating WikiPathways
unique(all.wp.terms$name)

####################

## Alternative source for disease terms for WikiPathways,
##  using same enrichment method used for PFOCR

gmt.wp.overlaps <- readRDS("gmt-wp-overlaps.RDS")
gmt.wp.overlaps.2 <- filter(gmt.wp.overlaps, wp.overlap.cnt >= 2)
all.wp.terms.2 <- gmt.wp.overlaps.2 %>%
  select(c(term,wpid)) %>%
  mutate(ontology = "Disease") %>%
  droplevels()
unique(all.wp.terms.2$term)

###############################################################################

## Identify disease terms used for PFOCR

gmt.pfocr.overlaps <- readRDS("gmt-pfocr-overlaps.RDS")
gmt.pfocr.overlaps.2 <- filter(gmt.pfocr.overlaps, pf.overlap.cnt >= 2)
all.pf.terms.2 <- gmt.pfocr.overlaps.2 %>%
  select(c(term,figid)) %>%
  mutate(ontology = "Disease") %>%
  droplevels()
unique(all.pf.terms.2$term)

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
wp.pf.terms.1 <- lapply(unique(all.wp.terms.2$name), function(p){
  hits <- agrep(p, unique(all.pf.terms.2$name), max = 1, ignore.case = T)
  if(length(hits) > 0){
  data.frame(wp.term = p, pf.term.index = hits) %>%
    mutate(pf.term.name =  unique(all.pf.terms.2$name)[pf.term.index]) %>%
    select(c(1,3))
  }
})
wp.pf.terms.1 <- bind_rows(wp.pf.terms.1)

wp.pf.terms.selected <- wp.pf.terms.1 %>%
  filter(!wp.term %in% c("cancer","SIDS", "ALS", "disease"))

## WP and Data
wp.data.terms.1 <- lapply(unique(all.wp.terms.2$term), function(p){
  hits <- agrep(p, unique(ds$dataset_annot), max = 1, ignore.case = T)
  if(length(hits) > 0){
    data.frame(wp.term = p, data.term.index = hits) %>%
      mutate(data.term.name =  unique(ds$dataset_annot)[data.term.index]) %>%
      select(c(1,3))
  }
})
wp.data.terms.1 <- bind_rows(wp.data.terms.1)

# 23 unique wp.terms overlap with 25 unique data.terms
wp.data.terms.selected <- wp.data.terms.1 %>%
  filter(!wp.term %in% c("cancer","Cancer", "SIDS", "ALS", "disease"))

## PFOCR and Data
pf.data.terms.1 <- lapply(unique(all.pf.terms.2$term), function(p){
  hits <- agrep(p, unique(ds$dataset_annot), max = 1, ignore.case = T)
  if(length(hits) > 0){
    data.frame(pf.term = p, data.term.index = hits) %>%
      mutate(data.term.name =  unique(ds$dataset_annot)[data.term.index]) %>%
      select(c(1,3))
  }
})
pf.data.terms.1 <- bind_rows(pf.data.terms.1)

# 26 unique pf.terms overlap with 29 unique data.terms
pf.data.terms.selected <- pf.data.terms.1 %>%
  filter(!pf.term %in% c("cancer","Cancer", "Myopathy")) %>%
  filter(!data.term.name %in% c("Nephroblastoma"))
