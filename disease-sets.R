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

###############################################################################

## Identify disease terms used for PFOCR
# retreive disease annotations 
enr.annots <- read.table("./enriched_annots.tsv", header=T, sep="\t", stringsAsFactors = F, quote = c("\""))

# transform into useful df
enr.annots.sub <- enr.annots %>%
  dplyr::select(c(1,2))%>%
  dplyr::filter(!is.na(.[2]))  
all.pf.terms <- enr.annots.sub %>%
  mutate(name = strsplit(!!as.name(names(enr.annots.sub)[2]), " | ", fixed = T),
         ontology = "Disease") %>% 
  unnest(name) %>%
  select(c(1,3,4))

# 151 unique disease terms annotating WikiPathways
unique(all.pf.terms$name)

###############################################################################

# Intersection of disease terms for WP and PFOCR