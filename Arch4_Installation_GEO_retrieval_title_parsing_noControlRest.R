# R script to download selected samples
# Copy code and run on a local machine to initiate download

#https://docs.google.com/document/d/1jm1sJtHgJkSV7hEjMC_VK7l7gzshMQzr4T7NdgY9rGE/edit#heading=h.skn982pchl79
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4930834/

# Check for dependencies and install if missing
packages <- c("rhdf5","pkgconfig")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  BiocManager::install("rhdf5")
}
#make ~/.R/Makevars  and add the following (remove # in the beginning of each line)
# LDFLAGS= -L/usr/local/clang4/lib
# FLIBS=-L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin16/6.3.0 -L/usr/local/gfortran/lib -lgfortran -lquadmath -lm
# CC=/usr/local/clang4/bin/clang
# CXX=/usr/local/gfortran/bin/g++
#   CXX1X=/usr/local/gfortran/bin/g++
#   CXX98=/usr/local/gfortran/bin/g++
#   CXX11=/usr/local/gfortran/bin/g++
#   CXX14=/usr/local/gfortran/bin/g++
#   CXX17=/usr/local/gfortran/bin/g++
#   
#   LLVM_LOC = /usr/local/opt/llvm
# CC=/usr/local/gfortran/bin/gcc -fopenmp
# CXX=/usr/local/gfortran/bin/g++ -fopenmp
# # -O3 should be faster than -O2 (default) level optimisation ..
# CFLAGS=-g -O3 -Wall -pedantic -std=gnu99 -mtune=native -pipe
# CXXFLAGS=-g -O3 -Wall -pedantic -std=c++11 -mtune=native -pipe
# LDFLAGS=-L/usr/local/opt/gettext/lib -L$(LLVM_LOC)/lib -Wl,-rpath,$(LLVM_LOC)/lib
# CPPFLAGS=-I/usr/local/opt/gettext/include -I$(LLVM_LOC)/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include

#add this path to your mac:  /usr/local/gfortran/bin
#If you have data.table error: ‘datatable.so’ not found
#https://github.com/Rdatatable/data.table/wiki/Installation#openmp-enabled-compiler-for-mac

# sample size
# platform: not array
# year
# organism


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("GEOquery", quietly = TRUE))
  install.packages("GEOquery")

library("rhdf5")
library("tools")
library("GEOquery")

setwd("~/Dropbox (Gladstone)/Pathway-Project/human_matrix_download.h5")
destination_file = "~/Dropbox (Gladstone)/Pathway-Project/human_matrix_download.h5"

gse = h5read(destination_file, "meta/Sample_series_id")
head(gse)
length(gse)
238522
length(unique(gse))
7909

unique_gse=unique(gse)


###gse issue example: "GSE27452Xx-xXGSE32303"
##removing
library(stringr)
str_count(unique_gse[3],"GSE")
GSE_counts=apply(unique_gse,1,function(x) str_count(x,"GSE"))
max(GSE_counts)
3
length(which(GSE_counts>1))
1607
unique_gse=unique_gse[-which(GSE_counts>1)]

dieases=c("Renal cell carcinoma","Alzheimer’s disease","Alzheimer","Thyroid cancer","Dilated cardiomyopathy",
          "Colorectal cancer","Prostate cancer","Huntington’s disease","Huntington","Acute Myeloid Leukemia",
          "Pancreatic cancer","Non-small cell lung cancer","Glioma","Parkinson’s disease","Parkinson",
          "Type II diabetes mellitus","Type 2 diabetes mellitus","Chronic Myeloid Leukemia","Endometrial cancer")



# NA       Renal cell carcinoma        Alzheimer’s disease                  Alzheimer             Thyroid cancer 
# 5795                         22                         10                         13                         13 
# Dilated cardiomyopathy          Colorectal cancer            Prostate cancer       Huntington’s disease                 Huntington 
# 4                         68                        133                          5                          9 
# Acute Myeloid Leukemia          Pancreatic cancer Non-small cell lung cancer                     Glioma        Parkinson’s disease 
# 70                         37                         27                         61                         12 
# Parkinson  Type II diabetes mellitus   Type 2 diabetes mellitus   Chronic Myeloid Leukemia         Endometrial cancer 
# 8                          1                          6                          4                          4 

#distinguish sample names



########## match additional terms
terms=read.table("87disease",header=FALSE)
terms=as.character(unlist((terms)))
terms=c(terms, dieases,"knock-out", "knock out", "knockout","siRNA","knockdown", "knock down", "knock-down", "crispr-i", "variant","pathway")
terms=unique(terms)
match_gse_terms<-function(gse_id){
  print(gse_id)
  
  
  
  tryCatch({
    
    data <- getGEO(gse_id)
    
  }, error = function(err) {
    

    return(NULL)
    
  }) # END tryCatch
  
  
  title=data[[1]]@experimentData@title
  abstract=data[[1]]@experimentData@abstract
  
  matched_disease_title=apply(as.matrix(terms),1,function(x) grep(x,title,ignore.case=TRUE) )
  matched_disease_abstract=apply(as.matrix(terms),1,function(x) grep(x,abstract,ignore.case=TRUE) )
  
  if(length(matched_disease_title)>0 | length(matched_disease_abstract)>0){
    if(length(matched_disease_title)>0){
      which_disease=which(matched_disease_title>0)
    }else{
      which_disease=which(matched_disease_abstract>0)
    }
    return(disease[which_disease]) 
  }else{
    return(NULL)
  }
  if(which_disease>0){print(disease[which_disease])}
  
}


match_gse_terms2<-function(description,abstract){

 
  matched_disease_title=apply(as.matrix(terms),1,function(x) grep(x,description,ignore.case=TRUE) )
  matched_disease_abstract=apply(as.matrix(terms),1,function(x) grep(x,abstract,ignore.case=TRUE) )
  
  if(length(matched_disease_title)>0 | length(matched_disease_abstract)>0){
    if(length(matched_disease_title)>0){
      which_disease=which(matched_disease_title>0)
    }else{
      which_disease=which(matched_disease_abstract>0)
    }
    return(terms[which_disease]) 
  }else{
    return(NULL)
  }
  if(which_disease>0){print(disease[which_disease])}
  
}


matched_terms=apply(unique_gse,1,match_gse_terms)

matched_terms_unlisted=as.matrix(unlist(lapply(matched_terms, `[[`, 1)))

length(which(matched_terms_unlisted>0))
1480

tb=table(matched_terms_unlisted)
names(tb)=c("NA",terms)
tb

#knockout    siRNA crispr-i=>knockdown  variant  pathway 
#   105      216        1      158     1000 



#code that retrieves samples
#https://github.com/ctlab/phantasus/blob/master/R/loadGEO.R#L708


#################sort titles



group_samples<-function(titles){

library(stringdist)
library(igraph)
library(raster)

controls=NULL
controls=c(controls,grep("control",titles,ignore.case = TRUE))
controls=c(controls,grep("wild",titles,ignore.case = TRUE))

controls=c(controls,which(regexpr("[^a-z]wt[^a-z]",titles,ignore.case = TRUE) >0))

titles[controls]

if(length(controls)==0){
  print("No control samples")
  return(list())
}

###case

distm=as.matrix(stringdist::stringdistmatrix(titles[-controls],method="lv"))
colnames(distm)=titles[-controls]
rownames(distm)=titles[-controls]
distm


# s.a <- strsplit(titles[1], "")[[1]]
# s.b <- strsplit(titles[2], "")[[1]]
# paste(s.a[s.a != s.b], collapse = "")


distm[which(distm>1)]=0

g  <- graph.adjacency(distm) 
#plot(g)
clu <- components(g)
gr=igraph::groups(clu)

gr=gr[which(lengths(gr)>1)]

if(length(gr)==0){
  print("No case samples")
  return(list())
}


sample_groups=list()

for(i in 1:length(gr)){
  #print(gr[[i]])
  sample_groups<-c(sample_groups,list(gr[[i]]))
  names(sample_groups)[i]=paste("case",i,sep="")
  last_case=i
}



###control

distm=as.matrix(stringdist::stringdistmatrix(titles[controls],method="lv"))
colnames(distm)=titles[controls]
rownames(distm)=titles[controls]
distm


distm[which(distm>1)]=0

g  <- graph.adjacency(distm) 
#plot(g)
clu <- components(g)
gr=igraph::groups(clu)

gr=gr[which(lengths(gr)>1)]

if(length(gr)==0){
  print("No control samples")
  return(list())
}

for(i in 1:length(gr)){
  #print(gr[[i]])
  sample_groups<-c(sample_groups,list(gr[[i]]))
  names(sample_groups)[last_case+i]=paste("control",i,sep="")
}

return(sample_groups)
}

group_samples2<-function(titles){
  
  library(stringdist)
  library(igraph)
  library(raster)
  
  controls=NULL
  controls=c(controls,grep("control",titles,ignore.case = TRUE))
  controls=c(controls,grep("wild",titles,ignore.case = TRUE))
  
  controls=c(controls,which(regexpr("[^a-z]wt[^a-z]",titles,ignore.case = TRUE) >0))
  
  titles[controls]
  
  if(length(controls)==0){
    print("No control samples")
    return(list())
  }
  
  if(length(titles)-length(controls)==0){
    print("No case samples")
    return(list())
  }
  ###case
  
  titles_cases=titles[-controls]
  rep_keys=get_rep_keys(titles_cases)
  
  sample_groups=list()
  
  last_case=1
  
  if(length(rep_keys)>0){
    for(i in 1:length(rep_keys)){
      #print(gr[[i]])
      matched=grep(rep_keys[i],titles_cases)
      sample_groups<-c(sample_groups,list(titles_cases[matched]))
      titles_cases=titles_cases[-matched]
      names(sample_groups)[i]=paste("case",i,sep="")
      last_case=i
    }
  }

  if(length(titles_cases)>1){
    
    distm=as.matrix(stringdist::stringdistmatrix(titles_cases,method="lv"))
    colnames(distm)=titles_cases
    rownames(distm)=titles_cases
    distm
    
    
    # s.a <- strsplit(titles[1], "")[[1]]
    # s.b <- strsplit(titles[2], "")[[1]]
    # paste(s.a[s.a != s.b], collapse = "")
    
    
    distm[which(distm>1)]=0
    
    g  <- graph.adjacency(distm) 
    #plot(g)
    clu <- components(g)
    gr=igraph::groups(clu)
    
    gr=gr[which(lengths(gr)>1)]
    
    if(length(gr)==0&last_case==1){
      print("No case samples")
      return(list())
    }
      
      for(i in (last_case+1):length(gr)){
        #print(gr[[i]])
        sample_groups<-c(sample_groups,list(gr[[i]]))
        names(sample_groups)[i]=paste("case",i,sep="")
        
      }
    
    
  }
 
  
  
  
  

  
  
  
  ###control
  
  titles_controls=titles[controls]
  
  rep_keys=get_rep_keys(titles_controls)
  
  sample_groups=list()
  
  last_control=1
  
  if(length(rep_keys)>0){
    for(i in 1:length(rep_keys)){
      #print(gr[[i]])
      matched=grep(rep_keys[i],titles_controls)
      sample_groups<-c(sample_groups,list(titles_controls[matched]))
      titles_controls=titles_controls[-matched]
      names(sample_groups)[i]=paste("control",i,sep="")
      last_control=i
    }
  }
  
  if(length(titles_controls)>1){
    
    distm=as.matrix(stringdist::stringdistmatrix(titles_controls,method="lv"))
    colnames(distm)=titles_controls
    rownames(distm)=titles_controls
    distm
    
    
    # s.a <- strsplit(titles[1], "")[[1]]
    # s.b <- strsplit(titles[2], "")[[1]]
    # paste(s.a[s.a != s.b], collapse = "")
    
    
    distm[which(distm>1)]=0
    
    g  <- graph.adjacency(distm) 
    #plot(g)
    clu <- components(g)
    gr=igraph::groups(clu)
    
    gr=gr[which(lengths(gr)>1)]
    
    if(length(gr)==0&last_control==1){
      print("No control samples")
      return(list())
    }
    
    for(i in (last_control+1):length(gr)){
      #print(gr[[i]])
      sample_groups<-c(sample_groups,list(gr[[i]]))
      names(sample_groups)[i]=paste("control",i,sep="")
    
    }
    
    
  }
  
  
  return(sample_groups)
}

get_rep_keys<-function(query){
  
inferCondition <- sub("(?![-+}\\]\\)])(\\W|_)*((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(r)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*\\d+(\\W|_)*", "", query, ignore.case = TRUE, perl = TRUE)
  
 inferCondition2 <- sub("(?![-+}\\]\\)])(\\W|_)*((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*[^a-z]((rep.*)|(r)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*\\d+(\\_)*", "", query, ignore.case = TRUE, perl = TRUE)

  if(length(inferCondition)<inferCondition2){inferCondition=inferCondition2}
  
  tb=as.data.frame(inferCondition) %>% group_by_all() %>% summarise(Unique_Elements = n()) 
  uniques=tb[which(tb[,2]>1),1]
  
  return(as.character(as.matrix(uniques)))
}

get_rep_keys2<-function(query){
  
   #inferCondition <- sub("(?![-+}\\]\\)])(\\W|_)*((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(r)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*\\d+(\\W|_)*", "", query, ignore.case = TRUE, perl = TRUE)
   m=regexpr("(?![-+}\\]\\)])(\\W|_)*((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(r)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*\\d+(\\W|_)*", query, ignore.case = TRUE, perl = TRUE)
   inferCondition=apply(as.matrix(1:length(query)),1,function(x) substr(query[x],1,m[x]-1))
   
  #inferCondition2 <- sub("(?![-+}\\]\\)])(\\W|_)*((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*[^a-z]((rep.*)|(r)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*\\d+(\\_)*", "", query, ignore.case = TRUE, perl = TRUE)
  m=regexpr("(?![-+}\\]\\)])(\\W|_)*((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*[^a-z]((rep.*)|(r)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*\\d+(\\_)*", query, ignore.case = TRUE, perl = TRUE)
  inferCondition2=apply(as.matrix(1:length(query)),1,function(x) substr(query[x],1,m[x]-1))
  
  if(length(inferCondition)<length(inferCondition2)){inferCondition=inferCondition2}
  
  tb=as.data.frame(inferCondition) %>% group_by_all() %>% summarise(Unique_Elements = n()) 
  uniques=tb[which(tb[,2]>1),1]
  
  return(as.character(as.matrix(uniques)))
}

get_rep_keys3<-function(query){
  
  #inferCondition <- sub("(?![-+}\\]\\)])(\\W|_)*((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(r)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*\\d+(\\W|_)*", "", query, ignore.case = TRUE, perl = TRUE)
  m=regexpr("(?![-+}\\]\\)])(\\W|_)*((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*\\d+(\\W|_)*", query, ignore.case = TRUE, perl = TRUE)
  inferCondition=apply(as.matrix(1:length(query)),1,function(x) substr(query[x],1,m[x]-1))
  nulls=which(inferCondition=="")
  if(length(nulls)>0){inferCondition=inferCondition[-nulls]}
  
  #inferCondition2 <- sub("(?![-+}\\]\\)])(\\W|_)*((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*[^a-z]((rep.*)|(r)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*\\d+(\\_)*", "", query, ignore.case = TRUE, perl = TRUE)
  m=regexpr("(?![-+}\\]\\)])(\\W|_)*((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*[^a-zA-Z]((rep.*)|(r)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*\\d+(\\_)*", query, ignore.case = TRUE, perl = TRUE)
  inferCondition2=apply(as.matrix(1:length(query)),1,function(x) substr(query[x],1,m[x]-1))
  nulls=which(inferCondition2=="")
  if(length(nulls)>0){inferCondition2=inferCondition2[-nulls]}
  
  if(length(inferCondition)==0 & length(inferCondition2) ==0){return(NULL)}
  
  if(length(inferCondition)<length(inferCondition2)){inferCondition=inferCondition2}
  
  tb=as.data.frame(inferCondition) %>% group_by_all() %>% summarise(Unique_Elements = n()) 
  uniques=tb[which(tb[,2]>1),1]
  
  return(as.character(as.matrix(uniques)))
}

group_samples3<-function(gse,titles){
  
  library(stringdist)
  library(igraph)
  library(raster)
  
  controls=NULL
  controls=c(controls,grep("control",titles,ignore.case = TRUE))
  controls=c(controls,grep("wild",titles,ignore.case = TRUE))
  controls=c(controls,grep("untreated",titles,ignore.case = TRUE))
  controls=c(controls,which(regexpr("[^a-z]wt[^a-z]",titles,ignore.case = TRUE) >0))
  pre_post_flag=0
  if( length(which(regexpr("[^a-z]pre[^a-z]",titles,ignore.case = TRUE) >0))>0 & length(which(regexpr("[^a-z]post[^a-z]",titles,ignore.case = TRUE) >0))>0 ){
    pre_post_flag=1
  }
  
  if(pre_post_flag){
    controls=c(controls,which(regexpr("[^a-z]pre[^a-z]",titles,ignore.case = TRUE) >0))
  }
  titles[controls]
  
  if(length(controls)==0){
    print("No control samples")
    return(list())
  }else if(length(controls)==1){
    print("Only one control sample")
    return(list())
  }
  
  if(length(titles)-length(controls)==0){
    print("No case samples")
    return(list())
  }else if(length(titles)-length(controls)==1){
    print("Only one case sample")
    return(list())
  }
  ###case
  
  titles_cases=titles[-controls]
  rep_keys=get_rep_keys3(titles_cases)
  
  sample_groups=list()
  
  last_case=0
  
  if(length(rep_keys)>0){
   
    if(rep_keys!=""){
      rep_keys=rep_keys[order(nchar(rep_keys), rep_keys,decreasing=TRUE)]
      if(length(which(nchar(rep_keys)==0))){
        rep_keys=rep_keys[-which(nchar(rep_keys)==0)]
      }
      for(i in 1:length(rep_keys)){
        #print(gr[[i]])
        rep_keys[i]=gsub("\\+","\\\\+",rep_keys[i])
        rep_keys[i]=gsub("\\(","\\\\(",rep_keys[i])
        rep_keys[i]=gsub("\\?","\\\\?",rep_keys[i])
        rep_keys[i]=gsub("\\*","\\\\*",rep_keys[i])
        rep_keys[i]=gsub("\\{","\\\\{",rep_keys[i])
        rep_keys[i]=gsub("\\[","\\\\[",rep_keys[i])
        rep_keys[i]=gsub("\\|","\\\\|",rep_keys[i])
        
        matched=grep(rep_keys[i],titles_cases)
        sample_groups<-c(sample_groups,list(titles_cases[matched]))
        titles_cases=titles_cases[-matched]
        names(sample_groups)[i]=paste("case",i,sep="")
        last_case=i
      }
    }

  }
  
  if(length(titles_cases)>1){
    
    distm=as.matrix(stringdist::stringdistmatrix(titles_cases,method="lv"))
    colnames(distm)=titles_cases
    rownames(distm)=titles_cases
    distm
    
    
    # s.a <- strsplit(titles[1], "")[[1]]
    # s.b <- strsplit(titles[2], "")[[1]]
    # paste(s.a[s.a != s.b], collapse = "")
    
    
    distm[which(distm>1)]=0
    
    g  <- graph.adjacency(distm) 
    #plot(g)
    clu <- components(g)
    gr=igraph::groups(clu)
    
    gr=gr[which(lengths(gr)>1)]
    
    if(length(gr)==0&last_case==0){
      print("No case samples")
      return(list())
    }else if(length(gr)>0){
      for(i in 1:length(gr)){
        #print(gr[[i]])
        sample_groups<-c(sample_groups,list(gr[[i]]))
        names(sample_groups)[length(sample_groups)]=paste("case",last_case+i,sep="")
        
      }
    }
    
   
    
    
  }
  
  
  
  
  
  
  
  
  
  ###control
  
  titles_controls=titles[controls]
  
  rep_keys=get_rep_keys3(titles_controls)
  

  
  last_control=0
  
  if(length(rep_keys)>0){
    
    if(rep_keys!=""){
      rep_keys=rep_keys[order(nchar(rep_keys), rep_keys,decreasing=TRUE)]
      if(length(which(nchar(rep_keys)==0))){
        rep_keys=rep_keys[-which(nchar(rep_keys)==0)]
      }
      
      for(i in 1:length(rep_keys)){
        #print(gr[[i]])
        rep_keys[i]=gsub("\\+","\\\\+",rep_keys[i])
        rep_keys[i]=gsub("\\(","\\\\(",rep_keys[i])
        rep_keys[i]=gsub("\\?","\\\\?",rep_keys[i])
        rep_keys[i]=gsub("\\*","\\\\*",rep_keys[i])
        rep_keys[i]=gsub("\\{","\\\\{",rep_keys[i])
        rep_keys[i]=gsub("\\[","\\\\[",rep_keys[i])
        rep_keys[i]=gsub("\\|","\\\\|",rep_keys[i])
        matched=grep(rep_keys[i],titles_controls)
        sample_groups<-c(sample_groups,list(titles_controls[matched]))
        titles_controls=titles_controls[-matched]
        names(sample_groups)[length(sample_groups)]=paste("control",i,sep="")
        last_control=i
      }
    }

  }
  
  if(length(titles_controls)>1){
    
    distm=as.matrix(stringdist::stringdistmatrix(titles_controls,method="lv"))
    colnames(distm)=titles_controls
    rownames(distm)=titles_controls
    distm
    
    
    # s.a <- strsplit(titles[1], "")[[1]]
    # s.b <- strsplit(titles[2], "")[[1]]
    # paste(s.a[s.a != s.b], collapse = "")
    
    
    distm[which(distm>1)]=0
    
    g  <- graph.adjacency(distm) 
    #plot(g)
    clu <- components(g)
    gr=igraph::groups(clu)
    
    gr=gr[which(lengths(gr)>1)]
    
    if(length(gr)==0&last_control==0){
      print("No control samples")
      return(list())
    }else if(length(gr)>0){
      for(i in 1:length(gr)){
        #print(gr[[i]])
        sample_groups<-c(sample_groups,list(gr[[i]]))
        names(sample_groups)[length(sample_groups)]=paste("control",last_control+i,sep="")
        
      }
    }
    
   
    
    
  }
  sample_groups<-c(sample_groups,list(GSE=gse))
  sample_groups<-c(sample_groups,list(titles=titles))
  sample_groups<-c(sample_groups,list(total_sample_count=length(titles)))
  sample_groups<-c(sample_groups,list(parsed_sample_count=length(unlist(sample_groups))-length(titles)-2))
  
  return(sample_groups)
}

group_samples4<-function(gse,titles){
  
  library(stringdist)
  library(igraph)
  library(raster)
  
  
  
  tryCatch({
    
    controls=NULL
    controls=c(controls,grep("control",titles,ignore.case = TRUE))
    controls=c(controls,grep("wild",titles,ignore.case = TRUE))
    controls=c(controls,grep("untreated",titles,ignore.case = TRUE))
    controls=c(controls,which(regexpr("[^a-z]wt[^a-z]",titles,ignore.case = TRUE) >0))
    pre_post_flag=0
    if( length(which(regexpr("[^a-z]pre[^a-z]",titles,ignore.case = TRUE) >0))>0 & length(which(regexpr("[^a-z]post[^a-z]",titles,ignore.case = TRUE) >0))>0 ){
      pre_post_flag=1
    }
    
    if(pre_post_flag){
      controls=c(controls,which(regexpr("[^a-z]pre[^a-z]",titles,ignore.case = TRUE) >0))
    }
    titles[controls]
    
    if(length(controls)==0){
      print("No control samples")
      return(list())
    }else if(length(controls)==1){
      print("Only one control sample")
      return(list())
    }
    
    if(length(titles)-length(controls)==0){
      print("No case samples")
      return(list())
    }else if(length(titles)-length(controls)==1){
      print("Only one case sample")
      return(list())
    }
    ###case
    
    titles_cases=titles[-controls]
    rep_keys=get_rep_keys3(titles_cases)
    
    sample_groups=list()
    
    last_case=0
    
    if(length(rep_keys)>0){
      
      if(rep_keys!=""){
        rep_keys=rep_keys[order(nchar(rep_keys), rep_keys,decreasing=TRUE)]
        if(length(which(nchar(rep_keys)==0))){
          rep_keys=rep_keys[-which(nchar(rep_keys)==0)]
        }
        for(i in 1:length(rep_keys)){
          #print(gr[[i]])
          rep_keys[i]=gsub("\\+","\\\\+",rep_keys[i])
          rep_keys[i]=gsub("\\(","\\\\(",rep_keys[i])
          rep_keys[i]=gsub("\\?","\\\\?",rep_keys[i])
          rep_keys[i]=gsub("\\*","\\\\*",rep_keys[i])
          rep_keys[i]=gsub("\\{","\\\\{",rep_keys[i])
          rep_keys[i]=gsub("\\[","\\\\[",rep_keys[i])
          rep_keys[i]=gsub("\\|","\\\\|",rep_keys[i])
          
          matched=grep(rep_keys[i],titles_cases)
          sample_groups<-c(sample_groups,list(titles_cases[matched]))
          titles_cases=titles_cases[-matched]
          names(sample_groups)[i]=paste("case",i,sep="")
          last_case=i
        }
      }
      
    }
    
    if(length(titles_cases)>1){
      
      distm=as.matrix(stringdist::stringdistmatrix(titles_cases,method="lv"))
      colnames(distm)=titles_cases
      rownames(distm)=titles_cases
      distm
      
      
      # s.a <- strsplit(titles[1], "")[[1]]
      # s.b <- strsplit(titles[2], "")[[1]]
      # paste(s.a[s.a != s.b], collapse = "")
      
      
      distm[which(distm>1)]=0
      
      g  <- graph.adjacency(distm) 
      #plot(g)
      clu <- components(g)
      gr=igraph::groups(clu)
      
      gr=gr[which(lengths(gr)>1)]
      
      if(length(gr)==0&last_case==0){
        print("No case samples")
        return(list())
      }else if(length(gr)>0){
        for(i in 1:length(gr)){
          #print(gr[[i]])
          sample_groups<-c(sample_groups,list(gr[[i]]))
          names(sample_groups)[length(sample_groups)]=paste("case",last_case+i,sep="")
          
        }
      }
      
      
      
      
    }
    
    
    
    
    
    
    
    
    
    ###control
    
    titles_controls=titles[controls]
    
    rep_keys=get_rep_keys3(titles_controls)
    
    
    
    last_control=0
    
    if(length(rep_keys)>0){
      
      if(rep_keys!=""){
        rep_keys=rep_keys[order(nchar(rep_keys), rep_keys,decreasing=TRUE)]
        if(length(which(nchar(rep_keys)==0))){
          rep_keys=rep_keys[-which(nchar(rep_keys)==0)]
        }
        
        for(i in 1:length(rep_keys)){
          #print(gr[[i]])
          rep_keys[i]=gsub("\\+","\\\\+",rep_keys[i])
          rep_keys[i]=gsub("\\(","\\\\(",rep_keys[i])
          rep_keys[i]=gsub("\\?","\\\\?",rep_keys[i])
          rep_keys[i]=gsub("\\*","\\\\*",rep_keys[i])
          rep_keys[i]=gsub("\\{","\\\\{",rep_keys[i])
          rep_keys[i]=gsub("\\[","\\\\[",rep_keys[i])
          rep_keys[i]=gsub("\\|","\\\\|",rep_keys[i])
          matched=grep(rep_keys[i],titles_controls)
          sample_groups<-c(sample_groups,list(titles_controls[matched]))
          titles_controls=titles_controls[-matched]
          names(sample_groups)[length(sample_groups)]=paste("control",i,sep="")
          last_control=i
        }
      }
      
    }
    
    if(length(titles_controls)>1){
      
      distm=as.matrix(stringdist::stringdistmatrix(titles_controls,method="lv"))
      colnames(distm)=titles_controls
      rownames(distm)=titles_controls
      distm
      
      
      # s.a <- strsplit(titles[1], "")[[1]]
      # s.b <- strsplit(titles[2], "")[[1]]
      # paste(s.a[s.a != s.b], collapse = "")
      
      
      distm[which(distm>1)]=0
      
      g  <- graph.adjacency(distm) 
      #plot(g)
      clu <- components(g)
      gr=igraph::groups(clu)
      
      gr=gr[which(lengths(gr)>1)]
      
      if(length(gr)==0&last_control==0){
        print("No control samples")
        return(list())
      }else if(length(gr)>0){
        for(i in 1:length(gr)){
          #print(gr[[i]])
          sample_groups<-c(sample_groups,list(gr[[i]]))
          names(sample_groups)[length(sample_groups)]=paste("control",last_control+i,sep="")
          
        }
      }
      
      
      
      
    }
    sample_groups<-c(sample_groups,list(GSE=gse))
    sample_groups<-c(sample_groups,list(titles=titles))
    sample_groups<-c(sample_groups,list(total_sample_count=length(titles)))
    sample_groups<-c(sample_groups,list(parsed_sample_count=length(unlist(sample_groups))-length(titles)-2))
    
    return(sample_groups)
    
  }, error = function(err) {
    
    
    
    return(list())
    
  }) # END tryCatch
  
  
  
}

group_samples5<-function(gse,titles){
  
  library(stringdist)
  library(igraph)
  library(raster)
  
  
  
  tryCatch({
    
    controls=NULL
    controls=c(controls,grep("control",titles,ignore.case = TRUE))
    controls=c(controls,grep("wild",titles,ignore.case = TRUE))
    controls=c(controls,grep("untreated",titles,ignore.case = TRUE))
    controls=c(controls,which(regexpr("[^a-z]wt[^a-z]",titles,ignore.case = TRUE) >0))
    controls=c(controls,which(regexpr("wt[^a-z]",titles,ignore.case = TRUE) >0))
    controls=c(controls,which(regexpr("[^a-z]wt$",titles,ignore.case = TRUE) >0))
    controls=c(controls,which(regexpr("[^a-z]healthy",titles,ignore.case = TRUE) >0))
    controls=c(controls,which(regexpr("^healthy",titles,ignore.case = TRUE) >0))
    controls=c(controls,which(regexpr("uninfected",titles,ignore.case = TRUE) >0))
    controls=c(controls,which(regexpr("^standard",titles,ignore.case = TRUE) >0))
    controls=c(controls,which(regexpr("[^a-z]standard",titles,ignore.case = TRUE) >0))
    controls=c(controls,which(regexpr("ctrl",titles,ignore.case = TRUE) >0))
    controls=c(controls,which(regexpr("[^a-z]normal",titles,ignore.case = TRUE) >0))
    controls=c(controls,which(regexpr("^normal",titles,ignore.case = TRUE) >0))
    
    
    
    
    pre_post_flag=0
    if( length(which(regexpr("[^a-z]pre[^a-z]",titles,ignore.case = TRUE) >0))>0 & length(which(regexpr("[^a-z]post[^a-z]",titles,ignore.case = TRUE) >0))>0 ){
      pre_post_flag=1
    }
    
    if(pre_post_flag){
      controls=c(controls,which(regexpr("[^a-z]pre[^a-z]",titles,ignore.case = TRUE) >0))
    }
    titles[controls]
    
    if(length(controls)==0){
      print("No control samples")
      return(list())
    }else if(length(controls)==1){
      print("Only one control sample")
      return(list())
    }
    
    if(length(titles)-length(controls)==0){
      print("No case samples")
      return(list())
    }else if(length(titles)-length(controls)==1){
      print("Only one case sample")
      return(list())
    }
    ###case
    
    titles_cases=titles[-controls]
    rep_keys=get_rep_keys3(titles_cases)
    
    sample_groups=list()
    
    last_case=0
    
    if(length(rep_keys)>0){
      
      if(rep_keys!=""){
        rep_keys=rep_keys[order(nchar(rep_keys), rep_keys,decreasing=TRUE)]
        if(length(which(nchar(rep_keys)==0))){
          rep_keys=rep_keys[-which(nchar(rep_keys)==0)]
        }
        for(i in 1:length(rep_keys)){
          #print(gr[[i]])
          rep_keys[i]=gsub("\\+","\\\\+",rep_keys[i])
          rep_keys[i]=gsub("\\(","\\\\(",rep_keys[i])
          rep_keys[i]=gsub("\\?","\\\\?",rep_keys[i])
          rep_keys[i]=gsub("\\*","\\\\*",rep_keys[i])
          rep_keys[i]=gsub("\\{","\\\\{",rep_keys[i])
          rep_keys[i]=gsub("\\[","\\\\[",rep_keys[i])
          rep_keys[i]=gsub("\\|","\\\\|",rep_keys[i])
          
          matched=grep(rep_keys[i],titles_cases)
          sample_groups<-c(sample_groups,list(titles_cases[matched]))
          titles_cases=titles_cases[-matched]
          names(sample_groups)[i]=paste("case",i,sep="")
          last_case=i
        }
      }
      
    }
    
    if(length(titles_cases)>1){
      
      distm=as.matrix(stringdist::stringdistmatrix(titles_cases,method="lv"))
      colnames(distm)=titles_cases
      rownames(distm)=titles_cases
      distm
      
      
      # s.a <- strsplit(titles[1], "")[[1]]
      # s.b <- strsplit(titles[2], "")[[1]]
      # paste(s.a[s.a != s.b], collapse = "")
      
      
      distm[which(distm>1)]=0
      
      g  <- graph.adjacency(distm) 
      #plot(g)
      clu <- components(g)
      gr=igraph::groups(clu)
      
      gr=gr[which(lengths(gr)>1)]
      
      if(length(gr)==0&last_case==0){
        print("No case samples")
        return(list())
      }else if(length(gr)>0){
        for(i in 1:length(gr)){
          #print(gr[[i]])
          sample_groups<-c(sample_groups,list(gr[[i]]))
          names(sample_groups)[length(sample_groups)]=paste("case",last_case+i,sep="")
          
        }
      }
      
      
      
      
    }
    
    
    
    
    
    
    
    
    
    ###control
    
    titles_controls=titles[controls]
    
    rep_keys=get_rep_keys3(titles_controls)
    
    
    
    last_control=0
    
    if(length(rep_keys)>0){
      
      if(rep_keys!=""){
        rep_keys=rep_keys[order(nchar(rep_keys), rep_keys,decreasing=TRUE)]
        if(length(which(nchar(rep_keys)==0))){
          rep_keys=rep_keys[-which(nchar(rep_keys)==0)]
        }
        
        for(i in 1:length(rep_keys)){
          #print(gr[[i]])
          rep_keys[i]=gsub("\\+","\\\\+",rep_keys[i])
          rep_keys[i]=gsub("\\(","\\\\(",rep_keys[i])
          rep_keys[i]=gsub("\\?","\\\\?",rep_keys[i])
          rep_keys[i]=gsub("\\*","\\\\*",rep_keys[i])
          rep_keys[i]=gsub("\\{","\\\\{",rep_keys[i])
          rep_keys[i]=gsub("\\[","\\\\[",rep_keys[i])
          rep_keys[i]=gsub("\\|","\\\\|",rep_keys[i])
          matched=grep(rep_keys[i],titles_controls)
          sample_groups<-c(sample_groups,list(titles_controls[matched]))
          titles_controls=titles_controls[-matched]
          names(sample_groups)[length(sample_groups)]=paste("control",i,sep="")
          last_control=i
        }
      }
      
    }
    
    if(length(titles_controls)>1){
      
      distm=as.matrix(stringdist::stringdistmatrix(titles_controls,method="lv"))
      colnames(distm)=titles_controls
      rownames(distm)=titles_controls
      distm
      
      
      # s.a <- strsplit(titles[1], "")[[1]]
      # s.b <- strsplit(titles[2], "")[[1]]
      # paste(s.a[s.a != s.b], collapse = "")
      
      
      distm[which(distm>1)]=0
      
      g  <- graph.adjacency(distm) 
      #plot(g)
      clu <- components(g)
      gr=igraph::groups(clu)
      
      gr=gr[which(lengths(gr)>1)]
      
      if(length(gr)==0&last_control==0){
        print("No control samples")
        return(list())
      }else if(length(gr)>0){
        for(i in 1:length(gr)){
          #print(gr[[i]])
          sample_groups<-c(sample_groups,list(gr[[i]]))
          names(sample_groups)[length(sample_groups)]=paste("control",last_control+i,sep="")
          
        }
      }
      
      
      
      
    }
    sample_groups<-c(sample_groups,list(GSE=gse))
    sample_groups<-c(sample_groups,list(titles=titles))
    sample_groups<-c(sample_groups,list(total_sample_count=length(titles)))
    sample_groups<-c(sample_groups,list(parsed_sample_count=length(unlist(sample_groups))-length(titles)-2))
    
    return(sample_groups)
    
  }, error = function(err) {
    
    
    
    return(list())
    
  }) # END tryCatch
  
  
  
}

group_samples6<-function(gse,titles){
  
  library(stringdist)
  library(igraph)
  library(raster)
  
  
  
  tryCatch({
    
   
    
    titles_cases=titles
    rep_keys=get_rep_keys3(titles_cases)
    
    sample_groups=list()
    
    last_case=0
    
    if(length(rep_keys)>0){
      
      if(rep_keys!=""){
        rep_keys=rep_keys[order(nchar(rep_keys), rep_keys,decreasing=TRUE)]
        if(length(which(nchar(rep_keys)==0))){
          rep_keys=rep_keys[-which(nchar(rep_keys)==0)]
        }
        
        for(i in 1:length(rep_keys)){
          #print(gr[[i]])
          rep_keys[i]=gsub("\\+","\\\\+",rep_keys[i])
          rep_keys[i]=gsub("\\(","\\\\(",rep_keys[i])
          rep_keys[i]=gsub("\\?","\\\\?",rep_keys[i])
          rep_keys[i]=gsub("\\*","\\\\*",rep_keys[i])
          rep_keys[i]=gsub("\\{","\\\\{",rep_keys[i])
          rep_keys[i]=gsub("\\[","\\\\[",rep_keys[i])
          rep_keys[i]=gsub("\\|","\\\\|",rep_keys[i])
          
          matched=grep(rep_keys[i],titles_cases)
          if(length(matched)>2){
            last_case=last_case+1
            sample_groups<-c(sample_groups,list(titles_cases[matched]))
            titles_cases=titles_cases[-matched]
            names(sample_groups)[last_case]=paste("group",last_case,sep="")
            
          }

        }
      }
      
    }
    
    if(length(titles_cases)>2){
      
      distm=as.matrix(stringdist::stringdistmatrix(titles_cases,method="lv"))
      colnames(distm)=titles_cases
      rownames(distm)=titles_cases
      distm
      

      distm[which(distm>1)]=0
      
      g  <- graph.adjacency(distm) 

      clu <- components(g)
      gr=igraph::groups(clu)
      
      gr=gr[which(lengths(gr)>2)]
      
      if(length(gr)==0&last_case==0){
        print("No enough samples")
        return(list())
      }else if(length(gr)>0){
        for(i in 1:length(gr)){
          #print(gr[[i]])
          sample_groups<-c(sample_groups,list(gr[[i]]))
          names(sample_groups)[length(sample_groups)]=paste("group",last_case+i,sep="")
          
        }
      }
      
      
      
      
    }
    
    if(length(sample_groups) >1){
      
      sample_groups<-c(sample_groups,list(failed_samples=setdiff(titles,unlist(sample_groups))))
      sample_groups<-c(sample_groups,list(GSE=gse))
      sample_groups<-c(sample_groups,list(matched_term=match_gse_terms(gse)))
      sample_groups<-c(sample_groups,list(titles=titles))
      sample_groups<-c(sample_groups,list(total_sample_count=length(titles)))
      sample_groups<-c(sample_groups,list(parsed_sample_count=length(unlist(sample_groups))-length(titles)-2))
      
      return(sample_groups)
    }else{
      print("Not enough groups")
      return(list())
    }
    
    
   

    
  }, error = function(err) {
    
    
    
    return(list())
    
  }) # END tryCatch
  
  
  
}

group_samples7<-function(gse,titles,description,abstract){
  
  library(stringdist)
  library(igraph)
  library(raster)
  library(stringr)
  
 
  
  tryCatch({
    
    titles_original=titles
    
    located=str_locate_all(titles_original,"[0-9]+")
    
    if( length( which(lengths(located)>0) ) <6 ){ # there should be at least 3 reps for at least two conditions
      print("Didn't meet the criteria:there should be at least 3 reps for at least two conditions ")
      return(list())
    }else{
      
      sample_groups_past=NULL
      mean_lengths_past=100
     
      while(1){
        
        titles=apply(as.matrix( 1:length(titles)),1,function(x){
        p1=substr(titles[x],1,located[[x]][dim(located[[x]])[1],1]-1);
        p2=substr(titles[x],located[[x]][dim(located[[x]])[1],2]+1,nchar(titles[x]));
        paste0(p1,p2)})
        
        sample_groups=list()
        
        distm=as.matrix(stringdist::stringdistmatrix(titles,method="lv"))
        colnames(distm)=titles_original
        rownames(distm)=titles_original
        distm
        
        
        distm[which(distm>0)]=100
        distm[which(distm==0)]=1
        distm[which(distm==100)]=0
        
        g  <- graph.adjacency(distm) 
        
        clu <- components(g)
        gr=igraph::groups(clu)
        
        gr=gr[which(lengths(gr)>2)]
        
        if(length(gr)<2 ){
          print("Not enough samples")
          return(list())
        }else if(length( which(lengths(gr)>2)) >1){
          for(i in 1:length(gr)){
            #print(gr[[i]])
            sample_groups<-c(sample_groups,list(gr[[i]]))
            names(sample_groups)[length(sample_groups)]=paste("group",i,sep="")
            
          }
        }else{ print("Not enough samples") }
        
        if(mean_lengths_past > mean(lengths(sample_groups))){ #splitted better or equal
          
          break

        }else{
          
          sample_groups_past=sample_groups
          mean_lengths_past=mean(lengths(sample_groups))
          
        }
        
        located=str_locate_all(titles,"[0-9]+")
        if( length( which(lengths(located)>0) ) <6 ){break}
    }
        

      
      if(length(sample_groups) >=2){
        
        sample_groups<-c(sample_groups,list(failed_samples=setdiff(titles_original,unlist(sample_groups))))
        sample_groups<-c(sample_groups,list(GSE=gse))
        sample_groups<-c(sample_groups,list(matched_term=match_gse_terms2(description, abstract)))
        sample_groups<-c(sample_groups,list(titles=titles_original))
        sample_groups<-c(sample_groups,list(total_sample_count=length(titles_original)))
        sample_groups<-c(sample_groups,list(failed_sample_count=length(sample_groups$failed_samples)))
        
        return(sample_groups)
      }else{
        print("Not enough groups")
        return(list())
      }
      
      
    }
    #end
 
    
    
    
    
  }, error = function(err) {
    
    
    
    return(list())
    
  }) # END tryCatch
  
  
  
}
#######################run the function
#result_1: successed datasets
#result_2: failed results
#results: successful results
#among 507 disease matched samples, 75 succeeded
start_t=Sys.time()
matched_count=0

result_1_7=list()
result_0_7=list()
results7=list()
chipseq=NULL
scRNAseq=NULL
for(i in 501:length(unique_gse)){
  print(i)
  
  error_flag=0
  
  tryCatch({
    
    data <- getGEO(unique_gse[i])
    titles=as.character(data[[1]]@phenoData@data$title)
    description=data[[1]]@experimentData@title
    abstract=data[[1]]@experimentData@abstract
   
  }, error = function(err) {
    
    
    
    error_flag=1
    
  }) # END tryCatch
  
  
  if(error_flag==1){
    result_0_7<-c(result_0_7,list(list(gse_id=unique_gse[i],titles="Failed read")))
  }else{
    if(length(grep("chipseq",paste(titles,description),ignore.case = TRUE))>0 | length(grep("chip-seq",paste(titles,description),ignore.case = TRUE))>0){
      chipseq=c(chipseq,i)
      next
      
    }
    if(length(grep("single-cell",paste(titles,description),ignore.case = TRUE))>0  | length(grep("single cell",paste(titles,description),ignore.case = TRUE))>0  | length(grep("scRNAseq",paste(titles,description),ignore.case = TRUE))>0 | length(grep("scRNA-seq",paste(titles,description),ignore.case = TRUE))>0){
      chipseq=c(chipseq,i)
      next
      
    }
    
    group_samples_result=group_samples7(unique_gse[i],titles,description,abstract)
    
    
    if(length(group_samples_result)>0){
      results7<-c(results7,list(group_samples_result))
      result_1_7<-rbind(result_1_7,c(i,unique_gse[i]))
      if( length(group_samples_result$matched_term)!=0){
        matched_count=matched_count+1
      }
      
    }else{
      
      result_0_7<-c(result_0_7,list(list(gse_id=unique_gse[i],titles=titles)))
    }
  }
  
  
  
  
}
end_t=Sys.time()

save.image("parsing_v7.RData")

false_positive=c(20)
false_negative=c(15)
histone=c(15)

result_1=NULL
result_0=list()
results=list()
chipseq=NULL
scRNAseq=NULL

for(i in 1:length(which(matched_diseases_unlisted>0))){
#for(i in 1:length(which(matched_diseases_unlisted>0))){
print(i)
  data <- getGEO(unique_gse[which(matched_diseases_unlisted>0)[i]])
  titles=as.character(data[[1]]@phenoData@data$title)
  titles
  description=as.character(data[[1]]@phenoData@data$description)
  
  if(length(grep("chipseq",paste(titles,description),ignore.case = TRUE))>0 | length(grep("chip-seq",paste(titles,description),ignore.case = TRUE))>0){
    chipseq=c(chipseq,i)
    next
    
  }
  if(length(grep("single-cell",paste(titles,description),ignore.case = TRUE))>0  | length(grep("single cell",paste(titles,description),ignore.case = TRUE))>0  | length(grep("scRNAseq",paste(titles,description),ignore.case = TRUE))>0 | length(grep("scRNA-seq",paste(titles,description),ignore.case = TRUE))>0){
    chipseq=c(chipseq,i)
    next
    
  }
  
  group_samples_result=group_samples3(unique_gse[which(matched_diseases_unlisted>0)[i]],titles)
  

  if(length(group_samples_result)>0){
    results<-c(results,list(group_samples_result))
    result_1<-rbind(result_1,c(i,unique_gse[which(matched_diseases_unlisted>0)[i]]))
    
  }else{
    
    result_0<-c(result_0,list(list(gse_id=unique_gse[which(matched_diseases_unlisted>0)[i]],titles=titles)))
  }
  
}


disease_result_1=result_1
disease_result_0=result_0
disease_results=results
