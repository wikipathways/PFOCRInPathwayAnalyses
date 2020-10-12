# R script to download selected samples
# Copy code and run on a local machine to initiate download

#https://docs.google.com/document/d/1jm1sJtHgJkSV7hEjMC_VK7l7gzshMQzr4T7NdgY9rGE/edit#heading=h.skn982pchl79
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4930834/

setwd("/Users/mingyoungshin/git/PFOCRInPathwayAnalyses")
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
terms=c(terms, dieases,"knock-out", "knock out", "knockout","siRNA","knockdown", "knock down", "knock-down", "crispr-i", "variant","pathway")#signaling
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


histone=NULL
matched_histone=NULL
duplicates=NULL
errors=NULL
##additional fitering
for(i in 1:length(results7)){
  print(i)
  print(results7[[i]]$GSE)
  
  if(i>1){
    if(length(which((results7[[i-1]]$titles==results7[[i]]$titles) ==TRUE))==length(results7[[i-1]]$titles)){
      duplicates=c(duplicates,i)
      next
    }
  }
  
  data=NULL
  tryCatch({
    
    data <- getGEO(results7[[i]]$GSE,getGPL = FALSE)
    
  }, error = function(err) {
    
    
  }) # END 
  
  if(length(data)==0){ errors=c(errors,i);next}
  titles=as.character(data[[1]]@phenoData@data$title)
  description=data[[1]]@experimentData@title
  abstract=data[[1]]@experimentData@abstract
  

  
  
  if(length(grep("histone",paste(titles,description),ignore.case = TRUE))>0 | length(grep("H[1-4]+[A-Z][0-9]+",paste(titles,description),ignore.case = TRUE))>0){
    histone=c(histone,i)
    if(length(grep("H[1-4]+[A-Z][0-9]+",paste(titles,description),ignore.case = TRUE))>0){
      matched_histone=c(matched_histone,paste(paste(titles,collapse=" "),description,collapse=" "))
    }else{
      matched_histone=c(matched_histone,"histone")
    }
    
    next
    
  }
  
  
}

histone=histone[-c(21,26)]
keywords=c(histone,duplicates)


results7_no_histone=list()
for(i in 1:length(results7)){
  if(!(i%in%keywords)){
    results7_no_histone<-c(results7_no_histone,list(results7[[i]]))
  }
}

#1792 ->1732 gses
groupC=apply(as.matrix(1:length(results7_no_histone)),1,function(x) length(results7_no_histone[[x]])-6 )
group_sampleC=apply(as.matrix(1:length(results7_no_histone)),1,function(x){groups=lengths(results7_no_histone[[x]]);max(lengths(results7_no_histone[[x]])[grep("group",names(groups))])}  )

filter1=which(groupC==2)#689
filter2=(intersect(filter1,which(group_sampleC<=40)))#671

results7_noHistone_1comparison_2groups=list()
for(i in 1:length(results7_no_histone)){
  if((i%in%filter2)){
    results7_noHistone_1comparison_2groups<-c(results7_noHistone_1comparison_2groups,list(results7_no_histone[[i]]))
  }
}

false_positive=c(20)
false_negative=c(15)


save.image("parsing_v7.2_results7_noHistone_1comparison_2groups.RData")


verifyItems <- function(list = NULL, range = NULL){
  if(is.null(list))
    stop("Must provide a list")
  
  if(is.null(range))
    range <- 1:length(list)
  
  if(!is.vector(list[range]))
    stop("Range must subset provided list.")
  
  sapply(list[range], function(x){
    print(x)
    r<-NULL
    while(!is.logical(r) ){
      r<-toupper(readline("Valid? (T/F): "))
      if (r %in% c("T","F"))
        r <- as.logical(r)
    }
    r
  })
}

fruits <- verifyItems(results7_noHistone_1comparison_2groups, 225:448)
write.table(t(fruits),file="255_448",col.names=FALSE,row.names=FALSE)
save.image("Min_parsed1.RData")


##### True False match
fruits <- verifyItems(results7_noHistone_1comparison_2groups, 225:448)

fruits2 <- verifyItems(results7_noHistone_1comparison_2groups[(225:448)[fruits]],1:199)

save.image("Min_groups_filtered_twice.RData")
second_set=results7_noHistone_1comparison_2groups[ ((225:448)[fruits])[fruits2] ]
#####

##### combine all indices
third_set<-read.csv("true_positive_449_671_results7_noHistone_1comparison_2groups.csv")
length(which(third_set$Curation=="T,N"))
third_set=results7_noHistone_1comparison_2groups[which(third_set$Curation=="T,N")]

load("true_positive_1_100.RData")
load("true_positive_101_224.RData")
second_set_1=results7_noHistone_1comparison_2groups[(1:100)[true_positive_1_100]]
second_set_2=results7_noHistone_1comparison_2groups[(101:224)[true_positive_101_224]]

merged_filtered_gses=c(second_set_1,second_set_2,second_set,third_set)

matched_count=0

result_1_7=list()
result_0_7=list()
results7=list()
chipseq=NULL
scRNAseq=NULL
for(i in 1:length(unique_gse)){
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

save.image("parsing_v7.2.RData")

groupC=apply(as.matrix(1:length(results7)),1,function(x) length(results7[[x]])-6 )
hist(groupC,200)
hist(groupC,200,ylim=c(1,15))
hist((groupC^2-groupC)/2,200)
hist((groupC^2-groupC)/2,200,ylim=c(1,15))

hist(groupC[which(groupC<10)],200)
hist(groupC[which(groupC<10)],200,ylim=c(1,15))
hist((groupC[which(groupC<10)]^2-groupC[which(groupC<10)])/2,200)
hist((groupC[which(groupC<10)]^2-groupC[which(groupC<10)])/2,200,ylim=c(1,15))

group_sampleC=apply(as.matrix(1:length(results7)),1,function(x){groups=lengths(results7[[x]]);max(lengths(results7[[x]])[grep("group",names(groups))])}  )
hist(group_sampleC,200)
hist(group_sampleC,200,ylim=c(1,15))

filter1=which(groupC==2)
length(intersect(filter1,which(group_sampleC<=40)))
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
