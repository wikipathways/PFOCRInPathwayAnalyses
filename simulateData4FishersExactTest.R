rm(list = ls())
setwd("/Users/reubenthomas/Dropbox (Gladstone)/PFOCRInPathwayAnalyses")

# solution of quadratic equation
result <- function(a,b,c){
  if(delta(a,b,c) > 0){ # first case D>0
    x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
    x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
    result = c(x_1,x_2)
    return(round(min(result)))
  }
  else if(delta(a,b,c) == 0){ # second case D=0
    x = -b/(2*a)
    return(round(x))
  }
  else {"There are no real roots."} # third case D<0
}

# Constructing delta
delta<-function(a,b,c){
  b^2-4*a*c
}


set.seed(1234)
##total no of genes
N=10000
##total no of DE genes
Ne=1000
##total no of genes in pathway
Np=c(10,20,50,100,500,1000,2000,5000, 7000, 9500)
##Odds ratio
lambda=2

##Coefficients of quadratic equation
pValues <- vector(mode = "numeric")
OddsRatio <- vector(mode = "numeric")
for(i in 1:length(Np)) {
  print(i)
  a <- (lambda - 1)
  b <- -((lambda -1)*(Ne + Np[i]) + N)
  c <- lambda*Ne*Np[i]
  
  
  
  l = result(a,b,c)
  o = (N - (Ne + Np[i]) + l)
  m = (Ne - l)
  n = (Np[i] - l)
  
  X = cbind(c(l,n), c(m,o))
  
  res <- fisher.test(X)
  OddsRatio[i] <- res$estimate
  pValues[i] <- res$p.value
}

plot_data <- data.frame(Ngenes=Np, OddsRatio=OddsRatio, pvalue=pValues)

require(ggplot2)
pdf("pathway_size_vs_fisher_exact_test_significance.pdf")
ggplot(plot_data, aes(x=Ngenes,y=OddsRatio)) + 
  geom_point() + 
  geom_line() + 
  xlab("Pathway size") + 
  ylab("Odds ratio") + 
  ylim(0,3) +
  theme(text = element_text(size = 10))

ggplot(plot_data, aes(x=Ngenes,y=-log10(pvalue))) + 
    geom_point() + 
    geom_line() + 
    xlab("Pathway size") + 
    ylab("-log10(pvalue)") +
  theme(text = element_text(size = 10))
dev.off()  

