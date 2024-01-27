library(diptest)

x<-seq(from=1,to=4,by=0.5)
seq1<-10^x # size of the samples that they're taking - 
seq2<-c(176,c(c(1:12)*2+176)) # heights - this is a fancy way of going up by twos


# so the pairs is the matrix of hieghts x sample sizes
pairs<-matrix(data=NA, nrow=length(seq1), ncol=length(seq2))
colnames(pairs)<-seq2
row.names(pairs)<-seq1



for(k in 1:length(seq1)){
  
  for(h in 1: length(seq2)){
    
    g<-seq1[k] # g is the sample size
    w<-seq2[h] # w is the height of the males
    
    Pvalue<-0
    
    iterations<-1000
    for(x in 1:iterations){
    
    m<-rnorm(g,w,7) # male sample
    f<-rnorm(g,162,7) # female sample
    
    c<-c(m,f) # this is the combined sample (male and female)
    
    dtest<-dip.test(c,simulate.p.value=TRUE, B=10000) # so they are getting the dip test p value the same way
    pval<-dtest$p.value
    
    Pvalue<-Pvalue+pval # sum of current p values
    
    }
    
    pairs[as.character(g),as.character(w)]<-Pvalue/iterations # this is getting the average p vale
    
    
    
  }
  
}


a<-dim(pairs)[1]*dim(pairs)[2]
output<-matrix(data=NA, ncol=3, nrow=a)
colnames(output)<-c("S.size", "M.Height", "P.Val")


for(y in 1:length(output[,1])){
  
      r<-ceiling(y/length(pairs[1,]))
      u<-(y-(ceiling(y/length(pairs[1,]))-1)*length(pairs[1,]))
  
      value<-pairs[r,u]
      
      output[y,1]<-as.numeric(row.names(pairs)[r])
      output[y,2]<-as.numeric(colnames(pairs)[u])
      output[y,3]<-as.numeric(value)
    
}



sig<-subset(output, output[,3]<=0.05, select=c(1:2))
n.sig<-subset(output,output[,3]>0.05, select=c(1:2))

tiff(filename="diptest.plot.tiff", width=150, height=86, units="mm", res=600)
par(mar=c(5,5,2,3), cex=0.65, lwd=0.5)
plot(output[,2],output[,1], log="y", ylab=expression("Total sample size (1:1 male:female)"),xlab=expression(paste("Average (", mu, ") male height (cm)")),yaxt="n",xaxt="n", col="WHITE", pch=".", cex=1, family="sans", ps=8)
points(n.sig[,2],n.sig[,1], col="RED", pch=16)
points(sig[,2],sig[,1], col="BLACK", pch=16)
axis(side=1, labels=(seq2), at=(seq2), lwd=0.5, lwd.ticks=0.5, cex=1)
axis(side=2, labels=(trunc(seq1)*2), at=(seq1), lwd=0.5, lwd.ticks = 0.5, cex=1)
dev.off()

