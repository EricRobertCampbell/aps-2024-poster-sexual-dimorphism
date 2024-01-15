#Alligator analyses
#Male curve from Wilkinson and Rhodes (1997)
x=runif(1000,min=0, max=50)
Lm=3.79*(1-0.94*exp(-0.0695*x))
Lmerrors=rnorm(1000,mean=Lm,sd=0.0589*log(Lm)+0.0816) #SD of errors increases logarithmically, following Wilkinson and Rhodes (1997)
plot(x,Lmerrors, xlab="Age (yr)", ylab="Total length (m)", pch=16, col="lightskyblue", ylim=c(0,4))
curve(3.79*(1-0.94*exp(-0.0695*x)), lwd=3, col="red", add=T)
curve1<-data.frame(cbind(x,Lmerrors))# added this to create one object with all the data for subsetting below

#Female curve from Wilkinson and Rhodes (1997)
x2=runif(1000,min=0, max=50)
Lf=2.78*(1-0.91*exp(-0.0926*x2))
Lferrors=rnorm(1000,mean=Lf,sd=0.0332*log(Lf)+0.046)
points(x2,Lferrors, pch=16, col="goldenrod")
curve(2.78*(1-0.91*exp(-0.0926*x)), lwd=3, col="red", add=T)
curve2<-data.frame(cbind(x2,Lferrors))# one object

# subsampling code is below
sampling_results<-matrix(ncol=2,nrow=1000)
colnames(sampling_results)<-c("T value","p value")

# sampling from both sizes of the line (10) 1000 times
for(i in 1:1000){ # for each iteration in 1000 iterations
  #for each curve subset down to those less than x = 10
  samped<-sample(curve1$Lmerrors,size=30,replace=FALSE)# sample size 30 with no replacement
  samped2<-sample(curve2$Lferrors,size=30,replace=FALSE)
  sig_test<-t.test(samped,samped2)
  sampling_results[i,1]<-sig_test$statistic
  sampling_results[i,2]<-sig_test$p.value
  
  print(i)
  }


#for sampling true population structure
#for males
ages1=runif(1000,min=0, max=50)
comp1<-1295.2*exp(-1*(ages1+87.185)^2/(2*719.76))
subsamples1<-sample(ages1,1000,replace=TRUE,prob=comp1)
Lm=3.79*(1-0.94*exp(-0.0695*subsamples1))
Lmerrors=rnorm(1000,mean=Lm,sd=0.0589*log(Lm)+0.0816)
plot(subsamples1,Lmerrors, xlab="Age (yr)", ylab="Total length (m)", pch=16, col="lightskyblue", ylim=c(0,4))
curve(3.79*(1-0.94*exp(-0.0695*x)), lwd=3, col="red", add=T)
curve1<-data.frame(cbind(subsamples1,Lmerrors))# added this to create one object with all the data for subsetting below

#for females
ages2=runif(1000,min=0, max=50)
comp2<-1295.2*exp(-1*(ages2+87.185)^2/(2*719.76))
subsamples2<-sample(ages2,1000,replace=TRUE,prob=comp2)
Lf=2.78*(1-0.91*exp(-0.0926*subsamples2))
Lferrors=rnorm(1000,mean=Lf,sd=0.0332*log(Lf)+0.046)
points(subsamples2,Lferrors, pch=16, col="goldenrod")
curve(2.78*(1-0.91*exp(-0.0926*x)), lwd=3, col="red", add=T)
curve2<-data.frame(cbind(subsamples2,Lferrors))# one object

# subsampling code is below
sampling_results<-matrix(ncol=2,nrow=1000)
colnames(sampling_results)<-c("T value","p value")

# sampling from both sizes of the line (10) 1000 times
for(i in 1:1000){ # for each iteration in 1000 iterations
  #for each curve subset down to those less than x = 10
  samped<-sample(curve1$Lmerrors,size=30,replace=FALSE)# sample size 30 with no replacement
  samped2<-sample(curve2$Lferrors,size=30,replace=FALSE)
  sig_test<-t.test(samped,samped2)
  sampling_results[i,1]<-sig_test$statistic
  sampling_results[i,2]<-sig_test$p.value
 
  print(i)
  }

#for taphonomic bias (sampling >60 kg, >2.45 m)
sampling_results_greaterthan2.45<-matrix(ncol=2,nrow=1000)
colnames(sampling_results_greaterthan2.45)<-c("T value","p value")
# sampling from both sizes of the line (10) 1000 times
for(i in 1:1000){ # for each iteration in 1000 iterations
  curve1_temp<-subset(curve1,curve1$x>2.45)
  samped<-sample(curve1_temp$Lmerrors,size=30,replace=FALSE)
  curve2_temp<-subset(curve2,curve2$x>2.45)
  samped2<-sample(curve2_temp$Lferrors,size=30,replace=FALSE)
  sig_test<-t.test(samped,samped2)
  sampling_results_greaterthan2.45[i,1]<-sig_test$statistic
  sampling_results_greaterthan2.45[i,2]<-sig_test$p.value

  print(i)  
  }

#Rhea analyses
#Male curve from Navarro et al. (2005)
x=runif(1000,min=0, max=10.5)
Mm=0.3544+(28.7063*exp(-1*exp(-10.85*(x-0.37369863)))) #equation of Navarro et al. (2005) modified to include starting size
Mmerrors=rnorm(1000,mean=Mm,sd=0.425*log(Mm)+0.474) #fit errors with logarithmically inceasing st dev
plot(x,Mmerrors, xlab="Age (years)", ylab="Mass (kg)", pch=16, col="lightskyblue", xlim=c(0,10.5))
curve(0.3544+(28.7063*exp(-1*exp(-10.85*(x-0.37369863)))), lwd=3, col="red", add=T)
curve1<-data.frame(cbind(x,Mmerrors))# added this to create one object with all the data for subsetting below

#Female curve from Navarro et al. (2005)
x2=runif(1000,min=0, max=10.5)
Mf=0.3544+(22.5074*exp(-1*exp(-12.12*(x2-0.32191781))))
Mferrors=rnorm(1000,mean=Mf,sd=0.6793*log(Mf)+0.7378)
points(x2,Mferrors, pch=16, col="goldenrod")
curve(0.3544+(22.5074*exp(-1*exp(-12.12*(x-0.32191781)))), lwd=3, col="red", add=T)
curve2<-data.frame(cbind(x2,Mferrors))# one object

# Iterative t-testing
sampling_results<-matrix(ncol=2,nrow=1000)
colnames(sampling_results)<-c("T value","p value")

for(i in 1:1000){ # for each iteration in 1000 iterations
  samped<-sample(curve1$Mmerrors,size=30,replace=FALSE)# sample size 30 with no replacement
  samped2<-sample(curve2$Mferrors,size=30,replace=FALSE)
  sig_test<-t.test(samped,samped2)
  sampling_results[i,1]<-sig_test$statistic
  sampling_results[i,2]<-sig_test$p.value
  
  print(i)
  }

