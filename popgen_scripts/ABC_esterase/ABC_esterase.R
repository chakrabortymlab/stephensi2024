#Script used for ABC inference of time since esterase duplication sweep
#in the Trivandrum population
#script adapted from Ormond et al. 2015
#Author: Alex Samano 2024


#bottle-growth model parameters inferred by dadi
#effective population size
Ne<-1757732

#load required R packages
library("pls"); 
library ("abc"); 
library("MASS")


#1) Run the simulations for ABC inference
#the length of the sequence
L=40000
#the sample size 
samplesize=36
#the position of the selected mutation halfway along the sequence
Sp=0.5
#the starting frequency assuming a model of de novo mutation
fs=1/Ne
#Condition on theta, rho and the number of segregating sites
#theta
mu=1.36e-9
theta=4*Ne*mu*L
#rho
r= 1e-7
rho=4*Ne*r*L
#number of segregating sites
S=753

#run the simulations
nsim <-500000 

msstats_sims=c()
for (i in 1:nsim)  {
  selec_coef=0.15 #See Methods for calculation of selection coefficient
  alpha_homo=2*Ne*selec_coef  #set alpha for the homozygote and heterozygote genotypes using the selection coefficient
  alpha_hetero=(alpha_homo)/2 
  b <- runif(1,-5,-3) #draw Ts from a log uniform prior 
  Ts<-10^b
  # write these parameters to a file for input to the msms command line
  write.table(data.frame("Ts"=format(Ts,scientific=F),"alpha_homo"=format(alpha_homo,scientific=F),"alpha_hetero"=format(alpha_hetero,scientific=F)),row.names=FALSE, col.names=FALSE, quote=FALSE, file="msms_input.txt")
  count<-0
  count1<-0
  repeat  {
    #run the msms simulation
    system(paste("java -jar /usr/local/bin/msms/lib/msms.jar ",samplesize," 1 -N 1757732 -t ",theta," -r ",rho," -s ",S," -I 1 ",samplesize," -en ",Tb," 1 ",nuB," -eg ",Tb," 1 ",alpha," -Sp 0.5 -SI ",Ts," 1 ",fs," -SFC -SAA ",alpha_homo," -SAa ",alpha_hetero," -oTrace < msms_input.txt > ms_full_output.txt",sep=''))
    #write the SNP and frequency tracing output from the msms simulation to a file
    #separate the frequency tracing output from the standard msms output and calculate the frequency of the selected mutation in the population
    system("python3 get_freq.py ms_full_output.txt ms_freq_trace.txt")
    freq=read.table("/asteph_popgen/ABC_ansteph/ms_freq_trace.txt",header=F)
    freq <- (tail(freq,1)[[3]])
    count1<-count1 + 1
    
    # if the final frequency is above 80% calculate statistics from the msms output using mssstats 
    if (freq > 0.80) {
      count=count+1
      system("python3 get_ms.py ms_full_output.txt ms_standard.txt")
      ms_stats<-system("cat ms_standard.txt | /usr/local/bin/libsequence/examples/msstats > run_msstats.txt", intern=T)
      stats <-read.table(file="run_msstats.txt", header = F)
      
      # write the statistics to a file and move to next random draw of input variables
      all_stats <- cbind(freq,selec_coef,Ts,stats)
    }
  }
  msstats_sims <- rbind(msstats_sims,all_stats)
}

colnames(msstats_sims)<-c("freq","selec_coeff","Ts", "NumPoly","ThetaW","ThetaPi","ThetaH","TajD","FuLiD","FuLiF","FuLiDStar","FuLiFStar")
#write the statistics from the simulations with the freq, s, and Ts
write.table(msstats_sims, "sim500K_esterase.txt", sep ="\t",quote = F, row.names = F,)


#Read in the file with the result of 500K simulations
msstats<-read.table(file="sim500K_esterase.txt", header=TRUE) 

#2) Calculate PLS components for the statistics 

#calculate PLS loadings based on a subset of the first 10,000 simulations (script adapted from Wegmann et al. 2010).  
#The PLS calculation takes a few minutes.
a<- msstats[1:10000,]
stat<-a[,6:12]; 
param<-as.data.frame(a[,3]); #the number of segregating sites S is excluded from the calculation as we have conditioned on S in the msms simulations

#standardize the parameters and the statistics
mymeanparam <- c()
mysdparam <- c()
for(i in 1:length(param)){
  mymeanparam <- c(mymeanparam, mean(param[,i])); mysdparam <- c(mysdparam, sd(param[,i]));
  param[,i]<-(param[,i]-mean(param[,i]))/sd(param[,i]);}
myMax<-c(); myMin<-c(); lambda<-c(); myGM<-c();
for(i in 1:length(stat)){
  myMax<-c(myMax, max(stat[,i])); myMin<-c(myMin, min(stat[,i]));
  stat[,i]<-1+(stat[,i] -myMin[i])/(myMax[i]-myMin[ i]);
}
# apply Box-Cox transformation to normalize the statistics prior to PLS calculation
for(i in 1:length(stat)){
  d<-cbind(stat[,i], param);
  mylm<-lm(as.formula(d), data=d);
  myboxcox<-boxcox(mylm, lambda=seq(-20,100,1/10), interp=T, eps=1/50);
  lambda<-c(lambda, myboxcox$x[myboxcox$y==max(myboxcox$y)]);
  myGM<-c(myGM, mean(exp(log(stat[,i]))));}
#standardize the Box-Coxed statistics
myBCMeans<-c(); myBCSDs<-c();
for(i in 1:length(stat)){
  stat[,i]<-(stat[,i] ^lambda[i] - 1)/(lambda[i]*myGM[i] ^(lambda[i]-1));
  myBCSDs<-c(myBCSDs, sd(stat[,i]));
  myBCMeans<-c(myBCMeans, mean(stat[,i]));
  stat[,i]<-(stat[,i] -myBCMeans[i])/myBCSDs[i];
}
#Calculate the PLS components
myPlsr<-plsr(as.matrix(param)~as.matrix(stat), scale=F, validation='LOO');
#write the PLS information to a file
myPlsrDataFrame<-data.frame(comp1=myPlsr$loadings[,1]);
for(i in 2:7){myPlsrDataFrame<-cbind(myPlsrDataFrame, myPlsr$loadings[,i])}
write.table(cbind(names(stat), myMax, myMin, lambda, myGM, myBCMeans, myBCSDs,myPlsrDataFrame), file="plssummary_non_eqm.txt", sep="\t", col.names=F,row.names=F, quote=F);
# the reduction in RMSE for each parameter from using PLS can be visualised
plot(RMSEP(myPlsr));

#convert the statistics from the simulations into PLS components
b<-msstats
scores<-b[,6:12];
param<-as.data.frame(b[,3]); 


#standardize the statistics 
for(i in 1:length(scores)){
  scores[,i]<- 1+(scores[,i] -myMin[i])/(myMax[i]-myMin[ i]);
}
#apply the Box Cox transformation and normalize 
for(i in 1:length(scores)){
  scores[,i]<-(scores[,i] ^lambda[i] - 1)/(lambda[i]*myGM[i] ^(lambda[i]-1));
  scores[,i]<-(scores[,i] -myBCMeans[i])/myBCSDs[i];
}
#convert into PLS components (here the top 2 components are used)
scores <- (as.matrix(scores))%*%(as.matrix(myPlsrDataFrame[,1:2])) 


#3) Apply the method to the esterase polymorphism data from Trivandrum
#read in stats for 20Kb sequences flanking the esterase dup
test <-read.table("tvm_esterase_dup_flank_snps.msstats", header = FALSE)
#remove S and invariant statistics and convert the statistics using the PLS loadings (same steps as above)
test <- test[,3:9]
for(i in 1:length(test)){
  test[,i]<- 1+(test[,i] -myMin[i])/(myMax[i]-myMin[ i]);
}
for(i in 1:length(test)){
  test[,i]<-(test[,i] ^lambda[i] - 1)/(lambda[i]*myGM[i] ^(lambda[i]-1));
  test[,i]<-(test[,i] -myBCMeans[i])/myBCSDs[i];
}
newtest <- (as.matrix(test))%*%(as.matrix(myPlsrDataFrame[,1:2]))

# run the ABC inference
sim<-abc(target = c(newtest[1,c(1,2)]), param=param, sumstat=scores[,c(1,2)], method ="rejection", tol=0.05, transf="none")
post_Ts=(sim$unadj.values[,1])

#write post distribution to file
posteriorDist<-as.data.frame(post_Ts)
posteriorDist$generations<-4*Ne*posteriorDist$post_Ts
write.table(posteriorDist, "esterase_dup_postdist.txt", sep ="\t",quote = F, row.names = F,)

#get Ts from the mode of posterior probability distribution
posteriorDist<-read.table("/asteph_popgen/ABC_ansteph/esterase_dup_postdist.txt", header=T)
density_estimate <- density(posteriorDist$generations)
Ts <- density_estimate$x[which.max(density_estimate$y)]

#get confidence interval 
conf<-quantile(posteriorDist$post_Ts,probs=c(0.025,0.975))
low_est<-conf[1]
high_est<-conf[2]
low_est_gens<-4*Ne*low_est
high_est_gens<-4*Ne*high_est

#density plot of posterior distribution
ggplot(posteriorDist, aes(x = generations)) +
  geom_density(fill = "#de2d26", alpha = 0.5, size=0.1) +
  labs(title = "Posterior Distribution for Age of Esterase Duplication Sweep", x = "Generations Since Selective Sweep", y = "Density")+
  theme_bw()+guides(color="none")+
  scale_x_continuous(limits = c(0, max(posteriorDist$generations)),breaks=c(0,seq(50,max(posteriorDist$generations), by=50)))+
  geom_vline(aes(xintercept = Ts), size=0.5)+
  geom_vline(aes(xintercept = high_est_gens), size=0.5,linetype="dashed")+ 
  geom_vline(aes(xintercept = low_est_gens), size=0.5, linetype="dashed")



