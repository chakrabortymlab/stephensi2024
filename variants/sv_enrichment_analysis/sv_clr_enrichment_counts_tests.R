#compare shuffled distribution to actual value

#actual number of SVs over 0.25 allele frequency, associated with CLR peak
bang_obs<-204
mang_obs<-147
kochi_obs<-252
tvm_obs<-261
lkd_obs<-255

#function to get pvalue, from North et al. 2002
get_pval<- function(distribution, observed){
  extremes<-subset(distribution,Count >=observed)
  pval<-(1+nrow(extremes))/100001
  return(pval)
}


#bang
bang_shuffled <-read.table("bang_100Kshuff_svs_nearpeak.txt", header=T)
bang<-ggplot(bang_shuffled, aes(x = Count)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(title = "Bangalore", x = "# SVs Associated with CLR peak", y = "Frequency")+
  geom_vline(aes(xintercept =bang_obs, color="red"),size=1)+ theme_bw()+guides(color="none")

bang_pval<- get_pval(bang_shuffled, bang_obs)
#bang pval 9.9999e-06


#mang
mang_shuffled <-read.table("mang_100Kshuff_svs_nearpeak.txt", header=T)
mang<-ggplot(mang_shuffled, aes(x = Count)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(title = "Mangalore", x = "# SVs Associated with CLR peak", y = "Frequency")+
  geom_vline(aes(xintercept =mang_obs, color="red"),size=1)+ theme_bw()+guides(color="none")

mang_pval<- get_pval(mang_shuffled, mang_obs)
#mang pval 0.0001399986

#kochi
kochi_shuffled <-read.table("kochi_100Kshuff_svs_nearpeak.txt", header=T)
kochi<-ggplot(kochi_shuffled, aes(x = Count)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(title = "Kochi", x = "# SVs Associated with CLR peak", y = "Frequency")+
  geom_vline(aes(xintercept =kochi_obs, color="red"),size=1)+ theme_bw()+guides(color="none")

kochi_pval<- get_pval(kochi_shuffled, kochi_obs)
#kochi pval 0.04501955

#tvm
tvm_shuffled <-read.table("tvm_100Kshuff_svs_nearpeak.txt", header=T)
tvm<-ggplot(tvm_shuffled, aes(x = Count)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(title = "Trivandrum", x = "# SVs Associated with CLR peak", y = "Frequency")+
  geom_vline(aes(xintercept =tvm_obs, color="red"),size=1)+ theme_bw()+guides(color="none")

tvm_pval<- get_pval(tvm_shuffled, tvm_obs)
#tvm pval 0.03476965



#lkd
lkd_shuffled <-read.table("/Users/Alex/Desktop/ChakLab/asteph_sv/fivepops115/sv_enrichment/lkd_100Kshuff_svs_nearpeak.txt", header=T)
lkd<-ggplot(lkd_shuffled, aes(x = Count)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(title = "Lakshadweep", x = "# SVs Associated with CLR peak", y = "Frequency")+
  geom_vline(aes(xintercept =lkd_obs, color="red"),size=1)+ theme_bw()+guides(color="none")

lkd_pval<- get_pval(lkd_shuffled, lkd_obs)
#lkd pval 0.04473955


grid.arrange(bang,mang,kochi,lkd,tvm,nrow=2)
