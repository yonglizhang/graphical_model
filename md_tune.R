
options(scipen=999)
vvv<-round(abs(rnorm(1)*100000000))

library(MASS)

.libPaths("/ibrix/home8/yongli/RCODE")
library(lattice)
.libPaths("/ibrix/home8/yongli/RCODE")
library(Matrix)
.libPaths("/ibrix/home8/yongli/RCODE")
library(QUIC)
.libPaths("/ibrix/home8/yongli/RCODE")
library(ncvreg)
.libPaths("/ibrix/home8/yongli/RCODE")
library(lars)
.libPaths("/ibrix/home8/yongli/RCODE")
library(glmnet)
.libPaths("/ibrix/home8/yongli/RCODE")
library(igraph)
.libPaths("/ibrix/home8/yongli/RCODE")
library(glasso)
#.libPaths("/ibrix/home8/yongli/RCODE")
#library("LambertW")

#source("tlp_matrix.R")
#source("C:/Users/yongli/Dropbox/Graphical/finnetwork/tlp_matrix.R")


ppmm<-20;



#md<-read.csv("md.csv", header = TRUE)
md<-read.csv("C:/Users/yongli/Dropbox/Graphical/finnetwork/md.csv", header = TRUE)
dim(md)
colnames(md)[-1]
p<-ncol(md)-1; p
#sel10<-c(0,1:10,21:30,41:50,61:70,81:90,101:110,121:130,141:150,161:170,181:190)
#md<-md[,sel10+1]
#dim(md)

#md <-subset(md,(as.Date(Date)<=as.Date("2009-08-31"))&(as.Date(Date)>=as.Date("2008-03-25") )  )
#dim(md)

############################
md_before<-subset(md,(as.Date(Date)>=as.Date("2005-01-01"))&(as.Date(Date)<as.Date("2005-12-31")))
sum(apply(is.na(md_before),1,sum))
md_before<-as.matrix(md_before[,-1]); 
logrn_before<-((log(md_before[-1,])-log(md_before[-nrow(md_before),]))/0.01)[-1,]
dim(logrn_before)
logrn_before<-scale(logrn_before, center = TRUE, scale = F)
apply(logrn_before,2,sd)
dim(logrn_before)
sum(logrn_before)


##########################
md_middle<-subset(md,(as.Date(Date)>=as.Date("2008-07-01"))&(as.Date(Date)<=as.Date("2009-06-30"))    )
md_middle<-as.matrix(md_middle[,-1]); 
logrn_middle<-((log(md_middle[-1,])-log(md_middle[-nrow(md_middle),]))/0.01 )[-1,]
dim(logrn_middle)
logrn_middle<-scale(logrn_middle, center = TRUE, scale = F)
apply(logrn_middle,2,sd)
dim(logrn_middle)
sum(logrn_middle)



##############################
md_after<-subset(md,(as.Date(Date)>=as.Date("2010-01-01"))&(as.Date(Date)<=as.Date("2010-12-31"))     )
md_after<-as.matrix(md_after[,-1]); logrn_after<-((log(md_after[-1,])-log(md_after[-nrow(md_after),]))/0.01)[-1,]
dim(logrn_after)
ogrn_after<-scale(logrn_after, center = TRUE, scale = F)
apply(logrn_after,2,sd)
dim(logrn_after)
sum(logrn_after)


logrn_all<-rbind(logrn_before,logrn_middle,logrn_after)
temp<-list(logrn_before,logrn_middle,logrn_after)


################################################################################
################################################################################

temp_matrix_t<-logrn_all

stocknames<-colnames(logrn_all)
stocknames

S_bar_m<-((cov(logrn_all,use="pairwise.complete.obs")) )
cov_all<-S_bar_m


cov_before<-(cov(logrn_before,use="pairwise.complete.obs"))
cor_before<-cov2cor(cov_before)
sum(abs(cor_before))
mean(diag(cov_before))



cov_middle<-(cov(logrn_middle,use="pairwise.complete.obs"))
cor_middle<-cov2cor(cov_middle)
sum(abs(cor_middle))
mean(diag(cov_middle))


cov_after<-(cov(logrn_after,use="pairwise.complete.obs"))
cor_after<-cov2cor(cov_after)
sum(abs(cor_after))
mean(diag(cov_after))



S_bar <-cbind(cov_before,cov_middle,cov_after)


L <- num_of_matrix <- 3

nn <- c(nrow(logrn_before), nrow(logrn_middle), nrow(logrn_after))
nn


################################################################################
################################################################################
################################################################################
################################################################################
Lambda_quic_ho<-mean(diag(cov_all))*sqrt(4*log(p)/sum(nn));
Lambda_quic_ha<-c( mean(diag(cov_before))*sqrt( 4*log(p)/nn[1]),mean(diag(cov_middle))*sqrt(4*log(p)/nn[2]),mean(diag(cov_after))*sqrt(4*log(p)/nn[3])  ); 

Lambda_quic_ho
Lambda_quic_ha





################################################################################
################################################################################
################################################################################
################################################################################



Lambda_quic_ho
Lambda_quic_ha<-rep(0,L)



#######################################
XX<-logrn_all
lak<-50; lamax<-mean(diag(cov(XX))); lamin<-mean(diag(cov(XX)))/sqrt(nrow(XX)); kf<-2;  II<-1;
la_seq<-exp(seq(log(lamin), log(lamax),length=lak))
la_seq
CV_m1<-matrix(0,kf*II,lak); CV_m2<-matrix(0,kf*II,lak);

for (ii in 1:II)
    { 
    am<-split(sample(1:nrow(XX),nrow(XX),replace=F),1:kf);

    for (jj in 1:kf)
                  {
                  val_ind<-am[[jj]];

                  
                  S_bar1 <-cov(XX[-val_ind,]) 
                  V_bar1 <-cov(XX[val_ind,])   
                  for (kk in 1:lak)
                                  {                                            

                                  Omega1<-glasso(S_bar1, penalize.diagonal = FALSE, rho =la_seq[kk])$wi

                                  
                                  CV_m1[jj+(ii-1)*kf,kk]<- -log(det(Omega1))+sum(diag(V_bar1 %*% Omega1))
                                  }
                  
                  S_bar1<-NULL; V_bar1<-NULL; 
                  }
    }

CV_mean1<-apply(CV_m1, 2,mean)
CV_mean1
la1<-(la_seq[CV_mean1==min(CV_mean1)],na.rm=T)
la1        
Lambda_quic_ho<-la1
Lambda_quic_ho
              
              
              
              

for (i in 1:L){

              #i<-2
              XX<-temp[[i]]
              lak<-50; lamax<-mean(diag(cov(XX))); lamin<-mean(diag(cov(XX)))/sqrt(nrow(XX)); kf<-5;  II<-1;
              la_seq<-exp(seq(log(lamin), log(lamax),length=lak))
              la_seq
              CV_m1<-matrix(0,kf*II,lak); CV_m2<-matrix(0,kf*II,lak);
              
              for (ii in 1:II)
                            { 
                            am<-split(sample(1:nrow(XX),nrow(XX),replace=F),1:kf);
                        
                            for (jj in 1:kf)
                                          {
                                          val_ind<-am[[jj]];
                        
                                          
                                          S_bar1 <-cov(XX[-val_ind,]) 
                                          V_bar1 <-cov(XX[val_ind,])   
                                          for (kk in 1:lak)
                                                          {                                            
               
                                                          Omega1<-glasso(S_bar1, penalize.diagonal = FALSE, rho =la_seq[kk])$wi
                                                          
                                                          #Omega2<-(TLP_matrix(S_bar1,Lambda=la_seq[kk],weighted=FALSE)$X.tlp)[,,1]
                                                          
                                                          CV_m1[jj+(ii-1)*kf,kk]<- -log(det(Omega1))+sum(diag(V_bar1 %*% Omega1))
                                                          #CV_m2[jj+(ii-1)*kf,kk]<- -log(det(Omega2))+sum(diag(V_bar1 %*% Omega2))
                                              
                                                          }
                                          
                                          S_bar1<-NULL; V_bar1<-NULL; 
                                          }
                            }
              
              CV_mean1<-apply(CV_m1, 2,mean)
              CV_mean1
              la1<-min(la_seq[CV_mean1==min(CV_mean1)],na.rm=T)        
              la1
              Lambda_quic_ha[i]<-la1
              
              }
Lambda_quic_ha              
              
              
              
