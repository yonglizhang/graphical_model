
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


ppmm<-1000;



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
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################




Lambda_quic_ho<-mean(diag(cov_all))*sqrt(4*log(p)/sum(nn));
Lambda_quic_ha<-c(mean(diag(cov_before))*sqrt( 4*log(p)/nn[1]),mean(diag(cov_middle))*sqrt(4*log(p)/nn[2]),mean(diag(cov_after))*sqrt(4*log(p)/nn[3])  ); 

Lambda_quic_ho
Lambda_quic_ha

offd_sel<-diag(rep(1,p))==0
offd_sel<-col(offd_sel)<row(offd_sel)
sum(offd_sel)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

Omega1_ho<-glasso(cov_all, penalize.diagonal = FALSE, rho =Lambda_quic_ho)$wi

logli1<-rep(0,L); 

for (i in 1:L)
              {
              
              #i<-1
              S_bar1<-S_bar[,((i-1)*p+1):(i*p)];
              sum(S_bar1)
              
              
              logli1[i]<-nn[i]*log( det(Omega1_ho))-nn[i]*sum(diag(S_bar1 %*% Omega1_ho))


              }
             
logli_quic_ho<- sum(logli1)

################################
################################


Omega1_ha_seq<-list(L);
pcor_ha_seq<-list(L);

logli1_ha<-rep(0,L); 
                  
for (i in 1:L)
            {
            #i<-3
            S_bar1<-S_bar[,((i-1)*p+1):(i*p)]
            
            Omega1_ha<-glasso(S_bar1, penalize.diagonal = FALSE, rho =Lambda_quic_ha[i])$wi
            
            Omega1_ha_seq[[i]]<-Omega1_ha                         
            pcor_ha_seq[[i]]<-cov2cor(Omega1_ha)
            logli1_ha[i]<-nn[i]*log(det(Omega1_ha))-nn[i]*sum(diag(S_bar1%*% Omega1_ha))

            }
    
logli_quic_ha<- sum(logli1_ha)

logli_quic<-logli_quic_ha-logli_quic_ho
logli_quic


(sum(Omega1_ha_seq[[1]]!=0)-p)/2
(sum(Omega1_ha_seq[[2]]!=0)-p)/2
(sum(Omega1_ha_seq[[3]]!=0)-p)/2





################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
logli_quic_seq_bp<-rep(0,ppmm)



for (b in 1:ppmm)


{


bp<-sample(1:sum(nn),sum(nn),replace=T)
logrn_all_bp<-logrn_all[bp, ]


logrn_before_bp<-logrn_all_bp[1:nrow(logrn_before),]
logrn_middle_bp<-logrn_all_bp[(nrow(logrn_before)+1):(nrow(logrn_before)+nrow(logrn_middle)),]
logrn_after_bp<-logrn_all_bp[(nrow(logrn_before)+nrow(logrn_middle)+1):nrow(logrn_all),]

cov_before_bp<-(cov(logrn_before_bp,use="pairwise.complete.obs"))
cov_middle_bp<-(cov(logrn_middle_bp,use="pairwise.complete.obs"))
cov_after_bp<-(cov(logrn_after_bp,use="pairwise.complete.obs"))

S_bar_m_bp<-cov(logrn_all_bp,use="pairwise.complete.obs") 
cov_all_bp<-S_bar_m_bp


Lambda_quic_ho<-mean(diag(cov_all_bp))*sqrt(4*log(p)/sum(nn));
Lambda_quic_ha<-c(mean(diag(cov_before_bp))*sqrt( 4*log(p)/nn[1]),mean(diag(cov_middle_bp))*sqrt(4*log(p)/nn[2]),mean(diag(cov_after_bp))*sqrt(4*log(p)/nn[3])  ); 



Omega1_bp<-glasso(cov_all_bp, penalize.diagonal = FALSE, rho =Lambda_quic_ho)$wi
#Omega1
#Omega2<-(TLP_matrix(S_bar_m,Lambda=Lambda_tlp_ho,tau=1e-2,tol=1e-4,maxIter=1000, max.dc.iter=20,weighted=FALSE)$X.tlp)[,,1]

S_bar_bp <-cbind(cov_before_bp,cov_middle_bp,cov_after_bp )


logli1_bp<-rep(0,L);     

for (i in 1:L)
              {

              S_bar1_bp<-S_bar_bp[,((i-1)*p+1):(i*p)];
              
              
              logli1_bp[i]<-nn[i]*log( det(Omega1_bp))-nn[i]*sum(diag(S_bar1_bp %*% Omega1_bp))


              }             
logli_quic_ho_bp<- sum(logli1_bp)

    
Omega1<-NULL; Omega2<-NULL



Omega1_seq_bp<-list(L);
logli1_bp<-rep(0,L); 
                  
for (i in 1:L)
            {
            #i<-3
            S_bar1_bp<-S_bar_bp[,((i-1)*p+1):(i*p)]
            
            Omega1_bp<-glasso(S_bar1_bp, penalize.diagonal = FALSE, rho =Lambda_quic_ha[i])$wi
            
            Omega1_seq_bp[[i]]<-Omega1_bp                         
            
            logli1_bp[i]<-nn[i]*log(det(Omega1_bp))-nn[i]*sum(diag(S_bar1_bp%*% Omega1_bp))

            }
    
logli_quic_ha_bp<- sum(logli1_bp)

logli_quic_ha_bp



logli_quic_bp<-logli_quic_ha_bp-logli_quic_ho_bp
logli_quic_bp




logli_quic_seq_bp[b]<-logli_quic_bp
}

mean(logli_quic<logli_quic_seq_bp)

logli_quic
logli_quic_seq_bp


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


L<-2

logrn_all<-rbind(logrn_middle,logrn_after)



nn2<-c(nrow(logrn_middle),nrow(logrn_after))
S_bar_m<-(cov(logrn_all,use="pairwise.complete.obs"))
cov_all<-S_bar_m
S_bar<-cbind(cov_middle,cov_after)

Lambda_quic_ho<-Lambda_tlp_ho<-mean(abs(diag(cov_all)))*sqrt(4*log(p)/sum(nn));
Lambda_quic_ha<-Lambda_tlp_ha<-c(mean(abs(diag(cov_middle)))*sqrt(4*log(p)/nn2[1]), mean(abs(diag(cov_after)))*sqrt(2*log(p)/nn2[2])); 
Lambda_quic_ho
Lambda_quic_ha

################################################################################




#lamatrix<-matrix(Lambda_quic_ho,p,p) ;  diag(lamatrix)<-0
#Omega1<-glasso(cov_all, rho=lamatrix,  msg=0)$X
Omega1<-glasso(cov_all, penalize.diagonal = FALSE, rho =Lambda_quic_ho)$wi
#Omega1
#Omega2<-(TLP_matrix(S_bar_m,Lambda=Lambda_tlp_ho,tau=1e-2,tol=1e-4,maxIter=1000, max.dc.iter=20,weighted=FALSE)$X.tlp)[,,1]


logli1<-rep(0,L); 

for (i in 1:L)
              {
              
              #i<-1
              S_bar1<-S_bar[,((i-1)*p+1):(i*p)];
              sum(S_bar1)
              
              
              logli1[i]<-nn[i]*log( det(Omega1))-nn[i]*sum(diag(S_bar1 %*% Omega1))


              }
             
logli_quic_ho<- sum(logli1)

logli_quic_ho

    
Omega1<-NULL; Omega2<-NULL
################################
################################


Omega1_seq<-list(L);
logli1<-rep(0,L); 
                  
for (i in 1:L)
            {
            S_bar1<-S_bar[,((i-1)*p+1):(i*p)]

            
            Omega1<-glasso(S_bar1, penalize.diagonal = FALSE, rho =Lambda_quic_ha[i])$wi
            
            Omega1_seq[[i]]<-Omega1                         

            
            logli1[i]<-nn[i]*log(det(Omega1))-nn[i]*sum(diag(S_bar1%*% Omega1))

            
            S_bar1<-NULL;  Omega1<-NULL;   Omega2<-NULL;  
            }
    
logli_quic_ha<- sum(logli1)

logli_quic<-logli_quic_ha-logli_quic_ho
logli_quic

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
logli_quic_seq_bp<-rep(0,ppmm)



for (b in 1:ppmm)


{


bp<-sample(1:sum(nn),sum(nn),replace=T)
logrn_all_bp<-logrn_all[bp, ]


logrn_middle_bp<-logrn_all_bp[1:nrow(logrn_middle),]
logrn_after_bp<-logrn_all_bp[(nrow(logrn_middle)+1):(nrow(logrn_middle)+nrow(logrn_after)),]


cov_middle_bp<-(cov(logrn_middle_bp,use="pairwise.complete.obs"))
cov_after_bp<-(cov(logrn_after_bp,use="pairwise.complete.obs"))

S_bar_m_bp<-cov(logrn_all_bp,use="pairwise.complete.obs") 
cov_all_bp<-S_bar_m_bp


Lambda_quic_ho<-mean(abs(diag(cov_all_bp)))*sqrt(4*log(p)/sum(nn));
Lambda_quic_ha<-c(mean(abs(diag(cov_middle_bp)))*sqrt(4*log(p)/nn2[1]), mean(abs(diag(cov_after_bp)))*sqrt(2*log(p)/nn2[2])); 
Lambda_quic_ho
Lambda_quic_ha




Omega1_bp<-glasso(cov_all_bp, penalize.diagonal = FALSE, rho =Lambda_quic_ho)$wi
#Omega1
#Omega2<-(TLP_matrix(S_bar_m,Lambda=Lambda_tlp_ho,tau=1e-2,tol=1e-4,maxIter=1000, max.dc.iter=20,weighted=FALSE)$X.tlp)[,,1]

S_bar_bp <-cbind(cov_middle_bp,cov_after_bp )


logli1_bp<-rep(0,L);     

for (i in 1:L)
              {

              S_bar1_bp<-S_bar_bp[,((i-1)*p+1):(i*p)];
              
              
              logli1_bp[i]<-nn[i]*log( det(Omega1_bp))-nn[i]*sum(diag(S_bar1_bp %*% Omega1_bp))


              }             
logli_quic_ho_bp<- sum(logli1_bp)

    
Omega1<-NULL; Omega2<-NULL



Omega1_seq_bp<-list(L);
logli1_bp<-rep(0,L); 
                  
for (i in 1:L)
            {
            #i<-3
            S_bar1_bp<-S_bar_bp[,((i-1)*p+1):(i*p)]

            
            Omega1_bp<-glasso(S_bar1_bp, penalize.diagonal = FALSE, rho =Lambda_quic_ha[i])$wi
            
            Omega1_seq_bp[[i]]<-Omega1_bp                         

            
            logli1_bp[i]<-nn[i]*log(det(Omega1_bp))-nn[i]*sum(diag(S_bar1_bp%*% Omega1_bp))

            }
    
logli_quic_ha_bp<- sum(logli1_bp)


logli_quic_ha_bp
#logli_tlp_ha

Omega1<-NULL; Omega2<-NULL



logli_quic_bp<-logli_quic_ha_bp-logli_quic_ho_bp
logli_quic_bp



logli_quic_seq_bp[b]<-logli_quic_bp
}



mean(logli_quic<logli_quic_seq_bp)
