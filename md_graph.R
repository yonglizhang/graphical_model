
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
md<-read.csv("C:/Users/Yongli Zhang/Documents/Graphical/finnetwork/md.csv", header = TRUE)
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
Lambda_quic_ha<-c( mean(diag(cov_before))*sqrt( 4*log(p)/nn[1]),mean(diag(cov_middle))*sqrt(4*log(p)/nn[2]),mean(diag(cov_after))*sqrt(4*log(p)/nn[3])  ); 

Lambda_quic_ho
Lambda_quic_ha


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

Omega1_ho<-glasso(cov_all, penalize.diagonal = FALSE, rho =Lambda_quic_ho)$wi

logli1_ho<-rep(0,L); 

for (i in 1:L)
              {
              
              #i<-1
              S_bar1<-S_bar[,((i-1)*p+1):(i*p)];
              sum(S_bar1)
              
              
              logli1_ho[i]<-nn[i]*log( det(Omega1_ho))-nn[i]*sum(diag(S_bar1 %*% Omega1_ho))


              }
             
logli_quic_ho<- sum(logli1_ho)

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
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

colorlist<-c( "black" ,  "red", "green" , "blue",  "cyan" ,   "magenta", "yellow" ,"gray"   ,"orange","brown")
cr<-c("Pre-Crisis (1/1/2005-12/31/2005)","Crisis (7/1/2008-6/30/2009)","Post-Crisis (1/1/2010-12/31/2010)")





dev.off()
dev.off()
dev.off()
postscript(file=paste("C:/Users/yongli/Dropbox/graphical/partial_cor.ps",sep=""), width = 11, height = 8,
                     horizontal = T, paper = "special")




set.seed(14)
par(mfrow=c(1,3))
for (i in 1:3)
{

#i<-2
print(sum(abs((Omega1_ha_seq[[i]]))))

#g <- graph.adjacency(abs((Omega1_ha_seq[[i]]))!=0,mode="undirected",diag=F,weighted=T)
g <- graph.adjacency(abs((pcor_ha_seq[[i]]))>0.09,mode="undirected",diag=F,weighted=T)

E(g)$color<-"black"
V(g)$color<-colorlist[ceiling(as.numeric(V(g)) /20)]
V(g)$size<-4
plot.igraph(g, vertex.label=NA,edge.width=1, layout=layout.fruchterman.reingold(g, niter=10000))
 title(cr[i])
}


dev.off()
dev.off()
dev.off()





################################################################################
################################################################################

dev.off()







postscript(file=paste("C:/Users/yongli/Dropbox/graphical/partial_cor_finance.ps",sep=""), width = 11, height = 4,
                     horizontal = T, paper = "special",fonts=c("serif", "Palatino"))               
                     
set.seed(14)
xindex<-101:120
par(mfrow=c(1,3))
for (i in 1:3)
{

#i<-2

sector_cor<-abs(cov2cor(pcor_ha_seq[[i]]))[xindex,xindex]


print(sum(abs(cov2cor(Omega1_ha_seq[[i]]))))
g <- graph.adjacency(sector_cor>0.09,mode="undirected",diag=F,weighted=T)
E(g)$color<-"black"

#V(g)$size<-8/diag(sector_cor)
#V(g)$color<-colorlist[ceiling(as.numeric(V(g)) /20)]
V(g)$color<-"white"
#V(g)$label<-stocknames[xindex]


plot.igraph(g,vertex.label=stocknames[xindex],edge.width=1,vertex.label.family="serif",layout=layout.fruchterman.reingold(g, niter=10000))
title(cr[i])

}


dev.off()
dev.off()



################################################################################
################################################################################






postscript(file=paste("C:/Users/yongli/Dropbox/graphical/partial_cor_energy.ps",sep=""), width = 11, height = 4,
                     horizontal = T, paper = "special",fonts=c("serif", "Palatino"))               

xindex<-81:100
#xindex<-181:200                     
set.seed(14)

par(mfrow=c(1,3))
for (i in 1:3)
{

#i<-2

sector_cor<-abs(cov2cor(pcor_ha_seq[[i]]))[xindex,xindex]


print(sum(abs(cov2cor(Omega1_ha_seq[[i]]))))
g <- graph.adjacency(sector_cor>0.05,mode="undirected",diag=F,weighted=T)
E(g)$color<-"black"

#V(g)$size<-8/diag(sector_cor)
#V(g)$color<-colorlist[ceiling(as.numeric(V(g)) /20)]
V(g)$color<-"white"
#V(g)$label<-stocknames[xindex]


plot.igraph(g,vertex.label=stocknames[xindex],edge.width=1,vertex.label.family="serif",layout=layout.fruchterman.reingold(g, niter=10000))
title(cr[i])

}




dev.off()
dev.off()
dev.off()