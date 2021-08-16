for (zzz in 1:100){
                  unlink(".RData")
                  options(scipen=999)
                  vvv<-round(abs(rnorm(1)*100000000))
                  
                  library(glasso)
                  library(MASS)
                  
                  output<-NULL
                  for (p in c(5,10,20,30,40,50,100,200)) {
                  #p<-100 
                  #########################################################
                  #########################################################
                  
                  n <- 100;   L<- 4;  nn <- rep(n,L);  nns<-c(0,cumsum(nn))
                  rho.1<-0.5; rho.2<-0.5;  d1<-1; d2<-1;
                  ################ matrices generation ends ###############
                  
                  covmat_inverse1<-covmat_inverse2<-matrix(0,p,p)
                  
                  covmat_inverse1<-(rho.1)^abs(row(covmat_inverse1)-col(covmat_inverse1))
                  
                  
                  covmat_inverse2<-(rho.2)^abs(row(covmat_inverse1)-col(covmat_inverse1))
                  
                  covmat1 <- solve(covmat_inverse1)
                  covmat2 <- solve(covmat_inverse2)
                  
                  covmat<-vector("list",L)
                  covmat[[1]]<-covmat1; covmat[[3]]<-covmat1; covmat[[2]]<-covmat2; covmat[[4]]<-covmat2
                  ######	generating sample covariance matrices S_bar = [S_1, ... S_L]	###########
                  Lambda_quic_ho<-sqrt(log(p)/sum(nn))
                  Lambda_quic_ha<-rep(sqrt(log(p)/sum(nn)),L)
                  
                  S_ha <- vector("list",L); 
                  Sigma_ha<-vector("list",L); Omega_ha<-vector("list",L);      logli_ha<-rep(0,L); logli_ha_approx<-rep(0,L); 
                  Sigma_quic_ha<-vector("list",L); Omega_quic_ha<-vector("list",L); logli_quic_ha<-rep(0,L);   logli_quic_ha_approx<-rep(0,L);
                  
                  temp_matrix_t<-NULL
                  for (i in 1:L){
                                temp <- scale(mvrnorm(nn[i], rep(0,p), covmat[[i]]),scale=FALSE)
                                S_ha[[i]] <- t(temp)%*%(temp)/nn[i]
                                
                                temp_matrix_t<-rbind(temp_matrix_t,temp) 
                                mm<-glasso(S_ha[[i]], rho=matrix(Lambda_quic_ha[i],p,p))
                                Omega_quic_ha[[i]]<-mm$wi
                                Sigma_quic_ha[[i]]<-mm$w
                                
                                logli_quic_ha[i]<-nn[i]*log(det(Omega_quic_ha[[i]]))-nn[i]*sum(diag(S_ha[[i]]%*% Omega_quic_ha[[i]]))           
                                logli_quic_ha_approx[i]<- p-nn[i]*sum(diag(Sigma_quic_ha[[i]])) -nn[i]*sum(diag(S_ha[[i]]%*% Omega_quic_ha[[i]]))
                                
                                if (p<n) {Omega_ha[[i]]<-solve(S_ha[[i]]);Sigma_ha[[i]]<-S_ha[[i]]} else {Omega_ha[[i]]<-Omega_quic_ha[[i]];Sigma_ha[[i]]<-Sigma_quic_ha[[i]] }
                                
                                logli_ha[i]<-nn[i]*log(det(Omega_ha[[i]]))-nn[i]*sum(diag(S_ha[[i]]%*% Omega_ha[[i]]))
                                logli_ha_approx[i]<- p-nn[i]*sum(diag(Sigma_ha[[i]]))-nn[i]*sum(diag(S_ha[[i]]%*% Omega_ha[[i]]))
                                }
                  
                  S_ho<-t(temp_matrix_t)%*%(temp_matrix_t)/sum(nn) 
                  mmo<-glasso(S_ho, rho=matrix(Lambda_quic_ho,p,p))
                  Omega_quic_ho<-mmo$wi
                  Sigma_quic_ho<-mmo$w
                  
                  logli_quic_ho<-sum(nn)*log(det(Omega_quic_ho))-sum(nn)*sum(diag(S_ho%*%Omega_quic_ho))             
                  logli_quic_ho_approx<-p- sum(nn)*sum(diag(Sigma_quic_ho))-sum(nn)*sum(diag(S_ho%*%Omega_quic_ho))
                  
                  if (p<sum(nn)) {Omega_ho<-solve(S_ho); Sigma_ho<-(S_ho)} else {Omega_ho<-Omega_quic_ho; Sigma_ho<-Sigma_quic_ho}
                  logli_ho<- sum(nn)*log(det(Omega_ho))-sum(nn)*sum(diag(S_ho%*%Omega_ho))
                  
                  logli<-sum(logli_ha)-logli_ho
                  logli_quic<-sum(logli_quic_ha)-logli_quic_ho
                  logli_quic_approx<-sum(logli_quic_ha_approx)-logli_quic_ho_approx
                  ################################################################################
                  ################################################################################
                  ################################################################################
                  
                  ppmm<-sum(nn);
                  logli_quic_bt<- rep(0,ppmm);
                  logli_quic_bt_approx<- rep(0,ppmm);
                  
                  
                  for (dp in 1:ppmm){
                                    btindex<-sample(1:(sum(nn)), sum(nn), replace=TRUE); temp_matrix_bt<-temp_matrix_t[btindex,]                  
                                    #temp_matrix_bt<-temp_matrix_t*rnorm(sum(nn))#sample(c(-1,1),sum(nn),replace=T)
                                    
                                    S_ho_bt<-t(temp_matrix_bt)%*%(temp_matrix_bt)/sum(nn)
                                    mmbt<-glasso(S_ho_bt, rho=matrix(Lambda_quic_ho,p,p))
                                    Omega_quic_ho_bt<-mmbt$wi
                                    Sigma_quic_ho_bt<-mmbt$w
                                    
                                    
                                    logli_quic_ho_bt<- sum(nn)*log(det(Omega_quic_ho_bt))-sum(nn)*sum(diag(S_ho_bt%*% Omega_quic_ho_bt))
                                    logli_quic_ho_bt_approx<- p-sum(nn)*sum(diag(Sigma_quic_ho_bt))-sum(nn)*sum(diag(S_ho_bt%*% Omega_quic_ho_bt))
                                    
                                    #logli_quic_ho_bt<- sum(nn)*log(det(Omega_quic_ho))-sum(nn)*sum(diag(S_ho_bt%*% Omega_quic_ho))
                                    logli_quic_ha_bt<-rep(0,L); logli_quic_ha_bt_approx<-rep(0,L);  Omega_quic_ha_bt<-vector("list",L); Sigma_quic_ha_bt<-vector("list",L)
                                    
                                    for (i in 1:L){
                                                  temp_bt<-temp_matrix_bt[(nns[i]+1):nns[i+1], ];
                                                  S_ha_bt<-t(temp_bt)%*%(temp_bt)/nn[i];
                                                  mmbt<-glasso(S_ha_bt, rho=matrix(Lambda_quic_ha[i],p,p))
                                                  Omega_quic_ha_bt[[i]]<-mmbt$wi
                                                  Sigma_quic_ha_bt[[i]]<-mmbt$w
                                                  
                                                  logli_quic_ha_bt[i]<-nn[i]*log(det(Omega_quic_ha_bt[[i]]))-nn[i]*sum(diag(S_ha_bt%*% Omega_quic_ha_bt[[i]]))                                                                 
                                                  logli_quic_ha_bt_approx[i]<-p-nn[i]*sum(diag(Sigma_quic_ha_bt[[i]]))-nn[i]*sum(diag(S_ha_bt%*% Omega_quic_ha_bt[[i]]))
                               
                                                  }                  
                                    logli_quic_bt[dp]<-sum(logli_quic_ha_bt)-(logli_quic_ho_bt)
                                    logli_quic_bt_approx[dp]<-sum(logli_quic_ha_bt_approx)-(logli_quic_ho_bt_approx)
                                    }
                  pv0<-pchisq(logli_quic, (p^2+p)*(L-1)/2, ncp = 0, lower.tail = F, log.p = FALSE)                  
                  pv1<-pchisq(logli, (p^2+p)*(L-1)/2, ncp = 0, lower.tail = F, log.p = FALSE)
                  pv2<-mean(logli_quic_bt>logli_quic)
                  pv3<-mean(logli_quic_bt_approx>logli_quic_approx)
                                    
                  output<-rbind(output,round(c(n,p, L, rho.1,rho.2,d1,d2,logli_quic,pv0,logli,pv1,logli_quic,pv2,pv3),4))
                  print(output)                 
                  }
                  
                  filename<-paste("C:/Users/alvin/Dropbox/Graphical/sim/global1/set2/",vvv,".txt",sep="");
                  write.table(output, file=filename, row.names=FALSE,col.names=FALSE)
                  } 