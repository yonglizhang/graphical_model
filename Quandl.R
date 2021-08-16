
stocklist <- NULL

#data_files_2 <- list.files(path = "C:/Users/yongli/Dropbox/Graphical/realdata/joint_analysis_data/data_2")
data_files_2 <- list.files(path = "C:/Users/Yongli Zhang/Documents/Graphical/finnetwork/co/")
for (file in data_files_2){
                          #file_path <- paste("C:/Users/yongli/Dropbox/Graphical/realdata/joint_analysis_data/data_2/",file,sep="")
                        	file_path <- paste("C:/Users/Yongli Zhang/Documents/Graphical/finnetwork/co/",file,sep="")
                        	
                        	stocklist<-cbind(stocklist,as.character(read.csv(file_path)[,"Symbol"] [1:100]))
                        	#data_cd3cd28icam2 <- rbind(data_cd3cd28icam2,scale(log(read.csv(file_path))))
                          }	

stocklist


library(Quandl)
Quandl.api_key("6KxSEL4WfTnTSzg45Uty")


##mydata = Quandl("NSE/OIL", start_date="2012-01-01", end_date="2015-12-31")
##mydata[1,]
##dim(mydata)


#stocklist<-read.delim("C:/Users/yongli/Dropbox/Graphical/finnetwork/top15stocks.txt", header = FALSE)#[,1:10]
#dim(stocklist)
#stocklist
 md<-NULL
 st<-"2004-01-01"
 ed<-"2015-12-31"
 Date<- seq(as.Date(st),as.Date(ed),by=1)
 md<-data.frame(Date)
 md
 dim(md)
 
 
 
for (j in 1:10)  {


                  #j<-1
                  k<-0
                  for (i in 1:100)  {
                                    #i<-2
                                    a<-stocklist[i,j]
                                    print(a)
                                    mydata<-NULL
                                    
                                    x.inv <- try(mydata <- Quandl(paste("WIKI/",a,sep=""), start_date=as.Date(st), end_date=as.Date(ed)))
                                    if ('try-error' %in% class(x.inv))     {print(i);next}  
                                    
                                    print(c(k, j, i,nrow(mydata)))
                                    
                                    myt<-mydata[,c("Date","Adj. Close")]
                                    names(myt)<-c("Date",a)
                                    if (nrow(mydata)>=3000)  {
                                    k<-k+1; 
                                    md<-merge(md, myt,by.x="Date",by.y="Date",all.x=TRUE) }
                                    k
                                    
                                    if (k==20) {break() }
                                
                                    }
                                    print(dim(md))
                                    print(md[1,])
                                    print(apply(is.na(md),2,sum))
}

dim(md)

mdna<-md[apply(is.na(md),1,sum) <1,]
 
dim(mdna)
write.csv(mdna, file = "C:/Users/yongli/Dropbox/Graphical/finnetwork/md.csv",row.names = F)











################################################################################
mds<-read.csv("C:/Users/yongli/Dropbox/Graphical/finnetwork/md.csv", header = TRUE)
dim(mds)
mds[1,]

mds[1:10,]

mds_before<-subset(mds,as.Date(Date)<as.Date("2008-09-14") )
dim(mds_before)
mds_mid<-subset(mds,(as.Date(Date)>=as.Date("2008-09-15"))&as.Date(Date)<=as.Date("2009-03-06") )
mds_after<-subset(mds,as.Date(Date)>as.Date("2009-03-06") )



