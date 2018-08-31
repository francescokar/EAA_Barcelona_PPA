## Import dataset

library(GISTools)
setwd("D:/R-Code_Projects/EAA_Barcelona")
FB018B<-readShapePoly("FB018B_Surface")
RR_FB018B<-readShapePoints("RR_FB018B")

## Assess size and preservation proprtion of debris

max(RR_FB018B$Dimension)

tiff("Hist_Size_Debris.tif",width=4800,height=3600,units="px",res=800)
par(mar=c(4.5,4.5,3,1))
plot(hist(RR_FB018B$Dimension,breaks=30),col="#999999",xlim=c(0,800),xlab="size",
     ylab="frequency",main="Size of debris in FB018B")
dev.off()

min(log(RR_FB018B$Dimension))
max(log(RR_FB018B$Dimension))

tiff("Hist_LogSize_Debris.tif",width=4800,height=3600,units="px",res=800)
plot(hist(log(RR_FB018B$Dimension),breaks=30),col="#666666",xlim=c(-2,7),
     xlab="log-size",ylab="frequency",main="Log-size of debris in FB018B")
dev.off()

## Create plot

summary(RR_FB018B$Dimension)

tiff("RR_FB018B.tif",width=6000,height=4800,units="px",res=800)
par(mar=c(0,0,4,0))
plot(FB018B,lwd=1.8,col="#CCCCCC",main="FB018B: intra-site debris distribution
     (log-size plotted)",cex.main=1.2)
plot(RR_FB018B,cex=log(RR_FB018B$Dimension),pch=c(1,2)[RR_FB018B$Preservati],
     lwd=1,col="red",add=T)
legend(-300,530,legend=c("Complete","Fragment"),pch=c(1,2),cex=1.2,col="red")
map.scale(200,70,len=200,units="meters",ndivs=2,subdiv=1)
dev.off()


### POINT PATTERN ANALYSIS

library(spatstat)

# Global Homogeneous L-Function

FB018Bowin<-as.owin(FB018B)
RR_FB018B_ppp<-as.ppp(coordinates(RR_FB018B),FB018Bowin)
RR_FB018B_Lest<-envelope(RR_FB018B_ppp,fun=Lest,nsim=9999,correction="Ripley")

tiff("L-Function_Plot.tif",width=4800,height=4800,units="px",res=800)
par(mar=c(5,5,5,3))
plot(RR_FB018B_Lest,main="FB018B",xlab="distance(cm)")
dev.off()

# Local L-Function

LocalL5<-localL(RR_FB018B_ppp,correction="Ripley",rvalue=5)
LocalL15<-localL(RR_FB018B_ppp,correction="Ripley",rvalue=15)
LocalL30<-localL(RR_FB018B_ppp,correction="Ripley",rvalue=30)
localL40<-localL(RR_FB018B_ppp,correction="Ripley",rvalue=40)
RR_FB018B_localL<-list("r5"=LocalL5,"r15"=LocalL15,"r30"=LocalL30,
                       "r40"=localL40)

# Local L maps

tiff("RR_FB018B_localL5.tif",width=6000,height=4800,units="px",res=800)
par(mar=c(0,0,4,0))
plot(FB018B,lwd=1.8,col="#CCCCCC",main="FB018B: local spatial autocorrelation
     (r = 5)",cex.main=1.2)
plot(RR_FB018B,cex=(RR_FB018B_localL$r5)*0.1,pch=1,
     lwd=1,col="red",add=T)
map.scale(200,70,len=200,units="meters",ndivs=2,subdiv=1)
dev.off()

tiff("RR_FB018B_localL15.tif",width=6000,height=4800,units="px",res=800)
par(mar=c(0,0,4,0))
plot(FB018B,lwd=1.8,col="#CCCCCC",main="FB018B: local spatial autocorrelation
     (r = 15)",cex.main=1.2)
plot(RR_FB018B,cex=(RR_FB018B_localL$r15)*0.1,pch=1,
     lwd=1,col="blue",add=T)
map.scale(200,70,len=200,units="meters",ndivs=2,subdiv=1)
dev.off()

tiff("RR_FB018B_localL30.tif",width=6000,height=4800,units="px",res=800)
par(mar=c(0,0,4,0))
plot(FB018B,lwd=1.8,col="#CCCCCC",main="FB018B: local spatial autocorrelation
     (r = 30)",cex.main=1.2)
plot(RR_FB018B,cex=(RR_FB018B_localL$r30)*0.1,pch=1,
     lwd=1,col="green4",add=T)
map.scale(200,70,len=200,units="meters",ndivs=2,subdiv=1)
dev.off()

tiff("RR_FB018B_localL40.tif",width=6000,height=4800,units="px",res=800)
par(mar=c(0,0,4,0))
plot(FB018B,lwd=1.8,col="#CCCCCC",main="FB018B: local spatial autocorrelation
     (r = 40)",cex.main=1.2)
plot(RR_FB018B,cex=(RR_FB018B_localL$r40)*0.1,pch=1,
     lwd=1,col="purple",add=T)
map.scale(200,70,len=200,units="meters",ndivs=2,subdiv=1)
dev.off()

# Inhomogeneous Possion model

X<-c(-180.42)
Y<-c(184.14)
Fire_FB018B<-ppp(X,Y,window=FB018Bowin)
Dist_Fire_FB018B<-distmap(Fire_FB018B)
FB018B_ppm<-ppm(unmark(RR_FB018B_ppp),~D,covariates=list(D=Dist_Fire_FB018B))


# Testing Inhomogeneity

library("MuMIn")
FB018B_null<-ppm(unmark(RR_FB018B_ppp),~1)
AICc(FB018B_ppm,FB018B_null)
Weights(AICc(FB018B_ppm,FB018B_null))


# L-Function conditioned on the Inhomogeneous Poisson model

RR_FB018B_Linhom<-envelope(FB018B_ppm,Lest,nsim=9999,correction="Ripley")

tiff("Linhom_Plot.tif",width=4800,height=4800,units="px",res=800)
par(mar=c(5,5,5,3))
plot(RR_FB018B_Linhom,main="FB018B: Inhomogeneous Poisson",xlab="distance(cm)")
dev.off()


### SPATIAL DEPENDENCY

library(pgirmess)
library(spdep)

# Moran's I correlogram (distance classes and lags)

coord<-subset(coordinates(RR_FB018B),select=-c(coords.x3))
# for some reason (Autocad?) my shp contained a third coordinate...

### Method used in Carrer 2017 to create correlogram for distance classes
#RR_FB018B_correl<-correlog(coords=coord,z=log(RR_FB018B$Dimension),
#                           method="Moran")
#RR_FB018B_correl_tab<-as.data.frame(RR_FB018B_correl,row.names=NULL)
#tiff("Correlogram_pgirmess.tif",width=6000,height=4800,units="px",res=800)
#plot(RR_FB018B_correl_tab$dist.class,
#     RR_FB018B_correl_tab$coef,type="b",pch=0,cex=2,
#     xlim=c(0,290),ylim=c(-0.12,0.27),xlab="Distance (cm)",
#     ylab="Moran's I",main="FB018B: Moran's I Correlogram")
#par(new=T)
#plot(RR_FB018B_correl_tab$dist.class[2],RR_FB018B_correl_tab$coef[2],
#     pch=15,cex=2,xlim=c(0,290),ylim=c(-0.12,0.27),ann=F)
#dev.off()

RRnb<-knn2nb(knearneigh(coord,k=8,RANN=F))

tiff("RR_FB018B_Neighbours.tif",width=6000,height=4800,units="px",res=800)
par(mar=c(0,0,4,0))
plot(FB018B,lwd=1.8,col="#CCCCCC",
     main="Neighbour relationships between RR (k = 8)",cex.main=1.2)
plot(RRnb, coordinates(RR_FB018B),pch=19,cex=0.6,add=T)
map.scale(200,70,len=200,units="meters",ndivs=2,subdiv=1)
dev.off()

RR_FB018B_correl_knn<-sp.correlogram(neighbours=RRnb,var=log(RR_FB018B$Dimension),
                         order=5,method="I",style="W")

tiff("Correlogram_spdep.tif",width=6000,height=4800,units="px",res=800)
par(mar=c(5,5,5,3))
par(cex=1.3)
par(cex.main=1.2)
plot(RR_FB018B_correl_knn,main="Moran's I correlogram log-size of RR")
dev.off()


# Local Moran's I test of RR

RR_FB018B_lMoran<-localmoran(log(RR_FB018B$Dimension),
                             listw=nb2listw(RRnb,style="W"))
RRlm_I1<-(RR_FB018B_lMoran[,1]<=-1)+0
RRlm_I2<-(RR_FB018B_lMoran[,1]>-1 & RR_FB018B_lMoran[,1]<=0)*2
RRlm_I3<-(RR_FB018B_lMoran[,1]>0 & RR_FB018B_lMoran[,1]<=1)*3
RRlm_I4<-(RR_FB018B_lMoran[,1]>1)*4
RRlm_Iclass<-RRlm_I1+RRlm_I2+RRlm_I3+RRlm_I4
RRlm_pch<-c(15,17,18,19)
RRlm_pr<-(RR_FB018B_lMoran[,5]<=0.05)+0
RRlm_pr2<-(RR_FB018B_lMoran[,5]>0.05)*2
RRlm_Prclass<-RRlm_pr+RRlm_pr2
RRlm_col<-c("#31A354","#DE2D26")

tiff("Localmoran.tif",width=6000,height=4800,units="px",res=800)
par(mar=c(0,2,5,0))
plot(FB018B,lwd=1.9,col="#CCCCCC",main="FB018B: Local-Moran's I",cex.main=1.5)
plot(RR_FB018B,cex=1.5,pch=RRlm_pch[RRlm_Iclass],
     col=RRlm_col[RRlm_Prclass],add=T)
legend(-300,580,legend=c("Pr < 0.05","Pr > 0.05"),pch=15,pt.cex=1.2,
       col=c("#31A354","#DE2D26"))
legend(-300,500,legend=c("I < -1","I > -1 < 0","I > 0 < 1","I > 1"),pt.cex=1.2,
       pch=RRlm_pch)
map.scale(200,70,len=200,units="meters",ndivs=2,subdiv=1)
dev.off()

### LOCAL JOIN COUNT STATISTICS

# Provisional example

RR_FB018B_cat<-as.numeric(RR_FB018B$Preservati)-1

for(i in 1:length(RR_FB018B_cat)){
  xi<-(RR_FB018B_cat[1:i])
  xj<-RR_FB018B_cat
  w<-(RRweight[,1:i])
  bb<-w*xj
  ww<-w*(1-xj)}

sbb<-colSums(bb)
BB_RRLICD<-RR_FB018B_cat*sbb
sww<-colSums(ww)
WW_RRLICD<-(1-RR_FB018B_cat)*sww

summary(BB_RRLICD)
summary(WW_RRLICD)

# plot the LICD results
RR_BB<-(BB_RRLICD>=0.5032)*3
RR_WW<-(WW_RRLICD>=0.3614)*2
RR_BW<-(BB_RRLICD<0.5032 & WW_RRLICD<0.3614)+0
RR_BBWW<-RR_BB+RR_WW+RR_BW
RRJCS_pch<-c(4,16,1)
tiff("FB018B_LICD.tif",width=6000,height=4800,units="px",res=800)
par(mar=c(0,2,5,0))
plot(FB018B,lwd=1.8,col="#CCCCCC",main="FB018B: Local-JCS",cex.main=1.5)
plot(RR_FB018B,cex=1.5,pch=RRJCS_pch[RR_BBWW],add=T)
legend(-300,530,legend=c("Not significant","BB","WW"),pt.cex=1.2,
       pch=RRJCS_pch)
map.scale(200,70,len=200,units="meters",ndivs=2,subdiv=1)
dev.off()

# compare with Complete/Fragment distribution
tiff("FB018B_Preservation.tif",width=6000,height=4800,units="px",res=800)
par(mar=c(0,2,5,0))
plot(FB018B,lwd=1.8,col="#CCCCCC",main="FB018B: Preservation",cex.main=1.5)
plot(RR_FB018B,cex=1.5,pch=c(16,1)[RR_FB018B$Preservati],add=T)
legend(-300,530,legend=c("Complete","Fragment"),pt.cex=1.2,
       pch=c(16,1))
map.scale(200,70,len=200,units="meters",ndivs=2,subdiv=1)
dev.off()