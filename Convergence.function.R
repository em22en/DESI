#rm(list=ls(all=TRUE))
#Libraries
library(ggplot2)
library(reshape2)
library(plm)
library(car)
library(lmtest)
#########################################

##########################################
#     Convergence function	             #
##########################################
G=function(a)
	{
a=na.omit(a)
a_temp=as.matrix(a)
Countries=rownames(a)
periods=colnames(a)
T=length(a[1,])
C=length(a[,1])

meanA=matrix(rep(0),C,1)
	rownames(meanA)=Countries
	colnames(meanA)=c("mean")
	for(h in 1:C){meanA[h,1]=mean(a_temp[h,])}

varT=matrix(rep(0),T,1)
	rownames(varT)=periods
	colnames(varT)=c("Variance")
	for(k in 1:T){varT[k,1]=var(a_temp[,k])}

varMin=min(varT)
varMax=max(varT)
minV=rownames(varT)[varMin==varT]
maxV=rownames(varT)[varMax==varT]
varStat=list()
	varStat$max=c(maxV,round(varMax,4))
	varStat$min=c(minV,round(varMin,4))

Growth=matrix(rep(0),C,T-1)
	rownames(Growth)=Countries
	colnames(Growth)=periods[2:T]
	for(j in 1:C){
		for (i in 1:(T-1)){
			Growth[j,i]=(a[j,i+1]-a[j,i])/a[j,i]
					}
		 	}
meanG=matrix(rep(0),C,1)
	rownames(meanG)=Countries
	colnames(meanG)=c("Av_Growth")
for(l in 1:C){meanG[l,1]=mean(Growth[l,])}

#-----------FE regression-----------------#
NOM=substring(colnames(Growth),3)
NOM=as.numeric(NOM)
years=rep(0,C*(T-1))
	for(t in 1:(T-1))
		{
		init=((C*t)-(C-1))
		fin=(C*t)
		years[init:fin]=NOM[t]
		}
original=a_temp
original=original[,-T]
original_m=melt(original)

gdata_m=melt(Growth)
gdata_m=gdata_m[,-2]
gdata_m=cbind(gdata_m,years)

panel=cbind(gdata_m,original_m[,3])
colnames(panel)=c("country","growth","years","Y0")

fixed =plm(growth~Y0,data=panel,index=c("country","years"), model="within")
FE=summary(fixed)
FE_int=fixef(fixed)   
DW_FE=pdwtest(fixed)#Durbin Watson test serial correlation
WR=pwartest(fixed) #Wooldridge test serial correlation

#---------OLS Growth Equation-------------#
Y0=log(a_temp[,1])
Y0_2=(Y0)^2
YT=log(a_temp[,T])
Ydiff=(YT-Y0)/T
fit=lm(Ydiff ~ Y0)
fit2=lm(Ydiff ~ Y0+Y0_2)
GEq=summary(fit)
GEq2=summary(fit2)
#------------Plot Av Growth vs Y0---------#
initial=a[,1]
Yt=cbind(initial,meanG)
Yt=as.data.frame(Yt)
plot1=ggplot(Yt, aes(x=initial, y=Av_Growth))+
	geom_point(size=2,colour="blue")+
	geom_text(colour="red",label=Countries,size=3,hjust=-0.5,vjust=0)+
	labs(title = "Average growth vs Initial conditions")+
	geom_smooth(method='lm')
#------ln--Plot Total Growth vs Y0---------#
Zt=cbind(a_temp[,1],((a_temp[,T]-a_temp[,1])/a_temp[,1])*100) #Unidades normales para plot
Zt=as.data.frame(Zt)
	colnames(Zt)=c("Y0","Growth")
plot2=ggplot(Zt, aes(x=Y0, y=Growth))+
	geom_point(size=2,colour="blue")+
	geom_text(colour="red",label=Countries,size=3,hjust=-0.5,vjust=0)+
	labs(title = "Total Growth",y="Total Growth %")+
	geom_smooth(method='lm')
#----------Growth-plot---------------------#
colnames(gdata_m)=c("Country","Growth","time")
plot3=ggplot(gdata_m, aes(x=time, y=Growth))+
	geom_line(aes(color=Country),size=1)+
		theme(axis.text=element_text(size=10),
		axis.text.x = element_text(angle=90),
            axis.title=element_text(size=10,face="bold"),
		legend.position = "bottom",
		legend.key = element_rect(fill = "white",colour = "white")  
	 	)+
	guides(col = guide_legend(nrow = 3))+
	labs(title = "Growth over the time")
#----------------Plot-Variable------------------#
NOM1=substring(colnames(a_temp),3)
NOM1=as.numeric(NOM1)
years1=rep(0,(C*T))
	for(t in 1:(T))
		{
		init=((C*t)-(C-1))
		fin=(C*t)
		years1[init:fin]=NOM1[t]
		}
a_melt=melt(a_temp)
a_melt[,2]=years1
	colnames(a_melt)=c("Country","time","Value")

plot4=ggplot(a_melt, aes(x=time, y=Value))+
	geom_line(aes(color=Country),size=1)+
		theme(axis.text=element_text(size=10),
		axis.text.x = element_text(angle=90),
            axis.title=element_text(size=10,face="bold"),
		legend.position = "bottom",
		legend.key = element_rect(fill = "white",colour = "white")  
	 	)+
	guides(col = guide_legend(nrow = 3))+
	labs(title = "Time perspective")
#---------Variance Evolution Plot---------------#
varP=cbind(sqrt(varT),NOM1)
varP=as.data.frame(varP)
	names(varP)=c("Std_Dev","time")

plot5=ggplot(varP, aes(x=time, y=Std_Dev))+
	geom_line(size=1)+
	labs(title = "Standard Deviation Evolution")

#-----------Speed of Catching up----------------#
last=a_temp[,T]
nt=cbind(last,meanG)
nt=as.data.frame(nt)
ctarget_name=rownames(nt)[nt$last==max(nt$last)]
ctarget=max(nt$last)
gtarget=nt$Av_Growth[rownames(nt)==ctarget_name]
best=c(ctarget_name,round(gtarget,6),round(ctarget,6))
	names(best)=c("Country","Avg.Growth","Last term")

n=rep(0,C)
nt=cbind(nt,n)
nt$n=(log(ctarget)-log(nt$last))/(log(1+nt$Av_Growth)-log(1+gtarget))
nt= round(nt,3)
nt=nt[order(-nt[,3]), , drop = FALSE]

tmean=mean(nt$last)
tmean=round(tmean,4)
nmean=cbind(nt[,1:2],n)
nmean$n=(log(tmean)-log(nt$last))/(log(1+nt$Av_Growth)-log(1))
nmean=nmean[order(-nmean[,3]), , drop = FALSE]
nmean=subset(nmean, last <= tmean)
nmean=round(nmean,4)
#-----------Objects-----------------------#
return(list(a,			#1 clean BD as data.frame
		Growth,	 	#2 Growth of Country per period
		meanG,	 	#3 Mean Growth over time
		meanA,	 	#4 Mean of Country over the time
		varT,			#5 Year Variance (sigma convergence)
		varStat,		#6 Stats on variance evolution
		Countries,		#7 Names of countries
		periods,		#8 time periods
		T,		 	#9 Number of time periods
		C,		 	#10 Number of Countries
		a_temp,	 	#11 Clean BD as matrix (for calculations)
		panel,	 	#12 Panel data for FE regression
		fixed,	 	#13 FE regression
		FE,		 	#14 Summary FE regression
		FE_int,	 	#15 Intercepts FE regression
		DW_FE,	 	#16 Durbin Watson serial correlation test
		fit,			#17 LS Growth Equation fit
		GEq,			#18 Summary Growth Equation fit
		plot1,	 	#19 Plot Average Growth vs Y0
		plot2,	 	#20 Plot Growth equation
		plot3,	 	#21 Plot of Growth over the time
		plot4,	 	#22 Plot variable over the time
		plot5,	 	#23 Plot Variance evolution over the time
		best,			#24 Best country to immitate Catch-up
		nt,		 	#25 Catch up speed to best performance country
		tmean,		#26 Mean of last period
		nmean,		#27 Catch up to the mean
		gdata_m,	 	#28 Melted growth data
		GEq2			#29 Summary OLS quadratic model
)
	)
	}
#----------------------------------------#