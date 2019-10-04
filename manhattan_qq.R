
args<-commandArgs()

# insert data file here
datafile<-args[4]

data<-read.table(datafile,header=T,strings=F,fill=T)


####
####	MANHATTAN mtag_pvalLOT

tiff(paste(args[4],".tiff",sep=""),2*2**10,2**10,res=300,compression="lzw",pointsize=7)
a<-NULL
labpos<-NULL
pos<-NULL

# set maximise base positions for each chromosome
a[1]<-233627066;a[2]<-243188919;a[3]<-197946621;a[4]<-191043593;a[5]<-180876273;a[6]<-171051269;a[7]<-159128574;a[8]<-146303866;a[9]<-141111026;a[10]<-135523864;a[11]<-134946451;a[12]<-133841510;a[13]<-115109852;a[14]<-107289436;a[15]<-102520751;a[16]<-90292811;a[17]<-81194907;a[18]<-78017157;a[19]<-59118838;a[20]<-62965028;a[21]<-48118513;a[22]<-51243297

#generate plot positions
cum_a<-NULL
cum_b<-NULL
for(i in 1:22){cum_a[i]<-sum(a[1:i])}
for(i in 1:22){if(i>1){cum_b[i]<-cum_a[i-1]}else cum_b[i]<-0}
data$CHR.MAX<-cum_b[data$CHR]
data$plot.pos<-data$CHR.MAX+as.numeric(data$BP)
plot.max<-max(-log10(.001),max(-log10(data$mtag_pval))+0.25)

#set colours
#colours<-colours<-rep(c("darkblue","darkred","darkolivegreen","cornflowerblue","darkgoldenrod1"),6)
colours<-rep(c("tomato","cornflowerblue"),11)

# make plot, with different colour for each chromosome
for(j in 1:22){
	if(j==22){
		plot(frame=F,data$plot.pos[data$CHR==j & data$mtag_pval<0.05],-log10(data$mtag_pval[data$CHR==j & data$mtag_pval<0.05]),ylim=c(2,plot.max),col=colours[j],ylab="-log10(p-value)",xlim=c(0,max(data$plot.pos)),pch=20,xaxt='n',xlab="Chromosomes",cex=1)
		par(new=T)
	}
	else { 
		plot(frame=F,data$plot.pos[data$CHR==j & data$mtag_pval<0.05],-log10(data$mtag_pval[data$CHR==j & data$mtag_pval<0.05]),ylim=c(2,plot.max),col=colours[j],ylab="",xlim=c(0,max(data$plot.pos)),pch=20,xaxt='n',xlab="",cex=1,yaxt='n')
		par(new=T)
	}
}


# add chromosome identifier
for(i in 1:22){pos[i]<-(cum_a[i]+cum_b[i])/2;mtext(i,1,at=pos[i],cex=1,line=0)}

# plot genome-wide and suggestive significance lines
abline(h=-log10(5e-8),col="seashell3",lty=2,cex=.75)

dev.off()



####
####	QQ-mtag_pvalLOT
## obs <- readfile; p-values only
## read in your p-values,
## here I generated some
tiff(paste(args[4],"_QQ.tiff",sep=""),2**10,2**10,res=200,compression="lzw",pointsize=7)
obs <- -log10(data$mtag_pval)
N <- length(obs) ## number of p-values

## create the null distribution
## (-log10 of the uniform)
null <- -log10( ppoints(length(obs) ))
MAX <- max(c(obs,null))


## create the confidence intervals
c95 <- rep(0,N)
c05 <- rep(0,N)

## the jth order statistic from a
## uniform(0,1) sample
## has a beta(j,n-j+1) distribution
## (Casella & Berger, 2002,
## 2nd edition, pg 230, Duxbury)

for(i in 1:N){
c95[i] <- qbeta(0.95,i,N-i+1)
c05[i] <- qbeta(0.05,i,N-i+1)
}

## plot the two confidence lines
plot(null, -log(c95,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
axes=FALSE, xlab="", ylab="",lty=3,col="red",lwd=2)
par(new=T)
plot(null, -log(c05,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
axes=FALSE, xlab="", ylab="",lty=3,col="red",lwd=2)
## add the diagonal
abline(0,1,lwd=2)
par(new=T)

## add the qqplot
qqplot(null,obs, xlab="Expected -log10(p)",ylab="Observed -log10(p)",ylim=c(0,MAX),xlim=c(0,MAX), main="",pch=18)
dev.off()


