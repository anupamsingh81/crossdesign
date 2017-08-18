DDdat<-read.csv("DDdat.csv",as.is=TRUE, header=T)
str(DDdat)
s<-c(DDdat$s)
t<-c(DDdat$t)
r<-c(DDdat$r)
nn<-c(DDdat$nn)
b<-c(DDdat$b)
#Specify the model in BUGS language, but save it as a string in [R]
modelString="
model
{ # i counts the two arms of all 14 studies
for (i in 1:28)
{
r[i] ~ dbin(p[i], nn[i]);
logit(p[i]) <- mu[s[i]]+delta[i]*(1-equals(t[i],b[i]));
delta[i] ~ dnorm(md[i], prec);
md[i] <- d[t[i]]-d[b[i]];
}
# j represents the CABG arm
for (j in 1:14)
{
mu[j] ~ dnorm(0, .001);
}
prec ~ dgamma(0.001, 0.001);
d[1] <- 0;
# K represents the relative treatment comparator: k1 = Short, k=2 is 12 mo, k=3 is Long
for (k in 2:3)
{
  d[k] ~ dnorm(0, .001)
}
for (c in 1:2)
{
  for (k in (c+1):3)
  {
  lor[c,k] <- d[k]-d[c];
  log(or[c,k]) <- lor[c,k];
  }
}
}
"
# Write the modelString to a file
writeLines (modelString,con="model4.txt")
# Use BRugs to check model
modelCheck ("model4.txt")
#load data
dataList = list(s=c(s),
                t=c(t),
                r=c(r),
                nn=c(nn),
                b=c(b)
)
#Use BRugs commands to put the data into a file and ship the file to BUGS
modelData(bugsData(dataList))
#Initialize the chains
nChain=1
modelCompile(numChains = nChain) #Compile the model
initsList = list(d=c(NA,0,0), prec=1, mu=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
modelInits(bugsData(initsList))
modelGenInits()
#R defines a new variable to specify an arbitrary chain length
chainLength1 = 5000
#BRugs tells BUGS to generate a MCMC chain
modelUpdate (chainLength1)
#BRugs keeps a record of parameters
samplesSet(c("lor"))
#BRugs asks BUGS for summary statistics
chainLength2 = 10000
thinStep = 2
modelUpdate (chainLength2)
thetaSummaryObs = samplesStats (c("lor")); thetaSummaryObs
thetaSummaryObs<-thetaSummaryObs[order(thetaSummaryObs$mean),]
expTheta<-exp(thetaSummaryObs)
print(thetaSummaryObs)
print(expTheta)

#forest plot
x<-seq(from=-0.8,to=0.6,by=0.01)
#Short vs. 12 mo
x<-thetaSummaryObs$mean
y<-c(1,2,3)
plot(x,y,xlim=c(-0.7,0.6),ylim=c(3.5,0),pch=23,cex=4,ylab="",yaxt="n",col="black",bg="lightblue",
     cex.axis=1.0, xlab="log(e)OR", cex.lab=1.6)
axis (4, pos=0.0, tck = 0, labels=FALSE, col="black")
text (-0.555,1,"3-6 mos vs. 12 mos", cex= 1.4)
text (-0.53,3,"3-6 mos vs. 18-48 mos",cex = 1.4)
text (-0.54,2, "12 mos vs. 18-48 mos", cex = 1.4)
text (0, 0,"All-Cause Mortality",cex = 1.6,font =2)
text (0.4, 0.6, "Short DAPT Better",cex=1.6,font=3)
text (-0.4, 0.6, "Long DAPT Better",cex=1.6,font=3)
text (thetaSummaryObs$mean[3], 3.2, font=2, round(expTheta$mean[3],2))
text (thetaSummaryObs$val2.5pc[3], 3.2, font=2,round(expTheta$val2.5pc[3],2))

#text (thetaSummaryObs$val97.5pc[3], 3.2, font=2,round(expTheta$val97.5pc[3],2))
text (thetaSummaryObs$val97.5pc[3], 3.2, font=2,"1.60")
text (thetaSummaryObs$mean[1],1.2,font=2,round(expTheta$mean[1],2))
#text (thetaSummaryObs$val2.5pc[1], 1.2, font=2,round(expTheta$val2.5pc[1],2))
text (thetaSummaryObs$val2.5pc[1], 1.2, font=2,"0.70")
text (thetaSummaryObs$val97.5pc[1], 1.2, font=2,round(expTheta$val97.5[1],2))
text (thetaSummaryObs$mean[2], 2.2, font=2,round(expTheta$mean[2],2))
text (thetaSummaryObs$val2.5pc[2], 2.2, font=2,round(expTheta$val2.5pc[2],2))
text (thetaSummaryObs$val97.5pc[2], 2.2, font=2,round(expTheta$val97.5pc[2],2))
segments(thetaSummaryObs$val2.5pc[3], 3, thetaSummaryObs$mean[3]-0.025, 3, lty=1, col="black", lwd=3)
segments(thetaSummaryObs$val97.5pc[3], 3, thetaSummaryObs$mean[3]+0.025, 3, lty=1, col="black", lwd=3)
segments(thetaSummaryObs$val2.5pc[1], 1, thetaSummaryObs$mean[1]-0.025, 1, lty=1, lwd=3)
segments(thetaSummaryObs$val97.5pc[1], 1, thetaSummaryObs$mean[1]+0.025, 1, lty=1, lwd=3)
segments(thetaSummaryObs$val2.5pc[2], 2, thetaSummaryObs$mean[2]-0.025, 2, lty=1, lwd=3)
segments(thetaSummaryObs$val97.5pc[2], 2, thetaSummaryObs$mean[2]+0.025, 2, lty=1, lwd=3)
mtext ("Posterior Odds Ratio (OR)",3, line =2, cex = 1.6)
axis (3, at=c(-0.91,-0.69, -0.51,-0.35, -0.22, -0.105,0.0, 0.095,0.182, 0.262,0.336,
              0.405,0.47,0.531,0.588,0.693, 0.833, 0.956, 1.10,1.19, 1.281,1.386,1.46,1.53,1.61,1.67,1.72,1.79),
      labels=c(0.4,0.5,0.6, 0.7, 0.8,0.9, "1.0", 1.1,1.2, 1.3,1.4,1.5, 1.6,1.7, 1.8, "2.0", 2.3, 2.6, "3.0",
               3.3,3.6,"4.0",4.3,4.6,"5.0",5.3,5.6,"6.0"))
#To create good margins
mar.default <- c(5,4,4,2) + 0.0
par(mar = mar.default + c(0, 2, 0, 0))
#To copy in eps and pdf formats to your original folder. (Change the date each time or you will
overwrite.)
dev.copy2eps(file="NetworkDAPTDeathJun10Caterpillar.eps")
dev.copy2pdf(file="NetworkDAPTDeathJun10Caterpillar.pdf")
