library(BRugs)

#model
modelString=" 
model
{ # i counts the two arms of all 19 studies
  for (i in 1:38)
  {
    r[i] ~ dbin(p[i], n[i]);
    logit(p[i]) <- mu[s[i]]+delta[i]*(1-equals(t[i],b[i]));
    delta[i] ~ dnorm(md[i], prec);
    md[i] <- d[t[i]]-d[b[i]];
  }
  # j represents the CABG arm
  for (j in 1:19)
  {
    mu[j] ~ dnorm(0, .001);
  }
  prec ~ dgamma(0.001, 0.001);
  d[1] <- 0;
  # K represents the relative treatment comparator: k=2 is PCI, k=3 is MT
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
}"

# Write the modelString to a file 
writeLines (modelString,con="model2.txt") 
# Use BRugs to check model 
modelCheck ("model2.txt") 


dataList=list(s=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19),t=c(1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,3,1,3,1,3,1,3,1,3,1,3,1,3),r=c(15,15,4,1,5,2,20,26,7,9,12,3,18,20,25,2,19,21,20,5,8,22,12,7,3,10,59,46,16,6,5,4,16,12,2,2,61,93),n=c(348,357,53,52,101,100,300,300,67,67,142,107,542,542,238,49,154,157,245,96,135,135,190,97,48,43,1183,309,141,24,40,17,89,32,28,31,899,440),b=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))





initsList = list(d=c(NA,0,0), prec=1, mu=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))


modelData(bugsData(dataList))
#Initialize the chains
nChain=1
modelCompile(numChains = nChain) #Compile the model

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


