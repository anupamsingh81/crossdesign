
library(BRugs)
#Specify the model in BUGS language, but save it as a string in [R] 
modelString=" 
model 
{ 
# K1 is the number of trials; 
for (k in 1:18) 
{ 
# calculate odds ratios; 
or[k] <- ((r.multi[k]+0.5)/(n.multi[k]-r.multi[k]+0.5))/((r.culprit[k]+0.5)/(n.culprit[k]-r.culprit[k]+0.5)) 
logor[k] <- log(or[k]); 
varlogor[k] <- (1/(r.multi[k]+0.5))+(1/(n.multi[k]-r.multi[k]+0.5))+(1/(r.culprit[k]+0.5))+(1/(n.culprit[k]-r.culprit[k]+0.5)); 
invlogor[k] <- 1/varlogor[k]; 
logor[k] ~ dnorm(theta[k], invlogor[k]); 
or.est[k] <- exp(theta[k]); 
# study-type level random-effects distributions 
theta[k] ~ dnorm(mu.theta.study[study[k]], prec.theta.study[study[k]]); 
} 
# K2 is the number of study types 
for (l in 1:3) 
{ 
mu.theta.study[l] ~ dnorm(mu.theta, prec.theta); 
or.theta.study[l] <- exp(mu.theta.study[l]); 
prec.theta.study[l] <- 1/(tau.theta.study[l]*tau.theta.study[l]); 
# prior distribution for tau.theta.study based on HN[0.36^2], giving precision 7.72 ,based on 95% belief that the
#true underlying OR for a study of a particular type will be <4x or >1/4 the overall OR of that type
tau.theta.study[l] ~ dnorm(0, 7.72)I(0,); 
}  
# prior distribution for mu.theta based on log(500)/1.96 = 3.17 for N[0,10], giving precision 0.1 
mu.theta ~ dnorm(0, 0.1); 
# prior distribution for tau.theta based on HN[0.18^2], giving precision 30.86 , ,based on 95% belief that the
# the underlying risk ratio for a particular study type will <2× or >1⁄2 the overall population effect
tau.theta ~ dnorm(0, 30.86)I(0,); 
prec.theta <- 1/(tau.theta*tau.theta); 
# global summary odds ratio; 
or.theta <- exp(mu.theta); 
# K1 is the number of trials; 
# DATA list(K1=21, K2=3); 
# INITIAL VALUES list(mu.theta=0, tau.theta = 1); 
# BUGS model specification ends 
} 
" 
# Write the modelString to a file 
writeLines (modelString,con="model.txt") 
# Use BRugs to check model 
modelCheck ("model.txt") 
#load data 
dataList = list(n.multi=c(52, 65, 234, 150, 79, 503, 403, 26, 95, 147, 3134, 70, 217, 
                          419, 1108, 442, 77, 367), 
                n.culprit=c(17, 84, 231, 146, 79, 503, 2418, 354, 25, 156, 25802, 707, 1984, 
                            2118, 3833, 1467, 180, 706), 
                r.multi=c(1, 6, 12, 4, 19, 59, 41, 5, 9, 12, 246, 11, 27, 6, 81, 12, 2, 26), 
                r.culprit=c(0, 13, 16, 10, 13, 54, 164, 42, 2, 8, 1321, 57, 111, 72, 168, 40, 14, 
                            127), 
                study=c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3) 
) 
#Use BRugs commands to put the data into a file and ship the file to BUGS 
modelData(bugsData(dataList)) 
#Initialize the chains 
nChain=1 
modelCompile(numChains = nChain) 
#Compile the model 
initsList = list(mu.theta=0, tau.theta=1) 
modelInits(bugsData(initsList)) 
modelGenInits() 
#R defines a new variable to specify an arbitrary chain length 
chainLength1 = 5000 
#BRugs tells BUGS to generate a MCMC chain 
modelUpdate (chainLength1) 
#BRugs keeps a record of parameters 
samplesSet(c("mu.theta","prec.theta","or.theta","tau.theta")) 
#BRugs asks BUGS for summary statistics 
chainLength2 = 10000 
thinStep = 2 
modelUpdate (chainLength2) 
thetaSummary = samplesStats (c("mu.theta","prec.theta","or.theta","tau.theta"));  
print(thetaSummary) 


#
# http://www.di.fc.ul.pt/~jpn/r/bugs/part2.html (refernce)

# Extract samples
OR  <- samplesSample( "or.theta" )
summary(OR)