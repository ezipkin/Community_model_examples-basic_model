#The model code below is written for program R and uses the R2WinBUGS package 
#to run WinBUGS as well as the reshape package to format the occurrence data.

#It is designed to estimate static species-specific occupancy and detection, 
#constant across sampling locations using the community model. The model also
#estimates the total species richness N, using data augmentation.
#The data are found in the file "occ data.csv". 

#Read in the occurence data
data1 <- read.table("occ data.csv", header=TRUE,sep=",",na.strings=TRUE)
data1$Occ <- rep(1, dim(data1)[1])
#See the first ten lines of data
data1[1:10,]
#How many citings for each species
total.count = tapply(data1$Occ, data1$Species, sum)

#Find the number of unique species
uspecies = as.character(unique(data1$Species))
#n is the number of observed species
n=length(uspecies)

#Find the number of unique sampling locations
upoints = as.character(unique(data1$Point))
#J is the number of sampled points
J=length(upoints)

#Reshape the data using the R package "reshape"
library(reshape)

#The detection/non-detection data is reshaped into a three dimensional 
#array X where the first dimension, j, is the point; the second 
#dimension, k, is the rep; and the last dimension, i, is the species. 
junk.melt=melt(data1,id.var=c("Species", "Point", "Rep"), measure.var="Occ")
X=cast(junk.melt, Point ~ Rep ~ Species)

#Add in the missing lines with NAs
for (i in 1: dim(X)[3]) {
   b = which(X[,,i] > 0) 
   X[,,i][b] = 1  
   X[,,i][-b] = 0  
   X[,,i][1:36,4] = NA;  X[,,i][38:56,4] = NA;  
   X[,,i][59:61,4] = NA;  X[,,i][66:70,4] = NA;        
}

#Create all zero encounter histories to add to the detection array X 
#as part of the data augmentation to account for additional 
#species (beyond the n observed species). 

#nzeroes is the number of all zero encounter histories to be added
  nzeroes = 50
#X.zero is a matrix of zeroes, including the NAs for when a point has not been sampled  
  X.zero = matrix(0, nrow=70, ncol=4)
  X.zero[1:36,4] = NA;  X.zero[38:56,4] = NA;  
  X.zero[59:61,4] = NA;  X.zero[66:70,4] = NA;   
#Xaug is the augmented version of X.  The first n species were actually observed
#and the n+1 through nzeroes species are all zero encounter histories  
  Xaug <- array(0, dim=c(dim(X)[1],dim(X)[2],dim(X)[3]+nzeroes))
  Xaug[,,(dim(X)[3]+1):dim(Xaug)[3]] = rep(X.zero, nzeroes)
  dimnames(X)=NULL
  Xaug[,,1:dim(X)[3]] <-  X

#K is a vector of length J indicating the number of reps at each point j  
KK <- X.zero
a=which(KK==0); KK[a] <- 1
K=apply(KK,1,sum, na.rm=TRUE)
K=as.vector(K)

################

#Write the model code to a text file (used to run WinBUGS)
cat("
	model{

#Define prior distributions for community-level model parameters

#The parameter estimating the prob that a species is included in the "hyper community"
omega ~ dunif(0,1)

#Prior distribution for the community level hyper parameter on mean occupancy
u.mean ~ dunif(0,1)	
mu.u <- log(u.mean) - log(1-u.mean)

#Prior distribution for the community level hyper parameter on mean detection
v.mean ~ dunif(0,1)
mu.v <- log(v.mean) - log(1-v.mean)

#Prior distribution for the community level hyper parameter for the precision (1/var) 
#on occupancy and detection
tau.u ~ dgamma(0.1,0.1)  
tau.v ~ dgamma(0.1,0.1)

#Loop over all species i (including the n observed species and the nzeroes nonobserved species
# that may or may not be in the hyper community
for (i in 1:(n+nzeroes)) {

#Create priors for species i from the community level prior distributions
    
	#Binary variable indicating whether a species is in fact in the hyper community
    # will always be a 1 if at least one individual was detected
	w[i] ~ dbern(omega)
	
	#Prior distribution for intercept term of species occupancy 
    u[i] ~ dnorm(mu.u, tau.u)
	#Prior distribution for intercept term of species detection
    v[i] ~ dnorm(mu.v, tau.v)    

#Create a loop to estimate the Z matrix (true occurrence for species i 
#at point j) for all J sites.      
   for (j in 1:J) {
       logit(psi[j,i]) <- u[i] 
       
  mu.psi[j,i] <- psi[j,i]*w[i]
  Z[j,i] ~ dbern(mu.psi[j,i])

#Create a loop to estimate detection for species i at point k during 
#sampling period k at site j.      
     for (k in 1:K[j]) {  	
    	logit(p[j,k,i]) <-  v[i]	
       mu.p[j,k,i] <- p[j,k,i]*Z[j,i]
       X[j,k,i] ~ dbern(mu.p[j,k,i])
}   	}		}

#Sum all species observed (n) and unobserved species (n0) to find the 
#total estimated richness
n0 <- sum(w[(n+1):(n+nzeroes)])	
N <- n + n0

#Finish writing the text file into a document called basicmodel.txt
}
",file="basicmodel.txt")

#Load the R2Winbugs library
library(R2WinBUGS)

#Create the necessary arguments to run the bugs() command 
#Load all the data
sp.data = list(n=n, nzeroes=nzeroes, J=J, K=K, X=Xaug)

#Specify the parameters to be monitored
sp.params = list('u', 'v', 'mu.u', 'mu.v', 'tau.u', 'tau.v', 'omega', 'N')

#Specify the initial values
    sp.inits = function() {
    omegaGuess = runif(1, n/(n+nzeroes), 1)
    psi.meanGuess = runif(1, .25,1)
    list(omega=omegaGuess,w=c(rep(1, n), rbinom(nzeroes, size=1, prob=omegaGuess)),
               u=rnorm(n+nzeroes), v=rnorm(n+nzeroes),
               Z = matrix(rbinom((n+nzeroes)*J, size=1, prob=psi.meanGuess), 
		   nrow=J, ncol=(n+nzeroes))
               )
           }

#Run the model and call the results “fit”
fit = bugs(sp.data, sp.inits, sp.params, "basicmodel.txt", debug=TRUE, 
       #I usually run for a short amount to check if it's working        
		debug=TRUE, n.chains=2, n.iter=100, n.burnin=10, n.thin=1)
		#If it works, pound out the statement above and use the line below
		#If the model doesn't converge, increase n.iter and n.burnin
         #debug=TRUE, n.chains=3, n.iter=6000, n.burnin=4000, n.thin=3)


#The model code below is written in the R language and designed to be run 
#after the “covariate model code” for summary of the model results

#See a summary of the parameter estimates
fit$summary

#See baseline estimates of species-specific occupancy and detection in one of 
#the habitat types (CATO)
species.occ = fit$sims.list$u
species.det = fit$sims.list$v

#Show occupancy and detection estimates for only the observed species (1:n)
psi = plogis(species.occ[,1:n]) 
p   = plogis(species.det[,1:n]) 

occ.matrix <- cbind(apply(psi,2,mean),apply(psi,2,sd))
colnames(occ.matrix) = c("mean occupancy", "sd occupancy")
rownames(occ.matrix) = uspecies
det.matrix <- cbind(apply(p,2,mean),apply(p,2,sd))
colnames(det.matrix) = c("mean detection", "sd detection")

round(occ.matrix, digits=2)

#See estimates of total richness (N)
N = fit$sims.list$N
mean(N) 
summary(N) 
table(N)
plot(table(N))