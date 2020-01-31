########### Multi-scale occupancy model #############
#### Based on Kéry & Royle 2016 AHM book

rm(list=ls())

## Open library
library(jagsUI)

####----------------------------------------------------------------------------------------------------####
####----------------------------                                         -------------------------------####
####---------------------------- Bayesian Multi-scale Occupancy Analysis -------------------------------####
####----------------------------                                         -------------------------------####
####----------------------------------------------------------------------------------------------------####

# Define model in BUGS language
cat(file="MSocc-HM7.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    for(t in 1:n.cts){                #Intercepts for availability probability  
    int.theta[t] <- logit(theta0[t])
    theta0[t] ~ dunif(0,1) 
    }
    
    for(t in 1:n.days){               #Intercepts for detection probability
    int.p[t] <- logit(p0[t])
    p0[t] ~ dunif(0,1)     
    }
    
    #Slope of site covariate
    for(t in 1:4){
    beta.lpsi[t] ~ dnorm(0,0.01)
    }        
    
    
    #Slopes of ct covariates
    for(t in 1:5){        
    beta.ltheta[t] ~ dnorm(0,0.01)
    }
    
    #Slope of p covariate
    beta.lp ~ dnorm(0,0.01)
    
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi + beta.lpsi[1] * elev[i] + beta.lpsi[2] * forest[i] + 
                                      beta.lpsi[3] * field[i] + beta.lpsi[4] * human[i]
    
    for(j in 1:n.cts){
    
    #Occurrence at camera location j
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z[i] * theta[i,j]
    logit(theta[i,j]) <- int.theta[j] + 
    beta.ltheta[1] * equals(micro[i,j],2) +   # forest road compared to no feature (intercept)
    beta.ltheta[2] * equals(micro[i,j],3) +   # wildlife trail compared to no feature (intercept)
    beta.ltheta[3] * equals(micro[i,j],4) +   # cliff compared to no feature (intercept)
    beta.ltheta[4] * rugged[i,j] +            # slope for ruggedness
    beta.ltheta[5] * human.passage[i,j]       # slope for human passages
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p[k] + 
    beta.lp * walktest[i,j] 
    }
    }
    }

    # Derived parameters
    logit(theta.act) <- int.theta[1] #Site-use probability active camera on prob. scale
    logit(theta.hab) <- int.theta[2] #Site-use probability random habitat patch camera on prob. scale
    logit(theta.lnd) <- int.theta[3] #Site-use probability random landscape camera on prob. scale

    logit(micro.nof) <- int.theta[1]                  #Site-use probability active camera on no feature prob. scale
    logit(micro.for) <- int.theta[1] + beta.ltheta[1] #Site-use probability active camera on forest road prob. scale
    logit(micro.wlt) <- int.theta[1] + beta.ltheta[2] #Site-use probability active camera on wildlife trail prob. scale
    logit(micro.clf) <- int.theta[1] + beta.ltheta[3] #Site-use probability active camera on cliff prob. scale
    
    } #End Model
    ")

# Parameters monitored
params <- c("int.psi","int.theta","int.p","beta.lpsi","beta.ltheta","beta.lp","theta.act","theta.hab","theta.lnd",
            "micro.nof","micro.for","micro.wlt","micro.clf")

# MCMC settings
ni <- 50000; nt <- 10; nb <- 20000; nc <- 3
#ni <- 100; nt <- 1; nb <- 0; nc <- 3 #not run

####---------- Moose ----------####
# Load data
load("MSocc-MooseHMv20200130.RData")
attach(data) 
win.data <- list(y=y, n.site=dim(y)[1],n.cts=dim(y)[2],n.days=dim(y)[3],elev=elev,forest=forest,
                 field=field,human=human,micro=micro,rugged=rugged,walktest=walktest,human.passage=human.passage)
# Initial values
zst <- apply(y,1,max,na.rm=T)      #inits for presence (z)
ast <- apply(y,c(1,2),max,na.rm=T) #inits for availability (a)
inits <- function()list(z=zst, a=ast, psi0=0.5, beta.lpsi=c(0,0,0,0))
outJ.mooseHM2 <- jags(win.data, inits, params, "MSocc-HM7.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.mooseHM2,file="outJ-mooseHMv20200130.RData")
detach(data); rm("data")

####---------- Roe deer ----------####
# Load data
load("MSocc-Roe.deerHMv20200130.RData")
attach(data) 
win.data <- list(y=y, n.site=dim(y)[1],n.cts=dim(y)[2],n.days=dim(y)[3],elev=elev,forest=forest,
                 field=field,human=human,micro=micro,rugged=rugged,walktest=walktest,human.passage=human.passage)
# Initial values
zst <- apply(y,1,max,na.rm=T)      #inits for presence (z)
ast <- apply(y,c(1,2),max,na.rm=T) #inits for availability (a)
inits <- function()list(z=zst, a=ast, psi0=0.5, beta.lpsi=c(0,0,0,0))
outJ.roe.deerHM2 <- jags(win.data, inits, params, "MSocc-HM7.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.roe.deerHM2,file="outJ-roe-deerHMv20200130.RData")
detach(data); rm("data")

####---------- Red fox ----------####
# Load data
load("MSocc-Red.foxHMv20200130.RData")
attach(data) 
win.data <- list(y=y, n.site=dim(y)[1],n.cts=dim(y)[2],n.days=dim(y)[3],elev=elev,forest=forest,
                 field=field,human=human,micro=micro,rugged=rugged,walktest=walktest,human.passage=human.passage)
# Initial values
zst <- apply(y,1,max,na.rm=T)      #inits for presence (z)
ast <- apply(y,c(1,2),max,na.rm=T) #inits for availability (a)
inits <- function()list(z=zst, a=ast, psi0=0.5, beta.lpsi=c(0,0,0,0))
outJ.red.foxHM2 <- jags(win.data, inits, params, "MSocc-HM7.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.red.foxHM2,file="outJ-red-foxHMv20200130.RData")
detach(data); rm("data")

####---------- Mountain hare ----------####
# Load data
load("MSocc-HareHMv20200130.RData")
attach(data) 

win.data <- list(y=y, n.site=dim(y)[1],n.cts=dim(y)[2],n.days=dim(y)[3],elev=elev,forest=forest,
                 field=field,human=human,micro=micro,rugged=rugged,walktest=walktest,human.passage=human.passage)
# Initial values
zst <- apply(y,1,max,na.rm=T)      #inits for presence (z)
ast <- apply(y,c(1,2),max,na.rm=T) #inits for availability (a)
inits <- function()list(z=zst, a=ast, psi0=0.5, beta.lpsi=c(0,0,0,0))
outJ.hareHM2 <- jags(win.data, inits, params, "MSocc-HM7.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.hareHM2,file="outJ-hareHMv20200130.RData")
detach(data); rm("data")

####---------- Wolf ----------####
# Load data
load("MSocc-WolfHMv20200130.RData")
attach(data)
win.data <- list(y=y, n.site=dim(y)[1],n.cts=dim(y)[2],n.days=dim(y)[3],elev=elev,forest=forest,
                 field=field,human=human,micro=micro,rugged=rugged,walktest=walktest,human.passage=human.passage)
# Initial values
zst <- apply(y,1,max,na.rm=T)      #inits for presence (z)
ast <- apply(y,c(1,2),max,na.rm=T) #inits for availability (a)
inits <- function()list(z=zst, a=ast, psi0=0.5, beta.lpsi=c(0,0,0,0))
outJ.wolfHM2 <- jags(win.data, inits, params, "MSocc-HM7.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.wolfHM2,file="outJ-wolfHMv20200130.RData")
detach(data); rm("data")

####---------- Lynx ----------####
# Load data
load("MSocc-LynxHMv20200130.RData")
attach(data) 
win.data <- list(y=y, n.site=dim(y)[1],n.cts=dim(y)[2],n.days=dim(y)[3],elev=elev,forest=forest,
                 field=field,human=human,micro=micro,rugged=rugged,walktest=walktest,human.passage=human.passage)
# Initial values
zst <- apply(y,1,max,na.rm=T)      #inits for presence (z)
ast <- apply(y,c(1,2),max,na.rm=T) #inits for availability (a)
inits <- function()list(z=zst, a=ast, psi0=0.5, beta.lpsi=c(0,0,0,0))
outJ.lynxHM2 <- jags(win.data, inits, params, "MSocc-HM7.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.lynxHM2,file="outJ-lynxHMv20200130.RData")
detach(data); rm("data")

####---------- Wolverine ----------####
# Load data
load("MSocc-WolverineHMv20200130.RData")
attach(data) 
win.data <- list(y=y, n.site=dim(y)[1],n.cts=dim(y)[2],n.days=dim(y)[3],elev=elev,forest=forest,
                 field=field,human=human, micro=micro,rugged=rugged,walktest=walktest,human.passage=human.passage)
# Initial values
zst <- apply(y,1,max,na.rm=T)      #inits for presence (z)
ast <- apply(y,c(1,2),max,na.rm=T) #inits for availability (a)
inits <- function()list(z=zst, a=ast, psi0=0.5, beta.lpsi=c(0,0,0,0))
outJ.wolverineHM2 <- jags(win.data, inits, params, "MSocc-HM7.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.wolverineHM2,file="outJ-wolverineHMv20200130.RData")
detach(data); rm("data")

####---------- Red deer ----------####
# Load data
load("MSocc-Red.deerHMv20200130.RData")
attach(data) 
win.data <- list(y=y, n.site=dim(y)[1],n.cts=dim(y)[2],n.days=dim(y)[3],elev=elev,forest=forest,
                 field=field,human=human,micro=micro,rugged=rugged,walktest=walktest,human.passage=human.passage)
# Initial values
zst <- apply(y,1,max,na.rm=T)      #inits for presence (z)
ast <- apply(y,c(1,2),max,na.rm=T) #inits for availability (a)
inits <- function()list(z=zst, a=ast, psi0=0.5, beta.lpsi=c(0,0,0,0))
outJ.red.deerHM2 <- jags(win.data, inits, params, "MSocc-HM7.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.red.deerHM2,file="outJ-red-deerHMv20200130.RData")
detach(data); rm("data")

####---------- Badger ----------####
# Load data
load("MSocc-BadgerHMv20200130.RData")
attach(data) 
win.data <- list(y=y, n.site=dim(y)[1],n.cts=dim(y)[2],n.days=dim(y)[3],elev=elev,forest=forest,
                 field=field,human=human,micro=micro,rugged=rugged,walktest=walktest,human.passage=human.passage)
# Initial values
zst <- apply(y,1,max,na.rm=T)      #inits for presence (z)
ast <- apply(y,c(1,2),max,na.rm=T) #inits for availability (a)
inits <- function()list(z=zst, a=ast, psi0=0.5, beta.lpsi=c(0,0,0,0))
outJ.badgerHM2 <- jags(win.data, inits, params, "MSocc-HM7.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.badgerHM2,file="outJ-badgerHMv20200130.RData")
detach(data); rm("data")

####---------- Pine marten ----------####
# Load data
load("MSocc-Pine.martenHMv20200130.RData")
attach(data) 
win.data <- list(y=y, n.site=dim(y)[1],n.cts=dim(y)[2],n.days=dim(y)[3],elev=elev,forest=forest,
                 field=field,human=human,micro=micro,rugged=rugged,walktest=walktest,human.passage=human.passage)
# Initial values
zst <- apply(y,1,max,na.rm=T)      #inits for presence (z)
ast <- apply(y,c(1,2),max,na.rm=T) #inits for availability (a)
inits <- function()list(z=zst, a=ast, psi0=0.5, beta.lpsi=c(0,0,0,0))
outJ.pine.martenHM2 <- jags(win.data, inits, params, "MSocc-HM7.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.pine.martenHM2,file="outJ-pine-martenHMv20200130.RData")
detach(data); rm("data")