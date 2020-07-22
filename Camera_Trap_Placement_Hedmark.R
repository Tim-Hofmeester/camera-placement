####----------------------------------------------------------------------------------------------------####
####----------------------------                                         -------------------------------####
####----------------------------      Camera trap placement analysis     -------------------------------####
####----------------------------                                         -------------------------------####
####----------------------------------------------------------------------------------------------------####

rm(list=ls())

####----------------------------  Defining all models in BUGS language  --------------------------------####

#### The five full models including Bayesian inclusion parameters ####

## Multi-scale model for all three cameras ##

# Define model in BUGS language
cat(file="MSocc-HM-all.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    
    #Slope of site covariate
    for(t in 1:3){
    beta.lpsi[t] ~ dnorm(0,0.01)
    }        
    
    #Slopes of ct covariates
    for(t in 1:4){        
    beta.ltheta[t] ~ dnorm(0,0.01)
    }
    
    #Slope of p covariates
    for(t in 1:6){
    beta.lp[t] ~ dnorm(0,0.01)
    }

    #Inclusion parameter priors
    for(c in 1:11){
    w[c] ~ dbern(0.5)
    }
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi + 
    w[1] * beta.lpsi[1] * elev[i] + 
    w[2] * beta.lpsi[2] * forest.grid[i] + 
    w[3] * beta.lpsi[3] * field.grid[i] 
    
    for(j in 1:n.cts){
    
    #Occurrence at camera location j
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z[i] * theta[i,j]
    logit(theta[i,j]) <- int.theta + 
    w[4] * beta.ltheta[1] * target.site[i,j] +          # if habitat patch is targeted (1) or not (0)
    w[5] * beta.ltheta[2] * ruggedness[i,j] +           # slope for ruggedness
    w[6] * beta.ltheta[3] * forest.site[i,j] +          # forest cover 1km2 around camera
    w[7] * beta.ltheta[4] * field.site[i,j]             # field cover 1km2 around camera
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p + 
    w[8] * beta.lp[1] * target.cam[i,j] +          # if camera is targeted (1) or not (0)
    w[9] * beta.lp[2] * equals(micro[i,j],2) +     # forest road compared to no feature (intercept)
    w[9] * beta.lp[3] * equals(micro[i,j],3) +     # wildlife trail compared to no feature (intercept)
    w[9] * beta.lp[4] * equals(micro[i,j],4) +     # cliff compared to no feature (intercept)
    w[10]* beta.lp[5] * cam.height[i,j] +          # height of the camera from the ground
    w[11]* beta.lp[6] * walktest[i,j]
    }
    }
    }
    } #End Model
    ")

## Multi-scale model for all targeted and landscape cameras ##

# Define model in BUGS language
cat(file="MSocc-HM-TL.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    #Slope of site covariate
    for(t in 1:3){
    beta.lpsi[t] ~ dnorm(0,0.01)
    }        
    
    #Slopes of ct covariates
    for(t in 1:4){        
    beta.ltheta[t] ~ dnorm(0,0.01)
    }
    
    #Slope of p covariates
    for(t in 1:6){
    beta.lp[t] ~ dnorm(0,0.01)
    }

    #Inclusion parameter priors
    for(c in 1:11){
    w[c] ~ dbern(0.5)
    }
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi + 
    w[1] * beta.lpsi[1] * elev[i] + 
    w[2] * beta.lpsi[2] * forest.grid[i] + 
    w[3] * beta.lpsi[3] * field.grid[i] 
    
    for(j in 1:n.cts){
    
    #Occurrence at camera location j
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z[i] * theta[i,j]
    logit(theta[i,j]) <- int.theta + 
    w[4] * beta.ltheta[1] * target.site[i,j] +            # if camera is targeted (1) or not (0)
    w[5] * beta.ltheta[2] * ruggedness[i,j] +            # slope for ruggedness
    w[6] * beta.ltheta[3] * forest.site[i,j] +          # forest cover 1km2 around camera
    w[7] * beta.ltheta[4] * field.site[i,j]             # field cover 1km2 around camera

    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p + 
    w[8] * beta.lp[1] * target.cam[i,j] +
    w[9] * beta.lp[2] * equals(micro[i,j],2) +   # forest road compared to no feature (intercept)
    w[9] * beta.lp[3] * equals(micro[i,j],3) +   # wildlife trail compared to no feature (intercept)
    w[9] * beta.lp[4] * equals(micro[i,j],4) +   # cliff compared to no feature (intercept)
    w[10]* beta.lp[5] * cam.height[i,j] +
    w[11]* beta.lp[6] * walktest[i,j]
    }
    }
    }
    } #End Model
    ")

## Multi-scale model for all targeted and habitat patch cameras ##

# Define model in BUGS language
cat(file="MSocc-HM-TH.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    #Slope of site covariate
    for(t in 1:3){
    beta.lpsi[t] ~ dnorm(0,0.01)
    }        
    
    
    #Slopes of ct covariates
    for(t in 1:3){        
    beta.ltheta[t] ~ dnorm(0,0.01)
    }
    
    #Slope of p covariates
    for(t in 1:6){
    beta.lp[t] ~ dnorm(0,0.01)
    }

    #Inclusion parameter priors
    for(c in 1:10){
    w[c] ~ dbern(0.5)
    }
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi + 
    w[1] * beta.lpsi[1] * elev[i] + 
    w[2] * beta.lpsi[2] * forest.grid[i] + 
    w[3] * beta.lpsi[3] * field.grid[i] 
    
    for(j in 1:n.cts){
    
    #Occurrence at camera location j
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z[i] * theta[i,j]
    logit(theta[i,j]) <- int.theta + 
    w[4] * beta.ltheta[1] * ruggedness[i,j] +            # slope for ruggedness
    w[5] * beta.ltheta[2] * forest.site[i,j] +          # forest cover 1km2 around camera
    w[6] * beta.ltheta[3] * field.site[i,j]             # field cover 1km2 around camera

    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p + 
    w[7] * beta.lp[1] * target.cam[i,j] +
    w[8] * beta.lp[2] * equals(micro[i,j],2) +   # forest road compared to no feature (intercept)
    w[8] * beta.lp[3] * equals(micro[i,j],3) +   # wildlife trail compared to no feature (intercept)
    w[8] * beta.lp[4] * equals(micro[i,j],4) +   # cliff compared to no feature (intercept)
    w[9] * beta.lp[5] * cam.height[i,j] +
    w[10]* beta.lp[6] * walktest[i,j]
    }
    }
    }
    } #End Model
    ")

## Single-scale model for targeted cameras ##

# Define model in BUGS language
cat(file="occ-HM-T.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    #Slope of site covariate
    for(t in 1:6){
    beta.lpsi[t] ~ dnorm(0,0.01)
    }        
    
    #Slope of p covariates
    for(t in 1:5){
    beta.lp[t] ~ dnorm(0,0.01)
    }

    #Inclusion parameter priors
    for(c in 1:9){
    w[c] ~ dbern(0.5)
    }
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi + 
    w[1] * beta.lpsi[1] * elev[i] + 
    w[2] * beta.lpsi[2] * forest.grid[i] + 
    w[3] * beta.lpsi[3] * field.grid[i] + 
    w[4] * beta.lpsi[4] * ruggedness[i] +
    w[5] * beta.lpsi[5] * forest.site[i] +
    w[6] * beta.lpsi[6] * field.site[i]
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,k] ~ dbern(mu.y[i,k])
    mu.y[i,k] <- z[i] * p[i,k]
    logit(p[i,k]) <- int.p + 
    w[7] * beta.lp[1] * equals(micro[i],2) +   # forest road compared to no feature (intercept)
    w[7] * beta.lp[2] * equals(micro[i],3) +   # wildlife trail compared to no feature (intercept)
    w[7] * beta.lp[3] * equals(micro[i],4) +   # cliff compared to no feature (intercept)
    w[8] * beta.lp[4] * cam.height[i] +
    w[9] * beta.lp[5] * walktest[i]
    }
    }
    } #End Model
    ")

## Single-scale model for landscape cameras ##

# Define model in BUGS language
cat(file="occ-HM-L.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    #Slope of site covariate
    for(t in 1:6){
    beta.lpsi[t] ~ dnorm(0,0.01)
    }        
    
    #Slope of p covariates
    for(t in 1:4){
    beta.lp[t] ~ dnorm(0,0.01)
    }

    #Inclusion parameter priors
    for(c in 1:9){
    w[c] ~ dbern(0.5)
    }
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi + 
    w[1] * beta.lpsi[1] * elev[i] + 
    w[2] * beta.lpsi[2] * forest.grid[i] + 
    w[3] * beta.lpsi[3] * field.grid[i] + 
    w[4] * beta.lpsi[4] * ruggedness[i] +
    w[5] * beta.lpsi[5] * forest.site[i] +
    w[6] * beta.lpsi[6] * field.site[i]
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,k] ~ dbern(mu.y[i,k])
    mu.y[i,k] <- z[i] * p[i,k]
    logit(p[i,k]) <- int.p + 
    w[7] * beta.lp[1] * equals(micro[i],3) +   # wildlife trail compared to no feature (intercept)
    w[7] * beta.lp[2] * equals(micro[i],4) +    # cliff compared to no feature (intercept)
    w[8] * beta.lp[3] * cam.height[i] +
    w[9] * beta.lp[4] * walktest[i]
    }
    }
    } #End Model
    ")

#### The null models ####

## Multi-scale model ##

# model
cat(file="MSocc-HMnull.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi 
    
    for(j in 1:n.cts){
    
    #Occurrence at camera location j
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z[i] * theta[i,j]
    logit(theta[i,j]) <- int.theta 
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p
    }
    }
    }
    } #End Model
    ")

## Single-scale model ##

# model
cat(file="occ-HMnull.txt", 
    "
    model{
    
   # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,k] ~ dbern(mu.y[i,k])
    mu.y[i,k] <- z[i] * p[i,k]
    logit(p[i,k]) <- int.p 
    }
    }
    } #End Model
    ")

#### Final models all three cameras ####

## Roe deer ##

# model
cat(file="MSocc-HM-all-roedeer.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    
    #Slope of site covariate
    for(t in 1:2){
    beta.ltheta[t] ~ dnorm(0,0.01)
    }
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi 
    
    for(j in 1:n.cts){
    
    #Occurrence at camera location j
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z[i] * theta[i,j]
    logit(theta[i,j]) <- int.theta + 
    beta.ltheta[1] * target.site[i,j] +          # if camera is targeted (1) or not (0)
    beta.ltheta[2] * field.site[i,j]             # percentage of field in 1km2 surrounding camera
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p
    }
    }
    }
    } #End Model
    ")

## Mountain hare ##

# model
cat(file="MSocc-HM-all-hare.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    

    #Slope of p covariate
    beta.lp ~ dnorm(0,0.01)
    

    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi 
    
    for(j in 1:n.cts){
    
    #Occurrence at camera location j
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z[i] * theta[i,j]
    logit(theta[i,j]) <- int.theta 
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p + 
    beta.lp * target.cam[i,j]           # if camera was targeted (1) or not (0)
    }
    }
    }
    } #End Model
    ")

## Red fox ##

# model
cat(file="MSocc-HM-all-fox.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    
    #Slope of site covariate
    beta.ltheta ~ dnorm(0,0.01)

    #Slope of p covariates
    beta.lp ~ dnorm(0,0.01)
    
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi 
    
    for(j in 1:n.cts){
    
    #Occurrence at camera location j
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z[i] * theta[i,j]
    logit(theta[i,j]) <- int.theta + 
    beta.ltheta * target.site[i,j]            # if the cameras were in a targeted patch (1) or not (0)
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p + 
    beta.lp * target.cam[i,j]                 # if the cameras were in a targeted micro site (1) or not (0)
    }
    }
    }
    } #End Model
    ")

#### Final models all targeted + landscape cameras ####

## Roe deer ##

# model
cat(file="MSocc-HM-TL-roedeer.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
 
    
    #Slope for site covariate
    beta.ltheta ~ dnorm(0,0.01)
    
    #Slope for p covariates
    for(t in 1:2){
    beta.lp[t] ~ dnorm(0,0.01)
    }


    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi 
    
    for(j in 1:n.cts){
    
    #Occurrence at camera location j
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z[i] * theta[i,j]
    logit(theta[i,j]) <- int.theta +
    beta.ltheta * field.site[i,j]
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p +
    beta.lp[1] * target.cam[i,j] +
    beta.lp[2] * walktest[i,j]
    }
    }
    }
    } #End Model
    ")

## Mountain hare ##

# model
cat(file="MSocc-HM-TL-hare.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    
    #Slope for p covariate
    beta.lp ~ dnorm(0,0.01)

    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi 
    
    for(j in 1:n.cts){
    
    #Occurrence at camera location j
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z[i] * theta[i,j]
    logit(theta[i,j]) <- int.theta 
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p +
    beta.lp * target.cam[i,j]
    }
    }
    }
    } #End Model
    ")

## Red fox ##

# model
cat(file="MSocc-HM-TL-fox.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    
    #Slope for site covariate
    beta.ltheta ~ dnorm(0,0.01)

    #Slope for p covariate
    beta.lp ~ dnorm(0,0.01)

    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi 
    
    for(j in 1:n.cts){
    
    #Occurrence at camera location j
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z[i] * theta[i,j]
    logit(theta[i,j]) <- int.theta +
    beta.ltheta * target.site[i,j]
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p +
    beta.lp * walktest[i,j]
    }
    }
    }
    } #End Model
    ")


#### Final models all targeted + habitat patch cameras ####

## Roe deer ##

# model
cat(file="MSocc-HM-TH-roedeer.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    
    #Slope for site covariate
    beta.ltheta ~ dnorm(0,0.01)

    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi 
    
    for(j in 1:n.cts){
    
    #Occurrence at camera location j
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z[i] * theta[i,j]
    logit(theta[i,j]) <- int.theta +
    beta.ltheta * field.site[i,j]
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p 
    }
    }
    }
    } #End Model
    ")

## Mountain hare ##

# model
cat(file="MSocc-HM-TH-hare.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    
    #Slope for p covariate
    beta.lp ~ dnorm(0,0.01)
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi 
    
    for(j in 1:n.cts){
    
    #Occurrence at camera location j
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z[i] * theta[i,j]
    logit(theta[i,j]) <- int.theta 
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p +
    beta.lp * target.cam[i,j]
    }
    }
    }
    } #End Model
    ")

## Red fox ##

# model
cat(file="MSocc-HM-TH-fox.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    
    #Slope for p covariate
    beta.lp ~ dnorm(0,0.01)
   

    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi 
    
    for(j in 1:n.cts){
    
    #Occurrence at camera location j
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z[i] * theta[i,j]
    logit(theta[i,j]) <- int.theta 
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p +
    beta.lp[1] * target.cam[i,j] 
    }
    }
    }
    } #End Model
    ")


#### Final models targeted cameras ####

## Moose ##

# model
cat(file="occ-HM-T-moose.txt", 
    "
    model{
    
   # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    
    #Slope for occupancy covariate
    beta.lpsi ~ dnorm(0,0.01)

    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi +
    beta.lpsi * ruggedness[i]
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,k] ~ dbern(mu.y[i,k])
    mu.y[i,k] <- z[i] * p[i,k]
    logit(p[i,k]) <- int.p 
    }
    }
    } #End Model
    ")

## Roe deer ##

# model
cat(file="occ-HM-T-roedeer.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    #Slope for p covariate
    beta.lp ~ dnorm(0,0.01)
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi 

    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,k] ~ dbern(mu.y[i,k])
    mu.y[i,k] <- z[i] * p[i,k]
    logit(p[i,k]) <- int.p +
    beta.lp * walktest[i]
    }
    }
    } #End Model
    ")

## Red squirrel ##

# model
cat(file="occ-HM-T-squirrel.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    #Slopes for psi covariates
    for(t in 1:3){
    beta.lpsi[t] ~ dnorm(0,0.01)
    }

    #Slope for p covariate
    beta.lp ~ dnorm(0,0.01)
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi +
    beta.lpsi[1] * forest.grid[i] +
    beta.lpsi[2] * field.grid[i] +
    beta.lpsi[3] * forest.site[i]
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,k] ~ dbern(mu.y[i,k])
    mu.y[i,k] <- z[i] * p[i,k]
    logit(p[i,k]) <- int.p +
    beta.lp * walktest[i]
    }
    }
    } #End Model
    ")

## Wolf ##

# model
cat(file="occ-HM-T-wolf.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    # slopes for p covariate
    for(t in 1:3){
    beta.lp[t] ~ dnorm(0,0.01)
    }

    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,k] ~ dbern(mu.y[i,k])
    mu.y[i,k] <- z[i] * p[i,k]
    logit(p[i,k]) <- int.p +
    beta.lp[1] * equals(micro[i],2) +   # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i],4)    # cliff compared to no feature (intercept)
    }
    }
    } #End Model
    ")

## Badger ##

# model
cat(file="occ-HM-T-badger.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     
    
    #Slope for psi covariate
    beta.lpsi ~ dnorm(0,0.01)

    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi +
    beta.lpsi * elev[i]
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,k] ~ dbern(mu.y[i,k])
    mu.y[i,k] <- z[i] * p[i,k]
    logit(p[i,k]) <- int.p
    }
    }
    } #End Model
    ")

## Red fox ##

# model
cat(file="occ-HM-T-fox.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)     

    #Slope for p covariate
    beta.lp ~ dnorm(0,0.01)
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,k] ~ dbern(mu.y[i,k])
    mu.y[i,k] <- z[i] * p[i,k]
    logit(p[i,k]) <- int.p +
    beta.lp * walktest[i]
    }
    }
    } #End Model
    ")

#### Final models landscape cameras ####

## Mountain hare ##

# model
cat(file="occ-HM-L-hare.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1) 
    
    #Slope for p covariates
    for(t in 1:3){
    beta.lp[t] ~ dnorm(0,0.01)
    }
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi 
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,k] ~ dbern(mu.y[i,k])
    mu.y[i,k] <- z[i] * p[i,k]
    logit(p[i,k]) <- int.p +
    beta.lp[1] * equals(micro[i],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[2] * equals(micro[i],4) +    # cliff compared to no feature (intercept)
    beta.lp[3] * walktest[i]
    }
    }
    } #End Model
    ")

## Red fox ##

# model
cat(file="occ-HM-L-fox.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1) 
    
    #Slope for occupancy covariate
    for(t in 1:2){
    beta.lpsi[t] ~ dnorm(0,0.01)
    }

    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi +
    beta.lpsi[1] * field.grid[i] + 
    beta.lpsi[2] * ruggedness[i]
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,k] ~ dbern(mu.y[i,k])
    mu.y[i,k] <- z[i] * p[i,k]
    logit(p[i,k]) <- int.p
    }
    }
    } #End Model
    ")

####----------------------------  Running the models  --------------------------------####

## Open library
library(jagsUI)

# Parameters monitored
params <- c("int.psi","int.theta","int.p","beta.lpsi","beta.ltheta","beta.lp","w")

# MCMC settings
## models with Bayesian inclusion parameters
ni <- 100000; nt <- 10; nb <- 50000; nc <- 3
## models without Bayesian inclusion parameters
ni1 <- 50000; nt1 <- 10; nb1 <- 20000; nc1 <- 3

#### Moose ####

# Load data
load("MSocc-MooseHM.RData")
attach(data) 

## All cameras ##
win.data <- list(y=y, n.site=dim(y)[1],n.cts=dim(y)[2],n.days=dim(y)[3],
                 elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                 target.site=target.site,ruggedness=ruggedness,forest.site=forest.site,field.site=field.site,
                 target.cam=target.cam,micro=micro,cam.height=cam.height,walktest=walktest)
# Initial values
zst <- apply(y,1,max,na.rm=T)      #inits for presence (z)
ast <- apply(y,c(1,2),max,na.rm=T) #inits for availability (a)
inits <- function()list(z=zst, a=ast, psi0=0.5, beta.lpsi=c(0,0,0))
inits.1 <- function()list(z=zst, a=ast, psi0=0.5)
# full model
outJ.mooseHM.all <- jags(win.data, inits, params, "MSocc-HM-all.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.mooseHM.all,file="outJ-mooseHM-all.RData")
# null model
outJ.mooseHM.all.nm <- jags(win.data, inits.1, params, "MSocc-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.mooseHM.all.nm,file="outJ-mooseHM-all-nm.RData")

## T+L cameras ##
y1 <- y[,c(1,3),]
win.data1 <- list(y=y1, n.site=dim(y1)[1],n.cts=dim(y1)[2],n.days=dim(y1)[3],
                 elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                 target.site=target.site[,c(1,3)],ruggedness=ruggedness[,c(1,3)],
                 forest.site=forest.site[,c(1,3)],field.site=field.site[,c(1,3)],
                 target.cam=target.cam[,c(1,3)],micro=micro[,c(1,3)],cam.height=cam.height[,c(1,3)],walktest=walktest[,c(1,3)])
# Initial values
zst1 <- apply(y1,1,max,na.rm=T)      #inits for presence (z)
ast1 <- apply(y1,c(1,2),max,na.rm=T) #inits for availability (a)
inits1 <- function()list(z=zst1, a=ast1, psi0=0.5, beta.lpsi=c(0,0,0))
inits1.1 <- function()list(z=zst1, a=ast1, psi0=0.5)
# full model
outJ.mooseHM.TL <- jags(win.data1, inits1, params, "MSocc-HM-TL.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.mooseHM.TL,file="outJ-mooseHM-TL.RData")
# null model
outJ.mooseHM.TL.nm <- jags(win.data1, inits1.1, params, "MSocc-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.mooseHM.TL.nm,file="outJ-mooseHM-TL-nm.RData")

## T+H cameras ##
y2 <- y[,c(1:2),]
win.data2 <- list(y=y2, n.site=dim(y2)[1],n.cts=dim(y2)[2],n.days=dim(y2)[3],
                 elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                 target.site=target.site,ruggedness=ruggedness,forest.site=forest.site,field.site=field.site,
                 target.cam=target.cam,micro=micro,cam.height=cam.height,walktest=walktest)
# Initial values
zst2 <- apply(y2,1,max,na.rm=T)      #inits for presence (z)
ast2 <- apply(y2,c(1,2),max,na.rm=T) #inits for availability (a)
inits2 <- function()list(z=zst2, a=ast2, psi0=0.5, beta.lpsi=c(0,0,0))
inits2.1 <- function()list(z=zst2, a=ast2, psi0=0.5)
# full model
outJ.mooseHM.TH <- jags(win.data2, inits2, params, "MSocc-HM-TH.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.mooseHM.TH,file="outJ-mooseHM-TH.RData")
# null model
outJ.mooseHM.TH.nm <- jags(win.data2, inits2.1, params, "MSocc-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.mooseHM.TH.nm,file="outJ-mooseHM-TH-nm.RData")

## T cameras ##
y3 <- y[,1,]
win.data3 <- list(y=y3, n.site=dim(y3)[1],n.days=dim(y3)[2],
                 elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                 ruggedness=ruggedness[,1],forest.site=forest.site[,1],field.site=field.site[,1],
                 micro=micro[,1],cam.height=cam.height[,1],walktest=walktest[,1])
# Initial values
zst3 <- apply(y3,1,max,na.rm=T)      #inits for presence (z)
inits3 <- function()list(z=zst3, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits3.1 <- function()list(z=zst3, psi0=0.5)
inits3.2 <- function()list(z=zst3, psi0=0.5, beta.lpsi=0)
# full model
outJ.mooseHM.T <- jags(win.data3, inits3, params, "occ-HM-T.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.mooseHM.T,file="outJ-mooseHM-T.RData")
# null model
outJ.mooseHM.T.nm <- jags(win.data3, inits3.1, params, "occ-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.mooseHM.T.nm,file="outJ-mooseHM-T-nm.RData")
# final model
outJ.mooseHM.T.fm <- jags(win.data3, inits3.2, params, "occ-HM-T-moose.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.mooseHM.T.fm,file="outJ-mooseHM-T-fm.RData")

## L cameras ##
y4 <- y[,3,]
win.data4 <- list(y=y4, n.site=dim(y4)[1],n.days=dim(y4)[2],
                 elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                 ruggedness=ruggedness[,3],forest.site=forest.site[,3],field.site=field.site[,3],
                 micro=micro[,3],cam.height=cam.height[,3],walktest=walktest[,3])
# Initial values
zst4 <- apply(y4,1,max,na.rm=T)      #inits for presence (z)
inits4 <- function()list(z=zst4, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits4.1 <- function()list(z=zst4, psi0=0.5)
# full model
outJ.mooseHM.L <- jags(win.data4, inits4, params, "occ-HM-L.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.mooseHM.L,file="outJ-mooseHM-L.RData")
# null model
outJ.mooseHM.L.nm <- jags(win.data4, inits4.1, params, "occ-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.mooseHM.L.nm,file="outJ-mooseHM-L-nm.RData")

## H cameras ##
y5 <- y[,2,]
win.data5 <- list(y=y5, n.site=dim(y5)[1],n.days=dim(y5)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,2],forest.site=forest.site[,2],field.site=field.site[,2],
                  micro=micro[,2],cam.height=cam.height[,2],walktest=walktest[,2])
# Initial values
zst5 <- apply(y5,1,max,na.rm=T)      #inits for presence (z)
inits5 <- function()list(z=zst5, psi0=0.5)
# null model
outJ.mooseHM.H.nm <- jags(win.data5, inits5, params, "occ-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.mooseHM.H.nm,file="outJ-mooseHM-H-nm.RData")

detach(data); rm("data")

#### Roe deer ####

# Load data
load("MSocc-Roe.deerHM.RData")
attach(data) 

## All cameras ##
win.data <- list(y=y, n.site=dim(y)[1],n.cts=dim(y)[2],n.days=dim(y)[3],
                 elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                 target.site=target.site,ruggedness=ruggedness,forest.site=forest.site,field.site=field.site,
                 target.cam=target.cam,micro=micro,cam.height=cam.height,walktest=walktest)
# Initial values
zst <- apply(y,1,max,na.rm=T)      #inits for presence (z)
ast <- apply(y,c(1,2),max,na.rm=T) #inits for availability (a)
inits <- function()list(z=zst, a=ast, psi0=0.5, beta.lpsi=c(0,0,0))
inits.1 <- function()list(z=zst, a=ast, psi0=0.5)
# full model
outJ.roedeerHM.all <- jags(win.data, inits, params, "MSocc-HM-all.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.roedeerHM.all,file="outJ-roedeerHM-all.RData")
# null model
outJ.roedeerHM.all.nm <- jags(win.data, inits.1, params, "MSocc-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.roedeerHM.all.nm,file="outJ-roedeerHM-all-nm.RData")
# final model
outJ.roedeerHM.all.fm <- jags(win.data, inits.1, params, "MSocc-HM-all-roedeer.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.roedeerHM.all.fm,file="outJ-roedeerHM-all-fm.RData")

## T+L cameras ##
y1 <- y[,c(1,3),]
win.data1 <- list(y=y1, n.site=dim(y1)[1],n.cts=dim(y1)[2],n.days=dim(y1)[3],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  target.site=target.site[,c(1,3)],ruggedness=ruggedness[,c(1,3)],
                  forest.site=forest.site[,c(1,3)],field.site=field.site[,c(1,3)],
                  target.cam=target.cam[,c(1,3)],micro=micro[,c(1,3)],cam.height=cam.height[,c(1,3)],walktest=walktest[,c(1,3)])
# Initial values
zst1 <- apply(y1,1,max,na.rm=T)      #inits for presence (z)
ast1 <- apply(y1,c(1,2),max,na.rm=T) #inits for availability (a)
inits1 <- function()list(z=zst1, a=ast1, psi0=0.5, beta.lpsi=c(0,0,0))
inits1.1 <- function()list(z=zst1, a=ast1, psi0=0.5)
# full model
outJ.roedeerHM.TL <- jags(win.data1, inits1, params, "MSocc-HM-TL.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.roedeerHM.TL,file="outJ-roedeerHM-TL.RData")
# null model
outJ.roedeerHM.TL.nm <- jags(win.data1, inits1.1, params, "MSocc-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.roedeerHM.TL.nm,file="outJ-roedeerHM-TL-nm.RData")
# final model
outJ.roedeerHM.TL.fm <- jags(win.data1, inits1.1, params, "MSocc-HM-TL-roedeer.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.roedeerHM.TL.fm,file="outJ-roedeerHM-TL-fm.RData")

## T+H cameras ##
y2 <- y[,c(1:2),]
win.data2 <- list(y=y2, n.site=dim(y2)[1],n.cts=dim(y2)[2],n.days=dim(y2)[3],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  target.site=target.site,ruggedness=ruggedness,forest.site=forest.site,field.site=field.site,
                  target.cam=target.cam,micro=micro,cam.height=cam.height,walktest=walktest)
# Initial values
zst2 <- apply(y2,1,max,na.rm=T)      #inits for presence (z)
ast2 <- apply(y2,c(1,2),max,na.rm=T) #inits for availability (a)
inits2 <- function()list(z=zst2, a=ast2, psi0=0.5, beta.lpsi=c(0,0,0))
inits2.1 <- function()list(z=zst2, a=ast2, psi0=0.5)
# full model
outJ.roedeerHM.TH <- jags(win.data2, inits2, params, "MSocc-HM-TH.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.roedeerHM.TH,file="outJ-roedeerHM-TH.RData")
# null model
outJ.roedeerHM.TH.nm <- jags(win.data2, inits2.1, params, "MSocc-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.roedeerHM.TH.nm,file="outJ-roedeerHM-TH-nm.RData")
# final model
outJ.roedeerHM.TH.fm <- jags(win.data2, inits2.1, params, "MSocc-HM-TH-roedeer.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.roedeerHM.TH.fm,file="outJ-roedeerHM-TH-fm.RData")

## T cameras ##
y3 <- y[,1,]
win.data3 <- list(y=y3, n.site=dim(y3)[1],n.days=dim(y3)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,1],forest.site=forest.site[,1],field.site=field.site[,1],
                  micro=micro[,1],cam.height=cam.height[,1],walktest=walktest[,1])
# Initial values
zst3 <- apply(y3,1,max,na.rm=T)      #inits for presence (z)
inits3 <- function()list(z=zst3, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits3.1 <- function()list(z=zst3, psi0=0.5)
# full model
outJ.roedeerHM.T <- jags(win.data3, inits3, params, "occ-HM-T.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.roedeerHM.T,file="outJ-roedeerHM-T.RData")
# null model
outJ.roedeerHM.T.nm <- jags(win.data3, inits3.1, params, "occ-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.roedeerHM.T.nm,file="outJ-roedeerHM-T-nm.RData")
# full model
outJ.roedeerHM.T.fm <- jags(win.data3, inits3.1, params, "occ-HM-T-roedeer.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.roedeerHM.T.fm,file="outJ-roedeerHM-T-fm.RData")

## L cameras ##
y4 <- y[,3,]
win.data4 <- list(y=y4, n.site=dim(y4)[1],n.days=dim(y4)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,3],forest.site=forest.site[,3],field.site=field.site[,3],
                  micro=micro[,3],cam.height=cam.height[,3],walktest=walktest[,3])
# Initial values
zst4 <- apply(y4,1,max,na.rm=T)      #inits for presence (z)
inits4 <- function()list(z=zst4, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits4.1 <- function()list(z=zst4, psi0=0.5)
# full model
outJ.roedeerHM.L <- jags(win.data4, inits4, params, "occ-HM-L.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.roedeerHM.L,file="outJ-roedeerHM-L.RData")
# null model
outJ.roedeerHM.L.nm <- jags(win.data4, inits4.1, params, "occ-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.roedeerHM.L.nm,file="outJ-roedeerHM-L-nm.RData")

## H cameras ##
y5 <- y[,2,]
win.data5 <- list(y=y5, n.site=dim(y5)[1],n.days=dim(y5)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,2],forest.site=forest.site[,2],field.site=field.site[,2],
                  micro=micro[,2],cam.height=cam.height[,2],walktest=walktest[,2])
# Initial values
zst5 <- apply(y5,1,max,na.rm=T)      #inits for presence (z)
inits5 <- function()list(z=zst5, psi0=0.5)
# null model
outJ.roedeerHM.H.nm <- jags(win.data5, inits5, params, "occ-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.roedeerHM.H.nm,file="outJ-roedeerHM-H-nm.RData")

detach(data); rm("data")

#### Mountain hare ####

# Load data
load("MSocc-HareHM.RData")
attach(data) 

## All cameras ##
win.data <- list(y=y, n.site=dim(y)[1],n.cts=dim(y)[2],n.days=dim(y)[3],
                 elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                 target.site=target.site,ruggedness=ruggedness,forest.site=forest.site,field.site=field.site,
                 target.cam=target.cam,micro=micro,cam.height=cam.height,walktest=walktest)
# Initial values
zst <- apply(y,1,max,na.rm=T)      #inits for presence (z)
ast <- apply(y,c(1,2),max,na.rm=T) #inits for availability (a)
inits <- function()list(z=zst, a=ast, psi0=0.5, beta.lpsi=c(0,0,0))
inits.1 <- function()list(z=zst, a=ast, psi0=0.5)
# full model
outJ.hareHM.all <- jags(win.data, inits, params, "MSocc-HM-all.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.hareHM.all,file="outJ-hareHM-all.RData")
# null model
outJ.hareHM.all.nm <- jags(win.data, inits.1, params, "MSocc-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.hareHM.all.nm,file="outJ-hareHM-all-nm.RData")
# final model
outJ.hareHM.all.fm <- jags(win.data, inits.1, params, "MSocc-HM-all-hare.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.hareHM.all.fm,file="outJ-hareHM-all-fm.RData")

## T+L cameras ##
y1 <- y[,c(1,3),]
win.data1 <- list(y=y1, n.site=dim(y1)[1],n.cts=dim(y1)[2],n.days=dim(y1)[3],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  target.site=target.site[,c(1,3)],ruggedness=ruggedness[,c(1,3)],
                  forest.site=forest.site[,c(1,3)],field.site=field.site[,c(1,3)],
                  target.cam=target.cam[,c(1,3)],micro=micro[,c(1,3)],cam.height=cam.height[,c(1,3)],walktest=walktest[,c(1,3)])
# Initial values
zst1 <- apply(y1,1,max,na.rm=T)      #inits for presence (z)
ast1 <- apply(y1,c(1,2),max,na.rm=T) #inits for availability (a)
inits1 <- function()list(z=zst1, a=ast1, psi0=0.5, beta.lpsi=c(0,0,0))
inits1.1 <- function()list(z=zst1, a=ast1, psi0=0.5)
# full model
outJ.hareHM.TL <- jags(win.data1, inits1, params, "MSocc-HM-TL.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.hareHM.TL,file="outJ-hareHM-TL.RData")
# null model
outJ.hareHM.TL.nm <- jags(win.data1, inits1.1, params, "MSocc-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.hareHM.TL.nm,file="outJ-hareHM-TL-nm.RData")
# final model
outJ.hareHM.TL.fm <- jags(win.data1, inits1.1, params, "MSocc-HM-TL-hare.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.hareHM.TL.fm,file="outJ-hareHM-TL-fm.RData")

## T+H cameras ##
y2 <- y[,c(1:2),]
win.data2 <- list(y=y2, n.site=dim(y2)[1],n.cts=dim(y2)[2],n.days=dim(y2)[3],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  target.site=target.site,ruggedness=ruggedness,forest.site=forest.site,field.site=field.site,
                  target.cam=target.cam,micro=micro,cam.height=cam.height,walktest=walktest)
# Initial values
zst2 <- apply(y2,1,max,na.rm=T)      #inits for presence (z)
ast2 <- apply(y2,c(1,2),max,na.rm=T) #inits for availability (a)
inits2 <- function()list(z=zst2, a=ast2, psi0=0.5, beta.lpsi=c(0,0,0))
inits2.1 <- function()list(z=zst2, a=ast2, psi0=0.5)
# full model
outJ.hareHM.TH <- jags(win.data2, inits2, params, "MSocc-HM-TH.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.hareHM.TH,file="outJ-hareHM-TH.RData")
# null model
outJ.hareHM.TH.nm <- jags(win.data2, inits2.1, params, "MSocc-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.hareHM.TH.nm,file="outJ-hareHM-TH-nm.RData")
# final model
outJ.hareHM.TH.fm <- jags(win.data2, inits2.1, params, "MSocc-HM-TH-hare.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.hareHM.TH.fm,file="outJ-hareHM-TH-fm.RData")

## T cameras ##
y3 <- y[,1,]
win.data3 <- list(y=y3, n.site=dim(y3)[1],n.days=dim(y3)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,1],forest.site=forest.site[,1],field.site=field.site[,1],
                  micro=micro[,1],cam.height=cam.height[,1],walktest=walktest[,1])
# Initial values
zst3 <- apply(y3,1,max,na.rm=T)      #inits for presence (z)
inits3 <- function()list(z=zst3, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits3.1 <- function()list(z=zst3, psi0=0.5)
# full model
outJ.hareHM.T <- jags(win.data3, inits3, params, "occ-HM-T.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.hareHM.T,file="outJ-hareHM-T.RData")
# null model
outJ.hareHM.T.nm <- jags(win.data3, inits3.1, params, "occ-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.hareHM.T.nm,file="outJ-hareHM-T-nm.RData")

## L cameras ##
y4 <- y[,3,]
win.data4 <- list(y=y4, n.site=dim(y4)[1],n.days=dim(y4)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,3],forest.site=forest.site[,3],field.site=field.site[,3],
                  micro=micro[,3],cam.height=cam.height[,3],walktest=walktest[,3])
# Initial values
zst4 <- apply(y4,1,max,na.rm=T)      #inits for presence (z)
inits4 <- function()list(z=zst4, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits4.1 <- function()list(z=zst4, psi0=0.5)
# full model
outJ.hareHM.L <- jags(win.data4, inits4, params, "occ-HM-L.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.hareHM.L,file="outJ-hareHM-L.RData")
# null model
outJ.hareHM.L.nm <- jags(win.data4, inits4.1, params, "occ-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.hareHM.L.nm,file="outJ-hareHM-L-nm.RData")
# full model
outJ.hareHM.L.fm <- jags(win.data4, inits4.1, params, "occ-HM-L-hare.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.hareHM.L.fm,file="outJ-hareHM-L-fm.RData")

## H cameras ##
y5 <- y[,2,]
win.data5 <- list(y=y5, n.site=dim(y5)[1],n.days=dim(y5)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,2],forest.site=forest.site[,2],field.site=field.site[,2],
                  micro=micro[,2],cam.height=cam.height[,2],walktest=walktest[,2])
# Initial values
zst5 <- apply(y5,1,max,na.rm=T)      #inits for presence (z)
inits5 <- function()list(z=zst5, psi0=0.5)
# null model
outJ.hareHM.H.nm <- jags(win.data5, inits5, params, "occ-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.hareHM.H.nm,file="outJ-hareHM-H-nm.RData")

detach(data); rm("data")

#### Red fox ####

# Load data
load("MSocc-Red.foxHM.RData")
attach(data) 

## All cameras ##
win.data <- list(y=y, n.site=dim(y)[1],n.cts=dim(y)[2],n.days=dim(y)[3],
                 elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                 target.site=target.site,ruggedness=ruggedness,forest.site=forest.site,field.site=field.site,
                 target.cam=target.cam,micro=micro,cam.height=cam.height,walktest=walktest)
# Initial values
zst <- apply(y,1,max,na.rm=T)      #inits for presence (z)
ast <- apply(y,c(1,2),max,na.rm=T) #inits for availability (a)
inits <- function()list(z=zst, a=ast, psi0=0.5, beta.lpsi=c(0,0,0))
inits.1 <- function()list(z=zst, a=ast, psi0=0.5)
# full model
outJ.foxHM.all <- jags(win.data, inits, params, "MSocc-HM-all.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.foxHM.all,file="outJ-foxHM-all.RData")
# null model
outJ.foxHM.all.nm <- jags(win.data, inits.1, params, "MSocc-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.foxHM.all.nm,file="outJ-foxHM-all-nm.RData")
# final model
outJ.foxHM.all.fm <- jags(win.data, inits.1, params, "MSocc-HM-all-fox.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.foxHM.all.fm,file="outJ-foxHM-all-fm.RData")

## T+L cameras ##
y1 <- y[,c(1,3),]
win.data1 <- list(y=y1, n.site=dim(y1)[1],n.cts=dim(y1)[2],n.days=dim(y1)[3],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  target.site=target.site[,c(1,3)],ruggedness=ruggedness[,c(1,3)],
                  forest.site=forest.site[,c(1,3)],field.site=field.site[,c(1,3)],
                  target.cam=target.cam[,c(1,3)],micro=micro[,c(1,3)],cam.height=cam.height[,c(1,3)],walktest=walktest[,c(1,3)])
# Initial values
zst1 <- apply(y1,1,max,na.rm=T)      #inits for presence (z)
ast1 <- apply(y1,c(1,2),max,na.rm=T) #inits for availability (a)
inits1 <- function()list(z=zst1, a=ast1, psi0=0.5, beta.lpsi=c(0,0,0))
inits1.1 <- function()list(z=zst1, a=ast1, psi0=0.5)
# full model
outJ.foxHM.TL <- jags(win.data1, inits1, params, "MSocc-HM-TL.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.foxHM.TL,file="outJ-foxHM-TL.RData")
# null model
outJ.foxHM.TL.nm <- jags(win.data1, inits1.1, params, "MSocc-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.foxHM.TL.nm,file="outJ-foxHM-TL-nm.RData")
# final model
outJ.foxHM.TL.fm <- jags(win.data1, inits1.1, params, "MSocc-HM-TL-fox.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.foxHM.TL.fm,file="outJ-foxHM-TL-fm.RData")

## T+H cameras ##
y2 <- y[,c(1:2),]
win.data2 <- list(y=y2, n.site=dim(y2)[1],n.cts=dim(y2)[2],n.days=dim(y2)[3],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  target.site=target.site,ruggedness=ruggedness,forest.site=forest.site,field.site=field.site,
                  target.cam=target.cam,micro=micro,cam.height=cam.height,walktest=walktest)
# Initial values
zst2 <- apply(y2,1,max,na.rm=T)      #inits for presence (z)
ast2 <- apply(y2,c(1,2),max,na.rm=T) #inits for availability (a)
inits2 <- function()list(z=zst2, a=ast2, psi0=0.5, beta.lpsi=c(0,0,0))
inits2.1 <- function()list(z=zst2, a=ast2, psi0=0.5)
# full model
outJ.foxHM.TH <- jags(win.data2, inits2, params, "MSocc-HM-TH.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.foxHM.TH,file="outJ-foxHM-TH.RData")
# null model
outJ.foxHM.TH.nm <- jags(win.data2, inits2.1, params, "MSocc-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.foxHM.TH.nm,file="outJ-foxHM-TH-nm.RData")
# final model
outJ.foxHM.TH.fm <- jags(win.data2, inits2.1, params, "MSocc-HM-TH-fox.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.foxHM.TH.fm,file="outJ-foxHM-TH-fm.RData")

## T cameras ##
y3 <- y[,1,]
win.data3 <- list(y=y3, n.site=dim(y3)[1],n.days=dim(y3)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,1],forest.site=forest.site[,1],field.site=field.site[,1],
                  micro=micro[,1],cam.height=cam.height[,1],walktest=walktest[,1])
# Initial values
zst3 <- apply(y3,1,max,na.rm=T)      #inits for presence (z)
inits3 <- function()list(z=zst3, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits3.1 <- function()list(z=zst3, psi0=0.5)
# full model
outJ.foxHM.T <- jags(win.data3, inits3, params, "occ-HM-T.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.foxHM.T,file="outJ-foxHM-T.RData")
# null model
outJ.foxHM.T.nm <- jags(win.data3, inits3.1, params, "occ-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.foxHM.T.nm,file="outJ-foxHM-T-nm.RData")
# final model
outJ.foxHM.T.fm <- jags(win.data3, inits3.1, params, "occ-HM-T-fox.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.foxHM.T.fm,file="outJ-foxHM-T-fm.RData")

## L cameras ##
y4 <- y[,3,]
win.data4 <- list(y=y4, n.site=dim(y4)[1],n.days=dim(y4)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,3],forest.site=forest.site[,3],field.site=field.site[,3],
                  micro=micro[,3],cam.height=cam.height[,3],walktest=walktest[,3])
# Initial values
zst4 <- apply(y4,1,max,na.rm=T)      #inits for presence (z)
inits4 <- function()list(z=zst4, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits4.1 <- function()list(z=zst4, psi0=0.5)
inits4.2 <- function()list(z=zst4, psi0=0.5, beta.lpsi=c(0,0))
# full model
outJ.foxHM.L <- jags(win.data4, inits4, params, "occ-HM-L.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.foxHM.L,file="outJ-foxHM-L.RData")
# null model
outJ.foxHM.L.nm <- jags(win.data4, inits4.1, params, "occ-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.foxHM.L.nm,file="outJ-foxHM-L-nm.RData")
# full model
outJ.foxHM.L.fm <- jags(win.data4, inits4.2, params, "occ-HM-L-fox.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.foxHM.L.fm,file="outJ-foxHM-L-fm.RData")

## H cameras ##
y5 <- y[,2,]
win.data5 <- list(y=y5, n.site=dim(y5)[1],n.days=dim(y5)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,2],forest.site=forest.site[,2],field.site=field.site[,2],
                  micro=micro[,2],cam.height=cam.height[,2],walktest=walktest[,2])
# Initial values
zst5 <- apply(y5,1,max,na.rm=T)      #inits for presence (z)
inits5 <- function()list(z=zst5, psi0=0.5)
# null model
outJ.foxHM.H.nm <- jags(win.data5, inits5, params, "occ-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.foxHM.H.nm,file="outJ-foxHM-H-nm.RData")

detach(data); rm("data")

#### Wolf ####

# Load data
load("MSocc-WolfHM.RData")
attach(data) 

## T cameras ##
y1 <- y[,1,]
win.data <- list(y=y1, n.site=dim(y1)[1],n.days=dim(y1)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,1],forest.site=forest.site[,1],field.site=field.site[,1],
                  micro=micro[,1],cam.height=cam.height[,1],walktest=walktest[,1])
# Initial values
zst <- apply(y1,1,max,na.rm=T)      #inits for presence (z)
inits <- function()list(z=zst, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits.1 <- function()list(z=zst, psi0=0.5)
# full model
outJ.wolfHM.T <- jags(win.data, inits, params, "occ-HM-T.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.wolfHM.T,file="outJ-wolfHM-T.RData")
# null model
outJ.wolfHM.T.nm <- jags(win.data, inits.1, params, "occ-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.wolfHM.T.nm,file="outJ-wolfHM-T-nm.RData")
# final model
outJ.wolfHM.T.fm <- jags(win.data, inits.1, params, "occ-HM-T-wolf.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.wolfHM.T.fm,file="outJ-wolfHM-T-fm.RData")

detach(data); rm("data")

#### Lynx ####

# Load data
load("MSocc-LynxHM.RData")
attach(data) 

## T cameras ##
y1 <- y[,1,]
win.data <- list(y=y1, n.site=dim(y1)[1],n.days=dim(y1)[2],
                 elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                 ruggedness=ruggedness[,1],forest.site=forest.site[,1],field.site=field.site[,1],
                 micro=micro[,1],cam.height=cam.height[,1],walktest=walktest[,1])
# Initial values
zst <- apply(y1,1,max,na.rm=T)      #inits for presence (z)
inits <- function()list(z=zst, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits.1 <- function()list(z=zst, psi0=0.5)
# full model
outJ.lynxHM.T <- jags(win.data, inits, params, "occ-HM-T.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.lynxHM.T,file="outJ-lynxHM-T.RData")
# null model
outJ.lynxHM.T.nm <- jags(win.data, inits.1, params, "occ-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.lynxHM.T.nm,file="outJ-lynxHM-T-nm.RData")

detach(data); rm("data")

#### Badger ####

# Load data
load("MSocc-BadgerHM.RData")
attach(data) 

## T cameras ##
y1 <- y[,1,]
win.data <- list(y=y1, n.site=dim(y1)[1],n.days=dim(y1)[2],
                 elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                 ruggedness=ruggedness[,1],forest.site=forest.site[,1],field.site=field.site[,1],
                 micro=micro[,1],cam.height=cam.height[,1],walktest=walktest[,1])
# Initial values
zst <- apply(y1,1,max,na.rm=T)      #inits for presence (z)
inits <- function()list(z=zst, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits.1 <- function()list(z=zst, psi0=0.5)
inits.2 <- function()list(z=zst, psi0=0.5, beta.lpsi=0)
# full model
outJ.badgerHM.T <- jags(win.data, inits, params, "occ-HM-T.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.badgerHM.T,file="outJ-badgerHM-T.RData")
# null model
outJ.badgerHM.T.nm <- jags(win.data, inits.1, params, "occ-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.badgerHM.T.nm,file="outJ-badgerHM-T-nm.RData")
# final model
outJ.badgerHM.T.fm <- jags(win.data, inits.2, params, "occ-HM-T-badger.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.badgerHM.T.fm,file="outJ-badgerHM-T-fm.RData")

detach(data); rm("data")

#### Red squirrel ####

# Load data
load("MSocc-Red.squirrelHM.RData")
attach(data) 

## T cameras ##
y1 <- y[,1,]
win.data <- list(y=y1, n.site=dim(y1)[1],n.days=dim(y1)[2],
                 elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                 ruggedness=ruggedness[,1],forest.site=forest.site[,1],field.site=field.site[,1],
                 micro=micro[,1],cam.height=cam.height[,1],walktest=walktest[,1])
# Initial values
zst <- apply(y1,1,max,na.rm=T)      #inits for presence (z)
inits <- function()list(z=zst, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits.1 <- function()list(z=zst, psi0=0.5)
inits.2 <- function()list(z=zst, psi0=0.5, beta.lpsi=c(0,0,0))
# full model
outJ.squirrelHM.T <- jags(win.data, inits, params, "occ-HM-T.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.squirrelHM.T,file="outJ-squirrelHM-T.RData")
# null model
outJ.squirrelHM.T.nm <- jags(win.data, inits.1, params, "occ-HMnull.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.squirrelHM.T.nm,file="outJ-squirrelHM-T-nm.RData")
# final model
outJ.squirrelHM.T.fm <- jags(win.data, inits.2, params, "occ-HM-T-squirrel.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.squirrelHM.T.fm,file="outJ-squirrelHM-T-fm.RData")

detach(data); rm("data")