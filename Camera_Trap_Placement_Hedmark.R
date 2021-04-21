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
    
    for(t in 1:3){
    int.p[t] <- logit(p0[t])          #Intercepts for detection probability
    p0[t] ~ dunif(0,1)     
    }
    
    #Slope of site covariate
    for(t in 1:3){
    beta.lpsi[t] ~ dnorm(0,0.01)
    }        
    
    #Slopes of ct covariates
    for(t in 1:4){        
    beta.ltheta[t] ~ dnorm(0,0.01)
    }
    
    #Slope of p covariates
    for(t in 1:5){
    beta.lp[t] ~ dnorm(0,0.01)
    }

    #Inclusion parameter priors
    for(c in 1:7){
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
    logit(p[i,j,k]) <- int.p[j] + 
    beta.lp[1] * equals(micro[i,j],2) +     # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i,j],3) +     # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i,j],4) +     # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i,j] +          # height of the camera from the ground
    beta.lp[5] * walktest[i,j]
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
    
    for(t in 1:2){
    int.p[t] <- logit(p0[t])          #Intercepts for detection probability
    p0[t] ~ dunif(0,1)
    }
    
    #Slope of site covariate
    for(t in 1:3){
    beta.lpsi[t] ~ dnorm(0,0.01)
    }        
    
    #Slopes of ct covariates
    for(t in 1:4){        
    beta.ltheta[t] ~ dnorm(0,0.01)
    }
    
    #Slope of p covariates
    for(t in 1:5){
    beta.lp[t] ~ dnorm(0,0.01)
    }

    #Inclusion parameter priors
    for(c in 1:7){
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
    w[4] * beta.ltheta[1] * target.site[i,j] +          # if camera is targeted (1) or not (0)
    w[5] * beta.ltheta[2] * ruggedness[i,j] +           # slope for ruggedness
    w[6] * beta.ltheta[3] * forest.site[i,j] +          # forest cover 1km2 around camera
    w[7] * beta.ltheta[4] * field.site[i,j]             # field cover 1km2 around camera

    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p[j] + 
    beta.lp[1] * equals(micro[i,j],2) +   # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i,j],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i,j],4) +   # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i,j] +
    beta.lp[5] * walktest[i,j]
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
    
    for(t in 1:2){
    int.p[t] <- logit(p0[t])          #Intercepts for detection probability
    p0[t] ~ dunif(0,1)
    }
    
    #Slope of site covariate
    for(t in 1:3){
    beta.lpsi[t] ~ dnorm(0,0.01)
    }        
    
    
    #Slopes of ct covariates
    for(t in 1:3){        
    beta.ltheta[t] ~ dnorm(0,0.01)
    }
    
    #Slope of p covariates
    for(t in 1:5){
    beta.lp[t] ~ dnorm(0,0.01)
    }

    #Inclusion parameter priors
    for(c in 1:6){
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
    logit(p[i,j,k]) <- int.p[j] + 
    beta.lp[1] * equals(micro[i,j],2) +   # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i,j],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i,j],4) +   # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i,j] +
    beta.lp[5] * walktest[i,j]
    }
    }
    }
    } #End Model
    ")

## Multi-scale model for all habitat patch and lanscape cameras ##

cat(file="MSocc-HM-HL.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    for(t in 1:2){
    int.p[t] <- logit(p0[t])          #Intercepts for detection probability
    p0[t] ~ dunif(0,1)
    }
    
    #Slope of site covariate
    for(t in 1:3){
    beta.lpsi[t] ~ dnorm(0,0.01)
    }        
    
    #Slopes of ct covariates
    for(t in 1:4){        
    beta.ltheta[t] ~ dnorm(0,0.01)
    }
    
    #Slope of p covariates
    for(t in 1:4){
    beta.lp[t] ~ dnorm(0,0.01)
    }

    #Inclusion parameter priors
    for(c in 1:7){
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
    w[5] * beta.ltheta[2] * ruggedness[i,j] +             # slope for ruggedness
    w[6] * beta.ltheta[3] * forest.site[i,j] +            # forest cover 1km2 around camera
    w[7] * beta.ltheta[4] * field.site[i,j]               # field cover 1km2 around camera

    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p[j] + 
    beta.lp[1] * equals(micro[i,j],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[2] * equals(micro[i,j],4) +   # cliff compared to no feature (intercept)
    beta.lp[3] * cam.height[i,j] +
    beta.lp[4] * walktest[i,j]
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
    for(c in 1:6){
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
    beta.lp[1] * equals(micro[i],2) +   # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i],4) +   # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i] +
    beta.lp[5] * walktest[i]
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
    for(c in 1:6){
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
    beta.lp[1] * equals(micro[i],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[2] * equals(micro[i],4) +    # cliff compared to no feature (intercept)
    beta.lp[3] * cam.height[i] +
    beta.lp[4] * walktest[i]
    }
    }
    } #End Model
    ")

#### Null models ####

## Three-camera null model ##

# model
cat(file="MSocc-HM-all-null.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    for(t in 1:3){
    int.p[t] <- logit(p0[t])          #Intercepts for detection probability
    p0[t] ~ dunif(0,1)     
    }
    
    #Slope of p covariates
    for(t in 1:5){
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
    logit(theta[i,j]) <- int.theta
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p[j] + 
    beta.lp[1] * equals(micro[i,j],2) +     # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i,j],3) +     # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i,j],4) +     # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i,j] +          # height of the camera from the ground
    beta.lp[5] * walktest[i,j]
    }
    }
    }
    } #End Model
    ")

## Two-camera null model ##

# model
cat(file="MSocc-HM-two-null.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    for(t in 1:2){
    int.p[t] <- logit(p0[t])          #Intercepts for detection probability
    p0[t] ~ dunif(0,1)     
    }
    
    #Slope of p covariates
    for(t in 1:5){
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
    logit(theta[i,j]) <- int.theta
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p[j] + 
    beta.lp[1] * equals(micro[i,j],2) +     # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i,j],3) +     # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i,j],4) +     # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i,j] +          # height of the camera from the ground
    beta.lp[5] * walktest[i,j]
    }
    }
    }
    } #End Model
    ")

## Single targeted-camera null model ##

# model
cat(file="occ-HM-T-null.txt", 
    "
    model{
    
   # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)
    
    #Slope of p covariates
    for(t in 1:5){
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
    beta.lp[3] * equals(micro[i],4) +   # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i] +
    beta.lp[5] * walktest[i]
    }
    }
    } #End Model
    ")

## Single landscape-camera null model ##

# model
cat(file="occ-HM-L-null.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1) 
    
    #Slope for p covariates
    for(t in 1:4){
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
    beta.lp[3] * cam.height[i] +
    beta.lp[4] * walktest[i]
    }
    }
    } #End Model
    ")

#### Final models all three cameras ####

## Moose ##

# model
cat(file="MSocc-HM-all-moose.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    for(t in 1:3){
    int.p[t] <- logit(p0[t])          #Intercepts for detection probability
    p0[t] ~ dunif(0,1)     
    }
    
    #Slope of theta covariate
    beta.ltheta ~ dnorm(0,0.01)
    
    #Slope of p covariates
    for(t in 1:5){
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
    beta.ltheta * field.site[i,j]            # field cover 1km2 around camera
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p[j] + 
    beta.lp[1] * equals(micro[i,j],2) +     # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i,j],3) +     # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i,j],4) +     # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i,j] +          # height of the camera from the ground
    beta.lp[5] * walktest[i,j]
    }
    }
    }
    } #End Model
    ")

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
    
    for(t in 1:3){
    int.p[t] <- logit(p0[t])          #Intercepts for detection probability
    p0[t] ~ dunif(0,1)     
    }
    
    #Slope of theta covariate
    beta.ltheta ~ dnorm(0,0.01)
    
    #Slope of p covariates
    for(t in 1:5){
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
    beta.ltheta * field.site[i,j]            # field cover 1km2 around camera
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p[j] + 
    beta.lp[1] * equals(micro[i,j],2) +     # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i,j],3) +     # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i,j],4) +     # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i,j] +          # height of the camera from the ground
    beta.lp[5] * walktest[i,j]
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
    
    for(t in 1:3){
    int.p[t] <- logit(p0[t])                #Intercepts for detection probability
    p0[t] ~ dunif(0,1)     
    }
    
    #Slope of theta covariate
    beta.ltheta ~ dnorm(0,0.01)
    
    #Slope of p covariates
    for(t in 1:5){
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
    beta.ltheta * target.site[i,j]            # if the cameras were in a targeted patch (1) or not (0)
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p[j] + 
    beta.lp[1] * equals(micro[i,j],2) +     # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i,j],3) +     # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i,j],4) +     # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i,j] +          # height of the camera from the ground
    beta.lp[5] * walktest[i,j]
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
    
    for(t in 1:2){
    int.p[t] <- logit(p0[t])          #Intercepts for detection probability
    p0[t] ~ dunif(0,1)     
    }
    
    #Slope of theta covariate
    beta.ltheta ~ dnorm(0,0.01)
    
    #Slope of p covariates
    for(t in 1:5){
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
    logit(p[i,j,k]) <- int.p[j] +
    beta.lp[1] * equals(micro[i,j],2) +   # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i,j],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i,j],4) +   # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i,j] +
    beta.lp[5] * walktest[i,j]    
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
    
    for(t in 1:2){
    int.p[t] <- logit(p0[t])          #Intercepts for detection probability
    p0[t] ~ dunif(0,1)     
    }
    
    #Slope of theta covariate
    beta.ltheta ~ dnorm(0,0.01)
    
    #Slope of p covariates
    for(t in 1:5){
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
    beta.ltheta * target.site[i,j]
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- int.p[j] +
    beta.lp[1] * equals(micro[i,j],2) +   # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i,j],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i,j],4) +   # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i,j] +
    beta.lp[5] * walktest[i,j]    
    }
    }
    }
    } #End Model
    ")


#### Final models all targeted + habitat patch cameras ####

## Moose ##

# model
cat(file="MSocc-HM-TH-moose.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    for(t in 1:2){
    int.p[t] <- logit(p0[t])          #Intercepts for detection probability
    p0[t] ~ dunif(0,1)     
    }
    
    #Slope of theta covariate
    beta.ltheta ~ dnorm(0,0.01)
    
    #Slope of p covariates
    for(t in 1:5){
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
    logit(p[i,j,k]) <- int.p[j] +
    beta.lp[1] * equals(micro[i,j],2) +   # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i,j],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i,j],4) +   # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i,j] +
    beta.lp[5] * walktest[i,j]
    }
    }
    }
    } #End Model
    ")

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
    
    for(t in 1:2){
    int.p[t] <- logit(p0[t])                #Intercepts for detection probability
    p0[t] ~ dunif(0,1)     
    }
    
    #Slope of theta covariate
    beta.ltheta ~ dnorm(0,0.01)
    
    #Slope of p covariates
    for(t in 1:5){
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
    logit(p[i,j,k]) <- int.p[j] +
    beta.lp[1] * equals(micro[i,j],2) +   # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i,j],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i,j],4) +   # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i,j] +
    beta.lp[5] * walktest[i,j]
    }
    }
    }
    } #End Model
    ")

#### Final models all habitat patch + landscape cameras ####

## Roe deer ##

# model
cat(file="MSocc-HM-HL-roedeer.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.theta <- logit(theta0)        #Intercepts for availability probability
    theta0 ~ dunif(0,1) 
    
    for(t in 1:2){
    int.p[t] <- logit(p0[t])          #Intercepts for detection probability
    p0[t] ~ dunif(0,1)     
    }
    
    #Slope of theta covariate
    beta.ltheta ~ dnorm(0,0.01)
    
    #Slope of p covariates
    for(t in 1:4){
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
    logit(p[i,j,k]) <- int.p[j] +
    beta.lp[1] * equals(micro[i,j],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[2] * equals(micro[i,j],4) +   # cliff compared to no feature (intercept)
    beta.lp[3] * cam.height[i,j] +
    beta.lp[4] * walktest[i,j]    
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
    
    #Slope of site covariate
    beta.lpsi ~ dnorm(0,0.01)
    
    #Slope of p covariates
    for(t in 1:5){
    beta.lp[t] ~ dnorm(0,0.01)
    }
    
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
    logit(p[i,k]) <- int.p +
    beta.lp[1] * equals(micro[i],2) +   # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i],4) +   # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i] +
    beta.lp[5] * walktest[i]
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
    for(t in 1:5){
    beta.lp[t] ~ dnorm(0,0.01)
    }
    
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
    beta.lp[1] * equals(micro[i],2) +   # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i],4) +   # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i] +
    beta.lp[5] * walktest[i]
    }
    }
    } #End Model
    ")

## Lynx ##

# model
cat(file="MSocc-HM-T-lynx.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1)
    
    #Slope of site covariate
    beta.lpsi ~ dnorm(0,0.01)
    
    #Slope of p covariates
    for(t in 1:5){
    beta.lp[t] ~ dnorm(0,0.01)
    }
    
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
    logit(p[i,k]) <- int.p +
    beta.lp[1] * equals(micro[i],2) +   # forest road compared to no feature (intercept)
    beta.lp[2] * equals(micro[i],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[3] * equals(micro[i],4) +   # cliff compared to no feature (intercept)
    beta.lp[4] * cam.height[i] +
    beta.lp[5] * walktest[i]
    }
    }
    } #End Model
    ")

#### Final models landscape cameras ####

## Moose ##

# model
cat(file="MSocc-HM-L-moose.txt", 
    "
    model{
    
    # Priors and model for params
    
    int.psi <- logit(psi0)            #Intercept for occupancy probability
    psi0 ~ dunif(0,1)        
    
    int.p <- logit(p0)                #Intercepts for detection probability
    p0 ~ dunif(0,1) 
    
    #Slope for psi covariate
    beta.lpsi ~ dnorm(0,0.01)
    
    #Slope for p covariates
    for(t in 1:4){
    beta.lp[t] ~ dnorm(0,0.01)
    }
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- int.psi +
    beta.lpsi * forest.site[i]
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,k] ~ dbern(mu.y[i,k])
    mu.y[i,k] <- z[i] * p[i,k]
    logit(p[i,k]) <- int.p +
    beta.lp[1] * equals(micro[i],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[2] * equals(micro[i],4) +   # cliff compared to no feature (intercept)
    beta.lp[3] * cam.height[i] +
    beta.lp[4] * walktest[i]
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
    
    #Slope for p covariates
    for(t in 1:4){
    beta.lp[t] ~ dnorm(0,0.01)
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
    logit(p[i,k]) <- int.p +
    beta.lp[1] * equals(micro[i],3) +   # wildlife trail compared to no feature (intercept)
    beta.lp[2] * equals(micro[i],4) +    # cliff compared to no feature (intercept)
    beta.lp[3] * cam.height[i] +
    beta.lp[4] * walktest[i]
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
# final model
outJ.mooseHM.all.fm <- jags(win.data, inits.1, params, "MSocc-HM-all-moose.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.mooseHM.all.fm,file="outJ-mooseHM-all-fm.RData")

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
# final model
outJ.mooseHM.TL.fm <- jags(win.data1, inits1.1, params, "MSocc-HM-two-null.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.mooseHM.TL.fm,file="outJ-mooseHM-TL-fm.RData")

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
# final model
outJ.mooseHM.TH.fm <- jags(win.data2, inits2.1, params, "MSocc-HM-TH-moose.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.mooseHM.TH.fm,file="outJ-mooseHM-TH-fm.RData")

## H+L cameras ##
y3 <- y[,c(2:3),]
win.data3 <- list(y=y3, n.site=dim(y3)[1],n.cts=dim(y3)[2],n.days=dim(y3)[3],
                 elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                 target.site=target.site,ruggedness=ruggedness,forest.site=forest.site,field.site=field.site,
                 target.cam=target.cam,micro=micro,cam.height=cam.height,walktest=walktest)
# Initial values
zst3 <- apply(y3,1,max,na.rm=T)      #inits for presence (z)
ast3 <- apply(y3,c(1,2),max,na.rm=T) #inits for availability (a)
inits3 <- function()list(z=zst3, a=ast3, psi0=0.5, beta.lpsi=c(0,0,0))
inits3.1 <- function()list(z=zst3, a=ast3, psi0=0.5)
# full model
outJ.mooseHM.HL <- jags(win.data3, inits3, params, "MSocc-HM-HL.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.mooseHM.HL,file="outJ-mooseHM-HL.RData")
# final model
outJ.mooseHM.HL.fm <- jags(win.data3, inits3.1, params, "MSocc-HM-two-null.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.mooseHM.HL.fm,file="outJ-mooseHM-HL-fm.RData")

## T cameras ##
y4 <- y[,1,]
win.data4 <- list(y=y4, n.site=dim(y4)[1],n.days=dim(y4)[2],
                 elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                 ruggedness=ruggedness[,1],forest.site=forest.site[,1],field.site=field.site[,1],
                 micro=micro[,1],cam.height=cam.height[,1],walktest=walktest[,1])
# Initial values
zst4 <- apply(y4,1,max,na.rm=T)      #inits for presence (z)
inits4 <- function()list(z=zst4, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits4.1 <- function()list(z=zst4, psi0=0.5, beta.lpsi=0)
# full model
outJ.mooseHM.T <- jags(win.data4, inits4, params, "occ-HM-T.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.mooseHM.T,file="outJ-mooseHM-T.RData")
# final model
outJ.mooseHM.T.fm <- jags(win.data4, inits4.1, params, "occ-HM-T-moose.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.mooseHM.T.fm,file="outJ-mooseHM-T-fm.RData")

## L cameras ##
y5 <- y[,3,]
win.data5 <- list(y=y5, n.site=dim(y5)[1],n.days=dim(y5)[2],
                 elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                 ruggedness=ruggedness[,3],forest.site=forest.site[,3],field.site=field.site[,3],
                 micro=micro[,3],cam.height=cam.height[,3],walktest=walktest[,3])
# Initial values
zst5 <- apply(y5,1,max,na.rm=T)      #inits for presence (z)
inits5 <- function()list(z=zst5, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits5.1 <- function()list(z=zst5, psi0=0.5, beta.lpsi=0)
# full model
outJ.mooseHM.L <- jags(win.data5, inits5, params, "occ-HM-L.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.mooseHM.L,file="outJ-mooseHM-L.RData")
# final model
outJ.mooseHM.L.fm <- jags(win.data5, inits5.1, params, "occ-HM-L-moose.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.mooseHM.L.fm,file="outJ-mooseHM-L-fm.RData")

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
# final model
outJ.roedeerHM.TH.fm <- jags(win.data2, inits2.1, params, "MSocc-HM-TH-roedeer.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.roedeerHM.TH.fm,file="outJ-roedeerHM-TH-fm.RData")

## H+L cameras ##
y3 <- y[,c(2:3),]
win.data3 <- list(y=y3, n.site=dim(y3)[1],n.cts=dim(y3)[2],n.days=dim(y3)[3],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  target.site=target.site,ruggedness=ruggedness,forest.site=forest.site,field.site=field.site,
                  target.cam=target.cam,micro=micro,cam.height=cam.height,walktest=walktest)
# Initial values
zst3 <- apply(y3,1,max,na.rm=T)      #inits for presence (z)
ast3 <- apply(y3,c(1,2),max,na.rm=T) #inits for availability (a)
inits3 <- function()list(z=zst3, a=ast3, psi0=0.5, beta.lpsi=c(0,0,0))
inits3.1 <- function()list(z=zst3, a=ast3, psi0=0.5)
# full model
outJ.roedeerHM.HL <- jags(win.data3, inits3, params, "MSocc-HM-HL.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.roedeerHM.HL,file="outJ-roedeerHM-HL.RData")
# final model
outJ.roedeerHM.HL.fm <- jags(win.data3, inits3.1, params, "MSocc-HM-HL-roedeer.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.roedeerHM.HL.fm,file="outJ-roedeerHM-HL-fm.RData")

## T cameras ##
y4 <- y[,1,]
win.data4 <- list(y=y4, n.site=dim(y4)[1],n.days=dim(y4)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,1],forest.site=forest.site[,1],field.site=field.site[,1],
                  micro=micro[,1],cam.height=cam.height[,1],walktest=walktest[,1])
# Initial values
zst4 <- apply(y4,1,max,na.rm=T)      #inits for presence (z)
inits4 <- function()list(z=zst4, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits4.1 <- function()list(z=zst4, psi0=0.5)
# full model
outJ.roedeerHM.T <- jags(win.data4, inits4, params, "occ-HM-T.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.roedeerHM.T,file="outJ-roedeerHM-T.RData")
# full model
outJ.roedeerHM.T.fm <- jags(win.data4, inits4.1, params, "occ-HM-T-null.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.roedeerHM.T.fm,file="outJ-roedeerHM-T-fm.RData")

## L cameras ##
y5 <- y[,3,]
win.data5 <- list(y=y5, n.site=dim(y5)[1],n.days=dim(y5)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,3],forest.site=forest.site[,3],field.site=field.site[,3],
                  micro=micro[,3],cam.height=cam.height[,3],walktest=walktest[,3])
# Initial values
zst5 <- apply(y5,1,max,na.rm=T)      #inits for presence (z)
inits5 <- function()list(z=zst5, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits5.1 <- function()list(z=zst5, psi0=0.5)
# full model
outJ.roedeerHM.L <- jags(win.data5, inits5, params, "occ-HM-L.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.roedeerHM.L,file="outJ-roedeerHM-L.RData")
# final model
outJ.roedeerHM.L.fm <- jags(win.data5, inits5.1, params, "occ-HM-L-null.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.roedeerHM.L.fm,file="outJ-roedeerHM-L-fm.RData")

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
# final model
outJ.hareHM.all.fm <- jags(win.data, inits.1, params, "MSocc-HM-all-null.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
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
# final model
outJ.hareHM.TL.fm <- jags(win.data1, inits1.1, params, "MSocc-HM-two-null.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
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
# final model
outJ.hareHM.TH.fm <- jags(win.data2, inits2.1, params, "MSocc-HM-two-null.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.hareHM.TH.fm,file="outJ-hareHM-TH-fm.RData")

## H+L cameras ##
y3 <- y[,c(2:3),]
win.data3 <- list(y=y3, n.site=dim(y3)[1],n.cts=dim(y3)[2],n.days=dim(y3)[3],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  target.site=target.site,ruggedness=ruggedness,forest.site=forest.site,field.site=field.site,
                  target.cam=target.cam,micro=micro,cam.height=cam.height,walktest=walktest)
# Initial values
zst3 <- apply(y3,1,max,na.rm=T)      #inits for presence (z)
ast3 <- apply(y3,c(1,2),max,na.rm=T) #inits for availability (a)
inits3 <- function()list(z=zst3, a=ast3, psi0=0.5, beta.lpsi=c(0,0,0))
inits3.1 <- function()list(z=zst3, a=ast3, psi0=0.5)
# full model
outJ.hareHM.HL <- jags(win.data3, inits3, params, "MSocc-HM-TH.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.hareHM.HL,file="outJ-hareHM-HL.RData")
# final model
outJ.hareHM.HL.fm <- jags(win.data3, inits3.1, params, "MSocc-HM-two-null.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.hareHM.HL.fm,file="outJ-hareHM-HL-fm.RData")

## T cameras ##
y4 <- y[,1,]
win.data4 <- list(y=y4, n.site=dim(y4)[1],n.days=dim(y4)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,1],forest.site=forest.site[,1],field.site=field.site[,1],
                  micro=micro[,1],cam.height=cam.height[,1],walktest=walktest[,1])
# Initial values
zst4 <- apply(y4,1,max,na.rm=T)      #inits for presence (z)
inits4 <- function()list(z=zst4, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits4.1 <- function()list(z=zst4, psi0=0.5)
# full model
outJ.hareHM.T <- jags(win.data4, inits4, params, "occ-HM-T.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.hareHM.T,file="outJ-hareHM-T.RData")
# final model
outJ.hareHM.T.fm <- jags(win.data4, inits4.1, params, "occ-HM-null.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.hareHM.T.fm,file="outJ-hareHM-T-fm.RData")

## L cameras ##
y5 <- y[,3,]
win.data5 <- list(y=y5, n.site=dim(y5)[1],n.days=dim(y5)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,3],forest.site=forest.site[,3],field.site=field.site[,3],
                  micro=micro[,3],cam.height=cam.height[,3],walktest=walktest[,3])
# Initial values
zst5 <- apply(y5,1,max,na.rm=T)      #inits for presence (z)
inits5 <- function()list(z=zst5, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits5.1 <- function()list(z=zst5, psi0=0.5)
# full model
outJ.hareHM.L <- jags(win.data5, inits5, params, "occ-HM-L.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.hareHM.L,file="outJ-hareHM-L.RData")
# full model
outJ.hareHM.L.fm <- jags(win.data5, inits5.1, params, "occ-HM-L-null.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.hareHM.L.fm,file="outJ-hareHM-L-fm.RData")

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
# final model
outJ.foxHM.TH.fm <- jags(win.data2, inits2.1, params, "MSocc-HM-two-null.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.foxHM.TH.fm,file="outJ-foxHM-TH-fm.RData")

## H+L cameras ##
y3 <- y[,c(2:3),]
win.data3 <- list(y=y3, n.site=dim(y3)[1],n.cts=dim(y3)[2],n.days=dim(y3)[3],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  target.site=target.site,ruggedness=ruggedness,forest.site=forest.site,field.site=field.site,
                  target.cam=target.cam,micro=micro,cam.height=cam.height,walktest=walktest)
# Initial values
zst3 <- apply(y3,1,max,na.rm=T)      #inits for presence (z)
ast3 <- apply(y3,c(1,2),max,na.rm=T) #inits for availability (a)
inits3 <- function()list(z=zst3, a=ast3, psi0=0.5, beta.lpsi=c(0,0,0))
inits3.1 <- function()list(z=zst3, a=ast3, psi0=0.5)
# full model
outJ.foxHM.HL <- jags(win.data3, inits3, params, "MSocc-HM-HL.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.foxHM.HL,file="outJ-foxHM-HL.RData")
# final model
outJ.foxHM.HL.fm <- jags(win.data3, inits3.1, params, "MSocc-HM-two-null.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.foxHM.HL.fm,file="outJ-foxHM-HL-fm.RData")

## T cameras ##
y4 <- y[,1,]
win.data4 <- list(y=y4, n.site=dim(y4)[1],n.days=dim(y4)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,1],forest.site=forest.site[,1],field.site=field.site[,1],
                  micro=micro[,1],cam.height=cam.height[,1],walktest=walktest[,1])
# Initial values
zst4 <- apply(y4,1,max,na.rm=T)      #inits for presence (z)
inits4 <- function()list(z=zst4, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits4.1 <- function()list(z=zst4, psi0=0.5)
# full model
outJ.foxHM.T <- jags(win.data4, inits4, params, "occ-HM-T.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.foxHM.T,file="outJ-foxHM-T.RData")
# final model
outJ.foxHM.T.fm <- jags(win.data4, inits4.1, params, "occ-HM-T-null.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.foxHM.T.fm,file="outJ-foxHM-T-fm.RData")

## L cameras ##
y5 <- y[,3,]
win.data5 <- list(y=y5, n.site=dim(y5)[1],n.days=dim(y5)[2],
                  elev=elev,forest.grid=forest.grid,field.grid=field.grid,
                  ruggedness=ruggedness[,3],forest.site=forest.site[,3],field.site=field.site[,3],
                  micro=micro[,3],cam.height=cam.height[,3],walktest=walktest[,3])
# Initial values
zst5 <- apply(y5,1,max,na.rm=T)      #inits for presence (z)
inits5 <- function()list(z=zst5, psi0=0.5, beta.lpsi=c(0,0,0,0,0,0))
inits5.1 <- function()list(z=zst5, psi0=0.5, beta.lpsi=c(0,0))
# full model
outJ.foxHM.L <- jags(win.data5, inits5, params, "occ-HM-L.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.foxHM.L,file="outJ-foxHM-L.RData")
# full model
outJ.foxHM.L.fm <- jags(win.data5, inits5.1, params, "occ-HM-L-fox.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.foxHM.L.fm,file="outJ-foxHM-L-fm.RData")

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
# final model
outJ.wolfHM.T.fm <- jags(win.data, inits.1, params, "occ-HM-T-null.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
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
inits.1 <- function()list(z=zst, psi0=0.5, beta.lpsi=0)
# full model
outJ.lynxHM.T <- jags(win.data, inits, params, "occ-HM-T.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.lynxHM.T,file="outJ-lynxHM-T.RData")
# final model
outJ.lynxHM.T.fm <- jags(win.data, inits.1, params, "occ-HM-T-lynx.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.lynxHM.T.fm,file="outJ-lynxHM-T-fm.RData")

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
# full model
outJ.badgerHM.T <- jags(win.data, inits, params, "occ-HM-T.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.badgerHM.T,file="outJ-badgerHM-T.RData")
# final model
outJ.badgerHM.T.fm <- jags(win.data, inits.1, params, "occ-HM-T-null.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
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
inits.1 <- function()list(z=zst, psi0=0.5, beta.lpsi=c(0,0,0))
# full model
outJ.squirrelHM.T <- jags(win.data, inits, params, "occ-HM-T.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
save(outJ.squirrelHM.T,file="outJ-squirrelHM-T.RData")
# final model
outJ.squirrelHM.T.fm <- jags(win.data, inits.1, params, "occ-HM-T-squirrel.txt", n.chains=nc1, n.thin=nt1, n.iter=ni1, n.burnin=nb1, parallel=T)
save(outJ.squirrelHM.T.fm,file="outJ-squirrelHM-T-fm.RData")

detach(data); rm("data")
