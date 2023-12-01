################################################################################
#### Nimble model code #########################################################

## Sollmann, Adenot, Spakovszky, Windt, Mattsson (submitted). Accounting for 
##   observation biases associated with counts of young when estimating fecundity: 
##   case study on the arboreal-nesting red kite (Milvus milvus). PCJ.

### Model code used in estimation simulation described in above referenced manu-
##    script
## For definition of nodes, see simulation script (Parameter estimation simulations.R)

##Model for Q1: paired data only - how well do we estimate p.obs?

modelQ1<-nimbleCode({
  
  l.p.nest[1]<-0
  for (j in 2:ncat.true){
    l.p.nest[j]~dnorm(0, sd=2)
  }
  p.nest[1:ncat.true]<-exp(l.p.nest[1:ncat.true])/sum(exp(l.p.nest[1:ncat.true]))
  ntrue[1:ncat.true]~dmulti(p.nest[1:ncat.true], sum.n)
  
  for (j in 1:ncat.true){
    
    #first as reference category
    l.p.obs[j,lims[j,1]]<-0
    for (k in (lims[j,1]+1):lims[j,2]){
    l.p.obs[j,k]~dnorm(0,sd=2)
    }
    p.obs[j,1:ncat.obs]<-get.pobs(lims=lims[j,1:2],
                                  l.p.obs=l.p.obs[j,1:ncat.obs],
                                  ncat=ncat.obs)

    obs[j,lims[j,1]:lims[j,2]]~dmulti(p.obs[j,lims[j,1]:lims[j,2]], 
                                      ntrue[j])
  }
  lam.hat<-sum((1:ncat.true)*p.nest[1:ncat.true])
})



##Model for Q2: combining paired with ground counts when both have same pi
modelQ2<-nimbleCode({
  l.p.nest[1]<-0
  for (j in 2:ncat.true){
    l.p.nest[j]~dnorm(0, sd=5)
  }
  p.nest[1:ncat.true]<-exp(l.p.nest[1:ncat.true])/sum(exp(l.p.nest[1:ncat.true]))
  ntrue[1:ncat.true]~dmulti(p.nest[1:ncat.true], sum.n)
  
  for (j in 1:ncat.true){
    #first as reference category
    l.p.obs[j,lims[j,1]]<-0
    for (k in (lims[j,1]+1):lims[j,2]){
      l.p.obs[j,k]~dnorm(0,sd=2)
    }
    p.obs[j,1:ncat.obs]<-get.pobs(lims=lims[j,1:2],
                                  l.p.obs=l.p.obs[j,1:ncat.obs],
                                  ncat=ncat.obs)
    
    obs[j,lims[j,1]:lims[j,2]]~dmulti(p.obs[j,lims[j,1]:lims[j,2]], ntrue[j])
  }
  
  for (i in 1:n){
    y[i]~dcat(p.obs[truth[i],1:ncat.obs])
    truth[i]~dcat(p.nest[1:ncat.true])
  }
  
  #lam - to be compared to output from model3
  lam.hat<-sum((1:ncat.true)*p.nest[1:ncat.true])
})


##Model for Q3: paired data from lam with ground counts only from lam.new
modelQ3<-nimbleCode({
  
  l.p.nest[1]<-0
  for (j in 2:ncat.true){
    l.p.nest[j]~dnorm(0, sd=2)
  }
  p.nest[1:ncat.true]<-exp(l.p.nest[1:ncat.true])/sum(exp(l.p.nest[1:ncat.true]))
  
  #model for double sampled nests
  for (j in 1:ncat.true){
    #first as reference category
    l.p.obs[j,lims[j,1]]<-0
    for (k in (lims[j,1]+1):lims[j,2]){
      l.p.obs[j,k]~dnorm(0,sd=2)
    }
    p.obs[j,1:ncat.obs]<-get.pobs(lims=lims[j,1:2],
                                  l.p.obs=l.p.obs[j,1:ncat.obs],
                                  ncat=ncat.obs)
    obs[j,lims[j,1]:lims[j,2]]~dmulti(p.obs[j,lims[j,1]:lims[j,2]], ntrue[j])
  }
  
  #model for ground counts only
  for (i in 1:n){
    y[i]~dcat(p.obs[truth[i],1:ncat.obs])
    truth[i]~dcat(p.nest[1:ncat.true])
  }
  lam.hat<-sum((1:ncat.true)*p.nest[1:ncat.true])
})


##Check effect of prior on p.nest: only estimating pi with log-linear implementation
##so that vague prior can be used
model0f<-nimbleCode({
  
  l.p.nest[1]<-0
  for (j in 2:ncat.true){
    l.p.nest[j]~dnorm(0, sd=10)
  }
  p.nest[1:ncat.true]<-exp(l.p.nest[1:ncat.true])/sum(exp(l.p.nest[1:ncat.true]))
  
  ntrue[1:ncat.true]~dmulti(p.nest[1:ncat.true], sum.n)
  lam.hat<-sum((1:ncat.true)*p.nest[1:ncat.true])
})


#################################################################################################################

##function to get observation probabilities when some cell probabilities in the p matrix
##are fixed to 0
get.pobs<-nimbleFunction(
  run=function(lims=double(1),
               l.p.obs=double(1),
               ncat=double(0)){
    returnType(double(1))
    
    realp<-exp(l.p.obs[lims[1]:lims[2]])/
           sum(exp(l.p.obs[lims[1]:lims[2]]))
    pvec<-rep(0, ncat)
    pvec[lims[1]:lims[2]]<-realp
    return(pvec)
  }
  )
