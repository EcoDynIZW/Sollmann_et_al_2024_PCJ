################################################################################
################################################################################

## Sollmann, Adenot, Spakovszky, Windt, Mattsson (submitted). Accounting for 
##   observation biases associated with counts of young when estimating fecundity: 
##   case study on the arboreal-nesting red kite (Milvus milvus). PCJ.

### Script to run parameter estimation simulation (Q1-Q3) described in above 
### referenced manuscript

library(nimble)
library(MCMCvis)

##write function to calculate posterior mode
get_mode<-function(x){
  dd<-density(x)
  md<-dd$x[dd$y == max(dd$y)]
  if (length(md)>1) md<-sample(md,1)
  return(md)}

##Source Nimble model code
source('R/Nimble model code.R')

##load empirical data of nest counts
##columns: true states; rows: observed states

datmat<-readRDS('data-raw/datmat.rds')
##observed counts
dm.vec<-rep(as.numeric(rownames(datmat)), apply(datmat, 1, sum))
##true counts - vectors don't match in order!!
dm.vec.true<-rep(as.numeric(colnames(datmat)), apply(datmat, 2, sum))
# 
# ##ML estimate of cell probabilities, pi|z
# p.mult.mean<-t(apply(datmat, 1, function(x)x/sum(x)))

##empirical true repro rate
lamr<-mean(dm.vec.true)

##distribution of #nestlings (pi) in sample
p.nest.samp<-table(dm.vec.true)/length(dm.vec.true)

##empirical classification matrix p(y|z)
pobsr<-t(apply(datmat, 2, function(x)x/sum(x)))


## First scenario dimension: reproductive output

##Create p.nest (pi) including z=4 and representing different scenarios of avg 
##  repro success

##general rule for pi: z=2 is most common, z=4 is rare
##baseline lam: 2.4, compare to higher and lower: 2 and 2.8

## baseline
p.nest<-c(0.1, 0.45, 0.4, 0.05)
lam<-sum((1:4)*p.nest)

## low
p.nest.low<-c(0.3, 0.45, 0.2,0.05 )
lam.low<-sum((1:4)*p.nest.low)

## high
p.nest.high<-c(0.05, 0.2, 0.65,0.1 )
lam.high<-sum((1:4)*p.nest.high)

## combine in single matrix
pnest.mat<-rbind(p.nest, p.nest.low, p.nest.high)

## number of true state categories, Z
ncat.true<-length(p.nest)


#### create p.obs (classification matrix, p)

## rules: at most one too many or too few
## p(correct count) goes down with true brood size
## more likely to under than overcount
## p(correct) >p(incorrect) for all states
## All rules based on empirical data

##rows: true state - column: observed state

p.obs<-matrix(c(0.95, 0.05,0,0,
                0.2, 0.75, 0.05, 0,
                0, 0.3, 0.6, 0.1,
                0,0,0.4, 0.6),
              nrow=4, ncol=4, byrow=T)

## number of observed state categories, K
ncat.obs<-ncol(p.obs)

## calculate observed lam (ie, lam if detection bias is ignored)
## note: lam = lambda, mean number of nestlings

lam.obs<-sum(apply(p.obs,1,function(x)sum(x*(1:4)))*p.nest)
(lam.obs-lam)/lam #7% negative bias

lam.obs.low<-sum(apply(p.obs,1,function(x)sum(x*(1:4)))*p.nest.low)
(lam.obs.low-lam.low)/lam.low #6% negative bias

lam.obs.high<-sum(apply(p.obs,1,function(x)sum(x*(1:4)))*p.nest.high)
(lam.obs.high-lam.high)/lam.high #7% negative bias


##Second scenario dimension: sample size
nnests<-c(10, 25, 50, 100, 250)

### combine n and pnest into scenarios
scen<-expand.grid(1:length(nnests), 1:nrow(pnest.mat))
nscen<-nrow(scen)


##now, create paired observation data sets for all scenarios
obs<-truth<-list()
niter<-100
  
for (sc in 1:nscen){
  
obs[[sc]]<-array(NA, c(ncat.true, ncat.obs ,niter))
truth[[sc]]<-array(NA, c(ncat.true, niter))

  for (iter in 1:niter){
    #generate #nests in each #nestling category, z
    truth[[sc]][,iter]<-as.vector(rmultinom(1, nnests[scen[sc,1]], 
                                            pnest.mat[scen[sc,2],]))
    #from these nests, generate observations based on p.obs (p(y|z)), y
    for (i in 1:ncat.true){
      obs[[sc]][i,,iter]<-rmultinom(1, truth[[sc]][i,iter], p.obs[i,])
    }
  }
}

## set parameters to monitor in Nimble for all questions/scenarios
## only monitor p.obs that are estimated
params<-c('lam.hat',paste('p.nest[',1:4, ']', sep=''), 
          paste('p.obs[', c('1, 1', 
                            '2, 1',
                            '1, 2',
                            '2, 2',
                            '3, 2',
                            '2, 3', 
                            '3, 3',
                            '4, 3',
                            '3, 4',
                            '4, 4'), ']',sep=''),
          paste('l.p.obs[',c(
            '1, 2',
            '2, 2',
            '2, 3', 
            '3, 3',
            '3, 4',
            '4, 4'), ']', sep=''),
          paste('l.p.nest[',2:4, ']', sep=''))

## provide cell ranges of p.obs (p) that are non-0
## if you set to 1,4 for all, you don't fix any cells to 0
lims<-matrix(c(1,2,
               1,3,
               2,4,
               3,4), 4, 2, byrow=T)

###############################################################################
##Q1: How much paired (observed + true) data do we need to estimate p.obs well

## create initial values and set up data
## Do so for a single simulated data set to run model once
## then, model can be re-run without compiling again

## number of nests in each true state z
ntrue<-apply(obs[[1]][,,2], 1, sum)

## data: observations for each true state, ntrue, total sample size
nimDat<-list(obs=obs[[1]][,,2], ntrue=ntrue, sum.n=sum(ntrue))

## constants: number of true and observation categories; indicator which cells
##            in p are not fixed to 0
nimConst<-list(ncat.obs=ncat.obs,
               ncat.true=ncat.true,
               lims=lims)

## inits log(p)
l.p.obs.in<-matrix(0, 4, 4)

inits<-function(){list(l.p.nest=rep(0,4),
                       l.p.obs=l.p.obs.in)}

## compile model and algorithm in Nimble once
model <- nimbleModel(modelQ1, constants = nimConst, data=nimDat,check = FALSE,
                     inits=inits())
cmodel <- compileNimble(model)       
conf.mcmc<-configureMCMC(model, monitors = params, thin=1)
mcmc <- buildMCMC(conf.mcmc)
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

## set up table to hold parameter estimates
p.hat<-array(NA, c(niter, length(params), 8, nscen))

## loop over scenarios
for(sc in 1:15){
  ## loop over iterations
  for (iter in 1:niter){
    ##set up data and replace in model object
    obs.new<-obs[[sc]][,,iter]
    ntrue.new<-apply(obs.new, 1, sum)
    cmodel$setData(obs=obs.new, ntrue=ntrue.new, sum.n=sum(ntrue.new))
    
    ## run model
    samp <- runMCMC(cmcmc, niter = 5000, nburnin = 2500, nchains=3,
                    progressBar = FALSE, inits = inits) 
    
    ## save parameters
    inn<-pmatch(params, colnames(samp[[1]]))
    samp2<-lapply(samp, function(x)x[,inn])
    summ<-MCMCsummary(samp2, func=get_mode)
    
    p.hat[iter,,,sc]<-as.matrix(summ)
    
    #check traceplot for ever 10th iteration
    #requires creating output folder
    if (iter %in% seq(10, niter, 10)){
      MCMCtrace(samp2, filename=paste('Parameter estimation/Sim results/Traceplot.',
                                     sc,'.', iter,'.pdf', sep=''))
    }
  }
  print(sc)
}

##save data, model output - adjust output folder as needed
saveRDS(list(obs=obs, p.hat=p.hat), 
        'Parameter estimation/Sim results loglinear pfixed/Data.Q1.rds')




################################################################################
##Q2: How much to we gain applying p.obs (p) from a small(ish) sample to a larger
##    sample of imperfect counts only when p.nest (pi) is the same

## Only consider scenario 2 (n=25) for paired data
## Maintain paired data from Q1 and simulate new imperfect counts only data

## For new data generation, don't use p.nest, but use average
## proportions from those data sets from Q1 where all true state categories
## were represented (see manuscript for details and rationale)

## read in data, results Q1 (adjust path as necessary)
ttt<-readRDS('Parameter estimation/Sim results loglinear pfixed/Data.Q1.rds')
obs<-ttt$obs #data
p.hat<-ttt$p.hat #estimates

##filter out 'incomplete z' data sets and non-converged iterations

keep2<-matrix(NA, niter, nscen)
for (j in 1:nscen){
  sub<-obs[[j]]
  nt<-apply(sub, c(1,3), sum)
  keep2[,j]<-apply(nt,2,function(x) all(x!=0))
}

##check convergence 
keep<-matrix(NA, niter, nscen)
for (j in 1:nscen){
  sub<-p.hat[,,6,j]
  maxR<-apply(sub, 1, max)
  keep[,j]<-maxR<=1.2
}

##set 'bad data' to FALSE
keep[!keep2]<-FALSE

##calculate pnest (pi) from sample of kept datasets
pi.hat<-array(NA, c(nscen, ncat.true, 2))
dimnames(pi.hat)[[3]]<-c('All data', 'Kept data')

for (i in 1: nscen){
  xx<-apply(obs[[i]], c(1,3), sum)
  n<-nnests[scen[i,1]]
  pi.hat[i,,1]<-apply(xx,1,sum)/(n*100)
  #remove iterations not kept for model results
  pix<-xx[,keep[,i]]
  pi.hat[i,,2]<-apply(pix,1,sum)/(n*sum(keep[,i]))
}

##create new pnest.mat using pi.hat for scenarios with n=25

pnest.mat.sample<-pnest.mat
pnest.mat.sample<- pi.hat[c(2,7,12),,2]

##set sample size scenarios for new (imperfect count only) data
n.single<-nnests

## create new data
new.truth<-new.obs<-list()

for (sc in 1:nscen){
  
new.truth[[sc]]<-array(NA, c(ncat.true,niter))
new.obs[[sc]]<-array(NA, c(ncat.true, ncat.obs,niter))

  for (iter in 1:niter){
    
    new.truth[[sc]][,iter]<-as.vector(rmultinom(1, n.single[scen[sc,1]], 
                                                pnest.mat.sample[scen[sc,2],]))
    for (i in 1:ncat.true){
      new.obs[[sc]][i,,iter]<-rmultinom(1, new.truth[[sc]][i,iter], p.obs[i,])
    }
  }#iteration loop
}#scenario loop


##create matrix to hold model output
outmat<-array(NA, c(niter,length(params),8,nscen))

## loop over scenarios
for (sc in 1:15){

  ## to select correct paired data, determine paired data scenarios from Q1 with n=25
  if(sc %in% (1:5)){
    get.scen.double<-2}
  if(sc %in% (6:10)){
    get.scen.double<-7}
  if(sc %in% (11:15)){
    get.scen.double<-12}
  
## set up data, inits and run once to compile; 
## here, do this once per scenario due to changing dimensions
  
y<-rep(1:4, apply(new.obs[[sc]][,,1], 2, sum))
n<-length(y)

nimDat<-list(y=y, obs=obs[[get.scen.double]][,,1], 
             ntrue=apply(obs[[get.scen.double]][,,1], 1, sum),
             sum.n=sum(obs[[get.scen.double]][,,1]))
nimConst<-list(n=n, ncat.obs=ncat.obs,
               ncat.true=ncat.true,
               lims=lims)

## create initial values for missing true states for new data
## make sure it conforms to limitations of data
u<-rbinom(n,1,0.5)
truth.in<-u*(y-1)+ (1-u)*(y+1)
truth.in[truth.in<1]<-1
truth.in[truth.in>4]<-4

inits<-function(){list(truth=truth.in,
                       l.p.nest=rep(0,4),
                       l.p.obs=l.p.obs.in)}

model <- nimbleModel(modelQ2, constants = nimConst, data=nimDat, 
                     inits=inits(), check = FALSE)
cmodel <- compileNimble(model)       
conf.mcmc<-configureMCMC(model, monitors = params, thin=1)
mcmc <- buildMCMC(conf.mcmc)
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

## loop over simulation iterations
for (iter in 1:niter){
  
  ## to speed things up, skip over paired data sets with missing true states
  if( any(apply(obs[[get.scen.double]][,,iter],1,sum) == 0) )next
  
  ## set up data and replace in model object
  y.new<-rep(1:4, apply(new.obs[[sc]][,,iter], 2, sum))
  obs.new<-obs[[get.scen.double]][,,iter]
  ntrue.new<-apply(obs.new, 1, sum)
  cmodel$setData(y=y.new, obs=obs.new, ntrue=ntrue.new, sum.n=sum(ntrue.new))

  ## create new inits (because here, they depend on data)
  u<-rbinom(n,1,0.5)
  truth.in<-u*(y.new-1)+ (1-u)*(y.new+1)
  truth.in[truth.in<1]<-1
  truth.in[truth.in>4]<-4
  
  inits<-function(){list(truth=truth.in,
                         l.p.nest=rep(0,4),
                         l.p.obs=l.p.obs.in)}
  samp <- runMCMC(cmcmc, niter = 5000, nburnin = 2500, nchains=3, 
                  inits = inits, progressBar = FALSE) 
  inn<-pmatch(params, colnames(samp[[1]]))
  samp2<-lapply(samp, function(x)x[,inn])
  summ<-MCMCsummary(samp2, func=get_mode)
  outmat[iter,,,sc]<-as.matrix(summ)
  
  #check traceplot for ever 10th iteration
  #will overwrite traceplots from Q1!
  if (iter %in% seq(10, niter, 10)){
    MCMCtrace(samp2, filename=paste('Parameter estimation/Sim results loglinear pfixed/Traceplot.',
                                   sc,'.', iter,'.pdf', sep=''))
  }
}
print(sc)
}

##save new data and model output
saveRDS(list(new.obs=new.obs, outmat=outmat), 
        'Parameter estimation/Sim results loglinear pfixed/Data.Q2.rds')



################################################################################
##Q3: When can we use estimates of p.obs (p) to estimate lam in imperfect count
##    only data set when p.nest (pi) is NOT the same

## for this, we use paired obs with lambda=2.4 (Q1) with imperfect obs only
## with lambda=2 or lambda=2.8 (newly generated)


## extract imperfect only data from those generated under Q1 (read in if not already 
## here from Q1 or Q2)
new.obs.low<-obs[[8]] #observations for n=50 for low lambda
new.obs.high<-obs[[13]] #observations for n=50 for high lambda

## these two data sets (in loop below) are combined with all sample sizes of paired 
## data for lam=2.4

nscen2<-5 #5 scenarios each for lam.low and lam.high

#create arrays to hold model output
outmat.high<-outmat.low<-array(NA, c(dim(outmat)[1:3], nscen2))

##set up inits and compile data to run model once

  y<-rep(1:4, apply(new.obs.low[,,1], 2, sum))
  n<-length(y)
  
  nimDat<-list(y=y, obs=obs[[1]][,,1], ntrue=apply(obs[[1]][,,1], 1, sum))
  nimConst<-list(n=n, 
                 ncat.obs=ncat.obs, 
                 ncat.true=ncat.true,
                 lims=lims)
  
  u<-rbinom(n,1,0.5)
  truth.in<-u*(y-1)+ (1-u)*(y+1)
  truth.in[truth.in<1]<-1
  truth.in[truth.in>4]<-4
  
  inits<-function(){list(truth=truth.in,
                         l.p.nest=rep(0,4),
                         l.p.obs=l.p.obs.in)}
  
  model <- nimbleModel(modelQ3, constants = nimConst, data=nimDat, 
                       inits=inits(), check = FALSE)
  cmodel <- compileNimble(model)       
  conf.mcmc<-configureMCMC(model, monitors = params, thin=1)
  mcmc <- buildMCMC(conf.mcmc)
  cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)
  
  ## loop over scenarios
  for (sc in 1:nscen2){
    ## loop over iterations
  for (iter in 1:niter){
    ##we combine the n=50 new imperfect count only data with different sized
    ##paired data sets
    y.new<-rep(1:4, apply(new.obs.low[,,iter], 2, sum))
    obs.new<-obs[[sc]][,,iter]

    ##skip iterations with missing truths in paired data
    if(any(apply(obs.new,1,sum)==0)) next

    ##create new inits
    u<-rbinom(n,1,0.5)
    truth.in<-u*(y.new-1)+ (1-u)*(y.new+1)
    truth.in[truth.in<1]<-1
    truth.in[truth.in>4]<-4

    inits<-function(){list(truth=truth.in,
                           l.p.nest=rep(0,4),
                           l.p.obs=l.p.obs.in)}

    cmodel$setData(y=y.new, obs=obs.new, ntrue=apply(obs.new, 1, sum))

    ##run these longer...
    samp <- runMCMC(cmcmc, niter = 10000, nburnin = 5000, nchains=3,
                    inits = inits, progressBar = FALSE)
    inn<-pmatch(params, colnames(samp[[1]]))
    samp2<-lapply(samp, function(x)x[,inn])
    summ<-MCMCsummary(samp2, func=get_mode)
    outmat.low[iter,,,sc]<-as.matrix(summ)
    
    ##traceplots for every 10th iteration
    if (iter %in% seq(10, niter, 10)){
    MCMCtrace(samp2, filename=paste('Parameter estimation/Sim results loglinear pfixed/TraceplotLOW.',
                                    sc,'.', iter,'.pdf', sep=''))
    }
    
    ##remove data, model output, to avoid possible errors
    rm(y.new); rm(samp); rm(summ); rm(inits); rm(samp2)

    ## same for high lam
    y.new<-rep(1:4, apply(new.obs.high[,,iter], 2, sum))
    obs.new<-obs[[sc]][,,iter] 
    u<-rbinom(n,1,0.5)
    truth.in<-u*(y.new-1)+ (1-u)*(y.new+1)
    truth.in[truth.in<1]<-1
    truth.in[truth.in>4]<-4
    
    inits<-function(){list(truth=truth.in,
                           l.p.nest=rep(0,4),
                           l.p.obs=l.p.obs.in)}
    
    cmodel$setData(y=y.new, obs=obs.new, ntrue=apply(obs.new, 1, sum))
    
    samp <- runMCMC(cmcmc, niter = 10000, nburnin = 5000, nchains=3, 
                    inits = inits, progressBar = FALSE) 
    inn<-pmatch(params, colnames(samp[[1]]))
    samp2<-lapply(samp, function(x)x[,inn])
    summ<-MCMCsummary(samp2, func=get_mode)
    outmat.high[iter,,,sc]<-as.matrix(summ)
    #traceplots every 10th iteration
    if (iter %in% seq(10, niter, 10)){
      MCMCtrace(samp2, filename=paste('Parameter estimation/Sim results loglinear pfixed/TraceplotHIGH.',
                                      sc,'.', iter,'.pdf', sep=''))
    }
    
  }
  print(sc)
  }
  
  #save model output (no new data generated)
  saveRDS(list(outmat.low=outmat.low, outmat.high=outmat.high), 
          'Parameter estimation/Sim results loglinear pfixed/Data.Q3.rds')
 