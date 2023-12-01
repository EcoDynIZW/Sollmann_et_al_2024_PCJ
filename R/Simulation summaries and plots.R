################################################################################
#### summarize estimation simulations ##########################################

## Sollmann, Wendt, Mattsson (submitted). Accounting for observation biases 
##   associated with counts of young when estimating fecundity: case study on 
##   the arboreal-nesting red kite (Milvus milvus). PCJ.

### Script to summarize output from parameter estimation simulation (Q1-Q3) 
### described in above referenced manuscript and produce figures

##DISCLAIMER: This is most likely not the cleanest or most efficient way to make
##            summarize results and make plots. Code provided only for transpa-
##            rency, not as template for similar analyses. I am on the 'what 
##            works is good enough' variety...


library(writexl)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(MASS) 
library(reshape2) 
library(reshape)
library(glue)


###############################################################################
### Input parameters, labels and other assorted pieces of information #########

##pi for lam=2.4
p.nest<-c(0.1, 0.45, 0.4, 0.05)
lam<-sum((1:4)*p.nest)
##pi for lam=2
p.nest.low<-c(0.3, 0.45, 0.2,0.05 )
lam.low<-sum((1:4)*p.nest.low)
##pi for lam=2.8
p.nest.high<-c(0.05, 0.2, 0.65,0.1 )
lam.high<-sum((1:4)*p.nest.high)

#number of possible categories for pi
ncat.true<-length(p.nest)

#number of nests surveyed
nnests<-c(10, 25, 50, 100, 250)
n.single<-nnests

#true observation matrix; rows = truth, columns = observations
p.obs<-matrix(c(0.95, 0.05,0,0,
                0.2, 0.75, 0.05, 0,
                0, 0.3, 0.6, 0.1,
                0,0,0.4, 0.6),
              nrow=4, ncol=4, byrow=T)

#number of observation categories
ncat.obs<-ncol(p.obs)

##some labels and parameter orders for some of the plots
labs<-c('y=1|z=1', 'y=2|z=1', 
        'y=1|z=2', 'y=2|z=2','y=3|z=2', 
        'y=2|z=3','y=3|z=3','y=4|z=3',
        'y=3|z=4', 'y=4|z=4')
ord<-c(6,8,7,9,11,10,12,14,13,15)
ord.lp<-c(16:21)


#make vector of true.p in correct order (by truth)
pb0<-as.vector(t(p.obs))
#turn into matrix for bias calculation
true.p<-matrix(c(pb0[pb0>0]), 100, 10, byrow=TRUE)

##calculate input l.p.obs - not needed for plots
lpo.in<-matrix(0, 4,4)
lpo.in[1,2]<-diff(log(p.obs[1,1:2]))
lpo.in[2,2:3]<- log(p.obs[2,2:3])-log(p.obs[2,1])
lpo.in[3,3:4]<-log(p.obs[3,3:4])-log(p.obs[3,2])
lpo.in[4,4]<-diff(log(p.obs[4,3:4]))

##vector only of estimated cells in lpo.in
lpo.in[p.obs==0]<-NA
#full log prob vector
lpvec<-na.omit(as.vector(t(lpo.in)))
#only those estimated
lpvec.hat<-lpvec[lpvec!=0]


###calculate lambda you'd observe if you didn't account for observation bias
lam.obs<-sum(apply(p.obs,1,function(x)sum(x*(1:4)))*p.nest)
(lam.obs-lam)/lam #7% negative bias

lam.obs.low<-sum(apply(p.obs,1,function(x)sum(x*(1:4)))*p.nest.low)
(lam.obs.low-lam.low)/lam.low #6% negative bias

lam.obs.high<-sum(apply(p.obs,1,function(x)sum(x*(1:4)))*p.nest.high)
(lam.obs.high-lam.high)/lam.high #7% negative bias


##make matrix with true lambda values, for bias calculation
truelam<-rep(c(lam, lam.low, lam.high), each=5)
tm<-matrix(truelam, 100, 15, byrow=T)


### combine nnest and pnest into scenarios
pnest.mat<-rbind(p.nest, p.nest.low, p.nest.high)
scen<-expand.grid(1:length(nnests), 1:nrow(pnest.mat))
nscen<-nrow(scen)

#calculate log(pnest) - not needed for plots
l.pn<-rbind(log(p.nest), log(p.nest.low), log(p.nest.high))
l.pn.mat<-t(apply(l.pn,1,function(x)x-x[1]))

#simulation settings and names of parameters in model output
niter<-100

#lam.hat: estimate of lambda based on estimates of pi
#p.nest: name for pi in model
#p.obs: name for p in model
#l.p.obs: log observation probabilities
#l.p.nest: log pi
#for both log(p) and log(pi), one category is fixed to 0
#as reference, so not all categories are monitored

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

#provide cell ranges of p.obs that are non-0
#if you set to 1,4 for all, you don't fix any cells to 0
lims<-matrix(c(1,2,
               1,3,
               2,4,
               3,4), 4, 2, byrow=T)


###############################################################################
#################### Summaries, plots Q1 ######################################

## NOTE: all 'write' commands are commented out to avoid substituting existing 
##       files; uncomment if you want to write output files

##read in results, data from Q1 simulations - adjust path as necessary!
ttt<-readRDS('Parameter estimation/Sim results loglinear pfixed/Data.Q1.rds')
obs<-ttt$obs
#obs is a list of length 15, corresponding to the 15 scenarios
#dependent on n and p.nest - encoded in scen

##each element of obs is a 3d array: truth by observation by iteration
##Eg, obs[[1]][,,1] is one data set for scenario 1 (n=10, lam=2.4); if you sum over rows, 
##you get true number of nests in each category, if you sum over columns,
##you get observed number of nests in each category

p.hat<-ttt$p.hat
#p.hat is array with model estimates
#iteration by parameter by summary statistic by scenario

##adding dimension names to make navigation easier
dimnames(p.hat)[[2]]<-params
dimnames(p.hat)[[3]]<-c('mean', 'sd', '2.5%', '50%', '97.5%', 'Rhat','n.eff',
                        'mode')
dimnames(p.hat)[[4]]<-paste('n=',nnests[scen[,1]], '; lam=', truelam, sep='')

##check for missing truths
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

##set 'bad data' (ie, w/ missing truths) to FALSE
keep[!keep2]<-FALSE

#set estimates, SDs of parameters in models with imperfect convergence
#to NA so they don't appear in plots and summary stats
for (i in 1:24){
  p.hat[,i,1,][!keep]<-NA
  p.hat[,i,2,][!keep]<-NA
}

##make matrix with iteration numbers

niter.Q1<-data.frame(Scenario=dimnames(p.hat)[[4]],
                     `Missing truths`=100-apply(keep2, 2, sum),
                     `Not converged` = apply(keep2, 2, sum)-apply(keep,2,sum),
                     `Total kept` = apply(keep,2,sum),
                     check.names = FALSE)
#write_xlsx(niter.Q1, 'Parameter estimation/Tables for ms/Niter.Q1.xlsx')

##calculate mean lam, p.nest in kept data
lhat<-matrix(NA, nscen,2)
pi.hat<-array(NA, c(nscen, ncat.true, 2))

colnames(lhat)<-dimnames(pi.hat)[[3]]<-c('All data', 'Kept data')

for (i in 1: nscen){
    xx<-apply(obs[[i]], c(1,3), sum)
    n<-nnests[scen[i,1]]
    lamx<-apply(xx, 2, function(x)sum(x*(1:4))/n)
    lhat[i,1]<-mean(lamx)
    #pix<-apply(xx, 2, function(x)x/n)
    pi.hat[i,,1]<-apply(xx,1,sum)/(n*100)
    #remove iterations not kept for model results
    lamx[!keep[,i]]<-NA
    lhat[i,2]<-mean(lamx, na.rm=TRUE)
    pix<-xx[,keep[,i]]
    pi.hat[i,,2]<-apply(pix,1,sum)/(n*sum(keep[,i]))
}

lam.sample<-data.frame(Scenario=dimnames(p.hat)[[4]],
                       `True lambda`=tm[1,],
                       apply(lhat, 2, function(x)round(x, dig=2)),
                       `Bias truth`=round(((lhat[,2]-tm[1,])/tm[1,])*100, dig=2),
                       `Bias all data`=round(((lhat[,2]-lhat[,1])/lhat[,1])*100, dig=2))
#write_xlsx(lam.sample, 'Parameter estimation/Tables for ms/Sample bias.xlsx')

## calculate summary stats: average bias, RMSE and avg CV (for lam only)
## only calculate for back-transformed probabilities
pars<-dimnames(p.hat)[[2]][1:15]
npar<-length(pars)

##make matrix of true input values, param by scenario, to match model output
tinp.o<-matrix(c(rep(c(lam, p.nest, pb0[pb0>0]), 5),
               rep(c(lam.low, p.nest.low, pb0[pb0>0]), 5),
               rep(c(lam.high, p.nest.high, pb0[pb0>0]), 5)), nrow=15, ncol=15)

##replace with observed lam for n=10 and n=25
##lambdas are in row 1
tinp<-tinp.o
tinp[1,1]<-lhat[1,2]
tinp[1,6]<-lhat[6,2]
tinp[1,11]<-lhat[11,2]
tinp[1,2]<-lhat[2,2]
tinp[1,7]<-lhat[7,2]
tinp[1,12]<-lhat[12,2]

Bias.Q1<-RMSE.Q1<-CV.Q1<-matrix(NA, npar, nscen)

#note: this works even though p.hat has more parameters because it's
#the first 15 we're interested in
for (i in 1:npar){
  for (j in 1:nscen){
    Bias.Q1[i,j]<-median(p.hat[,i,1,j]-tinp[i,j], na.rm=TRUE)#absolute bias
    RMSE.Q1[i,j]<-sum(((p.hat[,i,1,j]-tinp[i,j])^2)/niter.Q1$`Total kept`[j], 
                   na.rm=TRUE)
    #only for lam, relative bias, CV
    if (i == 1){
    Bias.Q1[i,j]<-median(((p.hat[,i,1,j]-tinp[i,j])/tinp[i,j])*100, na.rm=TRUE)
    CV.Q1[i,j]<-median((p.hat[,i,2,j]/p.hat[,i,1,j])*100, na.rm=TRUE)
    }
  }
}
Bias.Q1<-apply(Bias.Q1, 2, function(x)round(x, dig=3) )
CV.Q1<-apply(CV.Q1, 2, function(x)round(x, dig=2) )
RMSE.Q1<-apply(RMSE.Q1, 2, function(x)round(x, dig=3) )

bias1.df<-data.frame(Parameter=pars,Bias.Q1)
cv1.df<-data.frame(Parameter=pars,CV.Q1)
rmse1.df<-data.frame(Parameter=pars,RMSE.Q1)
colnames(bias1.df)<-colnames(cv1.df)<-colnames(rmse1.df)<-c('Parameter',
                                                          dimnames(p.hat)[[4]])
# write_xlsx(list(Bias.Q1 = bias1.df, 
#                 CV.Q1 = cv1.df,
#                 RMSE.Q1 = rmse1.df), 'Parameter estimation/Tables for ms/Summary stats Q1.xlsx')


################################################################################
############################### plots Q1 #######################################

##Figure 1 main ms: Estimates of detection probabilities for lam=2.4 (Q1)

##use scenarios 1-5

valuep1<-NULL
n.inp1<-NULL
parmp<-c('y=1|z=1', 'y=2|z=1', 
         'y=1|z=2', 'y=2|z=2','y=3|z=2', 
         'y=2|z=3','y=3|z=3','y=4|z=3',
         'y=3|z=4', 'y=4|z=4')
idx.inp1<-NULL


for (sc in 1:5){
  
  plt.mat2<-p.hat[,ord,1,sc]-true.p
  
  valuep1<-c(valuep1, as.vector(plt.mat2))
  n.inp1<-c(n.inp1, rep(nnests[sc], 1000))
  idx.inp1<-c(idx.inp1, rep(parmp, each=100))
  
}
bdfp.q1<-data.frame(value=valuep1, parameter=idx.inp1,
                    n=n.inp1)
bdfp.q1$n<-factor(bdfp.q1$n, levels=as.character(nnests))
bdfp.q1$parameter<-factor(bdfp.q1$parameter,
                          levels=parmp)

v_line <- data.frame(
  yintercept = 0
)

q1.p<-ggplot(bdfp.q1, aes(x=parameter, y=value)) +
  geom_boxplot(fatten=0.5)+
  labs(x = NULL, y = 'Error') +
  theme_bw() +
  geom_hline(data=v_line,aes(yintercept = yintercept),
             linetype='dashed', color = 'orangered', linewidth=0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~forcats::fct_inorder(glue("n = {n}")), 
             nrow=5, ncol=1, strip.position = 'right')

## Save in appropriate size to fit on page in-text
jpeg('Figure 1 Bias p Q1.jpg',
     width = 16, 
     height = 19, 
     units='cm',
     res=600)
q1.p
dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Figure S1: relative bias and CV of lambda for all scenarios

### tm has true value of lambda for each scenario
##substitute with observed lambda for n=10 and 25
tm.Q1<-tm
tm.Q1[,1]<-tinp[1,1]
tm.Q1[,2]<-tinp[1,2]
tm.Q1[,6]<-tinp[1,6]
tm.Q1[,7]<-tinp[1,7]
tm.Q1[,11]<-tinp[1,11]
tm.Q1[,12]<-tinp[1,12]

##first, bias in lam
df.lam<-data.frame(((p.hat[,1,1,]-tm.Q1)/tm.Q1)*100)
colnames(df.lam)<-c(paste( rep(1:3, each=5), nnests, sep='.' ))
meltQ1 <- melt(df.lam, id = NULL) 

coll2<-rep('2.4', nrow(meltQ1))
coll2[as.character(meltQ1$variable) %in% 
        c("2.10","2.25","2.50","2.100","2.250")]  <-'2.0'
coll2[as.character(meltQ1$variable) %in% 
        c("3.10","3.25","3.50","3.100","3.250")]  <-'2.8'
coll2<-factor(coll2,levels=c('2.4', '2.0', '2.8'))

pQ1b<-ggplot(meltQ1, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=coll2), fatten = 0.5)+
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', linewidth=0.5) +
  labs(x = '# paired counts', y = expression(paste("Relative error in  ", lambda, sep=''))) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  #scale_x_discrete(labels=rep(nnests, 3 ) ) +
  geom_text(x=0.7, y=24, label='A)', 
            hjust='left',vjust='bottom',  size=4, color='black')+
  scale_fill_manual(name=expression(lambda),
                    #labels=c("2.4", "2.0", "2.8"),
                    values = c("2.4"='white',
                               "2.0"='lightgrey',
                               "2.8"='darkgrey'))


##second, CV of lam

df.lam<-data.frame((p.hat[,1,2,]/p.hat[,1,1,])*100)
colnames(df.lam)<-c(paste( rep(1:3, each=5), nnests, sep='.' ))

meltQ1b <- melt(df.lam, id = NULL) 

coll2<-rep('2.4', nrow(meltQ1b))
coll2[as.character(meltQ1b$variable) %in% 
        c("2.10","2.25","2.50","2.100","2.250")]  <-'2.0'
coll2[as.character(meltQ1b$variable) %in% 
        c("3.10","3.25","3.50","3.100","3.250")]  <-'2.8'
coll2<-factor(coll2,levels=c('2.4', '2.0', '2.8'))

pQ1c<-ggplot(meltQ1b, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=coll2), fatten = 0.5)+
  labs(x = '# paired counts', y = expression(paste("CV of  ", lambda, sep=''))) +
  theme_bw() +
  scale_x_discrete(labels=rep(nnests, 3 ) ) +
  geom_text(x=0.7, y=14.2, label='B)', 
            hjust='left',vjust='bottom',  size=4, color='black')+
  scale_fill_manual(name=expression(paste(lambda, ': ', sep='')),
                    #labels=c("2.4", "2.0", "2.8"),
                    values = c("2.4"='white',
                               "2.0"='lightgrey',
                               "2.8"='darkgrey'))+
  
  theme(legend.direction = "horizontal", legend.box = "vertical",
        legend.position='bottom')# + 


##stack on top of each other
plotfull3<-plot_grid(pQ1b, pQ1c, align = 'v', ncol=1 , rel_heights = c(0.85,1))

## change size as appropriate
jpeg('Figure S1 Bias CV lambda Q1.jpg', width = 16, 
    height = 16, unit='cm', res=600)
plotfull3
dev.off()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Figure S2: bias of pnest

bvec<-NULL
scenario<-NULL
idx<-NULL

for (sc in 1:nscen){
  true.pnest<-matrix(pnest.mat[scen[sc,2],], 100, 4, byrow=TRUE)
  
  ##again, for n=10 and n=25, replace input value with mean from kept sample
  ##values in pi.hat[,,2]
  if (sc %in% c(1,2, 6,7, 11,12)){
    true.pnest<-matrix(pi.hat[sc,,2], 100, 4, byrow=TRUE)
  }
  
  bvec<-c(bvec, as.vector(p.hat[,2:5,1,sc]-true.pnest))
  scenario<-c(scenario, rep(sc, 400))
  idx<-rep(c('pi[1]', 'pi[2]', 'pi[3]', 'pi[4]'), each=100)
}

bdf<-data.frame(value=bvec, scenario=as.factor(scenario), 
                parameter=as.factor(idx))

coll3<-rep('2.4', length(scenario))
coll3[as.numeric(scenario) %in% (6:10)]<-'2.0'
coll3[as.numeric(scenario) %in% (11:15)]<-'2.8'
coll3<-factor(coll3, levels = c('2.4', '2.0', '2.8'))

#prnam<-paste('pi[', 1:4, ']', sep='')
#names(prnam)<-c('1', '2', '3', '4')

q1.pi<-ggplot(bdf, aes(x=scenario, y=value)) +
  geom_boxplot(aes(fill=coll3), fatten=0.5)+
  labs(x = '# paired counts', y = 'Error') +
  theme_bw() +
  scale_x_discrete(labels=rep(nnests, 3 ) ) +
  # geom_text(x=0.7, y=14.2, label='B)', 
  #           hjust='left',vjust='bottom',  size=4, color='black')+
  scale_fill_manual(name=expression(paste(lambda, ': ', sep='')),
                    #labels=c("2.4", "2.0", "2.8"),
                    values = c("2.4"='white',
                               "2.0"='lightgrey',
                               "2.8"='darkgrey'))+
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', linewidth=0.5) +
  theme(legend.position = 'right')+
  facet_grid(rows = vars(parameter),
             labeller = label_parsed)

jpeg('Figure S2 Bias pi Q1.jpg', width = 16, 
    height = 18, unit='cm', res=600)
q1.pi
dev.off()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Figure S3: bias of p for lam = 2 and lam = 2.8

##scenarios 6-10 and 11-15

valuep1b<-NULL
n.inp1b<-NULL
idx.inp1b<-NULL
lam.inp1b<-NULL

for (sc in 6:10){
  
  plt.mat2<-cbind(p.hat[,ord,1,sc]-true.p, #low
                  p.hat[,ord,1,sc+5]-true.p) #high
  
  valuep1b<-c(valuep1b, as.vector(plt.mat2))
  n.inp1b<-c(n.inp1b, rep(nnests[sc-5], 2000))
  idx.inp1b<-c(idx.inp1b, rep(rep(parmp, each=100), 2))
  lam.inp1b<-c(lam.inp1b, rep(c('2.0', '2.8'), each = 1000))
  
}

bdfp.q1b<-data.frame(value=valuep1b, parameter=idx.inp1b,
                    n=n.inp1b, lam=lam.inp1b)
bdfp.q1b$n<-factor(bdfp.q1b$n, levels=as.character(nnests))
bdfp.q1b$parameter<-factor(bdfp.q1b$parameter,
                          levels=parmp)

q1.pb<-ggplot(bdfp.q1b, aes(x=parameter, y=value)) +
  geom_boxplot(fatten=0.5)+
  labs(x = NULL, y = 'Error') +
  theme_bw() +
  geom_hline(data=v_line,aes(yintercept = yintercept),
             linetype='dashed', color = 'orangered', linewidth=0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #theme(legend.position = 'bottom')+
  # facet_grid(n~lam)
  # all this extra cajigger is to keep factor level order as it should be
  facet_grid(forcats::fct_inorder(glue("'n = {n}'")) ~ glue('lambda*" = {lam}"'),
             labeller = label_parsed)


jpeg('Figure S3 Bias p lam2&2.8 Q1.jpg', width = 16, 
    height = 19, unit='cm', res=600)
q1.pb
dev.off()


################################################################################
#####################Results, plots Q2 #########################################

#### summaries

##read in results, addl data for Q2
tt<-readRDS('Parameter estimation/Sim results loglinear pfixed/Data.Q2.rds')

new.obs<-tt$new.obs
##new.obs is list of new observations, of which only ground counts (ie, sum over
##columns) were used in the model, combined with the n=25 scenarios from Q1

outmat<-tt$outmat
##matrix holding model estimates from Q2 analysis
#iteration by parameter by summary statistic by scenario

##adding dimension names to make navigation easier
dimnames(outmat)[[2]]<-params
dimnames(outmat)[[3]]<-c('mean', 'sd', '2.5%', '50%', '97.5%', 'Rhat','n.eff',
                        'mode')
dimnames(outmat)[[4]]<-paste('n.ground=',nnests[scen[,1]], '; lam=', truelam, sep='')

##remove non-converged iterations 
##iterations for which paired data were 'bad' (ie, did not contain all
##true categries) already are set to NA

#get number of iterations kept after filtering out 'bad' datasets
nkeep<-NULL
for (j in 1:nscen){
  sub<-outmat[,1,6,j]
  nkeep[j]<-sum(!is.na(sub))
}

keep<-matrix(NA, niter, nscen)
for (j in 1:nscen){
  sub<-outmat[,,6,j]
  maxR<-apply(sub, 1, max)
  keep[,j]<-maxR<=1.2
}
apply(keep,2,sum, na.rm=TRUE)
keep[is.na(keep)]<-FALSE

##make new results matrix in which estimates of non-converged runs are set to NA
##so they don't show up in the plots
outmat2<-outmat
for(i in 1:24){
  outmat2[,i,1,][!keep]<-NA
  outmat2[,i,2,][!keep]<-NA
}

##make matrix with iteration numbers
niter.Q2<-data.frame(Scenario.ground=dimnames(p.hat)[[4]],
                     Scenario.paired = paste('n=25; lam=',rep(c('2.4', '2', '2.8'), each=5), 
                                             sep=''),
                     `Missing truths`=100-nkeep,
                     `Not converged` = nkeep-apply(keep,2,sum, na.rm=TRUE),
                     `Total kept` = apply(keep,2,sum, na.rm=TRUE),
                     check.names = FALSE)
#write_xlsx(niter.Q2, 'Parameter estimation/Tables for ms/Niter.Q2.xlsx')

## calculate summary stats: average bias, RMSE and avg CV (for lam only)
## only calculate for back-transformed probabilities

##make matrix of true input values, param by scenario, to match model output
##use observed lam, pnest from kept paired data
#lhat, pi.hat
tinp.o2<-matrix(c(rep(c(lhat[2,2], pi.hat[2,,2], pb0[pb0>0]), 5),
                 rep(c(lhat[7,2], pi.hat[7,,2], pb0[pb0>0]), 5),
                 rep(c(lhat[12,2], pi.hat[12,,2], pb0[pb0>0]), 5)), nrow=15, ncol=15)


Bias.Q2<-RMSE.Q2<-CV.Q2<-matrix(NA, npar, nscen)

#note: this works even though p.hat has more parameters because it's
#the first 15 we're interested in
for (i in 1:npar){
  for (j in 1:nscen){
    Bias.Q2[i,j]<-median(outmat2[,i,1,j]-tinp.o2[i,j], na.rm=TRUE)#absolute bias
    RMSE.Q2[i,j]<-sum(((outmat2[,i,1,j]-tinp.o2[i,j])^2)/niter.Q2$`Total kept`[j], 
                      na.rm=TRUE)
    #only for lam, relative bias
    if (i == 1){
      Bias.Q2[i,j]<-median(((outmat2[,i,1,j]-tinp.o2[i,j])/tinp.o2[i,j])*100, na.rm=TRUE)
      CV.Q2[i,j]<-median((outmat2[,i,2,j]/outmat2[,i,1,j])*100, na.rm=TRUE)
    }
  }
}
Bias.Q2<-apply(Bias.Q2, 2, function(x)round(x, dig=3) )
CV.Q2<-apply(CV.Q2, 2, function(x)round(x, dig=2) )
RMSE.Q2<-apply(RMSE.Q2, 2, function(x)round(x, dig=3) )

bias2.df<-data.frame(Parameter=pars,Bias.Q2)
cv2.df<-data.frame(Parameter=pars,CV.Q2)
rmse2.df<-data.frame(Parameter=pars,RMSE.Q2)

scennames<-gsub('n', 'n(ground)',niter.Q2$Scenario.ground, fixed=TRUE)
colnames(bias2.df)<-colnames(cv2.df)<-colnames(rmse2.df)<-c('Parameter',
                                                            scennames)
# write_xlsx(list(Bias.Q2 = bias2.df, 
#                 CV.Q2 = cv2.df,
#                 RMSE.Q2 = rmse2.df), 'Parameter estimation/Tables for ms/Summary stats Q2.xlsx')


##FOR PLOTS
##combine estimates for respective Q1, n=25 with estimates of Q2
plt.mat<-cbind( p.hat[,1,1,2], outmat2[,1,1,1:5],
                p.hat[,1,1,7], outmat2[,1,1,6:10],
                p.hat[,1,1,12], outmat2[,1,1,11:15])

##make matrix of true lambda to calculate bias
##adjusted to realized lambda in kept sample

lam.samp.vec<-c(lam.sample[lam.sample$Scenario == 'n=25; lam=2.4', 'Kept.data'],
                lam.sample[lam.sample$Scenario == 'n=25; lam=2', 'Kept.data'],
                lam.sample[lam.sample$Scenario == 'n=25; lam=2.8', 'Kept.data'])

lam.mat<-matrix(rep(lam.samp.vec, each=600), 
                nrow(plt.mat), ncol(plt.mat))

##for CV, make matrix of SD of estimates
plt.mat2<-cbind( p.hat[,1,2,2], outmat2[,1,2,1:5],
                 p.hat[,1,2,7], outmat2[,1,2,6:10],
                 p.hat[,1,2,12], outmat2[,1,2,11:15])

################################################################################
#### Plots Q2 ##################################################################

### Figure 2: bias and CV of lambda

##bias lam
df.lam<-data.frame(((plt.mat-lam.mat)/lam.mat)*100)
colnames(df.lam)<-c(paste( rep(1:3, each=6), c(0,nnests), sep='.' ))
meltQ2 <- melt(df.lam, id = NULL) 

coll2<-rep('2.4', nrow(meltQ2))
coll2[as.character(meltQ2$variable) %in% 
        c("2.0","2.10","2.25","2.50","2.100","2.250")]  <-'2.0'
coll2[as.character(meltQ2$variable) %in% 
        c("3.0","3.10","3.25","3.50","3.100","3.250")]  <-'2.8'
coll2<-factor(coll2,levels=c('2.4', '2.0', '2.8'))

pQ2a<-ggplot(meltQ2, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=coll2))+
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', linewidth=0.5)+
  labs(x = '# uncertain counts', y = expression(paste("Relative error in  ", lambda, sep=''))) +
  theme_bw()+
  theme(legend.position = 'none')+
  coord_cartesian(ylim=c(-20,20)) +
  scale_x_discrete(labels=rep(c(0,nnests), 3 ) )+
  geom_text(x=0.7, y=20, label='A)', 
            hjust='left',vjust='bottom',  size=4, color='black')+
  scale_fill_manual(name=expression(lambda),
                    #labels=c("2.4", "2.0", "2.8"),
                    values = c("2.4"='white',
                               "2.0"='lightgrey',
                               "2.8"='darkgrey'))+
  geom_segment(aes(x=0.5, y=-7, xend=6.5, yend=-7), color = 'blue')+
  geom_segment(aes(x=6.5, y=-6, xend=12.5, yend=-6), color = 'blue')+
  geom_segment(aes(x=12.5, y=-7, xend=18.5, yend=-7), color = 'blue')


##CV lam
df.lam<-data.frame((plt.mat2/plt.mat)*100)
colnames(df.lam)<-c(paste( rep(1:3, each=6), c(0,nnests), sep='.' ))
meltQ2b <- melt(df.lam, id = NULL) 

pQ2b<-ggplot(meltQ2b, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=coll2))+
  labs(x = '# uncertain counts', y = expression(paste("CV of  ", lambda, sep=''))) +
  theme_bw()+
  theme(legend.direction = "horizontal", legend.box = "vertical",
        legend.position='bottom') +
  scale_x_discrete(labels=rep(c(0,nnests), 3 ) ) +
  geom_text(x=0.7, y=9.8, label='B)', 
            hjust='left',vjust='bottom',  size=4, color='black') +
  scale_fill_manual(name=expression(paste(lambda, ': ', sep='')),
                    #labels=c("2.4", "2.0", "2.8"),
                    values = c("2.4"='white',
                               "2.0"='lightgrey',
                               "2.8"='darkgrey'))


##stack on top of each other
plotfullq2<-plot_grid(pQ2a, pQ2b, align = 'v', ncol=1 , rel_heights = c(0.85,1))

##this works for in-text figure
jpeg('Figure 2 bias CV lam Q2.jpg', width = 16, 
     height = 19, unit='cm', res=600)
plotfullq2
dev.off()

### Appendix figure S4:
### PLOT OF absolute bias in PNEST; SEPARATELY FOR EACH SCENARIO

bvec<-NULL
scenario<-NULL
idx<-NULL

for (sc in 1:nscen){
  
  if(sc %in% (1:5)){
    get.scen.double<-2
    sc.new<-sc+1}
  if(sc %in% (6:10)){
    get.scen.double<-7
    sc.new<-sc+2}
  if(sc %in% (11:15)){
    get.scen.double<-12
    sc.new<-sc+3}
  

  #true.pnest<-matrix(pnest.mat[scen[sc,2],], 100, 4, byrow=TRUE)

  ##for bias calculation use mean from kept sample from scenario 
  ## for paired counts
  ##values in pi.hat[,,2]
    true.pnest<-matrix(pi.hat[get.scen.double,,2], 100, 4, byrow=TRUE)
  
  bvec<-c(bvec, as.vector(outmat2[,2:5,1,sc]-true.pnest))
  scenario<-c(scenario, rep(sc.new, 400))
  idx<-rep(c('pi[1]', 'pi[2]', 'pi[3]', 'pi[4]'), each=100)
}
bdf<-data.frame(value=bvec, scenario=scenario,
                parameter=idx)

##add scenarios with n=0 (1=2.4, 7=2.0, 13=2.8)

omg<-rbind(p.hat[,2:5,1,2]-matrix(pi.hat[2,,2], 100, 4, byrow=TRUE),
      p.hat[,2:5,1,7]-matrix(pi.hat[7,,2], 100, 4, byrow=TRUE),
      p.hat[,2:5,1,12]-matrix(pi.hat[12,,2], 100, 4, byrow=TRUE)) 
omg<-data.frame(omg)
colnames(omg)<-c('pi[1]', 'pi[2]', 'pi[3]', 'pi[4]')
omg$scenario<-c(rep(1, 100), rep(7, 100), rep(13, 100))
omg.melt<-melt(omg, id = 'scenario') 
omg.melt2<-data.frame(value=omg.melt$value,
                      scenario=omg.melt$scenario,
                      parameter=omg.melt$variable)

bdf.full<-rbind(bdf, omg.melt2)
bdf.full$scenario<-as.factor(bdf.full$scenario)
bdf.full$parameter<-as.factor(bdf.full$parameter)

coll3<-rep('2.4', nrow(bdf.full))
coll3[as.numeric(bdf.full$scenario) %in% (7:12)]<-'2.0'
coll3[as.numeric(bdf.full$scenario) %in% (13:18)]<-'2.8'
coll3<-factor(coll3, levels = c('2.4', '2.0', '2.8'))

q2.pi<-ggplot(bdf.full, aes(x=scenario, y=value)) +
  geom_boxplot(aes(fill=coll3), fatten=0.5)+
  labs(x = '# uncertain counts', y = 'Error') +
  theme_bw() +
  scale_x_discrete(labels=rep(c(0,nnests), 3 ) ) +
  scale_fill_manual(name=expression(paste(lambda, ': ', sep='')),
                    #labels=c("2.4", "2.0", "2.8"),
                    values = c("2.4"='white',
                               "2.0"='lightgrey',
                               "2.8"='darkgrey'))+
  geom_hline(data=v_line,aes(yintercept = yintercept),
             linetype='dashed', color = 'orangered', linewidth=0.5) +
  theme(legend.position = 'bottom')+
  facet_grid(rows = vars(parameter),
             labeller = label_parsed)

jpeg('Figure S4 Bias pi Q2.jpg', width = 16, 
    height = 19, unit='cm', res=600)
q2.pi
dev.off()


### Figure S5: example figure of p, to compare with figure 1

##compile estimates of p from scenarios with addl ground counts
##for lam=2.4 only as an example (scenarios 1-5)
pvec<-NULL
n.id<-NULL
idx.p<-NULL

for (sc in 1:5){
  
  if(sc %in% (1:5)){
    get.scen.double<-2}
  # if(sc %in% (6:10)){
  #   get.scen.double<-7
  #   sc.new<-sc+2}
  # if(sc %in% (11:15)){
  #   get.scen.double<-12
  #   sc.new<-sc+3}
  
  pvec<-c(pvec, as.vector((outmat2[,ord,1,sc]-true.p)))
  n.id<-c(n.id, rep(nnests[scen[sc,1]], 1000))
  idx.p<-c(idx.p, rep(paste('p[', 1:10, ']', sep=''), each=100))
}
bdf.p<-data.frame(value=pvec, n=n.id,
                parameter=idx.p)

##get estimates of p from scenario 2 with no extra ground counts
omg.p<-p.hat[,ord,1,2]-true.p 

omg.p<-data.frame(omg.p)
colnames(omg.p)<-paste('p[', 1:10, ']', sep='')

omg.p$n<-c(rep(0, 100))

omgp.melt<-melt(omg.p, id = 'n') 

omgp.melt2<-data.frame(value=omgp.melt$value,
                      n=omgp.melt$n,
                      parameter=omgp.melt$variable)

bdfp.full<-rbind(bdf.p, omgp.melt2)
bdfp.full$n<-as.factor(bdfp.full$n)
bdfp.full$parameter<-factor(bdfp.full$parameter, 
                            levels = paste('p[', 1:10, ']', sep=''))

# New facet label names
nlabs <- paste ('# uncertain counts: ',c(0,nnests), sep='')
names(nlabs) <-as.character(c(0,nnests))

q2.p<-ggplot(bdfp.full, aes(x=parameter, y=value)) +
  geom_boxplot(fatten=0.5)+
  labs(x = NULL, y = 'Error') +
  theme_bw() +
  scale_x_discrete(labels=c('y=1|z=1', 'y=2|z=1', 
                            'y=1|z=2', 'y=2|z=2','y=3|z=2', 
                            'y=2|z=3','y=3|z=3','y=4|z=3',
                            'y=3|z=4', 'y=4|z=4')) +
  geom_hline(data=v_line,aes(yintercept = yintercept),
             linetype='dashed', color = 'orangered', linewidth=0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~n, nrow=3, ncol=2, labeller=labeller(n=nlabs))

jpeg('Figure S5 Bias p lam2.4 Q2.jpg', width = 16, 
    height = 19, unit='cm', res=600)
q2.p
dev.off()



################################################################################
######### Summaries, plots Q3 ##################################################

#### Summaries

##read in results from Q3 - no new data needed 
##n ground is always 50, n paired changes

ttq3<-readRDS('Parameter estimation/Sim results loglinear pfixed/Data.Q3.rds')
outmat.low<-ttq3$outmat.low
##outmat.low holds estimates for the new lambda=2 (paired observations for lam=2.4)
outmat.high<-ttq3$outmat.high
##outmat.high holds estimates for the new lambda=2.8 (paired observations for lam=2.4)

##adding dimension names to make navigation easier
dimnames(outmat.low)[[2]]<-dimnames(outmat.high)[[2]]<-params
dimnames(outmat.low)[[3]]<-dimnames(outmat.high)[[3]]<-
                          c('mean', 'sd', '2.5%', '50%', '97.5%', 'Rhat','n.eff',
                         'mode')
dimnames(outmat.low)[[4]]<-dimnames(outmat.high)[[4]]<-paste('n.paired=',nnests, sep='' )

#number of scenarios in this question
nscen2<-5

##remove non-converged iterations
keep.high<-matrix(NA, niter, nscen2)
for (j in 1:nscen2){
  sub<-outmat.high[,,6,j]
  maxR<-apply(sub, 1, max)
  keep.high[,j]<-maxR<=1.2
}
apply(keep.high,2,sum, na.rm=TRUE)
#for consistency, set NAs to false
keep.high[is.na(keep.high)]<-FALSE

keep.low<-matrix(NA, niter, nscen2)
for (j in 1:nscen2){
  sub<-outmat.low[,,6,j]
  maxR<-apply(sub, 1, max)
  keep.low[,j]<-maxR<=1.2
}
apply(keep.low,2,sum, na.rm=TRUE)
keep.low[is.na(keep.low)]<-FALSE


##make matrix with iteration numbers

niter.Q3<-data.frame(Scenario.ground=c(rep('n=50; lam=2.8', 5),
                                       rep('n=50; lambda=2', 5)),
                     Scenario.paired=rep(dimnames(p.hat)[[4]][1:5],2),
                     `Missing truths`=rep(100-apply(keep2, 2, sum)[1:5],2),
                     `Not converged` = c(apply(keep2, 2, sum)[1:5]-
                                           apply(keep.high,2,sum, na.rm=TRUE),
                                         apply(keep2, 2, sum)[1:5]-
                                           apply(keep.low,2,sum, na.rm=TRUE)),
                     `Total kept` = c(apply(keep.high,2,sum, na.rm=TRUE), 
                                      apply(keep.low,2,sum, na.rm=TRUE)),
                     check.names = FALSE)
# write_xlsx(niter.Q3, 'Parameter estimation/Tables for ms/Niter.Q3.xlsx')


##make new results matrices in which non-converged iterations are set to NA
outmat.low2<-outmat.low
outmat.high2<-outmat.high
for (i in 1:24){
outmat.low2[,i,1,][!keep.low]<-NA
outmat.low2[,i,2,][!keep.low]<-NA
outmat.high2[,i,1,][!keep.high]<-NA
outmat.high2[,i,2,][!keep.high]<-NA
}


##calculate summary results
Bias.Q3.h<-RMSE.Q3.h<-CV.Q3.h<-
  Bias.Q3.l<-RMSE.Q3.l<-CV.Q3.l<-matrix(NA, npar, nscen2)

#get input values high and low
tinp.h<-tinp.o[,11:15]
tinp.l<-tinp.o[,6:10]

#note: this works even though outmat has more parameters because it's
#the first 15 we're interested in
for (i in 1:npar){
  for (j in 1:nscen2){
    Bias.Q3.h[i,j]<-median(outmat.high2[,i,1,j]-tinp.h[i,j], na.rm=TRUE)#absolute bias
    RMSE.Q3.h[i,j]<-sum(((outmat.high2[,i,1,j]-tinp.h[i,j])^2)/niter.Q3$`Total kept`[j], 
                      na.rm=TRUE)
    
    Bias.Q3.l[i,j]<-median(outmat.low2[,i,1,j]-tinp.l[i,j], na.rm=TRUE)#absolute bias
    RMSE.Q3.l[i,j]<-sum(((outmat.low2[,i,1,j]-tinp.l[i,j])^2)/niter.Q3$`Total kept`[j+5], 
                        na.rm=TRUE)
    
    #only for lam, relative bias
    if (i == 1){
      Bias.Q3.h[i,j]<-median(((outmat.high2[,i,1,j]-tinp.h[i,j])/tinp.h[i,j])*100, na.rm=TRUE)
      CV.Q3.h[i,j]<-median((outmat.high2[,i,2,j]/outmat.high2[,i,1,j])*100, na.rm=TRUE)
      
      Bias.Q3.l[i,j]<-median(((outmat.low2[,i,1,j]-tinp.l[i,j])/tinp.l[i,j])*100, na.rm=TRUE)
      CV.Q3.l[i,j]<-median((outmat.low2[,i,2,j]/outmat.low2[,i,1,j])*100, na.rm=TRUE)
    }
  }
}
Bias.Q3.h<-round(Bias.Q3.h, dig=3); Bias.Q3.l<-round(Bias.Q3.l, dig=3) 
CV.Q3.h<-round(CV.Q3.h, dig=2); CV.Q3.l<-round(CV.Q3.l, dig=2)
RMSE.Q3.h<-round(RMSE.Q3.h, dig=3); RMSE.Q3.l<-round(RMSE.Q3.l, dig=3)

bias3h.df<-data.frame(Parameter=pars,Bias.Q3.h)
bias3l.df<-data.frame(Parameter=pars,Bias.Q3.l)

cv3h.df<-data.frame(Parameter=pars,CV.Q3.h)
cv3l.df<-data.frame(Parameter=pars,CV.Q3.l)

rmse3h.df<-data.frame(Parameter=pars,RMSE.Q3.h)
rmse3l.df<-data.frame(Parameter=pars,RMSE.Q3.l)

colnames(bias3h.df)<-colnames(cv3h.df)<-colnames(rmse3h.df)<-
    colnames(bias3l.df)<-colnames(cv3l.df)<-colnames(rmse3l.df)<-c('Parameter',
                                                               paste('n(paired)=', 
                                                               nnests, sep=''))
# write_xlsx(list(Bias.Q3.lam2.8 = bias3h.df, 
#                 CV.Q3.lam2.8 = cv3h.df,
#                 RMSE.Q3.lam2.8 = rmse3h.df,
#                 Bias.Q3.lam2 = bias3l.df, 
#                 CV.Q3.lam2 = cv3l.df,
#                 RMSE.Q3.lam2 = rmse3l.df), 'Parameter estimation/Tables for ms/Summary stats Q3.xlsx')


################################################################################
############################### Plots Q3 #######################################

##Figure 3 main ms: Estimates of lambda for all scenarios, Q3 (what else??)

##combine results from high and low lambda

plt.mat<-data.frame(cbind(((outmat.low2[,1,1,]-lam.low)/
                             lam.low)*100,
                          ((outmat.high2[,1,1,]-lam.high)/
                             lam.high)*100))
colnames(plt.mat)<-paste(rep(c('low', 'high'), each=5), rep(nnests, 2), sep='.')

meltQ3 <- melt(plt.mat, id = NULL) 

level<-factor(rep(c('2.0', '2.8'), each=500), levels=c('2.0', '2.8'))

pQ3<-ggplot(meltQ3, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=level))+
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', linewidth=0.5)+
  labs(x = '# paired counts', y = expression(paste("Relative error in  ", lambda, sep=''))) +
  theme_bw()+
  theme(legend.position = c(0.935,0.875),
        legend.box.background = element_rect(colour = "black"))+
  coord_cartesian(ylim=c(-20,20)) +
  scale_x_discrete(labels=rep(nnests, 2) )+
  scale_fill_manual(name=expression(lambda),
                    #labels=c("2.4", "2.0", "2.8"),
                    values = c("2.0"='lightgrey',
                               "2.8"='darkgrey'))+
  geom_segment(aes(x=0.5, y=-6, xend=5.5, yend=-6), color = 'blue')+
  geom_segment(aes(x=5.5, y=-7, xend=10.5, yend=-7), color = 'blue')

##this works for in-text figure
jpeg('Figure 3 Bias lam Q3.jpg', width = 16, 
     height = 10, unit='cm', res=600)
pQ3
dev.off()



### Supplementary figures

##Figure S6: Bias in pi

value<-NULL
lam.in<-NULL
n.in<-NULL
parm<-paste('pi[', 1:4, ']', sep='')
idx.in<-NULL

for (sc in 1:nscen2){
  
  tm.low<-matrix(p.nest.low, 100, 4, byrow=T)
  tm.high<-matrix(p.nest.high, 100, 4, byrow=T)
  plt.mat2<-cbind((outmat.low2[,2:5,1,sc]-tm.low),
                  (outmat.high2[,2:5,1,sc]-tm.high) )
  value<-c(value, as.vector(plt.mat2))
  lam.in<-c(lam.in, rep(c('2.0', '2.8'), each = 400))
  n.in<-c(n.in, rep(nnests[sc], 800))
  idx.in<-c(idx.in, rep(rep(parm, each=100), 2))
    
}

bdf.q3<-data.frame(value=value, parameter=idx.in,
                   n=n.in, lam=lam.in)
bdf.q3$n<-factor(bdf.q3$n, levels=as.character(nnests))

##new labels
lam.labs <- paste('lambda', ' = ', c('2.0', '2.8'), sep='')
names(lam.labs) <- c("2.0", "2.8")

q3.pi<-ggplot(bdf.q3, aes(x=n, y=value)) +
  geom_boxplot(fatten=0.5)+
  labs(x = '# paired counts', y = 'Error') +
  theme_bw() +
  #scale_x_discrete(labels=rep(c(0,nnests), 3 ) ) +
  # scale_fill_manual(name=expression(paste(lambda, ': ', sep='')),
  #                   #labels=c("2.4", "2.0", "2.8"),
  #                   values = c("2.4"='white',
  #                              "2.0"='lightgrey',
  #                              "2.8"='darkgrey'))+
  geom_hline(data=v_line,aes(yintercept = yintercept),
             linetype='dashed', color = 'orangered', linewidth=0.5) +
  #theme(legend.position = 'bottom')+
  # facet_grid(parameter~lam,
  #            labeller = label_parsed)
  facet_grid(parameter ~ glue('lambda*" = {lam}"'), 
             labeller = label_parsed)

jpeg('Figure S6 Bias pi Q3.jpg', width = 16, 
    height = 19, unit='cm', res=600)
q3.pi
dev.off()


##Figure S7: Bias in p

valuep<-NULL
lam.inp<-NULL
n.inp<-NULL
parmp<-c('y=1|z=1', 'y=2|z=1', 
         'y=1|z=2', 'y=2|z=2','y=3|z=2', 
         'y=2|z=3','y=3|z=3','y=4|z=3',
         'y=3|z=4', 'y=4|z=4')
idx.inp<-NULL


for (sc in 1:nscen2){
  
  plt.mat2<-cbind((outmat.low2[,ord,1,sc]-true.p),
                  (outmat.high2[,ord,1,sc]-true.p))
  
  valuep<-c(valuep, as.vector(plt.mat2))
  lam.inp<-c(lam.inp, rep(c('2.0', '2.8'), each = 1000))
  n.inp<-c(n.inp, rep(nnests[sc], 2000))
  idx.inp<-c(idx.inp, rep(rep(parmp, each=100), 2))
  
}

bdfp.q3<-data.frame(value=valuep, parameter=idx.inp,
                   n=n.inp, lam=lam.inp)
bdfp.q3$n<-factor(bdfp.q3$n, levels=as.character(nnests))
bdfp.q3$parameter<-factor(bdfp.q3$parameter,
                          levels=parmp)

q3.p<-ggplot(bdfp.q3, aes(x=parameter, y=value)) +
  geom_boxplot(fatten=0.5)+
  labs(x = NULL, y = 'Error') +
  theme_bw() +
  geom_hline(data=v_line,aes(yintercept = yintercept),
             linetype='dashed', color = 'orangered', linewidth=0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #theme(legend.position = 'bottom')+
  # facet_grid(n~lam)
  # all this extra cajigger is to keep factor level order as it should be
  facet_grid(forcats::fct_inorder(glue("'n = {n}'")) ~ glue('lambda*" = {lam}"'),
             labeller = label_parsed)

jpeg('Figure S7 Bias p Q3.jpg', width = 16, 
    height = 19, unit='cm', res=600)
q3.p
dev.off()



##############################################################################
############# Preparation for projections ####################################

##make table with pnest estimates for population projections

#1 three true pnest
#2 three expected observed pnest (based on pi and p)
#3 mean estimates of p.nest from Q3, all sample sizes, low, high
#  Note: not all under 3 are used for projections 

pnest.scenarios<-matrix(NA, nrow=16, ncol=4)

#true proportions of #nestlings
pnest.scenarios[1:3,]<-pnest.mat

#expected proportions of observed #nestlings
py<-matrix(NA,3,4)
for( jj in 1:3){
  for (i in 1:4){
    py[jj,i]<-0
    for(ii in 1:4){
      py[jj,i]<-py[jj,i]+pnest.mat[jj,ii]*p.obs[ii,i]
    }
  }
}
pnest.scenarios[4:6,]<-py

#Q3, low all sample sizes, high all sample sizes
for (jj in 1:5){
  pnest.scenarios[6+jj,]<-apply(outmat.low2[,2:5,1,jj],2,mean, na.rm=TRUE)
}
for (jj in 1:5){
  pnest.scenarios[11+jj,]<-apply(outmat.high2[,2:5,1,jj],2,mean, na.rm=TRUE)
}

rownames(pnest.scenarios)<-c(paste('True.lam.', c(2.4, 2, 2.8), sep=''),
                             paste('Expected.obs.lam.', c(2.4, 2, 2.8), sep=''),
                             paste('Q3.lam.2.ndouble.',nnests, sep=''),
                             paste('Q3.lam.2.8.ndouble.',nnests, sep=''))

#saveRDS(pnest.scenarios, 'Parameter estimation/Pnest.scenarios.rds')
