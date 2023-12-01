################################################################################
### super simple population projection #########################################

## Sollmann, Wendt, Mattsson (submitted). Accounting for observation biases 
##   associated with counts of young when estimating fecundity: case study on 
##   the arboreal-nesting red kite (Milvus milvus). PCJ.

### Script to run population projection described in above referenced manuscript

library(writexl)
library(ggplot2)
library(cowplot)

## function to generate number of female nestlings from number of breeders
## and pi
breed.to.nest<-function(n.breeders, pnest){
  
  ##get number of nestlings from multinomial for all breeders
  nnest<-rmultinom(1, sum(n.breeders), pnest)
  
  ##generate females only form 50:50 sex ratio, for all nestlings combined
  nfemale<-rbinom(1, sum((1:4)*nnest), 0.5)
  
  ##return total number of female nestlings
  return(sum(nfemale))
}

################################################################################
### Projection input values based on Sergio et al. 2020, EcoApps 

### All input values and sources are summarized in Appendix 1: Table S5 

###read in values for pi (true input, observed, estimated)
### the file is on the repo but can also be recreated with the script
### 'Simulation summaries and plots.R'

##order of values in pnest.scenarios:
#true medium-low-high (3 rows);
##expected observed medium-low-high (3 rows);
#mean estimated with 50 ground counts LOW (5 double sample sizes)
#mean estimated with 50 ground counts HIGH (5 double sample sizes)

pnest.scenarios<-readRDS('data-raw/Pnest.scenarios.rds')
nscen<-nrow(pnest.scenarios)
niter=1000

##set starting (female) population size, time frame
N0<-100
Tt<-50

#starting age distribution
p0<-c(rep(0.35/2,2), rep(0.16/4, 4), 0.49)

##proportion reproducing in each age class
prepro<-c(0,0,0.9, 0.90,0.90,1,1)

###############################################################################
#### SIMULATION LOOP ##########################################################

##track population size by year and iteration for all scenarios
N.traj<-list()
set.seed(123)

##first, loop over population trajectories (increasing or decreasing)
for (xx in 1:2){
  ##set survival according to population trajectory
  if(xx==1){surv<-c(0.55, 0.8, rep(0.86, 5))}else{
    surv<-c(0.456, rep(0.706, 4), rep(0.874,2))
  }

  #N holds adult population size (ie, excluding nestlings)
  N<-array(NA, c(Tt, niter, nscen))
  #list keeps track of abundance per (adult) age class
  N.age.l<-list()

  ##then, loop over scenarios for pi (pnest)
  for (sc in 1:nscen){

    #select pi
    pnest<-pnest.scenarios[sc,]

    #create array w age specific abundances
    N.age<-array(NA, c(length(p0), Tt, niter))

    #loop over iterations
    for (iter in 1:niter){
      #create ages at time 1
      age0<-sample(1:length(p0), N0, p0, replace = TRUE)
      N.age0<-table(age0)
      if(length(N.age0)<length(p0)){
        ppp<-rep(0, length(p0))
        ppp[as.numeric(names(N.age0))]<-N.age0
        N.age0<-ppp
      }

      #matrix to hold number of breeders in each age class for all years
      breeders<-matrix(NA,length(N.age0), Tt )
      nestlings<-NULL

      #number of breeders yr 1
      breeders[,1]<-rbinom(length(N.age0),N.age0, prepro)

      #total number of nestlings produced yr1
      nestlings[1]<-breed.to.nest(breeders[,1], pnest)

      #fill in N for yr 1
      N.age[,1, iter]<-N.age0
      N[1,iter,sc]<-sum(N.age[,1, iter])

      ##loop over years >1
      for (t in 2:Tt){

        #### survival ####

        ##each age class is survival from previous year's earlier age class
        ##adults are sum of surviving 6-yrs and >6-yrs
        ##1-y are surviving nestlings
        N.age[1,t, iter]<-rbinom(1, nestlings[t-1], 0.456)

        for (j in 2:6){
          N.age[j,t, iter]<-rbinom(1, N.age[j-1,t-1, iter], surv[j-1])
        }

        N.age[7,t, iter]<-rbinom(1,N.age[7,t-1, iter],surv[7]) +
                          rbinom(1, N.age[6,t-1, iter], surv[6])

        #### breeding ####
        breeders[,t]<-rbinom(length(N.age0),N.age[,t, iter], prepro)
        nestlings[t]<-breed.to.nest(breeders[,t], pnest)

        #fill in N
        N[t,iter,sc]<-sum(N.age[,t, iter])

      }#year loop
    }#iteration loop
    #save age structure just in case...
    N.age.l[[sc]]<-N.age
  }#end scenario loop
  N.traj[[xx]]<-N #save all abundances for each trajectory
} #end trajectory type loop

##calculate summary stats of trajectories (mean, median, 2.5th, 97.5th perc.)
meanNtraj<-lapply(N.traj, function(x){
  meanN<-array(NA,c(nscen, Tt, 4))
  dimnames(meanN)[[1]]<-paste('Scenario', 1:nscen, sep='')
  dimnames(meanN)[[2]]<-paste('Year', 1:Tt, sep='')
  dimnames(meanN)[[3]]<-c('Mean', 'Median', '2.5%', '97.5%')
for(j in 1:nscen){
  meanN[j,,1]<-apply(x[,,j],1,mean)
  meanN[j,,2]<-apply(x[,,j],1,median)
  meanN[j,,3]<-apply(x[,,j],1,function(y){quantile(y, probs=0.025)})
  meanN[j,,4]<-apply(x[,,j],1,function(y){quantile(y, probs=0.975)})
}
  return(meanN)
  })

##adjust path as necessary
saveRDS(meanNtraj,'Population projection/Trajectory summaries.rds')
saveRDS(N.traj, 'Population projection/Full trajectories.rds')



################################################################################
#### plots and post-processing #################################################

####read in full and avg population trajectories if no longer in workspace
##
N.traj<-readRDS('Population projection/Full trajectories.rds')
meanNtraj<-readRDS('Population projection/Trajectory summaries.rds')

##list of length 2
##first element is growing populations, second is declining population
##each element is an array, year (1:50) by iteration (1:1000) by scenario (1:16)
## scenarios correspond to row names in table holding all input values for pi
rownames(pnest.scenarios)


#### Calculate difference in trajectory at time T=50
#### estimated (Q3) and observed against mean truth
#### for growing and declining pop for high and low lambda

##diff: array iteration by true lam (2.4, 2, 2.8) by growing/declining pop
diff<-array(NA, c(1000, 4, 2))
dimnames(diff)<-list(NULL, 
                     c('High.increase','High.decrease',
                       'Low.increase', 'Low.decrease'),
                     c('Obs.vs.true', 'Est.vs.true'))
##growing pop
## lam=2.8
  ## est vs true
diff[,1,2]<-(N.traj[[1]][50,,14]-meanNtraj[[1]][3,50,1])/meanNtraj[[1]][3,50,1]
    ## obs vs true
diff[,1,1]<-(N.traj[[1]][50,,6]-meanNtraj[[1]][3,50,1])/meanNtraj[[1]][3,50,1]

## lam=2.0
  ## est vs true
diff[,3,2]<-(N.traj[[1]][50,,9]-meanNtraj[[1]][2,50,1])/meanNtraj[[1]][2,50,1]
  ## obs vs true
diff[,3,1]<-(N.traj[[1]][50,,5]-meanNtraj[[1]][2,50,1])/meanNtraj[[1]][2,50,1]


##declining pop
## lam=2.8
## est vs true
diff[,2,2]<-(N.traj[[2]][50,,14]-meanNtraj[[2]][3,50,1])/meanNtraj[[2]][3,50,1]
## obs vs true
diff[,2,1]<-(N.traj[[2]][50,,6]-meanNtraj[[2]][3,50,1])/meanNtraj[[2]][3,50,1]

## lam=2.0
## est vs true
diff[,4,2]<-(N.traj[[2]][50,,9]-meanNtraj[[2]][2,50,1])/meanNtraj[[2]][2,50,1]
## obs vs true
diff[,4,1]<-(N.traj[[2]][50,,5]-meanNtraj[[2]][2,50,1])/meanNtraj[[2]][2,50,1]

##summarize into tables (mean and 95% interval)

diff.mean<-round(apply(diff, 2:3, median), dig=4)*100
diff.lower<-round(apply(diff, 2:3, function(x)quantile(x, p=0.025)), dig=4)*100
diff.upper<-round(apply(diff, 2:3, function(x)quantile(x, p=0.975)), dig=4)*100

Table1<-matrix(NA, nrow(diff.mean), ncol(diff.mean))

for (i in 1:nrow(diff.mean)){
  for (j in 1:ncol(diff.mean)){

    Table1[i,j]<-paste(diff.mean[i,j], ' (', diff.lower[i,j], ' - ', 
                       diff.upper[i,j], ')',
      sep='')
}}
colnames(Table1)<-colnames(diff.mean)

# write_xlsx(data.frame(Population = rownames(diff.mean), 
#                          Table1), 'Table 1.xlsx')


##### Plots of trajectories 

### summaries of simulated trajectories (list of 2)
### for high survival trajectories (growing populations)
### use first element in meanNtraj
### for low survival (declining pop), use second element

###each element is an array, scenario by year by summary stat
###scenarios correspond to pnest.scenarios


##### Plot trajectories with estimated pi with confidence intervals
##### then, add true and observed trajectories on top

### First, growing population
##compile all trajectories in a data frame
tr.df<-data.frame(est.28=meanNtraj[[1]][14,,1],
                  l.e28 = meanNtraj[[1]][14,,3],
                  u.e28 = meanNtraj[[1]][14,,4],
                  obs.28 = meanNtraj[[1]][6,,1],
                  true.28 = meanNtraj[[1]][3,,1],
                  est.20 = meanNtraj[[1]][9,,1],
                  l.e20 = meanNtraj[[1]][9,,3],
                  u.e20 = meanNtraj[[1]][9,,4], 
                  obs.20 = meanNtraj[[1]][5,,1],
                  true.20 = meanNtraj[[1]][2,,1],
                  year = 1:Tt)


pi<-ggplot(tr.df, aes(y=est.28, x=year))+
  geom_line(aes(color = '2.8', linetype = 'Estimated'))+
  geom_line(aes(y=obs.28, x=year,color= '2.8', linetype = 'Observed'))+
  geom_line(aes(y=true.28, x=year,color= '2.8', linetype = 'Truth'))+
  
  geom_line(aes(y=est.20, x=year,color= '2.0', linetype = 'Estimated'))+
  geom_line(aes(y=obs.20, x=year,color= '2.0', linetype = 'Observed'))+
  geom_line(aes(y=true.20, x=year,color= '2.0', linetype = 'Truth'))+
  
  labs(x = 'Year',y = 'N') +
  theme_bw() +
  
  scale_color_manual(name=expression(lambda),
                     values=c('2.0'='seagreen', '2.8'="red"))+
  scale_linetype_manual(name='Value',
                        values = c('Estimated'='solid', 
                                   'Observed'='dashed',
                                   'Truth' = 'dotted')) +
  
  geom_ribbon(aes(ymin=l.e28, ymax=u.e28),
              alpha=0.3, color=NA, fill="red") +
  geom_ribbon(aes(ymin=l.e20, ymax=u.e20),
              alpha=0.3, color=NA, fill='seagreen') +
  
  theme(legend.position ='none') +
  annotate('text', x=0, y=8750, label='A)',
           hjust='left',vjust='top',  size=5, color='black')


### repeat for declining population

##compile all trajectories in a data frame
tr.df2<-data.frame(est.28=meanNtraj[[2]][14,,1],
                  l.e28 = meanNtraj[[2]][14,,3],
                  u.e28 = meanNtraj[[2]][14,,4],
                  obs.28 = meanNtraj[[2]][6,,1],
                  true.28 = meanNtraj[[2]][3,,1],
                  est.20 = meanNtraj[[2]][9,,1],
                  l.e20 = meanNtraj[[2]][9,,3],
                  u.e20 = meanNtraj[[2]][9,,4], 
                  obs.20 = meanNtraj[[2]][5,,1],
                  true.20 = meanNtraj[[2]][2,,1],
                  year = 1:Tt)

pd<-ggplot(tr.df2, aes(y=est.28, x=year))+
  geom_line(aes(color = '2.8', linetype = 'Estimated'))+
  geom_line(aes(y=obs.28, x=year,color= '2.8', linetype = 'Observed'))+
  geom_line(aes(y=true.28, x=year,color= '2.8', linetype = 'Truth'))+
  
  geom_line(aes(y=est.20, x=year,color= '2.0', linetype = 'Estimated'))+
  geom_line(aes(y=obs.20, x=year,color= '2.0', linetype = 'Observed'))+
  geom_line(aes(y=true.20, x=year,color= '2.0', linetype = 'Truth'))+
  
  labs(x = 'Year',y = 'N') +
  theme_bw() +
  
  scale_color_manual(name=expression(paste(lambda, ':', sep='')),
                     values=c('2.0'='seagreen', '2.8'="red"))+
  # scale_linetype_manual(name='Value:',
  #                       values = c('Estimated'='solid',
  #                                  'Observed'='dashed',
  #                                  'Truth' = 'dotted')) +
  
  geom_ribbon(aes(ymin=l.e28, ymax=u.e28),
              alpha=0.3, color=NA, fill="red") +
  geom_ribbon(aes(ymin=l.e20, ymax=u.e20),
              alpha=0.3, color=NA, fill='seagreen') +
  
  theme(legend.direction = "horizontal", legend.box = "vertical",
        legend.position=c('bottom') )+
  
  annotate('text', x=0, y=250, label='B)',
           hjust='left',vjust='top',  size=5, color='black') 
  # guides(shape = 'none')

## could not figure out how to suppress linetype legend so needs to be cut manually

##stack on top of each other
plotfull3<-plot_grid(pi, pd, align = 'v', ncol=1, rel_heights = c(1,1.5) )

##this works for Peer journal
jpeg('Figure 4 Projections.jpg', width = 16, 
    height = 14.5, unit='cm', res=600)
plotfull3
dev.off()

