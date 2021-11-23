# 08/05/2019
# Inference using random mixing SIR, static network SIR, and decoupled network SIR


#setwd("~/Documents/EpiNet/")
# run this under RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# OR run this when sourcing this file
# setwd(getSrcDirectory()[1])

source("./sim_inference.R")

library(dplyr)

parse_events2 <- function(G0, I0, events, model="SIR"){
  
  N = nrow(G0); epid = rep(0,N)
  adjmat =  G0; infected =  I0; epid[infected] = 1
  
  susceptible = which(epid==0); infected = which(epid==1)
  N.S = sum(epid==0); N.I = sum(epid==1)
  if(model == "SIR"){ N.R = 0 }
  Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
  M = c(nnzero(adjmat[susceptible,susceptible])/2, 
        nnzero(adjmat[susceptible,infected]), 
        nnzero(adjmat[infected,infected])/2)
  N.SI = M[2]
  
  # data storage
  n.E = sum(events$event == 1)
  n.R = sum(events$event == 2)
  C = rep(0,3); D = rep(0,3)
  
  # FOR NOW
  # assume the case of "quarantine"
  big.sums = numeric(8) 
  
  events = events[,-c(5:6)]
  
  # parse through the data
  t.cur = 0
  for(r in 1:nrow(events)){
    t.next = events$time[r]
    del.t = t.next - t.cur
    
    # update "big.sums"
    big.sums = big.sums + c(N.SI, N.I, Mmax-M, M) * del.t
    
    z =  events$event[r]
    if (z==1){
      # infection
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      epid[p1] = 1
      N.S = N.S - 1
      N.I = N.I + 1
      # update counts
      if(model=="SIR"){
        R.p1 = nnzero(adjmat[p1, which(epid==-1)])
        Mmax = c((N.S+N.R)*(N.S+N.R-1)/2, (N.S+N.R) * N.I, N.I*(N.I-1)/2)
        M[1] = M[1] - S.p1 - R.p1
        M[2] = M[2] + S.p1 + R.p1 - I.p1
        M[3] = M[3] + I.p1
        N.SI = N.SI + S.p1 - I.p1
      }else{
        Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
        M[1] = M[1] - S.p1
        M[2] = M[2] + S.p1 - I.p1
        M[3] = M[3] + I.p1
        N.SI = M[2]
      }
    }else if (z==2){
      # recovery
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      N.I = N.I - 1
      if(model=="SIR"){
        R.p1 = nnzero(adjmat[p1, which(epid==-1)])
        epid[p1] = -1
        N.R = N.R + 1
        M[1] = M[1] + S.p1 + R.p1
        M[2] = M[2] - S.p1 - R.p1 + I.p1
        M[3] = M[3] - I.p1
        Mmax = c((N.S+N.R)*(N.S+N.R-1)/2, (N.S+N.R) * N.I, N.I*(N.I-1)/2)
        N.SI = N.SI - S.p1
      }else{
        epid[p1] = 0
        N.S = N.S + 1
        M[1] = M[1] + S.p1
        M[2] = M[2] - S.p1 + I.p1
        M[3] = M[3] - I.p1
        Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
        N.SI = M[2]
      }
    }else{
      # some edge stuff
      p1 = events$per1[r]
      p2 = events$per2[r]
      ## NEED SOME KIND OF CODING HERE ##
      ## TO MARK RECONNECT/DISCONNECT  ##
      ##                               ##
      if(z %in% c(3:5)){
        # reconnection
        adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
        # bookkeeping
        if (epid[p1]!=1 & epid[p2]!=1){
          # S-S/S-R/R-R type
          C[1] = C[1] + 1
          M[1] = M[1] + 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          C[3] = C[3] + 1
          M[3] = M[3] + 1
        }else if (epid[p1]!=-1 & epid[p2]!=-1){
          # S-I type: neither is recovered
          C[2] = C[2] + 1
          M[2] = M[2] + 1
          N.SI = N.SI + 1
        }else{
          # R-I type
          C[2] = C[2] + 1
          M[2] = M[2] + 1
        }
      }else{
        # disconnection
        adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
        # bookkeeping
        if (epid[p1]!=1 & epid[p2]!=1){
          # S-S/S-R/R-R type
          D[1] = D[1] + 1
          M[1] = M[1] - 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          D[3] = D[3] + 1
          M[3] = M[3] - 1
        }else if (epid[p1]!=-1 & epid[p2]!=-1){
          # S-I type
          D[2] = D[2] + 1
          M[2] = M[2] - 1
          N.SI = N.SI - 1
        }else{
          # R-I type
          D[2] = D[2] + 1
          M[2] = M[2] - 1
        }
      }
      
    }
    # report progress
    # cat("Processing augmented dataset...", 
    #     r/nrow(events),"\r")
    t.cur = t.next
    
    # make inference once in a while
    
  }
  #cat("\nDone.\n")
  
  res = list(event.counts = c(n.E, n.R, C, D), big.sums = big.sums)
  
  return(res)
}


quick_check <- function(resmat){
  apply(resmat, 2, function(v) c(mean(v),sd(v), quantile(v, c(.025,.975))))
}



# run repetitions
REP = 50
event.thres = 400

betas = NULL

for(re in 1:REP){
  dat = stochastic_coevolve_infer2(N=50, bet=0.05, gam=0.12, model="SIR",
                                   quarantine = T, alpha.r = c(.005, 0.001, .005),
                                   alpha.d = c(.05, 0.1, .05),
                                   infer.interval = 1000, infer.verbose = F, 
                                   infer.plot = F, return.infer = F)
  while(nrow(dat$events) < event.thres){
    dat = stochastic_coevolve_infer2(N=50, bet=0.05, gam=0.12, model="SIR",
                                     quarantine = T, alpha.r = c(.005, 0.001, .005),
                                     alpha.d = c(.05, 0.1, .05),
                                     infer.interval = 1000, infer.verbose = F, 
                                     infer.plot = F, return.infer = F)
  }
  
  # a) do everything correctly
  PE = parse_events2(G0 = dat$G0, I0 = dat$I0, events = dat$events)
  beta1 = PE$event.counts[1]/PE$big.sums[1]
  
  # b) assume fixed network
  Events = dat$events %>% filter(event %in% c(1,2))
  PE = parse_events2(G0 = dat$G0, I0 = dat$I0, events = Events)
  beta2 = PE$event.counts[1]/PE$big.sums[1]
  
  # c) assume complete graph (random mixing)
  GG0 = dat$G0; GG0[,] = 1
  PE = parse_events2(G0 = GG0, I0 = dat$I0, events = Events)
  beta3 = PE$event.counts[1]/PE$big.sums[1]
  
  betas = rbind(betas, c(beta1, beta2, beta3))
}

betas = as.data.frame(betas)
names(betas) = c("dynamic", "static", "no.net")

quick_check(betas)


# 10/10/2019
# add histograms
compare_infer <- function(res, name, truth=0.05){
  
  xl = range(c(res, truth))
  hist(res, main=name, breaks=10,
      xlab = expression(hat(beta)), ylab=NULL,
      xlim=xl, col="lightblue", cex.main=3)
  abline(v=truth,col="red",lwd=3)
  
}

par(mar=c(5,2,2,2))

titles = c("dynamic net", "static net", "random mixing")

pdf("./net_assumptions.pdf")
for(i in 1:3){
  ti = titles[i]
  compare_infer(betas[,i], ti)
}
dev.off()

