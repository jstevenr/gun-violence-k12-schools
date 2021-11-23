# 06/02/2019
# Move all `sim_inference` running commands to this script
# So that things are less messy


source("./sim_inference.R")

# the prior settings
pr = data.frame(count = rep(1,4), avg = c(0.05, 0.1, 0.005, 0.05))


#################################
####### single experiment #######
#################################

# 1. MLEs

# a de-coupled case
EL = stochastic_coevolve_infer2(N=100, model="SIR", seed=53) 

# another de-coupled case
EL = stochastic_coevolve_infer2(N=100, bet=0.03, gam=0.15, model="SIR", 
                                init.p = .07, alpha.r = .003,
                                seed=83) 

# 06/02/2019
# generate a decoupled case for a coupled inference
set.seed(71)
EL = stochastic_coevolve_infer2(N=100, bet=0.03, gam=0.12, model="SIR", 
                                alpha.r = .003, alpha.d = .07) 
saveRDS(EL, "./EpiNet_decoup_2.rds")

# MLE on a coupled process
set.seed(71)
EL2 = stochastic_coevolve_infer2(N=100, bet=0.1, gam = 0.15, model = "SIR",
                                 quarantine=T,
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05), return.infer = F)
# save the dataset: looks interesting (almost everybody got sick once)
# Event type counts: 93 94 601 47 46 536 254 39
# saveRDS(EL2, "~/Documents/EpiNet_coupled_2.rds")

# another coupled process
set.seed(73)
EL2 = stochastic_coevolve_infer2(N=100, bet=0.05, gam = 0.15, model = "SIR",
                                 quarantine=T,
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05), return.infer = F)
# Event type counts: 60 61 1151 38 9 1049 208 15 
# saveRDS(EL2, "~/Documents/EpiNet_coupled_3.rds")




# 2. Bayesian: w/ Gamma priors
pr = data.frame(count = rep(1,4), avg = c(0.05, 0.1, 0.005, 0.05))

# a coupled case
set.seed(73)
EL3 = stochastic_coevolve_infer2(N=100, bet=0.05, gam = 0.15, model = "SIR",
                                 quarantine=T,
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05), 
                                 Bayes = T, priors = pr,
                                 return.infer = F, samples = 1000)

# another coupled case
set.seed(1003)
EL3 = stochastic_coevolve_infer2(N=100, bet=0.04, gam = 0.15, model = "SIR",
                                 quarantine=T,
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05), 
                                 Bayes = T, priors = pr,
                                 return.infer = F, samples = 1000)
# Event type counts: 35 36 1095 15 2 1035 108 4 
# saveRDS(EL3, "~/Documents/EpiNet_coupled_5.rds")

# then another case
EL3 = stochastic_coevolve_infer2(N=100, bet=0.03, gam = 0.12, model = "SIR",
                                 quarantine=T,
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05), 
                                 Bayes = T, priors = pr,
                                 return.infer = F, samples = 1000)
# Event type counts: 25 26 747 16 4 777 108 8 
# saveRDS(EL3, "~/Documents/EpiNet_coupled_6.rds")



# 3. Try "hubnet": one person connected to all, all others Erdos-Renyi
# 06/02/2019
hubnet = readRDS("./hubnet100.rds")

# MLE
set.seed(71)
pdf("./hubnet100_plots.pdf", width=9, height=6)
EL_hub1 = stochastic_coevolve_infer2(N=100, bet=0.03, gam = 0.12, model = "SIR",
                                     quarantine=T,
                                     alpha.r = c(.005, 0.001, .005),
                                     alpha.d = c(.05, 0.1, .05), return.infer = F)
# Event type counts:
#   21 22 1287 19 0 1297 91 6 

# Bayesian
EL_hub2 = stochastic_coevolve_infer2(N=100, bet=0.03, gam = 0.12, 
                                     model = "SIR", quarantine=T,
                                     alpha.r = c(.005, 0.001, .005),
                                     alpha.d = c(.05, 0.1, .05), 
                                     Bayes = T, priors = pr, infer.interval = 1000,
                                     return.infer = F, samples = 1000,
                                     init.net = hubnet)
# Event type counts:
#   44 45 1136 34 9 1097 170 15 
# save the data (and inference results)
saveRDS(EL_hub2, "./EpiNet_coupled_hubnet2.rds")

dev.off()



# 4. Try with larger social networks
# N = 500, for example

# MLEs
set.seed(73)
pdf("./N500_MLEs.pdf",height = 6, width = 9)
EL_big1 = stochastic_coevolve_infer2(N=500, bet=0.03, gam = 0.12, model = "SIR",
                                 quarantine=T,
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05), return.infer = T)
dev.off()
# MLEs with 50549 events in total:(which is good estimation)
#   0.03006962 0.1232847 
#   0.004999269 0.001043554 0.004855036 0.0499428 0.0994639 0.05106905 
# Event type counts: (pretty much everyone got sick once)
#   497 498 20770 1210 1745 17227 6855 1747 
saveRDS(EL_big1,"./EpiNet_coupled_big_1.rds")


# Bayesian
pdf("./N500_Bayes_2.pdf",height = 6, width = 9)
EL_big2 = stochastic_coevolve_infer2(N=500, bet=0.03, gam = 0.12, model = "SIR",
                                 quarantine=T,
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05), 
                                 Bayes = T, priors = pr,
                                 return.infer = F, samples = 2000,
                                 infer.interval = 10000)
dev.off()
# Posterior means and variances with 53466 events in total:
#   0.02883917 0.1237284 0.00500658 0.0009646837 0.005125877 0.04902018 0.09998641 0.05088475 
# 1.673436e-06 3.074038e-05 1.121113e-09 8.368837e-10 1.431859e-08 1.293245e-07 1.448042e-06 1.533012e-06 
# Event type counts:
#   496 497 22357 1111 1834 18580 6903 1688 


#################################
##### repeated experiments ######
#################################

# Bayesian inference
res1 = rep_stochastic_coevolve_infer(N=100, bet=0.03, gam = 0.12, model = "SIR",
                                     quarantine=T,
                                     alpha.r = c(.005, 0.001, .005),
                                     alpha.d = c(.05, 0.1, .05),
                                     seed=73, return.infer=T,
                                     Bayesian = T, priors = pr, 
                                     demo.sleep = F, rep= 4,
                                     event.thres = 1000,
                                     pdfpath = "./ex1_all",
                                     savedat = F)

res2 = rep_stochastic_coevolve_infer(N=100, bet=0.03, gam = 0.12, model = "SIR",
                                     quarantine=T,
                                     alpha.r = c(.005, 0.001, .005),
                                     alpha.d = c(.05, 0.1, .05),
                                     seed=83, return.infer=T,
                                     Bayesian = T, priors = pr, 
                                     demo.sleep = F, rep= 5,
                                     event.thres = 1000,
                                     pdfpath = "./ex2",
                                     savedat = T, fname = "./exdat2")

# MLEs
res3 = rep_stochastic_coevolve_infer(N=100, bet=0.03, gam = 0.12, model = "SIR",
                                     quarantine=T,
                                     alpha.r = c(.005, 0.001, .005),
                                     alpha.d = c(.05, 0.1, .05),
                                     seed=83, return.infer=T,
                                     Bayesian = F, priors = pr, 
                                     demo.sleep = F, rep= 5,
                                     event.thres = 1000,
                                     pdfpath = "./ex3",
                                     savedat = T, fname = "./ex3dat")

# Note ex3dat_5:
# Event type ounts: 47 48 621 35 13 573 189 17 
# MLE plots (w/ CIs) saved


# try the `hubnet`
# Bayesian
res_hub = rep_stochastic_coevolve_infer(N=100, bet=0.03, gam = 0.12, model = "SIR",
                                     quarantine=T, init.net = hubnet,
                                     alpha.r = c(.005, 0.001, .005),
                                     alpha.d = c(.05, 0.1, .05),
                                     seed=83, return.infer=T,
                                     Bayesian = T, priors = pr, 
                                     demo.sleep = F, rep= 4,
                                     event.thres = 1000,
                                     pdfpath = "./hubnet2",
                                     savedat = F)

# MLE
res_hub2 = rep_stochastic_coevolve_infer(N=100, bet=0.03, gam = 0.12, model = "SIR",
                                     quarantine=T,
                                     alpha.r = c(.005, 0.001, .005),
                                     alpha.d = c(.05, 0.1, .05),
                                     seed=83, return.infer=T,
                                     Bayesian = F, priors = pr, 
                                     demo.sleep = F, rep= 5,
                                     event.thres = 1000,
                                     pdfpath = "./hubnet_MLEs",
                                     savedat = T, fname = "./hubnet2")

