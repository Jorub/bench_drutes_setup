#run_all

# run drutes made benchmark
source("pso_full_vec.R")
# start mesh opti
source('drut_eval.R')

rndmint <- read.table("rndmint.dat", quote="\"", comment.char="")
seed=rndmint[42,1]
set.seed(seed)
write(seed,"rand_seed.out")

maxeval=8 # number of parallel oeprataions 
system(paste('./opti_setup.sh', maxeval))

xmin=c(0.05,1.5,0.3,10)
xmax=c(0.15,2.5,0.7,500)

PSO=PSO_fv_R(16,1,4,xmin=xmin,xmax=xmax,gen=200,maxeval=maxeval,printall=F)
write.csv(PSO,paste('bfnc_drut1D_PSO.csv',sep=''))

source("pso_full_vec_bn.R")
set.seed(seed)
write(seed,"rand_seedpso_bn.out")

xmin=c(0.05,1.5,0.3,10)
xmax=c(0.15,2.5,0.7,500)

PSO=PSO_fv_bn_R(80,1,4,xmin=xmin,xmax=xmax,gen=200,maxeval=maxeval,printall=F)
write.csv(PSO,paste('bfnc_drut1D_PSO_bn.csv',sep=''))

source("pso_full_vec_sce.R")
set.seed(seed)

maxeval=8 # number of parallel oeprataions 

xmin=c(0.05,1.5,0.3,10)
xmax=c(0.15,2.5,0.7,500)

PSO=PSO_fv_sce_R(80,2,4,xmin=xmin,xmax=xmax,gen=200,maxeval=maxeval,printall=F)
write.csv(PSO,paste('bfnc_drut1D_PSO_sce.csv',sep=''))

source("pso_full_vec_sce_bn.R")
set.seed(seed)
maxeval=8 # number of parallel oeprataions 


xmin=c(0.05,1.5,0.3,10)
xmax=c(0.15,2.5,0.7,500)

PSO=PSO_fv_sce_bn_R(80,2,4,xmin=xmin,xmax=xmax,gen=200,maxeval=maxeval,printall=F)
write.csv(PSO,paste('bfnc_drut1D_PSO_sce_bn.csv',sep=''))

source("tlbo_full_vec.R")
set.seed(seed)
maxeval=8 # number of parallel oeprataions 
xmin=c(0.05,1.5,0.3,10)
xmax=c(0.15,2.5,0.7,500)

TLBO=TLBO_fv_R(80,1,4,xmin=xmin,xmax=xmax,gen=200,maxeval=maxeval,printall=F,0.1)
write.csv(TLBO,paste('bfnc_drut1D_tlbo.csv',sep=''))

source("tlbo_full_vec_bn.R")
set.seed(seed)
maxeval=8 # number of parallel oeprataions 
xmin=c(0.05,1.5,0.3,10)
xmax=c(0.15,2.5,0.7,500)

TLBO=TLBO_bn_fv_R(80,1,4,xmin=xmin,xmax=xmax,gen=200,maxeval=maxeval,printall=F,0.1)
write.csv(TLBO,paste('bfnc_drut1D_tlbo_bn.csv',sep=''))

source("tlbo_full_vec_sce_bn.R")
set.seed(seed)
maxeval=8 # number of parallel oeprataions 
xmin=c(0.05,1.5,0.3,10)
xmax=c(0.15,2.5,0.7,500)

TLBO=TLBO_sce_bn_fv_R(80,2,4,xmin=xmin,xmax=xmax,gen=200,maxeval=maxeval,F,0.1)
write.csv(TLBO,paste('bfnc_drut1D_tlbo_sce_bn.csv',sep=''))

source("tlbo_full_vec_sce.R")
set.seed(seed)
maxeval=8 # number of parallel oeprataions 
xmin=c(0.05,1.5,0.3,10)
xmax=c(0.15,2.5,0.7,500)

TLBO=TLBO_fv_sce_R(80,2,4,xmin=xmin,xmax=xmax,gen=200,maxeval=maxeval,printall=F,0.1)
write.csv(TLBO,paste('bfnc_drut1D_tlbo_sce.csv',sep=''))
