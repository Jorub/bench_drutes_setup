
ref1=read.table("ref1.out", quote="\"")
timeref=ref1[,1]
ref1=ref1[,3]
ref2=read.table("ref2.out", quote="\"")
ref2=ref2[,3]
ref3=read.table("ref3.out", quote="\"")
ref3=ref3[,3]


eval_fun= function(ln_id){
  result=matrix(ncol=1,nrow=ln_id)
 
  for(i in 1:ln_id){ 
    t=0
    #while(t!=770){
      obs1=read.table(paste(i,"/out/obspt_RE_matrix-1.out",sep=''), quote="\"",comment.char="#", sep="")
      #timeobs1=obs1[,1]
      #t=length(timeobs1)
      #print(i)
   # }
    timeobs1=obs1[,1]
    print(length(timeobs1))
    obs1=obs1[,3]
    obs2=read.table(paste(i,"/out/obspt_RE_matrix-2.out",sep=''), quote="\"",comment.char="#", sep="")
    timeobs2=obs2[,1]
    obs2=obs2[,3]
    obs3=read.table(paste(i,"/out/obspt_RE_matrix-3.out",sep=''), quote="\"",comment.char="#", sep="")
    timeobs3=obs3[,1]
    obs3=obs3[,3]

    res1=RMSE(ref1,obs1,timeobs1)
    res2=RMSE(ref2,obs2,timeobs2)
    res3=RMSE(ref3,obs3,timeobs3)
    result[i,1]=1/3*sum(res1,res2,res3)
  }
  
  return(result)#
}

RMSE=function(real,sim, simtime){
  ln_time=length(timeref)-1
  newsim=approx(y=sim,x=simtime,xout=timeref[1:ln_time])
  result=sqrt(sum((newsim$y-real[1:ln_time])^2)/ln_time)
  return(result)
}

