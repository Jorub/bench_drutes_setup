# TLBO no shuffling, no bad neighbourhood approach
# Learning experience with teaching learning based optimization or LETLBO Zou et al. (2015)

TLBO_fv_R=function(class_size,classes,dim,xmin,xmax,gen,printall=T,maxeval,start_shuffle_prob=0.995){
  source("drut_eval.R")
  fit_func<-function(ln_id){
    result=eval_fun(ln_id)
    return(result) # this should be maxeval rows if combined opti is used
  }
  error<- 1e-8
  
  if(length(xmax)!=dim){
    stop('xmax: not enough boundaries defined. Dimension and length(xmax) differ')
  }
  if(length(xmin)!=dim){
    stop('xmin: not enough boundaries defined. Dimension and length(xmin) differ')
  }
  
  if(class_size<20){
    readline("Your class size is low (equivilent to population). There's a high risk to get stuck in while loops with low classsize due to strategy seperation of the population. You've been warned. Press enter if you want to continue but better restart with higher class size.") 
  }
  # Initialise learners 
  school_old=matrix(t(runif(class_size*dim,xmin,xmax)),nrow=class_size,byrow=T) 
  
  teacher=c()
  means=c()
  result_old=c()
  teacher_loc=matrix(ncol=dim,nrow=classes)
  mean_loc=matrix(ncol=dim,nrow=classes)
  # evaluate learners and define first teacher
  index=sort(rep(1:ceiling((class_size/maxeval)),maxeval))
  index=index[1:class_size]
  for(i in 1:ceiling(class_size/maxeval)){
    ln_id=(length(index[index==i]))
    pars_in=cbind(rep('p',ln_id),matrix(school_old[index==i,],nrow=ln_id))
    write(t(pars_in),'pars.in',append = F,ncol=dim+1)
    if(i>1){
      result_old=rbind(result_old,fit_func(ln_id))
    }else{
      result_old=fit_func(ln_id)
    }
  }
  teacher=min(result_old)
  pos=which.min(result_old)
  teacher_loc=school_old[pos,]
  mean_loc=colMeans(school_old)
  ranks=rank(result_old,ties.method = "random")
  result_new=matrix(ncol=classes,nrow=class_size)
  # optimisation algorithm update 
  k=1
  if(printall){
    results=list() # to be returned in the end
  }else{
    results=matrix(ncol=dim,nrow=1)
  }
  class_vec=1:class_size
  school_new=matrix(ncol=dim,nrow=class_size)
  reshuffle_prob=1
  stuck=F
  stuck_ind=0
  while((teacher>error) && (k<=gen) && !stuck){
    print(k)
    # two strategies in teachers phase, determined through a and b
    a=runif(class_size)
    b=runif(class_size)
    first=which(a<b)
    second=which(a>=b)
    while(length(first)<3&length(second)<3){
      a=runif(class_size)
      b=runif(class_size)
      first=which(a<b)
      second=which(a>=b)
    }
    # strategy a
    TF=round(1+runif(length(first)))
    teacher_loc_mat=matrix(rep(teacher_loc,length(first)),nrow=length(first),byrow=T)
    mean_loc_mat=matrix(rep(mean_loc,length(first)),nrow=length(first),byrow=T)           
    school_new[first,]=school_old[first,]+runif(length(first))*(teacher_loc_mat-TF*mean_loc_mat)
    # strategy b
    other=sample(class_vec[second])
    while(any(other==class_vec[second])){
      other=sample(class_vec[second])
    }
    # which is better
    other_better=which(result_old[second]>result_old[other])
    # make id based on result 
    secnd_id=second
    secnd_id[other_better]=other[other_better]
    teacher_loc_mat=matrix(rep(teacher_loc,length(second)),nrow=length(second),byrow=T)
    school_new[second,]=school_old[second,]+runif(length(second))*(teacher_loc_mat-school_old[secnd_id,])
    # parallel evaluation
    for(i in 1:ceiling(class_size/maxeval)){
      ln_id=(length(index[index==i]))
      pars_in=cbind(rep('p',ln_id),matrix(school_new[index==i,],nrow=ln_id))
      write(t(pars_in),'pars.in',append = F,ncol=dim+1)
      if(i>1){
        result_new=rbind(result_new,fit_func(ln_id))
      }else{
        result_new=fit_func(ln_id)
      }
    }
    #
    newbetter=which(result_new<result_old)
    result_old[newbetter]=result_new[newbetter]
    school_old[newbetter,]=school_new[newbetter,]
    ranks=rank(result_old,ties.method = 'random')
    # end of teachers phase
    ############################
    #
    #############################
    
    # beginning of learners' phase 
    # two strategies in learners' phase, determined through a and b
    a=runif(class_size)
    b=runif(class_size)
    first=which(a<b)
    second=which(a>=b)
    while(length(first)<4&length(second)<4){
      a=runif(class_size)
      b=runif(class_size)
      first=which(a<b)
      second=which(a>=b)
    }
    # first strategy, comparison to one other student
    other=sample(class_vec[first])
    while(any(other==class_vec[first])){
      other=sample(class_vec[first])
    }
    # which is better
    other_better=which(result_old[first]>result_old[other])
    # make ids based on result for substraction 
    # when other is better it's school_old[other,]-school_old[first,] and vice versa
    other_id=other
    other_id[other_better]=first[other_better]
    #
    frst_id=first
    frst_id[other_better]=other[other_better]
    school_new[first,]=school_old[first,]+runif(length(first))*(school_old[frst_id,]-school_old[other_id,])
    # 2nd strategy: comparing two other students
    other_one=sample(class_vec[second])
    while(any(other_one==class_vec[second])){
      other_one=sample(class_vec[second])
    }
    other_two=sample(class_vec[second])
    while(any(other_two==class_vec[second]&other_two==other_one)){
      other_two=sample(class_vec[second])
    }
    id_other=which(result_old[other_one]<result_old[other_two])
    # make ids based on result for substraction 
    # when other is better it's school_old[other,]-school_old[first,] and vice versa
    one_better=other_two
    one_better[id_other]=other_one[id_other]
    two_better=other_one
    two_better[id_other]=other_two[id_other]
    #
    school_new[second,]=school_old[second,]+runif(length(second))*(school_old[one_better,]-school_old[two_better,])

    # parallel evaluation
    for(i in 1:ceiling(class_size/maxeval)){
      ln_id=(length(index[index==i]))
      pars_in=cbind(rep('p',ln_id),matrix(school_new[index==i,],nrow=ln_id))
      write(t(pars_in),'pars.in',append = F,ncol=dim+1)
      if(i>1){
        result_new=rbind(result_new,fit_func(ln_id))
      }else{
        result_new=fit_func(ln_id)
      }
    }

    # replacing school with better solution
    newbetter=which(result_new<result_old)

    result_old[newbetter]=result_new[newbetter]
    school_old[newbetter,]=school_new[newbetter,]
    ranks=rank(result_old,ties.method = 'random')

    #calculating new teacher
    if((teacher-min(result_old))<0.01){
      stuck_ind=stuck_ind+1
    }else{
      stuck_ind=0
    }
    if(stuck_ind>=100){
      stuck=T
    }
    teacher=min(result_old)
    pos=which.min(result_old)
    teacher_loc=school_old[pos,]
    ranks=rank(result_old,ties.method = 'random')
    mean_loc=colMeans(school_old)
    if(printall){
      results[[k]]=cbind(school_old,result_old)
    }else{
      if(k==1){
        results=cbind(t(teacher_loc),teacher)
      }
      else{
        resulttemp=cbind(t(teacher_loc),teacher)
        results=rbind(results,resulttemp)
      }
    }
    k=k+1
    reshuffle_prob=reshuffle_prob*start_shuffle_prob
    a=runif(1)
    # reinitilize worst half if prob reached
    if(a>reshuffle_prob){
      worst=length(ranks[which(ranks>floor(class_size/2))])
      school_old[which(ranks>floor(class_size/2)),]=matrix(t(runif(worst*dim,xmin,xmax)),nrow=worst,byrow=T) 
      reshuffle_prob=start_shuffle_prob
    }
    
  }
  if(printall){
    results[[k]]=cbind(t(teacher_loc),teacher)
  }else{
    resulttemp=cbind(t(teacher_loc),teacher)
    results=rbind(results,resulttemp)
  }

  return(results)
}

