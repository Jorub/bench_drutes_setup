# TLBO with shuffling mechanism and bad neighbourhood
# Learning experience with teaching learning based optimization or LETLBO Zou et al. (2015)
TLBO_sce_bn_fv_R=function(class_size,classes,dim,xmin,xmax,gen,printall=T,maxeval,start_shuffle_prob=0.995){
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
  
  if(classes<2){stop('You are using a shuffling mechanism. Your number of classes should at least be two')}
  
  if(class_size<20*classes){
    readline("Your class size is low (equivilent to population). There's a high risk to get stuck in while loops with low class size due to strategy seperation of the population. You've been warned. Press enter if you want to continue but better restart with higher class size.") 
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
  ## groups according to ranks
  ## groups
  group_complex=cut(ranks,breaks=seq(0,class_size,length=classes+1),labels=1:classes)# length of pop group 1 contains lowest ranks
  rank_classes=tapply(result_old,group_complex,rank) # makes a list of size classes
  lbests=tapply(result_old,group_complex,min)
  lbests_loc=school_old[match(lbests,result_old,nomatch=FALSE),]
  result_new=matrix(ncol=classes,nrow=class_size)
  #bad nhood
  bad_hood=matrix(nrow=ceiling(length(ranks)/2),ncol=dim*2)
  centr=school_old[which(ranks>floor(class_size/2)),]
  bad_ranks=ranks[which(ranks>floor(class_size/2))]
  centr=centr[rev(sort(bad_ranks,index.return=T)$ix),]
  sph=seq(from=1,to=0.3/(ceiling(class_size/2)),length=ceiling(class_size/2))
  hood=as.matrix(sph)%*%t(as.matrix(xmax))
  bad_hood=cbind(centr-hood,centr+hood)
  modeall=T
  switch=0
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
    teachers_loc_mat=matrix(rep(teacher_loc,class_size),nrow=class_size,byrow=T)
    if(modeall){
      go_to_loc=teachers_loc_mat # global best
      switch=switch+1
      if(switch>10){
        modeall=F
        switch=0
        group_complex=cut(ranks,breaks=seq(0,class_size,length=classes+1),labels=1:classes)# length of pop group 1 contains lowest ranks
        rank_classes=tapply(result_old,group_complex,rank) # makes a list of size classes
      }
    }else{
      go_to_loc=lbests_loc[group_complex,] # group_complex
      switch=switch+1
      if(switch>10){modeall=T;switch=0}
    }
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
    mean_loc_mat=matrix(rep(mean_loc,length(first)),nrow=length(first),byrow=T)           
    school_new[first,]=school_old[first,]+runif(length(first))*(go_to_loc[first,]-TF*mean_loc_mat)
    # strategy b
    other=c(length(second))
    other_cplx=tapply(class_vec[second],group_complex[second],sample)
    for(i in 1:classes){
      other[group_complex[second]==i]=other_cplx[[i]]
    }

    while(any(other==class_vec[second])){
      other_cplx=tapply(class_vec[second],group_complex[second],sample)
      for(i in 1:classes){
        other[group_complex[second]==i]=other_cplx[[i]]
      }
      
    }
    # which is better
    other_better=which(result_old[second]>result_old[other])
    # make id based on result 
    secnd_id=second
    secnd_id[other_better]=other[other_better]
    school_new[second,]=school_old[second,]+runif(length(second))*(go_to_loc[second,]-school_old[secnd_id,])
    # test if in bad neighbourhood
    in_bad_hood=rep(FALSE,class_size)
    for(n in 1:(round(length(ranks)/2))){
      bool_vec=apply(school_new,1,function(x) sum(x>=bad_hood[n,1:dim]&x<=bad_hood[n,(dim+1):(2*dim)])==dim)
      in_bad_hood[bool_vec]=TRUE
    }
    ln_bn=length(in_bad_hood[in_bad_hood])
    if(ln_bn>0){
      rndm_gb=matrix(runif(dim*ln_bn),nrow=ln_bn,byrow=T)
      school_new[in_bad_hood,]=school_new[in_bad_hood,]-rndm_gb*(go_to_loc[in_bad_hood,]-school_new[in_bad_hood,])
    }
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
    ### new bad neighbourhood
    centr=school_old[which(ranks>floor(class_size/2)),]
    bad_ranks=ranks[which(ranks>floor(class_size/2))]
    centr=centr[rev(sort(bad_ranks,index.return=T)$ix),]
    sph=seq(from=0.1+(1-0.1)*(gen-k)/gen,to=0.3/(ceiling(class_size/2))*(gen-k)/gen,length=ceiling(class_size/2))
    hood=as.matrix(sph)%*%t(as.matrix(xmax))
    bad_hood=cbind(centr-hood,centr+hood)
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
    other=c(length(first))
    other_cplx=tapply(class_vec[first],group_complex[first],sample)
    for(i in 1:classes){
      other[group_complex[first]==i]=other_cplx[[i]]
    }
    while(any(other==class_vec[first])){
      other_cplx=tapply(class_vec[first],group_complex[first],sample)
      for(i in 1:classes){
        other[group_complex[first]==i]=other_cplx[[i]]
      }
      
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
    other_one=c(length(second))
    other_cplx=tapply(class_vec[second],group_complex[second],sample)
    for(i in 1:classes){
      other_one[group_complex[second]==i]=other_cplx[[i]]
    }
    while(any(other_one==class_vec[second])){
      other_cplx=tapply(class_vec[second],group_complex[second],sample)
      for(i in 1:classes){
        other_one[group_complex[second]==i]=other_cplx[[i]]
      }
      
    }
    other_two=c()
    other_cplx=tapply(class_vec[second],group_complex[second],sample)
    for(i in 1:classes){
      other_two[group_complex[second]==i]=other_cplx[[i]]
    }
    while(any(other_two==class_vec[second]|other_two==other_one)){
      other_cplx=tapply(class_vec[second],group_complex[second],sample)
      for(i in 1:classes){
        other_two[group_complex[second]==i]=other_cplx[[i]]
      }
      
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
    # test if in bad neighbourhood
    in_bad_hood=rep(FALSE,class_size)
    for(n in 1:(round(length(ranks)/2))){
      bool_vec=apply(school_new,1,function(x) sum(x>=bad_hood[n,1:dim]&x<=bad_hood[n,(dim+1):(2*dim)])==dim)
      in_bad_hood[bool_vec]=TRUE
    }
    ln_bn=length(in_bad_hood[in_bad_hood])
    if(ln_bn>0){
      rndm_gb=matrix(runif(dim*ln_bn),nrow=ln_bn,byrow=T)
      school_new[in_bad_hood,]=school_new[in_bad_hood,]-rndm_gb*(go_to_loc[in_bad_hood,]-school_new[in_bad_hood,])
    
    }
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
    ### new bad neighbourhood
    
    centr=school_old[which(ranks>floor(class_size/2)),]
    bad_ranks=ranks[which(ranks>floor(class_size/2))]
    centr=centr[rev(sort(bad_ranks,index.return=T)$ix),]
    hood=as.matrix(sph)%*%t(as.matrix(xmax))
    bad_hood=cbind(centr-hood,centr+hood)
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
    # updating lbests if applicable
    if(!modeall){
      lbests_new=tapply(result_old,group_complex,min)
      lbests_loc_new=school_old[match(lbests_new,result_old,nomatch=FALSE),]
      #replacement id
      rep_id=which(lbests_new<lbests)
      lbests=replace(lbests,rep_id,lbests_new[rep_id])
      lbests_loc[rep_id,]=lbests_loc_new[rep_id,]
    }
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
    reshuffle_prob=reshuffle_prob*start_shuffle_prob
    a=runif(1)
    # reinitilize worst half if prob reached
    if(a>reshuffle_prob){
      worst=length(ranks[which(ranks>floor(class_size/2))])
      school_old[which(ranks>floor(class_size/2)),]=matrix(t(runif(worst*dim,xmin,xmax)),nrow=worst,byrow=T) 
      reshuffle_prob=start_shuffle_prob
      }
      
    k=k+1
  }
  if(printall){
    results[[k]]=cbind(t(teacher_loc),teacher)
  }else{
    resulttemp=cbind(t(teacher_loc),teacher)
    results=rbind(results,resulttemp)
  }
  return(results)
}
