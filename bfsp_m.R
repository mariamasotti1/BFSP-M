library(fields);library(FastGP);library(ggplot2);library(sparseMVN);library(extraDistr);library(Matrix);library(Rfast);library(zoo)

ProgressBar = function(cur, total,...){
  n_equals = min(floor((cur/total) * (options("width")$width - 20)),
                 options("width")$width - 20)
  n_space = options("width")$width - n_equals - 20
  spinny_thing = rep(c("-", "\\", "|", "/"), each = 1)[(n_equals %% 4) + 1]
  if(cur > total){
    status = "  Done   \n\n"
    spinny_thing = "+"
  }else{
    status = "  Loading"
  }
  
  cat(c("\r[",
            rep("=", n_equals),
            rep("-", n_space),
            "]  ", spinny_thing, status), sep = ""
          )
}


prior_birth=function(prop_param,param,mu_out,x,y,prop_X,X,xlim,ylim,N){
  XB<-prop_X%*%prop_param[5:(5+2*N)]
  bp<-cbind(XB*sin(a2c(x,y,prop_param[1:2]))+prop_param[1],XB*cos(a2c(x,y,prop_param[1:2]))+prop_param[2])
  prop_cluster<-+(d2c(x,y,prop_param[1:2])<=XB)
  if(any(bp[,1]<xlim[1])|any(bp[,1]>xlim[2])|any(bp[,2]<ylim[1])|any(bp[,2]>ylim[2])){return(-Inf)}
  cluster<-prop_cluster
  for(n in which(complete.cases(param))){
    cluster<-cluster+(d2c(x,y,param[n,1:2])<=X[[n]]%*%param[n,5:15])
  }
  if(any(cluster>1)){return(-Inf)}
  
  return(sum(dnorm(prop_param[6:15],mean=c(rep(0,10)),sd=c(rep(1,10)),log=T))+dhnorm(prop_param[3]-mu_out, sigma = 10, log = TRUE)+dgamma(1/prop_param[4],1,1,log=T))
}

prior_death=function(prop_param,param,mu_out,x,y,prop_X,X,xlim,ylim,N){
  XB<-prop_X%*%prop_param[5:(5+2*N)]
  bp<-cbind(XB*sin(a2c(x,y,prop_param[1:2]))+prop_param[1],XB*cos(a2c(x,y,prop_param[1:2]))+prop_param[2])
  if(any(bp[,1]<xlim[1])|any(bp[,1]>xlim[2])|any(bp[,2]<ylim[1])|any(bp[,2]>ylim[2])){return(-Inf)}
  cluster<-rep(0,length(x))
  for(n in which(complete.cases(param))){
    cluster<-cluster+(d2c(x,y,param[n,1:2])<=X[[n]]%*%param[n,5:15])
  }
  if(any(cluster>1)){return(-Inf)}
  
dunif(prop_param[5],min=0,max=.5,log=T)+sum(dnorm(prop_param[6:15],mean=c(rep(0,10)),sd=c(rep(1,10)),log=T))+dhnorm(prop_param[3]-mu_out, sigma = 10, log = TRUE)+dgamma(1/prop_param[4],1,1,log=T)
}

proposal_birth=function(prop_param,param0,clust,mu_out){
  -log(length(which(colSums(clust,na.rm=T)==0)))+dhnorm(prop_param[3]-param0[1],10,log=T)+dgamma(prop_param[4],1,1,log=T)+dunif(prop_param[5],0,.5,log=T)
}

a2c<-function(x,y,center){
  shift.x<-x-center[1]
  shift.y<-y-center[2]
  a2c<-atan2(shift.x,shift.y)
  a2c<-ifelse(a2c<0,a2c+2*pi,a2c)
  return(a2c)
}

d2c<-function(x,y,center){
  d2c<-sqrt((x-center[1])^2+(y-center[2])^2)
  return(d2c)
}


likelihood=function(param,param0,value,x,y,Sigma,X,N,num_lesion,spatial){
  if(num_lesion>0){
    index<-which(complete.cases(param))
    par<-matrix(param[index,],nrow=num_lesion,ncol=(5+2*N))
    cluster<-matrix(nrow=num_lesion,ncol=length(x))
    eX<-X[c(index)]
    for(n in 1:num_lesion){
      cluster[n,]<-ifelse(d2c(x,y,par[n,1:2])<=eX[[n]]%*%par[n,5:(5+2*N)],1,0)
    }
  }else(cluster<-array(0,dim=c(1,length(x))))

  log.lik<-sapply(0:num_lesion,calc.lik,value=value,x=x,y=y,Sigma,param0=param0,param=par,cluster=cluster,spatial=spatial)
  return(log.lik)
}

calc.lik<-function(n,value,x,y,Sigma,param0,param,cluster,spatial){
  
 
    dmatrix<-Sigma[[1]]
    wend<-Sigma[[2]]
    if(n==0){
      ind0<-which(colSums(cluster)==0)
      x0<-x[ind0]
      y0<-y[ind0]
      center0<-c(mean(x0),mean(y0))
      
      index<-list()
      index[[1]]<-ind0[which(x0<center0[1]&y0<center0[2])]
      index[[2]]<-ind0[which(x0<center0[1]&y0>=center0[2])]
      index[[3]]<-ind0[which(x0>=center0[1]&y0<center0[2])]
      index[[4]]<-ind0[which(x0>=center0[1]&y0>=center0[2])]
      log.lik<-0
      for(i in 1:4){
        ind<-index[[i]]
        Y0<-value[ind]
        n0<-length(Y0)
        if(spatial==F){cov<-as(as(as(param0[2]*diag(n0), "dMatrix"), "symmetricMatrix"), "CsparseMatrix")}else{
        cov<-as(as(as(param0[3]*exp(-param0[4]*dmatrix[ind,ind])*wend[ind,ind]+param0[2]*diag(n0), "dMatrix"), "symmetricMatrix"), "CsparseMatrix")
        }
        if(n0>0){log.lik<-log.lik+ dmvn.sparse(t(Y0), rep(param0[1],length(Y0)), Cholesky(cov), prec = F, log = TRUE)}else{log.lik<-log.lik+0}
      }
    }else{
      ind1<-which(cluster[n,]==1)
      Y1<-value[ind1]
      n1<-length(Y1)
      if(n1<2) return(-Inf)
      if(spatial==F){cov<-as(as(as(param[n,4]*param0[2]*diag(n1), "dMatrix"), "symmetricMatrix"), "CsparseMatrix")}else{
      cov<-as(as(as(param[n,4]*exp(-param0[4]*dmatrix[ind1,ind1])*wend[ind1,ind1]+param0[2]*diag(n1), "dMatrix"), "symmetricMatrix"), "CsparseMatrix")
      }
      log.lik<- dmvn.sparse(t(Y1), rep(param[n,3],length(Y1)), Cholesky(cov), prec = F, log = TRUE)
      
    }
  

  return(as.numeric(log.lik))
}




prior_param=function(param,param0,x,y,X,xlim,ylim,N,num_lesion){
  cluster<-c(rep(0,length(x)))
  for(n in which(complete.cases(param))){
    XB<-X[[n]]%*%param[n,5:(5+2*N)]
    bp<-cbind(x=XB*sin(a2c(x,y,param[n,1:2]))+param[n,1],y=XB*cos(a2c(x,y,param[n,1:2]))+param[n,2])
    cluster<-cluster+(d2c(x,y,param[n,1:2])<=X[[n]]%*%param[n,5:15])
    if(any(bp[,1]<xlim[1])|any(bp[,1]>xlim[2])|any(bp[,2]<ylim[1])|any(bp[,2]>ylim[2])){return(-Inf)}
  }
  if(any(cluster>1)){return(-Inf)}
  par<-matrix(param[complete.cases(param),],nrow=num_lesion)
  prior=sum(dgamma(1/param[,4],1,1,log=T),na.rm=T)+sum(dhnorm(param[,3]-param0[1], sigma = 10, log = TRUE),na.rm=T)+sum(dunif(param[,5],min=0,max=.5,log=T),na.rm=T)+sum(dnorm(param[,6:15],mean=0,sd=1,log=T),na.rm=T)
  
  
  return(prior)
}

param_update=function(param,param0,la,gam,i,value,x,y,ll,Sigma,X,xlim,ylim,N,num_lesion,spatial){
  for(n in which(complete.cases(param))){
    
    if(i<=50) {
      param_star<-propose_new(param[n,3:4],(.01/2)*diag(2))
    }else{
      param_star<-.99*propose_new(param[n,3:4],exp(la[[n]][1])*gam[[n]][[1]])+.01*propose_new(param[n,3:4],(.01/2)*diag(2))}
    parameters_star<-param
    parameters_star[n,3:4]<-param_star
    prior_ratio=prior_param(parameters_star,param0,x,y,X,xlim,ylim,N,num_lesion)-prior_param(param,param0,x,y,X,xlim,ylim,N,num_lesion)
    if(is.na(prior_ratio)){prior_ratio=-Inf}
    if(prior_ratio>-Inf){
      ll_star<-likelihood(parameters_star,param0,value,x,y,Sigma,X,N,num_lesion,spatial)
      like_ratio=sum(ll_star)-sum(ll)
      alpha=exp(prior_ratio+like_ratio)
      if(runif(1)<alpha){
        param[n,3:4]<-param_star
        ll<-ll_star
      } 
    }
    
    if(i<=50) {
      param_star<-propose_new(param[n,5:10],(.01/6)*diag(6))
    }else{
      param_star<-.99*propose_new(param[n,5:10],exp(la[[n]][2])*gam[[n]][[2]])+.01*propose_new(param[n,5:10],(.01/6)*diag(6))}
    parameters_star<-param
    parameters_star[n,5:10]<-param_star
    prior_ratio=prior_param(parameters_star,param0,x,y,X,xlim,ylim,N,num_lesion)-prior_param(param,param0,x,y,X,xlim,ylim,N,num_lesion)
    if(is.na(prior_ratio)){prior_ratio=-Inf}
    if(prior_ratio>-Inf){
      ll_star<-likelihood(parameters_star,param0,value,x,y,Sigma,X,N,num_lesion,spatial)
      like_ratio=sum(ll_star)-sum(ll)
      alpha=exp(prior_ratio+like_ratio)
      if(runif(1)<alpha){
        param[n,5:10]<-param_star
        ll<-ll_star
      }
    }
    
    for(p in 2:(ceiling(2*N/5))){
      index<-(5*p+1):min((5*p+5),length(param[n,]))
      if(i<=50) {
        param_star<-propose_new(param[n,index],(.01/length(index))*diag(length(index)))
      }else{
        param_star<-.99*propose_new(param[n,index],exp(la[[n]][p+1])*gam[[n]][[p+1]])+.01*propose_new(param[n,index],(.01/length(index))*diag(length(index)))}
      parameters_star<-param
      parameters_star[n,index]<-param_star
      prior_ratio=prior_param(parameters_star,param0,x,y,X,xlim,ylim,N,num_lesion)-prior_param(param,param0,x,y,X,xlim,ylim,N,num_lesion)
      if(is.na(prior_ratio)){prior_ratio=-Inf}
      if(prior_ratio>-Inf){
        ll_star<-likelihood(parameters_star,param0,value,x,y,Sigma,X,N,num_lesion,spatial)
        like_ratio=sum(ll_star)-sum(ll)
        alpha=exp(prior_ratio+like_ratio)
        if(runif(1)<alpha){
          param[n,index]<-param_star
          ll<-ll_star
        }
      }
    }
    
  }
  return(list(param,ll))
}


birth_func=function(param0,param,num_lesion,value,x,y,ll,X,Sigma,xlim,ylim,N,P,spatial){
  if(num_lesion>0){
    clust<-array(dim=c(10,length(x)))
    for(c in which(complete.cases(param))){
      clust[c,]<-+(d2c(x,y,param[c,1:2])<=X[[c]]%*%param[c,5:15])
    }
    q<-sample(which(colSums(clust,na.rm=T)==0),1) 
  }else{q<-sample(1:length(x),1)}

  prop_size<-runif(1,0,.5)
  prop_param<-c(x[q],y[q],param0[1]+rhnorm(1,10),rgamma(1,1,1),prop_size,rep(0,10))
  

  prop_X<-cbind2(cbind2(1,sin(outer(a2c(x,y,prop_param[1:2]),(1:N)))),cos(outer(a2c(x,y,prop_param[1:2]),(1:N))))
  prior_ratio<-prior_birth(prop_param,param,param0[1],x,y,prop_X,X,xlim,ylim,N)

  
  if(prior_ratio==-Inf){return(list(param,num_lesion,ll,X))}
  if(num_lesion==0){
    pmm<-log(P[num_lesion+2,num_lesion+1])-log(P[num_lesion+1,num_lesion+2])
  }else{pmm<-log(P[num_lesion+1,num_lesion])-log(P[num_lesion,num_lesion+1])}
  
  if(num_lesion==0){
    clust<-array(0,dim=c(1,length(x)))
  }else{clust<-assign_clust(param,x,y,X)}
  
  proposal_ratio<--proposal_birth(prop_param,param0,clust,param0[1])+pmm
  
  param_star<-param
  
  param_star[min(which(!complete.cases(param))),]<-prop_param
  X_star<-X
  X_star[[min(which(!complete.cases(param)))]]<-prop_X
  ll_star<-likelihood(param_star,param0,value,x,y,Sigma,X_star,N,num_lesion+1,spatial)
  like_ratio<-sum(ll_star)-sum(ll)
  alpha=exp(prior_ratio+proposal_ratio+like_ratio)
  #alpha=exp(like_ratio)
  if(is.na(alpha)){alpha=0}
  if(runif(1)<alpha){
    num_lesion<-num_lesion+1
  
    X[[min(which(complete.cases(param)==F))]]<-prop_X
    
    param[min(which(complete.cases(param)==F)),]<-prop_param
    ll<-ll_star
    
  }
  return(list(param,num_lesion,ll,X))
}



death_func=function(param,param0,num_lesion,value,x,y,ll,X,Sigma,xlim,ylim,N,P,spatial){
  delete=sample(which(complete.cases(param)),1)
  param_star<-param
  param_star[delete,]<-NA
  X_star<-X
  X_star[[delete]]<-NA
  ll_star<-likelihood(param_star,param0,value,x,y,Sigma,X_star,N,num_lesion-1,spatial)
  
  prior_ratio<--prior_death(param[delete,],param,param0[1],x,y,X[[delete]],X,xlim,ylim,N)
  if(prior_ratio==-Inf){list(param,num_lesion,ll,X)}
  like_ratio<-sum(ll_star)-sum(ll)
  if(num_lesion==1){
    pmm<-log(P[num_lesion+1-1,num_lesion+1])-log(P[num_lesion+1,num_lesion+1-1])
  }else{pmm<-log(P[num_lesion-1,num_lesion])-log(P[num_lesion,num_lesion-1])}
  
  proposal_ratio<-pmm-log(num_lesion)
  alpha=exp(prior_ratio+proposal_ratio+like_ratio)
  #alpha=exp(like_ratio)
  if(is.na(alpha)){alpha=0}
  if(runif(1)<alpha){
    num_lesion<-num_lesion-1
   
    param[delete,]<-NA
    X[[delete]]<-NA
    ll<-ll_star

  }
  return(list(param,num_lesion,ll,X))
}

merge_function=function(param,param0,num_lesion,value,x,y,ll,X,Sigma,xlim,ylim,N,P,spatial){
 
  ind<-which(complete.cases(param))
  merge<-ind[sample(length(ind),2)]
  
  X1<-X[[merge[1]]]
  X2<-X[[merge[2]]]
  XB1<-X1%*%param[merge[1],5:15]
  XB2<-X2%*%param[merge[2],5:15]
  cluster<-+(d2c(x,y,param[merge[1],1:2])<XB1|d2c(x,y,param[merge[2],1:2])<XB2)
  cstar<-c(mean(x[cluster==1]),mean(y[cluster==1]))
  bp1<-cbind.data.frame(x=XB1*sin(a2c(x,y,param[merge[1],1:2]))+param[merge[1],1],y=XB1*cos(a2c(x,y,param[merge[1],1:2]))+param[merge[1],2])
  bp2<-cbind.data.frame(x=XB2*sin(a2c(x,y,param[merge[2],1:2]))+param[merge[2],1],y=XB2*cos(a2c(x,y,param[merge[2],1:2]))+param[merge[2],2])

  bpstar<-rbind(bp1,bp2)
  bpstar<-bpstar[d2c(bpstar[,1],bpstar[,2],cstar)>median(d2c(bpstar[,1],bpstar[,2],cstar)),]
  theta_star<-a2c(bpstar[,1],bpstar[,2],cstar)
  f_star<-d2c(bpstar[,1],bpstar[,2],cstar)
 
  beta_star<-try.solve(solve(t(X_str)%*%X_str)%*%t(X_str)%*%f_star)
  if(any(is.na(beta_star))){return(list(param,num_lesion,ll,X))}
  adj1<-c(rnorm(1,0,.01))
  adj<-rep(0,11)
  adj[sample(1:11,1)]<-adj1
  beta_star1<-beta_star+ adj

  mu_star<-rnorm(1,mean(param[merge,3]),1)
  var_star<-rnorm(1,mean(param[merge,4]),1)

  param_star<-param
  param_star[merge[1],]<-c(cstar,mu_star,var_star,beta_star1)
  param_star[merge[2],]<-NA
  X_star<-X
  X_star[[merge[1]]]<-cbind2(cbind2(1,sin(outer(a2c(x,y,cstar),(1:5)))),cos(outer(a2c(x,y,cstar),(1:5))))
  X_star[[merge[2]]]<-NA
  prior_ratio=prior_param(param_star,param0,x,y,X_star,xlim,ylim,N,num_lesion-1)-prior_param(param,param0,x,y,X,xlim,ylim,N,num_lesion)
  if(prior_ratio==-Inf|is.na(prior_ratio)){return(list(param,num_lesion,ll,X))}
  pmm<-log(P[num_lesion-1,num_lesion])-log(P[num_lesion,num_lesion-1])
  proposal_ratio<-pmm+log(choose(num_lesion,2))-dnorm(mu_star,mean(param[merge,3]),1,log=T)-dnorm(var_star,mean(param[merge,4]),1,log=T)-dnorm(adj1,0,.01,log=T)-log(1/11)
  ll_star<-likelihood(param_star,param0,value,x,y,Sigma,X_star,N,num_lesion-1,spatial)
  like_ratio=sum(ll_star)-sum(ll)
  alpha=exp(prior_ratio+proposal_ratio+like_ratio)
  #alpha=exp(like_ratio)
  if(is.na(alpha)){alpha=0}
  if(runif(1)<alpha){
    num_lesion<-num_lesion-1
 
    param<-param_star
    X<-X_star
    ll<-ll_star
   
  }
  return(list(param,num_lesion,ll,X))
}

try.solve <- function(code, silent = FALSE) {
  tryCatch(code, error = function(c) {
    if (!silent) {return(NA)}
    else{code}})}


split_function=function(param,param0,num_lesion,value,x,y,ll,X,Sigma,xlim,ylim,N,P,spatial){

  ind<-which(complete.cases(param))
  split<-ind[sample(length(ind),1)]
 
  basis<-cbind2(cbind2(1,sin(outer(seq(0,2*pi,length.out=200),(1:5)))),cos(outer(seq(0,2*pi,length.out=200),(1:5))))
  XB<-basis%*%param[split,5:15]
  xz <- as.zoo(XB)
  rxz <- rollapply(xz, 3, function(x) which.min(x)==2)
  minima<-index(rxz)[coredata(rxz)]
  if(length(minima)<2){return(list(param,num_lesion,ll,X))}
  mins=sample(minima,2)

  b<-cbind.data.frame(x=XB*sin(seq(0,2*pi,length.out=200))+param[split,1],y=XB*cos(seq(0,2*pi,length.out=200))+param[split,2])
  line<-lm(b[mins,2]~b[mins,1])$coef

  basis<-cbind2(cbind2(1,sin(outer(a2c(x,y,param[split,1:2]),(1:5)))),cos(outer(a2c(x,y,param[split,1:2]),(1:5))))
  XB<-basis%*%param[split,5:15]
  cluster<-+(d2c(x,y,param[split,1:2])<XB)
  new.clust<-ifelse(cluster==0,0,ifelse(line[1]+line[2]*x<=y,1,2))
  if(length(which(new.clust==1))==0|length(which(new.clust==2))==0){return(list(param,num_lesion,ll,X))}
 
  c1<-c(mean(x[new.clust==1]),mean(y[new.clust==1]))
  c2<-c(mean(x[new.clust==2]),mean(y[new.clust==2]))
  
  bp1<-b[which(line[1]+line[2]*b[,1]+.01<=b[,2]),]
  bp2<-b[which(line[1]+line[2]*b[,1]-.01>b[,2]),]
  
  line.pts<-cbind(x=seq(min(b[mins,1]),max(b[mins,1]),length.out = 100),y=line[1]+line[2]*seq(min(b[mins,1]),max(b[mins,1]),length.out = 100))
  
  line.pts1<-line.pts
  line.pts1[,2]<-line.pts1[,2]+.01
  line.pts2<-line.pts
  line.pts2[,2]<-line.pts2[,2]-.01
  bp1<-rbind(bp1,line.pts1)
  bp2<-rbind(bp2,line.pts2)
  

  theta_star<-a2c(bp1[,1],bp1[,2],c1)
  f_star<-d2c(bp1[,1],bp1[,2],c1)
  X_star<-cbind2(cbind2(1,sin(outer(theta_star,(1:N)))),cos(outer(theta_star,(1:N))))

  beta_star11<-try.solve(solve(t(X_star)%*%X_star)%*%t(X_star)%*%f_star)
  adj1<-c(rnorm(1,0,.01))
  adj<-rep(0,11)
  adj[sample(1:11,1)]<-adj1
  beta_star1<-beta_star11+ adj

  theta_star<-a2c(bp2[,1],bp2[,2],c2)
  f_star<-d2c(bp2[,1],bp2[,2],c2)
  X_star<-cbind2(cbind2(1,sin(outer(theta_star,(1:N)))),cos(outer(theta_star,(1:N))))

  beta_star22<-try.solve(solve(t(X_star)%*%X_star)%*%t(X_star)%*%f_star)
  adj2<-c(rnorm(1,0,.01))
  adj<-rep(0,11)
  adj[sample(1:11,1)]<-adj2
  beta_star2<-beta_star22+ adj

  if(any(is.na(beta_star1))|any(is.na(beta_star2))){return(list(param,num_lesion,ll,X))}

  mu1<-rnorm(1,param[split,3],1)
  mu2<-rnorm(1,param[split,3],1)
  v1<-rnorm(1,param[split,4],1)
  v2<-rnorm(1,param[split,4],1)

  param_star<-param
  param_star[split,]<-c(c1,mu1,v1,beta_star1)
  param_star[min(which(!complete.cases(param))),]<-c(c2,mu2,v2,beta_star2)
  X_star<-X
  X_star[[split]]<-cbind2(cbind2(1,sin(outer(a2c(x,y,c1),(1:5)))),cos(outer(a2c(x,y,c1),(1:5))))
  X_star[[min(which(!complete.cases(param)))]]<-cbind2(cbind2(1,sin(outer(a2c(x,y,c2),(1:5)))),cos(outer(a2c(x,y,c2),(1:5))))
  prior_ratio=prior_param(param_star,param0,x,y,X_star,xlim,ylim,N,num_lesion+1)-prior_param(param,param0,x,y,X,xlim,ylim,N,num_lesion)
  
  if(prior_ratio==-Inf|is.na(prior_ratio)){return(list(param,num_lesion,ll,X))}
  pmm<-log(P[num_lesion+1,num_lesion])-log(P[num_lesion,num_lesion+1])
  proposal_ratio=-(-log(num_lesion)+choose(length(minima),2)+dnorm(v1,param[split,4],1,log=T)+dnorm(v2,param[split,4],1,log=T)+dnorm(mu1,param[split,3],1,log=T)+dnorm(mu2,param[split,3],1,log=T)+dnorm(adj1,0,.01,log=T)+dnorm(adj2,0,.01,log=T)+log(1/11)+log(1/11))+pmm
  ll_star<-likelihood(param_star,param0,value,x,y,Sigma,X_star,N,num_lesion+1,spatial)
  like_ratio=sum(ll_star)-sum(ll)
  alpha=exp(prior_ratio+proposal_ratio+like_ratio)
  #alpha=exp(like_ratio)
  if(is.na(alpha)){alpha=0}
  if(runif(1)<alpha){
    num_lesion<-num_lesion+1
    param<-param_star
    X<-X_star
    ll<-ll_star
    
  }
  return(list(param,num_lesion,ll,X))
}


assign_clust=function(param,x,y,X){
  center<-param[,1:2]
  clust<-array(dim=c(10,length(x)))
  for(c in which(complete.cases(param))){
    clust[c,]<-+(d2c(x,y,center[c,])<=X[[c]]%*%param[c,5:15])
  }
  return(clust)
}

propose_new<-function(param,cov){
  prop<-rcpp_rmvnorm(1,cov,param)
  return(prop)
}

prior_param0=function(param0){
  prior=(dnorm(param0[1],0,100,log=T)+sum(dgamma(1/param0[2:3],.1,.1,log=T)))+dgamma(param0[4],3,.5,log=T)
  return(prior)
}


param0_update=function(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial){
  if(i<=50) {
    param0_star<-c(propose_new(param0,(.01/4)*diag(4)))
  }else{
    param0_star<-c(.99*propose_new(param0,exp(la0)*gam0)+.01*propose_new(param0,(.01/4)*diag(4)))}
  
  prior_ratio=prior_param0(param0_star)-prior_param0(param0)
  if(prior_ratio>-Inf){ 
    ll_star<-likelihood(param,param0_star,value,x,y,Sigma,X,N,num_lesion,spatial)
    like_ratio=sum(ll_star)-sum(ll)
    alpha=exp(prior_ratio+like_ratio)
    if(runif(1)<alpha){
      param0<-param0_star
      ll<-ll_star
    } 
  }
  return(list(param0,ll))
}


bd_bfsp<-function(x,y,value,iterations=10000,burn=2000,max_lesion=10,spatial=T){
  x<-range01(x)
  y<-range01(y)
  value<-scale(value)
  min_lesion=0
  loc<-cbind(x,y)
  dmatrix<-as.matrix(dist(loc))
  wend<-Wendland(dmatrix,theta=.05,dimension=2,k=2)
  Sigma<-list(dmatrix,wend)
  N<-5
  param<-matrix(ncol=15,nrow=max_lesion)
  param0<-c(-1,1,1,1)
  mu0<-param0
  gam0<-diag(4)
  la0<--1
  la<-list()
  for(n in 1:max_lesion){la[[n]]<-c(-1,-1,-1)}
  gam<-list()
  for(n in 1:max_lesion){gam[[n]]<-list(diag(2),diag(6),diag(5))}
  mu<-array(dim=c(max_lesion,15),data=0)
  param0_keeps<-array(dim=c(iterations,4))
  param_keeps<-list()

  for(n in 1:max_lesion){
    param_keeps[[n]]<-array(dim=c(iterations,15))
  }
  num_lesion_keeps<-array(dim=c(iterations,1))
  ll_keeps<-array(dim=c(iterations,1))
  prior_keeps<-array(dim=c(iterations,1))
  num_lesion=0
  xlim<-c(-.01,1.01)
  ylim<-c(-.01,1.01)
  P<-matrix(nrow=max_lesion,ncol=max_lesion,data=0)
  diag(P)<-1/3
  d <- row(P) - col(P)
  P[d==-1]<-1/3
  P[d==1]<-1/3
  P[1,1:2]<-.5
  P[max_lesion,(max_lesion-1):max_lesion]<-.5

  
  X <- vector(mode = "list", length = max_lesion)
  if(num_lesion>0){
    for(n in 1:num_lesion){
      X[[n]]<-cbind2(cbind2(1,sin(outer(a2c(x,y,param[n,1:2]),(1:N)))),cos(outer(a2c(x,y,param[n,1:2]),(1:N))))
    } 
  }

  
  ll<-likelihood(param,param0,value,x,y,Sigma,X,N,num_lesion,spatial)
  B<-list()

  start.time <- Sys.time()
  for(i in 1:iterations){
    
    ProgressBar(i+1,iterations)
    
    for(n in which(complete.cases(param))){ 
      XB<-X[[n]]%*%param[n,5:(5+2*N)]
      clust<-+(d2c(x,y,param[n,1:2])<=XB)
      center_star=c(mean(x[clust==1]),mean(y[clust==1]))
      if((abs(center_star[1]-param[n,1])>.01|abs(center_star[2]-param[n,2])>.01)){
        bp<-cbind(XB*sin(a2c(x,y,param[n,1:2]))+param[n,1],XB*cos(a2c(x,y,param[n,1:2]))+param[n,2])
        theta_star<-a2c(bp[,1],bp[,2],center_star)
        f_star<-d2c(bp[,1],bp[,2],center_star)
        X_star<-cbind2(cbind2(1,sin(outer(theta_star,(1:N)))),cos(outer(theta_star,(1:N))))
        beta_star<-try.solve(solve(t(X_star)%*%X_star)%*%t(X_star)%*%f_star)
        if(any(!is.na(beta_star))){
          if(beta_star[1]<0){
            beta_star[1]<-0
          }
          param_star<-param
          param_star[n,]<-c(center_star,param[n,3:4],beta_star)
          X_star<-X
          X_star[[n]]<-cbind2(cbind2(1,sin(outer(a2c(x,y,param_star[n,1:2]),(1:N)))),cos(outer(a2c(x,y,param_star[n,1:2]),(1:N))))
          prior<-prior_param(param_star,param0,x,y,X_star,xlim,ylim,N,num_lesion)
          if(prior>-Inf){
            ll_star<-likelihood( param_star,param0,value,x,y,Sigma,X_star,N,num_lesion,spatial)
            if(sum(ll_star)>-Inf){
              param<-param_star
              ll<-ll_star
              X<-X_star
            }
          }
        }
      }
    }
    

    if(num_lesion==1&num_lesion>min_lesion){
      step=sample(0:3,1,prob=c(rep(1/4,4)))
      if(step==0){ #birth
        birth<-birth_func(param0,param,num_lesion,value,x,y,ll,X,Sigma,xlim,ylim,N,P,spatial)
        param<-birth[[1]]
        num_lesion<-birth[[2]]
        ll<-birth[[3]]
        X<-birth[[4]]

        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]
        param_move<-param_update(param,param0,la,gam,i,value,x,y,ll,Sigma,X,xlim,ylim,N,num_lesion,spatial)
        param<-param_move[[1]]
        ll<-param_move[[2]]
      }
      if(step==1){ #wiggle
        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]
        param_move<-param_update(param,param0,la,gam,i,value,x,y,ll,Sigma,X,xlim,ylim,N,num_lesion,spatial)
        param<-param_move[[1]]
        ll<-param_move[[2]]
      }
      if(step==2){#split
        split<-split_function(param,param0,num_lesion,value,x,y,ll,X,Sigma,xlim,ylim,N,P,spatial)
        param<-split[[1]]
        num_lesion<-split[[2]]
        ll<-split[[3]]
        X<-split[[4]]
  
        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]
        param_move<-param_update(param,param0,la,gam,i,value,x,y,ll,Sigma,X,xlim,ylim,N,num_lesion,spatial)
        param<-param_move[[1]]
        ll<-param_move[[2]]
      }
      if(step==3){ #death
        death<-death_func(param,param0,num_lesion,value,x,y,ll,X,Sigma,xlim,ylim,N,P,spatial)
        param<-death[[1]]
        num_lesion<-death[[2]]
        ll<-death[[3]]
        X<-death[[4]]
    
        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]
      }
      
    }else if(num_lesion==1&num_lesion==min_lesion){
      step=sample(0:2,1,prob=c(rep(1/3,3)))
      if(step==0){ #birth
        birth<-birth_func(param0,param,num_lesion,value,x,y,ll,X,Sigma,xlim,ylim,N,P,spatial)
        param<-birth[[1]]
        num_lesion<-birth[[2]]
        ll<-birth[[3]]
        X<-birth[[4]]
       
        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]
        param_move<-param_update(param,param0,la,gam,i,value,x,y,ll,Sigma,X,xlim,ylim,N,num_lesion,spatial)
        param<-param_move[[1]]
        ll<-param_move[[2]]
      }
      if(step==1){ #wiggle
        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]
        param_move<-param_update(param,param0,la,gam,i,value,x,y,ll,Sigma,X,xlim,ylim,N,num_lesion,spatial)
        param<-param_move[[1]]
        ll<-param_move[[2]]
      }
      if(step==2){#split
        split<-split_function(param,param0,num_lesion,value,x,y,ll,X,Sigma,xlim,ylim,N,P,spatial)
        param<-split[[1]]
        num_lesion<-split[[2]]
        ll<-split[[3]]
        X<-split[[4]]
       
        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]
        param_move<-param_update(param,param0,la,gam,i,value,x,y,ll,Sigma,X,xlim,ylim,N,num_lesion,spatial)
        param<-param_move[[1]]
        ll<-param_move[[2]]
      }
    }else if(1<num_lesion&num_lesion<max_lesion&num_lesion>min_lesion){
      step=sample(0:4, 1,prob=c(rep(1/5,5)))
      
      if(step==0){ #birth
        birth<-birth_func(param0,param,num_lesion,value,x,y,ll,X,Sigma,xlim,ylim,N,P,spatial)
        param<-birth[[1]]
        num_lesion<-birth[[2]]
        ll<-birth[[3]]
        X<-birth[[4]]
        
        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]
        param_move<-param_update(param,param0,la,gam,i,value,x,y,ll,Sigma,X,xlim,ylim,N,num_lesion,spatial)
        param<-param_move[[1]]
        ll<-param_move[[2]]
        
      }
      if(step==1){ #wiggle
        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]
        param_move<-param_update(param,param0,la,gam,i,value,x,y,ll,Sigma,X,xlim,ylim,N,num_lesion,spatial)
        param<-param_move[[1]]
        ll<-param_move[[2]]
      }
      if(step==2){ #death
        death<-death_func(param,param0,num_lesion,value,x,y,ll,X,Sigma,xlim,ylim,N,P,spatial)
        param<-death[[1]]
        num_lesion<-death[[2]]
        ll<-death[[3]]
        X<-death[[4]]

        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]
        param_move<-param_update(param,param0,la,gam,i,value,x,y,ll,Sigma,X,xlim,ylim,N,num_lesion,spatial)
        param<-param_move[[1]]
        ll<-param_move[[2]]
      }
      if(step==3){ #split
        split<-split_function(param,param0,num_lesion,value,x,y,ll,X,Sigma,xlim,ylim,N,P,spatial)
        param<-split[[1]]
        num_lesion<-split[[2]]
        ll<-split[[3]]
        X<-split[[4]]
    
        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]
        param_move<-param_update(param,param0,la,gam,i,value,x,y,ll,Sigma,X,xlim,ylim,N,num_lesion,spatial)
        param<-param_move[[1]]
        ll<-param_move[[2]]
      }
      if(step==4){#merge
        merge<-merge_function(param,param0,num_lesion,value,x,y,ll,X,Sigma,xlim,ylim,N,P,spatial)
        param<-merge[[1]]
        num_lesion<-merge[[2]]
        ll<-merge[[3]]
        X<-merge[[4]]
       
        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]
        param_move<-param_update(param,param0,la,gam,i,value,x,y,ll,Sigma,X,xlim,ylim,N,num_lesion,spatial)
        param<-param_move[[1]]
        ll<-param_move[[2]]
      }
    }else if(num_lesion==max_lesion){
      step=sample(0:2, 1,prob=c(1/3,1/3,1/3))

      if(step==0){ #death
        death<-death_func(param,param0,num_lesion,value,x,y,ll,X,Sigma,xlim,ylim,N,P,spatial)
        param<-death[[1]]
        num_lesion<-death[[2]]
        ll<-death[[3]]
        X<-death[[4]]
   
        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]
        param_move<-param_update(param,param0,la,gam,i,value,x,y,ll,Sigma,X,xlim,ylim,N,num_lesion,spatial)
        param<-param_move[[1]]
        ll<-param_move[[2]]
      }
      if(step==1){ #wiggle
        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]
        param_move<-param_update(param,param0,la,gam,i,value,x,y,ll,Sigma,X,xlim,ylim,N,num_lesion,spatial)
        param<-param_move[[1]]
        ll<-param_move[[2]]
      }
      if(step==2){#merge
        merge<-merge_function(param,param0,num_lesion,value,x,y,ll,X,Sigma,xlim,ylim,N,P,spatial)
        param<-merge[[1]]
        num_lesion<-merge[[2]]
        ll<-merge[[3]]
        X<-merge[[4]]
      
        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]
        param_move<-param_update(param,param0,la,gam,i,value,x,y,ll,Sigma,X,xlim,ylim,N,num_lesion,spatial)
        param<-param_move[[1]]
        ll<-param_move[[2]]
      }
    }else if(num_lesion==min_lesion&min_lesion==0){
      step=sample(0:1,1,prob=c(rep(1/2,2)))
      if(step==0){ #birth
        birth<-birth_func(param0,param,num_lesion,value,x,y,ll,X,Sigma,xlim,ylim,N,P,spatial)
        param<-birth[[1]]
        num_lesion<-birth[[2]]
        ll<-birth[[3]]
        X<-birth[[4]]
  
        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]
      
      }
      if(step==1){ #wiggle
        param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,N,num_lesion,spatial)
        param0<-param0_move[[1]]
        ll<-param0_move[[2]]

      }
    }


    if(i>=50){
      mu0<-mu0+(1/(i+1))*(param0-mu0)
      gam0<-gam0+(1/(i+1))*((param0-mu0)%*%t(param0-mu0)-gam0)
      for(n in which(complete.cases(param))){
        mu[n,]<-mu[n,]+(1/(i+1))*(param[n,]-mu[n,])
        gam[[n]][[1]]<-gam[[n]][[1]]+(1/(i+1))*((param[n,3:4]-mu[n,3:4])%*%t(param[n,3:4]-mu[n,3:4])-gam[[n]][[1]])
        gam[[n]][[2]]<-gam[[n]][[2]]+(1/(i+1))*((param[n,5:10]-mu[n,5:10])%*%t(param[n,5:10]-mu[n,5:10])-gam[[n]][[2]])
        gam[[n]][[3]]<-gam[[n]][[3]]+(1/(i+1))*((param[n,11:15]-mu[n,11:15])%*%t(param[n,11:15]-mu[n,11:15])-gam[[n]][[3]])
      }
    }
  
    #order by center
    order<-order(param[,1],param[,2])
    param<-param[order,]  
    X<-X[order]
    mu<-mu[order,]  
    gam<-gam[order]  
    la<-la[order]  
   

    param0_keeps[i,]<-param0
    num_lesion_keeps[i,]<-num_lesion
    ll_keeps[i,]<-sum(ll)
    prior_keeps[i,]<-prior_param(param,param0,x,y,X,xlim,ylim,N,num_lesion)
 
    for(n in which(complete.cases(param))){
      param_keeps[[n]][i,]<-param[n,]
      
    }


    if(i%%20==0 & i>=100){
      last20<-param0_keeps[(i-20):i,]
      accept<-1-mean(duplicated(last20))
      if(accept<.4) la0<-la0-.1
      if(accept>.4) la0<-la0+.1
      for(n in which(complete.cases(param))){
        last20<-param_keeps[[n]][(i-20):i,3:4]
        accept<-1-mean(duplicated(round(last20,8)))
        if(accept<.2) la[[n]][1]<-la[[n]][1]-.1
        if(accept>.2) la[[n]][1]<-la[[n]][1]+.1
        
        last20<-param_keeps[[n]][(i-20):i,5:10]
        accept<-1-mean(duplicated(round(last20,8)))
        if(accept<.2) la[[n]][2]<-la[[n]][2]-.1
        if(accept>.2) la[[n]][2]<-la[[n]][2]+.1
        
        last20<-param_keeps[[n]][(i-20):i,11:15]
        accept<-1-mean(duplicated(round(last20,8)))
        if(accept<.2) la[[n]][3]<-la[[n]][3]-.1
        if(accept>.2) la[[n]][3]<-la[[n]][3]+.1
      }
    }
    
  }
   for(n in 1:max_lesion){
    param_keeps[[n]]<-param_keeps[[n]][-c(1:burn),]
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  value<-value*attr(value, 'scaled:scale') + attr(value, 'scaled:center')
  list1<-list(param_keeps,param0_keeps[-c(1:burn),],x,y,value,num_lesion_keeps[-c(1:burn),],ll_keeps[-c(1:burn),],prior_keeps[-c(1:burn),],time.taken)
  
  list1
}

post_process_bfspm<-function(results,prob=.5){
  rez<-result
  num_lesion<-rez[[6]]
  param_keeps<-rez[[1]]
  x<-rez[[3]]
  y<-rez[[4]]
  value<-rez[[5]]
  centroids<-lapply(param_keeps, "[", , 1:2)
  boundary_params<-lapply(param_keeps, "[", , c(1:2,5:15))

  centroids<-do.call(rbind, centroids)
  boundary_params<-do.call(rbind, boundary_params)

  centroids<-centroids[complete.cases(centroids),]
  boundary_params<-boundary_params[complete.cases(boundary_params),]

  centroids<-round(centroids,1)
  centroids<-paste0(centroids[,1],",",centroids[,2])
  centriods<-as.data.frame(centroids)
  
  table_centroids<-table(centroids)
  table_centroids<-table(centroids)[table(centroids)/dim(param_keeps[[1]])[1]>prob]
  param_avg<-array(dim=c(length(table_centroids),13))
  cluster<-rep(0,length(x))
  cluster_uncert<-rep(0,length(x))
  if(length(table_centroids)>0){
    for(i in 1:length(table_centroids)){
      centroid<-names(table_centroids)[i]
      param_avg[i,]<-colMeans(boundary_params[paste0(round(boundary_params[,1],1),",",round(boundary_params[,2],1))==centroid,])
      
    }
    for(i in 1:length(table_centroids)){
      X<-cbind2(cbind2(1,sin(outer(a2c(x,y,param_avg[i,1:2]),(1:N)))),cos(outer(a2c(x,y,param_avg[i,1:2]),(1:N))))
      beta<-param_avg[i,3:13]
      clust<-+(d2c(x,y,param_avg[i,1:2])<=X%*%beta)
      
      clust_uncert<-clust*table_centroids[i]/dim(param_keeps[[1]])[1]
      
      cluster<-cluster+i*clust
      cluster_uncert<-cluster_uncert+clust_uncert
      
    }
  }

  cluster_uncert<-replace(cluster_uncert, cluster_uncert>1, 1)
  cluster<-factor(cluster,levels=unique(cluster)[order(unique(cluster))],labels=0:(length(unique(cluster))-1))
  data<-cbind.data.frame(x,y,value,cluster,cluster_uncert)
  data<-data%>%group_by(cluster)%>%mutate(average_voxel_value=mean(value))
  data
}

