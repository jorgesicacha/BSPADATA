# Fitting homoscedastic SEM models

hom_sem=function(formula,data,W,nsim,burn,step,prior,initial,kernel="normal",seed=0){

  y_n <- as.character(formula[[2]])
  X0 <- as.character(formula[[3]])[-1]
  X1 <- as.character(do.call("c",sapply(X0, function(x){strsplit(x,"\\+")})))
  X_n <- gsub(" ","",X1)

  y <- data[,which(names(data)==y_n)]
  X <- as.matrix(data[,which(names(data)%in%X_n)])

  b_pri <- prior$b_pri
  B_pri <- prior$B_pri
  r_pri <- prior$r_pri
  lambda_pri <- prior$lambda_pri

  beta_0 <- initial$beta_0
  sigma2_0 <- initial$sigma2_0
  lambda_0 <- initial$lambda_0

  output <- hom_sem_int(y,X,W,nsim,burn,step,b_pri,B_pri,r_pri,lambda_pri,beta_0,sigma2_0,lambda_0,kernel="normal",seed=seed)

  return(output)
}

#' Hello
#'
#' @keywords internal
#'
hom_sem_int=function(y,X,W,nsim,burn,step,b_pri,B_pri,r_pri,lambda_pri,beta_0,sigma2_0,lambda_0,kernel="normal",seed=0)
{
  set.seed(seed)
  rowst=function(x){
    x1=c()
    x1=(x)/sum(x)
  }
  ########Lectura de la informaci?n
  y=as.matrix(y)
  if (is.null(X) | is.null(y) ){
    stop("No data")
  }
  if(burn>nsim | burn<0){
    stop("Burn must be between 0 and nsim")
  }
  if(nsim<=0){
    stop("There must be more than 0 simulations")
  }
  if(step<0 | step > nsim){
    stop("Jump length must not be lesser than 0 or greater than nsim")
  }
  if(class(W)=="nb"){
    matstand=nb2mat(W)
    mat0=nb2listw(W,style="B")
    mat=listw2mat(mat0)
  }
  else{
    if(class(W)=="listw"){
      mat=listw2mat(W)
      matstand=apply(mat,2,rowst)
      matstand=t(matstand)
    }
    else{
      if(sum(rowSums(W))==nrow(X))
      {
        matstand=W
        mat=matrix(nrow=nrow(X),ncol=nrow(X))
        for(i in 1:nrow(mat)){
          for(j in 1:ncol(mat)){
            if(matstand[i,j]==0){mat[i,j]=0}
            else{mat[i,j]=1/matstand[i,j]}
          }
        }
      }
      else{
        mat=W
        matstand=apply(mat,2,rowst)
        matstand=t(matstand)
      }
    }
  }
  dpost <- function(betas,sigma2,lambda) {
    A=diag(nrow(X))-lambda*matstand
    k=t(A%*%(y-X%*%betas))%*%(A%*%(y-X%*%betas))
    fc.y=k
    fc.beta=t(b_pri - betas)%*%solve(B_pri)%*%(b_pri-betas)
    fc.sigma2=(lambda_pri^(r_pri))*(sigma2)^(-r_pri-1)*exp(-lambda_pri/sigma2)/gamma(r_pri)
    logdp <- (-nrow(X)/2)*log(sigma2) + log(det(A)) -0.5*fc.y/sigma2 - 0.5*fc.beta + fc.sigma2
    return(logdp)
  }
  dproposal <- function(lambda) {
    a=(1/sigma2.now)*t(y-X%*%betas.now)%*%t(matstand)%*%matstand%*%(y-X%*%betas.now)
    b=(1/sigma2.now)*t(y-X%*%betas.now)%*%(matstand)%*%(y-X%*%betas.now)
    dmvnorm(lambda,b/a,1/sqrt(a),log = TRUE)
  }

  ind=rep(0,nsim)
  beta.mcmc=matrix(NA,nrow=nsim,ncol=ncol(X))
  sigma2.mcmc=c()
  lambda.mcmc=c()
  logV_DIC=c()
  Sigma_0=(sigma2_0)*diag(nrow(X))
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  if(kernel=="uniform"){
    for(i in 1:nsim){
      if(i==1){
        Sigma=Sigma_0
        Lambda=lambda_0
      }
      else{
        Sigma=sigma2.now*diag(nrow(X))
      }
      A=diag(nrow(X))-Lambda*matstand
      B_pos=solve(solve(B_pri)+t(X)%*%t(A)%*%solve(Sigma)%*%A%*%X)
      b_pos=B_pos%*%(solve(B_pri)%*%b_pri+t(X)%*%t(A)%*%solve(Sigma)%*%A%*%y)
      #Beta a posteriori condicional
      betas.now=c(rmvnorm(1,b_pos,B_pos))
      #A posteriori condicional completa para gammas
      r_pos=nrow(X)/2+r_pri
      A=diag(nrow(X))-Lambda*matstand
      k=t(A%*%(y-X%*%(betas.now)))%*%(A%*%(y-X%*%(betas.now)))
      lambda_pos=(k+2*lambda_pri)/2
      sigma2.now=rigamma(1,r_pos, lambda_pos)
      #A posteriori condicional completa para Lambda
      lambda.now=runif(1,1/abs(min(eigen(mat)$values)),1)
      p1=dpost(betas.now,sigma2.now,lambda.now)
      p2=dpost(betas.now,sigma2.now,Lambda)
      T.val=min(1,p1/p2)
      u<-runif(1)
      if (u <=T.val) {
        Lambda <- lambda.now
        ind[i] = 1
      }
      beta.mcmc[i,]<-betas.now
      sigma2.mcmc[i]<-sigma2.now
      lambda.mcmc[i]<-lambda.now
      Sigma=diag(sigma2.mcmc[i],nrow(X))
      detS=det(Sigma)
      detB=det(diag(nrow(X))-lambda.mcmc[i]*matstand)
      Yg=(diag(nrow(X))-lambda.mcmc[i]*matstand)%*%(y-X%*%beta.mcmc[i,])
      logV_DIC[i]=(-(nrow(X)/2)*log(pi))+log(detB)-0.5*log(detS)-0.5*t(Yg)%*%solve(Sigma)%*%Yg
      Sys.sleep(0.000000001)
      # update progress bar
      setTxtProgressBar(pb, i)
    }
  }
  if(kernel=="normal"){
    for(i in 1:nsim){
      #A posteriori condicional completa para Betas
      if(i==1){
        Sigma=Sigma_0
        Lambda=lambda_0
      }
      else{
        Sigma=sigma2.now*diag(nrow(X))
      }

      A=diag(nrow(X))-Lambda*matstand
      B_pos=solve(solve(B_pri)+t(X)%*%t(A)%*%diag(1/diag(Sigma))%*%A%*%X)
      b_pos=B_pos%*%(solve(B_pri)%*%b_pri+t(X)%*%t(A)%*%diag(1/diag(Sigma))%*%A%*%y)
      #Beta a posteriori condicional
      betas.now=c(rmvnorm(1,b_pos,B_pos))

      #A posteriori condicional completa para gammas
      r_pos=nrow(X)/2+r_pri
      A=diag(nrow(X))-Lambda*matstand
      k=t(A%*%(y-X%*%(betas.now)))%*%(A%*%(y-X%*%(betas.now)))
      lambda_pos=(k+2*lambda_pri)/2
      sigma2.now=rigamma(1,r_pos, lambda_pos)

      eigenvals <- eigen(matstand)$values
      lowlim <- -1/(max(abs(eigenvals[eigenvals<0])))

      a=(1/sigma2.now)*t(y-X%*%betas.now)%*%t(matstand)%*%matstand%*%(y-X%*%betas.now)
      b=(1/sigma2.now)*t(y-X%*%betas.now)%*%(matstand)%*%(y-X%*%betas.now)
      lambda.now=rnorm(1,b/a,1/sqrt(a))
      while(lambda.now>1 || lambda.now< lowlim){
        lambda.now <- rnorm(1,b/a,1/sqrt(a))
      }

      p1=dpost(betas.now,sigma2.now,lambda.now)
      p2=dpost(betas.now,sigma2.now,Lambda)
      q1=dproposal(lambda.now)
      q2=dproposal(Lambda)
      met.a <- ifelse(p1>p2,log(p1-p2),-log(p2-p1))
      met.b <- ifelse(q1>q2,log(q1-q2),-log(q2-q1))
      T.val=min(0,met.a+met.b)
      u<-runif(1)
      if (u <=exp(T.val)) {
        Lambda <- lambda.now
        ind[i] = 1
      }
      beta.mcmc[i,]<-betas.now
      sigma2.mcmc[i]<-sigma2.now
      lambda.mcmc[i]<-lambda.now
      Sigma=diag(sigma2.mcmc[i],nrow(X))
      detS=det(Sigma)
      detB=det(diag(nrow(X))-lambda.mcmc[i]*matstand)
      Yg=(diag(nrow(X))-lambda.mcmc[i]*matstand)%*%(y-X%*%beta.mcmc[i,])
      logV_DIC[i]=(-(nrow(X)/2)*log(pi))+log(detB)-0.5*log(detS)-0.5*t(Yg)%*%solve(Sigma)%*%Yg
      Sys.sleep(0.000000001)
      # update progress bar
      setTxtProgressBar(pb, i)

    }
  }
  beta.mcmc_1=beta.mcmc[(burn+1):nsim,]
  sigma2.mcmc_1=sigma2.mcmc[(burn+1):nsim]
  lambda.mcmc_1=lambda.mcmc[(burn+1):nsim]
  beta.mcmc_2=matrix(NA,nrow=(nsim-burn+1)/step,3)
  sigma2.mcmc_2=c()
  lambda.mcmc_2=c()
  for (i in 1:(nsim-burn+1))
  {
    if(i%%step==0)
    {
      beta.mcmc_2[i/step,]=beta.mcmc_1[i,]
      sigma2.mcmc_2[i/step]=sigma2.mcmc_1[i]
      lambda.mcmc_2[i/step]=lambda.mcmc_1[i]
    }
  }

  Bestimado = colMeans(beta.mcmc_2)
  Sigma2est = mean(sigma2.mcmc_2)
  lambda.mcmc_3=lambda.mcmc_2[lambda.mcmc_2<=1]
  Lambdaest=mean(lambda.mcmc_3)
  DesvBeta <- apply(beta.mcmc_2,2,sd)
  DesvSigma2 <- sd(sigma2.mcmc_2)
  DesvLambda<-sd(lambda.mcmc_3)
  Betaquant <- t(apply(beta.mcmc_2,2,function(x){quantile(x,c(0.025,0.5,0.975))}))
  Sigma2quant <- quantile(sigma2.mcmc_2,c(0.025,0.5,0.975))
  Lambdaquant <- quantile(lambda.mcmc_3,c(0.025,0.5,0.975))
  AccRate<-sum(ind)/nsim
  Sigma=diag(Sigma2est,nrow(X))
  detS=det(Sigma)
  detB=det(diag(nrow(X))-Lambdaest*matstand)
  Yg=(diag(nrow(X))-Lambdaest*matstand)%*%(y-X%*%Bestimado)
  Veros=detB*((detS)^(-0.5))*exp(-0.5*t(Yg)%*%solve(Sigma)%*%Yg)
  logV=(-(nrow(X)/2)*log(2*pi))+log(detB)-0.5*log(detS)-0.5*t(Yg)%*%solve(Sigma)%*%Yg
  p=ncol(X)+2
  BIC=-2*logV+p*log(nrow(X))
  logV_DIC=logV_DIC[is.nan(logV_DIC)==FALSE]
  Dbar=mean(-2*logV_DIC)
  logV1_DIC=(-(nrow(X)/2)*log(2*pi))+log(detB)-0.5*log(detS)-0.5*t(Yg)%*%solve(Sigma)%*%Yg
  Dev=-2*logV1_DIC
  DIC=2*Dbar+Dev

  summary = data.frame( mean=c(Bestimado,Sigma2est,Lambdaest),
                        sd = c(DesvBeta,DesvSigma2,DesvLambda),
                        q0.025=c(Betaquant[,1],Sigma2quant[1],Lambdaquant[1]),
                        q0.5=c(Betaquant[,2],Sigma2quant[2],Lambdaquant[2]),
                        q0.975=c(Betaquant[,3],Sigma2quant[3],Lambdaquant[3]))

  rownames(summary) = c("x0","x1","x2","sigma2","lambda")

  return(list(summary=summary,Acceptance_Rate=AccRate,Criteria=list(BIC=BIC,DIC=DIC),chains=mcmc(data.frame(beta_chain=beta.mcmc,sigma2_chain=sigma2.mcmc,lambda_chain=lambda.mcmc),thin = 1)))
}

