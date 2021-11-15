# Fitting heteroscedastic SEM models
#

hetero_sem=function(formulamean,formulavar,data,W,nsim,burn,step,prior,initial,kernel="normal",seed=0){
  ## Mean model ##
  y_n_mean <- as.character(formulamean[[2]])
  X0_mean <- as.character(formulamean[[3]])[-1]
  X1_mean <- as.character(do.call("c",sapply(X0_mean, function(x){strsplit(x,"\\+")})))
  X_n_mean <- gsub(" ","",X1_mean)

  y_mean <- data[,which(names(data)==y_n_mean)]
  X_mean <- as.matrix(data[,which(names(data)%in%X_n_mean)])

  ## Variance model ##

  X0_var <- as.character(formulavar)[[2]]
  X1_var <- as.character(do.call("c",sapply(X0_var, function(x){strsplit(x,"\\+")})))
  X_n_var <- gsub(" ","",X1_var)

  X_var <- as.matrix(data[,which(names(data)%in%X_n_var)])

  ## Prior information ##
  b_pri <- prior$b_pri
  B_pri <- prior$B_pri
  g_pri <- prior$g_pri
  G_pri <- prior$G_pri

  ## Initial values ##
  beta_0 <- initial$beta_0
  gammas_0 <- initial$gamma_0
  lambda_0 <- initial$lambda_0

  output <- hetero_sem_int(y_mean,X_mean,X_var,W,nsim,burn,step,b_pri,B_pri,g_pri,G_pri,beta_0,gammas_0,lambda_0,kernel="normal",seed=0)

  return(output)
}

#' Hello
#'
#' @keywords internal
#'
hetero_sem_int=function(y,X,Z,W,nsim,burn,step,b_pri,B_pri,g_pri,G_pri,beta_0,gammas_0,lambda_0,kernel="normal",seed=0)
{
  set.seed(seed)
  ########Lectura de la informaci?n
  rowst=function(x){
    x1=c()
    x1=(x)/sum(x)
  }
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
  dpost <- function(betas,gammas,lambda) {
    A=diag(nrow(X))-lambda*matstand
    Sigma=diag(c(exp((Z)%*%gammas)))
    solvSigma=diag(1/c(exp(Z%*%gammas)))
    k=t(A%*%(y-X%*%(betas)))%*%solvSigma%*%(A%*%(y-X%*%(betas)))
    fc.y=k
    fc.beta=t(betas - b_pri)%*%solve(B_pri)%*%(betas-b_pri)
    fc.gamma=t(gammas - g_pri)%*%solve(G_pri)%*%(gammas-g_pri)
    # dp <- det(Sigma)^(-1/2)*det(A)*exp(-0.5*(fc.y + fc.beta+fc.gamma))
    # dp
    logdp <- (-1/2)*log(det(Sigma)) + log(det(A)) -0.5*fc.y - 0.5*fc.beta - 0.5*fc.gamma
    return(logdp)
  }

  #Generacion de valores para las distribuciones propuestas
  r.proposal_gamma=function(Gammas){
    a.now=Z%*%Gammas
    A=diag(nrow(X))-Lambda*matstand
    b.now=A%*%(y-X%*%betas.now)
    y.now=a.now+(b.now^2/exp(a.now))-1
    G_pos=solve(solve(G_pri)+0.5*t(Z)%*%Z)
    g_pos=G_pos%*%(solve(G_pri)%*%g_pri+0.5*(t(Z)%*%y.now))
    gammas.pro=rmvnorm(1,g_pos,G_pos)
    gammas.pro
  }

  dproposal_gamma<-function(gammas.now, gammas.old){
    a.now=Z%*%gammas.old
    A=diag(nrow(X))-Lambda*matstand
    b.now=A%*%(y-X%*%betas.now)
    y.now=a.now+(b.now^2/exp(a.now))-1
    G_pos=solve(solve(G_pri)+0.5*t(Z)%*%Z)
    g_pos=G_pos%*%(solve(G_pri)%*%g_pri+0.5*(t(Z)%*%y.now))
    dmvnorm(gammas.now,g_pos,G_pos,log = TRUE)
  }

  dproposal_lambda<-function(lambda){
    Sigma=diag(c(exp(Z%*%Gammas)))
    solvSigma=diag(1/c(exp(Z%*%Gammas)))
    a=t(y-X%*%betas.now)%*%t(matstand)%*%solvSigma%*%matstand%*%(y-X%*%betas.now)
    b=t(y-X%*%betas.now)%*%t(matstand)%*%solvSigma%*%(y-X%*%betas.now)
    dnorm(lambda,b/a,1/sqrt(a),log = TRUE)
  }

  #Algoritmo Metropolis Hastings
  beta.mcmc=matrix(NA,nrow=nsim,ncol(X))
  gamma.mcmc=matrix(NA,nrow=nsim,ncol(Z))
  lambda.mcmc=c()
  ind1=rep(0,nsim)
  ind2=rep(0,nsim)
  logV_DIC=c()
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  if(kernel=="uniform"){
    for(i in 1:nsim){
      #Valores a posteriori condicional
      if(i==1){
        Lambda=lambda_0
        Gammas=gammas_0
        Sigma=diag(c(exp(Z%*%Gammas)))
      }
      else{
        Sigma=diag(c(exp(Z%*%Gammas)))
      }
      A=diag(nrow(X))-Lambda*matstand
      solvSigma=diag(1/c(exp(Z%*%Gammas)))
      B_pos=solve(solve(B_pri)+t(A%*%X)%*%solvSigma%*%A%*%X)
      b_pos=B_pos%*%(solve(B_pri)%*%b_pri+t(A%*%X)%*%solvSigma%*%A%*%y)
      #Beta a posteriori condicional
      betas.now=c(rmvnorm(1,b_pos,B_pos))
      #A posteriori condicional completa para Sigma2
      gammas.now=c(r.proposal_gamma(Gammas))
      q1.1=dproposal_gamma(gammas.now,Gammas)
      q2.1=dproposal_gamma(Gammas,gammas.now)
      p1.1=dpost(betas.now,gammas.now,Lambda)
      p2.1=dpost(betas.now,Gammas,Lambda)
      T.val=min(1,(p1.1/p2.1)*(q1.1/q2.1))
      u<-runif(1)
      if(p2.1==0){T.val=0}
      if(q2.1==0){T.val=0}
      if (u <=T.val) {
        Gammas= gammas.now
        ind1[i] = 1
      }
      #A posteriori condicional completa para Lambda
      lambda.now=runif(1,1/abs(min(eigen(mat)$values)),1)
      p1.2=dpost(betas.now,Gammas,lambda.now)
      p2.2=dpost(betas.now,Gammas,Lambda)
      T.val2=min(1,p1.2/p2.2)
      u<-runif(1)
      if(p2.2==0){T.val2=0}
      if (u <=T.val2) {
        Lambda <- lambda.now
        ind2[i] = 1
      }
      beta.mcmc[i,]<-betas.now
      gamma.mcmc[i,]<-gammas.now
      lambda.mcmc[i]<-lambda.now
      Sigma=diag(c(exp(Z%*%gamma.mcmc[i,])))
      detS=det(Sigma)
      detB=det(diag(nrow(X))-lambda.mcmc[i]*matstand)
      Yg=(diag(nrow(X))-lambda.mcmc[i]*matstand)%*%(y-X%*%beta.mcmc[i,])
      logV_DIC[i]=(-(nrow(X)/2)*log(pi))+log(detB)-0.5*log(detS)-0.5*t(Yg)%*%diag(1/c(exp(Z%*%gamma.mcmc[i,])))%*%Yg
      Sys.sleep(0.000000001)
      # update progress bar
      setTxtProgressBar(pb, i)
    }
  }
  if(kernel=="normal"){
    for(i in 1:nsim){
      #Valores a posteriori condicional
      if(i==1){
        Gammas=gammas_0
        Sigma=diag(c(exp(Z%*%Gammas)))
        Lambda=lambda_0
      }
      else{
        Sigma=diag(c(exp(Z%*%Gammas)))
      }
      A=diag(nrow(X))-Lambda*matstand
      solvSigma=diag(1/c(exp(Z%*%Gammas)))
      B_pos=solve(solve(B_pri)+t(A%*%X)%*%solvSigma%*%A%*%X)
      b_pos=B_pos%*%(solve(B_pri)%*%b_pri+t(A%*%X)%*%solvSigma%*%A%*%y)
      betas.now=c(rmvnorm(1,b_pos,B_pos))
      ###Propuesta de gammas
      gammas.now=c(r.proposal_gamma(Gammas))
      q1.1=dproposal_gamma(gammas.now,Gammas)
      q2.1=dproposal_gamma(Gammas,gammas.now)
      p1.1=dpost(betas.now,gammas.now,Lambda)
      p2.1=dpost(betas.now,Gammas,Lambda)
      met.a1 <- ifelse(p1.1>p2.1,log(p1.1-p2.1),-log(p2.1-p1.1))
      met.b1 <- ifelse(q1.1>q2.1,log(q1.1-q2.1),-log(q2.1-q1.1))
      T.val1=min(0,met.a1+met.b1)
      #T.val=min(1,(p1.1/p2.1)*(q1.1/q2.1))
      u<-runif(1)
      # if(p2.1==0){T.val=0}
      # if(q2.1==0){T.val=0}
      if (u <=exp(T.val1)) {
        Gammas= gammas.now
        ind1[i] = 1
      }
      ###Propuesta de Lambda
      Sigma=diag(c(exp(Z%*%Gammas)))
      solvSigma=diag(1/c(exp(Z%*%Gammas)))
      a=t(y-X%*%betas.now)%*%t(matstand)%*%solvSigma%*%matstand%*%(y-X%*%betas.now)
      b=t(y-X%*%betas.now)%*%t(matstand)%*%solvSigma%*%(y-X%*%betas.now)
      eigenvals <- eigen(matstand)$values
      lowlim <- -1/(max(abs(eigenvals[eigenvals<0])))
      lambda.now=rnorm(1,b/a,1/sqrt(a))
      while(lambda.now>1 || lambda.now< lowlim){
        lambda.now <- rnorm(1,b/a,1/sqrt(a))
      }
      p1.2=dpost(betas.now,Gammas,lambda.now)
      p2.2=dpost(betas.now,Gammas,Lambda)
      q1.2=dproposal_lambda(lambda.now)
      q2.2=dproposal_lambda(Lambda)
      met.a2 <- ifelse(p1.2>p2.2,log(p1.2-p2.2),-log(p2.2-p1.2))
      met.b2 <- ifelse(q1.2>q2.2,log(q1.2-q2.2),-log(q2.2-q1.2))
      T.val2=min(0,met.a2+met.b2)
      #T.val2=min(1,p1.2/p2.2)
      u<-runif(1)
      #if(p2.2==0){T.val2=0}
      if (u <=exp(T.val2)) {
        Lambda <- lambda.now
        ind2[i] = 1
      }
      beta.mcmc[i,]<-betas.now
      gamma.mcmc[i,]<-gammas.now
      lambda.mcmc[i]<-lambda.now
      Sigma=diag(c(exp(Z%*%gamma.mcmc[i,])))
      detS=det(Sigma)
      detB=det(diag(nrow(X))-lambda.mcmc[i]*matstand)
      Yg=(diag(nrow(X))-lambda.mcmc[i]*matstand)%*%(y-X%*%beta.mcmc[i,])
      logV_DIC[i]=(-(nrow(X)/2)*log(pi))+log(detB)-0.5*log(detS)-0.5*t(Yg)%*%diag(1/c(exp(Z%*%gamma.mcmc[i,])))%*%Yg
      Sys.sleep(0.000000001)
      # update progress bar
      setTxtProgressBar(pb, i)
    }
  }
  beta.mcmc_1=beta.mcmc[(burn+1):nsim,]
  gamma.mcmc_1=gamma.mcmc[(burn+1):nsim,]
  lambda.mcmc_1=lambda.mcmc[(burn+1):nsim]
  beta.mcmc_2=matrix(NA,nrow=(nsim-burn+1)/step,ncol(X))
  gamma.mcmc_2=matrix(NA,nrow=(nsim-burn+1)/step,ncol(Z))
  lambda.mcmc_2=c()
  for (i in 1:(nsim-burn+1))
  {
    if(i%%step==0)
    {
      beta.mcmc_2[i/step,]=beta.mcmc_1[i,]
      gamma.mcmc_2[i/step,]=gamma.mcmc_1[i,]
      lambda.mcmc_2[i/step]=lambda.mcmc_1[i]
    }
  }

  Bestimado = colMeans(beta.mcmc_2)
  Gammaest = colMeans(gamma.mcmc_2)
  lambda.mcmc_3=lambda.mcmc_2[lambda.mcmc_2<=1]
  Lambdaest=mean(lambda.mcmc_3)
  DesvBeta <- apply(beta.mcmc_2,2,sd)
  DesvGamma <- apply(gamma.mcmc_2,2,sd)
  DesvLambda<-sd(lambda.mcmc_3)
  Betaquant <- t(apply(beta.mcmc_2,2,function(x){quantile(x,c(0.025,0.5,0.975))}))
  Gammaquant <-  t(apply(gamma.mcmc_2,2,function(x){quantile(x,c(0.025,0.5,0.975))}))
  Lambdaquant <- quantile(lambda.mcmc_3,c(0.025,0.5,0.975))
  AccRate1<-sum(ind1)/nsim
  AccRate2<-sum(ind2)/nsim
  Sigma1=diag(c(exp(Z%*%Gammaest)))
  detS=det(Sigma1)
  detA=det(diag(nrow(X))-Lambdaest*matstand)
  Yg=(diag(nrow(X))-Lambdaest*matstand)%*%(y-X%*%Bestimado)
  Veros=detA*((detS)^(-0.5))*exp(-0.5*t(Yg)%*%solve(Sigma1)%*%Yg)
  p=ncol(X)+ncol(Z)+1
  BIC=-2*log(Veros)+p*log(nrow(X))
  logV_DIC=logV_DIC[is.nan(logV_DIC)==FALSE]
  Dbar=mean(-2*logV_DIC)
  logV1_DIC=(-(nrow(X)/2)*log(pi))+log(detA)-0.5*log(detS)-0.5*t(Yg)%*%solve(Sigma1)%*%Yg
  Dev=-2*logV1_DIC
  DIC=2*Dbar+Dev

  summary = data.frame( mean=c(Bestimado,Gammaest,Lambdaest),
                        sd = c(DesvBeta,DesvGamma,DesvLambda),
                        q0.025=c(Betaquant[,1],Gammaquant[,1],Lambdaquant[1]),
                        q0.5=c(Betaquant[,2],Gammaquant[,2],Lambdaquant[2]),
                        q0.975=c(Betaquant[,3],Gammaquant[,3],Lambdaquant[3]))

  rownames(summary) = c("x0","x1","x2","z0","z1","z2","lambda")
  #rownames(summary) = c("x0","x1","x2","z0","z1","lambda")

  return(list(summary=summary,Acceptance_Rates=list(Gamma_AccRate=AccRate1,Lambda_AccRate=AccRate2),Criteria=list(BIC=BIC,DIC=DIC),chains=mcmc(data.frame(beta_chain=beta.mcmc,gamma_chain=gamma.mcmc,lambda_chain=lambda.mcmc),thin = 1)))
}
