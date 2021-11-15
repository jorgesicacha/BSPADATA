# Fitting homoscedastic General model

hom_general=function(formula,data,W1,W2=NULL,nsim,burn,step,prior,initial,kernel="normal",mateq=TRUE,impacts=TRUE,seed=0){

  y_n <- as.character(formula[[2]])
  X0 <- as.character(formula[[3]])[-1]
  X1 <- as.character(do.call("c",sapply(X0, function(x){strsplit(x,"\\+")})))
  X_n <- gsub(" ","",X1)

  y <- data[,which(names(data)==y_n)]
  X <- as.matrix(data[,which(names(data)%in%X_n)])

  b_pri <- prior$b_pri
  B_pri <- prior$B_pri
  r_pri <- prior$r_pri
  rho_pri <- prior$rho_pri
  lambda_pri <- prior$lambda_pri

  beta_0 <- initial$beta_0
  sigma2_0 <- initial$sigma2_0
  rho_0 <- initial$rho_0
  lambda_0 <- initial$lambda_0

  output <- hom_general_int(y,X,W1,W2,nsim,burn,step,b_pri,B_pri,r_pri,lambda_pri,beta_0,sigma2_0,rho_0,lambda_0,kernel="normal",mateq = TRUE,seed=seed,impacts=impacts)

  return(output)
}


#' Hello
#'
#' @keywords internal
#'
hom_general_int=function(y,X,W1,W2=NULL,nsim,burn,step,b_pri,B_pri,r_pri,lambda_pri,beta_0,sigma2_0,rho_0,lambda_0,kernel="normal",seed=0,mateq=TRUE,impacts=TRUE)
{
  set.seed(seed)
  rowst=function(x){
    x1=c()
    x1=(x)/sum(x)}
  ########Lectura de la informaci?n
  y=as.matrix(y)
  if (is.null(X) | is.null(y) ){
    stop("No data")
  }
  #if(burn>1 | burn<0){
  #stop("Burn must be between 0 and 1. It is a proportion")
  #}
  if(nsim<=0){
    stop("There must be more than 0 simulations")
  }
  if(step<0 | step > nsim){
    stop("Jump length must not be lesser than 0 or greater than nsim")
  }
  rowst=function(x){
    x1=c()
    x1=(x)/sum(x)}

  if(mateq==TRUE){
    if(class(W1)=="nb"){
      matstand=nb2mat(W1)
      mat0=nb2listw(W1,style="B")
      mat=listw2mat(mat0)
      W2=W1
      matstand2=matstand
      mat2=mat
    }
    else{
      if(class(W1)=="listw"){
        mat=listw2mat(W1)
        matstand=apply(mat,2,rowst)
        matstand=t(matstand)
        W2=W1
        matstand2=matstand
        mat2=mat
      }
      else{
        if(sum(rowSums(W1))==nrow(X))
        {
          matstand=W1
          mat=matrix(nrow=nrow(X),ncol=nrow(X))
          for(i in 1:nrow(mat)){
            for(j in 1:ncol(mat)){
              if(matstand[i,j]==0){mat[i,j]=0}
              else{mat[i,j]=1/matstand[i,j]}
            }
          }
        }
        else{
          mat=W1
          matstand=apply(mat,2,rowst)
          matstand=t(matstand)
        }
        W2=W1
        matstand2=matstand
      }
    }
  }
  else{
    if(class(W1)=="nb"){
      matstand=nb2mat(W1)
      mat0=nb2listw(W1,style="B")
      mat=listw2mat(mat0)
    }
    else{
      if(class(W1)=="listw"){
        mat=listw2mat(W1)
        matstand=apply(mat,2,rowst)
        matstand=t(matstand)
      }
      else{
        if(sum(rowSums(W1))==nrow(X))
        {
          matstand=W1
          mat=matrix(nrow=nrow(X),ncol=nrow(X))
          for(i in 1:nrow(mat)){
            for(j in 1:ncol(mat)){
              if(matstand[i,j]==0){mat[i,j]=0}
              else{mat[i,j]=1}
            }
          }
        }
        else{
          mat=W1
          matstand=apply(mat,2,rowst)
          matstand=t(matstand)
        }
      }
    }
    if(class(W2)=="nb"){
      matstand2=nb2mat(W2)
      mat02=nb2listw(W2,style="B")
      mat2=listw2mat(mat02)
    }
    else{
      if(class(W2)=="listw"){
        mat2=listw2mat(W2)
        matstand2=apply(mat2,2,rowst)
        matstand2=t(matstand2)
      }
      else{
        if(sum(rowSums(W2))==nrow(X))
        {
          matstand2=W2
          mat2=matrix(nrow=nrow(X),ncol=nrow(X))
          for(i in 1:nrow(mat2)){
            for(j in 1:ncol(mat2)){
              if(matstand2[i,j]==0){mat2[i,j]=0}
              else{mat2[i,j]=1}
            }
          }
        }
        else{
          mat2=W2
          matstand2=apply(mat2,2,rowst)
          matstand2=t(matstand2)
        }
      }
    }
  }


  dpost <- function(betas,sigma2,rho,lambda) {
    A=diag(nrow(X))-rho*matstand
    B=diag(nrow(X))-lambda*matstand2
    k=t(B%*%(A%*%y-X%*%betas))%*%(B%*%(A%*%y-X%*%betas))
    fc.y=k
    fc.beta=t(betas - b_pri)%*%solve(B_pri)%*%(betas-b_pri)
    fc.sigma2=(lambda_pri^(r_pri))*(sigma2)^(-r_pri-1)*exp(-lambda_pri/sigma2)/gamma(r_pri)
    logdp <- (-nrow(X)/2)*log(sigma2) + log(det(A)) + log(det(B)) -0.5*fc.y/sigma2 - 0.5*fc.beta + fc.sigma2
    return(logdp)

  }
  dproposalrho<-function(rho){
    B=diag(nrow(X))-Lambda*matstand2
    a=(1/sigma2.now)*t(B%*%matstand%*%y)%*%(B%*%matstand%*%y)
    b=(1/sigma2.now)*t(B%*%(y-X%*%betas.now))%*%(B%*%matstand%*%y)
    dnorm(rho,b/a,1/sqrt(a),log = TRUE)
  }

  dproposallambda<-function(lambda){
    A=diag(nrow(X))-Rho*matstand
    a=(1/sigma2.now)*t(matstand%*%(A%*%y-X%*%betas.now))%*%(matstand%*%(A%*%y-X%*%betas.now))
    b=(1/sigma2.now)*t(A%*%y-X%*%betas.now)%*%matstand%*%(A%*%y-X%*%betas.now)
    dnorm(lambda,b/a,1/sqrt(a),log = TRUE)
  }

  ind=matrix(0,nsim,2)
  beta.mcmc=matrix(NA,nrow=nsim,ncol=ncol(X))
  sigma2.mcmc=c()
  rho.mcmc=c()
  lambda.mcmc=c()
  logV_DIC=c()
  Sigma_0=(sigma2_0)*diag(nrow(X))
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  if(kernel=="uniform"){
    for(i in 1:nsim){
      if(i==1){
        Sigma=Sigma_0
        Rho=rho_0
        Lambda=lambda_0
      }
      else{
        Sigma=sigma2.now*diag(nrow(X))
      }
      A=diag(nrow(X))-Rho*matstand
      B=diag(nrow(X))-Lambda*matstand2
      B_pos=solve(solve(B_pri)+t(X)%*%t(B)%*%solve(Sigma)%*%B%*%X)
      b_pos=B_pos%*%(t(X)%*%t(B)%*%solve(Sigma)%*%B%*%A%*%y+solve(B_pri)%*%b_pri)
      #Beta a posteriori condicional
      betas.now=c(rmvnorm(1,b_pos,B_pos))
      #A posteriori condicional completa para gammas
      r_pos=nrow(X)/2+r_pri
      k=t(B%*%(A%*%y-X%*%betas.now))%*%(B%*%(A%*%y-X%*%betas.now))
      lambda_pos=(k+2*lambda_pri)/2
      sigma2.now=rigamma(1,r_pos, lambda_pos)

      #A posteriori condicional completa para Rho
      rho.now=runif(1,1/abs(min(eigen(mat)$values)),1)
      p1=dpost(betas.now,sigma2.now,rho.now,Lambda)
      p2=dpost(betas.now,sigma2.now,Rho,Lambda)
      T.val=min(1,p1/p2)
      u<-runif(1)
      if (u <=T.val) {
        Rho <- rho.now
        ind[i,1] = 1
      }
      #A posteriori condicional completa para Lambda
      lambda.now=runif(1,1/abs(min(eigen(mat2)$values)),1)
      p1=dpost(betas.now,sigma2.now,Rho,lambda.now)
      p2=dpost(betas.now,sigma2.now,Rho,Lambda)
      T.val=min(1,p1/p2)
      u<-runif(1)
      if (u <=T.val) {
        Lambda <- lambda.now
        ind[i,2] = 1
      }

      beta.mcmc[i,]<-betas.now
      sigma2.mcmc[i]<-sigma2.now
      rho.mcmc[i]<-rho.now
      lambda.mcmc[i]<-lambda.now
      Sigma=diag(sigma2.mcmc[i],nrow(X))
      detS=det(Sigma)
      detB=det(diag(nrow(X))-lambda.mcmc[i]*matstand2)
      detA=det(diag(nrow(X))-rho.mcmc[i]*matstand)
      A=(diag(nrow(X))-rho.mcmc[i]*matstand)
      Yg=(diag(nrow(X))-lambda.mcmc[i]*matstand2)%*%(A%*%y-X%*%beta.mcmc[i,])
      logV_DIC[i]=(-(nrow(X)/2)*log(pi))+log(detA)+log(detB)-0.5*log(detS)-0.5*t(Yg)%*%solve(Sigma)%*%Yg
      Sys.sleep(0.000000001)
      # update progress bar
      setTxtProgressBar(pb, i)
    }
  }

  if(kernel=="normal")
  {

    for(i in 1:nsim)
    {
      if(i==1){
        Sigma=Sigma_0
        Rho=rho_0
        Lambda=lambda_0
      }
      else{
        Sigma=sigma2.now*diag(nrow(X))
      }
      A=diag(nrow(X))-Rho*matstand
      B=diag(nrow(X))-Lambda*matstand2
      B_pos=solve(solve(B_pri)+t(X)%*%t(B)%*%diag(1/diag(Sigma))%*%B%*%X)
      b_pos=B_pos%*%(t(X)%*%t(B)%*%diag(1/diag(Sigma))%*%B%*%A%*%y+solve(B_pri)%*%b_pri)
      #Beta a posteriori condicional
      betas.now=c(rmvnorm(1,b_pos,B_pos))
      #A posteriori condicional completa para gammas
      r_pos=nrow(X)/2+r_pri
      k=t(B%*%(A%*%y-X%*%betas.now))%*%(B%*%(A%*%y-X%*%betas.now))
      lambda_pos=(k+2*lambda_pri)/2
      sigma2.now=rigamma(1,r_pos, lambda_pos)

      #A posteriori condicional completa para Rho
      a_rho=(1/sigma2.now)*t(B%*%matstand%*%y)%*%(B%*%matstand%*%y)
      b_rho=(1/sigma2.now)*t(B%*%(y-X%*%betas.now))%*%(B%*%matstand%*%y)
      rho.now=rnorm(1,b_rho/a_rho,1/sqrt(a_rho))
      while(det(diag(nrow(X))-rho.now*matstand)<=0 || rho.now <= -1 || rho.now>=1){
        rho.now <- rnorm(1,b_rho/a_rho,1/sqrt(a_rho))
      }
      # while(rho.now>1){
      #   rho.now=rnorm(1,b/a,1/sqrt(a))
      # }


      p1.1=dpost(betas.now,sigma2.now,rho.now,Lambda)
      p2.1=dpost(betas.now,sigma2.now,Rho,Lambda)
      q1.1=dproposalrho(rho.now)
      q2.1=dproposalrho(Rho)

      met.a1 <- ifelse(p1.1>p2.1,log(p1.1-p2.1),-log(p2.1-p1.1))
      met.b1 <- ifelse(q1.1>q2.1,log(q1.1-q2.1),-log(q2.1-q1.1))
      T.val1=min(0,met.a1+met.b1)

      u<-runif(1)
      if (u <=exp(T.val1)) {
        Rho <- rho.now
        ind[i,1] = 1
      }

      #A posteriori condicional completa para Lambda
      A=diag(nrow(X))-Rho*matstand
      a_lambda=(1/sigma2.now)*t(matstand2%*%(A%*%y-X%*%betas.now))%*%(matstand2%*%(A%*%y-X%*%betas.now))
      b_lambda=(1/sigma2.now)*t(matstand2%*%(A%*%y-X%*%betas.now))%*%(A%*%y-X%*%betas.now)
      lambda.now=rnorm(1,b_lambda/a_lambda,1/sqrt(a_lambda))
      while(det(diag(nrow(X))-lambda.now*matstand)<=0){
        lambda.now <- rnorm(1,b_lambda/a_lambda,1/sqrt(a_lambda))
      }
      # while(lambda.now>1){
      #   lambda.now=rnorm(1,b/a,1/sqrt(a))
      # }

      p1.2=dpost(betas.now,sigma2.now,Rho,lambda.now)
      p2.2=dpost(betas.now,sigma2.now,Rho,Lambda)
      q1.2=dproposallambda(lambda.now)
      q2.2=dproposallambda(Lambda)

      met.a2 <- ifelse(p1.2>p2.2,log(p1.2-p2.2),-log(p2.2-p1.2))
      met.b2 <- ifelse(q1.2>q2.2,log(q1.2-q2.2),-log(q2.2-q1.2))
      T.val2=min(0,met.a2+met.b2)

      u<-runif(1)
      if (u <=exp(T.val2)) {
        Lambda <- lambda.now
        ind[i,2] = 1
      }
      beta.mcmc[i,]<-betas.now
      sigma2.mcmc[i]<-sigma2.now
      rho.mcmc[i]<-rho.now
      lambda.mcmc[i]<-lambda.now
      Sigma=diag(sigma2.mcmc[i],nrow(X))
      detS=det(Sigma)
      detB=det(diag(nrow(X))-lambda.mcmc[i]*matstand2)
      detA=det(diag(nrow(X))-rho.mcmc[i]*matstand)
      A=(diag(nrow(X))-rho.mcmc[i]*matstand)
      Yg=(diag(nrow(X))-lambda.mcmc[i]*matstand2)%*%(A%*%y-X%*%beta.mcmc[i,])
      logV_DIC[i]=(-(nrow(X)/2)*log(pi))+log(detA)+log(detB)-0.5*log(detS)-0.5*t(Yg)%*%solve(Sigma)%*%Yg
      Sys.sleep(0.000000001)
      # update progress bar
      setTxtProgressBar(pb, i)
    }
  }

  beta.mcmc_1=beta.mcmc[(burn+1):nsim,]
  sigma2.mcmc_1=sigma2.mcmc[(burn+1):nsim]
  rho.mcmc_1=rho.mcmc[(burn+1):nsim]
  lambda.mcmc_1=lambda.mcmc[(burn+1):nsim]
  beta.mcmc_2=matrix(NA,nrow=(nsim-burn+1)/step,3)
  sigma2.mcmc_2=c()
  rho.mcmc_2=c()
  lambda.mcmc_2=c()
  for (i in 1:(nsim-burn+1))
  {
    if(i%%step==0)
    {
      beta.mcmc_2[i/step,]=beta.mcmc_1[i,]
      sigma2.mcmc_2[i/step]=sigma2.mcmc_1[i]
      rho.mcmc_2[i/step]=rho.mcmc_1[i]
      lambda.mcmc_2[i/step]=lambda.mcmc_1[i]
    }
  }

  Bestimado=colMeans(beta.mcmc_2)
  Sigma2est=mean(sigma2.mcmc_2)
  Bestimado=colMeans(beta.mcmc_2)
  Sigma2est=mean(sigma2.mcmc_2)
  rho.mcmc_3=rho.mcmc_2[rho.mcmc_2<=1]
  Rhoest=mean(rho.mcmc_3)
  lambda.mcmc_3=lambda.mcmc_2[lambda.mcmc_2<=1]
  Lambdaest=mean(lambda.mcmc_3)
  Betaquant <- t(apply(beta.mcmc_2,2,function(x){quantile(x,c(0.025,0.5,0.975))}))
  Sigma2quant <- quantile(sigma2.mcmc_2,c(0.025,0.5,0.975))
  Rhoquant <- quantile(rho.mcmc_3,c(0.025,0.5,0.975))
  Lambdaquant <- quantile(lambda.mcmc_3,c(0.025,0.5,0.975))
  DesvBeta <- apply(beta.mcmc_2,2,sd)
  DesvSigma2 <- sd(sigma2.mcmc_2)
  DesvRho<-sd(rho.mcmc_3)
  DesvLambda<-sd(lambda.mcmc_3)
  AccRate1<-sum(ind[,1])/nsim
  AccRate2<-sum(ind[,2])/nsim
  Sigma=diag(Sigma2est,nrow(X))
  detS=det(Sigma)
  detA=det(diag(nrow(X))-Rhoest*matstand)
  A=(diag(nrow(X))-Rhoest*matstand)
  detB=det(diag(nrow(X))-Lambdaest*matstand2)
  Yg=(diag(nrow(X))-Lambdaest*matstand2)%*%(A%*%y-X%*%Bestimado)
  Veros=detB*detA*((detS)^(-0.5))*exp(-0.5*t(Yg)%*%solve(Sigma)%*%Yg)
  logV=log(Veros)
  logV1=-(nrow(X)/2)*log(2*pi)+log(detA)+log(detB)-0.5*log(detS)-0.5*(1/Sigma2est)*t(Yg)%*%Yg
  p=ncol(X)+3
  BIC=-2*logV1+p*log(nrow(X))
  logV_DIC=logV_DIC[is.nan(logV_DIC)==FALSE]
  Dbar=mean(-2*logV_DIC)
  logV1_DIC=(-(nrow(X)/2)*log(2*pi))+log(detB)+log(detA)-0.5*log(detS)-0.5*t(Yg)%*%solve(Sigma)%*%Yg
  Dev=-2*logV1_DIC
  DIC=2*Dbar+Dev
  # yestimado=solve(diag(nrow(X))-Rhoest*matstand)%*%X%*%Bestimado
  # residuals=(y)-yestimado

  summary = data.frame( mean=c(Bestimado,Sigma2est,Rhoest,Lambdaest),
                        sd = c(DesvBeta,DesvSigma2,DesvRho,DesvLambda),
                        q0.025=c(Betaquant[,1],Sigma2quant[1],Rhoquant[1],Lambdaquant[1]),
                        q0.5=c(Betaquant[,2],Sigma2quant[2],Rhoquant[2],Lambdaquant[2]),
                        q0.975=c(Betaquant[,3],Sigma2quant[3],Rhoquant[3],Lambdaquant[3]))

  rownames(summary) = c("x0","x1","x2","sigma2","rho","lambda")


  if(impacts){
    n <- nrow(X)
    V <-  pblapply(1:(nsim-burn),function(x){solve(diag(n)-rho.mcmc_1[x]*matstand)})
    S <- lapply(1:ncol(X),function(y){
      pblapply(1:(nsim-burn), function(x){
        beta.mcmc_1[x,y]*V[[x]]},cl=2)})
    retain <- rep(c(rep(0,step-1),1),(nsim-burn)/step)
    S <- lapply(S,function(x){x[which(retain>0)]})
    impact.df <- lapply(S, function(x){
      lapply(x, function(y){
        impact.direct <- sum(diag(y))/n
        impact.total<- sum((y))/n
        impact.indirect <- impact.total - impact.direct
        return(data.frame(Direct=impact.direct,Total=impact.total,Indirect=impact.indirect))
      })})
    impacts <- t(sapply(impact.df, function(x){tmp <- do.call("rbind",x)
    colMeans(tmp)
    }))
    rownames(impacts) <- rownames(summary)[1:ncol(X)]
    impacts <- impacts[-1,]
  }

  return(list(summary=summary,Acceptance_Rates=list(Rho_AccRate=AccRate1,Lambda_AccRate=AccRate2),Criteria=list(BIC=BIC,DIC=DIC),chains=mcmc(data.frame(beta_chain=beta.mcmc,sigma2_chain=sigma2.mcmc,rho_chain=rho.mcmc,lambda_chain=lambda.mcmc),thin = 1),
              impacts=impacts))
}

