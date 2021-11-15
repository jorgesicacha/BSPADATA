\name{hom_sem}
\alias{hom_sem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bayesian fitting of Spatial Error Model (SEM) with homoscedastic normal error term.
%%  ~~function to do ... ~~
}
\description{
Performs the Bayesian fitting of Homoscedastic Spatial Error Model (SEM) with normal error term
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hom_sem(formula, data, W, nsim, burn, step,prior, initial, kernel = "normal",seed=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ Object of class formula with the description of the model to be fitted.
%%     ~~Describe \code{y} here~~
}
  \item{data}{ Data frame object with covariates of model
%%     ~~Describe \code{X} here~~
}
  \item{W}{ Object of class matrix, nb or listw related to Spatial Contiguity Matrix, Anselin(1988)
%%     ~~Describe \code{W} here~~
}
  \item{nsim}{
  A number that indicates the amount of iterations
%%     ~~Describe \code{nsim} here~~
}
  \item{burn}{
  A number that indicates the amount of iterations to be burn at the beginning of the chain
%%     ~~Describe \code{burn} here~~
}
  \item{step}{
%%     ~~Describe \code{step} here~~
A number that indicates the length between samples in chain that generate the point estimates for each parameter.
}
  \item{prior}{
  List with prior distributions of the parameters in the hom_sem model
  }
  \item{initial}{
  List with the initial values for the chains of each the parameters of the hom_sem model
  }
  \item{kernel}{
%%     ~~Describe \code{kernel} here~~
Distribution used in transition kernel to get samples of lambda, it can be  "uniform" or "normal"
}
\item{seed}{
Random seed for generating the samples of the posterior distributions.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
hom_sem is a function made in order to fit Spatial Error Model (SEM) with a normal homoscedatic disturbance term through MCMC methods as Metropolis-Hastings algorithm, under two proposals for trasition kernel to get samples of spatial error lag parameter, lambda.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
List with the following:
\item{summary}{Data frame with summary statistics of the marginal posterior distributions of the parameters of the model}
\item{Acceptance_Rate}{Acceptance rate for the samples of lambda}
\item{Criteria}{List with values of both the Bayesian Information Criterion (BIC) and the Deviance Information Criterion (DIC)}
\item{chains}{Object of class mcmc with the samples of the marginal posterior distribution of each of the parameters of the model}
}
\references{
%% ~put references to the literature/web site here ~
1. Cepeda C. E. (2001). Modelagem da variabilidade em modelos lineares generalizados. Unpublished Ph.D. tesis. Instituto de Matematicas.
Universidade Federal do Rio do Janeiro.

2.Cepeda, E. and Gamerman D. (2005). Bayesian Methodology for modeling parameters in the two-parameter exponential family. Estadistica 57, 93 105.

3.Cepeda C., E. and Gamerman D. (2001). Bayesian Modeling of Variance Heterogeneity in Normal Regression Models. Brazilian Journal of Probability and Statistics. 14, 207-221.

4.Luc Anselin, Spatial Econometrics: Methods and Models, Kluwer Academic, Boston, 1988.

5. D. Gamerman, Markov Chains Monte Carlo: Stochastic Simulation for bayesian Inference, Chapman and Hall, 1997.

6. James Le Sage and Kelley Pace, Introduction to Spatial Econometrics, Chapman & Hall/CRC, Boca Raton, 2009.
}
\author{
Jorge  Sicacha-Parada <jasicachap@unal.edu.co>, Edilberto Cepeda-Cuervo <ecepedac@unal.edu.co>
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(oldcol)
n=49
x0=rep(1,n)
x1=runif(n,0,400)
x2=runif(n,10,23)
X=cbind(x0,x1,x2)
sigma2=rep(45,n)
Sigma=diag(sigma2)
W=COL.nb
matstand=nb2mat(W)
A=diag(n)-0.85*matstand
mu=(18+0.026*x1-0.4*x2)
Sigma2=t(solve(A))\%*\%Sigma\%*\%solve(A)
y=t(rmvnorm(1,mu,Sigma2))
formula <- y ~ x0+x1+x2
data <- data.frame(y=y,x0=x0,x1=x1,x2=x2)
prior <- list(b_pri=rep(0,3),B_pri=diag(rep(1000,3)),r_pri=0.01,lambda_pri=0.01)
initial <- list(beta_0=rep(0,3),sigma2_0=90,lambda_0=0.5)

hom_sem(formula=formula,data=data,W=W,nsim=500,burn=25,step=5,prior=prior,initial=initial,kernel="normal")

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Bayesian}
\keyword{Metropolis Hastings}
\keyword{Spatial Econometric Models}
\keyword{SEM Model}
