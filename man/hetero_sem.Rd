\name{hetero_sem}
\alias{hetero_sem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bayesian fitting of Spatial Error Model (SEM) model with heteroscedastic normal error term.
%%  ~~function to do ... ~~
}
\description{
Performs the Bayesian fitting of Heterocedastic Spatial Error Model (SEM) model with normal error term
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hetero_sem(formulamean,formulavar,data,W,nsim,burn,step,prior,initial,kernel="normal",seed=0)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{formulamean}{ Object of class formula with the description of the model to be fitted for the mean.
%%     ~~Describe \code{y} here~~
}
  \item{formulavar}{ Object of class formula with the description of the model to be fitted for the variance.
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
hetero_sem is a function made in order to fit Spatial Error Model (SEM) with a normal heteroscedatic disturbance term through MCMC methods as Metropolis-Hastings algorithm, under two proposals for trasition kernel to get samples of spatial lag parameter, lambda, and aided by working variables approach to get samples of conditional posterior distribution of gamma vector.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
List with the following:
\item{summary}{Data frame with summary statistics of the marginal posterior distributions of the parameters of the model}
\item{Acceptance_Rate}{List with the acceptance rate for the samples of gamma and lambda }
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
x3=runif(n,0,10)
X=cbind(x0,x1,x2)
Z=cbind(x0,x1,x3)
gammas=c(-8,0.026,-0.4)
Sigma=diag(c(exp(Z\%*\%gammas)))
W=COL.nb
matstand=nb2mat(W)
A=diag(n)-0.75*matstand
mu=-35+0.35*x1-1.7*x2
Sigma2=t(solve(A))\%*\%Sigma\%*\%solve(A)
y=t(rmvnorm(1,mu,Sigma2))

formulamean <- y ~ x0 + x1 + x2
formulavar <- ~ x0 + x1 +x3
data <- data.frame(y=y,x0=x0,x1=x1,x2=x2,x3=x3)
prior = list(b_pri=rep(0,3),B_pri=diag(rep(1000,3)),g_pri=rep(0,3),
             G_pri=diag(rep(1000,3)))
initial = list(beta_0=rep(0,3),gamma_0=c(0,0,0),lambda_0=0.5)

hetero_sem(formulamean=formulamean,formulavar=formulavar,data=data,W=W,nsim=500,burn=25,step=5,prior = prior,
           initial = initial,kernel="normal")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Bayesian}
\keyword{Metropolis Hastings}
\keyword{Spatial Econometric Models}
\keyword{Heteroscedastic Error Term}
\keyword{Heteroscedastic SEM Model}
