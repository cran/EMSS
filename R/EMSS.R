#' EM type Estimation Methods for the Heckman's Sample Selection Model
#'
#' Some algorithms: \code{ECM}, \code{ECMnr} and \code{ECME} can be used to estimate parameters
#' in Heckman selection model and contain the advantages of the \code{EM} algorithm: easy
#' implementation and numerical stability. \code{"ECMnr"} stands for Expectation/Conditioncal
#' Maximization with Newton-Raphson, and \code{"ECME"} for Expectation/Conditional Maximization Either.
#'
#' @details
#' The dependent variable of the selection equation (specified by argument selection) must have exactly
#' two levels (e.g., 'FALSE' and 'TRUE', or '0' and '1'). The default argument method is "ECM" and the
#' default start values ("NULL") are obtained by two-step estimation of this model through the command
#' \code{selection} from the package \code{sampleSelection}. NA's are allowed in the data. These are
#' ignored if the corresponding outcome is unobserved, otherwise observations which contain NA
#' (either in selection or outcome) are changed to 0.
#'
#' @param response a formula for the response equation.
#' @param selection a formula for the selection equation.
#' @param data a data frame and data has to be included with the form of \code{data.frame}.
#' @param method a character indicating which method to be used. \code{ECM} stands for Expectation Conditional Maximization, and
#' \code{ECMnr} stands for Expectation Conditioncal Maximization with Newton-Raphson, and \code{ECME} for Expectation Conditional Maximization Either.
#' @param initial.param a vector, initial parameter values for the estimation. The length of the initial parameters has to
#' be same as the length of parameters, which are to be estimated.
#' @param eps a numerical error value for the end of the loop. A minimum value that can be arbitrarily set
#' to terminate the iteration of the function, in order to find the optimal parameter estimation.
#' @return \code{ECM} returns an object of class \code{"ECM"}.
#' The object class \code{"ECM"} is a list
#' containing the following components.
#' @return
#' \item{call}{a matched call.}
#' \item{estimate_response}{estimated regression coefficients for the response formula.}
#' \item{estimate_selection}{estimated regression coefficients for the sample selection formula.}
#' \item{estimate_sigma}{an estimated scale paramter for the bivariate normal distribution.}
#' \item{estimate_rho}{an estimated correlation coefficient for the bivariate normal distribution.}
#' \item{hessian_mat}{hessian matrix for parameters.}
#' \item{resp_leng}{the numbers of coefficients for the response formula}
#' \item{select_leng}{the numbers of coefficients for the selection formula}
#' \item{Q_value}{the vallue of the Q function for EM type algorithms}
#' \item{names_response}{names of regression coefficients for the reponse formula.}
#' \item{names_selection}{names of regression coefficients for the selection formula.}
#' @importFrom mvtnorm rmvnorm
#' @importFrom sampleSelection selection
#' @importFrom stats dnorm model.matrix model.response na.pass pnorm qchisq qnorm symnum complete.cases
#' @examples
#' data(Smoke, package = "EMSS")
#' ex1 <- EMSS(response = cigs_intervals ~ educ,
#'            selection = smoker ~ educ + age,
#'            data = Smoke)
#' print(ex1)
#'
#' @examples
#' data(Smoke, package = "EMSS")
#' ex2 <- EMSS(response = cigs_intervals ~ educ,
#'            selection =  smoker ~ educ + age,
#'            data = Smoke, method="ECMnr")
#' print(ex2)
#'
#' @examples
#' ## example using random numbers with exclusion restriction
#'
#' N <- 1000
#' errps <- mvtnorm::rmvnorm(N,c(0,0),matrix(c(1,0.5,0.5,1),2,2) )
#' xs <- runif(N)
#' ys <- xs+errps[,1]>0
#' xo <- runif(N)
#' yo <- (xo+errps[,2])*(ys>0)
#' ex3 <- EMSS(response = yo ~ xo,
#'            selection = ys ~ xs,
#'            initial.param = c(rep(0,4), 0.3, 0.6), method="ECMnr")
#' print(ex3)
#'
#' @references Heckman, J. (1979) Sample selection bias as a specication error. \emph{Econometrica}, 47, 153-161.
#' @references Toomet, O. and Henningsen, A. (2008) Sample selection models in R:Package sampleSelection. \emph{Journal of Statistical Software}, 27, 1-23.
#' @references Zhao,J., Kim, H.-J. and Kim, H.-M. (2020) New EM-type algorithms for the Heckman selection model. \emph{Computational Statistics and Data Analysis}, 146, https://doi.org/10.1016/j.csda.2020.106930.
#' @references Zhelonkin, M., Genton, M.G. and Ronchetti, E. (2016) Robust inference in sample selection models. \emph{Journal of the Royal Statistical Society Series B}, 78, 805-827.
#'
#' @section Background:
#' Heckman selection model is classic to deal with the data where the outcome is partially observed and
#' the missing part is not at random. Heckman (1979) developed \code{2-step} and maximum likelihood
#' estimation (\code{MLE}) to do the estimation for this selection model. And these two method are
#' described in R package \code{sampleSelection} by Toomet and Henningsen (2008). Zhelonkin et al. (2016)
#' developed robust 2-stage method which performs more robustly than the 2-step method to deal with the
#' data where outlying observations exist and \code{ssmrob} package is available. Zhao et al. (2020) extended
#' EM algorithm to more general cases resulting in three algorithms: ECM, ECM(NR), and ECME. They also own
#' EM algorithm's main advantages, namely, stability and ease of implementation.
#' @export
EMSS <- function(response, selection, data, method="ECM",
                 initial.param = NULL,
                 eps = 10^(-10))
{
  if(method!="ECM" & method!="ECMnr" & method!="ECME"){
    message("Warning: Wrong maximization method. It runs by default.")
    method <- "ECM"
  }

  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)

  m1 <- match(c("response", "data"), names(mf), 0)
  mf1 <- mf[c(1, m1)]
  mf1$drop.unused.levels <- TRUE
  mf1$na.action <- na.pass
  mf1[[1]] <- quote(model.frame)
  names(mf1)[2] <- "formula"
  mf1 <- eval(mf1, parent.frame())
  mt1 <- attr(mf1, "terms")

  m2 <- match(c("selection", "data"), names(mf), 0)
  mf2 <- mf[c(1, m2)]
  mf2$drop.unused.levels <- TRUE
  mf2$na.action <- na.pass
  names(mf2)[2] <- "formula"
  mf2[[1]] <- quote(model.frame)
  mf2 <- eval(mf2, parent.frame())
  mt2 <- attr(mf2, "terms")

  y1 <- model.response(mf1, "numeric")
  y2 <- model.response(mf2, "numeric")

  x <- model.matrix(mt1, mf1, contrasts.arg = NULL, xlev = NULL)
  w <- model.matrix(mt2, mf2, contrasts.arg = NULL, xlev = NULL)
  
  complete <- complete.cases(x) & complete.cases(w)
  
  if(sum(!complete) > 0){
    x <- x[complete,]
    w <- w[complete,]
    y1 <- y1[complete]
    y2 <- y2[complete]
    warning("The response data or the selection data is not complete, the missing observations are ignored.")
  }
  
  x.name <- names(as.data.frame(x))
  w.name <- names(as.data.frame(w))

  x <- t(x)
  w <- t(w)

  N <- length(y2)

  obserindex <- which(y2>0)
  N1 <- length(obserindex) ; N2<-N-N1
  xy1obser <- x[,obserindex]
  xy1miss <- x[,-obserindex]
  x <- cbind(xy1obser,xy1miss )
  #put all wi that yi1(where yi2>0) observed in front
  wy1obser <- w[,obserindex]
  wy1miss <- w[,-obserindex]
  w <- cbind(wy1obser,wy1miss)
  largX <- t(x); largW <- t(w)

  y1obser <- matrix(y1[obserindex],ncol=1)#y1|y2>0
  y1miss <- matrix(rep(0,N2),ncol=1) #y1|y2<=0
  y1 <- drop(rbind(y1obser,y1miss))

  #put all y2i that yi1(where yi2>0) observed in front
  y2obser <- matrix(y2[obserindex],ncol=1)
  y2miss <- matrix(y2[-obserindex],ncol=1)
  y2 <- drop(rbind(y2obser,y2miss))
  u <- c(rep(1,N1),rep(0,N2))
  udiag <- diag(u)

  p <- dim(x)[1]
  q <- dim(w)[1]

  #-------------------------------------------------------------------------------------------------
  #estimate by two-step method------------------------------------------------------
  #-------------------------------------------------------------------------------------------------
  # to get initial values for ECM algorithm
  if(is.null(initial.param)){
    twostepsele <- selection( y2 ~ largW[,-1], y1 ~ largX[,-1], data =
                                as.data.frame(cbind(y1,y2,largX[,-1],largW[,-1])), method = "2step" )
    #initial values from two step---------------------------------------------------
    gammakp1 <- matrix(twostepsele$coefficients[1:q] ,ncol=1)
    betakp1 <- matrix(twostepsele$coefficients[(q+1):(p+q)],ncol=1)
    sigmakp1 <- as.vector(twostepsele$coefficients[p+q+2])
    rhokp1 <- as.vector(twostepsele$coefficients[p+q+3])
  }

  else {     if(length(initial.param) != p+q+2) stop("The length of the initial parameters has to be same as the parameters")
    betakp1 <- matrix(initial.param[1:p], ncol = 1)
    gammakp1 <- matrix(initial.param[(p+1):(p+q)], ncol = 1)
    sigmakp1 <- as.vector(initial.param[p+q+1])
    rhokp1 <- as.vector(initial.param[p+q+2])
  }

  rhokp1[which(abs(rhokp1)>1)] <- 0.9 ##if initial value of rho is larger than 1, set it to 0.9
  # we need to give the warning or the description if rho is bigger than 0.9
  psikp1 <- sigmakp1^2*(1-rhokp1^2)
  psistarkp1 <- log(sigmakp1^2*(1-rhokp1^2))
  rhostarkp1 <- sigmakp1*rhokp1

  #-------------------------------------------------------------------------------------------------
  #ECM algorithom---------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------------------------

  if(method=="ECM"){
    est_result = ECM.fit(betakp1, gammakp1, sigmakp1, rhokp1,
                         psikp1, psistarkp1, rhostarkp1,
                         xy1obser, wy1obser, y1obser, xy1miss, wy1miss, y1miss,
                         largX, largW, N, eps, err=100)
  }
  else if(method=="ECMnr"){
    est_result = ECMnr.fit(betakp1, gammakp1, sigmakp1, rhokp1,
                           psikp1, psistarkp1, rhostarkp1,
                           xy1obser, wy1obser, y1obser, xy1miss, wy1miss, y1miss,
                           largX, largW, N, eps, err=100)
  }
  else if(method=="ECME"){
    est_result = ECME.fit(betakp1, gammakp1, sigmakp1, rhokp1,
                          psikp1, psistarkp1, rhostarkp1,
                          xy1obser, wy1obser, y1obser, xy1miss, wy1miss, y1miss,
                          largX, largW, w, y1, x, u, N1, udiag, eps, err=100)
  }

  betaest <- as.matrix(est_result$betakp1)
  gammaest <- as.matrix(est_result$gammakp1)
  sigmaest <- drop(est_result$sigmakp1)
  rhoest <- drop(est_result$rhokp1)

  ### Q value---------------------------------------------------------------------------------
  if(method == "ECM" | method == "ECMnr"){
    Qfun.val = Q_ECM(betaest, gammaest, sigmaest, rhoest,
                  xy1obser, wy1obser, y1obser,
                  xy1miss, wy1miss, y1miss, N)
  } else if (method == "ECME") {
    Qfun.val = Q_ECME(betaest, gammaest, sigmaest, rhoest,
                      xy1obser, wy1obser, y1obser,
                      xy1miss, wy1miss, y1miss,
                      y1, x, w, u)
  }

  ### Hessian matrix--------------------------------------------------------------------------
  termS<-rhoest/sqrt(1-rhoest^2)*(y1-crossprod(x,betaest) )/sigmaest+crossprod(w,gammaest)/sqrt(1-rhoest^2)
  lambdaS<-dnorm(termS)/pnorm(termS)
  #1 partial logelikelihood/parialbetaest partial betaest
  parbetabeta<-diag(as.vector((-1-rhoest^2*(termS*lambdaS+lambdaS^2)/(1-rhoest^2) ) *u))
  parbetabeta<-x%*%parbetabeta%*% t(x) /(sigmaest^2)

  #2 partial logelikelihood/parial beta partial gamma
  parbetagamma<-diag(as.vector((termS*lambdaS+lambdaS^2)*u))
  parbetagamma<- x%*%parbetagamma%*%t(w) *rhoest/(sigmaest*(1-rhoest^2))

  #3 partial logelikelihood/parial beta partial sigma
  ##{}term in parbetasigma
  parterm<- -2*(y1-crossprod(x,betaest) )/sigmaest+rhoest*lambdaS/sqrt(1-rhoest^2)
  parterm<- parterm-rhoest^2/(sigmaest*(1-rhoest^2) )*(y1-crossprod(x,betaest) )*(termS*lambdaS+lambdaS^2)
  parbetasigma<-x%*%(parterm*u)/(sigmaest^2)

  #4 partial logelikelihood/parial beta partial rho
  ##cut {}term in parbetarho to three parts
  parbetarho1<- -lambdaS*( (1-rhoest^2)^(-1/2)+rhoest^2*((1-rhoest^2)^(-3/2) ) )/sigmaest
  parbetarho2<- rhoest/(sigmaest*sqrt(1-rhoest^2)) *(termS*lambdaS+lambdaS^2)
  parbetarho3<- ( (1-rhoest^2)^(-1/2)+rhoest^2*((1-rhoest^2)^(-3/2) ) )*(y1-crossprod(x,betaest) )/sigmaest
  parbetarho3<- parbetarho3+rhoest*((1-rhoest^2)^(-3/2) )*crossprod(w,gammaest)
  parbetarho<-x%*%((parbetarho1+parbetarho2*parbetarho3)*u)

  #5 partial logelikelihood/parial gamma partial gamma
  lambdaMwgmma<-dnorm(-crossprod(w,gammaest) )/pnorm(-crossprod(w,gammaest) )
  pargammagamma0<- -u*(termS*lambdaS+lambdaS^2)/(1-rhoest^2)
  pargammagamma<-diag(as.vector(pargammagamma0+(1-u)*(crossprod(w,gammaest)*lambdaMwgmma-lambdaMwgmma^2 ) ))
  pargammagamma<-w%*%pargammagamma%*%t(w)

  #6 partial logelikelihood/parial gamma partial sigma
  #Medianterm
  pargammasigmaM<-(y1-crossprod(x,betaest) )*(termS*lambdaS+lambdaS^2)
  pargammasigma<-rhoest/(sigmaest^2*(1-rhoest^2) )*w%*%(pargammasigmaM*u)

  #7 partial logelikelihood/parial gamma partial rho
  #cut {}term to two parts
  pargammarho1<-rhoest*((1-rhoest^2)^(-3/2) )*lambdaS
  pargammarho2<-( (1-rhoest^2)^(-1/2)+rhoest^2*((1-rhoest^2)^(-3/2) ) )*(y1-crossprod(x,betaest) )/sigmaest
  pargammarho2<-pargammarho2 +rhoest*((1-rhoest^2)^(-3/2) )*crossprod(w,gammaest)
  pargammarho2<-(termS*lambdaS+lambdaS^2)*pargammarho2/sqrt(1-rhoest^2)
  pargammarho<- w%*%((pargammarho1-pargammarho2 )*u)

  #8 partial logelikelihood/parial sigma partial sigma
  #cut {}term to three parts
  parsigmasigma1<-1/(sigmaest^2)-3*(y1-crossprod(x,betaest) )^2/(sigmaest^4)
  parsigmasigma2<-lambdaS*rhoest/sqrt(1-rhoest^2)*2*(y1-crossprod(x,betaest) )/(sigmaest^3)
  parsigmasigma3<-rhoest^2/(1-rhoest^2)*(y1-crossprod(x,betaest) )^2/(sigmaest^4)*(termS*lambdaS+lambdaS^2)
  parsigmasigma<-crossprod(u,(parsigmasigma1+parsigmasigma2-parsigmasigma3 ) )


  #9 partial logelikelihood/parial sigma partial rho
  #cut {}term to two parts
  #common term
  commterm<-(y1-crossprod(x,betaest) )*( (1-rhoest^2)^(-1/2)+rhoest^2*((1-rhoest^2)^(-3/2) ) )
  parsigmarho1<- -commterm*lambdaS/(sigmaest^2)
  parsigmarho2<- rhoest/sqrt(1-rhoest^2)*(y1-crossprod(x,betaest) )/(sigmaest^2)*(termS*lambdaS+lambdaS^2)
  parsigmarho2<-parsigmarho2*(commterm/sigmaest+rhoest*crossprod(w,gammaest)*((1-rhoest^2)^(-3/2) )  )
  parsigmarho<-crossprod(u,( parsigmarho1+parsigmarho2) )

  #10 partial logelikelihood/parial rho partial rho
  #cut {}term to two parts
  parrhorho1<-((y1-crossprod(x,betaest) )/sigmaest*(3*rhoest*((1-rhoest^2)^(-3/2) )+3*rhoest^3*((1-rhoest^2)^(-5/2) ))
               +((1-rhoest^2)^(-3/2)+3*rhoest^2*(1-rhoest^2)^(-5/2) )*crossprod(w,gammaest) )*lambdaS
  parrhorho2<- (termS*lambdaS+lambdaS^2)*(commterm/sigmaest+rhoest*((1-rhoest^2)^(-3/2) )*crossprod(w,gammaest) )^2 # fix -1/2 to -1
  parrhorho<-crossprod(u, ( parrhorho1-parrhorho2) )

  HessianMa<-matrix(0,nrow=p+q+2,ncol=p+q+2)
  HessianMa[1:p,1:p]<-parbetabeta
  HessianMa[1:p,(p+1):(p+q)]<-parbetagamma
  HessianMa[(p+1):(p+q),1:p]<-t(parbetagamma)
  HessianMa[1:p,(p+q+1)]<-parbetasigma
  HessianMa[(p+q+1),1:p]<-t(parbetasigma)
  HessianMa[1:p,(p+q+2)]<-parbetarho
  HessianMa[(p+q+2),1:p]<-t(parbetarho)
  HessianMa[(p+1):(p+q),(p+1):(p+q)]<-pargammagamma
  HessianMa[(p+1):(p+q),(p+q+1)]<-pargammasigma
  HessianMa[(p+q+1),(p+1):(p+q)]<-t(pargammasigma)
  HessianMa[(p+1):(p+q),(p+q+2)]<-pargammarho
  HessianMa[(p+q+2),(p+1):(p+q)]<-t(pargammarho)
  HessianMa[(p+q+1),(p+q+1)]<-parsigmasigma
  HessianMa[(p+q+1),(p+q+2)]<-parsigmarho
  HessianMa[(p+q+2),(p+q+1)]<-parsigmarho
  HessianMa[(p+q+2),(p+q+2)]<-parrhorho


  ECMsele <- list( call = cl,
                   estimate_response = est_result$betakp1,
                   estimate_selection = est_result$gammakp1,
                   estimate_sigma = sigmaest,
                   estimate_rho = rhoest,
                   hessian_mat = HessianMa,
                   resp_leng = p,
                   select_leng = q,
                   Q_value = Qfun.val
  )

  ECMsele$names_response <- c( x.name )
  ECMsele$names_selection <- c( w.name )

  class(ECMsele)<- "EMSS"
  ECMsele
}

ECM.fit <- function(betakp1, gammakp1, sigmakp1, rhokp1,
                    psikp1, psistarkp1, rhostarkp1,
                    xy1obser, wy1obser, y1obser,
                    xy1miss, wy1miss, y1miss,
                    largX, largW, N, eps, err=100
){
  while(err>eps){

    betak<-betakp1;gammak<-gammakp1;sigmak<-as.vector(sigmakp1)
    rhok<-as.vector(rhokp1); psik<-as.vector(psikp1)

    #E step moments-------------------------------------------------------------------------
    #yi1 observed  crossprod t(wy1obser)%*%gammakt(xy1obser)%*%betak
    mui21k<- crossprod(wy1obser,gammak)+rhok*(y1obser-crossprod(xy1obser,betak))/sigmak
    #lambda function
    mu21rho<-mui21k/sqrt(1-rhok^2)
    denomina<-pnorm(mu21rho)
    denomina[which(denomina==0),1]<-pnorm(-29)
    lambdamu21rho<-dnorm(mu21rho)/denomina

    #yi1 observed  equation 6-------------------
    alpha2ok<-mui21k+sqrt(1-rhok^2)*lambdamu21rho
    v2ok<-1-rhok^2+mui21k^2+mui21k*sqrt(1-rhok^2)*lambdamu21rho

    #yi1 missing
    mui1k<-crossprod(xy1miss,betak)#t(xy1miss)%*%betak
    mui2k<-crossprod(wy1miss,gammak) #t(wy1miss)%*%gammak
    #lambda function
    denomina2<-pnorm(-mui2k)
    denomina2[which(denomina2==0),1]<-pnorm(-29)
    lambdammui2<-dnorm(-mui2k)/denomina2
    rhostark<-rhok*sigmak

    #yi1 missing equations 1-5-------------------
    alpha1mk<-mui1k-rhostark*lambdammui2
    alpha2mk<-mui2k-lambdammui2
    v1mk<-mui1k^2+sigmak^2-rhostark*lambdammui2*(2*mui1k-rhostark*mui2k)
    v2mk<-1+mui2k^2-mui2k*lambdammui2
    alpha12mk<-mui1k*(mui2k-lambdammui2)+rhostark


    #CM step update parameters---------------------------------------------------------------
    #update beta(k)-------------------------------------------------------------------------
    #calculate the term in {} firstly
    #t(largX)%*%largW%*%gammak
    betakp1par<-crossprod(largX,largW)%*%gammak-xy1obser%*%alpha2ok-xy1miss%*% alpha2mk
    betakp1<-solve(crossprod(largX,largX), rhostark*betakp1par+xy1obser%*%y1obser+xy1miss%*%alpha1mk )

    #update gamama(k)-----------------------------------------------------------------------
    #calculate the term in {} firstly
    gammakp1par<-crossprod(largW,largX)%*%betakp1-wy1obser%*%y1obser-wy1miss%*%alpha1mk

    gammakp1<-solve(crossprod(largW), (rhostark*gammakp1par/(psik+rhostark^2)
                                       +wy1obser%*%alpha2ok+wy1miss%*%alpha2mk) )

    #update psi(k)----------------------------------------------------------------------
    #calcualte circled1
    #circled10<-t(y1obser-t(xy1obser)%*%betakp1 )%*%( y1obser-t(xy1obser)%*%betakp1 )
    circled10<-crossprod(y1obser-crossprod(xy1obser,betakp1),(y1obser-crossprod(xy1obser,betakp1) ) )
    #reusable function--------
    xy1misMulbetakp1<-crossprod(xy1miss,betakp1)
    sumpar11<-sum(v1mk-2*xy1misMulbetakp1*alpha1mk
                  +xy1misMulbetakp1*(xy1misMulbetakp1) )
    circled1<-circled10+sumpar11

    #calculate circled2
    #t(wy1obser)%*%gammakp1
    sumterm21<-sum(v2ok-2*crossprod(wy1obser,gammakp1)*alpha2ok)
    sumterm22<-sum(v2mk-2*crossprod(wy1miss,gammakp1)*alpha2mk )
    #circled2<-t(gammakp1)%*%t(largW)%*%largW%*%gammakp1+ sumterm21+ sumterm22
    circled2<-crossprod(gammakp1,crossprod(largW))%*%gammakp1+ sumterm21+ sumterm22

    #calculate circled3
    sumterm31<-sum(y1obser*alpha2ok-crossprod(wy1obser,gammakp1)*y1obser
                   -crossprod(xy1obser,betakp1)*alpha2ok )

    sumterm32<-sum(alpha12mk-crossprod(wy1miss,gammakp1)*alpha1mk
                   -xy1misMulbetakp1*alpha2mk)
    #circled3<-t(betakp1)%*%t(largX)%*%largW%*%gammakp1+sumterm31+sumterm32
    circled3<-crossprod(betakp1,crossprod(largX,largW) )%*%gammakp1+sumterm31+sumterm32

    #update psi
    psikp1<- (circled1+rhostark^2*circled2-2*rhostark*circled3)/N

    #update rhostar(k)----------------------------------------------------------------------
    rhostarkp1<-circled3/circled2

    #calculate sigma and rho by psi and rhostar
    sigmakp1<-as.vector( sqrt( psikp1+rhostarkp1^2 ) )
    rhokp1<- rhostarkp1/sigmakp1

    err<-sum(abs(betak-betakp1),abs(gammak-gammakp1),abs(sigmak-sigmakp1),abs(rhok-rhokp1) )
  }
  est_params = list(betakp1 = betakp1,
                    gammakp1 = gammakp1,
                    sigmakp1 = sigmakp1,
                    rhokp1 = rhokp1)
  est_params
}

ECMnr.fit <- function(betakp1, gammakp1, sigmakp1, rhokp1,
                      psikp1, psistarkp1, rhostarkp1,
                      xy1obser, wy1obser, y1obser,
                      xy1miss, wy1miss, y1miss,
                      largX, largW, N, eps, err=100){
  while(err>eps){
    betak<-betakp1;gammak<-gammakp1;sigmak<-as.vector(sigmakp1)
    rhok<-as.vector(rhokp1); psistark<-as.vector(psistarkp1)

    #E step moments-------------------------------------------------------------------------
    #yi1 observed
    mui21k <- crossprod(wy1obser,gammak)+rhok*(y1obser-crossprod(xy1obser,betak))/sigmak
    #lambda function
    mu21rho <- mui21k/sqrt(1-rhok^2)
    denomina <- pnorm(mu21rho)
    denomina[which(denomina==0),1]<-pnorm(-29)
    lambdamu21rho <- dnorm(mu21rho)/denomina

    #yi1 observed  equation 6-------------------
    alpha2ok <- mui21k+sqrt(1-rhok^2)*lambdamu21rho
    v2ok <- 1-rhok^2+mui21k^2+mui21k*sqrt(1-rhok^2)*lambdamu21rho

    #yi1 missing
    mui1k<-crossprod(xy1miss,betak)
    mui2k<-crossprod(wy1miss,gammak)
    denomina2<-pnorm(-mui2k)
    denomina2[which(denomina2==0),1]<-pnorm(-29)
    lambdammui2<-dnorm(-mui2k)/denomina2
    rhostark<-rhok*sigmak

    #yi1 missing equations 1-5----------------
    alpha1mk<-mui1k-rhostark*lambdammui2
    alpha2mk<-mui2k-lambdammui2
    v1mk<-mui1k^2+sigmak^2-rhostark*lambdammui2*(2*mui1k-rhostark*mui2k)
    v2mk<-1+mui2k^2-mui2k*lambdammui2
    alpha12mk<-mui1k*(mui2k-lambdammui2)+rhostark


    #M step update parameters---------------------------------------------------------------
    #update beta(k)-------------------------------------------------------------------------
    #calculate the term in {} firstly
    betakp1par<-crossprod(largX,largW)%*%gammak-xy1obser%*%alpha2ok-xy1miss%*% alpha2mk
    betakp1<-solve(crossprod(largX),( rhostark*betakp1par+xy1obser%*%y1obser+xy1miss%*%alpha1mk ) )

    #update gamama(k)-----------------------------------------------------------------------
    #calculate the term in {} firstly
    gammakp1par<-crossprod(largW,largX)%*%betakp1-wy1obser%*%y1obser-wy1miss%*%alpha1mk
    gammakp1<-solve(crossprod(largW),( rhostark*gammakp1par/(exp(psistark)+rhostark^2)
                                       +wy1obser%*%alpha2ok+wy1miss%*%alpha2mk ) )

    #update psistar(k)----------------------------------------------------------------------
    #calcualte circled1
    #reusable function-----------------------------
    xy1misMulbetakp1<-crossprod(xy1miss,betakp1)

    circled10<-crossprod(y1obser-crossprod(xy1obser,betakp1), (y1obser-crossprod(xy1obser,betakp1)) )
    sumpar11<-sum(v1mk-2*xy1misMulbetakp1*alpha1mk
                  +xy1misMulbetakp1*xy1misMulbetakp1 )
    circled1<-circled10+sumpar11

    #calculate circled2
    #reusable functions-----------------------------
    wy1obMulgammakp1 <- crossprod(wy1obser,gammakp1)
    wy1misMulgammakp1 <- crossprod(wy1miss,gammakp1)
    sumterm21<-sum(v2ok-2*wy1obMulgammakp1*alpha2ok)
    sumterm22<-sum(v2mk-2*wy1misMulgammakp1 *alpha2mk )
    circled2<-crossprod(gammakp1,crossprod(largW))%*%gammakp1+ sumterm21+ sumterm22

    #calculate circled3
    sumterm31<-sum(y1obser*alpha2ok-wy1obMulgammakp1*y1obser
                   -crossprod(xy1obser,betakp1)*alpha2ok )

    sumterm32<-sum(alpha12mk-wy1misMulgammakp1 *alpha1mk
                   -xy1misMulbetakp1*alpha2mk)
    circled3<-crossprod(betakp1,crossprod(largX,largW))%*%gammakp1+sumterm31+sumterm32

    #update psistar
    psik<-exp(psistark)
    commterm<- -N*psik+circled1+rhostark^2*circled2-2*rhostark*circled3
    #first derivative of psistar in Q function
    Qpartial<-commterm/(2*psik)
    #second derivatice of  psistar in Q function
    Qparpartial<- -commterm/(2*psik)-N/2

    psistarkp1<-psistark-Qparpartial^(-1)*Qpartial

    #update rhostar(k)----------------------------------------------------------------------
    rhostarkp1<-circled3/circled2

    sigmakp1<-as.vector( sqrt( exp(psistarkp1)+rhostarkp1^2 ) )
    rhokp1<- rhostarkp1/sigmakp1

    err<-sum(abs(betak-betakp1),abs(gammak-gammakp1),abs(sigmak-sigmakp1),abs(rhok-rhokp1) )
  }
  est_params = list(betakp1 = betakp1,
                    gammakp1 = gammakp1,
                    sigmakp1 = sigmakp1,
                    rhokp1 = rhokp1)
  est_params
}

ECME.fit <- function(betakp1, gammakp1, sigmakp1, rhokp1,
                     psikp1, psistarkp1, rhostarkp1,
                     xy1obser, wy1obser, y1obser,
                     xy1miss, wy1miss, y1miss,
                     largX, largW, w, y1, x, u,
                     N1, udiag, eps, err=100){
  while(err>eps ){
    #E-step---------------------------------------------------------------
    gammak<-gammakp1; betak<-betakp1; psik<-psikp1; rhostark<-rhostarkp1
    sigmak<-sigmakp1; rhok<- rhokp1 # rhostark/sigmak

    #parameters in lemma 1
    A<- -crossprod(w,gammak)
    xi<-rhok*(y1-crossprod(x,betak))/sigmak
    sigmainTN<-sqrt(1-rhok^2)
    #(xi-A)/sigma
    ximAdisigma<- (xi-A)/sigmainTN
    #denomina part in EY and VY
    denomina<-pnorm(ximAdisigma)
    denomina[which(denomina==0),1]<-pnorm(-29)

    #(dnorm((xi-A)/sigma))/denomian
    parterm<-( dnorm(ximAdisigma) )/denomina

    alphak<-EY<-xi+ parterm*sigmainTN
    VY<-( 1-ximAdisigma*parterm -parterm^2)*sigmainTN^2
    deltak<-VY+(EY)^2

    #M-step-------------------------------------------------------------------------------------
    #update psik--------------------------------------------------------------------------------

    psikp1<-as.vector(crossprod(u,( (y1-crossprod(x,betak))^2
                                    -2*rhostark*(y1-crossprod(x,betak))*alphak+(rhostark)^2*deltak ) ) )/N1

    #update rhostark----------------------------------------------------------------------------
    rhostarkp1<- as.vector(crossprod(y1-crossprod(x,betak),udiag)%*%alphak/(crossprod(u,deltak) ))

    #update betak--------------------------------------------------------------------------------
    betakp1<- solve( x%*%diag(u)%*%t(x),( x%*%udiag%*%(y1-rhostarkp1*alphak) ))

    #update gammak------------------------------------------------------------------------------
    sigmakp1<-as.vector(sqrt( psikp1+rhostarkp1^2 ));rhokp1<-as.vector(rhostarkp1/sigmakp1)
    #Newton-Raphson
    Sterm<- (sigmakp1*crossprod(w,gammak)+rhokp1*(y1-crossprod(x,betakp1) ) )/sqrt(psikp1)
    #pdfDcdfS=dnorm(Sterm)/pnorm(Sterm)
    pdfDcdfS<-dnorm(Sterm)/pnorm(Sterm)
    #pnorm/cnorm for -t(w)%*%gammak
    pdfDcdfmwg<-dnorm(-crossprod(w,gammak)) /pnorm(-crossprod(w,gammak))
    PartiallogL<-w%*%( sigmakp1/sqrt(psikp1)*udiag%*%pdfDcdfS+diag(u-1)%*%pdfDcdfmwg )
    ParparlogL<- w%*%diag(as.vector(-sigmakp1^2*u*(Sterm*pdfDcdfS+pdfDcdfS^2)/psikp1
                                    +(1-u)*( (crossprod(w,gammakp1))*pdfDcdfmwg-pdfDcdfmwg^2) ) )%*%t(w)
    gammakp1<- gammak-solve(ParparlogL,PartiallogL)

    err<-sum(abs(betakp1-betak),abs(gammakp1-gammak),abs(sigmakp1-sigmak),abs(rhokp1-rhok) )
  }
  est_params = list(betakp1 = betakp1,
                    gammakp1 = gammakp1,
                    sigmakp1 = sigmakp1,
                    rhokp1 = rhokp1)
  est_params
}

Q_ECM <- function(betaest, gammaest, sigmaest, rhoest,
                  xy1obser, wy1obser, y1obser,
                  xy1miss, wy1miss, y1miss, N
){
  psiest<-sigmaest^2*(1-rhoest^2)
  rhostarest<-rhoest*sigmaest
  #parameters arise from expectation---------------------------------------------
  #yi1 observed
  mui21k<-crossprod(wy1obser,gammaest)+rhoest*(y1obser-crossprod(xy1obser,betaest))/sigmaest
  #lambda function
  mu21rho<-mui21k/sqrt(1-rhoest^2)
  denomina<-pnorm(mu21rho)
  denomina[which(denomina==0),1]<-pnorm(-29)
  lambdamu21rho<-dnorm(mu21rho)/denomina

  alpha2ok<-mui21k+sqrt(1-rhoest^2)*lambdamu21rho
  v2ok<-1-rhoest^2+mui21k^2+mui21k*sqrt(1-rhoest^2)*lambdamu21rho

  #yi1 missing
  mui1k<-crossprod(xy1miss,betaest)
  mui2k<-crossprod(wy1miss,gammaest)
  #lambda function
  lambdammui2<-dnorm(-mui2k)/pnorm(-mui2k)

  alpha1mk<-mui1k-rhostarest*lambdammui2
  alpha2mk<-mui2k-lambdammui2
  v1mk<-mui1k^2+sigmaest^2-rhostarest*lambdammui2*(2*mui1k-rhostarest*mui2k)
  v2mk<-1+mui2k^2-mui2k*lambdammui2
  alpha12mk<-mui1k*(mui2k-lambdammui2)+rhostarest
  #parameters arise from expectation end---------------------------------------------

  ECMQfun1<- -N*log(2*pi)-N*log(psiest)/2

  ECMQfun21<- crossprod( (y1obser-crossprod(xy1obser,betaest) ),(y1obser-crossprod(xy1obser,betaest)) )
  ECMQfun22<- sum(v1mk-2*alpha1mk*crossprod(xy1miss,betaest) +crossprod(xy1miss,betaest)*(crossprod(xy1miss,betaest) ) )
  ECMQfun2<- -(ECMQfun21+ECMQfun22)/(2*psiest)

  ECMQfun31<- sum(v2ok-2*alpha2ok*crossprod(wy1obser,gammaest)+crossprod(wy1obser,gammaest)*(crossprod(wy1obser,gammaest)) )
  ECMQfun32<- sum(v2mk-2*alpha2mk*crossprod(wy1miss,gammaest) +crossprod(wy1miss,gammaest)*(crossprod(wy1miss,gammaest) ) )
  ECMQfun3<- -(1+rhostarest^2/psiest)/2*(ECMQfun31+ECMQfun32)

  ECMQfun41<- sum(y1obser*alpha2ok-y1obser*(crossprod(wy1obser,gammaest) )-alpha2ok*(crossprod(xy1obser,betaest) )
                  +crossprod(xy1obser,betaest)*(crossprod(wy1obser,gammaest) ) )
  ECMQfun42<- sum(alpha12mk-alpha1mk*(crossprod(wy1miss,gammaest) )-alpha2mk*(crossprod(xy1miss,betaest) )
                  +crossprod(xy1miss,betaest)*(crossprod(wy1miss,gammaest) ) )
  ECMQfun4<- rhostarest/psiest*(ECMQfun41+ECMQfun42)

  Q.val <- sum(ECMQfun1,ECMQfun2,ECMQfun3,ECMQfun4)
  Q.val
}

Q_ECME <- function(betaest, gammaest, sigmaest, rhoest,
                   xy1obser, wy1obser, y1obser,
                   xy1miss, wy1miss, y1miss,
                   y1, x, w, u
){
  psiest<-sigmaest^2*(1-rhoest^2)
  rhostarest<-rhoest*sigmaest

  #parameters in lemma 1
  A<- -crossprod(w,gammaest)
  xi<-rhoest*(y1-crossprod(x,betaest) )/sigmaest
  sigmaestinTN<-sqrt(1-rhoest^2)
  #(xi-A)/sigmaest
  ximAdisigmaest<-(xi-A)/sigmaestinTN
  #denomina part in EY and VY
  denomina<-pnorm(ximAdisigmaest)
  denomina[which(denomina==0),1]<-pnorm(-29)

  #(dnorm((xi-A)/sigmaest))/denomian
  parterm<-( dnorm(ximAdisigmaest) )/denomina

  alpha<-EY<-xi+ parterm*sigmaestinTN
  VY<-( 1-ximAdisigmaest*parterm -parterm^2)*sigmaestinTN^2
  delta<-VY+(EY)^2

  ECMEQfun1<- -1/2*sum(u)*log(2*pi*psiest)
  ECMEQfun21<- (crossprod(u,(y1-crossprod(x,betaest) ) ) )^2
  ECMEQfun22<-  -2*crossprod(u,(rhostarest*(y1-crossprod(x,betaest) )*alpha) )
  ECMEQfun23<- crossprod(u,delta) *rhostarest^2
  ECMEQfun2<-  -1/(2*psiest)*sum(ECMEQfun21,ECMEQfun22,ECMEQfun23 )
  ECMEQfun3<-  -sum(u)*log(2*pi)/2-crossprod(u,delta)/2+crossprod((1-u),pnorm(-t(w)%*%gammaest,log.p=TRUE) )

  Q.val <- sum(ECMEQfun1,ECMEQfun2,ECMEQfun3)
  Q.val
}


#' @method print EMSS
#' @export
print.EMSS <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("\nResponse coefficients:\n")
  print(t(data.frame(Estimate = x$estimate_response,
                     row.names = x$names_response, check.names=FALSE)), digits = digits )
  cat("\nSelection coefficients:\n")
  print(t(data.frame(Estimate = x$estimate_selection,
                     row.names = x$names_selection, check.names=FALSE)), digits = digits )
  cat("\nSigma: ", x$estimate_sigma,"\n")
  cat("\nRho: ", x$estimate_rho)

  cat("\n\n")
}

#' Summarizing EM type Sample Selection Model Fits
#'
#' \code{summary} method for a class "EMSS".
#' @param object an object of class "EMSS" made by the function \code{EMSS}.
#' @param tidy a logical value stands for whether the summary format is in tidy format or not, if \code{TRUE}, the summary function will return a tidy format.
#' @param conf.int a logical value stands for whether the confidence interval is included in the tiny format or not. If \code{TRUE}, confidence intervals are included.
#' If \code{tidy = FALSE}, this parameter does not control anything.
#' @param conf.level a numeric value between 0 and 1 for controlling the significance level of confidence interval; default value is 0.95.
#' @param ... not used, but exists because of the compatibility.
#' @param x an object of class "summary.EMSS".
#' @param digits a numeric number of significant digits.
#' @method summary EMSS
#' @examples
#' # examples continued from EMSS
#' data(Smoke, package = "EMSS")
#' ex1 <- EMSS(response = cigs_intervals ~ educ,
#'            selection = smoker ~ educ + age,
#'            data = Smoke)
#' summary(ex1)
#'
#' @examples
#' data(Smoke, package = "EMSS")
#' ex2 <- EMSS(response = cigs_intervals ~ educ,
#'            selection =  smoker ~ educ + age,
#'            data = Smoke, method="ECMnr")
#' summary(ex2)
#'
#' @examples
#' ## example using random numbers with exclusion restriction
#'
#' N <- 1000
#' errps <- mvtnorm::rmvnorm(N,c(0,0),matrix(c(1,0.5,0.5,1),2,2) )
#' xs <- runif(N)
#' ys <- xs+errps[,1]>0
#' xo <- runif(N)
#' yo <- (xo+errps[,2])*(ys>0)
#'
#' ex3 <- EMSS(response = yo ~ xo,
#'            selection = ys ~ xs,
#'            initial.param = c(rep(0,4), 0.3, 0.6), method="ECMnr")
#' summary(ex3)
#'
#' @export
summary.EMSS <- function(object, tidy = FALSE,
                         conf.int = FALSE,
                         conf.level = 0.95, ...)
{
  p <- object$resp_leng
  q <- object$select_leng
  std.err.mat <- sqrt(diag(vcov.EMSS(object)))
  std_error_response <- std.err.mat[1:p]
  std_error_selection <- std.err.mat[(p+1):(p+q)]
  std_error_sigma <- std.err.mat[p+q+1]
  std_error_rho <- std.err.mat[p+q+2]
  z.value_response <- as.vector(object$estimate_response)/std_error_response
  z.value_selection <- as.vector(object$estimate_selection)/std_error_selection
  z.value_sigma <- as.vector(object$estimate_sigma)/std_error_sigma
  z.value_rho <- as.vector(object$estimate_rho)/std_error_rho

  out <- data.frame( estimate = object$estimate_response, std.error = std_error_response,
                     z.value = z.value_response, pval = 2*pnorm(abs(z.value_response), lower.tail = F) )
  sel <- data.frame( estimate = object$estimate_selection, std.error = std_error_selection,
                     z.value = z.value_selection, pval = 2*pnorm(abs(z.value_selection), lower.tail = F) )
  sigma <- data.frame( estimate = object$estimate_sigma, std.error = std_error_sigma,
                       z.value = z.value_sigma, pval = 2*pnorm(abs(z.value_sigma), lower.tail = F)  )
  rho <- data.frame( estimate = object$estimate_rho, std.error = std_error_rho,
                     z.value = z.value_rho, pval = 2*pnorm(abs(z.value_rho), lower.tail = F)  )
  if(tidy){
    row.names(out) <- paste0(row.names(out), "_outcome")
    row.names(sel) <- paste0(row.names(sel), "_selection")
    tidy_form <- rbind(out, sel, sigma, rho)
    if(conf.int){
      conf <- confint.EMSS(object, level = conf.level)
      confresult <- rbind(conf$response, conf$selection, conf$sigma, conf$rho)
      tidy_form <- cbind(tidy_form, confresult)
    }
  } else {
    tidy_form <- NULL
  }

  result <- list( call = object$call,
                  estimate_response = object$estimate_response,
                  estimate_selection = object$estimate_selection,
                  estimate_sigma = object$estimate_sigma,
                  estimate_rho = object$estimate_rho,
                  Q_value = object$Q_value,
                  std_error_response = std_error_response,
                  std_error_selection = std_error_selection,
                  std_error_sigma = std_error_sigma,
                  std_error_rho = std_error_rho,
                  z.value_response = z.value_response,
                  z.value_selection = z.value_selection,
                  z.value_sigma = z.value_sigma,
                  z.value_rho = z.value_rho,
                  out = out,
                  sel = sel,
                  sigma = sigma,
                  rho = rho,
                  tidy_form = tidy_form
  )

  result$names_response <- object$names_response
  result$names_selection <- object$names_selection
  result$names_sigma <- "sigma"
  result$names_rho <- "rho"

  class(result) <- "summary.EMSS"
  result
}

#' @rdname summary.EMSS
#' @method print summary.EMSS
#' @export
print.summary.EMSS <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  if(is.null(x$tidy_form)){
    cat("\nCall:\n")
    print(x$call)
    cat("\nQ-Value:", x$Q_value, "\n")
    cat("\nResponse equation:\n")
    out <- x$out
    Signif.out <- symnum(out$pval, corr = FALSE, na = FALSE,
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                         symbols = c("***", "**", "*", ".", " "))
    out <- cbind(out, format(Signif.out))
    row.names(out) <- x$names_response ; colnames(out) <- c("Estimate", "Std. Error", "Z Value", "Pr(>|Z|)", "")
    print(out, digits = digits)

    cat("\nSelection equation:\n")
    sel <- x$sel
    Signif.sel <- symnum(sel$pval, corr = FALSE, na = FALSE,
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                         symbols = c("***", "**", "*", ".", " "))
    sel <- cbind(sel, format(Signif.sel))
    row.names(sel) <- x$names_selection ; colnames(sel) <- c("Estimate", "Std. Error", "Z Value", "Pr(>|Z|)", "")
    print(sel, digits = digits)

    cat("---\n")


    cat("\nSigma:\n")
    sigma <- x$sigma
    Signif.sigma <- symnum(sigma$pval, corr = FALSE, na = FALSE,
                           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                           symbols = c("***", "**", "*", ".", " "))
    sigma <- cbind(sigma, format(Signif.sigma))
    row.names(sigma) <- "sigma" ; colnames(sigma) <- c("Estimate", "Std. Error","Z Value", "Pr(>|Z|)", "")
    print(sigma, digits = digits)

    cat("\nRho:\n")
    rho <- x$rho
    Signif.rho <- symnum(rho$pval, corr = FALSE, na = FALSE,
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                         symbols = c("***", "**", "*", ".", " "))
    rho <- cbind(rho, format(Signif.rho))
    row.names(rho) <- "rho" ; colnames(rho) <- c("Estimate", "Std. Error","Z Value", "Pr(>|Z|)", "")
    print(rho, digits = digits)

    cat("---\n")
    cat("Signif. codes:", 0, "'***'", 0.001, "'**'", 0.01, "'*'", 0.05, "'.'", 0.1, "' '", 1, "\n\n")
  } else {
    print(x$tidy_form)
  }
}

#' Getting Coefficients of EM type Sample Selection Model Fits
#'
#' \code{coef} method for a class "EMSS".
#' @param object an object of class "EMSS" made by the function \code{EMSS}.
#' @param only a character value for choosing specific variable's coefficients. Initial value is \code{NULL},
#'  which shows all variable's coefficients. If "response" is written, only coefficients for response variables
#'  will be returned, and if "selection" is written, only coefficients for selection variables will be returned.
#' @param ... not used, but exists because of the compatibility.
#' @return a numeric vector or a list, containing one set or two sets, is given.
#' @method coef EMSS
#' @examples
#' # examples continued from EMSS
#' data(Smoke, package = "EMSS")
#' ex1 <- EMSS(response = cigs_intervals ~ educ,
#'            selection = smoker ~ educ + age,
#'            data = Smoke)
#' coef(ex1)
#'
#' @examples
#' data(Smoke, package = "EMSS")
#' ex2 <- EMSS(response = cigs_intervals ~ educ,
#'            selection =  smoker ~ educ + age,
#'            data = Smoke, method="ECMnr")
#' coef(ex2)
#'
#' @examples
#' ## example using random numbers with exclusion restriction
#'
#' N <- 1000
#' errps <- mvtnorm::rmvnorm(N,c(0,0),matrix(c(1,0.5,0.5,1),2,2) )
#' xs <- runif(N)
#' ys <- xs+errps[,1]>0
#' xo <- runif(N)
#' yo <- (xo+errps[,2])*(ys>0)
#'
#' ex3 <- EMSS(response = yo ~ xo,
#'            selection = ys ~ xs,
#'            initial.param = c(rep(0,4), 0.3, 0.6), method="ECMnr")
#' coef(ex3)
#'
#' @export
coef.EMSS <- function(object, only = NULL, ...)
{
  if (is.null(only)){
    names_response <- object$names_response
    names_selection <- object$names_selection
    result <- list( response = drop(object$estimate_response),
                    selection = drop(object$estimate_selection),
                    sigma = drop(object$estimate_sigma),
                    rho = drop(object$estimate_rho)
    )
    names(result$response) <- names_response
    names(result$selection) <- names_selection
  } else if (only=="response"){
    names_response <- object$names_response
    result <- list( response = drop(object$estimate_response)
    )
    names(result$response) <- names_response
  } else if (only=="selection"){
    names_selection <- object$names_selection
    result <- list( selection <- drop(object$estimate_selection)
    )
    names(result$selection) <- names_selection
  } else {
    stop("'only' has to be defined properly")
  }

  result
}

#' Getting Confidence Intervals for Parameters of EM type Sample Selection Model Fits
#'
#' \code{confint} method for a class "EMSS".
#'
#' @param object an object of class "EMSS" made by the function \code{EMSS}.
#' @param parm not used, but exists because of the compatibility.
#' @param level a numeric value between 0 and 1 for controlling the significance level of confidence interval; default value is 0.95.
#' @param ... not used, but exists because of the compatibility.
#' @method confint EMSS
#' @examples
#' # examples continued from EMSS
#' data(Smoke, package = "EMSS")
#' ex1 <- EMSS(response = cigs_intervals ~ educ,
#'            selection = smoker ~ educ + age,
#'            data = Smoke)
#' confint(ex1)
#'
#' @examples
#' data(Smoke, package = "EMSS")
#' ex2 <- EMSS(response = cigs_intervals ~ educ,
#'            selection =  smoker ~ educ + age,
#'            data = Smoke, method="ECMnr")
#' confint(ex2)
#'
#' @examples
#' ## example using random numbers with exclusion restriction
#'
#' N <- 1000
#' errps <- mvtnorm::rmvnorm(N,c(0,0),matrix(c(1,0.5,0.5,1),2,2) )
#' xs <- runif(N)
#' ys <- xs+errps[,1]>0
#' xo <- runif(N)
#' yo <- (xo+errps[,2])*(ys>0)
#'
#' ex3 <- EMSS(response = yo ~ xo,
#'            selection = ys ~ xs,
#'            initial.param = c(rep(0,4), 0.3, 0.6), method="ECMnr")
#' confint(ex3)
#'
#' @export
confint.EMSS<- function(object, parm, level = 0.95, ...)
{
  coef <- coef.EMSS(object)
  cf_resp <- coef$response
  cf_select <- coef$selection
  pnames_resp <- names(cf_resp)
  pnames_select <- names(cf_select)
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  pct <- paste(as.character(a * 100), "%")
  ci_resp <- array( NA, dim = c(length(cf_resp), 2), dimnames = list(pnames_resp, pct) )
  ci_select <- array( NA, dim = c(length(cf_select), 2), dimnames = list(pnames_select, pct) )
  ci_sigma <- array( NA, dim = c(1, 2), dimnames = list("sigma", pct) )
  ci_rho <- array( NA, dim = c(1, 2), dimnames = list("rho", pct) )
  p <- object$resp_leng
  q <- object$select_leng
  std.err.mat <- sqrt(diag(vcov.EMSS(object)))
  std_error_response <- std.err.mat[1:p]
  std_error_selection <- std.err.mat[(p+1):(p+q)]
  std_error_sigma <- std.err.mat[p+q+1]
  std_error_rho <- std.err.mat[p+q+2]
  ses_resp <- std_error_response
  ses_select <- std_error_selection
  ses_sigma <- std_error_sigma
  ses_rho <- std_error_rho
  ci_resp[] <- cf_resp[pnames_resp] + ses_resp %o% fac
  ci_select[] <- cf_select[pnames_select] + ses_select %o% fac
  ci_sigma[] <- object$estimate_sigma + ses_sigma * fac
  ci_rho[] <- object$estimate_rho + ses_rho * fac
  q <- qchisq(level, 1)

  result <- list(level = level, response = ci_resp, selection = ci_select, sigma = ci_sigma, rho = ci_rho)

  result
}

#' Getting Variance-Covariance Matrix for Parameters of EM type Sample Selection Model Fits
#'
#' \code{vcov} method for a class "EMSS".
#' @param object an object of class "EMSS" made by the function \code{EMSS}.
#' @param ... not used, but exists because of the compatibility.
#' @method vcov EMSS
#' @examples
#' # examples continued from EMSS
#' data(Smoke, package = "EMSS")
#' ex1 <- EMSS(response = cigs_intervals ~ educ,
#'            selection = smoker ~ educ + age,
#'            data = Smoke)
#' vcov(ex1)
#'
#' @examples
#' data(Smoke, package = "EMSS")
#' ex2 <- EMSS(response = cigs_intervals ~ educ,
#'            selection =  smoker ~ educ + age,
#'            data = Smoke, method="ECMnr")
#' vcov(ex2)
#'
#' @examples
#' ## example using random numbers with exclusion restriction
#'
#' N <- 1000
#' errps <- mvtnorm::rmvnorm(N,c(0,0),matrix(c(1,0.5,0.5,1),2,2) )
#' xs <- runif(N)
#' ys <- xs+errps[,1]>0
#' xo <- runif(N)
#' yo <- (xo+errps[,2])*(ys>0)
#'
#' ex3 <- EMSS(response = yo ~ xo,
#'            selection = ys ~ xs,
#'            initial.param = c(rep(0,4), 0.3, 0.6), method="ECMnr")
#' vcov(ex3)
#'
#' @export
vcov.EMSS <- function(object, ...)
{
  result <- solve(  - object$hessian_mat )
  colnames(result) <- c(object$names_response, object$names_selection, "sigma", "rho")
  rownames(result) <- c(object$names_response, object$names_selection, "sigma", "rho")
  result
}

#' Survey Data on Smoking Behaviour
#'
#' The Data is the subset of the original data from Mullahy (1985) and Mullahy (1997). The dataset is from
#' Wooldridge (2009) used for researches on cross sectinal data studies.
#' The dataset is also available from \code{\link[sampleSelection]{Smoke}} from the package \code{sampleSelection}.
#'
#' @usage data(Smoke, package = "EMSS")
#' @format a dataframe with 807 observations and 8 variables as below:
#' \describe{
#' \item{educ}{years of schooling (numeric)}
#' \item{age}{age of respondents (numeric)}
#' \item{cigpric}{cigarette price(state), cents per pack (numeric)}
#' \item{income}{annual income in us dollar (numeric)}
#' \item{restaurn}{state smoking restrictions for restaurants exist or not (categorical)}
#' \item{smoker}{smoked at least once or not (categorical)}
#' \item{cigs_intervals}{number of cigarettes smoked per day, with interval boundaries: 0,5,10,20,50 (numeric)}
#' \item{cigs}{number of cigarettes smoked per day (numeric)}
#' }
#' @source Wooldridge's dataset is available on \url{https://ideas.repec.org/p/boc/bocins/smoke.html#biblio}.
#' @references Jeffrey, M. Wooldridge (2009) \emph{Introductory Econometrics: A modern approach}, Canada: South-Western Cengage Learning.
#' @references Mullahy, John (1985) \emph{Cigarette Smoking: Habits, Health Concerns, and Heterogeneous Unobservables in a Microeconometric Analysis of Consumer Demand},
#' Ph.D. dissertation, University of Virginia.
#' @references Mullahy, John (1997), Instrumental-Variable Estimation of Count Data Models: Applications to Models of Cigarette Smoking Behavior,
#' \emph{Review of Economics and Statistics}, 79, 596-593.
"Smoke"

