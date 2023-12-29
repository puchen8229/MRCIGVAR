#' Data generating process of VAR(p)
#'
#' This function will generate data from a stationary VAR(p) process and return a list containing data and parameters used in the VAR(p) process.
#'
#' @param n     : number of variables
#' @param p     : lag length
#' @param T     : number of observations
#' @param r_np  : an n x p matrix of roots of the characteristic polynomials of n independent AR(p)-processes. If not provided, it will be generated randomly.
#' @param A     : an n x n full rank matrix of transformation to generate correlated VAR(p) from the n independent AR(p)
#' @param B     : (n,n,p) array of the coefficients of the VAR(p) process. If B is not given, it will be calculated out of r_np and A.
#' @param Co    : (n,k+1) matrix of the coefficients of deterministic components in a VAR(p) process. For type="none" Co = 0*(1:n), for type="const" Co is an n vector, for type="exog0" Co is an (n,k) matrix, and for type="exog1" Co is an (n,1+k) matrix.
#' @param U     : residuals, if it is not NA it will be used as input to generate the VAR(p) process.
#' @param Sigma : an n x n covariance matrix of the VAR(p) residuals.
#' @param type  : types of deterministic components. "none", "const" "exog0" and "exog1" are four options
#' @param X     : a (T x k) matrix of exogenous variables.
#' @param mu    : an n vector of the expected mean of the VAR(p) process
#' @param Yo    : a p x n matrix of initial values of the VAR(p) process
#' @return      An object of VAR(p) containing the generated data, the used parameters and the exogenous variables. res = list(n,p,type,r_np,Phi,A,B,Co,Sigma,Y,X,resid,U,Y1,Yo,check)
#' @examples
#' res_d = VARData(n=2,p=2,T=100,type="const")
#' res_d = VARData(n=2,p=2,T=10,Co=c(1:2)*0,type="none")
#' res_d = VARData(n=2,p=2,T=10,Co=c(1:2),  type="const")
#' res_d = VARData(n=3,p=2,T=200,type="exog1",X=matrix(rnorm(400),200,2))
#' res_d = VARData(n=3,p=2,T=200,Co=matrix(c(0,0,0,1,2,3,3,2,1),3,3), type="exog0",
#' X=matrix(rnorm(400),200,2))
#'
#' res_d = VARData(n=2,p=2,T=100,type="const")
#' res_d = VARData(n=3,p=2,T=200,type="exog1",X=matrix(rnorm(400),200,2))
#' @export
VARData = function(n,p,T,r_np,A,B,Co,U,Sigma,type,X,mu,Yo)  {
### T     : number of observations
### n     : number of variables
### p     : lag length
### r_np  : n x p matrix of roots outside unit circle for an n dynamically independent ar-processes (Lag equations).
###         the inverse if the roots will tbe within the unit circle such that the process is stationary.
###         If not provided, it will be generated randomly.
### Sigma : The covariance matrix of the n dynamically independent processes
### A     : an n x n full rank matrix of transformation to generate correlated VAR(p) from the n independent AR(p)
###         (r_np,Sigma,A) if not provided, they will be generated randomly.
### B     : (n,n,p) array of the AR(p) process. If B is not given, it will be calculated out of r_np and A.
### Co    :  (n,k+1) vector of intercept of the VAR(p) process.
###          for "none" Co = O*(1:n), for const Co is a n vector, exog0 Co is a (n,k) matrix for exog1 Co is a (n,1+k) matrix  Depending on type it will be zeros for none,
### type  : deterministic component "none" and "const" "exog0" and "exog1" are four options
### U     : residuals, is not NA will be used to gnerate the data
### X     : (T x k) matrix of exogeneous variables.
### Yo    : (p x n) matrix of initial values of the process
### mu    : n vector of os the expected mean of the VAR(p) process.
### output:
### c("Y","U","Phi","r_np","A","B","Sigma","Y1","check")
###
### Y     : the simulated data via A transformation
### U     : the simulated innovations of VAR(p)
### Phi   : n x p matrix collecting the dynamically independent ar coefficients,
### r_np  : n x p the corresponding roots (in L) of the charachsteristic functions
### A     : the transformation matrix
### B     : the n x n x p arraies of the VAR(p) coefficients
### Sigma : the Covariance matrix of the simulated VAR(p)
### Y1    : the same simulated data calculated via B form
### check : maximum of the data to check stationarity
### type  : can assume "none", "const" "exog0" "exog1"
### X     : exogenous variables
### mu    : mean of the time series
### check : containing stationarity check
### Remarks: the dynamically independent n AR process are stationary. Hence the A transformed VAR(P)
###          is stationary, but the coefficients in B are restricted through the transformation.
###
  if (missing(r_np)) {
      r_np = NA
  }
  if (missing(A)) {
	A = NA
  }
  if (missing(B)) {
	B = NA
  }
  if (missing(Co))   {  Co = NA  }

  if (missing(U))    {  U = NA   }

  if (missing(Sigma)){	Sigma = NA  }
  if (missing(type)) {  type  = NA  }
  if (missing(X))    {   X    = NA  }

  if (missing(mu))   {  mu    = NA  }
  if (missing(Yo))   {  Yo    = NA  }

  if (anyNA(r_np)) {
     r_np = matrix(0,n,p)
     for (i in 1:n)  r_np[i,] <- 0.51/(stats::runif(p)-0.5)  #random number outside unit circle
  }
  d = dim(r_np)
  if (!(n==d[1]&p==d[2])) {print("dimension problem"); return("dimension")}

  Phi = NA
  if (anyNA(Phi))  {
  Phi = matrix(0,n,p)
  for (i in 1:n) Phi[i,] = Roots2coef(p,r_np[i,])
  }

  if (anyNA(A)) {
  A = matrix(stats::runif(n*n),n,n);
  A = A%*%t(A)
  A = eigen(A)[[2]]
  }

  if (anyNA(B)) {
      B = matrix(0,n,n*p); dim(B) = c(n,n,p)
      if (n==1)  {B = Phi; dim(B) = c(n,n,p) }
      if (n >1)  {
		for (i in 1:p) B[,,i]= A%*%diag(Phi[,i])%*%t(A)
  	}
  }


  if (anyNA(Sigma)) {
  	Sigmao = matrix(stats::rnorm(n*n),n,n);
  	Sigmao = Sigmao%*%t(Sigmao)
	Sigma = A%*%Sigmao%*%t(A)
  }  else {
      Sigmao = solve(A)%*%Sigma%*%solve(t(A))
  }


  if (anyNA(U)) {
  	Uh = rnormSIGMA(T,Sigmao)
  	U = Uh%*%t(A)
     	}  else  {
        Uh =  U %*% solve(t(A))
  }


  if (anyNA(Yo)) {
  	Yo = Uh
      } else {
      Uh[1:p,] = Yo
      Yo       = Uh
  }

  for (i in 1:n)          {
    for ( t in (p+1):T )  {
      for (L in 1:p)    Yo[t,i] = Yo[t,i] + Yo[t-L,i]*Phi[i,L]
    }
  }
  Y1 = Yo%*%t(A)
  Ct = U*0



  if (anyNA(type)) {  type = "const" }

  type_soll = Type(Co,X,type)
  if (!type_soll==type)  return(list("type mismatch",type_soll,type))

  if (type=="none")   Co = c(1:n)*0
  if (type=="const") {
  	if (anyNA(mu))     {  mu = matrix(stats::rnorm(n),n,1)}
 	if (anyNA(Co)) {
            	Co = mu
            	for (L in 1:p) Co = Co - B[,,L]%*%mu
      	}   else   {

            H = diag(n)
            for (L in 1:p) H = H - B[,,L]
            mu = solve(H)%*%Co
      }
      Ct = matrix(1,T,1)%*%t(Co)
  }

  if (type=="exog0" & anyNA(Co)) { k = ncol(as.matrix(X)); Co = matrix(stats::rnorm(n*(k+1)),n,k+1);Co[,1]=0}
  if (type=="exog1" & anyNA(Co)) { k = ncol(as.matrix(X)); Co = matrix(stats::rnorm(n*(k+1)),n,k+1)}


  if (!anyNA(X)) {
  	if (type=="exog0" )  Ct = X%*%t(Co[,-1])
  	if (type=="exog1" )  Ct = as.matrix(cbind((1:T)/(1:T),X))%*%t(Co)
  }   else   {
      X = NA
  }

  if (n >0)  {
      Y = U + Ct
	for ( tt in (p+1):T ) {
	for (L in 1:p)  Y[tt,] = Y[tt,]+Y[tt-L,]%*%t(B[,,L])
      }
  }
  check = max(abs(Y1))

  #C = (1:n)*0;
  #if (type == "const") {
  #	Y1 = Y1 + as.matrix((1:T)/(1:T))%*%t(mu)
  #    for (L in 1:p) C = C + B[,,L]%*%mu
  #}
  resid = Uh
  ### to unify the name in MODEL
  result = list(n,p,type,r_np,Phi,A,B,Co,Sigma,Y,X,resid,U,Y1,Yo,check)
  names(result)=c("n","p","type","r_np","Phi","A","B","Co","Sigma","Y","X","resid","U","Y1","Yo","check")
  return(result)
}




#' Check the consistency between type, Co, and EXOG
#'
#' This function will output type according to specification of Co and EXOG and type.
#' @param Co The coefficients of the deterministic components
#' @param EXOG The exogenous variables
#' @param type The type of the deterministic specification for comparison
#'
#' @return type
#' @export
#'
#' @examples
#' Type(Co=matrix(c(0,0,1,1),2,2),EXOG= c(1:10),type="exog0")
#' Type(Co=matrix(c(1,1,1,1),2,2),EXOG= c(1:10),type="exog1")
#' Type(Co=matrix(c(1,1,1,1,0,0),2,3),EXOG= c(1:10),type=NA)
#' Type(Co=c(1,1),EXOG= NA,type=NA)
#' Type(Co=c(0,0),EXOG= NA,type=NA)
#' VARData(n=2,p=2,T=100,Co=matrix(c(2,2,1,1),2,2),type="exog0",X=c(1:100))
#' VARData(n=2,p=2,T=100,Co=matrix(c(1,1,1,1),2,2),type="exog1",X=c(1:100))
Type <- function (Co, EXOG, type)
{
  if (anyNA(Co) & anyNA(EXOG) & anyNA(type))
    type = "const"
  if (anyNA(Co) & !anyNA(EXOG) & anyNA(type))
    type = "exog1"
  if (anyNA(Co) & !anyNA(EXOG) & !anyNA(type))
    type = type
  if (!anyNA(Co) & anyNA(EXOG)) {
    if (sum(Co == 0) == length(Co))
      type = "none"
    else type = "const"
  }
  if (!anyNA(Co) & !anyNA(EXOG)) {
    if (ncol(as.matrix(Co)) == ncol(as.matrix(EXOG)) + 1) {
      if (sum(as.matrix(Co)[, 1] == 0) == length(as.matrix(Co)[, 1]))
        type = "exog0"
      else type = "exog1"
    }
    else type = " Co column<> EXOG column "
  }
  return(type)
}

#' Estimation of VAR(p)
#'
#' This function estimates parameters of a specified VAR(p) model based on provided data.
#'
#' @param  res  :a VAR(p) object which is an output of VARData including at least: n, p, type, Y and optionally X and type.
#' @return res  :a VAR(p) containing estimated parameter values, AIC, BIC and LH.
#' @examples
#' res_d = VARData(n=2,p=2,T=100,type="const")
#' res_e = VARest(res=res_d);
#' IRF_CB = irf_VAR_CB(res=res_e,nstep=20, comb=NA, irf = "gen1", runs = 100, conf = c(0.05, 0.95))
#' IRF_list <-IRF_graph(IRF_CB)
#' res_e$Summary
#'
#' @export
VARest <- function (res)
{     ### add summary
  n = res$n
  p = res$p
  Y = as.matrix(res$Y)
  X = as.matrix(res$X)
  sigma = res$Sigma
  type = res$type
  T = dim(Y)[1]
  Co = res$Co*0
  Z = stats::embed(Y, p + 1)
  if (type == "none")
    Z = Z
  if (type == "const")
    Z = cbind(Z, matrix(1, T - p, 1))
  if (type == "exog0")
    Z = cbind(Z, X[(p + 1):T, ])
  if (type == "exog1")
    Z = cbind(cbind(Z, matrix(1, T - p, 1)), X[(p + 1):T,
                                               ])
  m = dim(Z)[2]
  LREG = stats::lm(Z[, 1:n] ~ 0 + Z[, (n + 1):m])
  bs = as.matrix(LREG$coefficients)[1:(n * p), ]

  if (type == "const") {
    cs = as.matrix(LREG$coefficients)[m - n, ]
    dim(cs) = c(n, 1)
  }
  if (type == "none") {
    cs = as.matrix(LREG$coefficients)[1, ] * 0
    dim(cs) = c(n, 1)
  }
  if (type == "exog1")
    cs = t(as.matrix(LREG$coefficients)[(n * p + 1):(m -
                                                       n), ])
  if (type == "exog0") {
    if (ncol(X) > 1)
      cs = cbind(matrix(0, n, 1), t(as.matrix(LREG$coefficients)[(n *
                                                                    p + 1):(m - n), ]))
    if (ncol(X) == 1)
      cs = cbind(matrix(0, n, 1), as.matrix(LREG$coefficients)[(n *
                                                                  p + 1):(m - n), ])
  }
  sigma[,] = t(LREG$residuals) %*% (LREG$residuals)/(dim(Z)[1] -
                                                       dim(Z)[2])
  resid = Y * 0
  resid[(p + 1):T, ] = LREG$residuals
  bs = t(bs)
  dim(bs) = c(n, n, p)
  B = bs
  attr(cs,"dimnames") <- NULL
  if (type == "exog0")
    Co = (cs)
  if (type == "exog1")
    Co = (cs)
  if (type == "const")
    Co = (cs)
  if (type == "none")
    Co   = (cs)
  res$B <- B
  res$Co <- Co
  res$Sigma <- sigma
  res$resid <- resid
  LH = -(T * n/2) * log(2 * pi) - (T * n/2) + (T/2) * log(det(solve(sigma)))
  AIC = 2 * n * (m - n) + n * (n + 1) - 2 * LH
  BIC = log(T) * (n * (m - n) + n * (n + 1)/2) - 2 * LH
  if (is.null(colnames(Y)))
    colnames(Y) = sprintf("Y%s", 1:n)
  if (is.null(colnames(X)))
    colnames(X) = sprintf("X%s", 1:ncol(X))
  if (n > 1) {
    if (res$type == "none")
      varp = vars::VAR(y = Y, p = p, type = "none")
    if (res$type == "const")
      varp = vars::VAR(y = Y, p = p, type = "const")
    if (res$type == "exog0")
      varp = vars::VAR(y = Y, p = p, type = "none",
                       exogen = X)
    if (res$type == "exog1")
      varp = vars::VAR(y = Y, p = p, type = "const",
                       exogen = X)
  }
  else {
    varp = LREG
  }
  estimation_result = summary(varp)
  Summary = list(estimation_result,LH,AIC,BIC)
  names(Summary) = c("Estimation_Result","LH_function_value","AIC","BIC")
  res$Summary <- Summary
  res$varp    <- varp

  return(res)
}



#' Impulse response function of VAR(p)
#'
#' This function generates impulse response functions of an estimated VAR(p)
#' @param res an object of the output of VARest
#' @param nstep length of the impulse response functions
#' @param comb a vector containing the weights of a combined policy actions
#' @param irf types of impulse response functions
#' @param runs number of simulation runs in generating the confidence interval of the impulse response functions
#' @param conf confidence levels
#' @param smat chock to response transformation matrix for the option irf="smat"
#' @param Xshks the number of exogenous variables
#'
#' @return an (n, n, nstep, 3) array containing the impulse response functions
#' @export
irf_VAR_CB <- function (res, nstep, comb, irf = c("gen", "chol", "chol1", "gen1", "comb1","irfX"), runs = 200, conf = c(0.05, 0.95),smat=NA,Xshks=NA)
{
  n = res$n
  p = res$p
  T = dim(res$Y)[1]
  A = res$A
  B = res$B
  Co = res$Co
  type = res$type
  X = res$X
  crk = res$crk
  neq = dim(B)[1]
  nvar = dim(B)[2]
  sigma = res$Sigma
  response <- array(0, dim = c(neq, nvar, nstep, length(conf) + 1))
  if (irf=="irfX") { Xshk=t(Co[,1+Xshks])}
  response[, , , 1] <- irf_B_sigma(B = B, sigma = sigma, nstep, comb, irf, G=NA, A0=NA, B0=NA, smat, Xshk)
  responseR <- array(0, dim = c(neq, nvar, nstep, runs))
  for (i in 1:runs) {
    Uo_run = rnormSIGMA(T, sigma)
    res_run = VARData(n, p, T, r_np = NA, A = NA, B, Co, U = Uo_run, Sigma = NA, type, X, mu = NA, Yo = NA)
    res_e = VARest(res = res_run)
    B_run = res_e$B
    Co_run = res_e$Co
    sigma_run = res_e$Sigma
    if (irf=="irfX") { Xshk=t(Co_run[,1+Xshks])}
    responseR[, , , i] <- irf_B_sigma(B = res_e$B, sigma = res_e$Sigma, nstep, comb, irf, G=NA, A0=NA, B0=NA, smat, Xshk)
  }
  for (tt in 1:(nstep)) {
    for (i in 1:neq) {
      for (j in 1:nvar) {
        response[i, j, tt, -1] = stats::quantile(responseR[i, j, tt, ], conf)
      }
    }
  }
  return(response)
}





