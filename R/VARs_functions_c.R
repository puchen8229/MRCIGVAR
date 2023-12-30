### Root2coef
### This is function which maps roots (in L) of the characteristic function of an AR process as inout
### to the AR coefficients
###
### input
###
### p:   nummber of lags
### r_p  (optional) A p-vector of roots outside the unit circle
###
### output
### a_p  p-vector of lag-coefficients
###
### the recursion follows from the equation:    (1-a[3,1]L^1-a[3,2]L^2-a[3,3]L^3)(L-r_4)
###                                           = (1-a[4,1]L^1-a[4,2]L^2-a[4,3]L^3-a[4,4]L^4)r_4
###
###
###
#' Roots2Coefficients
#'
#' Given the roots of a characteristic polynom in lags, output the coefficients of the corresponding AR process.
#' @param p The lag length
#' @param r_p A p-vector of roots outside the unit circle
#'
#' @return A vector of AR coefficients
#' @export
#'
#' @examples
#' @keywords internal
#' Roots2coef(3,c(1.1,1.2,1.3))
Roots2coef = function(p,r_p) {
   if (missing(r_p)) r_p <- 0.5/(stats::runif(p)-0.5)  #random number outside unit circle
   if (min(abs(r_p)) < 1) {
   	return("r_p is within the unit circle")
   	#print("r_p is within the unit circle")
   }
   a = matrix(0,p,p)
   a[1,1] = 1/r_p[1]
   if ( p>1 ) {
   for (i in 2:p) {
      for ( j in 1:i ) {
         if (j == 1)           a[i,j] = a[i-1,j] + 1/r_p[i]
         if ((j > 1) & (j<i))  a[i,j] = a[i-1,j] - a[i-1,j-1]/r_p[i]
         if (j == i)           a[i,j] = -a[i-1,i-1]/r_p[i]
      }
   }
   }

   #R = matrix(0,p,p)
   #for (i in 1: p)     {
   #   for (j in 1: p ) {
   #        R[i,j] = r_p[j]^i
   #   }
   #}
   #a[p,]%*%R[,]
   return(a[p,])
}


#' Multivariate normal random series
#'
#' This function will generate iid multivariate normal random time series.
#'
#' @param T Length of the generated time series
#' @param sigma An (n x n) covariance matrix of the normal series
#' @return T x n matrix of iid normal time series
#' @export
#' @keywords internal
rnormSIGMA = function(T,sigma) {
    # generate random numbers from iid multivariate normal distribution with covariance matrix Sigma
    n = dim(sigma)[1]
    U = stats::rnorm(T*n)
    dim(U) = c(T,n)
    U = U%*%chol(sigma)
    return(U)
}



#' Conditional normal random numbers
#'
#' This function generates random numbers from iid multivariate conditional normal distribution with covariance matrix Sigma, given i-th component has the value of v, this will be an (n-1) dimensional random number
#'
#' @param T Length of generated time series
#' @param sigma The (n x n) covariance matrix of the normal series
#' @param I Index of conditioning component
#' @param v The value of the conditioning component
#' @param switch A switch variable: switch = 1 gives the conditional random series and switch = 0 gives the expected values.
#'
#' @return A (T x (n-1)) matrix of iid conditional normal time series or the conditional expected values
#' @export
#' @keywords internal
rnormSIGMA_cond = function(T,sigma,I,v,switch) {
    # generate random numbers from iid multivariate conditional normal distribution with covariance matrix Sigma, given
    # i-th component has the value of v, this will be an (n-1) dimensional random number
      sigma_cond = as.matrix(sigma[-I,-I])-sigma[-I,I]%*%solve(sigma[I,I])%*%sigma[I,-I]
      mu_cond = sigma[-I,I]%*%solve(sigma[I,I])%*%v
      U = rnormSIGMA(T,sigma_cond)*switch+as.vector(c(1:T)/c(1:T))%*%t(mu_cond)
      return(U)
}


#'  This function selects randomly N elements out of a set of T elements
#'
#' @param N The number of elements to be selected
#' @param T The total number of elements
#' @export
#' @keywords internal
NoutofT = function(N,T) {
  unique(round(stats::runif(3*N)*(T-1)))[1:N]+1
}


#' Impulse response function of a vector autoregressive model
#'
#' This function generates impulse response functions for VAR,CIVAR,MRVAR MRCIVAR, also for GVAR, CIGVAR, MRGVAR and MRCIGVAR.
#' For the later four classes of models it also provides the functionalities to calculate the global, regional and country-specific shocks.
#' It also calculates global and regional responses and coordinated policy actions.
#'
#' @param B An (nxnxp) coefficients array of an n-dimensional VAR(p) model
#' @param sigma The covariance matrix of the VAR(p) model
#' @param nstep Number of steps of the impulse response function
#' @param comb An n-vector of weights of the coordinated policy actions
#' @param irf Type of impulse response function
#' @param G The matrix used in the permanent and transitory decomposition
#' @param A0 The matrix for A/B identification in VAR
#' @param B0 The matrix for A/B identification in VAR
#' @param smat An explicit decomposition matrix that defines a structural shock.
#' @param Xshk The shock matrix of a one unit exogenous shock
#'
#' @export
irf_B_sigma <- function (B, sigma, nstep, comb, irf = c("gen", "chol",
                                                        "chol1", "gen1", "genN1", "comb1",
                                                        "smat", "concerts1","PTdecomp","ABSVAR","irfX"), G = NA, A0=NA, B0=NA, smat=NA,Xshk=NA)
{
  neq <- dim(B)[1]
  nvar <- dim(B)[2]
  lags <- dim(B)[3]
  n = dim(sigma)[1]
  if (irf == "irfX") {
    smat = matrix(0,n,n)
    smat[1:nrow(Xshk),] = Xshk
  }

  if (irf == "smat") {
    smat = (smat)
  }
  if (irf == "chol") {
    smat = chol(sigma)
  }
  if (irf == "PTdecomp") {
    smat = t(solve(G)%*%t(chol(G %*% sigma %*% t(G))))
  }
  if (irf == "ABSVAR") {
    smat = t(solve(A0)%*%B0)
  }
  if (irf == "gen") {
    smat = t(sigma %*% diag(1/(sqrt(diag(sigma)))))
  }
  if (irf == "chol1") {
    smat = chol(sigma) %*% diag(1/diag(chol(sigma)))
  }
  if (irf == "gen1") {
    smat = t(sigma %*% diag(1/(diag(sigma))))
  }
  if (irf == "genN1") {
    smat = t(sigma %*% diag(-1/(diag(sigma))))
  }
  if (irf == "concerts1") {
    smat = t((sigma %*% diag(1/(diag(sigma)))) %*% comb)
  }
  if (irf == "concerts0") {
    smat = t((sigma %*% diag(1/(sqrt(diag(sigma))))) %*%
               comb)
  }
  if (irf == "concertc") {
    c = as.numeric(!(comb[, 1] == 0))
    smat = matrix(0, n, n)
    for (i in 1:n) {
      smat[i, ] = sigma[i, ] %*% diag(c) %*% INVI(sigma,
                                                  c, i) %*% comb
    }
    smat = t(smat)
  }
  if (irf == "comb") {
    DD = diag(t(comb) %*% sigma %*% comb)
    for (i in 1:length(DD)) {
      if (DD[i] > 0)
        DD[i] = sqrt(1/DD[i])
    }
    DD = diag(DD)
    smat = t(sigma %*% comb %*% DD)
  }
  if (irf == "comb1") {
    DD = diag(t(comb) %*% sigma %*% comb)
    for (i in 1:length(DD)) {
      if (DD[i] > 0)
        DD[i] = (1/DD[i])
    }
    DD = diag(DD)
    smat = t(sigma %*% comb %*% DD)
  }
  if (dim(smat)[2] != dim(B)[2])
    stop("B and smat conflict on # of variables")
  response <- array(0, dim = c(neq, nvar, nstep))
  response[, , 1] <- t(smat)
  for (it in 2:nstep) {
    for (ilag in 1:min(lags, it - 1)) response[, , it] <- response[, , it] + B[, , ilag] %*% response[, , it - ilag]
  }
  dimnames(response) <- list(dimnames(B)[[2]], dimnames(smat)[[1]], as.character(0:(nstep - 1)))
  return(response)
}





#' Plot impulse response functions
#'
#' @param IRF_CB An (n x n x L x 3) array of impulse response function with confidence bands
#' @param Names An n-vector of strings of the variable names
#' @param INames An n-vector of string of the impulse names
#' @param response An vector of impulse indices
#' @param impulse An vector of response indices
#' @param ncol Number of columns of impulse response functions in the plot
#'
#' @return An ggplot object of impulse response functions
#' @export
#' @keywords internal
IRF_graph <- function(IRF_CB=IRF_CB,Names = NA,INames=NA,response=c(1:n),impulse=c(1:n),ncol=n) {
  ### This function create a list of ggplot objects for the impulse response functions
  IRF_list = list()
  n = dim(IRF_CB)[1]
  if (anyNA(Names))   Names  <-  colnames(IRF_CB)
  if (is.null(Names)) Names  <-  paste0(rep("Y",n),c(1:n))
  if (anyNA(INames))  INames <-  Names

  k = 0
  for (i in response) for (j in impulse) {
    k =  k +1
    myData <- as.data.frame(cbind(c(0:(dim(IRF_CB)[3]-1)),IRF_CB[i,j,,]))
    V1=1;V2=1;V3=1;V4=1
    IRF_list[[k]] <-  ggplot2::ggplot(myData,
                                      ggplot2::aes(x=V1, y=V2, ymin=V3, ymax=V4)) +
      ggplot2::geom_hline(yintercept = 0, color="red") +
      ggplot2::geom_ribbon(fill="grey", alpha=0.5) +
      ggplot2::geom_line() +
      ggplot2::theme_light() +
      ggplot2::ggtitle(paste(Names[i],"response to", INames[j], "shock"))+
      ggplot2::ylab("")+
      ggplot2::xlab("") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust=0.5), axis.title.y = ggplot2::element_text(size=11))
  }

  do.call(gridExtra::grid.arrange,c(IRF_list,ncol=ncol))
  return(IRF_list)
}


#' Embedding a time series
#'
#' Embeds the time series y into a low-dimensional Euclidean space with column names.
#'
#' @param y The time series to be embedded
#' @param p Number of lags for embedding
#' @param prefix Prefix of the column names
#'
#' @return A matrix containing the embedded time series y.
#' @export
#' @keywords internal
Embed <- function(y=tseries::ts(c(1:20)),p=3,prefix="") {
  YY = stats::embed(y,p)
  if (!is.null(colnames(y))) {
    strc = colnames(y)
    strc = paste(prefix,strc,sep="")
    str0 = strc
    if (p > 1) {
      for (i in 1:(p-1)) {
        stri <- paste(str0, ".", sep="")
        stri <- paste(stri,as.character(i),sep="")
        strc <- cbind(strc,stri)
      }
    }
    colnames(YY) <- as.vector(strc)
  }
  return(YY)
}




#' Inverse of a partial covariance matrix
#'
#' @param sigma The input covariance matrix
#' @param c The weighting vector of a concerted policy action
#' @param i The index of responding variables
#'
#' @return A response matrix
#' @export
#' @keywords internal
INVI = function(sigma,c,i) {
  n = dim(sigma)[1]
  sigmaout = matrix(0,n,n)
  cc = as.numeric(c>0)*(1:n)
  sigmai = sigma
  for (j in 1:n)         {
    for (k in 1:n)  {
      if ((j==i)|(k==i)) sigmai[j,k] = 0
    }
  }
  sigmai[i,i] = sigma[i,i]
  invsigmai = solve(sigmai[cc,cc])
  sigmaout[cc,cc] = invsigmai
  return(sigmaout)
}


#' An auxiliary function for nonlinear estimation of the cointegration vectors
#'
#' @param x A vector containing the parameters in the cointegration vectors.
#' @param beta The initial value of the cointegration vectors
#' @param Z1 The lagged I(1) series of the model
#' @param St The indicator series of regime 1
#' @param NSt The indicator series of regime 2
#' @param Y0 The first difference series
#' @param Z2 The lagged difference series
#'
#' @export
#' @keywords internal
f <- function(x,beta,Z1,St,NSt,Y0,Z2) {
  dim(x) = dim(beta)
  CI = Z1%*%x
  CI1 = CI*St
  CI2 = CI*NSt
  residuals <-stats::lm(Y0~0+CI1+CI2+Z2)$residuals
  SSR = sum(diag(t(residuals)%*%residuals))
  return(SSR)
}



#' Help function for nonlinear optimization in constrained estimation of CIVAR
#'
#' @param x Optimization variables
#' @param beta Cointegration vectors
#' @param alpha Adjustment vectors
#' @param G A matrix specifying restrictions on alpha
#' @param H A matrix specifying restrictions on beta
#' @param phi Freely varying parameters in beta
#' @param psi Freely varying parameters in alpha
#' @param h A vector specifying restrictions in beta
#' @param Z1 I(1) data matrix
#' @param St Regime 1 indicator series
#' @param NSt Regime 2 indicator series
#' @param Y0 I(0) data matrix
#' @param Z2 I(0) data matrix
#'
#' @return Sum of squared residuals
#' @export
#' @keywords internal
f_constrained <- function(x, beta=beta,alpha=alpha,G=G,H=H,phi=phi,psi=psi,h=h, Z1, St, NSt, Y0, Z2) {
  ## this is a help function to incoporate restrictions on beta and regime specific alpha
  ## vec(alpha'_1) = G_1 psi_1,  vec(alpha'_2) = G_2 psi_2   , vec(beta) = H phi + h
  ##
  x1 = x[1:length(phi)]
  x2 = x[(1+length(phi)):(length(phi)+length(psi[[1]]))]
  x3 = x[(1+(length(phi)+length(psi[[1]]))):(length(phi)+length(psi[[1]])+length(psi[[2]]))]

  ### restriction on beta
  phi  = x1
  beta_v = H%*% phi + h
  dim(beta_v) = dim(beta)
  #?dim(x) = dim(beta)

  CI = Z1 %*% beta_v
  ### restrictions on alpha_1,alpha_2
  psi_1 = x2
  psi_2 = x3
  alpha_1 = G[[1]]%*%psi_1
  alpha_2 = G[[2]]%*%psi_2
  dim(alpha_1) = dim(alpha)
  dim(alpha_2) = dim(alpha)
  CI1 = (CI * St)%*%t(alpha_1)
  CI2 = (CI * NSt)%*%t(alpha_2)
  residuals <- stats::lm(Y0-CI1-CI2 ~ 0 + Z2)$residuals
  #SSR = sum(diag(t(residuals) %*% residuals))
  SSR <- det(t(residuals) %*% residuals)
  return(SSR)
}


#' Shift of Z2 regressors
#'
#' This function generates impulse response functions of an estimated CIGVAR(m,n,p) with confidence bands.
#'
#' @param  Z2 : Regressor's matrix to be shifted
#' @param  kz : number of common exogenous factors
#' @param  n1 : number of domestic variables
#' @param  p  : the lag length of domestic variables in Z2 which is also the lag length of exogenous common factors
#' @return a shifted Z2 if there are exogenous common factors. It returns Z2 without shift if there are no common factors.
#'
#' @export
#'
#' @keywords internal
ShiftZ2 <- function(Z2,kz,n1,p) {
  # kz : number of Commfakts
  # Z2 : Regressiors on the auxiliery regression
  # n1 : number of domestic variables in a country
  # p  : lag length of domestic regressors. It is the same as the lag length of common factors (p[1,1]-1)
  if ((kz>0)&(p>0)) {
    cc = (1:(kz*p))*0
    for (j in 1:p) {
      for (i in 1:kz) {
        cc[(j-1)*kz+i] = (j-1)*(n1+kz)+i
      }
    }
    Z2 = cbind(Z2[,-cc],Z2[,cc])
  }
  return(Z2)
}



#' Shift of Z2 regressors
#'
#' This function generates impulse response functions of an estimated CIGVAR(m,n,p) with confidence bands.
#'
#' @param  Z2 : Regressor's matrix to be shifted
#' @param  kz : number of common exogenous factors
#' @param  n1 : number of domestic variables
#' @param  p  : the lag length of domestic variables in Z2 which is also the lag length of exogenous common factors
#' @param  p2 : the lag length of foreign  variables in Z2
#' @return a shifted Z2 if there are exogenous common factors. It returns Z2 without shift if there are no common factors.
#'
#' @export
#'
#' @keywords internal
ShiftZ2m <- function(Z2,kz,n1,p,p2) {
  if ((kz>0)&(p>0)) {
    cc = (1:(kz*p))*0
    for (j in 1:p) {
      for (i in 1:kz) {
        cc[(j-1)*kz+i] = (j-1)*(n1+kz)+i
      }
    }
    ccc=c(cc,(p*(n1+kz)+p2*n1)+cc)
    Z2 = cbind(Z2[,-ccc],Z2[,ccc])
  }
  return(Z2)
}


