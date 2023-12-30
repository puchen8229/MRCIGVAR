#' Data generating process of MRCIVAR(n,p,S)
#'
#' This function generates data from an multi-regime cointegrated VAR(p) process with identical cointegration relations and different adjustment speeds and returns a MRCIVAR which is a list containing data and parameters used in the multi regime cointegrated VAR(p) process.
#'
#' @param n     : number of variables
#' @param S     : number of regimes
#' @param p     : an (S x 2) matrix. Each row of p contains the lag length for the corresponding regime and the number of exogenous variables for the regime.
#' @param T     : number of observations
#' @param SESVI : index of the switching variable, switching sv_t = Y\[t-1,SESVI\] > Y\[t-2,SESVI\]
#'
#'                (n,S,p,T,SESVI) are parameters that have to be provided.
#' @param TH    : an (S-1) vector of threshold values
#' @param Bo    : an (n,n,p,S) array of the coefficients  of MRVAR(n,p,S). If not given it will be generated.
#' @param Co    : an (n,k+1,S) array of the coefficients of the deterministic components. For type="none" Co = O*(1:n,1:S), for "const" Co is an n-vector for each regime, "exog0" Co is a (n,k+1,S) array with first column of zeros for each regime respectively,for "exog1" Co is an (n,1+k, S) array without zero restrictions.
#' @param Sigmao : an (n,n,S) array containing S covariance matrices of the residuals
#' @param Uo    : residuals, if it is not NA it will be used as input to generate the MRVAR(n,p,S) for each regime.
#' @param SV    : exogenous switching variable
#' @param type  : type of the deterministic components type = ("none","const","exog0","exog1")
#' @param X     : a (T x k x S) matrix of exogenous variables for each state. The second dimension can be filled with zeros to take into account that the the exogenous variables are not identical in each state.
#' @param mu    : an (n x S) matrix of the regime specific mean of the variables
#' @param Yo    : a (p, n, S) array of initial values of the process
#' @param Do    : a (T, n, S) array of extra exogenous components (not used with value zero)
#' @param d     : lag delay of the self-exiting switching
#' @param r     : number of unit roots in the system. (n - r) is then the cointegration rank.
#'
#'               (TH,Bo,Sigmao,Uo,SV) if not provided, they will be generated randomly.
#' @return An MRCIVAR object containing the generated data, the used parameters, and the exogenous variables.
#'
#' @examples
#' Sigma = 1:(4*4*2)
#' dim(Sigma) = c(4,4,2)
#' Sigma[,,1] = diag(4)
#' Sigma[,,2] = diag(4)
#' p=matrix(0,2,2)
#' p[,1] = c(3,3)
#' res_d = MRCIVARDatam(n=4,p=p,T=250,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="const",r=2)
#' colnames(res_d$Y) = c("w","p","I","Q")
#' max(abs(res_d$Y))
#' res_e = MRCIVARestm1(res=res_d)
#' res_e$Summary
#'
#' p=matrix(0,2,2)
#' p[,1] = c(3,3)
#' p[,2] = 1
#' T = 200
#' XX = matrix(rnorm(T),T,1)
#' XX = cbind(XX,XX); dim(XX) = c(T,1,2)
#' dim(XX)
#'
#' res_d = MRCIVARDatam(n=4,p=p,T=T,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="exog0",r=1,X=XX)
#' res_e = MRCIVARestm1(res=res_d)
#' res_e$Summary
#'
#' res_d = MRCIVARDatam(n=4,p=p,T=T,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="exog1",r=1,X=XX)
#' res_e = MRCIVARestm1(res=res_d)
#' res_e$Summary
#'
#' @export
MRCIVARDatam = function(n=2,p=matrix(2,2,2),T=100,S=2,SESVI,TH,Bo,Co,Sigmao,Uo,SV,type,X,mu,Yo,Do,d=1,r=1) {
### T     : number of observations
### n     : number of variables
### p     : lag length, an S x 2 vector of lag length for each state and the number of exogenous variables for each state
### S     : number of states of the underlying Markov Chain
###         (n,p,T,S,SESVI) are parameters which must be provided.
### SESVI : index of the switching variable in the endogenous variables Y for the case of self-excited threshold model.
### TH    : (S-1)-vector of threshold values
### Bo    : n x n x p x S array collecting the coefficients of VAR(p) in S different states
### Sigmao :  n x n x S array of the covariance matrix of the VAR(p) in S different states
### Uo    : an T x n x S  array of the temporally independent innovation processes
### SV    : exogenous switching variable
### d     : lag of the self-exiting
###         (TH,Bo,Sigmao,Uo,SV) if not provided, they will be generated randomly.
###
### output:
### c("Y","Uo","Bo","Sigmao","TH","St","sv","SESVI","r_npo","check")
###
### Y     : the simulated data via of the MRVAR(p)
### Uo    : T x n x S array of the simulated innovations in S different states
### Bo    : n x n x p x S array collecting the VAR(p) coefficients in S different states
### Sigmao : n x n x S array collecting the covariance matrices of the simulated MRVAR(p) in S different states
### TH    : S-1 vector of threshold values
### St    : simulated time path of states
### sv    : the switching variable ( it is SW in the case of exogeneous switching and it is Y[,SESVI] in the case of self-excited threshold model.
### r_npo : n x p matrix collecting the roots of the characteristic functions in L of the n dynamically independent processes.
### check : maximum of the data to check stationarity
###
### Remarks: The VAR(p) process is A transformed from n dynamically independent ar processes, such that the stationary is guaranteed.
###
###
  if (missing(Yo)) {
    Yo = NA
  }
  if (missing(Do)) {
    Do = NA
  }
  if (missing(d)) {
    d = NA
  }
  if (missing(TH)) {
    TH = NA
  }
  if (missing(Bo)) {
    Bo = NA
  }
  if (missing(Sigmao)) {
    Sigmao = NA
  }
  if (missing(Uo)) {
    Uo = NA
  }
  if (missing(type)) {
    type = NA
  }
  if (missing(Co)) {
    Co = NA
  }
  if (missing(mu)) {
    mu = NA
  }
  if (missing(SV)) {
    SV = NA
  }
  if (missing(X)) {
    X = NA
  }
  if (missing(Yo)) {
    Yo = NA
  }
  if (missing(Do)) {
    Do = NA
  }
  check = c(1:(S + 1)) * 0
  if (anyNA(Yo)) {
    Yo = NA
  }
  if (anyNA(Do)) {
    Do = NA
  }
  if (anyNA(d)) {
    d = 1
  }
  P = max(p[, 1], d)
  Pmax = max(p[, 1])
  Pmin = min(p[, 1])
  if (anyNA(TH)) {
    TH = matrix(stats::runif(S - 1), S - 1, 1)
    TH = TH[order(TH)]
  }
  if (anyNA(Bo)) {
    Bo = (1:(n * n * Pmax * S)) * 0
    dim(Bo) = c(n, n, Pmax, S)
    r_npo = t(matrix(rep((2:(n + 1)), Pmax), Pmax, n))
    for (i in 1:r) r_npo[i, 1] = 1
    res_d = VARData(n, Pmax, T = 100, r_np = r_npo)
    B1 = res_d$B
    CIB1 = B2CIB(B1)[[1]]
    repeat {
      res_d2 = VARData(n, Pmin, T = 100)
      B2 = B1 * 0
      B2[, , 1:Pmin] = res_d2$B
      CIB2 = B2CIB(B2)[[1]]
      values = eigen(CIB1[, , 1])$values
      VECTOR = eigen(CIB1[, , 1])$vectors
      if ((n - r) == 1)
        CIB2[, , 1] = VECTOR[, (r + 1):n] %*% as.matrix(-stats::runif(n -
                                                                 r)) %*% solve(VECTOR)[(r + 1):n, ]
      if ((n - r) > 1)
        CIB2[, , 1] = VECTOR[, (r + 1):n] %*% (diag(-stats::runif(n -
                                                             r))) %*% solve(VECTOR)[(r + 1):n, ]
      CBB = CIB3B(CIB2)
      if (max(abs(STAT(CBB))) == 1)
        break
    }
    if (p[1, 1] > p[2, 1]) {
      Bo[, , , 1] = B1
      Bo[, , , 2] = CBB
    }
    else {
      Bo[, , , 1] = CBB
      Bo[, , , 2] = B1
    }
  }
  else {
    r_npo = NA
  }
  if (anyNA(Sigmao)) {
    Sigmao = (1:(n * n * S)) * 0
    dim(Sigmao) = c(n, n, S)
    for (i in 1:S) {
      VARD = VARData(n, Pmax, T)
      Sigmao[, , i] = VARD$Sigma
    }
  }
  if (anyNA(Uo)) {
    Uo = (1:(T * n * S)) * 0
    dim(Uo) = c(T, n, S)
    for (i in 1:S) Uo[, , i] = rnormSIGMA(T, as.matrix(Sigmao[,
                                                              , i]))
  }
  Ct = Uo * 0
  if (anyNA(type) | type == "none") {
    type = "none"
    Co = (1:(S * n)) * 0
    dim(Co) = c(n, 1, S)
  }
  if (type == "exog0") {
    k = dim(X)[2]
    CC = stats::rnorm(n * (1 + k) * S)
    dim(CC) = c(n, k + 1, S)
    CC[, 1, ] = c(1:n) * 0
    if (anyNA(Co))
      Co = CC
    for (s in 1:S) {
      if (p[s, 2] < k)
        Co[, (p[s, 2] + 2):(k + 1), s] = Co[, (p[s, 2] + 2):(k + 1), s] * 0
      Ct[, , s] = X[, , s] %*% t(Co[, -1, s])
    }
  }
  if (type == "exog1") {
    k = dim(X)[2]
    CC = stats::rnorm(n * (1 + k) * S)
    dim(CC) = c(n, k + 1, S)
    if (anyNA(Co))
      Co = CC
    for (s in 1:S) {
      if (p[s, 2] < k)
        Co[, (p[s, 2] + 2):(k + 1), s] = Co[, (p[s, 2] +
                                                 2):(k + 1), s] * 0
      Ct[, , s] = cbind((1:T)/(1:T), X[, , s]) %*% t(Co[,
                                                        , s])
    }
  }
  if (type == "const") {
    if (anyNA(mu)) {
      mu = matrix(stats::rnorm(n * S), n, S)
      dim(mu) = c(n, 1, S)
    }
    if (anyNA(Co)) {
      Co = mu
      for (s in 1:S) {
        for (L in 1:Pmax) Co[, 1, s] = Co[, 1, s] - Bo[,
                                                       , L, s] %*% mu[, 1, s]
        Ct[, , s] = matrix(1, T, 1) %*% t(Co[, 1, s])
      }
    }
    else {
      mu = Co * NA
      for (s in 1:S) {
        H = diag(n)
        for (L in 1:Pmax) H = H - Bo[, , L, s]
        mu = NA
        Ct[, , s] = matrix(1, T, 1) %*% t(Co[, 1, s])
      }

    }
  }
  St = (1:T) * 0
  Y = as.matrix(Uo[, , 1]) * NA
  if (anyNA(Yo)) {
    Y[1:P, ] = Uo[1:P, , 1]
    Yo = Y[1:P, ]
  }
  else {
    Y[1:P, ] = Yo
  }
  if (anyNA(Do)) {
    Do = matrix(0, n, S)
  }
  else {
    Do = Do
  }
  if (anyNA(SV)) {
    sv = Y[, SESVI]
  }
  else sv = SV
  for (tt in (P + 1):T) {
    if (anyNA(SV)) {
      sv = Y[tt - 1, SESVI] > Y[tt - 2, SESVI]
    }
    else {
      sv = SV[tt]
    }
    st = sv
    s = sum(st > TH) + 1
    St[tt] = s
    Y[tt, ] = Uo[tt, , s] + Ct[tt, , s] + Do[, s]
    for (L in 1:Pmax) Y[tt, ] = Y[tt, ] + Y[tt - L, ] %*%
      t(Bo[, , L, s])
  }
  check[S + 1] = max(abs(Y))
  resid = NA
  crk = n - r
  result = list(Y, X, Uo, resid, Bo, Co, Sigmao, TH, St, sv,
                d, SESVI, r_npo, check, n, p, S, type, Yo, r, crk)
  names(result) = c("Y", "X", "Uo", "resid", "Bo", "Co", "Sigmao", "TH", "St", "sv", "d", "SESVI", "r_npo", "check", "n", "p", "S", "type", "Yo", "r", "crk")
  return(result)
}




#' Estimation of MRCIVAR models
#'
#'This function executes a Johansen test of cointegration ranks and estimates the parameters of a MRCIVAR(n,p,S) model.
#' @param res an MRCIVAR object of the output of MRCIVARDatam
#'
#' @return an estimated MRCIVAR object containing estimated parameters and test results
#' @export
#'
#' @examples
#' Sigma = 1:(4*4*2)
#' dim(Sigma) = c(4,4,2)
#' Sigma[,,1] = diag(4)
#' Sigma[,,2] = diag(4)
#' p=matrix(0,2,2)
#' p[,1] = c(3,3)
#'
#' res_d = MRCIVARDatam(n=4,p=p,T=250,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="const",r=2)
#' colnames(res_d$Y) = c("w","p","I","Q")
#' max(abs(res_d$Y))
#' res_e = MRCIVARestm1(res=res_d)
#' res_e$Summary
#'
MRCIVARestm1 <- function (res)
{
  TH = res$TH
  type = res$type
  p = res$p
  S = res$S
  n = res$n
  Y = res$Y ;  if (max(abs(Y)) > 1000000) return("Check data!!!");
  if (is.null(colnames(Y)))
    colnames(Y) = paste0(rep("Y", ncol(Y)), c(1:ncol(Y)))
  X = res$X
  T = dim(Y)[1]
  SV = res$SV
  sv = res$sv
  d = res$d
  SESVI = res$SESVI
  St = 2 - res$St
  P = max(p[, 1], d)
  Pmax = max(p[, 1])
  Co = res$Co
  r = res$r
  crk = res$crk
  ms = (1:S) * 0
  LH_AIC = 0
  LH_BIC = 0
  LH_P = 0
  LH_N = 0
  b = res$Bo * NA
  cc = res$Co * NA
  sigma = res$Sigmao * NA
  resid = (1:(T * n * S)) * 0
  dim(resid) = c(T, n, S)
  Tresid = resid[, , 1]
  if (type == "none" | type == "exog0")
    Model = "I"
  if (type == "const" | type == "exog1")
    Model = "III"

  if (!anyNA(X)) {
    if (is.null(colnames(X)))
      colnames(X) = paste0(rep("X", ncol(X)), c(1:ncol(X)))
  }


  if (anyNA(X))
    Xi = NA
  else Xi = as.matrix(X[, , 1])

  ORCIVARD = CIVARData(n, p = Pmax, T = T, Co = Co[, , 1],
                       type = type, X = Xi)
  ORCIVARD$Y = Y
  ORCIVAR_e = CIVARest(res = ORCIVARD)
  tst <- MRCVECMestm(y = Y, x = 0, s = St, model = Model, type = "eigen",
                     P = p, crk = n - r, q = 0.95, X = X)
  Tresid[(T - dim(tst$estimation$residuals)[1] + 1):T, ] = tst$estimation$residuals
  Sigma_one = t(Tresid) %*% Tresid/(T - dim(tst$estimation[[1]])[1])
  resid[, , 1] = Tresid * St
  resid[, , 2] = Tresid * (1 - St)
  sigma[, , 1] = t(resid[, , 1]) %*% resid[, , 1]/(sum(St) -
                                                     1 * (n * (p[2, 1] - 1) + crk))
  sigma[, , 2] = t(resid[, , 2]) %*% resid[, , 2]/(sum(1 -
                                                         St) - 1 * (n * (p[2, 1] - 1) + crk))
  CIVAREST = VECM2VARm(param = tst$VECM1[[1]], beta = tst$betaS,
                       p = c(crk, p[1, 1] - 1, p[2, 1] - 1), s = 1)
  Bo = CIVAREST[[1]]
  if (type == "const" | type == "exog1")
    Co = CIVAREST[[3]]
  if (type == "exog0") {
    Co[, 2:(1 + dim(CIVAREST[[3]])[2]), ] = CIVAREST[[3]]
  }
  LH_P = -(T * n/2) * log(2 * pi) - (T * n/2) + (sum(St))/2 *
    log(det(solve(sigma[, , 1]))) + (sum(1 - St))/2 * log(det(solve(sigma[,
                                                                          , 2])))
  LH_P = -(T * n/2) * log(2 * pi) - (T * n/2) + (sum(St))/2 *
    log(det(solve(sigma[, , 1]))) + (sum(1 - St))/2 * log(det(solve(sigma[,
                                                                          , 2])))
  LH_AIC = 2 * (n * (dim(tst$estimation[[1]])[1]) + n * (n +
                                                           1)/2) - 2 * LH_P
  LHH_AIC = 2 * (n * (n * (p[1, 1] - 1 + crk/2) + n * (n +
                                                         1)/2) + n * (n * (p[2, 1] - 1 + crk/2) + n * (n + 1)/2)) -
    2 * LH_P
  LH_BIC = log(T) * (n * (dim(tst$estimation[[1]])[1]) + n *
                       (n + 1)/2) - 2 * LH_P
  LHH_BIC = log(sum(St)) * n * (n * (p[1, 1] - 1 + crk/2) +
                                  n * (n + 1)/2) + log(sum(1 - St)) * n * (n * (p[2, 1] -
                                                                                  1 + crk/2) + n * (n + 1)/2) - 2 * LH_P
  TT = (p[1, 1] - 1) * n + (p[2, 1] - 1) * n + crk - dim(tst$estimation[[1]])[1]
  TTT = dim(tst$estimation[[1]])[1]
  LH_N = 2 * n * (dim(tst$estimation[[1]])[1]) + n * (n + 1)
  LH_TN = log(sum(St)) * n * (n * (p[1, 1] - 1 + crk/2) + n *
                                (n + 1)/2) + log(sum(1 - St)) * n * (n * (p[2, 1] - 1 +
                                                                            crk/2) + n * (n + 1)/2)
  LH_TN1 = log(T) * (n * (dim(tst$estimation[[1]])[1]) + n *
                       (n + 1)/2)
  lm_obj = tst$VECM1
  ORBIC = ORCIVAR_e$Summary$BIC
  ORAIC = ORCIVAR_e$Summary$AIC
  ORLH = ORCIVAR_e$Summary$LH
  ORLHN = ORCIVAR_e$Summary$LH_N
  res$Bo <- Bo
  res$Co <- Co
  res$Sigmao <- sigma
  res$St <- 2 - St
  res$resid <- resid
  est_result = tst$VECM1S
  Summary = list(est_result, tst$erg, LH_AIC, LH_BIC, LH_P,
                 LH_N, LHH_AIC, LHH_BIC, ORBIC, ORAIC, ORLH, ORLHN, LH_TN,
                 LH_TN1, -2 * LH_P)
  names(Summary) = c("Estimation_Result", "CRk_test",
                     "LH_AIC", "LH_BIC", "LH_P", "LH_N",
                     "LHH_AIC", "LHH_BIC", "ORBIC", "ORAIC",
                     "ORLH", "ORLHN", "LH_TN", "LH_TN1",
                     "n2LHP")
  res$Summary = Summary
  res$tst <- tst
  return(res)
}






#' Regime specific impulse response functions of MRCIVAR(n,p,S)
#'
#' This function calculates the regime specific impulse response functions with confidence bands, using Bo\[,,,s\] and Sigma\[,,s\] matrices of the estimated MRCIVAR.
#'
#' @param res_e an object of MRCIVAR as output of MRVARestm
#' @param nstep the length of impulse response function
#' @param irf types of the impulse response function c("gen","chol","chol1","gen1","comb1"), gen for GIRF, gen1 for GIRF with unit impulse, chol Cholezky decomposition, chol1 Cholezky decomposition with unit impulse, comb1 concerted action with a one unit impulse.
#' @param runs Number of simulation runs
#' @param comb a vector specify the concerted action in policy-simulation impulse response function
#' @param G The matrix used in the permanent and transitory decomposition
#' @param smat An explicit decomposition matrix that defines a structural shock.
#' @param conf a vector of the tail probabilities of the confidence interval.
#' @return A list of impulse response function in two regimes and the bootstrap parameters.
#' @examples
#'
#' n =10
#'
#' Sigma = 1:(n*n*2)
#' dim(Sigma) = c(n,n,2)
#' Sigma[,,1] = diag(n)
#' Sigma[,,2] = diag(n)
#' p=matrix(0,2,2)
#' p[,1] = c(3,3)
#' res_d = MRCIVARDatam(n=n,p=p,T=250,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="const",r=2)
#' #colnames(res_d$Y) = c("w","p","I","Q")
#' max(abs(res_d$Y))
#' res_e = MRCIVARestm1(res=res_d)
#' #res_e$Summary
#' if (! max(Mod(STAT(res_e$Bo[,,,1])),Mod(STAT(res_e$Bo[,,,2])) ) > 1.0001 ) {
#'
#'   IRF  = irf_MRCIVAR_CB(res_e,nstep=20,irf="gen1",runs=100,comb=NA,G=NA,conf=c(0.05,0.95))
#'   IRF_g1 <- IRF_graph(IRF[[1]])        # irf of regime 1
#'   IRF_g2<- IRF_graph(IRF[[2]])         # irf of regime 2
#'
#' }
#'
#' @export
irf_MRCIVAR_CB = function (res_e, nstep = 20, irf = c("gen", "gen1"),runs = 100, comb = NA, G = NA, smat=NA, conf = c(0.05, 0.95))
{
  IRF1 = irf_B_sigma(B = res_e$Bo[, , , 1], sigma = res_e$Sigma[,, 1], nstep = nstep, comb = NA, irf = irf, G = NA,smat=NA)
  IRF2 = irf_B_sigma(B = res_e$Bo[, , , 2], sigma = res_e$Sigma[,, 2], nstep = nstep, comb = NA, irf = irf, G = NA,smat=NA)
  n = res_e$n
  p = res_e$p
  T = dim(res_e$Y)[1]
  S = res_e$S
  SV = res_e$SV
  SESVI = res_e$SESVI
  TH = res_e$TH
  Sigmao = res_e$Sigmao
  Bo = res_e$Bo
  Co = res_e$Co
  type = res_e$type
  r = res_e$r
  Uo = res_e$Uo
  Y = res_e$Y
  X = res_e$X
  res_d = MRCIVARDatam(n=n, p=p, T=T, S=S, SESVI=SESVI, TH=TH, Bo=Bo, Co=Co, Sigmao=Sigmao, Uo=Uo, type=type, X=X, r=r)
  neq = res_e$n
  nvar = res_e$n
  response1 <- array(0, dim = c(neq, nvar, nstep, length(conf) +         1))
  response2 <- array(0, dim = c(neq, nvar, nstep, length(conf) +         1))
  response1[, , , 1] <- IRF1
  response2[, , , 1] <- IRF2
  responseR <- array(0, dim = c(neq, nvar, nstep, runs, S))
  BoColect = array(0, c(dim(Bo), runs))
  UoColect = array(0, c(dim(Uo), runs))
  CoColect = array(0, c(dim(Co), runs))
  YColect = array(0, c(dim(Y), runs))
  Uo_run = array(0, c(T, n, S))
  for (i in 1:runs) {
    for (s in 1:S) Uo_run[, , s] = rnormSIGMA(T, Sigmao[, , s])
    if (length(p) > 1) {
      res_run = MRCIVARDatam(n = n, p = p, T = T, S = S,SESVI = SESVI, TH = TH, Bo = Bo, Co = Co, Sigmao = Sigmao,Uo = Uo_run, type = type,X=X,r = r)
      res_e_run = MRCIVARestm1(res = res_run)
    }
    IRF1 = irf_B_sigma(B = res_e_run$Bo[, , , 1], sigma = res_e_run$Sigma[, , 1], nstep = nstep, comb = NA, irf = irf, G = NA)
    IRF2 = irf_B_sigma(B = res_e_run$Bo[, , , 2], sigma = res_e_run$Sigma[, , 2], nstep = nstep, comb = NA, irf = irf, G = NA)
    responseR[, , , i, 1] <- IRF1
    responseR[, , , i, 2] <- IRF2
    BoColect[, , , , i] <- res_e_run$Bo
    CoColect[, , , i] <- res_e_run$Co
    UoColect[, , , i] <- res_run$Uo
    YColect[, , i] <- res_run$Y
  }
  for (tt in 1:(nstep)) {
    for (i in 1:neq) {
      for (j in 1:nvar) {
        response1[i, j, tt, -1] = stats::quantile(responseR[i, j, tt, , 1], conf)
        response2[i, j, tt, -1] = stats::quantile(responseR[i, j, tt, , 2], conf)
      }
    }
  }
  return(list(response1, response2, BoColect, UoColect, YColect))
}



#' Transformation of parameters of MRVECM to paramters of MRCIVAR
#'
#' @param param The coefficient matrix of the estimated multi regime vector error correction model
#' @param beta The estimated cointegration vectors
#' @param p The lags of the MRCIVAR model
#' @param s Types of MRCIVAR models
#'
#' @return The parameter of the multi regime cointegrated VAR
#' @export
#' @keywords internal
#'
#'
VECM2VARm = function (param, beta, p = c(1, 2, 2, 2, 2), s = NA)
{
    m = dim(param)[2]
    VECB = t(param)
    if (anyNA(s)) {                                            ## One Regime
        if (!p[3] == 0) {                                        ## CIGVAR
            B = (1:(m * m * (p[2] + 1))) * 0
            dim(B) = c(m, m, (p[2] + 1))
            A = (1:(m * m * (p[3] + 1))) * 0
            dim(A) = c(m, m, (p[3] + 1))
            AB1 = VECB[, 1:p[1]] %*% t(beta)
            AB2 = VECB[, (1 + p[1]):(p[1] + p[1])] %*% t(beta)
            B[, , 1] = AB[, 1:m]
            A[, , 1] = AB[, (m + 1):(2 * m)]
            for (i in 2:(p[2] + 1)) B[, , i] = VECB[, p[1] *
                2 + ((i - 2) * m + 1):((i - 2) * m + m)]
            for (i in 2:(p[3] + 1)) A[, , i] = VECB[, (p[1] *
                2 + p[2] * m) + ((i - 2) * m + 1):((i - 2) *
                m + m)]
            B = CIB3B(B)
            A = CIA2A(A)
        }
        else {                                                   ##CIVAR
            B = (1:(m * m * (p[2] + 1))) * 0
            dim(B) = c(m, m, (p[2] + 1))
            AB = VECB[, 1:p[1]] %*% t(beta)
            B[, , 1] = AB[, 1:m]
            for (i in 2:(p[2] + 1)) B[, , i] = VECB[, p[1] +
                ((i - 2) * m + 1):((i - 2) * m + m)]
            B = CIB3B(B)
            A = NA
        }
        if (dim(param)[1] > p[1] + (p[2] + p[3]) * m)            ## whehter C is empty
            C = as.matrix(t(param)[, (p[1] + (p[2] + p[3]) * m + 1):dim(param)[1]])
        else C = NA
    }
    else {                                                    ## Multi Regime
        if (length(p) == 5) {                                    ## MRCIGVAR
            P = max(p[-1])
            B = (1:(m * m * (P + 1) * 2)) * 0
            dim(B) = c(m, m, (P + 1), 2)
            A = B
            AB = VECB[, 1:p[1]] %*% t(beta)
            B[, , 1, 1] = AB[, 1:m]
            A[, , 1, 1] = AB[, (m + 1):(2 * m)]
            B[, , 1, 2] = AB[, 1:m]
            A[, , 1, 2] = AB[, (m + 1):(2 * m)]
            for (i in 2:(p[2] + 1)) B[, , i, 1] = VECB[, p[1] +
                ((i - 2) * m + 1):((i - 2) * m + m)]
            for (i in 2:(p[3] + 1)) A[, , i, 1] = VECB[, (p[1] +
                p[2] * m) + ((i - 2) * m + 1):((i - 2) * m +
                m)]
            for (i in 2:(p[4] + 1)) B[, , i, 2] = VECB[, (p[1] +
                (p[2] + p[3]) * m) + ((i - 2) * m + 1):((i -
                2) * m + m)]
            for (i in 2:(p[5] + 1)) A[, , i, 2] = VECB[, (p[1] +
                (p[2] + p[3] + p[4]) * m) + ((i - 2) * m + 1):((i -
                2) * m + m)]
            B[, , 1:(p[2] + 1), 1] = CIB3B(B[, , 1:(p[2] + 1),
                1])
            B[, , 1:(p[4] + 1), 2] = CIB3B(B[, , 1:(p[4] + 1),
                2])
            A[, , 1:(p[3] + 1), 1] = CIA2A(A[, , 1:(p[3] + 1),
                1])
            A[, , 1:(p[5] + 1), 2] = CIA2A(A[, , 1:(p[5] + 1),
                2])
            if (dim(param)[1] > p[1] + (sum(p[-1])) * m) {
                C = (1:(m * 2)) * 0
                dim(C) = c(m, 1, 2)
                C[, 1, 1] = as.matrix(t(param)[, p[1] * 2 + sum(p[-1]) *
                  m + 1])
                C[, 1, 2] = as.matrix(t(param)[, p[1] * 2 + sum(p[-1]) *
                  m + 2])
            }
            else {
                C = NA
            }
        }
        if (length(p) == 3) {                                     ## MRCIVAR
            P = max(p[-1])
            B = (1:(m * m * (P + 1) * 2)) * 0
            dim(B) = c(m, m, (P + 1), 2)
            AB1 = VECB[, 1:p[1]] %*% t(beta)
            AB2 = VECB[, (1 + p[1]):(p[1] + p[1])] %*% t(beta)
            B[, , 1, 1] = AB1[, 1:m]
            B[, , 1, 2] = AB2[, 1:m]
            for (i in 2:(p[2] + 1)) B[, , i, 1] = VECB[, 2 *
                p[1] + ((i - 2) * m + 1):((i - 2) * m + m)]
            for (i in 2:(p[3] + 1)) B[, , i, 2] = VECB[, (2 *
                p[1] + p[2] * m) + ((i - 2) * m + 1):((i - 2) *
                m + m)]
            B[, , 1:(p[2] + 1), 1] = CIB3B(B[, , 1:(p[2] + 1),
                1])
            B[, , 1:(p[3] + 1), 2] = CIB3B(B[, , 1:(p[3] + 1),
                2])
            A = NA
            if (dim(param)[1] > p[1] * 2 + (sum(p[-1])) * m) {		## Whether C is empty  differnt adjustment speed
                km = (dim(param)[1]- (p[1] * 2 + (sum(p[-1])) * m))/2
                C = (1:(m * km* 2)) * 0
                dim(C) = c(m, km, 2)
                C[, 1:km, 1] = as.matrix(t(param)[, (p[1] * 2 + sum(p[-1]) * m + 1):(p[1] * 2 + sum(p[-1]) * m + km)])
                C[, 1:km, 2] = as.matrix(t(param)[,  (p[1] * 2 + sum(p[-1]) * m + km +1):(p[1] * 2 + sum(p[-1]) * m + km+km) ])
            }
            else {
                C = NA
            }
        }
    }
    return(list(B, A, C))
}


#' Test of restrictions in the cointegration space of a MRCIVAR model
#'
#' @param res An MRCIVAR object of the output of MRCIGVARDatam
#' @param H A matrix specifying restrictions on beta
#' @param h A vector specifying restrictions on beta
#' @param phi A vector of freely varying parameters in beta
#' @param G A matrix specifying restrictions on alpha
#' @param psi A vector of freely varying parameters in alpha
#'
#' @return A list containing estimated MRVECM under restrictions, restricted parameter, and likelihood test statistic
#' @export
#'
ABC_MRCIVARestm <- function(res=res,H=H,h=h,phi=phi,G=G,psi=psi) {
  ###
  ### This function runs a likelihood ratio test of linear restrictions on alpha and beta in a CIVAR model
  ###
  ###
  ###        vec(alpha'_1) = G_1[[1]] psi[[1]],  vec(alpha'_2) = G[[2]] psi[[2]]   , vec(beta) = H phi + h
  ###
  ###        example 1 (restrictions on alpha) test of exogeneity
  ###                        vec(alpha) is 45 x 1  vector  ( N = 9 crk = 5, defined by the model )
  ###                  G          is 45 x 40 matrix  ( the first variable is exogeneous, i.e. the 5 adjustment coefficient of the first variable are zero )
  ###                      psi        40 x 1 vector (free variaring parameters not appearing in the specification but implied by G)
  ###                  vec(beta)  is 45 x 1 vector   ( N = 0 crk = 5 )
  ###                  H          is 45 x 45 identity matrix
  ###                  phi        45 x 1 vector (free variaring parameters not appearing in the specification but implied by h =0
  ###                  h          45 x 1 zero matrix implying ver(beta) = phi  >> no restrictions on beta.
  ###                        (H is identity and h is zero vector implies only restrictions on alpha)
  ###
  ###        example 2 (restrictions on beta ) test of PPP
  ###                        vec(alpha) is 40 x 1  vector  ( N = 8 crk = 5 ) conditioal VECM  or VECMX model
  ###                  G          is 40 x 40 identity matrix, implying there is no restriction on alpha
  ###                      psi        40 x 1 vector (free variaring parameters not appearing in the specification but implied by the identity matrix G
  ###                  vec(beta)  is 45 x 1 vector   ( N = 0 crk = 5, defined by the model )
  ###                  H          is 45 x 2 matrix that picks out the elements under restricitons ( two colunms out of the identity matrix ) a zero row in H and the corresponding h implies zero-restrictions on beta.
  ###                             ones in a row of H and zero in the corresponding h implies non-restricted beta.
  ###                  phi        2  x 1  vector (free variaring parameters not appearing in the specification but implied by h =0
  ###                  h          45 x 1  non zero elements in this vector together with the zero elements in the corresponding row in H are the normalization conditions.
  ###                        (H is identity and h is zero vector implies only restrictions on alpha)
  ###
  ###
  ### (y,x,model = c("I","II","III","IV","V"), bnorm = c("1","2"),type = c("eigen", "trace"),p = 1, r = 2, q = 0.95, H=H,h=h,G=G)

  res_e = MRCIVARestm1(res=res)
  #### regime-specific constrained  ML
  betaS = res_e$tst$betaS
  alphaS = res_e$tst$alphaS
  BoR    = res$Bo*0
  CoR    = res$Co*0
  type   = res$type


  S = res$S
  n = res$n
  X = res$X
  T = dim(res$Y)[1]


  Z1    = res_e$tst$Z1
  Z2    = res_e$tst$Z2
  St    = res_e$tst$St
  NSt   = res_e$tst$NSt
  Y0    = res_e$tst$Y0
  Sigma = res_e$Sigmao
  R1    = res_e$tst$R1
  R0    = res_e$tst$R0
  code =  res_e$tst$NLmcode
  p    =  res$p
  crk  =  res$crk
  sigmaR = res$Sigmao * NA
  resid = (1:(T * n * S)) * 0
  dim(resid) = c(T, n, S)
  Tresid = resid[, , 1]

  Omega = 1.0/(nrow(res_e$tst$VECM1$residuals))*t(res_e$tst$VECM1$residuals)%*%res_e$tst$VECM1$residuals

  phi_ini  = as.vector(betaS)[(H%*%phi-h)==1]
  psi_ini1 = as.vector(t(alphaS[[1]]))[!G[[1]]%*%psi[[1]]==0]
  psi_ini2 = as.vector(t(alphaS[[2]]))[!G[[2]]%*%psi[[2]]==0]


  #f_constrained(x=c(phi,psi[[1]],psi[[2]]), beta=betaS,alpha=alphaS[[1]],G=G,H=H,phi=phi,psi=psi,h, Z1, St, NSt, Y0, Z2)
  if (crk > 0) {
    x       = c(phi_ini,psi_ini1,psi_ini2)
    XX      = stats::nlm(f_constrained,x, beta=betaS, alpha=alphaS[[1]], G, H, phi, psi,h, Z1, St, NSt, Y0, Z2,iterlim=300)
    xR      = XX$estimate
    #### not converge for bad initial value but converge for good initial value
    phir    = xR[1:length(phi)]
    psi_1r  = xR[(1+length(phi)):(length(phi)+length(psi[[1]]))]
    psi_2r  = xR[(1+(length(phi)+length(psi[[1]]))):(length(phi)+length(psi[[1]])+length(psi[[1]]))]
    betaR   = H%*%phi + h
    alpha_1 = G[[1]]%*%psi[[1]]
    alpha_2 = G[[2]]%*%psi[[2]]

    dim(betaR) = dim(betaS)
    CI = Z1 %*% betaR
    ### restrictions on alpha_1,alpha_2
    dim(alpha_1) = dim(t(alphaS[[1]]))
    alpha_1 = t(alpha_1)
    dim(alpha_2) = dim(t(alphaS[[1]]))
    alpha_2 = t(alpha_2)

    alphaR  = list(alpha_1,alpha_2)
    CI1 = (CI * St) %*%t(alpha_1)
    CI2 = (CI * NSt)%*%t(alpha_2)
    LM  <- stats::lm(Y0-CI1-CI2 ~ 0 + Z2)
    residuals <- LM$residuals
    CI10 = CI * St
    CI20 = CI * NSt
    VECMR <- stats::lm(Y0 ~ 0 + CI10 + CI20 + Z2)
    VECMR$residuals <-  residuals
    #dim(VECMR$coefficients)
    VECMR$coefficients[1:crk,]  = t(alpha_1)
    VECMR$coefficients[(crk+1):(2*crk),] = t(alpha_2)
    VECMR$coefficients[(2*crk+1):nrow(VECMR$coefficients),] = LM$coefficients
    CIVAREST = VECM2VARm(param = VECMR$coefficients, beta = betaR, p = c(crk, p[1, 1] - 1, p[2, 1] - 1), s = 1)
    BoR = CIVAREST[[1]]
    if (type == "const" | type == "exog1")
      CoR = CIVAREST[[3]]
    if (type == "exog0") {
      CoR[,2:(1+dim(CIVAREST[[3]])[2]),] = CIVAREST[[3]]
    }

    Tresid[(T - dim(VECMR$residuals)[1] + 1):T, ] = VECMR$residuals
    resid[ (T - dim(VECMR$residuals)[1] + 1):T, , 1] = Tresid[ (T - dim(VECMR$residuals)[1] + 1):T, ] * St
    resid[ (T - dim(VECMR$residuals)[1] + 1):T, , 2] = Tresid[ (T - dim(VECMR$residuals)[1] + 1):T, ] * (1 - St)
    sigmaR[, , 1] = t(resid[, , 1]) %*% resid[, , 1]/(sum(St) - 1 * (n * (p[2, 1] - 1) + crk))
    sigmaR[, , 2] = t(resid[, , 2]) %*% resid[, , 2]/(sum(1 - St) - 1 * (n * (p[2, 1] - 1) + crk))
    LR = sum(St)*log(det(sigmaR[,,1])) + sum(NSt)*log(det(sigmaR[,,2])) - sum(St)*log(det(Sigma[,,1])) -  sum(NSt)*log(det(Sigma[,,2]))
    dgf =  n*crk - crk - length(phi)
    p_value = 1 - stats::pchisq(LR,dgf)
    LH_P = -(T * n/2) * log(2 * pi) - (T * n/2) + (sum(St))/2 * log(det(solve(sigmaR[, , 1]))) + (sum(NSt))/2 * log(det(solve(sigmaR[,, 2])))

  }

  SepIni = RZSt2VECM(R0,R1,Y0,Z2,Z1,St,crk)

  #tst = AB_MRCIVARTest(R0,R1,G=G,H=H,h=h,alphaR1=alphaS[[1]],alphaR2=alphaS[[2]],betaR=betaS,alpha1=SepIni$alpha1,alpha2=SepIni$alpha2,beta1=SepIni$beta1,beta2=SepIni$beta2,OmegaR=Omega,Omega1=SepIni$Omega1,Omega2=SepIni$Omega2,IC=1,T1=T,St=St,NSt=NSt)
  tst  = NA
  test = list( VECMR,BoR,CoR,sigmaR,betaR,alphaR,LR,p_value, LH_P, XX$code,tst)
  names(test) = c("VECMR","BoR","CoR","sigmaR","betaR","alphaR","LR","p_value","LH_P","code","tst")
  return(test)
}



#### R1,R0,Z2,St,crk to estimate CIVAR/VECM

#' Calculation of restricted MRCIVAR parameters from data moment matrices
#'
#' @param R0 I(0) residuals of the auxiliary regression
#' @param R1 I(0) residuals of the auxiliary regression
#' @param Y0 The first difference series of the endogenous I(1) variables
#' @param Z2 The first difference series with lags
#' @param Z1 The level series
#' @param St The regime indicator series
#' @param crk The cointegration rank
#'
#' @return A list containing estimated constrained parameters as initial values for further iterations
#' @export
#' @keywords internal
#'
RZSt2VECM <- function(R0,R1,Y0,Z2,Z1,St,crk) {
    #This function estimate separate MRCIVAR
    NSt  <- 1-St
    R0_1 <- R0*St
    R1_1 <- R1*St
    Z2_1 <- (Z2*St)[,which(colSums(Z2*St)!=0)]
    Z1_1 <- Z1*St
    Y0_1 <- Y0*St
    M1 <- ncol(as.matrix(Y0))


	 R0_2 <- R0*NSt
	 R1_2 <- R1*NSt
	 Z2_2 <- (Z2*NSt)[,which(colSums(Z2*NSt)!=0)]
	 Z1_2 <- Z1*NSt
    	 Y0_2 <- Y0*NSt


    S00 <- crossprod(R0_1)/sum(St)
    S01 <- crossprod(R0_1, R1_1)/sum(St)
    S10 <- crossprod(R1_1, R0_1)/sum(St)
    S11 <- crossprod(R1_1)/sum(St)

    Ctemp <- chol(S11, pivot = TRUE)
    pivot <- attr(Ctemp, "pivot")
    oo <- order(pivot)
    C <- t(Ctemp[, oo])
    Cinv <- solve(C)
    S00inv <- solve(S00)
    valeigen <- eigen(Cinv %*% S10 %*% S00inv %*% S01 %*% t(Cinv))
    lambda <- valeigen$values
    e <- valeigen$vector
    V <- t(Cinv) %*% e
    Vorg <- V
    V <- sapply(1:M1, function(j) V[, j]/V[1, j])
    W <- S01 %*% V %*% solve(t(V) %*% S11 %*% V)
    PI <- S01 %*% solve(S11)
    DELTA  <- S00 - S01 %*% V %*% solve(t(V) %*% S11 %*% V) %*% t(V) %*% S10
    #GAMMA <- M02 %*% M22inv - PI %*% M12 %*% M22inv
    beta1  <- as.matrix(V[, 1:crk])
    CI_1   <- Z1_1%*%beta1
    VECM_1 <- stats::lm(Y0_1[which(St!=0),]~0+CI_1[which(St!=0),]+Z2_1[which(St!=0),])
    alpha1 = t(VECM_1$coefficients[1:crk,])
    if (dim(alpha1)[1]==1) alpha1 <-t(alpha1)
    Omega1 =  t(VECM_1$residuals)%*%(VECM_1$residuals)/sum(St)

    VECM_1S <- summaryCIVAR(lm1=VECM_1,sname ="Z2_1")


    S00 <- crossprod(R0_2)/sum(NSt)
    S01 <- crossprod(R0_2, R1_2)/sum(NSt)
    S10 <- crossprod(R1_2, R0_2)/sum(NSt)
    S11 <- crossprod(R1_2)/sum(NSt)

    Ctemp <- chol(S11, pivot = TRUE)
    pivot <- attr(Ctemp, "pivot")
    oo <- order(pivot)
    C <- t(Ctemp[, oo])
    Cinv <- solve(C)
    S00inv <- solve(S00)
    valeigen <- eigen(Cinv %*% S10 %*% S00inv %*% S01 %*% t(Cinv))
    lambda <- valeigen$values
    e <- valeigen$vector
    V <- t(Cinv) %*% e
    Vorg <- V
    V <- sapply(1:M1, function(j) V[, j]/V[1, j])
    W <- S01 %*% V %*% solve(t(V) %*% S11 %*% V)
    PI <- S01 %*% solve(S11)
    DELTA <- S00 - S01 %*% V %*% solve(t(V) %*% S11 %*% V) %*% t(V) %*% S10
    #GAMMA <- M02 %*% M22inv - PI %*% M12 %*% M22inv
    beta2 <- as.matrix(V[, 1:crk])
    CI_2 = Z1_2%*%beta2
    VECM_2 <- stats::lm(Y0_2[which(NSt!=0),]~0+CI_2[which(NSt!=0),]+Z2_2[which(NSt!=0),])
    alpha2 = t(VECM_2$coefficients[1:crk,])
    if (dim(alpha2)[1]==1) alpha2 <-t(alpha2)

    Omega2 =  t(VECM_2$residuals)%*%(VECM_2$residuals)/sum(NSt)

    VECM_2S <- summaryCIVAR(lm1=VECM_2,sname ="Z2_2")

    ret = list(beta1,alpha1,VECM_1,VECM_1S,Omega1,beta2,alpha2,VECM_2,VECM_2S,Omega2)
    names(ret) = c("beta1","alpha1","VECM_1","VECM_1S","Omega1","beta2","alpha2","VECM_2","VECM_2S","Omega2")
    return(ret)
}













#' Calculation of information criteria AIC and BIC for MRCIVAR models
#'
#' @param  res  : an MRCIVAR object generated from MRCIVARDatam or estimated from MRCIVARestm1.
#' @param  L_V  : a 2 components vector containing the maxima of the lags length of the two regimes, respectively.
#' @param  TH_V  : a vector containing possible threshold values for the selection.

#' @return A matrix with different lag specifications, threshold values, and the corresponding model selection criteria.
#' @examples
#'
#' Sigma = 1:(4*4*2)
#' dim(Sigma) = c(4,4,2)
#' Sigma[,,1] = diag(4)
#' Sigma[,,2] = diag(4)
#' p=matrix(0,2,2)
#' p[,1] = c(3,2)
#'
#' res_d = MRCIVARDatam(n=4,p=p,T=2610,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="const",r=1)
#' res_e = MRCIVARestm1(res=res_d)
#'
#' TH_v = c(0,0.1)
#' L_v = c(6,6)
#'
#' Selm = MRCIVAR_Selectm(res=res_e,L_V=L_v,TH_V=TH_v)
#' MRVAR_Select_Summary(Selm)
#'
#' @export
MRCIVAR_Selectm <- function(res=res,L_V=L_V,TH_V=TH_V) {
  res_dd = res
  p      = res_dd$p
  n      = res_dd$n
  T      = dim(res_dd$Y)[1]
  SESVI  = res_dd$SESVI
  type   = res_dd$type
  TH     = res_dd$TH
  X      = res_dd$X
  p[,1]  = max(L_V)
  S      = res_dd$S
  r      = res_dd$r
  #res_dd = MRCIVARData(n=n,p=p,T=T,S=S,SESVI=SESVI,TH=TH,Sigmao=NA,type=type,r=r)

  Criteria = matrix(0,(L_V[1]-1)*(L_V[2]-1)*length(TH_V),9)
  idx = 0

  for (l_d in 2: L_V[1] )             {
    for (l_f in 2:L_V[2] )           {
      for (l_th in 1:length(TH_V))  {
        idx = idx + 1
        res_dd$p[1,1] = l_d
        res_dd$p[2,1] = l_f
        res_dd$TH     = TH_V[l_th]

        #res_dd = MRCIVARDatam(n=n,p=p,T=T,S=S,SESVI=SESVI,TH=TH,Sigmao=NA,type=type,r=r)
        #res_dd$Y= res$Y

        # for MRCIVARestm
        res_s         = MRCIVARestm1(res=res_dd)
        Criteria[idx,] = c (l_d,l_f,TH_V[l_th],res_s$Summary$LH_AIC,res_s$Summary$LH_BIC,res_s$Summary$LHH_AIC,res_s$Summary$LHH_BIC,res_s$Summary$ORAIC,res_s$Summary$ORBIC)
        colnames(Criteria) = c("Lag_regime1","Lag_regime2" ,"threshold","AIC","BIC","LHH_AIC","LHH_BIC","ORAIC","ORBIC")
        # for MECIVARest
        #res_s         = MRCIVARest(res=res_dd)
        #Criteria[idx,] = c (l_d,l_f,TH_V[l_th],res_s$LH_AIC,res_s$LH_BIC,res_s$LH_P,res_s$LH_N,res_s$ORAIC,res_s$ORBIC)
        #colnames(Criteria) = c("Lag_regime1","Lag_regime2" ,"threshold","AIC","BIC","LH_P","LH_N","ORAIC","ORBIC")
      }
    }
  }
  return(Criteria)
}


#' Generalized impulse response functions of MRCIVAR(n,p,S) with regime migrations
#'
#' This function calculates the generalized impulse response functions of an estimated MRVAR(n,p,S).
#'
#' For a given shock vector SHCK:
#'
#' GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK))
#'
#' See H.H. Pesaran and Y. Shin (1998) Generalized impulse response analysis in linear multivariate models, Economics Letters, 58(1) p. 17-29.
#' and G. Koop, M. H. Pesaran, and S. M. Potter (1996), Impulse response analysis in nonlinear multivariate models, Journal of Econometrics, 74 (1996) 119-74.
#'
#' @param  res   : an MRCIVAR object containing the components  of an output of MRCIVARestm1.
#' @param  shock : an n vector containing the shocks as impulse.
#' @param  R     : the number runs to integrate out the random effects in order to obtain the means (see equation above).
#' @param  nstep : the length of the responses
#' @param  Omega_hist : a (P x n) matrix of initial values, from which the impulse response functions start. Omega_hist determines from which regime the impulse response functions start. For Omega_hist=NA, the impulse response functions will start from the most resent observations.
#' @param  resid_method : resid_method = c("resid", "parametric"), It generate the random residuals from residuals bootstrap or parametric bootstrap.
#' @return an (n x n x nstep+1) matrix of impulse response functions. The rows represent response the columns represent impulses.
#' @examples
#'
#' n =4
#'
#' Sigma = 1:(n*n*2)
#' dim(Sigma) = c(n,n,2)
#' Sigma[,,1] = diag(n)
#' Sigma[,,2] = diag(n)
#' p=matrix(0,2,2)
#' p[,1] = c(3,3)
#' res_d = MRCIVARDatam(n=n,p=p,T=250,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="const",r=2)
#' #colnames(res_d$Y) = c("w","p","I","Q")
#' max(abs(res_d$Y))
#' res_e = MRCIVARestm1(res=res_d)
#' res_e$Summary
#' Mod(STAT(res_e$B[,,,1]))
#' Mod(STAT(res_e$B[,,,2]))
#'
#' if (!((max(Mod(STAT(res_e$B[,,,1])))>1)|(max(Mod(STAT(res_e$B[,,,2])))>1)) ) {
#'   girf_MRCIVAR_RM(res=res_e,shock=c(1,1,1,1),R=100,nstep=10,Omega_hist=NA,resid_method="parametric")
#'   GIRF <- girf_MRCIVAR_RM_CB(res=res_e,shock=c(1,1,1,1),R=100,nstep=10,Omega_hist=NA,
#'   resid_method="parametric",conf_level=c(0.05,0.95),N=100)
#'   IRF_g<- IRF_graph(GIRF)
#' }
#' @export
#' @keywords internal
girf_MRCIVAR_RM <- function(res,shock,R,nstep,Omega_hist=NA,resid_method) {
  ####  this function generate the impulse response function of MRVAR with migration
  ####
  ####                  	GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK))
  ####
  ####
  #### res is the output of MRVARest
  n = res$n
  p = res$p
  S = res$S
  SESVI= res$SESVI; if (is.null(SESVI)) SESVI = NA
  TH= res$TH;    if (is.null(TH))       TH = NA
  Bo= res$Bo
  Co= res$Co
  Sigmao= res$Sigmao
  d= res$d
  type= res$type
  P           = max(p,d);
  X= res$X;  if (!anyNA(X))  X = X[1:(P+nstep+1),]
  Yo= res$Yo; if (is.null(Yo)) Yo = NA
  Uo          = res$Uo
  #res = res$res
  vars_obj = res$vars_obj
  TT          = dim(res$Y)[1]

  IRF = list()
  YR          = list()
  YS          = list()
  residR   <-  list()
  residS   <-  residR

  if (anyNA(Omega_hist)) Omega_hist = res$Y[(dim(res$Y)[1]-P+1):dim(res$Y)[1],,drop=FALSE]
  state <- as.numeric(Omega_hist[P-d+1,SESVI]>TH) + 1

  #R = 400
  #shock       = (1:n)/(1:n)
  DIMresid     = dim(Uo)
  if (length(DIMresid) == 3)  residI = Uo[1:(P+1+nstep),,]*0
  if (length(DIMresid) == 2)  residI = Uo[1:(P+1+nstep),]*0
  shock_mat = Sigmao[,,state]%*%solve(diag(diag(Sigmao[,,state])))%*%diag(shock)

  residI    = array(0,c(P+1+nstep,n,n,DIMresid[3]))

  Yo = Omega_hist
  YR = array(0,c(P+1+nstep,n,n))
  YS <- YR
  MYS <- YS
  MYR <- YS
  GIRF <- array(0,c(n,n,nstep+1))

  for (i in 1:R) {
    if ( resid_method=="resid" ) {
      if (length(DIMresid) == 2)  {
        residI  = Uo[NoutofT(P+nstep+1,DIMresid[1]),];
      }
      if (length(DIMresid) == 3)  {
        for (k in 1:n) {
          residI[,,k,]  = Uo[NoutofT(P+nstep+1,DIMresid[1]),,];
        }
      }
    }
    if ( resid_method=="parametric" ) {
      if (length(DIMresid) == 2)  {
        residI   = rnormSIGMA(P+nstep+1,Sigmao)
        residI[1:P,]= residI[1:P,]*0
      }
      if (length(DIMresid) == 3)  {
        for (k in 1:n) {
          for (s in 1:S)   residI[,,k,s] = rnormSIGMA(P+nstep+1,Sigmao[,,s])
        }
      }
    }

    residR[[i]] <- residI
    for (s in 1:S) residI[P+1,,,s] <-shock_mat
    residS[[i]] = residI

    for (k in 1:n) {
      YR[,,k] = MRCIVARDatam(n=n,p=p,T=(P+nstep+1),S=S,SESVI=SESVI,TH=TH,Bo=Bo,Co=Co,Sigmao=Sigmao,Uo=residR[[i]][,,k,],type=type,X=X,Yo=Yo,d=d)$Y
      YS[,,k] = MRCIVARDatam(n=n,p=p,T=(P+nstep+1),S=S,SESVI=SESVI,TH=TH,Bo=Bo,Co=Co,Sigmao=Sigmao,Uo=residS[[i]][,,k,],type=type,X=X,Yo=Yo,d=d)$Y
    }

    MYR = MYR + 1/R*YR
    MYS = MYS + 1/R*YS

  }  # end of R loop
  ############## integrating out the disturbances
  for (i in 1:n ) {
    for (j in 1:n) {
      for (k in 1:(nstep+1)) GIRF[i,j,k] = MYS[P+k,i,j]-MYR[P+k,i,j]
    }
  }
  return(GIRF)
}




#' Generalized impulse response functions of MRCIVAR(n,p,S) with regime migrations with a confidence interval
#'
#' This function calculates the generalized impulse response functions of an estimated MRVAR(n,p,S) with a bootstrapping confidence interval.
#'
#' For a given shock vector SHCK:
#'
#' GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK))
#'
#' See H.H. Pesaran and Y. Shin (1998) Generalized impulse response analysis in linear multivariate models, Economics Letters, 58(1) p. 17-29.
#' and G. Koop, M. H. Pesaran, and S. M. Potter (1996), Impulse response analysis in nonlinear multivariate models, Journal of Econometrics, 74 (1996) 119-74.
#'
#' @param  res   : an MRCIVAR object containing the components  of an output of MRCIVARestm1.
#' @param  shock : an n vector containing the shocks as impulse.
#' @param  R     : the number runs to integrate out the random effects in order to obtain the means (see equation above).
#' @param  nstep : the length of the responses
#' @param  Omega_hist : a (P x n) matrix of initial values, from which the impulse response functions start. Omega_hist determines from which regime the impulse response functions start. For Omega_hist=NA, the impulse response functions will start from the most resent observations.
#' @param  resid_method : resid_method = c("resid", "parametric"), It generate the random residuals from residuals bootstrap or parametric bootstrap.
#' @param  conf_level : a vecter contain the level of confidences
#' @param  N            : the number of bootstrap runs in containing the bootstrapped confidence intervals.
#' @return an (n, n, nstep+1,3) array containing the impulse response functions with lower and upper confidence bonds. The rows represent response, and the columns represent impulses.
#' @examples
#' n =4
#'
#' Sigma = 1:(n*n*2)
#' dim(Sigma) = c(n,n,2)
#' Sigma[,,1] = diag(n)
#' Sigma[,,2] = diag(n)
#' p=matrix(0,2,2)
#' p[,1] = c(3,3)
#' res_d = MRCIVARDatam(n=n,p=p,T=250,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="const",r=2)
#' #colnames(res_d$Y) = c("w","p","I","Q")
#' max(abs(res_d$Y))
#' res_e = MRCIVARestm1(res=res_d)
#' res_e$Summary
#' Mod(STAT(res_e$B[,,,1]))
#' Mod(STAT(res_e$B[,,,2]))
#'
#' #res_e$Summary
#' if (! max(Mod(STAT(res_e$Bo[,,,1])),Mod(STAT(res_e$Bo[,,,2])) ) > 1.0001 ) {
#'
#'  IRF  = irf_MRCIVAR_CB(res_e,nstep=20,irf="gen1",runs=10,comb=NA,G=NA,conf=c(0.05,0.95))
#'  IRF_g1 <- IRF_graph(IRF[[1]])        # irf of regime 1
#'  IRF_g2<- IRF_graph(IRF[[2]])         # irf of regime 2
#'
#' }
#'
#' @export
girf_MRCIVAR_RM_CB <- function (res, shock, R, nstep, Omega_hist=NA, resid_method = "parametric", conf_level, N)
{
  n = res$n
  p = res$p
  S = res$S
  SESVI = res$SESVI
  TH = res$TH
  Bo = res$Bo
  Co = res$Co
  Sigmao = res$Sigmao
  d = res$d
  type = res$type
  X = res$X
  d = res$d
  T = dim(res$Y)[1]
  P = max(p, d)
  GIRF = (1:(n * n * (nstep + 1) * N)) * 0
  dim(GIRF) = c(n, n, nstep + 1, N)
  GIRFBd = (1:(n * n * (nstep + 1) * (length(conf_level) + 1))) * 0
  dim(GIRFBd) = c(n, n, nstep + 1, length(conf_level) + 1)
  GIRFBd[, , , 1] = girf_MRCIVAR_RM(res, shock = shock, R = 2*R, nstep, Omega_hist=Omega_hist, resid_method)
  for (i in 1:N) {
    res_d = MRCIVARDatam(n = n, p = p, T = T, S = S, SESVI = SESVI, TH = TH, Bo = Bo, Co = Co, Sigmao, type = type, X = X, d = d)
    if (length(colnames(res_d$Y))==0) colnames(res_d$Y) = paste("Y",1:ncol(res_d$Y),sep="")
    RESS = MRCIVARestm1(res_d)
    RF3 = girf_MRCIVAR_RM(res = RESS, shock = shock, R = R, nstep, Omega_hist=Omega_hist, resid_method)
    GIRF[, , , i] = RF3
  }
  GIRF[,,,1] = GIRFBd[, , , 1]
  for (tt in 1:( nstep + 1)) {
    for (i in 1:n) {
      for (j in 1:n) {
        GIRFBd[i, j, tt, -1] = stats::quantile(GIRF[i, j, tt, ], conf_level)
      }
    }
  }
  return(GIRFBd)
}









