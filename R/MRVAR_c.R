#' Data generating process of MRVAR(n,p,S)
#'
#' This function generates data from an multi-regime stationary VAR(p) process and returns an MRVAR(n,p,S) object that is a list containing data and the parameters used in the MRVAR(n,p,S) process.
#'
#' @param n     : number of variables
#' @param S     : number of regimes
#' @param p     : an (S x 2) matrix. Each row of p specifies the lag length and the number of exogenous variables for the corresponding regime.
#' @param T     : number of observations
#' @param SESVI : index of the switching variable
#'
#'                (n,S,p,T,SESVI) must be provided.
#' @param TH    : an (S-1)-vector of threshold values
#' @param Bo    : an (n,n,p,S) array of the coefficients of the MRVAR(n,p,S) model. If Bo is not given it will be generated.
#' @param Co    : an (n,k+1,S) array of the coefficients of the deterministic components. For type="none" Co = O*(1:n,1:S), for "const" Co is an (1:n,1:S) array, for "exog0" Co is an (n,k+1,S) array with first column zero for each regime respectively,for "exog1" Co is an (n,1+k, S) array without zero restrictions.
#' @param Sigmao : an (n,n,S)   array containing S covariance matrices of residuals, each for one regime.
#' @param Uo    : residuals, if it is not NA it will be used as input to generate data of the MRVAR(n,p,S) process.
#' @param SV    : exogenous switching variable
#' @param type  : type of the deterministic components type = c("none","const","exog0","exog1")
#' @param X     : a (T x k x S) array of exogenous variables for each state. The second dimension can be filled with zeros to take into account that the the exogenous variables are not identical in each state.
#' @param mu    : an (n x S) matrix of the regime specific means of the variables
#' @param Yo    : a (p, n, S) array of initial values of the process
#' @param Do    : a (T, n, S) array of extra exogenous components (not used with value zero)
#' @param d     : lag delays of the self-exiting switching variable.
#'
#'              (TH,Bo,Co,Sigmao,Uo,SV,type) if not provided, they will be generated randomly.
#' @return      an MRVAR(n,p,S) object that is a list containing the generated data, the used parameters and the exogenous variables.
#'
#' @examples
#' p = matrix(c(3,3,0,0),2,2)
#' res_d = MRVARData(n=2,p=p,T=200,S=2,SESVI=1,type="const")
#' colnames(res_d$Y) = c("Y","P")
#' STAT(res_d$B[,,,1])
#' STAT(res_d$B[,,,2])
#' res_e = MRVARest(res=res_d)
#' res_e$Summary
#' @export
MRVARData = function(n,p,T,S,SESVI,TH,Bo,Co,Sigmao,Uo,SV,type,X,mu,Yo,Do,d) {
  ### T     : number of observations
  ### n     : number of variables
  ### p     : lag length, an S x 2 vector of lag length for each state and the number of exogenous variables for each state
  ### S     : number of states of the underlying Markov Chain
  ###         (n,p,T,S,SESVI) are parameters which must be provided.
  ### SESVI : index of the switching variable in the endogeneous variables Y for the case of self-excited threshhold model.
  ### TH    : (S-1)-vector of threshold values
  ### Bo    : n x n x p x S array collecting the coefficients of VAR(p) in S different states
  ### Sigmao: n x n x S array of the covariance matrix of the VAR(p) in S different states
  ### Uo    : an T x n x S  array of the temorally independent innovation processes
  ### SV    : exogeneous switching variable
  ### d     : lag of the self-exiting
  ###         (TH,Bo,Sigmao,Uo,SV) if not provided, they will be generated randomly.
  ###
  ### output:
  ### c("Y","Uo","Bo","Sigmao","TH","St","sv","SESVI","r_npo","check")
  ###
  ### Y     : the simulated data via of the MRVAR(p)
  ### Uo    : T x n x S array of the simulated innovations in S different states
  ### Bo    : n x n x p x S array collecting the VAR(p) coefficients in S different states
  ### Sigmao: n x n x S array collecting the covariance matrices of the simulated MRVAR(p) in S different states
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
  if (anyNA(TH)) {
    TH = matrix(stats::runif(S - 1), S - 1, 1)
    TH = TH[order(TH)]
  }




    if (anyNA(Bo)) {
      Bo = (1:(n * n * Pmax * S)) * 0
      dim(Bo) = c(n, n, Pmax, S)
      r_npo = c(1:(n * Pmax * S)) * 0
      dim(r_npo) = c(n, Pmax, S)
      for (i in 1:S) {
        VARD = VARData(n, p[i], T)
        Bo[, , 1:p[i, 1], i] = VARD$B
        r_npo[, 1:p[i, 1], i] = VARD$r_np
        check[i] = max(abs(VARD$Y))
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
  if (anyNA(type)|type=="none") {
    type = "none"
    Co = matrix(0,n,2)
    dim(Co) = c(n,1,2)
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
        Co[, (p[s, 2] + 2):(k + 1), s] = Co[, (p[s, 2] +
                                                 2):(k + 1), s] * 0
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
        for (L in 1:Pmax) Co[, 1, s] = Co[, 1, s] - Bo[, , L, s] %*% mu[, 1, s]
        Ct[, , s] = matrix(1, T, 1) %*% t(Co[, 1, s])
      }
    }
    else {
      mu = Co * NA
      for (s in 1:S) {
        H = diag(n)
        for (L in 1:Pmax) H = H - Bo[, , L, s]
        mu[, 1, s] = solve(H) %*% Co[, 1, s]
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
      sv = Y[, SESVI]
    }
    else {
      sv = SV
    }
    st = sv[tt - d]
    s = sum(st > TH) + 1
    St[tt] = s
    Y[tt, ] = Uo[tt, , s] + Ct[tt, , s] + Do[, s]
    for (L in 1:Pmax) Y[tt, ] = Y[tt, ] + Y[tt - L, ] %*%
      t(Bo[, , L, s])
  }
  check[S + 1] = max(abs(Y))
  resid = NA
  result = list(Y, X, Uo, resid, Bo, Co, Sigmao, TH, St, sv, d, SESVI, r_npo, check, n, p, S, type, Yo)
  names(result) = c("Y", "X", "Uo", "resid", "Bo", "Co", "Sigmao", "TH", "St","sv", "d", "SESVI", "r_npo","check", "n", "p", "S", "type","Yo")
  return(result)
}

#' Estimation of MRVAR(n,p,S)
#'
#' This function estimates the parameters of a specified MRVAR(n,p,S) based on provided data.
#'
#' @param  res  :an object of MRVAR that is an output of MRVARData including at least: n, S, p, type, Y, SESVI, TH, d, and optionally X.
#' @return an MRVAR object that is a list containing estimated parameters and some test statistics.
#'
#' @examples
#' p = matrix(c(3,3,0,0),2,2)
#' res_d = MRVARData(n=2,p=p,T=200,S=2,SESVI=1,type="const")
#' colnames(res_d$Y) = c("Y","P")
#' STAT(res_d$B[,,,1])
#' STAT(res_d$B[,,,2])
#' res_e = MRVARest(res=res_d)
#' res_e$Summary
#'
#' @export
MRVARest <- function (res)
{
  TH = res$TH
  type = res$type
  p = res$p
  S = res$S
  n = res$n
  Y = res$Y;  if (max(abs(Y)) > 1000000) return("Check data!!!");
  X = res$X
  T = dim(Y)[1]
  SV = res$SV
  sv = res$sv
  d = res$d
  SESVI = res$SESVI
  P = max(p[, 1], d)
  Pmax = max(p[, 1])
  ms = (1:S) * 0
  LH_AIC = 0
  LH_BIC = 0
  LLH_AIC = 0
  LLH_BIC = 0
  LH_P = 0
  LH_N = 0
  LH_NN = 0
  b = res$Bo * 0
  cc = res$Co * 0
  sigma = res$Sigmao * 0
  resid = (1:(T * n * S)) * 0
  dim(resid) = c(T, n, S)
  Z = stats::embed(Y, Pmax + 1)
  if (!anyNA(X)) {
    XX = JointX(X)
    if (type == "none") {
      ZO = stats::embed(Y, Pmax + 1)
      cc = matrix(0, n, S)
      dim(cc) = c(n, 1, S)
    }
    if (type == "const")
      ZO = cbind(Z, (1:(T - Pmax))/(1:(T - Pmax)))
    if (type == "exog0")
      ZO = cbind(Z, XX[(Pmax + 1):T, ])
    if (type == "exog1")
      ZO = cbind(Z, ((Pmax + 1):T)/((Pmax + 1):T), XX[(Pmax +
                                                         1):T, ])
    mso = dim(ZO)[2]
    LREGOR = stats::lm(ZO[, 1:n] ~ 0 + ZO[, (n + 1):mso])
  }
  if (!anyNA(X)) {
    XX = JointX(X)
    Co = matrix(1, n, dim(XX)[2] + 1)
    if (type == "exog0")
      Co[, 1] = 0
    VARDataHelp = VARData(n = n, p = Pmax, T = T, type = type,
                          Co = Co, X = XX)
  }
  else {
    if (type == "none")
      Co = matrix(0, n, 1)
    if (type == "const")
      Co = matrix(1, n, 1)
    VARDataHelp = VARData(n = n, p = Pmax, T = T, type = type,
                          Co = Co)
  }
  VARDataHelp$Y = res$Y
  ORBIC = VARest(res = VARDataHelp)$Summary$BIC
  ORAIC = VARest(VARDataHelp)$Summary$AIC
  m <- ncol(Z)
  mm = c(1:S) * 0
  vars_obj = list()
  St = (1:T) * 0
  if (is.null(SV)) {
    sv = Y[, SESVI]
  }
  else sv = SV
  for (tt in (P + 1):T) {
    svt = sv[tt - d]
    s = which.min(abs(TH - svt))
    if (svt > TH[s]) {
      s = s + 1
    }
    else {
      s = s
    }
    St[tt] = s
  }
  select = matrix(0, T - P, S)
  select2 = select
  for (s in 1:S) {
    select[, s] = (1:(T - P)) * (St[(P + 1):T] == s)
    select2[, s] = ((P + 1):T) * (St[(P + 1):T] == s)
    Z = stats::embed(Y, P + 1)[, 1:((p[s, 1] + 1) * n)]
    if (type == "none") {
      Z = Z
      cc = matrix(0, n, S)
      dim(cc) = c(n, 1, S)
    }
    if (type == "const")
      Z = cbind(Z, (1:(T - P))/(1:(T - P)))
    if (type == "exog0")
      Z = cbind(Z, X[(P + 1):T, 1:p[s, 2], s])
    if (type == "exog1")
      Z = cbind(Z, ((P + 1):T)/((P + 1):T), X[(P + 1):T,
                                              1:p[s, 2], s])
    ms[s] = dim(Z)[2]
    if ((dim(as.matrix(Z[select[, s], 1:n]))[1] - dim(as.matrix(Z[select[, s], (n + 1):ms[s]]))[2] - n > 0)) {
      LREG = stats::lm(Z[select[, s], 1:n] ~ 0 + Z[select[, s], (n + 1):ms[s]])
      bs = as.matrix(LREG$coefficients)[1:(n * p[s, 1]), ]
      if (type == "none")
        cs = as.matrix(LREG$coefficients)[1, ] * 0
      if (type == "const")
        cs = as.matrix(LREG$coefficients)[ms[s] - n,]
      if (type == "exog0")
        cs = as.matrix(LREG$coefficients)[(n * p[s,1] + 1):(ms[s] - n), ]
      if (type == "exog1")
        cs = as.matrix(LREG$coefficients)[(n * p[s,1] + 1):(ms[s] - n), ]
      sigmas = t(LREG$residuals) %*% (LREG$residuals)/(dim(as.matrix(Z[select[,s], 1:n]))[1] - dim(as.matrix(Z[select[, s],(n + 1):ms[s]]))[2])
      bs = t(bs)
      dim(bs) = c(n, n, p[s, 1])
      b[, , 1:p[s, 1], s] = bs
      if (type == "none")
        cc[, 1, s] = (1:n) * 0
      if (type == "const")
        cc[, 1, s] = cs
      if (type == "exog0") {
        cc[, c(2:(1 + p[s, 2])), s] = t(cs)
        cc[, 1, s] = (1:n) * 0
      }
      if (type == "exog1")
        cc[, c(1:(1 + p[s, 2])), s] = t(cs)
      sigma[, , s] = sigmas
      resid[select2[, s], , s] = LREG$residuals
    }  else return("check data for regime")

    if (!is.na(log(det(as.matrix(sigma[, , s])))))
      LH_AIC = LH_AIC + 2 * n * (dim(Z[, (n + 1):ms[s]])[2] +
                                   (n + 1)/2) + dim(as.matrix(Z[select[, s], 1:n]))[1] *
      log(det(as.matrix(sigma[, , s]))) + dim(as.matrix(Z[select[,
                                                                 s], 1:n]))[1] * n * (1 + log(2 * pi))
    if (!is.na(log(det(as.matrix(sigma[, , s])))))
      LH_BIC = LH_BIC + log(dim(as.matrix(Z[select[, s],
                                            1:n]))[1]) * n * (dim(Z[, (n + 1):ms[s]])[2] +
                                                                (n + 1)/2) + dim(as.matrix(Z[select[, s], 1:n]))[1] *
      log(det(as.matrix(sigma[, , s]))) + dim(as.matrix(Z[select[,
                                                                 s], 1:n]))[1] * n * (1 + log(2 * pi))
    if (!is.na(log(det(as.matrix(sigma[, , s])))))
      LH_P = LH_P + dim(as.matrix(Z[select[, s], 1:n]))[1] *
      log(det(as.matrix(sigma[, , s]))) + dim(as.matrix(Z[select[,
                                                                 s], 1:n]))[1] * n * (1 + log(2 * pi))
    if (!is.na(log(det(as.matrix(sigma[, , s])))))
      LH_N = LH_N + 2 * n * (dim(Z[, (n + 1):ms[s]])[2] +
                               (n + 1)/2)
    if (!is.na(log(det(as.matrix(sigma[, , s])))))
      LH_NN = LH_NN + log(dim(as.matrix(Z[select[, s],
                                          1:n]))[1]) * (n * (dim(Z[, (n + 1):ms[s]])[2] +
                                                               (n + 1)/2))
    if (!is.na(log(det(as.matrix(sigma[, , s])))))
      LLH_AIC = LH_P + LH_N
    if (!is.na(log(det(as.matrix(sigma[, , s])))))
      LLH_BIC = LH_P + LH_NN
    Ys = as.matrix(Z[select[, s], 1:n])
    Zs = as.matrix(Z[select[, s], (n + 1):ms[s]])
    YY = Y[1:(nrow(as.matrix(Z[select[, s], 1:n])) + p[s,1]), ]
    if ((res$type == "exog0") | (res$type == "exog1")) {
      XX = X[1:(nrow(Z[select[, s], 1:n]) + p[s, 1]), 1:p[s,2], s]
      if (length(colnames(XX))==0) colnames(XX) = paste("XX",1:ncol(XX),sep="")
    }
    if (n > 1) {
      if (length(colnames(YY))==0) colnames(YY) = paste("YY",1:ncol(YY),sep="")
      if (res$type == "none")    var_yy = vars::VAR(y = YY, p = p[s, 1], type = "none")
      if (res$type == "const")   var_yy = vars::VAR(y = YY, p = p[s, 1], type = "const")
      if (res$type == "exog0")   var_yy = vars::VAR(y = YY, p = p[s, 1], type = "none",  exogen = XX)
      if (res$type == "exog1")   var_yy = vars::VAR(y = YY, p = p[s, 1], type = "const", exogen = XX)

      for (i in 1:n) {
        Names_of_coefficits = names(var_yy$varresult[[i]]$coefficients)
        LREG = stats::lm(Ys[, i] ~ 0 + Zs)
        var_yy$varresult[[i]] <- LREG
        names(var_yy$varresult[[i]]$coefficients) <- Names_of_coefficits
      }
      var_yy$datamat[, ] <- Z[select[, s], ]
      if (p[s, 1] > 1) {
        var_yy$y[, ] <- t(cbind(t(Ys[1:p[s, 1], ]), t(Ys)))
      }
      else {
        var_yy$y[, ] <- t(cbind(as.matrix(Ys[1:p[s, 1], ]), t(Ys)))
      }
    }
    else {
      var_yy = stats::lm(Ys ~ 0 + Zs)
    }
    vars_obj[[s]] = var_yy
  }
  res$Bo <- b
  res$Co <- cc
  res$Sigmao <- sigma
  res$St <- St
  res$resid <- resid

  MRVAR_result = list()
  for (i in 1:S) {
    MRVAR_result[[i]] = summary(vars_obj[[i]])
  }

  Summary = list(MRVAR_result, LH_AIC, LH_BIC, LH_P, LH_N, LH_NN, LLH_AIC, LLH_BIC, ORBIC, ORAIC)
  names(Summary) = c("Estimation_Result", "LH_AIC", "LH_BIC", "LH_P", "LH_N", "LH_NN", "LLH_AIC", "LLH_BIC", "ORBIC", "ORAIC")
  res$Summary = Summary
  res$vars_obj= vars_obj
  return(res)
}

#' Join exogenous variable matrix for MRVAR
#'
#' @param X The data matrix of exogenous variables
#'
#' @return An array of containing exogenous data
#' @export
JointX = function(X) {
### X: the T x Nx x S array collecting the exogenous variables in MRVAR with state dependent exog variables
	S  = dim(X)[3]
      Nx = dim(X)[2]
      ss = 0
      sidx = 0
	for (s in 1:S) {
        if (sum(colSums(X[,,s] != 0) !=0)>ss ) {
           ss = sum(colSums(X[,,s] != 0) !=0)
           sidx = s
      }
      }
      XX = X[,1:ss,sidx]
      return(XX)
}




#' Generalized impulse response functions of MRVAR(n,p,S) with regime migrations
#'
#' This function calculates the generalized impulse response functions of an estimated MRVAR(n,p,S) for a given shock vector SHCK.
#'
#'
#'
#' ### GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK))
#'
#' See H.H. Pesaran and Y. Shin (1998) Generalized impulse response analysis in linear multivariate models, Economics Letters, 58(1) p. 17-29.
#' and G. Koop, M. H. Pesaran, and S. M. Potter (1996), Impulse response analysis in nonlinear multivariate models, Journal of Econometrics, 74 (1996) 119-74.
#'###########################################################################################################################################
#' @param  RES   : an MRVAR object that is an output of MRVARest.
#' @param  shock : an n vector containing the shocks as impulse.
#' @param  R     : the number runs to integrate out the random effects in order to obtain the means (see equation above).
#' @param  nstep : the length of the responses
#' @param  Omega_hist : a (P x n) matrix of initial values, from which the impulse response functions start. Omega_hist determines from which regime the impulse response functions start. For Omega_hist=NA, the impulse response functions will start from the most resent observations.
#' @param  resid_method : resid_method = c("resid", "parametric"), It generate the random residuals from residuals bootstrap or parametric bootstrap.
#' @return an (n x n x nstep+1) matrix of impulse response functions. The rows represent response the columns represent impulses.
#' @examples
#' p = matrix(c(2,1,0,0),2,2)
#' res_d = MRVARData(n=2,p=p,T=300,S=2,SESVI=1,type="none")
#' max(res_d$Y)
#' colnames(res_d$Y) = c("R","P")
#' res_e = MRVARest(res=res_d)
#' RF3 = girf_MRVAR_RM(RES=res_e,shock=c(1,1),R = 100,nstep=20,Omega_hist=NA,resid_method="parametric")
#' RF4 = girf_MRVAR_RM_CB(RES=res_e, shock=c(1,1), R=100, nstep=20, Omega_hist=NA,
#' resid_method = "parametric", conf_level=c(0.05,0.95), N=100)
#' IRF_list <-IRF_graph(RF4)
#'
#' @export
girf_MRVAR_RM <- function(RES,shock,R,nstep,Omega_hist=NA,resid_method) {
  ####  this function generate the impulse response function of MRVAR with migration
  ####
  ####                  	GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK))
  ####
  ####
  #### RES is the output of MRVARest
  n = RES$n
  p = RES$p
  S = RES$S
  SESVI= RES$SESVI; if (is.null(SESVI)) SESVI = NA
  TH= RES$TH;    if (is.null(TH))       TH = NA
  Bo= RES$Bo
  Co= RES$Co
  Sigmao= RES$Sigmao
  d= RES$d
  type= RES$type
  P           = max(p,d);
  X= RES$X;  if (!anyNA(X))  X = X[1:(P+nstep+1),]
  Yo= RES$Yo; if (is.null(Yo)) Yo = NA
  Uo          = RES$Uo
  #res = RES$res
  vars_obj = RES$vars_obj
  TT          = dim(RES$Y)[1]

  IRF = list()
  YR          = list()
  YS          = list()
  residR   <-  list()
  residS   <-  residR

  if (anyNA(Omega_hist)) Omega_hist = RES$Y[(dim(RES$Y)[1]-P+1):dim(RES$Y)[1],,drop=FALSE]
  state <- as.numeric(Omega_hist[P-d+1,SESVI]>TH) + 1

  #R = 400
  #shock       = (1:n)/(1:n)
  DIMresid    = dim(Uo)
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
    residI[P+1,,,state] <-shock_mat
    residS[[i]] = residI

    for (k in 1:n) {
      YR[,,k] = MRVARData(n=n,p=p,T=(P+nstep+1),S=S,SESVI=SESVI,TH=TH,Bo=Bo,Co=Co,Sigmao=Sigmao,Uo=residR[[i]][,,k,],type=type,X=X,Yo=Yo,d=d)$Y
      YS[,,k] = MRVARData(n=n,p=p,T=(P+nstep+1),S=S,SESVI=SESVI,TH=TH,Bo=Bo,Co=Co,Sigmao=Sigmao,Uo=residS[[i]][,,k,],type=type,X=X,Yo=Yo,d=d)$Y
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




#' Generalized impulse response functions of MRVAR(n,p,S) with regime migrations and confidence bands
#'
#' This function calculates the generalized impulse response functions of an estimated MRVAR(n,p,S) for given a shock vector. It provides also a bootstrapped confidence interval.
#'
#'
#'
#' ### GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK))
#'
#' See H.H. Pesaran and Y. Shin (1998) Generalized impulse response analysis in linear multivariate models, Economics Letters, 58(1) p. 17-29.
#' and G. Koop, M. H. Pesaran, and S. M. Potter (1996), Impulse response analysis in nonlinear multivariate models, Journal of Econometrics, 74 (1996) 119-74.
#' @param  RES   : an estimated MRVAR object.
#' @param  shock : an n-vector containing the shocks as impulse.
#' @param  R     : number of runs to integrate out the random effects in order to obtain the means (see equation above).
#' @param  nstep : length of the impuse response function
#' @param  Omega_hist : a (P x n) matrix of initial values, from which the impulse response functions start. Omega_hist determines from which regime the impulse response functions start.
#' @param  resid_method : resid_method = c("resid", "parametric"), It generate the random residuals from residuals bootstrap or parametric bootstrap.
#' @param  conf_level : a vector contain the level of confidences
#' @param  N            : the number of bootstrap runs in containing the bootstrapped confidence intervals.
#' @return an (n, n, nstep+1,3) array containing the impulse response functions with lower and upper confidence bonds. The rows represent response, and the columns represent impulses.
#' @examples
#' p = matrix(c(2,1,0,0),2,2)
#' res_d = MRVARData(n=2,p=p,T=300,S=2,SESVI=1,type="none")
#' max(res_d$Y)
#' colnames(res_d$Y) = c("R","P")
#' res_e = MRVARest(res=res_d)
#' RF3 = girf_MRVAR_RM(RES=res_e,shock=c(1,1),R = 200,nstep=20,Omega_hist=NA,resid_method="parametric")
#' RF4 = girf_MRVAR_RM_CB(RES=res_e, shock=c(1,1), R=200, nstep=20, Omega_hist=NA,
#' resid_method = "parametric", conf_level=c(0.05,0.95), N=100)
#' IRF_list <-IRF_graph(RF4)
#' @export
girf_MRVAR_RM_CB <- function (RES, shock, R, nstep, Omega_hist=NA, resid_method = "parametric", conf_level, N)
{
  n = RES$n
  p = RES$p
  S = RES$S
  SESVI = RES$SESVI
  TH = RES$TH
  Bo = RES$Bo
  Co = RES$Co
  Sigmao = RES$Sigmao
  d = RES$d
  type = RES$type
  X = RES$X
  T = dim(RES$Y)[1]
  P = max(p, d)
  GIRF = (1:(n * n * (nstep + 1) * N)) * 0
  dim(GIRF) = c(n, n, nstep + 1, N)
  GIRFBd = (1:(n * n * (nstep + 1) * (length(conf_level) + 1))) * 0
  dim(GIRFBd) = c(n, n, nstep + 1, length(conf_level) + 1)
  GIRFBd[, , , 1] = girf_MRVAR_RM(RES, shock = shock, R = 2*R, nstep, Omega_hist=Omega_hist, resid_method)
  for (i in 1:N) {
    res_d = MRVARData(n = n, p = p, T = T, S = S, SESVI = SESVI, TH = TH, Bo = Bo, Co = Co, Sigmao, type = type, X = X, d = d)
    if (length(colnames(res_d$Y))==0) colnames(res_d$Y) = paste("Y",1:ncol(res_d$Y),sep="")
    RESS = MRVARest(res_d)
    RF3 = girf_MRVAR_RM(RES = RESS, shock = shock, R = R, nstep, Omega_hist=Omega_hist, resid_method)
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






#' Initial U1 for DATA generating IRF
#'
#' @param B The coefficient matrix of a VAR object
#' @param sigma The covariance matrix of a VAR object
#' @param irf Types of impulse response functions
#' @param comb Types of concerted policy actions
#' @param G The matrix used in the permanent and transitory decomposition
#' @param smat An explicit decomposition matrix that defines a structural shock.
#'
#' @return The transformation matrix for shock decomposition
#' @export
U1=function(B,sigma,irf = c("gen", "chol", "chol1","gen1","genN1","comb", "comb1","smat","concerts1","concerts0","concertc"),comb,G=NA,smat=NA)
  {
    ##-----debug--------
    ##browser()
    ##------------------
    neq <- dim(B)[1]
    nvar <- dim(B)[2]
    lags <- dim(B)[3]
    n <- dim(sigma)[1]
    if (irf == "smat")     {  smat = (smat)    }

    if (irf == "chol") 	   {        smat = chol(sigma)    					    }
    if (irf == "PTdecomp") {        smat = chol(G%*%sigma%*%t(G))          		    }
    if (irf == "gen")      {        smat = t(sigma %*% diag(1/(sqrt(diag(sigma)))))     }
    if (irf == "chol1")    {        smat = chol(sigma) %*% diag(1/diag(chol(sigma)))    }
    if (irf == "gen1")     {        smat = t(sigma %*% diag(1/(diag(sigma))))           }
    if (irf == "genN1")    {        smat = t(sigma %*% diag(-1/(diag(sigma))))          }
    if ( irf == "concerts1") { smat  = t((sigma%*%diag(1/(diag(sigma))))%*%comb)        };
    if ( irf == "concerts0") { smat  = t((sigma%*%diag(1/(sqrt(diag(sigma)))))%*%comb)  };

    if ( irf == "concertc") {
      c     = as.numeric(!(comb[,1]==0));
      smat  = matrix(0,n,n);
      for (i in 1: n) { smat[i,] = sigma[i,]%*%diag(c)%*%INVI(sigma,c,i)%*%comb; }
      smat = t(smat);
    };

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
    if(dim(smat)[2] != dim(B)[2]) stop("B and smat conflict on # of variables")
    return(smat)
  }


#' Regime specific impulse response functions of MRVAR(n,p,S)
#'
#' This function calculates the regime specific impulse response functions of an estimated MRVAR(n,p,S), using Bo\[,,,s\] and Sigma\[,,s\] matrices of the estimated MRVAR(n,p,S).
#' @param RESS an output of MRVARest
#' @param nstep the length of impulse response function
#' @param comb a vector specify the weights used in GVAR models. For MRVAR its value is NA.
#' @param irf types of the impulse response functions. irf=c("gen","chol","chol1","gen1","comb1"), gen for generalized IRF with one standard deviation shocks, gen1 for generalized IRF with one unit impulse, chol for IRF with Cholezky decomposition of the covariance matrix, chol1 for Cholezky decomposition with one unit impulse, comb1 for concerted action with one unit impulse in GVAR.
#' @return an (n,n,nstep,2) array of IRF with columns representing the impulse and rows the responses.
#' @examples
#'
#' p = matrix(c(2,1,0,0),2,2)
#' res_d = MRVARData(n=2,p=p,T=300,S=2,SESVI=1)
#' max(abs(res_d$Y))
#' res_e  = MRVARest(res_d)
#' res_e$Summary
#'
#' IRF    = irf_MRVAR_NM(res_e,nstep=10,comb=NA,irf="gen")
#' IRF_CB = irf_MRVAR_NM_CB(res_e,nstep=10,comb=NA,irf="gen",runs=200,conf=c(0.05,0.90))
#' IRF_list1 <-IRF_graph(IRF_CB[,,,,1])
#' IRF_list2 <-IRF_graph(IRF_CB[,,,,2])
#'
#' @export
irf_MRVAR_NM <-function (RESS, nstep, comb, irf = c("gen", "chol", "chol1", "gen1", "comb1"))
{
  p = RESS$p
  Bo = RESS$Bo
  Sigmao = RESS$Sigmao
  neq = dim(Bo)[1]
  nvar = dim(Bo)[2]
  S = dim(Bo)[4]
  response <- array(0, dim = c(neq, nvar, nstep, S))
  BBo = list()
  for (s in 1:S) {
    BBo[[s]] = Bo[, , 1:p[s, 1], s]
    if (length(dim(BBo[[s]])) < 3)
      dim(BBo[[s]]) = c(dim(BBo[[s]]), 1)
  }
  for (s in 1:S) response[, , , s] <- irf_B_sigma(B=BBo[[s]], sigma=Sigmao[, , s], nstep, comb, irf = irf)
  return(response)
}



#' Regime specific impulse response functions of MRVAR(n,p,S)
#'
#' This function calculates the regime specific impulse response functions of an estimated MRVAR(n,p,S), using Bo\[,,,s\] and Sigma\[,,s\] of the estimated MRVAR.
#' @param RESS an output of MRVARest.
#' @param nstep the length of impulse response function
#' @param comb a vector specify the concerted action in policy-simulation impulse response function
#' @param irf types of the impulse response functions. irf=c("gen","chol","chol1","gen1","comb1"), gen for generalized IRF with one standard deviation shocks, gen1 for generalized IRF with one unit impulse, chol for IRF with Cholezky decomposition of the covariance matrix, chol1 for Cholezky decomposition with one unit impulse, comb1 for concerted action with one unit impulse in GVAR.
#' @param runs number of simulation runs
#' @param conf a two components vector containing the tail probabilities of the bootstrap confidence interval.
#' @return an (n,n,nstep,3,2) array containing the IRF of the two regimes. The IRF columns represent the impulse and the rows the responses.
#' @examples
#' p = matrix(c(2,1,0,0),2,2)
#' res_d = MRVARData(n=2,p=p,T=300,S=2,SESVI=1)
#' max(abs(res_d$Y))
#' res_e  = MRVARest(res_d)
#' res_e$Summary
#'
#' IRF    = irf_MRVAR_NM(res_e,nstep=10,comb=NA,irf="gen")
#' IRF_CB = irf_MRVAR_NM_CB(res_e,nstep=10,comb=NA,irf="gen",runs=200,conf=c(0.05,0.90))
#' IRF_list1 <-IRF_graph(IRF_CB[,,,,1])
#' IRF_list2 <-IRF_graph(IRF_CB[,,,,2])
#'
#' @export
irf_MRVAR_NM_CB <- function (RESS, nstep, comb, irf = c("gen", "chol", "chol1", "gen1", "comb1"), runs = 200, conf = c(0.05, 0.95))
{
  n = RESS$n
  p = RESS$p
  T = dim(RESS$Y)[1]
  S = RESS$S
  SESVI = RESS$SESVI
  TH = RESS$TH
  Bo = RESS$Bo
  Co = RESS$Co
  Sigmao = RESS$Sigmao
  SV = RESS$SV
  if (is.null(SV))
    SV = NA
  type = RESS$type
  X = RESS$X
  if (is.null(X))
    X = NA
  mu = RESS$mu
  if (is.null(mu))
    mu = NA
  Yo = RESS$Yo
  if (is.null(Yo))
    Yo = NA
  Do = RESS$Do
  if (is.null(Do))
    Do = NA
  Uo_run = RESS$Uo
  neq = dim(Bo)[1]
  nvar = dim(Bo)[2]
  U_run = Uo_run * 0
  response <- array(0, dim = c(neq, nvar, nstep, length(conf) +
                                 1, S))
  response[, , , 1, ] <- irf_MRVAR_NM(RESS, nstep, comb, irf)
  responseR <- array(0, dim = c(neq, nvar, nstep, runs, S))
  for (i in 1:runs) {
    for (s in 1:S) {
      U_run[, , s] = rnormSIGMA(T, Sigmao[, , s])
      B = Bo[, , 1:p[s, 1], s]
      if (length(dim(B)) == 2)
        dim(B) = c(dim(B), 1)
      res_run = VARData(n = n, p = p[s, 1], T = T, B = B,
                        Co = 0 * (1:n), Sigma = Sigmao[, , s], U = U_run[,
                                                                         , s], type = "none")
      res_e = VARest(res_run)
      responseR[, , , i, s] <- irf_B_sigma(res_e$B, res_e$Sigma, nstep, comb, irf)
    }
  }
  for (tt in 1:(nstep)) {
    for (i in 1:neq) {
      for (j in 1:nvar) {
        for (s in 1:S) response[i, j, tt, -1, s] = stats::quantile(responseR[i, j, tt, , s], conf)
      }
    }
  }
  return(response)
}



#' Regime specific impulse response functions of MRVAR(n,p,S)
#'
#' This function calculates the regime specific impulse response functions of an estimated MRVAR(n,p,S).
#' @param res_e an output of MRVARest.
#' @param nstep the length of impulse response function
#' @param comb a vector specify the concerted action in policy-simulation impulse response function
#' @param irf types of the impulse response functions. irf=c("gen","chol","chol1","gen1","comb1"), gen for generalized IRF with one standard deviation shocks, gen1 for generalized IRF with one unit impulse, chol for IRF with Cholezky decomposition of the covariance matrix, chol1 for Cholezky decomposition with one unit impulse, comb1 for concerted action with one unit impulse in GVAR.
#' @param Xshks number of selected exogenous shocks
#' @param runs number of simulation runs
#' @param conf a two components vector containing the tail probabilities of the bootstrap confidence interval.
#' @return an (n,n,nstep,3,2) array containing the IRF of the two regimes. The IRF columns represent the impulse and the rows the responses.
#' @examples
#' p = matrix(c(2,1,0,0),2,2)
#' res_d = MRVARData(n=2,p=p,T=300,S=2,SESVI=1)
#' max(abs(res_d$Y))
#' res_e  = MRVARest(res_d)
#' res_e$Summary
#'
#'
#' IRF_CB = irf_MRVAR_CB(res_e,nstep=20,comb=NA,irf="gen1",runs=100,conf=c(0.05,0.90))
#' IRF_list1 <-IRF_graph(IRF_CB[[1]])
#' IRF_list2 <-IRF_graph(IRF_CB[[2]])
#'
#' @export
irf_MRVAR_CB <- function (res_e, nstep, comb, irf = c("gen", "chol", "chol1", "gen1", "comb1", "irfX"), Xshks=NA, runs = 200, conf = c(0.05, 0.95))
{
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

  Xshk1 = NA
  Xshk2 = NA
  if (irf=="irfX") {
    Xshk1 = t(Co[,1+Xshks,1])
    Xshk2 = t(Co[,1+Xshks,2])
  }

  IRF1 = irf_B_sigma(B = res_e$Bo[, , , 1], sigma = res_e$Sigma[,, 1], nstep = nstep, comb = NA, irf = irf, G = NA, smat = NA,Xshk=Xshk1)
  IRF2 = irf_B_sigma(B = res_e$Bo[, , , 2], sigma = res_e$Sigma[,, 2], nstep = nstep, comb = NA, irf = irf, G = NA, smat = NA,Xshk=Xshk2)

  neq = res_e$n
  nvar = res_e$n
  response1 <- array(0, dim = c(neq, nvar, nstep, length(conf) + 1))
  response2 <- array(0, dim = c(neq, nvar, nstep, length(conf) + 1))
  response1[, , , 1] <- IRF1
  response2[, , , 1] <- IRF2
  responseR <- array(0, dim = c(neq, nvar, nstep, runs, S))
  BoColect = array(0, c(dim(Bo), runs))
  UoColect = array(0, c(dim(Uo), runs))
  CoColect = array(0, c(dim(Co), runs))
  YColect = array(0, c(dim(Y), runs))
  Uo_run = array(0, c(T, n, S))

  for (i in 1:runs) {
    for (s in 1:S) Uo_run[, , s] = rnormSIGMA(T, Sigmao[,, s])
    if (length(p) > 1) {
      res_run = MRVARData(n = n, p = p, T = T, S = S, SESVI = SESVI, TH = TH, Bo = Bo, Co = Co, Sigmao = Sigmao, Uo = Uo_run, type = type, X = X)
      res_e_run = MRVARest(res = res_run)
    }
    Co_run <- res_e_run$Co
    Xshk1_run = NA
    Xshk2_run = NA

    if (irf == "irfX") {
      Xshk1_run = t(Co_run[,1+Xshks,1])
      Xshk2_run = t(Co_run[,1+Xshks,2])
    }
    IRF1 = irf_B_sigma(B = res_e_run$Bo[, , , 1], sigma = res_e_run$Sigma[, , 1], nstep = nstep, comb = NA, irf = irf, G = NA,Xshk=Xshk1_run)
    IRF2 = irf_B_sigma(B = res_e_run$Bo[, , , 2], sigma = res_e_run$Sigma[, , 2], nstep = nstep, comb = NA, irf = irf, G = NA,Xshk=Xshk2_run)
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
        response1[i, j, tt, -1] = stats::quantile(responseR[i,
                                                            j, tt, , 1], conf)
        response2[i, j, tt, -1] = stats::quantile(responseR[i,
                                                            j, tt, , 2], conf)
      }
    }
  }
  return(list(response1, response2, BoColect, UoColect, YColect))
}



#' Calculation of information criteria for MRVAR(n,p,S) models.
#'
#' @param  res  : an MRVAR object generated from MRVARData or estimated from MRVARest.
#' @param  L_V  : a 2 components vector containing the maxima lags in the two regimes, respectively.
#' @param  TH_V  : a vector containing the possible threshold values over which the model selection criteria values will be calculated.
#' @return      a matrix with different lag specifications and threshold values as well as the corresponding information criterion values.
#' @examples
#'
#' res_d = MRVARData(n=4,p=matrix(c(2,1,2,2,0,0,0,0),4,2),T=800,S=2,SESVI=1)
#' max(res_d$Y)
#' colnames(res_d$Y) = c("P","Y","R","U")
#' res_e = MRVARest(res=res_d)
#' TH_v = c(0,0.0)
#' L_v = c(5,5)
#' Sel = MRVAR_Select(res=res_e,L_V=L_v,TH_V=TH_v)
#' MRVAR_Select_Summary(Sel)
#' @export
MRVAR_Select <- function (res, L_V = L_V, TH_V = TH_V)
{
  res_dd = res
  p = res_dd$p
  n = res_dd$n
  T = dim(res_dd$Y)[1]
  SESVI = res_dd$SESVI
  type = res_dd$type
  X = res_dd$X
  p[, 1] = max(L_V)
  S = res_dd$S
  res_dd = MRVARData(n = n, p = p, T = T, S = S, SESVI = SESVI,
                     type = type, X = X)
  res_dd$Y = res$Y
  Criteria = matrix(0, L_V[1] * L_V[2] * length(TH_V), 13)
  idx = 0
  for (l_d in 1:L_V[1]) {
    for (l_f in 1:L_V[2]) {
      for (l_th in 1:length(TH_V)) {
        idx = idx + 1
        res_dd$p[1, 1] = l_d
        res_dd$p[2, 1] = l_f
        res_dd$TH = TH_V[l_th]
        res_ss = MRVARest(res_dd)
        TT = res_ss$Summary$LH_N + res_ss$Summary$LH_P
        Criteria[idx, ] = c(l_d, l_f, TH_V[l_th], res_ss$Summary$LH_AIC,
                            res_ss$Summary$LH_BIC, res_ss$Summary$LH_P, res_ss$Summary$LH_N, res_ss$Summary$LH_NN,
                            res_ss$Summary$ORAIC, res_ss$Summary$ORBIC, TT, res_ss$Summary$LLH_AIC,
                            res_ss$Summary$LLH_BIC)
        colnames(Criteria) = c("Lag_regime1", "Lag_regime2",
                               "threshold", "AIC", "BIC",
                               "LH_P", "LH_N", "LH_NN",
                               "ORAIC", "ORBIC", "TT", "LLH_AIC",
                               "LLH_BIC")
      }
    }
  }
  return(Criteria)
}


#' Summary result of MRVAR model selection
#'
#' @param Sel An output of MRVAR_Select
#'
#' @return The optimal model according to BIC or AIC
#' @export
#'
MRVAR_Select_Summary = function(Sel) {
  result = list()
  result[[1]] = Sel[which.min(Sel[,5]),c(1,2,3,5)]
  result[[2]] = Sel[which.min(Sel[,4]),c(1,2,3,4)]
  return(result)
}
