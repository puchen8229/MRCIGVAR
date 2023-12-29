#' Data generating process of MRGVAR(m,n,p,S)
#'
#' This function generates data from an multi-regime stationary global VAR(p) process and returns an MRGVAR object that is a list containing the generated data and the used parameters.
#'
#' @param m     : number of variables in a country/unit
#' @param n     : number of countries/units
#' @param p     : an (n, 3, S) array, each raw specifies the lag length of the domestic variables, the lag length of the foreign variables and the number of exogenous variables for the respective regime.
#' @param T     : number of observations.
#' @param S     : number of regimes.
#' @param W     : an (n x n) weighting matrix. w_ij is the weight of country j in the foreign variables of i-th country diag(W)=0
#' @param SESVI : an n-vector of indices of the switching variables across n countries. Eg. SESVI = seq(1,m*n,m).
#' @param TH    : an (n, S-1) matrix of threshold values
#' @param Go    : an (mn,mn,p,S) array of the MRGVAR(m,n,p,S) coefficients. G is constructed from Bo, Ao and W.
#' @param Ao    : an (m, m, p, n, S) array collecting coefficients of foreign variables
#' @param Bo    : an (m, m, p, n, S) array collecting coefficients of domestic variables
#' @param Uo    : a (T, mn, S) array of the temporally independent innovation processes
#' @param Sigmao : an (mn, mn, S) array of the covariance matrix of MRGVAR(m,n,p,S)
#' @param SV    : exogenous switching variables
#' @param type	: types of deterministic component "const", "none", "exog0", and  "exog1" are 4 options
#' @param Co    : an (m , k+1, n, S) array collecting the coefficients of the deterministic components of the n countries for every regime.
#' @param X	    : a (T x k x n x S) matrix of exogenous variables.
#' @param Yo    : Initial values
#' @param d     : the time lag between signal and switching
#' @return      a MRGVAR object containing the generated data, the used parameters and the exogenous variables.
#' res_d = list(Y,X,Uo,resid,Go,GDC,Co,Sigmao,TH,St,SV,SESVI,Ao,Bo,check,type,m,n,p,S,W,SigmaS,Yo,d)
#' \itemize{
#'    \item Y     : a (T x nm)  matrix of simulated data via of the MRGVAR(m,n,p,S)
#'    \item X     : a (T x k x n x S) matrix of exogenous variables.
#'    \item Uo    : a (T, mn,S) array of the simulated innovations of the MRGVAR(m,n,p,S)
#'    \item C     : an (nm, (k+1),S) array containing the coefficients of the deterministic components.
#'    \item St    : a (T x n) matrix of the simulated time path of states/regimes
#'    \item check : maximum of the data for checking the stationarity
#' }
#'
#' @examples
#' ## case of n = 2, m = 2, S = 2     ## m: number of variables, n: number of countries
#' p = rep(1,12); dim(p) = c(2,3,2)
#' p[1,1,2] = 2; p[2,2,2]=2; p[,3,] = 0
#' TH = c(1:2)*0; dim(TH) = c(1,2)
#' res_d <- MRGVARData(m=2,n=2,p=p,TH=TH,T=100,S=2,SESVI=c(1,3),type="const")
#' max(res_d$Y)
#'
#' ### estimation of the MRGVAR model
#' colnames(res_d$Y) = c("P","Q","Pa","Qa")
#' res_e = MRGVARest(res=res_d)
#' res_e$Summary
#'
#' IRF_CB  = irf_MRGVAR_CB(res=res_e,nstep=10,comb=NA,state=c(1,1),irf="gen1",runs=20,
#' conf=c(0.05,0.95))
#' IRF_g = IRF_graph(IRF_CB[[1]],Names=c("P","Q","Pa","Qa"))    #IRF
#' #IRF_g = IRF_graph(IRF_CB[[2]])   # accumulated IRF
#' @export

MRGVARData=function(m,n,p,T,S,W=NA,SESVI=NA,TH=NA,Go=NA,Ao=NA,Bo=NA,Sigmao=NA,Uo=NA,SV=NA,type=NA,Co=NA,X=NA,Yo=NA,d=NA) {
  ### m     : number of variables in a country
  ### n     : number of countries
  ### p     : lag length n x 3 x S array collecting lag length for each country with respect to domestic and foreign variable for each state. The last column specifies the number of exogenous variables for each state.
  ### T     : number of observations
  ### S     : number of regimes
  ### SESVI : n vector of indexes of the switching variable in the endogenous variables Y for the case of self-excited threshold model.
  ###         (m,n,p,T,S,SESVI) are parameters which must be provided.
  ### TH    : n x S-1  matrix of threshold values of n countries.
  ### Go    : mn x mn x p x S GVAR(m,n,p) coefficients matrix of S different regimes
  ### Sigmao: mn x mn x S array of the covariance matrix of the GVAR(m,n,p) in S different regimes
  ### Uo    : a T x mn x S  array of the temporally independent innovation processes
  ### X	    : a T x k x n x S array of exogenous variables which may be common/different to all countries and different across all states
  ###         (TH,Go,Sigmao,Uo,SV) if not provided, they will be generated randomly.
  ### type  : deterministic component "const", "none", "exog0", and "exog1" are foure options
  ### Co    : if type = "const" mu is an m x n x S array of the intercepts of the time series in the different regimes
  ###
  ###
  ### output:
  ### c("Y","Uo","Go","Sigmao","TH","St","sv","SESVI","check")
  ###
  ### Y     : T x mn data matrix of the simulated data via of the MSGVAR(m,n,p,S)
  ### Uo    : an T x mn x S   array of the temporally independent innovation processes of the MSGVAR(m,n,p,S)
  ### Go    : mn x mn x p x S array collecting the MSGVAR(m,n,p,S) coefficients in S different states
  ### Sigmao: mn x mn x S array collecting the covariance matrices of the simulated MSGVAR(m,n,p,S) in S different states
  ### TH    : S-1 vector of thresholds
  ### St    : simulated time path of states/regimes
  ### SESVI : index if the switching variable
  ### check : maximum of the data to check stationarity
  ###
  ### Remarks: The states of each country at each time step is governed by the switching variables of each country and hence there is a large numbers
  ###          of possible combinations of regimes.
  ###          The coefficients matrix can be constructed from Go, given a regime combination. The stationarity of the MRGVAR(m,n,p)
  ###          at each time step is not guaranteed. But this is not relevant. (This is a open question.)

  check = c(1:(S + 1)) * 0
  if (missing(TH)) {
    TH = NA
  }
  if (missing(W)) {
    W = NA
  }
  if (missing(SV)) {
    SV = NA
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
  if (missing(Go)) {
    Go = NA
  }
  if (missing(Bo)) {
    Bo = NA
  }
  if (missing(X)) {
    X = NA
  }
  if (missing(SV)) {
    SV = NA
  }
  if (missing(d)) {
    d = NA
  }
  if (missing(Yo)) {
    Yo = NA
  }
  if (anyNA(d)) {
    d = 1
  }
  P = max(p[, 1:2, ], d)
  Pmax = max(p[, 1:2, ])
  if (anyNA(TH)) {
    TH = matrix(stats::runif(S - 1) * n, S - 1, n)
    for (i in 1:n) {
      THh = stats::runif(S - 1) - 0.5
      THh = THh[order(THh)]
      TH[, i] = THh
    }
  }
  if (anyNA(W)) {
    W = matrix(((1:(n * n))/(1:(n * n)) * 1/(n - 1)), n,
               n)
    for (i in 1:n) W[i, i] = 0
  }
  if (anyNA(Sigmao)) {
    Sigmao = (1:(m * n * m * n * S)) * 0
    dim(Sigmao) = c(m * n, m * n, S)
    for (s in 1:S) {
      GVARD = GVARData(m, n, p[, , s], T, W = W)
      Sigmao[, , s] = GVARD$Sigmao
    }
  }
  if (anyNA(Uo)) {
    Uo = (1:(T * n * m * S)) * 0
    dim(Uo) = c(T, m * n, S)
    for (s in 1:S) Uo[, , s] = rnormSIGMA(T, Sigmao[, , s])
  }
  Ct = Uo * 0
  if (anyNA(type)) {
    type = "none"
  }
  if (type == "none") {
    CDC = c(1:(m * n * S)) * 0
    dim(CDC) = c(m, 1, n, S)
    Co = CDC
    GDC = Co
    dim(GDC) = c(n * m, 1, S)
    Ct = Uo * 0
  }
  if (type == "const") {
    CDC = stats::rnorm(m * n * S)
    dim(CDC) = c(m, 1, n, S)
    if (anyNA(Co))
      Co = CDC
    GDC = Co
    dim(GDC) = c(n * m, 1, S)
    CoV = as.vector(Co)
    for (s in 1:S) {
      CoV = as.vector(Co[, , , s])
      for (i in 1:n) Ct[, , s] = matrix(1, T, 1) %*% t(CoV)
    }
  }
  if (type == "exog1") {
    k = dim(X)[2] + 1
    CDC = stats::rnorm(m * k * n * S)
    dim(CDC) = c(m, k, n, S)
    if (anyNA(Co))
      Co = CDC
    GDC = Co
    dim(GDC) = c(n * m, k, S)
    CoV = matrix(0, k, m * n)
    for (ss in 1:S) {
      CoV = matrix(0, k, m * n)
      for (j in 1:k) CoV[j, ] = as.vector(Co[, j, , ss])
      for (i in 1:n) Ct[, ((i - 1) * m + 1):((i - 1) *
                                               m + m), ss] = cbind(matrix(1, dim(X)[1], 1),
                                                                   as.matrix(X[, , i, ss])) %*% CoV[, ((i - 1) *
                                                                                                         m + 1):((i - 1) * m + m)]
    }
  }
  if (type == "exog0") {
    k = dim(X)[2] + 1
    CDC = stats::rnorm(m * k * n * S)
    dim(CDC) = c(m, k, n, S)
    CDC[, 1, , ] = 0
    if (anyNA(Co))
      Co = CDC
    GDC = Co
    dim(GDC) = c(n * m, k, S)
    CoV = matrix(0, k, m * n)
    for (ss in 1:S) {
      CoV = matrix(0, k, m * n)
      for (j in 1:k) CoV[j, ] = as.vector(Co[, j, , ss])
      for (i in 1:n) Ct[, ((i - 1) * m + 1):((i - 1) *
                                               m + m), ss] = cbind(matrix(1, dim(X)[1], 1),
                                                                   as.matrix(X[, , i, ss])) %*% CoV[, ((i - 1) *
                                                                                                         m + 1):((i - 1) * m + m)]
    }
  }


  kkk = 1
  repeat {
    kkk = kkk + 1
    if (anyNA(Go) & anyNA(Bo)) {
      Go = (1:(m * n * m * n * Pmax * S)) * 0
      dim(Go) = c(m * n, m * n, Pmax, S)
      Bo = (1:(m * m * Pmax * n * S)) * 0
      dim(Bo) = c(m, m, Pmax, n, S)
      Ao = Bo * 0
      for (s in 1:S) {
        GVARD = GVARData(m, n, p[, , s], T, W = W)
        Go[, , 1:max(p[, 1:2, s]), s] = GVARD$G
        check[s] = max(abs(GVARD$Y))
        Ao[, , 1:max(p[, 1:2, s]), , s] = GVARD$Ao
        Bo[, , 1:max(p[, 1:2, s]), , s] = GVARD$Bo
      }
    }
    if (!anyNA(Go) & anyNA(Bo)) {
      Bo = (1:(m * m * Pmax * n * S)) * 0
      dim(Bo) = c(m, m, Pmax, n, S)
      Ao = Bo * 0
      for (s in 1:S) {
        Bo[, , , , s] = GW2BoAo(Go[, , , s], W)$Bo
        Ao[, , , , s] = GW2BoAo(Go[, , , s], W)$Ao
      }
    }
    if (anyNA(Go) & !anyNA(Bo)) {
      Go = (1:(m * n * m * n * Pmax * S)) * 0
      dim(Go) = c(m * n, m * n, Pmax, S)
      for (s in 1:S) Go[, , , s] = BoAoW2G(Bo[, , , , s], Ao[,
                                                             , , , s], W, m, n, Pmax)
    }

    if (! max(Mod(STAT(Go[,,,1])),Mod(STAT(Go[,,,2]))) > 1) break
    if ( kkk>100 ) break
  }



  St = matrix(0, T, n)
  Y = Uo[, , 1]
  if (!anyNA(Yo))
    Y[1:P, ] = Yo
  resid = Y * 0
  Yo = Y[1:P, ]
  if (anyNA(SV)) {
    sv = Y[, SESVI]
    SV = NA
  }
  else {
    sv = SV
  }
  ss = matrix(0, n, 1)
  for (tt in (P + 1):T) {
    if (anyNA(SV)) {
      sv = Y[, SESVI]
      SV = NA
    }
    else {
      sv = SV
    }
    svt = sv[tt - d, ]
    Bt = c(1:(m * n * m * n * Pmax)) * 0
    dim(Bt) = c(n * m, n * m, Pmax)
    for (i in 1:n) {
      ss[i] = which.min(abs(TH[, i] - svt[i]))
      if (svt[i] > TH[ss[i], i]) {
        ss[i] = ss[i] + 1
      }
      Bt[(1 + (i - 1) * m):(i * m), , ] = Go[(1 + (i -
                                                     1) * m):(i * m), , , ss[i]]
      Y[tt, (1 + (i - 1) * m):(i * m)] = Uo[tt, (1 + (i -
                                                        1) * m):(i * m), ss[i]] + Ct[tt, (1 + (i - 1) *
                                                                                            m):(i * m), ss[i]]
      resid[tt, (1 + (i - 1) * m):(i * m)] = Uo[tt, (1 +
                                                       (i - 1) * m):(i * m), ss[i]]
    }
    for (L in 1:Pmax) Y[tt, ] = Y[tt, ] + Y[tt - L, ] %*%
      t(Bt[, , L])
    St[tt, ] = ss
  }
  check[S + 1] = max(abs(Y))
  SigmaS = NA
  result = list(Y, X, Uo, resid, Go, GDC, Co, Sigmao, TH, St,
                SV, SESVI, Ao, Bo, check, type, m, n, p, S, W, SigmaS,
                Yo, d, kkk)
  names(result) = c("Y", "X", "Uo", "resid",
                    "Go", "GDC", "Co", "Sigmao",
                    "TH", "St", "SV", "SESVI", "Ao",
                    "Bo", "check", "type", "m", "n",
                    "p", "S", "W", "SigmaS", "Yo",
                    "d","kkk")
  return(result)
}

#' Estimation of MRGVAR(m,n,p,S) models
#'
#' This function estimates the parameters of a specified MRGVAR(m,n,p,S) model based on provided data.
#'
#' @param  res  : an MRGVAR object which is an output of MRGVARData.
#' @return res  : an MRGVAR object with estimated parameters and test statistics.
#' @examples
#'
#' ## case of n = 2, m = 2, S = 2     ## m: number of variables, n: number of countries
#' p = rep(1,12); dim(p) = c(2,3,2)
#' p[1,1,2] = 2; p[2,2,2]=2; p[,3,] = 0
#' TH = c(1:2)*0; dim(TH) = c(1,2)
#' res_d <- MRGVARData(m=2,n=2,p=p,TH=TH,T=100,S=2,SESVI=c(1,3),type="const")
#' max(res_d$Y)
#'
#' ### estimation of the MRGVAR model
#' colnames(res_d$Y) = c("P","Q","Pa","Qa")
#' res_e = MRGVARest(res=res_d)
#' res_e$Summary
#'
#' IRF_CB  = irf_MRGVAR_CB(res=res_e,nstep=10,comb=NA,state=c(1,1),irf="gen1",runs=20,
#' conf=c(0.05,0.95))
#' IRF_g = IRF_graph(IRF_CB[[1]],Names=c("P","Q","Pa","Qa"))    #IRF
#' #IRF_g = IRF_graph(IRF_CB[[2]])   # accumulated IRF
#' @export
MRGVARest <- function (res)
{
  m = res$m
  n = res$n
  p = res$p
  Y = as.matrix(res$Y);  if (max(abs(Y)) > 1000000) return("Check data!!!");
  X = res$X
  W = res$W
  type = res$type
  TH = res$TH
  SESVI = res$SESVI
  Go = res$Go
  S = res$S
  Ao = res$Ao
  Bo = res$Bo
  Co = res$Co
  GDC = res$GDC
  Sigmao = res$Sigmao * 0
  d = res$d
  r_npo = res$r_npo
  P = max(p, d)
  Pmax = max(p[, 1:2, ])
  SigmaS = (1:(n * m * S * n * m * S)) * 0
  dim(SigmaS) = c(n * m * S, n * m * S)
  VAR_domestic = list()
  SESVIi = SESVI/SESVI + SESVI[1] - 1
  kmax = max(p[, 3, ])
  T = dim(Y)[1]
  FY = Y %*% t(W %x% diag(m))
  resid = (1:(T * m * n * S)) * 0
  dim(resid) = c(T, m, n, S)
  if (type == "none" | type == "const")
    k = 1
  if (type == "exog0" | type == "exog1")
    k = dim(X)[2] + 1
  CC = c(1:(m * (Pmax * m + k) * S)) * 0
  dim(CC) = c(m, m * Pmax + k, S)
  EMCVG = c(1:n) * 0
  for (i in 1:n) {
    Pi = max(p[i, 1:2, ])
    CCi = c(1:(m * (Pi * m + k) * S)) * 0
    dim(CCi) = c(m, m * Pi + k, S)
    Z = matrix(0, T, Pi * m)
    HilfeFYp = FY[, (m * (i - 1) + 1):(i * m)]
    if (is.null(colnames( HilfeFYp)))  colnames(HilfeFYp) = sprintf("FY%s", 1:m)
    ## FYp = embed(FY[, (m * (i - 1) + 1):(i * m)], (Pi + 1))
    FYp = Embed(HilfeFYp, (Pi + 1),prefix="")
    Z[(Pi + 1):T, ] = FYp[, (m + 1):(Pi * m + m)]
    HilfeZ          = FYp[, (m + 1):(Pi * m + m)]
    colnames_Zi     <- colnames(HilfeZ)            # for none and const


    kz = dim(Z)[2]
    if (type == "none") {
      Z = rep(Z, S)
      dim(Z) = c(T, kz, S)
      ZZ = array(0, c(T, kz, S))
      for (s in 1:S) ZZ[, 1:kz, s] = Z[, , s]
    }
    if (type == "const") {
      Z = rep(Z, S)
      dim(Z) = c(T, kz, S)
      ZZ = array(0, c(T, kz, S))
      for (s in 1:S) ZZ[, 1:kz, s] = Z[, , s]
    }
    if (type == "exog0") {
      Z = rep(Z, S)
      dim(Z) = c(T, kz, S)
      ZZ = array(0, c(T, kz + kmax, S))
      for (s in 1:S) ZZ[, 1:(kz + p[i, 3, s]), s] = cbind(Z[, , s], X[, 1:p[i, 3, s], i, s])
      colnames_Zi =  c(colnames_Zi, colnames(X)[1:p[i, 3,s]] )
    }
    if (type == "exog1") {
      Z = rep(Z, S)
      dim(Z) = c(T, kz, S)
      ZZ = array(0, c(T, kz + kmax, S))
      for (s in 1:S) ZZ[, 1:(kz + p[i, 3, s]), s] = cbind(Z[, , s], X[, 1:p[i, 3, s], i, s])

      colnames_Zi =  c(colnames_Zi, colnames(X)[1:p[i, 3,s]] )
    }
    if (type == "none" | type == "exog0")
      type_holder = "exog0"
    if (type == "const" | type == "exog1")
      type_holder = "exog1"
    XX = (1:(T * dim(ZZ)[2] * S)) * 0
    dim(XX) = c(T, dim(ZZ)[2], S)
    if (!anyNA(X))
      for (s in 1:S) XX[, 1:(m * p[i, 2, s] + p[i, 3, s]), s] = ZZ[, c(1:(m * p[i, 2, s]), (m * Pi + 1):(m * Pi + p[i, 3, s])), s]
    if (anyNA(X))
      for (s in 1:S) XX[, 1:(m * p[i, 2, s]), s] = ZZ[, c(1:(m * p[i, 2, s])), s]
    colnames(XX) = colnames_Zi
    pp = t(p[i, , ])
    pp[, 2] = pp[, 2] * m
    ppp = pp[, 1:2]
    ppp[, 2] = pp[, 2] + pp[, 3]
    res_MRVAR = MRVARData(n = m, p = ppp, T = T, S = S, Co = NA, SESVI = SESVIi[i], type = type_holder, X = XX)
    res_MRVAR$Y = Y[, (m * (i - 1) + 1):(i * m)]
    #res_MRVAR$X = XX
    res_MRVAR$SV = res_MRVAR$Y[, SESVIi[i]]
    res_MRVAR$TH = TH[, i]

    RR = MRVARest(res = res_MRVAR)
    VAR_domestic[[i]] = RR
    for (s in 1:S) {
      Bo[, , 1:p[i, 1, s], i, s] = RR$Bo[, , 1:p[i, 1, s], s]
      Ai = RR$Co[, 2:(m * p[i, 2, s] + 1), s]
      dim(Ai) = c(m, m, p[i, 2, s])
      Ao[, , 1:p[i, 2, s], i, s] = Ai
      if (type == "none")
        Co[, , i, s] = RR$Co[, 1, s]
      if (type == "const")
        Co[, , i, s] = RR$Co[, 1, s]
      if (type == "exog1")
        Co[, 1:(1 + p[i, 3, s]), i, s] = RR$Co[, c(1, (m * p[i, 2, s] + 2):(m * p[i, 2, s] + p[i, 3, s] + 1)), s]
      if (type == "exog0")
        Co[, 1:(1 + p[i, 3, s]), i, s] = RR$Co[, c(1, (m * p[i, 2, s] + 2):(m * p[i, 2, s] + p[i, 3, s] + 1)), s]
    }
    Sigmao[(m * (i - 1) + 1):(i * m), (m * (i - 1) + 1):(i * m), ] = RR$Sigmao
    resid[, , i, ] = RR$resid
    GDC = Co
    dim(GDC) = c(m * n, k, S)
  }
  for (s in 1:S) {
    for (i in 1:n) {
      for (j in 1:n) Sigmao[(m * (i - 1) + 1):(i * m),
                            (m * (j - 1) + 1):(j * m), s] = t(resid[, , i,
                                                                    s]) %*% (resid[, , j, s])/(dim(Z)[1] - dim(Z)[2])
    }
  }
  for (s in 1:S) {
    for (ss in 1:S) {
      for (i in 1:n) {
        for (j in 1:n) SigmaS[((s - 1) * n * m + (i - 1) * m + 1):((s - 1) * n * m + i * m), ((ss - 1) * n * m + (j - 1) * m + 1):((ss - 1) * n * m + j * m)] = t(resid[, , i, s]) %*% (resid[, , j, ss])/sum((!resid[, 1, i, s] == 0) * (!resid[, 1, j, ss] == 0))
      }
    }
  }
  for (s in 1:S) {
    BBo = Bo[, , , , s]
    AAo = Ao[, , , , s]
    dim(BBo) = c(m, m, Pmax, n)
    dim(AAo) = c(m, m, Pmax, n)
    Go[, , , s] = BoAoW2G(BBo, AAo, W, m, n, Pmax)
  }
  dim(resid) = c(T, m * n, S)
  res$Go = Go
  res$Sigmao = Sigmao
  res$Ao = Ao
  res$Bo = Bo
  res$Co = Co
  res$GDC = GDC
  res$SigmaS = SigmaS
  res$resid = resid
  res$VAR_domestic = VAR_domestic
  est_result = list()
  for (i in 1:n) {
    for (s in 1:S) {
      est_result[[(i-1)*S+s]] =summary(VAR_domestic[[i]]$vars_obj[[s]])
    }
  }
  Summary = list(est_result,CC,r_npo)
  names(Summary) = c("Estimation_Result","C","r_npo")
  res$Summary = Summary
  return(res)
}


#' Calculation of auto regression parameters of a GVAR model from its country equations
#'
#' @param Bo Parmeter matrices of the domestic variables
#' @param Ao Parameter matrices of the foreign variables
#' @param W The weighting matrix of the GVAR model
#' @param m Number of variables of a country
#' @param n Number of countries in the GVAR model
#' @param p Lag specifications
#' @param state A vector specifying state of the GVAR parameters
#'
#' @return The auto regression paramters of the GVAR model
#' @export
BoAoWs2Gs = function(Bo,Ao,W,m,n,p,state) {
  Bo_s=Bo[,,,,1]*0;dim(Bo_s) = c(m,m,p,n)
  Ao_s=Ao[,,,,1]*0;dim(Ao_s) = c(m,m,p,n)
  for (i in 1:n ) {
    Bo_s[,,,i] = Bo[,,,i,state[i]]
    Ao_s[,,,i] = Ao[,,,i,state[i]]
  }
  Go_s = BoAoW2G(Bo_s,Ao_s,W,m,n,p)
  return(Go_s)
}



#' Calculate a state dependent covariance matrix from GVAR residuals
#'
#' @param res An estimated GVAR object
#' @param StateT A vector of selected states of each country
#'
#' @return The covariance matrix of the selected states
#' @export
SigmaNPD = function(res,StateT) {
      m = res$m
      p = res$p
      if (length(p)>1 ) pp = max(p[,1:2,]) else pp = p
      n = length(StateT);
      dim(StateT) = c(1,n);
      resid = res$resid[,,1]*0
      for (i in 1:n)  {
          resid[,(1+(i-1)*m):(m*i)]=res$resid[,(1+(i-1)*m):(m*i),StateT[i]]
      }

      sigmaT = diag(n*m)
	for (i in 1:n)        {
      	for (j in 1:n)  {
                NN = sum((abs(resid[,((i-1)*m+1):(i*m)])>0)*(abs(resid[,((j-1)*m+1):(j*m)])>0))/m-pp*m
                if (NN>5)  sigmaT[((i-1)*m+1):(i*m),((j-1)*m+1):(j*m)] = t(resid[,((i-1)*m+1):(i*m)])%*%resid[,((j-1)*m+1):(j*m)]/NN
                   else    sigmaT[((i-1)*m+1):(i*m),((j-1)*m+1):(j*m)] = t(res$resid[,((i-1)*m+1):(i*m),1]+res$resid[,((i-1)*m+1):(i*m),2])%*%(res$resid[,((j-1)*m+1):(j*m),1]+res$resid[,((j-1)*m+1):(j*m),2])/dim(resid)[1]
      }
      }
      sigmanpd = as.matrix(Matrix::nearPD(sigmaT,conv.tol = 1e-10)[[1]])
      return(sigmanpd)
}




#' Regime specific impulse response functions of MRGVAR(n,p,S)
#'
#' This function calculates the regime specific impulse response functions of an estimated MRGVAR(n,p,S).
#' Using the estimated G\[,,,s\] and Sigma\[,,s\] matrices of the MRGVAR, this function calculated the regime specific impulse response functions.
#' @param res a list of estimated MRGVAR as output of MRGVARest
#' @param nstep the length of impulse response function
#' @param comb a vector specify the concerted action in policy-simulation impulse response function
#' @param state an n vector specifying the specific state for each country.
#' @param  irf  : types of the impulse response irf=c("gen","chol","chol1","gen1","comb1")
#' "gen" for generalized impulse response with one standard deviation impulses, "gen1" for GIRF with one unit impulses, "chol" Cholezky decomposition, "chol1" Cholezky decomposition with unit impulses, "comb1" concerted action with unit impulse.
#' @param G For permanent and transitory decomposition
#' @param smat For explicit structural decomposition of the correlated shocks
#' @param sigmaNPDS the state-dependent covariance matrix
#' @return a list containing the impulse response functions and the accumulated impulse response function, and the boostrap parameters as well.
#' @examples
#'
#' ## case of n = 2, m = 2, S = 2     ## m: number of variables, n: number of countries
#' p = rep(1,12); dim(p) = c(2,3,2)
#' p[1,1,2] = 2; p[2,2,2]=2; p[,3,] = 0
#' TH = c(1:2)*0; dim(TH) = c(1,2)
#' res_d <- MRGVARData(m=2,n=2,p=p,TH=TH,T=100,S=2,SESVI=c(1,3),type="const")
#' max(res_d$Y)
#'
#' ### estimation of the MRGVAR model
#' colnames(res_d$Y) = c("P","Q","Pa","Qa")
#' res_e = MRGVARest(res=res_d)
#' res_e$Summary
#'
#' IRF  = irf_MRGVAR(res=res_e,nstep=10,comb=NA,state=c(1,1),irf="gen1")
#'
#' @export
#'
irf_MRGVAR = function(res=res,state=state,nstep=nstep,comb=comb,irf = c("gen", "chol", "chol1","gen1","genN1", "comb1","smat"),G=NA,smat=NA,sigmaNPDS=NA) {
      if (missing(G)) 		G = NA
      if (missing(sigmaNPDS))	sigmaNPDS=NA
      if (missing(smat))      smat = NA
      neq 	= dim(res$Go)[1]
	    nvar	= dim(res$Go)[2]
      m 	= res$m
	    n 	= res$n
      p 	= res$p
      if (length(p)>1) pp = max(p[,1:2,]) else pp = p
	    Bo    = res$Bo
	    Ao 	= res$Ao
      W	= res$W
      B     = BoAoWs2Gs(Bo,Ao,W,m,n,pp,state)
      #if (anyNA(sigmaNPDS))	sigma = SigmaNPD(res,state)  else  sigma = sigmaNPDS
      if (anyNA(sigmaNPDS))	sigma = SigmaNPDSelectR(res,state)   else  sigma = sigmaNPDS

      response <- array(0,dim=c(neq,nvar,nstep));
      response <- irf_B_sigma(B,sigma,nstep,comb,irf=irf,G=G,smat=smat)
	return(response)
}



#' Generalized impulse response functions of MRGVAR(m,n,p,S) with regime migrations
#'
#' This function calculates the generalized impulse response functions of an estimated MRGVAR(n,p,S) for given a shock vector.
#'
#'                   GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK))
#'
#' @param  res   : an MRGVAR object containing the components of the output of MRGVARData or MRGVARest.
#' @param  shock : an mn-vector containing the shocks as impulse.
#' @param  R     : the number of runs to integrate out the random effects in order to obtain the means (see equation above).
#' @param  nstep : the length of the responses
#' @param  Omega_hist : the initial values from which the simulation runs start.For Omega_hist=NA the most recent values are taken as the initial values. For Omega_hist=0, the initial values are zeros.
#' @param  resid_method : resid_method = c("resid", "parametric"), It generates the random residuals from residuals bootstrap or parametric bootstrap.
#' @return an (mn x mn x nstep+1) matrix of impulse response functions. The rows represent response the columns represent impulses.
#' @examples
#' ## case of n = 2, m = 2, S = 2     ## m: number of variables, n: number of countries
#' p = rep(1,12); dim(p) = c(2,3,2)
#' p[1,1,2] = 2; p[2,2,2]=2; p[,3,] = 0
#' TH = c(1:2)*0; dim(TH) = c(1,2)
#' res_d <- MRGVARData(m=2,n=2,p=p,TH=TH,T=100,S=2,SESVI=c(1,3),type="const")
#' max(res_d$Y)
#'
#' ### estimation of the MRGVAR model
#' colnames(res_d$Y) = c("P","Q","Pa","Qa")
#' res_e = MRGVARest(res=res_d)
#' res_e$Summary
#'
#' IRF_CB  = irf_MRGVAR_CB(res=res_e,nstep=10,comb=NA,state=c(1,1),irf="gen1",runs=20,
#' conf=c(0.05,0.95))
#' IRF_g = IRF_graph(IRF_CB[[1]],Names=c("P","Q","Pa","Qa"))    #IRF
#' IRF_g = IRF_graph(IRF_CB[[2]])   # accumulated IRF
#'
#' GIRF    = girf_MRGVAR_RM(res=res_e,shock=c(1,1,1,1),R=100,nstep=10,Omega_hist=NA,
#' resid_method='parametric')
#' GIRF_CB = girf_MRGVAR_RM_CB(res=res_e,shock=c(1,1,1,1),R=100,nstep=10,Omega_hist=NA,
#' resid_method='parametric',conf=c(0.05,0.95),N=10)
#' IRF_g   = IRF_graph(GIRF_CB,Names=c("P","Q","Pa","Qa"))    #IRF
#' @export
girf_MRGVAR_RM <- function(res,shock,R,nstep,Omega_hist,resid_method) {
  ####  this function generate the impulse response function of MRVAR with regime migration
  ####
  ####                  GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK))
  ####
  ####
  #### RES is the output of MRVAREST1
  n 		  = res$n
  m 		  = res$m
  p 		  = res$p
  d		    = res$d
  S 		  = res$S
  SESVI	  = res$SESVI
  TH		  = res$TH
  Go      = res$Go
  Bo		  = res$Bo
  Ao      = res$Ao
  Co		  = res$Co
  Sigmao  = res$Sigmao
  W       = res$W
  type	  = res$type
  SV      = res$SV
  P       = max(d,p)
  X 		  = res$X; if (!anyNA(X)) X = X[1:(P+nstep+1),];
  Yo		  = res$Yo
  Uo      = res$Uo
  TT      = dim(res$Y)[1]
  d       = res$d
  IRF 	  = list()
  YR      = list()
  YS      = list()
  residR  <-  list()
  residS  <-  residR

  if (anyNA(Omega_hist)) {
    Omega_hist = res$Y[(dim(res$Y)[1]-P+1):dim(res$Y)[1],,drop=FALSE]
  }  else  {
    if (Omega_hist==0) Omega_hist = res$Y[(dim(res$Y)[1]-P+1):dim(res$Y)[1],,drop=FALSE]*0
  }

  state <- as.numeric(Omega_hist[P-d+1,SESVI]>TH) + 1

  #R 		= 400
  #shock       = (1:n)/(1:n)

  DIMresid    = dim(Uo)
  if (length(DIMresid) == 3)  residI = Uo[1:(P+1+nstep),,]*0
  if (length(DIMresid) == 2)  residI = Uo[1:(P+1+nstep),]*0
  shock_mat = matrix(0,m*n,m*n)
  Sigma_s = SigmaNPD(res,state)

  shock_mat = Sigma_s%*%solve(diag(diag(Sigma_s)))%*%diag(shock)

  residI    = array(0,c(P+1+nstep,m*n,m*n,DIMresid[3]))

  residI    = array(0,c(P+1+nstep,m*n,m*n,DIMresid[3]))

  Yo = Omega_hist
  YR = array(0,c(P+1+nstep,m*n,m*n))
  YS <- YR
  MYS <- YS
  MYR <- YS
  GIRF <- array(0,c(m*n,m*n,nstep+1))



  for (i in 1:R) {
    if ( resid_method=="resid" ) {

      if (length(DIMresid) == 3)  {
        for (k in 1:(m*n)) {
          residI[,,k,]  = Uo[NoutofT(P+nstep+1,DIMresid[1]),,];
        }
      }
    }
    if ( resid_method=="parametric" ) {

      if (length(DIMresid) == 3)  {
        for (k in 1:(m*n)) {
          for (s in 1:S)   residI[,,k,s] = rnormSIGMA(P+nstep+1,Sigmao[,,s])
        }
      }
    }
    #residI[P+1,,,] = 0
    residR[[i]] <- residI
    for (s in 1:S) {residI[P+1,,,s] <-shock_mat}    # such that residI[P+1,,k,] corresponds to the k-th column shock in shock_mat for both regimes (budui dan ok)
    residS[[i]] = residI

    for (k in 1:DIMresid[2]) {
      YR[,,k] = MRGVARData(m=m,n=n,p=p,T=(P+nstep+1),S=S,W=W,SESVI=SESVI,TH=TH,Go=Go,Ao=Ao,Bo=Bo,Sigmao=Sigmao,Uo=residR[[i]][,,k,],SV=SV,type=type,Co=Co,X=X,Yo=Yo,d=d)$Y
      YS[,,k] = MRGVARData(m=m,n=n,p=p,T=(P+nstep+1),S=S,W=W,SESVI=SESVI,TH=TH,Go=Go,Ao=Ao,Bo=Bo,Sigmao=Sigmao,Uo=residS[[i]][,,k,],SV=SV,type=type,Co=Co,X=X,Yo=Yo,d=d)$Y
      ## to delete YS[[i]] = MRGVARData(m=m,n=n,p=p,T=(P+nstep+1),S=S,W=W,SESVI=SESVI,TH=TH,Go=Go,Ao=Ao,Bo=Bo,Sigmao=Sigmao,Uo=residS[[i]],SV=SV,type=type,Co=Co,X=X,Yo=Yo*0,d=d)
    }

    MYR = MYR + 1/R*YR
    MYS = MYS + 1/R*YS

  }  # end of R loop
  ############## integrating out the disturbances
  for (i in 1:(m*n) ) {
    for (j in 1:(m*n)) {
      for (k in 1:(nstep+1)) GIRF[i,j,k] = MYS[P+k,i,j]-MYR[P+k,i,j]
    }
  }
  return(list(GIRF,MYS[(P+1):(P+1+nstep),,],MYR[(P+1):(P+1+nstep),,]))
}


#' Generalized impulse response functions of MRGVAR(m,p,n,S,W,TH) with regime migrations
#'
#' This function calculates the generalized impulse response functions of an estimated MRVAR(n,p,S) for given a shock vector.
#'
#'                   GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK))
#'
#' It also generates the bootstrapped confidence intervals.
#'
#' @param  res   : a MRGVAR object containing the components of the output of MRGVARData or MRGVARest.
#' @param  shock : an mn-vector containing the shocks as impulse.
#' @param  R     : the number of runs to integrate out the random effects in order to obtain the means (see equation above).
#' @param  nstep : the length of the responses
#' @param  Omega_hist : the initial values from which the simulation runs of impulse and response functions start
#' @param  resid_method : resid_method = c("resid", "parametric"), It generate random residuals either from residuals bootstrap or parametric bootstrap.
#' @param  conf  : a two component vector containing the tail probabilities of the bootstrap confidence interval.
#' @param  N     : number of bootstrapping runs
#' @return an (n x n x nstep+1 x 3) array containing of impulse response functions with lower and upper confidence bonds. The rows represent response the columns represent impulses.
#' @examples
#' ## case of n = 2, m = 2, S = 2     ## m: number of variables, n: number of countries
#' p = rep(1,12); dim(p) = c(2,3,2)
#' p[1,1,2] = 2; p[2,2,2]=2; p[,3,] = 0
#' TH = c(1:2)*0; dim(TH) = c(1,2)
#' res_d <- MRGVARData(m=2,n=2,p=p,TH=TH,T=100,S=2,SESVI=c(1,3),type="const")
#' max(res_d$Y)
#'
#' ### estimation of the MRGVAR model
#' colnames(res_d$Y) = c("P","Q","Pa","Qa")
#' res_e = MRGVARest(res=res_d)
#' res_e$Summary
#'
#' Mod(STAT(res_e$Go[,,,1]))
#' Mod(STAT(res_e$Go[,,,2]))
#'
#' GIRF    = girf_MRGVAR_RM(res=res_e,shock=c(1,1,1,1),R=100,nstep=10,Omega_hist=NA,
#' resid_method='parametric')
#' GIRF_CB = girf_MRGVAR_RM_CB(res=res_e,shock=c(1,1,1,1),R=100,nstep=10,Omega_hist=NA,
#' resid_method='parametric',conf=c(0.05,0.95),N=10)
#' GIRF_g  = IRF_graph(GIRF_CB,Names=c("P","Q","Pa","Qa"))
#'
#' @export
girf_MRGVAR_RM_CB <- function(res,shock,R,nstep,Omega_hist=NA,resid_method="parametric",conf,N) {
  ##### this is to generate bootstrap GIRF for MRVAR with regime migrations
  #####
  n 		= res$n
  m 		= res$m
  p 		= res$p
  S 		= res$S
  SESVI		= res$SESVI
  TH		= res$TH
  Go          = res$Go
  Bo		= res$Bo
  Ao          = res$Ao
  Co		= res$Co
  Sigmao	= res$Sigmao
  W           = res$W
  type		= res$type
  SV          = res$SV
  X 		  = res$X;
  Yo		  = res$Yo
  d           = res$d
  Uo          = res$Uo
  T           = dim(res$Y)[1]
  P           = max(d,p)
  GIRF 		= (1:(m*n*n*m*(nstep+1)*N))*0; dim(GIRF) = c(m*n,m*n,nstep+1,N)

  GIRFBd	= (1:(m*n*n*m*(nstep+1)*(length(conf)+1)))*0; dim(GIRFBd) = c(m*n,m*n,nstep+1,length(conf)+1)

  GIRFBd[,,,1]= girf_MRGVAR_RM(res,shock,R=2*R,nstep,Omega_hist,resid_method)[[1]]

  for (i in 1:N) {
    #res_run = MRGVARData(m=m,n=n,p=p,T=T,S=S,W=W,SESVI=SESVI,TH=TH,Go=Go,Ao=Ao,Bo=Bo,Sigmao=Sigmao,Uo=NA,SV=SV,type=type,Co=Co,X=X,Yo=Yo,d=d)
    #res_e   = MRGVARest(res_run)
    res_run = MRGVARDataR(res)
    if (length(colnames(res_run$Y))==0) colnames(res_run$Y) = paste("res_runY",1:ncol(res_run$Y),sep="")
    res_erun   = MRGVARest(res_run)
    RF3     = girf_MRGVAR_RM(res=res_erun,shock,R,nstep,Omega_hist,resid_method)[[1]]
    GIRF[,,,i]  = RF3
  }


  for (tt in 1:(nstep+1) ) {
    for (i in 1:(n*m))           {
      for (j in 1:(n*m))     {
        GIRFBd[i,j,tt,-1] = stats::quantile(GIRF[i,j,tt,], conf)

      }
    }
  }
  return(GIRFBd)
}


#' Calculation of information criteria for a MRGVAR model
#'
#' @param res  : a MRGVAR object obtained from MRGVARData or estimated from MRGVARest.
#' @param I    : index of the country under investigation.
#' @param L_V  : a four components vector containing the maxims of the domestic lag and the foreign lag for each regime, respectively.
#' @param TH_V : a vector containing possible threshold values.
#' @return     a matrix with different lag specifications and the corresponding values of the model selection criteria.
#' @examples
#'
#'
## case of n = 5, m =3, S = 2 #m: number of variables, n: number of countries
#' p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' p[,2,] = 1; p[1,1,1] = 1; p[3,1,1] = 1; p[2,1,2] = 1
#'
#' TH = c(1:5)*0; dim(TH) = c(1,5)
#' res_d <- MRGVARData(m=3,n=5,p=p,TH =TH,T=400,S=2,SESVI=((1:5)*3-2))
#' max(abs(res_d$Y))		# to make sure it is not explosive
#' colnames(res_d$Y) = c("Y11","Y21","Y31","Y12","Y22","Y32","Y13","Y23","Y33","Y14",
#' "Y24","Y34","Y15","Y25","Y35")
#' res_e = MRGVARest(res=res_d)
#'
#' ### four numbers for the maxima lag length in country I:
#' ### regime 1: domestic foreign regime 2: domestic and foreign
#' L_v  = c(3,3,3,3)
#' ### a vector containing possible threshold values
#' TH_v = c(-0.1, -0.05, 0,0.05, 0.1  )
#' CC = MRGVAR_Select(res=res_d,I=1,L_V=L_v,TH_V=TH_v)
#'  CCC = CC[[1]]
#'
#' CCC[which.min(CCC[,9]),]
#' CCC[which.min(CCC[,18]),]
#'
#' CC = MRGVAR_Select(res=res_d,I=2,L_V=L_v,TH_V=TH_v)
#' CCC = CC[[1]]
#'
#' CCC[which.min(CCC[,9]),]
#' CCC[which.min(CCC[,18]),]
#'
#' @export
#'
MRGVAR_Select <- function(res,I,L_V,TH_V) {
  ##    res   object of MRGVARData
  ##    I     index of the equation under investigation, I =1 we investigate the first equation, I = 10 we investigate the 10th equation.
  ##    L_V   maximum of the lags of regimes (4 4 4 4) Regime 1 domestic lags 4   foreign lags 4   regime 2 lags 4  4
  ##    TH_V  vector of threshold values that will be estimated for each chosen lag combinations
  m= res$m
  n= res$n
  p= res$p
  Y= as.matrix(res$Y)
  X     = res$X
  W = res$W
  type  = res$type
  TH    = res$TH
  SESVI = res$SESVI
  Go= res$Go
  S     = res$S
  Ao= res$Ao
  Bo    = res$Bo
  Co    = res$Co
  GDC   = res$GDC
  Sigmao= res$Sigmao*0
  d     = res$d
  P     = max(p,d)
  Pmax  = max(p[,1:2,])
  SigmaS = (1:(n*m*S*n*m*S))*0
  dim(SigmaS) = c(n*m*S,n*m*S)
  VAR_domestic = list()

  foreignLagSame = 1


  k1 = p[I,3,1]
  k2 = p[I,3,2]
  kmax  = max(p[,3,])

  T= dim(Y)[1]
  FY= Y%*%t(W%x%diag(m))
  resid = (1:(T*m*n*S))*0; dim(resid) = c(T,m,n,S)
  ##PS  = (1:(T*n*S))*0;   dim(PS)    = c(T,n,S)
  ###   CC is the single country coefficients of determinsitic components i.e. the foreign variables and the deterministic variables
  ###   It is m x (m*p+k+1) where m*p is the number of foreign variables and k is the number of deterministic variables. with zeros for "none"

  P     = max(L_V,d)
  P_candidate = array(0,c(L_V[1],L_V[2],L_V[3],L_V[4]))
  Criteria   = matrix(0,prod(L_V),4+length(TH_V))
  Criteria_b = matrix(0,prod(L_V),4+length(TH_V))
  Criteria_c = matrix(0,prod(L_V),4+length(TH_V))
  Criteria_d = matrix(0,prod(L_V),4+length(TH_V))

  if (type=="none" |type=="const")  k = 1;
  if (type=="exog0"|type=="exog1")  k = dim(X)[2]+1;
  CC  =  c(1:(m*(Pmax*m+k)*S))*0
  dim(CC) = c(m,m*Pmax+k,S)

  #res_MSVAR = MSVARData(n=m,p=p,T=T,S=S,Co=CC,TM=TM,type="exog",X=Z)
  #res_MSVAR = MSVARData(n=m,p=p,T=T,S=S)

  EMCVG       = c(1:n)*0

  Pppp = matrix(0,3,2)

  for (h in 1:length(TH_V)) {
    idx = 0
    for (l_1 in 1:L_V[1])             {
      for (l_2 in 1:L_V[2])       {
        for (l_3 in 1:L_V[3])  {
          for (l_4 in 1:L_V[4]){
            if ( foreignLagSame == 1) l_4=l_2
            Pppp[,1] = c(l_1,l_2,k1);
            Pppp[,2] = c(l_3,l_4,k2);

            idx = idx+1;
            for (i in I:I) {
              p[i,,]      = Pppp
              Pi          = max(p[i,1:2,])
              CCi  = c(1:(m*(Pi*m+k)*S))*0
              dim(CCi) = c(m,m*Pi+k,S)
              Z         = matrix(0,T,Pi*m)
              FYp    = stats::embed(FY[,(m*(i-1)+1):(i*m)],(Pi+1))
              Z[(Pi+1):T,]= FYp[,(m+1):(Pi*m+m)]
              kz          = dim(Z)[2]
              if (type=="none")   {Z = rep(Z,S); dim(Z) = c(T,kz,S); ZZ = array(0,c(T,kz,S)); for (s in 1:S) ZZ[,1:kz,s]=Z[,,s]}
              if (type=="const")  {Z = rep(Z,S); dim(Z) = c(T,kz,S); ZZ = array(0,c(T,kz,S)); for (s in 1:S) ZZ[,1:kz,s]=Z[,,s]}
              if (type=="exog0")  {Z = rep(Z,S); dim(Z) = c(T,kz,S); ZZ = array(0,c(T,kz+kmax,S)); for (s in 1:S) ZZ[,1:(kz+p[i,3,s]),s] = cbind(Z[,,s],X[,1:p[i,3,s],i,s])}
              if (type=="exog1")  {Z = rep(Z,S); dim(Z) = c(T,kz,S); ZZ = array(0,c(T,kz+kmax,S)); for (s in 1:S) ZZ[,1:(kz+p[i,3,s]),s] = cbind(Z[,,s],X[,1:p[i,3,s],i,s])}
              if (type=="none"|type=="exog0")  type_holder = "exog0"
              if (type=="const"|type=="exog1") type_holder = "exog1"
              XX = (1:(T*dim(ZZ)[2]*S))*0; dim(XX) = c(T,dim(ZZ)[2],S)
              if (!anyNA(X)) for (s in 1:S) XX[,1:(m*p[i,2,s]+p[i,3,s]),s] = ZZ[,c(1:(m*p[i,2,s]),(m*Pi+1):(m*Pi+p[i,3,s])),s]
              if (anyNA(X))  for (s in 1:S) XX[,1:(m*p[i,2,s]),s] = ZZ[,c(1:(m*p[i,2,s])),s]

              ### put data into a MRVARData object in order to replace the data and order parameters
              pp = t(p[i,,]); pp[,2] = pp[,2]*m; ppp = pp[,1:2]; ppp[,2] = pp[,2]+pp[,3]
              if (SESVI[i]/m > SESVI[i]%/%m)   SESVIi = SESVI[i]- (SESVI[i]%/%m)*m  else SESVIi = m

              if (max(ppp[,1])>=max(ppp[,2])%/%m) {
                res_MRVAR =  MRVARData(n=m,p=ppp,T=T,S=S,Co=NA,SESVI=SESVIi,type=type_holder,X=XX)
                res_MRVAR$Y       =  Y[,(m*(i-1)+1):(i*m)]
                res_MRVAR$X       =  XX
                res_MRVAR$SV      =  res_MRVAR$Y[,SESVIi]
                res_MRVAR$TH      =  TH_V[h]

              }  else  {

                startT = max(ppp[,2])%/%m - max(ppp[,1])+1
                res_MRVAR =  MRVARData(n=m,p=ppp,T=(T-startT+1),S=S,Co=NA,SESVI=SESVIi,type=type_holder,X=XX[startT:T,,])
                res_MRVAR$Y       =  Y[startT:T,(m*(i-1)+1):(i*m)]
                res_MRVAR$X       =  XX[startT:T,,]
                res_MRVAR$SV      =  res_MRVAR$Y[,SESVIi]
                res_MRVAR$TH      =  TH_V[h]
              }
              RR=  MRVARest(res=res_MRVAR)
              if (length(RR)>1) {
                Criteria[idx,c(c(1:4),4+h)]     = c(l_1,l_2,l_3,l_4,RR$Summary$LH_AIC)
                Criteria_b[idx,c(c(1:4),4+h)]   = c(l_1,l_2,l_3,l_4,RR$Summary$LH_BIC)
                Criteria_c[idx,c(c(1:4),4+h)]   = c(l_1,l_2,l_3,l_4,RR$Summary$LH_P)
                Criteria_d[idx,c(c(1:4),4+h)]   = c(l_1,l_2,l_3,l_4,RR$Summary$LH_N)
              }


              Criterion = cbind(Criteria,Criteria_b,Criteria_c,Criteria_d)

            }
          }
        }
      }
    }
  }
  return(list(Criterion,RR))
}

#' Regime specific impulse response functions of MRGVAR(n,p,S)
#'
#' This function calculates the regime specific impulse response functions of an estimated MRGVAR(n,p,S).
#' Using the estimated G\[,,,s\] and Sigma\[,,s\] matrices of the MRGVAR, this function calculated the regime specific impulse response functions.
#' @param res a list of estimated MRGVAR as output of MRGVARest
#' @param state an n vector specifying the specific state for each country.
#' @param nstep the length of impulse response function
#' @param comb a vector specify the concerted action in policy-simulation impulse response function
#' @param  irf  : types of the impulse response irf=c("gen","chol","chol1","gen1","comb1")
#' "gen" for generalized impulse response with one standard deviation impulses, "gen1" for GIRF with one unit impulses, "chol" Cholezky decomposition, "chol1" Cholezky decomposition with unit impulses, "comb1" concerted action with unit impulse.
#' @param G For permanent and transitory decomposition
#' @param smat For an explicit structural decomposition of the correlated shocks
#' @param sigmaNPDS the state-dependent covariance matrix
#' @param runs number of bootstrapping runs
#' @param conf A vector containing confidence levels
#' @param NT number of impulse response scenarios in a simulation run
#' @return a list of bootstrap result. The first component contains the impulse response functions with confidence bands. It is an (mn,mn,nstep,3) array where the IRF columns represent the impulse and rows represent the responses.
#' @examples
#' ## case of n = 2, m = 2, S = 2     ## m: number of variables, n: number of countries
#' p = rep(1,12); dim(p) = c(2,3,2)
#' p[1,1,2] = 2; p[2,2,2]=2; p[,3,] = 0
#' TH = c(1:2)*0; dim(TH) = c(1,2)
#' res_d <- MRGVARData(m=2,n=2,p=p,TH=TH,T=100,S=2,SESVI=c(1,3),type="const")
#' max(res_d$Y)
#'
#' ### estimation of the MRGVAR model
#' colnames(res_d$Y) = c("P","Q","Pa","Qa")
#' res_e = MRGVARest(res=res_d)
#' res_e$Summary
#'
#' IRF_CB  = irf_MRGVAR_CB(res=res_e,nstep=10,comb=NA,state=c(1,1),irf="gen1",runs=20,
#' conf=c(0.05,0.95))
#' IRF_g = IRF_graph(IRF_CB[[1]],Names=c("P","Q","Pa","Qa"))    #IRF
#' #IRF_g = IRF_graph(IRF_CB[[2]])   # accumulated IRF
#'
#' @export
irf_MRGVAR_CB = function (res, state = c(2, 1), nstep, comb, irf = c("gen", "chol", "chol1", "gen1", "comb1"), G=NA,smat=NA,sigmaNPDS=NA,runs = 200, conf = c(0.05, 0.95),NT = 1)
{
    m = res$m
    n = res$n
    p = res$p
    if (length(p) > 1)
        pp = max(p[, 1:2, ])
    else pp = p
    T = dim(res$Y)[1] * NT
    W = res$W
    S = res$S
    SV = res$SV
    SESVI = res$SESVI
    TH = res$TH
    Ao = res$Ao
    Bo = res$Bo
    Co = res$Co
    Go = res$Go
    GoColect = array(0, c(dim(Go), runs))
    BoColect = array(0, c(dim(Bo), runs))
    AoColect = array(0, c(dim(Ao), runs))
    UoColect = array(0, c(T, m * n, S, runs))
    YColect = array(0, c(T, m * n, runs))
    type = res$type
    X = res$X
    mu = res$mu
    B = BoAoWs2Gs(Bo, Ao, W, m, n, pp, state)
    neq = dim(B)[1]
    nvar = dim(B)[2]
    sigma = res$Sigmao
    response <- array(0, dim = c(neq, nvar, nstep, length(conf) + 1))
    accresponse <- response
    sigmaNPDS = SigmaNPD(res, state)
    response[, , , 1] <- irf_MRGVAR(res=res, nstep=nstep, comb=comb, state=state, irf=irf)
    accresponse[, , , 1] <- response[, , , 1]
    for (tt in 2:nstep) {  accresponse[, ,tt, 1] = accresponse[, ,tt-1, 1] + response[, ,tt, 1] }

    responseR <- array(0, dim = c(neq, nvar, nstep, runs))
    accresponseR <- responseR
    Uo_run = array(0, c(T, n * m, S))
    for (i in 1:runs) {
        for (s in 1:S) Uo_run[, , s] = rnormSIGMA(T, res$Sigmao[,, s])
        if (length(p) > 1) {
            res_run = MRGVARDataR(res)
            if (length(colnames(res_run$Y))==0) colnames(res_run$Y) = paste("res_runY",1:ncol(res_run$Y),sep="")
            res_e   = MRGVARest(res_run)
        }

        responseR[, , , i] <- irf_MRGVAR(res=res_e,  state=state,nstep=nstep, comb=comb, irf=irf, G=G,smat=smat,sigmaNPDS=sigmaNPDS)
        accresponseR[, , , i] <-  responseR[, , , i]
        for (tt in 2:nstep) {  accresponseR[, ,tt, i] = accresponseR[, ,tt-1, i] + responseR[, ,tt, i] }
        GoColect[, , , , i] <- res_e$Go
        BoColect[, , , , , i] <- res_e$Bo
        AoColect[, , , , , i] <- res_e$Ao
        #UoColect[, , , i] <- res_run$Uo
        #YColect[, , i] <- res_run$Y
    }
    responseR[, , , 1]    = response[, , , 1]
    accresponseR[, , , 1] = accresponse[, , , 1]

    for (tt in 1:(nstep)) {
        for (i in 1:neq) {
            for (j in 1:nvar) {
                response[i, j, tt, -1]    = stats::quantile(responseR[i, j, tt, ], conf)
 		    accresponse[i, j, tt, -1] = stats::quantile(accresponseR[i, j, tt, ], conf)
            }
        }
    }
    return(list(response, accresponse, GoColect, BoColect, AoColect, UoColect, YColect,responseR))
}

#' Data generating process of MRGVARDataR(res)
#'
#' This function will generate data from a MRGVAR object. It will generate enough data for estimation purpose.
#'
#' @param res     : an output of MRGVARest
#' @return	: an MRGVAR object.
#'
#' @export
MRGVARDataR=function(res) {
### res_e : an estimated MRGVAR model that is an output of MRGVARest
### T     : number of observations
### Remarks: MRGVARDataR is used in bootstrapping to generate sufficient observations such that in the regime of rare occurrence also contains
###          sufficient observations for estimation purposed
###
###  MINH:   minimum number of observations of a nation in the rare occurrence regime.

    m = res$m
    n = res$n
    p = res$p
    if (length(p) > 1)
        pp = max(p[, 1:2, ])
    else pp = p
    T = dim(res$Y)[1]
    W = res$W
    S = res$S
    SV = res$SV
    SESVI = res$SESVI
    TH = res$TH
    Ao = res$Ao
    Bo = res$Bo
    Co = res$Co
    Go = res$Go
    type = res$type
    X = res$X
    mu = res$mu
    MINH = T
    repeat {
	 Uo_run = array(0, c(T, n * m, S))
       for (s in 1:S) Uo_run[, , s] = rnormSIGMA(T, res$Sigmao[,,s])       #### zhe ge you wenti
    	 res_run = MRGVARData(m, n, p, T, S, W, SESVI, TH, res$Go, Ao = NA, Bo = NA, Sigmao = NA, Uo = Uo_run, SV, type, Co, X)
       for (i in 1:n) { MINH = min(MINH,sum(res_run$Y[,SESVI[i]]>TH[1,i]),sum(res_run$Y[,SESVI[i]]<TH[1,i]))}
       if (MINH > (max(p)*(m*2)+2*m)) break
          else {T=T+T; print(c(MINH,T)); MINH=T}
    }
  return(res_run)
}


#' Calculation of the cumulative impulse response function from a impulse response function
#'
#' @param IRF An impulse response function
#'
#' @return The cumulative impulse response function
#' @export
#'
ACCIRFconfR = function(IRF) {
      ACCirf = IRF
        dm = dim(IRF)
      for (t in 2:dm[3])          {
                 ACCirf[,,t,] =    ACCirf[,,t-1,] + IRF[,,t,]
      }
      return(ACCirf)
}



#' Select a state-dependent covariance matrix from a GVAR model
#'
#' @param res An estimated GVAR object
#' @param StateT A vector of selected regimes of each country
#'
#' @return A positive definite covariance matrix of the selected states each country
#' @export
SigmaNPDSelectR = function(res,StateT) {
      n = length(StateT);
      dim(StateT) = c(1,n);
      Ranking = 1:n;
      SigmaS = res$SigmaS
      N = dim(SigmaS)[1]/(2*n)
      sigmaT = diag(n*N)
        for (i in 1:n)        {
        for (j in 1:n)  {
            sigmaT[((i-1)*N+1):(i*N),((j-1)*N+1):(j*N)] = SigmaS[((StateT[1,i]-1)*n*N+(i-1)*N+1):((StateT[1,i]-1)*n*N+i*N),((StateT[1,j]-1)*n*N+(j-1)*N+1):((StateT[1,j]-1)*n*N+j*N)]
      }
      }
      #if (anyNA(sigmaT)) { sigmanpd = SigmaRKSelect(Model,StateT,Ranking)[[1]] }
      #else {
      #    sigmanpd = as.matrix(Matrix::nearPD(sigmaT,conv.tol = 1e-10)[[1]])
      #}
      sigmanpd = as.matrix(Matrix::nearPD(sigmaT,conv.tol = 1e-10)[[1]])
      return(sigmanpd)
}


#' Global Responses
#'
#'This function calculates global or regional responses from a set of bootstrapped impulse response functions and a weighting matrix for the aggregation of the global responses.
#'
#' @param IRF_CB an output of irf_MRGVAR_CB
#' @param comb_all a weighting matrix for the aggregation of the global responses
#'
#' @return a list containing the global impulse response functions and the accumulated global impulse response functions.
#' @export
#'
irf_GloabalResponse_CB <- function(IRF_CB,comb_all) {

  bootarray <-IRF_CB[[8]]
  DIM <- dim(bootarray)
  #bootarray = IRFC[[8]]
  glresp <-dim(comb_all)[2]
  nN     <-DIM[1]
  nstep  <-DIM[3]
  nrun   <-DIM[4]
  globalbootarray <-array(0,dim=c(glresp,nN,nstep,nrun))

  for (m in 1:glresp)     {
    for (j in 1:(nN))    {
      for (k in 1:nstep ) {
        for (q in  1:nrun) {
          globalbootarray[m,j,k,q] = sum(bootarray[,j,k,q]*comb_all[,m])

        }}}}

  globalbootconf = c(1:(3*glresp*nN*nstep))*0
  dim(globalbootconf) = c(glresp,nN,nstep,3)
  globalbootconf[,,,1] = globalbootarray[,,,1]

  for (i in 1:glresp)       {
    for (j in 1:(nN))     {
      for (k in 1:nstep ) {
        globalbootconf[i,j,k,2:3] = stats::quantile(globalbootarray[i,j,k,],c(.05,.95))
      }}}
  dim(globalbootconf)

  ACCbootconf3N = ACCIRFconfR(globalbootconf)
  return(list(globalbootconf,ACCbootconf3N))
}



