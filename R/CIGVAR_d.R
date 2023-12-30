#' Data generating process of CIGVAR(m,n,p)
#'
#' This function generates data from a cointegrated global VAR process and returns an CIGVAR(m,n,p) object that is a list containing data and parameters used in the CIGVAR(m,n,p) process.
#'
#' @param m     : number of variables in each country/unit
#' @param n     : number of countries/units
#' @param p     : an (n x 3) matrix, each raw specifies the lag length of the domestic variables, the foreign variables and the number of the exogenous variables.
#' @param T     : number of observations.
#'
#'              (m,n,p,T) are parameters which must be provided.
#' @param W     : an (n x n) weighting matrix. w_ij is the weight of foreign country j in the foreign variables of i-th country diag(W)=0
#' @param r_npo : an (m, p, n) array collecting the roots of the characteristic functions in the lag operator of the country VAR. The number of ones in each ith row of r_npo is the number of unit roots i-th country/unit.
#' @param Ao    : an (m, m, p, n) array collecting the  coefficients of foreign variables.
#' @param Bo    : an (m, m, p, n) array collecting the  coefficients of domestic variables.
#' @param Co    : an (m , k+1, n) array collecting the coefficients of the deterministic components of the n countries.
#' @param Uo    : an (T x mn) matrix of the temporally independent innovation processes
#' @param Sigmao : (mn x mn) matrix of the covariance matrix of the CIGVAR(m,n,p)
#' (W,r_npo,Ao,Bo,Uo,Sigmao) if not provided, they will be generated randomly. The default assumption is one unit root in one country. Hence m-1 cointegration relations in each country.
#'
#' @param type	: deterministic component "const", "none", "exog0", and "exog1" are four options
#' @param X	: a (T x k x n) array of exogenous stationary variables.
#' @param mu    : if type = "const" mu has the same dimension as Co. is an muV is nm vector of the means of the time series in the system
#' @param DFYflag	: DFYflag = 0 implies foreign variables are not in the cointegration space. DFYflag = 1 allows foreign variables enter the cointegration relation.
#' @param crk   : n vector containing the cointegration rank in each country/unit. crk is to specified for estimation of parameters.
#' @param Ncommfakt   : number of common exogenous stochastic factors
#' @param A           : a transformation matrix that adds the common exogenous factors to the system.
#' @param uz          : innovations of the exogenous common factors
#' @return a CIGVAR object containing the generated data, the parameters used and the exogenous variables.
#'
#' @examples
#' n = 5
#' p = (1:15)*0; dim(p) = c(5,3)
#' p[,1] = 3; p[,2]=3;
#' res_d = CIGVARData(m=4,n=5,p=p,T=100,type="const")
#' max(abs(res_d$Y))
#' plot(ts(res_d$Y[,1:10]))
#' res_d$r_npo
#' STAT(res_d$G)
#' res_e = CIGVARest(res_d)
#' res_e$Summary
#' STAT(res_d$G)
#' plot(ts(res_d$Y[,1:10]))
#'
#' @export
CIGVARData <- function(m, n, p, T, W = NA, r_npo = NA, Ao = NA, Bo = NA, Co = NA, Uo = NA, Sigmao = NA, type = NA, X = NA, mu = NA, DFYflag = NA,  crk = NA, Ncommfakt=1,A=NA,uz=NA)
{

  if (missing(W)) {
    W = NA
  }
  if (missing(r_npo)) {
    maxPd = max(p[, 1])
    r_npo = c(1:(m * maxPd * n))/c(1:(m * maxPd * n)) * 4
    dim(r_npo) = c(m, maxPd, n)
    r_npo[1, 1, ] = 1
  }

  if (missing(Bo)) {
    Bo = NA
  }
  if (missing(Bo)) {
    Ao = NA
  }
  if (missing(Co)) {
    Co = NA
  }

  if (missing(Uo)) {
    Uo = NA
  }
  if (missing(Sigmao)) {
    Sigmao = NA
  }


  if (missing(type)) {
    type = NA
  }
  if (missing(X)) {
    X = NA
  }

  if (missing(mu)) {
    mu = NA
  }
  if (missing(DFYflag)) {
    DFYflag = 0
  }

  if (missing(crk)) {
    crk = (1:n)/(1:n)
  }
  if (missing(Ncommfakt)) {
    Ncommfakt = 0
  }


  alpha = list()
  beta = list()

  Pmax = max(p[, 1:2])
  P = max(p)
  if (!anyNA(X))
    k = dim(X)[2]
  if (anyNA(Bo)) {
    Bo = (1:(m * m * Pmax * n)) * 0
    dim(Bo) = c(m, m, Pmax, n)
    for (i in 1:n) {
      r_np = c(1:(m * p[i, 1])) * 0
      dim(r_np) = c(m, p[i, 1])
      r_np = r_npo[, 1:p[i, 1], i]
      VARD = VARData(m, p[i, 1], T, r_np = r_np)
      Bo[, , 1:p[i, 1], i] = VARD$B
      alpha[[i]] = B2CIB(VARD$B)[[2]]
      beta[[i]] = B2CIB(VARD$B)[[3]]
    }
  }
  if (anyNA(Sigmao)) {
    Sigmao = matrix(0, m * n, m * n)
    VARD = VARData(m * n, p[1, 1], T)
    Sigmao = VARD$Sigma
  }
  if (anyNA(Uo)) {
    Uo = rnormSIGMA(T, Sigmao)
  }
  if (anyNA(W)) {
    W = matrix(((1:(n * n))/(1:(n * n)) * 1/(n - 1)), n,
               n)
    for (i in 1:n) W[i, i] = 0
  }
  if (anyNA(Ao)) {
    Ao = (1:(m * m * Pmax * n)) * 0
    dim(Ao) = c(m, m, Pmax, n)
    for (i in 1:n) {
      if (p[i, 2] < 2)
        Ao = Ao
      if (p[i, 2] >= 2) {
        VARD = VARData(m, p = (p[i, 2] - 1), T)
        BB = VARD$B/1
        Ao[, , 1, i] = BB[, , 1] + alpha[[i]] %*% t(beta[[i]]) *DFYflag
        Ao[, , p[i, 2], i] = -BB[, , p[i, 2] - 1]
        if ((p[i, 2] - 1) >= 2) {
          for (L in 2:(p[i, 2] - 1)) Ao[, , L, i] = BB[,
                                                       , L] - BB[, , L - 1]
        }
      }
    }
  }
  if (anyNA(type)) {
    type = "none"
  }
  if (type == "none") {
    Co = matrix(0, m, n)
    dim(Co) = c(m, 1, n)
    mu = matrix(0, m, n)
    dim(mu) = c(m, 1, n)
  }
  G = (1:(n * m * n * m * Pmax)) * 0
  dim(G) = c(n * m, n * m, Pmax)
  for (i in 1:n) {
    for (j in 1:n) {
      for (L in 1:Pmax) G[(1 + (i - 1) * m):(i * m), (1 +
                                                        (j - 1) * m):(j * m), L] = Ao[, , L, i] * W[i,
                                                                                                    j]
    }
  }
  for (i in 1:n) {
    for (L in 1:Pmax) G[(1 + (i - 1) * m):(i * m), (1 + (i -
                                                           1) * m):(i * m), L] = Bo[, , L, i]
  }
  Ct = Uo * 0
  if (type == "const") {
    if (anyNA(mu)) {
      mu = matrix(stats::rnorm(n * m), m, n)
      dim(mu) = c(m, 1, n)
    }
    if (anyNA(Co)) {
      Co = mu
      muV = as.vector(mu)
      CoV = muV
      for (L in 1:Pmax) CoV = CoV - G[, , L] %*% muV
      Co = CoV
      dim(Co) = c(m, 1, n)
    }
    else {
      H = diag(n * m)
      for (L in 1:Pmax) H = H - G[, , L]
      CoV = as.vector(Co)
      if (min(abs(eigen(H)$values)) > 0.001) {
        muV = solve(H) %*% CoV
        mu = muV
        dim(mu) = c(m, 1, n)
      }
      else {
        mu = NA
      }
    }
    Ct = matrix(1, T, 1) %*% t(CoV)
  }
  if (type == "exog0") {
    if (anyNA(Co)) {
      Co = matrix(stats::rnorm(m * n * (k + 1)), m * (k +
                                                        1), n)
      dim(Co) = c(m, k + 1, n)
      Co[, 1, ] = (1:m) * 0
      for (i in 1:n) if (p[i, 3] < k)
        Co[, (p[i, 3] + 2):(k + 1), i] = Co[, (p[i, 3] +
                                                 2):(k + 1), i] * 0
    }
    DMCo = dim(Co[, -1, ])
    if (length(DMCo) < 3)
      DMCo = c(dim(Co)[1], 1, dim(Co)[3])
    if (!(DMCo[2] == dim(X)[2]) | !(DMCo[1] == m) | !(DMCo[3] ==
                                                      n)) {
      print("dimension problem")
      return("dimension")
    }
    CoV = matrix(0, dim(X)[2], m * n)
    for (i in 1:dim(X)[2]) CoV[i, ] = as.vector(Co[, 1 +
                                                     i, ])
    Ct = matrix(0, dim(X)[1], m * n)
    for (i in 1:n) Ct[, ((i - 1) * m + 1):((i - 1) * m +
                                             m)] = as.matrix(X[, , i]) %*% CoV[, ((i - 1) * m +
                                                                                    1):((i - 1) * m + m)]
    mu = NA
  }
  else {
    if (type == "exog1") {
      if (anyNA(Co)) {
        Co = matrix(stats::rnorm(m * n * (k + 1)), m *
                      (k + 1), n)
        dim(Co) = c(m, k + 1, n)
        for (i in 1:n) if (p[i, 3] < k)
          Co[, (p[i, 3] + 2):(k + 1), i] = Co[, (p[i,
                                                   3] + 2):(k + 1), i] * 0
      }
      DMCo = dim(Co[, -1, ])
      if (!(DMCo[2] == dim(X)[2]) | !(DMCo[1] == m) | !(DMCo[3] ==
                                                        n)) {
        print("dimension problem")
        return("dimension")
      }
      CoV = matrix(0, dim(X)[2] + 1, m * n)
      for (i in 1:(1 + dim(X)[2])) CoV[i, ] = as.vector(Co[,
                                                           i, ])
      Ct = matrix(0, dim(X)[1], m * n)
      for (i in 1:n) Ct[, ((i - 1) * m + 1):((i - 1) *
                                               m + m)] = cbind(matrix(1, dim(X)[1], 1), X[,
                                                                                          , i]) %*% CoV[, ((i - 1) * m + 1):((i - 1) *
                                                                                                                               m + m)]
      mu = NA
    }
    else {
      X = NA
    }
  }
  Y = Uo + Ct
  for (t in (P + 1):T) {
    for (L in 1:Pmax) Y[t, ] = Y[t, ] + Y[t - L, ] %*% t(G[,
                                                           , L])
  }
  check = max(abs(Y))

  if (Ncommfakt>0) {

    uz = rnormSIGMA(T,diag(Ncommfakt))
    Z = uz
    for (i in 3:T) Z[i,] = Z[i-1,]*1.4-Z[i-2,]*0.4+uz[i,]

    if (anyNA(A)) {
      A = diag(m*n+Ncommfakt)
      #for (i in 1:n) {A[(i-1)*m+Ncommfakt+(m-(Ncommfakt-1)):m ] = diag(Ncommfakt)*50*max(abs(Y[,i*m]))/max(abs(Z))}
      A[(Ncommfakt+1):(Ncommfakt+n*m),1:Ncommfakt]  = 0.3
    }

    YY <- cbind(Z,Y);

    Y = (YY%*%t(A))[,Ncommfakt+(1:(m*n))]
    dim_GC = dim(G); dim_GC[1:2] = dim_GC[1:2]+Ncommfakt
    GC = c(1:prod(dim_GC))*0; dim(GC) = dim_GC
    GC[(Ncommfakt+1):(Ncommfakt+n*m),(Ncommfakt+1):(Ncommfakt+n*m),] = G
    GC[1:Ncommfakt,1:Ncommfakt,1] = 1.4*diag(Ncommfakt); GC[1:Ncommfakt,1:Ncommfakt,2] = -0.4*diag(Ncommfakt)  # assuming unit root common factor(s)
    for (i in 1:Pmax) GC[,,i] = A%*%GC[,,i]%*%solve(A)

    Uo = (cbind(uz,Uo)%*%t(A))[,Ncommfakt+(1:(m*n))]
    Do = (1:(m*n*Pmax*Ncommfakt))*0; dim(Do) = c(m,Ncommfakt,Pmax,n)

  }  else  {
    Z  = NA
    uz = NA
    A  = NA
    GC = NA
    Do = NA
  }

  result = list(Y, X, Uo, G, Sigmao, r_npo, Ao, Bo, Co, W, m, n, p, mu, check, type, crk, Ct,DFYflag,GC,uz,Z,A,Ncommfakt,Do)
  names(result) = c("Y", "X", "Uo", "G", "Sigmao", "r_npo", "Ao", "Bo", "Co", "W", "m", "n", "p", "mu", "check","type", "crk","Ct","DFYflag","GC","uz","Z","A","Ncommfakt","Do")
  return(result)
}


#' Estimation of CIGVAR(res)
#'
#' This function estimates the parameters of a CIGVAR(m,n,p) object based on provided data.
#' It runs a VECM estimation country/unit by country/unit under a given cointegration rank and pieces the results together to obtain a CIGVAR.
#' @param  res  a CIGVAR object of an output of CIGVARData including at least values of m,n,p,type,Y,crk and optionally X.
#' @return a CIGVAR object with estimated parameter values, AIC, BIC and LH
#' @examples
#'
#' n = 5
#' p = (1:15)*0; dim(p) = c(5,3)
#' p[,1] = 3; p[,2]=3;
#' r_npo = (1:(3 * 3 * 5))*0; dim(r_npo) = c(3,3,5)
#' r_npo[,,1] = matrix(c(1,1,3,2,3,3,3,3,3),3,3)
#' r_npo[,,2] = matrix(c(1,2,3,2,3,3,3,3,3),3,3)
#' r_npo[,,3] = matrix(c(1,2,3,2,3,3,3,3,3),3,3)
#' r_npo[,,4] = matrix(c(1,2,3,2,3,3,3,3,3),3,3)
#' r_npo[,,5] = matrix(c(1,2,3,2,3,3,3,3,3),3,3)
#'
#' res_d = CIGVARData(m=3,n=5,p=p,T=500,r_npo=r_npo,type="const",DFYflag=0,Ncommfakt=1)
#' max(res_d$Y)
#' STAT(res_d$G)
#' plot(ts(res_d$Y[,1:10]))
#' res_e = CIGVARest(res_d)
#' res_e$Summary
#' res_e$Summary$CRK_Test
#' ## adjust the model specification according to the cointegration rank test
#' res_d$crk = c(1,2,2,2,2)
#' res_e = CIGVARest(res_d)
#' res_e$Summary
#' ## unit roots in the CIGVAR model
#' STAT(res_e$G)
#' @export
CIGVARest <- function (res)  {
  m = res$m
  n = res$n
  p = res$p
  Ncommfakt = res$Ncommfakt
  Z = res$Z
  Y = res$Y; if (max(abs(Y)) > 10000000) return("Check data!!!");
  X = res$X
  W = res$W
  crk = res$crk
  type = res$type
  Bo = res$Bo
  Ao = res$Ao
  Co = res$Co
  Do = res$Do
  r_npo = res$r_npo
  Ncommtrend = res$Ncommtrend
  DFYflag = res$DFYflag
  Pmax = max(p[, 1:2])
  k = max(p[, 3])
  AIC = c(1:n) * 0
  BIC = c(1:n) * 0
  LH = c(1:n) * 0
  resid = Y * 0
  Tresid = resid[, 1:m]
  VECM_domestic  = list()
  VECM_domesticS = list()
  CRKtst = list()
  T = dim(Y)[1]
  FY = Y %*% t(W %x% diag(m))
  Bo = Bo = (1:(m * m * Pmax * n)) * 0
  dim(Bo) = c(m, m, Pmax, n)
  Bo = Bo * 0
  Ao = Bo * 0
  Co = Co * 0
  Do = Do * 0
  if (type == "none")
    C = matrix(0, m * n, 1) * 0
  if (type == "const")
    C = matrix(0, m * n, 1) * 0
  if (type == "exog0")
    C = matrix(0, m * n, dim(X)[2] + 1) * 0
  if (type == "exog1")
    C = matrix(0, m * n, dim(X)[2] + 1) * 0
  for (i in 1:n) {
    y = Y[, ((i - 1) * m + 1):((i - 1) * m + m)]
    x = FY[, ((i - 1) * m + 1):((i - 1) * m + m)]
    ### add colnames
    HilfeYi  = Y[, ((i - 1) * m + 1):((i - 1) * m + m)]
    HilfeFYi = FY[, ((i - 1) * m + 1):((i - 1) * m + m)]
    if (is.null(colnames(HilfeYi))) colnames(HilfeYi)   = paste0(rep("Y",ncol(HilfeYi)),c(1:ncol(HilfeYi)))
    if (is.null(colnames(HilfeFYi))) colnames(HilfeFYi) = paste0(rep("FY",ncol(HilfeYi)),c(1:ncol(HilfeFYi)))
    y = HilfeYi
    x = HilfeFYi

    if (type == "none")
      Model = "I"
    if (type == "const")
      Model = "III"

    if (type == "exog0") {
      Model = "I"
    }
    if (type == "exog1") {
      Model = "III"
    }
    P = matrix(0, 2, 3)
    P[1, 1:2] = p[i, 1:2]
    P[2, 1:2] = p[i, 1:2]

    if (!anyNA(X)) {
      Xi = as.matrix(X[,,i])
      if (is.null(colnames(Xi))) colnames(Xi)   = paste0(rep("exog",ncol(Xi)),c(1:ncol(Xi)))
    }   else  { Xi = NA }

    #tst <- MRCVECMest2(y, x, model = Model, type = "eigen", P = P, crk = crk[i], q = 0.95, Dxflag = 1,X=Xi)
    tst <- MRCVECMest2(y, x, model = Model, type = "eigen", P = P, crk = crk[i], q = 0.95, Dxflag = DFYflag,X=Xi,CZ=as.matrix(Z))
    #tst <- MRCVECMest2a(y,x, model = Model, type = "eigen",    P = P, crk = crk[i], q = 0.95, Dxflag = DFYflag,X=Xi,CZ=as.matrix(Z))


    #BAC = VECM2VAR(param = tst[[2]][[1]], beta = tst$beta, q = c(crk[i], P[1, 1] - 1, P[1, 2] - 1),kz=0,Dxflag=DFYflag)
    BAC = VECM2VAR(param = tst[[2]][[1]], beta = tst$beta, q = c(crk[i], P[1, 1] - 1, P[1, 2] - 1),kz=Ncommfakt,Dxflag=DFYflag)

    Tresid[(T - dim(tst[[2]]$residuals)[1] + 1):T, ] = tst[[2]]$residuals
    Sigma_one = t(Tresid) %*% Tresid/(T - dim(tst[[2]][[1]])[1])
    LH_P = -(T * m/2) * log(2 * pi) - (T * m/2) + (T/2) *
      log(det(solve(Sigma_one)))
    LH_AIC = 2 * (n * (dim(tst[[2]][[1]])[1]) + n * (n +
                                                       1)/2) - 2 * LH_P
    LH_BIC = log(T) * (n * (dim(tst[[2]][[1]])[1]) + n *
                         (n + 1)/2) - 2 * LH_P
    LH_N = 2 * n * (dim(tst[[2]][[1]])[1]) + n * (n + 1)
    AIC[i] = LH_AIC
    BIC[i] = LH_BIC
    LH[i] = LH_P
    VECM_domestic[[i]] = tst[[2]]
    VECM_domesticS[[i]] <- summaryCIVAR(lm1=tst[[2]],sname ="Z2")


    CRKtst[[i]] = tst[[1]]
    Bo[, , 1:p[i, 1], i] = BAC[[1]]
    Ao[, , 1:p[i, 2], i] = BAC[[2]]
    if (!anyNA(Do)) Do[, , 1:p[i, 1], i] = BAC[[4]]
    if (type=="none") Co[,,i] = Co[,,i]
    if (type=="const"|type=="exog1") Co[,,i] = BAC[[3]]
    if (type=="exog0") Co[,-1,i] = BAC[[3]]
    resid[1:nrow(tst[[2]]$residuals), (m * (i - 1) + 1):(i *
                                                           m)] = tst[[2]]$residuals
  }
  Sigmao = t(resid) %*% resid/(T - m * (p[i, 1] + p[i, 2]) -  p[i, 3])

  Sigmao = as.matrix(Matrix::nearPD(Sigmao)[[1]])

  G = BoAoW2G(Bo, Ao, W, m, n, Pmax)
  Gs = diag(n * m)
  for (L in 1:Pmax) {
    Gs = Gs - G[, , L]
  }
  LH_g = -T * m * n/2 * log(2 * pi) - T * m * n/2 + T/2 * log(det(solve(Sigmao)))
  TwoN = 2 * (sum(p[, 1:2]) * m * m + sum(p[, 3]) * m + (type ==
                                                           "const") * n * m)
  AIC_g = TwoN - 2 * LH_g
  BIC_g = log(T) * TwoN/2 - 2 * LH_g
  res$G = G
  res$C = C
  res$Sigmao = Sigmao
  res$r_npo = NA
  res$Ao = Ao
  res$Bo = Bo
  res$Co = Co
  Summary = list(VECM_domesticS,CRKtst,LH,AIC,BIC,LH_g,AIC_g,BIC_g,r_npo,Ncommfakt)
  names(Summary) = c("Estimation_Result","CRK_Test","Country_LH_function_Value","Country_AIC","Country_BIC","LH_g","AIC_g","BIC_g","r_npo","Ncommfakt")
  res$VECM_domestic = VECM_domestic
  res$Summary = Summary
  return(res)
}





#' Common trends for CIGVAR models
#'
#' This function creates independent common stochastic trends and idiosyncratic trends.
#'
#' @param m number of variables
#' @param p lags
#' @param T number of observations
#' @param r_npo roots of the characteristic polynomial
#' @param Ncommtrend number of common trends
#' @param n number of countries
#'
#' @return A list containing the parameter matrix of the domestic variables, the cointegration vectors and the adjustment vectors
#' @export
#' @keywords internal
#'
VARB_commtrend <- function (m, p, T, r_npo, Ncommtrend, n)
{
  alpha = list()
  beta = list()
  if (missing(Ncommtrend)) {
    Ncommtrend = NA
  }
  if (Ncommtrend==0) Ncommtrend <- NA
  Pmax = max(p[, 1:2])
  Pmin = min(p[, 1:2][p[, 1:2] > 0])
  Bo = (1:(m * m * Pmax * n)) * 0
  dim(Bo) = c(m, m, Pmax, n)
  if (anyNA(Ncommtrend)) {
    for (i in 1:n) {
      r_np = c(1:(m * p[i, 1])) * 0
      dim(r_np) = c(m, p[i, 1])
      r_np = r_npo[, 1:p[i, 1], i]
      VARD = VARData(m, p[i, 1], T, r_np = r_np)
      Bo[, , 1:p[i, 1], i] = VARD$B
      alpha[[i]] = B2CIB(VARD$B)[[2]]
      beta[[i]] = B2CIB(VARD$B)[[3]]
    }
  }
  if (!anyNA(Ncommtrend)) {
    B = matrix(stats::rnorm(m * m), m, m)
    Nnocommtrend = m - Ncommtrend
    B[(Nnocommtrend + 1):m, ] = 0
    B = B + diag(m)
    r_np = c(1:(Ncommtrend * Pmin))/c(1:(Ncommtrend * Pmin)) *
      1.5
    dim(r_np) = c(Ncommtrend, Pmin)
    r_np[1:Ncommtrend, 1] = 1
    VAR_IONE = VARData(n = Ncommtrend, p = Pmin, T = 100,
                       r_np = r_np)
    for (i in 1:n) {
      B = matrix(stats::rnorm(m * m), m, m)
      B[(Nnocommtrend + 1):m, ] = 0
      B = B + diag(m)
      Nnocommtrend = m - Ncommtrend
      r_np = c(1:(Nnocommtrend * p[i, 1])) * 0
      dim(r_np) = c(Nnocommtrend, p[i, 1])
      r_np <- r_npo[(Ncommtrend + 1):m, 1:p[i, 1], i]
      dim(r_np) = c(Nnocommtrend, p[i, 1])
      VAR_IZero = VARData(Nnocommtrend, p[i, 1], T = 100,
                          r_np = r_np)
      for (j in 1:p[i, 1]) {
        Dblock = matrix(0, m, m)
        Dblock[1:Nnocommtrend, 1:Nnocommtrend] = VAR_IZero$B[,
                                                             , j]
        if (dim(VAR_IONE$B)[3] >= j)
          Dblock[(Nnocommtrend + 1):m, (Nnocommtrend +
                                          1):m] = VAR_IONE$B[, , j]
        Bo[, , j, i] = B %*% Dblock %*% solve(B)
      }
      alpha[[i]] = B2CIB(Bo[, , 1:p[i, 1], i])[[2]]
      beta[[i]] = B2CIB(Bo[, , 1:p[i, 1], i])[[3]]
    }
  }
  return(list(Bo, alpha, beta))
}































#' Calculation of information criteria for CIGVAR models
#'
#' This function calculates AIC and BIC criteria for a range of lag specifications of the domestic and foreign variables up to the maximum given in L_V.
#'
#' @param  res  : an object generated from CIGVARest.
#' @param  L_V  : a 2 components vector containing the maxima of the domestic lag and the foreign lag, respectively.
#' @return a matrix with different lag specifications and the corresponding information criteria.
#' @examples
#'
#'
#' n = 5
#' p = (1:15)*0; dim(p) = c(5,3)
#' p[,1] = 3; p[,2]=3;
#' p[2,1] = 2;  p[3,1] = 2
#' res_d = CIGVARData(m=2,n=5,p=p,T=5000,type="const",DFYflag=0)
#' res_d$r_npo
#' res_d$check
#'
#' plot(ts(res_d$Y[1:100,1:10]))
#'
#' res_e = CIGVARest(res_d)
#' res_e$Summary
#' res_e$crk = c(1,1,1,1,1)
#' L_V = c(4,4)
#' CIGVARSelecte = CIGVAR_Select(res=res_e,L_V=c(4,4))
#' CIGVARSelecte[which.min(CIGVARSelecte[,17]),]
#' CIGVARSelecte[which.min(CIGVARSelecte[,19]),]
#'
#' @export
CIGVAR_Select = function(res,L_V=L_V)  {
Criteria = matrix(0,n*(L_V[1]-1)*(L_V[2]-1),4+n*3)
idx = 0
res_dd   = res
idx = 0
n   = res$n
for (I in 1:n )           {
for (l_d in 2: L_V[1] )   {
for (l_f in 2: L_V[2] )   {

      idx = idx + 1
      res_dd$p[I,1] = l_d
      res_dd$p[I,2] = l_f


    	res_s = CIGVARest(res_dd)
      Criteria[idx,] = c(as.vector(res_dd$p),res_s$Summary$AIC_g,res_s$Summary$BIC_g,sum(res_s$Summary$AIC),sum(res_s$Summary$BIC))
      #colnames(Criteria) = c("l_d","l_f","AIC_g","BIC_g","AIC","BIC")
    }
 }
}
return(Criteria)
}


#' Impulse Response Functions of CIGVAR(m,n,p)
#'
#' This function generates impulse response functions of an estimated CIGVAR(m,n,p) with confidence bands.
#'
#' @param  res  : a CIGVAR object of an output of CIGVARest.
#' @param  nstep : length of the impulse response functions
#' @param  comb : an mn-vector specifying combined impulse such as global shocks, regional shocks, or concerted actions.
#' @param  irf  : types of the impulse response functions. irf=c("gen","chol","chol1","gen1","comb1"), gen for generalized IRF with one standard deviation shocks, gen1 for generalized IRF with one unit impulse, chol for IRF with Cholezky decomposition of the covariance matrix, chol1 for Cholezky decomposition with one unit impulse, comb1 for concerted action with one unit impulse.
#' @param  runs : number of bootstrap runs to generate the confidence bands
#' @param  conf : a two component vector containing the tail probabilities of the bootstrap confidence interval
#' @return a matrix of (mn,mn,nstep,3) as the IRF columns representing the impulse rows the responses
#' @examples
#'
#' n = 5
#' p = (1:15)*0; dim(p) = c(5,3)
#' p[,1] = 3; p[,2]=3;
#' p[2,1] = 2;  p[3,1] = 2
#' res_d = CIGVARData(m=2,n=5,p=p,T=500,type="const",DFYflag=0)
#' res_d$r_npo
#' res_d$check
#' res_e = CIGVARest(res_d)
#' res_e$Summary
#'
#' ### Impulse response function
#'
#' IRF_CB = irf_CIGVAR_CB(res=res_e,nstep=20,comb=NA,irf="gen1",runs=200,conf=c(0.05,0.95))
#'
#' dim(IRF_CB)
#' IRF_g = IRF_graph(IRF_CB,Names=NA,response=c(1,4),impulse=c(1,2,3,4), ncol=4)
#'
#'
#'
#'
#' @export
irf_CIGVAR_CB <- function (res, nstep, comb, irf = c("gen", "chol", "chol1", "gen1", "comb1"), runs = 200, conf = c(0.05, 0.95))
{
  m = res$m
  n = res$n
  p = res$p
  T = dim(res$Y)[1]
  W = res$W
  r_npo = res$r_npo
  Ao = res$Ao
  Bo = res$Bo
  Go = res$G
  Co = res$Co
  Sigmao = res$Sigmao
  type = res$type
  X = res$X
  mu = res$mu
  DFYflag = res$DFYflag
  crk = res$crk
  Ncommfakt = res$Ncommfakt
  A = res$A
  uz = res$uz

  B = res$G
  neq = dim(B)[1]
  nvar = dim(B)[2]
  sigma = res$Sigmao
  response <- array(0, dim = c(neq, nvar, nstep, length(conf) + 1))
  response[, , , 1] <- irf_GVAR(res, nstep, comb, irf)
  responseR <- array(0, dim = c(neq, nvar, nstep, runs))
  for (i in 1:runs) {
    Uo_run = rnormSIGMA(T, sigma)
    res_run = CIGVARData(m, n, p, T, W, r_npo = NA, Ao, Bo,Co, Uo = Uo_run, Sigmao = NA, type, X, mu,DFYflag,crk,Ncommfakt,A,uz)
    res_e = CIGVARest(res_run)
    B_run = res_e$G
    sigma_run = res_e$Sigmao
    responseR[, , , i] <- irf_GVAR(res_e, nstep, comb, irf)
  }
  for (tt in 1:(nstep)) {
    for (i in 1:neq) {
      for (j in 1:nvar) {
        Ao
        response[i, j, tt, -1] = stats::quantile(responseR[i,j, tt, ], conf)
      }
    }
  }
  return(response)
}






