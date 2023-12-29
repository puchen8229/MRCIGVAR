#' Data generating process of CIVAR(p)
#'
#' This function generates data from a cointegrated vector autoregressive process CIVAR(p) and returns CIVAR(p) object that is a list containing data and parameters used in the CIVAR(p) process.
#'
#' @param n     : number of variables
#' @param p     : lag length
#' @param T     : number of observations
#' @param r_np  : an (n x p) matrix of roots of the characteristic polynomials of n independent AR(p)-processes. An element 1 in a row implies a unit root process. If the matrix is not provided, it will be generated randomly with one unit root in the first row.
#' @param A     : an (n x n) full rank transformation matrix that is used to generate correlated a CIVAR(p) from the n independent AR(p) processes with unit roots.
#' @param B     : an (n,n,p) array containing the coefficients of the CIVAR(p) process. If B is not given, it will be calculated out of r_np and A.
#' @param Co    : an (n,k+1) matrix containing coefficients of deterministic components in a CIVAR(p) process. For type="none" Co = 0*(1:n), for type="const" Co is an n vector, for type="exog0" Co is a (n,k) matrix, and for type="exog1" Co is an (n,1+k) matrix.
#' @param C1    : an (n,1) matrix containing coefficients of the trend component. C1 = 0 if there is no trend in the data.
#' @param U     : residuals, if it is not NA it will be used as input to generate the CIVAR(p)
#' @param Sigma : an n x n covariance matrix of the CIVAR(p)
#' @param type  : types of deterministic components. "none","rconst","const","rtrend", "trend", "exog0" and "exog1" are 7 options
#'                \itemize{
#'                \item{none: }{No deterministic component in VECM (Case I)}
#'                \item{rconst: }{An intercept restricted in the cointegration space (Case II)}
#'                \item{const: }{An unrestricted intercept (Case III)}
#'                \item{rtrend: }{An trend component restricted in the cointegration space and an unrestricted intercept in VECM (Case IV)}
#'                \item{trend: }{An unrestricted trend component and an unretricted intercept in VECM (Case V)}
#'                \item{exog0: }{Model with exogeneous variable and no deterministic components}
#'                \item{exog1: }{Model with exogeneous variable and an unrestricted intercept}
#'                }
#' @param X     : a (T x k) matrix of exogenous variables.
#' @param mu    : an n vector of drifts of the CIVAR(p) process
#' @param Yo    : (p x n) matrix of initial values of the process
#' @param crk   : number of cointegration relations. It equals n-r, where r is the number of unit roots.
#' @return      An object of CIVAR(p) containing the generated data, the parameters used and the exogenous variables. res = list(n,p,type,r_np,Phi,A,B,Co,Sigma,Y,X,resid,U,Y1,Yo,check)
#' @examples
#' res_d = CIVARData(n=2,p=2,T=100,type="const")
#' res_d = CIVARData(n=2,p=2,T=10,Co=c(1:2)*0,type="none")
#' res_d = CIVARData(n=2,p=2,T=10,Co=c(1:2)*NA,  type="const")
#'
#' p = 3
#' n = 4
#' r_np = matrix(c(1,2,1.5,1.5,2.5,2.5,2,-1.5,2,-4,1.9,-2.1),4,3)
#'
#' res_d  = CIVARData(n=4,p=3,T=200,r_np=r_np,Co=matrix(0,n,1) ,type="none",crk=3)
#' res_e  = CIVARest(res=res_d)
#' sum(abs(res_e$B-res_d$B))
#' sum(abs(res_e$Co-res_d$Co))
#'
#'
#' res_d  = CIVARData(n=4,p=3,T=200)
#' res_e = CIVARest(res=res_d)
#' sum(abs(res_e$B-res_d$B))
#' sum(abs(res_e$Co-res_d$Co))
#' plot(ts(res_d$Y))
#' @export
CIVARData <- function (n, p, T, r_np, A, B, Co, C1,U, Sigma, type, X, mu, Yo, crk)
{
  if (missing(r_np)) {
    r_np = NA
  }
  if (missing(A)) {
    A = NA
  }
  if (missing(B)) {
    B = NA
  }
  if (missing(Co)) {
    Co = NA
  }
  if (missing(C1)) {
    C1 = NA
  }
  if (missing(U)) {
    U = NA
  }
  if (missing(Sigma)) {
    Sigma = NA
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
  if (missing(Yo)) {
    Yo = NA
  }
  if (missing(crk)) {
    crk = NA
  }
  if (anyNA(r_np)) {
    r_np = matrix(0, n, p)
    for (i in 1:n) r_np[i, ] <- 0.51/(stats::runif(p) - 0.5)
    if (anyNA(crk)) {
      r_np[1, 1] = 1
      crk = n - 1
    }
    else {
      for (j in 1:(n - crk)) r_np[j, 1] = 1
    }
  }
  d = dim(r_np)
  if (!(n == d[1] & p == d[2])) {
    print("dimension problem")
    return("dimension")
  }
  Phi = NA
  if (anyNA(Phi)) {
    Phi = matrix(0, n, p)
    for (i in 1:n) Phi[i, ] = Roots2coef(p, r_np[i, ])
  }
  if (anyNA(A)) {
    A = matrix(stats::runif(n * n), n, n)
    A = A %*% t(A)
    A = eigen(A)[[2]]
  }
  if (anyNA(B)) {
    B = matrix(0, n, n * p)
    dim(B) = c(n, n, p)
    if (n == 1) {
      B = Phi
      dim(B) = c(n, n, p)
    }
    if (n > 1) {
      for (i in 1:p) B[, , i] = A %*% diag(Phi[, i]) %*%
          t(A)
    }
  }
  if (anyNA(Sigma)) {
    Sigmao = matrix(stats::rnorm(n * n), n, n)
    Sigmao = Sigmao %*% t(Sigmao)
    Sigma = A %*% Sigmao %*% t(A)
  }
  else {
    Sigmao = solve(A) %*% Sigma %*% solve(t(A))
  }
  if (anyNA(U)) {
    Uh = rnormSIGMA(T, Sigmao)
    U = Uh %*% t(A)
  }
  else {
    Uh = U %*% solve(t(A))
  }
  if (anyNA(Yo)) {
    Yo = Uh
  }
  else {
    Uh[1:p, ] = Yo
    Yo = Uh
  }
  for (i in 1:n) {
    for (t in (p + 1):T) {
      for (L in 1:p) Yo[t, i] = Yo[t, i] + Yo[t - L, i] *
          Phi[i, L]
    }
  }
  Y1 = Yo %*% t(A)
  Ct = U * 0
  if (anyNA(type)) {
    type = "const"
  }
  #type_soll = Type(Co, X, type)
  #if (!type_soll == type)
  #  return(list("type mismatch", type_soll, type))
  if (!anyNA(Co) & (!nrow(as.matrix(Co)) == n))
    return("dimension mismatch")

  if ((type == "none") & (anyNA(Co)) & (anyNA(C1))) {
    Co = matrix(0, n, 1)
    C1 = matrix(0, n, 1)
  }
  if ((type == "rtrend") & (anyNA(Co)) & (anyNA(C1))) {
    if (anyNA(mu)) {
      mu = matrix(stats::rnorm(n), n, 1)
    }
    if (anyNA(Co)) {
      Co = mu
      for (L in 1:p) Co = Co - B[, , L] %*% mu
    }
    else {
      H = diag(n)
      for (L in 1:p) H = H - B[, , L]
      if (min(abs(eigen(H)$values)) > 1e-13)
        mu = solve(H) %*% Co
      else mu = NA
    }

    gamma0 <- rep(1,n)
    CIBo <- B2CIB(B)[[1]][,,1]
    C1 <- CIBo%*%gamma0
    Ct <- matrix(c(1:T),T,1)%*%t(C1) + matrix(1, T, 1) %*% t(Co)
  }

  if ((type == "trend") & (anyNA(C1)) & (anyNA(Co))) {
    if (anyNA(Co))  Co = matrix(2, n, 1)
    if (anyNA(C1))  C1 = matrix(1, n, 1)
    CIBo <- B2CIB(B)[[1]][,,1]
    Ct <- matrix(c(1:T),T,1)%*%t(C1) + matrix(1, T, 1) %*% t(Co)
  }
  if (type == "const") {
    if (anyNA(mu)) {
      mu = matrix(stats::rnorm(n), n, 1)
    }
    if (anyNA(Co)) {
      Co = mu
      for (L in 1:p) Co = Co - B[, , L] %*% mu
    }
    else {
      H = diag(n)
      for (L in 1:p) H = H - B[, , L]
      if (min(abs(eigen(H)$values)) > 1e-13)
        mu = solve(H) %*% Co
      else mu = NA
    }
    Ct = matrix(1, T, 1) %*% t(Co)
    C1 = matrix(0,n,1)
  }
  if (type == "rconst") {
    gamma0 = rep(1,n)
    CIBo <- B2CIB(B)[[1]][,,1]
    Co = CIBo%*%gamma0
    C1 = matrix(0,n,1)
    Ct = matrix(1, T, 1) %*% t(Co)
  }


  if (type == "exog0" & anyNA(Co)) {
    k = ncol(as.matrix(X))
    Co = matrix(stats::rnorm(n * (k + 1)), n, k + 1)
    Co[, 1] = 0
  }
  if (type == "exog1" & anyNA(Co)) {
    k = ncol(as.matrix(X))
    Co = matrix(stats::rnorm(n * (k + 1)), n, k + 1)
  }
  if (!anyNA(X)) {
    if (type == "exog0")
      Ct = X %*% t(Co[, -1])
    if (type == "exog1")
      Ct = as.matrix(cbind((1:T)/(1:T), X)) %*% t(Co)
  }
  else {
    X = NA
  }
  if (n > 0) {
    Y = U + Ct
    for (tt in (p + 1):T) {
      for (L in 1:p) Y[tt, ] = Y[tt, ] + Y[tt - L, ] %*%
          t(B[, , L])
    }
  }
  check = max(abs(Y1))
  resid = Uh
  result = list(n, p, type, r_np, Phi, A, B, Co, C1, Sigma, Y,X, resid, U, Y1, Yo, check, crk)
  names(result) = c("n", "p", "type", "r_np", "Phi", "A", "B", "Co","C1", "Sigma","Y", "X", "resid", "U", "Y1","Yo", "check", "crk")
  return(result)
}



#' Estimation of CIVAR(p)
#'
#' This function estimates parameters of a specified CIVAR(p) model based on provided data.
#'
#' @param  res  :an object of CIVAR(p) containing the components which are the output of CIVARData including at least: n, p, Y, crk, and optionally X and type.
#' @return res  an object of CIVAR(p) containing estimated parameter values, AIC, BIC, LH and the estimated VECM in regression format.
#' @examples
#' p = 3
#' n = 4
#' r_np = matrix(c(1,2,1.5,1.5,2.5,2.5,2,-1.5,2,-4,1.9,-2.1),4,3)
#'
#' res_d  = CIVARData(n=4,p=3,T=200,r_np=r_np,Co =(1:n)/(1:n)*0,type="none",crk=3)
#' res_e  = CIVARest(res=res_d)
#' res_d  = CIVARData(n=4,p=3,T=200)
#' B = res_d$B
#' plot(ts(res_d$Y))
#' res_d$Co
#' res_d$type
#' res_d$crk
#' res_e = CIVARest(res=res_d)
#'
#'
#'
#' @export
CIVARest <- function (res)
{
  n = res$n
  p = res$p
  crk = res$crk
  Y = as.matrix(res$Y)
  if (is.null(colnames(Y)))
    colnames(Y) = paste0(rep("Y", ncol(Y)), c(1:ncol(Y)))
  type = res$type
  T = dim(Y)[1]
  X = res$X
  if (!anyNA(X))
    if (is.null(colnames(X)))
      colnames(X) = paste0(rep("exog", ncol(X)),c(1:ncol(X)))
  if (type == "none" | type == "exog0")
    Model = "I"
  if (type == "const" | type == "exog1")
    Model = "III"
  if (type == "rconst")
    Model = "II"
  if (type == "rtrend")
    Model = "IV"
  if (type == "trend")
    Model = "V"
  P = matrix(0, 2, 3)
  P[, 1] = p
  Co = res$Co * 0
  C1 = res$C1 * 0
  tst <- MRCVECMest2(y = Y, x = 0, model = Model, type = "eigen",P = P, crk = crk, q = 0.95, Dxflag = 0, X = X)
  param = tst[[2]][[1]]
  CIVAREST = VECM2VAR(param = param, beta = tst$beta, q = c(crk, p - 1, 0))
  B = VECM2VAR(param = param, beta = tst$beta, q = c(crk, p - 1, 0))[[1]]
  C = VECM2VAR(param = param, beta = tst$beta, q = c(crk, p - 1, 0))[[3]]
  alpha <- as.matrix(t(param[1:crk, , drop = FALSE]))
  beta <- as.matrix(tst$beta)
  if (type == "rconst")
    Co = alpha %*% t(beta)[, n + 1]
  if (type == "rtrend")
    C1 = alpha %*% t(beta)[, n + 1]
  if (type == "trend") {
    C1 = param[crk + (p - 1) * n + 2, ]
    Co = param[crk + (p - 1) * n + 1, ]
  }
  LREG = tst[[2]]
  sigma = t(LREG$residuals) %*% (LREG$residuals)/(T)
  resid = Y * 0
  resid[(p + 1):T, ] = LREG$residuals
  if (type == "const" | type == "exog1")
    Co = C
  if (type == "exog0") {
    Co[, 2:(1 + dim(CIVAREST[[3]])[2])] = CIVAREST[[3]]
  }
  res$B <- B
  res$Co <- Co
  res$C1 <- C1
  res$Sigma <- sigma
  res$resid <- resid
  LH = -(T * n/2) * log(2 * pi) - (T * n/2) + ((T)/2) * log(det(solve(sigma)))
  AIC = 2 * n * (dim(tst[[2]][[1]])[1]) + n * (n + 1) - 2 *
    LH
  BIC = log(T) * (n * (dim(tst[[2]][[1]])[1]) + n * (n + 1)/2) -
    2 * LH
  LH_N = 2 * n * (dim(tst[[2]][[1]])[1]) + n * (n + 1)
  if (is.null(colnames(Y)))
    colnames(Y) = sprintf("Y%s", 1:n)
  estimation_result = summaryCIVAR(lm1 = LREG, sname = "Z2")
  Summary = list(estimation_result, LH, AIC, BIC, LH_N, tst$erg)
  names(Summary) = c("Estimation_Results", "LH_function_Value",
                     "AIC", "BIC", "LH_N", "Johansen_Test")
  res$tst = tst
  res$Summary = Summary
  return(res)
}






#' Estimation of a multi regime conditional vector error correction process.
#'
#' This function estimates the unknown parameters of a multi regime conditional VECM based on provided data.
#'
#' @param y	: data matrix of the endogenous variable
#' @param x	: x data matrix of the conditioning variables. If x is missing, it will estimate a unconditional VECM.
#' @param s     : the series of state variable of value (0,1) for the two regime case. Missing s implies single regime VECM.
#' @param model : It assumes one of the values in c("I","II","III","IV","V") corresponding to the five cases in Johansen test.
#' @param type  : c("eigen", "trace") corresponding to the trace test or the max eigenvalue test
#' @param crk   : cointegration rank
#' @param P     : 2 x 2 matrix containing the lag of the cointegrated VAR process. The first row is the lag of the first regime and the number of exogenous variables. If the second row is zero, this is a one regime VECM.
#' @param q     : The significance level
#' @param Dxflag : A flag indicating if the conditioning variables enter the cointegration space
#' @param X     : The conditioning variables
#' @param CZ    : Common factor variables in CIGVAR and MRCIGVAR
#' @return        A list containing: the result of the Johansen test, VECM in regression format, lambda, beta, PI, GAMMA, model, P, s
#'
#' @export
MRCVECMest2 <- function (y, x, s, model = c("I", "II", "III","IV", "V"), type = c("eigen", "trace"), crk = 2, P = matrix(2, 2, 2), q = 0.95, Dxflag = 0, X = NA,CZ=NA)
{
  y <- as.matrix(y)
  if (q != 0.9 && q != 0.95 && q != 0.99) {
    return("please correct significance level")
  }
  S = 2
  p <- as.integer(max(P))
  pmin <- as.integer(min(P))
  N1 <- ncol(as.matrix(y))
  NN1 <- ncol(as.matrix(x))
  N <- ncol(as.matrix(y))
  n = N
  if (N1 < crk) {
    return("y's dimension must be larger than crk")
  }
  if (missing(s)) {
    s = NA
  }
  if (missing(Dxflag)) {
    Dxflag = 0
  }
  if (!anyNA(s)) {
    St = s[(p + 1):length(s)]
    NSt = 1 - s[(p + 1):length(s)]
  }
  if (!anyNA(CZ)) {
    y = cbind(CZ,y);
    kz = ncol(CZ);
  } else {
    kz = 0
  }
  if (!sum(abs(x)) == 0) {
    z <- cbind(y, x)
    Zy <- Embed(diff(y), p, prefix = "d")
    Zx <- Embed(diff(x), p, prefix = "d")
    if (P[1, 1] < p) {
      Aa = P[1, 1] * N1 + (1:((p - P[1, 1]) * N1))
      ZyI <- Zy[, -Aa]
    }
    else ZyI = Zy
    if (P[1, 2] < p) {
      Ba = P[1, 2] * NN1 + (1:((p - P[1, 2]) * NN1))
      ZxI = Zx[, -Ba]
    }
    else ZxI = Zx
    if (P[2, 1] < p) {
      Aa = P[2, 1] * N1 + (1:((p - P[2, 1]) * N1))
      ZyII = Zy[, -Aa]
    }
    else ZyII = Zy
    if (P[2, 2] < p) {
      Bb = P[2, 2] * NN1 + (1:((p - P[2, 2]) * NN1))
      ZxII = Zx[, -Bb]
    }
    else ZxII = Zx
    Z = cbind(ZyI, ZxI)
    Z_2 = cbind(ZyII, ZxII)

    #Z2   <-   Z[, -c(1:N1, P[1, 1] * N1 + 1:NN1)]
    #ZS_2 <- Z_2[, -c(1:N1, P[2, 1] * N1 + 1:NN1)]
    Z2 <-      Z[, -c(1:(kz+N1), P[1, 1] * (kz+N1) + 1:NN1)]
    ZS_2 <-  Z_2[, -c(1:(kz+N1), P[2, 1] * (kz+N1) + 1:NN1)]

  }
  else {
    z = y
    Z = Embed(diff(y), p, prefix = "d")
    if (p > 1) {
      if (Dxflag == 0) {
        Z2 <- Z[, -c(1:N1)]
        ZS_2 <- Z2
        Z2 <- Z2[, 1:((P[1, 1] - 1) * n)]
        ZS_2 <- ZS_2[, 1:((P[2, 1] - 1) * n)]
      }
      else {
        Z2 <- Z[, -c(1:N1)]
        Zs_2 <- Z2
        Z2 <- Z2[, 1:((P[1, 1] - 1) * n)]
        ZS_2 <- ZS_2[, 1:((P[2, 1] - 1) * n)]
      }
    }
    if (p == 1) {
      Z2 <- Z[, -c(1:N1)]
      ZS_2 <- Z2
    }
  }
  M1 <- ncol(as.matrix(z))
  T <- nrow(as.matrix(z))
  MM1 = ncol(Z)
  Y0 <- Z[, kz+c(1:N1)]
  if (Dxflag==1)   Z1 <- z[-T, ][p:(T - 1), ]
  if (Dxflag==0)   Z1 <- z[-T,1:(kz+N1)][p:(T - 1), ]
  #Z1 <- z[-T, ][p:(T - 1), ]
  M1 <- ncol(as.matrix(Z1))
  if (!anyNA(s))
    Z2 = cbind(St * Z2, NSt * ZS_2)
  T1 <- nrow(as.matrix(Y0))
  TT <- nrow(as.matrix(y))
  lT = matrix(1, T1, 1)
  colnames(lT) <- "Const"
  Trend <- matrix((p + 1):TT, (TT - p), 1)
  colnames(Trend) <- "Trend"
  if (model == "I") {
    Y0 = Y0
    Z1 = Z1
    Z2 = Z2
    if (!anyNA(X) && !anyNA(s))
      Z2 = cbind(Z2, St * X[(nrow(X) - nrow(Z2) + 1):nrow(X),
                            , 1], NSt * X[(nrow(X) - nrow(Z2) + 1):nrow(X),
                                          , 2])
    if (!anyNA(X) && anyNA(s)) {
      Xi = as.matrix(X[(nrow(X) - nrow(Z2) + 1):nrow(X),
      ])
      colnames(Xi) <- colnames(X)
      Z2 = cbind(Z2, Xi)
    }
  }
  if (model == "II") {
    Y0 = Y0
    Z1 = cbind(Z1, lT)
    Z2 = Z2
  }
  if (model == "III") {
    Y0 = Y0
    Z1 = Z1
    if (!anyNA(X) && !anyNA(s))
      Z2 = cbind(Z2, St, St * X[(nrow(X) - nrow(Z2) + 1):nrow(X),
                                , 1], NSt, NSt * X[(nrow(X) - nrow(Z2) + 1):nrow(X),
                                                   , 2])
    if (anyNA(X) && !anyNA(s))
      Z2 = cbind(Z2, St, NSt)
    if (!anyNA(X) && anyNA(s)) {
      Xi = as.matrix(X[(nrow(X) - nrow(Z2) + 1):nrow(X),
      ])
      colnames(Xi) <- colnames(X)
      Z2 = cbind(Z2, lT, Xi)
    }
    if (anyNA(X) && anyNA(s))
      Z2 = cbind(Z2, lT)
  }
  if (model == "IV") {
    Y0 = Y0
    Z1 = cbind(Z1, Trend)
    Z2 = cbind(Z2, lT)
  }
  if (model == "V") {
    Y0 = Y0
    Z1 = cbind(Z1)
    Z2 = cbind(Z2, lT, Trend)
  }
  if (ncol(Z2)>0) {
    Z2 <- ShiftZ2(Z2,kz,N1,P[1,1]-1)
    M00 <- crossprod(Y0)/T1
    M11 <- crossprod(Z1)/T1
    M22 <- crossprod(Z2)/T1
    M01 <- crossprod(Y0, Z1)/T1
    M02 <- crossprod(Y0, Z2)/T1
    M10 <- crossprod(Z1, Y0)/T1
    M20 <- crossprod(Z2, Y0)/T1
    M12 <- crossprod(Z1, Z2)/T1
    M21 <- crossprod(Z2, Z1)/T1
    M22inv <- solve(M22)
    R0 <- Y0 - t(M02 %*% M22inv %*% t(Z2))
    R1 <- Z1 - t(M12 %*% M22inv %*% t(Z2))
  }
  if (ncol(Z2)==0) {
    R0 = Y0
    R1 = Z1
  }
  S00 <- crossprod(R0)/T1
  S01 <- crossprod(R0, R1)/T1
  S10 <- crossprod(R1, R0)/T1
  S11 <- crossprod(R1)/T1
  Ctemp <- chol(S11, pivot = TRUE)
  pivot <- attr(Ctemp, "pivot")
  oo <- order(pivot)
  C <- t(Ctemp[, oo])
  Cinv <- solve(C)
  S00inv <- solve(S00)
  valeigen <- eigen(Cinv %*% S10 %*% S00inv %*% S01 %*% t(Cinv))
  lambda <- Re(valeigen$values)
  e <- Re(valeigen$vector)
  V <- t(Cinv) %*% e
  Vorg <- V
  V <- sapply(1:M1, function(j) V[, j]/V[1, j])
  PI <- S01 %*% solve(S11)
  if (ncol(Z2)>0)    GAMMA <- M02 %*% M22inv - PI %*% M12 %*% M22inv
  if (ncol(Z2)==0)   GAMMA <- S00 - PI %*% S10

  beta <- (V[, 1:crk])
  if (crk > 0) {
    CI = Z1 %*% beta
    if (ncol(Z2)>0) {
      VECM <- stats::lm(Y0 ~ 0 + CI + Z2)
      VECMS <- summaryCIVAR(lm1 = VECM, sname = "Z2")
    }
    if (ncol(Z2)==0) {
      VECM <- stats::lm(Y0 ~ 0 + CI )
      VECMS <- summaryCIVAR(lm1 = VECM, sname = "")
    }

  }
  if (crk == 0) {
    CI = 0
    if (ncol(Z2)>0) {
      VECM <- stats::lm(Y0 ~ 0 + Z2)
    }
    if (ncol(Z2)==0) {
      VECM <- stats::lm(Y0 ~ 0)
    }
  }
  E = -T1 * log(1 - lambda)
  E = E[1:N]
  resultsvecm <- summary(VECM)
  Omega = t(VECM$residuals) %*% VECM$residuals/nrow(VECM$residuals)
  if (model == "I") {
    Tab <- Tab1
  }
  if (model == "II") {
    Tab <- Tab2
  }
  if (model == "III") {
    Tab <- Tab3
  }
  if (model == "IV") {
    Tab <- Tab4
  }
  if (model == "V") {
    Tab <- Tab5
  }
  b = c(1:12)
  for (i in 1:12) {
    b[i] = Tab[2 * (i - 1) + 1, 2 + M1 - N1]
  }
  a = c(1:12)
  for (i in 1:12) {
    a[i] = Tab[2 * i, 2 + M1 - N1]
  }
  if (type == "eigen") {
    critical_vals = b
    M = matrix(0, N, 1)
    j = 1
    rank = 0
    while (j <= N && E[j] > critical_vals[N + 1 - j]) {
      M[j, ] = M[j, ] + 1
      j = j + 1
      rank = rank + 1
    }
    erg <- cbind(E, critical_vals[N:1])
    colnames(erg) <- c("teststatistic", "critical_value")
    if (N > 1) {
      rownames(erg) <- c("crk <= 0 |", paste("crk <= ",
                                             1:(N - 1), " |", sep = ""))
    }
    if (N == 1) {
      rownames(erg) <- c("crk <= 0 |")
    }
    coint_rank <- paste("Johansen-Test (with maximum-eigenvalue-teststatistic) indicates",
                        rank, "cointegrating equation(s) at the", 1 -
                          q, "level")
    result <- new.env()
    result$erg = erg
    result$estimation = VECM
    result$VECMS = VECMS
    result$lambda = E
    result$z = z
    result$Z2 = Z2
    result$Z1 = Z1
    result$Y0 = Y0
    result$R1 = R1
    result$R0 = R0
    result$Omega = Omega
    result$beta = beta
    result$PI = PI
    result$GAMMA = GAMMA
    result$model = model
    result$P = P
    result$NN1 = NN1
    result$s = s
    rr <- as.list(result)
    return(rr)
  }
  else {
    type = "trace"
    critical_vals = a
    stat = matrix(0, N, 1)
    for (i in 1:N) {
      sum = 0
      for (j in i:N) {
        sum = sum + E[j]
      }
      stat[i] = sum
    }
    M = matrix(0, N, 1)
    j = 1
    rank = 0
    while (stat[j] > critical_vals[N + 1 - j] && j <= N) {
      M[j, ] = M[j, ] + 1
      j = j + 1
      rank = rank + 1
    }
    erg <- cbind(stat, critical_vals[N:1])
    colnames(erg) <- c("teststatistic", "critical_value")
    if (N > 1) {
      rownames(erg) <- c("crk <= 0 |", paste("crk <= ",
                                             1:(N - 1), " |", sep = ""))
    }
    if (N == 1) {
      rownames(erg) <- c("crk <= 0 |")
    }
    coint_rank <- paste("Johansen-Test (with trace-teststatistic) indicates",
                        rank, "cointegrating equation(s) at the", 1 -
                          q, "level")
    result <- new.env()
    result$erg = erg
    result$estimation = VECM
    result$VECMS = VECMS
    result$lambda = E
    result$z = z
    result$Z2 = Z2
    result$Z1 = Z1
    result$Y0 = Y0
    result$R1 = R1
    result$R0 = R0
    result$Omega = Omega
    result$beta = beta
    result$PI = PI
    result$GAMMA = GAMMA
    result$model = model
    result$P = P
    result$NN1 = NN1
    result$s = s
    rr <- as.list(result)
    return(rr)
  }
}


#' Estimation of a two regime conditional vector error correction model.
#'
#' This function estimates the unknown parameters of a two regime conditional VECM based on provided data. The cointegrating vectors are identical in the two regimes but the adjustment speeds are different in the two regimes.
#'
#' @param y	: data matrix of the endogenous variables
#' @param x	: data matrix of the conditioning variables. If x is missing, it will estimate an unconditional VECM.
#' @param s     : the series of state variable with values (0,1) for the two regimes. Missing s implies single regime VECM.
#' @param model : It assumes one of the values in c("I","II","III","IV","V") corresponding to the five cases in the Johansen test.
#' @param type  : c("eigen", "trace") corresponding to the trace test or the max eigenvalue test
#' @param crk   : cointegration rank
#' @param P     : a (2 x 2) matrix containing the lag of the cointegrated VAR process. The first row is the lag of the first regime and the number of exogenous variables. If the second row is zero, this is a one regime VECM.
#' @param q     : The significance level used
#' @param Dxflag : a flag indicating if the conditioning variables enter the cointegration space
#' @param X     : the exogenous stationary variables
#' @param CZ    : the exogenous non-stationary common factors
#' @return      : a list containing: the result of JH test, VECM in regression format, lambda, beta, PI, GAMMA, model, P, s
#' @export
MRCVECMestm <- function (y, x, s, model = c("I", "II", "III", "IV", "V"), type = c("eigen", "trace"), crk = 2, P = matrix(2, 2, 2), q = 0.95, Dxflag = 0, X = X, CZ=NA)
{
  y <- as.matrix(y)
  if (q != 0.9 && q != 0.95 && q != 0.99) {
    return("please correct significance level")

  }
  S = 2
  p <- as.integer(max(P))
  pmin <- as.integer(min(P))
  N1 <- ncol(as.matrix(y))
  NN1 <- ncol(as.matrix(x))
  N <- ncol(as.matrix(y))
  n = N
  if (N1 < crk) {
    return("y's dimension must be larger than crk")

  }
  if (missing(s)) {
    s = NA
  }
  if (missing(Dxflag)) {
    Dxflag = 0
  }
  if (!anyNA(s)) {
    St = s[(p + 1):length(s)]
    NSt = 1 - s[(p + 1):length(s)]
  }
  if (!anyNA(CZ)) {
    y = cbind(CZ,y);
    kz = ncol(CZ);
  } else {
    kz = 0
  }
  if (!sum(abs(x)) == 0) {
    z <- cbind(y, x)
    Zy <- Embed(diff(y), p, prefix = "d")
    Zx <- Embed(diff(x), p, prefix = "d")
    if (P[1, 1] < p) {
      Aa = P[1, 1] * N1 + (1:((p - P[1, 1]) * N1))
      Zy1 = Zy[, -Aa]
    }
    else Zy1 = Zy
    if (P[1, 2] < p) {
      Bb = P[1, 2] * NN1 + (1:((p - P[1, 2]) * NN1))
      Zx1 = Zx[, -Bb]
    }
    else Zx1 = Zx
    if (P[2, 1] < p) {
      Aa = P[2, 1] * N1 + (1:((p - P[2, 1]) * N1))
      Zy2 = Zy[, -Aa]
    }
    else Zy2 = Zy
    if (P[2, 2] < p) {
      Bb = P[2, 2] * NN1 + (1:((p - P[2, 2]) * NN1))
      Zx2 = Zx[, -Bb]
    }
    else Zx2 = Zx
    Z = cbind(Zy1, Zx1)
    Z_2 = cbind(Zy2, Zx2)

    Z2 <-     Z[, -c(1:(kz+N1), P[1, 1] * (kz+N1) + 1:NN1)]
    ZS_2 <- Z_2[, -c(1:(kz+N1), P[2, 1] * (kz+N1) + 1:NN1)]


  }
  else {
    z = y
    Z = Embed(diff(y), p, prefix = "d")
    if (Dxflag == 0) {
      Z2 <- Z[, -c(1:N1)]
      ZS_2 <- Z2
      Z2 <- Z2[, 1:((P[1, 1] - 1) * n)]
      ZS_2 <- ZS_2[, 1:((P[2, 1] - 1) * n)]
    }
    else {
      Z2 <- Z[, -c(1:N1)]
      Zs_2 <- Z2
      Z2 <- Z2[, 1:((P[1, 1] - 1) * n)]
      ZS_2 <- ZS_2[, 1:((P[2, 1] - 1) * n)]
    }
  }

  T <- nrow(as.matrix(z))
  MM1 = ncol(Z)
  Y0 <- Z[, kz+c(1:N1)]
  if (Dxflag==1)   Z1 <- z[-T, ][p:(T - 1), ]
  if (Dxflag==0)   Z1 <- z[-T,1:(kz+N1)][p:(T - 1), ]
  #Z1 <- z[-T, ][p:(T - 1), ]
  #M1 <- ncol(as.matrix(z))
  M1 <- ncol(as.matrix(Z1))
  if (!anyNA(s))
    Z2 = cbind(St * Z2, NSt * ZS_2)
  T1 <- nrow(as.matrix(Y0))
  lT = (1:T1)/(1:T1)
  Trend <- matrix(1:T1, T1, 1)
  if (model == "I") {
    Y0 = Y0
    Z1 = Z1
    Z2 = Z2
    if (!anyNA(X) & !anyNA(s)) {
      StX  = as.matrix(St * X[(nrow(X) - nrow(Z2) + 1):nrow(X), , 1]); colnames(StX)  = paste0(colnames(X),"1")
      NStX = as.matrix(NSt *X[(nrow(X) - nrow(Z2) + 1):nrow(X), , 2]); colnames(NStX) = paste0(colnames(X),"2")
      Z2 = cbind(Z2,StX,NStX)
      #Z2 = cbind(Z2, St * X[(nrow(X) - nrow(Z2) + 1):nrow(X), , 1], NSt * X[(nrow(X) - nrow(Z2) + 1):nrow(X),, 2])
    }
    if (!anyNA(X) & anyNA(s))
      Z2 = cbind(Z2, X[(nrow(X) - nrow(Z2) + 1):nrow(X),
      ])
  }
  if (model == "II") {
    Y0 = Y0
    Z1 = cbind(T1, Z1)
    Z2 = Z2
  }
  if (model == "III") {
    Y0 = Y0
    Z1 = Z1
    if (!anyNA(s) & !anyNA(X)) {
      StX  = as.matrix(St * X[(nrow(X) - nrow(Z2) + 1):nrow(X), , 1]); colnames(StX)  = paste0(colnames(X),"1")
      NStX = as.matrix(NSt *X[(nrow(X) - nrow(Z2) + 1):nrow(X), , 2]); colnames(NStX) = paste0(colnames(X),"2")


      Z2 = cbind(Z2,St,StX,NSt,NStX)
      #Z2 = cbind(Z2, St, St * X[(nrow(X) - nrow(Z2) + 1):nrow(X), , 1,drop=FALSE], NSt, NSt * X[(nrow(X) - nrow(Z2) + 1):nrow(X), , 2,drop=FALSE])

    }
    if (!anyNA(s) & anyNA(X))
      Z2 = cbind(Z2, St, NSt)
    if (anyNA(s) & !anyNA(X))
      Z2 = cbind(Z2, rep(1, T1), X[(nrow(X) - nrow(Z2) +
                                      1):nrow(X), ])
    if (anyNA(s) & anyNA(X))
      Z2 = cbind(Z2, rep(1, T1))
  }
  if (model == "IV") {
    Y0 = Y0
    Z1 = cbind(Trend, Z1)
    Z2 = cbind(Z2, lT)
  }
  if (model == "V") {
    Y0 = Y0
    Z1 = cbind(Z1)
    Z2 = cbind(Z2, lT, Trend)
  }
  if (ncol(Z2)>0) {
    Z2 <-  ShiftZ2m(Z2,kz,N1,max(P[,1:2])-1,min(P[,1:2])-1)

    M00 <- crossprod(Y0)/T1
    M11 <- crossprod(Z1)/T1
    M22 <- crossprod(Z2)/T1
    M01 <- crossprod(Y0, Z1)/T1
    M02 <- crossprod(Y0, Z2)/T1
    M10 <- crossprod(Z1, Y0)/T1
    M20 <- crossprod(Z2, Y0)/T1
    M12 <- crossprod(Z1, Z2)/T1
    M21 <- crossprod(Z2, Z1)/T1
    M22inv <- solve(M22)
    R0 <- Y0 - t(M02 %*% M22inv %*% t(Z2))
    R1 <- Z1 - t(M12 %*% M22inv %*% t(Z2))
  }

  if (ncol(Z2) == 0) {
    R0 = Y0
    R1 = Z1
  }

  S00 <- crossprod(R0)/T1
  S01 <- crossprod(R0, R1)/T1
  S10 <- crossprod(R1, R0)/T1
  S11 <- crossprod(R1)/T1
  Ctemp <- chol(S11, pivot = TRUE)
  pivot <- attr(Ctemp, "pivot")
  oo <- order(pivot)
  C <- t(Ctemp[, oo])
  Cinv <- solve(C)
  S00inv <- solve(S00)
  valeigen <- eigen(Cinv %*% S10 %*% S00inv %*% S01 %*% t(Cinv))
  lambda <- Re(valeigen$values)
  e <- Re(valeigen$vector)
  V <- t(Cinv) %*% e
  Vorg <- V
  V <- sapply(1:M1, function(j) V[, j]/V[1, j])
  PI <- S01 %*% solve(S11)
  GAMMA <- M02 %*% M22inv - PI %*% M12 %*% M22inv
  beta <- as.matrix(V[, 1:crk])
  beta0 <- as.matrix(V[, 1:crk])
  if (crk > 0) {
    NLm = stats::nlm(f, as.vector(beta0), beta, Z1, St, NSt, Y0, Z2)
    if ( NLm$code>2 )  betaS <- beta0  else   betaS <- NLm$estimate

    NLmcode = NLm$code
    dim(betaS) = dim(beta)
    for (i in 1:ncol(as.matrix(betaS))) betaS[, i] = betaS[, i]/(betaS[1, i])
    CI = Z1 %*% betaS
    estimation <- stats::lm(Y0 ~ 0 + CI + Z2)
    CI1 = CI * St
    CI2 = CI * NSt
    VECM1 <- stats::lm(Y0 ~ 0 + CI1 + CI2 + Z2)
    VECM1S <- summaryCIVAR(lm1 = VECM1, sname = "Z2")
    alphaS_1 = t(VECM1$coefficients[1:crk, ])
    alphaS_2 = t(VECM1$coefficients[(1 + crk):(2 * crk),])
    if (dim(alphaS_1)[1] == 1)
      alphaS_1 <- t(alphaS_1)
    if (dim(alphaS_2)[1] == 1)
      alphaS_2 <- t(alphaS_2)
  }
  if (crk == 0) {
    CI = 0
    estimation <- stats::lm(Y0 ~ 0 + Z2)
  }
  E = -T1 * log(1 - lambda)
  E = E[1:N]
  resultsvecm <- summary(estimation)
  if (model == "I") {
    Tab <- Tab1
  }
  if (model == "II") {
    Tab <- Tab2
  }
  if (model == "III") {
    Tab <- Tab3
  }
  if (model == "IV") {
    Tab <- Tab4
  }
  if (model == "V") {
    Tab <- Tab5
  }
  b = c(1:12)
  for (i in 1:12) {
    b[i] = Tab[2 * (i - 1) + 1, 2 + M1 - N1]
  }
  a = c(1:12)
  for (i in 1:12) {
    a[i] = Tab[2 * i, 2 + M1 - N1]
  }
  if (type == "eigen") {
    critical_vals = b
    M = matrix(0, N, 1)
    j = 1
    rank = 0
    while (j <= N && E[j] > critical_vals[N + 1 - j]) {
      M[j, ] = M[j, ] + 1
      j = j + 1
      rank = rank + 1
    }
    erg <- cbind(E, critical_vals[N:1])
    colnames(erg) <- c("teststatistic", "critical_value")
    if (N > 1) {
      rownames(erg) <- c("crk <= 0 |", paste("crk <= ",
                                             1:(N - 1), " |", sep = ""))
    }
    if (N == 1) {
      rownames(erg) <- c("crk <= 0 |")
    }
    coint_rank <- paste("Johansen-Test (with maximum-eigenvalue-teststatistic) indicates",
                        rank, "cointegrating equation(s) at the", 1 -
                          q, "level")
    result <- new.env()
    result$erg = erg
    result$estimation = estimation
    result$VECM1 = VECM1
    result$VECM1S = VECM1S
    result$lambda = E
    result$z = z
    result$Z2 = Z2
    result$Z1 = Z1
    result$St = St
    result$NSt = NSt
    result$R0 = R0
    result$R1 = R1
    result$Y0 = Y0
    result$beta = beta
    result$betaS = betaS
    result$alphaS = alphaS = list(alphaS_1, alphaS_2)
    result$NLmcode = NLmcode
    result$PI = PI
    result$GAMMA = GAMMA
    result$model = model
    result$P = P
    result$NN1 = NN1
    result$s = s
    rr <- as.list(result)
    return(rr)
  }
  else {
    type = "trace"
    critical_vals = a
    stat = matrix(0, N, 1)
    for (i in 1:N) {
      sum = 0
      for (j in i:N) {
        sum = sum + E[j]
      }
      stat[i] = sum
    }
    M = matrix(0, N, 1)
    j = 1
    rank = 0
    while (stat[j] > critical_vals[N + 1 - j] && j <= N) {
      M[j, ] = M[j, ] + 1
      j = j + 1
      rank = rank + 1
    }
    erg <- cbind(stat, critical_vals[N:1])
    colnames(erg) <- c("teststatistic", "critical_value")
    if (N > 1) {
      rownames(erg) <- c("crk <= 0 |", paste("crk <= ",
                                             1:(N - 1), " |", sep = ""))
    }
    if (N == 1) {
      rownames(erg) <- c("crk <= 0 |")
    }
    coint_rank <- paste("Johansen-Test (with trace-teststatistic) indicates",
                        rank, "cointegrating equation(s) at the", 1 -
                          q, "level")
    result <- new.env()
    result$erg = erg
    result$estimation = estimation
    result$VECM1 = VECM1
    result$VECM1S = VECM1S
    result$lambda = E
    result$z = z
    result$Z2 = Z2
    result$Z1 = Z1
    result$St = St
    result$NSt = NSt
    result$R0 = R0
    result$R1 = R1
    result$Y0 = Y0
    result$beta = beta
    result$betaS = betaS
    result$alphaS = alphaS = list(alphaS_1, alphaS_2)
    result$NLmcode = NLmcode
    result$PI = PI
    result$GAMMA = GAMMA
    result$model = model
    result$P = P
    result$NN1 = NN1
    result$s = s
    rr <- as.list(result)
    return(rr)
  }
}


#' This function replaces the variable names in an lm output.

#' @param lm1 The coefficient matrix of the output of lm
#' @param sname The elements in column names to be replaced
#' @export
summaryCIVAR <- function(lm1=lm1,sname ="Z2\\[,4:9\\]") {
  #### this function replace the Matrix name and leave the colnames in the output
  summ =summary(lm1)
  n = length(summ)
  if (n>1) {
    for (i in 1:n) {
      LL <- dimnames(summ[[i]]$coefficients)[[1]]
      dimnames(summ[[i]]$coefficients)[[1]] <- sub(sname, "", LL)
    }
  }  else {
    LL <- dimnames(summ$coefficients)[[1]]
    dimnames(summ$coefficients)[[1]] <- sub(sname, "", LL)

  }
  return(summ)
}



#' Test on restrictions in the cointegration space
#'
#' This function estimates constrained VECM using the iteration procedure proposed in Boswijk and Doornik (2003)
#'
#' @param R0    : T x n matrix containing the residuals of the auxiliary regression of dy_t on Z2_t
#' @param R1    : T x n matrix containing the residuals of the auxiliary regression of y_(t-1) on Z2_t
#' @param G     : The restriction matrix on alpha
#' @param H     : The restriction matrix on beta
#' @param h     : The restriction vector on beta ( free-varying parameters of beta )
#' @param alpha : n x crk matrix containing initial values of the adjustment parameters
#' @param beta  : n x crk matrix containing initial values of the cointegration vectors
#' @param Omega : n x n matrix containing initials of the covariance matrix
#' @param Y0    : T x n matrix of dY_t
#' @param Z2    : T x (n(p-1)+#det) matrix containing the lagged stationary regressors of the auxiliary regression
#' @param Z1    : T x n matrix containing the lagged cointegrated VAR process.
#'
#'
#' @section Details:
#' This function runs a likelihood ratio test of linear restrictions on \eqn{\alpha} and \eqn{\beta} in a CIVAR model in the following form:
#'		\deqn{vec(\alpha') = G \psi,  vec(\beta) = H\phi + h}
#'
#'        example 1 (restrictions on alpha) test of exogeneity
#'                  vec(alpha) is 45 x 1  vector  ( n = 9 crk = 5, defined by the model )
#'                  G          is 45 x 40 matrix  ( the first variable is exogenous, i.e. the 5 adjustment coefficient of the first variable are zero )
#'                  psi        40 x 1 vector (free varying parameters not appearing in the specification but implied by G)
#'                  vec(beta)  is 45 x 1 vector   ( N = 0 crk = 5 )
#'                  H          is 45 x 45 identity matrix
#'                  phi        45 x 1 vector (free varying parameters not appearing in the specification but implied by h =0
#'                  h          45 x 1 zero matrix implying ver(beta) = phi  >> no restrictions on beta.
#'                             (H is identity and h is zero vector implies only restrictions on alpha)
#'
#'        example 2 (restrictions on beta ) test of PPP
#'                  vec(alpha) is 40 x 1  vector  ( N = 8 crk = 5 ) conditional VECM  or VECMX model
#'                  G          is 40 x 40 identity matrix, implying there is no restriction on alpha
#'                  psi        40 x 1 vector (free varying parameters not appearing in the specification but implied by the identity matrix G
#'                  vec(beta)  is 45 x 1 vector   ( N = 8 crk = 5, defined by the model )
#'                  H          is 45 x 2 matrix that picks out the elements under restrictions ( two columns out of the identity matrix ) a zero row in H and the corresponding h implies zero-restrictions on beta.
#'                             ones in a row of H and zero in the corresponding h implies non-restricted beta.
#'                  phi        2  x 1  vector (free varying parameters not appearing in the specification but implied by h =0
#'                  h          45 x 1  non zero elements in this vector together with the zero elements in the corresponding row in H are the normalization conditions.
#'                             (H is identity and h is zero vector implies only restrictions on alpha)
#'
#' @return        a list containing (VECMRS, beter, alphar, LSKOEFR, LR, and error)
#'
#' @export
AB_CIVARTest <- function(R0,R1,G,H,h,alpha,beta,Omega,Y0,Z1,Z2) {
  ### this function runs Doornik Boswijk Iteration to estimate restricted alpha and beta
  ### It also estimates the separate MRCIVAR and compare it with the restricted including the model with identical beta.
  ### Add output of of separate Models.
  n   = dim(alpha)[1]
  crk = dim(alpha)[2]
  T1 = nrow(R1)
  S00 <- crossprod(R0)/T1
  S01 <- crossprod(R0, R1)/T1
  S10 <- crossprod(R1, R0)/T1
  S11 <- crossprod(R1)/T1

  alphai = alpha
  betai = beta
  Omegai = Omega
  VS10 = S10
  dim(VS10) <- c(length(S10),1)
  error = 1
  tol   = 0.0000001
  vecpip = solve(S11)%*%S10;dim(vecpip) = c(length(vecpip),1)

  phii = solve(t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%H)%*%
    (t(H)%*%kronecker(t(alphai)%*%solve(Omegai),diag(nrow(beta)))%*%VS10
     - t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%h)

  betai <-H%*%phii+h ; dim(betai) <- c(nrow(beta),crk)

  gammai = solve(t(G)%*%kronecker(solve(Omegai),t(betai)%*%S11%*%betai)%*%G)%*%
    t(G)%*%kronecker(solve(Omegai),t(betai))%*%VS10

  talphai <- G%*%gammai;dim(talphai) <- c(crk,nrow(alpha))
  alphai = t(talphai)

  Omegai = S00 - alphai%*%t(betai)%*%S10-S01%*%betai%*%t(alphai)+alphai%*%t(betai)%*%S11%*%betai%*%t(alphai)


  while  ( error > tol ) {

    phir   = phii
    gammar = gammai
    Omegar = Omegai

    phii = solve(t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%H)%*%
      (t(H)%*%kronecker(t(alphai)%*%solve(Omegai),diag(nrow(beta)))%*%VS10
       - t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%h)

    betai <-H%*%phii+h ; dim(betai) <- c(nrow(beta),crk)

    gammai = solve(t(G)%*%kronecker(solve(Omegai),t(betai)%*%S11%*%betai)%*%G)%*%
      t(G)%*%kronecker(solve(Omegai),t(betai))%*%VS10

    talphai <- G%*%gammai;dim(talphai) <- c(crk,nrow(alpha))

    alphai = t(talphai)

    Omegai = S00 - alphai%*%t(betai)%*%S10-S01%*%betai%*%t(alphai)+alphai%*%t(betai)%*%S11%*%betai%*%t(alphai)
    error = max( abs(log(det(Omegai))-log(det(Omegar))) )

  }
  LR =  T1*(log(det(Omegar))-log(det(Omega)))


  CI = Z1%*%betai

  if ((crk>0)&(length(Z2)>0)) {
    VECMR<-stats::lm(Y0-Z1%*%betai%*%t(alphai)~0+Z2)
    VECMRS <- summaryCIVAR(lm1=VECMR,sname ="Z2")
    LSKOEFR = VECMR$coefficients
  }
  if ((crk>0)&(length(Z2)==0)) {
    VECMR  <- NULL
    LSKOEFR =  NULL
  }
  if ((crk == 0)&(length(Z2)>0)) {
    VECMR <-stats::lm(Y0~0+Z2)
    LSKOEFR = VECMR$coefficients
  }
  if ((crk == 0)&(length(Z2)==0)) {
    VECMR<-NULL
    LSKOEFR<-NULL
  }

  result <- new.env()
  result$VECMR      = VECMR
  result$VECMRS     = VECMRS
  result$betar      = betai
  result$alphar     = alphai
  result$LSKOEFR    = LSKOEFR
  result$LR         = LR
  result$error      = error
  rr<-as.list(result)
  return(rr)
}




#' Test of restrictions in the cointegration space and estimation of restricted VECM
#'
#' This function estimates constrained VECM using the iteration procedure proposed in Boswijk and Doornik (2003)
#'
#' @param res an estimated CIVAR object with without restrictions
#' @param H the restriction matrix on beta
#' @param h the restriction vector on beta ( free-varying parameters of beta )
#' @param phi the restriction vector on beta
#' @param G the restriction matrix on alpha
#' @param Dxflag A flag that indicates whether X enters the cointegration space.
#'
#' @section Details:
#' This function runs a likelihood ratio test of linear restrictions on \eqn{\alpha} and \eqn{\beta} in a CIVAR model in the following form:
#'		\deqn{vec(\alpha') = G \psi,  vec(\beta) = H\phi + h}
#'
#'        example 1 (restrictions on alpha) test of exogeneity
#'                  vec(alpha) is 45 x 1  vector  ( n = 9 crk = 5, defined by the model )
#'                  G          is 45 x 40 matrix  ( the first variable is exogenous, i.e. the 5 adjustment coefficient of the first variable are zero )
#'                  psi        40 x 1 vector (free varying parameters not appearing in the specification but implied by G)
#'                  vec(beta)  is 45 x 1 vector   ( N = 0 crk = 5 )
#'                  H          is 45 x 45 identity matrix
#'                  phi        45 x 1 vector (free varying parameters not appearing in the specification but implied by h =0
#'                  h          45 x 1 zero matrix implying ver(beta) = phi  >> no restrictions on beta.
#'                             (H is identity and h is zero vector implies only restrictions on alpha)
#'
#'        example 2 (restrictions on beta ) test of PPP
#'                  vec(alpha) is 40 x 1  vector  ( N = 8 crk = 5 ) conditional VECM  or VECMX model
#'                  G          is 40 x 40 identity matrix, implying there is no restriction on alpha
#'                  psi        40 x 1 vector (free varying parameters not appearing in the specification but implied by the identity matrix G
#'                  vec(beta)  is 45 x 1 vector   ( N = 8 crk = 5, defined by the model )
#'                  H          is 45 x 2 matrix that picks out the elements under restrictions ( two columns out of the identity matrix ) a zero row in H and the corresponding h implies zero-restrictions on beta.
#'                             ones in a row of H and zero in the corresponding h implies non-restricted beta.
#'                  phi        2  x 1  vector (free varying parameters not appearing in the specification but implied by h =0
#'                  h          45 x 1  non zero elements in this vector together with the zero elements in the corresponding row in H are the normalization conditions.
#'                             (H is identity and h is zero vector implies only restrictions on alpha)
#'
#' @return a list containing (VECMRS, beter, alphae, LSKOEFR, LR and error)
#' @examples
#' X = matrix(stats::rnorm(2*207),207,2)
#' colnames(X) = c("ex1","ex2")
#' res_d  = CIVARData(n=6,p=2,T=207,crk=2,type="exog1", X=as.matrix(X))
#' colnames(res_d$Y) =  c("w","p","U_l","r","yn","y")
#' res_e =  CIVARest(res=res_d)
#' res_e$tst$erg
#' res_e$tst$beta
#' res_e$tst$VECMS
#' n = 6;  crk = 2
#' G = diag(n*crk); G0 = G; psi = matrix(1,n*crk,1); psi0=psi;   ### No restrictions on alpha
#' H = diag(n*crk); H2 = H[,-seq(1,n*crk,n)];
#' h2 = matrix(0,n*crk,1);
#' h2[seq(1,n*crk,n),1] = rep(1,crk);
#' phi2 = matrix(1,(n-1)*crk,1)                               ### no restrictions but normalization
#'
#' G0%*%psi0; H2%*%phi2+h2
#'
#' ABtest = CIVARTest(res=res_d,H=H2,h=h2,phi=phi2,G=G0,Dxflag=0)
#' ABtest$betar
#' ABtest$alphar
#' ABtest$VECMR$coefficients
#' ABtest$LR
#' #1-pchisq(ABtest$LR,2)   ### the fourth is the used restrictions
#'
#'
#' res_d  = CIVARData(n=7,p=2,T=207,crk=2,type="const")
#' colnames(res_d$Y) =  c("w","p","U_l","r","yn","y","fsi")
#' res_e =  CIVARest(res=res_d)
#' res_e$tst$erg
#' res_e$tst$beta
#' res_e$tst$VECMS
#'
#' n = 7;  crk = 2
#' G = diag(n*crk); G0 = G; psi = matrix(1,n*crk,1); psi0=psi;   ### No restrictions on alpha
#'
#' H = diag(n*crk); H2 = H[,-c(1,8)]; h2 = matrix(0,n*crk,1);
#' h2[c(1,8),1] = c(1,1);  phi2 = matrix(1,12,1)                 ### normalization
#'
#' G0%*%psi0; H2%*%phi2+h2
#'
#' ABtest = CIVARTest(res=res_d,H=H2,h=h2,phi=phi2,G=G0,Dxflag=1)
#' ABtest$LR
#' ABtest$betar
#' ABtest$alphar
#' #1-pchisq(ABtest$LR,4)
#'
#'
#' @export
CIVARTest <- function(res=res,H=H,h=h,phi=phi,G=G,Dxflag=0) {
###
### This function runs a likelihood ratio test of linear restrictions on alpha and beta in a CIVAR model
###
###             vec(alpha') = G psi , vec(beta) = H phi + h
###
###        example 1 (restrictions on alpha) test of exogeneity
###                  vec(alpha) is 45 x 1  vector  ( N = 9 crk = 5, defined by the model )
###                  G          is 45 x 40 matrix  ( the first variable is exogeneous, i.e. the 5 adjustment coefficient of the first variable are zero )
###                  psi        40 x 1 vector (free variaring parameters not appearing in the specification but implied by G)
###                  vec(beta)  is 45 x 1 vector   ( N = 0 crk = 5 )
###                  H          is 45 x 45 identity matrix
###                  phi        45 x 1 vector (free variaring parameters not appearing in the specification but implied by h =0
###                  h          45 x 1 zero matrix implying ver(beta) = phi  >> no restrictions on beta.
###                            (H is identity and h is zero vector implies only restrictions on alpha)
###
###        example 2 (restrictions on beta ) test of PPP
###                  vec(alpha) is 40 x 1  vector  ( N = 8 crk = 5 ) conditioal VECM  or VECMX model
###                  G          is 40 x 40 identity matrix, implying there is no restriction on alpha
###                  psi        40 x 1 vector (free variaring parameters not appearing in the specification but implied by the identity matrix G
###                  vec(beta)  is 45 x 1 vector   ( N = 8 crk = 5, defined by the model )
###                  H          is 45 x 2 matrix that picks out the elements under restricitons ( two colunms out of the identity matrix ) a zero row in H and the corresponding h implies zero-restrictions on beta.
###                             ones in a row of H and zero in the corresponding h implies non-restricted beta.
###                  phi        2  x 1  vector (free variaring parameters not appearing in the specification but implied by h =0
###                  h          45 x 1  non zero elements in this vector together with the zero elements in the corresponding row in H are the normalization conditions.
###                             (H is identity and h is zero vector implies only restrictions on alpha)
###
###
###
if (Dxflag==0) {
  res_e = CIVARest(res=res)
  crk   = res_e$crk
  R0    = res_e$tst$R0
  R1    = res_e$tst$R1
  alpha = t(res_e$tst$estimation$coefficients[1:crk,])
  beta  = res_e$tst$beta
  Omega = res_e$tst$Omega
  Y0 = res_e$tst$Y0
  Z1 = res_e$tst$Z1
  Z2 = res_e$tst$Z2
  test <- AB_CIVARTest(R0,R1,G,H,h,alpha,beta,Omega,Y0,Z1,Z2)
}   else  {
y = res$Y
if (is.na(res$X)) x = 0 else x = res$X

if (res$type == "const") model = "III"
if (res$type == "none")  model = "I"
bnorm = "1"
type  = "eigen"
p     =  res$p
r     =  res$crk
### N = 7, crk = 5, alphe beta are N  x crk  = 7 x 5   H = 35x35 G 35x35
###  G = diag(35); h = matrix(0,35,1); h[c(1,8,15,22,29)] = c(1,1,1,1,1)  ; H = diag(35); H2 = H[,-c(1,8,15,22,29)]
test<-abLRtest2(y,x,model = model, bnorm = bnorm,type = type,p = p, r=r, q = 0.95, H=H,h=h,G=G)
}
return(test)
}


#' Estimation of conditional VECM
#'
#' This function runs a test in cointegration space and provides an estimation of restricted alpha and beta.
#'
#' @param y the endogenous variables
#' @param x the exogenous variables
#' @param model Types of the deterministic specification
#' @param bnorm normalization of the cointegration vectors
#' @param type types of returns
#' @param p Lags of the CIVAR model
#' @param r the cointegration rank
#' @param q the significance level
#' @param H restriction matrix of beta
#' @param h restriction vector of beta
#' @param G restriction matrix of alpha
#' @export
abLRtest2 =function (y,x,model = c("I","II","III","IV","V"), bnorm = c("1","2"),type = c("eigen", "trace"),p = 1, r = 2, q = 0.95, H=H,h=h,G=G)   {
###
### This function runs a likelihood ratio test of linear restrictions on alpha and beta in a CIVAR model
###
###		vec(alpha) = G psi , vec(beta) = H phi + h
###
###        example 1 (restrictions on alpha) test of exogeneity
###			   vec(alpha) is 45 x 1  vector  ( N = 9 crk = 5, defined by the model )
###                  G          is 45 x 40 matrix  ( the first variable is exogeneous, i.e. the 5 adjustment coefficient of the first variable are zero )
###		         psi        40 x 1 vector (free variaring parameters not appearing in the specification but implied by G)
###                  vec(beta)  is 45 x 1 vector   ( N = 0 crk = 5 )
###                  H          is 45 x 45 identity matrix
###                  phi        45 x 1 vector (free variaring parameters not appearing in the specification but implied by h =0
###                  h          45 x 1 zero matrix implying ver(beta) = phi  >> no restrictions on beta.
###			   (H is identity and h is zero vector implies only restrictions on alpha)
###
###        example 2 (restrictions on beta ) test of PPP
###  			   vec(alpha) is 40 x 1  vector  ( N = 8 crk = 5 ) conditioal VECM  or VECMX model
###                  G          is 40 x 40 identity matrix, implying there is no restriction on alpha
###		         psi        40 x 1 vector (free variaring parameters not appearing in the specification but implied by the identity matrix G
###                  vec(beta)  is 45 x 1 vector   ( N = 0 crk = 5, defined by the model )
###                  H          is 45 x 2 matrix that picks out the elements under restricitons ( two colunms out of the identity matrix ) a zero row in H and the corresponding h implies zero-restrictions on beta.
###                             ones in a row of H and zero in the corresponding h implies non-restricted beta.
###                  phi        2  x 1  vector (free variaring parameters not appearing in the specification but implied by h =0
###                  h          45 x 1  non zero elements in this vector together with the zero elements in the corresponding row in H are the normalization conditions.
###			   (H is identity and h is zero vector implies only restrictions on alpha)
###
###

    y <- as.matrix(y)
    model <- match.arg(model)
    type <- match.arg(type)
    #ret <- match.arg(ret)
    #ctable <- match.arg(ctable)
    if (q != 0.9 && q != 0.95 && q != 0.99) {
        return("please correct significance level")
    }
    ret = "statistics"

    p <- as.integer(p)
    N1 <- ncol(as.matrix(y))
    #N  <- ncol(as.matrix(y))
    if (N1<r) {
        return("y's dimension must be larger than r")
    }

    if  ( length(x) == 1 ) {  z <- y }
    else {
       z <- cbind(y,x)
    }

    M1 <- ncol(as.matrix(z))
    T <- nrow(as.matrix(z))
    ## Z <- embed(diff(z), p)
    Z <- Embed(diff(z), p)
    Y0 <- Z[, c(1:N1)]                # \Delta Y
    #X0 <- Z[,c((N1+1):M1)]           # \Delta X
    Z1 <- z[-T, ][p:(T - 1), ]        #      Z_1
    if ((p==1)&(length(x)==1)) {
       Z2 = NULL
       NOLAG = 1
    }
    else {
       Z2 <- Z[, c((N1+1):(M1*p))]       # \Delta Z_
    }
    T1 <- nrow(as.matrix(Y0))
    #lT = (1:T1)/(1:T1)
    lT = matrix(1,T1,1); colnames(lT) <- "Const"
    #Trend <- matrix((p + 1):TT, (TT - p), 1); colnames(Trend) <- "Trend"
    Trend <- matrix(1:T1, T1, 1); colnames(Trend) <- "Trend"
    if ( model == "I" ) {
    Y0 = Y0
    Z1 = Z1
    Z2 = Z2
    }
    if ( model == "II" ) {
    Y0 = Y0
    Z1 = cbind(lT,Z1)
    Z2 = Z2
    }
    if ( model == "III" ) {
    Y0 = Y0
    Z1 = Z1
    Z2 = cbind(Z2,lT)
    }
    if ( model == "IV" ) {
    Y0 = Y0
    Z1 = cbind(Trend,Z1)
    Z2 = cbind(Z2,lT)
    }
    if ( model == "V" ) {
    Y0 = Y0
    Z1 = cbind(Z1)
    Z2 = cbind(Z2,lT,Trend)
    }

    if (length(Z2)==0) {
        R0 <- Y0
        R1 <- Z1
    }
    else {

        M00 <- crossprod(Y0)/T1
        M11 <- crossprod(Z1)/T1
        M22 <- crossprod(Z2)/T1
        M01 <- crossprod(Y0, Z1)/T1
        M02 <- crossprod(Y0, Z2)/T1
        M10 <- crossprod(Z1, Y0)/T1
        M20 <- crossprod(Z2, Y0)/T1
        M12 <- crossprod(Z1, Z2)/T1
        M21 <- crossprod(Z2, Z1)/T1
        M22inv <- solve(M22)
        R0 <- Y0 - t(M02 %*% M22inv %*% t(Z2))
        R1 <- Z1 - t(M12 %*% M22inv %*% t(Z2))
    }
    S00 <- crossprod(R0)/T1
    S01 <- crossprod(R0, R1)/T1
    S10 <- crossprod(R1, R0)/T1
    S11 <- crossprod(R1)/T1
    Ctemp <- chol(S11, pivot = TRUE)
    pivot <- attr(Ctemp, "pivot")
    oo <- order(pivot)
    C <- t(Ctemp[, oo])
    Cinv <- solve(C)
    S00inv <- solve(S00)
    valeigen <- eigen(Cinv %*% S10 %*% S00inv %*% S01 %*% t(Cinv))

    #valeigen <- eigen(S11inv %*% S10 %*% S00inv %*% S01 )
    lambda <- valeigen$values
    e      <- valeigen$vector

    V <- t(Cinv) %*% e
    Vorg <- V
    V <- sapply(1:M1, function(j) V[, j]/V[1, j])
    W <- S01 %*% V %*% solve(t(V) %*% S11 %*% V)
    PI <- S01 %*% solve(S11)
    DELTA <- S00 - S01 %*% V %*% solve(t(V) %*% S11 %*% V) %*% t(V) %*% S10
    if (length(Z2)==0) {
       GAMMA = NULL
    }
    else {
       GAMMA <- M02 %*% M22inv - PI %*% M12 %*% M22inv
    }
    if (bnorm == "1") {
        beta   <- as.matrix(V[,1:r])
    }
    if (bnorm == "2") {
        beta = as.matrix(V[,1:r]%*%solve(V[1:r,1:r]))
    }

    if ((r>0)&(length(Z2)>0)) {
       CI = Z1%*%beta
       VECM<-stats::lm(Y0~0+CI+Z2)
       if (r==1)  { alpha = as.matrix(VECM$coefficients[1:1,]) }
       if (r > 1) { alpha = t(as.matrix(VECM$coefficients[1:r,])) }
       LSKOEF = stats::lm(Y0-Z1%*%beta%*%t(alpha)~0+Z2)$coefficients

    }
    if ((r>0)&(length(Z2)==0)) {
       CI = Z1%*%beta
       VECM<-stats::lm(Y0~0+CI)
       LSKOEF = NULL
    }
    if ((r == 0)&(length(Z2)>0)) {
       CI = 0
       VECM<-stats::lm(Y0~0+Z2)
       LSKOEF = stats::lm(Y0~0+Z2)$coefficients
    }
    if ((r == 0)&(length(Z2)==0)) {
       CI = 0
       VECM<-NULL
       LSKOEF = NULL
    }
    E = -T1 * log(1 - lambda)
    E = E[1:N1]
    LOGL = -T1/2.0*(log(det(S00))+sum(log((1-lambda)[1:r])))-T1*N1/2.0*log(2*pi)-T1*N1/2.0
    #LOGL = -T1/2.0*(log(det(S00))+sum(log((1-lambda)[1:N1])))-T1*N1/2.0*log(2*pi)-T1*N1/2.0

    Omega = 1.0/(nrow(VECM$residuals))*t(VECM$residuals)%*%VECM$residuals
    # T1*(log(det(Omega))-log(det(S00))-sum(log(1-lambda)[1:r]))   # correct

    ##### Estimation of VECM
    if (r==1)  { alpha = as.matrix(VECM$coefficients[1:1,]) }
    if (r > 1) { alpha = t(as.matrix(VECM$coefficients[1:r,])) }
    resultsvecm <- summary(VECM)

    if (model == "I") {
      Tab <- Tab1
    }
    if (model == "II") {
      Tab <- Tab2
    }
    if (model == "III") {
      Tab <- Tab3
    }
    if (model == "IV") {
      Tab <- Tab4
    }
    if (model == "V") {
      Tab <- Tab5
    }


    b = c(1:12)
    for (i in 1:12)   { b[i] = Tab[2*(i-1)+1,2+M1-N1] }
    a = c(1:12)
    for (i in 1:12)   { a[i] = Tab[2*i,2+M1-N1] }

    if (type == "eigen") {
        critical_vals = b
        M = matrix(0, N1, 1)
        j = 1
        rank = 0
        while (j <= N1 && E[j] > critical_vals[N1 + 1 - j]) {
            M[j, ] = M[j, ] + 1
            j = j + 1
            rank = rank + 1
        }
        if (ret == "test") {
            return(M)
        }
        erge <- cbind(E, critical_vals[N1:1])
        colnames(erge) <- c("teststatistic", "critical_value")

        if ( N1 > 1 ) {
           rownames(erge) <- c("r <= 0 |", paste("r <= ", 1:(N1 - 1),  " |", sep = ""))
        }
        if ( N1 == 1 ) {
           rownames(erge) <- c("r <= 0 |" )
        }

        coint_rank <- paste("Johansen-Test (with maximum-eigenvalue-teststatistic) indicates",
            rank, "cointegrating equation(s) at the", 1 - q,
            "level")
    }
    else {
        type = "trace"
        critical_vals = a
        stat = matrix(0, N1, 1)
        for (i in 1:N1) {
            sum = 0
            for (j in i:N1) {
                sum = sum + E[j]
            }
            stat[i] = sum
        }
        M = matrix(0, N1, 1)
        j = 1
        rank = 0
        while (stat[j] > critical_vals[N1 + 1 - j] && j <= N1) {
            M[j, ] = M[j, ] + 1
            j = j + 1
            rank = rank + 1
        }
        if (ret == "test") {
            return(M)
        }
        ergt <- cbind(stat, critical_vals[N1:1])
        colnames(ergt) <- c("teststatistic", "critical_value")

        if ( N1 > 1 ) {
           rownames(ergt) <- c("r <= 0 |", paste("r <= ", 1:(N1 - 1),  " |", sep = ""))
        }
        if ( N1 == 1 ) {
           rownames(ergt) <- c("r <= 0 |")
        }

        coint_rank <- paste("Johansen-Test (with trace-teststatistic) indicates",
        rank, "cointegrating equation(s) at the", 1 - q, "level")
    }
######### Calculation of restricted MLE for cointegrated system
   alphai = alpha
   betai = beta
   Omegai = Omega
   VS10 = S10
   dim(VS10) <- c(length(S10),1)
   error = 1
   tol   = 0.0001
   vecpip = solve(S11)%*%S10;dim(vecpip) = c(length(vecpip),1)

       phii = solve(t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%H)%*%
          (t(H)%*%kronecker(t(alphai)%*%solve(Omegai),diag(nrow(beta)))%*%VS10
           - t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%h)

       betai <-H%*%phii+h ; dim(betai) <- c(nrow(beta),r)

       gammai = solve(t(G)%*%kronecker(solve(Omegai),t(betai)%*%S11%*%betai)%*%G)%*%
          t(G)%*%kronecker(solve(Omegai),t(betai))%*%VS10

       #phii = solve(t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%H)%*%
       #   t(H)%*%kronecker(t(alphai)%*%solve(Omegai),S11)%*%(vecpip-kronecker(alphai,diag(nrow(beta)))%*%h)

       #betai <-H%*%phii+h ; dim(betai) <- c(nrow(beta),r)

       #gammai = solve(t(G)%*%kronecker(solve(Omegai),t(betai)%*%S11%*%betai)%*%G)%*%
          t(G)%*%kronecker(solve(Omegai),t(betai)%*%S11)%*%vecpip

       talphai <- G%*%gammai;dim(talphai) <- c(r,N1)
       alphai = t(talphai)

       Omegai = S00 - alphai%*%t(betai)%*%S10-S01%*%betai%*%t(alphai)+alphai%*%t(betai)%*%S11%*%betai%*%t(alphai)


   while  ( error > tol ) {

         phir   = phii
         gammar = gammai
         Omegar = Omegai

       phii = solve(t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%H)%*%
          (t(H)%*%kronecker(t(alphai)%*%solve(Omegai),diag(nrow(beta)))%*%VS10
           - t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%h)

       betai <-H%*%phii+h ; dim(betai) <- c(nrow(beta),r)

       gammai = solve(t(G)%*%kronecker(solve(Omegai),t(betai)%*%S11%*%betai)%*%G)%*%
          t(G)%*%kronecker(solve(Omegai),t(betai))%*%VS10

       talphai <- G%*%gammai;dim(talphai) <- c(r,N1)
       alphai = t(talphai)

       Omegai = S00 - alphai%*%t(betai)%*%S10-S01%*%betai%*%t(alphai)+alphai%*%t(betai)%*%S11%*%betai%*%t(alphai)
       error = max( abs(log(det(Omegai))-log(det(Omegar))) )

       #phii = solve(t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%H)%*%
       #   t(H)%*%kronecker(t(alphai)%*%solve(Omegai),S11)%*%(vecpip-kronecker(alphai,diag(nrow(betai)))%*%h)

       #betai <-H%*%phii+h ; dim(betai) <- c(nrow(beta),r)

       #gammai = solve(t(G)%*%kronecker(solve(Omegai),t(betai)%*%S11%*%betai)%*%G)%*%
       #   t(G)%*%kronecker(solve(Omegai),t(betai)%*%S11)%*%vecpip

       #talphai <- G%*%gammai;dim(talphai) <- c(r,N1)
       #alphai = t(talphai)

       #Omegai = S00 - alphai%*%t(betai)%*%S10-S01%*%betai%*%t(alphai)+alphai%*%t(betai)%*%S11%*%betai%*%t(alphai)
       #error = max( abs(log(det(Omegai))-log(det(Omegar))) )

       log(det(Omegai))
   }
   LR = T1*(log(det(Omegar))-log(det(S00))-sum(log(1-lambda)[1:r]))

   #CI = Z1%*%beta
   #stats::lm( Z1%*%t(alphai%*%t(betai))

    if ((r>0)&(length(Z2)>0)) {
       VECMR<-stats::lm(Y0-Z1%*%betai%*%t(alphai)~0+Z2)
       VECMRS <- summaryCIVAR(lm1=VECMR,sname ="Z2")
       LSKOEFR = VECMR$coefficients
    }
    if ((r>0)&(length(Z2)==0)) {
       VECMR  <- NULL
       LSKOEFR =  NULL
    }
    if ((r == 0)&(length(Z2)>0)) {
       VECMR <-stats::lm(Y0~0+Z2)
       LSKOEFR = VECMR$coefficients
    }
    if ((r == 0)&(length(Z2)==0)) {
       VECMR<-NULL
       LSKOEFR<-NULL
    }

        result <- new.env()
        if (type == "trace")  {result$ergt = ergt}
        if (type == "eigen")  {result$erge = erge}
        result$estimation = VECM
  	  result$VECMRS     = VECMRS
        result$lambda     = E
        result$z          = z
        result$Z2         = Z2
        result$beta       = beta
        result$alpha      = alpha
        result$PI         = PI
        result$GAMMA      = GAMMA
        result$model      = model
        result$betar      = betai
        result$alphar     = alphai
        result$LSKOEFR    = LSKOEFR
        result$LSKOEF     = LSKOEF
        result$LR         = LR
        result$error      = error
        rr<-as.list(result)
        return(rr)
}

######################## copy from previous





#' Transformation of the estimated parameters of a VECM into the parameters of the corresponging VAR
#'
#' @param param The parameter of an estimated VECM
#' @param beta The cointegration vectors
#' @param q A vector specifying different type of VECMs
#' @param s The indicator of multi regime VECMs
#' @param kz The number of exogenous common factors
#' @param Dxflag Indicator whether the foreign variables enter the cointegration space
#'
#' @return A list of (A, B, C) with B the auto regression parameter matrix, C the parameter of the deterministic components, and A the parameter matrix of foreign variables for CIGVAR.
#' @export
VECM2VAR <- function(param, beta, q = c(1, 2, 2, 2, 2),s=NA,kz=0,Dxflag=0)
{
  m = dim(param)[2]
  VECB = t(param)
  if (missing(s))  s = NA
  if (anyNA(s)) {
    if (!q[3] == 0) {   # for GVAR with foreign variables
      B = (1:(m * m * (q[2] + 1))) * 0
      dim(B) = c(m, m, (q[2] + 1))
      A = (1:(m * m * (q[3] + 1))) * 0
      D = NA
      dim(A) = c(m, m, (q[3] + 1))
      AB = VECB[, 1:q[1]] %*% t(beta)
      B[, , 1] = AB[, kz+(1:m)]
      if (Dxflag == 1 ) A[, , 1] = AB[,kz+(m + 1):(2 * m)]
      if (Dxflag == 0 ) A[, , 1] = A[, , 1]*0
      for (i in 2:(q[2] + 1)) B[, , i] = VECB[, q[1] + ((i - 2) * m + 1):((i - 2) * m + m)]
      for (i in 2:(q[3] + 1)) A[, , i] = VECB[, (q[1] + q[2] * m) + ((i - 2) * m + 1):((i - 2) * m + m)]

      B = CIB3B(B)
      A = CIA2A(A)
      if ( kz > 0 ) {
        D = (1:(m*kz*(q[2]+1)))*0; dim(D) =c(m,kz,q[2]+1)
        D[,,1] = AB[,1:kz]
        for (i in 2:(q[2] + 1)) D[, , i] = VECB[, dim(VECB)[2]-q[2]*kz + (((i - 2) * kz + 1):((i - 2) * kz + kz))]
        D = CIA2A(D)
      }
    }
    else {
      B = (1:(m * m * (q[2] + 1))) * 0
      dim(B) = c(m, m, (q[2] + 1))
      AB = VECB[, 1:q[1]] %*% t(beta)
      B[, , 1] = AB[, 1:m]
      if (q[2] > 0) {
        for (i in 2:(q[2] + 1)) B[, , i] = VECB[, q[1] +
                                                  ((i - 2) * m + 1):((i - 2) * m + m)]
        B <- CIB3B(CIB = B)
      }
      if (q[2] == 0) {
        B[, , 1] = B[, , 1] + diag(m)
      }
      A = NA
    }
    if (dim(param)[1] > q[1] + (q[2] + q[3]) * m + q[2]*kz) {
      C = as.matrix(t(param)[, (q[1] + (q[2] + q[3]) * m + 1):(dim(param)[1]-q[2]*kz) ])
    }
    else C = NA
  }
  else {
    if (length(q) == 5) {
      P = max(q[-1])
      B = (1:(m * m * (P + 1) * 2)) * 0
      dim(B) = c(m, m, (P + 1), 2)
      A = B
      AB = VECB[, 1:q[1]] %*% t(beta)
      B[, , 1, 1] = AB[, 1:m]
      A[, , 1, 1] = AB[, (m + 1):(2 * m)]
      B[, , 1, 2] = AB[, 1:m]
      A[, , 1, 2] = AB[, (m + 1):(2 * m)]
      for (i in 2:(q[2] + 1)) B[, , i, 1] = VECB[, q[1] +
                                                   ((i - 2) * m + 1):((i - 2) * m + m)]
      for (i in 2:(q[3] + 1)) A[, , i, 1] = VECB[, (q[1] +
                                                      q[2] * m) + ((i - 2) * m + 1):((i - 2) * m +
                                                                                       m)]
      for (i in 2:(q[4] + 1)) B[, , i, 2] = VECB[, (q[1] +
                                                      (q[2] + q[3]) * m) + ((i - 2) * m + 1):((i -
                                                                                                 2) * m + m)]
      for (i in 2:(q[5] + 1)) A[, , i, 2] = VECB[, (q[1] +
                                                      (q[2] + q[3] + q[4]) * m) + ((i - 2) * m + 1):((i -
                                                                                                        2) * m + m)]
      B[, , 1:(q[2] + 1), 1] = CIB3B(B[, , 1:(q[2] + 1),
                                       1])
      B[, , 1:(q[4] + 1), 2] = CIB3B(B[, , 1:(q[4] + 1),
                                       2])
      A[, , 1:(q[3] + 1), 1] = CIA2A(A[, , 1:(q[3] + 1),
                                       1])
      A[, , 1:(q[5] + 1), 2] = CIA2A(A[, , 1:(q[5] + 1),
                                       2])
      if (dim(param)[1] > q[1] + (sum(q[-1])) * m) {
        C = (1:(m * 2)) * 0
        dim(C) = c(m, 1, 2)
        C[, 1, 1] = as.matrix(t(param)[, q[1] + sum(q[-1]) *
                                         m + 1])
        C[, 1, 2] = as.matrix(t(param)[, q[1] + sum(q[-1]) *
                                         m + 2])
      }
      else {
        C = NA
      }
    }
    if (length(q) == 3) {
      P = max(q[-1])
      B = (1:(m * m * (P + 1) * 2)) * 0
      dim(B) = c(m, m, (P + 1), 2)
      AB = VECB[, 1:q[1]] %*% t(beta)
      B[, , 1, 1] = AB[, 1:m]
      B[, , 1, 2] = AB[, 1:m]
      for (i in 2:(q[2] + 1)) B[, , i, 1] = VECB[, q[1] +
                                                   ((i - 2) * m + 1):((i - 2) * m + m)]
      for (i in 2:(q[3] + 1)) B[, , i, 2] = VECB[, (q[1] +
                                                      q[2] * m) + ((i - 2) * m + 1):((i - 2) * m +
                                                                                       m)]
      B[, , 1:(q[2] + 1), 1] = CIB3B(B[, , 1:(q[2] + 1),
                                       1])
      B[, , 1:(q[3] + 1), 2] = CIB3B(B[, , 1:(q[3] + 1),
                                       2])
      A = NA
      if (dim(param)[1] > q[1] + (sum(q[-1])) * m) {
        C = (1:(m * 2)) * 0
        dim(C) = c(m, 1, 2)
        C[, 1, 1] = as.matrix(t(param)[, q[1] + sum(q[-1]) *
                                         m + 1])
        C[, 1, 2] = as.matrix(t(param)[, q[1] + sum(q[-1]) *
                                         m + 2])
      }
      else {
        C = NA
      }
    }
  }
  return(list(B, A, C, D))
}





#' Transformation of the coefficient matrix of a CIGVAR to the coefficient matrix of the corresponding VECM.
#'
#' @param  B  an (m,m,L) array containing the coefficients of the level VAR.
#' @return      A list containing three components.
#' \itemize{
#'    \item CIB        : the coefficients matrix of the error correction form
#'    \item alpha	     : the adjustment coefficients
#'    \item beta       : the cointegration vectors
#'    \item sm         : the discrepancy between the selected alpha beta and CIB.
#' }
#' @export
B2CIB <- function (B)
{
  CIB = B * 0
  L = dim(B)[3]
  n = dim(B)[1]
  for (i in L:1) for (j in L:i) CIB[, , i] = CIB[, , i] - B[, , j]
  CIB[, , 1] = -diag(n) - CIB[, , 1]
  P = eigen(CIB[, , 1])$vectors
  lambda = eigen(CIB[, , 1])$values
  P=P[,sort(Mod(lambda),decreasing = FALSE,index.return=TRUE)$ix]
  lambda <- lambda[sort(Mod(lambda),decreasing = FALSE,index.return=TRUE)$ix]

  k = 0
  for (i in 1:(n - 1)) {
    sm = sum(abs(CIB[, , 1] - P[, i:n] %*% diag(lambda[i:n]) %*% solve(P)[i:n, ]))
    if (sm < 1e-06)
      k = i
  }
  sm = sum(abs(CIB[, , 1] - as.matrix(P[, n:n]) %*% lambda[n:n] %*% t(as.matrix(solve(P)[n:n, ]))))
  if (sm < 1e-06)
    k = n
  if (k < n) {
    alpha = P[, k:n] %*% diag(lambda[k:n])
    beta = t(solve(P)[k:n, ])
    sm = sum(abs(CIB[, , 1] - alpha %*% t(beta)))
  }
  if (k == n) {
    alpha = as.matrix(P[, n:n] * lambda[n:n])
    beta = as.matrix(solve(P)[n:n, ])
    sm = sum(abs(CIB[, , 1] - alpha %*% t(beta)))
  }
  return(list(CIB, alpha, beta, sm))
}






#' Transformation of coefficient matrix of a VECM to the coefficients matrix of the corresponding CIGVAR in level.
#'
#' @param  tst  : an output of CIGVARest
#' @return      A list containing three components.
#' \itemize{
#'     \item B          : the coefficient matrices of the domestic variables
#'     \item A	        : the coefficient matrices of the foreign variables
#'     \item C          : the coefficient matrices of the deterministic components
#' }
#' @export
CIB2B = function(tst) {
  P  = tst$P
  r  = ncol(as.matrix(tst$beta))
  m  = ncol(tst[[2]]$coefficients)
  C  = matrix(0,m,1)
  C1 = C
  C2 = C
  mx = tst$NN1
  if (r==1) PI = as.matrix(tst[[2]]$coefficients[1:r,])%*%t(tst$beta)
  if (r>1 ) PI = t(tst[[2]]$coefficients[1:r,])%*%t(tst$beta)
  if (tst$model=="II"|tst$model=="IV") C1 = as.matrix(PI[,m+mx+1])
  PI = PI[,1:(m+mx)]
  BB = t(tst[[2]]$coefficients[(r+1):(r+(P[1]-1)*m),])
  BA = t(tst[[2]]$coefficients[(r+(P[1]-1)*m+1):(r+(P[1]-1)*m+(P[2]-1)*mx) ,])
  #BBAA  = t(tst[[2]]$coefficients[(r+1:(r+(P[1]-1)*m+(P[2]-1)*mx),])

  if (tst$model=="III"|tst$model=="IV") C2 = as.matrix(tst[[2]]$coefficients[nrow(tst[[2]]$coefficients),])
  C  = C1 + C2
  B  = (1:(m*m*P[1]))*0; dim(B) = c(m,m,P[1])
  A  = (1:(m*mx*P[2]))*0; dim(A) = c(m,mx,P[2])

  if ( P[1]==1 )  B[,,1] = PI[,1:m]+diag(m)
  if ( P[1]==2 ) {B[,,1] = PI[,1:m] + BB[,1:m] +diag(m); B[,,2] = - BB[,1:m]}

  if ( P[1]>2 ) {
    B[,,1] = PI[,1:m] + BB[,1:m] + diag(m)
    for (i in 2:(P[1]-1)) {
      B[,,i]    = - BB[,((i-2)*m+1):((i-2)*m+m)] + BB[,((i-1)*m+1):((i-1)*m+m)];
    }
    B[,,P[1]] = - BB[,((P[1]-2)*m+1):((P[1]-2)*m+m)]
  }


  if ( P[2]==1 )  A[,,1] = PI[,(m+1):(m+mx)]
  if ( P[2]==2 ) {A[,,1] = PI[,(m+1):(m+mx)] + BA[,1:mx]; A[,,2] = - BA[,1:mx]}

  if ( P[2]>2 ) {
    A[,,1] = PI[,(m+1):(m+mx)] + BA[,1:mx]
    for (i in 2:(P[2]-1)) {
      A[,,i]    = - BA[,((i-2)*mx+1):((i-2)*mx+mx)] + BA[,((i-1)*mx+1):((i-1)*mx+mx)];
    }
    A[,,P[2]] = - BA[,((P[2]-2)*mx+1):((P[2]-2)*mx+mx)]
  }

  return(list(B,A,C))
}


#' Transformation of the coefficients matrix of a VECM to the coefficients of the corresponding VAR in level
#'
#' @param CIB an (m,m,L) array of the coefficients matrices of a VECM with lag L-1.
#' @return B (m,m,L) array the coefficients of a VAR with lag L
#' @export
CIB3B = function(CIB) {
  BB = CIB*0
  L = dim(CIB)[3]
  n = dim(BB)[1]
  if (L==1) BB = CIB + diag(4)
  if (L==2) {BB[,,1] = CIB[,,1]+CIB[,,2]+diag(n); BB[,,2] = - CIB[,,2]}
  if (L>2)  {
    BB[,,1] = CIB[,,1]+CIB[,,2]+diag(n)
    for (i in 2:(L-1))  BB[,,i] = CIB[,,i+1]-CIB[,,i]
    BB[,,L] = - CIB[,,L]
  }
  return(BB)
}







#' Transformation of VECM coefficient matrix to VAR coefficient matrix
#'
#' @param CIA a coefficient matrix of VECM
#'
#' @return the coefficient matrix of the corresponding VAR in level
#' @export
#'
CIA2A = function(CIA) {
	AA = CIA*0
	L = dim(CIA)[3]
      n = dim(AA)[1]
	if (L==1) AA = CIA
	if (L==2) {AA[,,1] = CIA[,,1]+CIA[,,2]; AA[,,2] = - CIA[,,2]}
	if (L>2)  {
     		AA[,,1] = CIA[,,1]+CIA[,,2]
     		for (i in 2:(L-1))  AA[,,i] = CIA[,,i+1]-CIA[,,i]
     		AA[,,L] = - CIA[,,L]
	}
return(AA)
}



#' Root of the characteristic polynomial in Lag
#'
#' @param G a coefficient matrix of dimension (n x n x p) of an VAR(p) object.
#'
#' @return Roots of the characteristic polynomial of VAR(p)
#' @export
STAT = function(G) {
   L = dim(G)[3]
   m = dim(G)[1]
   PP=matrix(0,m*L,m*L)
   PP[1:m,] = G[,,]
   if (L>1 ) {
     for ( i in 1:(L-1) ) PP[(i*m+1):(i*m+m),((i-1)*m+1):((i-1)*m+m)] = diag(m)
     eigen(PP)$values
   }  else {
     eigen(G[,,1])$values
   }
}


#' Impulse response of an estimated CIVAR.
#'
#' This function generates the impulse response functions of an estimated CIVAR model
#'
#' @param res	: a CIVAR object such as an output of CIVARest.
#' @param nstep	: the length of the impulse response functions.
#' @param comb  : a weighting matrix specifying the weights used in the impulse response functions of a global VAR. Its default value is NA for CIVAR(p).
#' @param irf   : types of the generated impulse response functions.
#' @param G     : the transformation matrix for PTdecomp
#' @param A0    : the transformation matrix for AB identification
#' @param B0    : the transformation matrix for AB identification
#' @param Xshks : the number of selected exogenous shocks
#' @param runs  : number of runs used in the the calculation of the bootstrap confidence interval.
#' @param conf  : a two component vector containing the tail probabilities of the bootstrap confidence interval.
#' @return       an array of dimension (n, n, nstep, 3).
#' @examples
#'
#' res_d = CIVARData(n=4,p=2,T=84,Co=matrix(c(1,1,1,1),4,1)*0,type="none",crk=1)
#' res_e = CIVARest(res=res_d)
#' res_e$Summary
#'
#' IRF_CB = irf_CIVAR_CB(res=res_e, nstep=20, comb=NA, irf = "gen1", runs = 20, conf = c(0.05, 0.95))
#' IRF_g = IRF_graph(IRF_CB)
#'
#' IRF_CB     = irf_CIVAR_CB(res=res_e, nstep=30, comb=NA, irf = "PTdecomp", G = NA, A0=NA,B0=NA,
#'  runs = 20, conf = c(0.05, 0.95))
#' IRF_g = IRF_graph(IRF_CB)
#' # The first three shocks have permanent effects,
#' # while the fourth shock does not have a permanent effect.
#'
#'
#'
#' T=100
#' X=matrix(rnorm(2*T),T,2)
#' res_d = VARData(n=3,p=2,T=T,Co=matrix(c(0,0,0,1,2,3,3,2,1),3,3), type="exog0",X=X) ;
#' res_e = VARest(res=res_d);
#' res_d$Co
#' res_e$Summary
#' IRF_CB = irf_VAR_CB(res=res_e,nstep=20, comb=NA, irf = "irfX", Xshks=c(1:2),
#' runs = 100, conf = c(0.05, 0.95))
#' IRF_list <-IRF_graph(IRF_CB,Names =c("Y1","Y2","Y3"),INames=c("X1","X2"),
#' response = c(1:3), impulse = c(1:3), n = 3)
#'
#' @export
irf_CIVAR_CB   <- function (res, nstep, comb, irf = c("gen", "chol", "chol1", "gen1", "comb1","PTdecomp","ABSVAR","irfX"), G=NA, A0=NA, B0=NA,Xshks=NA,runs = 200, conf = c(0.05, 0.95))
{
  n     = res$n
  p     = res$p
  T     = dim(res$Y)[1]
  A     = res$A
  B     = res$B
  Co    = res$Co
  C1    = res$C1
  type  = res$type
  X     = res$X
  crk   = res$crk
  neq   = dim(B)[1]
  nvar  = dim(B)[2]
  sigma = res$Sigma
  Uo    = res$resid
  G = NA
  if (irf=="PTdecomp") {
    crk = res$crk

    alpha  =  (as.matrix(res$tst$estimation$coefficients[1:crk,]))
    if (dim(alpha)[1]<dim(alpha)[2]) {
      alpha = t(alpha)
    }
    G     = t(cbind(alpha_perpf(alpha),alpha))
  }
  if (irf == "irfX") { Xshk = t(Co[,1+Xshks]) }
  response <- array(0, dim = c(neq, nvar, nstep, length(conf) + 1))
  response[, , , 1] <- irf_B_sigma(B = B, sigma = sigma, nstep, comb, irf, G=G,A0=A0,B0=B0,smat=NA,Xshk=Xshk)
  responseR         <- array(0, dim = c(neq, nvar, nstep, runs))




  for (i in 1:runs) {
    Uo_run = rnormSIGMA(T, sigma)
    #Index = c(1:T)
    #resample  = sample(Index)
    #Uo_run    = Uo[resample,]
    res_drun  = CIVARData(n, p, T, r_np = NA, A = NA, B, Co, C1, U = Uo_run, Sigma = NA, type, X, mu = NA, Yo = NA, crk)
    res_erun  = CIVARest(res_drun)
    B_run     = res_erun$B
    Co_run    = res_erun$Co
    sigma_run = res_erun$Sigma

    G_run = NA
    if (irf=="PTdecomp") {
      crk = res$crk

      alpha_run  =  (as.matrix(res_erun$tst$estimation$coefficients[1:crk,]))
      if (dim(alpha_run)[1]<dim(alpha_run)[2]) {
        alpha_run = t(alpha_run)
      }
      G_run     = t(cbind(alpha_perpf(alpha_run),alpha_run))
    }

    if ( irf =="ABSVAR") {
      x0 = stats::rnorm(n*(n+1)/2)
      ABD <- ABSVAR(x0,A0,B0,Sigma=sigma_run)
      A0_run = ABD[[1]]
      B0_run = ABD[[2]]
    }   else {
      A0_run = NA
      B0_run = NA
    }
    if (irf == "irfX") { Xshk = t(Co_run[,1+Xshks]) }
    responseR[, , , i] <- irf_B_sigma(B = B_run, sigma = sigma_run, nstep, comb, irf, G=G_run,A0=A0_run,B0=B0_run,smat=NA,Xshk=Xshk)
  }
  for (tt in 1:(nstep)) {
    for (i in 1:neq) {
      for (j in 1:nvar) {
        response[i, j, tt, -1] = stats::quantile(responseR[i,j, tt, ], conf)

      }
    }
  }
  #response[, , , 1] = responseR[,,,5]
  return(response)
}


#' Estimation of CCIVAR
#'
#' This function estimates parameters of a specified conditional cointegrated VAR based on provided data.
#'
#' @param res a CCIVAR object that can be an output of CCIVARData containing at least n1,n2, Y, X, and crk.
#'
#' @return a CCIVAR object with estimated parameter values, AIC, BIC, conditional VECM in a regression format.
#'
#' @examples
#' T = 100
#' res_d <- CCIVARData(n1=4,n2=3,crk=3,p=3,T=T,type="const")
#' res_e <- CCIVARest(res=res_d)
#' res_e$Summary
#'
#' @export
CCIVARest <- function (res)
{
  n = res$n1
  p = res$p
  crk = res$crk
  Y = as.matrix(res$Y)
  if (is.null(colnames(Y)))
    colnames(Y) = paste0(rep("Y", ncol(Y)), c(1:ncol(Y)))
  type = res$type
  T = dim(Y)[1]
  X = res$X
  if (is.null(colnames(X)))
    colnames(X) = paste0(rep("X", ncol(X)), c(1:ncol(X)))
  if (!anyNA(X)) if (is.null(colnames(X))) colnames(X) = paste0(rep("exog", ncol(X)), c(1:ncol(X)))
  if (type == "none" | type == "exog0")     Model = "I"
  if (type == "const" | type == "exog1")    Model = "III"
  if (type == "trend")                      Model = "IV"
  P = matrix(0, 2, 3)
  P[, 1] = p
  Co = res$Co * 0
  #      cvecm(y=Y,x=X,model="II",type = "eigen", crk=2, p = 3, q = 0.95)
  tst <- cvecm(y=Y,x=X,model=Model,type = "eigen", crk=crk, p = p, q = 0.95)
  param = tst[[2]][[1]]
  if (type=="trend") beta = tst$beta[-1,] else beta = tst$beta
  CIVAREST = CVECM2CVAR(param = param, beta = beta, p = c(crk,p - 1, p-1),N2 = ncol(X))
  B = CIVAREST[[1]]
  Bx = CIVAREST[[2]]
  C = CIVAREST[[3]]
  alpha <- t(param[1:crk,]); dim(alpha) = c(n,crk)
  #if (type=="trend")   C1 = rowSums(alpha)  else C1 = NA
  if (type=="trend")   C1 = alpha%*%t(as.matrix(tst$beta))[,1]  else C1 = NA
  LREG = tst[[2]]
  sigma = t(LREG$residuals) %*% (LREG$residuals)/(T)
  resid = Y * 0
  resid[(p + 1):T, ] = LREG$residuals
  if (type == "const" | type == "exog1")
    Co = C
  if (type == "exog0") {
    Co[, 2:(1 + dim(CIVAREST[[3]])[2])] = CIVAREST[[3]]
  }
  res$By <- B
  res$Bx <- Bx
  res$Co <- Co
  res$C1 <- C1
  res$Sigma <- sigma
  res$resid <- resid
  LH = -(T * n/2) * log(2 * pi) - (T * n/2) + ((T)/2) * log(det(solve(sigma)))
  AIC = 2 * n * (dim(tst[[2]][[1]])[1]) + n * (n + 1) - 2 *
    LH
  BIC = log(T) * (n * (dim(tst[[2]][[1]])[1]) + n * (n + 1)/2) -
    2 * LH
  LH_N = 2 * n * (dim(tst[[2]][[1]])[1]) + n * (n + 1)
  if (is.null(colnames(Y)))
    colnames(Y) = sprintf("Y%s", 1:n)
  estimation_result = summaryCIVAR(lm1 = LREG, sname = "Z2")
  Summary = list(estimation_result, LH, AIC, BIC, LH_N, tst$erg)
  names(Summary) = c("Estimation_Results", "LH_function_Value","AIC", "BIC", "LH_N", "Johansen_Test")
  res$tst = tst
  res$Summary = Summary
  return(res)
}




#' Transformation of the estimated parameters of a conditional VECM into the parameters of the corresponding conditional VAR
#'
#' @param param estimated parameters of a conditional VECM
#' @param beta the estimated cointegration vecters
#' @param p a vector specifying different types of conditional VECM
#' @param s indicator variable of different regimes
#' @param N2 dimension of the conditioning variables
#'
#' @return A list of (A, B, C) with B the conditional VAR parameter matrix, C the parameter of the deterministic components, and A the parameter matrix of foreign variables for CIGVAR.
#' @export
#'
CVECM2CVAR <- function (param, beta, p = c(1, 2, 2, 2, 2), s = NA, N2)
{
  m = dim(param)[2]
  VECB = t(param)
  if (anyNA(s)) {                                #### one regime
    if (!p[3] == 0) {
      B = (1:(m * m * (p[2] + 1))) * 0
      dim(B) = c(m, m, (p[2] + 1))
      A = (1:(m * N2 * (p[3] + 1))) * 0
      dim(A) = c(m, N2, (p[3] + 1))
      AB = VECB[, 1:p[1]] %*% t(beta)
      B[, , 1] = AB[, 1:m]
      A[, , 1] = AB[, (m + 1):(m+N2)]
      for (i in 2:(p[2] + 1)) B[, , i] = VECB[, p[1] +  ((i - 2) * m + 1):((i - 2) * m + m)]
      for (i in 2:(p[3] + 1)) A[, , i] = VECB[, (p[1] +  p[2] * m) + ((i - 2) * N2 + 1):((i - 2) * N2 + N2)]
      B <- CIB3B(B)
      A <- CIA2A(A)
    }
    else {
      B = (1:(m * m * (p[2] + 1))) * 0
      dim(B) = c(m, m, (p[2] + 1))
      AB = VECB[, 1:p[1]] %*% t(beta)
      B[, , 1] = AB[, 1:m]
      if (p[2] > 0)
        for (i in 2:(p[2] + 1)) B[, , i] = VECB[, p[1] +
                                                  ((i - 2) * m + 1):((i - 2) * m + m)]
      B <- CIB3B(CIB = B)
      A = NA
    }
    if (dim(param)[1] > p[1] + (p[2]*m + p[3]*N2) ) {
      C = as.matrix(t(param)[, (p[1] + (p[2]*m + p[3]*N2) + 1):dim(param)[1]])
    }      else C = NA
  }
  else {                                #### MRCIGVAR (crk, lag_domestic_1,lag_foreign_1,lag_domestic_2,lag_foreign_2)
    if (length(p) == 5) {
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
        C[, 1, 1] = as.matrix(t(param)[, p[1] + sum(p[-1]) *
                                         m + 1])
        C[, 1, 2] = as.matrix(t(param)[, p[1] + sum(p[-1]) *
                                         m + 2])
      }
      else {
        C = NA
      }
    }
    if (length(p) == 3) {                       ### MRCIVAR  crk, lag of regime 1,  lag of regime two
      P = max(p[-1])
      B = (1:(m * m * (P + 1) * 2)) * 0
      dim(B) = c(m, m, (P + 1), 2)
      AB = VECB[, 1:p[1]] %*% t(beta)
      B[, , 1, 1] = AB[, 1:m]
      B[, , 1, 2] = AB[, 1:m]
      for (i in 2:(p[2] + 1)) B[, , i, 1] = VECB[, p[1] +
                                                   ((i - 2) * m + 1):((i - 2) * m + m)]
      for (i in 2:(p[3] + 1)) B[, , i, 2] = VECB[, (p[1] +
                                                      p[2] * m) + ((i - 2) * m + 1):((i - 2) * m +
                                                                                       m)]
      B[, , 1:(p[2] + 1), 1] = CIB3B(B[, , 1:(p[2] + 1),
                                       1])
      B[, , 1:(p[3] + 1), 2] = CIB3B(B[, , 1:(p[3] + 1),
                                       2])
      A = NA
      if (dim(param)[1] > p[1] + (sum(p[-1])) * m) {
        C = (1:(m * 2)) * 0
        dim(C) = c(m, 1, 2)
        C[, 1, 1] = as.matrix(t(param)[, p[1] + sum(p[-1]) *
                                         m + 1])
        C[, 1, 2] = as.matrix(t(param)[, p[1] + sum(p[-1]) *
                                         m + 2])
      }
      else {
        C = NA
      }
    }
  }
  return(list(B, A, C))
}







#' Re-indexing the order of the endogenous and exogenous variables
#'
#' @param N1 dimension of the endogenous variables
#' @param N2 dimension of the exogenous variables
#' @param p lag length
#'
#' @return a new index
#' @export
#'
Redex = function(N1,N2,p) {
  M     = N1 + N2
  index = c(1:((p-1)*M))
  Rindex = index*0
  for (i in 1:(p-1)) {
    Rindex[((i-1)*N1+1):((i-1)*N1+N1)]          =  index[((i-1)*M+1):((i-1)*M+N1)]
    Rindex[N1*(p-1)+((i-1)*N2+1):((i-1)*N2+N2)] =  index[((i-1)*M+1+N1):((i-1)*M+N1+N2)]
  }
  Rindex
}







#' Conditional cointegrated VAR
#'
#' The function generates data from a conditional cointegrated VAR(p)
#'
#' @param n1 dimension of the conditional cointegrated process
#' @param n2 dimension of the conditioning variables
#' @param crk the cointegration rank
#' @param p lag
#' @param T number of observations
#' @param type types of the deterministic component. "none" and "const" are two options.
#' @param Bc Coefficient matrix of the joint cointegration process
#'
#' @return a CCIVAR object which is a list of (n1,n2,p,type,r_np,By,Bx,Cy,Sigma,Y,X,resid,U,check,crk,Bc,Cc)
#' @export
#'
#' @examples
#'
#' T = 100
#' res_d <- CCIVARData(n1=4,n2=3,crk=3,p=3,T=T,type="const")
#' res_e <- CCIVARest(res=res_d)
#' res_e$Summary
#'
#'
CCIVARData<- function(n1,n2,crk,p,T,type,Bc=NA) {
  n = n1+n2
  if (anyNA(Bc)) {
    res_d = CIVARData(n=n,p=p,T=T,type=type,crk=crk)
    B     = res_d$B
    BAB = B2CIB(B)
    alpha <- Re(BAB[[2]]); alpha[(n1+1):(n1+n2),]<- 0
    beta  <- Re(BAB[[3]])
    CIB   <- BAB[[1]]
    CIB[,,1] = alpha%*%t(beta)
    Bc = CIB3B(CIB)
  }
  red_d = CIVARData(n=n,p=p,T=T,B=Bc,Co = rep(1,n),type=type,crk=crk)
  r_np <- red_d$r_np
  Bd <- red_d$B[1:n1,,]
  By = Bd[,1:n1,]
  Bx = Bd[,(n1+1):(n1+n2),]
  Cy = red_d$Co[1:n1]
  Sigma = red_d$Sigma[1:n1,1:n1]
  Y = red_d$Y[,1:n1]
  X = red_d$Y[,(n1+1):(n1+n2)]
  resid = red_d$resid[,1:n1]
  check = max(abs(red_d$Y))
  U     = red_d$U[,1:n1]
  Bc    <- red_d$B
  Cc    <- red_d$Co
  result <- list(n1,n2,p,type,r_np,By,Bx,Cy,Sigma,Y,X,resid,U,check,crk,Bc,Cc)
  names(result) = c("n1","n2", "p", "type", "r_np", "By", "Bx", "Cy", "Sigma", "Y", "X", "resid","U","check","crk","Bc","Cc")
  return(result)
}









#' Estimation of a conditional vector error correction process.
#'
#' This function estimates the unknown parameters of a multi regime conditional VECM based on provided data.
#'
#'
#' @param y data matrix of the endogenous variables
#' @param x data matrix of the exogenous variables
#' @param model type of the specification of the deterministic components. "none" and "const" are two options.
#' @param type =c("eigen", "trace")
#' @param crk cointegration rank
#' @param p lag length in level
#' @param q significance level
#'
#' @return a list of estimated parameters and test statistics
#' @export
cvecm <- function (y,x,model = c("I","II","III","IV","V"),type = c("eigen", "trace"), crk=2, p = 1, q = 0.95)
{
  y <- as.matrix(y)
  if (q != 0.9 && q != 0.95 && q != 0.99) {
    return("please correct significance level")
  }

  p <- as.integer(p)
  N1 <- ncol(as.matrix(y))
  N  <- ncol(as.matrix(y))
  if (N1<crk) {
    return("y's dimension must be larger than crk")
  }

  if  ( length(x) == 1 ) {  z <- y }
  else {
    z <- cbind(y,x)
  }

  M1 <- ncol(as.matrix(z))
  T <- nrow(as.matrix(z))
  #Z <- embed(diff(z), p)
  Z <- Embed(diff(z), p)

  Y0 <- Z[, c(1:N1)]                # \Delta Y
  #X0 <- Z[,c((N1+1):M1)]           # \Delta X
  Z1 <- z[-T, ][p:(T - 1), ]        #      Z_1
  if ((p==1)&(length(x)==1)) {
    Z2 = NULL
    NOLAG = 1
  }
  else {
    Z2 <- Z[, c((M1+1):(M1*p))]               # \Delta Z_
    Z2 <- Z2[, Redex(N1=N1,N2=ncol(x),p)]       # \Delta Z_
  }
  T1 <- nrow(as.matrix(Y0))
  #lT = (1:T1)/(1:T1)
  #Trend <- matrix(1:T1, T1, 1)
  lT = matrix(1,T1,1); colnames(lT) <- "Const"
  Trend <- matrix(1, T1, 1); colnames(Trend) <- "Trend"

  if ( model == "I" ) {
    Y0 = Y0
    Z1 = Z1
    Z2 = Z2
  }
  if ( model == "II" ) {
    Y0 = Y0
    Z1 = cbind(lT,Z1)
    Z2 = Z2
  }
  if ( model == "III" ) {
    Y0 = Y0
    Z1 = Z1
    Z2 = cbind(Z2,lT)
  }
  if ( model == "IV" ) {
    Y0 = Y0
    Z1 = cbind(Trend,Z1)
    Z2 = cbind(Z2,lT)
  }
  if ( model == "V" ) {
    Y0 = Y0
    Z1 = cbind(Z1)
    Z2 = cbind(Z2,lT,Trend)
  }

  if (length(Z2)==0) {
    R0 <- Y0
    R1 <- Z1
  }
  else {

    M00 <- crossprod(Y0)/T1
    M11 <- crossprod(Z1)/T1
    M22 <- crossprod(Z2)/T1
    M01 <- crossprod(Y0, Z1)/T1
    M02 <- crossprod(Y0, Z2)/T1
    M10 <- crossprod(Z1, Y0)/T1
    M20 <- crossprod(Z2, Y0)/T1
    M12 <- crossprod(Z1, Z2)/T1
    M21 <- crossprod(Z2, Z1)/T1
    M22inv <- solve(M22)
    R0 <- Y0 - t(M02 %*% M22inv %*% t(Z2))
    R1 <- Z1 - t(M12 %*% M22inv %*% t(Z2))
  }
  S00 <- crossprod(R0)/T1
  S01 <- crossprod(R0, R1)/T1
  S10 <- crossprod(R1, R0)/T1
  S11 <- crossprod(R1)/T1
  Ctemp <- chol(S11, pivot = TRUE)
  pivot <- attr(Ctemp, "pivot")
  oo <- order(pivot)
  C <- t(Ctemp[, oo])
  Cinv <- solve(C)
  S00inv <- solve(S00)
  valeigen <- eigen(Cinv %*% S10 %*% S00inv %*% S01 %*% t(Cinv))

  #valeigen <- eigen(S11inv %*% S10 %*% S00inv %*% S01 )
  lambda <- Re(valeigen$values)
  e      <- Re(valeigen$vector)

  V <- t(Cinv) %*% e
  Vorg <- V
  V <- sapply(1:M1, function(j) V[, j]/V[1, j])
  #W <- S01 %*% V %*% solve(t(V) %*% S11 %*% V)
  PI <- S01 %*% solve(S11)
  #DELTA <- S00 - S01 %*% V %*% solve(t(V) %*% S11 %*% V) %*% t(V) %*% S10
  if (length(Z2)==0) {
    GAMMA = NULL
  }
  else {
    GAMMA <- M02 %*% M22inv - PI %*% M12 %*% M22inv
  }
  beta   <- V[,1:crk]

  if ((crk>0)&(length(Z2)>0)) {
    CI = Z1%*%beta
    VECM<-stats::lm(Y0~0+CI+Z2)
  }
  if ((crk>0)&(length(Z2)==0)) {
    CI = Z1%*%beta
    VECM<-stats::lm(Y0~0+CI)
  }
  if ((crk == 0)&(length(Z2)>0)) {
    CI = 0
    VECM<-stats::lm(Y0~0+Z2)
  }
  if ((crk == 0)&(length(Z2)==0)) {
    CI = 0
    VECM<-NULL
  }
  E = -T1 * log(1 - lambda)
  E = E[1:N]
  ##### Estimation of VECM

  resultsvecm <- summary(VECM)

  if (model == "I")   {    Tab <- Tab1  }
  if (model == "II")  {    Tab <- Tab2  }
  if (model == "III") {    Tab <- Tab3  }
  if (model == "IV")  {    Tab <- Tab4  }
  if (model == "V")   {    Tab <- Tab5  }


  b = c(1:12)
  for (i in 1:12)   { b[i] = Tab[2*(i-1)+1,2+M1-N1] }
  a = c(1:12)
  for (i in 1:12)   { a[i] = Tab[2*i,2+M1-N1] }

  if (type == "eigen") {
    critical_vals = b
    M = matrix(0, N, 1)
    j = 1
    rank = 0
    while (j <= N && E[j] > critical_vals[N + 1 - j]) {
      M[j, ] = M[j, ] + 1
      j = j + 1
      rank = rank + 1
    }
    erg <- cbind(E, critical_vals[N:1])
    colnames(erg) <- c("teststatistic", "critical_value")

    if ( N > 1 ) {
      rownames(erg) <- c("crk <= 0 |", paste("crk <= ", 1:(N - 1),  " |", sep = ""))
    }
    if ( N == 1 ) {
      rownames(erg) <- c("crk <= 0 |" )
    }

    coint_rank <- paste("Johansen-Test (with maximum-eigenvalue-teststatistic) indicates",
                        rank, "cointegrating equation(s) at the", 1 - q,
                        "level")
    #print(coint_rank)
    result <- new.env()
    result$erg = erg
    result$estimation = VECM
    result$lambda     = E
    result$z          = z
    result$Z2         = Z2
    result$beta       = beta
    result$PI         = PI
    result$GAMMA      = GAMMA
    result$model      = model
    rr<-as.list(result)
    return(rr)
  }
  else {
    type = "trace"
    critical_vals = a
    stat = matrix(0, N, 1)
    for (i in 1:N) {
      sum = 0
      for (j in i:N) {
        sum = sum + E[j]
      }
      stat[i] = sum
    }
    M = matrix(0, N, 1)
    j = 1
    rank = 0
    while (stat[j] > critical_vals[N + 1 - j] && j <= N) {
      M[j, ] = M[j, ] + 1
      j = j + 1
      rank = rank + 1
    }
    erg <- cbind(stat, critical_vals[N:1])
    colnames(erg) <- c("teststatistic", "critical_value")

    if ( N > 1 ) {
      rownames(erg) <- c("crk <= 0 |", paste("crk <= ", 1:(N - 1),  " |", sep = ""))
    }
    if ( N == 1 ) {
      rownames(erg) <- c("crk <= 0 |")
    }

    coint_rank <- paste("Johansen-Test (with trace-teststatistic) indicates",
                        rank, "cointegrating equation(s) at the", 1 - q,
                        "level")
    #print(coint_rank)
    result <- new.env()
    result$erg = erg
    result$estimation = VECM
    result$lambda     = E
    result$z          = z
    result$Z2         = Z2
    result$beta       = beta
    result$PI         = PI
    result$GAMMA      = GAMMA
    result$model      = model
    rr<-as.list(result)
    return(rr)
  }
}

#' Mixed VECM with I(0) and I(1) variables
#'
#' This function generates and estimates a mixed VECM with I(0) and I(1) variables
#'
#' @param n	: dimension of the mixed CIVAR
#' @param p	: lag of the mixed CIVAR
#' @param T : number of observations or length of the generated data
#' @param r : number of unit roots in the joint CIVAR
#' @param k  : number of I(0) variables
#' @param type : types of the deterministic components in the conditional CIVAR
#' @param Bo : coefficient matrix of the mixed CIVAR
#' @param Y : data of the mixed CIVAR
#' @param X : data of the exogenous variables
#' @param D : transformation matrix to mix the I(0) and I(1) components
#' @param Go : pre-loaded selection matrix
#' @param B : coefficient matrix of the mixed CIVAR
#' @param Sigma : covariance matrix of the residuals
#'
#' @return  a list contains estimation and test results of the mixed VECM
#' @examples
#'
#' #RR <- MIxCIVARData(n=9,p=2,T=209,r=5,k=2,type="const",Bo=NA,Y=NA,D=NA)
#' ## DGP of I(1) and I(0) mixed VECM
#' #plot(ts(RR$Y))
#'
#' #### testing the mixed VECM via testing the restrictions on beta
#' #res_d <- CIVARData(n=9,p=2,T=209,type="const",crk=4)
#' #res_e = CIVARest(res=res_d)
#' #res_e$Summary
#' #n = 9; crk = 4; k = 2; r = 5
#' #CC  <- c(8,9,17,18)
#' #GG  <- c(19:25,28:34)
#' #G = diag(n*crk); psi=matrix(1,n*crk,1)
#' #### this implies there is no restrictions on the adjustment coefficients alpha
#' #H = diag(n*crk);              H2 = H[,-c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n),CC,GG)]
#' #### only normalization
#' #h = matrix(0,n*crk,1);         h[c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n)),1] <- 1
#' #phi = matrix(1,ncol(H2),1)
#' ### check consistency of the restrictions    G%*%psi; H2%*%phi + h
#' #
#' #G%*%psi; H2%*%phi + h
#' #
#' #ABtest = CIVARTest(res=res_d,H=H2,h=h,phi=phi,G=G,Dxflag=0)
#' #ABtest$betar
#' #ABtest$alphar
#' #ABtest$VECMR$coefficients
#' #ABtest$LR
#' #1-pchisq(ABtest$LR,14)   ### The Ho of last two are I(0) can be rejected
#' #RR$GABtest$LR
#' #
#' @export
#'
MIxCIVARData = function(n,p,T,r,k,type,Bo=NA,Y=NA,X=NA,D=NA,Go=NA,B=NA,Sigma = NA) {

  if (anyNA(Bo)) {
    repeat {
      res_d  = VARData(n=n,p=p,T=T,type=type)
      Bo <- res_d$B
      Bo[,1:r,p] = 0
      if (max(Mod(STAT(Bo))) < 0.9) break
    }
  }
  res_d  = VARData(n=n,p=p,T=T,B=Bo,type=type)
  U = res_d$U
  Z = as.matrix(res_d$Y[,1:r]); for (i in 2:T) Z[i,] = Z[i-1,]+Z[i,]
  X <- res_d$Y[,(r+1):n]
  #plot(ts(cbind(Z,X)))

  if (anyNA(D)) {
    D = matrix(stats::rnorm(n*n),n,n)
    #solve(D)
    #D[1:(n-k),(n-k+1):n] <-0
    D[(n-k+1):n,1:(n-k)] <-0
    #D[(n-k+1):n,1:r] <-0
  }
  YS = cbind(Z,X)%*%t(D)
  #plot(ts(YS))
  crk <- n-r
  ########### restricted CIVAR

  B1 <- Bo
  B1[1:r,1:r,1]    =  Bo[1:r,1:r,1]  +  diag(r)
  B1[1:n,1:r,p]    = -Bo[1:n,1:r,p-1]
  if (p>2) for (i in 2:(p-1)) B1[1:n,1:r,i]  = Bo[1:n,1:r,i] - Bo[1:n,1:r,i-1]
  #STAT(B1)

  if (anyNA(B)) {
    B2 <- B1*0
    for (i in 1:p) B2[,,i] <- D%*%B1[,,i]%*%solve(D)
  } else B2 <- B

  #STAT(B2)

  #Co = seq(1,n,1)/seq(1,n,1)
  #U[1:2,] <- 0
  #res_d  = VARData(n=n,p=p,T=T,B=Bo,Co=Co*0,type="const",U=U,Yo=U[1:p,]*0)
  #ZZ <- res_d$Y[,1:r]; Z = ZZ*0; for (i in 3:T) Z[i,] = Z[i-1,]+ZZ[i,];X <- res_d$Y[,(r+1):n]
  #res_dd = VARData(n=n,p=p,T=T,B=B1,Co=Co*0,type="const",U=U,Yo=Yo)
  #res_ddd= VARData(n=n,p=p,T=T,B=B2,Co=Co*0,type="const",U=U%*%t(D),Yo=U[1:p,]*0)
  #Y =  cbind(Z,X)%*%t(D)
  #res_ddd$Y - Y

  ##############
  res_d <- CIVARData(n=n,p=p,T=T,B=B2,type=type,crk=crk,Sigma=Sigma)
  #plot(ts(res_d$Y))
  if (!anyNA(Y)) res_d$Y <-Y
  #res_d$Y <- Y
  #plot(ts(res_d$Y))
  res_e = CIVARest(res=res_d)

  G = diag(n*crk); psi=matrix(1,n*crk,1)          #### this implies there is no restrictions on the adjustment coefficient alpha
  H = diag(n*crk);              H2 = H[,-c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n))]   #### only normalization
  h = matrix(0,n*crk,1);        h[c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n)),1] <- 1 ;
  phi = matrix(1,(n-1)*crk,1)

  ### second set of normalization beta' = (I,beta_2') not done
  G = diag(n*crk); psi=matrix(1,n*crk,1)          #### this implies there is no restrictions on the adjustment coefficient alpha
  H = diag(n*crk);              H2 = H[,-c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n))]   #### only normalization
  h = matrix(0,n*crk,1);        h[c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n)),1] <- 1 ;
  phi = matrix(1,(n-1)*crk,1)

  ### check restrictions: G%*%psi  ; H2%*%phi+h
  ### c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n))

  ABtest = CIVARTest(res=res_d,H=H2,h=h,phi=phi,G=G,Dxflag=0)


  ABtest$betar
  ABtest$alphar
  ABtest$VECMR$coefficients
  ABtest$LR
  1-stats::pchisq(ABtest$LR,1)   ### the fourth is the used restrictions


  ##################
  CCo = seq(n-k,(n-r-k)*n,n)
  CC  <-CCo+1
  if (k > 1) for (i in 2:k) CC = cbind(CC,CCo+i)
  CC = as.vector(t(CC))
  GGo = 1:(n-k)
  GG <- (n-r-k)*n + GGo
  if (k>1) for (i in 2:k) GG = cbind(GG,GGo+(n-r-k+i-1)*n)
  GG <- as.vector(t(GG))
  GG0 = 1:r
  GG1 <- (n-r-k)*n + GG0
  #GG1 <- c(2,3)
  if (k>1) for (i in 2:k) GG1 = cbind(GG1,GG0+(n-r-k+i-1)*n)
  GG1 <- as.vector(GG1)
  G = diag(n*crk); psi=matrix(1,n*crk,1)          #### this implies there is no restrictions on the adjustment coefficient alpha
  H = diag(n*crk);              H2 = H[,-c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n),GG1)]   #### normalization + I(0) B31=0
  h = matrix(0,n*crk,1);        h[c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n)),1] <- 1 ;
  phi = matrix(1,(n-1)*crk-k*r,1)

  ### check restrictions: G%*%psi  ; H2%*%phi+h
  ### c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n))

  AABtest = CIVARTest(res=res_d,H=H2,h=h,phi=phi,G=G,Dxflag=0)

  AABtest$betar
  AABtest$alphar
  AABtest$VECMR$coefficients
  AABtest$LR
  1-stats::pchisq(AABtest$LR,k*(n-r))   ### the fourth is the used restrictions
  #plot(ts(res_d$Y))


  #################


  CCo = seq(n-k,(n-r-k)*n,n)
  CC  <-CCo+1
  if (k > 1) for (i in 2:k) CC = cbind(CC,CCo+i)
  CC = as.vector(t(CC))
  GGo = 1:(n-k)
  GG <- (n-r-k)*n + GGo
  if (k>1) for (i in 2:k) GG = cbind(GG,GGo+(n-r-k+i-1)*n)
  GG <- as.vector(t(GG))
  ## GG <- c(14,15,16,17,18)
  G = diag(n*crk); psi=matrix(1,n*crk,1)          #### this implies there is no restrictions on the adjustment coefficient alpha
  H = diag(n*crk);              H2 = H[,-c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n),GG)]   #### normalization + I(0) B31=0,B32=0
  h = matrix(0,n*crk,1);        h[c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n)),1] <- 1 ;
  phi = matrix(1,(n-1)*crk-k*(n-k),1)
  ### H2 = H[,-c(1,3,8,5,6,7)]; h[3]=-1; phi = matrix(1,2,1)
  ### check restrictions: G%*%psi  ; H2%*%phi+h
  ### c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n))

  ABtest = CIVARTest(res=res_d,H=H2,h=h,phi=phi,G=G,Dxflag=0)

  ABtest$betar
  ABtest$alphar
  ABtest$VECMR$coefficients
  ABtest$LR
  1-stats::pchisq(ABtest$LR,k*r)   ### the fourth is the used restrictions the degree of freedom is n*h-h*h







  #########################

  G = diag(n*crk); psi=matrix(1,n*crk,1)          #### this implies there is no restrictions on the adjustment coefficient alpha
  H = diag(n*crk);              H2 = H[,-c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n),CC,GG)]   #### normalization I(0) B31=0, B32=0, B23=0
  h = matrix(0,n*crk,1);        h[c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n)),1] <- 1 ;
  phi = matrix(1,(n-1)*crk-k*(n-r-k)-k*(n-k),1)

  ### H2 = H[,-c(1,3,4,8,5,6,7)]; h[3]=-1; phi = matrix(1,1,1)

  ### check restrictions: G%*%psi  ; H2%*%phi+h
  ### c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n))
  nG = 0
  if (!anyNA(Go)) { G = Go; nG = dim(Go)[1]-dim(Go)[2] }
  GABtest = CIVARTest(res=res_d,H=H2,h=h,phi=phi,G=G,Dxflag=0)

  GABtest$betar
  GABtest$alphar
  GABtest$VECMR$coefficients
  GABtest$LR
  1-stats::pchisq(GABtest$LR,k*r+nG)   ### the fourth is the used restrictions
  # collect restricted estimates
  BB3 <- res_d$B*0
  #BB3[,,1] <- GABtest$alphar%*%t(GABtest$betar)
  #AA  <- t(GABtest$LSKOEFR[1:((p-1)*n),])
  #dim(AA) <- c(n,n,(p-1))
  #BB3[,,2:p] <- AA
  #B3 <-CIB3B(BB3)
  #res_e$B <- B3
  #res_e$Sigma <- t(GABtest$VECMR$residuals)%*%GABtest$VECMR$residuals/T
  #for (i in 1:n) {
  #    res_e$tst$VECMS[[i]]$coefficients[(n-r+1):nrow(res_e$tst$VECMS[[i]]$coefficients),] <- GABtest$VECMRS[[i]]$coefficients;
  #    res_e$tst$VECMS[[i]]$coefficients[1:(n-r),1] <- GABtest$alpha[i,];
  #}
  BB3[,,1] <- GABtest$alphar%*%t(GABtest$betar)
  AA  <- t(GABtest$LSKOEFR[1:((p-1)*n),])
  dim(AA) <- c(n,n,(p-1))
  BB3[,,2:p] <- AA
  B3 <-CIB3B(BB3)
  res_e$B <- B3
  ##### transforming back to B

  #####
  res_e$Sigma <- t(GABtest$VECMR$residuals)%*%GABtest$VECMR$residuals/T
  for (i in 1:n) {
    res_e$tst$VECMS[[i]]$coefficients[(n-r+1):nrow(res_e$tst$VECMS[[i]]$coefficients),] <- GABtest$VECMRS[[i]]$coefficients;
    res_e$tst$VECMS[[i]]$coefficients[1:(n-r),1] <- GABtest$alpha[i,];
  }
  RR<-list(D,Bo,res_d$Y,AABtest, ABtest, GABtest,res_e,res_d,1-stats::pchisq(AABtest$LR,1),1-stats::pchisq(ABtest$LR,k*r), 1-stats::pchisq(GABtest$LR,k*r),k)
  names(RR) = c("D","Bo","Y","AABtest", "ABtest", "GABtest","res_e","res_d","p-value_A","p-value_E", "p-value_W","k")
  return(RR)
}

#' This function calculates the orthogonal complementary space
#'
#' @param alpha : a matrix with independent columns
#'
#' @return  the orthogonal complementary matrix
#' @export
#'
alpha_perpf = function(alpha=(alpha)) {

  alpha <- as.matrix(alpha)
  n = dim(alpha)[1]
  nx = dim(alpha)[2]
  nz = n-nx
  c = diag(n)[,(nz+1):(nz+nx)]
  c_perp = diag(n)[,1:nz]
  alpha_perp = (diag(n)-c%*%solve(t(alpha)%*%c)%*%t(alpha))%*%c_perp
  t(alpha)%*%alpha_perp
  return(alpha_perp)
}




#' This function solves AB_SVAR from the reduced form and return the A, B matrices and the sum of squared errors
#'
#'
#' @param x 	  : difference between the reduced form and the AB form
#' @param A0	  : A matrix in an AB_SVAR model
#' @param B0    : B matrix in an AB_SVAR model
#' @param Sigma : Covariance matrix of the reduced form
#'
#' @return  difference
#' @export
#'
xABF <-  function(x,A0,B0,Sigma) {
  ### this function solves A and B from the reduced form Sigma for a AB-SVAR
  n = dim(A0)[1]
  if (sum((A0==0)|(A0==1))+sum((B0==0)|(B0==1))!= 2*n*n-n*(n+1)/2) return("not exactly identified")
  if (length(x)!=n*(n+1)/2)  return("dimension of x is incorrect")
  A = A0
  B = B0
  VA = as.vector(A)
  VB = as.vector(B)

  VA[!((VA==0)|(VA==1))] = x[1:length(VA[!((VA==0)|(VA==1))])]

  VB[!((VB==0)|(VB==1))] = x[(length(VA[!((VA==0)|(VA==1))])+1):(length(VA[!((VA==0)|(VA==1))])+length(VB[!((VB==0)|(VB==1))]))]

  A = matrix(VA,n,n)
  B = matrix(VB,n,n)
  return(sum((A%*%Sigma%*%t(A)-B%*%t(B))^2))
}



#' This function solves AB_SVAR from the reduced form and return the A, B matrices and the sum of squared errors
#'
#'
#' @param x0	  : difference between the reduced form and the AB form
#' @param A0	  : A matrix in an AB_SVAR model
#' @param B0    : B matrix in an AB_SVAR model
#' @param Sigma : Covariance matrix of the reduced form
#'
#' @return  a list contains A B and the difference
#' @export
#'
ABSVAR = function(x0,A0,B0,Sigma) {
  ## This function solve AB_SVAR from the reduced form and return the A, B matrices and the sum of squared errors
  n =dim(A0)[1]
  ABdiff <- stats::nlm(xABF,x0,A0,B0,Sigma,gradtol=0.000000000000001 )
  x = ABdiff$estimate
  A = A0
  B = B0
  VA = as.vector(A)
  VB = as.vector(B)
  VA[!((VA==0)|(VA==1))] = x[1:length(VA[!((VA==0)|(VA==1))])]
  VB[!((VB==0)|(VB==1))] = (x[(length(VA[!((VA==0)|(VA==1))])+1):(length(VA[!((VA==0)|(VA==1))])+length(VB[!((VB==0)|(VB==1))]))])
  A = matrix(VA,n,n)
  B = abs(matrix(VB,n,n))
  diff <-sum((A%*%Sigma%*%t(A)-B%*%t(B))^2)
  ## identification
  diag(solve(A)%*%B)
  return(list(A,B,diff))
}



