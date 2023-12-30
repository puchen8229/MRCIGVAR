#' Data generating process of GVAR(m,n,p)
#'
#' This function generates data from a stationary global vector auto regressive process and return a GVAR(m,n,p) object containing data and parameters used in the GVAR(m,n,p) process.
#'
#' @param m     : number of variables
#' @param n     : number of countries/units
#' @param p     : an (n x 3) matrix in which each row contains the lag length of the domestic variables, the lag length of the foreign variables, and the number of exogenous variables.
#' @param T     : number of observations
#'
#'                (m,n,p,T) are parameters which must be provided.
#' @param W     : an (n x n) weighting matrix. w_ij is the weight of foreign country j in the foreign variables of i-th country diag(W)=0
#' @param r_npo : an (m, p, n) array collecting the roots of the characteristic polynomials in Lags for each of the m dimensional  variable across n countries.
#' @param Ao    : an (m, m, p, n) array collecting the off-diagonal block of coefficients that are coefficients of the foreign variables.
#' @param Bo    : an (m, m, p, n) array collecting the coefficients of the domestic variables.
#' @param Co    : an (m, k+1, n) array collecting the coefficients of the deterministic components of the n countries.
#' @param Uo    : an (T x mn) matrix of the temporally independent innovations
#' @param Sigmao : an (mn x mn) matrix of the covariance matrix of the GVAR(m,n,p)
#'                (W,r_npo,Ao,Bo,Uo,Sigmao) if not provided, they will be generated randomly.
#'
#' @param type  : types of deterministic components: "const", "none", "exog0", and "exog1" are the options.
#' @param X     : (T x k x n ) array of exogenous variables.
#' @param mu    : if type = "const" mu has the same dimension as Co. It contains the means of the time series in the system.
#' @return a GVAR object that is a list containing the generated data, the used parameters, and the inputted of exogenous variables.
#'
#' @export
#'
#' @examples
#'
#' n = 5
#' p = (1:15)*0; dim(p) = c(5,3)
#' p[,1] = 2; p[,2]=1;
#' res_d = GVARData(m=2,n=5,p=p,T=100,type="const")
#' max(res_d$Y)
#' dim(res_d$Y)
#' res_e = GVARest(res = res_d)
#' res_e$Summary
#'
#' X1 = matrix(1,200,1)
#' X2 = matrix(rnorm(200),200,1)
#' X3 = matrix(rnorm(200),200,1)
#' X4 = matrix(rnorm(200),200,1)
#' X  = cbind(X1,X2,X3,X4)
#' dim(X) = c(200,1,4)

#' n = 4
#' p = (1:12)*0; dim(p) = c(4,3);p[,1] = 2; p[,2]=1;   p[,3]=1; p[2,2]=2;
#' p    ## country-wise lag specification
#'
#' res_d = GVARData(m=2,n=4,p=p,T=200,type="exog0",X=X)
#' res_e = GVARest(res = res_d)
#' res_e$Summary
#'
#' IRF_CB = irf_GVAR_CB(res_e,nstep=10,comb=NA,irf="gen",runs=200,conf=c(0.05,0.95))
#' dim(IRF_CB)
#' IRF_g = IRF_graph(IRF_CB,Names=NA,response=c(1,4),impulse=c(1,2,3,4), ncol=4)
#'
#' @export
GVARData <- function (m, n, p, T, W = NA, r_npo = NA, Ao = NA, Bo = NA, Co = NA, Uo = NA, Sigmao = NA, type = NA, X = NA, mu = NA)
{
  if (missing(Bo)) {
    Bo = NA
  }
  if (missing(Sigmao)) {
    Sigmao = NA
  }
  if (missing(Uo)) {
    Uo = NA
  }
  if (missing(W)) {
    W = NA
  }
  if (missing(Ao)) {
    Ao = NA
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
  Pmax = max(p[, 1:2])
  P = max(p)
  if (!anyNA(X))
    k = dim(X)[2]
  if (anyNA(Bo)) {
    Bo = (1:(m * m * Pmax * n)) * 0
    dim(Bo) = c(m, m, Pmax, n)
    r_npo = c(1:(m * Pmax * n)) * 0
    dim(r_npo) = c(m, Pmax, n)
    for (i in 1:n) {
      VARD = VARData(m, p[i, 1], T)
      Bo[, , 1:p[i, 1], i] = VARD$B
      r_npo[, 1:p[i, 1], i] = VARD$r_np
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
      for (L in 1:p[i, 2]) Ao[, , L, i] = matrix((stats::runif(m * m) - 0.5)/10, m, m)
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
      for (L in 1:Pmax) G[(1 + (i - 1) * m):(i * m), (1 + (j - 1) * m):(j * m), L] = Ao[, , L, i] * W[i,
                                                                                                      j]
    }
  }
  for (i in 1:n) {
    for (L in 1:Pmax) G[(1 + (i - 1) * m):(i * m), (1 + (i - 1) * m):(i * m), L] = Bo[, , L, i]
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
    }
    else {
      H = diag(n * m)
      for (L in 1:Pmax) H = H - G[, , L]
      CoV = as.vector(Co)
      muV = solve(H) %*% CoV
      mu = muV
      dim(mu) = c(m, 1, n)
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
    for (i in 1:n) Ct[, ((i - 1) * m + 1):((i - 1) * m + m)] = as.matrix(X[, , i]) %*% CoV[, ((i - 1) * m +
                                                                                                1):((i - 1) * m + m)]
    mu = NA
  }
  else {
    if (type == "exog1") {
      if (anyNA(Co)) {
        Co = matrix(stats::rnorm(m * n * (k + 1)), m * (k + 1), n)
        dim(Co) = c(m, k + 1, n)
        for (i in 1:n) if (p[i, 3] < k)
          Co[, (p[i, 3] + 2):(k + 1), i] = Co[, (p[i, 3] + 2):(k + 1), i] * 0
      }
      DMCo = dim(Co[, -1, ])
      CoV = matrix(0, dim(X)[2] + 1, m * n)
      for (i in 1:(1 + dim(X)[2])) CoV[i, ] = as.vector(Co[,
                                                           i, ])
      Ct = matrix(0, dim(X)[1], m * n)
      for (i in 1:n) Ct[, ((i - 1) * m + 1):((i - 1) * m + m)] = cbind(matrix(1, dim(X)[1], 1), X[, , i]) %*% CoV[, ((i - 1) * m + 1):((i - 1) * m + m)]
      mu = NA
    }
    else {
      X = NA
    }
  }
  Y = Uo + Ct
  for (t in (P + 1):T) {
    for (L in 1:Pmax) Y[t, ] = Y[t, ] + Y[t - L, ] %*% t(G[, , L])
  }
  check = max(abs(Y))
  result = list(Y, X, Uo, G, Sigmao, r_npo, Ao, Bo, Co, W,
                m, n, p, mu, check, type)
  names(result) = c("Y", "X", "Uo", "G", "Sigmao", "r_npo", "Ao", "Bo", "Co", "W", "m", "n", "p", "mu", "check", "type")
  return(result)
}

#' Estimation of GVAR(m,n,p)
#'
#' This function estimates the parameters of a specified GVAR(m,n,p) model based on provided data.
#'
#' @param  res  : an GVAR object that is an output of GVARData including at least: m,n,p,type,Y and optionally X.
#' @return res  : an GVAR object with estimated parameter values, AIC, BIC, AIC_g, BIC_g and LH, where AIC and BIC are n-vectors of the country equations' AIC and BIC and AIG_g and BIC_g are the GVAR information criteria respectively.
#' @examples
#' n = 5
#' p = (1:15)*0; dim(p) = c(5,3)
#' p[,1] = 2; p[,2]=1;
#' res_d = GVARData(m=2,n=5,p=p,T=100,type="const")
#' max(res_d$Y)
#' dim(res_d$Y)
#' res_e = GVARest(res = res_d)
#' res_e$Summary
#'
#' X1 = matrix(1,200,1)
#' X2 = matrix(rnorm(200),200,1)
#' X3 = matrix(rnorm(200),200,1)
#' X4 = matrix(rnorm(200),200,1)
#' X  = cbind(X1,X2,X3,X4)
#' dim(X) = c(200,1,4)

#' n = 4
#' p = (1:12)*0; dim(p) = c(4,3);p[,1] = 2; p[,2]=1;   p[,3]=1; p[2,2]=2;
#' p
#'
#' res_d = GVARData(m=2,n=4,p=p,T=200,type="exog0",X=X)
#' res_e = GVARest(res = res_d)
#' res_e$Summary
#'
#' IRF_CB = irf_GVAR_CB(res_e,nstep=10,comb=NA,irf="gen",runs=200,conf=c(0.05,0.95))
#' dim(IRF_CB)
#' IRF_g = IRF_graph(IRF_CB,Names=NA,response=c(1,4),impulse=c(1,2,3,4), ncol=4)
#'
#'
#' @export
GVARest <- function (res)  {
  m = res$m
  n = res$n
  p = res$p
  Y = res$Y
  X = res$X
  W = res$W
  type = res$type
  Bo = res$Bo
  Ao = res$Ao

  Co = res$Co
  Pmax = max(p[, 1:2])
  k = max(p[, 3])
  AIC = c(1:n) * 0
  BIC = c(1:n) * 0
  LH = c(1:n) * 0
  resid = Y * 0
  VAR_domestic = list()
  T = dim(Y)[1]
  #if (is.null(colnames(Y))) colnames(Y) = sprintf("Y%s", 1:n)
  FY = Y %*% t(W %x% diag(m));
  #colnames(FY) =  paste("F",colnames(Y),sep="")
  Bo = Bo = (1:(m * m * Pmax * n)) * 0
  dim(Bo) = c(m, m, Pmax, n)
  Ao = Bo * 0
  if (type == "none")
    C = matrix(0, m * n, 1) * 0
  if (type == "const")
    C = matrix(0, m * n, 1) * 0
  if (type == "exog0")
    C = matrix(0, m * n, dim(X)[2] + 1) * 0
  if (type == "exog1")
    C = matrix(0, m * n, dim(X)[2] + 1) * 0
  for (i in 1:n) {
    res_VAR = VARData(n = m, p = p[i, 1], T = T)
    Z = matrix(0, T, p[i, 2] * m)
    HilfeFY = FY[, (m * (i - 1) + 1):(i * m)]
    if (is.null(colnames( HilfeFY)))  colnames(HilfeFY) = sprintf("FY%s", 1:m)
    ##FYp = embed(FY[, (m * (i - 1) + 1):(i * m)], (Pmax + 1))
    FYp = Embed(HilfeFY, (Pmax + 1),prefix="")

    Z[(Pmax + 1):T, ] = FYp[, (m + 1):(p[i, 2] * m + m)]
    HilfeZ            = FYp[, (m + 1):(p[i, 2] * m + m)]
    colnames(Z) <- colnames(HilfeZ)

    if (type == "none")
      Z = Z
    if (type == "const")
      Z = Z
    if ((type == "exog0")|(type == "exog1")) {
      Xi <- as.matrix(X[, 1:p[i, 3], i])
      if (is.null(colnames(Xi)))  colnames(Xi) = sprintf("exog%s", 1:ncol(Xi))
      Z = cbind(Z, Xi)
    }
    ki = dim(Z)[2]
    res_VAR$Y = Y[, (m * (i - 1) + 1):(i * m)]
    res_VAR$X = Z
    if ((type == "none") | (type == "exog0")) {
      res_VAR$type = "exog0"
    }
    if ((type == "const") | (type == "exog1")) {
      res_VAR$type = "exog1"
    }
    RR = VARest(res_VAR)
    AIC[i] = RR$Summary$AIC
    BIC[i] = RR$Summary$BIC
    LH[i] = RR$Summary$LH
    VAR_domestic[[i]] = RR$varp
    Bo[, , 1:p[i, 1], i] = RR$B
    if (type == "none")
      Ai = RR$Co[, -1]
    if (type == "const")
      Ai = RR$Co[, -1]
    if (type == "exog0")
      Ai = RR$Co[, c(2:(m * p[i, 2] + 1))]
    if (type == "exog1")
      Ai = RR$Co[, c(2:(m * p[i, 2] + 1))]
    dim(Ai) = c(m, m, p[i, 2])
    Ao[, , 1:p[i, 2], i] = Ai
    if (type == "none") {
      C[(m * (i - 1) + 1):(i * m), ] = RR$Co[, 1]
      Co[, , i] = RR$Co[, 1]
    }
    if (type == "const") {
      C[(m * (i - 1) + 1):(i * m), ] = RR$Co[, 1]
      Co[, , i] = RR$Co[, 1]
    }
    if (type == "exog0") {
      C[(m * (i - 1) + 1):(i * m), 1:(1 + p[i, 3])] = RR$Co[,
                                                            c(1, (m * p[i, 2] + 2):(m * p[i, 2] + 1 + p[i,
                                                                                                        3]))]
      Co[, 1:(p[i, 3] + 1), i] = RR$Co[, c(1, (m * p[i,
                                                     2] + 2):(m * p[i, 2] + 1 + p[i, 3]))]
    }
    if (type == "exog1") {
      C[(m * (i - 1) + 1):(i * m), ] = RR$Co[, c(1, (m *
                                                       p[i, 2] + 2):(m * p[i, 2] + 1 + p[i, 3]))]
      Co[, , i] = RR$Co[, c(1, (m * p[i, 2] + 2):(m * p[i,
                                                        2] + 1 + p[i, 3]))]
    }
    resid[, (m * (i - 1) + 1):(i * m)] = RR$resid
  }
  Sigmao = t(resid) %*% resid/(T - m * (p[i, 1] + p[i, 2]) - p[i, 3])
  Sigmao = as.matrix(Matrix::nearPD(Sigmao)[[1]])

  G = BoAoW2G(Bo, Ao, W, m, n, Pmax)
  Gs = diag(n * m)
  for (L in 1:Pmax) {
    Gs = Gs - G[, , L]
  }
  LH_g = -T * m * n/2 * log(2 * pi) - T * m * n/2 + T/2 * log(det(solve(Sigmao)))
  TwoN = sum(AIC) + 2 * sum(LH) - n * m * (m + 1) + (n * m) *
    (n * m + 1)
  AIC_g = TwoN - 2 * LH_g
  BIC_g = log(T) * TwoN/2 - 2 * LH_g
  res$G = G
  res$C = C
  res$Sigmao = Sigmao
  res$r_npo = NA
  res$Ao = Ao
  res$Bo = Bo
  res$Co = Co
  res$mu = solve(Gs) %*% C
  res$VAR_domestic = VAR_domestic
  est_result <- list()
  res$resid = resid
  for (i in 1: n)  est_result[[i]]= summary(VAR_domestic[[i]])
  Summary = list(est_result, LH, AIC, BIC, LH_g, AIC_g, BIC_g)
  names(Summary) = c("Estimation_Result", "Country_LH_function_Value", "Country_AIC", "Country_BIC", "LH_g", "AIC_g", "BIC_g")
  res$Summary = Summary
  return(res)
}



#' Transformation of country GVAR parameters to the Global VAR parameters
#'
#' @param Bo parameter matrix of domestic variables
#' @param Ao parameter matrix of foreign variables
#' @param W the weighting vector of the GVAR model
#' @param m the number of variables
#' @param n the number of countries
#' @param p lags
#'
#' @return An (mn, mn, p ) array of the GVAR(m,n,p) coefficients.
#' @export
#' @keywords internal
#'
BoAoW2G = function(Bo,Ao,W,m,n,p) {
  dim(Bo) = c(m,m,p,n)
  dim(Ao) = c(m,m,p,n)
  G = (1:(n*m*n*m*p))*0
  dim(G) = c(n*m,n*m,p)
  for (i in 1:n) {
      for (j in 1:n) {
         for (L in 1:p)    G[(1+(i-1)*m):(i*m),(1+(j-1)*m):(j*m),L] = Ao[,,L,i]*W[i,j]
      }
  }
  for (i in 1:n) {
	for (L in 1:p)       G[(1+(i-1)*m):(i*m),(1+(i-1)*m):(i*m),L]  = Bo[,,L,i]
  }
  dim(G) = c(n*m,n*m,p)
return(G)
}







#' Transformation of a global VAR parameter matrix to the country GVAR parameter matrices
#'
#' @param G The GVAR parameter matrix
#' @param W The weighting matrix of the GVAR model
#'
#' @return A list containing the parameter matrices of domestric variables and foreign variables
#' @export
#' @keywords internal
#'
GW2BoAo = function(G,W) {
  n = dim(W)[1]
  m = dim(G)[1]/n

  if (is.na(dim(G)[3])) dim(G) = c(m*n,m*n,1)
  p = dim(G)[3]
  Bo = (1:(m*m*p*n))*0
  dim(Bo) = c(m,m,p,n)
  Ao = Bo
  for (i in 1:n) {
      for (j in 1:n) {
         for (L in 1:p)   if (!j==i) Ao[,,L,i] = G[(1+(i-1)*m):(i*m),(1+(j-1)*m):(j*m),L]/W[i,j]
      }
  }
  for (i in 1:n) {
	for (L in 1:p)       Bo[,,L,i] = G[(1+(i-1)*m):(i*m),(1+(i-1)*m):(i*m),L]
  }
  result = list(Bo,Ao)
  names(result) = c("Bo","Ao")
return(result)
}


#' Impulse Response Functions of GVAR
#'
#' This function generates impulse response functions of an estimated GVAR
#'
#' @param  res  : an output of GVARest
#' @param  nstep : length of the impulse response functions
#' @param  comb : an mn vector specifying combined impulse such as global shocks, regional shocks, or concerted actions.
#' @param  irf  : types of the impulse response functions. irf=c("gen","chol","chol1","gen1","comb1"), gen for generalized IRF with one standard deviation shocks, gen1 for generalized IRF with one unit impulse, chol for IRF with Cholezky decomposition of the covariance matrix, chol1 for Cholezky decomposition with one unit impulse, comb1 for concerted action with one unit impulse.
#' @return an (mn,mn,nstep) array containing the IRF with columns representing the impulses rows representing the responses.
#' @examples
#'
#'
#' X1 = matrix(1,200,1)
#' X2 = matrix(rnorm(200),200,1)
#' X3 = matrix(rnorm(200),200,1)
#' X4 = matrix(rnorm(200),200,1)
#' X  = cbind(X1,X2,X3,X4)
#' dim(X) = c(200,1,4)
#'
#' n = 4
#' p = (1:12)*0; dim(p) = c(4,3);p[,1] = 2; p[,2]=2; p[,3]=1;
#' res_d = GVARData(m=2,n=4,p=p,T=200,type="exog0",X=X)
#' res_e = GVARest(res=res_d)
#' IRF = irf_GVAR(res_e,nstep=10,comb=NA,irf="gen")
#' IRF_CB = irf_GVAR_CB(res_e,nstep=10,comb=NA,irf="gen",runs=200,conf=c(0.05,0.95))
#' dim(IRF_CB)
#' IRF_g = IRF_graph(IRF_CB,Names=NA,response=c(1,4),impulse=c(1,2,3,4), ncol=4)
#'
#' @export
#' @keywords internal
#'
irf_GVAR = function(res,nstep,comb,irf=c("gen","chol","chol1","gen1","comb1")) {
	B 	= res$G
      neq 	= dim(B)[1]
	nvar	= dim(B)[2]
	sigma = res$Sigma
      response <- array(0,dim=c(neq,nvar,nstep));
      response <- irf_B_sigma(B,sigma,nstep,comb,irf=irf)
	return(response)
}



#' Impulse Response Functions of GVAR
#'
#' This function generates impulse response functions of an estimated GVAR with confidence bands
#'
#' @param  res  : an object of GVAR that is a list of the output of GVARest
#' @param  nstep : length of the impulse response functions
#' @param  comb : an mn vector specifying combined impulse such as global shocks, regional shocks, or concerted actions.
#' @param  irf  : types of the impulse response functions. irf=c("gen","chol","chol1","gen1","comb1"), gen for generalized IRF with one standard deviation shocks, gen1 for generalized IRF with one unit impulse, chol for IRF with Cholezky decomposition of the covariance matrix, chol1 for Cholezky decomposition with one unit impulse, comb1 for concerted action with one unit impulse.
#' @param  runs : number of bootstrap runs to generate the confidence bands
#' @param  conf : a two components vector of the tail probabilities of the confidence interval.
#' @return An (mn,mn,nstep,3) array of IRF with columns representing the impulse rows the responses.
#' @examples
#' n = 5
#' p = (1:15)*0; dim(p) = c(5,3)
#' p[,1] = 2; p[,2]=1;
#' res_d = GVARData(m=2,n=5,p=p,T=100,type="const")
#' max(res_d$Y)
#' dim(res_d$Y)
#' res_e = GVARest(res = res_d)
#' res_e$Summary
#'
#' IRF_CB = irf_GVAR_CB(res_e,nstep=10,comb=NA,irf="gen",runs=200,conf=c(0.05,0.95))
#' dim(IRF_CB)
#' IRF_g = IRF_graph(IRF_CB,Names=NA,response=c(1,4),impulse=c(1,2,3,4), ncol=4)
#'
#' @export
irf_GVAR_CB =function(res,nstep,comb,irf=c("gen","chol","chol1","gen1","comb1"),runs=200,conf=c(0.05,0.95)) {
        m = res$m
        n = res$n
      p = res$p
      T = dim(res$Y)[1]
      W = res$W
      Ao= res$Ao
      Bo= res$Bo
      Go= res$G
      Co= res$Co
      type=res$type
      X   = res$X
      mu  = res$mu

        B       = res$G
      neq       = dim(B)[1]
        nvar    = dim(B)[2]
        sigma = res$Sigmao
      response <- array(0,dim=c(neq,nvar,nstep,length(conf)+1))
     #response[,,,1] <- impulsdtrf_g(B,sigma,nstep,comb,irf=irf)
      response[,,,1] <- irf_GVAR(res,nstep,comb,irf)
      responseR <- array(0,dim=c(neq,nvar,nstep,runs))
      for (i in 1:runs) {
        Uo_run = rnormSIGMA(T,sigma)
            res_run = GVARData(m,n,p,T,W,r_npo=NA,Ao,Bo,Co,Uo=Uo_run,Sigmao=NA,type,X,mu)
            res_e   = GVARest(res_run)
            B_run   = res_e$G
            sigma_run = res_e$Sigmao
                #responseR[,,,i] <- irf_B_sigma(B_run,sigma_run,nstep,comb,irf=irf)
            responseR[,,,i]  <- irf_GVAR(res_e,nstep,comb,irf)
        }
        responseR[,,,1] = response[,,,1]
      for (tt in 1:(nstep) ) {
        for (i in 1:neq)           {
                for (j in 1:nvar)     {Ao
                        response[i,j,tt,-1] = stats::quantile(responseR[i,j,tt,], conf)

        }
        }
        }
        return(response)
}





############# GVAR_Selection ##################

#' GVAR lag selection
#'
#' Calculation of the information criteria of a GVAR country models for a given range of maximum lags.
#'
#'
#' @param  res  : a GVAR object obtained from GVARData or GVARest.
#' @param  L_V  : a two components vector containing the maxima of the domestic lag and the foreign lag, respectively.
#' @param  I    : Index of the country under investigation.

#' @return A matrix with different lag specifications and values of the model information criteria.
#' @examples
#'
#' n = 4
#' p = (1:12)*0; dim(p) = c(4,3);p[,1] = 2; p[,2]=1; p[2:3,2] = 2
#' res_d = GVARData(m=2,n=4,p=p,T=4000,type="const")
#'
#'
#' I = 3
#' L_V = c(4,4)
#' res_d$p
#' GVARSelect = GVAR_Select(res=res_d,L_V=c(4,4),I=2)
#' GVARSelect[which.min(GVARSelect[,3]),]
#' GVARSelect[which.min(GVARSelect[,4]),]
#'
#' @export
GVAR_Select = function(res=res,L_V=L_V,I=I)  {
m = res$m
n = res$n
Y = as.matrix(res$Y)
T = dim(Y)[1]
type = res$type
XXX  = res$X

Wmat = res$W
Wnmat = kronecker(Wmat,diag(m))
FYI   = Y%*%Wnmat[,((I-1)*m+1):(I*m)]
Yi = as.matrix(Y[,((I-1)*m+1):(I*m)])
Criteria = matrix(0,L_V[1]*L_V[2],4)
idx = 0

for (l_d in 1: L_V[1] )   {
   for (l_f in 1:L_V[2] ) {
      idx = idx + 1
	X = stats::embed(FYI,l_f+1)[,-(1:m)]
	XX = matrix(0,T,dim(X)[2])
	XX[(l_f+1):T,] = X
      if (type=="none")  { type_holder = "exog0"; Co = (1:(m*(1+l_f*m))); dim(Co) = c(m,l_f*m+1); Co[,1]=0}
      if (type=="const") { type_holder = "exog1"; Co = (1:(m*(1+l_f*m))); dim(Co) = c(m,l_f*m+1); Co[,1]=1}


	if (l_d>=l_f) {
            if (type=="exog0") { type_holder = type; XX = cbind(XX,XXX[,,I]); Co = matrix(1,m,1+dim(XX)[2]); Co[,1]=0}
	      if (type=="exog1") { type_holder = type; XX = cbind(XX,XXX[,,I]); Co = matrix(1,m,1+dim(XX)[2]); Co[,1]=1}
		res_dd = VARData(n=m,p=l_d,T=T,Co=Co,type=type_holder,X=as.matrix(XX))
		res_dd$Y = Yi
	}   else  {
            XX = XX[(l_f-l_d+1):T,];
	   	if (type=="exog0") { type_holder = type; XX = cbind(XX,XXX[(l_f-l_d+1):T,,I]);Co = matrix(1,m,1+dim(XX)[2]); Co[,1]=0 }
	   	if (type=="exog1") { type_holder = type; XX = cbind(XX,XXX[(l_f-l_d+1):T,,I]);Co = matrix(1,m,1+dim(XX)[2]); Co[,1]=1 }
		res_dd = VARData(n=m,p=l_d,T=(T-(l_f-l_d)),Co=Co,type=type_holder,X=XX)
		res_dd$Y = Yi[(l_f-l_d+1):T,]
	}
	res_e = VARest(res=res_dd)
	Criteria[idx,] = c (l_d,l_f,res_e$Summary$AIC,res_e$Summary$BIC)
    }
}
colnames(Criteria) = c("Domestic Lags","Foreign Lags","AIC","BIC")
return(Criteria)
}


#' Transformation of a weighting vector to a weighting matrix
#'
#' @param W the weighting vector
#'
#' @return the corresponding weighting matrix
#' @export
#' @keywords internal
#'
W2Wmat = function(W) {
 n = length(W)
 Wmat = matrix(0,n,n)
 for (i in 1:n ) {
    Wmat[i,-i] = W[-i]/sum(W[-i])
 }
 return(Wmat)
}

#' Transformation of the weighting vector to a weighting matrix
#'
#' @param SW The weighting vector
#' @param n The number of variables
#' @param N The number of countries
#' @param K The number
#'
#' @return The weighting matrix
#' @export
#' @keywords internal
#'
SW2comb = function(SW,n,N,K)
{
  comb = matrix(0,n*N,N*n)
  SW = SW/sum(SW)
  for (i in 1:n) {
    comb[(i-1)*N+K,1] = SW[i]
  }
  return(comb)
}

