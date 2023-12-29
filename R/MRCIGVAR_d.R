#' Data generating process of MRCIGVAR(m,n,p,S,r)
#'
#' This function generates data from an multi-regime cointegrated global VAR process and returns an MRCIGVAR(m,n,p,S) object that is a list containing the data and the parameters used in the MRCIGVAR(m,n,p,S) process.
#'
#' @param m     : number of variables in each country/unit
#' @param n     : number of countries/units
#' @param S     : number of regimes
#' @param p     : an (n, 3, S) array, each raw contains the lag length of the domestic variables, the lag length of the foreign variables, and the number of the exogenous variables of the corresponding country in the respective regimes.
#' @param T     : number of observations.
#' @param W     : an (n x n) weighting matrix. w_ij is the weight of country j in the foreign variables of i-th country diag(W)=0
#' @param SESVI : an n-vector of indices of the switching variables across n countries. E.g. SESVI = seq(1,m*n,m).
#' @param TH    : an (n x S-1) matrix of threshold values
#' @param Go    : an (mn,mn,p,S) array containing coefficients of MRCIGVAR. Go is constructed from Bo, Ao and W.
#' @param Ao    : an (m, m, p, n, S) array containing the coefficients of foreign variables
#' @param Bo    : an (m, m, p, n, S) array containing the coefficients of domestic variables.
#' @param Co    : an (m , k+1, n, S) array containing the coefficients of the deterministic components of the n countries.
#' @param Uo    : a (T, mn, S) array of the temporally independent innovations
#' @param Sigmao : an (mn, mn, S) array of the covariance matrix of MRCIGVAR(m,n,p,S)
#' @param SV    : exogenous switching variables
#' @param type	: types of the deterministic component. "const", "none", "exog0", and "exog1" are 4 options
#' @param X	    : a (T x k x n x S) array of exogenous stationary variables.
#' @param Yo    : initial values
#' @param d     : the time lag between signal and regime-switching
#' @param r     : an n-vector containing the number of unit root process in each country
#' @param r_np  : an (m, Pmax, n, S) array containing the roots of the characteristic polynomials of each country.
#' @param Ncommtrend : number of common stochastic trends in the MRCIGVAR.
#' @param DFYflag : indicator whether the foreign variables enter the cointegration space.
#' @param A       : transformation matrix for common exogenous non-stationary stochastic factors.
#' @param Ncommfakt : number of the common exogenous non-stationary stochastic factors.
#' @param uz      : innovations of the common exogenous non-stationary stochastic factors.
#' @return      an MRGVAR object containing the generated data, the used parameters and the exogenous variables.
#' \itemize{
#'    \item Y     : a (T x mn)  matrix of simulated data
#'    \item X     : a (T x k x n x S) array of exogenous variables
#'    \item Uo    : a (T, mn,S) array of the simulated innovations of the MRGVAR(m,n,p,S)
#'    \item C     : an (nm, (k+1),S) array containing the coefficients of the deterministic components.
#'    \item St    : a (T x n) matrix of the simulated time path of states/regimes
#'    \item check : maximum of the data for checking the stationarity
#' }
#'
#'
#' @examples
#' m = 3
#' n = 5
#'  p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' p = p[1:n,,]; p[,1,] = 3; p[,2,] = 2
#'
#' p[2,2,] = 3
#' TH = c(1:n)*0; dim(TH) = c(1,n)
#' SESVI=rep(1,4,7,10,13)
#' r  = rep(1,n)
#' i = 1
#' S = 2
#'
#' T= 200
#' XX = matrix(rnorm(T*10*2),T*10,2)
#' dim(XX) = c(T,2,5,2)
#' p[,3,]=2
#'
#' res_d <- MRCIGVARData(m=3,n=5,p=p,TH=TH,T=T,S=2, SESVI=c(1,4,7,10,13),type="exog0",r=rep(1,5),
#' DFYflag=0,Ncommfakt=1,X=XX)    ## m: number of variables, n: number of countries
#'
#' max(abs(res_d$Y))
#' plot(ts(res_d$Y[,1:10]))
#' STAT(res_d$Go[,,,2])
#' STAT(res_d$GC[,,,1]);STAT(res_d$GC[,,,2]);
#' STAT(res_d$Go[,,,1])
#' max(abs(res_d$Y))
#'
#' Ao=res_d$Ao;Bo=res_d$Bo;Co=res_d$Co;
#'
#' res_d <- MRCIGVARData(m=3,n=5,p=p,TH=TH,T=T,S=2, SESVI=c(1,4,7,10,13),type="exog0",Ao=Ao,Bo=Bo,
#' Co=Co,r = rep(1,5),DFYflag = 0,Ncommfakt=1,X=XX) ##m:number of variables, n:number of countries
#' res_e  = MRCIGVARest(res=res_d)
#'
#'
#' @export
MRCIGVARData <- function(m,n,p,T,S,W=NA,SESVI=NA,TH=NA,Go=NA,Ao=NA,Bo=NA,Sigmao=NA,Uo=NA,SV=NA,type=NA,Co=NA,X=NA,Yo=NA,d=NA,r=NA,r_np=NA,Ncommtrend=NA,DFYflag=0,A=NA,Ncommfakt=1,uz=NA) {
  ### m     : number of variables in a country
  ### n     : number of countries
  ### p     : lag length n x 3 x S array collecting lag length for each country with respect to domestic and foreigen variable for each state. The last column specifies the number of exogenous variabls for each state.
  ### T     : number of observations
  ### S     : number of regimes
  ### SESVI : n vector of indeces of the switching variable in the endogeneous variables Y for the case of self-excited threshhold model.
  ###         (m,n,p,T,S,SESVI) are parameters which must be provided.
  ### TH    : n x S-1  matrix of threshold values of n countries.
  ### Go    : mn x mn x p x S GVAR(m,n,p) coefficients matrix of S different regimes
  ### Sigmao: mn x mn x S array of the covariance matrix of the GVAR(m,n,p) in S different regimes
  ### Uo    : a T x mn x S  array of the temorally independent innovation processes
  ### X       : a T x k x n x S array of exogenous veriables which may be common/different to all countries and different across all states
  ###         (TH,Go,Sigmao,Uo,SV) if not provided, they will be generated randomly.
  ### type  : deterministic component "const" and "none" "exog0" "exog1" are foure options
  ### Co    : if type = "const" mu is an m x n x S array of the intercepts of the time series in the different regimes
  ###
  ###
  ### output:
  ### c("Y","Uo","Go","Sigmao","TH","St","sv","SESVI","check")
  ###
  ### Y     : T x mn data matrix of the simulated data via of the MSGVAR(m,n,p,S)
  ### Uo    : an T x mn x S   array of the temorally independent innovation processes of the MSGVAR(m,n,p,S)
  ### Go    : mn x mn x p x S array collecting the MSGVAR(m,n,p,S) coefficients in S different states
  ### Sigmao: mn x mn x S array collecting the covariance matrices of the simulated MSGVAR(m,n,p,S) in S different states
  ### TH    : S-1 vector of thresholds
  ### St    : simulated time path of states/regimes
  ### SESVI : index if the switching variable
  ### check : maximum of the data to check stationarity
  ###
  ### Remarks: The states of each country at each time step is governed by the switching variables of each country and hence there is a large numbers
  ###          of possible combinations of gegimes.
  ###          The coefficients matrix can be constructed from Go, given a regime combination. The stationarity of the GVAR(m,n,p)
  ###          at each time step is not guaranteed. But this is not relevant. (This is a open question.)
  ###
  ###
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
  if (missing(Ao)) {
    Ao = NA
  }
  if (missing(X)) {
    X = NA
  }
  if (missing(d)) {
    d = NA
  }
  if (missing(Yo)) {
    Yo = NA
  }
  if (missing(r)) {
    r = NA
  }
  if (missing(r_np)) {
    r_np = NA
  }
  if (missing(Ncommtrend)) {
    Ncommtrend = NA
  }
  if (missing(DFYflag)) {
    DFYflag = 0
  }
  if (missing(Ncommfakt)) {
    Ncommfakt  = 0
  }



  if (anyNA(d)) {
    d = 1
  }
  if (anyNA(r)) {
    r = rep(1, n)
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
      GVARD = GVARData(m, n, p[, , s], T)
      Sigmao[, , s] = GVARD$Sigmao
    }
  }
  if (anyNA(r_np)) {
    r_np = (1:(m * Pmax * n * S))/(1:(m * Pmax * n * S)) *
      (3.1 + abs(stats::rnorm(m * Pmax * n * S)))
    dim(r_np) = c(m, Pmax, n, S)
    r_np[1, , , ] = 1
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
    for (s in 1:S) {
      for (i in 1:n) {
        if (p[i, 3, s] < k - 1)
          Co[, (p[i, 3, s] + 1):k, i, s] = 0
        for (j in 1:k) CoV[j, ] = as.vector(Co[, j, ,
                                               s])
        Ct[, , s] = cbind(matrix(1, dim(X)[1], 1), X[,
                                                     , i, s]) %*% CoV
      }
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
    for (s in 1:S) {
      for (i in 1:n) {
        if (p[i, 3, s] < k - 1)
          Co[, (p[i, 3, s] + 2):k, i, s] = 0
        for (j in 1:k) CoV[j, ] = as.vector(Co[, j, ,
                                               s])
        Ct[, , s] = cbind(matrix(1, dim(X)[1], 1), X[,
                                                     , i, s]) %*% CoV
      }
    }
  }

  kkk    <- 0
  CIVARB <- NA
  ### Go = NA; Bo = NA; Ao = NA
  if (anyNA(Go) & anyNA(Bo)) {
    repeat {
      CIVARB = (1:(2 * m * 2 * m * Pmax * S * n)) * 0
      dim(CIVARB) = c(2 * m, 2 * m, Pmax, n, S)
      Go = (1:(m * n * m * n * Pmax * S)) * 0
      dim(Go) = c(m * n, m * n, Pmax, S)
      Bo = (1:(m * m * Pmax * n * S)) * 0
      dim(Bo) = c(m, m, Pmax, n, S)
      Bo = VARBS_commtrend(m, p, T, r, Ncommtrend, n, S, r_npo = r_np)[[1]]
      Ao = Bo2Ao(m, p, T, r, Ncommtrend, n, S, r_npo = r_np,DFYflag=DFYflag)
      BoAo_h = BoAo2TBoAo(Bo, Ao, W, S)
      Bo = BoAo_h[[1]]
      Ao = BoAo_h[[2]]
      for (i in 1:n) {
        Bi = matrix(stats::rnorm(m * m), m, m)
      }
      for (s in 1:S) Go[, , , s] = BoAoW2G(Bo[, , , , s], Ao[, , , , s], W, m, n, Pmax)

      if (!max(Mod(STAT(Go[,,,1])),Mod(STAT(Go[,,,2])))>1.000000001 ) break
      kkk = k + 1
      if (kkk > 100) break
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
    for (s in 1:S) Go[, , , s] = BoAoW2G(Bo[, , , , s], Ao[, , , , s], W, m, n, Pmax)
  }


  St = matrix(0, T, n)
  Y = Uo[, , 1]
  if (!anyNA(Yo))
    Y[1:P, ] = Yo
  resid = Y * 0
  Yo = Y[1:P, ]
  if (anyNA(SV)) {
    sv = Y[, SESVI] * 0
    SV = NA
  }
  else {
    sv = SV
  }
  ss = matrix(0, n, 1)
  for (tt in (P + 1):T) {
    if (anyNA(SV)) {
      sv[2:T, ] = Y[2:T, SESVI] - Y[1:(T - 1), SESVI]
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
    dim_GC = dim(Go); dim_GC[1:2] = dim_GC[1:2]+Ncommfakt
    GC = c(1:prod(dim_GC))*0; dim(GC) = dim_GC
    GC[(Ncommfakt+1):(Ncommfakt+n*m),(Ncommfakt+1):(Ncommfakt+n*m),,] = Go
    GC[1:Ncommfakt,1:Ncommfakt,1,] = 1.4*diag(Ncommfakt); GC[1:Ncommfakt,1:Ncommfakt,2,] = -0.4*diag(Ncommfakt)  # assuming unit root common factor(s)
    for (i in 1:Pmax) GC[,,i,1] = A%*%GC[,,i,1]%*%solve(A)
    for (i in 1:Pmax) GC[,,i,2] = A%*%GC[,,i,2]%*%solve(A)
    Uo[,,1] = (cbind(uz,Uo[,,1])%*%t(A))[,Ncommfakt+(1:(m*n))]
    Uo[,,2] = (cbind(uz,Uo[,,2])%*%t(A))[,Ncommfakt+(1:(m*n))]
    Do = (1:(m*n*Pmax*Ncommfakt*2))*0; dim(Do) = c(m,Ncommfakt,Pmax,n,2)

  }  else  {
    Z  = NA
    uz = NA
    A  = NA
    GC = NA
    Do = NA
  }


  result = list(Y, X, Uo, resid, Go, GDC, Co, Sigmao, TH, St,SV, SESVI, Ao, Bo, check, type, m, n, p, S, W, SigmaS, Yo, d, r, CIVARB, kkk,DFYflag,Z,uz,A,GC,Do)
  names(result) = c("Y", "X", "Uo", "resid","Go", "GDC", "Co", "Sigmao","TH", "St", "SV", "SESVI", "Ao","Bo", "check", "type", "m", "n","p", "S", "W", "SigmaS", "Yo",
                    "d", "r", "CIVARB", "kkk","DFYflag","Z","uz","A","GC","Do")
  return(result)
}






#' Generating parameter matrices with common trends across countries
#'
#' @param m Number of variables of a country
#' @param p Lags
#' @param T Number of observations
#' @param r Number of stochastic trends
#' @param Ncommtrend Number of common trends
#' @param n Number of countries in the GVAR model
#' @param S Number of regimes
#' @param r_npo A matrix specifying the  roots of the characteristic polynomial in lags
#'
#' @return A list contain VAR parameters B, alpha, and beta
#' @export
VARBS_commtrend = function(m,p,T,r, Ncommtrend,n,S,r_npo) {
  alpha   = list()
  beta    = list()
  VAR_Ione = list()
  if (missing(Ncommtrend))  {Ncommtrend=NA}
  Pmax = max(p[,1:2,])
  Pmin = min(p[,1:2,][p[,1:2,]>0])
  Bo = (1:(m * m * Pmax * n * S)) * 0
  dim(Bo) = c(m, m, Pmax, n, S)


  for (i in 1:n) {
      Pmini             = min(p[i,1:2,][p[i,1:2,]>0])
	r_npi 		= c(1:(r[i]*Pmini))/c(1:(r[i]*Pmini))*1.5
      dim(r_npi) 		= c(r[i],Pmini)
      r_npi[1:r[i],1]  	= 1
	VAR_Ione[[i]]  	= VARData(n=r[i],p=Pmini,T=100,r_np=r_npi)
  }



  if (!anyNA(Ncommtrend)) {
     	r_npc  = c(1:(Ncommtrend*Pmin))/c(1:(Ncommtrend*Pmin))*1.5
      dim(r_npc) = c(Ncommtrend,Pmin)
      r_npc[1:Ncommtrend,1]  = 1
      VAR_IoneC = VARData(n=Ncommtrend,p=Pmin,T=100,r_np=r_npc)
   	for (i in 1:n) {
		VAR_Ione[[i]]$B[1:Ncommtrend,1:Ncommtrend,]=0
            VAR_Ione[[i]]$B[1:Ncommtrend,1:Ncommtrend,1:Pmin] = VAR_IoneC$B
	}
  }

      for (i in 1:n) {
 		NIzero = m - r[i]
         	#B =  matrix(rnorm(m*m),m,m)                 ### jizhu 0 nadiao
      	#B[(NIzero+1):m,] = 0
      	#B = B + diag(m)
            B = diag(m)
            #B = eigen(t(B)%*%B)$vectors
            alphai = list()
            betai = list()
            for (s in 1:S) {
      		r_np = c(1:(NIzero * p[i,1,s]))*0
            	dim(r_np) = c(NIzero,p[i,1,s])
            	r_np <- r_npo[(r[i]+1):m,1:p[i,1,s],i,s]
            	dim(r_np) = c(NIzero,p[i,1,s])
                  #A = matrix(rnorm(NIzero*NIzero),NIzero,NIzero); A = eigen(t(A)%*%A)[[2]]
	            VAR_IZero = VARData(NIzero, p[i, 1,s],T=100,r_np=r_np,)

            	for (j in 1: p[i,1,s]) {
               		Dblock    = matrix(0,m,m)
               		Dblock[1:NIzero,1:NIzero]     = VAR_IZero$B[,,j]
               		if (dim(VAR_Ione[[i]]$B)[3]>=j) Dblock[(NIzero+1):m,(NIzero+1):m] = VAR_Ione[[i]]$B[,,j]
               		Bo[, , j, i,s] = B%*%Dblock%*%solve(B)
      		}
            	#alpha[[i]] = B2CIB(VARD$B)[[2]]
            	#beta[[i]]  = B2CIB(VARD$B)[[3]]
			B2CIBab  	= B2CIB(Bo[, , 1:p[i, 1,s], i,s])
            	alphai[[s]] = B2CIBab[[2]]
            	betai[[s]]  = B2CIBab[[3]]
      	}
            alpha[[i]] = alphai
		beta[[i]]  = betai
  	}



  return(list(Bo,alpha,beta))
}


################################################################################################
#' Transformation from Bo with common trends to Ao
#'
#' This function calculates Ao coefficients from Bo with common trends.
#'
#' @param m Number of variables
#' @param p Lags
#' @param T Number of observations
#' @param r Number of stochastic trends
#' @param Ncommtrend Number of the common trends
#' @param n Number of countries
#' @param S Number of regimes
#' @param r_npo The matrix containing the roots of the characteristic polynomial
#' @param DFYflag switching indicator whether the foreign variables enter the cointegration space
#'
#' @return Parameter matrix of the foreign variables in the country equation of the GVAR model
#' @export
Bo2Ao <- function(m,p,T,r, Ncommtrend,n,S,r_npo,DFYflag=NA) {
  if ( anyNA(DFYflag) ) DFYflag = 0
  Boab  = VARBS_commtrend(m,p,T,r,Ncommtrend,n,S,r_npo)
  Pmax = max(p[,1:2,])
  beta  = Boab[[3]]
  alpha = Boab[[2]]
  Ao = (1:(m * m * Pmax * n * S)) * 0
  dim(Ao) = c(m, m, Pmax, n, S)
  for (i in 1:n) {
    for (s in 1:S) {
      if (p[i,2,s] < 2) Ao = Ao
      if (p[i,2,s] >= 2) {
        VARD = VARData(m, p=(p[i, 2,s]-1),T)
        BB = VARD$B
        Ao[,,1,i,s] = BB[,,1]+alpha[[i]][[s]]%*%t(beta[[i]][[s]])*DFYflag
        Ao[,,p[i,2,s],i,s] = -BB[,,p[i,2,s]-1]
        if ((p[i, 2,s]-1)>=2) { for (L in 2:(p[i, 2,s]-1)) Ao[, ,L, i,s] = BB[,,L]-BB[,,L-1]}
      }
    }
  }

  return(Ao)
}


#' A full rank linear transformation of initial Bo and Ao matrices
#'
#' This function calculates a full rank linear transformation of initial Bo and Ao matrices with independent common trends and idiosyncratic trends.
#'
#' @param Bo The initial Bo matrix
#' @param Ao The initial Ao matrix
#' @param W  A weighting matrix of the GVAR model
#' @param S Number of regimes
#'
#' @return Transformed Bo and Ao matrices
#' @export
BoAo2TBoAo = function(Bo,Ao,W,S) {
  TBo = Bo*0
  TAo = Ao*0
  n = dim(W)[1]
  m = dim(Bo)[1]
  L = dim(Bo)[3]
  BD = matrix(0,m*n,m*n)
  B =  matrix(stats::rnorm(m*m),m,m)
  for (i in 1:n) {
    BD[((i-1)*m+1):((i-1)*m+m),((i-1)*m+1):((i-1)*m+m)] = B
  }
  for (j in 1:L) {
    for ( s in 1: S ) {

      BoDj =  matrix(0,m*n,m*n)
      AoDj =  matrix(0,m*n,m*n)

      for (i in 1:n)  {
        BoDj[((i-1)*m+1):((i-1)*m+m),((i-1)*m+1):((i-1)*m+m)]=Bo[,,j,i,s]
        AoDj[((i-1)*m+1):((i-1)*m+m),((i-1)*m+1):((i-1)*m+m)]=Ao[,,j,i,s]
      }
      B_h = BD%*%BoDj%*%solve(BD)
      A_h = BD%*%AoDj%*%kronecker(W,diag(m))%*%solve(BD)%*%solve(kronecker(W,diag(m)))
      for (i in 1:n)  {
        TBo[,,j,i,s] = B_h[((i-1)*m+1):((i-1)*m+m),((i-1)*m+1):((i-1)*m+m)]
        TAo[,,j,i,s] = A_h[((i-1)*m+1):((i-1)*m+m),((i-1)*m+1):((i-1)*m+m)]
      }
    }
  }
  return(list(TBo,TAo))
}






#' Estimation of MRCIGVAR(res)
#'
#' This function estimates parameters of a multi regime cointegrated global VAR(p) model. The adjustment speeds to the cointegration relations are assumed to be different in different regimes.
#'
#' @param res   : an MRCIVAR object of the output of MRCIGVARData
#' @return      : an MRCIVAR object with estimated parameters and test statistics.
#' @examples
#' m = 2
#' n = 3
#' p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' p = p[1:n,,]; p[,1,] = 3; p[,2,] = 2
#'
#' TH = c(1:n)*0; dim(TH) = c(1,n)
#' SESVI=rep(1,3,5)
#' r  = rep(1,n)
#'
#' ## case of n = 3, m = 2, S = 2     ## m: number of variables, n: number of countries
#'
#' res_d <- MRCIGVARData(m=2,n=3,p=p,TH=TH,T=300,S=2, SESVI=c(1,3,5),r=rep(1,3),Ncommtrend=1)
#' # Ao
#' max(abs(res_d$Y))
#' plot(ts(res_d$Y))
#' STAT(res_d$Go[,,,2])
#' STAT(res_d$Go[,,,1])
#' max(abs(res_d$Y))
#' res_d$type
#' res_d$W
#' #   Ao =  res_d$Ao;  Bo =  res_d$Bo;  Co =  res_d$Co;
#' #   res_d <- MRCIGVARData(m=2,n=3,p=p,TH=TH,T=300,S=2, SESVI=c(1,3,5),Ao=Ao,Bo=Bo,Co=Co)
#'
#' max(res_d$Y)
#' res_e  = MRCIGVARest(res=res_d)
#' res_e$Summary
#'
#' @export
MRCIGVARest <- function(res) {
  ### m     : number of variables in a country
  ### n     : number of countries
  ### p     : lag length n x 3 x S array collecting lag length for each country with respect to domestic and foreign variable for each state. The last column specifies the number of exogenous variables for each state.
  ### T     : number of observations
  ### S     : number of regimes
  ### SESVI : n vector of indeces of the switching variable in the endogeneous variables Y for the case of self-excited threshold model.
  ###         (m,n,p,T,S,SESVI) are parameters which must be provided.
  ### TH    : n x S-1  matrix of threshold values of n countries.
  ### Go    : mn x mn x p x S GVAR(m,n,p) coefficients matrix of S different regimes
  ### Sigmao: mn x mn x S array of the covariance matrix of the GVAR(m,n,p) in S different regimes
  ### Uo    : a T x mn x S  array of the temorally independent innovation processes
  ### X	    : a T x k x n x S array of exogenous veriables which may be common/different to all countries and different across all states
  ###         (TH,Go,Sigmao,Uo,SV) if not provided, they will be generated randomly.
  ### type  : deterministic component "const" and "none" "exog0" "exog1" are foure options
  ### Co    : if type = "const" mu is an m x n x S array of the intercepts of the time series in the different regimes
  ###
  ###
  ### output:
  ### c("Y","Uo","Go","Sigmao","TH","St","sv","SESVI","check")
  ###
  ### Y     : T x mn data matrix of the simulated data via of the MSGVAR(m,n,p,S)
  ### Uo    : an T x mn x S   array of the temorally independent innovation processes of the MSGVAR(m,n,p,S)
  ### Go    : mn x mn x p x S array collecting the MSGVAR(m,n,p,S) coefficients in S different states
  ### Sigmao: mn x mn x S array collecting the covariance matrices of the simulated MSGVAR(m,n,p,S) in S different states
  ### TH    : S-1 vector of thresholds
  ### St    : simulated time path of states/regimes
  ### SESVI : index if the switching variable
  ### check : maximum of the data to check stationarity
  ###
  ### Remarks: The states of each country at each time step is governed by the switching variables of each country and hence there is a large numbers
  ###          of possible combinations of gegimes.
  ###          The coefficients matrix can be constructed from Go, given a regime combination. The stationarity of the GVAR(m,n,p)
  ###          at each time step is not guaranteed. But this is not relevant. (This is a open question.)
  ###
  ###
  m = res$m
  n = res$n
  p = res$p
  Y = as.matrix(res$Y);  if (max(abs(Y)) > 10000000) return("Check data!!!");
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
  P = max(p[, 1:2, ], d)
  Pmax = max(p[, 1:2, ])
  St = res$St
  r = res$r
  Uo = res$Uo
  Z = res$Z
  if (is.null(Z)) Z   = NA
  if (!anyNA(Z))  {
    Z = as.matrix(Z)
    kz = ncol(Z)
  } else kz  = 0
  Do = res$Do
  DFYflag = res$DFYflag
  SigmaS = (1:(n * m * S * n * m * S)) * 0
  dim(SigmaS) = c(n * m * S, n * m * S)
  VECM_domestic = list()
  VECM_domesticS = list()
  CRK_tst = list()
  T = nrow(Y)
  check = c(1:(S + 1)) * 0
  LH_P = (1:n) * 0
  LH_AIC = (1:n) * 0
  LH_BIC = (1:n) * 0
  ORLHP = (1:n) * 0
  ORAIC = (1:n) * 0
  ORBIC = (1:n) * 0
  Ct = Uo * 0

  WMATn = W %x% diag(m)
  FY = Y %*% WMATn
  Go = (1:(m * n * m * n * Pmax * S)) * 0
  dim(Go) = c(m * n, m * n, Pmax, S)
  Bo = (1:(m * m * Pmax * n * S)) * 0
  dim(Bo) = c(m, m, Pmax, n, S)
  Ao = Bo * 0
  if (type == "none"|type == "exog0" )
    Model = "I"
  if (type == "const"|type =="exog1" )
    Model = "III"
  for (i in 1:n) {
    pj = t(p[i, , ])
    pjOR = pj * 0
    pjOR[1, 1] = max(pj[, 1])
    pjOR[1, 2] = max(pj[, 2])
    y = Y[, ((i - 1) * m + 1):(i * m)]
    x = FY[, ((i - 1) * m + 1):(i * m)]

    HilfeYi  = Y[, ((i - 1) * m + 1):((i - 1) * m + m)]
    HilfeFYi = FY[, ((i - 1) * m + 1):((i - 1) * m + m)]
    if (is.null(colnames(HilfeYi))) colnames(HilfeYi)   = paste0(rep("Y",ncol(HilfeYi)),c(1:ncol(HilfeYi)))
    if (is.null(colnames(HilfeFYi))) colnames(HilfeFYi) = paste0(rep("FY",ncol(HilfeYi)),c(1:ncol(HilfeFYi)))
    y = HilfeYi
    x = HilfeFYi
    crk = m - r[i]
    Sti = 2 - St[, i]
    Pmaxi = max(p[i, 1, ])
    if (!anyNA(X)) {
      if (is.null(colnames(X))) colnames(X)   = paste0(rep("X",ncol(X)),c(1:ncol(X)))
      ORXi = X[,,i,1]; dim(ORXi) = dim(X)[1:2]
      MRXi = X[,,i,];  #dim(MRXi) = c(dim(X)[1:2],2)
    } else {
      MRXi = NA
      ORXi = NA
    }
    #ORtst <- MRCVECMest2(y, x, model = Model, type = "trace",P = pjOR, crk = crk, q = 0.95, Dxflag = 0,X=ORXi)
    ORtst <- MRCVECMest2(y,x, model = Model, type = "trace", P = pjOR, crk = crk, q = 0.95, Dxflag = DFYflag, X=ORXi, CZ=as.matrix(Z))
    nrowresidOR = nrow(ORtst[[2]]$residuals)
    UOR = ORtst[[2]]$residuals
    SigmaOR = t(UOR) %*% UOR/(T - P * m + crk - ncol(x) + (Model == "III"))
    ORLHP[i] = -(T * m/2) * log(2 * 3.1415) - (T * m/2) + T/2 * log(det(solve(SigmaOR)))
    ORAIC[i] = 2 * m * (m * pjOR[1, 1] - m + crk + (m + 1)/2) -  2 * ORLHP[i]
    ORBIC[i] = log(T) * m * (m * pjOR[1, 1] - m + crk + (m +  1)/2) - 2 * ORLHP[i]
    #tst <- MRCVECMestm(y, x, s = Sti, model = Model, type = "trace", P = pj, crk = crk, q = 0.95,X=MRXi)
    tst <- MRCVECMestm(y, x, s = Sti, model = Model, type = "trace", P = pj, crk = crk, q = 0.95,Dxflag = DFYflag, X=MRXi,CZ=as.matrix(Z))
    VECM_domestic[[i]]     = tst$VECM1
    VECM_domesticS[[i]]    = tst$VECM1S
    CRK_tst[[i]] = tst$erg
    nrowresid = nrow(tst$VECM1$residuals)
    Uo[(T - nrowresid + 1):T, ((i - 1) * m + 1):(i * m),
       1] = (2 - St[(T - nrowresid + 1):T, i]) * tst$VECM1$residuals
    Uo[(T - nrowresid + 1):T, ((i - 1) * m + 1):(i * m),
       2] = (St[(T - nrowresid + 1):T, i] - 1) * tst$VECM1$residuals
    sigma1 = t(Uo[(T - nrowresid + 1):T, ((i - 1) * m + 1):(i *
                                                              m), 1]) %*% Uo[(T - nrowresid + 1):T, ((i - 1) *
                                                                                                       m + 1):(i * m), 1]/(sum(2 - St) - (m * sum(p[i, 1:2,
                                                                                                                                                    1]) + m - crk))
    sigma2 = t(Uo[(T - nrowresid + 1):T, ((i - 1) * m + 1):(i *
                                                              m), 2]) %*% Uo[(T - nrowresid + 1):T, ((i - 1) *
                                                                                                       m + 1):(i * m), 2]/(sum(St - 1) - (m * sum(p[i, 1:2,
                                                                                                                                                    2]) + m - crk))
    LH_P[i] = -(T * m/2) * log(2 * 3.1415) - (T * m/2) + (sum(2 -
                                                                St))/2 * log(det(solve(sigma1))) + (sum(St - 1))/2 *
      log(det(solve(sigma2)))
    LH_AIC[i] = 2 * (m * (sum(p[i, 1:2, 1]) - m + crk) +
                       m * (m + 1)/2) + 2 * (m * (sum(p[i, 1:2, 2]) - m +
                                                    crk) + m * (m + 1)/2) - 2 * LH_P[i]
    LH_BIC[i] = log(sum(2 - St)) * (m * (sum(p[i, 1:2, 1]) -
                                           m + crk) + m * (m + 1)/2) + log(sum(St - 1)) * (m * (sum(p[i, 1:2, 2]) - m + crk) + m * (m + 1)/2) - 2 * LH_P[i]


    ConnectPoint = 2*crk+(sum(pj[,1:2])-4)*m
    Diff =  (dim(tst$VECM1[[1]])[1] -  (2*crk+(sum(pj[,1:2])-4)*m))/2

    if ( Diff>0 ) {
      ADD1 = 1:Diff
      ADD2 = Diff + (1:Diff)
      if ((!type=="none")&(kz>0)) {
        CLength = Diff-kz*(pj[1,1]-1)
        ADD1 = c(1:CLength,2*CLength+1:(kz*(pj[1,1]-1)))
        ADD2 = c(CLength+1:CLength,2*CLength+kz*(pj[1,1]-1)+1:(kz*(pj[1,1]-1)))
      }
    }   else {
      ADD1 = NULL
      ADD2 = NULL

    }


    AB1 = VECM2VAR(param = tst$VECM1[[1]][c(1:crk,2*crk +(1: (m * (pj[1,1] - 1 + pj[1, 2] - 1))),ConnectPoint+ADD1), ], beta = tst$betaS, q = c(crk,(pj[1, 1] - 1), (pj[1, 2] - 1)),kz=kz,Dxflag=DFYflag)
    AB2 = VECM2VAR(param = tst$VECM1[[1]][c(crk+(1:crk), (crk*2 + m * (pj[1, 1] - 1 + pj[1, 2] - 1)) + 1:(m * (pj[2,1] + pj[2, 2] - 2)),ConnectPoint+ADD2), ], beta = tst$betaS, q = c(crk,(pj[2, 1] - 1), (pj[2, 2] - 1)),kz=kz,Dxflag=DFYflag)


    #AB1 = VECM2VAR(param = tst$VECM1[[1]][c(1:crk,2*crk +(1: (m * (pj[1,1] - 1 + pj[1, 2] - 1)))), ], beta = tst$betaS, q = c(crk,(pj[1, 1] - 1), (pj[1, 2] - 1)),Dxflag=DFYflag)
    #AB2 = VECM2VAR(param = tst$VECM1[[1]][c(crk+(1:crk), (crk*2 + m * (pj[1, 1] - 1 + pj[1, 2] - 1)) + 1:(m * (pj[2,1] + pj[2, 2] - 2))), ], beta = tst$betaS, q = c(crk,(pj[2, 1] - 1), (pj[2, 2] - 1)),Dxflag=DFYflag)
    Bo[, , 1:pj[1, 1], i, 1] = AB1[[1]]
    Ao[, , 1:pj[1, 2], i, 1] = AB1[[2]]
    Bo[, , 1:pj[2, 1], i, 2] = AB2[[1]]
    Ao[, , 1:pj[2, 2], i, 2] = AB2[[2]]
    if (!anyNA(Do)) {
      Do[, , 1:pj[1 ,1], i, 1] = AB1[[4]]
      Do[, , 1:pj[1 ,1], i, 2] = AB2[[4]]
    }


    if (type == "const" | type == "exog1") {
      dd = (dim(tst$VECM1[[1]])[1]-(crk*2 + (sum(pj[, -3]) - 4) * m)-2*kz*(pj[1,1]-1))/2
      Co[, , i, 1] = t(tst$VECM1[[1]][crk*2 + (sum(pj[, -3]) - 4) * m + 1:dd, ])
      Co[, , i, 2] = t(tst$VECM1[[1]][crk*2 + (sum(pj[, -3]) - 4) * m + dd +1:dd,])
      Co[, , i, 1] = AB1[[3]]
      Co[, , i, 2] = AB2[[3]]

    }
    if (type == "exog0") {
      dd = (dim(tst$VECM1[[1]])[1]-(crk*2 + (sum(pj[, -3]) - 4) * m)-2*kz*(pj[1,1]-1))/2
      Co[,-1, i, 1] =  t(tst$VECM1[[1]][crk*2 + (sum(pj[, -3]) - 4) * m + 1:dd, ])
      Co[,-1, i, 2] =  t(tst$VECM1[[1]][crk*2 + (sum(pj[, -3]) - 4) * m + 1:dd, ])
      Co[,-1, i, 1] =  AB1[[3]]
      Co[,-1, i, 2] =  AB2[[3]]

    }
  }
  for (s in 1:S) Go[, , , s] = BoAoW2G(Bo[, , , , s], Ao[,
                                                         , , , s], W, m, n, Pmax)
  for (s in 1:S) Sigmao[, , s] = t(Uo[, , s]) %*% Uo[, , s]/sum(St[,
                                                                   i] == s)
  for (s in 1:S) {
    for (i in 1:n) {
      for (j in 1:n) Sigmao[(m * (i - 1) + 1):(i * m),
                            (m * (j - 1) + 1):(j * m), s] = t(Uo[, (m * (i -
                                                                           1) + 1):(i * m), s]) %*% Uo[, (m * (j - 1) +
                                                                                                            1):(j * m), s]/sum((St[, i] == s) * (St[, j] ==
                                                                                                                                                   s))
    }
  }
  for (s in 1:S) {
    for (ss in 1:S) {
      for (i in 1:n) {
        for (j in 1:n) SigmaS[((s - 1) * n * m + (i -
                                                    1) * m + 1):((s - 1) * n * m + i * m), ((ss -
                                                                                               1) * n * m + (j - 1) * m + 1):((ss - 1) * n *
                                                                                                                                m + j * m)] = t(Uo[, (m * (i - 1) + 1):(i *
                                                                                                                                                                          m), s]) %*% Uo[, (m * (j - 1) + 1):(j * m),
                                                                                                                                                                                         ss]/sum((St[, i] == s) * (St[, j] == ss))
      }
    }
  }
  res$Go = Go
  res$Sigmao = Sigmao
  res$Ao = Ao
  res$Bo = Bo
  res$Co = Co
  res$Do = Do
  res$GDC = GDC
  res$SigmaS = SigmaS
  res$resid = Uo
  res$VECM_domestic = VECM_domestic
  ORAIC  = sum(ORAIC)
  ORBIC  = sum(ORBIC)
  ORLHP  = sum(ORLHP)
  LH_BIC = sum(LH_BIC)
  LH_AIC = sum(LH_AIC)
  LH_P   = sum(LH_P)
  Summary = list(VECM_domesticS,CRK_tst,r_npo,ORLHP,ORAIC,ORBIC,LH_P,LH_AIC,LH_BIC)
  names(Summary) = c("Estimation_Result","CRK_Test","r_npo","One_Regeime_LH","One_Regime_AIC","One_Regime_BIC","LH_function_Value","AIC","BIC")
  res$Summary = Summary
  return(res)
}




#' Calculation of information criteria AIC and BIC for an estimated MRCIGVAR model
#'
#' @param  res  : an MRCIGVAR object generated from MRCIGVARData or estimated from MRCIGVARest.
#' @param  L_V  : a 2 components vector containing the maxima of the domestic lag length and the foreign lag length respectively.
#' @param  TH_V  : a vector containing possible threshold values .
#' @return      a matrix with different lag specifications,  threshold values and the corresponding information criteria.
#' @examples
#'
#' m = 2                  ## m: number of variables, n: number of countries
#' n = 4
#' p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' p = p[1:4,,]; p[,1:2,2] = 3
#'
#' TH = c(1:4)*0; dim(TH) = c(1,4)
#'
#' res_d <- MRCIGVARData(m=2,n=4,p=p,TH=TH,T=400,S=2, SESVI=c(1,3,5,7),r=rep(1,4))
#'
#' res_em = MRCIGVARest(res=res_d)
#'
#' L_v = c(3,3)
#' TH_v = c(0,0.1)
#' res_d$p
#'
#' MRCIGVARSelect = MRCIGVAR_Select(res=res_d,L_V=L_v,TH_V=TH_v)
#' dim(MRCIGVARSelect)
#' MRCIGVARSelect[which.min(MRCIGVARSelect[,28]),]
#'
#' @export
MRCIGVAR_Select = function(res,L_V=L_V,TH_V=TH_V)  {
  Criteria = matrix(0,n*(L_V[1]-1)*(L_V[2]-1)*2*length(TH_V),length(as.vector(res$p))+6)
  idx = 0
  res_dd   = res
  idx = 0
  n   = res$n

  for (i in 1:n )           	{
    for (l_d in 2: L_V[1] )   	{
      for (l_f in 2: L_V[2] )   	{
        for (s  in 1:2       )    	{
          for (l_th in 1: length(TH_V) ){
            idx = idx + 1
            res_dd$p[i,1,s] = l_d
            res_dd$p[i,2,s] = l_f
            res_dd$TH       = TH_V[l_th]

            res_s = MRCIGVARest(res_dd)
            Criteria[idx,] = c(c(as.vector(res_dd$p),TH_V[l_th],s),res_s$Summary$AIC,res_s$Summary$BIC,res_s$Summary$One_Regime_AIC,res_s$Summary$One_Regime_BIC)
            #colnames(Criteria) = c("l_d","l_f","AIC","BIC","ORAIC","ORBIC")
          }
        }
      }
    }
  }
  return(Criteria)
}


#' Regime specific impulse response functions of an MRCIGVAR(m,n,p,S) model
#'
#' This function calculates the regime specific impulse response functions of an estimated MRCIGVAR(n,p,S,r). Using G\[,,,s\] and Sigma\[,,s\] matrices of the estimated MRCIGVAR, this function can produce impulse response functions for any possible combinations of states.
#' @param res a MRCIGVAR object that can be an output of MRCIGVARData, MRCIGVARest, or MRCIGVARest.
#' @param nstep the length of impulse response function
#' @param comb a vector specify the concerted action in policy-simulation impulse response function
#' @param state an n-vector specifying the specific state for each country.
#' @param irf types of the impulse response function c("gen","chol","chol1","gen1","comb1")
#' @param runs the number of bootstrapping runs
#' @param conf A vector specifying the confidence intervals
#' @param NT The number of impulse response scenarios
#' @return an list of an (mn,mn,nstep,3,S) array of the impulse response functions and test statistics. In the impulse response array, columns representing the impulse and rows the responses.
#' @examples
#' m = 2
#' n = 3
#' p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' p = p[1:n,,]; p[,1,] = 3; p[,2,] = 2
#'
#' TH = c(1:n)*0; dim(TH) = c(1,n)
#' SESVI=rep(1,3,5)
#' r  = rep(1,n)
#'
#' ## case of n = 3, m = 2, S = 2    m: number of variables, n: number of countries
#'
#'
#' res_d <- MRCIGVARData(m=2,n=3,p=p,TH=TH,T=100,S=2, SESVI=c(1,3,5),r=rep(1,3),Ncommtrend=1)
#' max(abs(res_d$Y))
#' STAT(res_d$Go[,,,2])
#' STAT(res_d$Go[,,,1])
#' res_e  = MRCIGVARest(res=res_d)
#'
#' res_e$Summary
#'
#' IRF_CB = irf_MRCIGVAR_CB(res=res_e, nstep=10, comb=NA, state = c(2,1,1), irf = "gen1",runs=20,
#'  conf = c(0.05, 0.95), NT = 1)
#'
#' str(IRF_CB)
#'
#' IRF_g = IRF_graph(IRF_CB[[1]])
#'
#' @export
#'
irf_MRCIGVAR_CB <- function (res, nstep, comb, state = c(2, 1), irf = c("gen", "chol", "chol1", "gen1", "comb1"), runs = 200, conf = c(0.05, 0.95), NT = 1)
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
  Sigmao = res$Sigmao
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
  response <- array(0, dim = c(neq, nvar, nstep, length(conf) +
                                 1))
  sigmaNPDS = SigmaNPD(res, state)
  sigmaNPD1 = SigmaNPD(res, rep(1, n))
  sigmaNPD2 = SigmaNPD(res, rep(2, n))
  Sigmao_sim = Sigmao * 0
  Sigmao_sim[, , 1] = sigmaNPD1
  Sigmao_sim[, , 2] = sigmaNPD2
  response[, , , 1] <- irf_MRCIGVAR(res, nstep, comb, state, irf)
  responseR <- array(0, dim = c(neq, nvar, nstep, runs))
  Uo_run = array(0, c(T, n * m, S))
  for (i in 1:runs) {
    for (s in 1:S) Uo_run[, , s] = rnormSIGMA(T, sigmaNPDS)
    if (length(p) > 1) {
      res_run = MRCIGVARData(m, n, p, T, S, W, SESVI, TH, res$Go, Sigmao = Sigmao_sim, SV, type = "none")
      res_e = MRCIGVARest(res_run)
    }
    if (length(p) == 1) {
      res_run = MRGVARData(m, n, p, T, S, W, SESVI, TH, res$Go, Ao = NA, Bo = NA, Sigmao = NA, Uo = Uo_run, SV, type, res$Co, X)
      res_e = MRCIGVARest(res_run)
    }
    responseR[, , , i] <- irf_MRGVAR(res_e,state, nstep, comb, irf,G=NA,smat=NA,sigmaNPDS=SigmaNPD(res_e, state))
    GoColect[, , , , i] <- res_e$Go
    BoColect[, , , , , i] <- res_e$Bo
    AoColect[, , , , , i] <- res_e$Ao
    UoColect[, , , i] <- res_run$Uo
    YColect[, , i] <- res_run$Y
  }
  for (tt in 1:(nstep)) {
    for (i in 1:neq) {
      for (j in 1:nvar) {
        response[i, j, tt, -1] = stats::quantile(responseR[i, j, tt, ], conf)
      }
    }
  }
  return(list(response, GoColect, BoColect, AoColect, UoColect, YColect))
}


#' Regime specific impulse response functions of MRCIGVAR(n,p,S)
#'
#' This function calculates the regime specific impulse response functions of an estimated MRCIGVAR(n,p,S).
#' Using the estimated G\[,,,s\] and Sigma\[,,s\] matrices of the MRGVAR, this function calculated the regime speicfic impulse response functions.
#' @param res a list of estimated MRCIGVAR as output of MRCIGVARest or MRCIGVARest
#' @param nstep the length of impulse response function
#' @param comb a vector specify the concerted action in policy-simulation impulse response function
#' @param state an n vector specifying the speciic state for each country.
#' @param  irf  : types of the impulse response irf=c("gen","chol","chol1","gen1","comb1")
#' @param sigmaNPDS the state-dependent covariance matrix
#' @return a list containing the impulse response functions and the accumulated impulse response function, and the boostrap parameters as well.
#' @examples
#' m = 2
#' n = 3
#' p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' p = p[1:n,,]; p[,1,] = 3; p[,2,] = 2
#'
#' TH = c(1:n)*0; dim(TH) = c(1,n)
#' SESVI=rep(1,3,5)
#' r  = rep(1,n)
#'
#' ## case of n = 3, m = 2, S = 2    m: number of variables, n: number of countries
#'
#'
#' res_d <- MRCIGVARData(m=2,n=3,p=p,TH=TH,T=300,S=2, SESVI=c(1,3,5),r=rep(1,3),Ncommtrend=1)
#' max(abs(res_d$Y))
#' STAT(res_d$Go[,,,2])
#' STAT(res_d$Go[,,,1])
#' res_e  = MRCIGVARest(res=res_d)
#'
#' res_e$Summary
#'
#' IRF_CB = irf_MRCIGVAR_CB(res=res_e, nstep=10, comb=NA, state = c(2,1,1), irf = "gen1",
#' runs = 20, conf = c(0.05, 0.95), NT = 1)
#'
#' str(IRF_CB)
#'
#' IRF_g = IRF_graph(IRF_CB[[1]])
#' IRFF2 = irf_MRCIGVAR(res=res_e,nstep=10,comb=NA,state=c(2,2,2),irf='gen')
#' plot(IRFF2[1,1,],type="l")
#'
#'
#' @export
#'
irf_MRCIGVAR = function(res,nstep,comb,state=state,irf=c("gen","chol","chol1","gen1","comb1"),sigmaNPDS=NA) {
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
      if (anyNA(sigmaNPDS))	sigma = SigmaNPD(res,state)  else  sigma = sigmaNPDS
      response <- array(0,dim=c(neq,nvar,nstep));
      response <- irf_B_sigma(B,sigma,nstep,comb,irf=irf)
	return(response)
}


#' Generalized impulse response functions of MRCIGVAR(m,n,p,S) with regime migrations
#'
#' This function calculates the generalized impulse response functions of an estimated MRCIGVAR(m,n,p,S) for given a shock vector and initial values.
#'
#'                   GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK))
#'
#' @param  res   : an MRCIGVAR object containing the components of an output of MRCIGVARest.
#' @param  shock : an mn-vector containing the shocks as impulse.
#' @param  R     : the number of runs to integrate out the random effects in order to obtain the means (see equation above).
#' @param  nstep : the length of the responses
#' @param  Omega_hist : the initial values from which the simulation runs start.For Omega_hist=NA the most recent values are taken as the initial values. For Omega_hist=0, the initial values are zeros.
#' @param  resid_method : resid_method = c("resid", "parametric"), It generates the random residuals from residuals bootstrap or parametric bootstrap.
#' @return an (mn x mn x nstep+1) matrix of impulse response functions. The rows represent response the columns represent impulses.
#' @examples
#' m = 2
#' n = 3
#' p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' p = p[1:n,,]; p[,1,] = 3; p[,2,] = 2
#'
#' TH = c(1:n)*0; dim(TH) = c(1,n)
#' SESVI=rep(1,3,5)
#'r  = rep(1,n)
#'
#' ## case of n = 3, m = 2, S = 2    m: number of variables, n: number of countries
#' res_d <- MRCIGVARData(m=2,n=3,p=p,TH=TH,T=200,S=2, SESVI=c(1,3,5),r=rep(1,3),Ncommtrend=1)
#' max(abs(res_d$Y))
#' STAT(res_d$Go[,,,2])
#' STAT(res_d$Go[,,,1])
#' res_e  = MRCIGVARest(res=res_d)
#'
#' STAT(res_e$Go[,,,2])
#' STAT(res_e$Go[,,,1])
#'
#' plot(ts(res_d$Y))
#' #res_e$Summary
#'
#' if (!((max(Mod(STAT(res_e$Go[,,,1])))>1)|(max(Mod(STAT(res_e$Go[,,,2])))>1)) ) {
#'   GIRF <- girf_MRCIGVAR_RM(res=res_e,shock=c(1,1,1,1,1,1),R=100,nstep=10,Omega_hist=NA,
#'   resid_method="parametric")
#'   GIRF_CB <- girf_MRCIGVAR_RM_CB(res=res_e,shock=c(1,1,1,1,1,1),R=100,nstep=10,Omega_hist=NA,
#'   resid_method="parametric",conf=c(0.05,0.95),N=100)
#'   IRF_g = IRF_graph(GIRF_CB)
#' }
#'
#' @export
girf_MRCIGVAR_RM <- function(res,shock,R,nstep,Omega_hist=NA,resid_method) {
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

  if (anyNA(Omega_hist)) Omega_hist = res$Y[(dim(res$Y)[1]-P+1):dim(res$Y)[1],,drop=FALSE]
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
          for (s in 1:S)   residI[,,k,s] = rnormSIGMA(P+nstep+1,Matrix::nearPD(Sigmao[,,s])[[1]])
        }
      }
    }

    residR[[i]] <- residI
    for (s in 1:S) {residI[P+1,,,s] <-shock_mat}    # such that residI[P+1,,k,] corresponds to the k-th column shock in shock_mat for both regimes (budui dan ok)
    residS[[i]] = residI

    for (k in 1:DIMresid[2]) {
      YR[,,k] = MRCIGVARData(m=m,n=n,p=p,T=(P+nstep+1),S=S,W=W,SESVI=SESVI,TH=TH,Go=Go,Ao=Ao,Bo=Bo,Sigmao=Sigmao,Uo=residR[[i]][,,k,],SV=SV,type=type,Co=Co,X=X,Yo=Yo,d=d)$Y
      YS[,,k] = MRCIGVARData(m=m,n=n,p=p,T=(P+nstep+1),S=S,W=W,SESVI=SESVI,TH=TH,Go=Go,Ao=Ao,Bo=Bo,Sigmao=Sigmao,Uo=residS[[i]][,,k,],SV=SV,type=type,Co=Co,X=X,Yo=Yo,d=d)$Y
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
  return(GIRF)
}

#' Generalized impulse response functions of MRCIGVAR(m,p,n,S,W,TH) with regime migrations
#'
#' This function calculates the generalized impulse response functions of an estimated MRVAR(n,p,S) for given a shock vector and initial values.
#'
#'                   GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK))
#'
#' It also generates the bootstrapped confidence intervals.
#'
#' @param  res   : a MRCIGVAR object containing the components of an output of MRCIGVARest.
#' @param  shock : an mn-vector containing the shocks as impulse.
#' @param  R     : the number of runs to integrate out the random effects in order to obtain the means (see equation above).
#' @param  nstep : the length of the responses
#' @param  Omega_hist : the initial values from which the simulation runs of impulse and response functions start
#' @param  resid_method : resid_method = c("resid", "parametric"), It generate random residuals either from residuals bootstrap or parametric bootstrap.
#' @param  conf  : a two component vector containing the tail probabilities of the bootstrap confidence interval.
#' @param  N     : number of bootstrapping runs
#' @return an (n x n x nstep+1 x 3) array containing of impulse response functions with lower and upper confidence bonds. The rows represent response the columns represent impulses.
#' @examples
#' m = 2
#' n = 3
#' p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' p = p[1:n,,]; p[,1,] = 3; p[,2,] = 2
#'
#' TH = c(1:n)*0; dim(TH) = c(1,n)
#' SESVI=rep(1,3,5)
#'r  = rep(1,n)
#'
#' ## case of n = 3, m = 2, S = 2    m: number of variables, n: number of countries
#' res_d <- MRCIGVARData(m=2,n=3,p=p,TH=TH,T=200,S=2, SESVI=c(1,3,5),r=rep(1,3),Ncommtrend=1)
#' max(abs(res_d$Y))
#' STAT(res_d$Go[,,,2])
#' STAT(res_d$Go[,,,1])
#' res_e  = MRCIGVARest(res=res_d)
#'
#' STAT(res_e$Go[,,,2])
#' STAT(res_e$Go[,,,1])
#'
#' plot(ts(res_d$Y))
#' #res_e$Summary
#'
#' if (!((max(Mod(STAT(res_e$Go[,,,1])))>1)|(max(Mod(STAT(res_e$Go[,,,2])))>1)) ) {
#'   GIRF <- girf_MRCIGVAR_RM(res=res_e,shock=c(1,1,1,1,1,1),R=100,nstep=10,Omega_hist=NA,
#'   resid_method="parametric")
#'   GIRF_CB <- girf_MRCIGVAR_RM_CB(res=res_e,shock=c(1,1,1,1,1,1),R=100,nstep=10,Omega_hist=NA,
#'   resid_method="parametric",conf=c(0.05,0.95),N=100)
#'   IRF_g = IRF_graph(GIRF_CB)
#' }
#'
#' @export
girf_MRCIGVAR_RM_CB <- function(res,shock,R,nstep,Omega_hist=NA,resid_method="parametric",conf,N) {
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

  GIRFBd[,,,1]= girf_MRCIGVAR_RM(res,shock,R,nstep,Omega_hist,resid_method)

  for (i in 1:N) {
    #res_run = MRGVARData(m=m,n=n,p=p,T=T,S=S,W=W,SESVI=SESVI,TH=TH,Go=Go,Ao=Ao,Bo=Bo,Sigmao=Sigmao,Uo=NA,SV=SV,type=type,Co=Co,X=X,Yo=Yo,d=d)
    #res_e   = MRGVARest(res_run)
    res_run  = MRCIGVARDataR(res)
    if (length(colnames(res_run$Y))==0) colnames(res_run$Y) = paste("res_runY",1:ncol(res_run$Y),sep="")
    res_erun   = MRCIGVARest(res_run)
    RF3     = girf_MRCIGVAR_RM(res=res_erun,shock,R,nstep,Omega_hist,resid_method)
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



#' Data generating process of MRCIGVARDataR(res)
#'
#' This function will generate data from a MRCIGVAR object. It will generate enough data for estimation purpose.
#'
#' @param res     : an output of MRCIGVARest
#' @return	: an MRCIGVAR object.
#'
#' @export
MRCIGVARDataR <- function(res) {
  ### res_e : an estimated MRGVAR model that is an output of MRCIGVARest
  ### T     : number of observations
  ### Remarks: MRCIGVARDataR is used in bootstrapping to generate sufficient observations such that in the regime of rare occurrence also contains
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
  Sigmao = res$Sigmao
  type = res$type
  X = res$X
  mu = res$mu
  MINH = T
  repeat {
    Uo_run = array(0, c(T, n * m, S))
    for (s in 1:S) Uo_run[, , s] = rnormSIGMA(T, Matrix::nearPD(Sigmao[,,s])[[1]])       #### zhe ge you wenti
    res_run = MRCIGVARData(m, n, p, T, S, W, SESVI, TH, res$Go, Ao = NA, Bo = NA, Sigmao = NA, Uo = Uo_run, SV, type, Co, X)
    for (i in 1:n) { MINH = min(MINH,sum(diff(res_run$Y)[,SESVI[i]]>TH[1,i]),sum(diff(res_run$Y)[,SESVI[i]]<TH[1,i]))}
    if (MINH > (max(p)*(m*2)+2*m)) break
    else {T=T+T; print(c(MINH,T)); MINH=T}
  }
  return(res_run)
}




