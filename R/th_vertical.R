

#' Title: estimating function for three causes vertical model
#'
#' @param xtime:observed time
#' @param Delta:censoring indicator 
#' @param D1:cause 2 indicator 
#' @param D2:cause 3 inicator 
#' @param X:covariates of total hazards model 
#' @param Z:covariates of relative hazards model 
#' @param Bt:b-spline basis for estimating function in relative hazards model 
#' @param W:threshold variable
#' @param Q:grouping indicator
#' @param ibeta:initial values of covariates cofficients in total hazards model
#' @param igamma1:initial values of covariates cofficients in cause 2 of relative hazards model
#' @param igamma2:initial values of covariates cofficients in cause 3 of relative hazards model 
#' @param ikappa1:initial values of b-spline cofficients  in cause 2 of relative hazards model 
#' @param ikappa2:initial values of b-spline cofficients  in cause 3 of relative hazards model 
#' @param iU:initial values of frailty 
#' @param izetah:initial value of threshold in total hazards model 
#' @param izetal:initial value of threshold in relative hazards model 
#' @param irho:initial value of frailty cofficients in relative hazards model 
#' @param h1:bandwidth of  the S-shaped approximation function in threshold estimation of total hazards model
#' @param h2:bandwidth of  the S-shaped approximation function in threshold estimation of relative hazards model
#' @param itheta:initial value of variance estimation of frailties
#' @param lr:step size reduction multiple
#' @param mite:The maximum number of iterations for coefficient estimation and frailties estimation 
#' @param miter:The maximum number of iterations for variance estimation  and coefficient estimation of frailties 
#' @param eps: 
#'
#' @returns
#' @export
#'
#' @examples
th_vertical <- function(xtime, Delta, D1, D2, X, Z, Bt, W, Q, ibeta, igamma1, igamma2, ikappa1, ikappa2, iU, izetah, izetal, irho, h1, h2, itheta, lr = 1, mite = 20L, miter = 20L, eps = 1E-3) {
  .Call(`_VerticalCR_th_vertical`, xtime, Delta, D1, D2, X, Z, Bt, W, Q, ibeta, igamma1, igamma2, ikappa1, ikappa2, iU, izetah, izetal, irho, h1, h2, itheta, lr, mite, miter, eps)
}