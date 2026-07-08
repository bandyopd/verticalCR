

#' Title
#'
#' @param xtime 
#' @param Delta 
#' @param D 
#' @param X 
#' @param Z 
#' @param Bt 
#' @param W 
#' @param Q 
#' @param ibeta 
#' @param igamma 
#' @param ikappa 
#' @param iU 
#' @param izetah 
#' @param izetal 
#' @param irho 
#' @param h1 
#' @param h2 
#' @param itheta 
#' @param lr 
#' @param mite 
#' @param miter 
#' @param eps 
#'
#' @returns
#' @export
#'
#' @examples
tw_vertical <- function(xtime, Delta, D, X, Z, Bt, W, Q, ibeta, igamma,
                        ikappa, iU, izetah, izetal, irho, h1, h2, itheta,
                        lr = 1, mite = 20L, miter = 20L, eps = 1E-3,
                        rho_bound = 1.0, zeta_lr_pi = 1.0) {
  .Call(`_VerticalCR_tw_vertical`, xtime, Delta, D, X, Z, Bt, W, Q, ibeta,
        igamma, ikappa, iU, izetah, izetal, irho, h1, h2, itheta, lr, mite,
        miter, eps, rho_bound, zeta_lr_pi)
}

