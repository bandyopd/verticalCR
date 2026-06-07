// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//[[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

arma::mat safe_inv(const arma::mat& x, const char* name) {
  arma::mat out;
  arma::mat xsym = 0.5 * (x + x.t());
  if (!arma::inv_sympd(out, xsym)) {
    if (!arma::inv(out, x)) {
      Rcpp::warning("%s is singular or numerically unstable; using Moore-Penrose pseudo-inverse. "
                    "Standard errors may be unreliable.", name);
      out = arma::pinv(x);
    }
  }
  return out;
}

double safe_exp(double x) {
  if (x > 35.0) return std::exp(35.0);
  if (x < -35.0) return std::exp(-35.0);
  return std::exp(x);
}

arma::vec damp_step(const arma::vec& step, double max_abs_step = 1.0) {
  double m = arma::max(arma::abs(step));
  if (m > max_abs_step && m > 0.0) {
    return step * (max_abs_step / m);
  }
  return step;
}

void check_finite_vec(const arma::vec& x, const char* name) {
  if (!x.is_finite()) Rcpp::stop("%s contains NA, NaN, or Inf.", name);
}

void check_finite_mat(const arma::mat& x, const char* name) {
  if (!x.is_finite()) Rcpp::stop("%s contains NA, NaN, or Inf.", name);
}

void validate_tw_vertical_inputs(const arma::vec& xtime,
                                 const arma::vec& Delta,
                                 const arma::vec& D,
                                 const arma::mat& X,
                                 const arma::mat& Z,
                                 const arma::mat& Bt,
                                 const arma::vec& W,
                                 const arma::mat& Q,
                                 const arma::vec& ibeta,
                                 const arma::vec& igamma,
                                 const arma::vec& ikappa,
                                 const arma::vec& iU,
                                 const arma::vec& izetah,
                                 const arma::vec& izetal,
                                 const arma::vec& irho,
                                 double h1,
                                 double h2,
                                 double itheta,
                                 int mite,
                                 int miter,
                                 double eps,
                                 double rho_bound) {
  const arma::uword n = xtime.n_elem;
  if (n == 0) Rcpp::stop("xtime must contain at least one observation.");
  if (Delta.n_elem != n || D.n_elem != n || W.n_elem != n) {
    Rcpp::stop("xtime, Delta, D, and W must have the same length.");
  }
  if (X.n_rows != n || Z.n_rows != n || Bt.n_rows != n || Q.n_rows != n) {
    Rcpp::stop("X, Z, Bt, and Q must have the same number of rows as xtime.");
  }
  if (Q.n_cols == 0) {
    Rcpp::stop("Q must have at least one column. Use N - 1 columns for sum-to-zero frailties.");
  }
  if (ibeta.n_elem != 2 * X.n_cols + 1) {
    Rcpp::stop("ibeta must have length 2 * ncol(X) + 1.");
  }
  if (igamma.n_elem != 2 * Z.n_cols + 1) {
    Rcpp::stop("igamma must have length 2 * ncol(Z) + 1.");
  }
  if (ikappa.n_elem != Bt.n_cols) {
    Rcpp::stop("ikappa must have length ncol(Bt).");
  }
  if (iU.n_elem != Q.n_cols) {
    Rcpp::stop("iU must have length ncol(Q).");
  }
  if (izetah.n_elem != 1 || izetal.n_elem != 1) {
    Rcpp::stop("izetah and izetal must be length-one numeric vectors.");
  }
  if (irho.n_elem != 1) Rcpp::stop("irho must have length 1 for the two-cause model.");
  if (!(h1 > 0.0) || !(h2 > 0.0)) Rcpp::stop("h1 and h2 must be positive.");
  if (!(itheta > 0.0)) Rcpp::stop("itheta must be positive.");
  if (mite < 1 || miter < 1) Rcpp::stop("mite and miter must be positive integers.");
  if (!(eps > 0.0)) Rcpp::stop("eps must be positive.");
  if (!R_finite(rho_bound) || rho_bound <= 0) Rcpp::stop("rho_bound must be a positive finite number.");

  check_finite_vec(xtime, "xtime");
  check_finite_vec(Delta, "Delta");
  check_finite_vec(D, "D");
  check_finite_vec(W, "W");
  check_finite_mat(X, "X");
  check_finite_mat(Z, "Z");
  check_finite_mat(Bt, "Bt");
  check_finite_mat(Q, "Q");
  check_finite_vec(ibeta, "ibeta");
  check_finite_vec(igamma, "igamma");
  check_finite_vec(ikappa, "ikappa");
  check_finite_vec(iU, "iU");
  check_finite_vec(irho, "irho");

  if (arma::any(xtime < 0)) Rcpp::stop("xtime must be non-negative.");
  if (arma::any(Delta < 0) || arma::any(Delta > 1) ||
      arma::any(D < 0) || arma::any(D > 1)) {
    Rcpp::stop("Delta and D must be 0/1 indicators.");
  }
  if (arma::any(D > Delta + 1e-8)) {
    Rcpp::stop("D must be zero for censored observations; reference-cause events use Delta=1 and D=0.");
  }
}

// [[Rcpp::export]]
List tw_updateb(arma::vec xtime,
                arma::vec Delta,
                arma::mat X,
                arma::mat Q,
                arma::vec ibeta,
                arma::vec iu,
                arma::vec izetah,
                arma::vec phih,
                arma::vec phih2,
                double lr);

// [[Rcpp::export]]
List tw_updateg(arma::vec Delta,
                arma::vec D,
                arma::mat Q,
                arma::mat ZB,
                arma::vec iu,
                arma::vec igammak,
                arma::vec phil,
                arma::vec izetal,
                arma::vec irho,
                double zeta_lr_pi,
                double lr);

// [[Rcpp::export]]
List tw_updateu(arma::vec xtime,
                arma::vec Delta,
                arma::vec D,
                arma::mat X,
                arma::mat Q,
                arma::mat ZB,
                arma::vec ibeta,
                arma::vec iu,
                arma::vec igammak,
                arma::vec irho,
                double itheta,
                double lr);

// [[Rcpp::export]]
List tw_verticalite2(arma::vec Delta,
                     arma::vec D,
                     arma::vec ZBg,
                     arma::vec qu,
                     arma::vec irho,
                     double lr,
                     int mite,
                     double eps,
                     double rho_bound
);
// [[Rcpp::export]]
List tw_verticalite3(arma::vec xtime,
                     arma::vec Delta,
                     arma::vec D,
                     arma::mat X,
                     arma::mat ZB,
                     arma::vec qu,
                     arma::mat Q,
                     arma::vec ibeta,
                     arma::vec iU,
                     arma::vec igammak,
                     arma::vec phih,
                     arma::vec phih2,
                     arma::vec phil,
                     arma::vec phil2,
                     arma::mat ZBz,
                     arma::vec irho,
                     double itheta);

// [[Rcpp::export]]
List tw_vertical(arma::vec xtime,
                 arma::vec Delta,
                 arma::vec D,
                 arma::mat X,
                 arma::mat Z,
                 arma::mat Bt,
                 arma::vec W,
                 arma::mat Q,
                 arma::vec ibeta,
                 arma::vec igamma,
                 arma::vec ikappa,
                 arma::vec iU,
                 arma::vec izetah,
                 arma::vec izetal,
                 arma::vec irho,
                 double h1,
                 double h2,
                 double itheta,
                 double lr=1,
                 int mite=20,
                 int miter=20,
                 double eps=1E-3,
                 double rho_bound=1.0,
                 double zeta_lr_pi=1.0
){
  validate_tw_vertical_inputs(xtime, Delta, D, X, Z, Bt, W, Q,
                              ibeta, igamma, ikappa, iU, izetah, izetal,
                              irho, h1, h2, itheta, mite, miter, eps,
                              rho_bound);
  if (!R_finite(zeta_lr_pi) || zeta_lr_pi <= 0.0) {
    Rcpp::stop("zeta_lr_pi must be a positive finite number.");
  }
  int i;
  int nuser= xtime.size();//观测时间数目
  int ncolX=X.n_cols;//特定原因生存模型协变量维度
  int ncolZ=Z.n_cols;//相对风险模型协变量维度
  int ncolBt=Bt.n_cols;//相对风险模型非参数部分B样条维度
  int ncolQ=Q.n_cols;//随机效应维度
  int ncenter=ncolQ+1;//sum-to-zero contrast uses N-1 coefficients for N centers
  double cgs=100,cgsr=100;//迭代差初始值
  int ite=1,iter=1;//迭代初始值
  arma::vec iUc(ncolQ);//设置随机效应变量
  arma::vec se_rho(1);//设置1维随机效应系数变量
  arma::vec ZBg(nuser);//设置汇总相对风险中协变量和非参部分的变量
 
  arma::vec qu(nuser);//设置随机效应值的变量
  arma::vec phih(nuser);//特定原因生存模型S曲线-正态密度的变量
  arma::vec phil(nuser);//相对风险模型S曲线-正态密度的变量
  arma::vec phih2(nuser);//second derivative of total threshold predictor wrt zeta_lambda
  arma::vec phil2(nuser);//second derivative of relative threshold predictor wrt zeta_pi
  arma::vec dPhil(nuser);//first derivative of smoothed threshold function wrt zeta_pi
  
  
  arma::vec Phih(nuser);//特定原因生存模型S曲线-正态分布的变量
  arma::vec Phil(nuser);//相对风险模型S曲线-正态分布的变量
  arma::vec igammak(2*ncolZ+1+ncolBt);//相对风险模型case1随机效应、回归系数和B样条系数变量
  
  
  double etheta;//随机效应方差
  arma::mat Xt(nuser,2*ncolX+1);//特定原因生存模型协变量和分组变量的整合
  arma::mat ZB(nuser,2*ncolZ+1+ncolBt);//相对风险模型协变量、分组变量和B样条的整合
  arma::mat ZBz(nuser,2*ncolZ+1+ncolBt, arma::fill::zeros);//d ZB / d zeta_pi
 arma::mat T33(ncolQ,ncolQ);//随机效应协方差
  
  arma::vec ur(ncolQ+2*ncolX+1+2*ncolZ+1+2+ncolBt);//目标函数的得分函数
  arma::vec sxtime=sort(xtime);//生存时间排序
  Delta=Delta(sort_index(xtime));//将删失指标根据生存时间排序
  D=D(sort_index(xtime));//将病因指标根据生存时间排序
  
  W=W(sort_index(xtime));//将阈值变量根据生存时间排序
  X=X.rows(sort_index(xtime));//将特定病因生存模型协变量根据生存时间排序
  Z=Z.rows(sort_index(xtime));//将相对风险模型协变量根据生存时间排序
  Q=Q.rows(sort_index(xtime));//将分组变量根据生存时间排序
  Bt=Bt.rows(sort_index(xtime));//将B样条根据生存时间排序
  iter=1;
  // Xt=join_rows(X,ones(nuser),X);
  // ZB=join_rows(join_rows(Z,ones(nuser),Z),Bt);
  igammak=join_cols(igamma,ikappa);
  
  while(iter<miter && cgsr>eps){
    ite=1;
    cgs=100;
    while(ite<mite && cgs>eps){
      //计算特定病因生存模型关于阈值的一阶导
      Xt=join_rows(X,ones(nuser),X);
      ZB=join_rows(join_rows(Z,ones(nuser),Z),Bt);
      phih=normpdf((W-ones(nuser)*izetah)/h1)*(-1/h1)%(ones(nuser)*ibeta(ncolX)+X*ibeta.subvec(ncolX+1,ncolX*2));
      phih2=-((W-ones(nuser)*izetah)/h1)%normpdf((W-ones(nuser)*izetah)/h1)/(h1*h1)%
        (ones(nuser)*ibeta(ncolX)+X*ibeta.subvec(ncolX+1,ncolX*2));
      //计算相对风险模型关于阈值的一阶导
      dPhil=normpdf((W-ones(nuser)*izetal)/h2)*(-1/h2);
      phil=dPhil%(ones(nuser)*igammak[ncolZ]+Z*igammak.subvec(ncolZ+1,ncolZ*2));
      phil2=-((W-ones(nuser)*izetal)/h2)%normpdf((W-ones(nuser)*izetal)/h2)/(h2*h2)%
        (ones(nuser)*igammak[ncolZ]+Z*igammak.subvec(ncolZ+1,ncolZ*2));
      
      //计算特定病因生存模型阈值部分S曲线-正态分布取值
      Phih=normcdf((W-ones(nuser)*izetah)/h1);
      //计算相对风险模型阈值部分S曲线-正态分布取值
      Phil=normcdf((W-ones(nuser)*izetal)/h2);
      ZBz.zeros();
      //合并特定病因生存模型中的分组变量和协变量
      
      for(i =0;i<=ncolX;i++ ){
        Xt.col(i+ncolX)=Xt.col(i+ncolX)%Phih;}
      //合并相对风险模型中的分组变量、协变量和B样条

      for(i =0;i<=ncolZ;i++){
        ZBz.col(ncolZ+i)=ZB.col(ncolZ+i)%dPhil;
        ZB.col(ncolZ+i)=ZB.col(ncolZ+i)%Phil;
      }
      //合并相对风险模型中的随机效应、回归系数和B样条系数


      List Lu=tw_updateu(sxtime,
                         Delta,
                         D,
                         Xt,
                         Q,
                         ZB,
                         ibeta,
                         iU,
                         igammak,
                         irho,
                         itheta,
                         lr);
      //计算给定参数下目标函数的一阶向量和二阶矩阵


      arma::vec u_temp = Lu["u_new"];

      List Lb= tw_updateb( sxtime,
                           Delta,
                           Xt,
                           Q,
                           ibeta,
                           u_temp,
                           izetah,
                           phih,
                           phih2,
                           lr
      );
      arma::vec beta_temp = Lb["beta_new"];
      arma::vec zetah_temp = Lb["zetah_new"];

      List Lg= tw_updateg( Delta,
                           D,
                           Q,
                           ZB,
                           u_temp,
                           igammak,
                           phil,
                           izetal,
                           irho,
                           zeta_lr_pi,
                           lr);
      arma::vec gammak_temp = Lg["gammak_new"];
      arma::vec zetal_temp = Lg["zetal_new"];
      //arma::mat imatg=Lg["imatg"];
      cgs=max(abs(ibeta-beta_temp));

      iU=u_temp.clamp(-5,5);
      ibeta=beta_temp.clamp(-5,5);
      izetah=zetah_temp.clamp(0.1,0.9);
      igammak=gammak_temp.clamp(-5,5);
      izetal=zetal_temp.clamp(0.1,0.9);

      ite++;
    }
  //   
    Xt=join_rows(X,ones(nuser),X);
    ZB=join_rows(join_rows(Z,ones(nuser),Z),Bt);
    phih=normpdf((W-ones(nuser)*izetah)/h1)*(-1/h1)%(ones(nuser)*ibeta(ncolX)+X*ibeta.subvec(ncolX+1,ncolX*2));
    phih2=-((W-ones(nuser)*izetah)/h1)%normpdf((W-ones(nuser)*izetah)/h1)/(h1*h1)%
      (ones(nuser)*ibeta(ncolX)+X*ibeta.subvec(ncolX+1,ncolX*2));
    dPhil=normpdf((W-ones(nuser)*izetal)/h2)*(-1/h2);
    phil=dPhil%(ones(nuser)*igammak[ncolZ]+Z*igammak.subvec(ncolZ+1,ncolZ*2));
    phil2=-((W-ones(nuser)*izetal)/h2)%normpdf((W-ones(nuser)*izetal)/h2)/(h2*h2)%
      (ones(nuser)*igammak[ncolZ]+Z*igammak.subvec(ncolZ+1,ncolZ*2));
    Phih=normcdf((W-ones(nuser)*izetah)/h1);
    Phil=normcdf((W-ones(nuser)*izetal)/h2);
    for(i =0;i<=ncolX;i++ ){
      Xt.col(i+ncolX)=Xt.col(i+ncolX)%Phih;}
    ZBz.zeros();
    for(i =0;i<=ncolZ;i++){
      ZBz.col(ncolZ+i)=ZB.col(ncolZ+i)%dPhil;
      ZB.col(ncolZ+i)=ZB.col(ncolZ+i)%Phil;
    }
    ZBg= ZB*igammak;

    qu= Q*iU;
    List R2=tw_verticalite2(Delta,D,ZBg,qu,irho,lr,mite,eps,rho_bound);
    arma::vec rho_temp=R2["est"];
    List L4= tw_verticalite3(sxtime, Delta, D, Xt,ZB,qu, Q,ibeta,iU, igammak,
                             phih,
                             phih2,
                             phil,
                             phil2,
                             ZBz,
                             rho_temp,
                             itheta);
    arma::mat imat_temp=L4["imatri"];

  T33=imat_temp.submat(0,0,ncolQ-1,ncolQ-1);
  etheta=sum(T33.diag(0)+iU%iU)/(ncenter);
  cgsr=max(abs(rho_temp-irho));

  irho=rho_temp.clamp(-rho_bound,rho_bound);
  itheta=etheta;
  iter++;
  }

  arma::mat T33c=T33*T33;
  List R4= tw_verticalite3(sxtime, Delta, D, Xt,ZB,qu, Q,ibeta,iU, igammak,
                           phih,
                           phih2,
                           phil,
                           phil2,
                           ZBz,
                           irho,
                           itheta);
  arma::mat info_base=R4["imatr"];
  arma::mat info_full(info_base.n_rows + 1, info_base.n_cols + 1, arma::fill::zeros);
  info_full.submat(0, 0, info_base.n_rows - 1, info_base.n_cols - 1) = info_base;
  for (i=0; i<ncolQ; i++) {
    info_full(i, info_base.n_cols) = -iU[i] / (itheta * itheta);
    info_full(info_base.n_cols, i) = info_full(i, info_base.n_cols);
  }
  info_full(info_base.n_cols, info_base.n_cols) = ncenter / (2.0 * itheta * itheta);
  arma::mat temp4=safe_inv(info_full, "Full information matrix including theta");
  arma::vec info_diag=temp4.diag(0);
  if (arma::any(info_diag < -1e-8)) {
    Rcpp::warning("The inverse information matrix has negative diagonal entries; "
                  "reported standard errors use absolute values and may be unreliable.");
  }
  arma::vec se_para=sqrt(abs(info_diag));
  T33=info_base.submat(0,0,ncolQ-1,ncolQ-1);
  T33c=T33*T33;
  double se_theta_hl=sqrt(2*itheta*itheta/(ncenter-2*sum(T33.diag(0))/itheta+sum(T33c.diag(0))/(itheta*itheta)));
  double se_theta=se_para[info_base.n_cols];

  iUc=iU;



  bool final_converged = (cgsr <= 1e-3 && cgs <= eps);
  if (!final_converged) {
    Rcpp::warning("Algorithm stopped before full convergence; inspect last_inner_change and last_outer_change.");
  }
  if (itheta <= 1.01e-4) {
    Rcpp::warning("theta is close to its lower bound; frailty variance may be weakly identified.");
  }
  if (std::abs(irho(0)) >= rho_bound - 1e-6) {
    Rcpp::warning("rho reached its stability bound; relative-cause frailty association may be weakly identified.");
  }
  if (izetah(0) <= 0.101 || izetah(0) >= 0.899 || izetal(0) <= 0.101 || izetal(0) >= 0.899) {
    Rcpp::warning("A threshold estimate is near the search boundary; consider checking W support or using a wider interior search interval.");
  }

  List R=List::create(Named("beta")=ibeta,
                      Named("se_beta")=se_para.subvec(ncolQ,ncolQ+2*ncolX),
                      Named("zetah")=izetah,
                      Named("se_zetah")=se_para(ncolQ+2*ncolX+1),
                      Named("gamma")=igammak.subvec(0,(2*ncolZ)),
                      Named("se_gamma")=se_para.subvec(ncolQ+2*ncolX+2,ncolQ+2*ncolX+1+2*ncolZ+1),
                      Named("kappa")=igammak.subvec((2*ncolZ+1),(2*ncolZ+ncolBt)),
                      Named("se_kappa")=se_para.subvec(ncolQ+2*ncolX+2+2*ncolZ+1,ncolQ+2*ncolX+1+2*ncolZ+1+ncolBt),
                      Named("zetal")=izetal,
                      Named("se_zetal")=se_para(ncolQ+2*ncolX+1+2*ncolZ+1+ncolBt+1),
                      Named("rho")=irho(0),
                      Named("se_rho1")=se_para(ncolQ+2*ncolX+1+2*ncolZ+1+ncolBt+1+1),
                      Named("theta")=itheta,
                      Named("se_theta")=se_theta,
                      Named("se_theta_hl")=se_theta_hl,
                      Named("theta_in_joint_information")=true,
                      Named("v")=iU,
                      Named("U")=iUc,
                      Named("U_se")=se_para.subvec(0,ncolQ-1),
                      Named("iter")=iter,
                      Named("converged")=final_converged,
                      Named("last_inner_change")=cgs,
                      Named("last_outer_change")=cgsr,
                      Named("parameter_layout")=List::create(
                        Named("beta")="c(beta1, beta2_threshold_shift, beta3_interactions)",
                        Named("gamma")="c(gamma1, gamma2_threshold_shift, gamma3_interactions)",
                        Named("cause_coding")="D is the non-reference cause; Delta=1,D=0 is the reference cause"
                      )
  );


  // List R=List::create(Named("beta")=1);
  return R;
}
//计算给定参数下目标函数的一阶得分向量和二阶信息矩阵

List tw_updateb(arma::vec xtime,
                arma::vec Delta,
                arma::mat X,
                arma::mat Q,
                arma::vec ibeta,
                arma::vec iu,
                arma::vec izetah,
                arma::vec phih,
                arma::vec phih2,
                double lr){
  int nuser= xtime.size();
  int person,i,j;
  int ncolX=X.n_cols;
  int lenU=Q.n_cols;
  arma::mat cmatb=zeros(ncolX+1,ncolX+1);
  arma::mat imatb=zeros(ncolX+1,ncolX+1);
  arma::vec ub=zeros(ncolX+1);
  arma::vec ab=zeros(ncolX+1);
  arma::vec betaz_new(ncolX+1);
  
  
  double denom=0;
  double risk=0;
  double dtime;
  
  int nrisk=0;
  double xbetau;
  double temp2;
  double ab2=0;
  double dead_phih2=0;
  int ndead;
  
  
  for (person=nuser-1; person>=0; ) {
    dtime = xtime[person];
    ndead =0; /*number of deaths at this time point */
    dead_phih2=0;
    
    
    
    while(person >=0 && xtime[person]==dtime) {
      /* walk through the this set of tied times */
      nrisk++;
      xbetau = 0;    /* form the term beta*x (vector mult) */
      
      for (i=0; i<ncolX; i++)
        xbetau +=  ibeta[i]*X(person,i);
      for (i=0; i<lenU; i++)
        xbetau += Q(person,i)*iu[i];
      risk = safe_exp(xbetau);
      if (Delta[person] ==1) {
        ndead++;
        
        
        for (i=0; i<ncolX; i++) {
          ub[i] += X(person,i);
          
        }
        ub[ncolX]+=phih[person];
        dead_phih2+=phih2[person];
      }
      
      denom += risk;
      /* a contains weighted sums of x, cmat sums of squares */
      
      for (i=0; i<ncolX; i++) {
        ab[i] += risk*X(person,i);
        
        for (j=0; j<=i; j++)
          cmatb(j,i) += risk*X(person,i)*X(person,j);
      }
      ab[ncolX]+=risk*phih[person];
      ab2+=risk*phih2[person];
      for (j=0; j<ncolX; j++)
        cmatb(j,ncolX) += risk*(phih[person]*X(person,j));
      
      cmatb(ncolX,ncolX)+=risk*(phih[person]*phih[person]);
      
      
      
      person--;
    }
    if(ndead>0){
      for (i=0; i<ncolX; i++) {
        temp2= ab[i]/ denom;  /* mean */
      ub[i] -=   ndead*temp2;
      for (j=0; j<=i; j++)
        imatb(j,i) += ndead*(cmatb(j,i) - temp2*ab[j])/denom;
      }
      
      temp2= ab[ncolX]/ denom;
      for(j=0;j<ncolX;j++)
        imatb(j,ncolX)+=ndead*(cmatb(j,ncolX) - temp2*ab[j])/denom;
      
      ub[ncolX]-=ndead*temp2;
      
      imatb(ncolX,ncolX)+= ndead*(cmatb(ncolX,ncolX) - temp2*ab[ncolX])/denom;
      
      
    }
  }
  
  
  
  for(i=0;i<(ncolX+1);i++){
    for(j=0;j<i;j++){
      imatb(i,j)=imatb(j,i);
    }
  }
  arma::vec beta_step = safe_inv(imatb, "Information matrix for beta/zetah update") * ub;
  beta_step = damp_step(beta_step, 0.5);
  betaz_new=join_cols(ibeta,izetah)+beta_step*lr;
  
  
  List L=List::create(Named("beta_new")=betaz_new.subvec(0,(ncolX-1)),Named("zetah_new")=betaz_new[ncolX],Named("imatb")=imatb,Named("ub")=ub);
  
  return L;
}
List tw_updateg(arma::vec Delta,
                arma::vec D,
                arma::mat Q,
                arma::mat ZB,
                arma::vec iu,
                arma::vec igammak,
                arma::vec phil,
                arma::vec izetal,
                arma::vec irho,
                double zeta_lr_pi,
                double lr){
  int nuser= Delta.size();
  int person,i,j;
  int ncolZ=ZB.n_cols;
  int lenU=Q.n_cols;
  arma::mat imatg=zeros(ncolZ+1,ncolZ+1);
  arma::vec ug=zeros(ncolZ+1);
  arma::vec gamkapz_new(ncolZ+1);
  
  double ezgamma=0;
  
  double zgamma;
  
  for (person=nuser-1; person>=0; ) {
    /* walk through the this set of tied times */
    
    zgamma=0;
    for(i=0;i<ncolZ;i++){
      zgamma+=igammak[i]*ZB(person,i);
      
    }
    
    for(i=0;i<lenU;i++){
      zgamma+=Q(person,i)*irho[0]*iu[i];
      
    }
    
    ezgamma=safe_exp(zgamma);
    for(i=0; i<(ncolZ);i++){
      ug[i]+=D[person]*ZB(person,i)-Delta[person]*ezgamma*ZB(person,i)/(1+ezgamma);
       for(j=0; j<=i;j++){
        imatg(j,i)+=Delta[person]*ezgamma*ZB(person,i)*ZB(person,j)/pow(1+ezgamma,2);
             }
      }
    
    ug[ncolZ]+=D[person]*phil[person]-Delta[person]*(ezgamma*phil[person])/(1+ezgamma);
    for(i=0; i<(ncolZ);i++){
      
      imatg(i,ncolZ)+=Delta[person]*(ezgamma*phil[person]*ZB(person,i))/pow(1+ezgamma,2);
      
    }
    
    imatg(ncolZ,ncolZ)+=Delta[person]*((ezgamma*phil[person]*phil[person])/pow(1+ezgamma,2));
    
    person--;
  }
  for(i=0;i<(ncolZ+1);i++){
    for(j=0;j<i;j++){
      imatg(i,j)=imatg(j,i);
    }
  }
  arma::vec gamma_step = safe_inv(imatg, "Information matrix for gamma/zetal update") * ug;
  gamma_step = damp_step(gamma_step, 0.5);
  gamkapz_new=join_cols(igammak,izetal);
  gamkapz_new.subvec(0, ncolZ - 1) += gamma_step.subvec(0, ncolZ - 1) * lr;
  gamkapz_new[ncolZ] += gamma_step[ncolZ] * lr * zeta_lr_pi;
  
  
  
  List L=List::create(Named("gammak_new")=gamkapz_new.subvec(0,(ncolZ-1)),Named("zetal_new")=gamkapz_new[ncolZ],Named("imatg")=imatg);
  //  List L=List::create(Named("imatg")=imatg);
  return L;
}

List tw_updateu(arma::vec xtime,
                arma::vec Delta,
                arma::vec D,
                arma::mat X,
                arma::mat Q,
                arma::mat ZB,
                arma::vec ibeta,
                arma::vec iu,
                arma::vec igammak,
                arma::vec irho,
                double itheta,
                double lr){
  int nuser= xtime.size();
  int person,i,j;
  int ncolX=X.n_cols;
  int ncolZ=ZB.n_cols;
  int lenU=Q.n_cols;
  arma::mat cmatu=zeros(lenU,lenU);
  arma::mat imatu=zeros(lenU,lenU);
  arma::vec uu=zeros(lenU);
  arma::vec au=zeros(lenU);
  arma::vec u_new(lenU);
  
  
  double denom=0;
  double risk=0;
  double ezgamma;
  double dtime;
  
  int nrisk=0;
  double xbetau;
  double zgamma;
  double temp2;
  int ndead;
  
  
  for (person=nuser-1; person>=0; ) {
    dtime = xtime[person];
    ndead =0; /*number of deaths at this time point */
    
    
    
    while(person >=0 && xtime[person]==dtime) {
      /* walk through the this set of tied times */
      nrisk++;
      xbetau = 0;    /* form the term beta*x (vector mult) */
      
      for (i=0; i<ncolX; i++)
        xbetau +=  ibeta[i]*X(person,i);
      for (i=0; i<lenU; i++)
        xbetau += Q(person,i)*iu[i];
      risk = safe_exp(xbetau);
      if (Delta[person] ==1) {
        ndead++;
        
        for (i=0; i<lenU; i++) {
          uu[i] += Q(person,i);
          
        }
        
      }
      
      denom += risk;
      /* a contains weighted sums of x, cmat sums of squares */
      
      for (i=0; i<lenU; i++) {
        au[i] += risk*Q(person,i);
        
        for (j=0; j<=i; j++)
          cmatu(j,i) += risk*Q(person,i)*Q(person,j);
      }
      
      
      zgamma=0;
      for(i=0;i<ncolZ;i++){
        zgamma+=igammak[i]*ZB(person,i);
      }
      for(i=0;i<lenU;i++){
        zgamma+=Q(person,i)*irho[0]*iu[i];
      }
      
      ezgamma=safe_exp(zgamma);
      
      for(j=0; j<lenU;j++){
        uu[j]+=D[person]*Q(person,j)*irho[0]-Delta[person]*(Q(person,j)*irho[0]*ezgamma)/(1+ezgamma);
        for(i=j;i<lenU;i++){
          imatu(j,i)+=Delta[person]*(Q(person,i)*irho[0]*Q(person,j)*irho[0]*ezgamma)/pow(1+ezgamma,2);
        }
      }
      
      person--;
    }
    if(ndead>0){
      
      for (i=0; i<lenU; i++) {
        temp2= au[i]/ denom;  /* mean */
      uu[i] -=   ndead*temp2;
      for (j=0; j<=i; j++)
        imatu(j,i) += ndead*(cmatu(j,i) - temp2*au[j])/denom;
      }
      
    }
  }
  
  
  
  
  
  for(i=0;i<(lenU);i++){
    for(j=0;j<i;j++){
      imatu(i,j)=imatu(j,i);
    }
  }
  
  uu-=iu/itheta;
  imatu.diag(0)+=ones(lenU)/itheta;
  arma::vec u_step = safe_inv(imatu, "Information matrix for frailty update") * uu;
  u_step = damp_step(u_step, 0.5);
  u_new=iu+u_step*lr;
  
  
  List L=List::create(Named("u_new")=u_new,Named("uu")=uu,Named("imatu")=imatu);
  
  return L;
}
List tw_verticalite2(arma::vec Delta,
                     arma::vec D,
                     arma::vec ZBg,
                     arma::vec qu,
                     arma::vec irho,
                     double lr=1,
                     int mite=20,
                     double eps=1e-3,
                     double rho_bound=1.0){
  
  int nuser= Delta.size();
  int person;
  arma::mat imat(1,1);
  arma::mat imat_ivs(1,1);
  arma::vec u(1);
  double ezgamma;
  double zgamma;
  int ite=1;
  double cgs=100;
  arma::vec erho(1);
  
  while(ite<mite&&cgs>eps){
    u=zeros(1);
    imat=zeros(1,1);
    for (person=nuser-1; person>=0; person--) {
      
      zgamma=ZBg(person)+qu(person)*irho(0);
      
      ezgamma=safe_exp(zgamma);
      
      u[0]+=D(person)*qu(person)-Delta[person]*ezgamma*qu(person)/(1+ezgamma);
      
      imat(0,0)+=Delta(person)*ezgamma*qu(person)*qu(person)/pow(1+ezgamma,2);
      
    }
    
    imat_ivs=safe_inv(imat, "Information matrix for rho update");
    erho=irho+damp_step(imat_ivs*u*lr, 0.2);
    erho=erho.clamp(-rho_bound,rho_bound);
    cgs=max(abs(erho-irho));
    ite++;
    irho=erho;
  }
  List L=List::create(Named("est")=erho,Named("imat")=imat,Named("u")=u);
  return L;
}

List tw_verticalite3(arma::vec xtime,
                     arma::vec Delta,
                     arma::vec D,
                     arma::mat X,
                     arma::mat ZB,
                     arma::vec qu,
                     arma::mat Q,
                     arma::vec ibeta,
                     arma::vec iU,
                     arma::vec igammak,
                     arma::vec phih,
                     arma::vec phih2,
                     arma::vec phil,
                     arma::vec phil2,
                     arma::mat ZBz,
                     arma::vec irho,
                     double itheta){
  int nuser= xtime.size();
  int person,i,j;
  double loglik=0;
  int ncolX=X.n_cols;
  int ncolZ=ZB.n_cols;
  int lenU=Q.n_cols;
  int zeta_idx=ncolX+1+ncolZ+lenU;
  arma::mat cmat(ncolX+lenU+1,ncolX+lenU+1);
  arma::mat imat(ncolX+lenU+ncolZ+3,ncolX+lenU+ncolZ+3);
  
  arma::vec a(lenU+ncolX+1);
  arma::mat Xq=join_rows(Q,X);
  arma::mat q=Q*irho[0];
  double denom=0;
  double risk=0;
  double ezgamma=0;
  double p_rel=0;
  double resid=0;
  double dtime;
  
  int nrisk=0;
  double xbetau;
  double zgamma;
  
  double temp2;
  double a2=0;
  double dead_phih2=0;
  int ndead;
  for(i=0; i<=ncolX+lenU+ncolZ+2;i++){
    for(j=0;j<=ncolX+lenU+ncolZ+2;j++)
      imat(i,j)=0;
  }
  
  
  for(i=0; i<=ncolX+lenU;i++){
    a[i]=0;
    for(j=0;j<=ncolX+lenU;j++){
      cmat(i,j)=0;
    }
  }
  
  
  for (person=nuser-1; person>=0; ) {
    dtime = xtime[person];
    ndead =0; /*number of deaths at this time point */
    dead_phih2=0;
      
      
      
      while(person >=0 && xtime[person]==dtime) {
        /* walk through the this set of tied times */
        nrisk++;
        xbetau = 0;    /* form the term beta*z (vector mult) */
        
        for (i=0; i<ncolX; i++)
          xbetau += ibeta[i]*X(person,i);
        for (i=0; i<lenU; i++)
          xbetau += iU[i]*Q(person,i);
        
        risk = safe_exp(xbetau);
        if (Delta[person] ==1) {
          ndead++;
          loglik += xbetau;
          dead_phih2+=phih2[person];
          
          
        }
        
        denom += risk;
        /* a contains weighted sums of x, cmat sums of squares */
        
        for (i=0; i<ncolX+lenU; i++) {
          a[i] += risk*Xq(person,i);
          
          for (j=0; j<=i; j++)
            cmat(i,j) += risk*Xq(person,i)*Xq(person,j);
        }
        a[ncolX+lenU]+=risk*phih[person];
        a2+=risk*phih2[person];
        for (j=0; j<ncolX+lenU; j++)
          cmat(ncolX+lenU,j) += risk*(phih[person]*Xq(person,j));
        
        cmat(ncolX+lenU,ncolX+lenU)+=risk*(phih[person]*phih[person]);
        zgamma=0;
        for(i=0;i<ncolZ;i++){
          zgamma+=igammak[i]*ZB(person,i);
          
        }
        for(i=0;i<lenU;i++){
          zgamma+=iU[i]*Q(person,i)*irho[0];
        }
        ezgamma=safe_exp(zgamma);
        p_rel=ezgamma/(1+ezgamma);
        resid=D[person]-Delta[person]*p_rel;
     
        loglik+=D[person]*zgamma-Delta[person]*log(1+ezgamma);
        for(i=0; i<(ncolZ);i++){
          for(j=0; j<=i;j++){
            imat(ncolX+lenU+1+j,ncolX+lenU+1+i)+=Delta[person]*ezgamma*ZB(person,i)*ZB(person,j)/pow(1+ezgamma,2);
           
          }
        }
        for(j=0; j<lenU;j++){
          
          for(i=j;i<lenU;i++){
            imat(j,i)+=Delta[person]*(q(person,j)*q(person,i)*ezgamma)/pow(1+ezgamma,2);
          }
          for(i=0;i<ncolZ;i++){
            imat(j,ncolX+lenU+1+i)+=Delta[person]*(q(person,j)*ZB(person,i)*ezgamma)/pow(1+ezgamma,2);
            
          }
          imat(j,ncolX+lenU+1+ncolZ)+=Delta[person]*(q(person,j)*phil[person]*ezgamma)/pow(1+ezgamma,2);
        }
        
        
        for(i=0; i<(ncolZ);i++){
          
          imat(ncolX+lenU+1+i,zeta_idx)+=
            Delta[person]*(ezgamma*phil[person]*ZB(person,i))/pow(1+ezgamma,2);
          
        }
        
        imat(zeta_idx,zeta_idx)+=
          Delta[person]*((ezgamma*phil[person]*phil[person])/pow(1+ezgamma,2));
        
        for(i=0;i<lenU;i++){
          imat(i,ncolX+lenU+2+ncolZ)+=Delta[person]*(Q(person,i)*ezgamma*(1+ezgamma)+ezgamma*qu[person]*q(person,i))/pow(1+ezgamma,2);
          
        }
        
        
        for(i=0;i<ncolZ;i++){
          imat(ncolX+lenU+1+i,ncolX+lenU+2+ncolZ)+=Delta[person]*(ezgamma*qu[person]*ZB(person,i)/pow(1+ezgamma,2));
         
        }
        
        imat(ncolX+2+ncolZ+lenU,ncolX+2+ncolZ+lenU)+=Delta[person]*ezgamma*qu(person)*qu(person)/pow(1+ezgamma,2);
       
        
        
        person--;
      }
      if(ndead>0){
        loglik -=  ndead*log(denom);
        for (i=0; i<ncolX+lenU; i++) {
          temp2= a[i]/ denom;  /* mean */
        
        for (j=0; j<=i; j++)
          imat(j,i) += ndead*(cmat(i,j) - temp2*a[j])/denom;
        }
        
        temp2= a[ncolX+lenU]/ denom;
        for(j=0;j<ncolX+lenU;j++)
          imat(j,ncolX+lenU)+=ndead*(cmat(ncolX+lenU,j) - temp2*a[j])/denom;
        
        
        imat(ncolX+lenU,ncolX+lenU)+= ndead*(cmat(ncolX+lenU,ncolX+lenU) - temp2*a[ncolX+lenU])/denom;
      }
  }
  
  for(i=1;i<(ncolX+ncolZ+2+lenU);i++){
    for(j=0;j<i;j++){
      imat(i,j)=imat(j,i);
    }
  }
  for(i=0;i<lenU;i++){
    imat(i,i)+=1/itheta;}
  List L=List::create(Named("imatri")=safe_inv(imat, "Information matrix in tw_verticalite3"),Named("imatr")=imat,Named("loglik")=loglik);
  // List L=List::create(Named("imatr")=imat,Named("loglik")=loglik);
  return L;
}
/*double obf(arma::vec xtime,
 arma::vec Delta,
 arma::vec D,
 arma::mat Xq,
 arma::vec W,
 arma::mat ZBq,
 arma::vec ibetau,
 arma::vec igammaku,double itheta,int lenU){
 int nuser= xtime.size();
 int person,i;
 double loglik=0;
 int ncolX=Xq.n_cols;
 int ncolZ=ZBq.n_cols;
 double denom;
 double risk=0;
 double ezgamma=0;
 double dtime;
 
 int nrisk=0;
 double xbetau;
 double zgamma;
 
 
 int ndead;
 
 
 for (person=nuser-1; person>=0; ) {
 dtime = xtime[person];
 ndead =0; //number of deaths at this time point 
 
 
 
 while(person >=0 && xtime[person]==dtime) {
 // walk through the this set of tied times 
 nrisk++;
 xbetau = 0;    // form the term beta*z (vector mult) 
 
 for (i=0; i<ncolX; i++)
 xbetau += ibetau[i]*Xq(person,i);
 
 risk = exp(xbetau);
 if (Delta[person] ==1) {
 ndead++;
 loglik += xbetau;
 }
 
 denom += risk;
 // a contains weighted sums of x, cmat sums of squares 
 zgamma=0;
 for(i=0;i<ncolZ;i++){
 zgamma+=igammaku[i]*ZBq(person,i);
 }
 ezgamma=exp(zgamma);
 loglik+=D[person]*zgamma-Delta[person]*log(1+ezgamma);
 
 person--;
 }
 if(ndead>0){
 loglik -=  ndead*log(denom);}
 }
 for(i=0;i<lenU;i++){
 loglik-=ibetau(i)*ibetau(i)/itheta/2;}
 
 
 return loglik;
 }
 arma::vec len_serch(arma::vec est_ini,arma::vec d, double f0,arma::vec df,arma::vec xtime,
 arma::vec Delta,
 arma::vec D,
 arma::mat X,
 arma::vec W,
 arma::mat Z,
 arma::mat Bt,
 arma::mat Q,
 double itheta,
 arma::vec irho,
 int lenU){
 double alpha=1;
 //arma::vec alpha=ones(10);
 double gamma=0.4;
 double sigma=0.5;
 
 int nsize=est_ini.size();
 int nuser=xtime.size();
 int ncolX=X.size();
 int ncolZ=Z.size();
 int ncolBt=Bt.size();
 double h=log(nuser)/sqrt(nuser)*stddev(W);
 arma::vec est_new(nsize);
 arma::vec Phih(nuser);
 arma::vec Phil(nuser);
 arma::vec ibetau(lenU+2*ncolX+1);
 arma::vec igammaku(lenU+2*ncolZ+1+ncolBt);
 int i,j;
 double newlik,dif_lik;
 arma::vec newlik1=ones(10);
 arma::vec iU(lenU);
 arma::vec ibeta(2*ncolX+1);
 arma::vec igamma(2*ncolZ+1);
 arma::vec ikappa(ncolBt);
 double izetah,izetal;
 arma::mat Xq(nuser,2*ncolX+1+lenU);
 arma::mat ZBq(nuser,2*ncolZ+1+lenU+ncolBt);
 double thre=0;
 for(i=0;i<nsize;i++){
 thre-=df(i)*d(i)*gamma;
 }
 j=0;
 while(j>=0 && j<10){
 est_new=alpha*d+est_ini;
 iU=est_new.subvec(0,lenU-1);
 ibeta=est_new.subvec(lenU,lenU+2*ncolX);
 izetah=est_new(lenU+2*ncolX+1);
 igamma=est_new.subvec(lenU+2*ncolX+2,lenU+2*ncolX+1+2*ncolZ+1);
 ikappa=est_new.subvec(lenU+2*ncolX+2+2*ncolZ+1,lenU+2*ncolX+1+2*ncolZ+1+ncolBt);
 izetal=est_new(lenU+2*ncolX+1+2*ncolZ+1+ncolBt+1);
 Phih=normcdf((W-ones(nuser)*izetah)/h);
 Phil=normcdf((W-ones(nuser)*izetal)/h);
 Xq=join_rows(Q,X,ones(nuser),X);
 for(i =0;i<=ncolX;i++ ){
 Xq.col(lenU+i+ncolX)=Xq.col(lenU+i+ncolX)%Phih;}
 ZBq=join_rows(Q*irho(0),join_rows(Z,ones(nuser),Z),Bt);
 for(i =0;i<=ncolZ;i++ )
 ZBq.col(lenU+ncolZ+i)=ZBq.col(lenU+ncolZ+i)%Phil;
 ibetau=join_cols(iU,ibeta);
 igammaku=join_cols(iU,igamma,ikappa);
 newlik=obf( xtime,
 Delta,
 D,
 Xq,
 W,
 ZBq,
 ibetau,
 igammaku,itheta,lenU);
 dif_lik=-newlik+f0;
 newlik1(j)=newlik;
 
 thre*=alpha;
 if(dif_lik<=thre){
 j=20;}
 else{
 alpha=sigma*alpha;
 j++;}
 
 }
 //return alpha;
 return newlik1;
 }
 */

