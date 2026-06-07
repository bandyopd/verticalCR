// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//[[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


List th_updateb(arma::vec xtime,
                arma::vec Delta,
                arma::mat X,
                arma::mat Q,
                arma::vec ibeta,
                arma::vec iu,
                arma::vec izetah,
                arma::vec phih,
                double lr);


List th_updateg(arma::vec Delta,
                arma::vec D1,
                arma::vec D2,
                arma::mat Q,
                arma::mat ZB,
                arma::vec iu,
                arma::vec igammak1,
                arma::vec igammak2,
                arma::vec phil1,
                arma::vec phil2,
                arma::vec dphil1,
                arma::vec dphil2,
                arma::vec izetal,
                arma::vec irho,
                double lr);


List th_updateu(arma::vec xtime,
                arma::vec Delta,
                arma::vec D1,
                arma::vec D2,
                arma::mat X,
                arma::mat Q,
                arma::mat ZB,
                arma::vec ibeta,
                arma::vec iu,
                arma::vec igammak1,
                arma::vec igammak2,
                arma::vec irho,
                double itheta,
                  double lr);


List th_verticalite2(arma::vec Delta,
                     arma::vec D1,
                     arma::vec D2,
                     arma::vec ZBg1,
                     arma::vec ZBg2,
                     arma::vec qu,
                     arma::vec irho,
                     double lr,
                     int mite,
                     double eps
);

List th_verticalite3(arma::vec xtime,
                     arma::vec Delta,
                     arma::vec D1,
                     arma::vec D2,
                     arma::mat X,
                     arma::mat ZB,
                     arma::vec qu,
                     arma::mat Q,
                     arma::vec ibeta,
                     arma::vec iU,
                     arma::vec igammak1,
                     arma::vec igammak2,
                     arma::vec phih,
                     arma::vec phil1,
                     arma::vec phil2,
                     arma::vec irho,
                     double itheta);

// [[Rcpp::export]]
List th_vertical(arma::vec xtime,
                 arma::vec Delta,
                 arma::vec D1,
                 arma::vec D2,
                 arma::mat X,
                 arma::mat Z,
                 arma::mat Bt,
                 arma::vec W,
                 arma::mat Q,
                 arma::vec ibeta,
                 arma::vec igamma1,
                 arma::vec igamma2,
                 arma::vec ikappa1,
                 arma::vec ikappa2,
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
                 double eps=1E-3
){
  int i;
  int nuser= xtime.size();//观测时间数目
  int ncolX=X.n_cols;//特定原因生存模型协变量维度
  int ncolZ=Z.n_cols;//相对风险模型协变量维度
  int ncolBt=Bt.n_cols;//相对风险模型非参数部分B样条维度
  int ncolQ=Q.n_cols;//随机效应维度
  double cgs=100,cgsr=100;//迭代差初始值
  int ite=1,iter=1;//迭代初始值
  arma::vec iUc(ncolQ+1);//设置随机效应变量
  arma::vec se_rho(1);//设置1维随机效应系数变量
  arma::vec ZBg1(nuser);//设置汇总相对风险中协变量和非参部分的变量
  arma::vec ZBg2(nuser);//设置汇总相对风险中协变量和非参部分的变量
  arma::vec qu(nuser);//设置随机效应值的变量
  arma::vec phih(nuser);//特定原因生存模型S曲线-正态密度的变量
  arma::vec phil1(nuser);//相对风险模型S曲线-正态密度的变量
  arma::vec phil2(nuser);//相对风险模型S曲线-正态密度的变量
  arma::vec dphil1(nuser);//相对风险模型S曲线-正态密度的变量
  arma::vec dphil2(nuser);//相对风险模型S曲线-正态密度的变量
  
 

  
  arma::vec Phih(nuser);//特定原因生存模型S曲线-正态分布的变量
  arma::vec Phil(nuser);//相对风险模型S曲线-正态分布的变量
  arma::vec igammak1(2*ncolZ+1+ncolBt);//相对风险模型case1随机效应、回归系数和B样条系数变量
  arma::vec igammak2(2*ncolZ+1+ncolBt);//相对风险模型case2随机效应、回归系数和B样条系数变量

  double etheta;//随机效应方差
  arma::mat Xt(nuser,2*ncolX+1);//特定原因生存模型协变量和分组变量的整合
  arma::mat ZB(nuser,2*ncolZ+1+ncolBt);//相对风险模型协变量、分组变量和B样条的整合
  arma::mat imatr(ncolQ+2*ncolX+1+4*ncolZ+2+2+ncolBt,ncolQ+2*ncolX+1+2+4*ncolZ+2+ncolBt);//目标函数2阶导
  arma::mat imatr_ivs(ncolQ+2*ncolX+1+4*ncolZ+2+2+ncolBt,ncolQ+2*ncolX+1+2+4*ncolZ+2+ncolBt);//目标函数2阶导的逆
  arma::mat T33(ncolQ,ncolQ);//随机效应协方差
  
  
  arma::vec ur(ncolQ+2*ncolX+1+2*ncolZ+1+2+ncolBt);//目标函数的得分函数
  arma::vec sxtime=sort(xtime);//生存时间排序
  Delta=Delta(sort_index(xtime));//将删失指标根据生存时间排序
  D1=D1(sort_index(xtime));//将病因指标根据生存时间排序
  D2=D2(sort_index(xtime));//将病因指标根据生存时间排序
  W=W(sort_index(xtime));//将阈值变量根据生存时间排序
  X=X.rows(sort_index(xtime));//将特定病因生存模型协变量根据生存时间排序
  Z=Z.rows(sort_index(xtime));//将相对风险模型协变量根据生存时间排序
  Q=Q.rows(sort_index(xtime));//将分组变量根据生存时间排序
  Bt=Bt.rows(sort_index(xtime));//将B样条根据生存时间排序
   iter=1;
    // Xt=join_rows(X,ones(nuser),X);
    // ZB=join_rows(join_rows(Z,ones(nuser),Z),Bt);
   igammak1=join_cols(igamma1,ikappa1);
   igammak2=join_cols(igamma2,ikappa2);
 while(iter<miter && cgsr>1E-3){
    ite=1;
    cgs=100;
 while(ite<mite && cgs>eps){
      //计算特定病因生存模型关于阈值的一阶导
      Xt=join_rows(X,ones(nuser),X);
      ZB=join_rows(join_rows(Z,ones(nuser),Z),Bt);
      phih=normpdf((W-ones(nuser)*izetah)/h1)*(-1/h1)%(ones(nuser)*ibeta(ncolX)+X*ibeta.subvec(ncolX+1,ncolX*2));
      //计算相对风险模型关于阈值的一阶导
      phil1=normpdf((W-ones(nuser)*izetal)/h2)*(-1/h2)%(ones(nuser)*igammak1[ncolZ]+Z*igammak1.subvec(ncolZ+1,ncolZ*2));
      phil2=normpdf((W-ones(nuser)*izetal)/h2)*(-1/h2)%(ones(nuser)*igammak2[ncolZ]+Z*igammak2.subvec(ncolZ+1,ncolZ*2));
    
      //计算特定病因生存模型阈值部分S曲线-正态分布取值
      Phih=normcdf((W-ones(nuser)*izetah)/h1);
      //计算相对风险模型阈值部分S曲线-正态分布取值
      Phil=normcdf((W-ones(nuser)*izetal)/h2);
//       //合并特定病因生存模型中的分组变量和协变量
     
     for(i=0;i<nuser;i++){
       dphil1[i]=(W[i]-izetal[0])/pow(h2,2)*phil1[i];
       dphil2[i]=(W[i]-izetal[0])/pow(h2,2)*phil2[i];
     }
      for(i =0;i<=ncolX;i++ ){
        
        Xt.col(i+ncolX)=Xt.col(i+ncolX)%Phih;}

      //合并相对风险模型中的分组变量、协变量和B样条

      for(i =0;i<=ncolZ;i++){
        
        ZB.col(ncolZ+i)=ZB.col(ncolZ+i)%Phil;

        }

//合并相对风险模型中的随机效应、回归系数和B样条系数


List Lu=th_updateu(sxtime,
                  Delta,
                  D1,
                  D2,
                  Xt,
                  Q,
                  ZB,
                  ibeta,
                  iU,
                  igammak1,
                  igammak2,
                  irho,
                  itheta,
                  lr);
//计算给定参数下目标函数的一阶向量和二阶矩阵


arma::vec u_temp = Lu["u_new"];

      List Lb= th_updateb( sxtime,
                               Delta,
                               Xt,
                               Q,
                               ibeta,
                               u_temp,
                               izetah,
                               phih,
                               lr
      );
      arma::vec beta_temp = Lb["beta_new"];
      arma::vec zetah_temp = Lb["zetah_new"];

    List Lg= th_updateg( Delta,
                     D1,
                     D2,
                    Q,
                    ZB,
                     u_temp,
                     igammak1,
                     igammak2,
                     phil1,
                     phil2,
                     dphil1,
                     dphil2,
                     izetal,
                     irho,
                     lr);
   arma::vec gammak1_temp = Lg["gammak1_new"];
   arma::vec gammak2_temp = Lg["gammak2_new"];
   arma::vec zetal_temp = Lg["zetal_new"];
   arma::mat imatg=Lg["imatg"];
   cgs=max(abs(ibeta-beta_temp));

   iU=u_temp.clamp(-5,5);
   ibeta=beta_temp.clamp(-5,5);
   izetah=zetah_temp.clamp(0.1,0.9);
   igammak1=gammak1_temp.clamp(-5,5);
   igammak2=gammak2_temp.clamp(-5,5);
   izetal=zetal_temp.clamp(0.1,0.9);

   ite++;
 }

 ZBg1= ZB*igammak1;
 ZBg2= ZB*igammak2;
 qu= Q*iU;
 List R2=th_verticalite2(Delta,D1,D2,ZBg1,ZBg2,qu,irho,lr,mite=20,eps);
 arma::vec rho_temp=R2["est"];
  List L4= th_verticalite3(sxtime, Delta, D1,D2, Xt,ZB,qu, Q,ibeta,iU, igammak1,igammak2,
                           phih,
                           phil1,
                           phil2,
                           rho_temp,
                           itheta);
  arma::mat imat_temp=L4["imatri"];

  T33=imat_temp.submat(0,0,ncolQ-1,ncolQ-1);
  //etheta=(sum(T33.diag(0))+sum((iUc-ones(ncolQ+1)*mean(iUc))%(iUc-ones(ncolQ+1)*mean(iUc))))/(ncolQ+1);
  etheta=sum(T33.diag(0)+iU%iU)/(ncolQ);
  cgsr=max(abs(rho_temp-irho));

  irho=rho_temp;
  itheta=etheta;
  iter++;
}

arma::mat T33c=T33*T33;
List R4= th_verticalite3(sxtime, Delta, D1,D2, Xt,ZB,qu, Q,ibeta,iU, igammak1,igammak2,
                         phih,
                         phil1,
                         phil2,
                         irho,
                         itheta);
arma::mat temp4=R4["imatri"];
arma::vec se_para=sqrt(abs(temp4.diag(0)));
double se_theta=sqrt(2*itheta*itheta/(ncolQ-2*sum(T33.diag(0))/itheta+sum(T33c.diag(0))/(itheta*itheta)));

iUc=join_cols(zeros(1),iU);
iUc-=mean(iUc);



  List R=List::create(Named("beta")=ibeta,
                      Named("se_beta")=se_para.subvec(ncolQ,ncolQ+2*ncolX),
                      Named("zetah")=izetah,
                      Named("se_zetah")=se_para(ncolQ+2*ncolX+1),
                      Named("gamma1")=igammak1.subvec(0,(2*ncolZ)),
                      Named("se_gamma1")=se_para.subvec(ncolQ+2*ncolX+2,ncolQ+2*ncolX+1+2*ncolZ+1),
                      Named("kappa1")=igammak1.subvec((2*ncolZ+1),(2*ncolZ+ncolBt)),
                      Named("se_kappa1")=se_para.subvec(ncolQ+2*ncolX+2+2*ncolZ+1,ncolQ+2*ncolX+1+2*ncolZ+1+ncolBt),
                      Named("gamma2")=igammak2.subvec(0,(2*ncolZ)),
                      Named("se_gamma2")=se_para.subvec(ncolQ+2*ncolX+1+2*ncolZ+1+ncolBt+1,ncolQ+2*ncolX+1+4*ncolZ+2+ncolBt),
                      Named("kappa2")=igammak2.subvec((2*ncolZ+1),(2*ncolZ+ncolBt)),
                      Named("se_kappa2")=se_para.subvec(ncolQ+2*ncolX+1+4*ncolZ+2+ncolBt+1,ncolQ+2*ncolX+1+4*ncolZ+2+2*ncolBt),
                      Named("zetal")=izetal,
                      Named("se_zetal")=se_para(ncolQ+2*ncolX+1+4*ncolZ+2+2*ncolBt+1),
                      Named("rho1")=irho(0),
                      Named("se_rho1")=se_para(ncolQ+2*ncolX+1+4*ncolZ+2+2*ncolBt+1+1),
                      Named("rho2")=irho(1),
                      Named("se_rho2")=se_para(ncolQ+2*ncolX+1+4*ncolZ+2+2*ncolBt+1+2),
                      Named("theta")=itheta,
                      Named("se_theta")=se_theta,
                      Named("U")=iUc,
                      Named("U_se")=join_cols(zeros(1),se_para.subvec(0,ncolQ-1)),
                      Named("iter")=iter
  );

  
  // List R=List::create(Named("beta")=1);
  return R;
}
//计算给定参数下目标函数的一阶得分向量和二阶信息矩阵

List th_updateb(arma::vec xtime,
                     arma::vec Delta,
                     arma::mat X,
                     arma::mat Q,
                     arma::vec ibeta,
                     arma::vec iu,
                     arma::vec izetah,
                     arma::vec phih,
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
      risk = exp(xbetau);
      if (Delta[person] ==1) {
        ndead++;
        
        
        for (i=0; i<ncolX; i++) {
          ub[i] += X(person,i);
          
        }
        ub[ncolX]+=phih[person];
      }
      
      denom += risk;
      /* a contains weighted sums of x, cmat sums of squares */
      
      for (i=0; i<ncolX; i++) {
        ab[i] += risk*X(person,i);
        
        for (j=0; j<=i; j++)
          cmatb(j,i) += risk*(X(person,i)*X(person,j));
      }
      ab[ncolX]+=risk*phih[person];
      for (j=0; j<ncolX; j++){
        cmatb(j,ncolX) += risk*(phih[person]*X(person,j));
      // imatb(j,ncolX)+=-Delta(person)*bphih(person,j);
      }
      
      cmatb(ncolX,ncolX)+=risk*(phih[person]*phih[person]);
      // imatb(ncolX,ncolX)+=-Delta(person)*dphih(person);
      
     
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
  if(max(inv(imatb)*ub)>10){
    betaz_new=join_cols(ibeta,izetah)+inv(imatb)*ub*lr;
  }else{
    betaz_new=join_cols(ibeta,izetah)+inv(imatb)*ub;
  }
  
  
  List L=List::create(Named("beta_new")=betaz_new.subvec(0,(ncolX-1)),Named("zetah_new")=betaz_new[ncolX],Named("imatb")=imatb,Named("ub")=ub);
  
  return L;
}
List th_updateg(arma::vec Delta,
                arma::vec D1,
                arma::vec D2,
                arma::mat Q,
                arma::mat ZB,
                arma::vec iu,
                arma::vec igammak1,
                arma::vec igammak2,
                arma::vec phil1,
                arma::vec phil2,
                arma::vec dphil1,
                arma::vec dphil2,
                arma::vec izetal,
                arma::vec irho,
                double lr){
  int nuser= Delta.size();
  int person,i,j;
  int ncolZ=ZB.n_cols;
  int lenU=Q.n_cols;
  arma::mat imatg=zeros(2*ncolZ+1,2*ncolZ+1);
  arma::vec ug=zeros(2*ncolZ+1);
  arma::vec gamkapz_new(2*ncolZ+1);
 
  double ezgamma1=0;
  double ezgamma2=0;

  double zgamma1;
  double zgamma2;

  for (person=nuser-1; person>=0; ) {
        /* walk through the this set of tied times */
        
        zgamma1=0;
        zgamma2=0;
        for(i=0;i<ncolZ;i++){
          zgamma1+=igammak1[i]*ZB(person,i);
          zgamma2+=igammak2[i]*ZB(person,i);
        }
        
        for(i=0;i<lenU;i++){
          zgamma1+=Q(person,i)*irho[0]*iu[i];
          zgamma2+=Q(person,i)*irho[1]*iu[i];
        }
        
        ezgamma1=exp(zgamma1);
        ezgamma2=exp(zgamma2);
        for(i=0; i<(ncolZ);i++){
          ug[i]+=D1[person]*ZB(person,i)-Delta[person]*ezgamma1*ZB(person,i)/(1+ezgamma1+ezgamma2);
          ug[ncolZ+i]+=D2[person]*ZB(person,i)-Delta[person]*ezgamma2*ZB(person,i)/(1+ezgamma1+ezgamma2);
          
          for(j=0; j<=i;j++){
            imatg(j,i)+=Delta[person]*ezgamma1*(1+ezgamma2)*ZB(person,i)*ZB(person,j)/pow(1+ezgamma1+ezgamma2,2);
            imatg(ncolZ+j,ncolZ+i)+=Delta[person]*ezgamma2*(1+ezgamma1)*ZB(person,i)*ZB(person,j)/pow(1+ezgamma1+ezgamma2,2);
            
          }
          for(j=0;j<ncolZ;j++){
            imatg(j,ncolZ+i)+=-Delta[person]*ezgamma1*ezgamma2*ZB(person,i)*ZB(person,j)/pow(1+ezgamma1+ezgamma2,2);
          }
          
        }
        
        ug[2*ncolZ]+=D1[person]*phil1[person]+D2[person]*phil2[person]-Delta[person]*(ezgamma1*phil1[person]+ezgamma2*phil2[person])/(1+ezgamma1+ezgamma2);
        for(i=0; i<(ncolZ);i++){

          imatg(i,2*ncolZ)+=Delta[person]*(ezgamma1*(1+ezgamma2)*phil1[person]*ZB(person,i)-ezgamma1*ezgamma2*phil1[person]*ZB(person,i))/pow(1+ezgamma1+ezgamma2,2);
          imatg(ncolZ+i,2*ncolZ)+=Delta[person]*(ezgamma2*(1+ezgamma1)*phil2[person]*ZB(person,i)-ezgamma1*ezgamma2*phil2[person]*ZB(person,i))/pow(1+ezgamma1+ezgamma2,2);
         }
        
        //imatg(2*ncolZ,2*ncolZ)+=Delta[person]*(ezgamma1*phil1[person]*phil1[person]+ezgamma2*phil2[person]*phil2[person]+ezgamma1*ezgamma2*pow(phil1[person]-phil2[person],2))/pow(1+ezgamma1+ezgamma2,2);
        imatg(2*ncolZ,2*ncolZ)+=Delta[person]*(ezgamma1*phil1[person]*phil1[person]+ezgamma2*phil2[person]*phil2[person]+ezgamma1*ezgamma2*pow(phil1[person]-phil2[person],2))/pow(1+ezgamma1+ezgamma2,2);
        person--;
  }
  for(i=0;i<(2*ncolZ+1);i++){
    for(j=0;j<i;j++){
      imatg(i,j)=imatg(j,i);
    }
  }
  if(max(inv(imatg)*ug)>10){
    gamkapz_new=join_cols(igammak1,igammak2,izetal)+inv(imatg)*ug*lr;
  }else{
    gamkapz_new=join_cols(igammak1,igammak2,izetal)+inv(imatg)*ug;
  }
  
  
   
   List L=List::create(Named("gammak1_new")=gamkapz_new.subvec(0,(ncolZ-1)),Named("gammak2_new")=gamkapz_new.subvec(ncolZ,(2*ncolZ-1)),Named("zetal_new")=gamkapz_new[2*ncolZ],Named("imatg")=imatg);
   // List L=List::create(Named("imatg")=imatg);
  return L;
}

List th_updateu(arma::vec xtime,
                arma::vec Delta,
                arma::vec D1,
                arma::vec D2,
                arma::mat X,
                arma::mat Q,
                arma::mat ZB,
                arma::vec ibeta,
                arma::vec iu,
                arma::vec igammak1,
                arma::vec igammak2,
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
  double ezgamma1;
  double ezgamma2;
  double dtime;
  
  int nrisk=0;
  double xbetau;
  double zgamma1;
  double zgamma2;
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
          risk = exp(xbetau);
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
          
          
          zgamma1=0;
          zgamma2=0;
          for(i=0;i<ncolZ;i++){
            zgamma1+=igammak1[i]*ZB(person,i);
            zgamma2+=igammak2[i]*ZB(person,i);
          }
          for(i=0;i<lenU;i++){
            zgamma1+=Q(person,i)*irho[0]*iu[i];
            zgamma2+=Q(person,i)*irho[1]*iu[i];
          }
          
          ezgamma1=exp(zgamma1);
          ezgamma2=exp(zgamma2);
         
          for(j=0; j<lenU;j++){
            uu[j]+=D1[person]*Q(person,j)*irho[0]+D2[person]*Q(person,j)*irho[1]-Delta[person]*(Q(person,j)*irho[0]*ezgamma1+Q(person,j)*irho[1]*ezgamma2)/(1+ezgamma1+ezgamma2);
            for(i=j;i<lenU;i++){
              imatu(j,i)+=Delta[person]*(Q(person,i)*irho[0]*Q(person,j)*irho[0]*ezgamma1*(1+ezgamma2)+Q(person,i)*irho[1]*Q(person,j)*irho[1]*ezgamma2*(1+ezgamma1)-2*ezgamma1*ezgamma2*Q(person,i)*irho[0]*Q(person,j)*irho[1])/pow(1+ezgamma1+ezgamma2,2);
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
  if(max(inv(imatu)*uu)>10){
    u_new=iu+inv(imatu)*uu*lr;
  }else{
    u_new=iu+inv(imatu)*uu;
  }
  
  
  u_new-=ones(lenU)*mean(u_new);
  List L=List::create(Named("u_new")=u_new,Named("uu")=uu,Named("imatu")=imatu);
  
  return L;
}
List th_verticalite2(arma::vec Delta,
                     arma::vec D1,
                     arma::vec D2,
                     arma::vec ZBg1,
                     arma::vec ZBg2,
                     arma::vec qu,
                     arma::vec irho,
                     double lr=1,
                     int mite=20,
                     double eps=1e-3){
  
  int nuser= Delta.size();
  int person;
  arma::mat imat(2,2);
  arma::mat imat_ivs(2,2);
  arma::vec u(2);
  double ezgamma1;
  double ezgamma2;
  double zgamma1;
  double zgamma2;
  int ite=1;
  double cgs=100;
  arma::vec erho(2);
  
  while(ite<mite&&cgs>eps){
    u=zeros(2);
    imat=zeros(2,2);
    for (person=nuser-1; person>=0; person--) {
      
      zgamma1=ZBg1(person)+qu(person)*irho(0);
      zgamma2=ZBg2(person)+qu(person)*irho(1);
      ezgamma1=exp(zgamma1);
      ezgamma2=exp(zgamma2);
      
      u[0]+=D1(person)*qu(person)-Delta[person]*ezgamma1*qu(person)/(1+ezgamma1+ezgamma2);
      u[1]+=D2(person)*qu(person)-Delta[person]*ezgamma2*qu(person)/(1+ezgamma1+ezgamma2);
      imat(0,0)+=Delta(person)*ezgamma1*(1+ezgamma2)*qu(person)*qu(person)/pow(1+ezgamma1+ezgamma2,2);
      imat(1,0)+=-Delta(person)*ezgamma1*ezgamma2*qu(person)*qu(person)/pow(1+ezgamma1+ezgamma2,2);
      imat(0,1)+=-Delta(person)*ezgamma1*ezgamma2*qu(person)*qu(person)/pow(1+ezgamma1+ezgamma2,2);
      imat(1,1)+=Delta(person)*(1+ezgamma1)*ezgamma2*qu(person)*qu(person)/pow(1+ezgamma1+ezgamma2,2);
    }
    
    imat_ivs=inv(imat);
    if(max(imat_ivs*u)>10){
      erho=irho+imat_ivs*u*lr;
    }else{
      erho=irho+imat_ivs*u;
    }
    cgs=max(abs(erho-irho));
    ite++;
    irho=erho;
  }
  List L=List::create(Named("est")=erho,Named("imat")=imat,Named("u")=u);
  return L;
}

List th_verticalite3(arma::vec xtime,
                     arma::vec Delta,
                     arma::vec D1,
                     arma::vec D2,
                     arma::mat X,
                     arma::mat ZB,
                     arma::vec qu,
                     arma::mat Q,
                     arma::vec ibeta,
                     arma::vec iU,
                     arma::vec igammak1,
                     arma::vec igammak2,
                     arma::vec phih,
                     arma::vec phil1,
                     arma::vec phil2,
                     arma::vec irho,
                     double itheta){
  int nuser= xtime.size();
  int person,i,j;
  double loglik=0;
  int ncolX=X.n_cols;
  int ncolZ=ZB.n_cols;
  int lenU=Q.n_cols;
  arma::mat cmat(ncolX+lenU+1,ncolX+lenU+1);
  arma::mat imat(ncolX+lenU+2*ncolZ+4,ncolX+lenU+2*ncolZ+4);
  
  arma::vec a(lenU+ncolX+1);
  arma::mat Xq=join_rows(Q,X);
  arma::mat q1=Q*irho[0];
  arma::mat q2=Q*irho[1];
  double denom=0;
  double risk=0;
  double ezgamma1=0;
  double ezgamma2=0;
  double dtime;
  
  int nrisk=0;
  double xbetau;
  double zgamma1;
  double zgamma2;
  
  double temp2;
  int ndead;
  for(i=0; i<=ncolX+lenU+2*ncolZ+3;i++){
    for(j=0;j<=ncolX+lenU+2*ncolZ+3;j++)
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
      
      
      
      while(person >=0 && xtime[person]==dtime) {
        /* walk through the this set of tied times */
        nrisk++;
        xbetau = 0;    /* form the term beta*z (vector mult) */
        
        for (i=0; i<ncolX; i++)
          xbetau += ibeta[i]*X(person,i);
        for (i=0; i<lenU; i++)
          xbetau += iU[i]*Q(person,i);
        
        risk = exp(xbetau);
        if (Delta[person] ==1) {
          ndead++;
          loglik += xbetau;
          
          
        }
        
        denom += risk;
        /* a contains weighted sums of x, cmat sums of squares */
        
        for (i=0; i<ncolX+lenU; i++) {
          a[i] += risk*Xq(person,i);
          
          for (j=0; j<=i; j++)
            cmat(i,j) += risk*Xq(person,i)*Xq(person,j);
        }
        a[ncolX+lenU]+=risk*phih[person];
        for (j=0; j<ncolX+lenU; j++)
          cmat(ncolX+lenU,j) +=  risk*phih[person]*Xq(person,j);
       // for(j=0;j<ncolX;j++){
       //   cmat(ncolX+lenU,j+lenU) +=  risk*bphih(person,j);
       //   imat(lenU+j,ncolX+lenU)+=-Delta(person)*bphih(person,j);}
       
        cmat(ncolX+lenU,ncolX+lenU)+=risk*phih[person]*phih[person];
      // imat(ncolX+lenU,ncolX+lenU)+=-Delta(person)*dphih(person);
        zgamma1=zgamma2=0;
        for(i=0;i<ncolZ;i++){
          zgamma1+=igammak1[i]*ZB(person,i);
          zgamma2+=igammak2[i]*ZB(person,i);
        }
        for(i=0;i<lenU;i++){
          zgamma1+=iU[i]*Q(person,i)*irho[0];
          zgamma2+=iU[i]*Q(person,i)*irho[1];;
        }
        ezgamma1=exp(zgamma1);
        ezgamma2=exp(zgamma2);
        loglik+=D1[person]*zgamma1+D2[person]*zgamma2-Delta[person]*log(1+ezgamma1+ezgamma2);
        for(i=0; i<(ncolZ);i++){
          for(j=0; j<=i;j++){
            imat(ncolX+lenU+1+j,ncolX+lenU+1+i)+=Delta[person]*ezgamma1*(1+ezgamma2)*ZB(person,i)*ZB(person,j)/pow(1+ezgamma1+ezgamma2,2);
            imat(ncolX+lenU+1+ncolZ+j,ncolX+lenU+ncolZ+1+i)+=Delta[person]*ezgamma2*(1+ezgamma1)*ZB(person,i)*ZB(person,j)/pow(1+ezgamma1+ezgamma2,2);
            
          }
          for(j=0;j<ncolZ;j++){
            imat(ncolX+lenU+1+j,ncolX+1+lenU+ncolZ+i)+=-Delta[person]*ezgamma1*ezgamma2*ZB(person,i)*ZB(person,j)/pow(1+ezgamma1+ezgamma2,2);
          }
          
        }
        for(j=0; j<lenU;j++){
          
          for(i=j;i<lenU;i++){
            imat(j,i)+=Delta[person]*(q1(person,j)*q1(person,i)*ezgamma1*(1+ezgamma2)+q2(person,j)*q2(person,i)*ezgamma2*(1+ezgamma1)-2*ezgamma1*ezgamma2*q1(person,j)*q2(person,i))/pow(1+ezgamma1+ezgamma2,2);
          }
          for(i=0;i<ncolZ;i++){
            imat(j,ncolX+lenU+1+i)+=Delta[person]*(q1(person,j)*ZB(person,i)*ezgamma1*(1+ezgamma2)-q2(person,j)*ZB(person,i)*ezgamma1*ezgamma2)/pow(1+ezgamma1+ezgamma2,2);
            imat(j,ncolX+lenU+1+ncolZ+i)+=Delta[person]*(q2(person,j)*ZB(person,i)*ezgamma2*(1+ezgamma1)-q1(person,j)*ZB(person,i)*ezgamma1*ezgamma2)/pow(1+ezgamma1+ezgamma2,2);
          }
          imat(j,ncolX+lenU+1+2*ncolZ)+=Delta[person]*(q1(person,j)*phil1[person]*ezgamma1*(1+ezgamma2)+q2(person,j)*phil2[person]*ezgamma2*(1+ezgamma1)-(q1(person,j)*phil2[person]+q2(person,j)*phil1[person])*ezgamma1*ezgamma2)/pow(1+ezgamma1+ezgamma2,2);
        }
        
        
        for(i=0; i<(ncolZ);i++){
         
          imat(ncolX+lenU+1+i,ncolX+1+2*ncolZ+lenU)+=Delta[person]*(ezgamma1*(1+ezgamma2)*phil1[person]*ZB(person,i)-ezgamma1*ezgamma2*phil1[person]*ZB(person,i))/pow(1+ezgamma1+ezgamma2,2);
          imat(ncolX+lenU+ncolZ+1+i,ncolX+1+2*ncolZ+lenU)+=Delta[person]*(ezgamma2*(1+ezgamma1)*phil2[person]*ZB(person,i)-ezgamma1*ezgamma2*phil2[person]*ZB(person,i))/pow(1+ezgamma1+ezgamma2,2);
        }
        
        imat(ncolX+1+2*ncolZ+lenU,ncolX+1+2*ncolZ+lenU)+=Delta[person]*((ezgamma1*phil1[person]*phil1[person]+ezgamma2*phil2[person]*phil2[person]+ezgamma1*ezgamma2*pow(phil1[person]-phil2[person],2))/pow(1+ezgamma1+ezgamma2,2));
        
        for(i=0;i<lenU;i++){
          imat(i,ncolX+lenU+2+2*ncolZ)+=-D1[person]*Q(person,i)+Delta[person]*(Q(person,i)*ezgamma1*(1+ezgamma1+ezgamma2)+ezgamma1*(1+ezgamma2)*qu[person]*q1(person,i)-qu[person]*q2(person,i)*ezgamma1*ezgamma2)/pow(1+ezgamma1+ezgamma2,2);
          imat(i,ncolX+lenU+3+2*ncolZ)+=-D2[person]*Q(person,i)+Delta[person]*(Q(person,i)*ezgamma2*(1+ezgamma1+ezgamma2)+ezgamma2*(1+ezgamma1)*qu[person]*q2(person,i)-qu[person]*q1(person,i)*ezgamma1*ezgamma2)/pow(1+ezgamma1+ezgamma2,2);
        }
        
        
        for(i=0;i<ncolZ;i++){
          imat(ncolX+lenU+1+i,ncolX+lenU+2+2*ncolZ)+=Delta[person]*(ezgamma1*(1+ezgamma2)*qu[person]*ZB(person,i)/pow(1+ezgamma1+ezgamma2,2));
          imat(ncolX+lenU+1+ncolZ+i,ncolX+lenU+2+2*ncolZ)+=-Delta[person]*(ezgamma2*ezgamma1)*qu[person]*ZB(person,i)/pow(1+ezgamma1+ezgamma2,2);
          imat(ncolX+lenU+1+i,ncolX+lenU+3+2*ncolZ)+=-Delta[person]*(ezgamma1*ezgamma2*qu[person]*ZB(person,i)/pow(1+ezgamma1+ezgamma2,2));
          imat(ncolX+lenU+1+ncolZ+i,ncolX+lenU+3+2*ncolZ)+=Delta[person]*(ezgamma2*(1+ezgamma1)*qu[person]*ZB(person,i)/pow(1+ezgamma1+ezgamma2,2));
        }
        
        imat(ncolX+2+2*ncolZ+lenU,ncolX+2+2*ncolZ+lenU)+=Delta[person]*ezgamma1*(1+ezgamma2)*qu(person)*qu(person)/pow(1+ezgamma1+ezgamma2,2);
        imat(ncolX+2+2*ncolZ+lenU,ncolX+3+2*ncolZ+lenU)+=-Delta[person]*ezgamma1*ezgamma2*qu(person)*qu(person)/pow(1+ezgamma1+ezgamma2,2);
        imat(ncolX+3+2*ncolZ+lenU,ncolX+3+2*ncolZ+lenU)+=Delta[person]*ezgamma2*(1+ezgamma1)*qu(person)*qu(person)/pow(1+ezgamma1+ezgamma2,2);
        
        
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
  
  for(i=1;i<(ncolX+2*ncolZ+4+lenU);i++){
    for(j=0;j<i;j++){
      imat(i,j)=imat(j,i);
    }
  }
  for(i=0;i<lenU;i++){
    imat(i,i)+=1/itheta;}
  List L=List::create(Named("imatri")=inv(imat),Named("imatr")=imat,Named("loglik")=loglik);
  // List L=List::create(Named("imatr")=imat,Named("loglik")=loglik);
  return L;
}


