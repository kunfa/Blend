#include<RcppArmadillo.h>
#include<stdio.h>
#include<vector>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;
//using namespace R;

// [[Rcpp::export()]]

Rcpp::List RBV_SI(arma::vec y, arma::mat e, arma::mat g, arma:: mat w, arma:: mat C,int maxSteps, int q,arma::vec hatBeta, arma:: mat hatEta, arma::vec hatAlpha, arma::vec hatata, double hatTau, arma::vec hatV, arma::vec hatSg1,arma::vec hatSg2, double hatSg3, double hatSg4,arma::mat invSigAlpha0, double hatPi3, double hatPi4, double hatEtaSq1, double hatEtaSq2,double hatEtaSq3, double hatEtaSq4, double xi1, double xi2, double r1,double r2,double r3, double r4,double a, double b, double sh1, double sh0,int progress)
{
  unsigned int n = g.n_rows,L = q-1, m = g.n_cols,p = w.n_cols,q1 = e.n_cols;
  arma::mat gsAlpha(maxSteps, q1),
  gsBeta(maxSteps,m),
  gsAta(maxSteps,3),
  gseta(maxSteps,p),
  gsV(maxSteps, n),
  gsSg1(maxSteps, m),
  gsSg2(maxSteps, m)
    ;
  
  arma::vec gsEtaSq1(maxSteps),
  gsEtaSq2(maxSteps),
  gsTau(maxSteps),
  gsEtaSq3(maxSteps),
  gsEtaSq4(maxSteps),
  gsPi3(maxSteps),
  gsPi4(maxSteps),
  gsSg3(maxSteps),
  gsSg4(maxSteps)
  ;
  arma::mat temp1;
  
  double meanb;
  double varb;
  double XgXgoV1,RXgoV1,tBgBg1, CgCgoV1,RCgoV1,meanc, varc, muS3, muS4,lj_3,lg_4,u3,u4;
  arma::vec meane;
  arma::mat varcove;
  arma::vec muV, muS1, REoV(q1), meanAlpha, res,meanc1,RCgoV2(2);
  double muS2;
  double lambV, xi1Sq = std::pow(xi1, 2), xi2Sq = std::pow(xi2, 2);
  arma::mat XgXgoV2(L,L), varAlpha, tEEoV(q1,q1),CgCgoV2(2,2),varcovc,temp2;
  arma::rowvec RXgoV2(L);
  
  for (int t = 0; t < maxSteps; t++) {
    
    
    // alpha|
    res = y - g*hatBeta-w*arma::vectorise(hatEta)-C*hatata-xi1*hatV;
    tEEoV = (e.each_col()/hatV).t() * e;
    REoV = arma::sum(e.each_col()% (res/hatV), 0).t();
    varAlpha = arma::inv_sympd(tEEoV*hatTau/xi2Sq+invSigAlpha0);
    meanAlpha = varAlpha* REoV * hatTau / xi2Sq;
    hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
    res -= e * hatAlpha;
    gsAlpha.row(t) = hatAlpha.t();
    
    // ata|
    
    res += C.col(0) * hatata(0);
    CgCgoV1 = arma::as_scalar((C.col(0)/ hatV).t() * C.col(0));
    varc = 1/(CgCgoV1 *hatTau/xi2Sq + 1/hatSg3);
    RCgoV1 = arma::sum(C.col(0) % res/ hatV)*hatTau / xi2Sq;
    meanc = varc* RCgoV1;
    double lj_temp_3 = std::sqrt(hatSg3)*std::exp(-0.5*varc*pow(RCgoV1,2))/std::sqrt(varc);
    lj_3 = hatPi3/(hatPi3+(1-hatPi3)*lj_temp_3);
    u3 = R::runif(0, 1);
    if(u3<lj_3){
      hatata(0) = R::rnorm(meanc, sqrt(varc));
    }else{
      hatata(0) = 0;
    }
    res -= C.col(0) * hatata(0);
    res += C.cols(1,2) * hatata.subvec(1,2);
    CgCgoV2 = (C.cols(1,2).each_col()/hatV).t()*C.cols(1,2);
    temp2 = CgCgoV2*hatTau/xi2Sq;
    temp2.diag() += 1/hatSg4;
    varcovc = arma::inv(temp2);
    RCgoV2 = arma::sum(C.cols(1,2).each_col()%(res/hatV),0).t();
    RCgoV2 *= hatTau/xi2Sq;
    meanc1 = varcovc*RCgoV2;
    double lg_temp_4 = arma::as_scalar(arma::exp(-0.5*(RCgoV2.t()*varcovc*RCgoV2)))*std::sqrt(arma::det(temp2))*std::pow(hatSg4, 2/2);
    lg_4 = hatPi4/(hatPi4+(1-hatPi4)*lg_temp_4);
    u4 = R::runif(0, 1);
    if(u4<lg_4){
      hatata.subvec(1,2) = mvrnormCpp(meanc1, varcovc);
    }else{
      hatata.subvec(1,2).zeros();
    }
    res -= C.cols(1,2) * hatata.subvec(1,2);
    
    gsAta.row(t) = hatata.t();
    
    //v|
    res += xi1*hatV;
    lambV = hatTau*xi1Sq/xi2Sq + 2*hatTau;
    muV = arma::sqrt((xi1Sq+2*xi2Sq) / arma::square(res));
    for(unsigned int i=0;i<n;i++){
      bool flag = true;
      while(flag){
        hatV(i) = 1/rinvGauss(muV(i), lambV);
        if(hatV(i)<=0 || std::isinf(hatV(i)) || std::isnan(hatV(i))){
          if(progress != 0){
            Rcpp::Rcout << "hatV(i) <= 0 or nan or inf" << std::endl; 
            Rcpp::checkUserInterrupt();
          }
        }else{
          flag = false;
        }
      }
    }
    res -= xi1*hatV;
    gsV.row(t) = hatV.t();
    
    //s1|
    muS1 = std::sqrt(hatEtaSq1)/ arma::abs(hatBeta);
    for(unsigned int j = 0; j<m; j++){
      bool flag = true;
      while(flag){
        hatSg1(j) = 1/rinvGauss(muS1(j), hatEtaSq1);
        if(hatSg1(j)<=0 || std::isinf(hatSg1(j)) || std::isnan(hatSg1(j))){
          if(progress != 0) Rcpp::Rcout << "hatSg1(j): " << hatSg1(j) << std::endl; 
          Rcpp::checkUserInterrupt();
        }else{
          flag = false;
        }
      }
    }
    gsSg1.row(t) = hatSg1.t();
    
    
    //s2|
    for(unsigned int j = 0; j<m; j++){
      muS2 = std::sqrt(hatEtaSq2)/ (arma::norm(hatEta.col(j)));
      bool flag = true;
      while(flag){
        hatSg2(j) = 1/rinvGauss(muS2, hatEtaSq2);
        if(hatSg2(j)<=0 || std::isinf(hatSg2(j)) || std::isnan(hatSg2(j))){
          if(progress != 0){
            Rcpp::Rcout << "hatSg2(j) = " << hatSg2(j) << " mu: " << muS2 << " lamb: " << hatEtaSq2 << std::endl;
            Rcpp::checkUserInterrupt();
          }
        }else{
          flag = false;
        }
      }
    }
    gsSg2.row(t) = hatSg2.t();
    
    //s3|
    
    muS3 = std::sqrt(hatEtaSq3/ pow(hatata(0),2));
    
    if(hatata(0) == 0){
      hatSg3 = R::rexp(2/hatEtaSq3);
    }else{
      bool flag = true;
      while(flag){
        hatSg3 = 1/rinvGauss(muS3, hatEtaSq3);
        if(hatSg3<=0 || std::isinf(hatSg3) || std::isnan(hatSg3)){
          if(progress != 0){
            Rcpp::Rcout << "hatSg3： " << hatSg3 << std::endl; 
            Rcpp::checkUserInterrupt();
          }
        }else{
          flag = false;
        }
      }
    }
    
    gsSg3(t) = hatSg3;
    
    //s4|
    arma::vec ata;
    ata = hatata.subvec(1,2);
    tBgBg1 = arma::sum(arma::square(ata));
    if(tBgBg1 == 0){
      hatSg4 = R::rgamma((2+1)/2, 2/hatEtaSq4);
    }else{
      muS4 = std::sqrt(hatEtaSq4/tBgBg1);
      bool flag = true;
      while(flag){
        hatSg4 = 1/rinvGauss(muS4, hatEtaSq4);
        if(hatSg4<=0 || std::isinf(hatSg4) || std::isnan(hatSg4)){
          if(progress != 0){
            Rcpp::Rcout << "hatSg4: " << hatSg4 << " muS4: " << muS4 << " hatEtaSq4: " << hatEtaSq4 << std::endl; 
            Rcpp::checkUserInterrupt();
          }
        }else{
          flag = false;
        }
      }
    }
    
    gsSg4(t) = hatSg4;
    
    // Beta|
    
    for(unsigned int j=0; j<m; j++){
      res += g.col(j) * hatBeta(j);
      XgXgoV1 = arma::as_scalar((g.col(j)/hatV).t() * g.col(j));
      varb = 1/(XgXgoV1*hatTau/xi2Sq + 1/hatSg1(j));
      
      RXgoV1 = arma::sum(g.col(j) % (res/hatV));
      meanb = varb * RXgoV1 * hatTau / xi2Sq;
      hatBeta(j) = R::rnorm(meanb, sqrt(varb));
      res -= g.col(j) * hatBeta(j);
    }
    
    gsBeta.row(t) = hatBeta.t();
    
    
    // eta|
    
    for(unsigned int j=0; j<m; j++){
      res += w.cols((j*L), (j*L+L-1))* hatEta.col(j);
      XgXgoV2 = (w.cols((j*L), (j*L+L-1)).each_col()/hatV).t() * w.cols((j*L), (j*L+L-1));
      temp1 = XgXgoV2*hatTau/xi2Sq;
      temp1.diag() += 1/hatSg2(j);
      varcove = arma::inv(temp1);
      RXgoV2 = arma::sum(w.cols((j*L), (j*L+L-1)).each_col()% (res/hatV), 0);
      meane = varcove * RXgoV2.t() * hatTau / xi2Sq;
      hatEta.col(j) = mvrnormCpp(meane, varcove);
      res -= w.cols((j*L), (j*L+L-1)) * hatEta.col(j);
    }
    gseta.row(t) = arma::vectorise(hatEta).t();
    
    //etasq1;
    
    double shape2 = m+1;
    double rate2 = arma::accu(hatSg1)/2 + r1;
    hatEtaSq1 = R::rgamma(shape2, 1/rate2);
    gsEtaSq1(t) = hatEtaSq1;
    
    //etasq2;
    
    double shape21 = (m+m*L)/2+1;
    double rate21 = arma::accu(hatSg2)/2 + r2;
    hatEtaSq2 = R::rgamma(shape21, 1/rate21);
    gsEtaSq2(t) = hatEtaSq2;
    
    
    //etasq3;
    
    double shape3 = 1+1;
    double rate3 = hatSg3/2 + r3;
    hatEtaSq3 = R::rgamma(shape3, 1/rate3);
    gsEtaSq3(t) = hatEtaSq3;
    
    //etasq4;
    
    double shape4 = (1+1*2)/2+1;
    double rate4 = hatSg4/2 + r4;
    hatEtaSq4 = R::rgamma(shape4, 1/rate4);
    gsEtaSq4(t) = hatEtaSq4;
    
    //tau|
    
    double shape = a + 3*n/2;
    double ResSqoV;
    ResSqoV = arma::accu(arma::square(res)/hatV);
    double rate = b + arma::accu(hatV)+ResSqoV/(2*xi2Sq);
    hatTau = R::rgamma(shape, 1/rate);
    gsTau(t) = hatTau;
    
    //pi3|
    double shapep3 = sh1 + (hatata(0) != 0);
    double shapep4 = sh0 + (hatata(0) == 0);
    hatPi3 = R::rbeta(shapep3, shapep4);
    gsPi3(t) = hatPi3;
    
    //pi4|
    double shape13 = sh1 + (tBgBg1 != 0);
    double shape24 = sh0 + (tBgBg1 == 0);
    hatPi4 = R::rbeta(shape13, shape24);
    gsPi4(t) = hatPi4;
    
    
  }
  
  return Rcpp::List::create(
    Rcpp::Named("GS.alpha") = gsAlpha,
    Rcpp::Named("GS.beta") = gsBeta,
    Rcpp::Named("GS.eta") = gseta,
    Rcpp::Named("GS.v") = gsV,
    Rcpp::Named("GS.s1") = gsSg1,
    Rcpp::Named("GS.s2") = gsSg2,
    Rcpp::Named("GS.eta21.sq") = gsEtaSq1,
    Rcpp::Named("GS.eta22.sq") = gsEtaSq2,
    Rcpp::Named("GS.tau") = gsTau
  
  );
  
}