#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List RBV (arma::mat xx, arma::vec y, arma::mat W, arma::mat C,unsigned int s, unsigned int q, int maxSteps, arma::vec hatAlpha, arma::mat hatBeta, arma::vec hatata,double hatTau, arma::vec hatV, arma::vec hatSg, double hatSg3, double hatSg4,arma::mat invSigAlpha0, double hatPi3, double hatPi4,double hatEtaSq,double hatEtaSq3, double hatEtaSq4, double xi1, double xi2, double r,double r3, double r4, double a, double b, double sh1, double sh0,int progress)
{
  unsigned int L=q, n = xx.n_rows, q1 = W.n_cols;
  arma::mat gsAlpha(maxSteps, q1),
  gsBeta(maxSteps, s*L),
  gsV(maxSteps, n),
  gsAta(maxSteps,3),
  gsSg(maxSteps, s);
  
  arma::vec gsEtaSq(maxSteps),
  gsTau(maxSteps),
  gsEtaSq3(maxSteps),
  gsEtaSq4(maxSteps),
  gsPi3(maxSteps),
  gsPi4(maxSteps),
  gsSg3(maxSteps),
  gsSg4(maxSteps),
  
  gsMSE(maxSteps);
  
  arma::mat varAlpha, tWWoV(q1,q1), XgXgoV(L,L), temp, varG, temp2,CgCgoV2(2,2),varcovc;
  arma::vec res, RWoV(q1), meanAlpha, muV, meanG, meanc1,RCgoV2(2);
  arma::rowvec RXgoV(L);
  double lambV, xi1Sq = std::pow(xi1, 2), xi2Sq = std::pow(xi2, 2), ResSqoV, muS,tBgBg1, CgCgoV1,RCgoV1,meanc, varc, muS3, muS4,lj_3,lg_4,u3,u4;
  
  std::vector<arma::mat> Xg(s);
  for(unsigned int j=0; j<s; j++){
    Xg[j] = xx.cols((j*L), (j*L+L-1));
  }
  
  for (int k = 0; k < maxSteps; k++) {
    // Rcpp::Rcout << "alpha" << std::endl;
    res = y - xx * arma::vectorise(hatBeta) - C*hatata-xi1*hatV;
    tWWoV = (W.each_col()/hatV).t() * W;
    RWoV = arma::sum(W.each_col()% (res/hatV), 0).t();
    varAlpha = arma::inv(tWWoV*hatTau/xi2Sq + invSigAlpha0);
    meanAlpha = varAlpha * RWoV * hatTau / xi2Sq;
    hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
    res -= W * hatAlpha;
    gsAlpha.row(k) = hatAlpha.t();
    
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
    
    gsAta.row(k) = hatata.t();
    
    // Rcpp::Rcout << "v" << std::endl;
    res += xi1*hatV;
    lambV = hatTau*xi1Sq/xi2Sq + 2*hatTau;
    muV = arma::sqrt((xi1Sq+2*xi2Sq) / arma::square(res));
    for(unsigned int i = 0; i<n; i++){
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
    gsV.row(k) = hatV.t();
    
    
    // Rcpp::Rcout << "S" << std::endl;
    for(unsigned int j = 0; j<s; j++){
      muS = std::sqrt(hatEtaSq)/ (arma::norm(hatBeta.col(j)));
      bool flag = true;
      while(flag){
        hatSg(j) = 1/rinvGauss(muS, hatEtaSq);
        if(hatSg(j)<=0 || std::isinf(hatSg(j)) || std::isnan(hatSg(j))){
          if(progress != 0){
            Rcpp::Rcout << "hatSg(j) = " << hatSg(j) << " mu: " << muS << " lamb: " << hatEtaSq << std::endl;
            Rcpp::checkUserInterrupt();
          }
        }else{
          flag = false;
        }
      }
    }
    gsSg.row(k) = hatSg.t();
    
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
    
    gsSg3(k) = hatSg3;
    
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
    
    gsSg4(k) = hatSg4;
    
    // res = y - xx * arma::vectorise(hatBeta) - xi1*hatV - W * hatAlpha; //$#@!
    // Rcpp::Rcout << "beta" << std::endl;
    for(unsigned int j=0; j<s; j++){
      res += Xg[j] * hatBeta.col(j);
      XgXgoV = (Xg[j].each_col()/hatV).t() * Xg[j];
      temp = XgXgoV*hatTau/xi2Sq;
      temp.diag() += 1/hatSg(j);
      varG = arma::inv(temp);
      
      RXgoV = arma::sum(Xg[j].each_col()% (res/hatV), 0);
      meanG = varG * RXgoV.t() * hatTau / xi2Sq;
      hatBeta.col(j) = mvrnormCpp(meanG, varG);
      res -= Xg[j] * hatBeta.col(j);
    }
    gsBeta.row(k) = arma::vectorise(hatBeta).t();
    
    
    // res = y - xx * arma::vectorise(hatBeta) - xi1*hatV - W * hatAlpha; //$#@!
    // Rcpp::Rcout << "tau" << std::endl;
    double shape = a + 3*n/2;
    ResSqoV = arma::accu(arma::square(res)/hatV);
    double rate = b + arma::accu(hatV) + ResSqoV/(2*xi2Sq);
    hatTau = R::rgamma(shape, 1/rate);
    gsTau(k) = hatTau;
    
    
    // Rcpp::Rcout << "eta2Sq" << std::endl;
    double shape2 = (s+s*L)/2+1;
    double rate2 = arma::accu(hatSg)/2 + r;
    hatEtaSq = R::rgamma(shape2, 1/rate2);
    gsEtaSq(k) = hatEtaSq;
    
    //etasq3;
    
    double shape3 = 1+1;
    double rate3 = hatSg3/2 + r3;
    hatEtaSq3 = R::rgamma(shape3, 1/rate3);
    gsEtaSq3(k) = hatEtaSq3;
    
    //etasq4;
    
    double shape4 = (1+1*2)/2+1;
    double rate4 = hatSg4/2 + r4;
    hatEtaSq4 = R::rgamma(shape4, 1/rate4);
    gsEtaSq4(k) = hatEtaSq4;
    
    //pi3|
    double shapep3 = sh1 + (hatata(0) != 0);
    double shapep4 = sh0 + (hatata(0) == 0);
    hatPi3 = R::rbeta(shapep3, shapep4);
    gsPi3(k) = hatPi3;
    
    //pi4|
    double shape13 = sh1 + (tBgBg1 != 0);
    double shape24 = sh0 + (tBgBg1 == 0);
    hatPi4 = R::rbeta(shape13, shape24);
    gsPi4(k) = hatPi4;
    
    
    gsMSE(k) = arma::mean(arma::abs(res));
    if(k % 100 == 0){
      Rcpp::checkUserInterrupt();
    }
    if(progress != 0 && k % progress == 0){
      Rcpp::Rcout << "\nIter." << k << "  MAD: " << gsMSE(k) << std::endl;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("GS.alpha") = gsAlpha,
                            Rcpp::Named("GS.beta") = gsBeta,
                            Rcpp::Named("GS.tau") = gsTau,
                            Rcpp::Named("GS.v") = gsV,
                            Rcpp::Named("GS.s") = gsSg,
                            Rcpp::Named("GS.eta2.sq") = gsEtaSq,
                            Rcpp::Named("GS.mad") = gsMSE);
}