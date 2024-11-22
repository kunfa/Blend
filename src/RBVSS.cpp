#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List RBVSS (arma::mat xx, arma::vec y, arma::mat W, arma:: mat C, unsigned int s, unsigned int q, int maxSteps, arma::vec hatAlpha, arma::mat hatBeta, arma::vec hatata, double hatTau, arma::vec hatV, arma::vec hatSg, double hatSg3, double hatSg4, arma::mat invSigAlpha0, double hatPi, double hatPi3, double hatPi4,double hatEtaSq,double hatEtaSq3, double hatEtaSq4, double xi1, double xi2, double r, double r3, double r4, double a, double b, double sh1, double sh0, int progress)
{
  unsigned int n = xx.n_rows, q1 = W.n_cols;
  arma::mat gsAlpha(maxSteps, q1),
  gsBeta(maxSteps, s*q),
  gsV(maxSteps, n),
  gsAta(maxSteps,3),
  gsRstRs(maxSteps, s),
  gsSg(maxSteps, s),
  gsLg(maxSteps, s),
  gsSg3(maxSteps, 1),
  gsSg4(maxSteps, 1) 
  ;
  
  arma::vec gsEtaSq(maxSteps),
  gsEtaSq3(maxSteps),
  gsEtaSq4(maxSteps),
  gsPi(maxSteps),
  gsPi3(maxSteps),
  gsPi4(maxSteps),
  gsTau(maxSteps),
  gsMSE(maxSteps),
  idgene(s);
  
  idgene.zeros();
  
  arma::mat varAlpha, tWWoV(q1,q1), XgXgoV(q,q), temp, temp1, varG, varcovc, CgCgoV2(2,2),temp2;
  arma::vec res, RWoV(q1), RXgoV(q), meanAlpha, muV, meanG, tBgBg, meanc1, RCgoV2(2);
  double lambV, xi1Sq = std::pow(xi1, 2), xi2Sq = std::pow(xi2, 2), ResSqoV, muS, lg, u,meanc,varc, tBgBg1, CgCgoV1,RCgoV1,muS3, muS4, lj_3, lg_4, u3, u4;
  
  for (int t = 0; t < maxSteps; t++) {
    // Rcpp::Rcout << "alpha" << std::endl;
    res = y - xx * arma::vectorise(hatBeta) - C*hatata- xi1*hatV;
    tWWoV = (W.each_col()/hatV).t() * W;
    RWoV = arma::sum(W.each_col()% (res/hatV), 0).t();
    varAlpha = arma::inv_sympd(tWWoV*hatTau/xi2Sq + invSigAlpha0);
    meanAlpha = varAlpha * RWoV * hatTau / xi2Sq;
    hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
    res -= W * hatAlpha;
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
    gsV.row(t) = hatV.t();
    
    
    // Rcpp::Rcout << "S" << std::endl;
    tBgBg = arma::sum(arma::square(hatBeta), 0).t();
    for(unsigned int j = 0; j<s; j++){
      if(tBgBg(j) == 0){
        hatSg(j) = R::rgamma((q+1)/2, 2/hatEtaSq);
      }else{
        muS = std::sqrt(hatEtaSq/tBgBg(j));
        bool flag = true;
        while(flag){
          hatSg(j) = 1/rinvGauss(muS, hatEtaSq);
          if(hatSg(j)<=0 || std::isinf(hatSg(j)) || std::isnan(hatSg(j))){
            if(progress != 0){
              Rcpp::Rcout << "hatSg(j): " << hatSg(j) << " muS: " << muS << " hatEtaSq: " << hatEtaSq << std::endl; 
              Rcpp::checkUserInterrupt();
            }
          }else{
            flag = false;
          }
        }
      }
    }
    gsSg.row(t) = hatSg.t();
    gsRstRs.row(t) = tBgBg.t();
    
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
    
    gsSg3.row(t) = hatSg3;
    
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
    
    gsSg4.row(t) = hatSg4;
    
    // Rcpp::Rcout << "beta" << std::endl;
    for(unsigned int j=0; j<s; j++){
      res += xx.cols(j*q, j*q+q-1) * hatBeta.col(j);
      
      XgXgoV = (xx.cols(j*q, j*q+q-1).each_col()/hatV).t() * xx.cols(j*q, j*q+q-1);
      temp = XgXgoV*hatTau/xi2Sq;
      temp.diag() += 1/hatSg(j);
      varG = arma::inv(temp);
      
      RXgoV = arma::sum(xx.cols(j*q, j*q+q-1).each_col()% (res/hatV), 0).t();
      RXgoV *=  hatTau/xi2Sq;
      meanG = varG * RXgoV;
      
      double lg_temp = arma::as_scalar(arma::exp(-0.5*(RXgoV.t()*varG*RXgoV)))*std::sqrt(arma::det(temp))*std::pow(hatSg(j), q/2);
      lg = hatPi/(hatPi+(1-hatPi)*lg_temp);
      u = R::runif(0, 1);
      if(u<lg){
        hatBeta.col(j) = mvrnormCpp(meanG, varG);
         idgene(j)++;
      }else{
        hatBeta.col(j).zeros();
      }
      res -= xx.cols(j*q, j*q+q-1) * hatBeta.col(j);
    }
    gsBeta.row(t) = arma::vectorise(hatBeta).t();
    
    // Rcpp::Rcout << "pi" << std::endl;
    double shape1 = sh1 + arma::accu(tBgBg != 0);
    double shape2 = sh0 + arma::accu(tBgBg == 0);
    hatPi = R::rbeta(shape1, shape2);
    gsPi(t) = hatPi;
    
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
   
    
    // Rcpp::Rcout << "tau" << std::endl;
    ResSqoV = arma::accu(arma::square(res)/hatV);
    double rate = b + arma::accu(hatV) + ResSqoV/(2*xi2Sq);
    hatTau = R::rgamma(a+3*n/2, 1/rate);
    gsTau(t) = hatTau;
    
    // Rcpp::Rcout << "eta2Sq" << std::endl;
    double rate2 = arma::accu(hatSg)/2 + r;
    hatEtaSq = R::rgamma((s+s*q)/2+1, 1/rate2);
    gsEtaSq(t) = hatEtaSq;
    
    //etasq3;
    
    double shape3 = 1+1;
    double rate3 = arma::accu(hatSg3)/2 + r3;
    hatEtaSq3 = R::rgamma(shape3, 1/rate3);
    gsEtaSq3(t) = hatEtaSq3;
    
    //etasq4;
    
    double shape4 = (1+1*2)/2+1;
    double rate4 = arma::accu(hatSg4)/2 + r4;
    hatEtaSq4 = R::rgamma(shape4, 1/rate4);
    gsEtaSq4(t) = hatEtaSq4;
    
    gsMSE(t) = arma::mean(arma::abs(res));
    if(t % 100 == 0){
      Rcpp::checkUserInterrupt();
    }
    if(progress != 0 && t % progress == 0){
      Rcpp::Rcout << "\nIter." << t << "  MAD: " << gsMSE(t) << std::endl;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("GS.alpha") = gsAlpha,
                            Rcpp::Named("GS.ata") = gsAta,
                            Rcpp::Named("GS.beta") = gsBeta,
                            Rcpp::Named("GS.tau") = gsTau,
                            Rcpp::Named("GS.v") = gsV,
                            Rcpp::Named("GS.s") = gsSg,
                            Rcpp::Named("GS.s3") = gsSg3,
                            Rcpp::Named("GS.tRsRs") = gsRstRs,
                            Rcpp::Named("GS.s4") = gsSg4,
                            Rcpp::Named("GS.pi") = gsPi,
                            Rcpp::Named("GS.pi3") = gsPi3,
                            Rcpp::Named("GS.pi4") = gsPi4,
                            Rcpp::Named("GS.eta2.sq") = gsEtaSq,
                            Rcpp::Named("GS.eta23.sq") = gsEtaSq3,
                            Rcpp::Named("GS.eta24.sq") = gsEtaSq4,
                            Rcpp::Named("gsLg") = gsLg
                            );
}