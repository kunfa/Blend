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

Rcpp::List RBVSS_SI(arma::vec y, arma::mat e, arma::mat g, arma:: mat w, arma:: mat C, int maxSteps, unsigned int q, arma::vec hatBeta, arma:: mat hatEta, arma::vec hatAlpha, arma:: vec hatata, double hatTau, arma::vec hatV, arma::vec hatSg1,arma::vec hatSg2, double hatSg3, double hatSg4, arma::mat invSigAlpha0,double hatPi1,double hatPi2, double hatPi3, double hatPi4, double hatEtaSq1, double hatEtaSq2, double hatEtaSq3, double hatEtaSq4, double xi1, double xi2, double r1,double r2, double r3, double r4,double a, double b,double sh1, double sh0, int progress)
{
  unsigned int n = g.n_rows, L = q-1, m = g.n_cols, p = w.n_cols,q1 = e.n_cols;
  arma::mat gsAlpha(maxSteps, q1),
  gsBeta(maxSteps,m),
  gseta(maxSteps,p),
  gsV(maxSteps, n),
  gsAta(maxSteps,3),
  gsRstRs(maxSteps, m),
  gsSg1(maxSteps, m),
  gsLg1(maxSteps, m),
  gsLg2(maxSteps, m),
  gsSg2(maxSteps, m),
  gsSg3(maxSteps, 1),
  gsSg4(maxSteps, 1)
    ;
  
  arma::vec gsEtaSq1(maxSteps),
  gsEtaSq2(maxSteps),
  gsEtaSq3(maxSteps),
  gsEtaSq4(maxSteps),
  gsTau(maxSteps),
  gsPi1(maxSteps),
  gsPi3(maxSteps),
  gsPi4(maxSteps),
  gsPi2(maxSteps),
  idgene(m),
  idgene1(m)
    ;
  idgene.zeros();
  idgene1.zeros();
  arma::mat temp,temp1, temp2;
  
  double meanb;
  double varb;
  double meanc;
  double varc;
  double XgXgoV1,RXgoV1,tBgBg1, CgCgoV1,RCgoV1;
  arma::vec meane;
  arma::mat varcove,matRStar;
  arma::vec meanc1;
  arma::mat varcovc;
  arma::vec muV, muS1,tBgBg, REoV(q1), meanAlpha, res;
  double muS2, muS3, muS4;
  double lambV, xi1Sq = std::pow(xi1, 2), xi2Sq = std::pow(xi2, 2),lj_1, lg_2, lj_3, lg_4,u1,u2, u3, u4;
  arma::mat XgXgoV2(L,L), varAlpha,varAlpha1, tEEoV(q1,q1), CgCgoV2(2,2);
  arma::vec RXgoV2(L),RCgoV2(2);
  
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
      if(hatBeta(j) == 0){
        hatSg1(j) = R::rexp(2/hatEtaSq1);
      }else{
        bool flag = true;
        while(flag){
          hatSg1(j) = 1/rinvGauss(muS1(j), hatEtaSq1);
          if(hatSg1(j)<=0 || std::isinf(hatSg1(j)) || std::isnan(hatSg1(j))){
            if(progress != 0){
              Rcpp::Rcout << "hatSg1(j)： " << hatSg1(j) << std::endl; 
              Rcpp::checkUserInterrupt();
            }
          }else{
            flag = false;
          }
        }
      }
      
    }
    gsSg1.row(t) = hatSg1.t();
    
    
    //s2|
    tBgBg = arma::sum(arma::square(hatEta), 0).t();
    for(unsigned int j = 0; j<m; j++){
      if(tBgBg(j) == 0){
        hatSg2(j) = R::rgamma((L+1)/2, 2/hatEtaSq2);
      }else{
        muS2 = std::sqrt(hatEtaSq2/tBgBg(j));
        bool flag = true;
        while(flag){
          hatSg2(j) = 1/rinvGauss(muS2, hatEtaSq2);
          if(hatSg2(j)<=0 || std::isinf(hatSg2(j)) || std::isnan(hatSg2(j))){
            if(progress != 0){
              Rcpp::Rcout << "hatSg2(j): " << hatSg2(j) << " muS2: " << muS2 << " hatEtaSq2: " << hatEtaSq2 << std::endl; 
              Rcpp::checkUserInterrupt();
            }
          }else{
            flag = false;
          }
        }
      }
    }
    gsSg2.row(t) = hatSg2.t();
    
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
    
    // Beta|
    
    for(unsigned int j=0;j<m;j++){
      res += g.col(j) * hatBeta(j);
      XgXgoV1 = arma::as_scalar((g.col(j)/ hatV).t() * g.col(j));
      varb = 1/(XgXgoV1 *hatTau/xi2Sq + 1/hatSg1(j));
      RXgoV1 = arma::sum(g.col(j) % res/ hatV)*hatTau / xi2Sq;
      meanb = varb* RXgoV1;
      double lj_temp_1 = std::sqrt(hatSg1(j))*std::exp(-0.5*varb*pow(RXgoV1,2))/std::sqrt(varb);
      lj_1 = hatPi1/(hatPi1+(1-hatPi1)*lj_temp_1);
      gsLg1(t, j) = lj_1;
      u1 = R::runif(0, 1);
      if(u1<lj_1){
        hatBeta(j) = R::rnorm(meanb, sqrt(varb));
        idgene(j)++;
      }else{
        hatBeta(j) = 0;
      }
      res -= g.col(j) * hatBeta(j);
    }
    
    gsBeta.row(t) = hatBeta.t();
    
    
    // eta|
    
    for(unsigned int j=0;j<m;j++){
      res += w.cols(j*L, j*L+L-1) * hatEta.col(j);
      XgXgoV2 = (w.cols((j*L),(j*L+L-1)).each_col()/hatV).t()*w.cols((j*L),(j*L+L-1));
      temp = XgXgoV2*hatTau/xi2Sq;
      temp.diag() += 1/hatSg2(j);
      varcove = arma::inv(temp);
      RXgoV2 = arma::sum(w.cols((j*L),(j*L+L-1)).each_col()%(res/hatV),0).t();
      RXgoV2 *= hatTau/xi2Sq;
      meane = varcove*RXgoV2;
      double lg_temp_2 = arma::as_scalar(arma::exp(-0.5*(RXgoV2.t()*varcove*RXgoV2)))*std::sqrt(arma::det(temp))*std::pow(hatSg2(j), L/2);
      lg_2 = hatPi2/(hatPi2+(1-hatPi2)*lg_temp_2);
      gsLg2(t, j) = lg_2;
      u2 = R::runif(0, 1);
      if(u2<lg_2){
        hatEta.col(j) = mvrnormCpp(meane, varcove);
        idgene1(j)++;
      }else{
        hatEta.col(j).zeros();
      }
      res -= w.cols(j*L, j*L+L-1) * hatEta.col(j);
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
    double rate3 = arma::accu(hatSg3)/2 + r3;
    hatEtaSq3 = R::rgamma(shape3, 1/rate3);
    gsEtaSq3(t) = hatEtaSq3;
    
    //etasq4;
    
    double shape4 = (1+1*2)/2+1;
    double rate4 = arma::accu(hatSg4)/2 + r4;
    hatEtaSq4 = R::rgamma(shape4, 1/rate4);
    gsEtaSq4(t) = hatEtaSq4;
    

    //tau|
    
    double shape = a + 3*n/2;
    double ResSqoV;
    ResSqoV = arma::accu(arma::square(res)/hatV);
    double rate = b + arma::accu(hatV)+ResSqoV/(2*xi2Sq);
    hatTau = R::rgamma(shape, 1/rate);
    gsTau(t) = hatTau;
    
    //pi1|
    double shapep1 = sh1 + arma::accu(hatBeta != 0);
    double shapep2 = sh0 + arma::accu(hatBeta == 0);
    hatPi1 = R::rbeta(shapep1, shapep2);
    gsPi1(t) = hatPi1;
    
    //pi2|
    double shape12 = sh1 + arma::accu(tBgBg != 0);
    double shape22 = sh0 + arma::accu(tBgBg == 0);
    hatPi2 = R::rbeta(shape12, shape22);
    gsPi2(t) = hatPi2;
    
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
    Rcpp::Named("GS.ata") = gsAta,
    Rcpp::Named("GS.v") = gsV,
    Rcpp::Named("GS.tRsRs") = gsRstRs,
    Rcpp::Named("GS.s1") = gsSg1,
    Rcpp::Named("GS.s2") = gsSg2,
    Rcpp::Named("GS.s3") = gsSg3,
    Rcpp::Named("GS.s4") = gsSg4,
    Rcpp::Named("GS.eta21.sq") = gsEtaSq1,
    Rcpp::Named("GS.eta22.sq") = gsEtaSq2,
    Rcpp::Named("GS.eta23.sq") = gsEtaSq3,
    Rcpp::Named("GS.eta24.sq") = gsEtaSq4,
    Rcpp::Named("GS.tau") = gsTau,
    Rcpp::Named("GS.pi1") = gsPi1,
    Rcpp::Named("GS.pi2") = gsPi2,
    Rcpp::Named("GS.pi3") = gsPi3
      );
  
}