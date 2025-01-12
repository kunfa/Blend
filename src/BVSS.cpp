#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include<vector>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BVSS (arma::vec y, arma::mat W, arma::mat xx,arma::mat C, unsigned int s, unsigned int q, int maxSteps, arma::vec hatAlpha, arma::vec hatBeta, arma::vec hatata, arma::vec hatInvTauSqStar,double hatInvTauSq1, double hatInvTauSq2, arma::mat invSigAlpha0, double hatPiStar, double hatPi, double hatPi1,double hatLambdaSqStar,double hatLambdaSq1, double hatLambdaSq2, double hatSigmaSq, double aStar, double bStar,double a1, double b1, double a2, double b2, double alpha, double gamma, double sh1, double sh0, double mu1, double nu1, double mu2, double nu2,int progress)
{
  unsigned int L = q, n = xx.n_rows,q1 = W.n_cols;
  arma::mat gsAlpha(maxSteps, q1),
  gsBeta(maxSteps, s*q),
  gsRstRs(maxSteps, s),
  gsInvTauSqStar(maxSteps, s),
  gsAta(maxSteps,3),
  gsLS(maxSteps, s);

  arma::vec gsLambdaStar(maxSteps),
  gsLambda1(maxSteps),
  gsLambda2(maxSteps),
  gsSigmaSq(maxSteps),
  gsInvTauSq1(maxSteps),
  gsInvTauSq2(maxSteps),
  gsPiStar(maxSteps),
  gsPi(maxSteps),
  gsPi1(maxSteps),
  idgene(s),
  gsMSE(maxSteps);
  idgene.zeros();
  arma::mat Br = xx;
  arma::mat tWW = W.t()*W;

  arma::mat Xr, varM, varAlpha, varRs, tempS, matRStar,temp1, vara1;
  arma::vec res, BrjtRes, meanAlpha, meanRs, tRsRs, repInvTauSq, muInvTauSqStar,CrjtRes, meana1;
  double lS, u, l,l1,lInvTauSqStar,lInvTauSq1, lInvTauSq2, muInvTauSq1, tR1R1, muInvTauSq2;

  std::vector<arma::mat> tBrBr(s);
  for(unsigned int j=0; j<s; j++){
    Xr = Br.cols((j*L), (j*L+L-1));
    tBrBr[j] = Xr.t()*Xr;
  }
  arma:: mat tCrCr,Cr;
  Cr = C.cols(1,2);
  tCrCr = Cr.t()*Cr;

  for (int t = 0; t < maxSteps; t++) {
    // alpha|
    varAlpha = arma::inv(tWW/hatSigmaSq + invSigAlpha0);
    res = y - Br * hatBeta - C*hatata;
    meanAlpha = varAlpha * (W.t() * res/hatSigmaSq);
    hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
    res -= W * hatAlpha;
    gsAlpha.row(t) = hatAlpha.t();

    // ata|

    double tBB, temp, vara, BjtRes, meana;
    tBB = arma::sum(arma::square(C.col(0)));
    temp = 1/(tBB + hatInvTauSq1);
    vara = hatSigmaSq * temp;
    res += C.col(0) * hatata(0);
    BjtRes = arma::as_scalar(C.col(0).t() * res);
    meana = temp * BjtRes;
    l = hatPi/(hatPi + (1-hatPi)*sqrt(hatInvTauSq1*temp)*exp(0.5/hatSigmaSq*temp*pow(BjtRes,2)));
    u = R::runif(0, 1);
    if(u<l){
      hatata(0) = 0;
    }else{
      hatata(0) = R::rnorm(meana, sqrt(vara));
    }
    res -= C.col(0) * hatata(0);

    temp1 = tCrCr;
    temp1.diag() += hatInvTauSq2;
    temp1 = arma::inv(temp1);
    vara1 = hatSigmaSq * temp1;
    res += C.cols(1,2) * hatata.subvec(1,2);
    CrjtRes = C.cols(1,2).t() * res;
    meana1 = temp1 * CrjtRes;
    l1 = arma::as_scalar(hatPi1/(hatPi1+(1-hatPi1)*pow(hatInvTauSq2,(2/2))*sqrt(arma::det(temp1))*exp(0.5/hatSigmaSq*(CrjtRes.t()*temp1*CrjtRes))));
    u = R::runif(0, 1);
    if(u<l1){
      hatata.subvec(1,2).zeros();
    }else{
      hatata.subvec(1,2) = mvrnormCpp(meana1, vara1);
    }
    res -= C.cols(1,2) * hatata.subvec(1,2);
    gsAta.row(t) = hatata.t();

    // beta|
    for(unsigned int j=0; j<s; j++){
      tempS = tBrBr[j];
      tempS.diag() += hatInvTauSqStar(j);
      tempS = arma::inv(tempS);
      varRs = hatSigmaSq * tempS;
      res += Br.cols((j*L), (j*L+L-1)) * hatBeta.subvec((j*L), (j*L+L-1));
      BrjtRes = Br.cols((j*L), (j*L+L-1)).t() * res;
      meanRs = tempS * BrjtRes;
      lS = arma::as_scalar(hatPiStar/(hatPiStar+(1-hatPiStar)*std::pow(hatInvTauSqStar(j),(L/2))*std::sqrt(arma::det(tempS))*arma::exp(0.5/hatSigmaSq*(BrjtRes.t()*tempS*BrjtRes))));
      u = R::runif(0, 1);
      if(u<lS){
        hatBeta.subvec((j*L), (j*L+L-1)).zeros();
      }else{
        hatBeta.subvec((j*L), (j*L+L-1)) = mvrnormCpp(meanRs, varRs);
        idgene(j)++;
      }
      res -= Br.cols((j*L), (j*L+L-1)) * hatBeta.subvec((j*L), (j*L+L-1));
    }
    gsBeta.row(t) = hatBeta.t();



    // invTAUsq.star|
    lInvTauSqStar = L * hatLambdaSqStar;
    matRStar = arma::reshape(hatBeta, L, s);
    tRsRs = arma::sum(arma::square(matRStar), 0).t();
    muInvTauSqStar = arma::sqrt(L * hatLambdaSqStar * hatSigmaSq / tRsRs);
    for(unsigned int j = 0; j<s; j++){
      if(tRsRs(j) == 0){
        hatInvTauSqStar(j) = 1/R::rgamma((L+1)/2, 2/lInvTauSqStar);
      }else{
        hatInvTauSqStar(j) = rinvgaussian(muInvTauSqStar(j), lInvTauSqStar);
      }
    }
    gsInvTauSqStar.row(t) = hatInvTauSqStar.t();
    gsRstRs.row(t) = tRsRs.t();

    // invTAUsq.0|lambda, r0
    lInvTauSq1 = hatLambdaSq1;
    muInvTauSq1 = sqrt(hatLambdaSq1 * hatSigmaSq / pow(hatata(0),2));
    if(hatata(0) == 0){
      hatInvTauSq1 = 1/R::rgamma(1, 2/lInvTauSq1);
    }else{
      hatInvTauSq1 = rinvgaussian(muInvTauSq1, lInvTauSq1);
    }
    gsInvTauSq1(t) = hatInvTauSq1;

    // invTAUsq.star|lambda.star, r.star
    lInvTauSq2 = 2 * hatLambdaSq2;
    arma::vec ata;
    ata = hatata.subvec(1,2);
    tR1R1 = arma::sum(arma::square(ata));
    muInvTauSq2 = std::sqrt(2 * hatLambdaSq2 * hatSigmaSq / tR1R1);
    if(tR1R1 == 0){
      hatInvTauSq2 = 1/R::rgamma((2+1)/2, 2/lInvTauSq2);
    }else{
      hatInvTauSq2 = rinvgaussian(muInvTauSq2, lInvTauSq2);
    }
    gsInvTauSq2(t) = hatInvTauSq2;

    // sigma.sq|
    double shapeSig = alpha + n/2 + L*arma::accu(tRsRs != 0)/2+2*(tR1R1 != 0)/2+(hatata(0) != 0)/2 ;
    repInvTauSq = arma::vectorise(arma::repelem(hatInvTauSqStar.t(), L, 1), 0);
    double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) +
                                  arma::accu(square(hatBeta) % repInvTauSq))+pow(hatata(0),2)*hatInvTauSq1+arma::accu(square(ata) * hatInvTauSq2);
    hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
    gsSigmaSq(t) = hatSigmaSq;

    // lambda.star|
    double shapeS = aStar + s*(L+1)/2;
    double rateS = bStar + L*arma::accu(1/hatInvTauSqStar)/2;
    hatLambdaSqStar = R::rgamma(shapeS, 1/rateS);
    gsLambdaStar(t) = hatLambdaSqStar;

    // lambda0|invTAUsq.0
    double shape1 = a1 + 1;
    double rate1 = b1 + (1/hatInvTauSq1)/2;
    hatLambdaSq1 = R::rgamma(shape1, 1/rate1);
    gsLambda1(t) = hatLambdaSq1;


    // lambda.star|invTAUsq.star
    double shape2 = a2 + 1*(2+1)/2;
    double rate2 = b2 + 2*(1/hatInvTauSq2)/2;
    hatLambdaSq2 = R::rgamma(shape2, 1/rate2);
    gsLambda2(t) = hatLambdaSq2;


    // pi.star|
    double shape1_s = sh0 + arma::accu(tRsRs == 0);
    double shape2_s = sh1 + arma::accu(tRsRs != 0);
    hatPiStar = R::rbeta(shape1_s, shape2_s);
    gsPiStar(t) = hatPiStar;

    // pi.0|
    double shape1_1 = mu1 + (hatata(0) == 0);
    double shape2_1 = nu1 + (hatata(0) != 0);
    hatPi = R::rbeta(shape1_1, shape2_1);
    gsPi(t) = hatPi;


    // pi.star|
    double shape1_2 = mu2 + (tR1R1 == 0);
    double shape2_2 = nu2 + (tR1R1 != 0);
    hatPi1 = R::rbeta(shape1_2, shape2_2);
    gsPi1(t) = hatPi1;

    gsMSE(t) = arma::mean(arma::square(res));
    if(t % 100 == 0){
      Rcpp::checkUserInterrupt();
    }
    if(progress != 0 && t % progress == 0){
      Rcpp::Rcout << "\nIter." << t << "  mse: " << gsMSE(t) << std::endl;
    }
  }
  return Rcpp::List::create(	Rcpp::Named("GS.alpha") = gsAlpha,
                             Rcpp::Named("GS.ata") = gsAta,
                             Rcpp::Named("GS.beta") = gsBeta,
                             Rcpp::Named("GS.tRsRs") = gsRstRs,
                             Rcpp::Named("GS.invTAUsq") = gsInvTauSqStar,
                             Rcpp::Named("GS.pi") = gsPiStar,
                             Rcpp::Named("GS.lambda.sq") = gsLambdaStar,
                             Rcpp::Named("GS.sigma.sq") = gsSigmaSq,
                             Rcpp::Named("idgene")=idgene,
                             Rcpp::Named("GS.lS") = gsLS);
}
