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
Rcpp::List BVSS_SI (arma::vec y, arma::mat e, arma::mat g, arma:: mat w, arma:: mat C, unsigned int s, unsigned int q, int maxSteps, arma::vec hatM, arma::vec hatR0, arma::vec hatata, arma::vec hatRStar, arma::vec hatInvSigM0, arma::vec hatInvTauSq0, arma::vec hatInvTauSqStar, double hatInvTauSq1, double hatInvTauSq2, double hatPi0, double hatPiStar,double hatPi, double hatPi1, double hatLambdaSq0, double hatLambdaSqStar, double hatLambdaSq1, double hatLambdaSq2, double hatSigmaSq, double a0, double b0, double aStar, double bStar, double a1, double b1, double a2, double b2, double alpha, double gamma, double mu0, double muStar, double mu1, double mu2, double nu0, double nuStar, double nu1, double nu2, int progress)
{
  unsigned int L = q-1, n = g.n_rows,q1 = e.n_cols;
  arma::mat gsM(maxSteps, q1),
  gsR0(maxSteps, s),
  gsRStar(maxSteps, s*(q-1)),
  gsRstRs(maxSteps, s),
  gsInvTauSq0(maxSteps, s),
  gsAta(maxSteps,3),
  gsInvTauSqStar(maxSteps, s),

  gsL0(maxSteps, s),
  gsLS(maxSteps, s);

  arma::vec gsLambda0(maxSteps),
  gsLambdaStar(maxSteps),
  gsLambda1(maxSteps),
  gsLambda2(maxSteps),
  gsSigmaSq(maxSteps),
  gsPi0(maxSteps),
  gsInvTauSq1(maxSteps),
  gsInvTauSq2(maxSteps),
  gsPiStar(maxSteps),
  gsPi(maxSteps),
  gsPi1(maxSteps),
  idgene(s),
  idgene1(s)
;
  idgene.zeros();
  idgene1.zeros();
  arma::mat tBmBm = e.t()*e, tB0B0 = g.t()*g;
  arma::vec tB0B0Diag = tB0B0.diag();
  arma::mat invSigM0 = arma::diagmat(hatInvSigM0);

  arma::mat Xr, varM, varRs, tempS, matRStar, temp1, vara1;
  arma::vec res, meanM, BrjtRes, meanRs, tRsRs, repInvTauSq, residual, muInvTauSq0, muInvTauSqStar, CrjtRes, meana1 ; // mu_m, mu_alpha,
  double temp0, varR0, B0jtRes, meanR0, l0, lS, l,l1,u, lInvTauSq0, lInvTauSqStar, lInvTauSq1, lInvTauSq2, muInvTauSq1, tR1R1, muInvTauSq2;

  std::vector<arma::mat> tBrBr(s);
  for(unsigned int j=0; j<s; j++){
    Xr = w.cols((j*L), (j*L+L-1));
    tBrBr[j] = Xr.t()*Xr;
  }
  arma:: mat tCrCr,Cr;
  Cr = C.cols(1,2);
  tCrCr = Cr.t()*Cr;

  for (int t = 0; t < maxSteps; t++) {
    // m|y, r0, r.star
    varM = arma::inv(tBmBm/hatSigmaSq + invSigM0);
    res = y - (g * hatR0 + w * hatRStar + C*hatata);
    meanM = varM * (e.t() * res/hatSigmaSq);
    hatM = mvrnormCpp(meanM, varM);
    res -= e * hatM;
    gsM.row(t) = hatM.t();

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

    for(unsigned int j=0; j<s; j++){
      temp0 = 1/(tB0B0Diag(j) + hatInvTauSq0(j));
      varR0 = hatSigmaSq * temp0;
      res += g.col(j) * hatR0(j);
      B0jtRes = arma::as_scalar(g.col(j).t() * res);
      meanR0 = temp0 * B0jtRes;
      l0 = hatPi0/(hatPi0 + (1-hatPi0)*sqrt(hatInvTauSq0(j)*temp0)*exp(0.5/hatSigmaSq*temp0*pow(B0jtRes,2)));
      u = R::runif(0, 1);
      if(u<l0){
        hatR0(j) = 0;
      }else{
        hatR0(j) = R::rnorm(meanR0, sqrt(varR0));
        idgene(j)++;
      }
      res -= g.col(j) * hatR0(j);
    }
    gsR0.row(t) = hatR0.t();

    for(unsigned int j=0; j<s; j++){
      tempS = tBrBr[j];
      tempS.diag() += hatInvTauSqStar(j);
      tempS = arma::inv(tempS);
      varRs = hatSigmaSq * tempS;
      res += w.cols((j*L), (j*L+L-1)) * hatRStar.subvec((j*L), (j*L+L-1));
      BrjtRes = w.cols((j*L), (j*L+L-1)).t() * res;
      meanRs = tempS * BrjtRes;
      lS = arma::as_scalar(hatPiStar/(hatPiStar+(1-hatPiStar)*pow(hatInvTauSqStar(j),(L/2))*sqrt(arma::det(tempS))*exp(0.5/hatSigmaSq*(BrjtRes.t()*tempS*BrjtRes))));
      u = R::runif(0, 1);
      if(u<lS){
        hatRStar.subvec((j*L), (j*L+L-1)).zeros();
      }else{
        hatRStar.subvec((j*L), (j*L+L-1)) = mvrnormCpp(meanRs, varRs);
        idgene1(j)++;
      }
      res -= w.cols((j*L), (j*L+L-1)) * hatRStar.subvec((j*L), (j*L+L-1));
    }

    gsRStar.row(t) = hatRStar.t();

    // invTAUsq.0|lambda, r0
    lInvTauSq0 = hatLambdaSq0;
    muInvTauSq0 = sqrt(hatLambdaSq0 * hatSigmaSq / square(hatR0));
    for(unsigned int j = 0; j < s; j++){
      if(hatR0(j) == 0){
        hatInvTauSq0(j) = 1/R::rgamma(1, 2/lInvTauSq0);
      }else{
        hatInvTauSq0(j) = rinvgaussian(muInvTauSq0(j), lInvTauSq0);
      }
    }
    gsInvTauSq0.row(t) = hatInvTauSq0.t();


    // invTAUsq.star|lambda.star, r.star
    lInvTauSqStar = L * hatLambdaSqStar;
    matRStar = arma::reshape(hatRStar, L, s);
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


    // lambda0|invTAUsq.0
    double shape = a0 + s;
    double rate = b0 + arma::accu(1/hatInvTauSq0)/2;
    hatLambdaSq0 = R::rgamma(shape, 1/rate);
    gsLambda0(t) = hatLambdaSq0;


    // lambda.star|invTAUsq.star
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

    // sigma.sq|
    double shapeSig = alpha + n/2 + L*arma::accu(tRsRs != 0)/2 + arma::accu(hatR0 != 0)/2+ 2*(tR1R1 != 0)/2+(hatata(0) != 0)/2 ;
    repInvTauSq = arma::vectorise(arma::repelem(hatInvTauSqStar.t(), L, 1), 0);
    double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) +
                                  arma::accu(square(hatR0) % hatInvTauSq0) +
                                  arma::accu(square(hatRStar) % repInvTauSq))+
                                  pow(hatata(0),2)*hatInvTauSq1+arma::accu(square(ata) * hatInvTauSq2);
    hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
    gsSigmaSq(t) = hatSigmaSq;


    // pi.0|
    double shape1_0 = mu0 + arma::accu(hatR0 == 0);
    double shape2_0 = nu0 + arma::accu(hatR0 != 0);
    hatPi0 = R::rbeta(shape1_0, shape2_0);
    gsPi0(t) = hatPi0;


    // pi.star|
    double shape1_s = muStar + arma::accu(tRsRs == 0);
    double shape2_s = nuStar + arma::accu(tRsRs != 0);
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

    if(progress != 0 && t % progress == 0){
      Rcpp::checkUserInterrupt();
      Rcpp::Rcout << "Iteration: " << t << std::endl;
      Rcpp::Rcout << "  mse    : " << arma::accu(arma::square(res))/n << std::endl;
      Rcpp::Rcout << "  sigmaSq: " << hatSigmaSq << std::endl;
      // Rcpp::Rcout << "  pi0    : " << hatPi0 << std::endl;
      // Rcpp::Rcout << "  piStar : " << hatPiStar << std::endl;
      // Rcpp::Rcout << "  piZeta : " << hatPiZeta << std::endl;
    }
  }
  return Rcpp::List::create(Rcpp::Named("GS.alpha") = gsM,
                            Rcpp::Named("GS.ata") = gsAta,
                            Rcpp::Named("GS.beta") = gsR0,
                            Rcpp::Named("GS.eta") = gsRStar,
                            Rcpp::Named("GS.tRsRs") = gsRstRs,
                            Rcpp::Named("GS.invTAUsq.1") = gsInvTauSq1,
                            Rcpp::Named("GS.invTAUsq.2") = gsInvTauSq2,
                            Rcpp::Named("GS.invTAUsq.star") = gsInvTauSqStar,
                            Rcpp::Named("GS.pi.star") = gsPiStar,
                            Rcpp::Named("GS.lambda.sq.0") = gsLambda0,
                            Rcpp::Named("GS.lambda.sq.1") = gsLambda1,
                            Rcpp::Named("GS.lambda.sq.2") = gsLambda2,
                            Rcpp::Named("GS.lambda.sq.star") = gsLambdaStar,
                            Rcpp::Named("GS.sigma.sq") = gsSigmaSq,
                            Rcpp::Named("GS.l0") = gsL0,
                            Rcpp::Named("GS.lS") = gsLS,
                            Rcpp::Named("idgene")=idgene,
                            Rcpp::Named("idgene1")=idgene1
                           );
}
