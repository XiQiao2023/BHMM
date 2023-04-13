#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;

// [[Rcpp::export]]
//Define my multivariate normal
arma::vec mymv(arma::mat mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat v = arma::randn(ncols);
  return mu + arma::chol(sigma) * v;
}

// [[Rcpp::export]]
List Gibbs_Posterior(arma::mat Group, arma::mat M, arma::mat T, 
                     arma::mat X, arma::mat Y,
                     int iteration, int burn, int chain, 
                     int p, int n2, int n, int c) {
  
  List posterior(chain);
  
  for (int l = 0; l < chain; ++l) {
    int iter = iteration + burn;
    
    // Mediator Model Parameters
    arma::mat PIP_alphaT(iter, p);
    arma::mat pos_alphaT(iter, p);
    
    // Mediator Model Initial Value
    arma::mat alpha0(n2, p);
    arma::mat alphaX(p, c);
    arma::vec Alpha0 (p);
    
    arma::mat Salpha0(p, p);
    for (int i = 0; i < p; i++) {
      Salpha0(i, i) = 1.0;
    }
    
    arma::mat Mintercept(n, p);
    arma::vec PalphaT(p);
    arma::vec PIPalphaT(p);
    arma::vec IalphaT(p);
    arma::vec alphaT(p);
    double SalphaT1 = 1;
    double SalphaT0 = 1;
    
    arma::mat sigmaM(p, p);
    
    for (int i = 0; i < p; i++) {
      sigmaM(i, i) = 1.0;
    }
    
    for (int i = 0; i < iter; ++i) {
      
      arma::mat inv_Salpha0 = arma::inv(Salpha0);
      arma::mat inv_sigmaM = arma::inv(sigmaM);
      
      //sampling alpha0
      for (int j = 0; j < n2; ++j) {
        arma::vec alpha0j = arma::conv_to<arma::vec>::from(Group == j+1);
        arma::mat va0 = arma::inv(inv_Salpha0 + arma::accu(arma::pow(alpha0j,2)) * inv_sigmaM);
        arma::mat sa0 = (M - T * (alphaT.t()) - X * (alphaX.t())).t() * alpha0j;
        arma::mat ea0 = va0 * (inv_Salpha0 * Alpha0 + inv_sigmaM * sa0);
        alpha0.row(j) = arma::trans(mymv(ea0, va0));
      }
      
      for (int s = 0; s < n; ++s) {
        int group_index = Group(s, 0);
        Mintercept.row(s) = alpha0.row(group_index-1);
      }
      
      //sampling Alpha0
      arma::mat vA0 = Salpha0 / n2;
      arma::mat eA0 = vA0  * inv_Salpha0 * (alpha0.t() * arma::ones(n2));
      Alpha0 = mymv(eA0, vA0);
      
      //sampling Salpha0
      arma::mat ssalpha0 = (alpha0 - arma::ones(n2) * (arma::trans(Alpha0))).t() * (alpha0 -  arma::ones(n2) * (arma::trans(Alpha0)));
      arma::mat diagp;
      Salpha0 = arma::iwishrnd(diagp.eye(p,p) + ssalpha0, p + 2 + n2);
      
      //sampling Indicator of alphaT
      for (int m = 0; m < p; ++m) {
        double a = R::dnorm(alphaT(m), 0, std::sqrt(SalphaT1),false) * PalphaT(m);
        double b = R::dnorm(alphaT(m), 0, std::sqrt(SalphaT0),false) * (1-PalphaT(m));
        PIPalphaT(m) =  a / (a+b);
        IalphaT(m) = R::rbinom(1, a / (a+b));
        PalphaT(m) = R::rbeta(1 +  IalphaT(m), 2 -  IalphaT(m));
      }
      
      //sampling alphaT
      arma::mat vaT0 = arma::diagmat(IalphaT) * SalphaT1 + arma::diagmat(1 - IalphaT) * SalphaT0;
      arma::mat vaT = arma::inv(arma::inv(vaT0) + arma::accu(arma::pow(T,2)) * inv_sigmaM);
      arma::mat saT = (M - Mintercept - X * (alphaX.t())).t() * T;
      arma::mat eaT = vaT * inv_sigmaM * saT; 
      alphaT = mymv(eaT, vaT);
      
      //sampling alphaX
      for (int e = 0; e < c; ++e) {
        arma::mat vaX = sigmaM / arma::accu(arma::pow(X.col(e), 2));
        arma::mat saX = (M - Mintercept - T * (alphaT.t()) - X * (alphaX.t()) + X.col(e) * (alphaX.col(e).t())) .t() * X.col(e);
        arma::mat eaX = vaX * inv_sigmaM * saX;
        alphaX.col(e) = mymv(eaX, vaX);
      }
      
      //sampling SalphaT1 and SalphaT0
      SalphaT1 = 1/(R::rgamma(arma::accu(IalphaT) / 2 + 2, 1/(arma::accu(arma::pow(IalphaT % alphaT,  2)) / 2 + 1)));
      SalphaT0 = 1/(R::rgamma((p - arma::accu(IalphaT)) / 2 + 2, 1/(arma::accu(arma::pow((1-IalphaT) % alphaT,  2)) / 2 + 10E-4)));
      
      //sampling sigmaM
      arma::mat MSSR = (M - Mintercept - T * (alphaT.t())- X * (alphaX.t())).t() * (M - Mintercept - T * (alphaT.t())- X * (alphaX.t()));
      sigmaM = arma::iwishrnd(diagp.eye(p,p) + MSSR, p + 2 + n);
      
      //Update Value
      PIP_alphaT.row(i) = arma::trans(PIPalphaT);
      pos_alphaT.row(i) = arma::trans(alphaT);
    }
    
    // Delete burning and get PIP
    pos_alphaT = pos_alphaT.rows(burn, pos_alphaT.n_rows - 1);
    arma::vec PIP_a = arma::mean(PIP_alphaT.rows(burn, PIP_alphaT.n_rows - 1),0).t();

    //Outcome Model Parameters 
    arma::mat PIP_betaM(iter, p);
    arma::mat pos_betaM (iter, p);
    
    //Outcome Model Initial Value
    arma::vec beta0(n2);
    double Beta0 = 0;
    double Sbeta0 = 1;
    arma::vec Yintercept(n);
      
    double betaT = 0;
    arma::vec PbetaM(p);
    arma::vec PIPbetaM(p);
    arma::vec IbetaM(p);
    arma::vec betaM(p);
    double SbetaM0 = 1;
    double SbetaM1 = 1;
    
    arma::vec betaX(c);
    double sigmaY = 1;
    
    for (int i = 0; i < iter; ++i) {
      //sampling beta0
      for (int j = 0; j < n2; ++j) {
        arma::vec beta0j = arma::conv_to<arma::vec>::from(Group == j+1);
        double vb0 = 1/((1 / Sbeta0 + arma::accu(pow(beta0j, 2)) / sigmaY));
        double sb0 = arma::accu(beta0j % (Y - T * betaT - M * betaM - X * betaX));
        double eb0 = vb0 * (Beta0  / Sbeta0 + sb0 / sigmaY);
        beta0(j) = R::rnorm(eb0, std::sqrt(vb0));
      }
      
      for (int s = 0; s < n; ++s) {
        int group_index = Group(s, 0);
        Yintercept(s) = beta0(group_index-1);
        }
      
      //sampling Beta0
      double vB0 = Sbeta0/n2;
      double eB0 = vB0 * n2 * arma::mean(beta0) / Sbeta0;
      Beta0 = R::rnorm(eB0, std::sqrt(vB0));
      
      //sampling Sbeta0
      double ssbeta0 = arma::accu(pow(beta0 - Beta0,  2));
      Sbeta0 = 1/(R::rgamma(n2/2, 2/ssbeta0));
      
      //sampling betaT
      double vbT = 1/(arma::accu(pow(T,2)) / sigmaY);
      double sbT = arma::accu(T % (Y - Yintercept - M*betaM - X*betaX));
      double ebT = vbT * sbT / sigmaY;
      betaT = R::rnorm(ebT, std::sqrt(vbT));
      
      for (int m = 0; m < p; ++m) {
        //sampling Indicator of betaM
        double a = R::dnorm(betaM[m], 0, std::sqrt(SbetaM1),false) * PbetaM(m);
        double b = R::dnorm(betaM[m], 0, std::sqrt(SbetaM0),false) * (1-PbetaM(m));
        PIPbetaM(m) = a / (a + b);
        IbetaM(m) = R::rbinom(1, a / (a + b));
        PbetaM(m) = R::rbeta(1 +  IbetaM(m), 2 -IbetaM(m));
        
        //sampling betaM
        double vbM0 = IbetaM(m) * SbetaM1 + (1 - IbetaM(m)) * SbetaM0;
        double vbM = 1 / (1 / vbM0 + arma::accu(pow(M.col(m), 2))/ sigmaY);
        double sbM = arma::accu(M.col(m) % (Y - Yintercept - T*betaT - X*betaX - M*betaM + M.col(m)*betaM(m)));
        double ebM = vbM * sbM / sigmaY;
        betaM(m) = R::rnorm(ebM, std::sqrt(vbM));
      }
      
      //sampling SbetaM1 and SbetaM0
      SbetaM1 = 1/(R::rgamma(arma::accu(IbetaM) / 2 + 2, 1/(arma::accu(pow(IbetaM % betaM, 2)) / 2  + 1)));
      SbetaM0 = 1/(R::rgamma((p - arma::accu(IbetaM)) / 2 + 2, 1/(arma::accu(pow((1 - IbetaM) % betaM,  2)) / 2  + 10E-4)));
      
      //sampling betaX
      for (int e = 0; e < c; ++e) {
        double vbX = sigmaY / arma::accu(pow(X.col(e), 2));
        double sbX =  arma::accu(X.col(e) % (Y - Yintercept - T*betaT - M*betaM - X*betaX + X.col(e)*betaX(e)));
        double ebX = vbX * sbX / sigmaY;
        betaX(e) = R::rnorm(ebX, std::sqrt(vbX));
        }
      
      //sampling sigmaY
      double YSSR = arma::accu(pow((Y - Yintercept - T * betaT - M * betaM - X * betaX), 2));
      sigmaY = 1/ (R::rgamma(n/2 , 2/YSSR));
      
      //Outcome Model Updated Value
      PIP_betaM.row(i) = PIPbetaM.t();
      pos_betaM.row(i) = betaM.t();
      }
    
    // Delete burning and get PIP
    pos_betaM = pos_betaM.rows(burn, pos_betaM.n_rows - 1);
    arma::vec PIP_b = arma::mean(PIP_betaM.rows(burn, PIP_betaM.n_rows - 1),0).t();
    
    //return necessary posteriors
    posterior[l] = List::create(Named("PIPbetaM") = PIP_b, Named("betaM") = pos_betaM,
                                Named("PIPalphaT") = PIP_a, Named("alphaT") = pos_alphaT);
    
  }
  
  return posterior;
}



