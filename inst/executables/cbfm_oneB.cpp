//-----------------------------
// Community-level basis function model, using spatio-temporal basis functions
// Author: Francis KC Hui
//-----------------------------
// #include <math.h> 
#include <TMB.hpp> 

template<class Type>

Type objective_function<Type>::operator() () { 
     using namespace density;
     using namespace Eigen;
     
     // Declare data
     DATA_VECTOR(y); 
     DATA_VECTOR(Xbeta);
     DATA_VECTOR(dispparam);
     DATA_VECTOR(powerparam);
     DATA_VECTOR(offset); 
     DATA_VECTOR(estep_weights); 
     DATA_MATRIX(other_basis_effects_mat);
     DATA_INTEGER(spp_ind);
     
     DATA_SPARSE_MATRIX(B); 
     DATA_MATRIX(Sigmainv);
     DATA_MATRIX(Ginv);
     
     DATA_INTEGER(family); 
     DATA_IVECTOR(trial_size);
     
     int num_units = y.size();
     int num_basisfns = B.cols();
     int num_spp = Ginv.rows();
     ////DATA_MATRIX(traits);
     ////int num_traits = traits.cols();
          
     // Declare parameters
     PARAMETER_VECTOR(basis_effects); //Spp specific coefficients for basis functions
     //// PARAMETER_MATRIX(fourthcorner_coef); // fourth corner coefficient matrix of dimension num_x by num_traits
          
     matrix<Type> basis_effects_mat(num_spp, num_basisfns);
     basis_effects_mat.row(spp_ind - 1) = basis_effects;
     for(int j=0; j<num_spp; j++) {
          if(j != (spp_ind-1))
               basis_effects_mat.row(j) = other_basis_effects_mat.row(j);
          }
     
     // Distribution of random effects -- see multivariatenormalformatrix.pdf
     matrix<Type> M = Ginv * basis_effects_mat * Sigmainv * basis_effects_mat.transpose();
     Type nll = Type(0.5) * M.trace();
     
     
     // Data likelihood     
     vector<Type> eta = Xbeta + B * basis_effects + offset;
     // eta(i,j) += x(i,k1)*traits(j,k1)*fourthcorner_coef(k1,k1); 
     Type predvalue = Type(0.0);

     if(family == 1) { //beta
          for(int i=0; i<num_units; i++) { 
               nll -= dbeta(y(i), dispparam(0)*invlogit(eta(i)), dispparam(0)*(1-invlogit(eta(i))), true);
               }
          }
     if(family == 2) { //binomial (logit link)          
          for(int i=0; i<num_units; i++) { 
               predvalue = dbinom(y(i), Type(trial_size(i)), invlogit(eta(i)), true);
               if(predvalue < -10000) 
                    predvalue = -10000;
               nll -= predvalue; 
               }
          }
     if(family == 3) { //gamma(logit link)
          for(int i=0; i<num_units; i++) { 
               nll -= dgamma(y(i), 1/dispparam(0), exp(eta(i))*dispparam(0), true);
               }
          }
     if(family == 4) { //negative binomial
          Type predvalue = Type(0.0);

          for(int i=0; i<num_units; i++) { 
            //nll -= dnbinom2(y(i), exp(eta(i)), exp(eta(i)) + dispparam(0)*pow(exp(eta(i)),2), true);
            predvalue = 1/(1+dispparam(0)*exp(eta(i)));
            if(predvalue > 0.9999)
                predvalue = 0.9999;
            nll -= dnbinom(y(i), 1/dispparam(0), predvalue, true);
            //lik_val(i) = dnbinom(y(i), 1/dispparam(0), predvalue, true);
            }
          }
     if(family == 5) { //normal
          for(int i=0; i<num_units; i++) { 
               nll -= dnorm(y(i), eta(i), dispparam(0), true);
               }
          }
     if(family == 6) { //poisson
          for(int i=0; i<num_units; i++) { 
               nll -= dpois(y(i), exp(eta(i)), true);
               }
          }
     if(family == 7) { //tweedie
          for(int i=0; i<num_units; i++) { 
               nll -= dtweedie(y(i), exp(eta(i)), dispparam(0), powerparam(0), true);
               }
          }
     if(family == 8) { //zero-truncated negative binomial
          for(int i=0; i<num_units; i++) { 
               nll -= (dnbinom2(y(i), exp(eta(i)), exp(eta(i)) + dispparam(0)*pow(exp(eta(i)),2), true) - log(1-dnbinom2(Type(0.0), exp(eta(i)), exp(eta(i)) + dispparam(0)*pow(exp(eta(i)),2), false)));
               }
          }
     if(family == 9) { //zero-truncated poisson
          for(int i=0; i<num_units; i++) { 
               nll -= (dpois(y(i), exp(eta(i)), true) - log(1 - exp(-eta(i))));
               }
          }
     if(family == 10) { //zero-inflated poisson
          for(int i=0; i<num_units; i++) { 
               nll -= (Type(1)-estep_weights(i))*dpois(y(i), exp(eta(i)), true);
               }
          }
     if(family == 11) { //zero-inflated negative binomial
          for(int i=0; i<num_units; i++) { 
               nll -= (Type(1)-estep_weights(i))*dnbinom2(y(i), exp(eta(i)), exp(eta(i)) + dispparam(0)*pow(exp(eta(i)),2), true);
               }
          }
          
          
     ADREPORT(basis_effects);
     return nll;     
     }
     
