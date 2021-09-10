// #define TMB_LIB_INIT R_init_CBFM
#include <math.h> 
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
     DATA_MATRIX(other_basis_effects_mat_B1);
     DATA_MATRIX(other_basis_effects_mat_B2);
     DATA_INTEGER(spp_ind);
     
     DATA_SPARSE_MATRIX(B); 
     DATA_MATRIX(Sigmainv_B1);
     DATA_MATRIX(Ginv_B1);
     DATA_MATRIX(Sigmainv_B2);
     DATA_MATRIX(Ginv_B2);
     
     DATA_INTEGER(family); 
     DATA_IVECTOR(trial_size);

     int num_units = y.size();
     int num_basisfns_B1 = other_basis_effects_mat_B1.cols();
     int num_basisfns_B2 = other_basis_effects_mat_B2.cols();
     int num_spp = Ginv_B1.rows();
     ////DATA_MATRIX(traits);
     ////int num_traits = traits.cols();
          
     // Declare parameters
     PARAMETER_VECTOR(basis_effects); //Spp specific coefficients for basis functions
     //// PARAMETER_MATRIX(fourthcorner_coef); // fourth corner coefficient matrix of dimension num_x by num_traits
          
     // Distribution of coefficients for basis functions corresponding to B1 -- see multivariatenormalformatrix.pdf
     matrix<Type> basis_effects_mat_B1(num_spp, num_basisfns_B1);
     for(int k=0; k<num_basisfns_B1; k++) {
          basis_effects_mat_B1(spp_ind - 1, k) = basis_effects(k);
          }
     for(int j=0; j<num_spp; j++) {
          if(j != (spp_ind-1))
               basis_effects_mat_B1.row(j) = other_basis_effects_mat_B1.row(j);
          }
     matrix<Type> M_B1 = Ginv_B1 * basis_effects_mat_B1 * Sigmainv_B1 * basis_effects_mat_B1.transpose();
     Type nll = Type(0.5) * M_B1.trace();
     
     
     // Distribution of coefficients for basis functions corresponding to B2 -- see multivariatenormalformatrix.pdf
     matrix<Type> basis_effects_mat_B2(num_spp, num_basisfns_B2);
     for(int k=0; k<num_basisfns_B2; k++) {
          basis_effects_mat_B2(spp_ind - 1, k) = basis_effects(num_basisfns_B1 + k);
          }
     for(int j=0; j<num_spp; j++) {
          if(j != (spp_ind-1))
               basis_effects_mat_B2.row(j) = other_basis_effects_mat_B2.row(j);
          }

     matrix<Type> M_B2 = Ginv_B2 * basis_effects_mat_B2 * Sigmainv_B2 * basis_effects_mat_B2.transpose();
     nll += Type(0.5) * M_B2.trace();
     
     
     // Data likelihood     
     vector<Type> lik_val(num_units);
     vector<Type> eta = Xbeta + B * basis_effects + offset;
     // eta(i,j) += x(i,k1)*traits(j,k1)*fourthcorner_coef(k1,k1); 
     if(family == 1) { //beta
          for(int i=0; i<num_units; i++) { 
               nll -= dbeta(y(i), dispparam(0)*invlogit(eta(i)), dispparam(0)*(1-invlogit(eta(i))), true);
               }
          }
     if(family == 2) { //binomial (logit link)
          Type predvalue = Type(0.0);
          
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
     
