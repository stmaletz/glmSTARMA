/* 
-----------------------------------------------------------------------------
    File: family.cpp
    Purpose: Implementation of Family factory methods
    Author: Steffen Maletz
    Last modified: 2026-08-11
-----------------------------------------------------------------------------
*/

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include "glmstarma.h"
using namespace Rcpp;

/*
  Factory method for R Input including copula object
*/

Family * Family::create(const Rcpp::List &family, Rcpp::S4& copula_obj)
{
  std::string distribution = Rcpp::as<std::string>(family["distribution"]);
  if(distribution == "Poisson")
  {
    return Poisson::create(family, copula_obj);
  }
  else if(distribution == "Negative Binomial")
  {
    return NegativeBinomial::create(family, copula_obj);
  }
  else if(distribution == "Binomial")
  {
    return Binomial::create(family, copula_obj);
  }
  else if(distribution == "Quasibinomial")
  {
    return QuasiBinomial::create(family, copula_obj);
  }
  else if(distribution == "Quasipoisson")
  {
    return QuasiPoisson::create(family, copula_obj);
  }
  else if(distribution == "Gamma")
  {
    return Gamma::create(family, copula_obj);
  }
  else if(distribution == "Inverse Gaussian")
  {
    return InverseGauss::create(family, copula_obj);
  }
  else if(distribution == "Gaussian")
  {
    return Normal::create(family, copula_obj);
  }
  else 
  {
    throw std::invalid_argument("The desired distribution is currently not yet implemented.");
  }
}

/*
  Factory method for R Input without copula object
*/

Family * Family::create(const Rcpp::List &family){
  std::string distribution = Rcpp::as<std::string>(family["distribution"]);
  if(distribution == "Poisson")
  {
    return Poisson::create(family);
  }
  else if(distribution == "Negative Binomial")
  {
    return NegativeBinomial::create(family);
  }
  else if(distribution == "Binomial")
  {
    return Binomial::create(family);
  }
  else if(distribution == "Quasibinomial")
  {
    return QuasiBinomial::create(family);
  }
  else if(distribution == "Quasipoisson")
  {
    return QuasiPoisson::create(family);
  }
  else if(distribution == "Gamma")
  {
    return Gamma::create(family);
  }
  else if(distribution == "Inverse Gaussian")
  {
    return InverseGauss::create(family);
  }
  else if(distribution == "Gaussian")
  {
    return Normal::create(family);
  }
  else 
  {
    throw std::invalid_argument("The desired distribution is currently not yet implemented.");
  }
}

