/* 
-----------------------------------------------------------------------------
    File: family_gamma.cpp
    Purpose: Implementation of Gamma family for different link functions
    Author: Steffen Maletz
    Last modified: 2025-12-06
-----------------------------------------------------------------------------
*/

#include "glmstarma.h"


/*
  Random number generators for Gamma family
  These either use functions from R's stats package or use inversion method on copula samples
*/

const arma::vec Gamma::random_observation_independent(const arma::vec &expectation) const
{
  if(!this->const_dispersion && dispersion_matrix.n_cols != 1)
  {
    Rcpp::stop("Dispersion is not constant. Cannot generate new observations without new dispersion values.");
  }
  arma::vec observations(expectation.n_elem, arma::fill::zeros);
  double shap;
  for(unsigned int i = 0; i < expectation.n_elem; i++)
  {
    shap = this->const_dispersion ? shape : shape_matrix(i, 0);
    observations(i) = Rcpp::as<double>(Rcpp::rgamma(1, shap, expectation(i) * (1.0 / shap)));
  }
  return observations;
}

const arma::vec Gamma::random_observation(const arma::vec &expectation) const
{
    if(!this->const_dispersion && dispersion_matrix.n_cols != 1)
    {
      Rcpp::stop("Dispersion is not constant. Cannot generate new observations without new dispersion values.");
    }
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    arma::vec copula_values(expectation.n_elem);
    copula_values = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));
    double shap;
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        shap = this->const_dispersion ? shape : shape_matrix(i, 0);
        observations(i) = R::qgamma( copula_values(i), shap, expectation(i) * (1.0 / shap), true, false );
    }
    return observations;
}

const arma::vec Gamma::random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const
{
  arma::vec observations(expectation.n_elem, arma::fill::zeros);
  double shap;
  for(unsigned int i = 0; i < expectation.n_elem; i++)
  {
    shap = 1.0 / dispersion_(i);
    observations(i) = Rcpp::as<double>(Rcpp::rgamma(1, shap, expectation(i) * dispersion_(i)));
  }
  return observations;
}

const arma::vec Gamma::random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const
{
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    arma::vec copula_values(expectation.n_elem);
    copula_values = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));
    double shap;
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        shap = 1.0 / dispersion_(i);
        observations(i) = R::qgamma( copula_values(i), shap, expectation(i) * dispersion_(i), true, false );
    }
    return observations;
}


/*
  Log-likelihood for Gamma family
*/
const double Gamma::log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const 
{
  double shap = 1.0 / dispersion_;
  return R::dgamma(observation, shap, expectation * dispersion, true);
}


/*
    Deviance and Pearson residuals for Gamma family
*/
const double Gamma::deviance_residual(const double &observation, const double &expectation) const
{
    return 2.0 * ((observation - expectation) / expectation - std::log(observation / expectation ) );  
}

const double Gamma::pearson_residual(const double &observation, const double &expectation) const
{
    return std::pow(observation - expectation, 2) / std::pow(expectation, 2);
}


/*
  Additional methods for Gamma family
*/

const bool Gamma::valid_expectation(const arma::mat &x) const 
{
  return x.is_finite() && arma::all(arma::vectorise(x) > 0.0);
}

const arma::mat Gamma::variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const
{
    arma::mat response = this->inverse_link(link_values);
    return (response % response) % dispersion_;
}

/*
    Define log-model, i.e., g(mu) = log(mu)
    Observation transformation: h(y) = log(y)
    Link transformation: h(psi) = psi
*/

class LogGamma : public Gamma {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    LogGamma(const Rcpp::RObject &dispersion, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : Gamma(dispersion, true, false, copula_obj, "log"){};
    LogGamma(const Rcpp::RObject &dispersion) : Gamma(dispersion, true, false, "log"){};
    ~LogGamma() = default;
    virtual const bool valid_link(const arma::mat &x) const;
    virtual Family* clone() const
    {
      Rcpp::NumericMatrix disp_mat;
      Rcpp::NumericVector disp_vec;
      if(this->const_dispersion)
      {
        disp_vec = Rcpp::NumericVector::create(this->dispersion);
        if(use_dependence)
        { 
          return new LogGamma(disp_vec, this->copula_object);
        } else {
          return new LogGamma(disp_vec);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new LogGamma(disp_mat, this->copula_object);
        } else {
          return new LogGamma(disp_mat);
        }
      }
    };
};


const double LogGamma::inverse_link(const double x) const
{
  return exp(x);
}

const double LogGamma::link(const double x) const
{
  return std::log(x);
}

const double LogGamma::derivative_inverse_link(const double x) const
{
  return exp(x);
}


const double LogGamma::observation_trafo(const double x) const
{
  return log(x);
}

const double LogGamma::link_trafo(const double x) const
{
  return x;
}

const double LogGamma::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool LogGamma::valid_link(const arma::mat &x) const
{
  return true;
}

/*
  Define various link functions for Gamma family
*/



/*
    Define linear, i.e. identity-model, i.e., g(mu) = mu
    Observation transformation: h(y) = y
    Link transformation: h(psi) = psi
*/

class LinearGamma : public Gamma {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    LinearGamma(const Rcpp::RObject &dispersion, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : Gamma(dispersion, true, true, copula_obj, "identity"){};
    LinearGamma(const Rcpp::RObject &dispersion) : Gamma(dispersion, true, true, "identity"){};
    ~LinearGamma() = default;
    virtual const bool valid_link(const arma::mat &x) const;
    virtual Family* clone() const
    {
      Rcpp::NumericMatrix disp_mat;
      Rcpp::NumericVector disp_vec;
      if(this->const_dispersion)
      {
        disp_vec = Rcpp::NumericVector::create(this->dispersion);
        if(use_dependence)
        { 
          return new LinearGamma(disp_vec, this->copula_object);
        } else {
          return new LinearGamma(disp_vec);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new LinearGamma(disp_mat, this->copula_object);
        } else {
          return new LinearGamma(disp_mat);
        }
      }
    };
};


const double LinearGamma::inverse_link(const double x) const
{
  return x;
}

const double LinearGamma::link(const double x) const
{
  return x;
}

const double LinearGamma::derivative_inverse_link(const double x) const
{
  return 1.0;
}


const double LinearGamma::observation_trafo(const double x) const
{
  return x;
}

const double LinearGamma::link_trafo(const double x) const
{
  return x;
}

const double LinearGamma::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool LinearGamma::valid_link(const arma::mat &x) const
{
  return x.is_finite() && arma::all(arma::vectorise(x) > 0.0);
}


/*
    Define inverse model, i.e., g(mu) = 1/mu
    Observation transformation: h(y) = 1/y
    Link transformation: h(psi) = psi
*/

class InverseGamma : public Gamma {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    InverseGamma(const Rcpp::RObject &dispersion, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : Gamma(dispersion, true, true, copula_obj, "inverse"){};
    InverseGamma(const Rcpp::RObject &dispersion) : Gamma(dispersion, true, true, "inverse"){};
    ~InverseGamma() = default;
    virtual const bool valid_link(const arma::mat &x) const;
    virtual Family* clone() const
    {
      Rcpp::NumericMatrix disp_mat;
      Rcpp::NumericVector disp_vec;
      if(this->const_dispersion)
      {
        disp_vec = Rcpp::NumericVector::create(this->dispersion);
        if(use_dependence)
        { 
          return new InverseGamma(disp_vec, this->copula_object);
        } else {
          return new InverseGamma(disp_vec);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new InverseGamma(disp_mat, this->copula_object);
        } else {
          return new InverseGamma(disp_mat);
        }
      }
    };
};


const double InverseGamma::inverse_link(const double x) const
{
  return 1.0 / x;
}

const double InverseGamma::link(const double x) const
{
  return 1.0 / x;
}

const double InverseGamma::derivative_inverse_link(const double x) const
{
  return -1.0 / (x * x);
}


const double InverseGamma::observation_trafo(const double x) const
{
  return 1.0 / x;
}

const double InverseGamma::link_trafo(const double x) const
{
  return x;
}

const double InverseGamma::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool InverseGamma::valid_link(const arma::mat &x) const
{
  return x.is_finite() && arma::all(arma::vectorise(x) > 0.0);
}



/*
    Factory methods for Gamma family
*/


Gamma* Gamma::create(const Rcpp::List &family, Rcpp::S4& copula_obj)
{
  std::string link = Rcpp::as<std::string>(family["link"]);
  Rcpp::RObject dispersion = family["dispersion"];
  if(link == "log")
  {
    return new LogGamma(dispersion, copula_obj);
  } 
  else if(link == "identity") 
  {
    return new LinearGamma(dispersion, copula_obj);
  } 
  else if(link == "inverse") 
  {
    return new InverseGamma(dispersion, copula_obj);
  } 
  else
  {
    throw std::invalid_argument("The desired link is currently not yet implemented.");
  }
}

Gamma* Gamma::create(const Rcpp::List &family)
{
  std::string link = Rcpp::as<std::string>(family["link"]);
  Rcpp::RObject dispersion = family["dispersion"];
  if(link == "log")
  {
    return new LogGamma(dispersion);
  } 
  else if(link == "identity") 
  {
    return new LinearGamma(dispersion);
  } 
  else if(link == "inverse") 
  {
    return new InverseGamma(dispersion);
  } 
  else
  {
    throw std::invalid_argument("The desired link is currently not yet implemented.");
  }
}










