/* 
-----------------------------------------------------------------------------
    File: family_normal.cpp
    Purpose: Implementation of Normal family for different link functions
    Author: Steffen Maletz
    Last modified: 2025-12-06
-----------------------------------------------------------------------------
*/

#include "glmstarma.h"

/*
    Random number generators for Normal family
    Uses inversion method for dependence case
*/

const arma::vec Normal::random_observation_independent(const arma::vec &expectation) const
{
  if(!this->const_dispersion && dispersion_matrix.n_cols != 1)
  {
    Rcpp::stop("Dispersion is not constant. Cannot generate new observations without new dispersion values.");
  }
  arma::vec observations(expectation.n_elem, arma::fill::randn);
  if(this->const_dispersion)
  {
    observations *= std::sqrt(dispersion);
  } else {
    observations %= arma::sqrt(dispersion_matrix.col(0));
  }
  observations += expectation;
  return observations;
}

const arma::vec Normal::random_observation(const arma::vec &expectation) const
{
    if(!this->const_dispersion && dispersion_matrix.n_cols != 1)
    {
      Rcpp::stop("Dispersion is not constant. Cannot generate new observations without new dispersion values.");
    }
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    arma::vec copula_values(expectation.n_elem);
    copula_values = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));
    double disp;
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        disp = this->const_dispersion ? dispersion : dispersion_matrix(i, 0);
        observations(i) = R::qnorm( copula_values(i), expectation(i), std::sqrt(disp), true, false );
    }
    return observations;
}


const arma::vec Normal::random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const
{
  arma::vec observations(expectation.n_elem, arma::fill::randn);
  observations %= arma::sqrt(dispersion_);
  observations += expectation;
  return observations;
}

const arma::vec Normal::random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const
{
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    arma::vec copula_values(expectation.n_elem);
    copula_values = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));
    double disp;
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        disp = dispersion_(i);
        observations(i) = R::qnorm( copula_values(i), expectation(i), std::sqrt(disp), true, false );
    }
    return observations;
}


/*
  Log-Likelihood for Normal family
*/

const double Normal::log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const 
{
  return R::dnorm(observation, expectation, std::sqrt(dispersion_), true);
}

/*
  Deviance and Pearson residuals for Normal family
*/

const double Normal::deviance_residual(const double &observation, const double &expectation) const
{
    return std::pow(observation - expectation, 2);
}

const double Normal::pearson_residual(const double &observation, const double &expectation) const
{
    return std::pow(observation - expectation, 2);
}


/*
  Additional methods for Normal family
*/

const bool Normal::valid_expectation(const arma::mat &x) const 
{
  return x.is_finite();
}

const arma::mat Normal::variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const
{
    return dispersion_;
}




/*
    Define various link functions for the Normal family
*/




/*
    Define identity link, i.e. g(mu) = mu
    Observation transformation: h(y) = y
    Link transformation: h(psi) = psi
*/

class LinearNormal : public Normal {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    LinearNormal(const Rcpp::RObject &dispersion, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : Normal(dispersion, true, false, copula_obj, "identity"){};
    LinearNormal(const Rcpp::RObject &dispersion) : Normal(dispersion, true, false, "identity"){};
    ~LinearNormal() = default;
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
          return new LinearNormal(disp_vec, this->copula_object);
        } else {
          return new LinearNormal(disp_vec);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new LinearNormal(disp_mat, this->copula_object);
        } else {
          return new LinearNormal(disp_mat);
        }
      }
    };
};


const double LinearNormal::inverse_link(const double x) const
{
  return x;
}

const double LinearNormal::link(const double x) const
{
  return x;
}

const double LinearNormal::derivative_inverse_link(const double x) const
{
  return 1.0;
}


const double LinearNormal::observation_trafo(const double x) const
{
  return x;
}

const double LinearNormal::link_trafo(const double x) const
{
  return x;
}

const double LinearNormal::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool LinearNormal::valid_link(const arma::mat &x) const
{
  return true;
}


/*
    Define log-link
    Will only result in "good" results for non-negative data, or only few negative data points (close to zero)
    Link function: g(mu) = log(mu)
    Observation transformation: h(y) = log(|y|)
    Link transformation: h(psi) = psi
*/

class LogNormal : public Normal {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    LogNormal(const Rcpp::RObject &dispersion, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : Normal(dispersion, true, false, copula_obj, "log"){};
    LogNormal(const Rcpp::RObject &dispersion) : Normal(dispersion, true, false, "log"){};
    ~LogNormal() = default;
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
          return new LinearNormal(disp_vec, this->copula_object);
        } else {
          return new LinearNormal(disp_vec);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new LogNormal(disp_mat, this->copula_object);
        } else {
          return new LogNormal(disp_mat);
        }
      }
    };
};


const double LogNormal::inverse_link(const double x) const
{
  return exp(x);
}

const double LogNormal::link(const double x) const
{
  return std::log(x);
}

const double LogNormal::derivative_inverse_link(const double x) const
{
  return exp(x);
}


const double LogNormal::observation_trafo(const double x) const
{
    double y = x < 0 ? -x : x;
    y = (y == 0.0) ? 1e-8 : y;
    return log(y);
}

const double LogNormal::link_trafo(const double x) const
{
  return x;
}

const double LogNormal::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool LogNormal::valid_link(const arma::mat &x) const
{
  return true;
}

/*
  Inverse link for Normal family, i.e. g(mu) = 1/mu
    Observation transformation: h(y) = 1/y
    Link transformation: h(psi) = psi
  Only for positive data implemented
*/

class InverseNormal : public Normal {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    InverseNormal(const Rcpp::RObject &dispersion, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : Normal(dispersion, true, true, copula_obj, "inverse"){};
    InverseNormal(const Rcpp::RObject &dispersion) : Normal(dispersion, true, true, "inverse"){};
    ~InverseNormal() = default;
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
          return new InverseNormal(disp_vec, this->copula_object);
        } else {
          return new InverseNormal(disp_vec);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new InverseNormal(disp_mat, this->copula_object);
        } else {
          return new InverseNormal(disp_mat);
        }
      }
    };
};


const double InverseNormal::inverse_link(const double x) const
{
  return 1.0 / x;
}

const double InverseNormal::link(const double x) const
{
  return 1.0 / x;
}

const double InverseNormal::derivative_inverse_link(const double x) const
{
  return -1.0 / (x * x);
}


const double InverseNormal::observation_trafo(const double x) const
{
  return 1.0 / x;
}

const double InverseNormal::link_trafo(const double x) const
{
  return x;
}

const double InverseNormal::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool InverseNormal::valid_link(const arma::mat &x) const
{
  return true;
}


Normal* Normal::create(const Rcpp::List &family, Rcpp::S4& copula_obj)
{
  std::string link = Rcpp::as<std::string>(family["link"]);
  Rcpp::RObject dispersion_ = family["dispersion"];
  if(link == "log")
  {
    return new LogNormal(dispersion_, copula_obj);
  } 
  else if(link == "identity") 
  {
    return new LinearNormal(dispersion_, copula_obj);
  } 
  else if(link == "inverse") 
  {
    return new InverseNormal(dispersion_, copula_obj);
  } 
  else
  {
    throw std::invalid_argument("The desired link is currently not yet implemented.");
  }
}

Normal* Normal::create(const Rcpp::List &family)
{
  
  std::string link = Rcpp::as<std::string>(family["link"]);
  Rcpp::RObject dispersion_ = family["dispersion"];
  if(link == "log")
  {
    return new LogNormal(dispersion_);
  } 
  else if(link == "identity") 
  {
    return new LinearNormal(dispersion_);
  } 
  else if(link == "inverse") 
  {
    return new InverseNormal(dispersion_);
  } 
  else
  {
    throw std::invalid_argument("The desired link is currently not yet implemented.");
  }
}



