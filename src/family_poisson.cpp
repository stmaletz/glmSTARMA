/* 
-----------------------------------------------------------------------------
    File: family_poisson.cpp
    Purpose: Implementation of Poisson family for different link functions
    Author: Steffen Maletz
    Last modified: 2025-12-06
-----------------------------------------------------------------------------
*/

#include "glmstarma.h"

/*
    Random number generators for Poisson family
*/
const arma::vec Poisson::random_observation_independent(const arma::vec &expectation) const
{
  arma::vec observations(expectation.n_elem, arma::fill::zeros);
  for(unsigned int i = 0; i < expectation.n_elem; i++)
  {
    observations(i) = Rcpp::as<double>(Rcpp::rpois(1, expectation(i)));
  }
  return observations;
}

const arma::vec Poisson::random_observation(const arma::vec &expectation) const
{
  if(use_fast_sampling)
  {
    return random_observation_fast(expectation);
  }
  else
  {
    return random_observation_slow(expectation);
  }
}

const arma::vec Poisson::random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const
{
  return random_observation_independent(expectation);
}
const arma::vec Poisson::random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const
{
  return random_observation(expectation);
}



/*
    Poisson process based generation of (multivariate) Poisson random variables
*/

const arma::vec Poisson::random_observation_slow(const arma::vec &expectation) const
{
  arma::vec observations(expectation.n_elem, arma::fill::zeros);
  arma::uvec finished(expectation.n_elem, arma::fill::value(0));
  arma::vec copula_values(expectation.n_elem);
  arma::vec temp(expectation.n_elem, arma::fill::zeros);
  
  while (any(finished == 0))
  {
    copula_values = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));
    temp -= log(copula_values);
    finished = 1 * (temp > expectation);
    observations(arma::find(finished == 0)) += 1;
  }
  return observations;
}

/*
    Random number generation using inversion method on copula samples
*/

const arma::vec Poisson::random_observation_fast(const arma::vec &expectation) const
{
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    arma::vec copula_values(expectation.n_elem);
    copula_values = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        observations(i) = R::qpois( copula_values(i), expectation(i), true, false );
    }
    return observations;
}


/*
  Log-Likelihood for Poisson family
*/
const double Poisson::log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const 
{
  return R::dpois( observation, expectation, true );
}

/*
  Deviance and Pearson residuals for Poisson family
*/

const double Poisson::deviance_residual(const double &observation, const double &expectation) const
{
    if(observation == 0.0)
    {
        return 2.0 * expectation;
    } else {
        return 2.0 * (observation * std::log(observation / expectation) + expectation - observation);
    }
}

const double Poisson::pearson_residual(const double &observation, const double &expectation) const
{
    return std::pow(observation - expectation, 2) / expectation;
}


/*
  Additional methods for Poisson family
*/

const bool Poisson::valid_expectation(const arma::mat &x) const 
{
  return x.is_finite() && arma::all(arma::vectorise(x) > 0.0);
}

const arma::mat Poisson::variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const
{
  return this->inverse_link(link_values);
}



/*
    Define various link functions for the Poisson family
*/




/*
  Define linear model, i.e. the 'identity' link of the Poisson distribution.
  This corresponds to the linear PSTARMA model from Maletz et al. (2024)
  Link function: g(mu) = mu
  Observation transformation: h(y) = y
  Link transformation: h(psi) = psi
*/

class LinearPoisson : public Poisson {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    LinearPoisson(const bool &fast_sampling, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : Poisson(fast_sampling, true, true, copula_obj, "identity"){};
    LinearPoisson(const bool &fast_sampling) : Poisson(fast_sampling, true, true, "identity"){};
    ~LinearPoisson() = default;
    virtual const bool valid_link(const arma::mat &x) const;
    virtual Family* clone() const
    {
      if(use_dependence)
      {
        return new LinearPoisson(use_fast_sampling, this->copula_object);
      } else {
        return new LinearPoisson(use_fast_sampling);
      }
    };
};

const double LinearPoisson::inverse_link(const double x) const
{
  return x;
}

const double LinearPoisson::link(const double x) const
{
  return x;
}

const double LinearPoisson::derivative_inverse_link(const double x) const
{
  return 1.0;
}


const double LinearPoisson::observation_trafo(const double x) const
{
  return x;
}

const double LinearPoisson::link_trafo(const double x) const
{
  return x;
}

const double LinearPoisson::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool LinearPoisson::valid_link(const arma::mat &x) const
{
  return x.is_finite() && arma::all(arma::vectorise(x) > 0.0);
}


/*
  Define log-linear model, i.e. the 'log' link of the Poisson distribution.
  This results in the log-linear PSTARMA model from Maletz et al. (2024)
  Link function: g(mu) = log(mu)
  Observation transformation: h(y) = log(y + constant)
  Link transformation: h(psi) = psi
*/

class LogPoisson : public Poisson {
  protected:
    const double added_constant;
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    LogPoisson(const bool &fast_sampling, double constant, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : Poisson(fast_sampling, true, false, copula_obj, "log"), added_constant(constant){};
    LogPoisson(const bool &fast_sampling, double constant) : Poisson(fast_sampling, true, false, "log"), added_constant(constant){};
    ~LogPoisson() = default;
    virtual const bool valid_link(const arma::mat &x) const;
    virtual Family* clone() const
    {
      if(use_dependence)
      {
        return new LogPoisson(use_fast_sampling, this->added_constant, this->copula_object);
      } else {
        return new LogPoisson(use_fast_sampling, this->added_constant);
      }
    };
};

const double LogPoisson::inverse_link(const double x) const
{
  return exp(x);
}

const double LogPoisson::link(const double x) const
{
  return std::log(x);
}

const double LogPoisson::derivative_inverse_link(const double x) const
{
  return exp(x);
}


const double LogPoisson::observation_trafo(const double x) const
{
  return log(x + added_constant);
}

const double LogPoisson::link_trafo(const double x) const
{
  return x;
}

const double LogPoisson::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool LogPoisson::valid_link(const arma::mat &x) const
{
  return true;
}



/*
  Define square root model, i.e. the 'sqrt' link of the Poisson distribution.
  Link function: g(mu) = sqrt(mu)
  Observation transformation: h(y) = 2.0 * sqrt(y + 3/8) Anscombe transformation
  Link transformation: h(psi) = psi
*/

class SqrtPoisson : public Poisson {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    SqrtPoisson(const bool &fast_sampling, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : Poisson(fast_sampling, true, true, copula_obj, "sqrt"){};
    SqrtPoisson(const bool &fast_sampling) : Poisson(fast_sampling, true, true, "sqrt"){};
    ~SqrtPoisson() = default;
    virtual const bool valid_link(const arma::mat &x) const;
    virtual Family* clone() const
    {
      if(use_dependence)
      {
        return new SqrtPoisson(use_fast_sampling, this->copula_object);
      } else {
        return new SqrtPoisson(use_fast_sampling);
      }
    };
};


const double SqrtPoisson::inverse_link(const double x) const
{
  return pow(x, 2);
}

const double SqrtPoisson::link(const double x) const
{
  return std::sqrt(x);
}

const double SqrtPoisson::derivative_inverse_link(const double x) const
{
  return 2 * x;
}

const double SqrtPoisson::observation_trafo(const double x) const
{
  return std::sqrt(x + 3.0 / 8.0) * 2.0;
}

const double SqrtPoisson::link_trafo(const double x) const
{
  return x;
}

const double SqrtPoisson::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool SqrtPoisson::valid_link(const arma::mat &x) const
{
  return x.is_finite() && arma::all(arma::vectorise(x) > 0.0);
}

/*
  Define approximately linear model, i.e. the 'softplus' link of the Poisson distribution from the work of Jahn et al. (2023)
  Link function: g(mu) = log(1 + exp(mu/constant)) * constant
  Observation transformation: h(y) = y
  Link transformation: h(psi) = log(1 + exp(psi/constant)) * constant = mu
*/

class SoftPlusPoisson : public Poisson {
  protected:
    const double tuning_param;
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    SoftPlusPoisson(const bool &fast_sampling, double constant, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : Poisson(fast_sampling, false, false, copula_obj, "softplus"), tuning_param(constant){};
    SoftPlusPoisson(const bool &fast_sampling, double constant) : Poisson(fast_sampling, false, false, "softplus"), tuning_param(constant){};
    ~SoftPlusPoisson() = default;
    virtual const bool valid_link(const arma::mat &x) const;
    virtual Family* clone() const
    {
      if(use_dependence)
      {
        return new SoftPlusPoisson(use_fast_sampling, this->tuning_param, this->copula_object);
      } else {
        return new SoftPlusPoisson(use_fast_sampling, this->tuning_param);
      }
    };
};



const double SoftPlusPoisson::inverse_link(const double x) const
{
  return tuning_param * std::log(1.0 + exp(x / tuning_param));
}

const double SoftPlusPoisson::link(const double x) const
{
  return tuning_param * std::log(exp(x / tuning_param) - 1.0);
}

const double SoftPlusPoisson::derivative_inverse_link(const double x) const
{
  double val = exp(x / tuning_param);
  return val / (1.0 + val);
}


const double SoftPlusPoisson::observation_trafo(const double x) const
{
  return x;
}

const double SoftPlusPoisson::link_trafo(const double x) const
{
  return tuning_param * std::log(1.0 + exp(x / tuning_param));
}

const double SoftPlusPoisson::derivative_link_trafo(const double x) const
{
  double val = exp(x / tuning_param);
  return val / (1.0 + val);
}

const bool SoftPlusPoisson::valid_link(const arma::mat &x) const
{
  return true;
}


/*
  Factory methods for Poisson family
*/

Poisson* Poisson::create(const Rcpp::List &family, Rcpp::S4& copula_obj)
{
    std::string link = Rcpp::as<std::string>(family["link"]);
    Rcpp::LogicalVector fast_sample = family["fast_sampling"];
    bool fast_sampling = is_true ( all(fast_sample) );
    if(link == "log")
    {
      double constant = 1.0;
      if(family.containsElementNamed("const"))
      {
        Rcpp::NumericVector temp_const = family["const"];
        constant = temp_const[0];
      }
      return new LogPoisson(fast_sampling, constant, copula_obj);
    } 
    else if(link == "identity") 
    {
      return new LinearPoisson(fast_sampling, copula_obj);
    } 
    else if(link == "sqrt")
    {
      return new SqrtPoisson(fast_sampling, copula_obj);
    }
    else if(link == "softplus") 
    {
      double constant = 1.0;
      if(family.containsElementNamed("const"))
      {
        Rcpp::NumericVector temp_const = family["const"];
        constant = temp_const[0];
      }
      return new SoftPlusPoisson(fast_sampling, constant, copula_obj);
    } 
    else 
    {
      throw std::invalid_argument("The desired link is currently not yet implemented.");
    }
}
Poisson* Poisson::create(const Rcpp::List &family)
{
  
  std::string link = Rcpp::as<std::string>(family["link"]);
  if(link == "log")
  {
    double constant = 1.0;
    if(family.containsElementNamed("const"))
    {
      Rcpp::NumericVector temp_const = family["const"];
      constant = temp_const[0];
    }
    return new LogPoisson(true, constant);
  } 
  else if(link == "identity") 
  {
    return new LinearPoisson(true);
  } 
  else if(link == "sqrt") 
  {
    return new SqrtPoisson(true);
  } 
  else if(link == "softplus") 
  {
    double constant = 1.0;
    if(family.containsElementNamed("const"))
    {
      Rcpp::NumericVector temp_const = family["const"];
      constant = temp_const[0];
    }
    return new SoftPlusPoisson(true, constant);
  } 
  else
  {
    throw std::invalid_argument("The desired link is currently not yet implemented.");
  }
}
