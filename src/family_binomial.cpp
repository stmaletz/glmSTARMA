/* 
-----------------------------------------------------------------------------
    File: family_binomial.cpp
    Purpose: Implementation of Binomial family for different link functions
    Author: Steffen Maletz
    Last modified: 2025-12-06
-----------------------------------------------------------------------------
*/

#include "glmstarma.h"



/*
  Random number generators for Binomial family
  These either use functions from R's stats package or use inversion method on copula samples
*/

const arma::vec Binomial::random_observation_independent(const arma::vec &expectation) const
{
  arma::vec observations(expectation.n_elem, arma::fill::zeros);
  unsigned int n;
  for(unsigned int i = 0; i < expectation.n_elem; i++)
  {
    n = n_upper.n_elem > 1 ? n_upper(i) : n_upper(0);
    observations(i) = Rcpp::as<double>(Rcpp::rbinom(1, n, expectation(i) / n));
  }
  return observations;
}

const arma::vec Binomial::random_observation(const arma::vec &expectation) const
{
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    arma::vec copula_values = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));
    unsigned int n;
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        n = n_upper.n_elem > 1 ? n_upper(i) : n_upper(0);
        observations(i) = R::qbinom(copula_values(i), n, expectation(i) / n, true, false);
    }
    return observations;
}

const arma::vec Binomial::random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const
{
  arma::vec observations(expectation.n_elem, arma::fill::zeros);
  unsigned int n;
  for(unsigned int i = 0; i < expectation.n_elem; i++)
  {
    n = n_upper.n_elem > 1 ? n_upper(i) : n_upper(0);
    observations(i) = Rcpp::as<double>(Rcpp::rbinom(1, n, expectation(i) / n));
  }
  return observations;
}

const arma::vec Binomial::random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const
{
  arma::vec observations(expectation.n_elem, arma::fill::zeros);
  arma::vec copula_values = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));
  unsigned int n;
  for(unsigned int i = 0; i < expectation.n_elem; i++)
  {
      n = n_upper.n_elem > 1 ? n_upper(i) : n_upper(0);
      observations(i) = R::qbinom(copula_values(i), n, expectation(i) / n, true, false);
  }
  return observations;
}


/*
  Log-likelihood for Binomial family
*/

const double Binomial::log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const 
{
  double value = R::dbinom( observation, n_upper(index_n_upper), expectation / n_upper(index_n_upper), true );
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
  return value;
}


/*
  Deviance and Pearson residuals for Binomial family
*/
const double Binomial::deviance_residual(const double &observation, const double &expectation) const
{
    int n = n_upper(index_n_upper);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
    if(observation <= 0.0)
    {
        return 2.0 * n * std::log((n - observation) / (n - expectation));
    } else if(observation >= n)
    {
        return 2.0 * n * std::log(n / expectation);
    } else {
        return 2.0 * (observation * std::log(observation / expectation) + (n - observation) * std::log((n - observation) / (n - expectation)));
    }
}

const double Binomial::pearson_residual(const double &observation, const double &expectation) const
{
    int n = n_upper(index_n_upper);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem;
    double p = expectation / n;
    return std::pow(observation - expectation, 2) / (n * p * (1.0 - p));
}

/*
  Additional methods for Binomial family
*/

const bool Binomial::valid_expectation(const arma::mat &x) const 
{
  bool is_valid = x.is_finite() && arma::all(arma::vectorise(x) >= 0.0);
  if(is_valid)
  {
    if(n_upper.n_elem > 1)
    {
      return arma::all(arma::all(x <= arma::repmat(n_upper, 1, x.n_cols)));
    } else {
      return arma::all(arma::vectorise(x) <= n_upper(0));
    }
  }
  return false;
}

const arma::mat Binomial::variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const
{
    arma::mat fitted = this->inverse_link(link_values);
    if(n_upper.n_elem > 1)
    {
        fitted.each_col() /=  arma::conv_to<arma::vec>::from(n_upper);
    } else {
        fitted = fitted / n_upper(0);
    }
    arma::mat ones = arma::ones(fitted.n_rows, fitted.n_cols);
    fitted %= (ones - fitted);
    if(n_upper.n_elem > 1)
    {
        fitted.each_col() %= arma::conv_to<arma::vec>::from(n_upper);
    } else {
        fitted *= n_upper(0);
    }
    return fitted;
}


/*
  Define various link functions for Binomial family
*/

/*
  Define approximately linear model, i.e. the 'softclipping' link of the Binomial distribution from the work of Jahn et al. (2023)
  link function: g(mu) = log( (1 + exp(mu/constant)) / (1 + exp((mu - n)/constant)) ) * constant
  Observation transformation: h(y) = y/n
  Link transformation: h(psi) = (1 + exp(psi/constant)) / (1 + exp((psi - n)/constant)) * constant = mu
*/

class SoftClippingBinomial : public Binomial {
  protected:
    const double tuning_param;
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    SoftClippingBinomial(const arma::uvec &n, double constant, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : Binomial(n, false, false, copula_obj, "softclipping"), tuning_param(constant){};
    SoftClippingBinomial(const arma::uvec &n, double constant) : Binomial(n, false, false, "softclipping"), tuning_param(constant){};
    virtual ~SoftClippingBinomial() = default;
    virtual const bool valid_link(const arma::mat &x) const;
    virtual Family* clone() const
    {
      if(use_dependence)
      {
        return new SoftClippingBinomial(this->n_upper, this->tuning_param, this->copula_object);
      } else {
        return new SoftClippingBinomial(this->n_upper, this->tuning_param);
      }
    };
};


const double SoftClippingBinomial::inverse_link(const double x) const
{
    double upper = 1.0 + exp(x / tuning_param);
    double lower = 1.0 + exp((x - 1.0) / tuning_param);
    double value = n_upper(index_n_upper) * tuning_param * std::log(upper / lower);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
    return value;
}

const double SoftClippingBinomial::link(const double x) const
{
    unsigned int n = n_upper(index_n_upper);
    double upper = std::exp(x / (n * tuning_param)) - 1.0;
    double lower = 1.0 - std::exp((x - n) / (n * tuning_param));
    double value = tuning_param * std::log(upper / lower);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
    return value;
}

const double SoftClippingBinomial::derivative_inverse_link(const double x) const
{
  double val_first = exp(x / tuning_param) / (1.0 + exp(x / tuning_param));
  double val_second = exp((x - 1.0) / tuning_param) / (1.0 + exp((x - 1.0) / tuning_param));
  double value = n_upper(index_n_upper) * (val_first - val_second);
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem;
  return value;
}

const double SoftClippingBinomial::observation_trafo(const double x) const
{
  double value = x / n_upper(index_n_upper);
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem;
  return value;
}


const double SoftClippingBinomial::link_trafo(const double x) const
{
    double value = inverse_link(x) / n_upper(index_n_upper);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem;
    return value;
}

const double SoftClippingBinomial::derivative_link_trafo(const double x) const
{
  double value = derivative_inverse_link(x) / n_upper(index_n_upper);
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem;
  return value;
}

const bool SoftClippingBinomial::valid_link(const arma::mat &x) const
{
  return true;
}

/*
  Implementation of the Identity link for the Binomial family, i.e., g(mu) = mu
  Observation transformation: h(y) = y/n
  Link transformation: h(psi) = psi
*/


class IdentityBinomial : public Binomial {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    IdentityBinomial(const arma::uvec &n, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : Binomial(n, false, true, copula_obj, "identity"){};
    IdentityBinomial(const arma::uvec &n) : Binomial(n, false, true, "identity"){};
    virtual ~IdentityBinomial() = default;
    virtual const bool valid_link(const arma::mat &x) const;
    virtual Family* clone() const
    {
      if(use_dependence)
      {
        return new IdentityBinomial(this->n_upper, this->copula_object);
      } else {
        return new IdentityBinomial(this->n_upper);
      }
    };
};


const double IdentityBinomial::inverse_link(const double x) const
{
    double value = n_upper(index_n_upper) * x; // x is already in the range [0, 1]
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
    return value;
}

const double IdentityBinomial::link(const double x) const
{
    double value = x / n_upper(index_n_upper);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
    return value;
}

const double IdentityBinomial::derivative_inverse_link(const double x) const
{
  double value = n_upper(index_n_upper);
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
  return value;
}

const double IdentityBinomial::observation_trafo(const double x) const
{
    double value = x / n_upper(index_n_upper);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
    return value;
}


const double IdentityBinomial::link_trafo(const double x) const
{
    return x;
}

const double IdentityBinomial::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool IdentityBinomial::valid_link(const arma::mat &x) const
{
  bool is_valid = x.is_finite() && arma::all(arma::vectorise(x) >= 0.0);
  if(is_valid)
  {
    if(n_upper.n_elem > 1)
    {
      return arma::all(arma::all(x <= arma::repmat(n_upper, 1, x.n_cols)));
    } else {
      return arma::all(arma::vectorise(x) <= n_upper(0));
    }
  }
  return false;
}



/*
  Implementation of the Logit link for the Binomial family, i.e. g(mu) = log(mu / (n - mu))
  Observation transformation: h(y) = y/n
  Link transformation: h(psi) = exp(psi) / (1 + exp(psi)) = mu/n
*/

class LogitBinomial : public Binomial {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    LogitBinomial(const arma::uvec &n, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : Binomial(n, false, false, copula_obj, "logit"){};
    LogitBinomial(const arma::uvec &n) : Binomial(n, false, false, "logit"){};
    virtual ~LogitBinomial() = default;
    virtual Family* clone() const
    {
      if(use_dependence)
      {
        return new LogitBinomial(this->n_upper, this->copula_object);
      } else {
        return new LogitBinomial(this->n_upper);
      }
    };

    virtual const bool valid_link(const arma::mat &x) const;
};

const double LogitBinomial::inverse_link(const double x) const
{
  double value = n_upper(index_n_upper) / (1.0 + exp(-x));
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
  return value;
}

const double LogitBinomial::link(const double x) const
{
  double p = x / n_upper(index_n_upper);
  double value = log(p / (1.0 - p));
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
  return value;
}

const double LogitBinomial::derivative_inverse_link(const double x) const
{
  double value = n_upper(index_n_upper) / (1.0 + exp(-x));
  value *= 1.0 - 1.0 / (1 + exp(-x));
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
  return value;
}

const double LogitBinomial::link_trafo(const double x) const
{
  return 1.0 / (1.0 + exp(-x));
}

const double LogitBinomial::derivative_link_trafo(const double x) const
{
  double value = 1.0 / (1.0 + exp(-x));
  value *= 1.0 - 1.0 / (1 + exp(-x));
  return value;
}

const bool LogitBinomial::valid_link(const arma::mat &x) const
{
  return true;
}


const double LogitBinomial::observation_trafo(const double x) const
{
  double value = x / n_upper(index_n_upper);
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem;
  return value;
}


/*
  Implementation of the Probit link for the Binomial family, i.e., g(mu) = Phi^{-1}(mu/n)
  Observation transformation: h(y) = y/n
  Link transformation: h(psi) = Phi(psi) = mu/n
*/


class ProbitBinomial : public Binomial {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    ProbitBinomial(const arma::uvec &n, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : Binomial(n, false, false, copula_obj, "probit"){};
    ProbitBinomial(const arma::uvec &n) : Binomial(n, false, false, "probit"){};
    virtual ~ProbitBinomial() = default;
    virtual Family* clone() const
    {
      if(use_dependence)
      {
        return new ProbitBinomial(this->n_upper, this->copula_object);
      } else {
        return new ProbitBinomial(this->n_upper);
      }
    };
    virtual const bool valid_link(const arma::mat &x) const;
};

const double ProbitBinomial::inverse_link(const double x) const
{
    double value = n_upper(index_n_upper) * R::pnorm( x, 0.0, 1.0, true, false);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
    return value;
}

const double ProbitBinomial::link(const double x) const
{
    double value = R::qnorm( x / n_upper(index_n_upper), 0.0, 1.0, true, false);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
    return value;
}

const double ProbitBinomial::derivative_inverse_link(const double x) const
{
  double value = n_upper(index_n_upper) * R::dnorm( x, 0.0, 1.0, false);
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
  return value;
}

const double ProbitBinomial::link_trafo(const double x) const
{
    return R::pnorm( x, 0.0, 1.0, true, false);;
}

const double ProbitBinomial::derivative_link_trafo(const double x) const
{
  return R::dnorm( x, 0.0, 1.0, false);
}

const bool ProbitBinomial::valid_link(const arma::mat &x) const
{
  return true;
}


const double ProbitBinomial::observation_trafo(const double x) const
{
  double value = x / n_upper(index_n_upper);
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem;
  return value;
}



/*
  Factory methods for Binomial family
*/


Binomial* Binomial::create(const Rcpp::List &family, Rcpp::S4& copula_obj)
{
  std::string link = Rcpp::as<std::string>(family["link"]);
  Rcpp::IntegerVector size = family["size"];
  arma::uvec n = Rcpp::as<arma::uvec>(size);
  if(link == "identity")
  {
    return new IdentityBinomial(n, copula_obj);
  } else if(link == "softclipping")
  {
    Rcpp::NumericVector temp_const = family["const"];
    double constant = temp_const[0];
    return new SoftClippingBinomial(n, constant, copula_obj);
  }
  else if(link == "logit")
  {
    return new LogitBinomial(n, copula_obj);
  } else if(link == "probit")
  {
    return new ProbitBinomial(n, copula_obj);
  } else {
    throw std::invalid_argument("The desired link is currently not yet implemented.");
  }
}

Binomial* Binomial::create(const Rcpp::List &family)
{
  std::string link = Rcpp::as<std::string>(family["link"]);
  Rcpp::IntegerVector size = family["size"];
  arma::uvec n = Rcpp::as<arma::uvec>(size);
  if(link == "identity")
  {
    return new IdentityBinomial(n);
  } else if(link == "softclipping")
  {
    Rcpp::NumericVector temp_const = family["const"];
    double constant = temp_const[0];
    return new SoftClippingBinomial(n, constant);
  } else if(link == "logit") {
      return new LogitBinomial(n);
  } else if(link == "probit") {
      return new ProbitBinomial(n);
  } else {
      throw std::invalid_argument("The desired link is currently not yet implemented.");
  }
}
