/* 
-----------------------------------------------------------------------------
    File: family_quasibinomial.cpp
    Purpose: Implementation of QuasiBinomial family for different link functions
    Author: Steffen Maletz
    Last modified: 2025-12-06
-----------------------------------------------------------------------------
*/

#include "glmstarma.h"


/*
    Random number generators for Quasi-Binomial family
*/

const arma::vec QuasiBinomial::random_observation_independent(const arma::vec &expectation) const
{
    if(!this->const_dispersion && dispersion_matrix.n_cols != 1)
    {
        Rcpp::stop("Dispersion is not constant. Cannot generate new observations without new dispersion values.");
    }
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    unsigned int n;
    double disp;
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        n = n_upper.n_elem > 1 ? n_upper(i) : n_upper(0);
        disp = const_dispersion ? dispersion : dispersion_matrix(i, 0);
        observations(i) = rquasi(n, expectation(i) / n, disp);
    }
    return observations;
}


const arma::vec QuasiBinomial::random_observation(const arma::vec &expectation) const
{
    if(!this->const_dispersion && dispersion_matrix.n_cols != 1)
    {
        Rcpp::stop("Dispersion is not constant. Cannot generate new observations without new dispersion values.");
    }
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    unsigned int n_max = n_upper.n_elem > 1 ? arma::max(n_upper) : n_upper(0);
    arma::mat uniforms1(expectation.n_elem, n_max + 1);
    arma::mat uniforms2(expectation.n_elem, n_max);
    for(unsigned int i = 0; i < uniforms1.n_cols; i++)
    {
      uniforms1.col(i) = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));
    }
    uniforms1 = uniforms1.t();
    for(unsigned int i = 0; i < uniforms2.n_cols; i++)
    {
      uniforms2.col(i) = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));
    }
    uniforms2 = uniforms2.t();


    unsigned int n;
    double disp;
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        n = n_upper.n_elem > 1 ? n_upper(i) : n_upper(0);
        disp = const_dispersion ? dispersion : dispersion_matrix(i, 0);
        observations(i) = rquasi_dependent(n, expectation(i) / n, disp, uniforms1.col(i).head(n + 1), uniforms2.col(i).head(n));
    }
    return observations;
}

const arma::vec QuasiBinomial::random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const
{
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    unsigned int n;
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        n = n_upper.n_elem > 1 ? n_upper(i) : n_upper(0);
        observations(i) = rquasi(n, expectation(i) / n, dispersion_(i));
    }
    return observations;
}

const arma::vec QuasiBinomial::random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const
{
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    unsigned int n;
    unsigned int n_max = n_upper.n_elem > 1 ? arma::max(n_upper) : n_upper(0);
    arma::mat uniforms1(expectation.n_elem, n_max + 1);
    arma::mat uniforms2(expectation.n_elem, n_max);
    for(unsigned int i = 0; i < uniforms1.n_cols; i++)
    {
      uniforms1.col(i) = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));
    }
    uniforms1 = uniforms1.t();
    for(unsigned int i = 0; i < uniforms2.n_cols; i++)
    {
      uniforms2.col(i) = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));
    }
    uniforms2 = uniforms2.t();
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        n = n_upper.n_elem > 1 ? n_upper(i) : n_upper(0);
        observations(i) = rquasi_dependent(n, expectation(i) / n, dispersion_(i), uniforms1.col(i).head(n + 1), uniforms2.col(i).head(n));
    }
    return observations;
}


/*
  Quasi-Log-Likelihood for Quasi-Binomial family
*/

const double QuasiBinomial::log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const 
{
    int n = n_upper(index_n_upper);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
    double log_comb = std::lgamma(n + 1) - std::lgamma(observation + 1) - std::lgamma(n - observation + 1);
    double dev = 0.0;
    if (observation > 0)
    {
      log_comb += observation * std::log(observation / n);
      dev += observation * std::log( observation / expectation );
    }
    if (observation < n)
    {
      log_comb += (n - observation) * std::log(1.0 - observation / n);
      dev += (n - observation) * std::log( (n - observation) / (n - expectation) );
    }
    log_comb *= -2.0;
    dev *= 2.0 / dispersion_;
    dev +=  log_comb;
    dev += std::log(dispersion_);
    return -0.5 * dev;
}


/*
  Residuals for Quasi-Binomial family
*/

const double QuasiBinomial::deviance_residual(const double &observation, const double &expectation) const
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

const double QuasiBinomial::pearson_residual(const double &observation, const double &expectation) const
{
    int n = n_upper(index_n_upper);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem;
    double p = expectation / n;
    return std::pow(observation - expectation, 2) / (n * p * (1.0 - p));
}


/*
  Additional methods for Quasi-Binomial family
*/

const bool QuasiBinomial::valid_expectation(const arma::mat &x) const 
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


const arma::mat QuasiBinomial::variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const
{
    arma::mat fitted = this->inverse_link(link_values);
    arma::mat n_up; // (link_values.n_rows, link_values.n_cols);
    arma::vec n_upper_dbl = arma::conv_to<arma::vec>::from(n_upper);
    if(n_upper.n_elem > 1){
        n_up = link_values.n_cols == 1 ? arma::repmat(n_upper_dbl, link_values.n_rows / n_upper_dbl.n_elem, 1) : arma::repmat(n_upper_dbl, 1, link_values.n_cols);
        fitted /= n_up;
    } else {
      fitted = fitted / n_upper(0);
    }
    arma::mat ones = arma::ones(fitted.n_rows, fitted.n_cols);
    fitted %= (ones - fitted);
    if(n_upper.n_elem > 1)
    {
       fitted %= n_up; 
    } else {
        fitted *= n_upper(0);
    }
    fitted %= dispersion_;
    return fitted;
}



/*
    Define various link functions for the Quasi-Binomial family
*/


/*
  Define approximately linear model, i.e. the 'softclipping' link of the binomial distribution from the work of Jahn et al. (2023)
  Link function: g(mu) = log( (exp(mu/constant) - 1) / (1 - exp((mu - n)/constant)) ) * constant
  Observation transformation: h(y) = y
  Link transformation: h(psi) = log( (exp(psi/constant) - 1) / (1 - exp((psi - n)/constant)) ) * constant = mu
*/

class SoftClippingQuasiBinomial : public QuasiBinomial {
  protected:
    const double tuning_param;
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    SoftClippingQuasiBinomial(const Rcpp::RObject &dispersion, const arma::uvec &n, double constant, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : QuasiBinomial(n, dispersion, false, false, copula_obj, "softclipping"), tuning_param(constant){};
    SoftClippingQuasiBinomial(const Rcpp::RObject &dispersion, const arma::uvec &n, double constant) : QuasiBinomial(n, dispersion, false, false, "softclipping"), tuning_param(constant){};
    virtual ~SoftClippingQuasiBinomial() = default;
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
          return new SoftClippingQuasiBinomial(disp_vec, this->n_upper, this->tuning_param, this->copula_object);
        } else {
          return new SoftClippingQuasiBinomial(disp_vec, this->n_upper, this->tuning_param);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new SoftClippingQuasiBinomial(disp_mat, this->n_upper, this->tuning_param, this->copula_object);
        } else {
          return new SoftClippingQuasiBinomial(disp_mat, this->n_upper, this->tuning_param);
        }
      }
    };
};


const double SoftClippingQuasiBinomial::inverse_link(const double x) const
{
    double upper = 1.0 + exp(x / tuning_param);
    double lower = 1.0 + exp((x - 1.0) / tuning_param);
    double value = n_upper(index_n_upper) * tuning_param * std::log(upper / lower);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
    return value;
}

const double SoftClippingQuasiBinomial::link(const double x) const
{
    unsigned int n = n_upper(index_n_upper);
    double upper = std::exp(x / (n * tuning_param)) - 1.0;
    double lower = 1.0 - std::exp((x - n) / (n * tuning_param));
    double value = tuning_param * std::log(upper / lower);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
    return value;
}

const double SoftClippingQuasiBinomial::derivative_inverse_link(const double x) const
{
  double val_first = exp(x / tuning_param) / (1.0 + exp(x / tuning_param));
  double val_second = exp((x - 1.0) / tuning_param) / (1.0 + exp((x - 1.0) / tuning_param));
  double value = n_upper(index_n_upper) * (val_first - val_second);
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem;
  return value;
}

const double SoftClippingQuasiBinomial::observation_trafo(const double x) const
{
  double value = x / n_upper(index_n_upper);
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem;
  return value;
}


const double SoftClippingQuasiBinomial::link_trafo(const double x) const
{
    double value = inverse_link(x) / n_upper(index_n_upper);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem;
    return value;
}

const double SoftClippingQuasiBinomial::derivative_link_trafo(const double x) const
{
  double value = derivative_inverse_link(x) / n_upper(index_n_upper);
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem;
  return value;
}

const bool SoftClippingQuasiBinomial::valid_link(const arma::mat &x) const
{
  return true;
}

/*
  Define identity link for Quasi-Binomial family
  Might be problematic as it must be ensured that the expectation is always in [0, n]
  Probably only possible without covariates
  Link function: g(mu) = mu
  Observation transformation: h(y) = y
  Link transformation: h(psi) = psi
*/

class IdentityQuasiBinomial : public QuasiBinomial {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    IdentityQuasiBinomial(const Rcpp::RObject dispersion, const arma::uvec &n, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : QuasiBinomial(n, dispersion, false, true, copula_obj, "identity"){};
    IdentityQuasiBinomial(const Rcpp::RObject dispersion, const arma::uvec &n) : QuasiBinomial(n, dispersion, false, true, "identity"){};
    virtual ~IdentityQuasiBinomial() = default;
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
          return new IdentityQuasiBinomial(disp_vec, this->n_upper, this->copula_object);
        } else {
          return new IdentityQuasiBinomial(disp_vec, this->n_upper);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new IdentityQuasiBinomial(disp_mat, this->n_upper, this->copula_object);
        } else {
          return new IdentityQuasiBinomial(disp_mat, this->n_upper);
        }
      }
    };
};


const double IdentityQuasiBinomial::inverse_link(const double x) const
{
    double value = n_upper(index_n_upper) * x; // x is already in the range [0, 1]
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
    return value;
}

const double IdentityQuasiBinomial::link(const double x) const
{
    double value = x / n_upper(index_n_upper);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
    return value;
}

const double IdentityQuasiBinomial::derivative_inverse_link(const double x) const
{
  double value = n_upper(index_n_upper);
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
  return value;
}

const double IdentityQuasiBinomial::observation_trafo(const double x) const
{
    double value = x / n_upper(index_n_upper);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
    return value;
}


const double IdentityQuasiBinomial::link_trafo(const double x) const
{
    return x;
}

const double IdentityQuasiBinomial::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool IdentityQuasiBinomial::valid_link(const arma::mat &x) const
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
  Define logit link for Quasi-Binomial family
  Link function: g(mu) = logit(mu) = log(mu / (n - mu))
  Observation transformation: h(y) = y / n
  Link transformation: h(psi) = exp(psi) / (1 + exp(psi)) = mu / n
*/

class LogitQuasiBinomial : public QuasiBinomial {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    LogitQuasiBinomial(const Rcpp::RObject dispersion, const arma::uvec &n, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : QuasiBinomial(n, dispersion, false, false, copula_obj, "logit"){};
    LogitQuasiBinomial(const Rcpp::RObject dispersion, const arma::uvec &n) : QuasiBinomial(n, dispersion, false, false, "logit"){};
    virtual ~LogitQuasiBinomial() = default;
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
          return new LogitQuasiBinomial(disp_vec, this->n_upper, this->copula_object);
        } else {
          return new LogitQuasiBinomial(disp_vec, this->n_upper);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new LogitQuasiBinomial(disp_mat, this->n_upper, this->copula_object);
        } else {
          return new LogitQuasiBinomial(disp_mat, this->n_upper);
        }
      }
    };
};

const double LogitQuasiBinomial::inverse_link(const double x) const
{
  double value = n_upper(index_n_upper) / (1.0 + exp(-x));
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
  return value;
}

const double LogitQuasiBinomial::link(const double x) const
{
  double p = x / n_upper(index_n_upper);
  double value = log(p / (1.0 - p));
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
  return value;
}

const double LogitQuasiBinomial::derivative_inverse_link(const double x) const
{
  double value = n_upper(index_n_upper) / (1.0 + exp(-x));
  value *= 1.0 - 1.0 / (1 + exp(-x));
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
  return value;
}

const double LogitQuasiBinomial::link_trafo(const double x) const
{
  return 1.0 / (1.0 + exp(-x));
}

const double LogitQuasiBinomial::derivative_link_trafo(const double x) const
{
  double value = 1.0 / (1.0 + exp(-x));
  value *= 1.0 - 1.0 / (1 + exp(-x));
  return value;
}

const bool LogitQuasiBinomial::valid_link(const arma::mat &x) const
{
  return true;
}


const double LogitQuasiBinomial::observation_trafo(const double x) const
{
  double value = x / n_upper(index_n_upper);
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem;
  return value;
}


/*
  Define Probit link for Quasi-Binomial family
  Link function: g(mu) = Phi^{-1}(mu / n)
  Observation transformation: h(y) = y / n
  Link transformation: h(psi) = Phi(psi) = mu / n
*/

class ProbitQuasiBinomial : public QuasiBinomial {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    ProbitQuasiBinomial(const Rcpp::RObject dispersion, const arma::uvec &n, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : QuasiBinomial(n, dispersion, false, false, copula_obj, "probit"){};
    ProbitQuasiBinomial(const Rcpp::RObject dispersion, const arma::uvec &n) : QuasiBinomial(n, dispersion, false, false, "probit"){};
    virtual ~ProbitQuasiBinomial() = default;
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
          return new ProbitQuasiBinomial(disp_vec, this->n_upper, this->copula_object);
        } else {
          return new ProbitQuasiBinomial(disp_vec, this->n_upper);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new ProbitQuasiBinomial(disp_mat, this->n_upper, this->copula_object);
        } else {
          return new ProbitQuasiBinomial(disp_mat, this->n_upper);
        }
      }
    };
};

const double ProbitQuasiBinomial::inverse_link(const double x) const
{
    double value = n_upper(index_n_upper) * R::pnorm( x, 0.0, 1.0, true, false);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
    return value;
}

const double ProbitQuasiBinomial::link(const double x) const
{
    double value = R::qnorm( x / n_upper(index_n_upper), 0.0, 1.0, true, false);
    index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
    return value;
}

const double ProbitQuasiBinomial::derivative_inverse_link(const double x) const
{
  double value = n_upper(index_n_upper) * R::dnorm( x, 0.0, 1.0, false);
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem; // Update index for next call
  return value;
}

const double ProbitQuasiBinomial::link_trafo(const double x) const
{
    return R::pnorm( x, 0.0, 1.0, true, false);;
}

const double ProbitQuasiBinomial::derivative_link_trafo(const double x) const
{
  return R::dnorm( x, 0.0, 1.0, false);
}

const bool ProbitQuasiBinomial::valid_link(const arma::mat &x) const
{
  return true;
}


const double ProbitQuasiBinomial::observation_trafo(const double x) const
{
  double value = x / n_upper(index_n_upper);
  index_n_upper = (index_n_upper + 1) % n_upper.n_elem;
  return value;
}



/*
  Factory methods for Quasi-Binomial family
*/


QuasiBinomial* QuasiBinomial::create(const Rcpp::List &family, Rcpp::S4& copula_obj)
{
  
  std::string link = Rcpp::as<std::string>(family["link"]);
  Rcpp::IntegerVector size = family["size"];
  arma::uvec n = Rcpp::as<arma::uvec>(size);
  Rcpp::RObject dispersion = family["dispersion"];
  if(link == "identity")
  {
    return new IdentityQuasiBinomial(dispersion, n, copula_obj);
  } else if(link == "softclipping")
  {
    Rcpp::NumericVector temp_const = family["const"];
    double constant = temp_const[0];
    return new SoftClippingQuasiBinomial(dispersion, n, constant, copula_obj);
  }
  else if(link == "logit")
  {
    return new LogitQuasiBinomial(dispersion, n, copula_obj);
  } else if(link == "probit")
  {
    return new ProbitQuasiBinomial(dispersion, n, copula_obj);
  } else {
    throw std::invalid_argument("The desired link is currently not yet implemented.");
  }
}

QuasiBinomial* QuasiBinomial::create(const Rcpp::List &family)
{
  std::string link = Rcpp::as<std::string>(family["link"]);
  Rcpp::IntegerVector size = family["size"];
  arma::uvec n = Rcpp::as<arma::uvec>(size);
  Rcpp::RObject dispersion = family["dispersion"];
  if(link == "identity")
  {
    return new IdentityQuasiBinomial(dispersion, n);
  } else if(link == "softclipping")
  {
    Rcpp::NumericVector temp_const = family["const"];
    double constant = temp_const[0];
    return new SoftClippingQuasiBinomial(dispersion, n, constant);
  } else if(link == "logit") {
      return new LogitQuasiBinomial(dispersion, n);
  } else if(link == "probit") {
      return new ProbitQuasiBinomial(dispersion, n);
  } else {
      throw std::invalid_argument("The desired link is currently not yet implemented.");
  }
}
