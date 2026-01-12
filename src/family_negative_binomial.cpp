/* 
-----------------------------------------------------------------------------
    File: family_negative_binomial.cpp
    Purpose: Implementation of Negative Binomial family for different link functions
    Author: Steffen Maletz
    Last modified: 2025-12-06
-----------------------------------------------------------------------------
*/

#include "glmstarma.h"


/*
  Random number generators for Negative Binomial family
  These either use functions from R's stats package or use inversion method on copula samples
*/
const arma::vec NegativeBinomial::random_observation_independent(const arma::vec &expectation) const
{
  if(!this->const_dispersion && dispersion_matrix.n_cols != 1)
  {
    Rcpp::stop("Dispersion is not constant. Cannot generate new observations without new dispersion values.");
  }
  arma::vec observations(expectation.n_elem, arma::fill::zeros);
  double disp;
  for(unsigned int i = 0; i < expectation.n_elem; i++)
  {
    disp = this->const_dispersion ? dispersion : dispersion_matrix(i, 0);
    if(disp == 0.0)
    {
      observations(i) = Rcpp::as<double>(Rcpp::rpois(1, expectation(i)));
    } else {
      observations(i) = Rcpp::as<double>(Rcpp::rnbinom_mu(1, 1.0 / dispersion, expectation(i)));
    }
  }
  return observations;
}

const arma::vec NegativeBinomial::random_observation(const arma::vec &expectation) const
{
    if(!this->const_dispersion && dispersion_matrix.n_cols != 1)
    {
      Rcpp::stop("Dispersion is not constant. Cannot generate new observations without new dispersion values.");
    }
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    arma::vec copula_values(expectation.n_elem);
    double disp;
    copula_values = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        disp = this->const_dispersion ? dispersion : dispersion_matrix(i, 0);
        if(disp == 0.0)
        {
          observations(i) = R::qpois( copula_values(i), expectation(i), true, false );
        } else {
          observations(i) = R::qnbinom_mu( copula_values(i), disp, expectation(i), true, false );
        }
    }
    return observations;
}


const arma::vec NegativeBinomial::random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const
{
  arma::vec observations(expectation.n_elem, arma::fill::zeros);
  for(unsigned int i = 0; i < expectation.n_elem; i++)
  {
    if(dispersion_(i) == 0.0)
    {
      observations(i) = Rcpp::as<double>(Rcpp::rpois(1, expectation(i)));
    } else {
      observations(i) = Rcpp::as<double>(Rcpp::rnbinom_mu(1, 1.0 / dispersion_(i), expectation(i)));
    }
  }
  return observations;
}

const arma::vec NegativeBinomial::random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const
{
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    arma::vec copula_values(expectation.n_elem);
    copula_values = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        if(dispersion_(i) == 0.0)
        {
          observations(i) = R::qpois( copula_values(i), expectation(i), true, false );
        } else {
          observations(i) = R::qnbinom_mu( copula_values(i), dispersion_(i), expectation(i), true, false );
        }
    }
    return observations;
}

/*
  Log-likelihood for Negative Binomial family
*/

const double NegativeBinomial::log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const 
{
  if(dispersion_ == 0.0)
  {
    return R::dpois( observation, expectation, true );
  }
  return R::dnbinom_mu(observation, 1.0 / dispersion_, expectation, true);
}

/*
  Deviance and Pearson residuals for Negative Binomial family
*/

const double NegativeBinomial::deviance_residual(const double &observation, const double &expectation) const
{
  return 1.0; // Not intended to be used
}

const double NegativeBinomial::pearson_residual(const double &observation, const double &expectation) const
{
    return std::pow(observation - expectation, 2) / expectation;
}


/*
  Additional methods for Negative Binomial family
*/

const bool NegativeBinomial::valid_expectation(const arma::mat &x) const 
{
  return x.is_finite() && arma::all(arma::vectorise(x) > 0.0);
}

const arma::mat NegativeBinomial::variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const
{
    arma::mat response = this->inverse_link(link_values);
    return response + (response % response) % dispersion_;
}


/*
  Implementation of various link functions for Negative Binomial family
*/

/*
  Define linear, i.e. identity-model, i.e. link function: g(mu) = mu
  Observation transformation: h(y) = y
  Link transformation: h(psi) = psi
  This corresponds to the linear model from Maletz et al. (2024), but with negative binomial marginals.
*/

class LinearNegativeBinomial : public NegativeBinomial {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    LinearNegativeBinomial(const Rcpp::RObject &dispersion, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : NegativeBinomial(dispersion, true, true, copula_obj, "identity"){};
    LinearNegativeBinomial(const Rcpp::RObject &dispersion) : NegativeBinomial(dispersion, true, true, "identity"){};
    ~LinearNegativeBinomial() = default;
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
          return new LinearNegativeBinomial(disp_vec, this->copula_object);
        } else {
          return new LinearNegativeBinomial(disp_vec);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new LinearNegativeBinomial(disp_mat, this->copula_object);
        } else {
          return new LinearNegativeBinomial(disp_mat);
        }
      }
    };
};

const double LinearNegativeBinomial::inverse_link(const double x) const
{
  return x;
}

const double LinearNegativeBinomial::link(const double x) const
{
  return x;
}

const double LinearNegativeBinomial::derivative_inverse_link(const double x) const
{
  return 1.0;
}


const double LinearNegativeBinomial::observation_trafo(const double x) const
{
  return x;
}

const double LinearNegativeBinomial::link_trafo(const double x) const
{
  return x;
}

const double LinearNegativeBinomial::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool LinearNegativeBinomial::valid_link(const arma::mat &x) const
{
  return x.is_finite() && arma::all(arma::vectorise(x) > 0.0);
}


/*
  Define log-linear model, i.e. the 'log' link of the Negative Binomial distribution, i.e., link function: g(mu) = log(mu)
  Observation transformation: h(y) = log(y + 1)
  Link transformation: h(psi) = psi
*/

class LogNegativeBinomial : public NegativeBinomial {
  protected:
    const double added_constant;
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    LogNegativeBinomial(const Rcpp::RObject &dispersion, double constant, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : NegativeBinomial(dispersion, true, false, copula_obj, "log"), added_constant(constant){};
    LogNegativeBinomial(const Rcpp::RObject &dispersion, double constant) : NegativeBinomial(dispersion, true, false, "log"), added_constant(constant){};
    ~LogNegativeBinomial() = default;
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
          return new LogNegativeBinomial(disp_vec, added_constant, this->copula_object);
        } else {
          return new LogNegativeBinomial(disp_vec, added_constant);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new LogNegativeBinomial(disp_mat, added_constant, this->copula_object);
        } else {
          return new LogNegativeBinomial(disp_mat, added_constant);
        }
      }
    };
};

const double LogNegativeBinomial::inverse_link(const double x) const
{
  return exp(x);
}

const double LogNegativeBinomial::link(const double x) const
{
  return std::log(x);
}


const double LogNegativeBinomial::derivative_inverse_link(const double x) const
{
  return exp(x);
}


const double LogNegativeBinomial::observation_trafo(const double x) const
{
  return log(x + added_constant);
}

const double LogNegativeBinomial::link_trafo(const double x) const
{
  return x;
}

const double LogNegativeBinomial::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool LogNegativeBinomial::valid_link(const arma::mat &x) const
{
  return true;
}



/*
  Define square root model, i.e. the 'sqrt' link of the Negative Binomial distribution, i.e., link function: g(mu) = sqrt(mu)
  Observation transformation: h(y) = 2.0 * sqrt(y + 3/8) Anscombe transformation (of Poisson case, as approximation)
  Link transformation: h(psi) = psi
*/

class SqrtNegativeBinomial : public NegativeBinomial {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    SqrtNegativeBinomial(const Rcpp::RObject &dispersion, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : NegativeBinomial(dispersion, true, true, copula_obj, "sqrt"){};
    SqrtNegativeBinomial(const Rcpp::RObject &dispersion) : NegativeBinomial(dispersion, true, true, "sqrt"){};
    ~SqrtNegativeBinomial() = default;
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
          return new SqrtNegativeBinomial(disp_vec, this->copula_object);
        } else {
          return new SqrtNegativeBinomial(disp_vec);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new SqrtNegativeBinomial(disp_mat, this->copula_object);
        } else {
          return new SqrtNegativeBinomial(disp_mat);
        }
      }
    };
};

const double SqrtNegativeBinomial::inverse_link(const double x) const
{
  return pow(x, 2);
}

const double SqrtNegativeBinomial::link(const double x) const
{
  return std::sqrt(x);
}

const double SqrtNegativeBinomial::derivative_inverse_link(const double x) const
{
  return 2 * x;
}

const double SqrtNegativeBinomial::observation_trafo(const double x) const
{
  return std::sqrt(x + 3.0 / 8.0) * 2.0;
}

const double SqrtNegativeBinomial::link_trafo(const double x) const
{
  return x;
}

const double SqrtNegativeBinomial::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool SqrtNegativeBinomial::valid_link(const arma::mat &x) const
{
  return x.is_finite() && arma::all(arma::vectorise(x) > 0.0);
}

/*
  Define approximately linear model, i.e. the 'softplus' link of the Poisson distribution from the work of Jahn et al. (2023), but for Negative Binomial marginals,
  link function: g(mu) = log(1 + exp(mu/constant)) * constant
  Observation transformation: h(y) = y
  Link transformation: h(psi) = log(1 + exp(psi/constant)) * constant = mu
*/

class SoftPlusNegativeBinomial : public NegativeBinomial {
  protected:
    const double tuning_param;
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    SoftPlusNegativeBinomial(const Rcpp::RObject &dispersion, double constant, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : NegativeBinomial(dispersion, false, false, copula_obj, "softplus"), tuning_param(constant){};
    SoftPlusNegativeBinomial(const Rcpp::RObject &dispersion, double constant) : NegativeBinomial(dispersion, false, false, "softplus"), tuning_param(constant){};
    ~SoftPlusNegativeBinomial() = default;
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
          return new SoftPlusNegativeBinomial(disp_vec, tuning_param, this->copula_object);
        } else {
          return new SoftPlusNegativeBinomial(disp_vec, tuning_param);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new SoftPlusNegativeBinomial(disp_mat, tuning_param, this->copula_object);
        } else {
          return new SoftPlusNegativeBinomial(disp_mat, tuning_param);
        }
      }
    };
};


const double SoftPlusNegativeBinomial::inverse_link(const double x) const
{
  return tuning_param * std::log(1.0 + exp(x / tuning_param));
}

const double SoftPlusNegativeBinomial::link(const double x) const
{
  return tuning_param * std::log(exp(x / tuning_param) - 1.0);
}

const double SoftPlusNegativeBinomial::derivative_inverse_link(const double x) const
{
  double val = exp(x / tuning_param);
  return val / (1.0 + val);
}


const double SoftPlusNegativeBinomial::observation_trafo(const double x) const
{
  return x;
}

const double SoftPlusNegativeBinomial::link_trafo(const double x) const
{
  return tuning_param * std::log(1.0 + exp(x / tuning_param));
}

const double SoftPlusNegativeBinomial::derivative_link_trafo(const double x) const
{
  double val = exp(x / tuning_param);
  return val / (1.0 + val);
}

const bool SoftPlusNegativeBinomial::valid_link(const arma::mat &x) const
{
  return true;
}

/*
  Factory methods for Negative Binomial family
*/


NegativeBinomial* NegativeBinomial::create(const Rcpp::List &family, Rcpp::S4& copula_obj)
{
  Rcpp::CharacterVector link = family["link"];
    Rcpp::RObject dispersion = family["dispersion"];
    if(link[0] == "log")
    {
      double constant = 1.0;
      if(family.containsElementNamed("const"))
      {
        Rcpp::NumericVector temp_const = family["const"];
        constant = temp_const[0];
      }
      return new LogNegativeBinomial(dispersion, constant, copula_obj);
    } 
    else if(link[0] == "identity") 
    {
      return new LinearNegativeBinomial(dispersion, copula_obj);
    } 
    else if(link[0] == "sqrt")
    {
      return new SqrtNegativeBinomial(dispersion, copula_obj);
    }
    else if(link[0] == "softplus") 
    {
      double constant = 1.0;
      if(family.containsElementNamed("const"))
      {
        Rcpp::NumericVector temp_const = family["const"];
        constant = temp_const[0];
      }
      return new SoftPlusNegativeBinomial(dispersion, constant, copula_obj);
    } 
    else 
    {
      throw std::invalid_argument("The desired link is currently not yet implemented.");
    }
}

NegativeBinomial* NegativeBinomial::create(const Rcpp::List &family)
{
  Rcpp::CharacterVector link = family["link"];
  Rcpp::RObject dispersion = family["dispersion"];
  if(link[0] == "log")
  {
    double constant = 1.0;
    if(family.containsElementNamed("const"))
    {
      Rcpp::NumericVector temp_const = family["const"];
      constant = temp_const[0];
    }
    return new LogNegativeBinomial(dispersion, constant);
  } 
  else if(link[0] == "identity") 
  {
    return new LinearNegativeBinomial(dispersion);
  } 
  else if(link[0] == "sqrt") 
  {
    return new SqrtNegativeBinomial(dispersion);
  } 
  else if(link[0] == "softplus") 
  {
    double constant = 1.0;
    if(family.containsElementNamed("const"))
    {
      Rcpp::NumericVector temp_const = family["const"];
      constant = temp_const[0];
    }
    return new SoftPlusNegativeBinomial(dispersion, constant);
  } 
  else
  {
    throw std::invalid_argument("The desired link is currently not yet implemented.");
  }
}