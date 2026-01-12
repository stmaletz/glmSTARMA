/* 
-----------------------------------------------------------------------------
    File: family_quasipoisson.cpp
    Purpose: Implementation of QuasiPoisson family for different link functions
    Author: Steffen Maletz
    Last modified: 2025-12-06
-----------------------------------------------------------------------------
*/

#include "glmstarma.h"


/*
    Random number generators for Quasi-Poisson familys
*/

const arma::vec QuasiPoisson::random_observation_independent(const arma::vec &expectation) const
{
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    arma::vec uniforms(expectation.n_elem, arma::fill::randu);
    double disp;
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        disp = this->const_dispersion ? this->dispersion : this->dispersion_matrix(i, 0);
        observations(i) = (this->*random_generator)(expectation(i), disp, uniforms(i));
    }
    return observations;
}

const arma::vec QuasiPoisson::random_observation(const arma::vec &expectation) const
{
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    arma::vec copula_values(expectation.n_elem);
    copula_values = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));
    double disp;
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        disp = this->const_dispersion ? this->dispersion : this->dispersion_matrix(i, 0);
        observations(i) = (this->*random_generator)(expectation(i), disp, copula_values(i));
    }
    return observations;
}

const arma::vec QuasiPoisson::random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const
{
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    arma::vec uniforms(expectation.n_elem, arma::fill::randu);
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        observations(i) = (this->*random_generator)(expectation(i), dispersion_(i), uniforms(i));
    }
    return observations;
}


const arma::vec QuasiPoisson::random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const
{
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    arma::vec copula_values(expectation.n_elem);
    copula_values = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        observations(i) = (this->*random_generator)(expectation(i), dispersion_(i), copula_values(i));
    }
    return observations;
}


/*
  Quasi-Log-Likelihood for Quasi-Poisson family
*/


const double QuasiPoisson::log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const 
{
    double value = std::log(dispersion_) + 2.0 * observation + 2.0 * std::lgamma(observation + 1);
    if(observation == 0.0)
    {
        value += 2.0 * expectation / dispersion_;
    } else {
        value -= 2.0 * observation * std::log(observation);
        value += 2.0 * (observation * std::log(observation / expectation) - (observation - expectation)) / dispersion_;
    }
    return -0.5 * value;
}


/*
  Residuals for Quasi-Poisson family
*/


const double QuasiPoisson::deviance_residual(const double &observation, const double &expectation) const
{
    if(observation == 0.0)
    {
        return 2.0 * expectation;
    } else {
        return 2.0 * (observation * std::log(observation / expectation) + expectation - observation);
    }
}

const double QuasiPoisson::pearson_residual(const double &observation, const double &expectation) const
{
    return std::pow(observation - expectation, 2) / expectation;
}



/*
  Additional methods for Quasi-Poisson family
*/

const bool QuasiPoisson::valid_expectation(const arma::mat &x) const 
{
  return x.is_finite() && arma::all(arma::vectorise(x) > 0.0);
}

const arma::mat QuasiPoisson::variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const
{
    arma::mat result = this->inverse_link(link_values);
    return result % dispersion_;
}

/*
    Define various link functions for the Quasi-Poisson family
*/



/*
  Define linear model, i.e. the 'identity' link of the Poisson distribution.
  This results in the extended linear PSTARMA model
  Link function: g(mu) = mu
  Observation transformation: h(y) = y
  Link transformation: h(psi) = psi
*/
class LinearQuasiPoisson : public QuasiPoisson {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    LinearQuasiPoisson(const Rcpp::RObject &dispersion, const std::string& sampling_method, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : QuasiPoisson(sampling_method, dispersion, true, true, copula_obj, "identity"){};
    LinearQuasiPoisson(const Rcpp::RObject &dispersion, const std::string& sampling_method) : QuasiPoisson(sampling_method, dispersion, true, true, "identity"){};
    ~LinearQuasiPoisson() = default;
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
          return new LinearQuasiPoisson(disp_vec, sampling_method, this->copula_object);
        } else {
          return new LinearQuasiPoisson(disp_vec, sampling_method);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new LinearQuasiPoisson(disp_mat, sampling_method, this->copula_object);
        } else {
          return new LinearQuasiPoisson(disp_mat, sampling_method);
        }
      }
    };
};

const double LinearQuasiPoisson::inverse_link(const double x) const
{
  return x;
}

const double LinearQuasiPoisson::link(const double x) const
{
  return x;
}

const double LinearQuasiPoisson::derivative_inverse_link(const double x) const
{
  return 1.0;
}


const double LinearQuasiPoisson::observation_trafo(const double x) const
{
  return x;
}

const double LinearQuasiPoisson::link_trafo(const double x) const
{
  return x;
}

const double LinearQuasiPoisson::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool LinearQuasiPoisson::valid_link(const arma::mat &x) const
{
  return x.is_finite() && arma::all(arma::vectorise(x) > 0.0);
}


/*
  Define log-linear model, i.e. the 'log' link of the Poisson distribution.
  This results in the log-linear PSTARMA model
  Link function: g(mu) = log(mu)
  Observation transformation: h(y) = log(y + constant)
  Link transformation: h(psi) = psi
*/

class LogQuasiPoisson : public QuasiPoisson {
  protected:
    const double added_constant;
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    LogQuasiPoisson(const Rcpp::RObject &dispersion, const std::string& sampling_method, double constant, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : QuasiPoisson(sampling_method, dispersion, true, false, copula_obj, "log"), added_constant(constant){};
    LogQuasiPoisson(const Rcpp::RObject &dispersion, const std::string& sampling_method, double constant) : QuasiPoisson(sampling_method, dispersion, true, false, "log"), added_constant(constant){};
    ~LogQuasiPoisson() = default;
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
          return new LogQuasiPoisson(disp_vec, sampling_method, this->added_constant, this->copula_object);
        } else {
          return new LogQuasiPoisson(disp_vec, sampling_method, this->added_constant);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new LogQuasiPoisson(disp_mat, sampling_method, this->added_constant, this->copula_object);
        } else {
          return new LogQuasiPoisson(disp_mat, sampling_method, this->added_constant);
        }
      }
    };
};

const double LogQuasiPoisson::inverse_link(const double x) const
{
  return exp(x);
}

const double LogQuasiPoisson::link(const double x) const
{
  return std::log(x);
}

const double LogQuasiPoisson::derivative_inverse_link(const double x) const
{
  return exp(x);
}


const double LogQuasiPoisson::observation_trafo(const double x) const
{
  return log(x + added_constant);
}

const double LogQuasiPoisson::link_trafo(const double x) const
{
  return x;
}

const double LogQuasiPoisson::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool LogQuasiPoisson::valid_link(const arma::mat &x) const
{
  return true;
}



/*
  Define square root model, i.e. the 'sqrt' link of the Poisson distribution.
  Link function: g(mu) = sqrt(mu)
  Observation transformation: h(y) = 2.0 * sqrt(y + 3/8) Anscombe transformation (of Poisson for approximation)
  Link transformation: h(psi) = psi
*/

class SqrtQuasiPoisson : public QuasiPoisson {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    SqrtQuasiPoisson(const Rcpp::RObject &dispersion, const std::string& sampling_method, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : QuasiPoisson(sampling_method, dispersion, true, true, copula_obj, "sqrt"){};
    SqrtQuasiPoisson(const Rcpp::RObject &dispersion, const std::string& sampling_method) : QuasiPoisson(sampling_method, dispersion, true, true, "sqrt"){};
    ~SqrtQuasiPoisson() = default;
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
          return new SqrtQuasiPoisson(disp_vec, sampling_method, this->copula_object);
        } else {
          return new SqrtQuasiPoisson(disp_vec, sampling_method);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new SqrtQuasiPoisson(disp_mat, sampling_method, this->copula_object);
        } else {
          return new SqrtQuasiPoisson(disp_mat, sampling_method);
        }
      }
    };


};


const double SqrtQuasiPoisson::inverse_link(const double x) const
{
  return pow(x, 2);
}

const double SqrtQuasiPoisson::link(const double x) const
{
  return std::sqrt(x);
}

const double SqrtQuasiPoisson::derivative_inverse_link(const double x) const
{
  return 2 * x;
}

const double SqrtQuasiPoisson::observation_trafo(const double x) const
{
  return 2.0 * std::sqrt(x + 3.0/8.0);
}

const double SqrtQuasiPoisson::link_trafo(const double x) const
{
  return x;
}

const double SqrtQuasiPoisson::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool SqrtQuasiPoisson::valid_link(const arma::mat &x) const
{
  return x.is_finite() && arma::all(arma::vectorise(x) > 0.0);
}

/*
  Define approximately linear model, i.e. the 'softplus' link of the Poisson distribution from the work of Jahn et al. (2023)
  Link function: g(mu) = log(1 + exp(mu/constant)) * constant
  Observation transformation: h(y) = y
  Link transformation: h(psi) = log(1 + exp(psi/constant)) * constant = mu
*/

class SoftPlusQuasiPoisson : public QuasiPoisson {
  protected:
    const double tuning_param;
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    SoftPlusQuasiPoisson(const Rcpp::RObject &dispersion, const std::string& sampling_method, double constant, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : QuasiPoisson(sampling_method, dispersion, false, false, copula_obj, "softplus"), tuning_param(constant){};
    SoftPlusQuasiPoisson(const Rcpp::RObject &dispersion, const std::string& sampling_method, double constant) : QuasiPoisson(sampling_method, dispersion, false, false, "softplus"), tuning_param(constant){};
    ~SoftPlusQuasiPoisson() = default;
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
          return new SoftPlusQuasiPoisson(disp_vec, sampling_method, tuning_param, this->copula_object);
        } else {
          return new SoftPlusQuasiPoisson(disp_vec, sampling_method, tuning_param);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new SoftPlusQuasiPoisson(disp_mat, sampling_method, tuning_param, this->copula_object);
        } else {
          return new SoftPlusQuasiPoisson(disp_mat, sampling_method, tuning_param);
        }
      }
    };
};



const double SoftPlusQuasiPoisson::inverse_link(const double x) const
{
  return tuning_param * std::log(1.0 + exp(x / tuning_param));
}

const double SoftPlusQuasiPoisson::link(const double x) const
{
  return tuning_param * std::log(exp(x / tuning_param) - 1.0);
}

const double SoftPlusQuasiPoisson::derivative_inverse_link(const double x) const
{
  double val = exp(x / tuning_param);
  return val / (1.0 + val);
}


const double SoftPlusQuasiPoisson::observation_trafo(const double x) const
{
  return x;
}

const double SoftPlusQuasiPoisson::link_trafo(const double x) const
{
  return tuning_param * std::log(1.0 + exp(x / tuning_param));
}

const double SoftPlusQuasiPoisson::derivative_link_trafo(const double x) const
{
  double val = exp(x / tuning_param);
  return val / (1.0 + val);
}

const bool SoftPlusQuasiPoisson::valid_link(const arma::mat &x) const
{
  return true;
}



/*
  Factory methods for Quasi-Poisson family
*/
QuasiPoisson* QuasiPoisson::create(const Rcpp::List &family, Rcpp::S4& copula_obj)
{
    std::string link = Rcpp::as<std::string>(family["link"]);
    std::string sampling_method = Rcpp::as<std::string>(family["sampling_method"]);
    Rcpp::RObject dispersion = family["dispersion"];
    if(link == "log")
    {
      double constant = 1.0;
      if(family.containsElementNamed("const"))
      {
        Rcpp::NumericVector temp_const = family["const"];
        constant = temp_const[0];
      }
      return new LogQuasiPoisson(dispersion, sampling_method, constant, copula_obj);
    } 
    else if(link == "identity") 
    {
      return new LinearQuasiPoisson(dispersion, sampling_method, copula_obj);
    } 
    else if(link == "sqrt")
    {
      return new SqrtQuasiPoisson(dispersion, sampling_method, copula_obj);
    }
    else if(link == "softplus") 
    {
      double constant = 1.0;
      if(family.containsElementNamed("const"))
      {
        Rcpp::NumericVector temp_const = family["const"];
        constant = temp_const[0];
      }
      return new SoftPlusQuasiPoisson(dispersion, sampling_method, constant, copula_obj);
    } 
    else 
    {
      throw std::invalid_argument("The desired link is currently not yet implemented.");
    }
}
QuasiPoisson* QuasiPoisson::create(const Rcpp::List &family)
{
  std::string link = Rcpp::as<std::string>(family["link"]);
  Rcpp::RObject dispersion = family["dispersion"];
  std::string sampling_method = Rcpp::as<std::string>(family["sampling_method"]);
  if(link == "log")
  {
    double constant = 1.0;
    if(family.containsElementNamed("const"))
    {
      Rcpp::NumericVector temp_const = family["const"];
      constant = temp_const[0];
    }
    return new LogQuasiPoisson(dispersion, sampling_method, constant);
  } 
  else if(link == "identity") 
  {
    return new LinearQuasiPoisson(dispersion, sampling_method);
  } 
  else if(link == "sqrt") 
  {
    return new SqrtQuasiPoisson(dispersion, sampling_method);
  } 
  else if(link == "softplus") 
  {
    double constant = 1.0;
    if(family.containsElementNamed("const"))
    {
      Rcpp::NumericVector temp_const = family["const"];
      constant = temp_const[0];
    }
    return new SoftPlusQuasiPoisson(dispersion, sampling_method, constant);
  } 
  else
  {
    throw std::invalid_argument("The desired link is currently not yet implemented.");
  }
}
