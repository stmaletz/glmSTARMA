/* 
-----------------------------------------------------------------------------
    File: family_inverse_gauss.cpp
    Purpose: Implementation of Inverse Gaussian family for different link functions
    Author: Steffen Maletz
    Last modified: 2026-01-10
-----------------------------------------------------------------------------
*/

#include "glmstarma.h"


/*
  Random number generators for Inverse Gaussian family
  These implement the Michael, Schucany and Haas (1976) method,
  which is applied used with contemporary independent samples 
    or to samples where contemporaneous dependence is generated via a copula.
  See also:
  https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution
  Michael, John R.; Schucany, William R.; Haas, Roy W. (1976), "Generating Random Variates Using Transformations with Multiple Roots", 
    The American Statistician, 30 (2): 88–90, doi:10.1080/00031305.1976.10479147, JSTOR 2683801
  dispersion parameter is inverse of second parameter in Wikipedia notation
*/

const arma::vec InverseGauss::random_observation_independent(const arma::vec &expectation) const
{
    if(!this->const_dispersion && dispersion_matrix.n_cols != 1)
    {
      Rcpp::stop("Dispersion is not constant. Cannot generate new observations without new dispersion values.");
    }
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    arma::vec random_normal(expectation.n_elem, arma::fill::randn);
    random_normal = arma::square(random_normal);
    arma::vec transform = expectation;
    if(this->const_dispersion)
    {
      transform += (arma::square(expectation) % random_normal) * dispersion / (2.0);
      transform -= (expectation * dispersion) / (2.0) % arma::sqrt((4.0 / dispersion) * (expectation % random_normal) + arma::square(expectation) % arma::square(random_normal));
    } else {
      transform += (arma::square(expectation) % random_normal) % dispersion_matrix.col(0) / (2.0);
      transform -= (expectation % dispersion_matrix.col(0) / 2.0) % arma::sqrt((4.0 / dispersion_matrix.col(0)) % (expectation % random_normal) + arma::square(expectation) % arma::square(random_normal));
    }
    arma::vec accept(expectation.n_elem, arma::fill::randu);
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        observations(i) = accept(i) * (expectation(i) + transform(i)) <= expectation(i) ? transform(i) : expectation(i) * expectation(i) / transform(i);
    }
    return observations;
}

const arma::vec InverseGauss::random_observation(const arma::vec &expectation) const
{
    if(!this->const_dispersion && dispersion_matrix.n_cols != 1)
    {
      Rcpp::stop("Dispersion is not constant. Cannot generate new observations without new dispersion values.");
    }
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    arma::vec copula_values(expectation.n_elem);
    copula_values = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));

    arma::vec random_normal(expectation.n_elem, arma::fill::randn);
    random_normal = arma::square(random_normal);
    arma::vec transform = expectation;
    if(this->const_dispersion)
    {
      transform += (arma::square(expectation) % random_normal * dispersion) / (2.0);
      transform -= (expectation * dispersion / (2.0)) % arma::sqrt((4.0 / dispersion) * (expectation % random_normal) + arma::square(expectation) % arma::square(random_normal));
    } else {
      transform += (arma::square(expectation) % random_normal % dispersion_matrix.col(0)) / (2.0);
      transform -= (expectation % dispersion_matrix.col(0) / 2.0) % arma::sqrt((4.0 / dispersion_matrix.col(0)) % (expectation % random_normal) + arma::square(expectation) % arma::square(random_normal));
    }
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        observations(i) = copula_values(i) * (expectation(i) + transform(i)) <= expectation(i) ? transform(i) : expectation(i) * expectation(i) / transform(i);
    }
    return observations;
}

const arma::vec InverseGauss::random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const
{
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    arma::vec random_normal(expectation.n_elem, arma::fill::randn);
    random_normal = arma::square(random_normal);
    arma::vec transform = expectation;
    transform += (arma::square(expectation) % random_normal % dispersion_) / 2.0;
    transform -= ((expectation % dispersion_) / 2.0) % arma::sqrt((4.0 / dispersion_)% (expectation % random_normal) + arma::square(expectation) % arma::square(random_normal));
    arma::vec accept(expectation.n_elem, arma::fill::randu);
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        observations(i) = accept(i) * (expectation(i) + transform(i)) <= expectation(i) ? transform(i) : expectation(i) * expectation(i) / transform(i);
    }
    return observations;
}

const arma::vec InverseGauss::random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const
{
    arma::vec observations(expectation.n_elem, arma::fill::zeros);
    arma::vec copula_values(expectation.n_elem);
    copula_values = Rcpp::as<arma::vec>(this->copula_sample(1, this->copula_object));

    arma::vec random_normal(expectation.n_elem, arma::fill::randn);
    random_normal = arma::square(random_normal);
    arma::vec transform = expectation;
    transform += (arma::square(expectation) % random_normal % dispersion_) / 2.0;
    transform -= ((expectation % dispersion_) / 2.0) % arma::sqrt((4.0 / dispersion_) % (expectation % random_normal) + arma::square(expectation) % arma::square(random_normal));
    for(unsigned int i = 0; i < expectation.n_elem; i++)
    {
        observations(i) = copula_values(i) * (expectation(i) + transform(i)) <= expectation(i) ? transform(i) : expectation(i) * expectation(i) / transform(i);
    }
    return observations;
}

/*
  Log-likelihood for Inverse Gaussian family
*/

const double InverseGauss::log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const 
{
    double value = - (observation - expectation) * (observation - expectation) / (2.0 * expectation * expectation * observation * dispersion_);
    value -= 0.5 * std::log(2 * arma::datum::pi * dispersion_ * observation * observation * observation);
    return value;
}

/*
    Deviance and Pearson residuals for Inverse Gaussian family
*/

const double InverseGauss::deviance_residual(const double &observation, const double &expectation) const
{
    return std::pow(observation - expectation, 2) / (std::pow(expectation, 2) * observation);  
}

const double InverseGauss::pearson_residual(const double &observation, const double &expectation) const
{
    return std::pow(observation - expectation, 2) / std::pow(expectation, 3);
}


/*
  Additonal methods for Inverse Gaussian family
*/

const bool InverseGauss::valid_expectation(const arma::mat &x) const 
{
  return x.is_finite() && arma::all(arma::vectorise(x) > 0.0);
}

const arma::mat InverseGauss::variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const
{
    arma::mat response = this->inverse_link(link_values);
    return (response % response % response) * dispersion_;
}

/*
  Define various link functions for Inverse Gaussian family
*/


/*
    Define model with natural link, i.e., link function: g(mu) = 1/mu^2
    Observation transformation: h(y) = 1/y^2
    Link transformation: h(psi) = psi
*/

class NaturalInverseGauss : public InverseGauss {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    NaturalInverseGauss(const Rcpp::RObject &dispersion, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : InverseGauss(dispersion, true, true, copula_obj, "1/mu^2"){};
    NaturalInverseGauss(const Rcpp::RObject &dispersion) : InverseGauss(dispersion, true, true, "1/mu^2"){};
    ~NaturalInverseGauss() = default;
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
          return new NaturalInverseGauss(disp_vec, this->copula_object);
        } else {
          return new NaturalInverseGauss(disp_vec);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new NaturalInverseGauss(disp_mat, this->copula_object);
        } else {
          return new NaturalInverseGauss(disp_mat);
        }
      }
    };
};


const double NaturalInverseGauss::inverse_link(const double x) const
{
  return 1.0 / std::sqrt(x);
}

const double NaturalInverseGauss::link(const double x) const
{
  return 1.0 / (x * x);
}

const double NaturalInverseGauss::derivative_inverse_link(const double x) const
{
  return -0.5 / (x * std::sqrt(x));
}


const double NaturalInverseGauss::observation_trafo(const double x) const
{
  return 1.0 / (x * x);
}

const double NaturalInverseGauss::link_trafo(const double x) const
{
  return x;
}

const double NaturalInverseGauss::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool NaturalInverseGauss::valid_link(const arma::mat &x) const
{
  return x.is_finite() && arma::all(arma::vectorise(x) > 0.0);
}

/*
    Define inverse link, i.e., link function: g(mu) = 1/mu
    Observation transformation: h(y) = 1/y
    Link transformation: h(psi) = psi
*/

class InverseInverseGauss : public InverseGauss {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    InverseInverseGauss(const Rcpp::RObject &dispersion, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : InverseGauss(dispersion, true, true, copula_obj, "inverse"){};
    InverseInverseGauss(const Rcpp::RObject &dispersion) : InverseGauss(dispersion, true, true, "inverse"){};
    ~InverseInverseGauss() = default;
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
          return new InverseInverseGauss(disp_vec, this->copula_object);
        } else {
          return new InverseInverseGauss(disp_vec);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new InverseInverseGauss(disp_mat, this->copula_object);
        } else {
          return new InverseInverseGauss(disp_mat);
        }
      }
    };
};



const double InverseInverseGauss::inverse_link(const double x) const
{
  return 1.0 / x;
}

const double InverseInverseGauss::link(const double x) const
{
  return 1.0 / x;
}

const double InverseInverseGauss::derivative_inverse_link(const double x) const
{
  return -1.0 / (x * x);
}


const double InverseInverseGauss::observation_trafo(const double x) const
{
  return 1.0 / x;
}

const double InverseInverseGauss::link_trafo(const double x) const
{
  return x;
}

const double InverseInverseGauss::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool InverseInverseGauss::valid_link(const arma::mat &x) const
{
  return x.is_finite() && arma::all(arma::vectorise(x) > 0.0);
}

/*
    Define identity link, i.e., link function: g(mu) = mu
    Observation transformation: h(y) = y
    Link transformation: h(psi) = psi
*/

class LinearInverseGauss : public InverseGauss {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    LinearInverseGauss(const Rcpp::RObject &dispersion, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : InverseGauss(dispersion, true, true, copula_obj, "identity"){};
    LinearInverseGauss(const Rcpp::RObject &dispersion) : InverseGauss(dispersion, true, true, "identity"){};
    ~LinearInverseGauss() = default;
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
          return new LinearInverseGauss(disp_vec, this->copula_object);
        } else {
          return new LinearInverseGauss(disp_vec);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new LinearInverseGauss(disp_mat, this->copula_object);
        } else {
          return new LinearInverseGauss(disp_mat);
        }
      }
    };
};

const double LinearInverseGauss::inverse_link(const double x) const
{
  return x;
}

const double LinearInverseGauss::link(const double x) const
{
  return x;
}

const double LinearInverseGauss::derivative_inverse_link(const double x) const
{
  return 1.0;
}


const double LinearInverseGauss::observation_trafo(const double x) const
{
  return x;
}

const double LinearInverseGauss::link_trafo(const double x) const
{
  return x;
}

const double LinearInverseGauss::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool LinearInverseGauss::valid_link(const arma::mat &x) const
{
  return x.is_finite() && arma::all(arma::vectorise(x) > 0.0);
}


/*
    Define log link, i.e., link function: g(mu) = log(mu)
    Observation transformation: h(y) = log(y)
    Link transformation: h(psi) = psi
*/

class LogInverseGauss : public InverseGauss {
  protected:
    virtual const double inverse_link(const double x) const;
    virtual const double link(const double x) const;
    virtual const double derivative_inverse_link(const double x) const;
    virtual const double observation_trafo(const double x) const;
    virtual const double link_trafo(const double x) const;
    virtual const double derivative_link_trafo(const double x) const;
  public:
    LogInverseGauss(const Rcpp::RObject &dispersion, const Rcpp::Nullable<Rcpp::S4> &copula_obj) : InverseGauss(dispersion, true, false, copula_obj, "log"){};
    LogInverseGauss(const Rcpp::RObject &dispersion) : InverseGauss(dispersion, true, false, "log"){};
    ~LogInverseGauss() = default;
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
          return new LogInverseGauss(disp_vec, this->copula_object);
        } else {
          return new LogInverseGauss(disp_vec);
        }
      } else {
        disp_mat = Rcpp::wrap(this->dispersion_matrix);
        if(use_dependence)
        { 
          return new LogInverseGauss(disp_mat, this->copula_object);
        } else {
          return new LogInverseGauss(disp_mat);
        }
      }
    };
};

const double LogInverseGauss::inverse_link(const double x) const
{
  return exp(x);
}

const double LogInverseGauss::link(const double x) const
{
  return std::log(x);
}

const double LogInverseGauss::derivative_inverse_link(const double x) const
{
  return exp(x);
}


const double LogInverseGauss::observation_trafo(const double x) const
{
  return log(x);
}

const double LogInverseGauss::link_trafo(const double x) const
{
  return x;
}

const double LogInverseGauss::derivative_link_trafo(const double x) const
{
  return 1.0;
}

const bool LogInverseGauss::valid_link(const arma::mat &x) const
{
  return true;
}




/*
    InverseGauss factory methods
*/

InverseGauss* InverseGauss::create(const Rcpp::List &family, Rcpp::S4& copula_obj)
{
  std::string link = Rcpp::as<std::string>(family["link"]);
  Rcpp::RObject dispersion = family["dispersion"];
  if(link == "log")
  {
    return new LogInverseGauss(dispersion, copula_obj);
  } 
  else if(link == "identity") 
  {
    return new LinearInverseGauss(dispersion, copula_obj);
  } 
  else if(link == "inverse") 
  {
    return new InverseInverseGauss(dispersion, copula_obj);
  }
  else if(link == "1/mu^2") 
  {
    return new NaturalInverseGauss(dispersion, copula_obj);
  } 
  else
  {
    throw std::invalid_argument("The desired link is currently not yet implemented.");
  }
}


InverseGauss* InverseGauss::create(const Rcpp::List &family)
{
  std::string link = Rcpp::as<std::string>(family["link"]);
  Rcpp::RObject dispersion = family["dispersion"];
  if(link == "log")
  {
    return new LogInverseGauss(dispersion);
  } 
  else if(link == "identity") 
  {
    return new LinearInverseGauss(dispersion);
  } 
  else if(link == "inverse") 
  {
    return new InverseInverseGauss(dispersion);
  }
  else if(link == "1/mu^2") 
  {
    return new NaturalInverseGauss(dispersion);
  } 
  else
  {
    throw std::invalid_argument("The desired link is currently not yet implemented.");
  }
}






