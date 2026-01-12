/* 
-----------------------------------------------------------------------------
    File: family.h
    Purpose: Declaration of Family classes
    Author: Steffen Maletz
    Last modified: 2025-12-06
-----------------------------------------------------------------------------
*/

/* 
 Parts of the following algorithm are translated and adapted from
 the R package 'RNGforGDP', licensed under the GNU General Public License.
 The translated code is used in the QuasiPoisson family for random number generation,
 See lines 307-358 for the relevant functions.
 Original source: Function GenUniGpois in RNGforGDP package
 Original authors: Hesen Li, Ruizhe Chen, Hai Nguyen, Yu-Che Chung, Ran Gao, Hakan Demirtas
*/

#ifndef FAMILY_H
#define FAMILY_H


/*
    Family: Abstract base class for different types of distributions (families) used in the model
    Defines the vectorized interface for family-specific functions and factory method for creating Family objects
    Abstract methods:
      * random_observation_independent: Random number generation assuming contemporaneous independence
      * random_observation: Random number generation considering contemporaneous dependence via copula 
      * inverse_link: Inverse of the link function
      * link: Link function
      * derivative_link_trafo: Derivative of the function applied to transform past link values
      * derivative_inverse_link: Derivative of the inverse of the link function
      * observation_trafo: Transformation applied to observations
      * link_trafo: Transformation applied to past link values
      * log_likelihood: (Quasi-)Log-likelihood of a single observation
      * deviance_residual: Deviance residual for a single observation
      * pearson_residual: Pearson residual for a single observation
    Contains common attributes:
      * dispersion: Dispersion parameter (if constant)
      * dispersion_matrix: Dispersion parameters for each time point and spatial location (if not constant)
      * const_dispersion: Boolean indicating whether dispersion is constant
      * identity_link_trafo: Boolean indicating whether the transformation applied to past link values is the identity function. Used to optimize calculations.
      * only_positive_parameters: Boolean indicating whether the family only allows (non-negative) parameters.
      * family_in_R: Name of the family as used in R (e.g. "poisson", "gaussian", etc.)
      * link_name: Name of the link function as used in R (e.g. "log", "identity", etc.)
*/  
class Family
{
protected:
    Rcpp::S4 copula_object;
    bool use_dependence = false;
    const Rcpp::Environment pkg = Rcpp::Environment::namespace_env("copula");
    Rcpp::Function copula_sample = Rcpp::Environment("package:stats")["rpois"]; // placeholder, will be overwritten in constructor
    virtual const arma::vec random_observation_independent(const arma::vec &expectation) const = 0;
    virtual const arma::vec random_observation(const arma::vec &expectation) const = 0;
    virtual const arma::vec random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const = 0;
    virtual const arma::vec random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const = 0;
    virtual const double inverse_link(const double x) const = 0;
    virtual const double link(const double x) const = 0;
    virtual const double derivative_link_trafo(const double x) const = 0;
    virtual const double derivative_inverse_link(const double x) const = 0;
    virtual const double observation_trafo(const double x) const = 0;
    virtual const double link_trafo(const double x) const = 0;
    virtual const double log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const = 0;
    virtual const double deviance_residual(const double &observation, const double &expectation) const = 0;
    virtual const double pearson_residual(const double &observation, const double &expectation) const = 0;
public:
    double dispersion;
    arma::mat dispersion_matrix;
    bool const_dispersion;

    const bool identity_link_trafo;
    const bool only_positive_parameters;
    const std::string family_in_R; // Name of Family in R
    const std::string link_name; // Link-Name for R

    Family(Rcpp::RObject dispersion_, const bool identity_link_trafo, const bool &only_positive, const Rcpp::Nullable<Rcpp::S4> &copula_obj, const char* family_name, const char* link_name_r)  : identity_link_trafo(identity_link_trafo), only_positive_parameters(only_positive), family_in_R(family_name), link_name(link_name_r)
    {
        copula_sample = pkg["rCopula"];
        if(copula_obj.isNotNull())
        {
            this->use_dependence = true;
            this->copula_object = Rcpp::S4( copula_obj );
        }

        // If dispersion_ is a numeric vector of length 1, then constant, otherwise matrix
        if(Rf_isMatrix(dispersion_))
        {
            const_dispersion = false;
            Rcpp::NumericMatrix disp_mat(dispersion_);
            dispersion_matrix = Rcpp::as<arma::mat>(disp_mat);
        } else
        {
            Rcpp::NumericVector disp_vec(dispersion_);
            const_dispersion = disp_vec.size() == 1;
            if(const_dispersion)
            {
              dispersion = disp_vec[0];
            } else {
              dispersion_matrix = arma::mat(disp_vec.size(), 1);
              dispersion_matrix.col(0) = Rcpp::as<arma::vec>(disp_vec);
            }
        }
    };
    Family(Rcpp::RObject dispersion_, const bool identity_link_trafo, const bool &only_positive, const char* family_name, const char* link_name_r) : identity_link_trafo(identity_link_trafo), only_positive_parameters(only_positive), family_in_R(family_name), link_name(link_name_r)
    {
        copula_sample = pkg["rCopula"];

        if(Rf_isMatrix(dispersion_))
        {
            const_dispersion = false;
            Rcpp::NumericMatrix disp_mat(dispersion_);
            dispersion_matrix = Rcpp::as<arma::mat>(disp_mat);
        } else
        {
            Rcpp::NumericVector disp_vec(dispersion_);
            const_dispersion = disp_vec.size() == 1;
            if(const_dispersion)
            {
              dispersion = disp_vec[0];
            } else {
              dispersion_matrix = arma::mat(disp_vec.size(), 1);
              dispersion_matrix.col(0) = Rcpp::as<arma::vec>(disp_vec);
            }
        }
    };
    virtual ~Family() = default;
    virtual Family* clone() const = 0;
    
    static Family * create(const Rcpp::List &family, Rcpp::S4 &copula_obj);
    static Family * create(const Rcpp::List &family);

    const arma::mat inverse_link(const arma::mat &x) const
    {
        arma::mat result(x.n_rows, x.n_cols);
        for(unsigned int i = 0; i < x.n_rows; i++)
        {
            for(unsigned int j = 0; j < x.n_cols; j++)
            {
                result(i, j) = this->inverse_link( (double) x(i, j));
            }
        }
        return result;
    };
    const arma::mat link(const arma::mat &x) const
    {
        arma::mat result(x.n_rows, x.n_cols);
        for(unsigned int i = 0; i < x.n_rows; i++)
        {
            for(unsigned int j = 0; j < x.n_cols; j++)
            {
                result(i, j) = this->link( (double) x(i, j));
            }
        }
        return result;
    };
    const arma::mat derivative_inverse_link(const arma::mat &x) const
    {
        arma::mat result(x.n_rows, x.n_cols);
        for(unsigned int i = 0; i < x.n_rows; i++)
        {
            for(unsigned int j = 0; j < x.n_cols; j++)
            {
                result(i, j) = this->derivative_inverse_link( (double) x(i, j));
            }
        }
        return result;
    };
    const arma::mat observation_trafo(const arma::mat &x) const
    {
      arma::mat result(x.n_rows, x.n_cols);
        for(unsigned int i = 0; i < x.n_rows; i++)
        {
            for(unsigned int j = 0; j < x.n_cols; j++)
            {
                result(i, j) = this->observation_trafo( (double) x(i, j));
            }
        }
        return result;
    };
    const arma::mat link_trafo(const arma::mat &x) const
    {  
      arma::mat result(x.n_rows, x.n_cols);
        for(unsigned int i = 0; i < x.n_rows; i++)
        {
            for(unsigned int j = 0; j < x.n_cols; j++)
            {
                result(i, j) = this->link_trafo( (double) x(i, j));
            }
        }
        return result;
    };
    const arma::mat derivative_link_trafo(const arma::mat &x) const
    {
        arma::mat result(x.n_rows, x.n_cols);
        for(unsigned int i = 0; i < x.n_rows; i++)
        {
            for(unsigned int j = 0; j < x.n_cols; j++)
            {
                result(i, j) = this->derivative_link_trafo( (double) x(i, j));
            }
        }
        return result;
    };
    const arma::mat log_likelihood(const arma::mat &observations, const arma::mat &expectations, const arma::mat &dispersion_) const
    {
        arma::mat result(observations.n_rows, observations.n_cols);
        for(unsigned int i = 0; i < observations.n_rows; i++)
        {
            for(unsigned int j = 0; j < observations.n_cols; j++)
            {
                result(i, j) = this->log_likelihood( (double) observations(i, j), (double) expectations(i, j), (double) dispersion_(i, j)); 
            }
        }
        return result;
    };

    const arma::mat deviance_residual(const arma::mat &observations, const arma::mat &expectations) const
    {
      arma::mat result(observations.n_rows, observations.n_cols);
        for(unsigned int i = 0; i < observations.n_rows; i++)
        {
            for(unsigned int j = 0; j < observations.n_cols; j++)
            {
                result(i, j) = this->deviance_residual( (double) observations(i, j), (double) expectations(i, j)); 
            }
        }
        return result;
    };
    const arma::mat pearson_residual(const arma::mat &observations, const arma::mat &expectations) const
    {
      arma::mat result(observations.n_rows, observations.n_cols);
        for(unsigned int i = 0; i < observations.n_rows; i++)
        {
            for(unsigned int j = 0; j < observations.n_cols; j++)
            {
                result(i, j) = this->pearson_residual( (double) observations(i, j), (double) expectations(i, j)); 
            }
        }
        return result;
    };
    virtual const arma::mat variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const = 0;

    const arma::vec sample(const arma::vec &expectation) const
    {
        if(this->use_dependence){
            return random_observation(expectation);
        } else {
            return random_observation_independent(expectation);
        }
    };
    const arma::vec sample(const arma::vec &expectation, const arma::vec &dispersion_) const
    {
        if(this->use_dependence){
            return random_observation(expectation, dispersion_);
        } else {
            return random_observation_independent(expectation, dispersion_);
        }
    };
    virtual const bool valid_expectation(const arma::mat &x) const = 0;
    virtual const bool valid_link(const arma::mat &x) const = 0;
};


/*
    Poisson distribution family
    Additional attribute:
      * fast_sampling: Boolean indicating whether to use fast sampling method for random number generation
        If true, uses inversion method is used, otherwise the data is generated using the Poisson process
*/ 
class Poisson : public Family 
{
  protected:
    const bool use_fast_sampling;
    virtual const arma::vec random_observation_independent(const arma::vec &expectation) const;
    virtual const arma::vec random_observation(const arma::vec &expectation) const;
    virtual const arma::vec random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const;
    virtual const arma::vec random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const;
    virtual const arma::vec random_observation_fast(const arma::vec &expectation) const;
    virtual const arma::vec random_observation_slow(const arma::vec &expectation) const;
    virtual const double log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const;
    virtual const double deviance_residual(const double &observation, const double &expectation) const;
    virtual const double pearson_residual(const double &observation, const double &expectation) const;
  public:
    Poisson(const bool &fast_sampling, const bool &identity_link_trafo, const bool &only_positive, const Rcpp::Nullable<Rcpp::S4> &copula_obj, const char* link_name_r) : Family(Rcpp::RObject(Rcpp::NumericVector::create(1.0)), identity_link_trafo, only_positive, copula_obj, "poisson", link_name_r), use_fast_sampling(fast_sampling) {};
    Poisson(const bool &fast_sampling, const bool &identity_link_trafo, const bool &only_positive, const char* link_name_r) : Family(Rcpp::RObject(Rcpp::NumericVector::create(1.0)), identity_link_trafo, only_positive, "poisson", link_name_r), use_fast_sampling(fast_sampling) {};
    static Poisson* create(const Rcpp::List &family, Rcpp::S4& copula_obj);
    static Poisson* create(const Rcpp::List &family);
    virtual ~Poisson() = default;
    virtual const bool valid_expectation(const arma::mat &x) const;
    virtual const arma::mat variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const;
};

/*
    Family for Quasi-Poisson model (over- or underdispersed Poisson)
    Additional attribute:
      * sampling_method: Method used for random number generation ("branching", "negbin", "chop_down", "build_up")
        - "branching": Branching process method (only for overdispersion/no contemporaneous dependence possible)
        - "negbin": Use reparameterized Negative binomial distribution (only for overdispersion)
        - "chop_down": Generalized Poisson distribution using chop-down algorithm
        - "build_up": Generalized Poisson distribution using build-up algorithm
*/
class QuasiPoisson : public Family 
{
  private:
    int (QuasiPoisson::*random_generator)(const double&, const double&, const double&) const;
    // Translated code from RNGforGDP package with suitable transformation of dispersion parameter
    int rgpoisson_branching(const double &mu, const double &phi, const double &u) const {
      // Only for phi >= 1
      double theta = mu / std::sqrt(phi);
      double lambda = 1.0 - 1.0 / std::sqrt(phi);
      int total = R::rpois(theta);
      int current = total;
      while((current > 0) && (phi >= 1.0))
      {
        current = R::rpois(lambda * current);
        total += current;
      }
      return total;
    };
    // Translated code from RNGforGDP package with suitable transformation of dispersion parameter
    int rgpoisson_chop_down(const double &mu, const double &phi, const double &u) const {
      double theta = mu / std::sqrt(phi);
      double lambda = 1.0 - 1.0 / std::sqrt(phi);
      int total = 0;
      double px = std::exp(-theta);  // P(X=0)
      double u_temp = u;
      while (u_temp > px) {
        u_temp -= px;
        total += 1;
        double log_px = std::log(theta) +
          (total - 1) * std::log(theta + lambda * total) -
          (theta + lambda * total) -
          std::lgamma(total + 1);
        px = std::exp(log_px);
      }
      return total;
    };
    // Translated code from RNGforGDP package with suitable transformation of dispersion parameter
    int rgpoisson_build_up(const double &mu, const double &phi, const double &u) const {
      double theta = mu / std::sqrt(phi);
      double lambda = 1.0 - 1.0 / std::sqrt(phi);
      int total = 0;
      double t = std::exp(-theta);
      double px = t;
      double s = px;
      double log_px;
      while (u > s) {
        total += 1;
        log_px = std::log(theta) + 
          (total - 1) * std::log(theta + lambda * total) -
          (theta + lambda * total) -
          std::lgamma(total + 1);
        px = std::exp(log_px);
        s += px;
      }
      return total;
    };
    int rquasi_negbin(const double &mu, const double &phi, const double &u) const {
      // Only for phi >= 1
      int sample;
      if(phi > 1.0)
      {
        double size = mu / (phi - 1.0);  // r
        double prob = size / (size + mu);
        sample = R::qnbinom(u, size, prob, true, false);
      } else {
        sample = R::qpois(u, mu, true, false);
      }
      return sample;
    };
  protected:
    const std::string sampling_method;
    virtual const arma::vec random_observation_independent(const arma::vec &expectation) const;
    virtual const arma::vec random_observation(const arma::vec &expectation) const;
    virtual const arma::vec random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const;
    virtual const arma::vec random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const;
    virtual const double log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const;
    virtual const double deviance_residual(const double &observation, const double &expectation) const;
    virtual const double pearson_residual(const double &observation, const double &expectation) const;
  public:
    QuasiPoisson(const std::string& sampling_method, const Rcpp::RObject &dispersion, const bool &identity_link_trafo, const bool &only_positive, const Rcpp::Nullable<Rcpp::S4> &copula_obj, const char* link_name_r) : Family(dispersion, identity_link_trafo, only_positive, copula_obj, "quasipoisson", link_name_r), sampling_method(sampling_method){
      if (sampling_method == "branching") {
        random_generator = &QuasiPoisson::rgpoisson_branching;
      } else if (sampling_method == "negbin") {
        random_generator = &QuasiPoisson::rquasi_negbin;
      } else if (sampling_method == "chop_down") {
        random_generator = &QuasiPoisson::rgpoisson_chop_down;
      } else if (sampling_method == "build_up") {
        random_generator = &QuasiPoisson::rgpoisson_build_up;
      }
    };
    QuasiPoisson(const std::string& sampling_method, const Rcpp::RObject &dispersion, const bool &identity_link_trafo, const bool &only_positive, const char* link_name_r) : Family(dispersion, identity_link_trafo, only_positive, "quasipoisson", link_name_r), sampling_method(sampling_method) {
      if (sampling_method == "branching") {
        random_generator = &QuasiPoisson::rgpoisson_branching;
      } else if (sampling_method == "negbin") {
        random_generator = &QuasiPoisson::rquasi_negbin;
      } else if (sampling_method == "chop_down") {
        random_generator = &QuasiPoisson::rgpoisson_chop_down;
      } else if (sampling_method == "build_up") {
        random_generator = &QuasiPoisson::rgpoisson_build_up;
      }
    };
    static QuasiPoisson* create(const Rcpp::List &family, Rcpp::S4& copula_obj);
    static QuasiPoisson* create(const Rcpp::List &family);
    virtual ~QuasiPoisson() = default;
    virtual const bool valid_expectation(const arma::mat &x) const;
    virtual const arma::mat variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const;
};


/*
    Normal (Gaussian) distribution family
*/
class Normal : public Family 
{
  protected:
    virtual const arma::vec random_observation_independent(const arma::vec &expectation) const;
    virtual const arma::vec random_observation(const arma::vec &expectation) const;
    virtual const arma::vec random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const;
    virtual const arma::vec random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const;
    virtual const double log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const;
    virtual const double deviance_residual(const double &observation, const double &expectation) const;
    virtual const double pearson_residual(const double &observation, const double &expectation) const;
  public:
    Normal(const Rcpp::RObject &dispersion, const bool &identity_link_trafo, const bool &only_positive, const Rcpp::Nullable<Rcpp::S4> &copula_obj, const char* link_name_r) : Family(dispersion, identity_link_trafo, only_positive, copula_obj, "gaussian", link_name_r) {};
    Normal(const Rcpp::RObject &dispersion, const bool &identity_link_trafo, const bool &only_positive, const char* link_name_r) : Family(dispersion, identity_link_trafo, only_positive, "gaussian", link_name_r) {};
    static Normal* create(const Rcpp::List &family, Rcpp::S4& copula_obj);
    static Normal* create(const Rcpp::List &family);
    virtual ~Normal() = default;
    virtual const bool valid_expectation(const arma::mat &x) const;
    virtual const arma::mat variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const;
};


// Negative Binomial
class NegativeBinomial : public Family 
{
  protected:
    // const double dispersion;
    virtual const arma::vec random_observation_independent(const arma::vec &expectation) const;
    virtual const arma::vec random_observation(const arma::vec &expectation) const;
    virtual const arma::vec random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const;
    virtual const arma::vec random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const;
    virtual const double log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const;
    virtual const double deviance_residual(const double &observation, const double &expectation) const;
    virtual const double pearson_residual(const double &observation, const double &expectation) const;
  public:
    NegativeBinomial(const Rcpp::RObject &dispersion, const bool &identity_link_trafo, const bool &only_positive, const Rcpp::Nullable<Rcpp::S4> &copula_obj, const char* link_name_r) : Family(dispersion, identity_link_trafo, only_positive, copula_obj, "negative_binomial", link_name_r) {};
    NegativeBinomial(const Rcpp::RObject &dispersion, const bool &identity_link_trafo, const bool &only_positive, const char* link_name_r) : Family(dispersion, identity_link_trafo, only_positive, "negative_binomial", link_name_r) {};
    static NegativeBinomial* create(const Rcpp::List &family, Rcpp::S4& copula_obj);
    static NegativeBinomial* create(const Rcpp::List &family);
    virtual ~NegativeBinomial() = default;
    virtual const bool valid_expectation(const arma::mat &x) const;
    virtual const arma::mat variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const;
};


/*
    Inverse Gaussian distribution family
*/
class InverseGauss : public Family 
{
  protected:
    // const double dispersion;
    virtual const arma::vec random_observation_independent(const arma::vec &expectation) const;
    virtual const arma::vec random_observation(const arma::vec &expectation) const;
    virtual const arma::vec random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const;
    virtual const arma::vec random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const;
    virtual const double log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const;
    virtual const double deviance_residual(const double &observation, const double &expectation) const;
    virtual const double pearson_residual(const double &observation, const double &expectation) const;
  public:
    InverseGauss(const Rcpp::RObject &dispersion, const bool &identity_link_trafo, const bool &only_positive, const Rcpp::Nullable<Rcpp::S4> &copula_obj, const char* link_name_r) : Family(dispersion, identity_link_trafo, only_positive, copula_obj, "inverse_gaussian", link_name_r) {};
    InverseGauss(const Rcpp::RObject &dispersion, const bool &identity_link_trafo, const bool &only_positive, const char* link_name_r) : Family(dispersion, identity_link_trafo, only_positive, "inverse_gaussian", link_name_r) {};
    static InverseGauss* create(const Rcpp::List &family, Rcpp::S4& copula_obj);
    static InverseGauss* create(const Rcpp::List &family);
    virtual ~InverseGauss() = default;
    virtual const bool valid_expectation(const arma::mat &x) const;
    virtual const arma::mat variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const;
};


/*    
    Gamma distribution family
*/
class Gamma : public Family 
{
  protected:
    double shape;
    arma::mat shape_matrix;
    virtual const arma::vec random_observation_independent(const arma::vec &expectation) const;
    virtual const arma::vec random_observation(const arma::vec &expectation) const;
    virtual const arma::vec random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const;
    virtual const arma::vec random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const;
    virtual const double log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const;
    virtual const double deviance_residual(const double &observation, const double &expectation) const;
    virtual const double pearson_residual(const double &observation, const double &expectation) const;
  public:
    Gamma(const Rcpp::RObject &dispersion, const bool &identity_link_trafo, const bool &only_positive, const Rcpp::Nullable<Rcpp::S4> &copula_obj, const char* link_name_r) : Family(dispersion, identity_link_trafo, only_positive, copula_obj, "gamma", link_name_r){
      if(this->const_dispersion)
      {
        shape = 1.0 / this->dispersion;
      } else {
        shape_matrix = 1.0 / this->dispersion_matrix;
      }
    };
    Gamma(const Rcpp::RObject &dispersion, const bool &identity_link_trafo, const bool &only_positive, const char* link_name_r) : Family(dispersion, identity_link_trafo, only_positive, "gamma", link_name_r) {
      if(this->const_dispersion)
      {
        shape = 1.0 / this->dispersion;
      } else {
        shape_matrix = 1.0 / this->dispersion_matrix;
      }
    };
    static Gamma* create(const Rcpp::List &family, Rcpp::S4& copula_obj);
    static Gamma* create(const Rcpp::List &family);
    virtual ~Gamma() = default;
    virtual const bool valid_expectation(const arma::mat &x) const;
    virtual const arma::mat variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const;
};

/*
    Binomial distribution family
    Additional attribute:
      * n_upper: Number of trials (can be vector for varying number of trials)
      * index_n_upper: Index for n_upper, in case n_upper contains multiple values (one vor each location)
*/
class Binomial : public Family 
{
  protected:
    mutable unsigned int index_n_upper = 0; // Index fuer n_upper, falls n_upper mehrere Werte enthaelt
    const arma::uvec n_upper;
    virtual const arma::vec random_observation_independent(const arma::vec &expectation) const;
    virtual const arma::vec random_observation(const arma::vec &expectation) const;
    virtual const arma::vec random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const;
    virtual const arma::vec random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const;
    virtual const double log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const;
    virtual const double deviance_residual(const double &observation, const double &expectation) const;
    virtual const double pearson_residual(const double &observation, const double &expectation) const;
  public:
    Binomial(const arma::uvec &n, const bool &identity_link_trafo, const bool &only_positive, const Rcpp::Nullable<Rcpp::S4> &copula_obj, const char* link_name_r) : Family(Rcpp::NumericVector::create(1.0), identity_link_trafo, only_positive, copula_obj, "binomial", link_name_r), n_upper(n) {};
    Binomial(const arma::uvec &n, const bool &identity_link_trafo, const bool &only_positive, const char* link_name_r) : Family(Rcpp::NumericVector::create(1.0), identity_link_trafo, only_positive, "binomial", link_name_r), n_upper(n) {};
    static Binomial* create(const Rcpp::List &family, Rcpp::S4& copula_obj);
    static Binomial* create(const Rcpp::List &family);
    virtual ~Binomial() = default;
    virtual const bool valid_expectation(const arma::mat &x) const;
    virtual const arma::mat variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const;
};


/*
    Family for Quasi-Binomial model (over- or underdispersed Binomial)
    Additional attribute:
      * n_upper: Number of trials (can be vector for varying number of trials)
      * index_n_upper: Index for n_upper, in case n_upper contains multiple values (one vor each location)
    Additional method:
      * rquasi: Random number generator for Quasi-Binomial distribution (contemporaneous independence)
      * rquasi_dependent: Random number generator for Quasi-Binomial distribution based on uniforms (generated with contemporaneous dependence by copula)
        - For underdispersion a normal approximation is used
        - For overdispersion a correlated Bernoulli sequence is generated based and the number of successes is counted. 
          The correlation is calculated based on the dispersion parameter.
    Note: The dispersion parameter is limited to the interval (0, n]. It is not possible to have dispersion values larger than n (number of trials).
*/

class QuasiBinomial : public Family 
{
  private:
    const int rquasi(const unsigned int &n, const double &p, const double &dispersion) const {
      int number;
      if(dispersion >= 1.0){
      double rho = std::sqrt((dispersion - 1.0) / (n - 1));
      arma::uvec bernoulli_sequence = (arma::randu<arma::vec>(n + 1) < p);
      arma::uvec mixing = (arma::randu<arma::vec>(n) < rho);
      arma::uvec result = mixing * bernoulli_sequence(0) + (1 - mixing) % bernoulli_sequence.subvec(1, n);
      number = static_cast<int>(arma::accu(result));
    } else {
      double mu = n * p;
      double sd = std::sqrt(n * p * (1.0 - p) * dispersion);
      arma::vec value(1, arma::fill::randn);
      value *= sd;
      value += mu;
      value.clamp(0.0, n);
      value = arma::round(value);
      number = static_cast<int>(value(0));
    }
    return number;
  };
  const int rquasi_dependent(const unsigned int &n, const double &p, const double &dispersion, const arma::vec &uniforms1, const arma::vec &uniforms2) const {
      int number;
      if(dispersion >= 1.0){
      double rho = std::sqrt((dispersion - 1.0) / (n - 1));
      arma::uvec bernoulli_sequence = (uniforms1 < p);
      arma::uvec mixing = (uniforms2 < rho);
      arma::uvec result = mixing * bernoulli_sequence(0) + (1 - mixing) % bernoulli_sequence.subvec(1, n);
      number = static_cast<int>(arma::accu(result));
    } else {
      double mu = n * p;
      double sd = std::sqrt(n * p * (1.0 - p) * dispersion);
      arma::vec value(1);
      value(0) = R::qnorm(uniforms1(0), mu, sd, true, false);
      value.clamp(0.0, n);
      value = arma::round(value);
      number = static_cast<int>(value(0));
    }
    return number;
  };


  protected:
    mutable unsigned int index_n_upper = 0;
    const arma::uvec n_upper;
    virtual const arma::vec random_observation_independent(const arma::vec &expectation) const;
    virtual const arma::vec random_observation(const arma::vec &expectation) const;
    virtual const arma::vec random_observation_independent(const arma::vec &expectation, const arma::vec &dispersion_) const;
    virtual const arma::vec random_observation(const arma::vec &expectation, const arma::vec &dispersion_) const;
    virtual const double log_likelihood(const double &observation, const double &expectation, const double &dispersion_) const;
    virtual const double deviance_residual(const double &observation, const double &expectation) const;
    virtual const double pearson_residual(const double &observation, const double &expectation) const;
  public:
    QuasiBinomial(const arma::uvec &n, const Rcpp::RObject &dispersion, const bool &identity_link_trafo, const bool &only_positive, const Rcpp::Nullable<Rcpp::S4> &copula_obj, const char* link_name_r) : Family(dispersion, identity_link_trafo, only_positive, copula_obj, "binomial", link_name_r), n_upper(n) {};
    QuasiBinomial(const arma::uvec &n, const Rcpp::RObject &dispersion, const bool &identity_link_trafo, const bool &only_positive, const char* link_name_r) : Family(dispersion, identity_link_trafo, only_positive, "binomial", link_name_r), n_upper(n) {};
    static QuasiBinomial* create(const Rcpp::List &family, Rcpp::S4& copula_obj);
    static QuasiBinomial* create(const Rcpp::List &family);
    virtual ~QuasiBinomial() = default;
    virtual const bool valid_expectation(const arma::mat &x) const;
    virtual const arma::mat variance_fun(const arma::mat &link_values, const arma::mat &dispersion_) const;
};


#endif
