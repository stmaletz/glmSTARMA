/* 
-----------------------------------------------------------------------------
    File: simulation.cpp
    Purpose: Implementation of simulation functions for glmstarma and dglmstarma models
    Author: Steffen Maletz
    Last modified: 2026-01-12
-----------------------------------------------------------------------------
*/

#include "glmstarma.h"


/*
    Internal function to simulate data from a glmstarma process
    Used in the R wrapper function glmstarma_sim_cpp
    Inputs:
    * observations: Matrix to store the simulated observations (dim x (n_obs + burn_in))
    * link_values: Matrix to store the link values (dim x (n_obs + burn_in))
    * initial_link: Initial link values to start the simulation (dim x n_start)
    * fam: Family object
    * intercept: Intercept parameters
    * model_order: Model order information
    * wlist_ar: Neighborhood object for autoregressive parameters
    * wlist_ma: Neighborhood object for moving average parameters
    * covariates: Covariate list object
    * wlist_cov: Neighborhood object for covariates
    * burn_in: Number of burn-in observations
    Outputs:
    * None (observations and link_values matrices are modified in place)
*/


void glmstarma_sim_internal(arma::mat &observations, arma::mat &link_values, const arma::mat &initial_link, const Family * fam, 
                            const arma::vec &intercept, const Orders &model_order, Neighborhood* wlist_ar, Neighborhood* wlist_ma, 
                            const CovariateList &covariates, Neighborhood* wlist_cov, const unsigned int &burn_in)
{
    const unsigned int dim = observations.n_rows;
    arma::vec sample(dim);
    arma::vec temp_link(dim);
    arma::vec temp(dim);
    arma::vec temp_expectation(dim);
    arma::vec inter = (intercept.n_elem == 1) ? arma::vec(dim, arma::fill::value(intercept(0))) : intercept;
    arma::mat initial_mean = fam->inverse_link( initial_link );
    // Simulate initial observations
    for(unsigned int t = 0; t < initial_link.n_cols; t++)
    {
        sample = fam->sample(initial_mean.col(t));
        observations.col(t) = sample;
        link_values.col(t) = initial_link.col(t);
    }
    // Simulate rest of burn-in period ignoring covariate values
    for(unsigned int t = initial_link.n_cols; t < burn_in; t++)
    {
        Model::calculate_link_value_at(t, link_values, observations, model_order, fam, inter, &covariates, wlist_ar, wlist_ma, wlist_cov, true);
        temp_link = link_values.col(t);
        temp_expectation = fam->inverse_link(temp_link);
        sample = fam->sample(temp_expectation);
        observations.col(t) = sample;
    }
    // Simulate remaining observations
    for(unsigned int t = burn_in; t < observations.n_cols; t++)
    {
        Model::calculate_link_value_at(t, link_values, observations, model_order, fam, inter, &covariates, wlist_ar, wlist_ma, wlist_cov, false);
        temp_link = link_values.col(t);
        temp_expectation = fam->inverse_link(temp_link);
        sample = fam->sample(temp_expectation);
        observations.col(t) = sample;
    }
}


/*
    Exported Rcpp function to simulate data from a glmstarma process
    Inputs:
    * n_obs_in_time: Number of observation times to simulate
    * parameters: Parameters of the glmstarma process
    * model: Model specification of the glmstarma process
    * wlist_ar: List of neighborhood matrices for autoregressive part
    * wlist_ma: Optional list of neighborhood matrices for moving average part
    * family: Family object for the glmstarma process
    * covariates: List of covariates used to simulate the data
    * wlist_covariates: Optional list of neighborhood matrices for covariates
    * copula_obj: Optional copula object for contemporaneous dependencies
    * n_start: Number of burn-in observation times
    * control: List of control parameters for the simulation
    Outputs:
    * A list with the simulated observations, link values, model and parameter specifications
*/



// [[Rcpp::export]]
Rcpp::List glmstarma_sim_cpp(const Rcpp::IntegerVector& n_obs_in_time, const Rcpp::List& parameters, const Rcpp::List& model, const Rcpp::List& wlist_ar, 
                             const Rcpp::Nullable<Rcpp::List> wlist_ma, const Rcpp::List& family, const Rcpp::List& covariates, const Rcpp::Nullable<Rcpp::List> wlist_covariates,
                             const Rcpp::Nullable<Rcpp::S4> copula_obj, const Rcpp::IntegerVector& n_start, const Rcpp::List &control)
{
    // Get Control Arguments
    const bool use_sparsity = Rcpp::as<bool>(control["use_sparsity"]);
    const bool return_burn_in = Rcpp::as<bool>(control["return_burn_in"]);
    const double sparsity_threshold = Rcpp::as<double>(control["sparsity_threshold"]);
    const int number_of_observations_in_time = n_obs_in_time[0];
    const int burn_in = n_start[0];

    const unsigned int dim = Neighborhood::get_matrix_dim(wlist_ar[0]); 
    Family * fam = nullptr;
    if(copula_obj.isNotNull())
    {
        Rcpp::S4 copula_object(copula_obj);
        fam = Family::create(family, copula_object);
    } else {
        fam = Family::create(family);
    }
    arma::mat ts(dim, number_of_observations_in_time + burn_in);
    arma::mat link_values(dim, number_of_observations_in_time + burn_in);
    CovariateList covariate_list(covariates, number_of_observations_in_time, dim, burn_in, 0);
    Orders model_order(model, dim, number_of_observations_in_time);

    Neighborhood* W_ar = Neighborhood::create(wlist_ar, sparsity_threshold, use_sparsity);
    Neighborhood * W_ma = Neighborhood::create_neighbor(model_order.n_param_ma, wlist_ma, W_ar, dim, sparsity_threshold, use_sparsity);
    Neighborhood * W_covariates = Neighborhood::create_neighbor(model_order.n_param_ma, wlist_covariates, W_ar, dim, sparsity_threshold, use_sparsity);

    // Set parameters for simulation
    Parameter param(parameters, &model_order);
    param.set_param_vector(&model_order);
    W_ar->set_parameter_matrices(param.ar_params);
    W_ma->set_parameter_matrices(param.ma_params);
    W_covariates->set_parameter_matrices(param.cov_params);

    // Calculate unconditional (mean) of link values as initial mean values
    arma::vec denominator(1, arma::fill::zeros);
    denominator(0) = param.ma_params.is_empty() ? 0.0 : arma::accu(param.ma_params);
    denominator(0) += param.ar_params.is_empty() ? 0.0 : arma::accu(param.ar_params);
    denominator = 1.0 - denominator;
    arma::vec u_link = arma::mean(param.intercept) / denominator;
    double unconditional_link = u_link(0);
    arma::mat initial_links(dim, model_order.max_time_lag, arma::fill::value(unconditional_link));
    
    // Simulate data
    glmstarma_sim_internal(ts, link_values, initial_links, fam, param.intercept, model_order, W_ar, W_ma, covariate_list, W_covariates, burn_in);

    // Return values
    if(!return_burn_in){
        ts = ts.tail_cols(number_of_observations_in_time);
        link_values = link_values.tail_cols(number_of_observations_in_time);
    }
    Rcpp::List result = Rcpp::List::create(Rcpp::Named("observations") = ts, Rcpp::Named("link_values") = link_values, Rcpp::Named("model") = model, Rcpp::Named("parameters") = parameters);
    delete W_ar;
    delete W_ma;
    delete W_covariates;
    delete fam;
    return result;
}



/*
    Internal function to simulate data from a dglmstarma process
    Used in the R wrapper function dglmstarma_sim_cpp
    Inputs:
    * observations: Matrix to store the simulated observations (dim x (n_obs + burn_in))
    * link_values: Matrix to store the link values for the mean (dim x (n_obs + burn_in))
    * pseudo_observations: Matrix to store the pseudo-observations (dim x (n_obs + burn_in))
    * dispersion_values: Matrix to store the dispersion values (dim x (n_obs + burn_in))
    * initial_link: Initial link values to start the simulation (dim x n_start)
    * mean_fam: Family object for the mean
    * disp_fam: Family object for the dispersion
    * pseudo_obs_type: Type of pseudo-observations ("deviance" or "pearson")
    * intercept_mean: Intercept parameters for the mean
    * intercept_dispersion: Intercept parameters for the dispersion
    * orders_mean: Model order information for the mean
    * orders_dispersion: Model order information for the dispersion
    * wlist_ar: Neighborhood object for autoregressive parameters of the mean
    * wlist_ma: Neighborhood object for moving average parameters of the mean
    * wlist_cov: Neighborhood object for covariates of the mean
    * wlist_ar_dispersion: Neighborhood object for autoregressive parameters of the dispersion
    * wlist_ma_dispersion: Neighborhood object for moving average parameters of the dispersion
    * wlist_covariates_dispersion: Neighborhood object for covariates of the dispersion
    * covariates_mean: Covariate list object for the mean
    * covariates_dispersion: Covariate list object for the dispersion
    * burn_in: Number of burn-in observations
    Outputs:
    * None (observations, link_values, pseudo_observations and dispersion_values matrices are modified in place)
*/


void dglmstarma_sim_internal(arma::mat &observations, arma::mat &link_values, arma::mat &pseudo_observations, arma::mat &dispersion_values,
                             const arma::mat &initial_link, const Family * mean_fam, const Family * disp_fam, const std::string &pseudo_obs_type,
                             const arma::vec &intercept_mean, const arma::vec &intercept_dispersion, 
                             const Orders &orders_mean, const Orders &orders_dispersion, 
                             Neighborhood* wlist_ar, Neighborhood* wlist_ma, Neighborhood* wlist_cov,
                             Neighborhood* wlist_ar_dispersion, Neighborhood* wlist_ma_dispersion, Neighborhood* wlist_covariates_dispersion,
                             const CovariateList &covariates_mean, const CovariateList &covariates_dispersion, 
                             const unsigned int &burn_in)
{
    const unsigned int dim = observations.n_rows;
    const unsigned int max_time_lag = initial_link.n_cols;

    arma::vec sample(dim);
    arma::vec temp_link(dim);
    arma::vec temp_expectation(dim);
    arma::vec temp_dispersion(dim);

    arma::mat initial_mean = mean_fam->inverse_link( initial_link ); 
    dispersion_values.head_cols(max_time_lag).fill(1.0);
    arma::mat dispersion_link_values(dim, observations.n_cols);
    dispersion_link_values.head_cols(max_time_lag) = disp_fam->link(dispersion_values.head_cols(max_time_lag));

    // create intercept vectors
    arma::vec inter_mean = (intercept_mean.n_elem == 1) ? arma::vec(dim, arma::fill::value(intercept_mean(0))) : intercept_mean;
    arma::vec inter_dispersion = (intercept_dispersion.n_elem == 1) ? arma::vec(dim, arma::fill::value(intercept_dispersion(0))) : intercept_dispersion;
    
    // First values without model
    for(unsigned int t = 0; t < max_time_lag; t++)
    {
        temp_expectation = initial_mean.col(t);
        temp_dispersion = dispersion_values.col(t);
        sample = mean_fam->sample(temp_expectation, temp_dispersion);
        observations.col(t) = sample;
        link_values.col(t) = initial_link.col(t);
        pseudo_observations.col(t) = (pseudo_obs_type == "deviance") ? mean_fam->deviance_residual(sample, temp_expectation) : mean_fam->pearson_residual(sample, temp_expectation);
    }
    // Rest of the burn-in period, ignoring covariate effects
    for(unsigned int t = max_time_lag; t < burn_in; t++)
    {
        Model::calculate_link_value_at(t, link_values, observations, orders_mean, mean_fam, inter_mean, &covariates_mean, wlist_ar, wlist_ma, wlist_cov, true);
        temp_expectation = mean_fam->inverse_link(link_values.col(t));

        Model::calculate_link_value_at(t, dispersion_link_values, pseudo_observations, orders_dispersion, disp_fam, inter_dispersion, &covariates_dispersion, wlist_ar_dispersion, wlist_ma_dispersion, wlist_covariates_dispersion, true);
        temp_dispersion = disp_fam->inverse_link(dispersion_link_values.col(t));
        dispersion_values.col(t) = temp_dispersion;
        
        if(mean_fam->family_in_R == "negative_binomial")
        {
            temp_dispersion -= 1.0;
            temp_dispersion /= temp_expectation;
            temp_dispersion.clamp(0.0, arma::datum::inf);
        }

        sample = mean_fam->sample(temp_expectation, temp_dispersion);
        observations.col(t) = sample;
        pseudo_observations.col(t) = (pseudo_obs_type == "deviance") ? mean_fam->deviance_residual(sample, temp_expectation) : mean_fam->pearson_residual(sample, temp_expectation);
    }
    // Remaining observations
    for(unsigned int t = burn_in; t < observations.n_cols; t++)
    {
        Model::calculate_link_value_at(t, link_values, observations, orders_mean, mean_fam, inter_mean, &covariates_mean, wlist_ar, wlist_ma, wlist_cov, false);
        temp_expectation = mean_fam->inverse_link(link_values.col(t));

        Model::calculate_link_value_at(t, dispersion_link_values, pseudo_observations, orders_dispersion, disp_fam, inter_dispersion, &covariates_dispersion, wlist_ar_dispersion, wlist_ma_dispersion, wlist_covariates_dispersion, false);
        temp_dispersion = disp_fam->inverse_link(dispersion_link_values.col(t));
        dispersion_values.col(t) = temp_dispersion;

        if(mean_fam->family_in_R == "negative_binomial")
        {
            temp_dispersion -= 1.0;
            temp_dispersion /= temp_expectation;
            temp_dispersion.clamp(0.0, arma::datum::inf);
        }
        
        sample = mean_fam->sample(temp_expectation, temp_dispersion);
        observations.col(t) = sample;
        pseudo_observations.col(t) = (pseudo_obs_type == "deviance") ? mean_fam->deviance_residual(sample, temp_expectation) : mean_fam->pearson_residual(sample, temp_expectation);
    }
}
/*
    Exported Rcpp function to simulate data from a dglmstarma process
    Inputs:
    * n_obs_in_time: Number of observation times to simulate
    * parameters_mean: Parameters of the mean model
    * parameters_dispersion: Parameters of the dispersion model
    * mean_family: Family object for the mean
    * dispersion_family: Family object for the dispersion
    * mean_model: Model specification of the mean
    * dispersion_model: Model specification of the dispersion
    * pseudo_observations_type: Type of pseudo-observations ("deviance" or "pearson")
    * wlist: List of neighborhood matrices for autoregressive part of the mean
    * mean_covariates: List of covariates for the mean
    * dispersion_covariates: List of covariates for the dispersion
    * wlist_past_mean: Optional list of neighborhood matrices for moving average part of the mean
    * wlist_covariates: Optional list of neighborhood matrices for covariates of the mean
    * wlist_ar_dispersion: Optional list of neighborhood matrices for autoregressive part of the dispersion
    * wlist_past_mean_dispersion: Optional list of neighborhood matrices for moving average part of the dispersion
    * wlist_covariates_dispersion: Optional list of neighborhood matrices for covariates of the dispersion
    * copula_obj_mean: Optional copula object for contemporaneous dependencies in the mean
    * n_start: Number of burn-in observation times
    * control: List of control parameters for the simulation
    Outputs:
    * A list with the simulated observations, link values, pseudo-observations, dispersion values, model and parameter specifications
*/

// [[Rcpp::export]]
Rcpp::List dglmstarma_sim_cpp(const unsigned int& n_obs_in_time, const Rcpp::List& parameters_mean, const Rcpp::List &parameters_dispersion,
                              const Rcpp::List &mean_family, const Rcpp::List &dispersion_family,
                              const Rcpp::List& mean_model, const Rcpp::List& dispersion_model, const std::string &pseudo_observations_type, const Rcpp::List& wlist, 
                              const Rcpp::List &mean_covariates, const Rcpp::List &dispersion_covariates, const Rcpp::Nullable<Rcpp::List> wlist_past_mean,
                              const Rcpp::Nullable<Rcpp::List> wlist_covariates, const Rcpp::Nullable<Rcpp::List> wlist_ar_dispersion,
                              const Rcpp::Nullable<Rcpp::List> wlist_past_mean_dispersion, const Rcpp::Nullable<Rcpp::List> wlist_covariates_dispersion,
                              const Rcpp::List &control, const Rcpp::Nullable<Rcpp::S4> copula_obj_mean,
                              const unsigned int& n_start)
{
    const bool use_sparsity = Rcpp::as<bool>(control["use_sparsity"]);
    const bool return_burn_in = Rcpp::as<bool>(control["return_burn_in"]);
    const double sparsity_threshold = Rcpp::as<double>(control["sparsity_threshold"]);

    Family * mean_fam = nullptr;
    Family * disp_fam = nullptr;
    if(copula_obj_mean.isNotNull())
    {
        Rcpp::S4 copula_object(copula_obj_mean);
        mean_fam = Family::create(mean_family, copula_object);
    } else {
        mean_fam = Family::create(mean_family);
    }
    disp_fam = Family::create(dispersion_family);

    const unsigned int dim = Neighborhood::get_matrix_dim(wlist[0]);

    arma::mat ts(dim, n_obs_in_time + n_start);
    arma::mat link_values(dim, n_obs_in_time + n_start);
    arma::mat pseudo_observations(dim, n_obs_in_time + n_start);
    arma::mat dispersion_values(dim, n_obs_in_time + n_start);

    CovariateList covariates_mean(mean_covariates, n_obs_in_time, dim, n_start, 0);
    Orders orders_mean(mean_model, dim, n_obs_in_time);
    CovariateList covariates_dispersion(dispersion_covariates, n_obs_in_time, dim, n_start, 0);
    Orders orders_dispersion(dispersion_model, dim, n_obs_in_time);

    // Check if dispersion model is constant in time
    const bool const_dispersion_model = (orders_dispersion.n_param_ar == 0) && (orders_dispersion.n_param_ma == 0) && (orders_dispersion.n_param_cov == 0);

    Neighborhood * W_ar_mean = Neighborhood::create(wlist, sparsity_threshold, use_sparsity);
    Neighborhood * W_ma_mean = Neighborhood::create_neighbor(orders_mean.n_param_ma, wlist_past_mean, W_ar_mean, ts.n_rows, sparsity_threshold, use_sparsity);
    Neighborhood * W_covariates_mean = Neighborhood::create_neighbor(orders_mean.n_param_cov, wlist_covariates, W_ar_mean, ts.n_rows, sparsity_threshold, use_sparsity);
    Neighborhood * W_ar_dispersion = nullptr;
    Neighborhood * W_ma_dispersion = nullptr;
    Neighborhood * W_covariates_dispersion = nullptr;
    if(!const_dispersion_model)
    {
        W_ar_dispersion = Neighborhood::create_neighbor(orders_dispersion.n_param_ar, wlist_ar_dispersion, W_ar_mean, ts.n_rows, sparsity_threshold, use_sparsity);
        W_ma_dispersion = Neighborhood::create_neighbor(orders_dispersion.n_param_ma, wlist_past_mean_dispersion, W_ar_mean, ts.n_rows, sparsity_threshold, use_sparsity);
        W_covariates_dispersion = Neighborhood::create_neighbor(orders_dispersion.n_param_cov, wlist_covariates_dispersion, W_ar_mean, ts.n_rows, sparsity_threshold, use_sparsity);
    }

    // Set parameters for simulation
    Parameter param_mean(parameters_mean, orders_mean);
    Parameter param_dispersion(parameters_dispersion, orders_dispersion);
    param_mean.set_param_vector(orders_mean);
    W_ar_mean->set_parameter_matrices(param_mean.ar_params);
    W_ma_mean->set_parameter_matrices(param_mean.ma_params);
    W_covariates_mean->set_parameter_matrices(param_mean.cov_params);
    param_dispersion.set_param_vector(orders_dispersion);
    W_ar_dispersion->set_parameter_matrices(param_dispersion.ar_params);
    W_ma_dispersion->set_parameter_matrices(param_dispersion.ma_params);
    W_covariates_dispersion->set_parameter_matrices(param_dispersion.cov_params);

    // Calculate unconditional (mean) of link values as initial mean values
    arma::vec denominator(1, arma::fill::zeros);
    denominator(0) = param_mean.ma_params.is_empty() ? 0.0 : arma::accu(param_mean.ma_params);
    denominator(0) += param_mean.ar_params.is_empty() ? 0.0 : arma::accu(param_mean.ar_params);
    denominator = 1.0 - denominator;
    arma::vec u_link = arma::mean(param_mean.intercept) / denominator;
    double unconditional_link = u_link(0);
    const unsigned int max_lag = std::max(orders_mean.max_time_lag, orders_dispersion.max_time_lag);
    arma::mat initial_links(dim, max_lag, arma::fill::value(unconditional_link));

    // Simulate data
    dglmstarma_sim_internal(ts, link_values, pseudo_observations, dispersion_values, initial_links, mean_fam, disp_fam, pseudo_observations_type,
                             param_mean.intercept, param_dispersion.intercept, orders_mean, orders_dispersion,
                             W_ar_mean, W_ma_mean, W_covariates_mean,
                             W_ar_dispersion, W_ma_dispersion, W_covariates_dispersion,
                             covariates_mean, covariates_dispersion, n_start);

    // Return values
    if(!return_burn_in)
    {
        ts = ts.tail_cols(n_obs_in_time);
        link_values = link_values.tail_cols(n_obs_in_time);
        pseudo_observations = pseudo_observations.tail_cols(n_obs_in_time);
        dispersion_values = dispersion_values.tail_cols(n_obs_in_time);
    }
    Rcpp::List result = Rcpp::List::create(Rcpp::Named("observations") = ts, 
                                            Rcpp::Named("link_values") = link_values, 
                                            Rcpp::Named("pseudo_observations") = pseudo_observations,
                                            Rcpp::Named("dispersion_values") = dispersion_values,
                                            Rcpp::Named("mean_model") = mean_model, 
                                            Rcpp::Named("dispersion_model") = dispersion_model,
                                            Rcpp::Named("parameters_mean") = parameters_mean, 
                                            Rcpp::Named("parameters_dispersion") = parameters_dispersion);

    delete W_ar_mean;
    delete W_ma_mean;
    delete W_covariates_mean;
    delete W_ar_dispersion;
    delete W_ma_dispersion;
    delete W_covariates_dispersion;
    delete mean_fam;
    delete disp_fam;
    return result;
}
