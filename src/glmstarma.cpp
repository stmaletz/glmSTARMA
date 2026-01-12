/* 
-----------------------------------------------------------------------------
    File: glmstarma.cpp
    Purpose: Implementation of the main glmstarma fitting function called from R
    Author: Steffen Maletz
    Last modified: 2025-12-06
-----------------------------------------------------------------------------
*/


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include "glmstarma.h"


/*
    Function called in R to fit a glmstarma model, i.e., to estimate the parameters of the mean model.
    Input are the preprocessed data and model specifications from the `glmstarma` R function.
    Output is a List containing estimation results
    It has to be further processed in R.
*/


// [[Rcpp::export]]
Rcpp::List glmstarma_cpp(const arma::mat &ts, const Rcpp::List &covariate_list, 
    const Rcpp::List &model, const Rcpp::List &wlist_ar, 
    const Rcpp::Nullable<Rcpp::List> wlist_ma, const Rcpp::Nullable<Rcpp::List> wlist_covariates, 
    const Rcpp::List &family, const Rcpp::List &control)
{
    // Get Control Arguments
    const bool use_sparsity = Rcpp::as<bool>(control["use_sparsity"]);
    const double sparsity_threshold = Rcpp::as<double>(control["sparsity_threshold"]);
    const std::string dispersion_estimation_type = Rcpp::as<std::string>(control["dispersion_est_type"]);
    const bool estimate_dispersion = Rcpp::as<bool>(family["estimate_dispersion"]);

    Rcpp::RObject init_method_param = control["parameter_init"];
    Rcpp::RObject init_method_link = control["init_link"];

    // Set up the model:
    const unsigned int dim = ts.n_rows;

    CovariateList covariates(covariate_list, ts.n_cols, ts.n_rows, 0, 0);
    Orders orders(model, ts.n_rows, ts.n_cols);
    Family * fam = Family::create(family);

    Neighborhood * W_ar = Neighborhood::create(wlist_ar, sparsity_threshold, use_sparsity);
    Neighborhood * W_ma = Neighborhood::create_neighbor(orders.n_param_ma, wlist_ma, W_ar, dim, sparsity_threshold, use_sparsity);
    Neighborhood * W_covariates = Neighborhood::create_neighbor(orders.n_param_ma, wlist_covariates, W_ar, dim, sparsity_threshold, use_sparsity);

    // Initialize model
    arma::mat link_init = Model::init_link(init_method_link, orders, ts, fam);
    arma::vec start_param = Model::init_param(init_method_param, orders, ts, fam);

    arma::vec intercept(orders.n_param_intercept);
    arma::mat ar_params(orders.autoregressive_orders.n_rows, orders.autoregressive_orders.n_cols);
    arma::mat ma_params(orders.moving_average_orders.n_rows, orders.moving_average_orders.n_cols);
    arma::mat cov_params(orders.covariate_orders.n_rows, orders.covariate_orders.n_cols);
    Parameter::create_param_matrices(start_param, &orders, intercept, ar_params, ma_params, cov_params);

    if(orders.n_param_ma > 0)
    {
        W_ma->set_parameter_matrices(ma_params);
    }

    Design model_design(ts, &covariates, link_init, &orders, fam, W_ar, W_ma, W_covariates);
    FittingObject * fitting = FittingObject::create(ts, &model_design, &orders, &covariates, W_ma, fam, control);

    // Fit the model:
    arma::vec estimation = fitting->fit(start_param);
    Rcpp::List fitting_info = fitting->get_fitting_info();
    Rcpp::List algorithm_info = fitting->get_algorithm_info();
    delete fitting;
    fitting = nullptr;

    Parameter::create_param_matrices(estimation, &orders, intercept, ar_params, ma_params, cov_params);
    if(orders.n_param_ma > 0)
    {
        W_ma->set_parameter_matrices(ma_params);
    }
    
    // Prepare output:
    arma::vec fitted = model_design.update_design(estimation, &orders, fam, W_ma);
    double dispersion_est = 1.0;
    arma::vec ts_vec = arma::vectorise(ts.tail_cols(orders.n_obs_effective));
    if(estimate_dispersion)
    {
        fitted = fam->inverse_link(fitted);
        arma::vec residuals = dispersion_estimation_type == "deviance" ? fam->deviance_residual(ts_vec, fitted) : fam->pearson_residual(ts_vec, fitted);
        if(fam->family_in_R == "negative_binomial"){
            residuals -= 1;
            residuals /= fitted;
        }
        residuals /= (orders.n_obs_effective * dim - orders.n_param);
        dispersion_est = arma::accu(residuals);
        fam->dispersion = dispersion_est;
        fam->const_dispersion = true;
    }

    Parameter param(estimation, orders);
    Rcpp::List param_est = param.return_param_list();
    arma::vec score(estimation.n_elem);
    arma::mat information(estimation.n_elem, estimation.n_elem);
    
    double log_likelihood = Model::log_likelihood(estimation, score, information, ts_vec, &model_design, &orders, W_ma, fam, true, true);
    log_likelihood /= orders.n_obs_effective;
    log_likelihood *= ts.n_cols;
    arma::mat variance_estimation = Model::variance_estimation(estimation, ts_vec, &model_design, fam, information, W_ma, &orders);
    arma::mat var_temp = information * variance_estimation / orders.n_obs_effective;

    double aic = -2.0 * log_likelihood + 2.0 * orders.n_param;
    double bic = -2.0 * log_likelihood + std::log( static_cast<double>(orders.n_obs_effective * dim) ) * orders.n_param;
    double qic = -2.0 * log_likelihood + 2.0 * arma::trace(var_temp);
    arma::mat fitted_values = fam->inverse_link(model_design.link_vals);

    Rcpp::List results;
    results["covariates"] = covariate_list;
    results["model"] = model;
    results["family"] = family;
    if(estimate_dispersion)
    {
        results["dispersion"] = dispersion_est;
    }
    results["control"] = control;
    results["coefficients_list"] = param_est;
    results["coefficients"] = estimation;
    results["target_dim"] = dim;
    results["n_obs_effective"] = orders.n_obs_effective;
    results["max_time_lag"] = orders.max_time_lag;
    results["log_likelihood"] = log_likelihood;
    results["score"] = score;
    results["information"] = information;
    results["variance_estimation"] = variance_estimation;
    results["aic"] = aic;
    results["bic"] = bic;
    results["qic"] = qic;
    results["design_matrix"] = model_design.design_matrix;
    results["derivatives"] = model_design.derivative_link;
    results["fitted.values"] = fitted_values;
    results["link_values"] = model_design.link_vals;
    results["algorithm"] = algorithm_info;
    results["convergence"] = fitting_info;

    delete fam;
    delete W_ar;
    delete W_ma;
    delete W_covariates;

    return results;
}

