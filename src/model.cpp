/* 
-----------------------------------------------------------------------------
    File: model.cpp
    Purpose: Implementation of Model class
    Author: Steffen Maletz
    Last modified: 2025-12-13
-----------------------------------------------------------------------------
*/

#include "glmstarma.h"


/*
    Helper function to initialize the link values
    Possible methods:
    * Matrix with link values, supplied by user
    * "first_obs" -> Use first observations as link values
    * "mean" -> Use mean of transformed observations as link values
    * "transformed_mean" -> Use transformed mean of observations as link values
    Output:
    * Matrix with initial link values
*/

arma::mat Model::init_link(const Rcpp::RObject &method, const Orders &orders, const arma::mat &ts, const Family* fam)
{
    arma::mat linkvals(ts.n_rows, orders.max_time_lag);
    if(orders.n_param_ma > 0){
        if(Rf_isMatrix(method)){
            Rcpp::NumericMatrix init_values = Rcpp::as<Rcpp::NumericMatrix>( method );
            return Rcpp::as<arma::mat>( init_values );
        }
        Rcpp::CharacterVector init_method = Rcpp::as<Rcpp::CharacterVector>( method );
        arma::mat transformed_obs = fam->observation_trafo(ts);
        if(init_method[0] == "first_obs"){
            linkvals = transformed_obs.head_cols(orders.max_time_lag);
        } else if(init_method[0] == "mean" || init_method[0] == "transformed_mean"){
            arma::vec row_mean = init_method[0] == "mean" ? arma::mean(transformed_obs, 1) : fam->observation_trafo(arma::mean(ts, 1));
            for(unsigned int i = 0; i < orders.max_time_lag; i++)
            {
                linkvals.col(i) = row_mean;
            }
        }
    }
    return linkvals;
}


/*
    Helper function to initialize the parameters for parameter estimation
    Possible methods:
    * List with parameter matrices, supplied by user
    * Vector with parameters, supplied by user
    * "zero" -> All parameters set to zero (or small positive values if only positive parameters are allowed)
    * "random" -> All parameters set to random values (within bounds if only positive parameters are allowed)
    Output:
    * Vector with initial parameters
*/

arma::vec Model::init_param(const Rcpp::RObject &method, const Orders &orders, const arma::mat &ts, const Family* fam)
{
    arma::vec intercept;
    arma::mat autoregressive_params;
    arma::mat moving_average_params;
    arma::mat covariate_params;

    if(Rcpp::is<Rcpp::List>(method)){
        Rcpp::List initial_parameters = Rcpp::as<Rcpp::List>( method );
        Rcpp::NumericVector intercept_vector = initial_parameters["intercept"];
        intercept = Rcpp::as<arma::vec>(intercept_vector);

        if(orders.n_param_ar > 0){
            Rcpp::NumericMatrix temp_ar = initial_parameters["past_obs"];
            autoregressive_params = Rcpp::as<arma::mat>(temp_ar);
        }

        if(orders.n_param_ma > 0){
            Rcpp::NumericMatrix temp_ma = initial_parameters["past_mean"];
            moving_average_params = Rcpp::as<arma::mat>(temp_ma);
        }

        if(orders.n_param_cov > 0){
            Rcpp::NumericMatrix temp_cov = initial_parameters["covariates"];
            covariate_params = Rcpp::as<arma::mat>(temp_cov);
        }
    }  else if(Rcpp::is<Rcpp::NumericVector>(method)){
        Rcpp::NumericVector initial_parameters = Rcpp::as<Rcpp::NumericVector>( method );
        arma::vec params = Rcpp::as<arma::vec>(initial_parameters);
        return params;
    } else {
        Rcpp::CharacterVector init_method = Rcpp::as<Rcpp::CharacterVector>( method );
        // Create param-matrices:
        intercept = arma::vec(orders.n_param_intercept);
        if(orders.n_param_ar > 0){
            autoregressive_params = arma::mat(orders.autoregressive_orders.n_rows, orders.autoregressive_orders.n_cols);
        }
        if(orders.n_param_ma > 0){
            moving_average_params = arma::mat(orders.moving_average_orders.n_rows, orders.moving_average_orders.n_cols);
        }
        if(orders.n_param_cov > 0){
            covariate_params = arma::mat(orders.covariate_orders.n_rows, orders.covariate_orders.n_cols);
        }

        double param_sum = 0.0;
        if(init_method[0] == "zero"){
            if(fam->only_positive_parameters){
                if(orders.n_param_ar > 0){
                    autoregressive_params.elem(orders.which_autoregressive_orders).fill( 1.0 / (orders.n_param_ar + orders.n_param_ma + 10.0) );
                    param_sum += arma::accu(autoregressive_params);
                }
                if(orders.n_param_ma > 0){
                    moving_average_params.elem(orders.which_moving_average_orders).fill( 1.0 / (orders.n_param_ar + orders.n_param_ma + 10.0) );
                    param_sum += arma::accu(moving_average_params);
                }
                if(orders.n_param_cov > 0){
                    covariate_params.elem(orders.which_covariate_orders).fill( 1.0 );
                }
            }
        } else if(init_method[0] == "random"){
            if(fam->only_positive_parameters){
                if(orders.n_param_ar > 0){
                    autoregressive_params.elem(orders.which_autoregressive_orders) = arma::randu(orders.n_param_ar, arma::distr_param(0.0, 1.0 / (orders.n_param + orders.n_param_ma)));
                    param_sum += arma::accu(autoregressive_params);
                }
                if(orders.n_param_ma > 0){
                    moving_average_params.elem(orders.which_moving_average_orders) = arma::randu(orders.n_param_ma, arma::distr_param(0.0, 1.0 / (orders.n_param + orders.n_param_ma)));
                    param_sum += arma::accu(moving_average_params);
                }
                if(orders.n_param_cov > 0){
                    covariate_params.elem(orders.which_covariate_orders).randu();
                }
            } else {
                if(orders.n_param_ar > 0){
                    autoregressive_params.elem(orders.which_autoregressive_orders) = arma::randu(orders.n_param_ar, arma::distr_param(-1.0 / (orders.n_param + orders.n_param_ma), 1.0 / (orders.n_param + orders.n_param_ma)));
                    param_sum += arma::accu(autoregressive_params);
                }
                if(orders.n_param_ma > 0){
                    moving_average_params.elem(orders.which_moving_average_orders) = arma::randu(orders.n_param_ma, arma::distr_param(-1.0 / (orders.n_param + orders.n_param_ma), 1.0 / (orders.n_param + orders.n_param_ma)));
                    param_sum += arma::accu(moving_average_params);
                }
                if(orders.n_param_cov > 0){
                    covariate_params.elem(orders.which_covariate_orders) = arma::randu(orders.n_param_cov, arma::distr_param(-1.0, 1.0 ));
                }
            }
        }
        arma::mat transformed_obs = fam->observation_trafo( ts );
        if(orders.intercepts_equal){
            intercept(0) = (1.0 - param_sum) * arma::mean(arma::mean(transformed_obs));
        } else {
            intercept = (1.0 - param_sum) * arma::mean(transformed_obs, 1);
        }
    }

    arma::vec parameter(orders.n_param);
    if(orders.n_param_ma > 0){
        parameter.head(orders.n_param_ma) = moving_average_params.elem( orders.which_moving_average_orders );
    }
    parameter.subvec(orders.n_param_ma, orders.n_param_ma + orders.n_param_intercept - 1) = intercept;
    if(orders.n_param_ar > 0){
        parameter.subvec(orders.n_param_ma + orders.n_param_intercept, orders.n_param_ma + orders.n_param_intercept + orders.n_param_ar - 1) = autoregressive_params.elem( orders.which_autoregressive_orders );
    }
    if(orders.n_param_cov > 0){
        parameter.tail( orders.n_param_cov ) = covariate_params.elem( orders.which_covariate_orders );
    }
    return parameter;
}


/*
    Helper function to calculate the log-likelihood of the model, as well as the score vector and information matrix if necessary
    Output:
    * Log-likelihood value
    * Score vector (if calc_score = true)
    * Information matrix (if calc_information = true)
*/

double Model::log_likelihood(const arma::vec &x, arma::vec &score, arma::mat &information, const arma::vec &ts_vec, Design * design, const Orders * orders, Neighborhood * W_ma, const Family * family, bool calc_score, bool calc_information)
{
    arma::vec intercept(orders->n_param_intercept);
    arma::mat ar_params(orders->autoregressive_orders.n_rows, orders->autoregressive_orders.n_cols);
    arma::mat ma_params(orders->moving_average_orders.n_rows, orders->moving_average_orders.n_cols);
    arma::mat cov_params(orders->covariate_orders.n_rows, orders->covariate_orders.n_cols); 
    Parameter::create_param_matrices(x, *orders, intercept, ar_params, ma_params, cov_params);
    if(orders->n_param_ma > 0)
    {
        W_ma->set_parameter_matrices(ma_params);
    }
    arma::vec link_values = design->update_design(x, orders, family, W_ma);
    arma::vec fitted_values = family->inverse_link( link_values ) ; //
    arma::vec dispersion_(link_values.n_elem);
    if(family->const_dispersion)
    {
        dispersion_.fill(family->dispersion);
    } else {
        if(family->dispersion_matrix.n_cols > 1)
        {
            dispersion_ = arma::vectorise(family->dispersion_matrix);
        } else {
            dispersion_ = arma::repmat(family->dispersion_matrix.col(0), link_values.n_elem / family->dispersion_matrix.n_elem, 1);
        }
    }

    double log_like = arma::accu( family->log_likelihood( ts_vec, fitted_values, dispersion_ ) );

    if(calc_score || calc_information)
    {
        design->update_derivative(orders, W_ma, family);
        arma::vec residuals = ts_vec - fitted_values;
        arma::vec middle = family->derivative_inverse_link( link_values ) / family->variance_fun( link_values, dispersion_ );
        residuals = middle % residuals;
    
        if(calc_score)
        {
            score = (residuals.t() * design->derivative_link).t();
        }
        if(calc_information)
        {
            middle = middle % family->derivative_inverse_link( link_values );
            information = design->derivative_link.t() * (design->derivative_link.each_col() % middle);
        }
    }    
    return log_like;
};


/*
    Helper function to estimate the variance of the parameter estimates using a sandwich estimator
    Output:
    * Variance-covariance matrix of the parameter estimates
*/


arma::mat Model::variance_estimation(const arma::vec &x, const arma::vec &ts_vec, Design * design, const Family * family, const arma::mat &information, Neighborhood * W_ma, const Orders * orders)
{
    arma::vec link_values = design->design_matrix * x;
    arma::vec fitted_values = family->inverse_link( link_values ) ;
    arma::vec residuals = ts_vec - fitted_values;
    arma::vec dispersion_(link_values.n_elem);
    if(family->const_dispersion)
    {
        dispersion_.fill(family->dispersion);
    } else {
        if(family->dispersion_matrix.n_cols > 1)
        {
            dispersion_ = arma::vectorise(family->dispersion_matrix);
        } else {
            dispersion_ = arma::repmat(family->dispersion_matrix.col(0), link_values.n_elem / family->dispersion_matrix.n_elem, 1);
        }
    }
    residuals = family->derivative_inverse_link( link_values ) % residuals;
    residuals = residuals / family->variance_fun( link_values, dispersion_ );

    design->update_derivative(orders, W_ma, family);
    arma::mat derivatives_link = design->derivative_link;
    derivatives_link.each_col() %= residuals;
    arma::mat H(derivatives_link.n_cols, derivatives_link.n_cols);
    unsigned int dim = ts_vec.n_elem / orders->n_obs_effective;
    arma::vec score_temp(derivatives_link.n_cols);
    for(unsigned int t = 0; t < orders->n_obs_effective; t++)
    {
        score_temp = arma::sum(derivatives_link.rows(t * dim, (t + 1) * dim - 1), 0).t();
        H += score_temp * score_temp.t();
    }
    H /= orders->n_obs_effective;
    arma::mat variance = arma::inv_sympd(information / orders->n_obs_effective);
    return (variance * H * variance);
};


/*
    Helper function to calculate the link value at time t, i.e. applies the model equation for the linear process
    Does not return anything, but updates the link_values matrix at column t
*/

void Model::calculate_link_value_at(unsigned int &t, arma::mat &link_values, arma::mat &ts, const Orders &model_orders, const Family * fam, const arma::vec &intercept,
                                const CovariateList * covariates, const Neighborhood * W_ar, const Neighborhood * W_ma, const Neighborhood * W_covariates, const bool &ignore_covariates)
{
    // link_values.col(t) = intercept;
    arma::vec temp_link(ts.n_rows);
    arma::vec temp(ts.n_rows);

    for(unsigned int i = 0; i < model_orders.ma_time_lags.n_elem; i++)
    {
        temp = link_values.col(t - model_orders.ma_time_lags(i));
        temp = fam->link_trafo( temp );
        temp_link += W_ma->multiply_param_with_x(i, temp);
    }
    temp_link += intercept;
    for(unsigned int i = 0; i < model_orders.ar_time_lags.n_elem; i++)
    {
        temp = fam->observation_trafo( ts.col(t - model_orders.ar_time_lags(i)) );
        temp_link += W_ar->multiply_param_with_x(i, temp);
    }
    if(!ignore_covariates)
    {
        for(unsigned int i = 0; i < model_orders.covariate_orders.n_cols; i++)
        {
            temp = covariates->get_values_at(i, t);
            temp_link += W_covariates->multiply_param_with_x(i, temp);
        }
    }
    link_values.col(t) = temp_link;
}



