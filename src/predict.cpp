/* 
-----------------------------------------------------------------------------
    File: predict.cpp
    Purpose: Implementation of predict functions for glmstarma and dglmstarma models
    Author: Steffen Maletz
    Last modified: 2025-12-07
-----------------------------------------------------------------------------
*/

#include "glmstarma.h"


/*
    glmstarma_predict: Function to make predictions in the future for a glmstarma model
    Input:
    * n_ahead: How many time points to predict ahead
    * pred_type: Type of prediction, "sample", "response", "link"
    * ts: Time series to predict (automatically extracted from model object in R)
    * covariate_list: List of covariates (automatically extracted from model object in R)
    * model: Model to use (automatically extracted from model object in R)
    * wlist_ar: List of neighborhoods for autoregressive parameters (automatically extracted from model object in R)
    * wlist_ma: List of neighborhoods for moving average parameters (automatically extracted from model object in R)
    * wlist_covariates: List of neighborhoods for covariates (automatically extracted from model object in R)
    * family: marginal distribution family (automatically extracted from model object in R)
    * new_covariates: New covariate values for prediction. If no new values are provided, the last known values are used
    * new_obs: New observations. These are used as ground truth to predict observations more than one step ahead (moving window prediction)
    * fitted_values: Fitted values from the model (automatically extracted from model object in R)
    * control: Control parameters argument (automatically extracted from model object in R)
    * parameter_est: Estimated parameters from the model (automatically extracted from model object in R)
    * copula_obj: Optional, if a copula is used (automatically extracted from model object in R)
    Output:
    * Matrix of predictions with dimensions (dim x n_ahead)
*/


// [[Rcpp::export]]
arma::mat glmstarma_predict(const unsigned int &n_ahead, const std::string &pred_type, const arma::mat &ts, const Rcpp::List &covariate_list,
    const Rcpp::List &model, const Rcpp::List &wlist_ar, const Rcpp::Nullable<Rcpp::List> wlist_ma, const Rcpp::Nullable<Rcpp::List> wlist_covariates, 
    const Rcpp::List &family, const Rcpp::Nullable<Rcpp::List> &new_covariates, const Rcpp::Nullable<arma::mat> &new_obs, const arma::mat &fitted_values, const Rcpp::List &control,
    const Rcpp::List &parameter_est, Rcpp::Nullable<Rcpp::S4> copula_obj)
{
    const bool use_sparsity = Rcpp::as<bool>(control["use_sparsity"]);
    const double sparsity_threshold = Rcpp::as<double>(control["sparsity_threshold"]);

    const unsigned int n_obs = ts.n_cols + n_ahead;
    const unsigned int dim = ts.n_rows;
    CovariateList covariates(covariate_list, n_obs, dim, 0, 0);
    Orders orders(model, dim, n_obs);

    if(orders.n_param_cov > 0 && covariates.has_time_variant_covariates() && new_covariates.isNull())
    {
        Rcpp::warning("The model includes time-variant covariates, but no new covariate values were provided. The last known values will be used for prediction.");
    }
    if(orders.n_param_cov > 0 && !new_covariates.isNull())
    {
        Rcpp::List new_covariates_list(new_covariates);
        CovariateList newx(new_covariates_list, n_ahead, dim, 0, 0);
        covariates.append_covariate_values(newx);
    }

    Family * fam = nullptr;
    if(pred_type == "sample" && copula_obj.isNotNull())
    {
        Rcpp::S4 copula_object(copula_obj);
        fam = Family::create(family, copula_object);
    } else {
        fam = Family::create(family);
    }
    if(pred_type == "sample" && !fam->const_dispersion)
    {
        Rcpp::stop("Prediction by sampling is only supported for models with constant dispersion parameter.");
    }    

    Neighborhood * W_ar = Neighborhood::create(wlist_ar, sparsity_threshold, use_sparsity);
    Neighborhood * W_ma = Neighborhood::create_neighbor(orders.n_param_ma, wlist_ma, W_ar, dim, sparsity_threshold, use_sparsity);
    Neighborhood * W_covariates = Neighborhood::create_neighbor(orders.n_param_cov, wlist_covariates, W_ar, dim, sparsity_threshold, use_sparsity);

    arma::vec param_est = Model::init_param(parameter_est, orders, ts, fam);
    arma::vec intercept(orders.n_param_intercept);
    arma::mat ar_params(orders.autoregressive_orders.n_rows, orders.autoregressive_orders.n_cols);
    arma::mat ma_params(orders.moving_average_orders.n_rows, orders.moving_average_orders.n_cols);
    arma::mat cov_params(orders.covariate_orders.n_rows, orders.covariate_orders.n_cols);
    Parameter::create_param_matrices(param_est, &orders, intercept, ar_params, ma_params, cov_params);
    W_ar->set_parameter_matrices(ar_params);
    W_ma->set_parameter_matrices(ma_params);
    W_covariates->set_parameter_matrices(cov_params);

    arma::mat ts_extended(ts.n_rows, ts.n_cols + n_ahead, arma::fill::zeros);
    ts_extended.head_cols(ts.n_cols) = ts;
    arma::mat link_values_extended(ts.n_rows, ts.n_cols + n_ahead, arma::fill::zeros);
    link_values_extended.head_cols(ts.n_cols) = fam->link(fitted_values);
    arma::mat predictions(ts.n_rows, n_ahead, arma::fill::zeros);

    unsigned int n_new_obs = 0;
    arma::mat new_obs_;
    if(new_obs.isNotNull())
    {
        new_obs_ = Rcpp::as<arma::mat>(new_obs);
        n_new_obs = new_obs_.n_cols;
        ts_extended.cols(ts.n_cols, ts.n_cols + n_ahead - 1) = new_obs_;
    }

    for(unsigned int t = ts.n_cols; t < n_obs; t++)
    {
        Model::calculate_link_value_at(t, link_values_extended, ts_extended, orders, fam, intercept, &covariates, W_ar, W_ma, W_covariates, false);
        arma::vec prediction = link_values_extended.col(t);
        prediction = (pred_type != "link") ? fam->inverse_link(prediction) : prediction;
        prediction = (pred_type == "sample") ? fam->sample(prediction) : prediction; // TODO: Dispersion beachten
        predictions.col(t - ts.n_cols) = prediction;
        if(t >= ts.n_cols + n_new_obs)
        {
            ts_extended.col(t) = (pred_type == "link") ? fam->inverse_link(prediction) : prediction;
        }
    }

    delete W_ar;
    delete W_ma;
    delete W_covariates;
    delete fam;

    return predictions;
}




/*
    dglmstarma_predict: Function to make predictions in the future for a dglmstarma model
    Input:
    * n_ahead: How many time points to predict ahead
    * pred_type_mean: Type of prediction for the mean model, "sample", "response", "link"
    * pred_type_dispersion: Type of prediction for the dispersion model, "response", "link"
    * ts: Time series to predict (automatically extracted from model object in R)
    * mean_model: Mean model to use (automatically extracted from model object in R)
    * dispersion_model: Dispersion model to use (automatically extracted from model object in R)
    * pseudo_obs_type: Type of pseudo observations used (automatically extracted from model object in R)
    * mean_family: Marginal distribution family for the mean model (automatically extracted from model object in R)
    * dispersion_family: Marginal distribution family for the dispersion model (automatically extracted from model object in R)
    * wlist: List of neighborhoods for the mean model (automatically extracted from model object in R)
    * mean_covariates: List of covariates for the mean model (automatically extracted from model object in R)
    * dispersion_covariates: List of covariates for the dispersion model (automatically extracted from model object in R)
    * pseudo_observations: Pseudo observations to predict (automatically extracted from model object in R)
    * wlist_past_mean: List of neighborhoods for past mean values in the dispersion model (automatically extracted from model object in R)
    * wlist_covariates: List of neighborhoods for covariates in the mean model (automatically extracted from model object in R)
    * wlist_ar_dispersion: List of neighborhoods for autoregressive parameters in the dispersion model (automatically extracted from model object in R)
    * wlist_past_mean_dispersion: List of neighborhoods for past mean values in the dispersion model (automatically extracted from model object in R)
    * wlist_covariates_dispersion: List of neighborhoods for covariates in the dispersion model (automatically extracted from model object in R)
    * parameter_est_mean: Estimated parameters from the mean model (automatically extracted from model object in R)
    * parameter_est_dispersion: Estimated parameters from the dispersion model (automatically extracted from model object in R)
    * covariate_list_mean_new: New covariate values for the mean model prediction. If no new values are provided, the last known values are used
    * covariate_list_dispersion_new: New covariate values for the dispersion model prediction. If no new values are provided, the last known values are used
    * new_obs_mean: New observations for the mean model. These are used as ground truth to predict observations more than one step ahead (moving window prediction)
    * fitted_values_mean: Fitted values from the mean model (automatically extracted from model object in R)
    * fitted_values_dispersion: Fitted values from the dispersion model (automatically extracted from model object in R)
    * control: Control parameters argument (automatically extracted from model object in R)
    * copula_obj_mean: Optional, if a copula is used in the mean model (automatically extracted from model object in R)
    Output:
    List with:
    * Matrix of predictions for the mean model with dimensions (dim x n_ahead)
    * Matrix of predictions for the dispersion model with dimensions (dim x n_ahead)
*/


// [[Rcpp::export]]
Rcpp::List dglmstarma_predict(const unsigned int &n_ahead, const std::string &pred_type_mean, const std::string &pred_type_dispersion, 
    const arma::mat &ts, const Rcpp::List &mean_model, const Rcpp::List &dispersion_model, const std::string &pseudo_obs_type,
    const Rcpp::List &mean_family, const Rcpp::List &dispersion_family, const Rcpp::List &wlist,
    const Rcpp::List &mean_covariates, const Rcpp::List &dispersion_covariates,
    const arma::mat &pseudo_observations, const Rcpp::Nullable<Rcpp::List> wlist_past_mean,
    const Rcpp::Nullable<Rcpp::List> wlist_covariates, const Rcpp::Nullable<Rcpp::List> wlist_ar_dispersion,
    const Rcpp::Nullable<Rcpp::List> wlist_past_mean_dispersion, const Rcpp::Nullable<Rcpp::List> wlist_covariates_dispersion,
    const Rcpp::List &parameter_est_mean, const Rcpp::List &parameter_est_dispersion, const Rcpp::Nullable<Rcpp::List> &covariate_list_mean_new,
    const Rcpp::Nullable<Rcpp::List> &covariate_list_dispersion_new, const Rcpp::Nullable<arma::mat> new_obs_mean,  
    const arma::mat &fitted_values_mean, const arma::mat &fitted_values_dispersion, const Rcpp::List &control,
    const Rcpp::Nullable<Rcpp::S4> copula_obj_mean)
{
    const bool use_sparsity = Rcpp::as<bool>(control["use_sparsity"]);
    const double sparsity_threshold = Rcpp::as<double>(control["sparsity_threshold"]);

    const unsigned int n_obs_mean =  ts.n_cols + n_ahead;
    const unsigned int n_obs_dispersion = pseudo_observations.n_cols + n_ahead;
    const unsigned int dim = ts.n_rows;

    CovariateList covariates_mean(mean_covariates, n_obs_mean, dim, 0, 0);
    Orders orders_mean(mean_model, dim, n_obs_mean);
    CovariateList covariates_dispersion(dispersion_covariates, n_obs_dispersion, dim, 0, 0);
    Orders orders_dispersion(dispersion_model, dim, n_obs_dispersion);

    if(orders_mean.n_param_cov > 0 && covariates_mean.has_time_variant_covariates() && covariate_list_mean_new.isNull())
    {
        Rcpp::warning("The mean model includes time-variant covariates, but no new covariate values were provided. The last known values will be used for prediction.");
    }
    if(orders_dispersion.n_param_cov > 0 && covariates_dispersion.has_time_variant_covariates() && covariate_list_dispersion_new.isNull())
    {
        Rcpp::warning("The dispersion model includes time-variant covariates, but no new covariate values were provided. The last known values will be used for prediction.");
    }
    if(orders_mean.n_param_cov > 0 && covariate_list_mean_new.isNotNull())
    {
        Rcpp::List new_covariates_list(covariate_list_mean_new);
        CovariateList newx(new_covariates_list, n_ahead, dim, 0, 0);
        covariates_mean.append_covariate_values(newx);
    }
    if(orders_dispersion.n_param_cov > 0 && covariate_list_dispersion_new.isNotNull())
    {
        Rcpp::List new_covariates_list(covariate_list_dispersion_new);
        CovariateList newx(new_covariates_list, n_ahead, dim, 0, 0);
        covariates_dispersion.append_covariate_values(newx);
    }
    Family * fam_mean = nullptr;
    Family * fam_dispersion = Family::create(dispersion_family);
    if(pred_type_mean == "sample" && copula_obj_mean.isNotNull())
    {
        Rcpp::S4 copula_object(copula_obj_mean);
        fam_mean = Family::create(mean_family, copula_object);
    } else {
        fam_mean = Family::create(mean_family);
    }
    Neighborhood * W_ar_mean = Neighborhood::create(wlist, sparsity_threshold, use_sparsity);
    Neighborhood * W_ma_mean = Neighborhood::create_neighbor(orders_mean.n_param_ma, wlist_past_mean, W_ar_mean, dim, sparsity_threshold, use_sparsity);
    Neighborhood * W_covariates_mean = Neighborhood::create_neighbor(orders_mean.n_param_cov, wlist_covariates, W_ar_mean, dim, sparsity_threshold, use_sparsity);
    Neighborhood * W_ar_dispersion = Neighborhood::create_neighbor(orders_dispersion.n_param_ar, wlist_ar_dispersion, W_ar_mean, dim, sparsity_threshold, use_sparsity);
    Neighborhood * W_ma_dispersion = Neighborhood::create_neighbor(orders_dispersion.n_param_ma, wlist_past_mean_dispersion, W_ar_mean, dim, sparsity_threshold, use_sparsity);
    Neighborhood * W_covariates_dispersion = Neighborhood::create_neighbor(orders_dispersion.n_param_cov, wlist_covariates_dispersion, W_ar_mean, dim, sparsity_threshold, use_sparsity);

    arma::vec param_est_mean = Model::init_param(parameter_est_mean, orders_mean, ts, fam_mean);
    arma::vec intercept_mean(orders_mean.n_param_intercept);
    arma::mat ar_params_mean(orders_mean.autoregressive_orders.n_rows, orders_mean.autoregressive_orders.n_cols);
    arma::mat ma_params_mean(orders_mean.moving_average_orders.n_rows, orders_mean.moving_average_orders.n_cols);
    arma::mat cov_params_mean(orders_mean.covariate_orders.n_rows, orders_mean.covariate_orders.n_cols);
    Parameter::create_param_matrices(param_est_mean, &orders_mean, intercept_mean, ar_params_mean, ma_params_mean, cov_params_mean);
    W_ar_mean->set_parameter_matrices(ar_params_mean);
    W_ma_mean->set_parameter_matrices(ma_params_mean);
    W_covariates_mean->set_parameter_matrices(cov_params_mean);

    arma::vec param_est_dispersion = Model::init_param(parameter_est_dispersion, orders_dispersion, pseudo_observations, fam_dispersion);
    arma::vec intercept_dispersion(orders_dispersion.n_param_intercept);
    arma::mat ar_params_dispersion(orders_dispersion.autoregressive_orders.n_rows, orders_dispersion.autoregressive_orders.n_cols);
    arma::mat ma_params_dispersion(orders_dispersion.moving_average_orders.n_rows, orders_dispersion.moving_average_orders.n_cols);
    arma::mat cov_params_dispersion(orders_dispersion.covariate_orders.n_rows, orders_dispersion.covariate_orders.n_cols);
    Parameter::create_param_matrices(param_est_dispersion, &orders_dispersion, intercept_dispersion, ar_params_dispersion, ma_params_dispersion, cov_params_dispersion);
    W_ar_dispersion->set_parameter_matrices(ar_params_dispersion);
    W_ma_dispersion->set_parameter_matrices(ma_params_dispersion);
    W_covariates_dispersion->set_parameter_matrices(cov_params_dispersion);

    arma::mat ts_extended(ts.n_rows, ts.n_cols + n_ahead, arma::fill::zeros);
    ts_extended.head_cols(ts.n_cols) = ts;
    arma::mat pseudo_observations_extended(ts.n_rows, ts.n_cols + n_ahead, arma::fill::zeros);
    pseudo_observations_extended.cols(orders_mean.max_time_lag, orders_mean.max_time_lag + pseudo_observations.n_cols - 1) = pseudo_observations;

    arma::mat link_values_extended(ts.n_rows, ts.n_cols + n_ahead, arma::fill::zeros);
    link_values_extended.head_cols(fitted_values_mean.n_cols) = fam_mean->link(fitted_values_mean);
    arma::mat dispersion_link_values_extended(ts.n_rows, ts.n_cols + n_ahead, arma::fill::zeros);
    dispersion_link_values_extended.head_cols(fitted_values_dispersion.n_cols) = fam_dispersion->link(fitted_values_dispersion);

    arma::mat predictions_mean(ts.n_rows, n_ahead, arma::fill::zeros);
    arma::mat predictions_dispersion(ts.n_rows, n_ahead, arma::fill::zeros);
    unsigned int n_new_obs = 0;
    if(new_obs_mean.isNotNull())
    {
        arma::mat new_obs_mean_ = Rcpp::as<arma::mat>(new_obs_mean);
        n_new_obs = new_obs_mean_.n_cols;
        ts_extended.cols(ts.n_cols, ts.n_cols + n_new_obs - 1) = new_obs_mean_;
    }

    arma::vec temp_link(ts.n_rows);
    arma::vec temp_dispersion(ts.n_rows);
    arma::vec temp_prediction_mean(ts.n_rows);
    arma::vec temp_prediction_dispersion(ts.n_rows);

    arma::vec inter_mean(ts.n_rows);
    arma::vec inter_dispersion(ts.n_rows);
    if(intercept_mean.n_elem == 1){
        inter_mean.fill(intercept_mean(0));
    } else {
        inter_mean = intercept_mean;
    }
    if(intercept_dispersion.n_elem == 1){
        inter_dispersion.fill(intercept_dispersion(0));
    } else {
        inter_dispersion = intercept_dispersion;
    }

    for(unsigned int t = ts.n_cols; t < ts.n_cols + std::min(n_new_obs, n_ahead); t++)
    {
        Model::calculate_link_value_at(t, dispersion_link_values_extended, pseudo_observations_extended, orders_dispersion, fam_dispersion, inter_dispersion, &covariates_dispersion, W_ar_dispersion, W_ma_dispersion, W_covariates_dispersion, false);
        temp_prediction_dispersion = dispersion_link_values_extended.col(t);
        temp_dispersion = fam_dispersion->inverse_link(temp_prediction_dispersion);
        temp_prediction_dispersion = (pred_type_dispersion != "link") ? temp_dispersion : temp_prediction_dispersion;
        Model::calculate_link_value_at(t, link_values_extended, ts_extended, orders_mean, fam_mean, inter_mean, &covariates_mean, W_ar_mean, W_ma_mean, W_covariates_mean, false);
        temp_prediction_mean = link_values_extended.col(t);
        temp_prediction_mean = (pred_type_mean != "link") ? fam_mean->inverse_link(temp_prediction_mean) : temp_prediction_mean;
        temp_prediction_mean = (pred_type_mean == "sample") ? fam_mean->sample(temp_prediction_mean, temp_dispersion) : temp_prediction_mean; 
        predictions_mean.col(t - ts.n_cols) = temp_prediction_mean;
        predictions_dispersion.col(t - ts.n_cols) = temp_prediction_dispersion;
        pseudo_observations_extended.col(t) = (pseudo_obs_type == "deviance") ? fam_mean->deviance_residual(ts_extended.col(t), temp_prediction_mean) : fam_mean->pearson_residual(ts_extended.col(t), temp_prediction_mean);
    }
    if(n_ahead > n_new_obs)
    {
        for(unsigned int t = ts.n_cols + n_new_obs; t < ts.n_cols + n_ahead; t++)
        {
            Model::calculate_link_value_at(t, dispersion_link_values_extended, pseudo_observations_extended, orders_dispersion, fam_dispersion, inter_dispersion, &covariates_dispersion, W_ar_dispersion, W_ma_dispersion, W_covariates_dispersion, false);
            temp_prediction_dispersion = dispersion_link_values_extended.col(t);
            temp_dispersion = fam_dispersion->inverse_link(temp_prediction_dispersion);
            temp_prediction_dispersion = (pred_type_dispersion != "link") ? temp_dispersion : temp_prediction_dispersion;
            Model::calculate_link_value_at(t, link_values_extended, ts_extended, orders_mean, fam_mean, inter_mean, &covariates_mean, W_ar_mean, W_ma_mean, W_covariates_mean, false);
            temp_prediction_mean = link_values_extended.col(t);
            temp_prediction_mean = (pred_type_mean != "link") ? fam_mean->inverse_link(temp_prediction_mean) : temp_prediction_mean;
            temp_prediction_mean = (pred_type_mean == "sample") ? fam_mean->sample(temp_prediction_mean, temp_dispersion) : temp_prediction_mean; 
            predictions_dispersion.col(t - ts.n_cols) = temp_prediction_dispersion;
            predictions_mean.col(t - ts.n_cols) = temp_prediction_mean;
            ts_extended.col(t) = (pred_type_mean == "link") ? fam_mean->inverse_link(temp_prediction_mean) : temp_prediction_mean;
            pseudo_observations_extended.col(t) = (pred_type_mean != "sample") ? temp_dispersion : ((pseudo_obs_type == "deviance") ? fam_mean->deviance_residual(ts_extended.col(t), temp_prediction_mean) : fam_mean->pearson_residual(ts_extended.col(t), temp_prediction_mean));
        }
    }

    delete W_ar_mean;
    delete W_ma_mean;
    delete W_covariates_mean;
    delete W_ar_dispersion;
    delete W_ma_dispersion;
    delete W_covariates_dispersion;
    delete fam_mean;
    delete fam_dispersion;
    Rcpp::List result;
    result["mean"] = predictions_mean;
    result["dispersion"] = predictions_dispersion;
    return result;
}