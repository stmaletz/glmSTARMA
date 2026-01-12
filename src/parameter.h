/* 
-----------------------------------------------------------------------------
    File: parameter.h
    Purpose: Declarations for Parameter classes, to store model parameters
    Author: Steffen Maletz
    Last modified: 2026-01-11
-----------------------------------------------------------------------------
*/

#ifndef PARAMETER_H
#define PARAMETER_H

/*
    Parameter:
    Class to store model parameters as 
    * Vector of all (unknown) model parameters
    * Matrices for intercept, AR, MA, and covariate parameters
    Additional methods to set and get these representations
*/

class Parameter
{
    public:
    arma::vec param_vector;
    arma::vec intercept;
    arma::mat ar_params;
    arma::mat ma_params;
    arma::mat cov_params;

    Parameter(const arma::vec &param, const Orders &orders) : param_vector(param)
    {
        set_param_matrices(orders);
    };
    Parameter(const Rcpp::List &params, const Orders &orders) 
    {
        intercept = Rcpp::as<arma::vec>(params["intercept"]);
        if(params.containsElementNamed("past_obs")){
            ar_params = Rcpp::as<arma::mat>(params["past_obs"]);
        }
        if(params.containsElementNamed("past_mean")){
            ma_params = Rcpp::as<arma::mat>(params["past_mean"]);
        }
        if(params.containsElementNamed("covariates")){
            cov_params = Rcpp::as<arma::mat>(params["covariates"]);
        }
        set_param_vector(orders);
    }
    Parameter(const Rcpp::List &params) 
    {
        intercept = Rcpp::as<arma::vec>(params["intercept"]);
        if(params.containsElementNamed("past_obs")){
            ar_params = Rcpp::as<arma::mat>(params["past_obs"]);
        }
        if(params.containsElementNamed("past_mean")){
            ma_params = Rcpp::as<arma::mat>(params["past_mean"]);
        }
        if(params.containsElementNamed("covariates")){
            cov_params = Rcpp::as<arma::mat>(params["covariates"]);
        }
    }
    Parameter(const arma::vec &param) : param_vector(param) {};
    Parameter(Parameter* to_clone) : param_vector(to_clone->param_vector),
                                     intercept(to_clone->intercept),
                                     ar_params(to_clone->ar_params),
                                     ma_params(to_clone->ma_params),
                                     cov_params(to_clone->cov_params) {};
    ~Parameter() = default;
    void set_param_matrices(const Orders &orders){
        if(orders.n_param_ma > 0)
        {
            ma_params = arma::mat(orders.moving_average_orders.n_rows, orders.moving_average_orders.n_cols, arma::fill::zeros);
            ma_params.elem(orders.which_moving_average_orders) = param_vector.head(orders.n_param_ma);
        }

        if(orders.intercepts_equal){
            intercept = arma::vec(1, arma::fill::value(param_vector(orders.n_param_ma)));
        } else {
            intercept = param_vector.subvec(orders.n_param_ma, orders.n_param_ma + orders.n_param_intercept - 1);
        }

        if(orders.n_param_ar > 0){
            ar_params = arma::mat(orders.autoregressive_orders.n_rows, orders.autoregressive_orders.n_cols, arma::fill::zeros);
            ar_params.elem(orders.which_autoregressive_orders) = param_vector.subvec(orders.n_param_ma + orders.n_param_intercept, orders.n_param_ma + orders.n_param_intercept + orders.n_param_ar - 1);
        }

        if(orders.n_param_cov > 0){
            cov_params = arma::mat(orders.covariate_orders.n_rows, orders.covariate_orders.n_cols, arma::fill::zeros);
            cov_params.elem(orders.which_covariate_orders) = param_vector.subvec(orders.n_param_ma + orders.n_param_intercept + orders.n_param_ar, orders.n_param - 1);
        }
    };
    void add_parameters(const arma::vec &params_to_add){
        param_vector = arma::join_vert(param_vector, params_to_add);
    };
    void set_param_vector(const Orders &orders)
    {
        param_vector = arma::vec(orders.n_param);
        if(orders.n_param_ma > 0){
            param_vector.head(orders.n_param_ma) = arma::vectorise(ma_params.elem(orders.which_moving_average_orders));
        }
        param_vector.subvec(orders.n_param_ma, orders.n_param_ma + orders.n_param_intercept - 1) = intercept;
        if(orders.n_param_ar > 0){
            param_vector.subvec(orders.n_param_ma + orders.n_param_intercept, orders.n_param_ma + orders.n_param_intercept + orders.n_param_ar - 1) = arma::vectorise(ar_params.elem(orders.which_autoregressive_orders));
        }
        if(orders.n_param_cov > 0){
            param_vector.tail(orders.n_param_cov) = arma::vectorise(cov_params.elem(orders.which_covariate_orders));
        }
    };
    Parameter* clone(){
        return new Parameter(*this);
    };
    Rcpp::List return_param_list()
    {
        Rcpp::List param_list = Rcpp::List::create(Rcpp::Named("intercept") = intercept,
                                                   Rcpp::Named("past_obs") = ar_params,
                                                   Rcpp::Named("past_mean") = ma_params,
                                                   Rcpp::Named("covariates") = cov_params);
        return param_list;
    }


    static void create_param_matrices(const arma::vec &param, const Orders &orders, arma::vec &intercept, arma::mat &ar_params, arma::mat &ma_params, arma::mat &cov_params)
    {
        if(orders.n_param_ma > 0)
        {
            ma_params = arma::mat(orders.moving_average_orders.n_rows, orders.moving_average_orders.n_cols, arma::fill::zeros);
            ma_params.elem(orders.which_moving_average_orders) = param.head(orders.n_param_ma);
        }

        if(orders.intercepts_equal){
            intercept = arma::vec(1, arma::fill::value(param(orders.n_param_ma)));
        } else {
            intercept = param.subvec(orders.n_param_ma, orders.n_param_ma + orders.n_param_intercept - 1);
        }

        if(orders.n_param_ar > 0){
            ar_params = arma::mat(orders.autoregressive_orders.n_rows, orders.autoregressive_orders.n_cols, arma::fill::zeros);
            ar_params.elem(orders.which_autoregressive_orders) = param.subvec(orders.n_param_ma + orders.n_param_intercept, orders.n_param_ma + orders.n_param_intercept + orders.n_param_ar - 1);
        }

        if(orders.n_param_cov > 0){
            cov_params = arma::mat(orders.covariate_orders.n_rows, orders.covariate_orders.n_cols, arma::fill::zeros);
            cov_params.elem(orders.which_covariate_orders) = param.subvec(orders.n_param_ma + orders.n_param_intercept + orders.n_param_ar, orders.n_param - 1);
        }
    };
};


#endif