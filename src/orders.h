/* 
-----------------------------------------------------------------------------
    File: orders.h
    Purpose: Declarations for Orders class, to store model orders
    Author: Steffen Maletz
    Last modified: 2025-12-07
-----------------------------------------------------------------------------
*/

#ifndef ORDERS_H
#define ORDERS_H

class Orders
{
    private:
    public:
    // Model orders
    bool intercepts_equal;
    arma::umat autoregressive_orders;
    arma::uvec which_autoregressive_orders;
    arma::uvec ar_time_lags;
    arma::umat moving_average_orders;
    arma::uvec which_moving_average_orders;
    arma::uvec ma_time_lags;
    arma::umat covariate_orders;
    arma::uvec which_covariate_orders;

    // Additional information about model:
    unsigned int max_time_lag;
    unsigned int max_ma_lag;
    unsigned int n_param_intercept;
    unsigned int n_param_ar;
    unsigned int n_param_ma;
    unsigned int n_param_cov;
    unsigned int n_param;
    unsigned int n_obs_effective;
                                                                                
    Orders(const Rcpp::List &model_definition, const unsigned int &dim, const unsigned int &obs)
    {
        Rcpp::CharacterVector intercept_definition = model_definition["intercept"];
        intercepts_equal = (intercept_definition[0] == "homogeneous");
        n_param_intercept = intercepts_equal ? 1 : dim;
        n_param_ar = 0;
        n_param_ma = 0;
        n_param_cov = 0;
        n_param = n_param_intercept;
        max_time_lag = 0;
        max_ma_lag = 0;

        if(model_definition.containsElementNamed("past_obs")){
            Rcpp::IntegerMatrix ar_orders = model_definition["past_obs"];
            Rcpp::IntegerVector ar_lags = model_definition["past_obs_time_lags"];
            autoregressive_orders = Rcpp::as<arma::umat>( ar_orders );
            ar_time_lags = Rcpp::as<arma::uvec>( ar_lags );
            which_autoregressive_orders = arma::find(autoregressive_orders != 0);
            max_time_lag = ar_time_lags.max();
            n_param_ar = arma::accu(autoregressive_orders);
            n_param += n_param_ar;
        } 

        if(model_definition.containsElementNamed("past_mean")){
            Rcpp::IntegerMatrix ma_orders = model_definition["past_mean"];
            Rcpp::IntegerVector ma_lags = model_definition["past_mean_time_lags"];
            moving_average_orders = Rcpp::as<arma::umat>( ma_orders );
            ma_time_lags = Rcpp::as<arma::uvec>( ma_lags );
            which_moving_average_orders = arma::find(moving_average_orders != 0);
            max_ma_lag = ma_time_lags.max();
            max_time_lag = max_ma_lag < max_time_lag ? max_time_lag : max_ma_lag;
            n_param_ma = arma::accu(moving_average_orders);
            n_param += n_param_ma;
        }

        if(model_definition.containsElementNamed("covariates")){
            Rcpp::IntegerMatrix cov_orders = model_definition["covariates"];
            covariate_orders = Rcpp::as<arma::umat>( cov_orders );
            which_covariate_orders = arma::find(covariate_orders != 0);
            n_param_cov = arma::accu(covariate_orders);
            n_param += n_param_cov;
        }
        n_obs_effective = obs - max_time_lag;
    };
    void add_covariate_orders(arma::uvec &orders_to_add, arma::uvec &external_to_add){
        arma::umat cov_orders_to_add(orders_to_add.max(), orders_to_add.n_elem);
        for(unsigned int k = 0; k < orders_to_add.n_elem; k++){
            cov_orders_to_add.submat(0, k, orders_to_add(k) - 1, k).fill(1);
        }
        if(covariate_orders.n_rows < cov_orders_to_add.n_rows){
            arma::umat zero_rows(cov_orders_to_add.n_rows - covariate_orders.n_rows, covariate_orders.n_cols);
            covariate_orders = arma::join_cols(covariate_orders, zero_rows);
        }
        covariate_orders = arma::join_rows(covariate_orders, cov_orders_to_add);
        which_covariate_orders = arma::find(covariate_orders != 0);
        n_param_cov = arma::accu(covariate_orders);
        n_param += arma::accu(orders_to_add);
    }
    Orders(Orders* to_clone) : intercepts_equal(to_clone->intercepts_equal),
                               autoregressive_orders(to_clone->autoregressive_orders),
                               which_autoregressive_orders(to_clone->which_autoregressive_orders),
                               ar_time_lags(to_clone->ar_time_lags),
                               moving_average_orders(to_clone->moving_average_orders),
                               which_moving_average_orders(to_clone->which_moving_average_orders),
                               ma_time_lags(to_clone->ma_time_lags),
                               covariate_orders(to_clone->covariate_orders),
                               which_covariate_orders(to_clone->which_covariate_orders),
                               max_time_lag(to_clone->max_time_lag),
                               max_ma_lag(to_clone->max_ma_lag),
                               n_param_intercept(to_clone->n_param_intercept),
                               n_param_ar(to_clone->n_param_ar),
                               n_param_ma(to_clone->n_param_ma),
                               n_param_cov(to_clone->n_param_cov),
                               n_param(to_clone->n_param),
                               n_obs_effective(to_clone->n_obs_effective) {};
    Orders* clone(){
        return new Orders(*this);
    }
    ~Orders() = default;
};

#endif
