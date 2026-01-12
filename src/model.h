/* 
-----------------------------------------------------------------------------
    File: model.h
    Purpose: Declaration of Model class with static functions for model fitting
        and initialization
    Author: Steffen Maletz
    Last modified: 2025-12-06
-----------------------------------------------------------------------------
*/


#ifndef MODEL_H
#define MODEL_H

/*
    Class with static functions
    - for initialization of parameter estimation,
    - initialization of link values
    - calculation of log-likelihood
    - calculation of the score
    - calculation of the information
    - calculation of the variance of parameter estimation
*/

class Model
{
    public:
    static arma::mat init_link(const Rcpp::RObject &method, const Orders &orders, const arma::mat &ts, const Family* fam);
    static arma::vec init_param(const Rcpp::RObject &method, const Orders &orders, const arma::mat &ts, const Family* fam);
    static double log_likelihood(const arma::vec &x, arma::vec &score, arma::mat &information, const arma::vec &ts_vec, Design * design, const Orders * orders, Neighborhood * W_ma, const Family * family, bool calc_score, bool calc_information);
    static arma::mat variance_estimation(const arma::vec &x, const arma::vec &ts_vec, Design * design, const Family * family, const arma::mat &information, Neighborhood * W_ma, const Orders * orders);
    static void calculate_link_value_at(unsigned int &t, arma::mat &link_values, arma::mat &ts, const Orders &model_orders, const Family * fam, const arma::vec &intercept, const CovariateList * covariates, const Neighborhood * W_ar, const Neighborhood * W_ma, const Neighborhood * W_covariates, const bool &ignore_covariates);
};


#endif
