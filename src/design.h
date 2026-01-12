/* 
-----------------------------------------------------------------------------
    File: design.h
    Purpose: Declaration of Design class
    Author: Steffen Maletz
    Last modified: 2025-12-06
-----------------------------------------------------------------------------
*/

#ifndef DESIGN_H
#define DESIGN_H

// Bei Update des link_designs fixe Startwerte beachten

/*
    Class to store the design matrix and related link values of the model
    If regression on past values of the feedback is included, some columns of
        the design matrix need to be updated in each iteration of the fitting process.
    The columns corresponding to the intercept(s), past observations of the time series
        and covariates remain constant throughout the fitting process.
    Contains functions to update the design matrix and link values accordingly
    Static function is provided to update derivative matrices for score test calculations
*/

class Design
{
    public:
    arma::mat transformed_link_values;
    arma::mat link_vals;
    arma::mat design_matrix;
    arma::mat derivative_link;

    Design(const arma::mat &ts, const CovariateList * covariates, const arma::mat &link_values, 
           const Orders * orders, const Family * fam, Neighborhood * W_ar, Neighborhood * W_ma, Neighborhood * W_covariates)
    {
        arma::mat transformed_obs = fam->observation_trafo(ts);
        link_vals = arma::mat(ts.n_rows, orders->n_obs_effective + orders->max_time_lag, arma::fill::zeros);
        transformed_link_values = arma::mat(ts.n_rows, orders->n_obs_effective + orders->max_time_lag, arma::fill::zeros);
        design_matrix = arma::mat(ts.n_rows * orders->n_obs_effective, orders->n_param);
        // Intercept
        if(orders->intercepts_equal){
            design_matrix.col(orders->n_param_ma) = arma::ones<arma::vec>(design_matrix.n_rows);
        } else {
            arma::mat identity = arma::eye(ts.n_rows, ts.n_rows);
            for(int i = 0; i < orders->n_obs_effective; i++){
                design_matrix.submat(i * ts.n_rows, 0, (i + 1) * ts.n_rows - 1, ts.n_rows - 1) = identity;
            }
        }
        // Past Observations
        unsigned int col_index = orders->n_param_ma + orders->n_param_intercept;
        arma::mat temp_obs;
        if(orders->n_param_ar > 0)
        {
            for(unsigned int i = 0; i < orders->autoregressive_orders.n_cols; i++)
            {
                temp_obs = transformed_obs.cols(orders->max_time_lag - orders->ar_time_lags(i), orders->max_time_lag - orders->ar_time_lags(i) + orders->n_obs_effective - 1);
                arma::uvec orders_lag = arma::find( orders->autoregressive_orders.col(i) );
                for(unsigned int o : orders_lag){
                    design_matrix.col(col_index) = arma::vectorise( W_ar->multiply_with_x(o, temp_obs) );
                    col_index++;
                }
            }
        }
        // Covariates
        if(orders->n_param_cov > 0)
        {
            design_matrix.tail_cols(orders->n_param_cov) = covariates->create_design_matrix(orders->covariate_orders, W_covariates).tail_rows(orders->n_obs_effective * ts.n_rows);
        }
        derivative_link = design_matrix;
        // Initial link values + transformation
        transformed_link_values = arma::mat(ts.n_rows, ts.n_cols, arma::fill::zeros);
        transformed_link_values.head_cols( orders->max_time_lag ) = fam->link_trafo( link_values );
        link_vals.head_cols( orders->max_time_lag ) = link_values;
    };
    arma::vec update_design(const arma::vec &param, const Orders * orders, const Family * fam, Neighborhood * W_ma)
    {
        arma::vec link_temp(transformed_link_values.n_rows);
        for(unsigned int t = orders->max_time_lag; t < transformed_link_values.n_cols; t++)
        {
            unsigned int col_index = 0;
            unsigned int t_rows = t - orders->max_time_lag;
            for(unsigned int i = 0; i < orders->moving_average_orders.n_cols; i++)
            {
                arma::uvec orders_lag = arma::find( orders->moving_average_orders.col(i) );
                for(unsigned int o : orders_lag)
                {
                    design_matrix.submat(t_rows * transformed_link_values.n_rows, col_index, (t_rows + 1) * transformed_link_values.n_rows - 1, col_index) = W_ma->multiply_with_x(o, transformed_link_values.col(t - orders->ma_time_lags(i)));
                    col_index++;
                }
            }
            link_temp = design_matrix.rows(t_rows * transformed_link_values.n_rows, (t_rows + 1) * transformed_link_values.n_rows - 1) * param;
            link_vals.col(t) = link_temp;
            transformed_link_values.col(t) = fam->link_trafo( link_temp );
        }
        return design_matrix * param;
    };
    void update_derivative(const Orders * orders, Neighborhood * W_ma, const Family * family)
    {
        if(orders->n_param_ma > 0)
        {
            unsigned int n_obs = orders->n_obs_effective + orders->max_time_lag;
            unsigned int dim = transformed_link_values.n_rows;
            unsigned int t_lagged;
            unsigned int t_index;
            bool changed_temp_mat = false;

            derivative_link = design_matrix;
            arma::mat temp_mat(dim, orders->n_param);

            for(unsigned int t = orders->max_time_lag; t < n_obs; t++)
            {
                temp_mat.zeros();
                t_index = t - orders->max_time_lag;
                changed_temp_mat = false;
                for(unsigned int i = 0; i < orders->moving_average_orders.n_cols; i++)
                {
                    t_lagged = t - orders->ma_time_lags(i);
                    changed_temp_mat = changed_temp_mat || (t_lagged > orders->max_time_lag);
                    if(t_lagged > orders->max_time_lag)
                    {
                        t_lagged -= orders->max_time_lag;
                        if(family->identity_link_trafo)
                        {
                            temp_mat += W_ma->multiply_param_with_x(i, derivative_link.rows(t_lagged * dim, (t_lagged + 1) * dim - 1));
                        } else {
                            arma::mat temp = derivative_link.rows(t_lagged * dim, (t_lagged + 1) * dim - 1);
                            temp.each_col() %= family->derivative_link_trafo(link_vals.col(t_lagged));
                            temp_mat += W_ma->multiply_param_with_x(i, temp);
                        }
                    }
                }
                if(changed_temp_mat)
                {
                    derivative_link.rows(t_index * dim, (t_index + 1) * dim - 1) += temp_mat;
                }
            }
        }
    };
    static void update_derivative_extern(arma::mat &deriv, const Orders * orders, Neighborhood * W_ma, const Family * family, const unsigned int &dim, const arma::mat &link_vals)
    {
        if(orders->n_param_ma > 0)
        {
            unsigned int n_obs = orders->n_obs_effective + orders->max_time_lag;
            unsigned int t_lagged;
            unsigned int t_index;
            bool changed_temp_mat = false;
            arma::mat temp_mat(dim, deriv.n_cols);

            for(unsigned int t = orders->max_time_lag; t < n_obs; t++)
            {
                temp_mat.zeros();
                t_index = t - orders->max_time_lag;
                changed_temp_mat = false;
                for(unsigned int i = 0; i < orders->moving_average_orders.n_cols; i++)
                {
                    t_lagged = t - orders->ma_time_lags(i);
                    changed_temp_mat = changed_temp_mat || (t_lagged > orders->max_time_lag);
                    if(t_lagged > orders->max_time_lag)
                    {
                        t_lagged -= orders->max_time_lag;
                        if(family->identity_link_trafo)
                        {
                            temp_mat += W_ma->multiply_param_with_x(i, deriv.rows(t_lagged * dim, (t_lagged + 1) * dim - 1));
                        } else {
                            arma::mat temp = deriv.rows(t_lagged * dim, (t_lagged + 1) * dim - 1);
                            temp.each_col() %= family->derivative_link_trafo(link_vals.col(t_lagged));
                            temp_mat += W_ma->multiply_param_with_x(i, temp);
                        }
                    }
                }
                if(changed_temp_mat)
                {
                    deriv.rows(t_index * dim, (t_index + 1) * dim - 1) += temp_mat;
                }
            }
        }
    };
    ~Design() = default;
};

#endif

