/* 
-----------------------------------------------------------------------------
    File: covariates.cpp
    Purpose: Implementation of constructors and functions for CovariateList class
    Author: Steffen Maletz
    Last modified: 2025-12-12
-----------------------------------------------------------------------------
*/


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include "glmstarma.h"
using namespace Rcpp;

/*
    Constructor to create CovariateList from R input

    Parameters:
    * covariates: Rcpp::List of covariates specified in glmstarma, dglmstarma or the simulation functions
    * n_obs: Number of time points (calculated automatically from data)
    * dim: Dimension of the multivariate time series (calculated automatically from data)
    * burn_in: Number of burn-in observations (to use correct time points later)
    * shift_interventions: Sift time-points of interventions covariates by this ammount (needed for dispersion model if max lag of mean model is dropped)
    
    Details:
    This constructor processes the R list of covariates and creates the appropriate Covariate objects (DefaultCovariate, TimeConstantCovariate, SpatialConstantCovariate, Intervention) based on the attributes of each covariate.
*/
CovariateList::CovariateList(const Rcpp::List &covariates, const unsigned int &n_obs, const unsigned int &dim, const int &burn_in, const unsigned int &shift_interventions) : n_obs(n_obs), dim(dim), burn_in(burn_in)
{
    n_covariates = covariates.size();
    if(n_covariates > 0)
    {
        Rcpp::CharacterVector covariate_names_cv = covariates.names();
        covariate_names = Rcpp::as<std::vector<std::string>>(covariate_names_cv);
        if(covariate_names.size() == 0){
            covariate_names.resize(n_covariates);
            for(unsigned int i = 0; i < n_covariates; i++){
                covariate_names[i] = "X" + std::to_string(i + 1);
            }
        }
        arma::vec temp(n_covariates); // For finding time-constant covariates
        for(unsigned int i = 0; i < n_covariates; i++){
            Rcpp::RObject cov_temp = covariates[i];
            if(cov_temp.hasAttribute("const")){
                Rcpp::CharacterVector attribute = cov_temp.attr("const");
                if(attribute[0] == "time"){
                    Rcpp::NumericVector values = Rcpp::as<Rcpp::NumericVector>( cov_temp );
                    arma::vec arma_vals = Rcpp::as<arma::vec>( values );
                    list_of_covariates.push_back( new TimeConstantCovariate(arma_vals, burn_in) );
                    temp(i) = 1;
                } else if(attribute[0] == "space"){
                    Rcpp::NumericVector values = Rcpp::as<Rcpp::NumericVector>( cov_temp );
                    arma::vec arma_vals = Rcpp::as<arma::vec>( values );
                    if(arma_vals.n_elem < n_obs) {
                        Rcpp::warning("At least one covariate has fewer number of observations in time than needed. Last value will be repeated to match the number of observations.");
                        arma::vec extended_vals(n_obs, arma::fill::zeros);
                        extended_vals.head(arma_vals.n_elem) = arma_vals;
                        extended_vals.tail(n_obs - arma_vals.n_elem).fill(arma_vals(arma_vals.n_elem - 1));
                        arma_vals = extended_vals;
                    }
                    arma_vals = arma_vals.head(n_obs);
                    list_of_covariates.push_back( new SpatialConstantCovariate(arma_vals, dim, burn_in) );
            }    else {
                    Rcpp::stop("Covariate not supported");
                }
            } else {
                Rcpp::NumericMatrix values = covariates[i];
                arma::mat arma_vals = Rcpp::as<arma::mat>( values );
                if(arma_vals.n_cols < n_obs) {
                    Rcpp::warning("At least one covariate has fewer number of observations in time than needed. Last value will be repeated to match the number of observations.");
                    arma::mat extended_vals(dim, n_obs, arma::fill::zeros);
                    extended_vals.head_cols(arma_vals.n_cols) = arma_vals;
                    extended_vals.tail_cols(n_obs - arma_vals.n_cols) = arma::repmat(arma_vals.col(arma_vals.n_cols - 1), 1, n_obs - arma_vals.n_cols);
                    arma_vals = extended_vals;
                }
                arma_vals = arma_vals.head_cols(n_obs);
                list_of_covariates.push_back( new DefaultCovariate(arma_vals, burn_in) );
            }
        }
        time_variant_covariates = arma::find(temp == 0);
    }
}

/*
    Copy constructor for CovariateList

    Parameters:
    * to_clone: Pointer to the CovariateList object to clone

    Details:
    This constructor creates a deep copy of the provided CovariateList object, including cloning each individual Covariate object in the list.
*/
CovariateList::CovariateList(CovariateList* to_clone) : n_obs(to_clone->n_obs), dim(to_clone->dim), time_variant_covariates(to_clone->time_variant_covariates), burn_in(to_clone->burn_in), n_covariates(to_clone->n_covariates)
{
    for(unsigned int i = 0; i < n_covariates; i++){
        list_of_covariates.push_back( to_clone->list_of_covariates.at(i)->clone() );
    }
}


/*
    Constructor to create empty CovariateList

    Parameters:
    * n_obs: Number of time points (calculated automatically from data)
    * dim: Dimension of the multivariate time series (calculated automatically from data)

    Details:
    This constructor initializes an empty CovariateList with specified number of observations and dimensions.
*/
CovariateList::CovariateList(const unsigned int &n_obs, const unsigned int &dim) : n_obs(n_obs), dim(dim), burn_in(0), n_covariates(0)
{
    time_variant_covariates = arma::uvec();
}

/*
    Function to add a single Covariate to the CovariateList

    Parameters:
    * cov: Pointer to the Covariate object to add

    Details:
    This function appends the provided Covariate object to the list of covariates and updates the count of covariates.
*/
void CovariateList::add_covariate(Covariate* cov)
{
    list_of_covariates.push_back(cov);
    n_covariates++;
    if(cov->is_time_constant){
        time_variant_covariates.resize(n_covariates);
        time_variant_covariates(n_covariates - 1) = 1;
    }
}



/*
    Function to get values of k-th covariate at time t
*/
arma::vec CovariateList::get_values_at(const int &k, const int &t) const {
    return (list_of_covariates.at(k))->get_values_at(t);
}

/*
    Function to clone the CovariateList object
*/
CovariateList* CovariateList::clone() {
    return new CovariateList(this);
}

/*
    Function to delete the last 'to_delete' covariates from the CovariateList
*/
void CovariateList::delete_last_covariates(unsigned int &to_delete) {
    for (unsigned int i = 0; i < to_delete; i++) {
        delete list_of_covariates.back();
        list_of_covariates.pop_back();
    }
    n_covariates -= to_delete;
    time_variant_covariates = time_variant_covariates.head(n_covariates);
}

/*
    Function to create a design matrix from the CovariateList according to the model orders
*/
arma::mat CovariateList::create_design_matrix(const arma::umat& orders, Neighborhood * W_covariates) const
{
    unsigned int n_params = arma::accu(orders);
    arma::mat design_matrix(n_obs * dim, n_params);
    unsigned int col_index = 0;

    unsigned int start_row = 0;
    unsigned int end_row = 0;
    arma::vec cov_values;

    for(unsigned int k = 0; k < n_covariates; k++)
    {
        arma::uvec orders_lag = arma::find( orders.col(k) );
        for(unsigned int o : orders_lag)
        {
            for(unsigned int t = 0; t < n_obs; t++)
            {
                start_row = t * dim;
                end_row = (t + 1) * dim - 1;
                cov_values = get_values_at(k, t);
                design_matrix.submat(start_row, col_index, end_row, col_index) = W_covariates->multiply_with_x(o, cov_values);
            }
            col_index++;
        }
    }
    return design_matrix;
}


/*
    Function to create a design matrix from some covariates in the CovariateList according to the model orders

    Parameters:
    * covariate_indices: Indices of the covariates to be included in the design matrix
    * orders: Orders of the covariates as matrix. Same number of columns as the length of covariate_indices
    * W_covariates: Pointer to Neighborhood object containing the spatial weights for covariates
*/
arma::mat CovariateList::create_design_matrix(const arma::uvec &covariate_indices, const arma::umat &orders, Neighborhood * W_covariates) const
{
    unsigned int n_params = arma::accu(orders);
    arma::mat design_matrix(n_obs * dim, n_params);
    unsigned int col_index = 0;

    unsigned int start_row = 0;
    unsigned int end_row = 0;
    arma::vec cov_values;

    for(unsigned int k = 0; k < covariate_indices.n_elem; k++)
    {
        arma::uvec orders_lag = arma::find( orders.col(k) );
        for(unsigned int o : orders_lag)
        {
            for(unsigned int t = 0; t < n_obs; t++)
            {
                start_row = t * dim;
                end_row = (t + 1) * dim - 1;
                cov_values = get_values_at(covariate_indices(k), t);
                design_matrix.submat(start_row, col_index, end_row, col_index) = W_covariates->multiply_with_x(o, cov_values);
            }
            col_index++;
        }
    }
    return design_matrix;
}





