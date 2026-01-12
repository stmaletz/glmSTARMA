/* 
-----------------------------------------------------------------------------
    File: neighbors.cpp
    Purpose: Implementations for Neighborhood classes
    Author: Steffen Maletz
    Last modified: 2025-12-07
-----------------------------------------------------------------------------
*/

#include "glmstarma.h"


/*
    Functions for SparseNeighborhood class
*/

arma::mat SparseNeighborhood::multiply_with_x(const int &order, const arma::mat &x) const
{
    if(order < max_order && !is_empty){
        return wlist(order) * x;
    } else if(is_empty) {
        return x;
    } else {
        return arma::zeros(dim, x.n_cols);
    }
}

void SparseNeighborhood::set_parameter_matrices(const arma::mat &param) {
    if(!is_empty)
    {
        if(param_matrices.n_elem == 0){
            param_matrices = arma::field<arma::sp_mat>(param.n_cols);
        }
        for(unsigned int i = 0; i < param.n_cols; i++)
        {
            param_matrices(i) = arma::sp_mat(dim, dim);
            arma::uvec indices = arma::find(param.col(i) != 0);
            for(unsigned int l : indices){
                param_matrices(i) += param(l, i) * wlist(l);
            }
        }
    } else {
        this->param = arma::vec(param.n_cols);
        if(this->param.n_elem > 0)
        {
            this->param = param.row(0).t();
        }
    }
}

arma::mat SparseNeighborhood::multiply_param_with_x(const int &index, const arma::mat &x) const 
{
    if(is_empty){
        return param(index) * x;
    } else {
        return param_matrices(index) * x;
    }
}


/*
    Functions for DenseNeighborhood class
*/

arma::mat DenseNeighborhood::multiply_with_x(const int &order, const arma::mat &x) const 
{
    if(order < max_order)
    {
        return wlist.slice(order) * x;
    } else {
        return arma::zeros(dim, x.n_cols);
    }
}


void DenseNeighborhood::set_parameter_matrices(const arma::mat &param) {
    if(param_matrices.n_slices == 0){
        param_matrices = arma::cube(dim, dim, param.n_cols);
    } else {
        param_matrices.zeros();
    }
   for(unsigned int i = 0; i < param.n_cols; i++)
   {
        arma::uvec indices = arma::find(param.col(i) != 0);
        for(unsigned int l : indices){
            param_matrices.slice(i) += param(l, i) * wlist.slice(l);
        }
    }
}

arma::mat DenseNeighborhood::multiply_param_with_x(const int &index, const arma::mat &x) const 
{
    return param_matrices.slice(index) * x;
}



/*
    Helper functions for factory method
*/

std::size_t Neighborhood::count_nonzero(const SEXP obj)
{
    if (Rf_isS4(obj)) {
        // Sparse matrix (dgCMatrix) is an S4 object
        Rcpp::S4 mat(obj);
        if (!mat.is("dgCMatrix")) {
            Rcpp::stop("Sparse matrices must be of class 'dgCMatrix'.");
        }
        return Rcpp::as<Rcpp::NumericVector>(mat.slot("x")).size();
    }
    arma::mat M = Rcpp::as<arma::mat>(obj);
    return arma::accu(M != 0);
}


unsigned int Neighborhood::get_matrix_dim(const SEXP obj)
{
    if (Rf_isS4(obj)) {
        // dgCMatrix is an S4 object
        Rcpp::S4 mat(obj);
        if (!mat.is("dgCMatrix")) {
            Rcpp::stop("Sparse matrices must be of class 'dgCMatrix'.");
        }
        Rcpp::IntegerVector Dim = mat.slot("Dim");
        unsigned int res = Dim[0];
        return res;
    }
    Rcpp::NumericMatrix M(obj);
    return M.nrow();
}


/*
    Factory method to create Neighborhood objects
*/

Neighborhood* Neighborhood::create(const Rcpp::List &wlist, const double &sparsity_threshold, const bool &use_sparsity) {
    if(!use_sparsity){
        return new DenseNeighborhood(wlist);
    } else {
        const unsigned int K   = wlist.size();
        const unsigned int dim = get_matrix_dim(wlist[0]);
        std::size_t non_zero_total = 0;
        for (int i = 0; i < K; ++i) {
            non_zero_total += count_nonzero(wlist[i]);
        }
        const double total_entries = static_cast<double>(K) * dim * dim;
        const double density = static_cast<double>(non_zero_total) / total_entries;
        if(density <= sparsity_threshold){
            return new SparseNeighborhood(wlist);
        } else {
            return new DenseNeighborhood(wlist);
        }
    }
}



Neighborhood* Neighborhood::create_neighbor(const unsigned int &param_count, const Rcpp::Nullable<Rcpp::List>& wlist, Neighborhood* fallback_clone, const unsigned int &dim, const double &sparsity_threshold, const bool &use_sparsity)
{
    if (param_count > 0 && wlist.isNotNull()) 
    {
        Rcpp::List wlist_ = Rcpp::as<Rcpp::List>(wlist);
        return Neighborhood::create(wlist_, sparsity_threshold, use_sparsity);
    } else if (param_count > 0) 
    {
        return fallback_clone->clone();
    } else {
        return new SparseNeighborhood(dim);
    }
}



