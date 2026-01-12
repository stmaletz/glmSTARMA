/* 
-----------------------------------------------------------------------------
    File: neighbors.h
    Purpose: Declarations for Neighborhood classes, to store spatial weight matrices
    Author: Steffen Maletz
    Last modified: 2025-12-07
-----------------------------------------------------------------------------
*/


#ifndef NEIGHBORS_H
#define NEIGHBORS_H

/*
    Abstract base class for different types of Neighborhood structures
    Defines the interface for multiplying vectors with spatial weight matrices
    Attributes:
      * max_order: Maximum order of neighborhood matrices stored
      * is_empty: Boolean indicating whether no neighborhood structure is used
    Methods:
      * multiply_with_x: Multiplies the spatial weight matrix of given order with matrix x
      * multiply_param_with_x: Multiplies the spatial weight matrix corresponding to given parameter index with matrix x
      * set_parameter_matrices: Sets the spatial weight matrices according to given parameter matrix
      * clone: Creates a deep copy of the Neighborhood object
      * create: Static factory method to create Neighborhood instances based on R list input
      * create_neighbor: Static factory method to create individual Neighborhood instances for AR, MA, or covariate effects
*/


class Neighborhood {
    private:
    static std::size_t count_nonzero(const SEXP obj);
    protected:
    const unsigned int max_order;
    const bool is_empty;
    public:
    virtual ~Neighborhood() = default;
    Neighborhood(const unsigned int &max_order, const bool &empty) : max_order(max_order), is_empty(empty) {};
    virtual arma::mat multiply_with_x(const int &order, const arma::mat &x) const = 0;
    virtual arma::mat multiply_param_with_x(const int &index, const arma::mat &x) const = 0;
    virtual void set_parameter_matrices(const arma::mat &param) = 0;
    virtual Neighborhood* clone() = 0;
    static Neighborhood* create(const Rcpp::List &wlist, const double &sparsity_threshold, const bool &use_sparsity);
    static Neighborhood* create_neighbor(const unsigned int &param_count, const Rcpp::Nullable<Rcpp::List>& wlist, Neighborhood* fallback_clone, const unsigned int &dim, const double &sparsity_threshold, const bool &use_sparsity);
    static unsigned int get_matrix_dim(const SEXP obj);
};

/*
    Dense Neighborhood class
    Stores the spatial weight matrices as dense matrices (arma::mat)
*/

class DenseNeighborhood : public Neighborhood {
    private:
    arma::cube wlist;
    arma::cube param_matrices;
    unsigned int dim;

    static arma::mat as_dense_mat(const SEXP obj)
    {
        if (Rf_isS4(obj)) {
            // Sparse matrix (dgCMatrix) is an S4 object
            Rcpp::S4 mat(obj);
            if (!mat.is("dgCMatrix")) {
                    Rcpp::stop("Sparse matrices must be of class 'dgCMatrix'.");
            }
            return Rcpp::as<arma::mat>(mat);
        }
        return Rcpp::as<arma::mat>(obj);
    };
    public:
    DenseNeighborhood(const Rcpp::List &W) : Neighborhood(W.size(), false)
    {
        arma::mat first = as_dense_mat(W[0]);
        dim = first.n_rows;
        wlist = arma::cube(dim, dim, W.size());
        wlist.slice(0) = first;
        for (int i = 1; i < W.size(); i++) {
            wlist.slice(i) = as_dense_mat(W[i]);
        }
    };
    DenseNeighborhood(DenseNeighborhood* to_clone) : Neighborhood(to_clone->max_order, to_clone->is_empty), wlist(to_clone->wlist), param_matrices(to_clone->param_matrices), dim(to_clone->dim) {};
    DenseNeighborhood(const unsigned int &dim) : Neighborhood(0, true), dim(dim) {};
    Neighborhood* clone() override
    {
        return new DenseNeighborhood(this);
    };
    ~DenseNeighborhood() = default;
    arma::mat multiply_with_x(const int &order, const arma::mat &x) const override;
    void set_parameter_matrices(const arma::mat &param) override;
    arma::mat multiply_param_with_x(const int &index, const arma::mat &x) const override;
};


/*
    Sparse Neighborhood class
    Stores the spatial weight matrices as sparse matrices (arma::sp_mat)
*/


class SparseNeighborhood : public Neighborhood {
    private:
    arma::field<arma::sp_mat> wlist;
    arma::field<arma::sp_mat> param_matrices;
    unsigned int dim;
    arma::vec param;

    // Helper function to allow both dgCMatrix and dense matrices as input
    static arma::sp_mat as_sp_mat(const SEXP obj)
    {
        if (Rf_isS4(obj)) 
        {
            // dgCMatrix is an S4 object
            Rcpp::S4 mat(obj);
            if (!mat.is("dgCMatrix")) {
                Rcpp::stop("Sparse matrices must be of class 'dgCMatrix'.");
            }
            return Rcpp::as<arma::sp_mat>(mat);
        }
        // Otherwise the matrix is dense
        arma::mat dense = Rcpp::as<arma::mat>(obj);
        return arma::sp_mat(dense);
    };

    public:
    SparseNeighborhood(const Rcpp::List &W) : Neighborhood(W.size(), false)
    {
        wlist = arma::field<arma::sp_mat>(W.size());
        arma::sp_mat first = as_sp_mat(W[0]);
        dim = first.n_rows;
        wlist(0) = first;
        for(int i = 1; i < W.size(); i++)
        {
            wlist(i) = as_sp_mat(W[i]);
        }
    };
    SparseNeighborhood(SparseNeighborhood* to_clone) : Neighborhood(to_clone->max_order, to_clone->is_empty), wlist(to_clone->wlist), param_matrices(to_clone->param_matrices), dim(to_clone->dim), param(to_clone->param) {};
    SparseNeighborhood(const unsigned int &dim) : Neighborhood(0, true), dim(dim) {};
    Neighborhood* clone() override
    {
        return new SparseNeighborhood(this);
    };
    ~SparseNeighborhood() = default;
    arma::mat multiply_with_x(const int &order, const arma::mat &x) const override;
    void set_parameter_matrices(const arma::mat &param) override;
    arma::mat multiply_param_with_x(const int &index, const arma::mat &x) const override;
};

#endif
