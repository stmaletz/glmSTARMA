/* 
-----------------------------------------------------------------------------
    File: covariates.h
    Purpose: Declarations for Covariate and CovariateList classes
    Author: Steffen Maletz
    Last modified: 2025-12-12
-----------------------------------------------------------------------------
*/


#ifndef COVARIATES_H
#define COVARIATES_H

/*
    Abstract base class for different types of covariates
*/
class Covariate{
    protected:
    int burn_in;

    public:
    bool is_time_constant;
    arma::mat values;
    Covariate() {};
    Covariate(int burn_in, bool is_time_constant) {
        this->burn_in = burn_in;
        this->is_time_constant = is_time_constant;
    };
    Covariate(arma::mat covariate_values, int burn_in, bool is_time_constant) {
        if(covariate_values.is_rowvec()){
            values = covariate_values.t();
        } else {
            values = covariate_values;
        }
        this->burn_in = burn_in;
        this->is_time_constant = is_time_constant;
    };
    Covariate(const Covariate* to_clone) : burn_in(to_clone->burn_in), is_time_constant(to_clone->is_time_constant), values(to_clone->values) {};
    virtual ~Covariate() = default;
    virtual Covariate* clone() = 0;
    virtual arma::vec get_values_at(int t) = 0;
    virtual void append_values(arma::mat vals) = 0;
    virtual void delete_last_values(int to_delete) = 0;
};

/*
    Class for time-constant covariates, i.e. covariates that do not change over time but only over space
    Only the values at each location are stored
*/
class TimeConstantCovariate : public Covariate {
    public:
    TimeConstantCovariate(arma::mat covariate_values, int burn_in) : Covariate(covariate_values, burn_in, true) {};
    TimeConstantCovariate(const TimeConstantCovariate* to_clone) : Covariate(to_clone) {};
    ~TimeConstantCovariate() = default;
    arma::vec get_values_at(int t){
        if(t >= burn_in){
            return values.col(0);
        } else {
            return arma::zeros(values.n_rows);
        }
    }
    TimeConstantCovariate* clone(){
        return new TimeConstantCovariate(this);
    }
    void append_values(arma::mat vals){};
    void delete_last_values(int to_delete){};
};

/*
    Class for spatial-constant covariates, i.e. covariates that do not change over space but only over time
    Only the values at each time point are stored
*/
class SpatialConstantCovariate : public Covariate {
    private:
    unsigned int dim;
    arma::vec returned_value;
    public:
    SpatialConstantCovariate(arma::mat covariate_values, int dim, int burn_in) : Covariate(covariate_values, burn_in, false), dim(dim) {
        returned_value = arma::vec(dim, arma::fill::zeros);
    }
    SpatialConstantCovariate(const SpatialConstantCovariate* to_clone) : Covariate(to_clone), dim(to_clone->dim) {
        returned_value = arma::vec(dim, arma::fill::zeros);
    }
    ~SpatialConstantCovariate() = default;
    arma::vec get_values_at(int t){
        if(t >= burn_in){
            if(t - burn_in >= values.n_elem){
                returned_value.fill(values(values.n_elem - 1));
            } else {
                returned_value.fill(values(t - burn_in));
            }
        } else {
            returned_value.zeros();
        }
        return returned_value;
    }
    void append_values(arma::mat vals){
        values = arma::join_cols(values, vals);
    }
    SpatialConstantCovariate* clone(){
        return new SpatialConstantCovariate(this);
    }
    void delete_last_values(int to_delete){
        values = values.head_rows(values.n_elem - to_delete);
    }
};

/*
    Class for default covariates, i.e. covariates that change over both space and time
    These are submitted as a matrix with same dimensions as the time series that is modeled or simulated
*/
class DefaultCovariate : public Covariate {
    public:
    DefaultCovariate(arma::mat covariate_values, int burn_in) : Covariate(covariate_values, burn_in, false) {};
    DefaultCovariate(const DefaultCovariate* to_clone) : Covariate(to_clone) {};
    ~DefaultCovariate() = default;
    arma::vec get_values_at(int t){
        if(t >= burn_in){
            if(t - burn_in >= values.n_cols){
                return values.col(values.n_cols - 1);
            } else {
                return values.col(t - burn_in);
            }
        } else {
            return arma::zeros(values.n_rows);
        }
    }
    DefaultCovariate* clone(){
        return new DefaultCovariate(this);
    }
    void append_values(arma::mat vals){
        values = arma::join_rows(values, vals);
    }
    void delete_last_values(int to_delete){
        values = values.head_cols(values.n_cols - to_delete);
    }
};

/*
    Class to manage a list of covariates
*/
class CovariateList {
    private:
    unsigned int n_obs;
    const unsigned int dim;
    arma::uvec time_variant_covariates;
    const int burn_in;
    friend class CovariateList;
    public:
    std::vector<Covariate*> list_of_covariates;
    unsigned int n_covariates;
    std::vector<std::string> covariate_names;
    ~CovariateList()
    {
        for(Covariate* cov : list_of_covariates)
        {
            delete cov;
        }
    };
    CovariateList(const Rcpp::List &covariates, const unsigned int &n_obs, const unsigned int &dim, const int &burn_in, const unsigned int &shift_interventions);
    CovariateList(CovariateList* to_clone);
    CovariateList(const unsigned int &n_obs, const unsigned int &dim);
    arma::vec get_values_at(const int &k, const int &t) const;
    bool has_time_variant_covariates() const {
        return arma::accu(time_variant_covariates) > 0;
    }
    // Function to append additional covariate values for prediction
    void append_covariate_values(const CovariateList &other)
    {
        for(unsigned int i = 0; i < n_covariates; i++)
        {
            if(time_variant_covariates(i) != 0)
            {
                auto it = std::find(other.covariate_names.begin(), other.covariate_names.end(), covariate_names[i]);
                if(it == other.covariate_names.end())
                {
                    Rcpp::warning("Covariate " + covariate_names[i] + " was not found in the new covariate values.");
                } else {
                    int index = std::distance(other.covariate_names.begin(), it);
                    list_of_covariates[i]->append_values(other.list_of_covariates[index]->values);
                }
            }
        }
        n_obs += other.n_obs;
    }
    CovariateList* clone();
    void delete_last_covariates(unsigned int &to_delete);
    void add_covariate(Covariate* cov);
    arma::mat create_design_matrix(const arma::umat &orders, Neighborhood * W_covariates) const;
    arma::mat create_design_matrix(const arma::uvec &covariate_indices, const arma::umat &orders, Neighborhood * W_covariates) const;
};


#endif
