/* 
-----------------------------------------------------------------------------
    File: fitting.h
    Purpose: Declaration of FittingObject class
    Author: Steffen Maletz
    Last modified: 2025-12-06
-----------------------------------------------------------------------------
*/


#ifndef FITTING_H
#define FITTING_H


/*
    Abstract base class for different optimization methods used for fitting the model
    Defines the interface for fitting functions and contains common attributes
    Attributes:
      * method: Character vector indicating the optimization method used
      * algorithm: Character vector indicating the specific algorithm used within the optimization method
      * fcounter: Counter for number of function evaluations during fitting
      * grcounter: Counter for number of gradient evaluations during fitting
      * info_counter: Counter for number of information matrix evaluations during fitting
      * fitting_time: Time taken for the fitting process
      * start_vector: Starting values used for optimization
      * converged: Boolean indicating whether optimization converged
      * convergence_code: Integer code indicating convergence status
      * message: Character vector with message from optimization routine
    Methods:
      * fit: Pure virtual method to perform the fitting given starting values
      * get_algorithm_info: Returns a list with method and algorithm information
      * get_fitting_info: Returns a list with fitting statistics
      * create: Static factory method to create FittingObject instances based on control settings
*/

class FittingObject{
    protected:
    const Rcpp::CharacterVector method; // Set in Constructor
    const Rcpp::CharacterVector algorithm; // Set in Konstruktor
    unsigned int fcounter; // set in fit
    unsigned int grcounter; // set in fit
    unsigned int info_counter; // set in fit
    long long int fitting_time; // set in fit
    arma::vec start_vector; // set in fit
    bool converged; // set in fit
    int convergence_code;
    Rcpp::CharacterVector message; // set in fit
    public:
    FittingObject(const char* methode, const char* alg) : method(methode), algorithm(alg)
    {
        fcounter = 0;
        grcounter = 0;
        info_counter = 0;
    };
    FittingObject(FittingObject* to_clone) : method(to_clone->method), algorithm(to_clone->algorithm)
    {
        fcounter = to_clone->fcounter;
        grcounter = to_clone->grcounter;
        info_counter = to_clone->info_counter;
        fitting_time = to_clone->fitting_time;
        start_vector = to_clone->start_vector;
        converged = to_clone->converged;
        convergence_code = to_clone->convergence_code;
        message = to_clone->message;
    };
    virtual FittingObject * clone() = 0;
    virtual ~FittingObject() = default;
    virtual arma::vec fit(arma::vec start_value) = 0;
    Rcpp::List get_algorithm_info() const 
    {
        return Rcpp::List::create(Rcpp::Named("method") = method, Rcpp::Named("algorithm") = algorithm);
    };
    Rcpp::List get_fitting_info() const
    {
        return Rcpp::List::create(Rcpp::Named("start") = start_vector, Rcpp::Named("fncount") = fcounter, Rcpp::Named("grcount") = grcounter,
        Rcpp::Named("hecount") = info_counter, Rcpp::Named("fitting_time") = fitting_time,
        Rcpp::Named("convergence") = Rcpp::LogicalVector(converged), Rcpp::Named("message") = message);
    };
    const arma::vec get_start() const 
    {
        return start_vector;
    }
    const long long int get_fitting_time() const 
    {
        return fitting_time;
    }
    const int get_convergence_code() const 
    {
        return convergence_code;
    }
    const bool is_converged() const 
    {
        return converged;
    }
    const unsigned int get_fncount() const 
    {
        return fcounter;
    }
    const unsigned int get_grcount() const 
    {
        return grcounter;
    }
    const unsigned int get_hecount() const 
    {
        return info_counter;
    }
    static FittingObject * create(const arma::mat &ts, Design * model_design, Orders * orders, CovariateList * covariates, Neighborhood * W_ma, Family * fam, const Rcpp::List &control);
};


#endif
