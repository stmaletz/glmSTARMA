/* 
-----------------------------------------------------------------------------
    File: fitting.cpp
    Purpose: Implementation of FittingObject classes
    Author: Steffen Maletz
    Last modified: 2025-12-07
-----------------------------------------------------------------------------
*/

// [[Rcpp::depends(roptim)]]
// [[Rcpp::depends(nloptr)]]
#include "glmstarma.h"
#include <nloptrAPI.h>
#include <roptim.h>
#include <chrono>
using namespace roptim;


/*
    Help class for optimization with the nloptr library
*/

class optim_model : public Functor
{
    public:
    arma::vec ts_vec;
    const Orders * orders;
    Design * design;
    Family * family;
    Neighborhood * W_ma;

    arma::vec score;
    arma::mat information;

    optim_model(const arma::mat &ts, const Orders * model_orders, Design * design, Family * fam, Neighborhood * W_ma) : ts_vec(arma::vectorise(ts.tail_cols(model_orders->n_obs_effective))), orders(model_orders), design(design), family(fam), W_ma(W_ma), score(model_orders->n_param), information(model_orders->n_param, model_orders->n_param) {};
    ~optim_model()
    {
        orders = nullptr;
        design = nullptr;
        family = nullptr;
        W_ma = nullptr;
    }

    double operator()(const arma::vec &x) override {
        return -Model::log_likelihood(x, score, information, ts_vec, design, orders, W_ma, family, false, false);
    }
    void Gradient(const arma::vec &x, arma::vec &gr) override {
        Model::log_likelihood(x, gr, information, ts_vec, design, orders, W_ma, family, true, false);
        gr = -gr;
    }
    void Hessian(const arma::vec &x, arma::mat &he) override {
        Model::log_likelihood(x, score, he, ts_vec, design, orders, W_ma, family, true, true);
        he = -he;
    }
};

/*
    Class to wrap nloptr optimization in a FittingObject
*/

class nloptr : public FittingObject {
    nlopt_opt opt;
    optim_model to_optimize;
    arma::uvec * constrain_param;
    arma::uvec * constrain_param2;
    public:
    inline static int fcount = 0;
    inline static int gcount = 0;
    inline static int ccount = 0;
    ~nloptr(){
        nlopt_destroy(opt);
        delete constrain_param;
        delete constrain_param2;
    }
    nloptr(const arma::mat &ts, const Orders * model_orders, Design * design, Family * fam, Neighborhood * W_ma, const Rcpp::List control) : FittingObject("nloptr", "SLSQP"), to_optimize(ts, model_orders, design, fam, W_ma)
    {
        Rcpp::LogicalVector constrained = control["constrained"]; // Default TRUE
        Rcpp::NumericVector constraint_tol = control["constraint_tol"]; // Default 1e-8
        Rcpp::CharacterVector constrain_method = control["constrain_method"];
        Rcpp::NumericVector xtol = control["reltol"]; // Default sqrt(.Machine$double.eps)
        Rcpp::IntegerVector maxit = control["maxit"];

        opt = nlopt_create(NLOPT_LD_SLSQP, model_orders->n_param);
        nlopt_set_min_objective(opt, neg_log_likelihood, &to_optimize);

        if(fam->only_positive_parameters){
            double* lb = new double[model_orders->n_param];
            for(int i = 0; i < model_orders->n_param; i++){
                lb[i] = 0.0;
            }
            nlopt_set_lower_bounds(opt, lb);
        }
        
        constrain_param = new arma::uvec(model_orders->n_param);
        constrain_param2 = new arma::uvec(model_orders->n_param);
        if(model_orders->n_param_ma > 0){
            constrain_param->subvec(0, model_orders->n_param_intercept - 1).fill(1);
            constrain_param2->subvec(0, model_orders->n_param_intercept - 1).fill(1);
        }
        if(model_orders->n_param_ar > 0)
        {
            constrain_param->subvec(model_orders->n_param_ma + model_orders->n_param_intercept, model_orders->n_param_ma + model_orders->n_param_intercept + model_orders->n_param_ar - 1).fill(1);
        }

        if(constrained[0] & (fam->only_positive_parameters || (constrain_method[0] == "sum_of_absolutes"))){
            nlopt_add_inequality_constraint(opt, constraint_sum_of_absolutes, constrain_param, constraint_tol[0]);
        } else if(constrained[0] & (constrain_method[0] == "absolute_sum")){
            if(!fam->only_positive_parameters){
                // All autoregressive and moving average parameters in [−1, 1]
                // If non-negative parameters are needed, a constraint is already set above
                double* lb = new double[model_orders->n_param];
                double* ub = new double[model_orders->n_param];
                for(int i = 0; i < model_orders->n_param; i++){
                    lb[i] = constrain_param->at(i) > 0 ? -1.0 : -HUGE_VAL;
                    ub[i] = constrain_param->at(i) > 0 ? 1.0 : HUGE_VAL;
                }
                nlopt_set_lower_bounds(opt, lb);
                nlopt_set_upper_bounds(opt, ub);
            }
            nlopt_add_inequality_constraint(opt, constraint_absolute_sum_lower, constrain_param, constraint_tol[0]);
            nlopt_add_inequality_constraint(opt, constraint_absolute_sum_upper, constrain_param, constraint_tol[0]);
        } else if(constrained[0] & (constrain_method[0] == "soft") & (fam->link_name == "softplus" || fam->link_name == "softclipping")){
            // Only available for softplus and softclipping link functions
            nlopt_add_inequality_constraint(opt, constraint_soft, constrain_param, constraint_tol[0]);
            nlopt_add_inequality_constraint(opt, constraint_sum_of_absolutes, constrain_param2, constraint_tol[0]);
        }

        // Set stopping criteria
        nlopt_set_xtol_rel(opt, xtol[0]);
        nlopt_set_maxeval(opt, maxit[0]);
    }
    nloptr(nloptr* to_clone) : FittingObject(to_clone), to_optimize(to_clone->to_optimize)
    {
        opt = nlopt_copy(to_clone->opt);
        constrain_param = new arma::uvec(*to_clone->constrain_param);
    }
    FittingObject * clone() override {
        return new nloptr(this);
    }

    arma::vec fit(arma::vec start_value) override 
    {
        fcount = 0;
        gcount = 0;
        ccount = 0;
        double* start = new double[start_value.n_elem];
        for(unsigned int i = 0; i < start_value.n_elem; i++){
            start[i] = start_value(i);
        }
        double min_f = 0.0;
        start_vector = arma::vec(start_value);

        auto start_time = std::chrono::high_resolution_clock::now();
        nlopt_result res = nlopt_optimize(opt, start, &min_f);
        auto stop_time = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time);
        fitting_time = duration.count();

        fcounter = fcount;
        grcounter = gcount; 
        info_counter = ccount;
        convergence_code = res;
        converged = (res > 0 && res <= 4);
        message = Rcpp::CharacterVector("No message available");

        arma::vec res_vec(start_value.n_elem);
        for(int i = 0; i < start_value.n_elem; i++)
        {
            res_vec(i) = start[i];
        }
        delete[] start;
        return res_vec;
    }

    static double neg_log_likelihood(unsigned int n, const double *x, double *grad, void *func_data){
        fcount++;
        optim_model * model = (optim_model *) func_data;
        arma::vec param(x, n);
        double value;
        arma::vec score(n);
        if(grad)
        {
            gcount++;
            value = -Model::log_likelihood(param, score, model->information, model->ts_vec, model->design, model->orders, model->W_ma, model->family, true, false);
            for(unsigned int i = 0; i < n; i++){
                grad[i] = -score(i);
            }
        }
        else
        {
            value = -Model::log_likelihood(param, score, model->information, model->ts_vec, model->design, model->orders, model->W_ma, model->family, false, false);
        }
        return value;
    }

    static double constraint_sum_of_absolutes(unsigned int n, const double *x, double *grad, void *func_data){
        // sum_i |x_i| <= 1 (smooth approximation)
        ccount++;
        arma::uvec* y = static_cast<arma::uvec*>(func_data);
        arma::uvec temp;
        if(y)
        {
            temp = *y;
        }
        if(grad){
            for(int i = 0; i < n; i++)
            {
                grad[i] = temp(i) > 0 ? x[i] / std::sqrt( x[i] * x[i] + 1e-5) : 0.0;
            }
        }
        double value = 0.0;
        for(int i = 0; i < n; i++)
        {
            value += temp(i) * std::sqrt( x[i] * x[i] + 1e-5);
        }
        return (value - 1.0);
    }

    static double constraint_absolute_sum_lower(unsigned int n, const double *x, double *grad, void *func_data){
        // sum_i x_i >= -1
        ccount++;
        arma::uvec* y = static_cast<arma::uvec*>(func_data);
        arma::uvec temp;
        if(y)
        {
            temp = *y;
        }
        double value = 0.0;
        for(int i = 0; i < n; i++)
        {
            value += temp(i) * x[i];
        }
        if(grad){
            for(int i = 0; i < n; i++)
            {
                grad[i] = temp(i) > 0 ? -1.0 : 0.0;
            }
        }
        return (-value - 1.0);
    }

    static double constraint_absolute_sum_upper(unsigned int n, const double *x, double *grad, void *func_data){
        // sum_i x_i <= 1
        ccount++;
        arma::uvec* y = static_cast<arma::uvec*>(func_data);
        arma::uvec temp;
        if(y)
        {
            temp = *y;
        }
        double value = 0.0;
        for(int i = 0; i < n; i++)
        {
            value += temp(i) * x[i];
        }
        if(grad){
            for(int i = 0; i < n; i++)
            {
                grad[i] = temp(i) > 0 ? 1.0 : 0.0;
            }
        }
        return (value - 1.0);
    }

    // Constraint for softplus and softclipping
    static double constraint_soft(unsigned int n, const double *x, double *grad, void *func_data){
        // sum_i max(0, x_i) <= 1  with smooth Approximation
        ccount++;
        arma::uvec* y = static_cast<arma::uvec*>(func_data);
        arma::uvec mask;
        if(y) mask = *y;
        const double eps = 1e-8; // smoothing parameter
        double value = 0.0;

        // compute value = sum mask(i) * (x_i + sqrt(x_i^2 + eps)) / 2
        for(unsigned int i = 0; i < n; ++i){
            if(mask(i) > 0){
                double xi = x[i];
                double s = std::sqrt(xi*xi + eps);
                double term = 0.5 * (xi + s);
                value += term;
            }
        }

        if(grad){
            // derivative of 0.5*(x + sqrt(x^2+eps)) w.r.t x is 0.5*(1 + x / sqrt(x^2+eps))
            for(unsigned int i = 0; i < n; ++i){
                if(mask(i) > 0){
                    double xi = x[i];
                    double s = std::sqrt(xi*xi + eps);
                    grad[i] = 0.5 * (1.0 + xi / s);
                } else {
                    grad[i] = 0.0;
                }
            }
        }
        return (value - 1.0);
    }
};


/*
    Optimization with optim
    Uses the C implementation of the optim function in R via roptim
    * For models with exclusively positive parameters, L-BFGS-B is used (otherwise unconstrained)
    * For models with unconstrained parameters, BFGS is used (unconstrained)
*/

class optim : public FittingObject {
    Roptim<optim_model> opt;
    optim_model to_optimize;
    public:
    optim(const arma::mat &ts, const Orders * model_orders, Design * design, Family * fam, Neighborhood * W_ma, const Rcpp::List control) : FittingObject("optim", fam->only_positive_parameters ? "L-BFGS-B" : "BFGS"), to_optimize(ts, model_orders, design, fam, W_ma)
    {
        if(fam->only_positive_parameters)
        {
            opt.set_method("L-BFGS-B");
            opt.set_lower(arma::zeros(model_orders->n_param));
            Rcpp::IntegerVector lmm = control["lmm"];
            Rcpp::NumericVector factr = control["factr"];
            Rcpp::NumericVector pgtol = control["pgtol"];

            opt.control.lmm = lmm[0];
            opt.control.factr = factr[0];
            opt.control.pgtol = pgtol[0];
        } 
        else 
        {
            opt.set_method("BFGS");
            Rcpp::NumericVector abstol = control["abstol"];
            Rcpp::NumericVector reltol = control["reltol"];

            opt.control.abstol = abstol[0];
            opt.control.reltol = reltol[0];
        }
        Rcpp::IntegerVector trace = control["trace"];
        Rcpp::NumericVector fnscale = control["fnscale"];
        Rcpp::IntegerVector maxit = control["maxit"];

        opt.control.trace = trace[0];
        opt.control.fnscale = fnscale[0];
        opt.control.maxit = maxit[0];
        opt.set_hessian(false);
    }
    optim(const optim &to_clone) : FittingObject(to_clone), opt(to_clone.opt), to_optimize(to_clone.to_optimize) {};
    FittingObject * clone() override
    {
        return new optim(*this);
    }

    arma::vec fit(arma::vec start_value) override
    {
        start_vector = arma::vec(start_value);
        arma::vec start(start_value);

        auto start_time = std::chrono::high_resolution_clock::now();
        opt.minimize(to_optimize, start);
        auto stop_time = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time);
        fitting_time = duration.count();

        converged = (opt.convergence() == 0); // 0 means convergence in roptim
        convergence_code = opt.convergence();
        message = opt.message();
        fcounter = opt.fncount();
        grcounter = opt.grcount();
        info_counter = 0;

        return opt.par();
    };
};





/*
    Factory method for FittingObject
    * nloptr -> (constrained) optimization with the nloptr library
    * optim -> optimization with the C implementation of the optim function in R via roptim (unconstrained/box-constrained)
*/

FittingObject * FittingObject::create(const arma::mat &ts, Design * model_design, Orders * orders, CovariateList * covariates, Neighborhood * W_ma, Family * fam, const Rcpp::List &control)
{
    std::string method = Rcpp::as<std::string>(control["method"]);
    if(method == "nloptr"){
        return new nloptr(ts, orders, model_design, fam, W_ma, control);
    } 
    else {
        return new optim(ts, orders, model_design, fam, W_ma, control);
    }
}