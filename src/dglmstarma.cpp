/* 
-----------------------------------------------------------------------------
    File: dglmstarma.cpp
    Purpose: Implementation dglmstarma fitting function
    Author: Steffen Maletz
    Last modified: 2026-01-12
-----------------------------------------------------------------------------
*/

#include "glmstarma.h"


/*
    dglmstarma_cpp: Function called from R to fit a dglmstarma model
    Input:
    * ts: Time series data matrix
    * mean_model: List specifying the mean model
    * dispersion_model: List specifying the dispersion model
    * mean_family: List specifying the family for the mean model
    * dispersion_family: List specifying the family for the dispersion model
    * wlist: List of neighborhoods for the mean model
    * mean_covariates: List of covariates for the mean model
    * dispersion_covariates: List of covariates for the dispersion model
    * pseudo_observations: Type of pseudo observations to use
    * wlist_past_mean: List of neighborhoods for past mean values in the dispersion model
    * wlist_covariates: List of neighborhoods for covariates in the mean model
    * wlist_ar_dispersion: List of neighborhoods for autoregressive parameters in the dispersion model
    * wlist_past_mean_dispersion: List of neighborhoods for past mean values in the dispersion model
    * wlist_covariates_dispersion: List of neighborhoods for covariates in the dispersion model
    * control: List of control parameters for fitting
    Output:
    * List containing estimation results
*/


// [[Rcpp::export]]
Rcpp::List dglmstarma_cpp(const arma::mat &ts, const Rcpp::List &mean_model, const Rcpp::List &dispersion_model,
    const Rcpp::List &mean_family, const Rcpp::List &dispersion_family, const Rcpp::List &wlist,
    const Rcpp::List &mean_covariates, const Rcpp::List &dispersion_covariates,
    const Rcpp::CharacterVector &pseudo_observations, const Rcpp::Nullable<Rcpp::List> wlist_past_mean,
    const Rcpp::Nullable<Rcpp::List> wlist_covariates, const Rcpp::Nullable<Rcpp::List> wlist_ar_dispersion,
    const Rcpp::Nullable<Rcpp::List> wlist_past_mean_dispersion, const Rcpp::Nullable<Rcpp::List> wlist_covariates_dispersion,
    const Rcpp::List &control)
{
    // Control Arguments for model
    const bool use_sparsity = Rcpp::as<bool>(control["use_sparsity"]);
    const double sparsity_threshold = Rcpp::as<double>(control["sparsity_threshold"]);
    Rcpp::RObject init_method_param_mean = control["parameter_init"];
    Rcpp::RObject init_method_link = control["init_link"];
    Rcpp::RObject init_method_param_dispersion = control["parameter_init_dispersion"];
    Rcpp::RObject init_method_dispersion = control["init_dispersion"];

    // Control Arguments for iterated model fitting
    const bool print_progress = Rcpp::as<bool>(control["print_progress"]);
    const bool print_warnings = Rcpp::as<bool>(control["print_warnings"]);
    const bool use_backtracking = Rcpp::as<bool>(control["use_backtracking"]);

    const double alpha_shrink = Rcpp::as<double>(control["alpha_shrink"]);
    const double alpha_start = Rcpp::as<double>(control["alpha_start"]);
    const double min_alpha = Rcpp::as<double>(control["min_alpha"]);

    const double convergence_threshold = control["convergence_threshold"];
    const unsigned int max_fits = control["max_fits"];
    const bool drop_max_mean_lag = Rcpp::as<bool>(control["drop_max_mean_lag"]);
    const bool previous_param_as_start = Rcpp::as<bool>(control["previous_param_as_start"]);
    const bool use_fast_if_const_dispersion = Rcpp::as<bool>(control["use_fast_if_const_dispersion"]);
    // Control arguments for calculating pseudo observations
    const double lower_dispersion = control["lower_dispersion"];
    const double upper_dispersion = control["upper_dispersion"];
    const std::string pseudo_obs(pseudo_observations[0]);

    Rcpp::List control_mean = Rcpp::clone(control);
    control_mean["constrained"] = control["constrained_mean"];
    control_mean["constrain_method"] = control["constrain_method_mean"];
    Rcpp::List control_dispersion = Rcpp::clone(control);
    control_dispersion["constrained"] = control["constrained_dispersion"];
    control_dispersion["constrain_method"] = control["constrain_method_dispersion"];


    // Initalize model orders, covariates and weight matrices:
    Orders mean_orders(mean_model, ts.n_rows, ts.n_cols);
    const unsigned int dispersion_obs = drop_max_mean_lag ? mean_orders.n_obs_effective : ts.n_cols;
    Orders dispersion_orders(dispersion_model, ts.n_rows, dispersion_obs);

    CovariateList mean_covariates_list(mean_covariates, ts.n_cols, ts.n_rows, 0, 0);
    unsigned int shift_interventions = drop_max_mean_lag ? mean_orders.max_time_lag : 0;
    CovariateList dispersion_covariates_list(dispersion_covariates, dispersion_obs, ts.n_rows, 0, shift_interventions); // problematic if interventions are present
    const bool const_dispersion_model = (dispersion_orders.n_param_ar == 0) && (dispersion_orders.n_param_ma == 0) && (dispersion_orders.n_param_cov == 0);

    Family * mean_fam = Family::create(mean_family);
    Family * dispersion_fam = Family::create(dispersion_family);

    Neighborhood * W_ar_mean = Neighborhood::create(wlist, sparsity_threshold, use_sparsity);
    Neighborhood * W_ma_mean = Neighborhood::create_neighbor(mean_orders.n_param_ma, wlist_past_mean, W_ar_mean, ts.n_rows, sparsity_threshold, use_sparsity);
    Neighborhood * W_covariates_mean = Neighborhood::create_neighbor(mean_orders.n_param_cov, wlist_covariates, W_ar_mean, ts.n_rows, sparsity_threshold, use_sparsity);
    Neighborhood * W_ar_dispersion = nullptr;
    Neighborhood * W_ma_dispersion = nullptr;
    Neighborhood * W_covariates_dispersion = nullptr;
    if(!const_dispersion_model)
    {
        W_ar_dispersion = Neighborhood::create_neighbor(dispersion_orders.n_param_ar, wlist_ar_dispersion, W_ar_mean, ts.n_rows, sparsity_threshold, use_sparsity);
        W_ma_dispersion = Neighborhood::create_neighbor(dispersion_orders.n_param_ma, wlist_past_mean_dispersion, W_ar_mean, ts.n_rows, sparsity_threshold, use_sparsity);
        W_covariates_dispersion = Neighborhood::create_neighbor(dispersion_orders.n_param_cov, wlist_covariates_dispersion, W_ar_mean, ts.n_rows, sparsity_threshold, use_sparsity);
    }
    // Initialize parameters of mean model:
    arma::mat link_init = Model::init_link(init_method_link, mean_orders, ts, mean_fam);
    arma::vec start_param_mean = Model::init_param(init_method_param_mean, mean_orders, ts, mean_fam);
    arma::vec intercept_mean(mean_orders.n_param_intercept);
    arma::mat ar_params_mean(mean_orders.autoregressive_orders.n_rows, mean_orders.autoregressive_orders.n_cols);
    arma::mat ma_params_mean(mean_orders.moving_average_orders.n_rows, mean_orders.moving_average_orders.n_cols);
    arma::mat cov_params_mean(mean_orders.covariate_orders.n_rows, mean_orders.covariate_orders.n_cols);
    Parameter::create_param_matrices(start_param_mean, &mean_orders, intercept_mean, ar_params_mean, ma_params_mean, cov_params_mean);
    W_ma_mean->set_parameter_matrices(ma_params_mean);

    // Initialize variables for the dispersion model:
    arma::mat pseudo;
    arma::mat dispersion_init;
    arma::vec start_param_dispersion;
    arma::vec dispersion_params;
    arma::vec dispersion_params_previous;
    arma::vec intercept_dispersion(dispersion_orders.n_param_intercept);
    arma::mat ar_params_dispersion(dispersion_orders.autoregressive_orders.n_rows, dispersion_orders.autoregressive_orders.n_cols);
    arma::mat ma_params_dispersion(dispersion_orders.moving_average_orders.n_rows, dispersion_orders.moving_average_orders.n_cols);
    arma::mat cov_params_dispersion(dispersion_orders.covariate_orders.n_rows, dispersion_orders.covariate_orders.n_cols);
    FittingObject * dispersion_fitting = nullptr;
    arma::mat dispersion_vals_mat;
    arma::mat init_dispersion = ts.head_cols(mean_orders.max_time_lag);
    init_dispersion = arma::repmat(arma::mean(init_dispersion, 1), 1, init_dispersion.n_cols);

    // Variables for tracking the fitting process
    unsigned int fit_counter = 0;
    bool not_finished = true;
    bool convergence_issues = false;
    bool mean_converged = false;
    bool dispersion_converged = false;
    double change_mean = convergence_threshold + 1.0;
    double change_dispersion = convergence_threshold + 1.0;
    double change_likelihood = convergence_threshold + 1.0;
    double change_likelihood_dispersion = convergence_threshold + 1.0;
    double change_likelihood_mean = convergence_threshold + 1.0;
    arma::vec log_likelihoods_mean(max_fits + 1, arma::fill::zeros);
    arma::vec log_likelihoods_dispersion(max_fits + 1, arma::fill::zeros);
    arma::vec total_log_likelihoods(max_fits + 1, arma::fill::zeros);
    arma::vec fitting_times_mean(max_fits + 1, arma::fill::zeros);
    arma::vec fitting_times_dispersion(max_fits + 1, arma::fill::zeros);
    arma::uvec mean_convergence_indices(max_fits + 1);
    arma::uvec dispersion_convergence_indices(max_fits + 1);
    arma::uvec fncounts_mean(max_fits + 1, arma::fill::zeros);
    arma::uvec fncounts_dispersion(max_fits + 1, arma::fill::zeros);
    arma::uvec grcounts_mean(max_fits + 1, arma::fill::zeros);
    arma::uvec grcounts_dispersion(max_fits + 1, arma::fill::zeros);
    arma::uvec hecounts_mean(max_fits + 1, arma::fill::zeros);
    arma::uvec hecounts_dispersion(max_fits + 1, arma::fill::zeros);
    arma::ivec convergence_codes_mean(max_fits + 1, arma::fill::zeros);
    arma::ivec convergence_codes_dispersion(max_fits + 1, arma::fill::zeros);
    arma::mat params_history_mean(mean_orders.n_param, max_fits + 1, arma::fill::zeros);
    arma::mat params_history_dispersion(dispersion_orders.n_param, max_fits + 1, arma::fill::zeros);
    arma::mat start_params_mean(mean_orders.n_param, max_fits + 1, arma::fill::zeros);
    arma::mat start_params_dispersion(dispersion_orders.n_param, max_fits + 1, arma::fill::zeros);


    // Perform first fit of mean model assuming constant dispersion parameter:
    Design mean_design(ts, &mean_covariates_list, link_init, &mean_orders, mean_fam, W_ar_mean, W_ma_mean, W_covariates_mean);
    Design dispersion_design = mean_design; // Placeholder
    FittingObject * mean_fitting = FittingObject::create(ts, &mean_design, &mean_orders, &mean_covariates_list, W_ma_mean, mean_fam, control_mean);
    arma::vec mean_params = mean_fitting->fit(start_param_mean);
    arma::vec mean_params_previous = mean_params; // For convergence check
    convergence_issues = convergence_issues || !mean_fitting->is_converged();
    mean_convergence_indices(0) = mean_fitting->is_converged() ? 1 : 0;
    params_history_mean.col(0) = mean_params;
    start_params_mean.col(0) = start_param_mean;
    fitting_times_mean(0) = mean_fitting->get_fitting_time() / 1000.0; // Convert to seconds
    fncounts_mean(0) = mean_fitting->get_fncount();
    grcounts_mean(0) = mean_fitting->get_grcount();
    hecounts_mean(0) = mean_fitting->get_hecount();
    convergence_codes_mean(0) = mean_fitting->get_convergence_code();
    if(mean_orders.n_param_ma > 0)
    {
        W_ma_mean->set_parameter_matrices(ma_params_mean);
    }
    mean_design.update_design(mean_params, &mean_orders, mean_fam, W_ma_mean);
    arma::mat link_values = mean_design.link_vals;
    Parameter::create_param_matrices(mean_params, &mean_orders, intercept_mean, ar_params_mean, ma_params_mean, cov_params_mean);
    if(print_progress)
    {
        Rcpp::Rcout << "Mean fit 0 - Time required: " << mean_fitting->get_fitting_time() / 1000.0 << "s\n";
    }
    log_likelihoods_mean(0) = arma::accu( mean_fam->log_likelihood( ts.tail_cols(mean_orders.n_obs_effective), mean_fam->inverse_link(link_values.tail_cols(mean_orders.n_obs_effective)), arma::ones(ts.n_rows, mean_orders.n_obs_effective) ) );
  
    // Use these observations to calculate the total likelihood (without bias)
    const unsigned int n_obs_effective = ts.n_cols - (mean_orders.max_time_lag + dispersion_orders.max_time_lag);


    // arma::mat dispersion_vals_mat_old(ts.n_rows, mean_orders.n_obs_effective, arma::fill::ones);
    arma::vec proposed_params_mean;
    arma::vec proposed_params_dispersion;
    dispersion_params_previous = arma::vec(dispersion_orders.n_param, arma::fill::zeros);
    dispersion_params = arma::vec(dispersion_orders.n_param, arma::fill::zeros);
    // Start iterative procedure
    while(not_finished)
    {
        Rcpp::checkUserInterrupt();
        if(!dispersion_converged)
        {
            // Calculate pseudo observations for dispersion model
            pseudo = link_values;
            pseudo = mean_fam->inverse_link(pseudo);
            if(!drop_max_mean_lag)
            {
                pseudo.head_cols(mean_orders.max_time_lag) = init_dispersion;
            }
            pseudo = (pseudo_obs == "deviance") ? mean_fam->deviance_residual(ts, pseudo) : mean_fam->pearson_residual(ts, pseudo);
            if(drop_max_mean_lag) {
                pseudo = pseudo.tail_cols(mean_orders.n_obs_effective);
            }
            pseudo.clamp(lower_dispersion, upper_dispersion); // Clamping to avoid numerical issues

            // Estimate dispersion Model
            if(const_dispersion_model && use_fast_if_const_dispersion)
            {
                dispersion_vals_mat = (dispersion_orders.intercepts_equal) ? arma::mat(ts.n_rows, mean_orders.n_obs_effective, arma::fill::value(arma::accu(pseudo))) : arma::repmat(arma::sum(pseudo, 1), 1, mean_orders.n_obs_effective);
                dispersion_vals_mat /= (dispersion_orders.intercepts_equal) ? (ts.n_rows * mean_orders.n_obs_effective - mean_orders.n_param) : mean_orders.n_obs_effective - (mean_orders.n_param / ts.n_rows); // TODO: Check if this is correct
                dispersion_params = (dispersion_orders.intercepts_equal) ? arma::vec(1, arma::fill::value(dispersion_vals_mat(0, 0))) : dispersion_vals_mat.col(0);
                log_likelihoods_dispersion(fit_counter) = arma::accu(arma::log(dispersion_vals_mat) + pseudo.tail_cols(mean_orders.n_obs_effective) / dispersion_vals_mat);
            } else {
                dispersion_init = Model::init_link(init_method_dispersion, dispersion_orders, pseudo, dispersion_fam);
                if(previous_param_as_start && fit_counter > 0)
                {
                    start_param_dispersion = dispersion_params;
                } else {
                    start_param_dispersion = Model::init_param(init_method_param_dispersion, dispersion_orders, pseudo, dispersion_fam);
                }
                start_params_dispersion.col(fit_counter) = start_param_dispersion;
                Parameter::create_param_matrices(start_param_dispersion, &dispersion_orders, intercept_dispersion, ar_params_dispersion, ma_params_dispersion, cov_params_dispersion);
                if(dispersion_orders.n_param_ma > 0)
                {
                    W_ma_dispersion->set_parameter_matrices(ma_params_dispersion);
                }
                
                dispersion_design = Design(pseudo, &dispersion_covariates_list, dispersion_init, &dispersion_orders, dispersion_fam, W_ar_dispersion, W_ma_dispersion, W_covariates_dispersion);
                delete dispersion_fitting; // Free memory if it was already created
                dispersion_fitting = FittingObject::create(pseudo, &dispersion_design, &dispersion_orders, &dispersion_covariates_list, W_ma_dispersion, dispersion_fam, control_dispersion);
                
                proposed_params_dispersion = dispersion_fitting->fit(start_param_dispersion);
                convergence_issues = convergence_issues || !dispersion_fitting->is_converged();
                dispersion_convergence_indices(fit_counter) = dispersion_fitting->is_converged() ? 1 : 0;
                fncounts_dispersion(fit_counter) = dispersion_fitting->get_fncount();
                grcounts_dispersion(fit_counter) = dispersion_fitting->get_grcount();
                hecounts_dispersion(fit_counter) = dispersion_fitting->get_hecount();
                fitting_times_dispersion(fit_counter) = dispersion_fitting->get_fitting_time() / 1000.0; // Convert to seconds
                convergence_codes_dispersion(fit_counter) = dispersion_fitting->get_convergence_code();

                if(use_backtracking)
                {
                    double old_ll = (fit_counter > 1) ? total_log_likelihoods(fit_counter - 1) : -arma::datum::inf;
                    arma::vec old_params = dispersion_params_previous;

                    // Backtracking line search to ensure increase in likelihood
                    double current_alpha = alpha_start;
                    bool success = false;
                    arma::vec best_params_this_round = proposed_params_dispersion;

                    while(current_alpha >= min_alpha)
                    {
                        // weighted mixing: alpha * new + (1-alpha) * old
                        arma::vec candidate_params = current_alpha * proposed_params_dispersion + (1.0 - current_alpha) * dispersion_params;
                        
                        // Update Design & calculate Likelihood
                        Parameter::create_param_matrices(candidate_params, &dispersion_orders, intercept_dispersion, ar_params_dispersion, ma_params_dispersion, cov_params_dispersion);
                        if(dispersion_orders.n_param_ma > 0)
                        {
                            W_ma_dispersion->set_parameter_matrices(ma_params_dispersion);  
                        }
                        dispersion_design.update_design(candidate_params, &dispersion_orders, dispersion_fam, W_ma_dispersion);
                        arma::mat test_link = dispersion_design.link_vals;
                        test_link = dispersion_fam->inverse_link(test_link);
                        test_link = test_link.tail_cols(n_obs_effective);
                        if(mean_fam->family_in_R == "negative_binomial")
                        {
                            test_link -= 1.0;
                            test_link /= mean_fam->inverse_link(link_values.tail_cols(n_obs_effective));
                            test_link.clamp(0.0, arma::datum::inf);
                        }
                        double test_ll = arma::accu(mean_fam->log_likelihood(ts.tail_cols(n_obs_effective), mean_fam->inverse_link(link_values.tail_cols(n_obs_effective)), test_link));
                        if(test_ll > old_ll || fit_counter == 0) {
                            best_params_this_round = candidate_params;
                            success = true;
                            break; // Improvement found!
                        }
                        current_alpha *= alpha_shrink; // Reduce step size
                    }
                    // If no better point found, stay at old (prevents divergence)
                    if(!success) {
                        dispersion_params = old_params;
                        total_log_likelihoods(fit_counter) = old_ll;
                        if(print_progress && print_warnings) {
                            Rcpp::warning("No improvement found in dispersion parameter update, reverting to previous parameters.");
                        } 
                    } else {
                        dispersion_params = best_params_this_round;
                    }
                } else {
                    dispersion_params = proposed_params_dispersion;
                }
                // Update model with accepted parameters
                params_history_dispersion.col(fit_counter) = dispersion_params;
                Parameter::create_param_matrices(dispersion_params, &dispersion_orders, intercept_dispersion, ar_params_dispersion, ma_params_dispersion, cov_params_dispersion);
                if(dispersion_orders.n_param_ma > 0)
                {
                    W_ma_dispersion->set_parameter_matrices(ma_params_dispersion);  
                }
                dispersion_design.update_design(dispersion_params, &dispersion_orders, dispersion_fam, W_ma_dispersion);
                dispersion_vals_mat = dispersion_design.link_vals;
                dispersion_vals_mat = dispersion_fam->inverse_link(dispersion_vals_mat);
                log_likelihoods_dispersion(fit_counter) = arma::accu(arma::log(dispersion_vals_mat.tail_cols(dispersion_orders.n_obs_effective)) + pseudo.tail_cols(dispersion_orders.n_obs_effective) / dispersion_vals_mat.tail_cols(dispersion_orders.n_obs_effective));
            }
            
            dispersion_vals_mat = dispersion_vals_mat.tail_cols(mean_orders.n_obs_effective);
            if(mean_fam->family_in_R == "negative_binomial")
            {
                dispersion_vals_mat -= 1.0;
                dispersion_vals_mat /= mean_fam->inverse_link(link_values.tail_cols(mean_orders.n_obs_effective));
                dispersion_vals_mat.clamp(0.0, arma::datum::inf);
            }
            total_log_likelihoods(fit_counter) = arma::accu( mean_fam->log_likelihood(ts.tail_cols(n_obs_effective), mean_fam->inverse_link(link_values.tail_cols(n_obs_effective)), dispersion_vals_mat.tail_cols(n_obs_effective)));
            
            if(fit_counter > 0)
            {
                change_dispersion = arma::norm(dispersion_params - dispersion_params_previous, 2) / (arma::norm(dispersion_params_previous, 2) + 1.0);
                change_likelihood = std::abs(total_log_likelihoods(fit_counter) - total_log_likelihoods(fit_counter - 1)) / (std::abs(total_log_likelihoods(fit_counter - 1)) + 1.0);
                change_likelihood_dispersion = std::abs(log_likelihoods_dispersion(fit_counter) - log_likelihoods_dispersion(fit_counter - 1)) / (std::abs(log_likelihoods_dispersion(fit_counter - 1)) + 1.0);
            }

            if(print_progress)
            {
                if(const_dispersion_model && use_fast_if_const_dispersion)
                {
                    Rcpp::Rcout << "Dispersion fit " << fit_counter << " - Time required: NA s";
                } else {
                    Rcpp::Rcout << "Dispersion fit " << fit_counter << " - Time required: " << dispersion_fitting->get_fitting_time() / 1000.0 << "s";
                }
                (fit_counter > 0) ? Rcpp::Rcout << ", rel. change in parameters: " << change_dispersion << ", rel. change in ll: " << change_likelihood_dispersion << "\n" : Rcpp::Rcout << "\n";
            }
            dispersion_params_previous = dispersion_params;
            dispersion_converged = (change_dispersion < convergence_threshold) || (change_likelihood_dispersion < convergence_threshold);
        } else {
            params_history_dispersion.col(fit_counter) = dispersion_params;
        }



        // Check stopping criteria
        not_finished = !((mean_converged && dispersion_converged) || fit_counter >= max_fits || change_likelihood < convergence_threshold);
        if(!not_finished)
        {
            break;
        }

        Rcpp::checkUserInterrupt();
        // Refit Mean-Model if not yet converged
        if(!mean_converged)
        {
            fit_counter++;
            mean_fam->const_dispersion = false;
            mean_fam->dispersion_matrix = dispersion_vals_mat;
            mean_params = previous_param_as_start ? mean_params_previous : start_param_mean;
            start_params_mean.col(fit_counter) = mean_params;
            Parameter::create_param_matrices(mean_params, &mean_orders, intercept_mean, ar_params_mean, ma_params_mean, cov_params_mean);
            if(mean_orders.n_param_ma > 0)
            {
                W_ma_mean->set_parameter_matrices(ma_params_mean);
            }
            proposed_params_mean = mean_fitting->fit(mean_params);
            // mean_params = mean_fitting->fit(mean_params); 

            convergence_issues = convergence_issues || !mean_fitting->is_converged();
            mean_convergence_indices(fit_counter) = mean_fitting->is_converged() ? 1 : 0;
            fncounts_mean(fit_counter) = mean_fitting->get_fncount();
            grcounts_mean(fit_counter) = mean_fitting->get_grcount();
            hecounts_mean(fit_counter) = mean_fitting->get_hecount();
            fitting_times_mean(fit_counter) = mean_fitting->get_fitting_time() / 1000.0; // Convert to seconds
            convergence_codes_mean(fit_counter) = mean_fitting->get_convergence_code();

            if(use_backtracking)
            {
                double old_ll = (fit_counter > 0) ? total_log_likelihoods(fit_counter - 1) : -arma::datum::inf;
                arma::vec old_params = mean_params_previous;

                // Backtracking line search to ensure increase in likelihood
                double current_alpha = alpha_start;
                bool success = false;
                arma::vec best_params_this_round = proposed_params_mean;

                while(current_alpha >= min_alpha)
                {
                    // weighted mixing: alpha * new + (1-alpha) * old
                    arma::vec candidate_params = current_alpha * proposed_params_mean + (1.0 - current_alpha) * mean_params;
        
                    // Update Design & calculate test likelihood
                    Parameter::create_param_matrices(candidate_params, &mean_orders, intercept_mean, ar_params_mean, ma_params_mean, cov_params_mean);
                    if(mean_orders.n_param_ma > 0)
                    {
                        W_ma_mean->set_parameter_matrices(ma_params_mean);
                    }
                    mean_design.update_design(candidate_params, &mean_orders, mean_fam, W_ma_mean);
                    arma::mat test_link = mean_design.link_vals;
                    double test_ll = arma::accu(mean_fam->log_likelihood(ts.tail_cols(n_obs_effective), mean_fam->inverse_link(test_link.tail_cols(n_obs_effective)), dispersion_vals_mat.tail_cols(n_obs_effective)));
                    // Check if improvement occurred (or first iteration)
                    if(test_ll > old_ll || fit_counter == 1) {
                        best_params_this_round = candidate_params;
                        total_log_likelihoods(fit_counter) = test_ll;
                        //link_values = test_link;
                        success = true;
                        break; // Verbesserung gefunden!
                    }
                    current_alpha *= alpha_shrink; // Reduce step size
                }
                // If no better point was found, stay with the old one (prevents divergence)
                if(!success) {
                    mean_params = old_params;
                    total_log_likelihoods(fit_counter) = old_ll;
                    if(print_progress && print_warnings) {
                            Rcpp::warning("No improvement found in dispersion parameter update, reverting to previous parameters.");
                    } 
                } else {
                    mean_params = best_params_this_round;
                }
            } else {
                // Ohne Backtracking
                mean_params = proposed_params_mean;
            }

            params_history_mean.col(fit_counter) = mean_params;
            Parameter::create_param_matrices(mean_params, &mean_orders, intercept_mean, ar_params_mean, ma_params_mean, cov_params_mean);
            if(mean_orders.n_param_ma > 0)
            {
                W_ma_mean->set_parameter_matrices(ma_params_mean);
            }
            mean_design.update_design(mean_params, &mean_orders, mean_fam, W_ma_mean);
            link_values = mean_design.link_vals;
            change_mean = arma::norm(mean_params - mean_params_previous, 2) / (arma::norm(mean_params_previous, 2) + 1.0);
            log_likelihoods_mean(fit_counter) = arma::accu( mean_fam->log_likelihood( ts.tail_cols(mean_orders.n_obs_effective), mean_fam->inverse_link(link_values.tail_cols(mean_orders.n_obs_effective)), dispersion_vals_mat.tail_cols(mean_orders.n_obs_effective) ) );
            change_likelihood_mean = std::abs(log_likelihoods_mean(fit_counter) - log_likelihoods_mean(fit_counter - 1)) / (std::abs(log_likelihoods_mean(fit_counter - 1)) + 1.0);
            if(print_progress)
            {
                Rcpp::Rcout << "Mean fit " << fit_counter << " - Time required: " << mean_fitting->get_fitting_time() / 1000.0 << "s";
                (fit_counter > 0) ? Rcpp::Rcout << ", rel. change in parameters: " << change_mean << ", rel. change in ll: " << change_likelihood_mean << "\n" : Rcpp::Rcout << "\n";
            }
            mean_params_previous = mean_params;
            mean_converged = (change_mean < convergence_threshold) || (change_likelihood_mean < convergence_threshold);
        } 
    }

    if(convergence_issues && print_warnings)
    {
        Rcpp::warning("There were convergence issues in at least one of the model fits (mean or dispersion). Please check the results carefully and consider increasing 'maxit' or different model orders.");
    }
    if(change_likelihood > convergence_threshold && !(mean_converged && dispersion_converged) && print_warnings)
    {
        Rcpp::warning("Model did not converge after " + std::to_string(fit_counter) + " iterations. Use results with care and consider increasing the maximum number of iterations ('max_fits') or a different model specification.");
    }

    // Merge results and clean up

    // Information about algorithm and convergence
    Rcpp::List algorithm_info;
    algorithm_info["mean"] = mean_fitting->get_algorithm_info();
    if(const_dispersion_model)
    {
        algorithm_info["dispersion"] = Rcpp::List::create(
            Rcpp::Named("method") = "moment estimator", 
            Rcpp::Named("algorithm") = "none"
        );
    } else {
        algorithm_info["dispersion"] = dispersion_fitting->get_algorithm_info();
    }
    Rcpp::List convergence_info;
    convergence_info["start"] = Rcpp::List::create(
        Rcpp::Named("mean") = start_params_mean.head_cols(fit_counter),
        Rcpp::Named("dispersion") = start_params_dispersion.head_cols(fit_counter)
    );
    convergence_info["fncount"] = Rcpp::List::create(
        Rcpp::Named("mean") = fncounts_mean.head(fit_counter),
        Rcpp::Named("dispersion") = fncounts_dispersion.head(fit_counter)
    );
    convergence_info["grcount"] = Rcpp::List::create(
        Rcpp::Named("mean") = grcounts_mean.head(fit_counter),
        Rcpp::Named("dispersion") = grcounts_dispersion.head(fit_counter)
    );
    convergence_info["hecount"] = Rcpp::List::create(
        Rcpp::Named("mean") = hecounts_mean.head(fit_counter),
        Rcpp::Named("dispersion") = hecounts_dispersion.head(fit_counter)
    );
    convergence_info["fitting_time"] = Rcpp::List::create(
        Rcpp::Named("mean") = fitting_times_mean.head(fit_counter),
        Rcpp::Named("dispersion") = fitting_times_dispersion.head(fit_counter)
    );
    convergence_info["convergence"] = Rcpp::List::create(
        Rcpp::Named("mean") = mean_convergence_indices.head(fit_counter),
        Rcpp::Named("dispersion") = dispersion_convergence_indices.head(fit_counter)
    );
    convergence_info["convergence_codes"] = Rcpp::List::create(
        Rcpp::Named("mean") = convergence_codes_mean.head(fit_counter),
        Rcpp::Named("dispersion") = convergence_codes_dispersion.head(fit_counter)
    );
    convergence_info["outer_iterations"] = fit_counter;
    
    // Get estimations, score and variance:
    Parameter mean_parameters(mean_params, mean_orders);
    Rcpp::List param_est_mean = mean_parameters.return_param_list();
    Rcpp::List param_est_dispersion;
    if(const_dispersion_model)
    {
        param_est_dispersion = Rcpp::List::create(
            Rcpp::Named("intercept") = dispersion_params
        );
    } else {
        Parameter dispersion_parameters(dispersion_params, dispersion_orders);
        param_est_dispersion = dispersion_parameters.return_param_list();
    }

    arma::mat fitted_values_mean = mean_fam->inverse_link(arma::reshape(link_values, ts.n_rows, link_values.n_elem / ts.n_rows));
    arma::mat fitted_values_dispersion = dispersion_vals_mat;

    arma::vec score_mean(mean_params.n_elem);
    arma::mat information_mean(mean_params.n_elem, mean_params.n_elem);
    arma::vec ts_vec_mean = arma::vectorise(ts.tail_cols(mean_orders.n_obs_effective));
    double log_likelihood_mean = Model::log_likelihood(mean_params, score_mean, information_mean, ts_vec_mean, &mean_design, &mean_orders, W_ma_mean, mean_fam, true, true);
    log_likelihood_mean /= mean_orders.n_obs_effective;
    log_likelihood_mean *= ts.n_cols;
    arma::mat variance_estimation_mean = Model::variance_estimation(mean_params, ts_vec_mean, &mean_design, mean_fam, information_mean, W_ma_mean, &mean_orders);
    arma::mat var_temp_mean = information_mean * variance_estimation_mean / mean_orders.n_obs_effective;
    double aic_mean = -2.0 * log_likelihood_mean + 2.0 * mean_orders.n_param;
    double bic_mean = -2.0 * log_likelihood_mean + std::log(static_cast<double>(mean_orders.n_obs_effective * ts.n_rows)) * mean_orders.n_param;
    double qic_mean = -2.0 * log_likelihood_mean + 2.0 * arma::trace(var_temp_mean);
    double qic_total = 2.0 * arma::trace(var_temp_mean);

    Rcpp::List results_mean;
    results_mean["ts"] = ts;
    results_mean["covariates"] = mean_covariates;
    results_mean["model"] = mean_model;
    results_mean["family"] = mean_family;
    results_mean["coefficients_list"] = param_est_mean;
    results_mean["coefficients"] = mean_params;
    results_mean["n_obs_effective"] = mean_orders.n_obs_effective;
    results_mean["max_time_lag"] = mean_orders.max_time_lag;
    results_mean["log_likelihood"] = log_likelihood_mean;
    results_mean["score"] = score_mean;
    results_mean["information"] = information_mean;
    results_mean["variance_estimation"] = variance_estimation_mean;
    results_mean["aic"] = aic_mean;
    results_mean["bic"] = bic_mean;
    results_mean["qic"] = qic_mean;
    results_mean["design_matrix"] = mean_design.design_matrix;
    results_mean["derivatives"] = mean_design.derivative_link;
    results_mean["fitted.values"] = fitted_values_mean;
    results_mean["link_values"] = mean_design.link_vals;
    results_mean["param_history"] = params_history_mean.head_cols(fit_counter);
    results_mean["log_likelihood_history"] = log_likelihoods_mean.head(fit_counter);


    Rcpp::List results_dispersion;
    results_dispersion["ts"] = pseudo;
    results_dispersion["pseudo_type"] = pseudo_obs;
    results_dispersion["covariates"] = dispersion_covariates;
    results_dispersion["model"] = dispersion_model;
    results_dispersion["family"] = dispersion_family;
    results_dispersion["coefficients_list"] = param_est_dispersion;
    results_dispersion["coefficients"] = dispersion_params;
    results_dispersion["n_obs_effective"] = dispersion_orders.n_obs_effective;
    results_dispersion["max_time_lag"] = dispersion_orders.max_time_lag;

    if(!(const_dispersion_model && use_fast_if_const_dispersion))
    {
        arma::vec score_dispersion(dispersion_params.n_elem);
        arma::mat information_dispersion(dispersion_params.n_elem, dispersion_params.n_elem);
        arma::vec ts_vec_dispersion = arma::vectorise(pseudo.tail_cols(dispersion_orders.n_obs_effective));
        double log_likelihood_dispersion = Model::log_likelihood(dispersion_params, score_dispersion, information_dispersion, ts_vec_dispersion, &dispersion_design, &dispersion_orders, W_ma_dispersion, dispersion_fam, true, true);
        log_likelihood_dispersion /= dispersion_orders.n_obs_effective;
        log_likelihood_dispersion *= ts.n_cols;
        arma::mat variance_estimation_dispersion = Model::variance_estimation(dispersion_params, ts_vec_dispersion, &dispersion_design, dispersion_fam, information_dispersion, W_ma_dispersion, &dispersion_orders);
        arma::mat var_temp_dispersion = information_dispersion * variance_estimation_dispersion / dispersion_orders.n_obs_effective;
        double aic_dispersion = -2.0 * log_likelihood_dispersion + 2.0 * dispersion_orders.n_param;
        double bic_dispersion = -2.0 * log_likelihood_dispersion + std::log(static_cast<double>(dispersion_orders.n_obs_effective * ts.n_rows)) * dispersion_orders.n_param;
        double qic_dispersion = -2.0 * log_likelihood_dispersion + 2.0 * arma::trace(var_temp_dispersion);
        qic_total += 2.0 * arma::trace(var_temp_dispersion);
        results_dispersion["log_likelihood"] = log_likelihood_dispersion;
        results_dispersion["score"] = score_dispersion;
        results_dispersion["information"] = information_dispersion;
        results_dispersion["variance_estimation"] = variance_estimation_dispersion;
        results_dispersion["aic"] = aic_dispersion;
        results_dispersion["bic"] = bic_dispersion;
        results_dispersion["qic"] = qic_dispersion;
        results_dispersion["design_matrix"] = dispersion_design.design_matrix;
        results_dispersion["derivatives"] = dispersion_design.derivative_link;
        results_dispersion["link_values"] = dispersion_design.link_vals;
    }
    results_dispersion["fitted.values"] = fitted_values_dispersion;
    results_dispersion["dispersion_matrix"] = dispersion_vals_mat;
    results_dispersion["param_history"] = params_history_dispersion.head_cols(fit_counter);
    results_dispersion["log_likelihood_history"] = log_likelihoods_dispersion.head(fit_counter);

    // Free memory
    delete mean_fam;
    delete dispersion_fam;
    delete W_ar_mean;
    delete W_ma_mean;
    delete W_covariates_mean;
    delete W_ar_dispersion;
    delete W_ma_dispersion;
    delete W_covariates_dispersion;
    delete mean_fitting;
    delete dispersion_fitting;
    mean_fitting = nullptr; 
    dispersion_fitting = nullptr;

    // Scale average log-likelihood to number of observations
    double total_ll = total_log_likelihoods(fit_counter - 1);
    total_ll /= n_obs_effective;
    total_ll *= ts.n_cols;
    double aic_total = -2.0 * total_ll + 2.0 * (dispersion_orders.n_param + mean_orders.n_param);
    double bic_total = -2.0 * total_ll + std::log(static_cast<double>(ts.n_cols * ts.n_rows)) * (dispersion_orders.n_param + mean_orders.n_param);
    qic_total += -2.0 * total_ll;
    if(const_dispersion_model && use_fast_if_const_dispersion)
    {
        qic_total += 2.0 * dispersion_orders.n_param;
    } 

    Rcpp::List results;
    results["mean"] = results_mean;
    results["dispersion"] = results_dispersion;
    results["target_dim"] = ts.n_rows;
    results["algorithm_info"] = algorithm_info;
    results["convergence_info"] = convergence_info;
    results["total_log_likelihood_history"] = total_log_likelihoods.head(fit_counter);
    results["total_log_likelihood"] = total_ll;
    results["aic"] = aic_total;
    results["bic"] = bic_total;
    results["qic"] = qic_total;
    return results;
}
