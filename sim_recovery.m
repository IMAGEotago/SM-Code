function all_results = sim_recovery()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BLT Behavioural Data: Simulation and Parameter Recovery                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------- %
% Author: Sophie Morris                                                   %
% Created: 18/03/2026                                                     %
% ----------------------------------------------------------------------- %
% Description: This script runs simulation and parameter recovery for     %
% behavioural data from the Breathing Learning Task (BLT). A timestamped  %
% analysis plan is available on (https://github.com/IMAGEotago).          %
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set paths
addpath(genpath("C:/Users/sarmo/OneDrive/Desktop/VBA-toolbox-master"));
addpath(genpath("C:/Users/sarmo/OneDrive/Desktop/computational_modelling"));
load('pilot_analysis_uninformed.mat', 'pilot_results');

%% Configurations
% n_seeds may be adjusted for less 'runs' if needed
n_seeds = 10;
seeds = [42, 234, 1264, 2001, 5247, 4521, 7816, 2458, 441, 647];

%% Initialize parallel processing
use_parallel = true;  
n_workers = [];

if use_parallel
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        try
            if isempty(n_workers)
                poolobj = parpool('Processes');
            else
                poolobj = parpool('Processes', n_workers);
            end
            actual_workers = poolobj.NumWorkers;
        catch ME
            use_parallel = false;
            actual_workers = 1;
        end
    else
        actual_workers = poolobj.NumWorkers;
    end
end 

%% Define priors to test
priors_to_test = {};

% Starting priors are shown below. To test multiple configurations, specify
% e.g. [-0.75 -0.5 -0.25] for any of the 'lists' below.
mutheta_list = 0;
sigmatheta_list = 1;
muphi_list = 0;
sigmaphi_list = 1;

alpha_priors = repmat(struct('muTheta',[],'SigmaTheta',[]), 0, 1);

for mutheta = mutheta_list
    for sgtheta = sigmatheta_list
        alpha_priors(end+1) = struct('muTheta', mutheta, 'SigmaTheta', sgtheta);
    end
end

zeta_priors = repmat(struct('muPhi',[],'SigmaPhi',[]), 0, 1);

for muphi = muphi_list
    for sgphi = sigmaphi_list
        zeta_priors(end+1) = struct('muPhi', muphi, 'SigmaPhi', sgphi);
    end
end

idx = 1;

for a = 1:length(alpha_priors)
    for z = 1:length(zeta_priors)
        priors_to_test{idx} = struct(...
            'muTheta', alpha_priors(a).muTheta, ...
            'SigmaTheta', alpha_priors(a).SigmaTheta, ...
            'muPhi', zeta_priors(z).muPhi, ...
            'SigmaPhi', zeta_priors(z).SigmaPhi, ...
            'alpha_name', sprintf('α(%.2f,%.2f)', alpha_priors(a).muTheta, alpha_priors(a).SigmaTheta), ...
            'zeta_name', sprintf('ζ(%.2f,%.2f)', zeta_priors(z).muPhi, zeta_priors(z).SigmaPhi));
        idx = idx + 1;
    end
end 

n_priors = length(priors_to_test);

%% Pilot parameter ranges and noise levels
alpha_min = max(0.01, min(pilot_results.alphas) * 0.8); % 0.07
alpha_max = min(0.99, max(pilot_results.alphas) * 1.2);  % 0.77

zeta_min = max(0.1, min(pilot_results.zetas) * 0.8);  % floor at 0.1 
zeta_max = max(pilot_results.zetas) * 1.2;  % 10.14

% Below values may be adjusted to test varying levels of observation noise
sigma_values_test = [0.05; 0.15; 0.25; 0.35];

%% Benchmarking
pairings = load("C:/Users/sarmo/OneDrive/Desktop/computational_modelling/blt_pairings.mat");
pairings_for_benchmark = double(pairings.data(:))';
N_for_benchmark = length(pairings_for_benchmark);

% Save parameters for main loop
pairings_main = pairings_for_benchmark;
N_main = N_for_benchmark;
alpha_range = [alpha_min, alpha_max];
zeta_range = [zeta_min, zeta_max];

fprintf('\n== ESTIMATING RUNTIME ==\n');
fprintf('Running benchmark with 10 simulations...\n');

benchmark_times = zeros(10, 1);
benchmark_prior = priors_to_test{1};

for bench_idx = 1:10
    tic;
    run_single_two_param_recovery_with_prior(0.3, 2.0, 0.1, ...
        pairings_for_benchmark, N_for_benchmark, benchmark_prior);
    benchmark_times(bench_idx) = toc;
end

mean_fit_time = mean(benchmark_times);
total_fits = n_priors * 4 * 500 * n_seeds;
estimated_total_seconds = total_fits * mean_fit_time;
estimated_hours = estimated_total_seconds / 3600;
estimated_mins = estimated_total_seconds/ 60;
estimated_hours_parallel = estimated_hours / actual_workers;
estimated_mins_parallel = estimated_mins / actual_workers;

fprintf('Total fits needed: %s\n', num2str(total_fits, '%d'));
fprintf('\nEstimated total runtime:\n');
fprintf('  Parallel (%d workers): %.1f minutes\n', actual_workers, estimated_mins_parallel);
fprintf('  Parallel (%d workers): %.1f hours\n', actual_workers, estimated_hours_parallel);

%% Run across multiple seeds
all_seeds_results = cell(n_seeds, 1);
total_prior_fits = n_priors * n_seeds;
completed_prior_fits = 0;
start_time_total = tic;

fprintf('\n== STARTING MAIN ANALYSIS ==\n');

for seed_idx = 1:n_seeds
    current_seed_results = cell(n_priors, 1);
    
     if use_parallel
        parfor p = 1:n_priors
            current_seed_results{p} = run_two_param_recovery_with_prior(priors_to_test{p}, ...
                seeds(seed_idx), pairings_main, N_main, alpha_range, zeta_range, ...
                sigma_values_test);
        end
    else
        for p = 1:n_priors
            if mod(p, 20) == 0 || p == 1
                elapsed_total = toc(start_time_total);
                if completed_prior_fits > 0
                    avg_time_per_prior = elapsed_total / completed_prior_fits;
                    remaining_priors = total_prior_fits - completed_prior_fits;
                    estimated_remaining_hours = (remaining_priors * avg_time_per_prior) / 3600;
                    fprintf('  Prior %d/%d (%.1f%% complete, ETA: %.1f hours)...\n', ...
                        p, n_priors, (completed_prior_fits/total_prior_fits)*80, ...
                        estimated_remaining_hours);
                else
                    fprintf('  Prior %d/%d...\n', p, n_priors);
                end
            end
            
            current_seed_results{p} = run_two_param_recovery_with_prior(priors_to_test{p}, ...
                seeds(seed_idx), pairings_main, N_main, alpha_range, zeta_range, ...
                sigma_values_test);
            completed_prior_fits = completed_prior_fits + 1;
        end
    end

    all_seeds_results{seed_idx} = current_seed_results;

    if use_parallel
        completed_prior_fits = seed_idx * n_priors;
    end
    
    save('pilot_informed_recovery_INPROGRESS.mat', 'all_seeds_results', 'priors_to_test', ...
        'seeds', 'seed_idx', 'alpha_range', 'zeta_range', 'sigma_values_test', 'pilot_results');
end

%% Aggregate results
fprintf('\n== AGGREGATING RESULTS ==\n');
all_results = aggregate_across_seeds(all_seeds_results, priors_to_test);

%% Analysis sections
check_convergence_comprehensive(all_results);
best_priors = identify_best_priors_multiseed(all_results);
analyze_seed_stability(all_results, best_priors);
analyze_error_correlation_by_sigmaPhi(all_results);
create_correlation_plots(all_results, best_priors);

%% STOPPING RULE CHECK
fprintf('\n== STOPPING RULE CHECK ==\n');

summary_file = 'refinement_summary.mat';

best_idx = best_priors.best_prior_number;
best_result = all_results{best_idx};
sigma_vals = best_result.sigma_values;

fprintf('\nBest prior recovery correlations:\n');

for j = 1:length(sigma_vals)
    fprintf('σ = %.2f   alpha r = %.3f   zeta r = %.3f\n', sigma_vals(j), ...
        best_result.alpha_corr_mean(j), best_result.zeta_corr_mean(j));
end

fprintf('\nOverall mean recovery:\n');
fprintf('alpha r = %.3f\n', mean(best_result.alpha_corr_mean));
fprintf('zeta r  = %.3f\n', mean(best_result.zeta_corr_mean));

current_summary = struct();

current_summary.best_prior_index = best_idx;
current_summary.muTheta = best_result.alpha_prior_used.muTheta;
current_summary.SigmaTheta = best_result.alpha_prior_used.SigmaTheta;
current_summary.muPhi = best_result.zeta_prior_used.muPhi;
current_summary.SigmaPhi = best_result.zeta_prior_used.SigmaPhi;

current_summary.alpha_corr_mean = mean(best_result.alpha_corr_mean);
current_summary.zeta_corr_mean  = mean(best_result.zeta_corr_mean);
current_summary.error_corr_mean = mean(best_result.error_corr_mean);
current_summary.mean_abs_error_corr = mean(abs(best_result.error_corr_mean));

current_summary.alpha_corr_by_sigma = best_result.alpha_corr_mean;
current_summary.zeta_corr_by_sigma  = best_result.zeta_corr_mean;
current_summary.error_corr_by_sigma = best_result.error_corr_mean;

stop_refinement = false;

if exist(summary_file,'file')

    S = load(summary_file);
    previous_summary = S.summary;

    delta_alpha_corr = abs(current_summary.alpha_corr_mean - previous_summary.alpha_corr_mean);
    delta_zeta_corr  = abs(current_summary.zeta_corr_mean  - previous_summary.zeta_corr_mean);

    error_corr_ok = current_summary.mean_abs_error_corr <= 0.35;

    fprintf('Δ alpha recovery = %.4f\n', delta_alpha_corr);
    fprintf('Δ zeta recovery  = %.4f\n', delta_zeta_corr);
    fprintf('Mean |error corr| across noise levels = %.4f\n', ...
        current_summary.mean_abs_error_corr);

    if delta_alpha_corr < 0.02 && delta_zeta_corr < 0.02 && error_corr_ok
        stop_refinement = true;
        fprintf('Stopping rule met.\n');
    else
        fprintf('Stopping rule not met.\n');
    end

else
    fprintf('No previous comparison.\n');
end

% Save summary for next refinement round
summary = current_summary;
save(summary_file,'summary','stop_refinement');

fprintf('Summary saved for next refinement round.\n');

%% Save results
save('pilot_informed_recovery_results.mat', 'all_results', ...
    'all_seeds_results', 'priors_to_test', 'best_priors', 'seeds', ...
    'pilot_results', 'alpha_range', 'zeta_range', 'sigma_values_test');

if use_parallel
    delete(gcp('nocreate'));
end

end

%% Recovery function
function results = run_two_param_recovery_with_prior(prior_settings, seed, pairings, N, ...
    alpha_range, zeta_range, sigma_values)
    
    % Adjust n_subjects as needed - 100 were used here initially to
    % iteratively test prior configurations, and 500 were then used to
    % validate the final selected prior configuration
    n_subjects = 100; 
    n_sigmas = length(sigma_values);
    
    rng(seed);
     
    % Generate true parameters
    true_alphas = alpha_range(1) + (alpha_range(2) - alpha_range(1)) * rand(n_subjects, 1);
    true_zetas = zeta_range(1) + (zeta_range(2) - zeta_range(1)) * rand(n_subjects, 1);
    
    results = struct();
    results.true_alphas = true_alphas;
    results.true_zetas = true_zetas;
    results.sigma_values = sigma_values;
    results.n_subjects = n_subjects;
    results.seed = seed;
    results.recovered_alphas = zeros(n_subjects, n_sigmas);
    results.recovered_zetas = zeros(n_subjects, n_sigmas);
    results.alpha_errors = zeros(n_subjects, n_sigmas);
    results.zeta_errors = zeros(n_subjects, n_sigmas);
    results.model_evidence = zeros(n_subjects, n_sigmas);
    results.model_fits = zeros(n_subjects, n_sigmas);
    results.alpha_prior_used = struct('muTheta', prior_settings.muTheta, ...
        'SigmaTheta', prior_settings.SigmaTheta);
    results.zeta_prior_used = struct('muPhi', prior_settings.muPhi, ...
        'SigmaPhi', prior_settings.SigmaPhi);
    
    for s = 1:n_subjects
        for j = 1:n_sigmas
            [recovered_alpha, recovered_zeta, model_evidence, model_fit] = ...
                run_single_two_param_recovery_with_prior(true_alphas(s), ...
                true_zetas(s), sigma_values(j), pairings, N, prior_settings);
            
            results.recovered_alphas(s, j) = recovered_alpha;
            results.recovered_zetas(s, j) = recovered_zeta;
            results.alpha_errors(s, j) = abs(recovered_alpha - true_alphas(s));
            results.zeta_errors(s, j) = abs(recovered_zeta - true_zetas(s));
            results.model_evidence(s, j) = model_evidence;
            results.model_fits(s, j) = model_fit;
        end
    end
end

%% Single subject recovery
function [recovered_alpha, recovered_zeta, model_evidence, model_fit] = ...
    run_single_two_param_recovery_with_prior(true_alpha, true_zeta, ...
    sigma, pairings, N, prior_settings)
    
    theta_true = log(true_alpha / (1 - true_alpha));
    phi_true = log(true_zeta);
    
    opt_sim = struct();
    opt_sim.priors.muX0 = 0.5;
    opt_sim.priors.SigmaX0 = 0;
    
    sigma_precision = 1 / (sigma^2);
    [y_sim_raw, ~] = VBA_simulate(N, @f_contingency_rw_single_param, ...
        @g_sigmoid, theta_true, phi_true, pairings, Inf, sigma_precision, opt_sim);
    y_sim = max(0, min(1, y_sim_raw));
    
    % Set up VBA inversion
    dim.n = 1; 
    dim.n_theta = 1; 
    dim.n_phi = 1; 
    dim.p = 1; 
    dim.n_t = N;
    
    opt = struct();
    opt.verbose = 0; 
    opt.DisplayWin = 0;
    opt.GnFigs = 0;
    opt.TolFun = 1e-4;
    
    % Set priors
    opt.priors.muTheta = prior_settings.muTheta;
    opt.priors.SigmaTheta = prior_settings.SigmaTheta;
    opt.priors.muPhi = prior_settings.muPhi;
    opt.priors.SigmaPhi = prior_settings.SigmaPhi;
    opt.priors.muX0 = 0.5; 
    opt.priors.SigmaX0 = 0;
    
    % Hyperpriors on observation noise
    [a_sigma, b_sigma] = VBA_guessHyperpriors(y_sim, [0.01, 0.99]);
    opt.priors.a_sigma = a_sigma;
    opt.priors.b_sigma = b_sigma;
    
    % Run VBA inversion
    [post, out] = VBA_NLStateSpaceModel(y_sim, pairings, ...
        @f_contingency_rw_single_param, @g_sigmoid, dim, opt);
    
    % Extract recovered parameters
    recovered_alpha = 1/(1+exp(-post.muTheta));
    recovered_zeta = exp(post.muPhi);
    model_evidence = out.F;
    
    if isfield(out.suffStat, 'muX')
        model_fit = corr(y_sim', out.suffStat.muX');
    else
        model_fit = NaN;
    end
end

%% VBA model functions
function [xnext, dfdx, dfdtheta] = f_contingency_rw_single_param(x, theta, u, ~)
    alpha = 1/(1+exp(-theta(1)));
    prediction_error = u(1) - x;
    xnext = x + alpha * prediction_error;
    if nargout > 1
        dfdx = 1 - alpha;
        dfdtheta = prediction_error * alpha * (1 - alpha);
    end
end

function [gx, dgdx, dgdphi] = g_sigmoid(x, phi, ~, ~)
    zeta = exp(phi(1));
    gx = 1./(1 + exp(-zeta * (x-0.5)));
    if nargout > 1
        s = gx;
        dgdx = zeta * s .* (1 - s);
        dgdphi = (x-0.5) .* s .* (1 - s) * zeta;
    end
end

%% Aggregate results across seeds
function aggregated_results = aggregate_across_seeds(all_seeds_results, ~)
    n_seeds = length(all_seeds_results);
    n_priors = length(all_seeds_results{1});
    n_sigmas = length(all_seeds_results{1}{1}.sigma_values);
    n_subjects = all_seeds_results{1}{1}.n_subjects;
    
    aggregated_results = cell(n_priors, 1);
    
    for p = 1:n_priors
        agg = struct();
        agg.alpha_prior_used = all_seeds_results{1}{p}.alpha_prior_used;
        agg.zeta_prior_used = all_seeds_results{1}{p}.zeta_prior_used;
        agg.sigma_values = all_seeds_results{1}{p}.sigma_values;
        agg.n_subjects = n_subjects;
        agg.all_seeds = cell(n_seeds, 1);
        
        alpha_corrs = zeros(n_seeds, n_sigmas);
        zeta_corrs = zeros(n_seeds, n_sigmas);
        alpha_errors = zeros(n_seeds, n_sigmas);
        zeta_errors = zeros(n_seeds, n_sigmas);
        error_corrs = zeros(n_seeds, n_sigmas);
        
        for seed_idx = 1:n_seeds
            agg.all_seeds{seed_idx} = all_seeds_results{seed_idx}{p};
            
            for j = 1:n_sigmas
                alpha_corrs(seed_idx, j) = corr(all_seeds_results{seed_idx}{p}.true_alphas, ...
                    all_seeds_results{seed_idx}{p}.recovered_alphas(:, j));
                zeta_corrs(seed_idx, j) = corr(all_seeds_results{seed_idx}{p}.true_zetas, ...
                    all_seeds_results{seed_idx}{p}.recovered_zetas(:, j));

                alpha_errors(seed_idx, j) = ...
                    mean(all_seeds_results{seed_idx}{p}.alpha_errors(:, j));
                zeta_errors(seed_idx, j) = ...
                    mean(all_seeds_results{seed_idx}{p}.zeta_errors(:, j));
                
                alpha_err = all_seeds_results{seed_idx}{p}.true_alphas - ...
                    all_seeds_results{seed_idx}{p}.recovered_alphas(:, j);
                zeta_err = all_seeds_results{seed_idx}{p}.true_zetas - ...
                    all_seeds_results{seed_idx}{p}.recovered_zetas(:, j);

                error_corrs(seed_idx, j) = corr(alpha_err, zeta_err);
            end
        end
        
        agg.alpha_corr_mean = mean(alpha_corrs, 1);
        agg.alpha_corr_std = std(alpha_corrs, 0, 1);
        agg.zeta_corr_mean = mean(zeta_corrs, 1);
        agg.zeta_corr_std = std(zeta_corrs, 0, 1);
        agg.alpha_error_mean = mean(alpha_errors, 1);
        agg.alpha_error_std = std(alpha_errors, 0, 1);
        agg.zeta_error_mean = mean(zeta_errors, 1);
        agg.zeta_error_std = std(zeta_errors, 0, 1);
        agg.error_corr_mean = mean(error_corrs, 1);
        agg.error_corr_std = std(error_corrs, 0, 1);
        agg.alpha_corrs_all_seeds = alpha_corrs;
        agg.zeta_corrs_all_seeds = zeta_corrs;
        agg.error_corrs_all_seeds = error_corrs;
        
        aggregated_results{p} = agg;
    end
end

%% Identify best priors
% Priors are evaluated here based on recovery correlations, correlations
% between estimation errors for the two parameters, and stability across
% runs
function best_priors = identify_best_priors_multiseed(all_results)
    fprintf('\n== BEST PRIOR IDENTIFICATION ==\n');

    n_priors = length(all_results);

    alpha_corr_overall = zeros(n_priors, 1);
    zeta_corr_overall = zeros(n_priors, 1);
    error_corr_overall = zeros(n_priors, 1);
    stability_score = zeros(n_priors, 1);

    for p = 1:n_priors
        alpha_corr_overall(p) = mean(all_results{p}.alpha_corr_mean);
        zeta_corr_overall(p) = mean(all_results{p}.zeta_corr_mean);
        error_corr_overall(p) = mean(abs(all_results{p}.error_corr_mean));
        stability_score(p) = 1 / (1 + mean(all_results{p}.alpha_corr_std) + ...
            mean(all_results{p}.zeta_corr_std));
    end

    valid_priors = error_corr_overall <= 0.3;
    fprintf('Priors meeting mean |error corr| <= 0.3: %d/%d\n', sum(valid_priors), n_priors);

    if any(valid_priors)
        valid_idx = find(valid_priors);
        [~, best_local] = max(alpha_corr_overall(valid_priors));
        best_prior = valid_idx(best_local);

        fprintf(['Best valid prior #%d: muTheta=%.2f, SigmaTheta=%.2f, muPhi=%.2f, ...' ...
            'SigmaPhi=%.2f\n'], best_prior, all_results{best_prior}.alpha_prior_used.muTheta, ...
            all_results{best_prior}.alpha_prior_used.SigmaTheta, ...
            all_results{best_prior}.zeta_prior_used.muPhi, ...
            all_results{best_prior}.zeta_prior_used.SigmaPhi);

        fprintf('  mean alpha r = %.3f\n', alpha_corr_overall(best_prior));
        fprintf('  mean zeta r = %.3f\n', zeta_corr_overall(best_prior));
        fprintf('  mean |err corr|= %.3f\n', error_corr_overall(best_prior));

        % Rank all valid priors by alpha recovery
        [~, order] = sort(alpha_corr_overall(valid_priors), 'descend');
        ranked_idx = valid_idx(order);

    else
        
        [~, best_prior] = min(error_corr_overall);

        fprintf('Fallback prior #%d: muTheta=%.2f, SigmaTheta=%.2f, muPhi=%.2f, SigmaPhi=%.2f\n', ...
            best_prior, ...
            all_results{best_prior}.alpha_prior_used.muTheta, ...
            all_results{best_prior}.alpha_prior_used.SigmaTheta, ...
            all_results{best_prior}.zeta_prior_used.muPhi, ...
            all_results{best_prior}.zeta_prior_used.SigmaPhi);

        fprintf('  mean alpha r   = %.3f\n', alpha_corr_overall(best_prior));
        fprintf('  mean zeta r    = %.3f\n', zeta_corr_overall(best_prior));
        fprintf('  mean |err corr|= %.3f\n', error_corr_overall(best_prior));

        ranked_idx = best_prior;
    end

    best_priors.ranked_indices = ranked_idx;
    best_priors.alpha_corr = alpha_corr_overall;
    best_priors.zeta_corr = zeta_corr_overall;
    best_priors.error_corr = error_corr_overall;
    best_priors.stability = stability_score;
    best_priors.best_prior_number = best_prior;

end

%% Analyze stability across seeds
function analyze_seed_stability(all_results, best_priors)
    fprintf('\n== SEED STABILITY ANALYSIS ==\n');

    p = best_priors.best_prior_number;
    max_alpha_std = max(all_results{p}.alpha_corr_std);
    max_zeta_std = max(all_results{p}.zeta_corr_std);

    fprintf('Best prior #%d: max alpha SD across seeds = %.3f, max zeta SD across seeds = %.3f\n', ...
        p, max_alpha_std, max_zeta_std);

    if max(max_alpha_std, max_zeta_std) < 0.02
        fprintf('Seed stability: VERY STABLE\n');
    elseif max(max_alpha_std, max_zeta_std) < 0.05
        fprintf('Seed stability: STABLE\n');
    else
        fprintf('Seed stability: VARIABLE\n');
    end
end

%% Analyze error correlation by zeta prior width
function analyze_error_correlation_by_sigmaPhi(all_results)
    fprintf('\n== ERROR CORRELATION BY ZETA PRIOR WIDTH ==\n');

    n_priors = length(all_results);
    sigmaPhi_vals = zeros(n_priors, 1);
    error_corrs = zeros(n_priors, 1);

    for p = 1:n_priors
        sigmaPhi_vals(p) = all_results{p}.zeta_prior_used.SigmaPhi;
        error_corrs(p) = mean(abs(all_results{p}.error_corr_mean));
    end

    [best_error_corr, best_idx] = min(error_corrs);

    fprintf('Lowest absolute error correlation: prior #%d, SigmaPhi=%.2f, |error corr|=%.3f\n', ...
        best_idx, sigmaPhi_vals(best_idx), best_error_corr);
end

%% Create correlation plots
function create_correlation_plots(all_results, best_priors)
    fprintf('\n== CREATING CORRELATION PLOTS ==\n');
    
    best_prior_idx = best_priors.best_prior_number;
    n_sigmas = length(all_results{best_prior_idx}.sigma_values);
    best_results = all_results{best_prior_idx}.all_seeds{1};
    
    figure('Position', [50, 50, 1400, 800]);
    
    for j = 1:n_sigmas
        subplot(2, n_sigmas, j);
        scatter(best_results.true_alphas, best_results.recovered_alphas(:, j), ...
                50, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.6);
        hold on;
        plot([0 1], [0 1], 'k--', 'LineWidth', 2);
        
        r_alpha = all_results{best_prior_idx}.alpha_corr_mean(j);
        r_std = all_results{best_prior_idx}.alpha_corr_std(j);
        
        xlabel('True Alpha', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Recovered Alpha', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('σ=%.3f, r=%.3f±%.3f', all_results{best_prior_idx}.sigma_values(j), ...
            r_alpha, r_std), 'FontSize', 12, 'FontWeight', 'bold');
        grid on; axis equal; xlim([0 1]); ylim([0 1]);
        
        subplot(2, n_sigmas, j + n_sigmas);
        scatter(best_results.true_zetas, best_results.recovered_zetas(:, j), 50, ...
            [0.8 0.2 0.4], 'filled', 'MarkerFaceAlpha', 0.6);
        hold on;
        
        z_max = max([best_results.true_zetas; best_results.recovered_zetas(:, j)]) * 1.1;
        plot([0 z_max], [0 z_max], 'k--', 'LineWidth', 2);
        
        r_zeta = all_results{best_prior_idx}.zeta_corr_mean(j);
        r_zeta_std = all_results{best_prior_idx}.zeta_corr_std(j);
        
        xlabel('True Zeta', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Recovered Zeta', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('σ=%.3f, r=%.3f±%.3f', all_results{best_prior_idx}.sigma_values(j), ...
            r_zeta, r_zeta_std), 'FontSize', 12, 'FontWeight', 'bold');
        grid on; axis equal;
    end
    
    sgtitle(sprintf('Best Prior #%d: θ=(%.2f,%.2f), φ=(%.2f,%.2f)', ...
        best_prior_idx, ...
        all_results{best_prior_idx}.alpha_prior_used.muTheta, ...
        all_results{best_prior_idx}.alpha_prior_used.SigmaTheta, ...
        all_results{best_prior_idx}.zeta_prior_used.muPhi, ...
        all_results{best_prior_idx}.zeta_prior_used.SigmaPhi), ...
        'FontSize', 14, 'FontWeight', 'bold');
end


%% Comprehensive convergence check
function check_convergence_comprehensive(all_results)
    fprintf('\n========================================\n');
    fprintf('CONVERGENCE CHECK\n');
    fprintf('========================================\n\n');
    
    n_subjects_actual = length(all_results{1}.all_seeds{1}.true_alphas);
    
    % Test sizes may be adjusted here
    test_sizes = [20, 40, 60, 80, 100];
    test_sizes = test_sizes(test_sizes <= n_subjects_actual);
    
    test_priors = [1, round(length(all_results)/2), length(all_results)];
    n_iterations = 20;
    sigma_idx = 2;
    
    fprintf('Testing convergence across sample sizes...\n');
    fprintf('Test sizes: [%s]\n', num2str(test_sizes));
    fprintf('Using σ=%.2f\n\n', all_results{1}.sigma_values(sigma_idx));
    
    for p = test_priors
        fprintf('Prior #%d:\n', p);
        
        true_alphas = all_results{p}.all_seeds{1}.true_alphas;
        recovered_alphas = all_results{p}.all_seeds{1}.recovered_alphas(:, sigma_idx);
        
        for s_idx = 1:length(test_sizes)
            n_subs = test_sizes(s_idx);

            corrs = zeros(n_iterations, 1);

            for iter = 1:n_iterations
                sample_idx = randperm(n_subjects_actual, n_subs);
                corrs(iter) = corr(true_alphas(sample_idx), recovered_alphas(sample_idx));
            end
            
            mean_corr = mean(corrs);
            std_corr = std(corrs);
            
            fprintf('  N=%3d: r=%.3f ± %.3f ', n_subs, mean_corr, std_corr);
            
            if std_corr < 0.02
                fprintf('[STABLE]\n');
            else
                fprintf('[UNSTABLE]\n');
            end
        end
        fprintf('\n');
    end
end
