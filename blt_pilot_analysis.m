function pilot_results = blt_pilot_analysis()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BLT Pilot Data Analysis                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------- %
% Author: Sophie Morris                                                   %
% Created: 18/03/2026                                                     %
% ----------------------------------------------------------------------- %
% Description: This script analyses the behavioural BLT pilot data        %
% collected as part of the 'Breathing and Anxiety' study, conducted       %
% collaboratively with the Department of Psychology, School of Sport,     %
% Physical Education and Sport Science, and the School of Pharmacy at the %
% Universty of Otago. A timestamped analysis plan is available on         %
% (https://github.com/IMAGEotago).                                        %
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath("/Users/sophiemorris/Desktop/study/analysis/matlab/VBA-toolbox-master"));
addpath(genpath("/Users/sophiemorris/Desktop/study/analysis/matlab"));
pairings = load("/Users/sophiemorris/Desktop/study/analysis/matlab/blt_pairings.mat");
pilot_dir = "/Users/sophiemorris/Desktop/study/analysis/matlab/data-pilots";

%% Load pilot data
pilot_files = dir(fullfile(pilot_dir, '*.mat'));
n_pilots = length(pilot_files);
pilot_data = cell(n_pilots, 1);

for i = 1:n_pilots
    tmp = load(fullfile(pilot_dir, pilot_files(i).name));
    tmp.source_filename = pilot_files(i).name;
    pilot_data{i} = tmp;
end

%% Fit data
pilot_results = fit_pilots(pilot_data, pairings);

trim_percent = 20;

alpha_trimmed = trimmean(pilot_results.alphas, trim_percent);
alpha_median = median(pilot_results.alphas);
alpha_mean = mean(pilot_results.alphas);

theta_trimmed = trimmean(pilot_results.thetas, trim_percent);
theta_median = median(pilot_results.thetas);
theta_mean = mean(pilot_results.thetas);
theta_std_trimmed = std_trimmed(pilot_results.thetas, trim_percent);
theta_std_median = mad(pilot_results.thetas) * 1.4826;  % MAD (robust std)
theta_std_mean = std(pilot_results.thetas);

zeta_trimmed = trimmean(pilot_results.zetas, trim_percent);
zeta_median = median(pilot_results.zetas);
zeta_mean = mean(pilot_results.zetas);

phi_trimmed = trimmean(pilot_results.phis, trim_percent);
phi_median = median(pilot_results.phis);
phi_mean = mean(pilot_results.phis);
phi_std_trimmed = std_trimmed(pilot_results.phis, trim_percent);
phi_std_median = mad(pilot_results.phis) * 1.4826;  % MAD (robust std)
phi_std_mean = std(pilot_results.phis);

sigma_results = struct();
sigma_results.individual_sigmas = pilot_results.sigmas;
sigma_results.trimmed_mean = trimmean(pilot_results.sigmas, trim_percent);
sigma_results.median_sigma = median(pilot_results.sigmas);
sigma_results.mean_sigma = mean(pilot_results.sigmas);
sigma_results.std_sigma = std(pilot_results.sigmas);
sigma_results.min_sigma = min(pilot_results.sigmas);
sigma_results.max_sigma = max(pilot_results.sigmas);

% Save
pilot_results.statistics = struct();
pilot_results.statistics.alpha.trimmed_mean = alpha_trimmed;
pilot_results.statistics.alpha.median = alpha_median;
pilot_results.statistics.alpha.mean = alpha_mean;

pilot_results.statistics.theta.trimmed_mean = theta_trimmed;
pilot_results.statistics.theta.median = theta_median;
pilot_results.statistics.theta.mean = theta_mean;
pilot_results.statistics.theta.trimmed_sd = theta_std_trimmed;
pilot_results.statistics.theta.mad = theta_std_median;
pilot_results.statistics.theta.sd = theta_std_mean;

pilot_results.statistics.zeta.trimmed_mean = zeta_trimmed;
pilot_results.statistics.zeta.median = zeta_median;
pilot_results.statistics.zeta.mean = zeta_mean;

pilot_results.statistics.phi.trimmed_mean = phi_trimmed;
pilot_results.statistics.phi.median = phi_median;
pilot_results.statistics.phi.mean = phi_mean;
pilot_results.statistics.phi.trimmed_sd = phi_std_trimmed;
pilot_results.statistics.phi.mad = phi_std_median;
pilot_results.statistics.phi.sd = phi_std_mean;

pilot_results.statistics.sigma_results = sigma_results;

save('pilot_analysis.mat', 'pilot_results');

end

%% Functions
%% Function: Fit pilots
function results = fit_pilots(pilot_data, ~)
    n_pilots = length(pilot_data);
    
    results.alphas = [];
    results.zetas = []; 
    results.thetas = [];
    results.phis = [];
    results.R2 = [];
    results.alpha_CI_lower = [];
    results.alpha_CI_upper = [];
    results.excluded_subjects = [];
    results.subject_ids = [];
    results.sigmas = [];

    results.free_energy_per_subject = [];
    
    dim.n = 1;
    dim.n_theta = 1;
    dim.n_phi = 1;
    dim.p = 1;
    
    n_excluded_invalid = 0;
    n_excluded_learning = 0;
    
    alpha_threshold = 0.05; 
    
    for i = 1:n_pilots
        fprintf('Fitting pilot %d/%d... ', i, n_pilots);
        
        data = pilot_data{i}.data;
        params = pilot_data{i}.params;
        
        cue = double(params.cue(:));
        pairings = double(params.pairings(:));
        pred = (double(data.pred_answer(:)) / 100)';
        
        y = nan(size(pred));
        u = nan(size(pred));
        
        for n = 1:length(pred)
            if cue(n) == 1
                y(n) = pred(n);
                u(n) = pairings(n);
            elseif cue(n) == 2
                y(n) = 1 - pred(n);
                u(n) = pairings(n);
            end
        end
        
        invalid_trials = isnan(y) | y > 1 | y < 0;

        if sum(invalid_trials) > 0.2 * length(y)
            fprintf('EXCLUDED (>20%% invalid)\n');
            n_excluded_invalid = n_excluded_invalid + 1;
            results.excluded_subjects(end+1) = i;
            continue;
        end
        
        y(invalid_trials) = 0.5;
        
        if size(y, 1) > size(y, 2), y = y'; end
        if size(u, 1) > size(u, 2), u = u'; end
        
        dim.n_t = length(y);
        
        opt = struct();
        opt.verbose = 0; 
        opt.DisplayWin = 0;
        % Initially, theta and phi priors are centred on weakly-informative,
        % theoretically neutral values
        opt.priors.muTheta = 0;
        opt.priors.SigmaTheta = 1;
        opt.priors.muPhi = 0;
        opt.priors.SigmaPhi = 1;
        opt.priors.muX0 = 0.5; 
        opt.priors.SigmaX0 = 0;
        
        [a_sigma, b_sigma] = VBA_guessHyperpriors(y, [0.01, 0.99]);
        opt.priors.a_sigma = a_sigma;
        opt.priors.b_sigma = b_sigma;
        opt.isYout = invalid_trials;
        
        [post, out] = VBA_NLStateSpaceModel(y, u, @f_rw, @g_sigmoid, dim, opt);
        
        precision = post.a_sigma / post.b_sigma;
        sigma = sqrt(1 / precision);

        % Exclusion criterion
        theta_mean = post.muTheta;
        alpha_mean = 1/(1+exp(-theta_mean));
        
        zeta_mean = exp(post.muPhi);
        
        if alpha_mean < alpha_threshold
            fprintf('Excluded (α = %.2f, lower < %.2f)\n', alpha_mean, alpha_threshold);
            n_excluded_learning = n_excluded_learning + 1;
            results.excluded_subjects(end+1) = i;
            continue;
        end
        
        % Store results
        results.thetas(end+1) = post.muTheta;
        results.phis(end+1) = post.muPhi;
        results.alphas(end+1) = alpha_mean;
        results.zetas(end+1) = zeta_mean;
        results.R2(end+1) = out.fit.R2;
        results.sigmas(end+1) = sigma;

        if isfield(pilot_data{i}, 'source_filename')
            subj_id = pilot_data{i}.source_filename;
        elseif isfield(pilot_data{i}, 'id')
            subj_id = pilot_data{i}.id;
        else
            subj_id = sprintf('subj_%03d', i);
        end

        results.subject_ids{end+1} = subj_id;

        results.free_energy_per_subject(end+1) = out.F;

        fprintf('α=%.3f, ζ=%.3f, R²=%.3f\n', ...
            alpha_mean, zeta_mean, out.fit.R2);
    end
    
    fprintf('\nTotal pilots: %d\n', n_pilots);
    fprintf('Excluded (invalid trials): %d\n', n_excluded_invalid);
    fprintf('Excluded (α < %.2f): %d\n', alpha_threshold, n_excluded_learning);
    fprintf('Included in analysis: %d\n', length(results.alphas));    
end    

%% Function: Calculate robust SD
function trimmed_sd = std_trimmed(data, trim_percent)
    n = length(data);
    trim_count = floor(n * trim_percent / 200);
    
    sorted_data = sort(data);
    trimmed_data = sorted_data(trim_count+1:end-trim_count);
    
    trimmed_sd = std(trimmed_data);
end
