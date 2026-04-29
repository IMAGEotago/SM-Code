%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BLT Behavioural Data Analysis                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------- %
% Author: Sophie Morris                                                   %
% Created: 19/11/2025                                                     %
% ----------------------------------------------------------------------- %
% Description: This script analyses the behavioural BLT data collected as %
% part of the 'Breathing and Anxiety' study, conducted collaboratively    %
% with the Department of Psychology, School of Sport, Physical Education  %
% and Sport Science, and the School of Pharmacy at the University of      %
% Otago. A timestamped analysis plan is available on                      %
% (https://github.com/IMAGEotago).                                        %
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up data directory and parallel processing
addpath("/Users/sophiemorris/Desktop/study/analysis/matlab/VBA-toolbox-master");
data = "/Users/sophiemorris/Desktop/study/analysis/matlab/data";
data_files = dir(fullfile(data, '*.mat'));
n_datasets = length(data_files);

%% Preallocate storage
participant_ids = strings(0,1);
i_alpha = [];
i_phi = [];
i_sigmaphi = [];
i_zeta = [];
i_theta = [];
i_sigmatheta = [];
i_certainty = [];
i_modelfit_rw = [];
i_modelfit_null = [];
i_freeenergy_rw = [];
i_freeenergy_null = [];
i_aic_rw = [];
i_aic_null = [];
i_bic_rw = [];
i_bic_null = [];

%% Loop through datasets
for i = 1:n_datasets
    file_path = fullfile(data_files(i).folder, data_files(i).name);
    loaded = load(file_path);

    data = loaded.data;
    params = loaded.params;

    cue = double(params.cue(:));
    pairings = double(params.pairings(:));
    pred = (double(data.pred_answer(:)) / 100)';

    y = nan(size(pred));
    u = nan(size(pred));

    % Convert to contingency space
    for n = 1:length(pred)
        if cue(n) == 1
            y(n) = pred(n);
            u(n) = pairings(n);
        elseif cue(n) == 2
            y(n) = 1 - pred(n);
            u(n) = pairings(n);
        end
    end

    % Exclusions
    invalid_trials = isnan(y) | y > 1 | y < 0;
    if sum(invalid_trials) > 0.2 * length(y)
        fprintf('Excluded (>20%% invalid): %s\n', data_files(i).name);
        continue;
    end

    % Replace invalid with neutral response
    y(invalid_trials) = 0.5;

    % Calculate average certainty
    t_certainty = abs(y - 0.5);
    m_certainty = mean(t_certainty(~invalid_trials));

    if size(y,1) > size(y,2), y = y'; end
    if size(u,1) > size(u,2), u = u'; end

    % Base dimensions
    dim = struct;
    dim.n = 1;
    dim.n_phi = 1;
    dim.n_theta = 1;
    dim.p = 1;
    dim.n_t = length(y);

    % Base options
    options = struct;
    options.isYout = invalid_trials(:)';
    options.verbose = 1;
    options.DisplayWin = false;
    options.priors.muPhi = 1.9;
    options.priors.SigmaPhi = 0.6;
    options.priors.muX0 = 0.5;
    options.priors.SigmaX0 = 0;

    % Hyperpriors on observation noise
    [a_sigma, b_sigma] = VBA_guessHyperpriors(y, [0.01, 0.99]);
    options.priors.a_sigma = a_sigma;
    options.priors.b_sigma = b_sigma;

    %% RW model (learning)
    options_rw = options;
    options_rw.priors.muTheta = -0.3;
    options_rw.priors.SigmaTheta = 0.5;

    [post_rw, out_rw] = VBA_NLStateSpaceModel(y, u, @f_rw, @g_sigmoid, ...
        dim, options_rw);

    %% Null model (α = 0)
    options_null = options;
    options_null.priors.muTheta = -10;
    options_null.priors.SigmaTheta = 0;

    [post_null, out_null] = VBA_NLStateSpaceModel(y, u, @f_rw, ...
        @g_sigmoid, dim, options_null);

    %% Store
    % Extract participant ID
    [~, name, ~] = fileparts(data_files(i).name);
    token = regexp(name, 'sub-(\d{4,5}(?:-\d)?)', 'tokens');
    if isempty(token)
        fprintf('Could not parse ID from %s\n', name);
        continue
    end
    sub_id = token{1}{1};
    participant_ids(end+1, 1) = string(sub_id);

    i_certainty(end+1) = m_certainty;

    % RW parameters
    i_theta(end+1) = post_rw.muTheta(1);
    i_sigmatheta(end+1) = post_rw.SigmaTheta(1);
    i_alpha(end+1) = 1/(1 + exp(-post_rw.muTheta(1)));

    i_phi(end+1) = post_rw.muPhi(1);
    i_sigmaphi(end+1) = post_rw.SigmaPhi(1);
    i_zeta(end+1) = exp(post_rw.muPhi(1));

    % RW model fit and goodness of fit metrics
    if isfield(out_rw, 'fit') && isfield(out_rw.fit, 'R2')
        thisR2_rw = out_rw.fit.R2;
    else
        y_obs = out_rw.y(:);
        y_hat = out_rw.suffStat.gx(:);
        valid = ~isnan(y_obs) & ~options.isYout(:);
        y_obs = y_obs(valid);
        y_hat = y_hat(valid);

        SSE = sum((y_obs - y_hat).^2);
        SST = sum((y_obs - mean(y_obs)).^2);
        thisR2_rw = 1 - SSE ./ SST;
    end

    i_modelfit_rw(end+1) = thisR2_rw;

    i_freeenergy_rw(end+1) = out_rw.F(1);
    i_aic_rw(end+1) = out_rw.fit.AIC(1);
    i_bic_rw(end+1) = out_rw.fit.BIC(1);

    % Null model fit and goodness of fit metrics
    if isfield(out_null, 'fit') && isfield(out_null.fit, 'R2')
        thisR2_null = out_null.fit.R2;
    else
        y_obs = out_null.y(:);
        y_hat = out_null.suffStat.gx(:);
        valid = ~isnan(y_obs) & ~options.isYout(:);
        y_obs = y_obs(valid);
        y_hat = y_hat(valid);

        SSE = sum((y_obs - y_hat).^2);
        SST = sum((y_obs - mean(y_obs)).^2);
        thisR2_null = 1 - SSE ./ SST;
    end
    i_modelfit_null(end+1) = thisR2_null;

    i_freeenergy_null(end+1) = out_null.F(1);
    i_aic_null(end+1) = out_null.fit.AIC(1);
    i_bic_null(end+1) = out_null.fit.BIC(1);
end

%% Per-subject model comparison: RW vs no-learning null
% Variational free energy difference is the primary model-comparison metric
deltaF  = i_freeenergy_rw' - i_freeenergy_null';    % RW - Null
deltaR2 = i_modelfit_rw' - i_modelfit_null';        % RW - Null

% Use deltaF directly rather than exp(deltaF), to avoid numerical overflow
F_threshold = 3;   % moderate evidence; use 5 for stricter evidence
retain_F = deltaF > F_threshold;

decision_str = repmat("Retain", numel(deltaF), 1);
decision_str(~retain_F) = "Exclude";

%% Per-file results table
T_subjects = table(participant_ids, i_alpha(:), i_theta(:), ...
    i_sigmatheta(:), i_zeta(:), i_phi(:), i_sigmaphi(:), ...
    i_certainty(:), ...
    i_modelfit_rw(:), i_modelfit_null(:), deltaR2(:), ...
    i_freeenergy_rw(:), i_freeenergy_null(:), deltaF(:), ...
    decision_str, ...
    'VariableNames', {'PPID','Alpha','muTheta','SigmaTheta','Zeta', ...
    'muPhi','SigmaPhi','AvgCertainty', ...
    'ModelFit_RW','ModelFit_Null','DeltaR2', ...
    'FreeEnergy_RW','FreeEnergy_Null','DeltaFreeEnergy', ...
    'OverallDecision'});

%% Print excluded files
fprintf('\nFiles excluded by free-energy comparison (DeltaF <= %.1f):\n', F_threshold);

excluded_idx = find(~retain_F);
if isempty(excluded_idx)
    fprintf('  None\n');
else
    for k = 1:numel(excluded_idx)
        iSub = excluded_idx(k);
        fprintf('  %s  [DeltaF = %.3f, DeltaR2 = %.3f]\n', ...
            participant_ids(iSub), deltaF(iSub), deltaR2(iSub));
    end
    fprintf('Total excluded files: %d / %d\n', numel(excluded_idx), numel(participant_ids));
end

%% Keep only complete participant pairs where both files are retained
base_id = regexprep(string(T_subjects.PPID), '-\d+$', '');
[G, base_groups] = findgroups(base_id);

% Number of files available for each participant base ID
n_files_per_base = splitapply(@numel, base_id, G);

% Whether all files for that participant were retained
row_included = retain_F;
all_included_per_base = splitapply(@all, row_included, G);

% Keep only participants with exactly 2 files and both retained
keep_base = (n_files_per_base == 2) & all_included_per_base;

% Expand group-level decision back to row level
keep_rows = keep_base(G);

%% Optional reporting for dropped participants
missing_pair_base = base_groups(n_files_per_base ~= 2);
failed_pair_base  = base_groups((n_files_per_base == 2) & ~all_included_per_base);

fprintf('\nParticipants dropped because they did not have both files:\n');
if isempty(missing_pair_base)
    fprintf('  None\n');
else
    for i = 1:numel(missing_pair_base)
        fprintf('  %s\n', string(missing_pair_base(i)));
    end
end

fprintf('\nParticipants dropped because at least one file failed the DeltaF criterion:\n');
if isempty(failed_pair_base)
    fprintf('  None\n');
else
    for i = 1:numel(failed_pair_base)
        fprintf('  %s\n', string(failed_pair_base(i)));
    end
end

%% Split into all vs included-only tables
T_subjects_all = T_subjects;
T_subjects_included = T_subjects(keep_rows, :);

%% Summary table for ALL rows
total_files_all = height(T_subjects_all);
n_excluded_all  = sum(T_subjects_all.OverallDecision == "Exclude");
n_retained_all  = total_files_all - n_excluded_all;

summaryNames = ["Mean"; "SD"];

T_summary_all = table(summaryNames, ...
    [mean(T_subjects_all.Alpha); std(T_subjects_all.Alpha)], ...
    [mean(T_subjects_all.muTheta); std(T_subjects_all.muTheta)], ...
    [mean(T_subjects_all.SigmaTheta); std(T_subjects_all.SigmaTheta)], ...
    [mean(T_subjects_all.Zeta); std(T_subjects_all.Zeta)], ...
    [mean(T_subjects_all.muPhi); std(T_subjects_all.muPhi)], ...
    [mean(T_subjects_all.SigmaPhi); std(T_subjects_all.SigmaPhi)], ...
    [mean(T_subjects_all.AvgCertainty); std(T_subjects_all.AvgCertainty)], ...
    [mean(T_subjects_all.ModelFit_RW); std(T_subjects_all.ModelFit_RW)], ...
    [mean(T_subjects_all.ModelFit_Null); std(T_subjects_all.ModelFit_Null)], ...
    [mean(T_subjects_all.DeltaR2); std(T_subjects_all.DeltaR2)], ...
    [mean(T_subjects_all.FreeEnergy_RW); std(T_subjects_all.FreeEnergy_RW)], ...
    [mean(T_subjects_all.FreeEnergy_Null); std(T_subjects_all.FreeEnergy_Null)], ...
    [mean(T_subjects_all.DeltaFreeEnergy); std(T_subjects_all.DeltaFreeEnergy)], ...
    [n_excluded_all; NaN], ...
    [n_retained_all; NaN], ...
    [total_files_all; NaN], ...
    'VariableNames', {'Statistic','Alpha','muTheta','SigmaTheta', ...
    'Zeta','muPhi','SigmaPhi','AvgCertainty', ...
    'ModelFit_RW','ModelFit_Null','DeltaR2', ...
    'FreeEnergy_RW','FreeEnergy_Null','DeltaFreeEnergy', ...
    'N_Excluded','N_Retained','N_Total'});

%% Summary table for INCLUDED paired rows only
total_files_included = height(T_subjects_included);

if total_files_included == 0
    T_summary_included = table(summaryNames, ...
        [NaN; NaN], [NaN; NaN], [NaN; NaN], [NaN; NaN], ...
        [NaN; NaN], [NaN; NaN], [NaN; NaN], [NaN; NaN], ...
        [NaN; NaN], [NaN; NaN], [NaN; NaN], [NaN; NaN], ...
        [NaN; NaN], [0; NaN], [0; NaN], [0; NaN], ...
        'VariableNames', {'Statistic','Alpha','muTheta','SigmaTheta', ...
        'Zeta','muPhi','SigmaPhi','AvgCertainty', ...
        'ModelFit_RW','ModelFit_Null','DeltaR2', ...
        'FreeEnergy_RW','FreeEnergy_Null','DeltaFreeEnergy', ...
        'N_Excluded','N_Retained','N_Total'});
else
    T_summary_included = table(summaryNames, ...
        [mean(T_subjects_included.Alpha); std(T_subjects_included.Alpha)], ...
        [mean(T_subjects_included.muTheta); std(T_subjects_included.muTheta)], ...
        [mean(T_subjects_included.SigmaTheta); std(T_subjects_included.SigmaTheta)], ...
        [mean(T_subjects_included.Zeta); std(T_subjects_included.Zeta)], ...
        [mean(T_subjects_included.muPhi); std(T_subjects_included.muPhi)], ...
        [mean(T_subjects_included.SigmaPhi); std(T_subjects_included.SigmaPhi)], ...
        [mean(T_subjects_included.AvgCertainty); std(T_subjects_included.AvgCertainty)], ...
        [mean(T_subjects_included.ModelFit_RW); std(T_subjects_included.ModelFit_RW)], ...
        [mean(T_subjects_included.ModelFit_Null); std(T_subjects_included.ModelFit_Null)], ...
        [mean(T_subjects_included.DeltaR2); std(T_subjects_included.DeltaR2)], ...
        [mean(T_subjects_included.FreeEnergy_RW); std(T_subjects_included.FreeEnergy_RW)], ...
        [mean(T_subjects_included.FreeEnergy_Null); std(T_subjects_included.FreeEnergy_Null)], ...
        [mean(T_subjects_included.DeltaFreeEnergy); std(T_subjects_included.DeltaFreeEnergy)], ...
        [0; NaN], ...
        [total_files_included; NaN], ...
        [total_files_included; NaN], ...
        'VariableNames', {'Statistic','Alpha','muTheta','SigmaTheta', ...
        'Zeta','muPhi','SigmaPhi','AvgCertainty', ...
        'ModelFit_RW','ModelFit_Null','DeltaR2', ...
        'FreeEnergy_RW','FreeEnergy_Null','DeltaFreeEnergy', ...
        'N_Excluded','N_Retained','N_Total'});
end

%% Write to Excel
writetable(T_subjects_all, 'blt_beh_data_all.xlsx', 'Sheet', 'Subjects');
writetable(T_summary_all,  'blt_beh_data_all.xlsx', 'Sheet', 'Summary');

writetable(T_subjects_included, 'blt_beh_data_included.xlsx', 'Sheet', 'Subjects');
writetable(T_summary_included,  'blt_beh_data_included.xlsx', 'Sheet', 'Summary');

fprintf('\nSaved %d rows to blt_beh_data_all.xlsx\n', height(T_subjects_all));
fprintf('Saved %d rows to blt_beh_data_included.xlsx\n', height(T_subjects_included));