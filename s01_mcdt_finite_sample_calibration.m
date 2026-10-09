%% ========================================================================
% s01_mcdt_finite_sample_calibration.m
% ========================================================================
%
% SIMULATION of MCD Robust-Distance Calibration Under Multivariate t
%
% Purpose: 
% Reproduces the Monte Carlo calibration results reported in the paper
% "Finite-sample calibration of robust covariance estimation under heavy tails"
%
% Workflow:
%   1. Estimate determinant-bias correction factors.
%   2. Run split-sample Monte Carlo simulations.
%   3. Estimate empirical robust-distance cutoffs.
%   4. Evaluate empirical test size on an independent evaluation sample.
%   5. Save calibration tables and diagnostic figures.
%
% Main Input:
%   - Grid: n_vec, p_vec, nu_vec, alpha0_vec
%   - B: Monte Carlo replications (even; split into calibration/evaluation)
%
% Main output:
%   - CalibrationDatabase (CSV/MAT), diagnostic figures
%     mcdt_calibration.csv
%
% Requirements:
%   - MATLAB
%   - FSDA toolbox
%   - Statistics Toolbox
%
% See README.md for the simulation design and parameter definitions.
%
% Notes:
%   - Detailed algorithmic description and derivations moved to README.md.
%   - Dependencies: mcd, msdcutoff (and Statistics Toolbox for chi2rnd if used).
%   - Split-sample design:
%       B total replications are generated and split into two equal halves:
%        - calibration replications:  1:B_cal   (used to compute calibration cutoffs)
%        - evaluation replications:   B_cal+1:B (used only to evaluate empirical sizes)
%       All calibration quantiles (pooled or per-rep) are computed from the
%       calibration half only. Evaluation of empirical size is always done on
%       the independent evaluation half to avoid circularity.
%
% Full call chain:
%   s01_mcdt_finite_sample_calibration.m
%       -> estimateBiasDetFactors.m
%           -> loadBiasDetTable.m
%       -> buildBiasDetFun.m
%           -> buildBiasDetFun_exact
%               -> lookupBiasDetExact
%       -> runPairedFullGrid.m
%           -> runPairedSimConfig.m
%               -> onePairedReplication.m
%                   -> verifyCorrectedMD.m
%       -> resultColumnsToTable.m
%       -> writeResultRow.m
%       -> initialiseResultColumns.m
%       -> plotKappaRDDiagnostics.m
%       -> plotDeterminantBiasDiagnostics.m
%       -> plotSizeDiagnostics.m


clear; clc; close all   % cleaning

seed = 20260711;     % Reproducibility
useParfor   = true;  % Parallelization
projectRoot = pwd;   % Your project folder: results will be saved there

%% ========================================================================
% SIMULATION SETTINGS
% =========================================================================

% Monte Carlo replications with calibration / evaluation split
B = 20;    
if mod(B,2) ~= 0
    error('For the split-sample calibration diagnostics, B must be even. Current B = %d.', B);
end
B_cal  = B/2;
B_eval = B/2;
fprintf('Monte Carlo replications:  B = %d\n', B);
fprintf('Calibration replications:  B_cal = %d\n', B_cal);
fprintf('Evaluation  replications:  B_eval = %d\n\n', B_eval);

% Single confidence level (scalar). Run separate scripts for other conflev.
conflev = 0.95;

% Simulation grid
n_vec       = 50:10:500;       
p_vec       = [2,5];             
nu_vec      = [3,5,7];           
alpha0_vec  = 0.05:0.05:0.50;

% Quick test override
quickTest.n_vec      = [50:10:80];  % [50];
quickTest.p_vec      = [2,5];  % [2,5];
quickTest.nu_vec     = [3,5];  % [3,5,7];
quickTest.alpha0_vec = [0.10, 0.25, 0.5];  % [0.10, 0.25, 0.5];
quickTest.B          = [20];  % 100;

% Which pooled vs. avg comparison to show in diagnostics. 
% - 'avg' (default, Eq. 14–15) collects per-rep calibration Quantiles from 
%    each replication and averages them (q_emp_avgquantile). 
% - 'pooled' uses all calibration observations together and computes one
%    pooled quantile (q_emp_pooled_sq) from the combined calibration sample.
% - 'replicationwise' pairs each calibration replication with the matching
%    evaluation replication and uses the per-rep calibration quantile directly.
%
plotPooledComparison = 'avg'; % 'avg' | 'pooled' | 'replicationwise' | 'all' 

% Quantile method used to compute calibration cutoffs from the sample.
% 'inclusive' and 'exclusive' control endpoint handling, 'midpoint' uses
% midpoint interpolation, 'approximate' uses a fast T-Digest estimate,
% and 'HarrellDavis' uses a smooth weighted-order-statistics estimator.
%
q_method             = 'inclusive'; % 'inclusive' | 'exclusive' | 'midpoint' | 'approximate' | 'HarrellDavis'

% Workflow options
options.saveRawDatabase = true;
options.saveDiagnostics = true;

%% ========================================================================
% BIAS_DET SETTINGS
% ========================================================================

% Bias-determinant correction workflow:
%   'load'               -> load an existing correction table/interpolant
%   'estimate'           -> estimate the correction table, save it, and stop
%   'estimate-and-apply' -> estimate the correction table and use it immediately
%
biasDetMode = 'estimate-and-apply';    % 'load' | 'estimate' | 'estimate-and-apply'

% Parameters used only when biasDetMode = 'load'
BIASDET_TABLE = table();          % optional in-memory table
BIASDET_TABLE_FILE   = '';        % optional .mat / .csv file
BIASDET_TABLE_FOLDER = 'results_radius_xBiasDet';  % fallback folder for loading tables
BIASDET_TABLE_MODE   = 'exact';   % 'exact' | 'interpolate'
BIASDET_INTERP = [];              % optional interpolant object
BIASDET_INTERP_FILE  = '';        % optional .mat file containing interpolant
% see also biasDetOutFile below, used when biasDetMode = 'estimate' or 'estimate-and-apply'

% Determinant-bias constructions to estimate.
% 'mean_root_det' computes the mean of the per-replication determinant factors.
% 'root_det_mean_cov' computes the determinant factor from the mean covariance.
BIASDET_CONSTRUCTIONS = {'mean_root_det'};   % {'mean_root_det'} | {'root_det_mean_cov'} | both

% Optional verification of MD vs covariance-level correction.
% Checks that applying the covariance-level correction yields the same
% outlier flags as multiplying MDs (and checks numerical accuracy).
VERIFY_BIASDET_CORRECTION = false;

%% Add project folders to the path

addpath(genpath(fullfile(projectRoot, 'src_srep')));
addpath(genpath(fullfile(projectRoot, 'src_core')));
addpath(genpath(fullfile(projectRoot, 'src_plots')));

%% Storage and output files

saveResults = true;

if saveResults
    resultsRoot = fullfile(projectRoot, 'REPLICABILITY_RESULTS');
    mkdir(resultsRoot);

    ts = datestr(now, 'yyyy-mm-dd_HHMMSS');
    results_folder = fullfile(resultsRoot, ['calibration_' ts]);
    mkdir(results_folder);
end

outMatFile     = fullfile(results_folder, 'mcd_heavytails_RD_calibration_results.mat');
outCsvFile     = fullfile(results_folder, 'mcd_heavytails_RD_calibration.csv');
outDiagFile    = fullfile(results_folder, 'mcd_heavytails_RD_calibration_diagnostics.csv');
biasDetOutFile = fullfile(results_folder, 'BiasDet_table.mat');


%% Apply quick-test overrides
if ~isempty(quickTest.n_vec)
    n_vec = quickTest.n_vec;
end
if ~isempty(quickTest.p_vec)
    p_vec = quickTest.p_vec;
end
if ~isempty(quickTest.nu_vec)
    nu_vec = quickTest.nu_vec;
end
if ~isempty(quickTest.alpha0_vec)
    alpha0_vec = quickTest.alpha0_vec;
end
if ~isempty(quickTest.B)
    B = quickTest.B;
end

rng(seed,'twister');

%% ========================================================================
% Sanity check: verify that the MCD robust-distance cutoff matches the
% expected multivariate t quantile on a small random sample
% ========================================================================
n0 = 100;
p0 = p_vec(1);
nu0 = nu_vec(1);
conflev0 = conflev;
X0 = randn(n0,p0).*sqrt(nu0./chi2rnd(nu0,n0,1));
RAW0 = mcd(X0,...
    'modelT',nu0,...
    'bdp',alpha0_vec(1),...
    'conflev',conflev0,...
    'smallsamplecor',false,...
    'nsamp',500,...
    'refsteps',10,...
    'plots',0,...
    'msg',0,...
    'nocheck',1);
cutoffT0 = msdcutoff(conflev0,p0,nu0);
assert(numel(RAW0.outliers)==sum(RAW0.md>cutoffT0),...
    'Sanity check failed');
fprintf('Sanity check passed.\n\n');

%% ========================================================================
%% Build determinant correction functions
% ========================================================================
% Set up the correction-factor functions used to adjust the MCD distances.
% Depending on biasDetMode, the factors are either:
%   - loaded from an existing table/interpolant ('load'), or
%   - estimated from Monte Carlo and then applied ('estimate' / 'estimate-and-apply').
% The two functions correspond to the two determinant-bias constructions:
%   - mean_root_det
%   - root_det_mean_cov

switch biasDetMode
    case 'load'
        meanRootDetFactorFun = buildBiasDetFun(...
            'table', BIASDET_TABLE, BIASDET_TABLE_FILE, BIASDET_TABLE_FOLDER, BIASDET_TABLE_MODE, ...
            'mean_root_det', BIASDET_INTERP, BIASDET_INTERP_FILE);
        rootDetMeanCovFactorFun = buildBiasDetFun(...
            'table', BIASDET_TABLE, BIASDET_TABLE_FILE, BIASDET_TABLE_FOLDER, BIASDET_TABLE_MODE, ...
            'root_det_mean_cov', BIASDET_INTERP, BIASDET_INTERP_FILE);

    case {'estimate','estimate-and-apply'}
        fprintf('Estimating BiasDet table (this may be expensive)...\n');
        BiasDetTable = estimateBiasDetFactors(n_vec, p_vec, nu_vec, alpha0_vec, B, ...
            useParfor, biasDetOutFile, conflev);
        save(biasDetOutFile, 'BiasDetTable');
        fprintf('BiasDet table saved to %s\n', biasDetOutFile);

        meanRootDetFactorFun    = buildBiasDetFun('table', BiasDetTable, '', '', 'exact', 'mean_root_det', [], '');
        rootDetMeanCovFactorFun = buildBiasDetFun('table', BiasDetTable, '', '', 'exact', 'root_det_mean_cov', [], '');

        if strcmp(biasDetMode,'estimate')
            fprintf('BiasDet estimation complete. Re-run with biasDetMode=''load'' to apply corrections.\n');
            return
        end

    otherwise
        error('Unknown biasDetMode: %s', biasDetMode);
end

%% ========================================================================
% RUN MONTE CARLO SIMULATION
% ========================================================================
%
% Main simulation step.
% This calls:
%   s01_mcdt_finite_sample_calibration.m
%       -> runPairedFullGrid.m
%           -> runPairedSimConfig.m
%               -> onePairedReplication.m
%
% At this stage the script:
%   - loops over the full (n, p, nu, alpha0) grid,
%   - computes the asymptotic cutoff,
%   - applies the determinant-bias correction factors,
%   - runs the split-sample paired replications,
%   - collects the baseline and calibrated empirical sizes.
%

[ResultsBaseline, ResultsDetMeanRoot, ResultsDetRootMeanCov] = ...
    runPairedFullGrid(...
        n_vec,...
        p_vec,...
        nu_vec,...
        alpha0_vec,...
        conflev,...
        B,...
        meanRootDetFactorFun,...
        rootDetMeanCovFactorFun,...
        useParfor,...
        VERIFY_BIASDET_CORRECTION,...
        q_method);
disp(ResultsBaseline);
disp(ResultsDetMeanRoot);
disp(ResultsDetRootMeanCov);


%% ========================================================================
% Combine results and create calibration database
% ========================================================================
ResultsAll = [...
    addMethodColumn(ResultsBaseline,'asymptotic'); ...
    addMethodColumn(ResultsDetMeanRoot,'det_mean_root'); ...
    addMethodColumn(ResultsDetRootMeanCov,'det_root_mean_cov')];
CalibrationDatabase = ResultsAll(:,...
    {'method',...
    'n',...
    'p',...
    'nu',...
    'alpha0',...
    'conflev',...
    'q_asym_sq',...
    'q_emp_avgquantile_sq',...
    'q_emp_pooled_sq',...
    'kappa_RD',...
    'q_RD_calibrated_sq',...
    'q_asym',...
    'q_emp_avgquantile',...
    'q_emp_pooled',...
    'q_RD_calibrated',...
    'emp_size',...
    'nominal_size',...
    'sizeRatio',...
    'emp_size_calibrated',...
    'sizeRatio_calibrated',...
    'emp_size_pooled',...
    'sizeRatio_calibrated_pooled',...
    'emp_size_replicationwise',...
    'sizeRatio_calibrated_replicationwise',...
    'biasdet_found'});


%% ========================================================================
% Simulation metadata
% ========================================================================
SimulationInfo.B = B;
SimulationInfo.seed = seed;
SimulationInfo.date = datestr(now);
SimulationInfo.version = ...
    'sim_MCDt_correction_revised_patch1';
SimulationInfo.options = options;

%% ========================================================================
% Save results
% ========================================================================
if options.saveRawDatabase
    save(outMatFile,...
        'ResultsBaseline',...
        'ResultsDetMeanRoot',...
        'ResultsDetRootMeanCov',...
        'ResultsAll',...
        'CalibrationDatabase',...
        'SimulationInfo',...
        'n_vec',...
        'p_vec',...
        'nu_vec',...
        'alpha0_vec',...
        'conflev');
end

writetable(CalibrationDatabase,outCsvFile);

fprintf('\nResults saved in:\n%s\n\n',results_folder);

%% ========================================================================
% Diagnostic plots
% ========================================================================
if options.saveDiagnostics
    set(groot,'defaultTextInterpreter','latex');
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    fs = 14;
    set(groot,...
        'defaultAxesFontSize',fs,...
        'defaultTextFontSize',fs,...
        'defaultLegendFontSize',fs);

    % Section 4.3: kappa_RD (plotMode: 'avg'|'pooled'|'both')
    fh1 = plotKappaRDDiagnostics(CalibrationDatabase, p_vec(1), alpha0_vec(1), plotPooledComparison);

    % Section 4.1 & 4.2: determinant-bias diagnostics (plotMode supported but determinant plots unaffected)
    fh2 = plotDeterminantBiasDiagnostics(BiasDetTable.mean_root, BiasDetTable.root_mean_cov, ...
        p_vec(1), alpha0_vec(1));

    % Outlier-test size before/after calibration for the chosen correction method.
    % use CalibrationDatabase (contains emp_size_pooled / sizeRatio_calibrated_pooled)
    fh3 = plotSizeDiagnostics(CalibrationDatabase, p_vec(1), alpha0_vec(1), 'det_mean_root', plotPooledComparison); 
    figHandles = [fh1, fh2, fh3];


    desired = { ...
        struct('handle',fh1,'fname','f1_kappa_RD','suptitle','Distance-based finite-sample calibration'), ...
        struct('handle',fh2,'fname','f2_determinant_calibration','suptitle','Residual finite-sample scatter discrepancy and determinant-based calibration factor'), ...
        struct('handle',fh3,'fname','f3_test_size_calibration','suptitle','Empirical outlier-test size (relative to nominal)') ...
        };

    for k = 1:numel(desired)
        entry = desired{k};
        f = entry.handle;
        f.Name = entry.fname;
    
        % Force undocked window and set a fixed figure size / aspect ratio
        try
            f.WindowStyle = 'normal';        % undock (makes standalone window)
        catch
            % older MATLAB versions may not support WindowStyle on figure handle
            set(f,'WindowStyle','normal');
        end
        % Choose a sensible size: [left bottom width height] in inches
        % Adjust width/height to taste (here 10 by 6 inches -> 5:3 aspect)
        try
            set(f,'Units','inches','Position',[1 1 12 5]);
        catch
            set(f,'Units','pixels'); % fallback; don't change size if fails
        end
        drawnow; pause(0.1); % ensure graphics updated and undocked
    
        pngFile = fullfile(results_folder, [entry.fname '.png']);
        epsFile = fullfile(results_folder, [entry.fname '.eps']);
        figFile = fullfile(results_folder, [entry.fname '.fig']);
    
        % Export raster PNG at high resolution
        f.WindowStyle = 'normal';
        set(f,'Units','inches'); drawnow; % ,'Position',[1 1 15 5]
        print(f, '-dpng', '-r300', pngFile);
    
        % Export vector EPS (use ContentType 'vector' where supported)
        set(f,'Renderer','painters'); drawnow;
        print(f, '-depsc2', epsFile);   
    
        % Also save MATLAB figure for later interactive editing
        savefig(f, figFile);
    end

end

%% quick consistency check on squared-unsquared distances

row = CalibrationDatabase(1,:);

cutoff_sq = msdcutoff(row.conflev , row.p , row.nu);

assert(abs(cutoff_sq - row.q_asym_sq) < 1e-12, ...
    'q_asym_sq mismatch');

assert(isfinite(row.q_emp_pooled_sq), ...
    'q_emp_pooled_sq is not available.');

assert(abs(row.kappa_RD - ...
    (row.q_emp_avgquantile_sq / row.q_asym_sq)) < 1e-12, ...
    'kappa_RD mismatch');

fprintf('Consistency checks passed for row 1 (squared-distance scale confirmed).\n');

assert(abs(sqrt(row.q_emp_pooled_sq) - row.q_emp_pooled) < 1e-12, ...
    'q_emp_pooled (unsquared) mismatch');

assert(abs(sqrt(cutoff_sq) - row.q_asym) < 1e-12, ...
    'q_asym (unsquared) mismatch');

assert(abs(row.kappa_RD - ...
    (row.q_emp_avgquantile_sq / row.q_asym_sq)) < 1e-12, ...
    'kappa_RD mismatch');

fprintf('Consistency checks passed for row 1 (unsquared reporting matches computed values).\n');


%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function Results = addMethodColumn(Results,method)
Results.method = repmat({method},height(Results),1);
Results = movevars(Results,'method','Before','n');
end

