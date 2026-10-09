%% ========================================================================
% Initialise storage for one simulation result table
% ========================================================================
%
% Full call chain:
%   s01_mcdt_finite_sample_calibration.m
%       -> runPairedFullGrid.m
%           -> initialiseResultColumns.m
%
% Preallocates the columns used to store one full grid of simulation
% results. The fields include:
%   - grid parameters (n, p, nu, alpha0, conflev)
%   - asymptotic and calibrated cutoffs, on squared and unsquared scales
%   - empirical sizes and size ratios
%   - determinant-bias bookkeeping
%% ========================================================================


function cols = initialiseResultColumns(nRows)

cols.n = zeros(nRows,1);
cols.p = zeros(nRows,1);
cols.nu = zeros(nRows,1);
cols.alpha0 = zeros(nRows,1);
cols.conflev = zeros(nRows,1);

% Cutoffs and calibration quantiles on the squared-distance scale
%
cols.q_asym_sq = nan(nRows,1);                  % asymptotic cutoff (squared)
cols.q_emp_avgquantile_sq = nan(nRows,1);       % Eq.14-15: average of per-rep quantiles (squared)
cols.q_emp_pooled_sq = nan(nRows,1);            % pooled sample quantile (squared)
cols.kappa_RD = nan(nRows,1);                   % ratio using avg-quantile: q_emp_avgquantile_sq / q_asym_sq
cols.q_RD_calibrated_sq = nan(nRows,1);         % calibrated cutoff (squared)

% Same quantities reported on the unsquared-distance scale
%
cols.q_asym = nan(nRows,1);
cols.q_emp_avgquantile = nan(nRows,1);
cols.q_emp_pooled = nan(nRows,1);
cols.q_RD_calibrated = nan(nRows,1);

% Empirical sizes and size ratios
%
cols.emp_size = nan(nRows,1);
cols.nominal_size = nan(nRows,1);
cols.sizeRatio = nan(nRows,1);

% Average-quantile calibration Eq. (14)/(15) 
%
cols.emp_size_calibrated = nan(nRows,1);
cols.sizeRatio_calibrated = nan(nRows,1);

% Pooled calibration
%
cols.emp_size_pooled = nan(nRows,1);
cols.sizeRatio_calibrated_pooled = nan(nRows,1);

% Replication-wise calibration
%
cols.emp_size_replicationwise = nan(nRows,1);
cols.sizeRatio_calibrated_replicationwise = nan(nRows,1);

% Whether the bias-determinant correction was available
%
cols.biasdet_found = true(nRows,1);
end
