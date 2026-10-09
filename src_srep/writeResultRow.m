%% ========================================================================
% Write one simulation result row
% ========================================================================
%
% Full call chain:
%   s01_mcdt_finite_sample_calibration.m
%       -> runPairedFullGrid.m
%           -> writeResultRow.m
%
% Fills one preallocated row with the grid parameters, calibration
% quantiles, empirical sizes, and determinant-bias bookkeeping for a
% single (n, p, nu, alpha0) configuration.
%% ========================================================================

function cols = writeResultRow(...
    cols,...
    row,...
    n,...
    p,...
    nu,...
    alpha0,...
    conflev,...
    cutoff_sq,...
    calib,...
    evalInfo,...
    found)

cols.n(row) = n;
cols.p(row) = p;
cols.nu(row) = nu;
cols.alpha0(row) = alpha0;
cols.conflev(row) = conflev;

% ------------------------------------------------------------------------
% Asymptotic cutoff
% -------------------------------------------------------------------------

cols.q_asym_sq(row) = cutoff_sq;

% ------------------------------------------------------------------------
% Calibration quantiles
%
% Eq. (14):
%
% q_emp^2 = (1/B_cal) sum_b q_emp^{2,(b)}
%
% Here the q_emp^{2,(b)} are computed ONLY in the calibration sample.
%
% -------------------------------------------------------------------------

%  average of per-rep calibration quantiles (Eq.14)
cols.q_emp_avgquantile_sq(row) = calib.q_avg_sq;

%  pooled quantile computed on calibration half
cols.q_emp_pooled_sq(row)      = calib.q_pooled_sq;

% ------------------------------------------------------------------------
% Eq. (15): kappa_RD
%
% The paper's definition corresponds to the average-of-replication
% quantiles.
% -------------------------------------------------------------------------

cols.kappa_RD(row) = ...
    calib.q_avg_sq / cutoff_sq;

% ------------------------------------------------------------------------
% Eq. (16): calibrated cutoff
%
% q_cal^2 = kappa_RD q_asym^2
%
% Algebraically this equals q_emp_avgquantile^2.
% -------------------------------------------------------------------------

cols.q_RD_calibrated_sq(row) = ...
    cols.kappa_RD(row) * cutoff_sq;

% ------------------------------------------------------------------------
% Unsquared reporting
% -------------------------------------------------------------------------

cols.q_asym(row) = sqrt(cols.q_asym_sq(row));

cols.q_emp_avgquantile(row) = ...
    sqrt(cols.q_emp_avgquantile_sq(row));

cols.q_emp_pooled(row) = ...
    sqrt(cols.q_emp_pooled_sq(row));

cols.q_RD_calibrated(row) = ...
    sqrt(cols.q_RD_calibrated_sq(row));

% ------------------------------------------------------------------------
% Nominal and asymptotic empirical size
%
% IMPORTANT:
% size is now evaluated ONLY on the independent evaluation sample.
% -------------------------------------------------------------------------

%  empirical rejection prob evaluated on the independent evaluation half
cols.emp_size(row) = ...
    evalInfo.size_asym;

cols.nominal_size(row) = ...
    1-conflev;

%  emp_size_* / nominal_size
cols.sizeRatio(row) = ...
    cols.emp_size(row) / cols.nominal_size(row);

% ------------------------------------------------------------------------
% Calibrated empirical size: Eq. (14)/(15)
%
% The average-quantile calibration cutoff is determined from the
% calibration sample and evaluated on the independent evaluation sample.
% -------------------------------------------------------------------------

cols.emp_size_calibrated(row) = ...
    evalInfo.size_avg;

cols.sizeRatio_calibrated(row) = ...
    cols.emp_size_calibrated(row) / cols.nominal_size(row);

% ------------------------------------------------------------------------
% Additional diagnostics for pooled and replication-wise calibration
% -------------------------------------------------------------------------

cols.emp_size_pooled(row) = ...
    evalInfo.size_pooled;

cols.sizeRatio_calibrated_pooled(row) = ...
    evalInfo.size_pooled / cols.nominal_size(row);

cols.emp_size_replicationwise(row) = ...
    evalInfo.size_replicationwise;

cols.sizeRatio_calibrated_replicationwise(row) = ...
    evalInfo.size_replicationwise / cols.nominal_size(row);

% ------------------------------------------------------------------------
% Store whether determinant correction was available
% -------------------------------------------------------------------------

cols.biasdet_found(row) = found;
end
