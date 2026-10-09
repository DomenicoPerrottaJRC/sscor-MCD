function hf = plotSizeDiagnostics(CalibrationDatabase, p, alpha0Plot, method, plotMode)
%% ========================================================================
% plotSizeDiagnostics
% ========================================================================
%
% Plot empirical outlier-test size before and after finite-sample
% calibration, for a selected MCD correction method.
%
% SYNTAX
%   fh = plotSizeDiagnostics(CalibrationDatabase,p,alpha0Plot,method,plotMode)
%
% INPUTS
%   CalibrationDatabase
%       Table produced by the main simulation script. It contains the
%       empirical sizes, nominal size, and calibrated-size diagnostics
%       for the calibration constructions considered in the simulation.
%
%   p
%       Dimension of the multivariate observations to be plotted.
%
%   alpha0Plot
%       MCD breakdown-point parameter to be plotted.
%
%   method
%       MCD/determinant-correction method to be plotted, e.g.
%           'asymptotic'
%           'det_mean_root'
%           'det_root_mean_cov'
%
%   plotMode
%       Calibration construction to display:
%
%           'pooled'
%               Option A. The empirical calibration cutoff is obtained
%               from the pooled calibration distribution, i.e. all
%               observations from the calibration replications are
%               pooled before computing the empirical quantile.
%
%           'avg'
%               Option B. The empirical calibration cutoff is obtained
%               by first computing the conflev-quantile separately within
%               each calibration replication and then averaging these
%               replication-specific quantiles, corresponding to the
%               averaging construction in Equations (14)-(15).
%
%           'replicationwise'
%               Option C. The replication-specific empirical calibration
%               quantiles are retained separately and are applied to the
%               corresponding independent evaluation replications.
%
%           'all'
%               Plot the three calibration constructions together,
%               allowing a direct comparison of pooled, average, and
%               replication-wise calibration.
%
%       If omitted, 'avg' is used.
%
%
% PURPOSE
% -------
%
% The figure compares the empirical rejection probability of the MCD
% robust-distance outlier test with its nominal level before and after
% finite-sample calibration.
%
% The nominal test size is
%
%       alpha = 1 - conflev.
%
% The uncalibrated size is obtained by applying the asymptotic reference
% cutoff q_asym^2 to the simulated robust distances.
%
% For the calibrated diagnostics, the calibration cutoff is estimated
% from a set of Monte Carlo replications that is disjoint from the
% replications used for the final size evaluation.
%
% Therefore, the calibrated-size curves shown by this function represent
% an out-of-sample validation of the corresponding calibration
% construction. The observations used to estimate the calibration
% cutoff are not reused in the evaluation of the calibrated empirical
% size.
%
%
% THREE CALIBRATION CONSTRUCTIONS
% --------------------------------
%
% The function distinguishes three different ways of constructing the
% finite-sample calibration target.
%
% OPTION A: POOLED
%
% All squared robust distances from the calibration replications are
% pooled into a single empirical distribution. The calibration cutoff is
%
%       q_cal,pooled^2
%       = Q_conflev( {MD_ib^2 : b in B_cal, i=1,...,n} ).
%
% The corresponding size is evaluated independently on B_eval.
%
%
% OPTION B: AVERAGE
%
% For each calibration replication b, the empirical conflev-quantile is
% computed separately:
%
%       q_emp,b^2
%       = Q_conflev(MD_1b^2,...,MD_nb^2).
%
% The calibration target is then
%
%       q_cal,avg^2
%       = (1/B_cal) sum_b q_emp,b^2.
%
% This is the calibration construction corresponding to the
% average-of-replication-quantiles formulation in Equations (14)-(15).
%
% The resulting common cutoff is subsequently evaluated on the independent
% evaluation replications.
%
%
% OPTION C: REPLICATION-WISE
%
% The replication-specific calibration quantiles are retained:
%
%       q_cal,b^2 = q_emp,b^2,
%
% and are applied replication-wise to the corresponding independent
% evaluation replications according to the simulation design.
%
% This construction therefore does not collapse the calibration
% information into one pooled quantile or one averaged quantile.
%
%
% SIZE RATIO
% ----------
%
% The principal quantity plotted by this function is the empirical-size
% ratio
%
%       sizeRatio = empirical size / nominal size.
%
% Hence:
%
%       sizeRatio = 1
%
% corresponds to exact nominal calibration;
%
%       sizeRatio > 1
%
% indicates over-rejection (empirical size larger than nominal);
%
%       sizeRatio < 1
%
% indicates under-rejection (empirical size smaller than nominal).
%
% For the calibrated diagnostics, the corresponding quantities are
% computed from the independent evaluation sample:
%
%       sizeRatio_calibrated
%
%       sizeRatio_calibrated_pooled
%
%       sizeRatio_calibrated_replicationwise.
%
% The exact column used depends on the selected plotMode.
%
%
% INTERPRETATION
% --------------
%
% The purpose of the figure is not to establish that one calibration
% construction is universally preferable. Rather, it provides a direct
% finite-sample comparison of the size behaviour induced by the three
% alternative definitions of the empirical calibration target.
%
% In particular, differences between the three curves can arise because
% pooled calibration, averaging of replication-specific quantiles, and
% replication-wise calibration are mathematically different operations.
%
% A calibrated size ratio close to one indicates that the corresponding
% calibration procedure restores the empirical rejection probability
% towards its nominal value in the independent evaluation sample.
%
%
% IMPORTANT DISTINCTION
% ---------------------
%
% This function is a diagnostic plotting routine. It does not recompute
% the MCD estimator, empirical calibration quantiles, or calibration
% factors. Those quantities are produced by the simulation code and
% stored in CalibrationDatabase.
%
% The function therefore assumes that CalibrationDatabase has been
% generated by the current two-sample calibration/evaluation simulation
% design and that the relevant columns for the selected plotMode are
% already present.
%
% In particular, this function does NOT use the old pooledMD-based
% diagnostic calculation. The calibrated-size quantities plotted here
% must be those produced by the current simulation, where calibration
% and evaluation replications are separated.
%
%
% SCALE
% -----
%
% The underlying calibration is performed on squared robust distances
% (MD^2). The plotting routine displays the resulting empirical-size
% ratios and therefore does not alter the squared/unsquared distance
% convention.
%
%
% OUTPUT
% ------
%
%   fh
%       Handle to the generated MATLAB figure.
%
%% ========================================================================

%% ------------------------------------------------------------------------
% Defaults
% -------------------------------------------------------------------------

logscale = false;

if nargin < 3 || isempty(alpha0Plot)
    alpha0Plot = min(CalibrationDatabase.alpha0);
end

if nargin < 4 || isempty(method)
    method = 'det_mean_root';
end

if nargin < 5 || isempty(plotMode)
    plotMode = 'avg';
end

plotMode = lower(char(plotMode));

validModes = {'pooled','avg','replicationwise','all'};

if ~ismember(plotMode,validModes)
    error('plotSizeDiagnostics:invalidPlotMode', ...
        ['Unknown plotMode "%s". Use one of: ' ...
        '''pooled'', ''avg'', ''replicationwise'', ''all''.'], ...
        plotMode);
end

%% ------------------------------------------------------------------------
% Check required variables
% -------------------------------------------------------------------------

requiredVars = { ...
    'method', ...
    'n', ...
    'p', ...
    'nu', ...
    'alpha0', ...
    'sizeRatio', ...
    'q_RD_calibrated'};

% All three calibrated-size columns are needed for 'all'.
% For a single mode, only the relevant column is strictly required.

switch plotMode
    case 'avg'
        requiredVars{end+1} = 'sizeRatio_calibrated';

    case 'pooled'
        requiredVars{end+1} = 'sizeRatio_calibrated_pooled';

    case 'replicationwise'
        requiredVars{end+1} = 'sizeRatio_calibrated_replicationwise';

    case 'all'
        requiredVars = [requiredVars, ...
            {'sizeRatio_calibrated', ...
            'sizeRatio_calibrated_pooled', ...
            'sizeRatio_calibrated_replicationwise'}];
end

missingVars = setdiff(requiredVars,...
    CalibrationDatabase.Properties.VariableNames);

if ~isempty(missingVars)
    error('plotSizeDiagnostics:missingVariables', ...
        'CalibrationDatabase is missing: %s', ...
        strjoin(missingVars,', '));
end

%% ------------------------------------------------------------------------
% Select requested MCD correction construction
% -------------------------------------------------------------------------

T = CalibrationDatabase( ...
    CalibrationDatabase.p == p & ...
    CalibrationDatabase.alpha0 == alpha0Plot & ...
    strcmp(CalibrationDatabase.method,method), :);

if isempty(T)
    error(['No data found for p=%d, alpha0=%.4g, ' ...
        'method=%s.'], ...
        p, alpha0Plot, method);
end

nuGrid = unique(T.nu);
nCols = numel(nuGrid);

%% ------------------------------------------------------------------------
% Plotting defaults
% -------------------------------------------------------------------------

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

fs = 20;

set(groot,...
    'defaultAxesFontSize',fs,...
    'defaultTextFontSize',fs,...
    'defaultLegendFontSize',fs);

% Before calibration
colBefore = [0.40 0.40 0.40];

% Calibration curves
colPooled = [0.0000 0.4470 0.7410];
colAvg    = [0.8500 0.3250 0.0980];
colRep    = [0.4940 0.1840 0.5560];

lw = 2;
mkSize = 7;

%% ------------------------------------------------------------------------
% Determine which calibrated columns are plotted
% -------------------------------------------------------------------------

switch plotMode

    case 'pooled'

        calibratedVars = {'sizeRatio_calibrated_pooled'};
        calibratedLabels = {'after pooled calibration'};
        calibratedStyles = {'--'};
        calibratedMarkers = {'s'};
        calibratedColors = {colPooled};

    case 'avg'

        calibratedVars = {'sizeRatio_calibrated'};
        calibratedLabels = {'after average calibration'};
        calibratedStyles = {'--'};
        calibratedMarkers = {'s'};
        calibratedColors = {colAvg};

    case 'replicationwise'

        calibratedVars = {'sizeRatio_calibrated_replicationwise'};
        calibratedLabels = {'after replication-wise calibration'};
        calibratedStyles = {'--'};
        calibratedMarkers = {'d'};
        calibratedColors = {colRep};

    case 'all'

        calibratedVars = { ...
            'sizeRatio_calibrated_pooled', ...
            'sizeRatio_calibrated', ...
            'sizeRatio_calibrated_replicationwise'};

        calibratedLabels = { ...
            'after pooled calibration', ...  % (A)
            'after average calibration', ... % Eq.~(14)--(15) -  (B)
            'after replication-wise calibration'}; %(C)

        calibratedStyles = {'--','-.',':'};

        calibratedMarkers = {'s','o','d'};

        calibratedColors = { ...
            colPooled, ...
            colAvg, ...
            colRep};
end

%% ------------------------------------------------------------------------
% Determine common y-axis range
% -------------------------------------------------------------------------

allRatios = T.sizeRatio;
allRatios = [allRatios; 1];

for j = 1:numel(calibratedVars)
    allRatios = [allRatios; T.(calibratedVars{j})]; %#ok<AGROW>
end

allRatios = allRatios(isfinite(allRatios));

yMin = min(allRatios);
yMax = max(allRatios);

yPad = 0.05*(yMax-yMin);

if yPad == 0
    yPad = 0.05;
end

yLim = [yMin-yPad, yMax+yPad];

%% ------------------------------------------------------------------------
% Create figure
% -------------------------------------------------------------------------

if strcmp(plotMode,'all')

    modeTitle = 'A/B/C calibration comparison';

else

    modeTitle = calibratedLabels{1};

end

hf = figure( ...
    'Name', ...
    sprintf('Outlier-test size diagnostics: p=%d, alpha0=%.2f, %s', ...
    p,alpha0Plot,modeTitle));

%% ------------------------------------------------------------------------
% Plot each nu
% -------------------------------------------------------------------------

for iu = 1:nCols

    nu = nuGrid(iu);

    Tn = T(T.nu == nu,:);
    Tn = sortrows(Tn,'n');

    subplot(1,nCols,iu);
    hold on;

    % -------------------------------------------------------------
    % Before calibration
    % -------------------------------------------------------------

    plot(Tn.n,...
        Tn.sizeRatio,...
        'LineStyle','-',...
        'Marker','o',...
        'Color',colBefore,...
        'LineWidth',lw,...
        'MarkerSize',mkSize,...
        'DisplayName','before calibration');

    % -------------------------------------------------------------
    % Calibrated curves
    % -------------------------------------------------------------

    for j = 1:numel(calibratedVars)

        plot(Tn.n,...
            Tn.(calibratedVars{j}),...
            'LineStyle',calibratedStyles{j},...
            'Marker',calibratedMarkers{j},...
            'Color',calibratedColors{j},...
            'LineWidth',lw,...
            'MarkerSize',mkSize,...
            'DisplayName',calibratedLabels{j});

    end

    % -------------------------------------------------------------
    % Nominal size
    % -------------------------------------------------------------

    yline(1,'k:',...
        'HandleVisibility','off', 'LineWidth', 2.5);

    if logscale
        set(gca,'XScale','log');
    
        %Show an x-tick at every plotted sample size and use integer formatting
        xi = unique(Tn.n);
        xticks(xi);
        % Use plain numeric labels (avoid scientific notation for integers)
        if all(mod(xi,1)==0)
            xtickformat('%,.0f');
        end
    end

    ylim(yLim);

    xlabel('$n$','Interpreter','latex');

    % -------------------------------------------------------------
    % Display relevant calibrated cutoff
    % -------------------------------------------------------------
    %
    % q_RD_calibrated is the common stored cutoff associated with the
    % selected correction construction. It is retained here only as
    % a descriptive value in the title.

    [~,imax] = max(Tn.n);

    qCal = Tn.q_RD_calibrated(imax);

    title(sprintf('$p=%d$, $\\alpha_0=%.2f$, $\\nu=%d$ -- Cutoff for largest $n$: $\\sqrt{q_{cal}}=%.2f$', p,alpha0Plot,nu,qCal), 'Interpreter','latex');

    if iu == 1

        ylabel('empirical size / nominal size',...
            'Interpreter','latex');

        legend('Location','best',...
            'Interpreter','latex');

    end

end

%% ------------------------------------------------------------------------
% Console diagnostic
% -------------------------------------------------------------------------

fprintf('\nSize diagnostic summary\n');
fprintf('p = %d, alpha0 = %.4g, MCD method = %s\n',...
    p,alpha0Plot,method);
fprintf('calibration plot mode = %s\n\n',plotMode);

for iu = 1:nCols

    nu = nuGrid(iu);

    Tn = T(T.nu == nu,:);
    Tn = sortrows(Tn,'n');

    fprintf('nu = %g:\n',nu);

    fprintf('  before calibration : min = %.4f, max = %.4f\n',...
        min(Tn.sizeRatio),...
        max(Tn.sizeRatio));

    for j = 1:numel(calibratedVars)

        x = Tn.(calibratedVars{j});

        fprintf('  %-36s: min = %.4f, max = %.4f\n',...
            calibratedLabels{j},...
            min(x),...
            max(x));

    end

    fprintf('\n');

end

end
