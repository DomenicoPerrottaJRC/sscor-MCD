function fh = plotKappaRDDiagnostics(CalibrationDatabase, p, alpha0Plot, plotMode)
% plotKappaRDDiagnostics  Section 4.3 diagnostic: the distance-based
% finite-sample calibration factor kappa_RD (Eq. 16).
%
% Inputs:
%   CalibrationDatabase : table with kappa_RD computed from the
%                         average-of-per-rep quantiles (q_emp_avgquantile_sq)
%   p, alpha0Plot        : select rows to plot (alpha0Plot defaults to min available)
%   plotMode (optional)  : 'avg' | 'pooled' | 'all'  (default: 'avg')

logscale = false;

% optional: whether to show the second determinant construction (left false)
plot_det_root_mean_cov = false;

if nargin < 2
    error('Usage: plotKappaRDDiagnostics(CalibrationDatabase, p [, alpha0Plot [, plotMode]])');
end
if nargin < 3 || isempty(alpha0Plot)
    alpha0Plot = min(CalibrationDatabase.alpha0);
end

% plotMode handling: accept only canonical values (no aliases)
if nargin < 4 || isempty(plotMode)
    plotMode = 'avg';
end
if ~(ischar(plotMode) || isstring(plotMode))
    error('plotKappaRDDiagnostics:badPlotMode', 'plotMode must be a character vector or string.');
end
plotMode = lower(char(plotMode));
validModes = {'pooled','avg','all'};
if ~ismember(plotMode, validModes)
    error('plotKappaRDDiagnostics:badPlotMode', ...
        'Unknown plotMode "%s". Use one of: ''pooled'', ''avg'', ''all''.', plotMode);
end

% booleans used later
doPooled = ismember(plotMode, {'pooled','all'});
doAvg    = ismember(plotMode, {'avg','all'});


% ---- plotting defaults ----
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
fs = 20;
set(groot,'defaultAxesFontSize',fs,'defaultTextFontSize',fs,'defaultLegendFontSize',fs);

colBase = [0.40 0.40 0.40];      % grey
colMRD  = [0.00 0.4470 0.7410];  % blue
colRDM  = [0.8500 0.3250 0.0980];% orange
colPooled = [0 0.6 0];           % green for pooled series
lw = 2; mkSize = 7;

% ---- select data ----
T = CalibrationDatabase(CalibrationDatabase.p==p & CalibrationDatabase.alpha0==alpha0Plot,:);
% validate required columns for plotting
requiredCols = {'kappa_RD','q_asym_sq','q_emp_avgquantile_sq'};
missing = setdiff(requiredCols, T.Properties.VariableNames);
if ~isempty(missing)
    error('plotKappaRDDiagnostics:missingCols','CalibrationDatabase is missing required columns: %s', strjoin(missing,','));
end
% pooled columns required when pooled plotting requested
if doPooled && ~ismember('q_emp_pooled_sq', T.Properties.VariableNames)
    error('plotKappaRDDiagnostics:missingPooled','plotMode includes ''pooled'' but CalibrationDatabase lacks ''q_emp_pooled_sq''.');
end
% ensure method column exists for method filtering
if ~ismember('method', T.Properties.VariableNames)
    error('plotKappaRDDiagnostics:missingMethod','CalibrationDatabase must have a ''method'' column (strings/cell chars).');
end

if isempty(T)
    error('No data found for p=%d, alpha0=%.4g in CalibrationDatabase.', p, alpha0Plot);
end
nuGrid = unique(T.nu);
nCols = numel(nuGrid);

% shared y-range across all panels, always including the reference value 1
allKappa = T.kappa_RD;
% include pooled kappa in range if available and requested
if doPooled && ismember('q_emp_pooled_sq', T.Properties.VariableNames) && ismember('q_asym_sq', T.Properties.VariableNames)
    allKappa = [allKappa; T.q_emp_pooled_sq ./ T.q_asym_sq];
end
allKappa = [allKappa; 1];
yPad = 0.05*(max(allKappa)-min(allKappa));
if yPad == 0, yPad = 0.05; end
yLim = [min(allKappa)-yPad, max(allKappa)+yPad];

% ---- figure ----
fh = figure('Name', sprintf('kappa_RD diagnostics: p=%d, alpha0=%.2f', p, alpha0Plot));
for iu = 1:nCols
    nu = nuGrid(iu);
    subplot(1,nCols,iu); hold on

    % plot avg-based kappa curves when requested
    if doAvg
        plotOneMethod(T, nu, 'asymptotic',        colBase, '-',  'o', 'Asymptotic',                          lw, mkSize);
        plotOneMethod(T, nu, 'det_mean_root',     colMRD,  '--', 's', 'Mean root determinant',                lw, mkSize);
        if plot_det_root_mean_cov
            plotOneMethod(T, nu, 'det_root_mean_cov', colRDM,  '-.', 'd', 'Root determinant of mean covariance',  lw, mkSize);
        end
    end

    % pooled-based kappa comparison 
    if doPooled && ismember('q_emp_pooled_sq', T.Properties.VariableNames) && ismember('q_asym_sq', T.Properties.VariableNames)
        idxP = T.nu==nu;
        if any(idxP)
            nnP = T.n(idxP);
            kappa_pooled = T.q_emp_pooled_sq(idxP) ./ T.q_asym_sq(idxP);
            [nnP,ordP] = sort(nnP);
            kappa_pooled = kappa_pooled(ordP);
            plot(nnP, kappa_pooled, 'LineStyle',':', 'Marker','^', 'Color',colPooled, ...
                'LineWidth',1.5, 'MarkerSize',mkSize, 'DisplayName','pooled kappa');
        end
    end

    yline(1,'k:','HandleVisibility','off', 'LineWidth', 2.5);
    
    if logscale
        set(gca,'XScale','log');
    end

    ylim(yLim);
    xlabel('$n$');
    title(sprintf('$p=%d$, $\\alpha_0=%.2f$, $\\nu=%d$',p, alpha0Plot,nu), 'Interpreter','latex');
    if iu==1
        ylabel('$\kappa_{RD}$', 'Interpreter','latex');
        legend('Location','best')
    end
end

% ---- summary diagnostics ----
fprintf('\nDiagnostic summary (p=%d, alpha0=%.4g):\n', p, alpha0Plot);
methods = {'asymptotic','det_mean_root','det_root_mean_cov'};
for m = 1:numel(methods)
    Tm = T(strcmp(T.method,methods{m}),:);
    if isempty(Tm), continue; end
    fprintf('  %-20s : min=%.4f  max=%.4f  max|kappa-1|=%.4f  (rows=%d)\n', ...
        methods{m}, min(Tm.kappa_RD), max(Tm.kappa_RD), max(abs(Tm.kappa_RD-1)), height(Tm));
end
end

function plotOneMethod(T,nu,method,col,style,marker,label,lw,mkSize)
idx = strcmp(T.method,method) & T.nu==nu;
if ~any(idx), return; end
[nn,ord] = sort(T.n(idx));
k = T.kappa_RD(idx); k = k(ord);
plot(nn, k, 'LineStyle',style, 'Marker',marker, 'Color',col, ...
    'LineWidth',lw, 'MarkerSize',mkSize, 'DisplayName',label);
end
