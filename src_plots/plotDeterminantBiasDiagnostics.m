function hf = plotDeterminantBiasDiagnostics(T_mean_root, T_root_mean_cov, p, alpha0Plot)

% plotDeterminantBiasDiagnostics
%
% Inputs:
%   T_mean_root, T_root_mean_cov : long-format tables with columns
%       [alpha0 n p nu BiasDet]
%   p          : dimension to plot (scalar)
%   alpha0Plot : which alpha0 to plot (scalar; defaults to the smallest
%                value present in T_mean_root if omitted)

logscale = false;

plot_det_root_mean_cov = false; % leave second construction off by default

if nargin < 3
    error('Usage: plotDeterminantBiasDiagnostics(T_mean_root, T_root_mean_cov, p [, alpha0Plot])');
end

if nargin < 4 || isempty(alpha0Plot)
    alpha0Plot = min(T_mean_root.alpha0);
end

% ---- plotting defaults ----
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
fs = 20;
set(groot,'defaultAxesFontSize',fs,'defaultTextFontSize',fs,'defaultLegendFontSize',fs);

colMRD = [0.00 0.4470 0.7410];   % blue
colRDM = [0.8500 0.3250 0.0980]; % orange
lw = 2; mkSize = 7;

% ---- select data ----
idxMRall = T_mean_root.p==p & T_mean_root.alpha0==alpha0Plot;
idxRMall = T_root_mean_cov.p==p & T_root_mean_cov.alpha0==alpha0Plot;
nuGrid = unique(T_mean_root.nu(idxMRall));
nCols = numel(nuGrid);
if nCols == 0
    error('No data found for p=%d, alpha0=%.4g in T_mean_root.', p, alpha0Plot);
end

% shared y-range per row, across all nu columns (not across rows)
dAll     = [T_mean_root.BiasDet(idxMRall); T_root_mean_cov.BiasDet(idxRMall); 1];
padRow1  = 0.05*(max(dAll)-min(dAll)); if padRow1==0, padRow1=0.05; end
yLimRow1 = [min(dAll)-padRow1, max(dAll)+padRow1];

deltaAll = [1./T_mean_root.BiasDet(idxMRall); 1./T_root_mean_cov.BiasDet(idxRMall); 1];
padRow2  = 0.05*(max(deltaAll)-min(deltaAll)); if padRow2==0, padRow2=0.05; end
yLimRow2 = [min(deltaAll)-padRow2, max(deltaAll)+padRow2];

% ---- figure ----
hf = figure('Name', sprintf('Determinant bias diagnostics: p=%d, alpha0=%.2f', p, alpha0Plot));
for iu = 1:nCols
    nu = nuGrid(iu);

    idxMR = T_mean_root.p==p & T_mean_root.nu==nu & T_mean_root.alpha0==alpha0Plot;
    [nMR,ordMR] = sort(T_mean_root.n(idxMR));
    dMR = T_mean_root.BiasDet(idxMR); dMR = dMR(ordMR);

    if plot_det_root_mean_cov
        idxRM = T_root_mean_cov.p==p & T_root_mean_cov.nu==nu & T_root_mean_cov.alpha0==alpha0Plot;
        [nRM,ordRM] = sort(T_root_mean_cov.n(idxRM));
        dRM = T_root_mean_cov.BiasDet(idxRM); dRM = dRM(ordRM);
    end

    % ---- row 1: d-bar (Eq. 9) ----
    subplot(2,nCols,iu); hold on
    plot(nMR, dMR, 'LineStyle','--', 'Marker','s', 'Color',colMRD, ...
        'LineWidth',lw, 'MarkerSize',mkSize, 'DisplayName','mean root determinant');
    if plot_det_root_mean_cov
        plot(nRM, dRM, 'LineStyle','-.', 'Marker','d', 'Color',colRDM, ...
            'LineWidth',lw, 'MarkerSize',mkSize, 'DisplayName','root determinant of mean covariance');
    end
    yline(1,'k:','HandleVisibility','off');

    if logscale
        set(gca,'XScale','log');
    end

    ylim(yLimRow1);
    title(sprintf('$ p=%d $, $ \\alpha_0=%.2f $, $\\nu=%d$', p, alpha0Plot, nu), 'Interpreter', 'latex');
    if iu==1
        ylabel('$\bar d^{\alpha_0}_{p,n,\nu}$', 'Interpreter','latex');
        %legend('Location','best')
    end

    % ---- row 2: delta-hat (Eq. 10) ----
    subplot(2,nCols,nCols+iu); hold on
    plot(nMR, 1./dMR, 'LineStyle','--', 'Marker','s', 'Color',colMRD, 'LineWidth',lw, 'MarkerSize',mkSize);
    if plot_det_root_mean_cov
        plot(nRM, 1./dRM, 'LineStyle','-.', 'Marker','d', 'Color',colRDM, 'LineWidth',lw, 'MarkerSize',mkSize);
    end
    yline(1,'k:','HandleVisibility','off', 'LineWidth', 2.5);

    if logscale
        set(gca,'XScale','log');
    end

    ylim(yLimRow2);
    xlabel('$n$', 'Interpreter','latex');
    if iu==1
        ylabel('$\widehat{\delta}_n$', 'Interpreter','latex');
    end
end

end
