%% ========================================================================
% Plot LOESS interpolants for RD and delta
% ========================================================================
%
% Full call chain:
%   s02_mcdt_interpolation_analysis.m
%       -> plot_loess_interpolants.m
%
% Loads the saved dense LOESS lookup tables produced by the interpolation
% analysis and renders one surface for:
%   - Section 6: kappa_RD^delta
%   - Section 5.4: delta_hat
%
% Inputs:
%   - pSel, nuSel, conflevSel: surface slice to display
%   - InterpolationRoot: root folder containing the interpolation outputs
%
% Outputs:
%   - hf_RD: figure handle for the RD surface
%   - hf_delta: figure handle for the delta surface
%
% Example:
%{
    % The output root of the interpolation analysis.
    InterpolationRoot = './NEW/MCDt_all_interpolation_analysis';
    
    % Pick one surface to plot
    pSel = 2;
    nuSel = 5;
    conflevSel = 0.95; % for RD only

   [hf_RD , hf_delta] = plot_loess_interpolants(pSel,nuSel,conflevSel,InterpolationRoot)

%}
%% ========================================================================

function [hf_RD , hf_delta] = plotLoessInterpolants(pSel,nuSel,conflevSel,InterpolationRoot)

% Folder layout created by s02.
rdOut    = [InterpolationRoot filesep 'RD_interpolation'];
deltaOut = [InterpolationRoot filesep 'delta_interpolation'];

% Dense lookup tables saved by the interpolation analysis.
T_RD    = readtable(fullfile(rdOut,    'RD_loess_lookup.csv'));
T_delta = readtable(fullfile(deltaOut, 'delta_loess_lookup.csv'));

% RD surface: select one (p, nu, conflev) slice.
T   = T_RD;
sel = T.p == pSel & T.nu == nuSel & T.conflev == conflevSel;
S   = T(sel, :);

nVals = unique(S.n, 'sorted');
aVals = unique(S.alpha0, 'sorted');

Z = nan(numel(aVals), numel(nVals));
[~, in] = ismember(S.n, nVals);
[~, ia] = ismember(S.alpha0, aVals);
Z(sub2ind(size(Z), ia, in)) = S.loess;

hf_RD = figure;
setfont = 20;
surf(nVals, aVals, Z, 'EdgeColor', 'none')
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 0.5;
xlabel('$n$','Interpreter','Latex','FontSize',setfont)
ylabel('$\alpha_0$','Interpreter','Latex','FontSize',setfont)
zlabel('$\kappa_{\mathrm{RD}}^{\delta}(n,p,\nu,\alpha_0,c)$ LOESS surface','Interpreter','Latex','FontSize',setfont)
%title(sprintf('RD LOESS surface: p=%g, nu=%g, conflev=%g', pSel, nuSel, conflevSel))
%colorbar
view(45, 30)
grid on

saveas(gcf, fullfile(rdOut, sprintf('RD_loess_p%g_nu%g_conf%g.png', ...
    pSel, nuSel, conflevSel)));

% Delta surface: select one (p, nu) slice.
T   = T_delta;
sel = T.p == pSel & T.nu == nuSel ;
S   = T(sel, :);

nVals = unique(S.n, 'sorted');
aVals = unique(S.alpha0, 'sorted');

Z = nan(numel(aVals), numel(nVals));
[~, in] = ismember(S.n, nVals);
[~, ia] = ismember(S.alpha0, aVals);
Z(sub2ind(size(Z), ia, in)) = S.loess;

hf_delta = figure;
setfont = 20;
surf(nVals, aVals, Z, 'EdgeColor', 'none')
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 0.5;
xlabel('$n$','Interpreter','Latex','FontSize',setfont)
ylabel('$\alpha_0$','Interpreter','Latex','FontSize',setfont)
zlabel('$\widehat{\delta}(n,p,\alpha_0,\nu)$ LOESS surface','Interpreter','Latex','FontSize',setfont)
%title(sprintf('RD LOESS surface: p=%g, nu=%g, conflev=%g', pSel, nuSel, conflevSel))
%colorbar
view(45, 30)
grid on

saveas(gcf, fullfile(deltaOut, sprintf('delta_loess_p%g_nu%g.png', ...
    pSel, nuSel)))

end
