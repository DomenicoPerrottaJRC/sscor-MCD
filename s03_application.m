%% Application to a real dataset
%
%% Application to Real Dataset (Grass) - MCD-t Finite-Sample Calibration
%
% This script demonstrates the finite-sample calibration of MCD-t on the
% Grass dataset (Barabesi et al., 2026):
%   1. Loads Grass_contaminated.prn (n=330, p=3) and extracts the
%      uncontaminated sample X0 (n0=300, p0=3).
%   2. Computes the Student-t MCD on X0 (nu0=8, alpha0=0.25, conflev0=0.95).
%   3. Obtains calibration factors:
%      - a1/a2: Determinant scatter calibration factor delta_n (MC & LOESS)
%      - b1/b2: Robust-distance calibration factor kappa_RD^delta (MC & LOESS)
%      - c:     Theoretical asymptotic cutoff q^2_{nu,p}(c) via msdcutoff
%   4. Compares outlier classifications: shows that borderline observations
%      exceeding the asymptotic cutoff are correctly retained by the
%      calibrated threshold.
%   5. Produces diagnostic plots (Distance Index Plot & Scatterplot Matrix).


clear; clc; close all;

datatype = 1; % 1 = Grass ; 2 = t ; 3 = S&P

%% 1. Configuration and Paths
ds_path      = './datasets';
dsfile1      = 'Grass.prn';
dsfile2      = 'Grass_original_contaminated.prn';
interpFolder = './NEW/MCDt_all_interpolation_analysis_p3_nu8';

% Model parameters (Barabesi et al., 2026)
nu0      = 8;
alpha0   = 0.50;
conflev0 = 0.95;

%% 2a. Load Grass Dataset

if datatype == 1

    dsfile = fullfile(ds_path, dsfile1);
    if ~isfile(dsfile)
        error('Dataset file not found: %s', dsfile);
    end
    
    X0 = load(dsfile);
    n0 = size(X0, 1);
    p0 = size(X0, 2);

    fprintf('=== Grass Dataset Overview ===\n');
    fprintf('Uncontaminated sample size: n0 = %d, Dimension: p0 = %d\n', n0, p0);
    fprintf('Target model: Student-t with nu0 = %d, alpha0 = %.2f, conflev = %.2f\n\n', ...
        nu0, alpha0, conflev0);

    %{
        % just a check to see if the dataset is correct
        % Label the known 30 contaminated observations
        Xcont    = load(fullfile(ds_path, dsfile2));
        n1       = size(Xcont, 1);
        Xlabels1 = zeros(n1, 1);
        Xlabels1(n1 - 30 + 1 : end) = 1;
        Xcheck = Xcont(Xlabels1 == 0, :); % Extract uncontaminated data (n0 = 300)
        if norm(X0)-norm(Xcheck) > 10^(-6) , disp('check dataset'); end;
    %}

end

%% 2b. Generate independent standard Student-t observations

if datatype == 2
    p0  = 3;
    n0  = 100;
    rng(12345);
    X0 = trnd(nu0, n0, p0);
end

%% 2c. Financial from S&P
if datatype == 3
    dataset_SP500;
end

%% 3. Compute MCD on Uncontaminated Data
RAW0 = mcd(X0, ...
    'modelT', nu0, ...
    'bdp', alpha0, ...
    'conflev', conflev0, ...
    'smallsamplecor', false, ...
    'nsamp', 500, ...
    'refsteps', 10, ...
    'plots', 0, ...
    'msg', 0, ...
    'nocheck', 1);

% Extract raw squared robust distances
RD2_raw = RAW0.md;

%% 4. Obtain Calibration Factors and Cutoffs

% --- (c) Theoretical Asymptotic Cutoff (squared scale) ---
q2_asym = msdcutoff(conflev0, p0, nu0);

% --- Load calibration tables and LOESS models ---
deltaMatFile = fullfile(interpFolder, 'delta_interpolation', 'MCDt_delta_interpolation_analysis.mat');
rdMatFile    = fullfile(interpFolder, 'RD_interpolation',    'MCDt_RD_interpolation_analysis.mat');

% Default fallbacks if exact grid lookup is needed
delta_mc     = NaN;
delta_loess  = NaN;
kappa_mc     = NaN;
kappa_loess  = NaN;

% Retrieve delta (Section 5.4)
if isfile(deltaMatFile)
    deltaData = load(deltaMatFile);
    % (a1) Exact Monte Carlo grid lookup (if available)
    matchRow = deltaData.D.p == p0 & deltaData.D.nu == nu0 & ...
        deltaData.D.alpha0 == alpha0 & deltaData.D.n == n0;
    if any(matchRow)
        delta_mc = deltaData.D.delta_hat(find(matchRow, 1));
    end

    % (a2) LOESS evaluation
    % Find model corresponding to (p0, nu0)
    for m = 1:numel(deltaData.DeltaSummary.loessModels)
        mdl = deltaData.DeltaSummary.loessModels(m);
        if mdl.p == p0 && mdl.nu == nu0 && ~isempty(mdl.fit)
            delta_loess = mdl.fit(alpha0, n0);
            break;
        end
    end
end

% Retrieve kappa_RD^delta (Section 6)
if isfile(rdMatFile)
    rdData = load(rdMatFile);
    % (b1) Exact Monte Carlo grid lookup (if available)
    matchRow = rdData.R.p == p0 & rdData.R.nu == nu0 & ...
        rdData.R.conflev == conflev0 & ...
        rdData.R.alpha0 == alpha0 & rdData.R.n == n0;
    if any(matchRow)
        kappa_mc = rdData.R.kappa_RD(find(matchRow, 1));
    end

    % (b2) LOESS evaluation
    for m = 1:numel(rdData.RDSummary.loessModels)
        mdl = rdData.RDSummary.loessModels(m);
        if mdl.p == p0 && mdl.nu == nu0 && mdl.conflev == conflev0 && ~isempty(mdl.fit)
            kappa_loess = mdl.fit(alpha0, n0);
            break;
        end
    end
end

% Stop if the requested configuration is not available in the precomputed
% interpolation models.
if isnan(delta_loess)
    error(['Could not estimate delta_loess for (p=%g, nu=%g, alpha0=%.3f, n=%d). ' ...
        'This configuration is probably not available in the precomputed ' ...
        'delta interpolation table. Check p0, nu0, alpha0, and n0.'], ...
        p0, nu0, alpha0, n0);
end

if isnan(kappa_loess)
    error(['Could not estimate kappa_loess for (p=%g, nu=%g, conflev=%.3f, alpha0=%.3f, n=%d). ' ...
        'This configuration is probably not available in the precomputed ' ...
        'RD interpolation table. Check p0, nu0, conflev0, alpha0, and n0.'], ...
        p0, nu0, conflev0, alpha0, n0);
end

fprintf('=== Calibration Factors ===\n');
fprintf('(c)  Asymptotic squared cutoff:  q2_asym     = %.4f\n', q2_asym);
fprintf('(a1) Monte Carlo delta_n:         delta_mc    = %.4f\n', delta_mc);
fprintf('(a2) LOESS-interpolated delta_n:  delta_loess = %.4f\n', delta_loess);
fprintf('(b1) Monte Carlo kappa_RD^delta:  kappa_mc    = %.4f\n', kappa_mc);
fprintf('(b2) LOESS kappa_RD^delta:        kappa_loess = %.4f\n\n', kappa_loess);

%% 5. Objective A & B: Compare Corrected Distances and Calibrated Cutoffs

% Corrected squared distances: RD2_delta = RD2_raw / delta
RD2_delta_mc    = RD2_raw ./ delta_mc;
RD2_delta_loess = RD2_raw ./ delta_loess;

% Calibrated squared cutoffs: q2_cal = kappa_RD^delta * q2_asym (equation (24))
q2_cal_mc    = kappa_mc    * q2_asym;
q2_cal_loess = kappa_loess * q2_asym;

% Outlier identification
% 1. Asymptotic rule (Uncalibrated): RD2_raw > q2_asym
outliers_asym = find(RD2_raw >  q2_asym);
good_asym     = find(RD2_raw <= q2_asym);

% 2. Fully calibrated rule (using LOESS): RD2_delta_loess > q2_cal_loess
%    (Equivalently: RD2_raw > delta_loess * kappa_loess * q2_asym)
outliers_cal_loess = find(RD2_delta_loess > q2_cal_loess);
outliers_cal_mc    = find(RD2_delta_mc    > q2_cal_mc);
cal_outliers_loess = outliers_cal_loess;

% Borderline points (false positives under asymptotic cutoff)
borderline_idx_loess = setdiff(outliers_asym, outliers_cal_loess);
borderline_idx_mc    = setdiff(outliers_asym, outliers_cal_mc);

fprintf('=== Outlier Detection Comparison ===\n');
fprintf('Nominal size: %.2f%% (Expected false positives ~ %d)\n', (1-conflev0)*100, round((1-conflev0)*n0));
fprintf('Outliers flagged by Asymptotic cutoff: %d\n', numel(outliers_asym));
fprintf('Outliers flagged by Calibrated cutoff: %d\n', numel(outliers_cal_loess));
fprintf('Borderline observations reclassified as regular: %d (indices: %s)\n\n', ...
    numel(borderline_idx_loess), mat2str(borderline_idx_loess'));

%% groups and colors

% groups labels, to be used in the plots
group_labels = nan(n0, 1);
group_labels(good_asym) = 0;            % asym regular, unfilled
group_labels(cal_outliers_loess) = 1;   % calibrated final outliers
group_labels(borderline_idx_loess) = 2; % asymptotic-only outliers

regularColor = [0.2 0.2 0.2];
regularFace  = 'none';

calColor = [0.00 0.45 0.74];
calFace  = [0.00 0.45 0.74];

borderColor = [0.85 0.33 0.10];
borderFace  = [0.85 0.33 0.10];

%% 6. Diagnostic Plot 1: Robust Distance Index Plot
figure('Name', 'Grass_RD_Diagnostic', 'Position', [100, 100, 900, 500]);
hold on; grid on;

idx0 = find(group_labels == 0); % asym regular, unfilled
idx1 = find(group_labels == 1); % calibrated final outliers
idx2 = find(group_labels == 2); % asymptotic-only outliers

plot(idx0, RD2_delta_loess(idx0), 'o', ...
    'Color', regularColor, 'MarkerFaceColor', regularFace, ...
    'MarkerSize', 5, 'DisplayName', 'Regular units');

plot(idx1, RD2_delta_loess(idx1), 'o', ...
    'Color', calColor, 'MarkerFaceColor', calFace, ...
    'MarkerSize', 6, 'DisplayName', 'Calibrated outliers');

plot(idx2, RD2_delta_loess(idx2), 'o', ...
    'Color', borderColor, 'MarkerFaceColor', borderFace, ...
    'MarkerSize', 8, 'DisplayName', 'Borderline (asymptotic only)');

yline(q2_asym, '--r', 'LineWidth', 1.8, ...
    'DisplayName', 'Asymptotic Cutoff $q^2_{\nu,p}(c)$');
yline(q2_cal_loess, '-b', 'LineWidth', 2.0, ...
    'DisplayName', 'Calibrated Cutoff $q^{2,\delta}_{\mathrm{cal}}(c)$');

xlabel('Observation index $i$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Squared Robust Distance $RD_i^{2,\delta}$', 'Interpreter', 'latex', 'FontSize', 14);
title(sprintf('Grass Dataset ($n=%d, p=%d$): Asymptotic vs Calibrated Cutoffs', n0, p0), ...
    'Interpreter', 'latex', 'FontSize', 16);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 12);
ylim([0, max([RD2_delta_loess; q2_cal_loess]) * 1.15]);


%% 7. Diagnostic Plot 2: Scatterplot Matrix (spmplot)
% Same group convention as section 6:
%   0 = regular units, open gray
%   1 = calibrated outliers, filled blue
%   2 = borderline units, filled orange

groupOrder = [0 1 2];
[~, groupIdx] = ismember(group_labels, groupOrder);

groupNames = strings(n0, 1);
groupNames(group_labels == 0) = "Regular units";
groupNames(group_labels == 1) = "Calibrated outliers";
groupNames(group_labels == 2) = "Borderline (asymptotic only)";

plo = struct;
plo.sym = repmat({'o'}, 1, numel(groupOrder));
plo.doleg = 'on';

% [H1, AX1, BigAx1] = spmplot(X0, 'group', groupIdx, ...
%     'dispopt', 'box', 'tag', 'Grass', 'plo', plo);
[H1, AX1, BigAx1] = spmplot(X0, 'group', groupNames, ...
    'plo', plo, 'dispopt', 'box');


sgtitle('Grass Data: Outlier Calibration Diagnostics', ...
    'FontSize', 16, 'Interpreter', 'latex');

% Recolor the scatterplot matrix to match section 6 exactly.
pData = size(X0, 2);

if ismatrix(AX1) && size(AX1,1) >= pData && size(AX1,2) >= pData
    axGrid = AX1(1:pData, 1:pData);
    mask = ~eye(pData);
    axOff = axGrid(mask);
    [ii, jj] = ind2sub([pData, pData], find(mask));
else
    axOff = AX1(:);
    axOff = axOff(1:min(numel(axOff), pData * (pData - 1)));
    [ii, jj] = ind2sub([pData, pData], (1:numel(axOff)).');
    keep = ii ~= jj;
    axOff = axOff(keep);
    ii = ii(keep);
    jj = jj(keep);
end


for k = 1:numel(axOff)
    ax = axOff(k);
    i = ii(k);
    j = jj(k);

    hLines = findall(ax, 'Type', 'line');
    hLines = hLines(arrayfun(@(h) isprop(h, 'Marker') && ~strcmp(h.Marker, 'none'), hLines));

    if isempty(hLines)
        continue;
    end

    for h = 1:numel(hLines)
        coords = [hLines(h).XData(:), hLines(h).YData(:)];

        for g = 0:2
            expected = X0(group_labels == g, [j i]);

            if isempty(expected)
                continue;
            end

            if size(coords, 1) == size(expected, 1) && isequal(sortrows(coords), sortrows(expected))
                switch g
                    case 0
                        hLines(h).MarkerEdgeColor = regularColor;
                        hLines(h).MarkerFaceColor = regularFace;
                    case 1
                        hLines(h).MarkerEdgeColor = calColor;
                        hLines(h).MarkerFaceColor = calFace;
                    case 2
                        hLines(h).MarkerEdgeColor = borderColor;
                        hLines(h).MarkerFaceColor = borderFace;
                end
                break;
            end
        end
    end
end

%set(AX1, 'FontSize', 12);
% saveas(BigAx1,[ds_path filesep 'grass.fig'],'fig');
% saveas(BigAx1,[ds_path filesep 'grass.eps'],'epsc');




%% Notes

% This script demonstrates the effect of the correction on a real dataset
% discussed in Barabesi et al (2026). The script does the following:
% 1. Load the Grass.prn dataset, say X0, of size n0=300, p0=3.
% 2. Compute the mcd on Grass, using these options:
%    RAW0 = mcd(X0,...
%     'modelT',nu0,...
%     'bdp',alpha0,...
%     'conflev',conflev0,...
%     'smallsamplecor',false,...
%     'nsamp',500,...
%     'refsteps',10,...
%     'plots',0,...
%     'msg',0,...
%     'nocheck',1);
%    The parameters are chosen as in Barabesi et al (2026):
%    - nu0 = 8
%    - alpha0 = 0.25
%    - conflev0 = 0.95
% 3. Extract the squared robust distances: RD2 = RAW.md
% 4. Get:
%    a1. the montecarlo   Pison-like coefficient \delta_n
%    a2. the interpolated Pison-like coefficient \hat{\delta}_n
%    b1. the distance-based calibration factor \kappa_{\mathrm{RD}}^{\delta}
%    b2. the interpolated distance-based calibration factor \hat{\kappa}_{\mathrm{RD}}^{\delta}
%    c.  the asymptotic cutoff \kappa \equiv \eta_{\alpha_0,p}(\nu),
%        which is computed in matlab using
%        cutoffT0 = msdcutoff(conflev0,p0,nu0);
%
% OBJECTIVES:
%  A. Compare the distances corrected with the inverse of (a1) with the
%     quantile corrected with (b1).
%  B. Same as A, but for the interpolated quantities (a2) and (b2).
%
% EXPECTATIONS:
% show that some distances above the asymptotic quantile (c) are below the
% corrected quantile (b)
%
% Files location (under the NEW folder):
% - Results_correction_alpha0extended_2026-09-27_120522
% - MCDt_RD_analysis
% - MCDt_all_interpolation_analysis
