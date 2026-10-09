%% ========================================================================
% Paired simulation for one configuration
% ========================================================================
%
% Full call chain:
%   s01_mcdt_finite_sample_calibration.m
%       -> runPairedFullGrid.m
%           -> runPairedSimConfig.m
%
% Split-sample calibration:
%   - replications 1:B/2 are used as the calibration sample
%   - replications B/2+1:B are used as the evaluation sample
%
% Three calibration definitions are computed:
%   - pooled: one quantile from all calibration observations
%   - avg: Eq. (14), average of replication-specific calibration quantiles
%   - replication-wise: each evaluation replication uses its matching
%     calibration replication quantile
%
% This avoids using the same observations to both determine the cutoff and
% evaluate the empirical size.
%% ========================================================================

function [calibBase,...
    calibMeanRoot,...
    calibRootMeanCov,...
    evalBase,...
    evalMeanRoot,...
    evalRootMeanCov] = ...
    runPairedSimConfig(...
    n,...
    p,...
    nu,...
    alpha0,...
    B,...
    cutoff,...
    factorMeanRoot,...
    factorRootMeanCov,...
    useParfor,...
    verifyCorrection,...
    conflev,...
    q_method)


% ------------------------------------------------------------------------
% Split B into calibration and evaluation halves
% -------------------------------------------------------------------------
if mod(B,2) ~= 0
    error('B must be even for split-sample calibration.');
end

B_cal = B/2;
B_eval = B/2;

% ------------------------------------------------------------------------
% Storage
% -------------------------------------------------------------------------
mdBaseMat      = nan(n,B);
mdMeanRootMat  = nan(n,B);
mdRootMeanCovMat = nan(n,B);

countsBase     = nan(B,1);
countsMeanRoot = nan(B,1);
countsRootMeanCov = nan(B,1);

% ------------------------------------------------------------------------
% Run the paired Monte Carlo simulation
% -------------------------------------------------------------------------
if useParfor

    parfor b = 1:B

        [mdBaseMat(:,b),...
            countsBase(b),...
            mdMeanRootMat(:,b),...
            countsMeanRoot(b),...
            mdRootMeanCovMat(:,b),...
            countsRootMeanCov(b)] = ...
            onePairedReplication(...
            n,...
            p,...
            nu,...
            alpha0,...
            cutoff,...
            factorMeanRoot,...
            factorRootMeanCov,...
            verifyCorrection,...
            conflev);

    end

else

    for b = 1:B

        [mdBaseMat(:,b),...
            countsBase(b),...
            mdMeanRootMat(:,b),...
            countsMeanRoot(b),...
            mdRootMeanCovMat(:,b),...
            countsRootMeanCov(b)] = ...
            onePairedReplication(...
            n,...
            p,...
            nu,...
            alpha0,...
            cutoff,...
            factorMeanRoot,...
            factorRootMeanCov,...
            verifyCorrection,...
            conflev);

    end

end

% ------------------------------------------------------------------------
% Separate calibration and evaluation replications
% -------------------------------------------------------------------------

calibIdx = 1:B_cal;
evalIdx  = B_cal+1:B;

mdBaseCal = mdBaseMat(:,calibIdx);
mdBaseEval = mdBaseMat(:,evalIdx);

mdMeanRootCal = mdMeanRootMat(:,calibIdx);
mdMeanRootEval = mdMeanRootMat(:,evalIdx);

mdRootMeanCovCal = mdRootMeanCovMat(:,calibIdx);
mdRootMeanCovEval = mdRootMeanCovMat(:,evalIdx);

% ------------------------------------------------------------------------
% Replication-specific calibration quantiles
%
% These are the q_emp^{2,(b)} quantities entering Eq. (14).
% -------------------------------------------------------------------------

qBaseRepCal_sq = ...
    quantile(mdBaseCal,conflev,1, Method=q_method)';%, Method='inclusive'

qMeanRootRepCal_sq = ...
    quantile(mdMeanRootCal,conflev,1, Method=q_method)';%, Method='inclusive'

qRootMeanCovRepCal_sq = ...
    quantile(mdRootMeanCovCal,conflev,1, Method=q_method)';%, Method='inclusive'

% ------------------------------------------------------------------------
% OPTION B: Eq. (14)-(15), average of replication-specific quantiles
%
% q_emp^2 = (1/B_cal) sum_b q_emp^{2,(b)}
%
% IMPORTANT:
% the average is now formed using ONLY the calibration replications.
% -------------------------------------------------------------------------

qBaseAvgCal_sq = mean(qBaseRepCal_sq);
qMeanRootAvgCal_sq = mean(qMeanRootRepCal_sq);
qRootMeanCovAvgCal_sq = mean(qRootMeanCovRepCal_sq);

% ------------------------------------------------------------------------
% OPTION A: pooled calibration
%
% One quantile from all n*B_cal calibration observations.
% -------------------------------------------------------------------------

qBasePooledCal_sq = ...
    quantile(mdBaseCal(:),conflev, Method=q_method);%, Method='inclusive'

qMeanRootPooledCal_sq = ...
    quantile(mdMeanRootCal(:),conflev, Method=q_method);%, Method='inclusive'

qRootMeanCovPooledCal_sq = ...
    quantile(mdRootMeanCovCal(:),conflev, Method=q_method); %, Method='inclusive'

% ------------------------------------------------------------------------
% Store calibration information
% -------------------------------------------------------------------------

calibBase.q_rep_sq = qBaseRepCal_sq;
calibBase.q_avg_sq = qBaseAvgCal_sq;
calibBase.q_pooled_sq = qBasePooledCal_sq;

calibMeanRoot.q_rep_sq = qMeanRootRepCal_sq;
calibMeanRoot.q_avg_sq = qMeanRootAvgCal_sq;
calibMeanRoot.q_pooled_sq = qMeanRootPooledCal_sq;

calibRootMeanCov.q_rep_sq = qRootMeanCovRepCal_sq;
calibRootMeanCov.q_avg_sq = qRootMeanCovAvgCal_sq;
calibRootMeanCov.q_pooled_sq = qRootMeanCovPooledCal_sq;

% ------------------------------------------------------------------------
% OPTION C: replication-wise calibration
%
% Pair calibration replication b with evaluation replication b.
%
% For b = 1,...,B_eval:
%
%       cutoff_b = q_emp,cal^{2,(b)}
%
% and evaluate:
%
%       P(MD_eval,b^2 > cutoff_b)
%
% This is NOT circular because the observations defining cutoff_b and
% the observations testing cutoff_b come from disjoint replications.
% -------------------------------------------------------------------------

repSizeBase = nan(B_eval,1);
repSizeMeanRoot = nan(B_eval,1);
repSizeRootMeanCov = nan(B_eval,1);

for b = 1:B_eval

    repSizeBase(b) = ...
        mean(mdBaseEval(:,b) > qBaseRepCal_sq(b));

    repSizeMeanRoot(b) = ...
        mean(mdMeanRootEval(:,b) > qMeanRootRepCal_sq(b));

    repSizeRootMeanCov(b) = ...
        mean(mdRootMeanCovEval(:,b) > qRootMeanCovRepCal_sq(b));

end

evalBase.size_replicationwise = mean(repSizeBase);
evalMeanRoot.size_replicationwise = mean(repSizeMeanRoot);
evalRootMeanCov.size_replicationwise = mean(repSizeRootMeanCov);

evalBase.size_replicationwise_vec = repSizeBase;
evalMeanRoot.size_replicationwise_vec = repSizeMeanRoot;
evalRootMeanCov.size_replicationwise_vec = repSizeRootMeanCov;

% ------------------------------------------------------------------------
% OPTION A: pooled calibration evaluated on independent evaluation sample
% -------------------------------------------------------------------------

evalBase.size_pooled = ...
    mean(mdBaseEval(:) > qBasePooledCal_sq);

evalMeanRoot.size_pooled = ...
    mean(mdMeanRootEval(:) > qMeanRootPooledCal_sq);

evalRootMeanCov.size_pooled = ...
    mean(mdRootMeanCovEval(:) > qRootMeanCovPooledCal_sq);

% ------------------------------------------------------------------------
% OPTION B: Eq. (14) average calibration evaluated on independent
% evaluation sample
% -------------------------------------------------------------------------

evalBase.size_avg = ...
    mean(mdBaseEval(:) > qBaseAvgCal_sq);

evalMeanRoot.size_avg = ...
    mean(mdMeanRootEval(:) > qMeanRootAvgCal_sq);

evalRootMeanCov.size_avg = ...
    mean(mdRootMeanCovEval(:) > qRootMeanCovAvgCal_sq);

% ------------------------------------------------------------------------
% Baseline size on independent evaluation sample
% -------------------------------------------------------------------------

evalBase.size_asym = ...
    mean(mdBaseEval(:) > cutoff);

evalMeanRoot.size_asym = ...
    mean(mdMeanRootEval(:) > cutoff);

evalRootMeanCov.size_asym = ...
    mean(mdRootMeanCovEval(:) > cutoff);

end
