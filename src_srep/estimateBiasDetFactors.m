% =========================================================================
% Estimate determinant-based finite-sample correction factors
% =========================================================================
%
% Full call chain:
%   s01_mcdt_finite_sample_calibration.m
%       -> estimateBiasDetFactors.m
%
% Estimates the two determinant-based finite-sample correction factors used
% by the simulation:
%
%   1. mean of the replication-wise p-th roots of the determinant ratio
%   2. p-th root of the determinant ratio computed from the mean covariance
%
% The factors are estimated by Monte Carlo for each (alpha0, n, p, nu)
% combination, using the Student-t model and the chosen confidence level.
%
% The output is a struct with two tables:
%   - BiasDetTable.mean_root
%   - BiasDetTable.root_mean_cov
%
% If `saveFile` is provided, the struct is also saved to a MAT file.
% =========================================================================
    
function BiasDetTable = estimateBiasDetFactors(n_vec, p_vec, nu_vec, alpha0_vec, B_est, useParfor_est, saveFile, conflev)
% Runs MC to estimate BiasDet constructions per (alpha0,n,p,nu).
% Produces table with columns [alpha0 n p nu BiasDet_mean_root BiasDet_root_mean_cov]
rows = [];
for ip = 1:numel(p_vec)
    p = p_vec(ip);
    for iu = 1:numel(nu_vec)
        nu = nu_vec(iu);
        for ia = 1:numel(alpha0_vec)
            alpha0 = alpha0_vec(ia);
            for in = 1:numel(n_vec)
                n = n_vec(in);

                % initialize the Sigma^-1
                detMeanRootVec = nan(B_est,1);
                SigmaTrue    = (nu/(nu-2)) * eye(p);   % matches paper's simulation design (§4.1)
                SigmaTrueInv = inv(SigmaTrue);         % = (nu-2)/nu * eye(p), written generally for traceability
                SigmaAcc = zeros(p,p);

                if useParfor_est
                    % parfor needs sliced containers; collect local arrays
                    detTemp = nan(B_est,1);
                    SigmaCells = cell(B_est,1);
                    parfor b=1:B_est
                        X = randn(n,p).*sqrt(nu./chi2rnd(nu,n,1));
                        RAW = mcd(X,'modelT',nu,'bdp',alpha0,'conflev',conflev,'smallsamplecor',false,'nsamp',500,'refsteps',10,'plots',0,'msg',0,'nocheck',1);
                        % per-rep determinant of scaled covariance (assume underlying Sigma=I)
                        detSigma = det(RAW.cov * SigmaTrueInv);
                        detTemp(b) = detSigma^(1/p);
                        SigmaCells{b} = RAW.cov;
                    end
                    detMeanRootVec = detTemp;
                    for b=1:B_est
                        SigmaAcc = SigmaAcc + SigmaCells{b};
                    end
                else
                    for b=1:B_est
                        X = randn(n,p).*sqrt(nu./chi2rnd(nu,n,1));
                        RAW = mcd(X,'modelT',nu,'bdp',alpha0,'conflev',conflev,'smallsamplecor',false,'nsamp',500,'refsteps',10,'plots',0,'msg',0,'nocheck',1);
                        detSigma = det(RAW.cov * SigmaTrueInv);
                        detMeanRootVec(b) = detSigma^(1/p);
                        SigmaAcc = SigmaAcc + RAW.cov;
                    end
                end

                BiasDet_mean_root = mean(detMeanRootVec);
                SigmaMean = SigmaAcc / B_est;
                BiasDet_root_mean_cov = det(SigmaMean * SigmaTrueInv)^(1/p);

                rows = [rows; {alpha0, n, p, nu, BiasDet_mean_root, BiasDet_root_mean_cov}]; %#ok<AGROW>
                fprintf('Estimated BiasDet for n=%d p=%d nu=%d alpha0=%.3g: mean_root=%.6g root_mean_cov=%.6g\n', ...
                    n,p,nu,alpha0,BiasDet_mean_root,BiasDet_root_mean_cov);
            end
        end
    end
end

BiasDetTable = cell2table(rows, 'VariableNames', {'alpha0','n','p','nu','BiasDet_mean_root','BiasDet_root_mean_cov'});

% Convert to long table with single BiasDet column to match buildBiasDetFun expectations
alpha0_all = [];
n_all = [];
p_all = [];
nu_all = [];
biasdet_all = [];
for i=1:height(BiasDetTable)
    alpha0_all(end+1,1) = BiasDetTable.alpha0(i); %#ok<AGROW>
    n_all(end+1,1) = BiasDetTable.n(i); %#ok<AGROW>
    p_all(end+1,1) = BiasDetTable.p(i); %#ok<AGROW>
    nu_all(end+1,1) = BiasDetTable.nu(i); %#ok<AGROW>
    biasdet_all(end+1,1) = BiasDetTable.BiasDet_mean_root(i); %#ok<AGROW>
end
BiasDetTable_meanroot = table(alpha0_all,n_all,p_all,nu_all,biasdet_all, ...
    'VariableNames',{'alpha0','n','p','nu','BiasDet'});

% same for root_mean_cov, saved to a separate table variable inside a struct file
alpha0_all = [];
n_all = [];
p_all = [];
nu_all = [];
biasdet_all = [];
for i=1:height(BiasDetTable)
    alpha0_all(end+1,1) = BiasDetTable.alpha0(i); %#ok<AGROW>
    n_all(end+1,1) = BiasDetTable.n(i); %#ok<AGROW>
    p_all(end+1,1) = BiasDetTable.p(i); %#ok<AGROW>
    nu_all(end+1,1) = BiasDetTable.nu(i); %#ok<AGROW>
    biasdet_all(end+1,1) = BiasDetTable.BiasDet_root_mean_cov(i); %#ok<AGROW>
end
BiasDetTable_rootmeancov = table(alpha0_all,n_all,p_all,nu_all,biasdet_all, ...
    'VariableNames',{'alpha0','n','p','nu','BiasDet'});

% Pack both constructions into one MAT-friendly consolidated table-like struct
% To be compatible with loadBiasDetTable(...,estimator) we will save two separate tables
BiasDetTable = struct();
BiasDetTable.mean_root     = BiasDetTable_meanroot;
BiasDetTable.root_mean_cov = BiasDetTable_rootmeancov;
% Save the full struct to MAT
save(saveFile, 'BiasDetTable', '-v7.3');

% Save the determinant-bias estimation table
[saveFolder, saveBase] = fileparts(saveFile);

% Save the two tables to CSV
writetable(BiasDetTable_meanroot, ...
    fullfile(saveFolder, [saveBase '_mean_root.csv']));
writetable(BiasDetTable_rootmeancov, ...
    fullfile(saveFolder, [saveBase '_root_mean_cov.csv']));

end
