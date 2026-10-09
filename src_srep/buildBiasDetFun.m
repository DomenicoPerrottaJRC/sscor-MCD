%% ========================================================================
% Build BiasDet correction function handle
% ========================================================================
%
% Full call chain:
%   s01_mcdt_finite_sample_calibration.m
%       -> buildBiasDetFun.m
%
% Returns a function handle BiasDetFun(alpha0,n,p,nu) that provides the
% determinant-bias correction factor for the requested construction.
%
% The correction factor can be obtained in three ways:
%   - 'none': return factor 1
%   - 'table': load a correction table and use either exact lookup or
%      build an interpolant from the table
%   - 'interpolant': load a previously saved interpolant object directly
%
% The 'interpolant' option is useful when the interpolant has already been
% constructed offline and saved to a MAT file. It avoids rebuilding the
% interpolant from the raw table at every run.
%
% The returned handle is used later to scale the robust distances in the
% simulation.
% ========================================================================

function BiasDetFun = buildBiasDetFun(...
    source,...
    tableVar,...
    tableFile,...
    tableFolder,...
    tableMode,...
    estimator,...
    interpVar,...
    interpFile)

switch source
    case 'none'
        BiasDetFun = @(alpha0,n,p,nu) deal(1,true);

    case 'table'
        if isstruct(tableVar) && isfield(tableVar,'mean_root') && isfield(tableVar,'root_mean_cov')
            switch estimator
                case 'mean_root_det'
                    tableVar = tableVar.mean_root;
                case 'root_det_mean_cov'
                    tableVar = tableVar.root_mean_cov;
                otherwise
                    error('Unknown BiasDet estimator: %s', estimator);
            end
        end

        if isempty(tableVar) || (istable(tableVar) && height(tableVar)==0)
            if ~isempty(tableFile)
                [~,~,ext]=fileparts(tableFile);
                if strcmpi(ext,'.mat')
                    S=load(tableFile,'BiasDetTable');
                    tableVar=S.BiasDetTable;
                else
                    tableVar=readtable(tableFile);
                end
            elseif ~isempty(tableFolder)
                tableVar=loadBiasDetTable(tableFolder,estimator);
            else
                error('BIASDET_SOURCE=''table'' but no table available.');
            end
        end

        switch tableMode
            case 'exact'
                BiasDetFun = buildBiasDetFun_exact(tableVar);
            case 'interpolate'
                F = scatteredInterpolant(...
                    tableVar.alpha0,...
                    tableVar.n,...
                    tableVar.p,...
                    tableVar.nu,...
                    tableVar.BiasDet,...
                    'linear','nearest');
                BiasDetFun = @(alpha0,n,p,nu) deal(F(alpha0,n,p,nu), true);

            otherwise
                error('Unknown BiasDet table mode.');
        end

    case 'interpolant'
        if isempty(interpVar)
            if isempty(interpFile)
                error('Missing BiasDet interpolant.');
            end
            S = load(interpFile,'BiasDetInterpolant');
            interpVar = S.BiasDetInterpolant;
        end
        BiasDetFun = @(alpha0,n,p,nu) deal(interpVar(alpha0,n,p,nu), true);

    otherwise
        error('Unknown BIASDET_SOURCE.');
end
end


%% ========================================================================
% Load BiasDet simulations
% ========================================================================
%
% Full call chain:
%   s01_mcdt_finite_sample_calibration.m
%       -> buildBiasDetFun.m
%           -> loadBiasDetTable.m
%
% Loads the raw BiasDet simulation files from a folder and converts them
% into one long-format table with columns:
%   - alpha0
%   - n
%   - p
%   - nu
%   - BiasDet
%
% The `estimator` argument selects which stored bias-determinant field is
% extracted from the simulation output.
%% ========================================================================

function BiasDetTable = loadBiasDetTable(folder,estimator)
switch estimator
    case {'original','mean_root_det'}      % legacy name + current name
        field = 'bias_det';
    case {'new','root_det_mean_cov'}       % legacy name + current name
        field = 'bias_det_new';
    otherwise
        error('Unknown BiasDet estimator.');
end

files=dir(fullfile(folder,'F_outSIM_v*_nu*.mat'));
if isempty(files)
    error('No BiasDet simulation files found.');
end

alpha0_all=[];
n_all=[];
p_all=[];
nu_all=[];
bias_all=[];

for i=1:numel(files)
    filename=files(i).name;
    tok=regexp(filename,...
        'F_outSIM_v(\d+)_nu(\d+)',...
        'tokens',...
        'once');
    if isempty(tok)
        continue
    end
    p0=str2double(tok{1});
    nu0=str2double(tok{2});
    S=load(fullfile(files(i).folder,filename),'outSIM');
    if ~isfield(S,'outSIM')
        continue
    end
    outSIM=S.outSIM;
    for k=1:numel(outSIM)
        rec=outSIM(k);
        if ~isfield(rec,'alpha0') || ...
                ~isfield(rec,'n') || ...
                ~isfield(rec,field)
            continue
        end
        a=rec.alpha0(:);
        b=rec.(field)(:);
        if numel(a)~=numel(b)
            continue
        end
        nr=numel(a);
        alpha0_all=[alpha0_all;a];
        n_all=[n_all;repmat(rec.n,nr,1)];
        p_all=[p_all;repmat(p0,nr,1)];
        nu_all=[nu_all;repmat(nu0,nr,1)];
        bias_all=[bias_all;b];
    end
end

BiasDetTable=table(...
    alpha0_all,...
    n_all,...
    p_all,...
    nu_all,...
    bias_all,...
    'VariableNames',...
    {'alpha0','n','p','nu','BiasDet'});
end



%% ========================================================================
% Exact BiasDet lookup
% ========================================================================
%
% Full call chain:
%   s01_mcdt_finite_sample_calibration.m
%       -> buildBiasDetFun.m
%           -> buildBiasDetFun_exact
%               -> lookupBiasDetExact
%
% Returns a function handle that performs an exact table lookup for
% BiasDet(alpha0,n,p,nu). If the requested combination is missing, the
% factor defaults to 1 and `found` is false.
%% ========================================================================

function BiasDetFun=buildBiasDetFun_exact(tableVar)
tol=1e-6;
BiasDetFun=@(alpha0,n,p,nu) ...
    lookupBiasDetExact(...
    tableVar,...
    alpha0,...
    n,...
    p,...
    nu,...
    tol);
end

function [factor,found]=lookupBiasDetExact(...
    T,...
    alpha0,...
    n,...
    p,...
    nu,...
    tol)

idx = T.n==n & ...
    T.p==p & ...
    T.nu==nu & ...
    abs(T.alpha0 - alpha0) < tol;

if ~any(idx)
    warning(...
        'BiasDet value missing: using factor 1.');
    factor=1;
    found=false;
    return
end

id=find(idx,1);
factor = T.BiasDet(id);
found=true;
end
