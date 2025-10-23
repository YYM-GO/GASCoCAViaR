function [H,univariate,parameters,parametersallPOT] = dcc_fit_variance_POT(data2d,p,o,q,gjrType,startingVals)
% Fits TARCH models for use in DCC and related estimators
%
% USAGE:
%   [H,UNIVARIATE] = dcc_fit_variance(DATA,P,O,Q,GJRTYPE)
%
% INPUTS:
%   DATA    - A column of mean zero data
%   P       - K by 1 vector of positive, scalar integers representing the number of symmetric innovations
%   O       - K by 1 vector of non-negative scalar integers representing the number of asymmetric innovations (0
%                    for symmetric processes)
%   Q       - K by 1 vector of non-negative, scalar integers representing the number of lags of conditional
%                    variance (0 for ARCH)
%   GJRTYPE - K by 1 vector of model types:
%                    1 - Model evolves in absolute values
%                    2 - Model evolves in squares [DEFAULT]
%   STARTINGVALS - [OPTIONAL] K+sum(P)+sum(O)+sum(Q) vector of starting values
%
% OUTPUTS:
%   H          - T by K matrix of conditional variances
%   UNIVARIATE - A cell array of structures used to reconstruct the variance
%
% COMMENTS:
%
%  See also TARCH, DCC, CCC_MVGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 17/4/2012


if size(startingVals,2)>size(startingVals,1)
    startingVals = startingVals';
end

[T,k] = size(data2d);
H = zeros(T,k);
univariate = cell(k,1);
univariteOptions = optimset('fminunc');
univariteOptions.Display = 'none';
univariteOptions.LargeScale = 'off';
offset = 0;

for i=1:k
    if ~isempty(startingVals)
        count = 1+p(i)+o(i)+q(i);
        volStartingVals = startingVals(offset + (1:count));
        offset = offset + count;
    else
        volStartingVals = [];
    end
    %     [pr, ~, ~, ~, ~, scores, diagnostics] = tarch(data2d(:,i),p(i),o(i),q(i), [], gjrType(i), volStartingVals, univariteOptions);
    
    options0 = optimset('Display','off','TolCon',10^-3,'TolFun',10^-3,'TolX',10^-3,'Algorithm','interior-point');
    
    optPOT = @(params) Garchpot(params,data2d(:,i));
    [parametersall,~,~]  = fmincon(optPOT ,[cov(data2d(:,i))*0.05;0.03;0.9;0.5;0.1],[0,1,1,0,0],1,[],[],-0.5*ones(5,1),[0.5,0.5,0.99,0.9,0.9],[],options0);
    [~, ht,uL(i)] = Garchpot(parametersall,data2d(:,i));
    parameters = parametersall(1:end-2);
    parametersallPOT(:,i) = [parametersall;uL(i)];
    
    % Store output for later use
    univariate{i}.p = p(i);
    univariate{i}.o = o(i);
    univariate{i}.q = q(i);
    %     univariate{i}.fdata = diagnostics.fdata;
    %     univariate{i}.fIdata = diagnostics.fIdata;
    %     univariate{i}.back_cast = diagnostics.back_cast;
    univariate{i}.m = 1;
    univariate{i}.T = 1+T;
    univariate{i}.tarch_type = gjrType(i);
    univariate{i}.parameters = parameters;
    univariate{i}.ht = ht;
    %     univariate{i}.A = diagnostics.A;
    %     univariate{i}.scores = scores;
    H(:,i) = ht;
end
