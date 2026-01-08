function [BetaHat,RQ,VaR] = optcavsav27(THETA,xvec,yvec,xvar)

y=yvec;
T = size(yvec,1);

nInitialVectors = [10000, 4]; % Number of initial vector fed in the uniform random number generator COSAV  models.
nInitialCond = 30;
MaxFunEvals = 5000; % Parameters for the optimisation algorithm. Increase them in case the algorithm does not converge.
MaxIter     = 5000;
REP = 10;
options = optimset('LargeScale', 'off', 'HessUpdate', 'dfp','MaxFunEvals', MaxFunEvals, ...
    'display', 'off', 'MaxIter', MaxIter, 'TolFun', 1e-8, 'TolX', 1e-8);
rng(100,'philox')
VaR            = zeros(size(y));
Hit            = VaR;
DQinSample     = zeros(1, T);
DQoutOfSample  = DQinSample;
nSamples = 1;
empiricalQuantile = quantile(yvec(xvec<-xvar),THETA);
initialTargetVectors = unifrnd(-1, 1, nInitialVectors);
RQfval = zeros(nInitialVectors(1), 1);
for i = 1:nInitialVectors(1)
    RQfval(i) = CAVSAV27(xvec,yvec,initialTargetVectors(i,:),empiricalQuantile,THETA,xvar);
end
Results          = [RQfval, initialTargetVectors];
SortedResults    = sortrows(Results,1);
BestInitialCond  = SortedResults(1:nInitialCond,2:end);
for i = 1:size(BestInitialCond,1)
    flike = @(beta)CAVSAV27(xvec,yvec,beta,empiricalQuantile,THETA,xvar);
    [Beta(i,:), fval(i,1), exitflag(i,1)] = fmincon(flike,BestInitialCond(i,:),[],[],[],[],[-0.3*ones(1,3),0.7],[ones(1,4)],[],options);
    for it = 1:REP
        [Beta(i,:), fval(i,1), exitflag(i,1)] = fminunc(flike,Beta(i,:),options);
        [Beta(i,:), fval(i,1), exitflag(i,1)] = fminsearch(flike,Beta(i,:),options);
        if exitflag(i,1) == 1
            break
        end
    end
end
SortedFval  = sortrows([fval, Beta], 1);
BetaHat    = SortedFval(1, 2:end)';
[RQ,VaR] = CAVSAV27(xvec,yvec,BetaHat,empiricalQuantile,THETA,xvar);
end

