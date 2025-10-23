function VaRest = VaRestimation(data,THETA1,MODEL)
%Estimate VaR
[T,N] = size(data);
options0 = optimset('Display','off','TolCon',10^-3,'TolFun',10^-3,'TolX',10^-3,'Algorithm','interior-point');
if MODEL == 1    % GARCH-Gaussian
    [~, ~, ~, ~, HNormal] = multigarch_AJP(data,1,0,1,'GJRGARCH','NORMAL',[],[cov(data)*0.05;0.05;0.9]);
    VaRest = sqrt(HNormal)*norminv(THETA1);
else
    optPOT = @(params) Garchpot(params,data);
    paramestGPOT  = fmincon(optPOT ,[cov(data)*0.05;0.03;0.9;0.5;0.1],[0,1,1,0,0],1,[],[],-0.5*ones(5,1),[0.5,0.5,0.99,0.9,0.9],[],options0);
    [~, HPOT,uL1] = Garchpot(paramestGPOT,data);
    VaRest = sqrt(HPOT)*POTinv(uL1,paramestGPOT(4),paramestGPOT(5),THETA1);
end
end

