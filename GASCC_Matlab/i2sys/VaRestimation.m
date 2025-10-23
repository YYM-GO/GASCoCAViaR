function VaRest = VaRestimation(data,THETA1,MODEL)
%  Estimate VaR for Xt
% MODEL =1 :GARCH-Gaussian      =2:GARCH-POT
[T,N] = size(data);
options0 = optimset('Display','off','TolCon',10^-3,'TolFun',10^-3,'TolX',10^-3,'Algorithm','interior-point');
parfor i = 1:N
    if MODEL == 1           %GARCH-Gaussian
        [~, ~, ~, ~, HNormal] = multigarch_AJP(data(:,i),1,0,1,'GJRGARCH','NORMAL',[],[cov(data(:,i))*0.05;0.05;0.9]);
        VaRest(:,i) = sqrt(HNormal)*norminv(THETA1);
    else            %GARCH-POT
        optPOT = @(params) Garchpot(params,data(:,i));
        paramestGPOT  = fmincon(optPOT ,[cov(data(:,i))*0.05;0.03;0.9;0.5;0.1],[0,1,1,0,0],1,[],[],-0.5*ones(5,1),[0.5,0.5,0.99,0.9,0.9],[],options0);
        [~, HPOT,uL1] = Garchpot(paramestGPOT,data(:,i));
        VaRest(:,i) = sqrt(HPOT)*POTinv(uL1,paramestGPOT(4),paramestGPOT(5),THETA1);
    end
end
end

