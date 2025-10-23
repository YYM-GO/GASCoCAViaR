function [pvalue,LR] = uccovar(x,xvar,y,yvar,THETA2)
% backtesting CoVaR in Kupiec (1995) 
yhit = y(x<xvar);
yvarhit = yvar(x<xvar);
n = sum(x<xvar);
pii = sum(yhit<yvarhit);
p = pii/n;
LR = -2*log((1-THETA2)^(n-pii)*THETA2^pii/(1-p)^(n-pii)/p^pii);
pvalue = 1-chi2cdf(LR,1);
end

