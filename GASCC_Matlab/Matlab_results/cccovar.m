function [pvalue,LR] = cccovar(x,xvar,y,yvar,p)
% backtesting CoVaR in  Girardi & Ergün
yhit = y(x<xvar);
yvarhit = yvar(x<xvar);
n00=0; n01=0; n10=0; n11=0;
for i = 2:length(yhit)
    if yhit(i-1)<yvarhit(i-1)
        if yhit(i)<yvarhit(i)
            n11 = n11+1;
        else
            n10 = n10+1;
        end
    else
        if yhit(i)<yvarhit(i)
            n01 = n01+1;
        else
            n00 = n00+1;
        end
    end
end

p00 = n00/(n00+n01);
p01 = n01/(n00+n01);
p10 = n10/(n10+n11);
p11 = n11/(n10+n11);
p2 = (n01+n11)/(n00+n01+n10+n11);
if n11 == 0
    L2 = (1-p)^(n00+n10)*p^(n01+n11);
else
    L2 = (1-p2)^(n00+n10)*p2^(n01+n11);
end

L1 = (1-p01)^n00*p01^n01*(1-p11)^n10*p11^n11;

LR = -2*log(L2/L1);

pvalue = 1-chi2cdf(LR,1);

end

