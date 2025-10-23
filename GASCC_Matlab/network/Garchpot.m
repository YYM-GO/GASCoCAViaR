function [out,sigma2,uL] = Garchpot(params,data)
%likelihood function of GARCH POT
g1 = params(1);
g2 = params(2);
g3 = params(3);
betaL = params(4);
xiL = params(5);

T = length(data);
sigma2 = zeros(T,1);
sigma2(1) = var(data);
for t = 2:T
    sigma2(t) = g1+g2*data(t-1)^2+g3*sigma2(t-1);
end
res = data./sqrt(sigma2);
uL = quantile(res,0.1);

for t = 1:T
    if res(t)<uL
        LL(t) = -log(POTpdf2(res(t),uL,round(T*0.1),T,betaL,xiL)/sqrt(sigma2(t)));
    else
        LL(t) = -log(normpdf(res(t))/sqrt(sigma2(t)));
    end
end

out = sum(LL);
end


