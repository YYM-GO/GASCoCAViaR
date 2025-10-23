function [RQ,VaR] = CAVAS28(xvec,yvec,beta,empiricalQuantile,THETA,xvar)
VaRtest2(1) = -empiricalQuantile;
for i = 2:length(yvec)
    VaRtest2(i) = beta(1)+beta(2)*xvec(i-1)*(xvec(i-1)>0)-beta(3)*xvec(i-1)*(xvec(i-1)<0)+beta(4)*yvec(i-1)*(yvec(i-1)>0)-beta(5)*yvec(i-1)*(yvec(i-1)<0)+beta(6)*VaRtest2(i-1);
end
VaR =  VaRtest2;
VaR = VaR';
HitQ = -(xvec<-xvar).*((yvec<-VaR)-THETA).*(yvec+VaR);
RQ  = sum(HitQ);
if RQ == Inf | (RQ ~= RQ) | ~isreal(RQ)
    RQ = 1e+100;
end
end

