function [CoVoRforc,CoESforc,CoVaRin,CoESin,param,delta] = CoforcOOS(xvec,yvec,VaRin,VaRforc,THETA1,THETA2,CoMODEL,p,epsil,miter,xnew,ynew)
%calculate CoVaR and CoES OOS forecasts

[param,CoVaRi,CoESi,gammai,~] = GASGARCH1OOS(xvec,yvec,-VaRin,THETA1,THETA2,CoMODEL,p,epsil,miter);
beta = param(4:end);
CoVaR1 = -quantile(yvec((xvec<VaRin)),THETA2);
es = mean(yvec(xvec<VaRin));
if CoMODEL == 1
    [~,CoVaRnew] = CAVSAV27(xnew(2:end),ynew,beta,-CoVaR1,THETA2,VaRforc);
else
    [~,CoVaRnew] = CAVAS28(xnew(2:end),ynew,beta,-CoVaR1,THETA2,VaRforc);
end
CoVaRnew = -CoVaRnew;
delta0 = -THETA2*(1+exp(gammai(end)))*(-CoVaR1);
[~,~,delta] = denAL2new(ynew,param(1:3),gammai(end),delta0,THETA2,CoVaRnew,es,p);
[~,zt] = uandz(ynew,CoVaRnew,delta,THETA2);
txlv = find(xnew(2:end)<VaRforc);
if isnan(txlv) == 1
    txlv = 1;
end
[~,~,~,CoVoRforc,CoESforc0] = ElcALGAS6(param,ynew,zt,xnew(2:end),VaRforc,-CoVaR1,es,THETA1,THETA2,p,CoMODEL,gammai(end),delta0,txlv);
CoESforc = CoESforc0';
CoVaRin = -reshape(CoVaRi,[],1);
CoESin = -reshape(CoESi,[],1);
delta = delta';
end

