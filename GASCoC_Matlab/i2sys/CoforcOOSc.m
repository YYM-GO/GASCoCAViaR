function [CoVoRforc,CoESforc,CoVaRin,CoESin,param] = CoforcOOSc(xvec,yvec,VaRin,VaRforc,THETA1,THETA2,CoMODEL,epsil,miter,xnew,ynew)
%calculate CoVaR and CoES OOS forecasts

[param,CoVaRi,CoESi] = constOOS(xvec,yvec,-VaRin,THETA1,THETA2,CoMODEL,epsil,miter);
beta = param(2:end);
CoVaR1 = -quantile(yvec((xvec<VaRin)),THETA2);
es = mean(yvec(xvec<VaRin));
[~,~,delta] = denALconst(ynew,param(1),THETA2,-CoVaRi,es);
[~,zt] = uandz(ynew,-CoVaRi,delta,THETA2);
txlv = find(xnew(2:end)<VaRforc);
if isnan(txlv) == 1
    txlv = 1;
end
[~,~,~,CoVoRforc,CoESforc0] = ElcALGAS6const(param,ynew,zt,xnew(2:end),VaRforc,-CoVaR1,es,THETA1,THETA2,CoMODEL,txlv);
                
CoESforc = CoESforc0';
CoVaRin = -CoVaRi;
CoESin = -CoESi;

end

