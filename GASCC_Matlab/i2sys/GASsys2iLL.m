function ElAL = GASsys2iLL(params,xvec,yvec,xvar,THETA1,THETA2,CoMODEL,p)
empiricalQuantile = quantile(yvec(xvec<-xvar),THETA2);
if CoMODEL == 1
    [~,CoVaR] = CAVSAV27(xvec,yvec,params(4:end),empiricalQuantile,THETA2,xvar); 
else
    [~,CoVaR] = CAVAS28(xvec,yvec,params(4:end),empiricalQuantile,THETA2,xvar); 
end
txlv = find(xvec<-xvar);
gamma0 = log(mean(yvec(find(yvec<quantile(yvec(xvec<quantile(xvec,THETA1)),THETA2))))/quantile(yvec(xvec<quantile(xvec,THETA1)),THETA2)-1);
delta0 = -THETA2*(1+exp(gamma0))*(-CoVaR(1));
es = mean(yvec(find(xvec<-xvar)));
[f,gamma,delta] = denAL2new(yvec,params(1:3),gamma0,delta0,THETA2,-CoVaR,es,p);
[ut,zt] = uandz(yvec,-CoVaR,delta,THETA2);
ElAL = ElcALGAS5new(params,yvec,zt,xvec,-xvar,es,THETA1,THETA2,p,CoMODEL,gamma0,delta0,txlv);
end