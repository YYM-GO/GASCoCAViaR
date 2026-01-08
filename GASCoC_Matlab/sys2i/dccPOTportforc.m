function [DCCCoVaR,DCCCoES] = dccPOTportforc(xvec,yvec,xnew,ynew,THETA1,THETA2,tforc)
%calculate the VaR, CoVaR and CoES forcasts using DCC-GARCH model with
%Gaussian innovations.

rates = [xvec,yvec];
ratesnew = [xnew,ynew];
MUhat = [0,0]; % only GARCH, No ARMA
[parameters, ~ ,Ht,Qt,Rt] = dccPOTport(rates,[],1,0,1);
Htnew(:,:,1) = Ht(:,:,end);
Qtnew(:,:,1) = Qt(:,:,end);
Rtnew(:,:,1) = Rt(:,:,end);
for k =2:(tforc+1)
    Htnew(1,1,k) = parameters(1)+parameters(2)*xnew(k-1)^2+parameters(3)*Htnew(1,1,k-1);
    Htnew(2,2,k) = parameters(7)+parameters(8)*ynew(k-1)^2+parameters(9)*Htnew(2,2,k-1);
    e = ratesnew(k-1,:)'./diag(sqrt(Htnew(:,:,k-1)));
    Qtnew(:,:,k) = (1-parameters(end-1)-parameters(end))*[1,parameters(end-2);parameters(end-2),1]+parameters(end-1)*e*e'+parameters(end)*Qtnew(:,:,k-1);
    Rtnew(:,:,k) = [1,Qtnew(1,2,k)/(sqrt(Qtnew(1,1,k))*sqrt(Qtnew(2,2,k)));Qtnew(1,2,k)/(sqrt(Qtnew(1,1,k))*sqrt(Qtnew(2,2,k))),1];
    Htnew(:,:,k) = diag([sqrt(Htnew(1,1,k)),sqrt(Htnew(2,2,k))])*Rtnew(:,:,k)*diag([sqrt(Htnew(1,1,k)),sqrt(Htnew(2,2,k))]);
    DCCVaR0(k) = sqrt(Htnew(1,1,k))*POTinv(parameters(6),parameters(4),parameters(5),THETA1);
    DCCCoVaR0(k) = CoVaRfrocNormalPOT(THETA1,THETA2,Rtnew(1,2,k),diag((Htnew(:,:,k))),ones(tforc+1,1),MUhat,parameters(12),parameters(10),parameters(11)); warning off
    grid = [0.005:0.01:THETA2];
    for l = 1:length(grid)
        esgrid(l) = CoVaRfrocNormalPOT(THETA1,grid(l),Rtnew(1,2,k),diag((Htnew(:,:,k))),ones(tforc+1,1),MUhat,parameters(12),parameters(10),parameters(11)); warning off
    end
    DCCCoES0(k) = mean(esgrid);
end
DCCCoVaR = DCCCoVaR0(tforc+1);
DCCCoES = DCCCoES0(tforc+1);
end

