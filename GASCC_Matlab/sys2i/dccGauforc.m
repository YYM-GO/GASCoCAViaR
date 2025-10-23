function [DCCVaR,DCCCoVaR,DCCCoES,condvar2out,stdresidsIS2] = dccGauforc(xvec,yvec,xnew,ynew,THETA1,THETA2,stepst)
%calculate the VaR, CoVaR and CoES forcasts using DCC-GARCH model with
%Gaussian innovations.

rates = [xvec,yvec];
ratesnew = [xnew,ynew];
MUhat = [0,0]; % only GARCH, No ARMA
[parameters, ~ ,Ht,Qt,Rt] = dcc(rates,[],1,0,1);
stdresidsIS = rates./ [squeeze(sqrt(Ht(1,1,:))),squeeze(sqrt(Ht(2,2,:)))];
stdresidsIS2 = stdresidsIS(:,2);
Htnew(:,:,1) = Ht(:,:,end);
Qtnew(:,:,1) = Qt(:,:,end);
Rtnew(:,:,1) = Rt(:,:,end);
for k =2:(stepst+1)
    Htnew(1,1,k) = parameters(1)+parameters(2)*xnew(k-1)^2+parameters(3)*Htnew(1,1,k-1);
    Htnew(2,2,k) = parameters(4)+parameters(5)*ynew(k-1)^2+parameters(6)*Htnew(2,2,k-1);
    e = ratesnew(k-1,:)'./diag(sqrt(Htnew(:,:,k-1)));
    Qtnew(:,:,k) = (1-parameters(end-1)-parameters(end))*[1,parameters(end-2);parameters(end-2),1]+parameters(end-1)*e*e'+parameters(end)*Qtnew(:,:,k-1);
    Rtnew(:,:,k) = [1,Qtnew(1,2,k)/(sqrt(Qtnew(1,1,k))*sqrt(Qtnew(2,2,k)));Qtnew(1,2,k)/(sqrt(Qtnew(1,1,k))*sqrt(Qtnew(2,2,k))),1];
    Htnew(:,:,k) = diag([sqrt(Htnew(1,1,k)),sqrt(Htnew(2,2,k))])*Rtnew(:,:,k)*diag([sqrt(Htnew(1,1,k)),sqrt(Htnew(2,2,k))]);
    DCCVaR0(k) = sqrt(Htnew(1,1,k))*norminv(THETA1);
    %             DCCES0(k,i) = - sqrt(Htnew(1,1,k+1))*normpdf(norminv(THETA1))/(THETA1);
    DCCCoVaR0(k) = CoVaRfrocNormalN(THETA1,THETA2,Rtnew(1,2,k),diag(Htnew(:,:,k)),ones(stepst+1,1),MUhat);
    grid = [0.005:0.01:THETA2];
    for l = 1:length(grid)
        esgrid(l) = CoVaRfrocNormalN(THETA1,grid(l),Rtnew(1,2,k),diag(Htnew(:,:,k)),ones(stepst+1,1),MUhat); warning off
    end
    DCCCoES0(k) = mean(esgrid);
end
DCCVaR = DCCVaR0(2:end)';
DCCCoVaR = DCCCoVaR0(2:end)';
DCCCoES = DCCCoES0(2:end)';
condvar2out = squeeze(Htnew(2,2,2:end));
end

