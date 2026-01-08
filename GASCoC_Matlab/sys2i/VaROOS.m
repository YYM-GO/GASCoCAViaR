function [VaRin,VaRforc] = VaROOS(xvec,xnew,GMODEL,THETA1)
%calculate the VaR for insample and outsample period
%GMODEL: 1-GARCH-Gaussian   2-GARCH-POT   3-CAViaR-SAV   4-CAViaR-AS
if GMODEL == 1     %GARCH-Gaussian
    [paramG, ~, ~, ~, HNormal] = multigarch_AJP(xvec,1,0,1,'GJRGARCH','NORMAL',[],[cov(xvec)*0.05;0.05;0.9]);
    VaRin = sqrt(HNormal)*norminv(THETA1);
    hN(1) = HNormal(end);
    for k = 2:length(xnew)
        hN(k) = paramG(1)+paramG(2)*xnew(k-1)^2+paramG(3)*hN(k-1);
    end
    VaRforc = norminv(THETA1)*sqrt(hN(2:end))';
else if GMODEL == 2     %GARCH-POT
        options0 = optimset('Display','off','TolCon',10^-3,'TolFun',10^-3,'TolX',10^-3,'Algorithm','interior-point');
        optPOT = @(params) Garchpot(params,xvec);
        paramGPOT  = fmincon(optPOT ,[cov(xvec)*0.05;0.03;0.9;0.5;0.1],[0,1,1,0,0],1,[],[],-0.5*ones(5,1),[0.5,0.5,0.99,0.9,0.9],[],options0);
        [~, HPOT,uL1] = Garchpot(paramGPOT,xvec);
        VaRin = sqrt(HPOT)*POTinv(uL1,paramGPOT(4),paramGPOT(5),THETA1);
        hP(1) = HPOT(end);
        for k = 2:length(xnew)
            hP(k) = paramGPOT(1)+paramGPOT(2)*xnew(k-1)^2+paramGPOT(3)*hP(k-1);
        end
        VaRforc = sqrt(hP(2:end))'*POTinv(uL1,paramGPOT(4),paramGPOT(5),THETA1);
    elseif GMODEL == 3     %CAViaR-SAV
        [BetaHat,~,VaRin0,~] = CAViaR2(1,THETA1,xvec);
        VaRin = -VaRin0;
        VaRforc0 = RQobjectiveFunction(BetaHat, 2, 1, length(xnew), xnew, THETA1, VaRin(end));
        VaRforc = -VaRforc0(2:end,1);
    else               %CAViaR-AS
        [BetaHat,~,VaRin0,~] = CAViaR2(2,THETA1,xvec);
        VaRin = -VaRin0;
        VaRforc0 = RQobjectiveFunction(BetaHat, 2, 2, length(xnew), xnew, THETA1, VaRin(end));
        VaRforc = -VaRforc0(2:end,1);
    end
end
end
