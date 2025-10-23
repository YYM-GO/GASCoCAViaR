function [VaRin,VaRforc,best] = VaROOSchoose(xvec,xnew,THETA1)
%calculate the VaR for insample and outsample period
%GMODEL: 1-GARCH-Gaussian   2-GARCH-POT   3-CAViaR-SAV   4-CAViaR-AS
% 1-GARCH-Gaussian
[paramG, ~, ~, ~, HNormal] = multigarch_AJP(xvec,1,0,1,'GJRGARCH','NORMAL',[],[cov(xvec)*0.05;0.05;0.9]);
VaRinG = sqrt(HNormal)*norminv(THETA1);
% 2-GARCH-POT
options0 = optimset('Display','off','TolCon',10^-3,'TolFun',10^-3,'TolX',10^-3,'Algorithm','interior-point');
optPOT = @(params) Garchpot(params,xvec);
paramGPOT  = fmincon(optPOT ,[cov(xvec)*0.05;0.03;0.9;0.5;0.1],[0,1,1,0,0],1,[],[],-0.5*ones(5,1),[0.5,0.5,0.99,0.9,0.9],[],options0);
[~, HPOT,uL1] = Garchpot(paramGPOT,xvec);
VaRinP = sqrt(HPOT)*POTinv(uL1,paramGPOT(4),paramGPOT(5),THETA1);
% 3-CAViaR-SAV
[BetaHatS,~,VaRin0S,~] = CAViaR2(1,THETA1,xvec);
VaRinS = -VaRin0S;
% 4-CAViaR-AS
[BetaHatA,~,VaRin0A,~] = CAViaR2(2,THETA1,xvec);
VaRinA = -VaRin0A;
for t = 1:length(xvec)
    S1(t) = SVaR3(xvec(t),VaRinG(t),THETA1);
    S2(t) = SVaR3(xvec(t),VaRinP(t),THETA1);
    S3(t) = SVaR3(xvec(t),VaRinS(t),THETA1);
    S4(t) = SVaR3(xvec(t),VaRinA(t),THETA1);
end
all = [mean(S1),mean(S2),mean(S3),mean(S4)];
best = find(all == min(all));

if best == 1     %GARCH-Gaussian
    VaRin = VaRinG;
    hN(1) = HNormal(end);
    for k = 2:length(xnew)
        hN(k) = paramG(1)+paramG(2)*xnew(k-1)^2+paramG(3)*hN(k-1);
    end
    VaRforc = norminv(THETA1)*sqrt(hN(2:end))';
else if best == 2     %GARCH-POT
        VaRin = VaRinP;
        hP(1) = HPOT(end);
        for k = 2:length(xnew)
            hP(k) = paramGPOT(1)+paramGPOT(2)*xnew(k-1)^2+paramGPOT(3)*hP(k-1);
        end
        VaRforc = sqrt(hP(2:end))'*POTinv(uL1,paramGPOT(4),paramGPOT(5),THETA1);
        
    else if best == 3     %CAViaR-SAV
            VaRin = VaRinS;
            VaRforc0 = RQobjectiveFunction(BetaHatS, 2, 1, length(xnew), xnew, THETA1, VaRinS(end));
            VaRforc = -VaRforc0(2:end,1);
        else               %CAViaR-AS
            VaRin = VaRinA;
            VaRforc0 = RQobjectiveFunction(BetaHatA, 2, 2, length(xnew), xnew, THETA1, VaRinA(end));
            VaRforc = -VaRforc0(2:end,1);
        end
    end
end
