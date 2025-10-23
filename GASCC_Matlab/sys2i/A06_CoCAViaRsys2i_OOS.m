close
clear
clc


load shenwan       %Shenwan industry index
load hushen        %CSI 300




[T,N] = size(shenwan);
THETA1 = 0.05; THETA2 = 0.05; %probability level for VaR and CoVaR
wind = 1400;
Tout = T-wind;
stepst0 = 100;
nstepst = ceil(Tout/stepst0);
p = 0;      % the scale of the score, we can choose 0.
epsil = 1e-6;
TolX = 1e-6;
miter = 100;


%%                  GAS-CoCAViaR

GMODEL = 2;          %GMODEL: 1-GARCH-Gaussian   2-GARCH-POT
CoMODEL = 1;         %CoCAViaR: 1-CoSAV     2-CoAS             (Need to run all combinations: GMODEL+CoGMODEL = 1+1/ 1+2/ 2+1/ 2+2)

VaRforcn = [];  CoVoRforcn = [];  CoESforcn = []; deltan = [];  paramn = [];  param = []; deltaforc = [];
paramsys2ioos = []; VaRin = []; VaRforc = [];
CoVaRin = []; CoESin = [];
VaRsys2ioos = zeros(Tout,N);
CoVaRsys2ioos = zeros(Tout,N);
CoESsys2ioos = zeros(Tout,N);
deltasys2ioos = zeros(Tout,N);
corrmatsys2ioos = zeros(N,N,Tout);
condiexpsys2ioos = zeros(Tout,N);
ksen2 = (1-THETA2*2)/(THETA2*(1-THETA2));
sig2 = sqrt(2/(THETA2*(1-THETA2)));
for t = 1:nstepst
    txlv = []; w = []; wm = []; Sigma = [];
    tic
    if t ~= nstepst
        stepst = stepst0;
    else
        stepst = Tout-(nstepst-1)*stepst0;
    end
    xvec = hushen(stepst*(t-1)+1:stepst*(t-1)+wind);
    xvecall = [xvec;hushen(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst)];
    xnew = xvecall(end-stepst:end);
    [VaRin,VaRforc] = VaROOS(xvec,xnew,GMODEL,THETA1); warning off
    parfor i = 1:N
        yvec = shenwan(stepst*(t-1)+1:stepst*(t-1)+wind,i);
        yvecall = [yvec;shenwan(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst,i)];
        ynew = yvecall(end-stepst+1:end);
        [CoVoRforc,CoESforc,CoVaRin,CoESin,param,deltain,deltaforc,condiexp] = CoforcOOS(xvec,yvec,VaRin,VaRforc,THETA1,THETA2,CoMODEL,p,epsil,miter,xnew,ynew);
        VaRforcn(:,i) = VaRforc;
        CoVoRforcn(:,i) = CoVoRforc;
        CoESforcn(:,i) = CoESforc;
        paramn(:,i) = param;
        CoVaRinall(:,i) = CoVaRin;
        deltainall(:,i) = deltain;
        condiexpall(:,i) = condiexp;
        deltan(:,i) = deltaforc;
    end
    ymat =  shenwan(stepst*(t-1)+1:stepst*(t-1)+wind,:);
    txlv = find(xvec<VaRin);
    w1=max(mean((ymat(txlv,:)-CoVaRinall(txlv,:))./(deltainall(txlv,:)*ksen2),2),0.001);
    for k = 1:length(txlv)
        wm1(k,:) = (deltainall(txlv(k),:)*ksen2)*w1(k);
        Sigma0(k,:) = (ymat(txlv(k),:)-CoVaRinall(txlv(k),:)-wm1(k,:))/sqrt(w1(k));
        Sigma1(k,:) = Sigma0(k,:)./(deltainall(txlv(k),:).^2)./(sig2*ones(1,N).^2);
    end
    corrmatin =  corr(Sigma1);
    if t ~= nstepst
        VaRsys2ioos(stepst*(t-1)+1:stepst*t,:) = VaRforcn;
        CoVaRsys2ioos(stepst*(t-1)+1:stepst*t,:) = CoVoRforcn;
        CoESsys2ioos(stepst*(t-1)+1:stepst*t,:) = CoESforcn;
        deltasys2ioos(stepst*(t-1)+1:stepst*t,:) = deltan;
        corrmatsys2ioos(:,:,stepst*(t-1)+1:stepst*t) = repmat(corrmatin,1,1,stepst);
        condiexpsys2ioos(stepst*(t-1)+1:stepst*t,:) = repmat(condiexpall,stepst,1);
    else
        VaRsys2ioos(end-stepst+1:end,:) = VaRforcn;
        CoVaRsys2ioos(end-stepst+1:end,:) = CoVoRforcn;
        CoESsys2ioos(end-stepst+1:end,:) = CoESforcn;
        deltasys2ioos(end-stepst+1:end,:) = deltan;
        corrmatsys2ioos(:,:,end-stepst+1:end) = repmat(corrmatin,1,1,stepst);
        condiexpsys2ioos(end-stepst+1:end,:) = repmat(condiexpall,stepst,1);
    end
    [t,toc]
end



if GMODEL == 1
    if CoMODEL == 1
        VaRGauSAVooss2i = VaRsys2ioos;
        CoVaRGauSAVooss2i = CoVaRsys2ioos;
        CoESGauSAVooss2i = CoESsys2ioos;
        paramGauSAVooss2i = paramsys2ioos;
        deltaGauSAVooss2i = deltasys2ioos;
        corrmatGauSAVooss2i = corrmatsys2ioos;
        condiexpGauSAVooss2i = condiexpsys2ioos;
    else
        VaRGauASooss2i = VaRsys2ioos;
        CoVaRGauASooss2i = CoVaRsys2ioos;
        CoESGauASooss2i = CoESsys2ioos;
        paramGauASooss2i = paramsys2ioos;
        deltaGauASooss2i = deltasys2ioos;
        corrmatGauASooss2i = corrmatsys2ioos;
        condiexpGauASooss2i = condiexpsys2ioos;
    end
else
    if CoMODEL == 1
        VaRPOTSAVooss2i = VaRsys2ioos;
        CoVaRPOTSAVooss2i = CoVaRsys2ioos;
        CoESPOTSAVooss2i = CoESsys2ioos;
        paramPOTSAVooss2i = paramsys2ioos;
        deltaPOTSAVooss2i = deltasys2ioos;
        corrmatPOTSAVooss2i = corrmatsys2ioos;
        condiexpPOTSAVooss2i = condiexpsys2ioos;
    else
        VaRPOTASooss2i = VaRsys2ioos;
        CoVaRPOTASooss2i = CoVaRsys2ioos;
        CoESPOTASooss2i = CoESsys2ioos;
        paramPOTASooss2i = paramsys2ioos;
        deltaPOTASooss2i = deltasys2ioos;
        corrmatPOTASooss2i = corrmatsys2ioos;
        condiexpPOTASooss2i = condiexpsys2ioos;
    end
end



save VaRGauSAVooss2i.mat VaRGauSAVooss2i
save CoVaRGauSAVooss2i.mat CoVaRGauSAVooss2i
save CoESGauSAVooss2i.mat CoESGauSAVooss2i
save deltaGauSAVooss2i.mat deltaGauSAVooss2i
save corrmatGauSAVooss2i.mat corrmatGauSAVooss2i
save condiexpGauSAVooss2i.mat condiexpGauSAVooss2i


save VaRGauASooss2i.mat VaRGauASooss2i
save CoVaRGauASooss2i.mat CoVaRGauASooss2i
save CoESGauASooss2i.mat CoESGauASooss2i
save deltaGauASooss2i.mat deltaGauASooss2i
save corrmatGauASooss2i.mat corrmatGauASooss2i
save condiexpGauASooss2i.mat condiexpGauASooss2i


save VaRPOTSAVooss2i.mat VaRPOTSAVooss2i
save CoVaRPOTSAVooss2i.mat CoVaRPOTSAVooss2i
save CoESPOTSAVooss2i.mat CoESPOTSAVooss2i
save deltaPOTSAVooss2i.mat deltaPOTSAVooss2i
save corrmatPOTSAVooss2i.mat corrmatPOTSAVooss2i
save condiexpPOTSAVooss2i.mat condiexpPOTSAVooss2i


save VaRPOTASooss2i.mat VaRPOTASooss2i
save CoVaRPOTASooss2i.mat CoVaRPOTASooss2i
save CoESPOTASooss2i.mat CoESPOTASooss2i
save paramPOTASooss2i.mat paramPOTASooss2i
save deltaPOTASooss2i.mat deltaPOTASooss2i
save corrmatPOTASooss2i.mat corrmatPOTASooss2i
save condiexpPOTASooss2i.mat condiexpPOTASooss2i

%%               CoCAViaR-c
GMODEL = 2;          %GMODEL: 1-GARCH-Gaussian   2-GARCH-POT
CoMODEL = 1;         %CoCAViaR: 1-CoSAV     2-CoAS                    (Need to run all combinations: GMODEL+CoGMODEL = 1+1/ 1+2/ 2+1/ 2+2)

VaRforcnc = [];  CoVoRforcnc = [];  CoESforcnc = []; deltanc = [];  paramnc = [];  param = []; deltaforc = [];
paramsys2ioosc = [];  VaRin = []; VaRforc = [];
CoVaRin = []; CoESin = [];
VaRsys2ioosc = zeros(Tout,N);
CoVaRsys2ioosc = zeros(Tout,N);
CoESsys2ioosc = zeros(Tout,N);
deltasys2ioosc = zeros(Tout,N);
corrmatsys2ioosc = zeros(N,N,Tout);
condiexpsys2ioosc = zeros(Tout,N);
for t = 1:nstepst
    tic
    if t ~= nstepst
        stepst = stepst0;
    else
        stepst = Tout-(nstepst-1)*stepst0;
    end
    xvec = hushen(stepst*(t-1)+1:stepst*(t-1)+wind);
    xvecall = [xvec;hushen(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst)];
    xnew = xvecall(end-stepst:end);
    [VaRin,VaRforc] = VaROOS(xvec,xnew,GMODEL,THETA1); warning off    
    parfor i = 1:N
        yvec = shenwan(stepst*(t-1)+1:stepst*(t-1)+wind,i);
        yvecall = [yvec;shenwan(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst,i)];
        ynew = yvecall(end-stepst+1:end);
        [CoVoRforc,CoESforc,CoVaRin,CoESin,param,deltain,deltaforc,condiexp] = CoforcOOSc(xvec,yvec,VaRin,VaRforc,THETA1,THETA2,CoMODEL,epsil,miter,xnew,ynew);
        VaRforcnc(:,i) = VaRforc;
        CoVoRforcnc(:,i) = CoVoRforc;
        CoESforcnc(:,i) = CoESforc;
        paramnc(:,i) = param;
        CoVaRinallc(:,i) = CoVaRin;
        deltainallc(:,i) = deltain;
        condiexpallc(:,i) = condiexp;
        deltanc(:,i) = deltaforc;
    end
    ymat =  shenwan(stepst*(t-1)+1:stepst*(t-1)+wind,:);
    txlv = find(xvec<VaRin);
    w1=max(mean((ymat(txlv,:)-CoVaRinallc(txlv,:))./(deltainallc(txlv,:)*ksen2),2),0.001);
    for k = 1:length(txlv)
        wm1(k,:) = (deltainallc(txlv(k),:)*ksen2)*w1(k);
        Sigma0(k,:) = (ymat(txlv(k),:)-CoVaRinallc(txlv(k),:)-wm1(k,:))/sqrt(w1(k));
        Sigma1(k,:) = Sigma0(k,:)./(deltainallc(txlv(k),:).^2)./(sig2*ones(1,N).^2);
    end
    corrmatinc =  corr(Sigma1);
    if t ~= nstepst
        VaRsys2ioosc(stepst*(t-1)+1:stepst*t,:) = VaRforcnc;
        CoVaRsys2ioosc(stepst*(t-1)+1:stepst*t,:) = CoVoRforcnc;
        CoESsys2ioosc(stepst*(t-1)+1:stepst*t,:) = CoESforcnc;
        deltasys2ioosc(stepst*(t-1)+1:stepst*t,:) = deltanc;
        corrmatsys2ioosc(:,:,stepst*(t-1)+1:stepst*t) = repmat(corrmatinc,1,1,stepst);
        condiexpsys2ioosc(stepst*(t-1)+1:stepst*t,:) = repmat(condiexpallc,stepst,1);
    else
        VaRsys2ioosc(end-stepst+1:end,:) = VaRforcnc;
        CoVaRsys2ioosc(end-stepst+1:end,:) = CoVoRforcnc;
        CoESsys2ioosc(end-stepst+1:end,:) = CoESforcnc;
        deltasys2ioosc(end-stepst+1:end,:) = deltanc;
        corrmatsys2ioosc(:,:,end-stepst+1:end) = repmat(corrmatinc,1,1,stepst);
        condiexpsys2ioosc(end-stepst+1:end,:) = repmat(condiexpallc,stepst,1);
    end
    [t,toc]
end


if GMODEL == 1
    if CoMODEL == 1
        VaRGauSAVooss2ic = VaRsys2ioosc;
        CoVaRGauSAVooss2ic = CoVaRsys2ioosc;
        CoESGauSAVooss2ic = CoESsys2ioosc;
        paramGauSAVooss2ic = paramsys2ioosc;
        deltaGauSAVooss2ic = deltasys2ioosc;
        corrmatGauSAVooss2ic = corrmatsys2ioosc;
        condiexpGauSAVooss2ic = condiexpsys2ioosc;
    else
        VaRGauASooss2ic = VaRsys2ioosc;
        CoVaRGauASooss2ic = CoVaRsys2ioosc;
        CoESGauASooss2ic = CoESsys2ioosc;
        paramGauASooss2ic = paramsys2ioosc;
        deltaGauASooss2ic = deltasys2ioosc;
        corrmatGauASooss2ic = corrmatsys2ioosc;
        condiexpGauASooss2ic = condiexpsys2ioosc;
    end
else
    if CoMODEL == 1
        VaRPOTSAVooss2ic = VaRsys2ioosc;
        CoVaRPOTSAVooss2ic = CoVaRsys2ioosc;
        CoESPOTSAVooss2ic = CoESsys2ioosc;
        paramPOTSAVooss2ic = paramsys2ioosc;
        deltaPOTSAVooss2ic = deltasys2ioosc;
        corrmatPOTSAVooss2ic = corrmatsys2ioosc;
        condiexpPOTSAVooss2ic = condiexpsys2ioosc;
    else
        VaRPOTASooss2ic = VaRsys2ioosc;
        CoVaRPOTASooss2ic = CoVaRsys2ioosc;
        CoESPOTASooss2ic = CoESsys2ioosc;
        paramPOTASooss2ic = paramsys2ioosc;
        deltaPOTASooss2ic = deltasys2ioosc;
        corrmatPOTASooss2ic = corrmatsys2ioosc;
        condiexpPOTASooss2ic = condiexpsys2ioosc;
    end
end


save VaRGauSAVooss2ic.mat VaRGauSAVooss2ic
save CoVaRGauSAVooss2ic.mat CoVaRGauSAVooss2ic
save CoESGauSAVooss2ic.mat CoESGauSAVooss2ic
save deltaGauSAVooss2ic.mat deltaGauSAVooss2ic
save corrmatGauSAVooss2ic.mat corrmatGauSAVooss2ic
save condiexpGauSAVooss2ic.mat condiexpGauSAVooss2ic


save VaRGauASooss2ic.mat VaRGauASooss2ic
save CoVaRGauASooss2ic.mat CoVaRGauASooss2ic
save CoESGauASooss2ic.mat CoESGauASooss2ic
save deltaGauASooss2ic.mat deltaGauASooss2ic
save corrmatGauASooss2ic.mat corrmatGauASooss2ic
save condiexpGauASooss2ic.mat condiexpGauASooss2ic


save VaRPOTSAVooss2ic.mat VaRPOTSAVooss2ic
save CoVaRPOTSAVooss2ic.mat CoVaRPOTSAVooss2ic
save CoESPOTSAVooss2ic.mat CoESPOTSAVooss2ic
save deltaPOTSAVooss2ic.mat deltaPOTSAVooss2ic
save corrmatPOTSAVooss2ic.mat corrmatPOTSAVooss2ic
save condiexpPOTSAVooss2ic.mat condiexpPOTSAVooss2ic


save VaRPOTASooss2ic.mat VaRPOTASooss2ic
save CoVaRPOTASooss2ic.mat CoVaRPOTASooss2ic
save CoESPOTASooss2ic.mat CoESPOTASooss2ic
save deltaPOTASooss2ic.mat deltaPOTASooss2ic
save corrmatPOTASooss2ic.mat corrmatPOTASooss2ic
save condiexpPOTASooss2ic.mat condiexpPOTASooss2ic