close
clear
clc

load shenwan       % Shenwan industry index
load hushen        % CSI 300




[T,N] = size(shenwan);
THETA1 = 0.05; THETA2 = 0.05; %probability level for VaR and CoVaR
wind = 1400;
Tout = T-wind;
stepst0 = 100;
nstepst = ceil(Tout/stepst0);
p = 0;      
epsil = 1e-6;
TolX = 1e-6;
miter = 100;

%%             GAS-COCAViaR

GMODEL = 2;          %GMODEL: 1-GARCH-Gaussian   2-GARCH-POT
CoMODEL = 1;         %CoCAViaR: 1-CoSAV     2-CoAS         (Need to run all combinations: GMODEL+CoGMODEL = 1+1/ 1+2/ 2+1/ 2+2)
VaRforcn = [];  CoVoRforcn = [];  CoESforcn = []; deltan = [];  paramn = [];  param = []; deltaforc = [];
parami2sysoos = []; VaRin = []; VaRforc = [];
CoVaRin = []; CoESin = [];
VaRi2sysoos = zeros(Tout,N);
CoVaRi2sysoos = zeros(Tout,N);
CoESi2sysoos = zeros(Tout,N);
deltai2sysoos = zeros(Tout,N);
for t = 1:nstepst                    % about 40-70 seconds for each t
    tic
    if t ~= nstepst
        stepst = stepst0;
    else
        stepst = Tout-(nstepst-1)*stepst0;
    end
    yvec = hushen(stepst*(t-1)+1:stepst*(t-1)+wind);
    yvecall = [yvec;hushen(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst)];
    ynew = yvecall(end-stepst+1:end);
    parfor i = 1:N
        xvec = shenwan(stepst*(t-1)+1:stepst*(t-1)+wind,i);
        xvecall = [xvec;shenwan(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst,i)];
        xnew = xvecall(end-stepst:end);
        [VaRin,VaRforc] = VaROOS(xvec,xnew,GMODEL,THETA1); warning off
        [CoVoRforc,CoESforc,CoVaRin,CoESin,param,deltaforc] = CoforcOOS(xvec,yvec,VaRin,VaRforc,THETA1,THETA2,CoMODEL,p,epsil,miter,xnew,ynew);
        VaRforcn(:,i) = VaRforc;
        CoVoRforcn(:,i) = CoVoRforc;
        CoESforcn(:,i) = CoESforc;
        paramn(:,i) = param;
        deltan(:,i) = deltaforc;
        %           i
    end
    
    if t ~= nstepst
        VaRi2sysoos(stepst*(t-1)+1:stepst*t,:) = VaRforcn;
        CoVaRi2sysoos(stepst*(t-1)+1:stepst*t,:) = CoVoRforcn;
        CoESi2sysoos(stepst*(t-1)+1:stepst*t,:) = CoESforcn;
        deltai2sysoos(stepst*(t-1)+1:stepst*t,:) = deltan;
    else
        VaRi2sysoos(end-stepst+1:end,:) = VaRforcn;
        CoVaRi2sysoos(end-stepst+1:end,:) = CoVoRforcn;
        CoESi2sysoos(end-stepst+1:end,:) = CoESforcn;
        deltai2sysoos(end-stepst+1:end,:) = deltan;
    end
    [t,toc]
end



if GMODEL == 1
    if CoMODEL == 1
        VaRGauSAVoosi2s = VaRi2sysoos;
        CoVaRGauSAVoosi2s = CoVaRi2sysoos;
        CoESGauSAVoosi2s = CoESi2sysoos;
        paramGauSAVoosi2s = parami2sysoos;
    else
        VaRGauASoosi2s = VaRi2sysoos;
        CoVaRGauASoosi2s = CoVaRi2sysoos;
        CoESGauASoosi2s = CoESi2sysoos;
        paramGauASoosi2s = parami2sysoos;
    end
else  
    if CoMODEL == 1
        VaRPOTSAVoosi2s = VaRi2sysoos;
        CoVaRPOTSAVoosi2s = CoVaRi2sysoos;
        CoESPOTSAVoosi2s = CoESi2sysoos;
        paramPOTSAVoosi2s = parami2sysoos;
    else
        VaRPOTASoosi2s = VaRi2sysoos;
        CoVaRPOTASoosi2s = CoVaRi2sysoos;
        CoESPOTASoosi2s = CoESi2sysoos;
        paramPOTASoosi2s = parami2sysoos;
    end   
end



save VaRGauSAVoosi2s.mat VaRGauSAVoosi2s
save CoVaRGauSAVoosi2s.mat CoVaRGauSAVoosi2s
save CoESGauSAVoosi2s.mat CoESGauSAVoosi2s

save VaRGauASoosi2s.mat VaRGauASoosi2s
save CoVaRGauASoosi2s.mat CoVaRGauASoosi2s
save CoESGauASoosi2s.mat CoESGauASoosi2s

save VaRPOTSAVoosi2s.mat VaRPOTSAVoosi2s
save CoVaRPOTSAVoosi2s.mat CoVaRPOTSAVoosi2s
save CoESPOTSAVoosi2s.mat CoESPOTSAVoosi2s

save VaRPOTASoosi2s.mat VaRPOTASoosi2s
save CoVaRPOTASoosi2s.mat CoVaRPOTASoosi2s
save CoESPOTASoosi2s.mat CoESPOTASoosi2s

%%      constant multiplicative factor
GMODEL = 2;          %GMODEL: 1-GARCH-Gaussian   2-GARCH-POT
CoMODEL = 1;         %CoCAViaR: 1-CoSAV     2-CoAS            (Need to run all combinations: GMODEL+CoGMODEL = 1+1/ 1+2/ 2+1/ 2+2)


VaRforcnc = [];  CoVoRforcnc = [];  CoESforcnc = []; deltanc = [];  paramnc = [];  param = []; deltaforc = []; 
parami2sysoosc = []; CoVaRi2sysoosc = []; VaRi2sysoosc = [];  CoESi2sysoosc = []; VaRin = []; VaRforc = [];
CoVaRin = []; CoESin = [];
VaRi2sysoosc = zeros(Tout,N);
CoVaRi2sysoosc = zeros(Tout,N);
CoESi2sysoosc = zeros(Tout,N);

for t = 1:nstepst
    tic
    if t ~= nstepst
        stepst = stepst0;
    else
        stepst = Tout-(nstepst-1)*stepst0;
    end
    yvec = hushen(stepst*(t-1)+1:stepst*(t-1)+wind);
    yvecall = [yvec;hushen(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst)];
    ynew = yvecall(end-stepst+1:end);
    parfor i = 1:N
        xvec = shenwan(stepst*(t-1)+1:stepst*(t-1)+wind,i);
        xvecall = [xvec;shenwan(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst,i)];
        xnew = xvecall(end-stepst:end);
        [VaRin,VaRforc] = VaROOS(xvec,xnew,GMODEL,THETA1); warning off
        [CoVoRforc,CoESforc,CoVaRin,CoESin,param] = CoforcOOSc(xvec,yvec,VaRin,VaRforc,THETA1,THETA2,CoMODEL,epsil,miter,xnew,ynew);
        VaRforcnc(:,i) = VaRforc;
        CoVoRforcnc(:,i) = CoVoRforc;
        CoESforcnc(:,i) = CoESforc;
        paramnc(:,i) = param;
    end
    
    if t ~= nstepst
        VaRi2sysoosc(stepst*(t-1)+1:stepst*t,:) = VaRforcnc;
        CoVaRi2sysoosc(stepst*(t-1)+1:stepst*t,:) = CoVoRforcnc;
        CoESi2sysoosc(stepst*(t-1)+1:stepst*t,:) = CoESforcnc;
    else
        VaRi2sysoosc(end-stepst+1:end,:) = VaRforcnc;
        CoVaRi2sysoosc(end-stepst+1:end,:) = CoVoRforcnc;
        CoESi2sysoosc(end-stepst+1:end,:) = CoESforcnc;
    end
    parami2sysoosc(:,:,t) = paramnc;
    [t,toc]
end


if GMODEL == 1
    if CoMODEL == 1
        VaRGauSAVoosi2sc = VaRi2sysoosc;
        CoVaRGauSAVoosi2sc = CoVaRi2sysoosc;
        CoESGauSAVoosi2sc = CoESi2sysoosc;
        paramGauSAVoosi2sc = parami2sysoosc;
    else
        VaRGauASoosi2sc = VaRi2sysoosc;
        CoVaRGauASoosi2sc = CoVaRi2sysoosc;
        CoESGauASoosi2sc = CoESi2sysoosc;
        paramGauASoosi2sc = parami2sysoosc;
    end
else  GMODEL == 2
    if CoMODEL == 1
        VaRPOTSAVoosi2sc = VaRi2sysoosc;
        CoVaRPOTSAVoosi2sc = CoVaRi2sysoosc;
        CoESPOTSAVoosi2sc = CoESi2sysoosc;
        paramPOTSAVoosi2sc = parami2sysoosc;
    else
        VaRPOTASoosi2sc = VaRi2sysoosc;
        CoVaRPOTASoosi2sc = CoVaRi2sysoosc;
        CoESPOTASoosi2sc = CoESi2sysoosc;
        paramPOTASoosi2sc = parami2sysoosc;
    end
end

save VaRGauSAVoosi2sc.mat VaRGauSAVoosi2sc
save CoVaRGauSAVoosi2sc.mat CoVaRGauSAVoosi2sc
save CoESGauSAVoosi2sc.mat CoESGauSAVoosi2sc

save VaRGauASoosi2sc.mat VaRGauASoosi2sc
save CoVaRGauASoosi2sc.mat CoVaRGauASoosi2sc
save CoESGauASoosi2sc.mat CoESGauASoosi2sc

save VaRPOTSAVoosi2sc.mat VaRPOTSAVoosi2sc
save CoVaRPOTSAVoosi2sc.mat CoVaRPOTSAVoosi2sc
save CoESPOTSAVoosi2sc.mat CoESPOTSAVoosi2sc

save VaRPOTASoosi2sc.mat VaRPOTASoosi2sc
save CoVaRPOTASoosi2sc.mat CoVaRPOTASoosi2sc
save CoESPOTASoosi2sc.mat CoESPOTASoosi2sc