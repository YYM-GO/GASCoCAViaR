close
clear
clc
%    GARCH-DCC model the standardized residuals follows bavariate Gaussian distribution

load shenwan       % Shenwan industry index
load hushen        % CSI 300
[T,N] = size(shenwan);
THETA1 = 0.05; THETA2 = 0.05; %probability level for VaR and CoVaR
wind = 1400;
Tout = T-wind;
stepst0 = 100;
nstepst = ceil(Tout/stepst0);

DCCVaRi2sPOTOOS = nan(Tout,N);
DCCCoVaRi2sPOTOOS = nan(Tout,N);
DCCCoESi2sPOTOOS = nan(Tout,N);

for t = 1:nstepst
    tic
    if t ~= nstepst
        stepst = stepst0;
    else
        stepst = Tout-(nstepst-1)*stepst0;
    end
    parfor i = 1:N
        xvec = shenwan(stepst*(t-1)+1:stepst*(t-1)+wind,i);
        xvecall = [xvec;shenwan(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst,i)];
        yvec = hushen(stepst*(t-1)+1:stepst*(t-1)+wind);
        yvecall = [yvec;hushen(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst)];
        xnew = xvecall(end-stepst:end);
        ynew = yvecall(end-stepst:end);
        [DCCVaR,DCCCoVaR,DCCCoES,parameters] = dccPOTforc(xvec,yvec,xnew,ynew,THETA1,THETA2,stepst);
        DCCVaRall(:,i) = DCCVaR;
        DCCCoVaRall(:,i) = DCCCoVaR;
        DCCCoESall(:,i) = DCCCoES;
    end
    if t ~= nstepst
        DCCVaRi2sPOTOOS(stepst*(t-1)+1:stepst*t,:) = DCCVaRall;
        DCCCoVaRi2sPOTOOS(stepst*(t-1)+1:stepst*t,:) = DCCCoVaRall;
        DCCCoESi2sPOTOOS(stepst*(t-1)+1:stepst*t,:) = DCCCoESall;
    else
        DCCVaRi2sPOTOOS(end-stepst+1:end,:) = DCCVaRall;
        DCCCoVaRi2sPOTOOS(end-stepst+1:end,:) = DCCCoVaRall;
        DCCCoESi2sPOTOOS(end-stepst+1:end,:) = DCCCoESall;
    end
    [t,toc]
end


save DCCVaRi2sPOTOOS.mat DCCVaRi2sPOTOOS
save DCCCoVaRi2sPOTOOS.mat DCCCoVaRi2sPOTOOS
save DCCCoESi2sPOTOOS.mat DCCCoESi2sPOTOOS