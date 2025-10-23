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

DCCVaRs2iPOTOOS = nan(Tout,N);
DCCCoVaRs2iPOTOOS = nan(Tout,N);
DCCCoESs2iPOTOOS = nan(Tout,N);
DCCcondvars2iPOTOOS = nan(Tout,N);
DCCstdresidss2iPOTOOS = nan(N,N,Tout);

for t = 1:nstepst
    tic
    if t ~= nstepst
        stepst = stepst0;
    else
        stepst = Tout-(nstepst-1)*stepst0;
    end
    parfor i = 1:N
        yvec = shenwan(stepst*(t-1)+1:stepst*(t-1)+wind,i);
        yvecall = [yvec;shenwan(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst,i)];
        xvec = hushen(stepst*(t-1)+1:stepst*(t-1)+wind);
        xvecall = [xvec;hushen(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst)];
        xnew = xvecall(end-stepst:end);
        ynew = yvecall(end-stepst:end);
        [DCCVaR,DCCCoVaR,DCCCoES,condvar2out,stdresidsIS2] = dccPOTforc(xvec,yvec,xnew,ynew,THETA1,THETA2,stepst);
        DCCVaRall(:,i) = DCCVaR;
        DCCCoVaRall(:,i) = DCCCoVaR;
        DCCCoESall(:,i) = DCCCoES;
        DCCcondvarall(:,i) = condvar2out;
        DCCstdresidsall(:,i) = stdresidsIS2;
    end
    if t ~= nstepst
        DCCVaRs2iPOTOOS(stepst*(t-1)+1:stepst*t,:) = DCCVaRall;
        DCCCoVaRs2iPOTOOS(stepst*(t-1)+1:stepst*t,:) = DCCCoVaRall;
        DCCCoESs2iPOTOOS(stepst*(t-1)+1:stepst*t,:) = DCCCoESall;
        DCCcondvars2iPOTOOS(stepst*(t-1)+1:stepst*t,:) = DCCcondvarall;
        DCCstdresidss2iPOTOOS(:,:,stepst*(t-1)+1:stepst*t) = repmat(corr(DCCstdresidsall),1,1,stepst);
    else
        DCCVaRs2iPOTOOS(end-stepst+1:end,:) = DCCVaRall;
        DCCCoVaRs2iPOTOOS(end-stepst+1:end,:) = DCCCoVaRall;
        DCCCoESs2iPOTOOS(end-stepst+1:end,:) = DCCCoESall;
        DCCcondvars2iPOTOOS(end-stepst+1:end,:) = DCCcondvarall;
        DCCstdresidss2iPOTOOS(:,:,end-stepst+1:end) = repmat(corr(DCCstdresidsall),1,1,stepst);
    end
    [t,toc]
end


save DCCVaRs2iPOTOOS.mat DCCVaRs2iPOTOOS
save DCCCoVaRs2iPOTOOS.mat DCCCoVaRs2iPOTOOS
save DCCCoESs2iPOTOOS.mat DCCCoESs2iPOTOOS
save DCCcondvars2iPOTOOS.mat DCCcondvars2iPOTOOS
save DCCstdresidss2iPOTOOS.mat DCCstdresidss2iPOTOOS

