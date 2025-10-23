close
clear
clc


load shenwan       % Shenwan industry index
load hushen        % CSI 300
[T,N] = size(shenwan);
THETA1 = 0.05; THETA2 = 0.05; %probability level for VaR and CoVaR

DCCVaRi2sPOTfull = nan(T,N);
DCCESi2sPOTfull = nan(T,N);
DCCCoVaRi2sPOTfull = nan(T,N);
DCCCoESi2sPOTfull = nan(T,N);
parfor i = 1:N
    tic
    esgrid = nan(5,1);
    xvec = shenwan(:,i);
    yvec = hushen;
    resids = [xvec,yvec];
    MUhat = [0,0];
    [parameters, ~ ,Ht,~,Rt] = dccPOT(resids,[],1,0,1);
    stdresids = resids./sqrt([squeeze(Ht(1,1,:)),squeeze(Ht(2,2,:))]);
    for t = 1:T
        DCCVaRi2sPOTfull(t,i) = sqrt(Ht(1,1,t))*POTinv(parameters(6),parameters(4),parameters(5),THETA1);
        DCCESi2sPOTfull(t,i) = (DCCVaRi2sPOTfull(t,i)-parameters(6)*parameters(4)-parameters(5))/(1-parameters(4));
        DCCCoVaRi2sPOTfull(t,i) = CoVaRfrocNormalPOT(THETA1,THETA2,Rt(1,2,t),diag(Ht(:,:,t)),stdresids(:,2),MUhat,parameters(12),parameters(10),parameters(11));
        grid = [0.005:0.01:THETA2];
        for l = 1:length(grid)
            esgrid(l) = CoVaRfrocNormalPOT(THETA1,grid(l),Rt(1,2,t),diag(Ht(:,:,t)),stdresids(:,2),MUhat,parameters(12),parameters(10),parameters(11)); warning off
        end
        DCCCoESi2sPOTfull(t,i) = mean(esgrid);
    end
    [i,toc]
end
save DCCCoESi2sPOTfull.mat DCCCoESi2sPOTfull
save DCCCoVaRi2sPOTfull.mat DCCCoVaRi2sPOTfull
save DCCVaRi2sPOTfull.mat DCCVaRi2sPOTfull
