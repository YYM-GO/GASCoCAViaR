close
clear
clc


load shenwan       %Shenwan industry index
load hushen        %CSI 300
[T,N] = size(shenwan);


THETA1 = 0.05; THETA2 = 0.05; %probability level for VaR and CoVaR

DCCVaRs2iPOTfull = nan(T,N);
DCCESs2iPOTfull = nan(T,N);
DCCCoVaRs2iPOTfull = nan(T,N);
DCCCoESs2iPOTfull = nan(T,N);
parfor i = 1:N
    tic
    esgrid = nan(5,1);
    yvec = shenwan(:,i);
    xvec = hushen;
    resids = [xvec,yvec];
    MUhat = [0,0];
    [parameters, ~ ,Ht,~,Rt] = dccPOT(resids,[],1,0,1);
    stdresids = resids./sqrt([squeeze(Ht(1,1,:)),squeeze(Ht(2,2,:))]);
    for t = 1:T
        DCCVaRs2iPOTfull(t,i) = sqrt(Ht(1,1,t))*POTinv(parameters(6),parameters(4),parameters(5),THETA1);
        DCCESs2iPOTfull(t,i) = (DCCVaRs2iPOTfull(t,i)-parameters(6)*parameters(4)-parameters(5))/(1-parameters(4));
        DCCCoVaRs2iPOTfull(t,i) = CoVaRfrocNormalPOT(THETA1,THETA2,Rt(1,2,t),diag(Ht(:,:,t)),stdresids(:,2),MUhat,parameters(12),parameters(10),parameters(11));
        grid = [0.005:0.01:THETA2];
        for l = 1:length(grid)
            esgrid(l) = CoVaRfrocNormalPOT(THETA1,grid(l),Rt(1,2,t),diag(Ht(:,:,t)),stdresids(:,2),MUhat,parameters(12),parameters(10),parameters(11)); warning off
        end
        DCCCoESs2iPOTfull(t,i) = mean(esgrid);
    end
    [i,toc]
end

save DCCVaRs2iPOTfull.mat DCCVaRs2iPOTfull
save DCCCoVaRs2iPOTfull.mat DCCCoVaRs2iPOTfull
save DCCCoESs2iPOTfull.mat DCCCoESs2iPOTfull