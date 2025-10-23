close
clear
clc

load shenwan       %Shenwan industry index
load hushen        %CSI 300
[T,N] = size(shenwan);


THETA1 = 0.05; THETA2 = 0.05; %probability level for VaR and CoVaR
DCCVaRs2iGaufull = nan(T,N);
DCCESs2iGaufull = nan(T,N);
DCCCoVaRs2iGaufull = nan(T,N);
DCCCoESs2iGaufull = nan(T,N);
parfor i = 1:N
    tic
    esgrid = nan(5,1);
    yvec = shenwan(:,i);
    xvec = hushen;
    resids = [xvec,yvec];
    MUhat = [0,0]; % only GARCH, No ARMA
    [parameters, ~ ,Ht,~,Rt] = dcc(resids,[],1,0,1);
    stdresids = resids./sqrt([squeeze(Ht(1,1,:)),squeeze(Ht(2,2,:))]);
    for t = 1:T
        DCCVaRs2iGaufull(t,i) = sqrt(Ht(1,1,t))*norminv(THETA1);
        DCCESs2iGaufull(t,i) = - sqrt(Ht(1,1,t))*normpdf(norminv(THETA1))/(THETA1);
        DCCCoVaRs2iGaufull(t,i) = CoVaRfrocNormalN(THETA1,THETA2,Rt(1,2,t),diag(Ht(:,:,t)),stdresids(:,2),MUhat);
        grid = [0.005:0.01:THETA2];
        for l = 1:length(grid)
            esgrid(l) = CoVaRfrocNormalN(THETA1,grid(l),Rt(1,2,t),diag(Ht(:,:,t)),stdresids(:,2),MUhat); warning off
        end
        DCCCoESs2iGaufull(t,i) = mean(esgrid);
    end
    i
end

save DCCVaRs2iGaufull.mat DCCVaRs2iGaufull
save DCCCoVaRs2iGaufull.mat DCCCoVaRs2iGaufull
save DCCCoESs2iGaufull.mat DCCCoESs2iGaufull