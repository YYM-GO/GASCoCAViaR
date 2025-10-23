close
clear
clc


load shenwan       % Shenwan industry index
load hushen        % CSI 300
[T,N] = size(shenwan);
THETA1 = 0.05; THETA2 = 0.05; %probability level for VaR and CoVaR

DCCVaRi2sGaufull = nan(T,N);
DCCESi2sGaufull = nan(T,N);
DCCCoVaRi2sGaufull = nan(T,N);
DCCCoESi2sGaufull = nan(T,N);
parfor i = 1:N
    tic
    
    esgrid = nan(5,1);
    xvec = shenwan(:,i);
    yvec = hushen;
    resids = [xvec,yvec];
    MUhat = [0,0]; % only GARCH, No ARMA
    [parameters, ~ ,Ht,~,Rt] = dcc(resids,[],1,0,1);
    stdresids = resids./sqrt([squeeze(Ht(1,1,:)),squeeze(Ht(2,2,:))]);
    for t = 1:T
        DCCVaRi2sGaufull(t,i) = sqrt(Ht(1,1,t))*norminv(THETA1);
        DCCESi2sGaufull(t,i) = - sqrt(Ht(1,1,t))*normpdf(norminv(THETA1))/(THETA1);
        DCCCoVaRi2sGaufull(t,i) = CoVaRfrocNormalN(THETA1,THETA2,Rt(1,2,t),diag(Ht(:,:,t)),stdresids(:,2),MUhat);
        grid = [0.005:0.01:THETA2];
        for l = 1:length(grid)
            esgrid(l) = CoVaRfrocNormalN(THETA1,grid(l),Rt(1,2,t),diag(Ht(:,:,t)),stdresids(:,2),MUhat); warning off
        end
        DCCCoESi2sGaufull(t,i) = mean(esgrid);
    end
    i
end
save DCCCoESi2sGaufull.mat DCCCoESi2sGaufull
save DCCCoVaRi2sGaufull.mat DCCCoVaRi2sGaufull
save DCCVaRi2sGaufull.mat DCCVaRi2sGaufull
