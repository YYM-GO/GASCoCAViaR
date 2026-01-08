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

% copula models
tic
reps = 10000;
KK = [1.0001;(1.1:0.1:3)';(3.25:0.25:5)';(6:30)'];% values of kappa at which to evaluate the function
HESSgumbel = nan(length(KK),1);
parfor ii=1:length(KK)
    Usim = Gumbel_rnd(KK(ii),reps,1000*100*ii);  warning('off')
    Usim = 1-Usim;  % rotating the data
    scoreSIM = LLgrad_1('gumbelCLa',KK(ii),Usim);
    HESSgumbel(ii,1) = mean(scoreSIM.^2);
end
toc
tic
reps = 1000;
NN = [2.1;3;4;6;10;20;100];
RR = [-0.99;-0.75;-0.5;-0.25;-0.1;0;0.1;0.25;0.5;0.75;0.8;0.85;0.9;0.95;0.99];
HESSstudt = nan(length(RR),length(NN));
rng(122,'philox')
for rr=1:length(RR)
    for nn=1:length(NN)
        temp = mvtrnd([1,RR(rr);RR(rr),1],NN(nn),reps);
        temp = [temp;-temp];  % imposing symmetry
        temp = [temp;[temp(:,2),temp(:,1)]];  % imposing exchangeability
        U = tdis_cdf(temp(:,1),NN(nn));
        V = tdis_cdf(temp(:,2),NN(nn));
        scoreSIM = LLgrad_1('tcopulaCLa',[RR(rr);NN(nn)],[U,V]);
        HESSstudt(rr,nn) = mean(scoreSIM(:,1).^2);
    end
end
toc
tic
RRnorm = [-0.99:0.01:0.99999]';
HESSnorm = nan(length(RRnorm),1);
rng(1334,'philox')
for rr=1:length(RRnorm)
    temp = mvnrnd([0,0],[1,RRnorm(rr);RRnorm(rr),1],reps);
    temp = [temp;-temp];  % imposing symmetry
    temp = [temp;[temp(:,2),temp(:,1)]];  % imposing exchangeability
    U = normcdf(temp(:,1));
    V = normcdf(temp(:,2));
    scoreSIMnorm = LLgrad_1('NormalCopula_CLa',RRnorm(rr),[U,V]);
    HESSnorm(rr) = mean(scoreSIMnorm.^2);
end
toc
reps = 10000;
KKclayton = [-0.5:0.03:20]';
HESSclayton = nan(length(KKclayton),1);
parfor ii=1:length(KKclayton)
    Usim = clayton_rnd(KKclayton(ii),reps,1667+100*ii);  warning('off')
    scoreSIM = LLgrad_1('claytonCLa',KKclayton(ii),Usim);
    HESSclayton(ii) = mean(scoreSIM.^2);
end
copulaVaRall = []; copulaCoVaRall = [];copulaCoESall = []; convarall = []; stdresidsall = [];
CopulaVaRs2iGauOOS = nan(Tout,N,4);
CopulaCoVaRs2iGauOOS = nan(Tout,N,4);
CopulaCoESs2iGauOOS = nan(Tout,N,4);
Copulaconvars2iGauOOS = nan(Tout,N);
Copulacorrs2iGauOOS = nan(N,N,Tout);
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
        [copulaVaR,copulaCoVaR,copulaCoES,condvar2out,stdresidsIS2] = copulaGauforc(xvec,yvec,xnew,ynew,THETA1,THETA2,KK,HESSgumbel,NN,RR,HESSstudt,RRnorm,HESSnorm,KKclayton,HESSclayton,t,i);
        copulaVaRall(:,i,:) = copulaVaR;
        copulaCoVaRall(:,i,:) = copulaCoVaR;
        copulaCoESall(:,i,:) = copulaCoES;
        convarall(:,i) = condvar2out;
        stdresidsall(:,i) = stdresidsIS2;
    end
    if t ~= nstepst
        CopulaVaRs2iGauOOS(stepst*(t-1)+1:stepst*t,:,:) = copulaVaRall;
        CopulaCoVaRs2iGauOOS(stepst*(t-1)+1:stepst*t,:,:) = copulaCoVaRall;
        CopulaCoESs2iGauOOS(stepst*(t-1)+1:stepst*t,:,:) = copulaCoESall;
        Copulaconvars2iGauOOS(stepst*(t-1)+1:stepst*t,:) = convarall;
        Copulacorrs2iGauOOS(:,:,stepst*(t-1)+1:stepst*t) = repmat(corr(stdresidsall),1,1,stepst);
    else
        CopulaVaRs2iGauOOS(end-stepst+1:end,:,:) = copulaVaRall;
        CopulaCoVaRs2iGauOOS(end-stepst+1:end,:,:) = copulaCoVaRall;
        CopulaCoESs2iGauOOS(end-stepst+1:end,:,:) = copulaCoESall;
        Copulaconvars2iGauOOS(end-stepst+1:end,:) = convarall;
        Copulacorrs2iGauOOS(:,:,end-stepst+1:end) = repmat(corr(stdresidsall),1,1,stepst);
    end
    [t,toc]
end


save CopulaVaRs2iGauOOS.mat CopulaVaRs2iGauOOS
save CopulaCoVaRs2iGauOOS.mat CopulaCoVaRs2iGauOOS
save CopulaCoESs2iGauOOS.mat CopulaCoESs2iGauOOS
save Copulaconvars2iGauOOS.mat Copulaconvars2iGauOOS
save Copulacorrs2iGauOOS.mat Copulacorrs2iGauOOS
