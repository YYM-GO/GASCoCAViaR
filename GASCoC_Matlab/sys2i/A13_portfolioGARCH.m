close
clear
clc

load shenwan
load hushen
% load these .mat files from Matlab_results folder
load Copulaconvars2iGauOOS
load Copulacorrs2iGauOOS
load Copulaconvars2iPOTOOS
load Copulacorrs2iPOTOOS
load DCCcondvars2iGauOOS
load DCCstdresidss2iGauOOS
load DCCcondvars2iPOTOOS
load DCCstdresidss2iPOTOOS

% Then go back to i2sys folder


[T,N] = size(shenwan);
wind = 1400;
Tout = T-wind;
stepst = 100;
options = optimset( 'Display', 'off','TolFun',1e-10,'Tolcon',1e-10,'TolX',1e-10,'Algorithm', 'sqp');
THETA1 = 0.05; THETA2 = 0.05; %probability level for VaR and CoVaR

weightori = 1/N*ones(N,1);
weightlb = -1*ones(N,1);
weightub = 1*ones(N,1);
Aeq = ones(1,N);
Beq = 1;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GAS copula with Gaussian innovation

concov = Copulaconvars2iGauOOS;
corrmat = Copulacorrs2iGauOOS;
parfor i = 1:Tout
    weightopt = @(weight)mvobj(weight,concov(i,:),corrmat(:,:,i));
    weightresult = fmincon(weightopt, weightori,[],[],Aeq,Beq,weightlb,weightub,[],options);
    weightcopulaGauall(i,:) = reshape(weightresult,1,[]);
    HHIcopulaGau(i) = weightresult'*weightresult;
    %     i
end
copulaCoVaRGauport = zeros(Tout,4);   copulaCoESGauport = zeros(Tout,4);
tic
parfor i = 1:Tout    % This code will take a long time
    t = ceil(i/stepst);
    xvec = hushen(stepst*(t-1)+1:stepst*(t-1)+wind);
    xvecall = [xvec;hushen(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst)];
    yvec = shenwan(stepst*(t-1)+1:stepst*(t-1)+wind,:)*weightcopulaGauall(i,:)';
    yvecall = [yvec;shenwan(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst,:)*weightcopulaGauall(i,:)'];
    xnew = xvecall(end-stepst:end);
    ynew = yvecall(end-stepst:end);
    [copulaCoVaR,copulaCoES] = copulaGauportforc(xvec,yvec,xnew,ynew,THETA1,THETA2,KK,HESSgumbel,NN,RR,HESSstudt,RRnorm,HESSnorm,KKclayton,HESSclayton,i-stepst*(t-1),i);
    copulaCoVaRGauport(i,:) = copulaCoVaR;
    copulaCoESGauport(i,:) = copulaCoES;
    i
end
toc


save weightcopulaGauall.mat weightcopulaGauall
save copulaCoVaRGauport.mat copulaCoVaRGauport
save copulaCoESGauport.mat copulaCoESGauport
%%%%%%%%%%%%%%%%%%%%%% GAS Copula with POT innovation

concov = Copulaconvars2iPOTOOS;
corrmat = Copulacorrs2iPOTOOS;
parfor i = 1:Tout
    weightopt = @(weight)mvobj(weight,concov(i,:),corrmat(:,:,i));
    weightresult = fmincon(weightopt, weightori,[],[],Aeq,Beq,weightlb,weightub,[],options);
    weightcopulaPOTall(i,:) = reshape(weightresult,1,[]);
    HHIcopulaPOT(i) = weightresult'*weightresult;
    %     i
end
copulaCoVaRPOTport = zeros(Tout,4);  copulaCoESPOTport = zeros(Tout,4);

tic
parfor i = 1:Tout
    t = ceil(i/stepst);
    xvec = hushen(stepst*(t-1)+1:stepst*(t-1)+wind);
    xvecall = [xvec;hushen(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst)];
    yvec = shenwan(stepst*(t-1)+1:stepst*(t-1)+wind,:)*weightcopulaPOTall(i,:)';
    yvecall = [yvec;shenwan(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst,:)*weightcopulaPOTall(i,:)'];
    xnew = xvecall(end-stepst:end);
    ynew = yvecall(end-stepst:end);
    [copulaCoVaR,copulaCoES] = copulaPOTportforc(xvec,yvec,xnew,ynew,THETA1,THETA2,KK,HESSgumbel,NN,RR,HESSstudt,RRnorm,HESSnorm,KKclayton,HESSclayton,i-stepst*(t-1),i);
    copulaCoVaRPOTport(i,:) = copulaCoVaR;
    copulaCoESPOTport(i,:) = copulaCoES;
    i
end
toc

save weightcopulaPOTall.mat weightcopulaPOTall
save copulaCoVaRPOTport.mat copulaCoVaRPOTport
save copulaCoESPOTport.mat copulaCoESPOTport
%%%%%%%%%%%%%%%%%%%%%% DCC with Gaussian innovation

concov = DCCcondvars2iGauOOS;
corrmat = DCCstdresidss2iGauOOS;
parfor i = 1:Tout
    weightopt = @(weight)mvobj(weight,concov(i,:),corrmat(:,:,i));
    weightresult = fmincon(weightopt, weightori,[],[],Aeq,Beq,weightlb,weightub,[],options);
    weightDCCGauall(i,:) = reshape(weightresult,1,[]);
    HHIDCCGau(i) = weightresult'*weightresult;
    %      i
end
DCCCoVaRGauport = zeros(Tout,1);  DCCCoESGauport = zeros(Tout,1);
tic
parfor i = 1:Tout
    t = ceil(i/stepst);
    xvec = hushen(stepst*(t-1)+1:stepst*(t-1)+wind);
    xvecall = [xvec;hushen(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst)];
    yvec = shenwan(stepst*(t-1)+1:stepst*(t-1)+wind,:)*weightDCCGauall(i,:)';
    yvecall = [yvec;shenwan(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst,:)*weightDCCGauall(i,:)'];
    xnew = xvecall(end-stepst:end);
    ynew = yvecall(end-stepst:end);
    [DCCCoVaR,DCCCoES] = dccGauportforc(xvec,yvec,xnew,ynew,THETA1,THETA2,i-stepst*(t-1));
    DCCCoVaRGauport(i) = DCCCoVaR;
    DCCCoESGauport(i) = DCCCoES;
    i
end
toc

save weightDCCGauall.mat weightDCCGauall
save DCCCoVaRGauport.mat DCCCoVaRGauport
save DCCCoESGauport.mat DCCCoESGauport
%%%%%%%%%%%%%%%%%%%%%% DCC with POT innovation

concov = DCCcondvars2iPOTOOS;
corrmat = DCCstdresidss2iPOTOOS;
parfor i = 1:Tout
    weightopt = @(weight)mvobj(weight,concov(i,:),corrmat(:,:,i));
    weightresult = fmincon(weightopt, weightori,[],[],Aeq,Beq,weightlb,weightub,[],options);
    weightDCCPOTall(i,:) = reshape(weightresult,1,[]);
    HHIDCCPOT(i) = weightresult'*weightresult;
    %     i
end
DCCCoVaRPOTport = zeros(Tout,1);   DCCCoESPOTport = zeros(Tout,1);
tic
parfor i = 1:Tout
    t = ceil(i/stepst);
    xvec = hushen(stepst*(t-1)+1:stepst*(t-1)+wind);
    xvecall = [xvec;hushen(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst)];
    yvec = shenwan(stepst*(t-1)+1:stepst*(t-1)+wind,:)*weightDCCPOTall(i,:)';
    yvecall = [yvec;shenwan(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst,:)*weightDCCPOTall(i,:)'];
    xnew = xvecall(end-stepst:end);
    ynew = yvecall(end-stepst:end);
    [DCCCoVaR,DCCCoES] = dccPOTportforc(xvec,yvec,xnew,ynew,THETA1,THETA2,i-stepst*(t-1));warning off
    DCCCoVaRPOTport(i) = DCCCoVaR;
    DCCCoESPOTport(i) = DCCCoES;
    i
end
toc
save weightDCCPOTall.mat weightDCCPOTall
save DCCCoVaRPOTport.mat DCCCoVaRPOTport
save DCCCoESPOTport.mat DCCCoESPOTport
