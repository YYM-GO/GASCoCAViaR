load shenwan
load hushen

% load these .mat file from Matlab_results folder
load deltaGauSAVooss2ic
load corrmatGauSAVooss2ic
load deltaGauASooss2ic
load corrmatGauASooss2ic
load deltaPOTSAVooss2ic
load corrmatPOTSAVooss2ic
load deltaPOTASooss2ic
load corrmatPOTASooss2ic
load CoVaRGauSAVooss2ic
load CoVaRGauASooss2ic
load CoVaRPOTSAVooss2ic
load CoVaRPOTASooss2ic
load condiexpGauSAVooss2ic
load condiexpGauASooss2ic
load condiexpPOTSAVooss2ic
load condiexpPOTASooss2ic


[T,N] = size(shenwan);
wind = 1400;
Tout = T-wind;
stepst0 = 100;
nstepst = ceil(Tout/stepst0);

THETA1 = 0.05; THETA2 = 0.05; %probability level for VaR and CoVaR
ksen1 = (1-THETA1*2)/(THETA1*(1-THETA1));
ksen2 = (1-THETA2*2)/(THETA2*(1-THETA2));
sig1 = sqrt(2/(THETA1*(1-THETA1)));
sig2 = sqrt(2/(THETA2*(1-THETA2)));

options = optimset( 'Display', 'off','TolFun',1e-10,'Tolcon',1e-10,'TolX',1e-10,'Algorithm', 'sqp');

GMODEL = 2      %GMODEL = 1:GARCH-Gasiaaian. GMODEL = 2:GARCH-POT
CoMODEL = 1     % CoMODEL = 1: GAS-CoSAV.  CoMODEL = 2: GAS-CoAS        (Need to run all combinations: GMODEL+CoGMODEL = 1+1/ 1+2/ 2+1/ 2+2)

if GMODEL == 1
    if CoMODEL == 1
        covarmat = CoVaRGauSAVooss2ic;
        deltai = deltaGauSAVooss2ic;
        corrmat = corrmatGauSAVooss2ic;
        condexpmat = condiexpGauSAVooss2ic;
    else
        covarmat = CoVaRGauASooss2ic;
        deltai = deltaGauASooss2ic;
        corrmat = corrmatGauASooss2ic;
        condexpmat = condiexpGauASooss2ic;
    end
else
    if CoMODEL == 1
        covarmat = CoVaRPOTSAVooss2ic;
        deltai = deltaPOTSAVooss2ic;
        corrmat = corrmatPOTSAVooss2ic;
        condexpmat = condiexpPOTSAVooss2ic;
    else
        covarmat = CoVaRPOTASooss2ic;
        deltai = deltaPOTASooss2ic;
        corrmat = corrmatPOTASooss2ic;
        condexpmat = condiexpPOTASooss2ic;
    end
end

weightori = 1/N*ones(N,1);
weightlb = -1*ones(N,1);
weightub = 1*ones(N,1);
Aeq = ones(1,N);
Beq = 1;
parfor i = 1:Tout
    ksen = ksen2*ones(N,1);
    G = sig2*eye(N);
    weightopt = @(weight)smvobj(weight,deltai(i,:),sig2,corrmat(:,:,i));
    nonlinear = @(weight)const(weight,ksen,deltai(i,:),G,THETA2,covarmat(i,:),corrmat(:,:,i),condexpmat(i,:));
    weightresult = fmincon(weightopt, weightori,[],[],Aeq,Beq,weightlb,weightub,nonlinear,options);
    weightmalall(i,:) = reshape(weightresult,1,[]);
    covarport(i) = covarmat(i,:)*weightresult;
    D = diag(deltai(i,:));
    taustar(i) = 0.5*(1-weightresult'*D*ksen/sqrt(2*weightresult'*D*G*corrmat(:,:,i)*G*D*weightresult+(weightresult'*D*ksen)^2));
    deltport(i) = delstar(weightresult,ksen,deltai(i,:),G,corrmat(:,:,i));
    condexpport(i) = condexpmat(i,:)*weightresult;
    coesport(i) = condexpport(i)-deltport(i)/taustar(i);
    HHIal(i) = weightresult'*weightresult;
    % i
end
if GMODEL == 1
    if CoMODEL == 1
        covarportGauSAVc = covarport;
        coesportGauSAVc = coesport;
        weightGauSAVc = weightmalall;
        HHIGauSAVc = HHIal;
        taustarGauSAVc = taustar;
    else
        covarportGauASc = covarport;
        coesportGauASc = coesport;
        weightGauASc = weightmalall;
        HHIGauASc = HHIal;
        taustarGauASc = taustar;
    end
else
    if CoMODEL == 1
        covarportPOTSAVc = covarport;
        coesportPOTSAVc = coesport;
        weightPOTSAVc = weightmalall;
        HHIPOTSAVc = HHIal;
        taustarPOTSAVc = taustar;
    else
        covarportPOTASc = covarport;
        coesportPOTASc = coesport;
        weightPOTASc = weightmalall;
        HHIPOTASc = HHIal;
        taustarPOTASc = taustar;
    end
end

save covarportGauSAVc.mat covarportGauSAVc
save coesportGauSAVc.mat coesportGauSAVc
save weightGauSAVc.mat weightGauSAVc
save covarportGauASc.mat covarportGauASc
save coesportGauASc.mat coesportGauASc
save weightGauASc.mat weightGauASc
save covarportPOTSAVc.mat covarportPOTSAVc
save coesportPOTSAVc.mat coesportPOTSAVc
save weightPOTSAVc.mat weightPOTSAVc
save covarportPOTASc.mat covarportPOTASc
save coesportPOTASc.mat coesportPOTASc
save weightPOTASc.mat weightPOTASc