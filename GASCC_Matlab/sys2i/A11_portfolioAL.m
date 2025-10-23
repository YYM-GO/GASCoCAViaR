load shenwan
load hushen

% load these .mat file from Matlab_results folder
load deltaGauSAVooss2i
load corrmatGauSAVooss2i
load deltaGauASooss2i
load corrmatGauASooss2i
load deltaPOTSAVooss2i
load corrmatPOTSAVooss2i
load deltaPOTASooss2i
load corrmatPOTASooss2i
load CoVaRGauSAVooss2i
load CoVaRGauASooss2i
load CoVaRPOTSAVooss2i
load CoVaRPOTASooss2i
load condiexpGauSAVooss2i
load condiexpGauASooss2i
load condiexpPOTSAVooss2i
load condiexpPOTASooss2i



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
CoMODEL = 1     % CoMODEL = 1: GAS-CoSAV.  CoMODEL = 2: GAS-CoAS       (Need to run all combinations: GMODEL+CoGMODEL = 1+1/ 1+2/ 2+1/ 2+2)

if GMODEL == 1
    if CoMODEL == 1
        covarmat = CoVaRGauSAVooss2i;
        deltai = deltaGauSAVooss2i;
        corrmat = corrmatGauSAVooss2i;
        condexpmat = condiexpGauSAVooss2i;
    else
        covarmat = CoVaRGauASooss2i;
        deltai = deltaGauASooss2i;
        corrmat = corrmatGauASooss2i;
        condexpmat = condiexpGauASooss2i;
    end
else
    if CoMODEL == 1
        covarmat = CoVaRPOTSAVooss2i;
        deltai = deltaPOTSAVooss2i;
        corrmat = corrmatPOTSAVooss2i;
        condexpmat = condiexpPOTSAVooss2i;
    else
        covarmat = CoVaRPOTASooss2i;
        deltai = deltaPOTASooss2i;
        corrmat = corrmatPOTASooss2i;
        condexpmat = condiexpPOTASooss2i;
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
        covarportGauSAV = covarport;
        coesportGauSAV = coesport;
        weightGauSAV = weightmalall;
        HHIGauSAV = HHIal;
        taustarGauSAV = taustar;
    else
        covarportGauAS = covarport;
        coesportGauAS = coesport;
        weightGauAS = weightmalall;
        HHIGauAS = HHIal;
        taustarGauAS = taustar;
    end
else
    if CoMODEL == 1
        covarportPOTSAV = covarport;
        coesportPOTSAV = coesport;
        weightPOTSAV = weightmalall;
        HHIPOTSAV = HHIal;
        taustarPOTSAV = taustar;
    else
        covarportPOTAS = covarport;
        coesportPOTAS = coesport;
        weightPOTAS = weightmalall;
        HHIPOTAS = HHIal;
        taustarPOTAS = taustar;
    end
end

save covarportGauSAV.mat covarportGauSAV
save coesportGauSAV.mat coesportGauSAV
save weightGauSAV.mat weightGauSAV
save covarportGauAS.mat covarportGauAS
save coesportGauAS.mat coesportGauAS
save weightGauAS.mat weightGauAS
save covarportPOTSAV.mat covarportPOTSAV
save coesportPOTSAV.mat coesportPOTSAV
save weightPOTSAV.mat weightPOTSAV
save covarportPOTAS.mat covarportPOTAS
save coesportPOTAS.mat coesportPOTAS
save weightPOTAS.mat weightPOTAS

