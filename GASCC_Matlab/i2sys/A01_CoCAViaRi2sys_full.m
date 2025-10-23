close
clear
clc


load shenwan       % Shenwan industry index
load hushen        % CSI 300
[T,N] = size(shenwan);
THETA1 = 0.05; THETA2 = 0.05; %probability level for VaR and CoVaR

%% calculate VaR
MODEL = 1;    % GARCH-Gaussian for Xt
VaRNormali2s = VaRestimation(shenwan,THETA1,MODEL);
MODEL = 2;    % GARCH-POT for Xt
VaRPOTi2s = VaRestimation(shenwan,THETA1,MODEL);

%
p = 0;      % the scale of the score, we choose 0.
epsil = 1e-6;
TolX = 1e-6;
miter = 100;

%%        full sample

%%-------------------------------------------- GAS-CoCAViaR-----------each about 150-200 seconds
CoMODEL = 1;      % CoMODEL = 1: CoSAVmodel;   CoMODEL  = 2:CoAS model
tic           % GARCH-Gaussian+CoSAV
[paramfullGauSAVi2s,CoVaRfullGauSAVi2s,CoESfullGauSAVi2s] = GASGARCH1i2sysnew(shenwan,hushen,-VaRNormali2s,THETA1,THETA2,CoMODEL,p,epsil,miter);
toc


CoMODEL = 2;      % GARCH-Gaussian+CoAS
tic
[paramfullGauASi2s,CoVaRfullGauASi2s,CoESfullGauASi2s] = GASGARCH1i2sysnew(shenwan,hushen,-VaRNormali2s,THETA1,THETA2,CoMODEL,p,epsil,miter);
toc

CoMODEL = 1;       % GARCH-POT+CoSAV
tic
[paramfullPOTSAVi2s,CoVaRfullPOTSAVi2s,CoESfullPOTSAVi2s] = GASGARCH1i2sysnew(shenwan,hushen,-VaRPOTi2s,THETA1,THETA2,CoMODEL,p,epsil,miter);
toc


CoMODEL = 2;    % GARCH-POT+CoAS
tic
[paramfullPOTASi2s,CoVaRfullPOTASi2s,CoESfullPOTASi2s,gammafullPOTASi2s,deltafullPOTASi2s,ScorePOTASi2s] = GASGARCH1i2sysnew(shenwan,hushen,-VaRPOTi2s,THETA1,THETA2,CoMODEL,p,epsil,miter);
toc
%%% standard error
parfor i = 1:N
opt3 = @(params)GASsys2iLL(params,shenwan(:,i),hushen,-VaRPOTi2s(:,i),THETA1,THETA2,CoMODEL,p);
[~,se_S] = Compute_se_errors_ML(opt3,paramfullPOTASi2s(:,i),T);
se_SPOTi2sAS(:,i) = se_S;
end

paramfullPOTASi2sresult = paramfullPOTASi2s;
paramfullPOTASi2sresult(1,:) = paramfullPOTASi2s(1,:)*10;
paramfullPOTASi2sresult(2,:) = paramfullPOTASi2s(2,:)*10;
paramfullPOTASi2sresult(4,:) = paramfullPOTASi2s(4,:)*10;
paramfullPOTASi2sresult(5,:) = paramfullPOTASi2s(5,:)*10;
paramfullPOTASi2sresult(6,:) = paramfullPOTASi2s(6,:)*10;
paramfullPOTASi2sresult(7,:) = paramfullPOTASi2s(7,:)*10;
paramfullPOTASi2sresult(8,:) = paramfullPOTASi2s(8,:)*10;
se_SPOTi2sASresult = se_SPOTi2sAS;
se_SPOTi2sASresult(1,:) = se_SPOTi2sAS(1,:)*10;
se_SPOTi2sASresult(2,:) = se_SPOTi2sAS(2,:)*10;
se_SPOTi2sASresult(4,:) = se_SPOTi2sAS(4,:)*10;
se_SPOTi2sASresult(5,:) = se_SPOTi2sAS(5,:)*10;
se_SPOTi2sASresult(6,:) = se_SPOTi2sAS(6,:)*10;
se_SPOTi2sASresult(7,:) = se_SPOTi2sAS(7,:)*10;
se_SPOTi2sASresult(8,:) = se_SPOTi2sAS(8,:)*10;
paramfullPOTASi2sresult'
se_SPOTi2sASresult'

%--------------------------------------- CoCAViaR-c -----------------------


CoMODEL = 1;
tic
[paramfullGauSAVi2sc,CoVaRfullGauSAVi2sc,CoESfullGauSAVi2sc] = consti2sys(shenwan,hushen,-VaRNormali2s,THETA1,THETA2,CoMODEL,epsil,miter);
toc

CoMODEL = 2;
tic
[paramfullGauASi2sc,CoVaRfullGauASi2sc,CoESfullGauASi2sc] = consti2sys(shenwan,hushen,-VaRNormali2s,THETA1,THETA2,CoMODEL,epsil,miter);
toc

CoMODEL = 1;
tic
[paramfullPOTSAVi2sc,CoVaRfullPOTSAVi2sc,CoESfullPOTSAVi2sc] = consti2sys(shenwan,hushen,-VaRPOTi2s,THETA1,THETA2,CoMODEL,epsil,miter);

CoMODEL = 2;
tic
[paramfullPOTASi2sc,CoVaRfullPOTASi2sc,CoESfullPOTASi2sc] = consti2sys(shenwan,hushen,-VaRPOTi2s,THETA1,THETA2,CoMODEL,epsil,miter);
toc

save VaRNormali2s.mat VaRNormali2s
save VaRPOTi2s.mat VaRPOTi2s


save CoVaRfullGauSAVi2s.mat CoVaRfullGauSAVi2s
save CoESfullGauSAVi2s.mat CoESfullGauSAVi2s

save CoVaRfullGauASi2s.mat CoVaRfullGauASi2s
save CoESfullGauASi2s.mat CoESfullGauASi2s

save CoVaRfullPOTSAVi2s.mat CoVaRfullPOTSAVi2s
save CoESfullPOTSAVi2s.mat CoESfullPOTSAVi2s

save paramfullPOTASi2s.mat paramfullPOTASi2s
save CoVaRfullPOTASi2s.mat CoVaRfullPOTASi2s
save CoESfullPOTASi2s.mat CoESfullPOTASi2s
save ScorePOTASi2s.mat ScorePOTASi2s
save gammafullPOTASi2s.mat gammafullPOTASi2s
save paramfullPOTASi2sresult.mat paramfullPOTASi2sresult
save se_SPOTi2sASresult.mat se_SPOTi2sASresult


save CoVaRfullGauSAVi2sc.mat CoVaRfullGauSAVi2sc
save CoESfullGauSAVi2sc.mat CoESfullGauSAVi2sc

save CoVaRfullGauASi2sc.mat CoVaRfullGauASi2sc
save CoESfullGauASi2sc.mat CoESfullGauASi2sc

save CoVaRfullPOTSAVi2sc.mat CoVaRfullPOTSAVi2sc
save CoESfullPOTSAVi2sc.mat CoESfullPOTSAVi2sc

save CoVaRfullPOTASi2sc.mat CoVaRfullPOTASi2sc
save CoESfullPOTASi2sc.mat CoESfullPOTASi2sc

