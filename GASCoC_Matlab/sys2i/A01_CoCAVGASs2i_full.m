close
clear
clc


load shenwan       %Shenwan industry index
load hushen        %CSI 300

[T,N] = size(shenwan);
THETA1 = 0.05; THETA2 = 0.05; %probability level for VaR and CoVaR

%% calculate VaR

MODEL = 1; % GARCH-Gaussian
VaRNormals2i = VaRestimation(hushen,THETA1,MODEL);
MODEL = 2; % GARCH-POT
VaRPOTs2i = VaRestimation(hushen,THETA1,MODEL);


p = 0;      % the scale of the score, we can choose 0.
epsil = 1e-6;
TolX = 1e-6;
miter = 100;

%%        full sample
%%%%    ------------ parameter estimations -----------------------
%---------------------------------GAS-CoCAViaR ----------------------------

CoMODEL = 1;   % CoMODEL = 1: CoSAVmodel;   CoMODEL  = 2:CoAS model
tic       % GARCH-Gaussian+CoSAV
[paramfullGauSAVs2i,CoVaRfullGauSAVs2i,CoESfullGauSAVs2i] = GASGARCH1sys2inew(hushen,shenwan,-VaRNormals2i,THETA1,THETA2,CoMODEL,p,epsil,miter);
toc


CoMODEL = 2;      % GARCH-Gaussian+CoAS
tic
[paramfullGauASs2i,CoVaRfullGauASs2i,CoESfullGauASs2i] = GASGARCH1sys2inew(hushen,shenwan,-VaRNormals2i,THETA1,THETA2,CoMODEL,p,epsil,miter);
toc

CoMODEL = 1;     % GARCH-POT+CoSAV
tic
[paramfullPOTSAVs2i,CoVaRfullPOTSAVs2i,CoESfullPOTSAVs2i] = GASGARCH1sys2inew(hushen,shenwan,-VaRPOTs2i,THETA1,THETA2,CoMODEL,p,epsil,miter);
toc

CoMODEL = 2;      % GARCH-POT+CoAS
tic
[paramfullPOTASs2i,CoVaRfullPOTASs2i,CoESfullPOTASs2i,~,~,ScorePOTASs2i] = GASGARCH1sys2inew(hushen,shenwan,-VaRPOTs2i,THETA1,THETA2,CoMODEL,p,epsil,miter);
toc    



%-------------------------------------  CoCAViaR-c ------------------------
CoMODEL = 1;
tic
[paramfullGauSAVs2ic,CoVaRfullGauSAVs2ic,CoESfullGauSAVs2ic] = constsys2i(hushen,shenwan,-VaRNormals2i,THETA1,THETA2,CoMODEL,epsil,miter);
toc

CoMODEL = 2;
tic
[paramfullGauASs2ic,CoVaRfullGauASs2ic,CoESfullGauASs2ic] = constsys2i(hushen,shenwan,-VaRNormals2i,THETA1,THETA2,CoMODEL,epsil,miter);
toc

CoMODEL = 1;
tic
[paramfullPOTSAVs2ic,CoVaRfullPOTSAVs2ic,CoESfullPOTSAVs2ic] = constsys2i(hushen,shenwan,-VaRPOTs2i,THETA1,THETA2,CoMODEL,epsil,miter);
toc

CoMODEL = 2;
tic
[paramfullPOTASs2ic,CoVaRfullPOTASs2ic,CoESfullPOTASs2ic] = constsys2i(hushen,shenwan,-VaRPOTs2i,THETA1,THETA2,CoMODEL,epsil,miter);
toc



save VaRNormals2i.mat VaRNormals2i
save VaRPOTs2i.mat VaRPOTs2i

save CoVaRfullGauSAVs2i.mat CoVaRfullGauSAVs2i
save CoESfullGauSAVs2i.mat CoESfullGauSAVs2i

save CoVaRfullGauASs2i.mat CoVaRfullGauASs2i
save CoESfullGauASs2i.mat CoESfullGauASs2i

save CoVaRfullPOTSAVs2i.mat CoVaRfullPOTSAVs2i
save CoESfullPOTSAVs2i.mat CoESfullPOTSAVs2i

save CoVaRfullPOTASs2i.mat CoVaRfullPOTASs2i
save CoESfullPOTASs2i.mat CoESfullPOTASs2i
save ScorePOTASs2i.mat ScorePOTASs2i




save CoVaRfullGauSAVs2ic.mat CoVaRfullGauSAVs2ic
save CoESfullGauSAVs2ic.mat CoESfullGauSAVs2ic

save CoVaRfullGauASs2ic.mat CoVaRfullGauASs2ic
save CoESfullGauASs2ic.mat CoESfullGauASs2ic

save CoVaRfullPOTSAVs2ic.mat CoVaRfullPOTSAVs2ic
save CoESfullPOTSAVs2ic.mat CoESfullPOTSAVs2ic

save CoVaRfullPOTASs2ic.mat CoVaRfullPOTASs2ic
save CoESfullPOTASs2ic.mat CoESfullPOTASs2ic


