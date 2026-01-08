function [paramestall,CoVaRall,CoESall,delta,condiexp] = constOOS(dataX,dataY,VaRsector,THETA1,THETA2,CoMODEL,epsil,miter)
%%%%  This is the main code to estmae all the parameters in the CoCAViaR-c model.  %%%%

% dataX--conditioning return
% dataY--interested return
% THETA1/THETA2--prabability for conditioning return and interested return (in the CoVaR and CoES)
% CoMODEL--CoSAV or CoAS (1 or 2)
% p--scale parameter to calculate the score (0,1/2 or 1)
% epsil-- difference between two iterations of EM algorithm
% miter-- max number of iteration


T = length(dataX);
%
%%  initial values
% initial values of CoCAViaR parameters
if CoMODEL == 1
    [cobetaSAV,~,CoVaRSAV] = optcavsav27(THETA2,dataX,dataY,VaRsector);  warning off
else
    [cobetaAS,~,CoVaRAS] = optcavAS28(THETA2,dataX,dataY,VaRsector);  warning off
end
if CoMODEL == 1
    CoVaR = CoVaRSAV;
else
    CoVaR = CoVaRAS;
end
% initial values for gamma
data1 = dataX;
data2 = dataY;
CoVaRS = CoVaR;
gamma20 = log(mean(data2(data2<quantile(data2,THETA2)))/quantile(data2,THETA2)-1);
deltavec2 = -THETA2*(1+exp(gamma20))*(-CoVaRS);
gamma0 = gamma20;
delta0 = deltavec2(1);
%% calculate parameters using EM algorithm
CoVaREM = CoVaR;
if CoMODEL == 1
    cobeta = cobetaSAV;
else
    cobeta = cobetaAS;
end
options = optimset( 'Display', 'off','Algorithm', 'sqp');
CoESEM = zeros(T,1);
if CoMODEL == 1
    paramest2 = zeros(5,1);
else
    paramest2 = zeros(7,1);
end
es = mean(dataY(find(dataX<-VaRsector)));
gamma00 = gamma0;
CoVEM = CoVaREM;
paraminit = [gamma00;cobeta];
paraminit0 = paraminit;
dif = 1;
Lcnew = 1;
iter = 1;
while dif>epsil && iter<miter
    Lc = Lcnew;
    [f,gamma,delta] = denALconst(dataY,gamma00,THETA2,-CoVEM,es);
    % EM (E step) for AL density
    [ut,zt] = uandz(dataY,-CoVEM,delta,THETA2);
    %%% M-step
    txlv = find(dataX<-VaRsector); %find the time that X < VaR(X)
    opt3 = @(params)ElcALGAS5const(params,dataY,zt,dataX,-VaRsector,es,THETA1,THETA2,CoMODEL,gamma0,txlv);
    [paramest2,fval,exitflag]  = fmincon(opt3,paraminit,[],[],...
        [],[],[max(paraminit-abs(paraminit)*0.3,paraminit0-abs(paraminit0)*0.3)],...
        [min(max(paraminit+abs(paraminit)*0.3,paraminit0+abs(paraminit0)*0.3),0.999*ones(length(paraminit),1))],[],options);
    [ElAL,~,~,CoVEMnew,CoESEMnew] = ElcALGAS5const(paraminit,dataY,zt,dataX,-VaRsector,es,THETA1,THETA2,CoMODEL,gamma0,txlv);
    CoVEM = -CoVEMnew;
    CoESEM = -CoESEMnew;
    Lcnew = ElAL;
    paraminit = paramest2;
    dif = abs(Lcnew-Lc);
    iter = iter+1;
end
paramestall = paramest2;
CoVaRall = CoVEM;
CoESall = CoESEM;
condiexp = es;
delta = reshape(delta,[],1);
end

