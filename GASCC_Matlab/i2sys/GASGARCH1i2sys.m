function [paramestall,CoVaRall,CoESall,gammaall,deltaall,Scoreall] = GASGARCH1i2sys(dataX,dataY,VaRsector,THETA1,THETA2,CoMODEL,p,epsil,miter)
%%%%  This is the main code to estmae all the parameters in the  GAS-CoCAViaR model.  %%%%

% dataX--conditioning return
% dataY--interested return
% THETA1/THETA2--prabability for conditioning return and interested return (in the CoVaR and CoES)
% MODEL--SAV or AS (1 or 2)
% CoMODEL--CoSAV or CoAS (1 or 2)
% p--scale parameter to calculate the score (0,1/2 or 1)
% epsil-- difference between two iterations of EM algorithm
%


[T,N] = size(dataX);

%%  initial values
% initial values of CoCAViaR parameters

parfor i = 1:N
    if CoMODEL == 1
        [cobetaSAV(:,i),~,CoVaRSAV(:,i)] = optcavsav27(THETA2,dataX(:,i),dataY,VaRsector(:,i));
    else
        [cobetaAS(:,i),~,CoVaRAS(:,i)] = optcavAS28(THETA2,dataX(:,i),dataY,VaRsector(:,i));
    end
end

if CoMODEL == 1
    CoVaR = CoVaRSAV;
else
    CoVaR = CoVaRAS;
end



% initial values for GAS parameter 

gamma0 = zeros(1,N);
delta0 = zeros(1,N);
options = optimset( 'Display', 'off','Algorithm', 'sqp');

parfor i = 1:N
    data1 = dataX(:,i);
    data2 = dataY;
    VaRS = VaRsector(:,i);
    CoVaRS = CoVaR(:,i);
     gamma20 = log(mean(data2(find(data2<quantile(data2(data1<quantile(data1,THETA1)),THETA2))))/quantile(data2(data1<quantile(data1,THETA1)),THETA2)-1);
     deltavec2 = -THETA2*(1+exp(gamma20))*(-CoVaRS);
    [u2,z2] = uandz(data2,-CoVaRS,deltavec2,THETA2);
    try
        rng(100*i,'philox')
        opt2 = @(params)Elc20scale(params,data1,data2,-VaRS,-CoVaRS,THETA2,gamma20,z2,deltavec2(1),0);
        [paramest2,fval,exitflag]  = fmincon(opt2,[0.02;0.02;0.95]-rand(3,1)*0.02,[0,1,1],1,[],[],[-0.3;-0.3;0.8],[0.3;0.3;1.1],[],options);
    catch
        rng(200*i,'philox')
        opt2 = @(params)Elc20scaleconst(params,data1,data2,-VaRS,-CoVaRS,THETA2,gamma20,z2,deltavec2(1),0);
        [paramest2,fval,exitflag]  = fmincon(opt2,[0.02;0.02;0.95]-rand(3,1)*0.02,[0,1,1],1,[],[],[-0.3;-0.3;0.8],[0.5;0.5;1],[],options);
    end
    gasSpara(:,i) = paramest2;
    gamma0(:,i) = gamma20;
    delta0(:,i) = deltavec2(1);
end

%% calculate parameters using EM algorithm


CoVaREM = CoVaR;
if CoMODEL == 1
    cobeta = cobetaSAV;
    A = [0,1,1,zeros(1,4)];
else
    cobeta = cobetaAS;
    A = [0,1,1,zeros(1,6)];
end

parfor i = 1:N
    ut = zeros(1,T);
    zt = zeros(1,T);
    gamma = zeros(1,T+1);
    delta = zeros(1,T);
    CoESEM = zeros(T,1);
    scoreal = zeros(T,1);
    if CoMODEL == 1
        paramest2 = zeros(7,1);
    else
        paramest2 = zeros(9,1);
    end
     es = mean(dataY(find(dataX(:,i)<-VaRsector(:,i))));
    gamma00 = gamma0(:,i);
    delta00 = delta0(:,i);
    CoVEM = CoVaREM(:,i);
    paraminit = [gasSpara(:,i);cobeta(:,i)];
    paraminit0 = paraminit;
    dif = 1;
    Lcnew = 1;
    iter = 1;
    
    
    while dif>epsil && iter<miter
        Lc = Lcnew;
         [f,gamma,delta] = denAL5(dataY,paraminit(1:3),gamma00,delta00,THETA2,-CoVEM,es,p);
%         [f,gamma,delta] = denAL2(dataY,paraminit(1:3),gamma00,delta00,THETA2,-CoVEM,es,p);
        if isnan(sum(sum(gamma))) | isnan(sum(sum(f)))== 1
            [f,gamma,delta] = denAL2(dataY,paraminit(1:3),gamma00,delta00,THETA2,-CoVEM,es,p);
        end
        % EM (E step) for AL density
        
        [ut,zt] = uandz(dataY,-CoVEM,delta,THETA2);
        
        %%% M-step
        txlv = find(dataX(:,i)<-VaRsector(:,i)); %find the time that X < VaR(X)
        
        opt3 = @(params)ElcALGAS5(params,dataY,zt,dataX(:,i),-VaRsector(:,i),es,THETA1,THETA2,p,CoMODEL,gamma0(:,i),delta0(:,i),txlv);
        [paramest2,fval,exitflag]  = fmincon(opt3,paraminit,A,1,...
                          [],[],[-0.2;-0.2;0.85;max(paraminit(4:end)-abs(paraminit(4:end))*0.3,paraminit0(4:end)-abs(paraminit0(4:end))*0.3)],...
                          [0.3;0.3;0.995;min(max(paraminit(4:end)+abs(paraminit(4:end))*0.3,paraminit0(4:end)+abs(paraminit0(4:end))*0.3),0.995*ones(length(paraminit(4:end)),1))],[],options);

        [ElAL,~,~,CoVEMnew,CoESEMnew,scoreal] = ElcALGAS5(paraminit,dataY,zt,dataX(:,i),-VaRsector(:,i),es,THETA1,THETA2,p,CoMODEL,gamma0(:,i),delta0(:,i),txlv);
        CoVEM = -CoVEMnew;
        CoESEM = -CoESEMnew;
        Lcnew = ElAL;
        paraminit = paramest2
        dif = abs(Lcnew-Lc)
        iter = iter+1
    end
    paramestall(:,i) = paramest2;
    gammaall(:,i) = gamma';
    deltaall(:,i) = delta';
    CoVaRall(:,i) = CoVEM;
    CoESall(:,i) = CoESEM;
    Scoreall(:,i) = scoreal;
%      i
end

end


