function [copulaVaR,copulaCoVaR,copulaCoES,condvar2out,stdresidsIS2] = copulaGauforc(xvec,yvec,xnew,ynew,THETA1,THETA2,KK,HESSgumbel,NN,RR,HESSstudt,RRnorm,HESSnorm,KKclayton,HESSclayton,tseed,iseed)
%%calculate the VaR, CoVaR and CoES forcasts using GAS-Copula-GARCH model with Gaussian innovations.

rates = [xvec,yvec];
ratesnew = [xnew,ynew];
wind = length(xvec);
stepst = length(xnew)-1;
MUhat = [0,0]; % only GARCH, No ARMA
for k=1:2
    %  GARCH(1,1)
    parameters3 = multigarch_AJP(rates(:,k),1,0,1,'GJRGARCH','NORMAL',[],[cov(rates(:,k))*0.05;0.05;0.9]);
    
    if k==1
        GARCHparams1 = parameters3;
    elseif k==2
        GARCHparams2 = parameters3;
    end
end
[GARCHparams1IS, ~, ~, ~, hhat1IS] = multigarch_AJP(rates(:,1),1,0,1,'GJRGARCH','NORMAL',[],GARCHparams1); warning off
[GARCHparams2IS, ~, ~, ~, hhat2IS] = multigarch_AJP(rates(:,2),1,0,1,'GJRGARCH','NORMAL',[],GARCHparams2); warning off
stdresidsIS = rates./sqrt([hhat1IS,hhat2IS]);
stdresidsIS2 = stdresidsIS(:,2);

UhatIS = nan(length(xvec),2);
for mm=1:2
    UhatIS(:,mm) = empiricalCDF(stdresidsIS(:,mm));
end
KAPPAhatIS = nan(4,8);  % parameter (nan if num params<4) ; [Normal, Clayton, RGum, Studt, RGum-GAS, Studt-GAS, Gaussian-GAS, Clayton-GAS] ; [nonparam, skewt]
kappaseries = nan(length(xvec),8);
options = optimset('Display','off','TolCon',10^-3,'TolFun',10^-3,'TolX',10^-3);
% 1. Normal Copula
KAPPAhatIS(1,1) = corrcoef12(norminv(UhatIS));
% 2. Clayton's copula
lower = 0.0001;    warning off;
theta0 = 1;
KAPPAhatIS(1,2) = fmincon('claytonCL',theta0,[],[],[],[],lower,[],[],options,UhatIS);
% 3. Rotated Gumbel copula
lower = 1.1;
upper = 5;
theta0 = 2;
KAPPAhatIS(1,3) = fmincon('gumbelCL',theta0,[],[],[],[],lower,upper,[],options,1-UhatIS);
% 4. Student's t copula
lower = [-0.9 , 0.01 ];
upper = [ 0.9 , 0.45 ];
theta0 = [KAPPAhatIS(1,1);10];
KAPPAhatIS(1:2,4) = fmincon('tcopulaCL2',theta0,[],[],[],[],lower,upper,[],options,UhatIS);
% 5. Rotated Gumbel GAS copula
lower = [-0.1; -3 ; 0.6 ];
upper = [ 0.2 ; 0.5 ; 0.999 ];
theta0 = [log(KAPPAhatIS(1,3)-1)*(1-0.05-0.93);0.05;0.93];
KAPPAhatIS(1:3,5) = fmincon('Rotgumbel_GAS_CL',theta0,[],[],[],[],lower,upper,[],options,UhatIS,KAPPAhatIS(1,3),[KK,HESSgumbel(:,1)]);
[~,kappatISRgum] = Rotgumbel_GAS_CL(KAPPAhatIS(1:3,5),UhatIS,KAPPAhatIS(1,3),[KK,HESSgumbel(:,1)]);
kappaseries(:,5) = kappatISRgum;
% 6. Student's t copula（GAS）
%
lower = [-0.1 ; -0.5 ; 0.6 ; 0.01];
upper = [ 0.2 ; 0.5 ; 0.999 ; 0.45];
theta0 = [log( (0.9999+KAPPAhatIS(1,4))/(0.9999-KAPPAhatIS(1,4)) )*(1-0.05-0.8);0.05;0.9;KAPPAhatIS(2,4)];
%
KAPPAhatIS(1:4,6) = fmincon('tcopula_GAS_CL',theta0,[],[],[],[],lower,upper,[],options,UhatIS,KAPPAhatIS(1,4),RR,NN,HESSstudt);
[~,kappatIStgas] = tcopula_GAS_CL( KAPPAhatIS(1:4,6),UhatIS,KAPPAhatIS(1,4),RR,NN,HESSstudt);
kappaseries(:,6) = kappatIStgas;
% Gaussian copula (GAS)
lower = [-0.2 ; -3 ; 0.6 ];
upper = [ 0.7 ; 3 ; 0.999 ];
theta0 = [log( (0.9999+KAPPAhatIS(1,1))/(0.9999-KAPPAhatIS(1,1)) )*(1-0.05-0.8);0.05;0.9];
KAPPAhatIS(1:3,7) = fmincon('Normalcopula_GAS_CL',theta0,[],[],[],[],lower,upper,[],options,UhatIS,KAPPAhatIS(1,1),[RRnorm ,HESSnorm(:,1)]);
[~, kappatIStnormalgas] = Normalcopula_GAS_CL(KAPPAhatIS(1:3,7),UhatIS,KAPPAhatIS(1,1),[RRnorm ,HESSnorm(:,1)]);
kappaseries(:,7) = kappatIStnormalgas;
%Clayton copula (GAS)
lower = [-0.3 ; -0.3 ; 0.6];
upper = [ 0.2 ; 0.2 ; 0.999 ];
theta0 = [log(KAPPAhatIS(1,3)+1)*(1-0.05-0.93);0.05;0.95];
[KAPPAhatIS(1:3,8),p,q] = fmincon('clayton_GAS_CLm',theta0,[],[],[],[],lower,upper,[],options,UhatIS,KAPPAhatIS(1,2),[KKclayton,HESSclayton(:,1)]);
[~,kappatISclaytongas] =clayton_GAS_CL(KAPPAhatIS(1:3,8),UhatIS,KAPPAhatIS(1,2),[KKclayton,HESSclayton(:,1)]);
kappaseries(:,8) = kappatISclaytongas;
%%%%%%%%%%%%%%%%%       forecasting       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vout = zeros(stepst+1,2);
vout(1,:) = [hhat1IS(end),hhat2IS(end)];
for t = 2:(stepst+1)
    vout(t,1) = GARCHparams1IS(1)+GARCHparams1IS(2)*xnew(t-1)^2+GARCHparams1IS(3)*vout(t-1,1);
    vout(t,2) = GARCHparams2IS(1)+GARCHparams2IS(2)*ynew(t-1)^2+GARCHparams2IS(3)*vout(t-1,2);
end
condvar2out = vout(2:end,2);
stdresidsOS = ratesnew./sqrt(vout);
UhatOOS=nan(stepst+1,2);
for l=1:2
    UhatOOS(:,l) = min(max(empiricalCDF(stdresidsIS(:,l),stdresidsOS(:,l)),1/(wind+1)),(wind)/(wind+1));
end

VaRforccopula = repmat(MUhat(1) + sqrt(vout(2:end,1))*norminv(THETA1),1,4);
% ESforccopula = repmat(MUhat(1) - sqrt(vout(:,1))*normpdf(norminv(THETA1))/(THETA1),1,4);
kappatOOS = nan(stepst+1,4);
kappatOOS(1,:) = kappaseries(end,5:end);
kappatOOS(:,1) = Rotgumbel_GAS_n_step(KAPPAhatIS(1:3,5),kappaseries(end,5),UhatOOS,[KK,HESSgumbel(:,1)],stepst);  %kappa of gumbel copula
kappatOOS(:,2) = tcopula_GAS_n_step( KAPPAhatIS(1:4,6),kappaseries(end,6),UhatOOS,RR,NN,HESSstudt,stepst);   %rho of t copula
kappatOOS(:,3) = gaussian_GAS_n_step(KAPPAhatIS(1:3,7),kappaseries(end,7),UhatOOS,[RRnorm ,HESSnorm(:,1)],stepst);
kappatOOS(:,4) = clayton_GAS_n_step( KAPPAhatIS(1:3,8),kappaseries(end,8),UhatOOS,[KKclayton,HESSclayton(:,1)],stepst);
kappatOOS = kappatOOS(end-stepst+1:end,:);
CoVaRforc = nan(stepst,4);
for k = 1:stepst
    titer = 321*iseed+k*654+987*tseed;
    CoVaRforc(k,1) = CoVaRfrocRGumN(THETA1,THETA2,kappatOOS(k,1),vout(k+1,:),stdresidsOS,MUhat);
    CoVaRforc(k,2) = CoVaRfroctcopN(THETA1,THETA2,kappatOOS(k,2),1/KAPPAhatIS(4,6),vout(k+1,:),stdresidsOS,MUhat,titer);
    CoVaRforc(k,3) = CoVaRfrocNormalN(THETA1,THETA2,kappatOOS(k,3),vout(k+1,:),stdresidsOS,MUhat);
    CoVaRforc(k,4) = CoVaRfrocClaytonN(THETA1,THETA2,kappatOOS(k,4),vout(k+1,:),stdresidsOS,MUhat);
end

tr = 1;
while isnan(sum(kappatOOS(:,4)))==1 || sum(CoVaRforc(:,4)>0) > 0
    try
        reps = 10000;
        KKclayton = [-0.6:0.03:20]';
        HESSclayton = nan(length(KKclayton),1);
        for ii=1:length(KKclayton)
            Usim = clayton_rnd(KKclayton(ii),reps,100*ii+10*tr);  warning('off')
            scoreSIM = LLgrad_1('claytonCLa',KKclayton(ii),Usim);
            HESSclayton(ii,1) = mean(scoreSIM.^2);
        end
        lower = [-0.3 ; -0.3 ; 0.7];
        upper = [ 0.2 ; 0.2 ; 0.999 ];
        theta0 = [log(KAPPAhatIS(1,3)+1)*(1-0.05-0.93);0.05;0.95];
        [KAPPAhatIS(1:3,8),p,q] = fmincon('clayton_GAS_CLm',theta0,[],[],[],[],lower,upper,[],options,UhatIS,KAPPAhatIS(1,2),[KKclayton,HESSclayton(:,1)]);
        [~,kappatISclaytongas] =clayton_GAS_CL(KAPPAhatIS(1:3,8),UhatIS,KAPPAhatIS(1,2),[KKclayton,HESSclayton(:,1)]);
        kappaseries(:,8) = kappatISclaytongas;
        kappatOOS4 = nan(stepst+1,1);
        kappatOOS4 = clayton_GAS_n_step( KAPPAhatIS(1:3,8),kappaseries(end,8),UhatOOS,[KKclayton,HESSclayton(:,1)],stepst);
        kappatOOS(:,4) = kappatOOS4(end-stepst+1:end);
        for k = 1:stepst
            CoVaRforc(k,4) = CoVaRfrocClaytonN(THETA1,THETA2,kappatOOS(k,4),vout(k+1,:),stdresidsOS,MUhat);
        end
    end
    tr = tr+1;
end
grid = [0.005:0.01:0.05];
esgrid = nan(stepst,length(grid),4);

for k = 1:stepst
    for l = 1:length(grid)
        esgrid(k,l,1) = CoVaRfrocRGumN(THETA1,grid(l),kappatOOS(k,1),vout(k+1,:),stdresidsOS,MUhat); warning off
    end
end
for k = 1:stepst
    for l = 1:length(grid)
        titeres = 123*k+456*l+789*iseed+111*tseed;
        esgrid(k,l,2) = CoVaRfroctcopN(THETA1,grid(l),kappatOOS(k,2),1/KAPPAhatIS(4,6),vout(k+1,:),stdresidsOS,MUhat,titeres);warning off
    end
end

for k = 1:stepst
    for l = 1:length(grid)
        esgrid(k,l,3) = CoVaRfrocNormalN(THETA1,grid(l),kappatOOS(k,3),vout(k+1,:),stdresidsOS,MUhat);warning of
    end
end
for k = 1:stepst
    for l = 1:length(grid)
        esgrid(k,l,4) = CoVaRfrocClaytonN(THETA1,grid(l),kappatOOS(k,4),vout(k+1,:),stdresidsOS,MUhat);warning off
    end
end
CoESforc1matGarchN = mean(esgrid,2);
CoESforc1GarchN = nan(stepst,4);
for mm = 1:4
    CoESforc1GarchN(:,mm) = CoESforc1matGarchN(:,:,mm);
end
copulaVaR = VaRforccopula;
copulaCoVaR = CoVaRforc;
copulaCoES = CoESforc1GarchN;
end

