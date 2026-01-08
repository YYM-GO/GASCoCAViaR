close
clear
clc

load shenwan       %Shenwan industry index
load hushen        %CSI 300

[T,N] = size(shenwan);
THETA1 = 0.05; THETA2 = 0.05; %probability level for VaR and CoVaR

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
options = optimset('Display','off','TolCon',10^-4,'TolFun',10^-4,'TolX',10^-4,'Algorithm','interior-point');
e = nan(T,N,4);
f = nan(T,N,4);
o = nan(T,N,4);
m = nan(T,N,4);
tic
xvec = hushen;
MUhat = [0,0];
[GARCHparams1, likelihood1, ~, ~, hhat1] = multigarch_AJP(xvec,1,0,1,'GJRGARCH','NORMAL',[],[cov(xvec)*0.05;0.05;0.9]);
[GARCHparams1IS, ~, ~, ~, hhat1IS] = multigarch_AJP(xvec,1,0,1,'GJRGARCH','NORMAL',[],GARCHparams1);
for k = 1:T
    VaRforcgarch(k,:) = repmat(MUhat(1) + sqrt(hhat1IS(k))*norminv(THETA1),1,4);
    ESforcgarch(k,:) = repmat(MUhat(1) - sqrt(hhat1IS(k))*normpdf(norminv(THETA1))/(THETA1),1,4);
end
parfor i = 1:N
    tic
    yvec = shenwan(:,i);
    resids = [xvec,yvec];
    [GARCHparams2IS, ~, ~, ~, hhat2IS] = multigarch_AJP(yvec,1,0,1,'GJRGARCH','NORMAL',[],[cov(xvec)*0.05;0.05;0.9]);
    stdresidsIS = resids./sqrt([hhat1IS,hhat2IS]);
    GARCHparamsIS = [GARCHparams1IS,GARCHparams2IS];    
    UhatIS = nan(T,2);    
    for mm=1:2
        UhatIS(:,mm) = empiricalCDF(stdresidsIS(:,mm));
    end
    v=[hhat1IS,hhat2IS];
    e(:,i,:) = VaRforcgarch;
    f(:,i,:) = ESforcgarch;
    KAPPAhatIS = nan(4,8);  % parameter (nan if num params<4) ; [Normal, Clayton, RGum, Studt, RGum-GAS, Studt-GAS, Gaussian-GAS, Clayton-GAS] ; [nonparam, skewt]
    kappaseries = nan(T,8);
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
    lower = [-0.1 ; -0.5 ; 0.6 ; 0.01];
    upper = [ 0.2 ; 0.5 ; 0.999 ; 0.45];
    theta0 = [log( (0.9999+KAPPAhatIS(1,4))/(0.9999-KAPPAhatIS(1,4)) )*(1-0.05-0.8);0.05;0.9;KAPPAhatIS(2,4)];
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
    lower = [-1 ; -1 ; 0.6];
    upper = [ 1 ; 1 ; 0.999 ];
    theta0 = [log(KAPPAhatIS(1,3)+1)*(1-0.05-0.93);0.05;0.93];
    [KAPPAhatIS(1:3,8),p,q] = fmincon('clayton_GAS_CLm',theta0,[],[],[],[],lower,upper,[],options,UhatIS,KAPPAhatIS(1,2),[KKclayton,HESSclayton(:,1)]);
    [~,kappatISclaytongas] =clayton_GAS_CL(KAPPAhatIS(1:3,8),UhatIS,KAPPAhatIS(1,2),[KKclayton,HESSclayton(:,1)]);
    kappaseries(:,8) = kappatISclaytongas;
    %%%%%-------------------------- covar and coes -------------------------
    kappat = kappaseries(:,5:end);
    % one step forecast of CoVaR of dim 2 given dim 1
    CoVaRforc = zeros(T,4);
    for k = 1:T
        titer = 321*i+k*654;
        CoVaRforc(k,1) = CoVaRfrocRGumN(THETA1,THETA2,kappat(k,1),v(k,:),stdresidsIS,MUhat);
        CoVaRforc(k,2) = CoVaRfroctcopN(THETA1,THETA2,kappat(k,2),1/KAPPAhatIS(4,6),v(k,:),stdresidsIS,MUhat,titer);
        CoVaRforc(k,3) = CoVaRfrocNormalN(THETA1,THETA2,kappat(k,3),v(k,:),stdresidsIS,MUhat);
        CoVaRforc(k,4) = CoVaRfrocClaytonN(THETA1,THETA2,kappat(k,4),v(k,:),stdresidsIS,MUhat);
    end    
    grid = [0.005:0.01:THETA2];
    esgrid = nan(T,length(grid),4);  
    for k = 1:T
        for l = 1:length(grid)
            esgrid(k,l,1) = CoVaRfrocRGumN(THETA1,grid(l),kappat(k,1),v(k,:),stdresidsIS,MUhat); warning off
        end
    end   
    for k = 1:T
        for l = 1:length(grid)
            titeres = 123*k+456*l+789*i;
            esgrid(k,l,2) = CoVaRfroctcopN(THETA1,grid(l),kappat(k,2),1/KAPPAhatIS(4,6),v(k,:),stdresidsIS,MUhat,titeres);warning off
        end
    end    
    for k = 1:T
        for l = 1:length(grid)
            esgrid(k,l,3) = CoVaRfrocNormalN(THETA1,grid(l),kappat(k,3),v(k,:),stdresidsIS,MUhat);warning off
        end
    end
    for k = 1:T
        for l = 1:length(grid)
            esgrid(k,l,4) = CoVaRfrocClaytonN(THETA1,grid(l),kappat(k,4),v(k,:),stdresidsIS,MUhat);warning off
        end
    end
    CoESforc1matGarchN = mean(esgrid,2);
    CoESforc1GarchN = nan(T,4);
    for mm = 1:4
        CoESforc1GarchN(:,mm) = CoESforc1matGarchN(:,:,mm);
    end
    %----------------summary---------------------------------    
    o(:,i,:) = CoVaRforc;
    m(:,i,:) = CoESforc1GarchN;
    [i,toc]
end
toc

CopulaVaRs2iGaufull = e;
CopulaESs2iGaufull = f;
CopulaCoVaRs2iGaufull = o;
CopulaCoESs2iGaufull = m;


save CopulaVaRs2iGaufull.mat CopulaVaRs2iGaufull
save CopulaESs2iGaufull.mat CopulaESs2iGaufull
save CopulaCoVaRs2iGaufull.mat CopulaCoVaRs2iGaufull
save CopulaCoESs2iGaufull.mat CopulaCoESs2iGaufull

