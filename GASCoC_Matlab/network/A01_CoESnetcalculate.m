close
clear
clc


load shenwan
load hushen        %CSI 300


[T,N] = size(shenwan);
THETA1 = 0.05; THETA2 = 0.05; %probability level for VaR and CoVaR
wind = 1400;
Tout = T-wind;
stepst0 = 100;
nstepst = ceil(Tout/stepst0);
p = 0;      % the scale of the score, we can choose 0, 1/2 or 1. 0 is the fastest.
epsil = 1e-6;
TolX = 1e-6;
miter = 50;


CoMODEL = 2;
Covarnet = zeros(Tout,N,N);
Coesnet = zeros(Tout,N,N);

bar = waitbar(0,'Please wait...');
for t = 1:nstepst
    if t ~= nstepst
        stepst = stepst0;
    else
        stepst = Tout-(nstepst-1)*stepst0;
    end
    for i = 1:N
        tic
        xvec = shenwan(stepst*(t-1)+1:stepst*(t-1)+wind,i);
        xvecall = [xvec;shenwan(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst,i)];
        xnew = xvecall(end-stepst:end);
        [VaRin,VaRforc,bestall(t,i)] = VaROOSchoose(xvec,xnew,THETA1); warning off
        parfor j = 1:N
            yvec = shenwan(stepst*(t-1)+1:stepst*(t-1)+wind,j);
            yvecall = [yvec;shenwan(stepst*(t-1)+wind+1:stepst*(t-1)+wind+stepst,j)];
            ynew = yvecall(end-stepst+1:end);
            [CoVoRforc,CoESforc] = CoforcOOS(xvec,yvec,VaRin,VaRforc,THETA1,THETA2,CoMODEL,p,epsil,miter,xnew,ynew);
            CoVoRforcn(:,j) = CoVoRforc;
            CoESforcn(:,j) = CoESforc;
        end
        if t~= nstepst
            Covarnet(stepst*(t-1)+1:stepst*t,i,:) = CoVoRforcn;
            Coesnet(stepst*(t-1)+1:stepst*t,i,:) = CoESforcn;
        else
            Covarnet(end-stepst+1:end,i,:) = CoVoRforcn;
            Coesnet(end-stepst+1:end,i,:) = CoESforcn;
        end
        [t,i,toc]
    end
    waitbar(t/nstepst)
end

save Covarnet.mat Covarnet
save Coesnet.mat Coesnet

%%       Chen et al.(2019)
load marketvalue
marketvalue = marketvalue*1e-11;
mean(marketvalue)';

load shenwan
load hushen
load Covarnet
load Coesnet

Covar = Covarnet;  %T*N*N
CoES = Coesnet;   %T*N*N

[T,N] = size(shenwan);
wind = 1400;
Tout = T-wind;

THETA1 = 0.05; THETA2 = 0.05; %probability level for VaR and CoVaR
ksen1 = (1-THETA1*2)/(THETA1*(1-THETA1));
ksen2 = (1-THETA2*2)/(THETA2*(1-THETA2));
sig1 = sqrt(2/(THETA1*(1-THETA1)));
sig2 = sqrt(2/(THETA2*(1-THETA2)));

for t = 1:size(CoES,1)
    for i = 1:N
        Coesivec(:,i,t) = reshape(CoES(t,:,i),[],1);
    end
end

parfor t = 1:size(CoES,1)
    for i = 1:N
        for j = 1:N
            Xit = reshape(Coesivec(:,i,t),[],1);
            Yit = reshape(Coesivec(:,j,t),[],1);
            rhoijt(i,j,t) = corr(Xit,Yit);
        end
    end
%     t
end

Adjacency = nan(N,N,Tout);
wall = nan(N,N,Tout);
for t = 1:size(CoES,1)
    mmarketvalue = log(marketvalue(i,:));
    rhopos = []; rhoneg = []; indexpos = []; indexneg = []; phipos = []; phineg = []; Deltapos = []; Deltaneg = [];
    Lposdelta = []; Sposdelta = [];  Lnegdelta = []; Snegdelta = [];
    Lrhopos_rho = []; Lrhopos_location = [];Lrhoneg_rho = []; Lrhoneg_location = [];
    Srhopos_rho = []; Srhopos_location = [];Srhoneg_rho = []; Srhoneg_location = [];
    Lposrhopair = []; Lnegrhopair = [];  Sposrhopair = []; Snegrhopair = [];
    a0 = zeros(N,N);
    w0 = zeros(N,N);
    rhoall(:,t) = tritovec(rhoijt(:,:,t));  % make the upper triangle of rho to vector
    [rhot,indexrho] = sort(rhoall(:,t));
    n1(t) = length(find(rhot>0));
    n2(t) = N*(N-1)/2-n1(t);
    rhoneg = rhot(1:n2(t));
    rhopos = rhot(n2(t)+1:end);
    indexneg = indexrho(1:n2(t));
    indexpos = indexrho(n2(t)+1:end);
    phipos = normcdf(sqrt(N)*rhopos);
    phineg = normcdf(sqrt(N)*rhoneg);
    Deltapos = diff(phipos);
    Deltaneg = diff(phineg);
    [deltordpos,indexdeltordpos] = sort(Deltapos);
    [deltordneg,indexdeltordneg] = sort(Deltaneg);
    [thetapos(t),Lposdelta,Sposdelta,nsmallpos] = threshold2(Deltapos);
    [thetaneg(t),Lnegdelta,Snegdelta,nsmallneg] = threshold2((Deltaneg));
    
    Lrhopos_rho = rhopos(nsmallpos+1:end);
    Lrhopos_location = indexpos(nsmallpos+1:end);
    Lrhoneg_location = indexneg(1:(n2(t)-nsmallneg));
    Lrhoneg_rho = rhoneg(1:(n2(t)-nsmallneg));
    
    Srhopos_rho = rhopos(1:nsmallpos);
    Srhopos_location = indexpos(1:nsmallpos);
    Srhoneg_location = indexneg(end-nsmallneg+1:end);
    Srhoneg_rho = rhoneg(end-nsmallneg+1:end);
    
    Lposrhopair = findposition2(Lrhopos_location);
    Lnegrhopair = findposition2(Lrhoneg_location);
    Sposrhopair = findposition2(Srhopos_location);
    Snegrhopair = findposition2(Srhoneg_location);
    % find which pair is Lpos\Lneg\Spos\Sneg
    
    Lpos{t} = Lposrhopair;
    Spos{t} = Sposrhopair;
    Lneg{t} = Lnegrhopair;
    Sneg{t} = Snegrhopair;
    
    for k = 1:size(Lposrhopair,1)
        m = Lposrhopair(k,:);
        i = m(1);j = m(2);
        a0(i,j) = 1;
        a0(j,i) = 1;
    end
    for k = 1:size(Lnegrhopair,1)
        m = Lnegrhopair(k,:);
        i = m(1);j = m(2);
        a0(i,j) = -1;
        a0(j,i) = -1;
    end
    Adjacency(:,:,t) = a0;
    %  --------------   test   ----------------------------
    mSpos = length(Sposdelta);
    mSneg = length(Snegdelta);
    mLpos = length(Lposdelta);
    mLneg = length(Lnegdelta);
    w0 = a0.*abs(rhoijt(:,:,t));
    wall(:,:,t) = w0;
    % --------------score------------------------
    for i = 1:N
        sc1(t,i) = 2*a0(i,:)*mmarketvalue'*mmarketvalue(i);
    end
    S(t) = mmarketvalue*a0*mmarketvalue';
    rationeg(t) = n2(t)/(N*(N-1)/2);
    ratioLpos(t) = (n1(t)-nsmallpos)/(N*(N-1)/2);  %ND+ (10)
    ratioLneg(t) = (n2(t)-nsmallneg)/(N*(N-1)/2);  %ND-
    for j = 1:N
        rhoposvec(t,j) = mean(rhoijt(find(rhoijt(:,j,t)>0),j,t));
        rhonegvec(t,j) = mean(rhoijt(find(rhoijt(:,j,t)<0),j,t));
        rhoLipos(t,j) = mean(w0(find(a0(:,j)==1),j));
        rhoLineg(t,j) = mean((w0(find(a0(:,j)==-1),j)));
        rLipos(t,j) = length(find(a0(:,j)==1))/(N-1);    %DC+
        rLineg(t,j) = length(find(a0(:,j)==-1))/(N-1);   %DC-(11)
    end
    densitypos(t) = mean(w0(find(a0(:,:)==1)));
    densityneg(t) = -mean(w0(find(a0(:,:)==-1)));
    density1(t) = mean(mean(abs(w0)));
    density2(t) = mean(mean(abs(rhoall(:,t))));
    t
end
save S.mat S
save sc1.mat sc1
save Adjacency.mat Adjacency

