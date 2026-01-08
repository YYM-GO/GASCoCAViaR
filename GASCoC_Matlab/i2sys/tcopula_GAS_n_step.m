function [rhot,ft] = tcopula_GAS_n_step(theta,lag_rho,Umat,RR,NN,HESSstudt,stepst)
%

RBAR = 0.9999;  % can make this equal to 1, in theory, but issues occur very close to that boundary

w = theta(1);
a = theta(2);
b = theta(3);
nuinv = theta(4);
nu = 1/theta(4);

h = 0.00001;  % step size to use when computing score

rhot = nan(stepst+1,1);
ft = nan(stepst+1,1);
rhot(1) = lag_rho;
ft(1) = log( (RBAR+rhot(1))/(RBAR-rhot(1)) );  % used below.

% now computing the n-step ahead value of rho and f
for i = 2:(stepst+1)
It = interp2(RR',NN',squeeze(HESSstudt(:,:,1,1))',rhot(i-1),nu,'nearest')    ;  % hessian for rho

DELTAt = (-tcopulaCL([rhot(i-1)+h;nu],Umat(i-1,:))--tcopulaCL([rhot(i-1);nu],Umat(i-1,:)))/h    ;         % estimated score for rho. NOTE: my tcopulaCL2 code returns the *neg* log-like, so need to undo that here

drhodf = 2*RBAR*exp(-ft(i-1))/((1+exp(-ft(i-1)))^2);  % used below
Itildat = It / ( drhodf^2) ;                         % estimated hessian for f
DELTAtildat = DELTAt / (  drhodf  );                    % estimated score for f

Itildat = max(Itildat,1e-6);                            % imposing a min value here to stop the likelihood blowing up when Rho is very close to the boundary
DELTAtildat = max(min(DELTAtildat,1e4),-1e4);           % imposing that this is inside (-1e6,1e6)

ft(i) = w + b*ft(i-1) + a*DELTAtildat/sqrt(Itildat);

ft(i) = max(min(ft(i),100),-100);                     % imposing that this is inside (-100,100)
rhot(i) = RBAR*(1-exp(-ft(i)))/(1+exp(-ft(i)));
end

