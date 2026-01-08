function [kappat,ft] = clayton_GAS_n_step(theta,lag_k,Umat,hessINFO,stepst)
%

w = theta(1);
a = theta(2);
b = theta(3);

% will model  f[t] = w + b*f[t-1] + a*DELTA[t-1]*S[t-1]
% where f[t] = log(kappa[t]-1),  so that kappa[t] = 1+exp{ f[t] }, ensuring kappa is always above 1

h = 0.00001;  % step size to use when computing score




kappat = nan(stepst+1,1);
ft = nan(stepst+1,1);
kappat(1) = lag_k;
ft(1) = log(kappat(1) +1);

for i = 2:(stepst+1)
    It = interp1(hessINFO(:,1),hessINFO(:,2),kappat(i-1)) ;    % estimated hessian for kappa
DELTAt = (-claytonCL(kappat(i-1) +h,Umat(i-1,:))--claytonCL(kappat(i-1),Umat(i-1,:)))/h;         % estimated score for kappa. NOTE: my gumbelCL code returns the *neg* log-like, so need to undo that here

dkappadf = exp(kappat(i-1));  % used below
Itildat = It / ( dkappadf^2)       ;                  % estimated hessian for f

DELTAtildat = DELTAt / (  dkappadf  )    ;               % estimated score for f

ft(i) = w + b*ft(i-1)+ a*DELTAtildat/sqrt(Itildat);  % dkappadf drops out in this univar case (as Creal et al. note) but I keep it here for completeness
%   kappat(i) = max(min(-1+exp(ft(i)),20),0.0001);
     kappat(i) = -1+exp(ft(i));
end
