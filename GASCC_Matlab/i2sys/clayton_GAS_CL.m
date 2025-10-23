function [CL,kappat,ft] = clayton_GAS_CL(theta,data,kappa0,hessINFO)
%function CL = rotgumbel_GAS_CL(theta,data)
% The negative copula log-likelihood of a member of ROTATED Gumbel's family
% (from Joe (1997), p142), where the parameter varies through time according to the "GAS" model for Creal, Koopman and Lucas (2011, JAE)
% This code is modified according to Patton's code
% INPUTS:   theta = [w,a,b], the parameters of the GAS specification
%           data, a Tx2 matrix of Unif(0,1) random variables.
%           kappa0, a scalar, the value of the Rotated Gumbel copula to use as the starting value (perhaps the estimate from a constant version of this model)
%           hessINFO, a kx2 matrix containing the grid of values at which the hessian was evaluated, and the value of the hessian at those points
%
% OUTPUTS:  CL, a scalar, the negative log-likelihood
%
%

%theta'

if sum(isnan(theta))==0
    
    
    T = size(data,1);
    %
    
    w = theta(1);
    a = theta(2);
    b = theta(3);
    
    
    h = 0.00001;  % step size to use when computing score
    
    ft = nan(T,1);
    kappat = nan(T,1);
    kappat(1) = kappa0;
    ft(1) = log(kappat(1)+1);
    
    for tt=2:T
        It = interp1(hessINFO(:,1),hessINFO(:,2),kappat(tt-1));     % estimated hessian for kappa
        
        DELTAt = (-claytonCL(kappat(tt-1)+h,data(tt-1,:))--claytonCL(kappat(tt-1),data(tt-1,:)))/h;         % estimated score for kappa. NOTE: my gumbelCL code returns the *neg* log-like, so need to undo that here
        
        dkappadf = exp(kappat(tt-1));  % used below
        Itildat = It / ( dkappadf^2) ;                         % estimated hessian for f
        
        DELTAtildat = DELTAt / (  dkappadf  );                    % estimated score for f
        
        ft(tt) = w + b*ft(tt-1) + a*DELTAtildat/sqrt(Itildat);  % dkappadf drops out in this univar case (as Creal et al. note) but I keep it here for completeness
        k = max(min(-1+exp(ft(tt)),14.9),0.0001);
        kappat(tt) = k;
    end
    
    
    
    CL = log(data(:,1).*data(:,2)).*(1+kappat);
    CL = CL + (2+1./kappat).*log( max((data(:,1).^(-kappat)) + (data(:,2).^(-kappat)) -1,0.001));
    CL = log(1+kappat) - CL;
    CL = sum(CL);
    CL = -CL;
    
    %
    if isnan(CL) || isinf(CL)
        CL = 1e8;
    end
    
else
    CL = 1e7;
end