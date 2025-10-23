function [CL,rhot,ft] = Normalcopula_GAS_CL(theta,data,rho0,HESSnorm)

% The negative copula log-likelihood of the gaussian copula, 

%
% INPUTS:   theta = [w,a,b], the parameters of the GAS specification, and the *inverse* DoF parameter
%           data, a Tx2 matrix of Unif(0,1) random variables.
%           rho0, a scalar, the value of the correlation parameter to use as the starting value (perhaps the estimate from a constant version of this model)
%           RR, a k1x1 vector of values of rho at which the hessian was computed
%           NN, a k2x1 vector of values of nu at which the hessian was computed
%           HESSstudt, a k1 x k2 x 2 x 2 matrix, containing the 2x2 hessian for each combination of values of [rho,nu]
%


RBAR = 0.9999;  % can make this equal to 1, in theory, but issues occur very close to that boundary

%theta'

if sum(isnan(theta))==0
    
    T = size(data,1);
    
    w = theta(1);
    a = theta(2);
    b = theta(3);

    
    h = 0.0001;  % step size to use when computing score

    % generating the time series of rho
    ft = nan(T,1);
    rhot = nan(T,1);
    rhot(1) = rho0;
    ft(1) = log( (RBAR+rhot(1))/(RBAR-rhot(1)) );
%     count = zeros(1,2);
    for tt=2:T

        %rhot(tt-1)
%         hessINFO=HESSnorm
%         It = interp2(RR',squeeze(HESSnorm)',rhot(tt-1),'nearest')    ;  % hessian for rho
        It = interp1(HESSnorm(:,1),HESSnorm(:,2),rhot(tt-1));
        DELTAt = (-NormalCopula_CL(rhot(tt-1)+h,data(tt-1,:))--NormalCopula_CL(rhot(tt-1),data(tt-1,:)))/h    ;         % estimated score for rho. NOTE: my tcopulaCL2 code returns the *neg* log-like, so need to undo that here

        
        drhodf = 2*RBAR*exp(-ft(tt-1))/((1+exp(-ft(tt-1)))^2);  % used below
        Itildat = It / ( drhodf^2) ;                         % estimated hessian for f
        DELTAtildat = DELTAt / (  drhodf  );                    % estimated score for f

        
        % the search algorithm sometimes looks in parts of te parameter space where things blow up. so put in the following checks to stop this
        % happening. Those parts are never selected, so this is just to keep the search algorithm going, rather than stopping with a NaN of Inf error
        % message.
        % code below is used to count how many times the trunaction is needed. comment out for speed
%         if Itildat<=1e-6
%             count(1)=count(1)+1;
%         end
%         if abs(DELTAtildat)>=1e4
%             count(2)=count(2)+1;
%         end

        % next two lines, and the similar line for ft(tt) are useful when running simulations and want to automate estimation.
        % parameters will generally not be such that these bounds are binding (at least, they should not be) but during the numerical optimization they are
        % hit sometimes and this helps the code keep going
        Itildat = max(Itildat,1e-6);                            % imposing a min value here to stop the likelihood blowing up when Rho is very close to the boundary
        DELTAtildat = max(min(DELTAtildat,1e4),-1e4);           % imposing that this is inside (-1e6,1e6)

        ft(tt) = w + b*ft(tt-1) + a*DELTAtildat/sqrt(Itildat);
        
        ft(tt) = max(min(ft(tt),100),-100);                     % imposing that this is inside (-100,100)
        rhot(tt) = RBAR*(1-exp(-ft(tt)))/(1+exp(-ft(tt)));
        
%        [tt,rhot(tt),ft(tt),Itildat,DELTAtildat]
    end
%    count
    
    % transforming data to feed into likelihood
		
        x = norminv(data(:,1));
        y = norminv(data(:,2));
    

    % the norma; likelihood
    CL = -1*(2*(1-rhot.^2)).^(-1).*(x.^2+y.^2-2*rhot.*x.*y);
    CL = CL + 0.5*(x.^2+y.^2);
    CL = sum(CL) - size(x,1)/2*log(1-rhot.^2);
    CL = -CL;
   CL = sum(CL);
    if isnan(CL) || isinf(CL)
        CL = 1e8;
    end
    
else
    CL = 1e7;
end