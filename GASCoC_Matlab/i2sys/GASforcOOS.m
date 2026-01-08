function [CoVaRt,CoESt] = GASforcOOS(params,z,zt,x,xvar,empiricalQuantile,e,THETA2,p,CoMODEL,gamma00)
% expected log likelihood function of AL distribution 
% params: all parameters for the conditioning distribution, including CoCAViaR parameters and GAS parameters.
% z: observations
% zt: EM for AL : E[W^(-1)|Z]
% CoMODEL: 1 CoSAV, 2 CoAS
gasp = params(1:3);
CAVp = params(4:end);

T = length(z);
if CoMODEL == 1
    [~,CoVaRt] = CAVSAV27(x,z,CAVp,empiricalQuantile,THETA2,-xvar);
else
    [~,CoVaRt] = CAVAS28(x,z,CAVp,empiricalQuantile,THETA2,-xvar); 
end

CoVaRt = -CoVaRt;
CoVaRt(CoVaRt>0) = -0.4;
ksen = (1-THETA2*2)/(THETA2*(1-THETA2));
sig = sqrt(2/(THETA2*(1-THETA2)));
d = ksen^2/sig^2;
delta = zeros(size(zt,1),size(zt,2)+1);
gamma = zeros(size(zt,1),size(zt,2)+1);
gamma(:,1) = gamma00;

for t = 1:length(z)
    dt(t) = z(t)-CoVaRt(t);
    if abs(dt(t))<1e-2
        dt(t) = 1e-2;
    end
    for j = 1:1
        delta(j,t) = max(THETA2*(e-(1+exp(gamma(j,t)))*(CoVaRt(t))),0.01);
        m(j,t) = (dt(t))^2/(delta(j,t)*sig)^2;
        d = ksen^2/sig^2;
        p1(j,t) = -dt(t)*ksen/(delta(j,t)*sig)^2-1.5/delta(j,t)+0.5*(k32(sqrt((2+d)*m(j,t)))/k12(sqrt((2+d)*m(j,t)))+1)*sqrt((2+d)*m(j,t))/delta(j,t);
        p2(j,t) = -THETA2*exp(gamma(j,t))*CoVaRt(t);
        if p == 0
            s(j,t) = p1(j,t)*p2(j,t);
        else
            S1(j,t) = integral(@(z)Informdelta(z,CoVaRt(t),ksen,sig,delta(j,t)),-10,10)*(-THETA2*CoVaRt(t)*exp(gamma(j,t)))^2;
            s(j,t) = p1(j,t)*p2(j,t)/(S1(j,t))^p;
        end
        gamma(j,t+1) = min(max(gasp(1)+gasp(2)*s(j,t)+gasp(3)*gamma(j,t),-2),2);
        CoESt(t,j) = (1+exp(gamma(j,t)))*CoVaRt(t);
    end
end




end




