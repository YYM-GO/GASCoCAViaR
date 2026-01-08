function [El,gamma,delta,CoVaRt,CoESt] = ElcALGAS6(params,z,zt,x,xvar,empiricalQuantile,e,THETA1,THETA2,p,CoMODEL,gamma00,delta00,txlv)
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
CoVaRt = -abs(CoVaRt);
ksen = (1-THETA2*2)/(THETA2*(1-THETA2));
sig = sqrt(2/(THETA2*(1-THETA2)));
delta = zeros(size(zt,1),size(zt,2)+1);
gamma = zeros(size(zt,1),size(zt,2)+1);
dt = zeros(length(z),1);
p1 = zeros(length(z),1);
p2 = zeros(length(z),1);
s = zeros(length(z),1);
Elt0 = zeros(length(z),1);
CoESt = zeros(length(z),1);
gamma(:,1) = gamma00;
for t = 1:length(z)
    dt(t) = z(t)-CoVaRt(t);
    if abs(dt(t))<1e-2
        dt(t) = 1e-2;
    end
    delta(t) = max(THETA2*(e-(1+exp(gamma(t)))*(CoVaRt(t))),0.005);
    p1(t) = -1/delta(t)-dt(t)*ksen/(delta(t)*sig)^2+sqrt(2*sig^2+ksen^2)/((delta(t)*sig)^2)*abs(dt(t));
    p2(t) = -THETA2*exp(gamma(t))*CoVaRt(t);
    if p == 0
        s(t) = p1(t)*p2(t);
    else
        S1(t) = integral(@(z)Informdelta(z,CoVaRt(t),ksen,sig,delta(t)),-10,10)*(-THETA2*CoVaRt(t)*exp(gamma(t)))^2;
        s(t) = p1(t)*p2(t)/(S1(t))^p;
    end
    gamma(t+1) = min(max(gasp(1)+gasp(2)*s(t)+gasp(3)*gamma(t),-2),2);
    Elt0(t) = -log(delta(t))+dt(t)*ksen/delta(t)/sig^2-0.5*zt(t)*dt(t)^2/(delta(t)*sig)^2;
    CoESt(t) = (1+exp(gamma(t)))*CoVaRt(t);
end
El = sum(-Elt0(txlv));
end


