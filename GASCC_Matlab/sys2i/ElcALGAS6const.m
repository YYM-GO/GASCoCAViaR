function [El,gamma,delta,CoVaRt,CoESt] = ElcALGAS6const(params,z,zt,x,xvar,empiricalQuantile,e,THETA1,THETA2,CoMODEL,txlv)
% expected log likelihood function of AL distribution
% params: all parameters for the conditioning distribution.
% z: observations
% zt: EM for AL : E[W^(-1)|Z]
% CoMODEL: 1 CoSAV, 2 CoAS
gamma = params(1);
CAVp = params(2:end);
T = length(z);
if CoMODEL == 1
    [~,CoVaRt] = CAVSAV27(x,z,CAVp,empiricalQuantile,THETA2,-xvar);
else
    [~,CoVaRt] = CAVAS28(x,z,CAVp,empiricalQuantile,THETA2,-xvar);
end
CoVaRt = -max(CoVaRt,0.5);
ksen = (1-THETA2*2)/(THETA2*(1-THETA2));
sig = sqrt(2/(THETA2*(1-THETA2)));
for t = 1:length(z)
    dt(t) = z(t)-CoVaRt(t);
    if abs(dt(t))<1e-2
        dt(t) = 1e-2;
    end
    delta(t) = max(THETA2*(e-(1+exp(gamma))*(CoVaRt(t))),0.01);
    m(t) = (dt(t))^2/(delta(t)*sig)^2;
    Elt0(t) = -log(delta(t))+dt(t)*ksen/delta(t)/sig^2-0.5*zt(t)*dt(t)^2/(delta(t)*sig)^2;
    CoESt(t) = (1+exp(gamma))*CoVaRt(t);
    Elt(t) = -Elt0(:,t);
end
El = sum(Elt(txlv));
delta = reshape(delta,[],1);
end


