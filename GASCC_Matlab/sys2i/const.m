function [c,ceq] = const(weight,ksen,deltai0,G,THETA,covarvec,cor,condexp)
D = diag(deltai0);
taustar = 0.5*(1-weight'*D*ksen/sqrt(2*weight'*D*G*cor*G*D*weight+(weight'*D*ksen)^2));
deltastar1 = delstar(weight,ksen,deltai0,G,cor);
covar = covarvec*weight;
e = condexp*weight;
c = [e-deltastar1/THETA-covar];
ceq = taustar-THETA;
end

