function [El,gammavec,gammanew,delta] = Elc20scaleconst(params,xvec,yvec,xvar,yvar,THETA2,gamma0,z,delta0,p)
%calcualte params for GAS (initial value), score is calculated using log f
%
omega = params(1);
A = params(2);   %A
B = params(3);   %B
k = find(xvec<xvar);
ksen = (1-THETA2*2)/(THETA2*(1-THETA2));
sig = sqrt(2/(THETA2*(1-THETA2)));
delta(1) = delta0;
gamma(1) = gamma0;
for t = 1:size(xvec)
    dt(t) = yvec(t)-yvar(t);
    if abs(dt(t))<1e-2
        dt(t) = 1e-2;
    end
    delta(t) = max(THETA2*(mean(yvec(xvec<xvar))-(1+exp(gamma(t)))*(yvar(t))),0.005);
    m(t) = (dt(t))^2/(delta(t)*sig)^2;
    d = ksen^2/sig^2;
    p1(t) = -dt(t)*ksen/(delta(t)*sig)^2-1.5/delta(t)+0.5*(k32(sqrt((2+d)*m(t)))/k12(sqrt((2+d)*m(t)))+1)*sqrt((2+d)*m(t))/delta(t);
    p2(t) = -THETA2*exp(gamma(t))*yvar(t);
    if p == 0
        s(t) = p1(t)*p2(t);
    else
        S1 = integral(@(z)Informdelta(z,yvar(t),ksen,sig,delta(t)),-10,10)*(-THETA2*yvar(t)*exp(gamma(t)))^2;
        s(t) = p1(t)*p2(t)/(S1)^p;
    end
    gamma(t+1) = max(min(omega+A*s(t)+B*gamma(t),2),-2);
    Elt(t) = -log(delta(t))+dt(t)*ksen/delta(t)/sig^2-0.5*z(t)*dt(t)^2/(delta(t)*sig)^2;
end
El = -sum(Elt(k));
gammavec = gamma(1:size(xvec));
gammanew = gamma(length(xvec)+1);
if  ~isreal(El)
    El = 1e+100;
end
end

