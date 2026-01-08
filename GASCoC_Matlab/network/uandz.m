function [u,z] = uandz(ratevec,VaR,deltavec,THETA)
% calculate the weight in EM algorithm.
ksen = (1-THETA*2)/(THETA*(1-THETA));
sig = sqrt(2/(THETA*(1-THETA)));
for t = 1:length(ratevec)
    dt(t) = ratevec(t)-VaR(t);
    if abs(dt(t))<1e-2
        dt(t) = 1e-2;
    end
    m(t) = (dt(t))^2/(deltavec(t)*sig)^2;
    d = ksen^2/sig^2;
    u(t) = (m(t)/(2+d))^0.5*k32(sqrt((2+d)*m(t)))/k12(sqrt((2+d)*m(t)));
    z(t) = sqrt((2+d)/m(t))*k32(sqrt((2+d)*m(t)))/k12(sqrt((2+d)*m(t)))-1/m(t);
end
end

