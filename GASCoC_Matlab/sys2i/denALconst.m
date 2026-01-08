function [f,gamma,delta] = denALconst(z,gamma,THETA2,Q,es0)
%calculate the density function of AL distribution for all t=1,...,T  and l = 1,...,L
% z:observations
%p: scale of score- 0, 0,5 or 1
% es: expected value of Y given X<VaR(X)
ksen = (1-THETA2*2)/(THETA2*(1-THETA2));
sig = sqrt(2/(THETA2*(1-THETA2)));
d = ksen^2/sig^2;
for t = 1:size(z,1)
    delta(t) = max(THETA2*(es0-(1+exp(gamma))*Q(t)),0.01);
    m(t) = (z(t)-Q(t))^2/(delta(t)*sig)^2;
    f(t) = sqrt(2/pi)*exp((z(t)-Q(t))*ksen/(delta(t)*sig^2))/(delta(t)*sig)*(m(t)/(2+d))^0.25*k12(sqrt((2+d)*m(t)));
end
end



