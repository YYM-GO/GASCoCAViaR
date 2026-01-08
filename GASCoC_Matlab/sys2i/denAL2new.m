function [f,gamma,delta] = denAL2new(z,gasp,gamma0,delta0,THETA2,Q,es0,p)
%calculate the density function of AL distribution for all t=1,...,T  and l = 1,...,L
% z:observations
% gamma0/delta0: initial values
%p: scale of score- 0, 0,5 or 1
% es: expected value of Y given X<VaR(X)
ksen = (1-THETA2*2)/(THETA2*(1-THETA2));
sig = sqrt(2/(THETA2*(1-THETA2)));
d = ksen^2/sig^2;
delta = zeros(size(gasp,2),length(z));
gamma = zeros(size(gasp,2),length(z)+1);
delta(:,1) = delta0;
gamma(:,1) = gamma0;
m = zeros(size(z,1),1);
f = zeros(size(z,1),1);
p1 = zeros(size(z,1),1);
p2 = zeros(size(z,1),1);
s = zeros(size(z,1),1);
for t = 1:size(z,1)
    delta(t) = max(THETA2*(es0-(1+exp(gamma(t)))*Q(t)),0.01);
    m(t) = (z(t)-Q(t))^2/(delta(t)*sig)^2;
    f(t) = sqrt(2/pi)*exp((z(t)-Q(t))*ksen/(delta(t)*sig^2))/(delta(t)*sig)*(m(t)/(2+d))^0.25*k12(sqrt((2+d)*m(t)));
    p1(t) = -1/delta(t)-(z(t)-Q(t))*ksen/(delta(t)*sig)^2+sqrt(2*sig^2+ksen^2)/((delta(t)*sig)^2)*abs((z(t)-Q(t)));
    p2(t) = -THETA2*exp(gamma(t))*Q(t);
    if p == 0
        s(t) = p1(t)*p2(t);
    else
        S1(t) = integral(@(z)Informdelta(z,Q(t),ksen,sig,delta(t)),-10,10)*(-THETA2*Q(t)*exp(gamma(t)))^2;
        s(t) = p1(t)*p2(t)/(S1(t))^p;
    end
    gamma(t+1) = max(min(gasp(1)+gasp(2)*s(t)+gasp(3)*gamma(t),2),-2);    
end
end



