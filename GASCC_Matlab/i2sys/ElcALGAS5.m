function [El,gamma,delta,CoVaRt,CoESt,s] = ElcALGAS5(params,z,zt,x,xvar,es,THETA1,THETA2,p,CoMODEL,gamma00,delta00,txlv)
% expected log likelihood function of AL distribution (No regime switching, only GAS)
% params: all parameters for the conditioning distribution, including CoCAViaR parameters and GAS parameters.
% z: observations
% zt: EM for AL : E[W^(-1)|Z]
% CoMODEL: 1 SAV, 2 AS
gasp = params(1:3);
CAVp = params(4:end);
% empiricalQuantile = quantile(z(find(x<quantile(x,THETA1))),THETA2);
empiricalQuantile = quantile(z(find(x<xvar)),THETA2);
T = length(z);
if CoMODEL == 1
    [~,CoVaRt] = CAVSAV27(x,z,CAVp,empiricalQuantile,THETA2,-xvar); CoVaRt = -CoVaRt;
else
    [~,CoVaRt] = CAVAS28(x,z,CAVp,empiricalQuantile,THETA2,-xvar); CoVaRt = -CoVaRt;
end


% es = mean(z(find(x<quantile(x,THETA1))));
% gamma200 = log(mean(z(find(z<quantile(z(x<quantile(x,THETA1)),THETA2))))/quantile(z(x<quantile(x,THETA1)),THETA2)-1);
% % delta200 = -THETA2*(1+exp(gamma200))*CoVaRt(1);
% delta200 = THETA2*(mean(z(find(x<quantile(x,THETA1))))-(1+exp(gamma200))*CoVaRt(1));


ksen = (1-THETA2*2)/(THETA2*(1-THETA2));
sig = sqrt(2/(THETA2*(1-THETA2)));
d = ksen^2/sig^2;
delta = zeros(size(zt,1),size(zt,2)+1);
gamma = zeros(size(zt,1),size(zt,2)+1);


delta(:,1) = delta00;
gamma(:,1) = gamma00;

for t = 1:length(z)
    dt(t) = z(t)-CoVaRt(t);
    if abs(dt(t))<1e-2
        dt(t) = 1e-2;
    end
    for j = 1:1
        delta(j,t) = max(THETA2*(mean(z(x<xvar))-(1+exp(gamma(j,t)))*(CoVaRt(t))),0.005);
%         delta(j,t) = THETA2*(es-(1+exp(gamma(j,t)))*(CoVaRt(t)));
        m(j,t) = (dt(t))^2/(delta(j,t)*sig)^2;
        d = ksen^2/sig^2;
        %     delta(t) = max(-THETA2*(1+exp(gamma(t)))*VaR(t),0.03);
        p1(j,t) = -dt(t)*ksen/(delta(j,t)*sig)^2-1.5/delta(j,t)+0.5*(k32(sqrt((2+d)*m(j,t)))/k12(sqrt((2+d)*m(j,t)))+1)*sqrt((2+d)*m(j,t))/delta(j,t);
        p2(j,t) = -THETA2*exp(gamma(j,t))*CoVaRt(t);
        if p == 0
            s(j,t) = p1(j,t)*p2(j,t);
        else
            S1(j,t) = integral(@(z)Informdelta(z,CoVaRt(t),ksen,sig,delta(j,t)),-10,10)*(-THETA2*CoVaRt(t)*exp(gamma(j,t)))^2;
            s(j,t) = p1(j,t)*p2(j,t)/(S1(j,t))^p;
        end
        %     gamma(t+1) = max(min(omega+A*s(t)+B*gamma(t),1.5),-0.5);
        gamma(j,t+1) = min(gasp(1)+gasp(2)*s(j,t)+gasp(3)*gamma(j,t),3);
        
        %     Elt(t) = -log(delta(t))+dt(t)*ksen/delta(t)/sig^2-0.5*z(t)*dt(t)^2/(delta(t)*sig)^2-0.5*ksen^2*u(t)/sig^2;
        Elt0(j,t) = -log(delta(j,t))+dt(t)*ksen/delta(j,t)/sig^2-0.5*zt(t)*dt(t)^2/(delta(j,t)*sig)^2;
        CoESrt(t,j) = (1+exp(gamma(j,t)))*CoVaRt(t);
    end
    CoESt(t) = CoESrt(t,1);
    Elt(t) = -Elt0(:,t);
end

El = sum(Elt(txlv));
% if El == Inf | ~isreal(El)
%     El = 1e+100;
% end


end


