function out = POTinv(uL,beta,xi,tau)
%calculate the VaR of mixed POT (with Gaussian distribution)
if tau > 0.1
    out = norminv(tau);
else
    out = uL-((tau*10)^(-xi)-1)*beta/xi;
end

end

