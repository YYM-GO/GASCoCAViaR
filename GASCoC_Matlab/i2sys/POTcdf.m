function out = POTcdf(x,uL,beta,xi)
%CDF of POT
if x<uL
    out = 0.1*(1+xi*(uL-x)/beta)^(-1/xi);
else
    out = normcdf(x);
end
end

