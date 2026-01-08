function out = POTpdf2(x,uL,NuL,N,betaL,xiL)
% the distribution funciton of POT 
% ref:Portfolio optimization based on GARCH-EVT-Copula forecasting models

    out = NuL/N*(1+xiL*(uL-x)/betaL)^(-1/xiL-1)/betaL;
end