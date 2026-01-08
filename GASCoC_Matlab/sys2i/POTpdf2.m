function out = POTpdf2(x,uL,NuL,N,betaL,xiL)
% the distribution funciton of POT 
% ref:Portfolio optimization based on GARCH-EVT-Copula forecasting models

    out = NuL/N*(1+xiL*(uL-x)/betaL)^(-1/xiL-1)/betaL;
% else if x> uR
%         out = NuR/N*(1+xiR*(x-uR)/betaR)^(-1/xiR-1)/betaR;
%     else
%         out = normpdf(x);
%     end
% end




end