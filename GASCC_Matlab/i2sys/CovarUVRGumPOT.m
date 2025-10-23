function out = CovarUVRGumPOT(THETA1,THETA2,kappa,stdresidsIS,uL,beta,xi)
%calculate the covar of U<=THETA2|V<=THETA1
% see Measuring systemic risk in the European bankingsector: a copula CoVaR
% approach (14)

z1 = (-log(THETA1*THETA2)).^kappa-(-log(THETA1)).^kappa;
z2 = exp(-(z1).^(1/kappa));
out = POTinv(uL,beta,xi,z2);



end