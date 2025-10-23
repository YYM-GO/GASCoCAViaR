function out = CovarUVClaytonN(THETA1,THETA2,kappa,stdresidsIS)
%calculate the covar of U<=THETA2|V<=THETA1
% see Measuring systemic risk in the European bankingsector: a copula CoVaR
% approach (14)

z1 = ((THETA1*THETA2)^(-kappa)-THETA1^(-kappa))/kappa;
z2 = (z1*kappa+1)^(-1/kappa);
out = norminv(z2);

end