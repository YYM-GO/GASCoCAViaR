function CoVaR = CoVaRfrocNormalN(THETA1,THETA2,kappa,VOL2hatOS1,stdresidsIS,meanforc)

options = optimoptions('fsolve','Display','off','TolX',10^-3,'Maxiter',100);
FNor = @(v)CovarUVNormalCopula(v,kappa,THETA1,THETA2);
cavarUnorm = fsolve(FNor,0.001,options); warning off

covarstdresid = norminv(cavarUnorm);
CoVaR = covarstdresid*sqrt(VOL2hatOS1(2))+meanforc(2);



end
