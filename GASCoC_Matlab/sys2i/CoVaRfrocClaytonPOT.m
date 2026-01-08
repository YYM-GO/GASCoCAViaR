function CoVaRforc = CoVaRfrocClaytonPOT(THETA1,THETA2,kappa,VOL2hatOS,stdresidsIS,meanforc,uL,beta,xi)

covarstdresidual = CovarUVClaytonPOT(THETA1,THETA2,kappa,stdresidsIS,uL,beta,xi);

CoVaRforc = covarstdresidual*sqrt(VOL2hatOS(2))+meanforc(2);
end