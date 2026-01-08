function CoVaRforc = CoVaRfrocClaytonN(THETA1,THETA2,kappa,VOL2hatOS,stdresidsIS,meanforc)
%
covarstdresidual = CovarUVClaytonN(THETA1,THETA2,kappa,stdresidsIS);

CoVaRforc = covarstdresidual*sqrt(VOL2hatOS(2))+meanforc(2);

end