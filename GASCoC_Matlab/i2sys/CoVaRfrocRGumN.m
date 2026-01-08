function CoVaRforc = CoVaRfrocRGumN(THETA1,THETA2,kappa,VOL2hatOS,stdresidsIS,meanforc)

covarstdresidual = CovarUVRGumN(THETA1,THETA2,kappa,stdresidsIS);
     CoVaRforc = covarstdresidual*sqrt(VOL2hatOS(2))+meanforc(2);

end