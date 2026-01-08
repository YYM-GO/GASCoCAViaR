function CoVaR = CoVaRfroctcopN(THETA1,THETA20,kappa1,kappa2,VOL2hatOS1,stdresidsIS,meanforc,titer)

options1 = optimset('TolX',1e-5,'Display', 'off','Maxfuneval',50,'MaxTime',10);
Ftcop = @(v)CovarUVtCopula(v,kappa1,kappa2,THETA1,THETA20,titer);
cavarUt = fzero(Ftcop,[0.000001,0.999999],options1);
covarstdresid = norminv(cavarUt);
CoVaR = covarstdresid*sqrt(VOL2hatOS1(2))+meanforc(2);

end

