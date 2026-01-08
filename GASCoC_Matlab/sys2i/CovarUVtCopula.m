function out = CovarUVtCopula(V,RHO,NU,THETA10,THETA20,titer)



U = THETA10;


rng(100*titer,'philox')
nobs = 5000;




for tt = 1:1
    if NU(tt)>100		% then just use Normal copula
        out1(tt) = NormalCopula(U(tt),V(tt),RHO(tt));
    else
        x = tdis_inv(U(tt),NU(tt))*sqrt((NU(tt)-2)/NU(tt));  % need to adjust these as the bivartcdfmc.m is for *standardised* t random variables
        y = tdis_inv(V(tt),NU(tt))*sqrt((NU(tt)-2)/NU(tt));
        out1(tt) = bivartcdfmc(x,y,RHO(tt),NU(tt),nobs);
        
    end
    out=out1-THETA10*THETA20;
end

