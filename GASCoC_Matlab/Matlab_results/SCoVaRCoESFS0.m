function out = SCoVaRCoESFS0(x,y,v,c,e,q)

a = (y<=c)*(y-c)/(q*e)+c/e+log(-e)-1;
out = (x<v)*a;

end
