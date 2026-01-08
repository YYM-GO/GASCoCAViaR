function out = SCoVaRCoESN(x,y,v,c,e,q)

a = ((y<=c)-q)*c/(2*q*sqrt(-e))-((y<=c)*y/q-e)/(2*sqrt(-e))+sqrt(-e);
out = (x<v)*a;

end

