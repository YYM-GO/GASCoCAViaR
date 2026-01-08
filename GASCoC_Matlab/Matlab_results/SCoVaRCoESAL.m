function out = SCoVaRCoESAL(x,y,v,c,e,q)

 out = (x<v)*(-log((q-1)/e)-(y-c)*(q-(y<=c))/(q*e)+y/e);

end

