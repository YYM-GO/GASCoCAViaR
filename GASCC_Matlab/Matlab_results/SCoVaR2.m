function out = SCoVaR2(x,y,v,c,q)

out = (x>v)*(((y<=c)-q)*(c-y));
end