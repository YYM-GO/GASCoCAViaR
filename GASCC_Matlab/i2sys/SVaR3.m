function out = SVaR3(x,v,q)

out = ((x<=v)-q)*(v-x);
end

