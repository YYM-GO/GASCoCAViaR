function vecall = tritovec(M)
%lower triangle to vector
N = size(M,1);vecall = [];
for i = 1:(N-1)
    vec = M((i+1):N,i);
    vecall = [vecall;vec];
end
    

end

