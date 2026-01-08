function out = skewport(weight,rate)
%calculate the skewness of the portfolio
Tout = size(weight,1);
for i = 1:Tout
    rateall(i) = weight(i,:)*rate(i,:)';
end
out = skewness(rateall');
end

