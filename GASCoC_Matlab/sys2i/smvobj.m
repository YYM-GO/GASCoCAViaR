function output = smvobj(weight,deltai,sigma,cor)
D = diag(deltai);
G = sigma*eye(length(deltai));
output = weight'*D*G*cor*G*D*weight;

end

