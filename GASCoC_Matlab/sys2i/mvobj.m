function output = mvobj(weight,sigma,cor)

G = diag(sqrt(sigma));
output = weight'*G*cor*G*weight;

end

