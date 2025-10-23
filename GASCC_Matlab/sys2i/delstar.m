function out = delstar(weight,ksen,del,G,cor)
D = diag(del);
out = 0.5*weight'*D*G*cor*G*D*weight/sqrt(2*weight'*D*G*cor*G*D*weight+(weight'*D*ksen)^2);
end

