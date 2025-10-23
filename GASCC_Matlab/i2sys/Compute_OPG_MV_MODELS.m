function OPG_mat = Compute_OPG_MV_MODELS(f_handle, theta_vec0,T)
% compute the expected value of the outer product of the gradient
% associated wit the log-likelihood function (given by
% f_handle)


dH = 2.5e-10*max( abs(theta_vec0), 1e-2);
k = length(theta_vec0);
g_mat = zeros(T,k);

for i = 1:k
    
    theta_vec = theta_vec0;
    
    theta_vec(i) = theta_vec0(i) + dH(i);
    [LogLike] = f_handle(theta_vec);
    L_plus = -LogLike;
    
    theta_vec(i) = theta_vec0(i) - dH(i);
    [LogLike] = f_handle(theta_vec);
    L_min = -LogLike;
    
    g_mat(:,i) = (L_plus - L_min)/ (2 * dH(i));
    
end

OPG_mat = zeros(k,k);
for j = 1:T
    OPG_mat = OPG_mat + g_mat(j,:)'*g_mat(j,:);
end

OPG_mat = OPG_mat/T;

end