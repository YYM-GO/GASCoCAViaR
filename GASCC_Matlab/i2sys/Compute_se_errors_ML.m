function [se_H,se_S,se_G] = Compute_se_errors_ML(f_handle,theta_opt,T)

% step 1: compute hessian
H = hessian_2sided(f_handle,theta_opt);
inv_H    = inv(H);
se_H = sqrt(diag(inv_H));

% step 2 Outer product of gradient
OPG_mat = Compute_OPG_MV_MODELS(f_handle,theta_opt,T); warning off
S_G = inv(OPG_mat)/T; warning off
se_G = sqrt(diag(S_G));

% step 3; build sandwich covariance matrix
S_COV = T * inv_H * OPG_mat * inv_H;
se_S = sqrt(diag(S_COV));
end