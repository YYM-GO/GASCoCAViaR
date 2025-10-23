function [out1,out2] = clayton_rnd(kappa,n,state)
%function out1 = clayton_rnd(kappa,n)
% Generates n (bivariate) random numbers from the Clayton copula
% with parameter kappa.
%

if nargin==1
   n = 1;	
end
if size(kappa,1)==1
   kappa = kappa*ones(n,1);		% stretching vector
end

rng(state,'philox')
U = rand(n,1);
t = rand(n,1);		% interim variable
V = (1-U.^(-kappa)+(t.*(U.^(1+kappa))).^(-kappa./(1+kappa))).^(-1./kappa);
out1 = [U,V];
if nargout==2
   out2 = out1(:,2);
   out1 = out1(:,1);
end