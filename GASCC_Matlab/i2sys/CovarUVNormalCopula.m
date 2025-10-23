% Computes the value of the Normal copula at a specified point
% 
% INPUTS:	U,   a Tx1 vector (or a scalar) of F(X[t])
%				V,   a Tx1 vector (or a scalat) of G(Y[t])
%				RHO, a Tx1 vector (or a scalar) of correlation coefficients
%

function out1 = CovarUVNormalCopula(V,RHO,THETA1,THETA2)


U=THETA1;


X = norminv(U,0,1);
Y = norminv(V,0,1);

out1 = mvncdf([X,Y],[0,0],[1,RHO;RHO,1])-THETA1*THETA2;
