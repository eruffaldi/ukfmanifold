function [mu,Sigma] = se3dinv(mu0,Sigma0)

mu= inv(mu0);
Ai = se3adj(mu);
Sigma = Ai*Sigma0*Ai';
