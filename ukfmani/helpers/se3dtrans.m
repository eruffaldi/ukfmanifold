function [mu,Sigma] = se3dtrans(T,mu0,Sigma0)

mu = T*mu0;
AT = se3adj(T);
Sigma = AT*Sigma0*AT';