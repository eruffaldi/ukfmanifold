%
% Kalman Smoorther using the following convention
%
% x(t) is the state variable
% x(t+) f(x(t),A1,..An) is the process function that accepts x(t) and other
% inputs, eventually noise
% where Ai are time invariant variables with their distribution
%
% manifold mxa = mx + ma 
%
% mxa   = joint manifold (product manifold mandatory indexable)
% x_manis = indices of manifold (typically the first k elements)
% Xest  = estimated 
%   mean T,G_X
%   cov  T,a_X,a_X
% Xpred = prediction (optional) or built using f
%   mean T,G_X
%   cov  T,a_X,a_X
% A     = constant associated to the manifold spec 
%   mean 1,G_A
%   cov  1,a_A,a_A
%
% Emanuele Ruffaldi 2018
function [A_mean,A_cov,x0_mean,x0_cov] = manismoother(mxa, x_manis,f, Xpred_mean, Xpred_cov, Xest_mean, Xest_cov, A_mean, A_cov)

% compute Xpred if needed and store it
% loop backward 
%   compute intial value: x0_mean,x0_cov that are based 
%   accumulate contributions for A 
% use last value to update x0
% average the contributions for A 