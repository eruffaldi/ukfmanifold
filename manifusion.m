% given two estimates this operation performs the Covariance weighted
% fusion of the two estimates
%
% Emanuele Ruffaldi 2017
%
% TODO: implement the multiple fusion algorithms == weighted mean 
function [X,C] = manifusion(m,x0,x1,C0,C1)

% inv(inv(C0)+inv(C1))
C = C0 - C0/(C0+C1)*C0;

v = m.delta(x1,x0); % from x0 to x1

X = m.step(x0,C/C1*v(:)); % idem but weighted