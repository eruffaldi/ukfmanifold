% given X in SE3 unpacked form we would like to compute mean and covariance
% using manopt
%
% sympositivedefinitefactory(6,1)
% specialeuclideanfactory(3, 1)
function [m,S] = meancovSE3manopt(X)

%bishop page 93
%loglikelihood
% cost + -N/2 log det(Sigma) - sum (xi-mu)' inv(Sigma) (xi-mu)
%
% derivative for mu
%   sum inv(sigma) (x-mu)
%   => mu = 1/N sum (x)
%
% derivative for Sigma
%   -N/2 inv(A)' - ... TBC
%   => Sigma = 1/N sum (x-mu)' (x-mu)
%
%
% The maximization of (2.118) with respect to Σ is rather more involved. The simplest approach is to ignore the symmetry constraint and show that the resulting solution is symmetric as required. Alternative derivations of this result, which impose the symmetry and positive defi- niteness constraints explicitly, can be found in Magnus and Neudecker (1999). 
% es. 2.34 Using the results (C.21), (C.26), and (C.28) from Appendix C, show that the covariance matrix Σ that maximizes the log likelihood function (2.118) is given by the sample covariance (2.122). We note that the final result is necessarily symmetric and positive definite (provided the sample covariance is nonsingular).
%
% using: 
%   d/dA ln(|A|) = (inv(A))'
%   d/dA tr(A) = I
%   d/dx inv(F) = - inv(F) dF/dx inv(F)
%   d/dA tr(A B A') = A (B + B')
%
% Manifold loglikelihood
% cost + -N/2 log det(Sigma) - sum delta(xi,mu)' inv(Sigma) delta(xi-mu)
%
% delta(m1,m2) = log(m1*inv(m2)) if the manifold is Lie Group
%
% http://www.janmagnus.nl/misc/mdc2007-3rdedition
% https://docs.google.com/document/d/1QlzmQe1U1vEwN2ShqPTB3V2uG1tW0snGEosWZgbOb4Q/edit



problem1.M = specialeuclideanfactory(3, 1);

% zi = data
% x = mu
% sum_i (zi-x)
problem1.cost  = @(x) -x'*(A*x);
problem1.egrad = @(x) -2*A*x;
 
% Numerically check gradient consistency (optional).
checkgradient(problem1);
 
% Solve.
[x, xcost, info, options] = trustregions(problem1);



problem2.M = sympositivedefinitefactory(6, 1);

% zi = data
% x = mu

% 1/N sum_i (z-x)
problem2.cost  = @(x) -x'*(A*x);
problem2.egrad = @(x) -2*A*x;
 
% Numerically check gradient consistency (optional).
checkgradient(problem1);
 
% Solve.
[x, xcost, info, options] = trustregions(problem1);