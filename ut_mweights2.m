%   n     - Dimensionality of random variable
%   k     - 
%   alpha - Transformation parameter  (optional, default 0.5)
%   beta  - Transformation parameter  (optional, default 2)
%   kappa - Transformation parameter  (optional, default 3-n)
% Copyright (C) 2006 Simo Sarkka
%
% Modified Emanele Ruffaldi 2014
function [wei] = ut_weights2(k,alpha,beta,kappa)
if nargin < 2
    alpha = 0.5;
end
if nargin < 3
    beta = 2;
end
if nargin < 4
    kappa = 3-k;
end

% mean: lambda/(k+lambda)
% others: 1/(2*(k+lambda))
%
% 2*k / (2*(k+lambda)) + lambda /(k+lambda)
% (k / (k+lambda) + lambda/ (k+lambda))
% (k+lambda) / (k+lambda) == 1

c = alpha^2*(k+kappa);
lambda = c-k; % automatic lambda

WM = repmat(1/(2*(k+lambda)), 2*k+1,1); % except first
WM(1) = lambda / (k+lambda); % first different

WC = WM;
WC(1) = WC(1) + (1-alpha^2+beta); % actually never used

W0 = eye(length(WC)) - repmat(WM,1,length(WM));
W = W0*diag(WC)*W0';

% Original: 
% W0(c) = l/(n+l) + (1-aa+b)
% Wi(c) = Wi(c)

wei = [];
wei.c = sqrt(c);
wei.WM = WM;
wei.W = W;
