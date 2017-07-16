% sqrt using SVD
%
% Originally introduced in the demonstration https://github.com/eruffaldi/compare-mvn-transform
%
% Emanuele Ruffaldi 2016-2017
function Y = svdsqrt(X)

[~,S,V] = svd(X);
R = V;
S = sqrt(S);
Y = R*S;
%disp('svd')
%Y*Y'-X
%disp('cholcov')
%C=cholcov(X)';
%if isempty(C) == 0
 %   C*C'-X
%end

