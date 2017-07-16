% Emanuele Ruffaldi 2017 @ SSSA
function [m,S] = meancov(x)

m = mean(x);
S = cov(x);