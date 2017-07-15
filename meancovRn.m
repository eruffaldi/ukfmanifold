function [m,S] = meancov(x)

m = mean(x);
S = cov(x);