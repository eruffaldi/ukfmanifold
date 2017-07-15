% computes sigmapoints of the manidolf using
%
% Inputs:
%  model
%  mean:mu
%  covariance: S
%  sigma weights: sigmainfo
%
% Outputs:
%  sigmapoints Chi
%  deltsa vChi
function [Chi,vChi] = manisigmas(model,mu,S,sigmainfo)

% decompose the covariance
C = cholcov(S); % TODO use svd
k = size(C,1);

Chi = zeros(2*k+1,model.group); % flattened
Chi(1,:) = mu; % not weighted
vChi = zeros(2*k+1,model.alg); % delta

c = sigmainfo.c;
for I=1:k
    psi = c*C(I,:)';
    vChi(I+1,:) = psi;
    vChi(I+1+k,:) = -psi;
	Chi(I+1,:) = model.step(mu,psi);
	Chi(I+1+k,:) = model.step(mu,-psi);
end
