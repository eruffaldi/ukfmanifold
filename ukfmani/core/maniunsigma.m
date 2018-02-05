% recomposes the sigma points
%
% Input:
%  model 
%  sigma: Chiz
%  sigmainfo: the weights
%  the delta of input sigma: vChi (see manisigma)
%
% Emanuele Ruffaldi 2017 @ SSSA
function [mz,Czz,Cxz] = maniunsigma(model,Chiz,sigmainfo,vChi)

% estimates the mean in a weighted way
steps = 20;
N=size(Chiz,1);

v = zeros(size(Chiz,1),model.alg); % preallocated
mz = Chiz(1,:)'; % COL

% for lie group we make a little optimization using inv
if isfield(model,'log')
    % estimate mean but weighted of WM
    for k=1:steps
        imz = model.inv(mz);
        for i=1:N
            v(i,:) = model.log(model.prod(Chiz(i,:)',imz));
        end
        % update mz by weighted v
        mz = model.prod(model.exp(v'*sigmainfo.WM),mz);
    end
    
    % update v for computing covariance, skips the log
    imz = model.inv(mz);
    for i=1:N
        v(i,:) = model.log(model.prod(Chiz(i,:)',imz));
    end
else    
    % estimate mean but weighted of WM
    for k=1:steps
        for i=1:N
            % same as: se3_logdelta but with igk once
            v(i,:) = model.delta(Chiz(i,:),mz);
        end
        mz = model.step(mz,(v'*sigmainfo.WM)'); % [A,S] [S,1]
    end
    
    % update v for computing covariance
    for i=1:N
        v(i,:) = model.delta(Chiz(i,:),mz);
    end
end

Czz = v'*sigmainfo.W*v; % covariance ZZ

if nargin >= 4 && nargout > 2
    Cxz = vChi'*sigmainfo.W*v; % cross XZ
end