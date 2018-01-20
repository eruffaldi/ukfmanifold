% given a manifold definition and a set of manifold data (packed) computes
% the mean and the variance
%
% Inputs
%  model is the model of the manifold
%  X     is the input as a packed manifold (N x packedsize)
%
% Outputs
%  mz    is the mean in the group
%  S     is the covariance
%
% Note: for obtaining the sqrt of S use cholcov
%
% Emanuele Ruffaldi 2017
function [mz,S,nv] = manimeancov(model,X,steps)

if nargin < 3
    steps = 10;
end


% estimates the mean in a weighted way
N=size(X,1);

v = zeros(size(X,1),model.alg); % preallocated
mz = X(1,:); % first is picked as mean

% for lie group we make a little optimization using inv
if nargout > 2
    nv = zeros(steps,3);
end
if isfield(model,'log')
    % estimate mean but weighted of WM
    for k=1:steps
        imz = model.inv(mz);
        for i=1:N
            v(i,:) = model.log(model.prod(X(i,:),imz));
        end
        if nargout > 2
            p = sqrt(sum(v.^2,2));
            nv(k,:) = [norm(mean(v,1)),max(p),mean(p)];
        end
        % update mz by weighted v
        mz = model.prod(model.exp(mean(v,1)'),mz);
    end
    
    % update v for computing covariance, skips the log
    imz = model.inv(mz);
    for i=1:N
        v(i,:) = model.log(model.prod(X(i,:),imz));
    end
else    
    % estimate mean but weighted of WM
    for k=1:steps
        for i=1:N
            % same as: se3_logdelta but with igk once
            v(i,:) = model.delta(X(i,:),mz);
        end
        if nargout > 2
            p = sqrt(sum(v.^2,2)); % norm of each vector, sized as Nx1
            nv(k,:) = [norm(mean(v,1)),max(p),mean(p)];
        end
        mz = model.step(mz,mean(v,1)); % [A,S] [S,1]
    end
    
    % update v for computing covariance
    for i=1:N
        v(i,:) = model.delta(X(i,:),mz);
    end
end

S = 1/N*(v'*v); % covariance ZZ
