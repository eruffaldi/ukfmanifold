classdef AGaussian
    % A Gaussian in a given manifold
    properties
        model % belongs to AManifold
        amean   % G
        acov    % AxA
    end
    
    methods
        function obj = AGaussian(amanifold,amean,acov)
            obj.model = amanifold;
            obj.amean = amean;
            obj.acov = acov;
        end
        
        function n = length(obj)
            n = size(obj.amean,1);
        end
        
        function m = mean(obj)
            m = obj.amean;
        end
        
        function c = cov(obj)
            c = obj.acov;
        end
        
        % compute sigma points for each Gaussian 
        % Gaussian => [q,G] Sigma points
        function [Chi,vChi] = sigmas(obj,sigmainfo)
            % decompose the covariance
            model = obj.model;
            S = obj.acov;
            mu = obj.amean;
            C = sigmainfo.sqrt(S);
            if size(C,1) ~= size(S,1)
                C = eye(size(S,1));
            end
            k = size(C,1);
            Chi = zeros(2*k+1,model.group); % flattened
            Chi(1,:) = mu; % not weighted
            vChi = zeros(2*k+1,model.alg); % delta

            c = sigmainfo.c;
            for I=1:k
                psi = c*C(I,:)'; % COLUMN
                vChi(I+1,:) = psi;
                vChi(I+1+k,:) = -psi;
                Chi(I+1,:) = model.step(mu,psi);
                Chi(I+1+k,:) = model.step(mu,-psi);
            end            
        end

        % transform all the Gaussians using f and outputs new 
        function [r,Cxy] = transform(obj,sigmainfo,f,outmodel,xsigma)
            [Chi,vChi] = obj.sigmas(sigmainfo);
            cXs = obj.model.unpack(Chi);
            cZs = cell(size(cXs,1),outmodel.count);  % [S,Cz]
            for I=1:size(cXs,1)
                [cZs{I,:}] = f(cXs{I,:});
            end
            [r,Cxz] = 
            [zm,Czz,Cxz] = maniunsigma(mz,manipack(mz,cZs),wsigmax,vXs); % [Gz,Gx]            
        end
        
        function r = prod(obj,other)
        end
        
        function r = inv(obj)
        end
        
        % returns n samples that is [n,G]
        function r = sample(obj,n)
            
        end
    end
    
    methods(Static)
        % assembles 1 Gaussian from sigma points 
        function [r,Cxz] = unsigma(model,Chi,sigmainfo,vChi,steps)
            if nargin < 5
                steps = 20;
            end
            % estimates the mean in a weighted way
            N=size(Chi,1);

            v = zeros(size(Chi,1),model.alg); % preallocated
            mz = Chi(1,:)'; % COLUMN vector

            % for lie group we make a little optimization using inv
            if isfield(model,'log')
                % estimate mean but weighted of WM
                for k=1:steps
                    imz = model.inv(mz);
                    for i=1:N
                        v(i,:) = model.log(model.prod(Chi(i,:)',imz));
                    end
                    % update mz by weighted v
                    mz = model.prod(model.exp(v'*sigmainfo.WM),mz);
                end

                % update v for computing covariance, skips the log
                imz = model.inv(mz);
                for i=1:N
                    v(i,:) = model.log(model.prod(Chi(i,:)',imz));
                end
            else    
                % estimate mean but weighted of WM
                for k=1:steps
                    for i=1:N
                    % same as: se3_logdelta but with igk once
                    v(i,:) = model.delta(Chi(i,:)',mz);
                    end
                      mz = model.step(mz,(v'*sigmainfo.WM)); % [A,S] [S,1]
                end
                % update v for computing covariance
                for i=1:N
                    v(i,:) = model.delta(Chi(i,:)',mz);
                end
            end

            Czz = v'*sigmainfo.WC*v; % covariance ZZ - NOTE THAT WE USE DIFFERENTIAL
            if nargin >= 4 && nargout > 1
                Cxz = vChi'*sigmainfo.WC*v; % cross XZ - NOTE THAT WE USE DIFFERENTIAL
            end    
            r = AGaussian(model,mz,Czz);
        end

        % given [n,G] points estimates 1 Gaussian
        function r = estimate(model,X,steps)
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
            r = AGaussian(model,mz,S);
        end

        % given n Gaussians computes the fusion 
        function r = fusion(ma,pts)
        end

    end
end

