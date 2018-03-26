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
            if nargin == 1
                obj.amean = zeros(1,obj.model.group);
                obj.acov = zeros(obj.model.alg);
            else
                assert(all(size(amean) == [1,obj.model.group]),'mismatch mean group size in AGaussian');
                assert(all(size(acov) == [obj.model.alg,obj.model.alg]),'mismatch cov alg size in AGaussian');
                obj.amean = amean;
                obj.acov = acov;
            end
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
            if nargin == 1
                sigmainfo = obj.model.wsigma;
            end
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
                psi = c*C(I,:); % ROW as ALL ALG
                vChi(I+1,:) = psi;
                vChi(I+1+k,:) = -psi;
                Chi(I+1,:) = model.step(mu,psi);
                Chi(I+1+k,:) = model.step(mu,-psi);
            end            
        end
        
        % given the sigma points asscoated to the current object compute
        % the relative vector
        function vChi = vchifromchi(obj,Chi)
            model = obj.model;
            k = size(C,1);
            vChi = zeros(2*k+1,model.alg); % delta
            for I=1:k
                vChi(I+1,:) = model.delta(obj.mu,Chi(I+1,:));
                vChi(I+1+k,:) = model.delta(obj.mu,Chi(I+1+k,:));
            end            

        end
        
        % marginalize the current Gaussian using the provided target
        % manifold and the indices that extract data from mean (group) and
        % covariance (algebra)
        function x = marginalize(obj,mx,indices,aindices)
            x = AGaussian(mx,obj.amean(indices),obj.acov(aindices,aindices));
        end
        
        % given a multivariate manifold gaussian (xy) that is this object
        % 
        % assigns the mean and covariance of the partition defined by (indices,aindices) using the Gaussian x 
        function xy = assign(xy,x,indices,aindices)
            xy.amean(indices) = mean(x);
            xy.acov(aindices,aindices) = cov(x);
        end
        
        % applies the Kalman correction:
        %   obj is the predicted state
        %   z   is the estimated observation
        %   Cxz is the cross variance of xp and z from the h(x) computation
        %   zvalur is the value of the observation
        %
        % note that if z has an additive noise R it has to be applied
        % BEFORE this condition
        function [xc,deltax] = condition(obj,z,Cxz,zvalue)
            Czz = cov(z);
            Cxx_p = cov(obj);
            K = Cxz/Czz;
            Cxx_c = (Cxx_p - K * Czz * K');
            dz = z.model.delta(zvalue,mean(z));
            deltax = (K*dz')';
            mx_c = obj.model.step(mean(obj),deltax);
            
            xc = AGaussian(obj.model,mx_c,Cxx_c);
        end

        % transform all Gaussians using f and outputs new model
        % if model is not provided we assume that is the current model
        function [z,Cxz,Chiz] = transform(obj,f,outmodel)
            if nargin < 3
                outmodel = obj.model;
            end
            sigmainfo = obj.model.wsigma;
            [Chi,vChi] = obj.sigmas(sigmainfo); % Chi is [sigmas,packed group]
            assert(size(Chi,1) == size(vChi,1),'same group and vchi');
            assert(size(Chi,2) == obj.model.group,'correct model group size');
            assert(size(vChi,2) == obj.model.alg,'correct model alg size');
            cXs = obj.model.unpack(Chi); % [sigmas, unpacked group]
            assert(size(cXs,1) == size(Chi,1),'same rows');
            cZs = cell(size(cXs,1),outmodel.count);  % [unpacked z group,sigmas] cell
            for I=1:size(cXs,1) % for every sigma point
                [cZs{I,:}] = f(cXs{I,:}); 
            end
            % manipack wants [sigmas, unpacked z group]
            Chiz = manipack(outmodel,cZs);
            [zm,Czz,Cxz] = maniunsigma(outmodel,Chiz,sigmainfo,vChi); % [Gz,Gx]   
            z = AGaussian(outmodel,zm,Czz);
        end

        
        
        % transform all Gaussians using f and outputs new model
        % if model is not provided we assume that is the current model
        function [z,Cxz,Chiz] = transformSigma(obj,Chix,f,outmodel)
            if nargin < 4
                outmodel = obj.model;
            end
            sigmainfo = obj.model.wsigma;
            cXs = obj.model.unpack(Chix); % cXs is [unpacked group,sigmas] cell
            cZs = cell(size(cXs,2),outmodel.count);  % [unpacked z group,sigmas] cell
            for I=1:size(cXs,1) % for every sigma point
                [cZs{I,:}] = f(cXs{I,:}); 
            end
            % manipack wants [sigmas, unpacked z group]
            Chiz = manipack(outmodel,cZs');
            vChi = obj.vchifromchi(Chix);
            [zm,Czz,Cxz] = maniunsigma(outmodel,Chiz,sigmainfo,vChi); % [Gz,Gx]   
            z = AGaussian(outmodel,zm,Czz);
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

