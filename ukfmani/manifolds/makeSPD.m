% Symmetric Positive Definite
%
% Could be used for modeling covariance of Gaussians, but NOTE that a
% typical prior of MVN Gaussians is the Normal-Wishart Distribution defined
% by:
% - mu0
% - lambda
% - W
% - n dofs
%
% 
% 
% Emanuele Ruffaldi 2018 @ SSSA
function m = makeSPD(n)


m = [];
m.type = {'SPD',n};
m.delta = @spd_delta;
m.step = @spd_step;
m.transport = @spd_transport;
m.count = 1;
m.group = n*n; % as matrix
m.alg = n*(n+1)/2;
m.pack = @(x) x(:)';
m.unpack = @(x) reshape(x,n,n);
m.vpack = @vpack;
m.vunpack = @(v) vunpack(v,n);

function Y = symm(X)
    Y = 0.5*(X+X');
end

function tp = vpack(t)
    a = triu(t);
    tp = a(:)';
end

function t = vunpack(tp,n)
    t = zeros(n,n);  % constant
    u = triu(true(n,n)); % constant
    t(u(:)) = tp;
end

function tpy = spd_transport(X,tpx,Y)

% from manopt
tx = vunpack(tpx);
E = sqrtm((Y/X));
tpy = vpack(E*tx*E');

end

function t = spd_delta(X,Y)
t = vpack(symm(X*real(logm(X\Y))));
end

function Y = spd_step(X,tp)    
t = vunpack(tp);
Y = symm(X*real(expm(X\(t*eta))));

end
 
end