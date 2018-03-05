% Takes a manifold and produces a new one whose tangential space is a
% transformation of the original one, potentially projected
%
% This meta manifold works with a projection of the tangential space (e.g.
% single axis rotation). The projection is expressed as follows:
%   R matrix that rotates the tangential space
%   k output degrees of freedom (the first k over the overall incoming n
%   degrees)
%
% Note that the state is still in the original frame
%
% 
%
% Emanuele Ruffaldi 2017-2018 @ SSSA
function m = makeProject(m0,R,k)

n=m0.alg;
assert(k >= n && k > 0,'new space should have correct number of dofs');
assert(all(size(R)==n),'R should be n by n');
%assert R not singular
m = m0;
m.alg = k;
m.models = {m0};
m.type = {'Project',m0.type,R};
m.alginc = [0,k];
ostep=m0.step;
m.step = @(X,v) ostep(X,R\zeroextend(v,n,k));
m.delta =@(X,Y) firstk(R*m0.delta(X,Y),k);

% for Lie Group delta is: log(X*inv(Y))

if isfield(m0,'log')
    olog=m0.log;
oexp=m0.exp;

    m.log = @(X) firstk(R*olog(X),k);
    m.exp = @(v) oexp(R\zeroextend(v,n,k));
end

function v = firstk(v0,k)
v = v0(1:k);

function v = zeroextend(v0,n,k)
v=zeros(n,1);
v(1:k) = v0;
