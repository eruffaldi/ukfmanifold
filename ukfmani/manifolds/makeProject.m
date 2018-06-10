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
% m = makeProject(m0,R)
%   applis the linear transformation R
%
% m = makeProject(m0,R,k)
%   applies the transformation R and takes the first k tangent vectors
%
% The transformation can be buil using:
%   adjoint matrix (6x6) for SE(3) to constraint along direction
%   rotation matrix (3x3) for SO(3)
%   permutatio of axis (permVector)
% 
%
% Emanuele Ruffaldi 2017-2018 @ SSSA
function m = makeProject(m0,R,k)

n=m0.alg;
if nargin == 2
    k=n;
end
assert(k <= n && k > 0,'new space should have correct number of dofs');
assert(all(size(R)==n),'R should be n by n');
%assert R not singular
m = m0;
m.alg = k;
m.refmodels = {m0};
m.type = {'Project',m0.type,R};
m.alginc = [0,k];
ostep=m0.step;

m.islie = m0.islie;
m.unpack = @(x) projectunpack(m0,x);
m.pack = @(x) m0.pack(x);


if n == k
    m.step = @(X,v) ostep(X,(R\v(:))');
    m.delta =@(X,Y) (R*m0.delta(X,Y)')';
else
    m.step = @(X,v) ostep(X,(R\zeroextendcol(v,n,k))');
    m.delta =@(X,Y) firstk(R*m0.delta(X,Y)',k)';
end
m.transport = @(X,t,Y) xtransport(X,t,Y,m0.transport,n,k);
m.projectcov = @(C) projectcov(C,R,k);
m.unprojectcov = @(C) unprojectcov(C,R,k,n);
m.s = int_manisetup([],[],m);
    
% for Lie Group delta is: log(X*inv(Y))

if isfield(m0,'log')
    olog=m0.log;
    oexp=m0.exp;
    if n == k
        m.log = @(X) (R*olog(X))';
        m.exp = @(v) oexp((R\v(:))');
    else
        m.log = @(X) firstk((R*reshape(olog(X),[],1)),k)';
        m.exp = @(v) oexp((R\zeroextendrow(v,n,k))');
    end
end

function v = xtransport(X,t,Y,tra,n,k)

ve = tra(X,zeroextendcol(t,n,k,Y),Y);
v = ve(1:k); 

function r = projectunpack(m0,x)
r = m0.unpack(x);
if iscell(r)
    r = r{1};
end
function v = firstk(v0,k)
v = v0(1:k);

% extend as row
function v = zeroextendrow(v0,n,k)
v=zeros(1,n);
v(1:k) = v0;

% extend as col
function v = zeroextendcol(v0,n,k)
v=zeros(n,1);
v(1:k) = v0;

function CR = projectcov(C,R,k)

tmp = R*C/R;
CR = tmp(1:k,1:k);

function C = unprojectcov(CR,R,k,n)

C = zeros(n);
C(1:k,1:k) = R\CR*R;
