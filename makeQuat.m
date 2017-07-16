% SO3 as quaternion with (xyz w)
%
% Emanuele Ruffaldi 2017 @ SSSA
function m = makeQuat()


m = [];
m.inv = @(X) qconj(X);
m.prod = @(X,Y) qmult(X,Y);
m.delta = @mdelta;
m.step = @mstep;
m.meancov = @manimeancov; % default
m.count = 1;
m.group = 4; % as unitary quaternion
m.alg = 3;
m.pack = @(x) x;
m.unpack = @(x) x;

function q = mstep(X,y)

if norm(y) == 0
    q = X;
else
    q = qnorm(qmult(qomega2q(y),X));
end

function w = mdelta(X,Y)

w = getomega(qmult(X,qconj(Y))); % log(x*inv(y))

function w = getomega(q)

[v,phi] = qdecomp(q);
assert(isnan(phi) == 0)
w = v*phi;
w = w(:)';