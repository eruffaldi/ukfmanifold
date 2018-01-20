% quaternion constrained to single axis
%
function m = makeQuat1(axis)

axis = axis/norm(axis);

m = [];
m.axis = axis;
m.inv = @(X) qconj(X);
m.prod = @(X,Y) qmult(X,Y);
m.count = 1;
m.group = 4; % as unitary quaternion
m.alg = 1;
m.pack = @(x) x;
m.unpack = @(x) x;
m.meancov = @manimeancov; % default
m.step = @(X,y)mstep(m,X,y);
m.delta = @(X,Y) mdelta(m,X,Y);

function q = mstep(m,X,y)

if norm(y) == 0
    q = X;
else
    q = qnorm(qmult(qomega2q(m.axis*y),X));
end

function q = mdelta(m.X,Y)

w = getomega(qmult(X,qconj(Y))); % log(x*inv(y))
q = dot(w,m.axis); % project

function w = getomega(q)

[v,phi] = qdecomp(q);
assert(isnan(phi) == 0)
w = v*phi;
w = w(:)';