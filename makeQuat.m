% SO3 as quaternion
function m = makeQuat()


m = [];
m.inv = @(X) qconj(X);
m.prod = @(X,Y) qmult(X,Y);
m.delta = @(X,Y) getomega(qmult(X,qconj(Y))); % log(x*inv(y))
m.step =@(X,y) qmult(X,[y(:);0]');
m.meancov = @manimeancov; % default
m.count = 1;
m.group = 4; % as unitary quaternion
m.alg = 3;
m.pack = @(x) x;
m.unpack = @(x) x;

function w = getomega(q)

w = q(1:3)';