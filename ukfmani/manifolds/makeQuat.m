% SO3 as quaternion with (xyz w)
%
% Emanuele Ruffaldi 2017 @ SSSA
function m = makeQuat()


m = [];
m.type = {'Quat'};
m.inv = @(X) qconj(X);
m.prod = @(X,Y) qmult(X,Y);
m.delta = @qmdelta;
m.step = @qmstep;
m.meancov = @manimeancov; % default
m.count = 1;
m.group = 4; % as unitary quaternion
m.alg = 3;
m.pack = @(x) x;
m.unpack = @(x) x;
