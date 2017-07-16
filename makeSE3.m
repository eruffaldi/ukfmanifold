% SE3 to be tested, alternatively use makeCom(makeRot(),makeRn(3))
%
% Emanuele Ruffaldi 2017 @ SSSA
function m = makeSE3()

m1 = makeRot();

m = [];
m.inv = @(x) uflatten(inv(unflat(x));
m.prod = @(x,y) uflatten(unflat(x)*unflat(y));
m.log = @(x) [m1.log(getrot(x)); getpos(x)]; % not exact
m.exp = @(x) ubuild2(m1.exp(getomega(x)),x(4:6)); % not exact
m.delta = @(x,y) [m1.delta(getrot(x)); x(4:6)-y(4:6)];
m.group = 16;
m.alg = 6;
m.count = 1;
m.meancov = @manimeancov;
m.pack = @(x) x(:);
m.unpack = @(x) reshape(x,4,4);

function u = uflatten(x)

u = x(:);

function u = unflat(x)

u = reshape(x,4,4);

function u = ubuild2(rot,pos)

u = eye(4);
u(1:3,1:3) = rot;
u(1:3,4) = pos;

function u = getrot(x)

u = x(1:3,1:3);

function u = getpos(x)

u = x(1:3,4);

function u = getomega(x)

u = x(1:3);

function u = getvel(x)

u = x(4:6);